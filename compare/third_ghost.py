"""Third independent recomputation (exact sympy; imports nothing from engine_a/engine_b/classical code).
Inputs: the action, Q^GF propagator components and Fourier convention of classical/classical_spec.md sections 4,7,8;
the one-loop formula W_1 = hbar sum_n 1/(2n) str[(P (-1)^{|.|} d2I)^n] verified independently in toy3.py.
g = su(2), orthonormal basis e_a with [e_a,e_b] = eps_abc e_c and (e_a,e_b) = delta_ab (code only).

Computes, as coefficients of I^CT[eps] in units u = hbar c_adj/(16 pi^2):
 (1) graph 'AV-AV ghost loop'   -> B1, B2, B3     (engines: 1/24, 1/8, 1/16)
 (2) graph 'A*-c, lines A_V and c-A*' -> B7       (engines: -1)
 (3) the one-loop log coefficient of the ghost vertex int ([A,c],A^*) (degree 0), all graphs, zero external momenta.
"""
import itertools, json, sys
import sympy as sp
from grassmann3 import G

I_ = sp.I
k = sp.symbols('k0:4')
p = sp.symbols('p0:4')
ls = sp.symbols('l0:4')
uu, T = sp.symbols('u T', positive=True)
DIMG = 3
EPS = lambda a, b, c: sp.LeviCivita(a, b, c)

# ---------------- colour scalar c_adj: X -> sum_a [t_a,[t^a,X]] ----------------
def bracket(X, Y):
    return [sum(EPS(a, b, c) * X[a] * Y[b] for a in range(3) for b in range(3)) for c in range(3)]
Xs = sp.symbols('X0:3')
E = [[1 if i == j else 0 for i in range(3)] for j in range(3)]
cas = [0, 0, 0]
for a in range(3):
    inner = bracket(E[a], list(Xs))
    outer = bracket(E[a], inner)
    cas = [cas[i] + outer[i] for i in range(3)]
C_ADJ = sp.simplify(cas[0] / Xs[0])
assert all(sp.simplify(cas[i] - C_ADJ * Xs[i]) == 0 for i in range(3))

# ---------------- fluctuation index set: A_V (mu,a) even, c (a) odd, A*_V (mu,a) odd ----------------
COMPS = [('A', m, a) for a in range(3) for m in range(4)] + [('c', None, a) for a in range(3)] + \
        [('As', m, a) for a in range(3) for m in range(4)]
IDX = {c: i for i, c in enumerate(COMPS)}
PAR = [0 if c[0] == 'A' else 1 for c in COMPS]
N = len(COMPS)
# odd generator ids: x-modes 0..14, y-modes 15..29, background r_c = 30, r_* = 31
ODD_X = {c: i for i, c in enumerate([c for c in COMPS if c[0] != 'A'])}
ODD_Y = {c: 15 + i for i, c in enumerate([c for c in COMPS if c[0] != 'A'])}
RC, RS = 30, 31
XA = {c: sp.Symbol('xA_%d_%d' % (c[1], c[2])) for c in COMPS if c[0] == 'A'}
YA = {c: sp.Symbol('yA_%d_%d' % (c[1], c[2])) for c in COMPS if c[0] == 'A'}

def xvar(c):
    return G.c(XA[c]) if c[0] == 'A' else G.gen(ODD_X[c])
def yvar(c):
    return G.c(YA[c]) if c[0] == 'A' else G.gen(ODD_Y[c])
def dx(F, c):
    return F.deven(XA[c]) if c[0] == 'A' else F.dodd(ODD_X[c])
def dy(F, c):
    return F.deven(YA[c]) if c[0] == 'A' else F.dodd(ODD_Y[c])

# ---------------- cubic density terms (classical_spec section 8), slot order as printed ----------------
# term = (coef, [(field, mu, a, deriv_index_or_None)]*3)
TERMS = []
for m in range(4):
    for n in range(m + 1, 4):
        for a, b, c in itertools.product(range(3), repeat=3):
            e = EPS(b, c, a)
            if e == 0:
                continue
            # -1/2 (d_m A_n - d_n A_m , [A_m, A_n]) ; [A_m,A_n]^a = eps_bca A_m^b A_n^c
            TERMS.append((-sp.Rational(1, 2) * e, [('A', n, a, m), ('A', m, b, None), ('A', n, c, None)]))
            TERMS.append((sp.Rational(1, 2) * e, [('A', m, a, n), ('A', m, b, None), ('A', n, c, None)]))
for m in range(4):
    for a, b, d in itertools.product(range(3), repeat=3):
        e = EPS(a, b, d)
        if e != 0:   # ([A_m, c], A*_m) = eps_abd A_m^a c^b A*_m^d
            TERMS.append((e, [('A', m, a, None), ('c', None, b, None), ('As', m, d, None)]))

def hessian(bg, pbg, qL, qR):
    """(d2I)_{jl} = d_{x_j}( d_{y_l} I ), x-modes at incoming momentum qL, y-modes at qR, one background leg.
    bg: dict {(field, mu, a): G} background polarisation (with r for odd fields) at incoming momentum pbg."""
    H = {}
    mom = {'bg': pbg, 'x': qL, 'y': qR}
    for coef, slots in TERMS:
        for kinds in itertools.permutations(('bg', 'x', 'y')):
            fac = coef
            elems = []
            ok = True
            for (f, m, a, der), kd in zip(slots, kinds):
                comp = (f, m, a)
                if kd == 'bg':
                    if comp not in bg:
                        ok = False; break
                    elems.append(bg[comp])
                elif kd == 'x':
                    elems.append(xvar(comp)); cx = comp
                else:
                    elems.append(yvar(comp)); cy = comp
                if der is not None:
                    fac = fac * I_ * mom[kd][der]
            if not ok:
                continue
            prod = elems[0] * elems[1] * elems[2]
            val = dx(dy(prod, cy), cx) * fac
            if not val.iszero():
                key = (IDX[cx], IDX[cy])
                H[key] = H.get(key, G()) + val
    return H

def prop(kv, mask):
    """numerator of P(eps,infty)(k) = Num(k) e^{-eps k^2}/k^2, component X(x)Y = d_P(phi^X(k) phi^Y(-k))."""
    P = {}
    for a in range(3):
        if 'AA' in mask:
            for m in range(4):
                P[(IDX[('A', m, a)], IDX[('A', m, a)])] = sp.Integer(2)
        if 'gh' in mask:
            for m in range(4):
                P[(IDX[('c', None, a)], IDX[('As', m, a)])] = -I_ * kv[m]
                P[(IDX[('As', m, a)], IDX[('c', None, a)])] = -I_ * kv[m]
    return P

def PPiH(P, H):
    M = {}
    byrow = {}
    for (i, j), v in P.items():
        byrow.setdefault(j, []).append((i, v))
    for (j, l), h in H.items():
        for (i, v) in byrow.get(j, []):
            M[(i, l)] = M.get((i, l), G()) + h * (v * (-1) ** PAR[j])
    return M

def mul(A, B):
    out = {}
    brow = {}
    for (k_, l), v in B.items():
        brow.setdefault(k_, []).append((l, v))
    for (i, k_), v in A.items():
        for (l, w) in brow.get(k_, []):
            out[(i, l)] = out.get((i, l), G()) + v * w
    return out

def strace(M):
    tot = G()
    for (i, j), v in M.items():
        if i == j:
            tot = tot + v * (-1) ** PAR[i]
    return tot

# ---------------- loop-momentum integration (Schwinger, exact eps-singular parts) ----------------
def gauss_avg(expr):
    """<expr(l)> for weight e^{-T l^2}: <l_m^{2r}> = (2r-1)!!/(2T)^r, componentwise."""
    poly = sp.Poly(sp.expand(expr), *ls)
    tot = 0
    for mon, cf in poly.terms():
        val = cf
        for e in mon:
            if e % 2:
                val = 0; break
            val = val * sp.factorial2(e - 1) / (2 * T) ** (e // 2)
        tot += val
    return sp.expand(tot)

def bubble_sing(Num, Delta):
    """Num(k) for lines k (t1) and k - Delta (t2); returns (pole, log) with
    int d^4k/(2pi)^4 int_{[eps,inf)^2} e^{-t1 k^2 - t2 (k-Delta)^2} Num = (pole/eps + log*log eps + ...)/(16 pi^2)."""
    sub = {k[m]: ls[m] + uu * Delta[m] for m in range(4)}
    avg = gauss_avg(Num.subs(sub, simultaneous=True))          # sum_m C_m(u) T^{-m}
    w = uu * (1 - uu) * sum(x * x for x in Delta)
    expo = sum((-T * w) ** j / sp.factorial(j) for j in range(4))
    full = sp.expand(avg * expo)                                 # times T^{-2}/(16pi^2)
    Tm = sp.Poly(sp.expand(full * T**6), T)                      # powers shifted by 6
    logpart = 0; polepart = 0
    for (e,), cf in Tm.terms():
        power = e - 6 - 2        # integrand ~ T^power
        if power == -2:
            logpart += -sp.integrate(cf, (uu, 0, 1))
        elif power == -3:
            polepart += sp.integrate((1 - uu) * cf, (uu, sp.Rational(1, 2), 1)) + sp.integrate(uu * cf, (uu, 0, sp.Rational(1, 2)))
        elif power < -3:
            raise ValueError('unexpected power %d' % power)
    return sp.expand(polepart), sp.expand(logpart)

def triangle_log_zero_momentum(Num):
    """all three lines at loop momentum k; Num homogeneous of degree 2 -> T^{-3}; int_{[eps,inf)^3} T^{-3} = -(1/2) log eps."""
    avg = gauss_avg(Num.subs({k[m]: ls[m] for m in range(4)}, simultaneous=True))
    cf = sp.expand(avg * T)
    assert not cf.has(T), cf
    return sp.expand(-sp.Rational(1, 2) * cf)

def coeff_key(Gel, key):
    return Gel.coeff(key)

results = {'c_adj_su2': str(C_ADJ)}
U = C_ADJ   # result in units hbar c_adj/(16 pi^2): divide the 1/(16pi^2)-coefficient by c_adj

# ======== (1) AV-AV ghost loop ========
uA = {('A', m, a): G.c(sp.Symbol('u_%d_%d' % (m, a))) for m in range(4) for a in range(3)}
vA = {('A', m, a): G.c(sp.Symbol('v_%d_%d' % (m, a))) for m in range(4) for a in range(3)}
mp = [-x for x in p]
kmp = [k[m] - p[m] for m in range(4)]
kpp = [k[m] + p[m] for m in range(4)]
def two_point(bgL, pL, bgR, pR, mask, key):
    """(1/4)[str(P(k) Pi H_L P(k-pL) Pi H_R) + str(P(k) Pi H_R P(k-pR) Pi H_L)]; returns (pole, log)/(16pi^2)."""
    pole = 0; lg = 0
    for (b1, p1, b2, p2) in ((bgL, pL, bgR, pR), (bgR, pR, bgL, pL)):
        k1 = list(k); k2 = [k[m] - p1[m] for m in range(4)]
        H1 = hessian(b1, p1, [-x for x in k1], k2)
        H2 = hessian(b2, p2, [-x for x in k2], k1)
        M = mul(PPiH(prop(k1, mask), H1), PPiH(prop(k2, mask), H2))
        num = coeff_key(strace(M), key)
        a_, b_ = bubble_sing(sp.expand(num * sp.Rational(1, 4)), p1)
        pole += a_; lg += b_
    return sp.expand(pole), sp.expand(lg)

pole1, log1 = two_point(uA, list(p), vA, mp, {'gh'}, ())
# basis cross terms on A = u e^{ipx} + v e^{-ipx}
def cross_B(name):
    tot = 0
    for (x, px, y, py) in ((uA, p, vA, mp), (vA, mp, uA, p)):
        for a in range(3):
            if name == 'B1':
                for m in range(4):
                    for n in range(m + 1, 4):
                        Fx = I_ * px[m] * x[('A', n, a)].coeff(()) - I_ * px[n] * x[('A', m, a)].coeff(())
                        Fy = I_ * py[m] * y[('A', n, a)].coeff(()) - I_ * py[n] * y[('A', m, a)].coeff(())
                        tot += Fx * Fy
            elif name == 'B2':
                dx_ = sum(I_ * px[m] * x[('A', m, a)].coeff(()) for m in range(4))
                dy_ = sum(I_ * py[m] * y[('A', m, a)].coeff(()) for m in range(4))
                tot += dx_ * dy_
            elif name == 'B3':
                tot += sum(x[('A', m, a)].coeff(()) * y[('A', m, a)].coeff(()) for m in range(4))
    return sp.expand(tot)
# cross term of a quadratic functional F on u+v is F(u,v)+F(v,u): computed above as the sum over both assignments
def fit(target, names):
    cs = sp.symbols('C0:%d' % len(names))
    expr = sp.expand(target - sum(c * cross_B(nm) for c, nm in zip(cs, names)))
    syms = sorted(expr.free_symbols - set(cs), key=str)
    eqs = sp.Poly(expr, *syms).coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    assert sol, 'not in span'
    return {nm: sp.nsimplify(sol[0][c] / U) for c, nm in zip(cs, names)}
r1 = {'eps_inv': fit(pole1, ['B1', 'B2', 'B3']), 'log_eps': fit(log1, ['B1', 'B2', 'B3'])}
results['AV-AV ghost loop'] = {kk: {a: str(b) for a, b in vv.items()} for kk, vv in r1.items()}
print('AV-AV ghost loop', results['AV-AV ghost loop'], flush=True)

# ======== (2) A*-c two-point, lines A_V and c-A* ========
gam = [sp.Symbol('g_%d' % a) for a in range(3)]
alp = {(m, a): sp.Symbol('al_%d_%d' % (m, a)) for m in range(4) for a in range(3)}
bg_c = {('c', None, a): G.gen(RC, gam[a]) for a in range(3)}
bg_s = {('As', m, a): G.gen(RS, alp[(m, a)]) for m in range(4) for a in range(3)}
pole2, log2 = two_point(bg_c, list(p), bg_s, mp, {'AA', 'gh'}, (RC, RS))
# tree (dc, A*) on c = r_c gam e^{ipx}, A* = r_* al e^{-ipx}:  sum_m (i p_m r_c gam^a)(r_* al^{m a})
tree7 = G()
for m in range(4):
    for a in range(3):
        tree7 = tree7 + G.gen(RC, I_ * p[m] * gam[a]) * G.gen(RS, alp[(m, a)])
t7 = tree7.coeff((RC, RS))
b7 = sp.simplify(log2 / t7 / U)
assert sp.simplify(pole2) == 0
results['As-c graph (A_V, c-A*)'] = {'log_eps B7': str(b7), 'eps_inv': '0'}
print('B7', b7, flush=True)

# ======== (3) ghost vertex ([A,c],A*) at zero external momenta, all one-loop graphs ========
avec = {('A', m, a): G.c(sp.Symbol('a_%d_%d' % (m, a))) for m in range(4) for a in range(3)}
zero = [0, 0, 0, 0]
Hs = {'a': hessian(avec, zero, [-x for x in k], list(k)),
      'g': hessian(bg_c, zero, [-x for x in k], list(k)),
      's': hessian(bg_s, zero, [-x for x in k], list(k))}
Pk = prop(list(k), {'AA', 'gh'})
Ms = {n: PPiH(Pk, h) for n, h in Hs.items()}
tri = G()
for order in itertools.permutations('ags'):
    M = mul(mul(Ms[order[0]], Ms[order[1]]), Ms[order[2]])
    tri = tri + strace(M)
num3 = sp.expand(tri.coeff((RC, RS)) * sp.Rational(1, 6))      # 1/(2n), n = 3; all 6 orderings
log3 = triangle_log_zero_momentum(num3)
tree3 = G()
for m in range(4):
    for a, b, d in itertools.product(range(3), repeat=3):
        e = EPS(a, b, d)
        if e:
            tree3 = tree3 + avec[('A', m, a)] * G.gen(RC, gam[b]) * G.gen(RS, alp[(m, d)]) * e
t3 = tree3.coeff((RC, RS))
ratio3 = sp.simplify(log3 / t3 / U)
# per-graph split: graph (i) two ghost lines + one A line; graph (ii) one ghost line + two A lines
results['ghost vertex ([A,c],A*)'] = {'log_eps / tree': str(ratio3)}
print('ghost vertex log coefficient / tree', ratio3, flush=True)

# β-direction bookkeeping (see report): r(monomial) = beta (1 - #antifields) + sum of Euler weights
beta, wA, wc = sp.symbols('beta w_A w_c')
sol = sp.solve([beta + 2 * wA - 4, wc - wA - (-1), wc - ratio3], [beta, wA, wc], dict=True)[0]
results['coupling direction'] = {k_.name: str(v) for k_, v in sol.items()}
print('coupling-direction solve', sol)
json.dump(results, open('third_ghost.json', 'w'), indent=1)
