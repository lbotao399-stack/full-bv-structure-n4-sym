"""Independent finite-dimensional audit model of the Feynman-type Q^GF (DESIGN 'Gauge fixing' section).
One 'momentum mode' analogue of the minimal N=4 complex: c (1), A (3: two 'V' components + one 'W-perp' component),
A^* (3), c^* (1), lambda (2), lambda^* (2).  'd' = D (column, D^T D = h), 'd^*d' on A = M = h - D D^T on the V block, h on
the W-perp component, 'dslash' = G (antisymmetric, G^2 = -h).  Conventions copied only from DESIGN / ledger / classical_spec:
pairs (x even, xi odd) = (A,A^*), (lambda^*,lambda), (c^*,c), <x,xi> = +1, Delta_0 = sum d_x d_xi, {f,g} = Delta_0(fg) -
Delta_0(f) g - (-1)^{|f|} f Delta_0(g), K_0 = -sum (x(x)xi + xi(x)x), partial_phi = 1/2 sum phi^{ij} d_j d_i,
S_free = -1/4 (dA,dA) + 1/2 (lambda, dslash lambda) + (dc, A^*), Q := {S_free, -}.
Everything exact (sympy).  Output: toy_gf.json."""
import json, random
import sympy as sp
from superpoly import Alg

h = sp.Integer(25)
D = [sp.Integer(3), sp.Integer(4)]          # D^T D = 25 = h
G = sp.Matrix([[0, 5], [-5, 0]])            # G^2 = -25 = -h, antisymmetric
eps, hb = sp.symbols('epsilon hbar', positive=True)

# generator indices
names = ['c', 'A0', 'A1', 'A2', 'As0', 'As1', 'As2', 'cs', 'l0', 'l1', 'ls0', 'ls1']
par = [1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]
I = {n: k for k, n in enumerate(names)}
alg = Alg(par)
y = {n: alg.gen(I[n]) for n in names}
N = len(names)
pairs = [('A0', 'As0'), ('A1', 'As1'), ('A2', 'As2'), ('ls0', 'l0'), ('ls1', 'l1'), ('cs', 'c')]

M = sp.zeros(3, 3)
for a in range(2):
    for b in range(2):
        M[a, b] = (h if a == b else 0) - D[a] * D[b]
M[2, 2] = h
assert all(sum(M[a, b] * D[b] for b in range(2)) == 0 for a in range(3))

def Delta0(p):
    out = {}
    for x, xi in pairs:
        out = alg.add(out, alg.d(I[x], alg.d(I[xi], p)))
    return out

def br(f, g):
    pf = alg.parity(f) if f else 0
    return alg.add(Delta0(alg.mul(f, g)), alg.scale(-1, alg.mul(Delta0(f), g)), alg.scale(-(-1) ** pf, alg.mul(f, Delta0(g))))

A = ['A0', 'A1', 'A2']; As = ['As0', 'As1', 'As2']; L = ['l0', 'l1']; Ls = ['ls0', 'ls1']
quadA = {}
for a in range(3):
    for b in range(3):
        quadA = alg.add(quadA, alg.scale(M[a, b], alg.mul(y[A[a]], y[A[b]])))           # (dA,dA)
quadL = {}
for a in range(2):
    for b in range(2):
        quadL = alg.add(quadL, alg.scale(G[a, b], alg.mul(y[L[a]], y[L[b]])))            # (lambda, dslash lambda)
B7 = {}
for k in range(2):
    B7 = alg.add(B7, alg.scale(D[k], alg.mul(y['c'], y[As[k]])))                          # (dc, A^*)
S_free = alg.add(alg.scale(sp.Rational(-1, 4), quadA), alg.scale(sp.Rational(1, 2), quadL), B7)

res = {}
# ---- Q on generators and as a matrix on E:  Q y_j = sum_i Q[j,i] y_i,  Q e_i = sum_j Q[j,i] e_j
Qimg = [br(S_free, y[n]) for n in names]
Qm = sp.zeros(N, N)
for j in range(N):
    for (e, o), v in Qimg[j].items():
        mons = [k for k in range(N) if alg.gen(k) == {(e, o): 1}]
        assert len(mons) == 1, 'Q image not linear'
        Qm[j, mons[0]] = v
Qfun = lambda p: alg.vf(Qimg, p)
res['Q_on_E (target<-source: coefficient)'] = {f'{names[i]}->{names[j]}': str(Qm[j, i]) for i in range(N) for j in range(N) if Qm[j, i] != 0}
res['Q^2=0 on E'] = (Qm * Qm == sp.zeros(N, N))
# classical_spec sec 5 statements, toy form
spec5 = {
    'Qc=-dc[A]': all(Qm[I[A[k]], I['c']] == -D[k] for k in range(2)) and Qm[I['A2'], I['c']] == 0,
    'QA=-1/2 d*dA [A*]': all(Qm[I[As[j]], I[A[k]]] == -sp.Rational(1, 2) * M[j, k] for j in range(3) for k in range(3)),
    'QA*=d*(A*|V)[c*]': all(Qm[I['cs'], I[As[k]]] == D[k] for k in range(2)) and Qm[I['cs'], I['As2']] == 0,
    'Qlambda=+dslash lambda [lambda*]': all(Qm[I[Ls[i]], I[L[j]]] == G[i, j] for i in range(2) for j in range(2)),
}
res['classical_spec_sec5_signs_reproduced_by_Q={S_free,-}'] = spec5

# ---- Q^GF of DESIGN form with unknown coefficients; solve H = [Q,Q^GF] = h id
al, be, ga, de = sp.symbols('alpha beta gamma delta')
GF = sp.zeros(N, N)
for k in range(3):
    GF[I[A[k]], I[As[k]]] = al                       # A^* -> alpha A   (all components)
for k in range(2):
    GF[I['c'], I[A[k]]] = be * D[k]                  # A -> beta d^*A|_V  [c]
    GF[I[As[k]], I['cs']] = ga * D[k]                # c^* -> gamma d c^*  [A^*|_V]
for i in range(2):
    for j in range(2):
        GF[I[L[i]], I[Ls[j]]] = de * G[i, j]         # lambda^* -> delta dslash lambda^*  [lambda]
Hm = sp.expand(Qm * GF + GF * Qm)                    # both odd: graded commutator = anticommutator
sol = sp.solve(list(set(sp.expand(Hm - h * sp.eye(N)))), [al, be, ga, de], dict=True)
res['Q^GF coefficients solving H=h*id'] = [{str(k): str(v) for k, v in s.items()} for s in sol]
assert len(sol) == 1
GF = GF.subs(sol[0]); Hm = sp.expand(Qm * GF + GF * Qm)
res['H = h*id on all components'] = (Hm == h * sp.eye(N))
res['[Q^GF,H]=0'] = (sp.expand(GF * Hm - Hm * GF) == sp.zeros(N, N))
GF2 = sp.expand(GF * GF)
res['(Q^GF)^2 nonzero components'] = {f'{names[i]}->{names[j]}': str(GF2[j, i]) for i in range(N) for j in range(N) if GF2[j, i] != 0}

# ---- tensors in E (x) E
def op1(X, T, parX=1):   # (X (x) 1)
    out = {}
    for (i, j), c in T.items():
        for k in range(N):
            if X[k, i] != 0:
                out[(k, j)] = out.get((k, j), 0) + X[k, i] * c
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

def op2(X, T, parX=1):   # (1 (x) X), Koszul sign (-1)^{|X||e_i|}
    out = {}
    for (i, j), c in T.items():
        for k in range(N):
            if X[k, j] != 0:
                out[(i, k)] = out.get((i, k), 0) + (-1) ** (parX * par[i]) * X[k, j] * c
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

def tadd(*Ts):
    out = {}
    for T in Ts:
        for k, v in T.items():
            out[k] = out.get(k, 0) + v
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}

def tscale(c, T):
    return {k: sp.expand(c * v) for k, v in T.items() if sp.expand(c * v) != 0}

K0 = {}
for x, xi in pairs:
    K0[(I[x], I[xi])] = -1
    K0[(I[xi], I[x])] = -1
# Delta_0 = -partial_{K0}
rng = random.Random(7)
tests = [alg.random_poly(rng, 3, 6) for _ in range(6)]
res['Delta_0 = -partial_{K_0}'] = all(alg.add(Delta0(F), alg.dphi(K0, F)) == {} for F in tests)
res['(Q(x)1+1(x)Q)K_0 = 0'] = tadd(op1(Qm, K0), op2(Qm, K0)) == {}
res['(Q^GF(x)1-1(x)Q^GF)K_0 = 0'] = tadd(op1(GF, K0), tscale(-1, op2(GF, K0))) == {}
res['(Q^GF(x)1+1(x)Q^GF)K_0 = 0 (control)'] = tadd(op1(GF, K0), op2(GF, K0)) == {}

Keps = tscale(sp.exp(-eps * h), K0)
P = tscale(sp.exp(-eps * h) / h, op1(GF, K0))        # int_eps^inf (Q^GF(x)1) e^{-tH} K_0 dt, H = h id
res['P(eps,infty) components'] = {f'{names[i]}(x){names[j]}': str(v) for (i, j), v in P.items()}
res['P graded symmetric P^{ji}=(-1)^{|i||j|}P^{ij}'] = all(sp.expand(P.get((j, i), 0) - (-1) ** (par[i] * par[j]) * v) == 0 for (i, j), v in P.items())
QP = tadd(op1(Qm, P), op2(Qm, P))
res['(Q(x)1+1(x)Q)P(eps,infty) = K_eps'] = tadd(QP, tscale(-1, Keps)) == {}
res['P has even total parity (only even(x)even, odd(x)odd)'] = all(par[i] == par[j] for (i, j) in P)

# ---- operator identities on random functionals
def comm_dP_Q(F):
    return alg.add(alg.dphi(P, Qfun(F)), alg.scale(-1, Qfun(alg.dphi(P, F))))
res['[partial_P, Q] = partial_{QP} = partial_{K_eps} on random F'] = all(alg.add(comm_dP_Q(F), alg.scale(-1, alg.dphi(Keps, F))) == {} for F in tests)
uv_plus, uv_minus = True, True
for F in tests:
    lhs = Qfun(alg.expo(hb, P, F))
    rhs_p = alg.expo(hb, P, alg.add(Qfun(F), alg.scale(-hb, alg.dphi(Keps, F))))     # Q + hbar Delta_eps, Delta_eps=-d_{K_eps}
    rhs_m = alg.expo(hb, P, alg.add(Qfun(F), alg.scale(+hb, alg.dphi(Keps, F))))
    uv_plus = uv_plus and alg.add(lhs, alg.scale(-1, rhs_p)) == {}
    uv_minus = uv_minus and alg.add(lhs, alg.scale(-1, rhs_m)) == {}
res['UV-boundary identity Q e^{hbar d_P} = e^{hbar d_P}(Q + hbar Delta_eps), Delta_eps=-d_{K_eps}'] = uv_plus
res['(control) same with -hbar Delta_eps'] = uv_minus
res['[Delta_eps, Q] = 0'] = all(alg.add(alg.dphi(Keps, Qfun(F)), alg.scale(1, Qfun(alg.dphi(Keps, F)))) == {} or
                                alg.add(alg.dphi(Keps, Qfun(F)), alg.scale(-1, Qfun(alg.dphi(Keps, F)))) == {} for F in tests)
res['Q^2 = 0 on functionals'] = all(Qfun(Qfun(F)) == {} for F in tests)

# ---- what (Q^GF)^2 != 0 breaks: isotropy of Im Q^GF, C07 Lemma 6.5.1, Def 4.3.3(5)
pairing = sp.zeros(N, N)
for x, xi in pairs:
    pairing[I[x], I[xi]] = 1; pairing[I[xi], I[x]] = -1
iso = sp.expand(GF.T * pairing * GF)
res['Im Q^GF isotropic (<Q^GF a, Q^GF b> = 0 for all a,b)'] = (iso == sp.zeros(N, N))
res['nonzero <Q^GF e_i, Q^GF e_j>'] = {f'<GF {names[i]}, GF {names[j]}>': str(iso[i, j]) for i in range(N) for j in range(N) if iso[i, j] != 0}
Pop = GF / h          # P(0,infty) as an operator Q^GF H^{-1}
imGF = GF.columnspace()
lemma651 = all(sp.expand(Pop * Qm * v - v) == sp.zeros(N, 1) for v in imGF)
res['C07 Lemma 6.5.1: P(0,inf) Q = id on Im Q^GF'] = lemma651
res['P(0,inf)Q + Q P(0,inf) = id (homotopy, needs only H=[Q,Q^GF])'] = (sp.expand(Pop * Qm + Qm * Pop) == sp.eye(N))
inter = sp.Matrix.hstack(*Qm.columnspace()).rank() + sp.Matrix.hstack(*imGF).rank() - sp.Matrix.hstack(*(Qm.columnspace() + imGF)).rank()
res['dim(Im Q cap Im Q^GF) (Def 4.3.3(5) needs 0)'] = int(inter)

# ---- quadratic Q-exactness of the log counterterm (Assumption (A), quadratic level) and the coupling direction
NA = {}; NL = {}
for k in range(3):
    NA = alg.add(NA, alg.mul(y[A[k]], y[As[k]]))
for k in range(2):
    NL = alg.add(NL, alg.mul(y[Ls[k]], y[L[k]]))
NC = alg.mul(y['cs'], y['c'])
Psi = alg.add(alg.scale(2, NA), alg.scale(-4, NL), alg.scale(-1, NC))
target = alg.add(alg.scale(-1, quadA), alg.scale(4, quadL), alg.scale(-1, B7))   # -B1-B4 + 4 B6 - B7 (log coefficients)
QPsi = Qfun(Psi)
res['Q(2N_A - 4N_lam - N_c) = +(-B_A + 4 B_lam - B7)'] = alg.add(QPsi, alg.scale(-1, target)) == {}
res['Q(2N_A - 4N_lam - N_c) = -(-B_A + 4 B_lam - B7)'] = alg.add(QPsi, target) == {}
# general rescaling weights (a, l, g): Q(a N_A - l N_lam - g N_c) vs b_A B_A + b_l B_lam + b_7 B7
a_, l_, g_ = sp.symbols('a l g')
Qgen = alg.add(alg.scale(a_, Qfun(NA)), alg.scale(-l_, Qfun(NL)), alg.scale(-g_, Qfun(NC)))
res['Q(a N_A - l N_lam - g N_c) in basis (B_A, B_lam, B7)'] = {
    'B_A': str(sp.simplify(sum(v for k, v in Qgen.items() if k in quadA) / sum(v for k, v in quadA.items() if k in Qgen) if False else 0)),
}
# extract coefficients by matching one monomial of each basis element
def coeff_of(poly, basis):
    k0 = next(iter(basis)); return sp.simplify(poly.get(k0, 0) / basis[k0])
cA, cL, c7 = coeff_of(Qgen, quadA), coeff_of(Qgen, quadL), coeff_of(Qgen, B7)
chk = alg.add(Qgen, alg.scale(-cA, quadA), alg.scale(-cL, quadL), alg.scale(-c7, B7)) == {}
res['Q(a N_A - l N_lam - g N_c) = cA B_A + cL B_lam + c7 B7'] = {'cA': str(cA), 'cL': str(cL), 'c7': str(c7), 'exact': chk}
coup = alg.add(alg.scale(sp.Rational(-1, 4), quadA), alg.scale(sp.Rational(1, 2), quadL))  # S - (A* d_A* + ...)S, quadratic
Qc = Qfun(alg.add(alg.scale(sp.Rational(1, 2), NA), alg.scale(sp.Rational(-1, 2), NL), alg.scale(sp.Rational(-1, 2), NC)))
res['coupling direction (quadratic) = Q(1/2 N_A - 1/2 N_lam - 1/2 N_c)'] = alg.add(Qc, alg.scale(-1, coup)) == {}
res['coupling direction (quadratic) = -Q(1/2 N_A - 1/2 N_lam - 1/2 N_c)'] = alg.add(Qc, coup) == {}

json.dump(res, open('toy_gf.json', 'w'), indent=1, default=str)
for k, v in res.items():
    print(k, '=>', v)
