"""Audit 2b: independent recomputation of the (A_perp) sector of Ob(Delta) (Part D, thm:d-D2 (ii)).

Own Feynman rules from the manuscript only (ch01 vertex functions eq:c1-AAA and the AcA* table, ch03 Prop c3-P
propagator table, ch01 Def c1-bracket K_0, Example c2-dil S^(1)(Delta)), own Wick/Grassmann bookkeeping, exact integrals.
Units: hbar tr_g(ad ad)/(16 pi^2 (X,Y)); su(2) with (T_a,T_b) = delta_ab, [T_a,T_b] = eps_abc T_c (ratio tr(ad ad)/(,) = -2).

 (1) normalisation check: singular part of the Step-1 graph AW,AW:g3[VWW]^2 (A_perp self-energy bubble, lines A_V and
     A_perp): expected (compare/compare.json) eps^-1 -> -(1/4) int(A_perp,A_perp), log eps -> +1 int(dA_perp,dA_perp).
 (2) Q Pi_Psi|(A_perp,A_perp) = Q(class (i)) = (hbar/2) d_Psi(Q W_0) = -hbar x [special-edge composite triangle]:
     vertices AAA(phi_ext, phi_int, A_V), ([phi_ext,c],A*_perp), S^(1)(Delta); edges K_eps - K_0 (phi_int -- A*_perp),
     P_AA (A_V -- S.A), P_{A* c} (S.A*_V -- c); both assignments of the two external phi legs; exact eps -> 0 limit
     (Taylor expansion in the external momenta, angular average on S^3, exact radial integrals).
     Engine A (Q of its Pi_int) and engine B (Ob - 4 C_log) both give +(1/3) int(dA_perp,dA_perp), no x-part.
 Ob(Delta)|_perp = 4 x (Step-1 log coefficient of int(dA_perp,dA_perp)) + (2) = -4 + (2).
Canonical two-leg value: X1 X2 [alpha(p) delta(p+p2) + beta^nu(p) d/dp2_nu delta(p+p2)], leg X1 at p fixed, p2 = -p-kappa;
S_E = sum_nu (-i d/dkappa_nu)|_0 of the vertex 2 int e^{i kappa x}(d_nu A, A*) (injection kappa), S_w = 2 int(A,A*).
Writes perp_sector_from_scratch.json.  Pure sympy, single process.
"""
import itertools, json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
eps = sp.Symbol('epsilon', positive=True)
u = sp.Symbol('u', positive=True)            # |m|^2
t = sp.Symbol('t')                          # bookkeeping for powers of the external momentum p
sk = sp.Symbol('sk')                        # bookkeeping for powers of kappa
m = sp.Matrix(sp.symbols('m0:4', real=True))
kap = sp.Matrix(sp.symbols('k0:4', real=True))
PV = sp.Matrix([1, 2, 0, -1])               # generic external momentum, |p|^2 = 6
PSQ = 6
NV, NW = 4, 6
DIM = NV + NW
COL = 3
TMAX, SMAX = 3, 1

def trunc(expr):
    """drop t^a sk^b with a > TMAX or b > SMAX"""
    expr = sp.expand(expr)
    out = 0
    for term in sp.Add.make_args(expr):
        pd = term.as_powers_dict()
        a = pd.get(t, 0); b = pd.get(sk, 0)
        if a <= TMAX and b <= SMAX:
            out += term
    return out

def dot(a, b):
    return sum(a[i] * b[i] for i in range(4))

def exp_shift(a, nmax=5):
    """e^{-eps(2 m.a + a.a)} truncated; a is O(t) + O(sk)"""
    X = 2 * dot(m, a) + dot(a, a)
    return trunc(sum((-eps * X) ** n / sp.factorial(n) for n in range(nmax)))

def inv_shift(a, nmax=5):
    """u / |m + a|^2 = 1/(1 + (2 m.a + a.a)/u) truncated"""
    X = (2 * dot(m, a) + dot(a, a)) / u
    return trunc(sum((-X) ** n for n in range(nmax)))

def levi(a, b, c):
    return int(sp.LeviCivita(a, b, c))

def unitvec(i):
    v = [0] * DIM; v[i] = 1; return v

def dotV(e, k):
    return sum(e[i] * k[i] for i in range(NV))

def dot10(e1, e2):
    return sum(e1[i] * e2[i] for i in range(DIM))

def V_AAA(e1, k1, c1, e2, k2, c2, e3, k3, c3):
    """eq:c1-AAA, legs (A(k1),A(k2),A(k3)) all incoming, colours c_i: -(i/2)[(i_{k2-k1}A3,[A1,A2]) + cyclic]"""
    f = levi(c1, c2, c3)
    if f == 0:
        return 0
    d21 = [k2[i] - k1[i] for i in range(NV)]; d32 = [k3[i] - k2[i] for i in range(NV)]; d13 = [k1[i] - k3[i] for i in range(NV)]
    return -sp.I / 2 * f * (dotV(e3, d21) * dot10(e1, e2) + dotV(e1, d32) * dot10(e2, e3) + dotV(e2, d13) * dot10(e3, e1))

def V_AcAs(eA, cA, cc, es, cs):
    """([A,c],A*) vertex function, legs (A, c, A*), coefficient of r_c r_A*"""
    return dot10(eA, es) * levi(cA, cc, cs)

def koszul(parities, order):
    seq = list(order); sign = 1
    for i in range(len(seq)):
        for j in range(len(seq) - 1 - i):
            if seq[j] > seq[j + 1]:
                if parities[seq[j]] and parities[seq[j + 1]]:
                    sign = -sign
                seq[j], seq[j + 1] = seq[j + 1], seq[j]
    return sign

def sphere_avg(exps):
    if any(a % 2 for a in exps):
        return 0
    n = sum(exps) // 2
    num = sp.prod([sp.gamma(sp.Rational(a + 1, 2)) for a in exps])
    return sp.nsimplify(sp.simplify(num / (sp.gamma(sp.Rational(1, 2)) ** 4 * sp.gamma(n + 2))))

def angular_average(expr):
    expr = sp.expand(expr)
    poly = sp.Poly(expr, *list(m))
    out = 0
    for mon, coeff in poly.terms():
        avg = sphere_avg(mon)
        if avg != 0:
            out += coeff * avg * u ** (sum(mon) // 2)
    return sp.expand(out)

def split_term(term):
    """term = c * eps^j * u^b * exp(-N eps u) -> (c, j, b, N); c free of eps, u"""
    c = 1; j = 0; b = 0; N = 0
    for fac in sp.Mul.make_args(term):
        if fac.func == sp.exp:
            N += sp.simplify(-sp.expand(fac.args[0]) / (eps * u))
        elif fac == u:
            b += 1
        elif fac.is_Pow and fac.base == u:
            b += fac.exp
        elif fac == eps:
            j += 1
        elif fac.is_Pow and fac.base == eps:
            j += fac.exp
        else:
            c *= fac
    return c, j, b, N

def radial_exact(expr):
    """16 pi^2 int d^4m/(2pi)^4 (...) = int_0^inf du u (...), exact; terms grouped by (j, b) so that each group's
    u -> 0 behaviour is checked (combinations sum_N c_N e^{-N eps u} with s = 2 + b <= 0 must vanish to that order)."""
    groups = {}
    for term in sp.Add.make_args(sp.expand(expr)):
        c, j, b, N = split_term(term)
        groups.setdefault((j, b), {}).setdefault(N, 0)
        groups[(j, b)][N] += c
    total = 0
    for (j, b), byN in groups.items():
        byN = {N: sp.simplify(c) for N, c in byN.items() if sp.simplify(c) != 0}
        if not byN:
            continue
        s_ = 2 + b
        if s_ > 0:
            for N, c in byN.items():
                total += c * eps ** j * sp.gamma(s_) / (N * eps) ** s_
        elif s_ == 0:
            if sp.simplify(sum(byN.values())) != 0:
                raise ValueError('IR log divergence at u->0 in group (j,b)=%s: %s' % ((j, b), byN))
            for N, c in byN.items():
                total += -c * eps ** j * sp.log(N)
        else:
            # s = -1: needs sum c_N = 0 and sum c_N N = 0; analytic continuation Gamma(s) sum c_N (N eps)^{-s}
            if s_ == -1 and sp.simplify(sum(byN.values())) == 0 and sp.simplify(sum(N * c for N, c in byN.items())) == 0:
                for N, c in byN.items():
                    total += c * eps ** j * N * sp.log(N) * eps
                continue
            raise ValueError('IR divergent group (j,b)=%s: %s' % ((j, b), byN))
    return sp.expand(total)

# ------------------------------------------------------------------ (2) special-edge composite triangle
def triangle_integrand(assignment):
    """integrand (polynomial in m, u, t, sk times exponentials in eps u) of the special-edge triangle, summed over
    internal indices, with the Koszul sign; returns (Gw, [G_nu]) before the m-integration."""
    p = PV * t
    kv = kap * sk
    p2 = -p - kv
    q1, q2 = (p, p2) if assignment == 1 else (p2, p)
    uflav = unitvec(4)
    # loop routing: S.A leg momentum k5 = m (incoming at S); P_AA joins V1.A (k2 = -m) to S.A (k5 = m)
    k2 = -m; k1 = -q1 - k2; k5 = m; k6 = -k5 - kv; k4 = -k6; k3 = -k1
    assert sp.simplify(q2 + k3 + k4) == sp.zeros(4, 1)
    # factors with exact exponentials written around |m|^2 = u:
    #   Khat(k1) = -(e^{-eps |m - q1|^2} - 1) = 1 - e^{-eps u} E(-q1)
    #   P_AA(k2) = 2 e^{-eps u}/u
    #   P_{A* c}(k6) = -i k6_nu e^{-eps |k6|^2}/|k6|^2 = i (m + kv)_nu e^{-eps u} E(kv) inv(kv)/u
    Khat = 1 - sp.exp(-eps * u) * exp_shift(-q1)
    PAA = 2 * sp.exp(-eps * u) / u
    Ekap = exp_shift(kv); Ikap = inv_shift(kv)
    PAsc = [trunc(sp.I * (m[nu] + kv[nu]) * Ekap * Ikap) * sp.exp(-eps * u) / u for nu in range(4)]
    par = [0, 0, 0, 0, 1, 1, 0, 1]
    names = ['X_V1', 'phi_int', 'A_V1', 'X_V2', 'c', 'As_V2', 'A_S', 'As_S']
    order_ = [names.index(x) for x in ['phi_int', 'As_V2', 'A_V1', 'A_S', 'As_S', 'c', 'X_V1', 'X_V2']]
    sgn = koszul(par, order_)
    Gw = 0; Gnu = [0, 0, 0, 0]
    for I, M, ca, cb in itertools.product(range(NV, DIM), range(NV), range(COL), range(COL)):
        e_int = unitvec(I); eA = unitvec(M)
        v1 = V_AAA(uflav, list(q1), 0, e_int, list(k1), ca, eA, list(k2), cb)
        if v1 == 0:
            continue
        cc = cb                                     # P_{A* c} colour delta; S colour delta; P_AA colour delta
        v2 = V_AcAs(uflav, 0, cc, e_int, ca)        # K-hat delta: A*_V2 has flavour I, colour ca
        if v2 == 0:
            continue
        Np = M                                      # S delta between its A (comp M) and A* (comp N'); A*|_V only
        base = trunc(sgn * v1 * v2 * Khat * PAA * PAsc[Np])
        Gw += 2 * base
        for nu in range(4):
            Gnu[nu] += trunc(2 * sp.I * k5[nu] * base)
    return sp.expand(Gw), [sp.expand(g) for g in Gnu]

def coeff_t_sk(expr, a, b):
    expr = sp.expand(expr)
    return sp.expand(expr.coeff(t, a).coeff(sk, b)) if (a or b) else sp.expand(expr.subs({t: 0, sk: 0}))

def to_units_triangle(expr):
    """-hbar x (pattern) integrated, in units: (-1) * radial / (-2)"""
    ang = angular_average(expr)
    rad = radial_exact(ang)
    return sp.expand(-1 * rad / (-2))

def run_triangle():
    out = {}
    alpha_tot = {a: 0 for a in range(TMAX + 1)}; beta_tot = {a: [0] * 4 for a in range(TMAX + 1)}
    for asg in (1, 2):
        Gw, Gnu = triangle_integrand(asg)
        rec = {}
        for a in range(TMAX + 1):
            # alpha at order p^a: Gw at (t^a, sk^0) + sum_nu (-i) d/dkappa_nu of G_nu at (t^a, sk^1)
            aw = to_units_triangle(coeff_t_sk(Gw, a, 0))
            aE = 0
            for nu in range(4):
                g1 = coeff_t_sk(Gnu[nu], a, 1)
                aE += -sp.I * to_units_triangle(sp.diff(g1, kap[nu]))
            al = sp.simplify(aw + aE)
            be = [sp.simplify(-sp.I * to_units_triangle(coeff_t_sk(Gnu[nu], a, 0))) for nu in range(4)]
            rec['order p^%d' % a] = {'alpha': str(al), 'beta': [str(x) for x in be]}
            alpha_tot[a] += al
            beta_tot[a] = [beta_tot[a][i] + be[i] for i in range(4)]
        out['assignment %d' % asg] = rec
    out['total'] = {('order p^%d' % a): {'alpha': str(sp.simplify(alpha_tot[a])),
                                         'beta': [str(sp.simplify(x)) for x in beta_tot[a]]} for a in range(TMAX + 1)}
    out['coefficient of int(dA_perp,dA_perp) = alpha(order p^2) / (2|p|^2)'] = str(sp.simplify(alpha_tot[2] / (2 * PSQ)))
    return out, sp.simplify(alpha_tot[2] / (2 * PSQ))

# ------------------------------------------------------------------ (1) Step-1 bubble AW,AW:g3[VWW]^2
def bubble_singular():
    """multilinear value on (X1 phi(p), X2 phi(-p)) = hbar int V1 V2 P_phi P_A (Aut = 1); eps^-1 and log eps parts."""
    p = PV * t
    uflav = unitvec(4)
    # A line momentum: k2 = -m at V1, k4 = m at V2; phi line k1 = -p - k2 = m - p at V1, k3 = -k1
    k2 = -m; k1 = -p - k2; k3 = -k1; k4 = -k2
    Pphi = 2 * sp.exp(-eps * u) * exp_shift(-p) * inv_shift(-p) / u      # 2 e^{-eps|m-p|^2}/|m-p|^2
    PA = 2 * sp.exp(-eps * u) / u
    val = 0
    for I, M, ca, cb in itertools.product(range(NV, DIM), range(NV), range(COL), range(COL)):
        e_int = unitvec(I); eA = unitvec(M)
        v1 = V_AAA(uflav, list(p), 0, e_int, list(k1), ca, eA, list(k2), cb)
        v2 = V_AAA(uflav, list(-p), 0, e_int, list(k3), ca, eA, list(k4), cb)
        if v1 == 0 or v2 == 0:
            continue
        val += trunc(v1 * v2 * Pphi * PA)
    val = sp.expand(val.subs(sk, 0))
    res = {}
    for a in (0, 2):
        ang = angular_average(coeff_t_sk(val, a, 0))
        # UV singular part: terms c eps^j u^b e^{-N eps u}: int du u^{1+b} e^{-N eps u} = Gamma(s)(N eps)^{-s}, s = 2+b
        einv = 0; leps = 0
        for term in sp.Add.make_args(ang):
            c, j, b, N = split_term(term)
            s_ = 2 + b
            if s_ >= 1 and j - s_ == -1:
                einv += c * sp.gamma(s_) / N ** s_
            if s_ == 0 and j == 0:
                leps += -c
        res[a] = (sp.simplify(einv / (-2)), sp.simplify(leps / (-2)))
    return {'eps^-1 coefficient of the two-leg value, order p^0 (expected 2 x (-1/4) = -1/2)': str(res[0][0]),
            'log eps coefficient of the two-leg value, order p^2, divided by 2|p|^2 (expected +1)':
                str(sp.simplify(res[2][1] / (2 * PSQ))),
            'eps^-1 at order p^2 (expected 0)': str(res[2][0]), 'log eps at order p^0 (expected 0)': str(res[0][1])}

# ------------------------------------------------------------------ (0) sign of {S^(1)(Delta), -}
def sign_link():
    """{S_M, int(phi,phi)} with {f,g} = -C_{K_0}(f,g), K_0^{A* x A} = -1: the A* leg of S_M is contracted with an A leg
    of g; the S_M A leg replaces it with (2E + 2) acting: {S_M, int(phi,phi)} = -(-1) x 2 int(phi, (2E+2) phi)
    = 4 int(phi,E phi) + 4 int(phi,phi) = (-8 + 4) int(phi,phi) since int(phi, E phi) = -2 int(phi,phi) in 4d.
    D_M int(phi,phi) = 2 int((2E+2)phi, phi) = -4 int(phi,phi).  Both -4: {S_M,-} = +D_M."""
    return {'{S_M, int(phi,phi)} / int(phi,phi)': '-4', 'D_M int(phi,phi) / int(phi,phi)': '-4'}

if __name__ == '__main__':
    res = {}
    res['(0) sign link'] = sign_link()
    res['(1) Step-1 bubble AW,AW:g3[VWW]^2 singular parts (normalisation check)'] = bubble_singular()
    tri, coef = run_triangle()
    res['(2) Q Pi_Psi|(A_perp,A_perp) via the special-edge composite triangle'] = tri
    res['Ob(Delta)|(A_perp,A_perp) = 4 x (-1) + Q Pi_Psi (units, coefficient of int(dA_perp,dA_perp))'] = str(-4 + coef)
    with open(os.path.join(HERE, 'perp_sector_from_scratch.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
