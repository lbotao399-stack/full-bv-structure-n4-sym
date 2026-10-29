"""Hand-derived closed forms (refuter #2), evaluated without the graph engine.

Loop v0 -> v1 -> v2 -> v0, q0 = l (special edge, X = K_0 - K_eps), q1 = l + p1, q2 = l - p0; f = (X0,[X1,X2]).

(1) (c, A_V, A_V), gaugino triangle c lam lam^*[c] - X(lam^*-lam) - lam lam A[A_V#1] - P_lam - lam lam A[A_V#2] - P_lam:
    edge rule: Koszul sign -1, X numerator -|q0|^2, P_lam numerators -i B(q1), -i B(q2) (k = -q entering the
    first slot), colour sum eps_{X0 b g} eps_{X1 g' g} eps_{X2 b g'} = f, Clifford trace
    tr[A(a1)B(q1)A(a2)B(q2)] = 16[(a1,q1)(a2,q2) - (a1,a2)(q1,q2) + (a1,q2)(a2,q1)]:
        N(l) = -16 f |q0|^2 [(a1,q1)(a2,q2) - (a1,a2)(q1,q2) + (a1,q2)(a2,q1)].
    exchange partner (even legs): N with (p1,a1,X1) <-> (p2,a2,X2), f -> -f.
(2) (c, lam, lam), row c lam lam^*[c] - X(lam^*-lam) - lam lam A[lam#1] - P_{A_V} - lam lam A[lam#2] - P_lam:
    Koszul sign +1, X: -|q0|^2, P_AA: 2 delta (V), P_lam: -i B(q2), sum_{M in V} A(e_M) B(q) A(e_M) = -2 A(q),
    colour -f:      N(l) = 4 i f |q0|^2 (xi1, rho(q2) xi2).
    exchange partner (odd legs): value(1,2) -> -value(2,1).
eps^0 coefficient (Ch.4 Lemma c4-gauss to fourth order, scaled region tau_X in [0,1], tau_P in [1,oo)):
    c0 = 1/(16 pi^2) int dtau [ (1/(4T^3)) lap N(l0) - Q0 lap^2 N / (32 T^4) ],
    l0 = -(tau1 p1 - tau2 p0)/T,  Q0 = (tau0 tau1 |p1|^2 + tau1 tau2 |p2|^2 + tau2 tau0 |p0|^2)/T
    (tau0: X edge, tau1: edge q1, tau2: edge q2), masters J3 checked by quadrature (rb_masters_check.json).
Units: divide by tr(ad ad)/(X,Y) = -2 (su(2), (e_a,e_b) = delta_ab).
"""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
import numpy as np
from fractions import Fraction as Fr
from rb_loop import J3
from rb_alg import Amat, p10

l = sp.symbols('l0:4')
L32 = sp.Symbol('L')       # log(3/2)


def vec(v):
    return sp.Matrix([sp.Rational(int(x)) for x in v])


def dotp(u, v):
    return sum(u[i] * v[i] for i in range(4))


def lap(expr):
    return sum(sp.diff(expr, l[i], 2) for i in range(4))


def master(ex, m):
    """ex = (a, b, c) exponents of (tau0 [X], tau1, tau2), returns rational + coefficient*L."""
    r, g = J3(ex[0], ex[1], ex[2], m)
    return sp.Rational(r.numerator, r.denominator) + sp.Rational(g.numerator, g.denominator) * L32


t0, t1, t2, T = sp.symbols('t0 t1 t2 T', positive=True)


def eps0(N, p0, p1, p2):
    lv = sp.Matrix(l)
    L1 = sp.expand(lap(N))
    L2 = sp.expand(lap(L1))
    assert all(sp.diff(L2, l[i]) == 0 for i in range(4))
    # l0 = -(t1 p1 - t2 p0)/T expressed with x1 = t1/T, x2 = t2/T
    x1, x2 = sp.symbols('x1 x2')
    sub = {l[i]: -(x1 * p1[i] - x2 * p0[i]) for i in range(4)}
    A1 = sp.Poly(sp.expand(L1.subs(sub)), x1, x2)
    total = 0
    for (e1, e2), co in A1.terms():
        total += sp.Rational(1, 4) * co * master((0, e1, e2), 3 + e1 + e2)
    q = {(1, 1, 0): dotp(p1, p1), (0, 1, 1): dotp(p2, p2), (1, 0, 1): dotp(p0, p0)}
    for ex, co in q.items():
        total -= L2 * co / 32 * master(ex, 5)
    return sp.expand(total)


rng = np.random.default_rng(2026)
res = {'gaugino_triangle_cAVAV': [], 'row4_clamlam': []}
EPS = np.zeros((3, 3, 3))
for (i, j, k), sg in {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.items():
    EPS[i, j, k] = sg
rows = []
for trial in range(4):
    while True:
        P1 = rng.integers(-3, 4, 4); P2 = rng.integers(-3, 4, 4); P0 = -P1 - P2
        if P0.any() and P1.any() and P2.any():
            break
    p0, p1, p2 = vec(P0), vec(P1), vec(P2)
    a1, a2 = vec(rng.integers(-3, 4, 4)), vec(rng.integers(-3, 4, 4))
    X0, X1, X2 = [rng.integers(-3, 4, 3) for _ in range(3)]
    f = int(np.einsum('abc,a,b,c->', EPS, X0, X1, X2))
    lv = sp.Matrix(l)

    def Ngau(pa, pb, aa, ab, ff):
        q0 = lv; q1 = lv + pa; q2 = lv - p0
        h = dotp(aa, q1) * dotp(ab, q2) - dotp(aa, ab) * dotp(q1, q2) + dotp(aa, q2) * dotp(ab, q1)
        return sp.expand(-16 * ff * dotp(q0, q0) * h)

    # graph a (X next to leg 1) and its exchange partner (legs 1 <-> 2, f -> -f)
    ca = eps0(Ngau(p1, p2, a1, a2, f), p0, p1, p2)
    cb = eps0(Ngau(p2, p1, a2, a1, -f), p0, p2, p1)
    val = sp.expand((ca + cb) / (-2))            # units (the 1/16pi^2 is the unit's)
    # G1, G2 plane-wave values (Prop s2c-dims, checked against my Grassmann evaluator in rb_plugback)
    G1 = f * ((dotp(p2, p2) - dotp(p1, p1)) * dotp(a1, a2) - dotp(a1, p2) * dotp(a2, p2) + dotp(a1, p1) * dotp(a2, p1))
    G2 = f * (dotp(a1, p2) * dotp(a2, p2) - dotp(a1, p1) * dotp(a2, p1))
    rows.append((val, G1, G2))
x, y, u, w = sp.symbols('x y u w')
eqs = []
for val, G1, G2 in rows:
    rat = val.subs(L32, 0); lg = sp.expand(val).coeff(L32)
    eqs += [sp.Eq(rat, x * G1 + y * G2), sp.Eq(lg, u * G1 + w * G2)]
sol = sp.solve(eqs, [x, y, u, w], dict=True)
res['gaugino_triangle_cAVAV'] = {'solution(G1 rat, G2 rat, G1 log, G2 log)': {str(k): str(v) for k, v in sol[0].items()} if sol else 'INCONSISTENT',
                                 'claimed': '(10/9 + 8/3 L) G1 + 4/3 G2', 'n_equations': len(eqs)}
print(res['gaugino_triangle_cAVAV'], flush=True)

rows = []
for trial in range(4):
    while True:
        P1 = rng.integers(-3, 4, 4); P2 = rng.integers(-3, 4, 4); P0 = -P1 - P2
        if P0.any() and P1.any() and P2.any():
            break
    p0, p1, p2 = vec(P0), vec(P1), vec(P2)
    xi1 = rng.integers(-2, 3, 16); xi2 = rng.integers(-2, 3, 16)
    X0, X1, X2 = [rng.integers(-3, 4, 3) for _ in range(3)]
    f = int(np.einsum('abc,a,b,c->', EPS, X0, X1, X2))
    Am = [sp.Matrix(16, 16, lambda i, j: sp.nsimplify(complex(Amat(np.eye(10)[M])[i, j]).real) +
                    sp.I * sp.nsimplify(complex(Amat(np.eye(10)[M])[i, j]).imag)) for M in range(4)]
    def spin(xa, xb, q):            # (xa, rho(q) xb) = xa^T A(q) xb
        Aq = sum((q[M] * Am[M] for M in range(4)), sp.zeros(16, 16))
        return (sp.Matrix(list(xa)).T * Aq * sp.Matrix(list(xb)))[0, 0]
    lv = sp.Matrix(l)

    def Nrow4(xa, xb, ff):
        q0 = lv; q2 = lv - p0
        return sp.expand(4 * sp.I * ff * dotp(q0, q0) * spin(xa, xb, q2))

    def eps0_lin(N, p0_, p1_):
        # omega = 1: only (1/(4T^3)) lap N(l0)
        L1 = sp.expand(lap(N))
        x1, x2 = sp.symbols('x1 x2')
        sub = {l[i]: -(x1 * p1_[i] - x2 * p0_[i]) for i in range(4)}
        A1 = sp.Poly(sp.expand(L1.subs(sub)), x1, x2)
        tot = 0
        for (e1, e2), co in A1.terms():
            tot += sp.Rational(1, 4) * co * master((0, e1, e2), 3 + e1 + e2)
        return sp.expand(tot)

    ca = eps0_lin(Nrow4(xi1, xi2, f), p0, p1)
    ca_swapped = eps0_lin(Nrow4(xi2, xi1, -f), p0, p2)       # graph a with legs (2,1)
    val = sp.expand((ca - ca_swapped) / (-2))                 # partner = -(graph a with legs exchanged)
    G5 = -sp.I * f * spin(xi1, xi2, p0)
    rows.append((val, G5))
eqs = []
for val, G5 in rows:
    rat = val.subs(L32, 0); lg = sp.expand(val).coeff(L32)
    eqs += [sp.Eq(rat, x * G5), sp.Eq(lg, u * G5)]
sol = sp.solve(eqs, [x, u], dict=True)
res['row4_clamlam'] = {'solution(G5 rat, G5 log)': {str(k): str(v) for k, v in sol[0].items()} if sol else 'INCONSISTENT',
                       'claimed': '(1/2 - L) G5', 'n_equations': len(eqs)}
print(res['row4_clamlam'], flush=True)
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rb_hand.json'), 'w'), indent=1)
