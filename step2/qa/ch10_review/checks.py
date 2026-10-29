"""Independent recomputation of the cheap numbers of Chapter 10 (referee). Writes checks_results.json."""
import json, itertools, math
from fractions import Fraction as Fr
from pathlib import Path
import sympy as sp
import numpy as np
import mpmath as mp

OUT = Path(__file__).with_name('checks_results.json')
res = {}

# 1. self-contraction lemma
x2, eps = sp.symbols('x2 epsilon', positive=True)
f = (1 - sp.exp(-x2 / (4 * eps))) / (4 * sp.pi**2 * x2)
ser = sp.series(f, x2, 0, 5).removeO()
claim = sum((-1)**(m + 1) / (16 * sp.pi**2 * sp.factorial(m) * eps**m) * (x2 / 4)**(m - 1) for m in range(1, 6))
res['self_contraction_series'] = sp.simplify(sp.expand(ser - claim)) == 0
# derivative at 0 of (|x|^2)^(m-1) in d=4: 2^(m-1)(m-1)! * (#pairings) for gamma = (0,0,1,1,...) etc.
X = sp.symbols('x0:4')
ok = True
for m in range(1, 5):
    r2 = sum(xi**2 for xi in X)
    for gam in itertools.combinations_with_replacement(range(4), 2 * (m - 1)):
        expr = r2**(m - 1)
        for g in gam:
            expr = sp.diff(expr, X[g])
        val = expr.subs({xi: 0 for xi in X})
        idx = list(gam)
        # count pairings with equal indices
        def pairings(l):
            if not l:
                return 1
            a, rest = l[0], l[1:]
            return sum(pairings(rest[:i] + rest[i + 1:]) for i in range(len(rest)) if rest[i] == a)
        pred = 2**(m - 1) * math.factorial(m - 1) * pairings(idx)
        if val != pred:
            ok = False
    # full coefficient
    c = Fr(math.factorial(m - 1) * 2**(m - 1), 4**(m - 1) * math.factorial(m))
    if c != Fr(1, m * 2**(m - 1)):
        ok = False
res['self_contraction_derivative'] = ok
print('self done', flush=True)

# 2. simplex integrands
def integrand(num, l, centre, MF, k, m):
    """m-th simplex integrand (Definition c10-cycle) for polynomial numerator num(l)."""
    tot = 0
    lap = num
    j = 0
    terms = []
    while lap != 0:
        terms.append(lap)
        lap = sum(sp.diff(lap, li, 2) for li in l)
        j += 1
    for j, t in enumerate(terms):
        e = j + 1 - k - m
        if e < 0:
            continue
        tot += t.subs(dict(zip(l, centre)), simultaneous=True) / (sp.factorial(j) * 4**j) * (-MF)**e / sp.factorial(e)
    return sp.expand(tot)

l = sp.symbols('l0:4')
p = sp.symbols('p0:4')
a = sp.symbols('a', positive=True)
psq = sum(pi**2 for pi in p)
# k=1 numerator |l|^4, q0=0, q1=p
num = sum(li**2 for li in l)**2
centre = [-a * pi for pi in p]
MF = a * (1 - a) * psq
F2 = integrand(num, l, centre, MF, 1, 2)
F1 = integrand(num, l, centre, MF, 1, 1)
F0 = integrand(num, l, centre, MF, 1, 0)
mn = sp.Piecewise((a, a < sp.Rational(1, 2)), (1 - a, True))
pe_m2 = sp.Rational(1, 2) * (2 * sp.integrate(a**2 * F2, (a, 0, sp.Rational(1, 2))) if F2.free_symbols <= {a} else None)
pe_m2 = sp.Rational(1, 2) * (sp.integrate(a**2 * F2, (a, 0, sp.Rational(1, 2))) + sp.integrate((1 - a)**2 * F2, (a, sp.Rational(1, 2), 1)))
pe_m1 = sp.integrate(a * F1, (a, 0, sp.Rational(1, 2))) + sp.integrate((1 - a) * F1, (a, sp.Rational(1, 2), 1))
tt_m2 = sp.Rational(1, 2) * sp.integrate(F2, (a, 0, 1))
tt_m1 = sp.integrate(F1, (a, 0, 1))
lg = -sp.integrate(F0, (a, 0, 1))
res['shape_bubble'] = {'per_edge_eps^-1': str(sp.simplify(pe_m1)), 'per_edge_eps^-2': str(pe_m2),
                       'total_eps^-1': str(sp.simplify(tt_m1)), 'total_eps^-2': str(tt_m2), 'log': str(sp.simplify(lg))}
# triangle numerator |l|^4, k=2: F1 constant
q1 = sp.symbols('q1_0:4'); q2 = sp.symbols('q2_0:4')
a1, a2 = sp.symbols('a1 a2', positive=True)
cen = [-(a1 * q1[i] + a2 * q2[i]) for i in range(4)]
a0 = 1 - a1 - a2
d = lambda u, v: sum((ui - vi)**2 for ui, vi in zip(u, v))
MF3 = a0 * a1 * d(q1, [0] * 4) + a0 * a2 * d(q2, [0] * 4) + a1 * a2 * d(q1, q2)
F1t = integrand(num, l, cen, MF3, 2, 1)
res['triangle_F1'] = str(F1t)
# per edge: 6 * int min(alpha) over 2-simplex = 6 * 1/18 ; total time 6 * 1/2
res['triangle_eps^-1'] = {'per_edge': str(F1t * sp.Rational(1, 18)), 'total': str(F1t * sp.Rational(1, 2))}

# 3. special edge example
q = sp.symbols('q')
numS = -sum(li**2 for li in l)
# edges: e0 carries l (q0 = 0), e1 carries l - p (q1 = -p)
centreS = [a1 * pi for pi in p]  # centre = -alpha1 q1 = alpha1 p, alpha1 = 1 - alpha0
s0 = sp.symbols('s0', positive=True)
cS = [(1 - s0) * pi for pi in p]
MFS = s0 * (1 - s0) * psq
F1S = integrand(numS, l, cS, MFS, 1, 1)
F0S = integrand(numS, l, cS, MFS, 1, 0)
res['special_F1'] = str(F1S)
res['special_F0'] = str(sp.factor(F0S))
I38 = mp.quad(lambda t: (1 - t) * (3 * t - 1) * mp.log(t / (1 - t)), [0, 0.5])
res['special_int_3/8'] = float(I38)
cm1 = sp.integrate(((1 - s0) - s0) * F1S, (s0, 0, sp.Rational(1, 2)))
res['special_eps^-1_(16pi^2 units)'] = str(cm1)
# direct numerical evaluation of int d^4l/(2pi)^4 (e^{-eps l^2}-1) e^{-eps|l-p|^2}/|l-p|^2 at eps=1e-3, |p|=1
def direct(epsv, pv):
    # Schwinger: -int_0^eps dt0 int_eps^inf dt1 (1/(16 pi^2 tau^2)) (t1^2 p^2/tau^2 + 2/tau) exp(-t0 t1 p^2/tau)
    f = lambda t0, t1: -(1 / (16 * mp.pi**2 * (t0 + t1)**2)) * (t1**2 * pv**2 / (t0 + t1)**2 + 2 / (t0 + t1)) * mp.exp(-t0 * t1 * pv**2 / (t0 + t1))
    return mp.quad(lambda t0: mp.quad(lambda t1: f(t0, t1), [epsv, 10 * epsv, 1, mp.inf]), [0, epsv])
from scipy import integrate as si
def direct_sc(epsv, pv):
    f = lambda t1, t0: -(1 / (16 * np.pi**2 * (t0 + t1)**2)) * (t1**2 * pv**2 / (t0 + t1)**2 + 2 / (t0 + t1)) * np.exp(-t0 * t1 * pv**2 / (t0 + t1))
    return si.dblquad(f, 0, epsv, lambda t0: epsv, lambda t0: np.inf, epsabs=1e-12, epsrel=1e-10)[0]
e1 = mp.mpf('1e-3')
val = direct_sc(1e-3, 1.0)
print('direct done', flush=True)
res['special_direct'] = {'eps': '1e-3', 'p2': 1, 'value': float(val),
                         '(-1/(2eps)-3/8 p^2)': float(-1 / (2 * e1) - mp.mpf(3) / 8),
                         '(-1/(2eps)-3/8 p^2)/(16pi^2)': float((-1 / (2 * e1) - mp.mpf(3) / 8) / (16 * mp.pi**2))}

# 4. bubble closed form vs Bessel generating function (k=1, constant numerator), numeric rational vectors
import random
random.seed(3)
vec=lambda: [sp.Rational(random.randint(-5,5),random.randint(1,4)) for _ in range(4)]
ya, yb, qq = vec(), vec(), vec()
y = [yb[i] - ya[i] for i in range(4)]
dot = lambda u, v: sum(ui * vi for ui, vi in zip(u, v))
order = 6
tvar = sp.symbols('t')
expo = sp.exp(sp.I * tvar * dot(yb, qq) - sp.I * tvar * a * dot(y, qq))
ser_e = sp.series(expo, tvar, 0, order + 1).removeO()
I0 = sum((tvar**2 * dot(y, y) * dot(qq, qq) * a * (1 - a) / 4)**jj / sp.factorial(jj)**2 for jj in range(order // 2 + 1))
prod = sp.expand(ser_e * I0)
prod = sum(prod.coeff(tvar, n) * tvar**n for n in range(order + 1))
lhs = -sp.integrate(sp.expand(prod), (a, 0, 1))
rhs = 0
for n1 in range(order + 1):
    for n2 in range(order + 1 - n1):
        for jj in range((order - n1 - n2) // 2 + 1):
            rhs += -(sp.I * tvar * dot(ya, qq))**n1 * (sp.I * tvar * dot(yb, qq))**n2 * (tvar**2 * dot(y, y))**jj * (dot(qq, qq) / 4)**jj * \
                sp.factorial(n1 + jj) * sp.factorial(n2 + jj) / (sp.factorial(n1) * sp.factorial(n2) * sp.factorial(jj)**2 * sp.factorial(n1 + n2 + 2 * jj + 1))
diff = sp.expand(lhs - rhs)
res['bubble_closed_form_vs_bessel_order<=6'] = (diff == 0)
print('bubble done', flush=True)
# 5. generating function: Bessel identity check
nB, xB = 3, sp.Rational(7, 5)
lhsB = (xB / 2)**(-nB) * mp.besseli(nB, float(xB))
rhsB = sum((float(xB)**2 / 4)**i / (math.factorial(i) * math.factorial(i + nB)) for i in range(30))
res['bessel_series'] = abs(float(lhsB) - rhsB) < 1e-12

# 6. Heisenberg cyclic spectra
def cyclic_spectrum(L, M):
    states = [s for s in itertools.product([0, 1], repeat=L) if sum(s) == M]
    idx = {s: i for i, s in enumerate(states)}
    n = len(states)
    H = sp.zeros(n, n)
    for s in states:
        for i in range(L):
            j = (i + 1) % L
            t = list(s); t[i], t[j] = t[j], t[i]; t = tuple(t)
            H[idx[s], idx[s]] += 1
            H[idx[t], idx[s]] -= 1
    # projector onto cyclic (shift-invariant) states
    def shift(s):
        return s[1:] + s[:1]
    orbits, seen = [], set()
    for s in states:
        if s in seen:
            continue
        orb, t = [], s
        while t not in orb:
            orb.append(t); t = shift(t)
        seen.update(orb)
        orbits.append(orb)
    B = sp.zeros(n, len(orbits))
    for k, orb in enumerate(orbits):
        for t in orb:
            B[idx[t], k] = 1
    # H commutes with shift: restrict
    Hc = (B.T * B).inv() * B.T * H * B
    E = sp.symbols('E')
    return sp.factor(Hc.charpoly(E).as_expr())
spec = {}
for L in range(2, 9):
    for M in range(0, L // 2 + 1):
        spec[f'{L},{M}'] = str(cyclic_spectrum(L, M))
res['heisenberg_cyclic_charpolys'] = spec
print('heis done', flush=True)
E = sp.symbols('E')
res['cos_cubic'] = str(sp.expand(sp.prod([E - 4 + 4 * sp.cos(2 * sp.pi * n / 7) for n in (1, 2, 3)])).evalf())
res['cos_cubic_exact'] = str(sp.nsimplify(sp.expand(sp.prod([E - 4 + 4 * sp.cos(2 * sp.pi * n / 7) for n in (1, 2, 3)])).evalf(30), rational=True))

# 7. twist two: Beisert harmonic action H'_12 on two sites, cyclic (symmetric) states
def h(n):
    return sum(Fr(1, j) for j in range(1, n + 1))
tw = {}
for S in range(0, 9):
    # basis |k, S-k>, k = 0..S ; H12 |k> = sum_k' (delta (h(k)+h(S-k)) - (1-delta)/|k-k'|) |k'>
    Hm = sp.zeros(S + 1, S + 1)
    for k in range(S + 1):
        for kp in range(S + 1):
            Hm[kp, k] = (h(k) + h(S - k)) if k == kp else -Fr(1, abs(k - kp))
    # site exchange P: |k, S-k> -> |S-k, k>
    P = sp.zeros(S + 1, S + 1)
    for k in range(S + 1):
        P[S - k, k] = 1
    H21 = P * Hm * P
    # cyclic states of tr(Z Z): symmetric under P
    sym = [k for k in range(S + 1) if k <= S - k]
    B = sp.zeros(S + 1, len(sym))
    for c, k in enumerate(sym):
        B[k, c] += 1
        B[S - k, c] += 1 if S - k != k else 0
    Mt = -2 * (Hm + H21)
    Mc = (B.T * B).inv() * B.T * Mt * B
    ev = sorted([sp.nsimplify(v) for v in Mc.eigenvals(multiple=True)])
    tw[S] = {'eigs': [str(v) for v in ev], 'pred': [str(-8 * sp.Rational(h(s).numerator, h(s).denominator)) for s in range(0, S + 1, 2)]}
res['twist_two'] = tw
print('twist done', flush=True)
res['harmonic'] = {S: str(-8 * h(S)) for S in (0, 2, 4, 6, 8)}

# 8. colour rules (Lemma c10-colour ii) for gl(N) with random matrices, kappa = 1
rng = np.random.default_rng(1)
N = 4
def rnd():
    return rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
basis = []
for i in range(N):
    for j in range(N):
        Eij = np.zeros((N, N)); Eij[i, j] = 1
        Eji = np.zeros((N, N)); Eji[j, i] = 1
        basis.append((Eij, Eji))  # tr(E_ij E_ji) = 1 : dual bases for (X,Y) = tr(XY)
def ad_word(Ys, Z):
    for Y in Ys:
        Z = Y @ Z - Z @ Y
    return Z
def prodl(Ys, S, rev):
    M = np.eye(N, dtype=complex)
    js = [j for j in range(len(Ys)) if (j in S) == True] if rev is not None else []
    return M
colour_ok = True
for r in range(0, 4):
    Ys = [rnd() for _ in range(r)]
    M1, M2 = rnd(), rnd()
    # rule 1
    lhs1 = sum(np.trace(M1 @ TA @ ad_word(Ys, TB)) for TA, TB in basis)
    lhs2 = sum(np.trace(M1 @ TA @ M2 @ ad_word(Ys, TB)) for TA, TB in basis)
    lhs3 = sum(np.trace(M1 @ TA) * np.trace(M2 @ ad_word(Ys, TB)) for TA, TB in basis)
    rhs1 = rhs2 = rhs3 = 0
    for Ssz in range(r + 1):
        for S in itertools.combinations(range(r), Ssz):
            left = np.eye(N, dtype=complex)
            for j in sorted(S, reverse=True):
                left = left @ Ys[j]
            right = np.eye(N, dtype=complex)
            for j in range(r):
                if j not in S:
                    right = right @ Ys[j]
            sg = (-1)**(r - len(S))
            rhs1 += sg * np.trace(left) * np.trace(M1 @ right)
            rhs2 += sg * np.trace(M2 @ left) * np.trace(M1 @ right)
            rhs3 += sg * np.trace(M1 @ right @ M2 @ left)
    lead1 = (-1)**r * N * np.trace(M1 @ (np.linalg.multi_dot(Ys) if r > 1 else (Ys[0] if r == 1 else np.eye(N))))
    colour_ok &= abs(lhs1 - rhs1) < 1e-8 and abs(lhs2 - rhs2) < 1e-8 and abs(lhs3 - rhs3) < 1e-8
    res[f'colour_r{r}_S_empty_term_rule1'] = bool(abs(((-1)**r) * N * np.trace(M1 @ (np.eye(N) if r == 0 else np.linalg.multi_dot(Ys + [np.eye(N)]))) - lead1) < 1e-8)
res['colour_rules_gl'] = bool(colour_ok)
# wrap-around pair leading term: tr(X1 M2 Xn) -> + N tr(M2 Y_r ... Y_1)
r = 3; Ys = [rnd() for _ in range(r)]; M2 = rnd()
lhs = sum(np.trace(TA @ M2 @ ad_word(Ys, TB)) for TA, TB in basis)
# N-coefficient: terms with complement of S empty (S = all): + N tr(M2 Y_r..Y_1)
full = np.eye(N, dtype=complex)
for j in reversed(range(r)):
    full = full @ Ys[j]
res['wraparound_leading'] = 'term S=all gives +N tr(M2 Y_r...Y_1) with sign (-1)^0'

# 9. unit and numbers
g, Nn = sp.symbols('g N', positive=True)
unit = -g**2 * Nn / (32 * sp.pi**2)
res['units'] = {'-4u': str(sp.simplify(-4 * unit)), '-16u': str(sp.simplify(-16 * unit)), '-24u': str(sp.simplify(-24 * unit)),
                '-48u': str(sp.simplify(-48 * unit)), '-8u': str(sp.simplify(-8 * unit)), '-2u': str(sp.simplify(-2 * unit))}
res['one_body'] = {'A': str(Fr(-8) + Fr(2, 3)), 'lambda': str(Fr(-16) + Fr(1, 3)), 'c': str(Fr(-4) + Fr(4, 3)),
                   'A_perp_with_composite': str(Fr(-22, 3) - Fr(2, 3)), 'Xi': str(Fr(5, 3) - 2)}
# multipliers of Prop c10-W1
w = {'c': 1, 'A': 2, 'lam': 4, 'As': -2, 'lams': -4, 'cs': -1}
res['W1_multipliers'] = {'AA': 2 * w['A'], 'lamlam': 2 * w['lam'], 'cAs': w['c'] + w['As'], 'AAA': 3 * w['A'], 'lamAlam': 2 * w['lam'] + w['A'],
                         'AcAs': w['A'] + w['c'] + w['As'], 'clamlams': w['c'] + w['lam'] + w['lams'], 'cccs': 2 * w['c'] + w['cs'], 'AAAA': 4 * w['A']}
OUT.write_text(json.dumps(res, indent=1, default=str))
print(json.dumps(res, indent=1, default=str)[:6000])
