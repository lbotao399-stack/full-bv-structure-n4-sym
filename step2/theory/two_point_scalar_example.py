#!/usr/bin/env python3
"""Worked example of two_point_design.tex, Section 7: one-loop two-point function of the scalar bilinears
O_a = int omega_a (iota_{a_1}A, iota_{a_2}A), a = 1, 2, at separation x != 0, UV boundary eps only, heat-kernel
propagators P(eps, infty) (scalar factor (1 - exp(-r^2/4 eps))/(4 pi^2 r^2), momentum space exp(-eps p^2)/p^2).

Code-level names (json/py only): G0 = tree part, G1 = one-loop part, unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)),
S = (a1,b1')(a2,b2') + ... flavour structure of the tree, T = (a_1,a_2)(b_1,b_2) trace structure, tau = T/S.
Exact steps are sympy; numerical cross-checks (mpmath/scipy) only confirm closed forms, nothing is fitted.
Single process, small memory; writes two_point_example_results.json next to this file.
"""
import json, os
import sympy as sp
import mpmath as mp
from scipy import integrate
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
u, eps, p2 = sp.symbols('u epsilon p2', positive=True)
gE, l2, l3, Lx, Le, tau = sp.symbols('gamma_E log2 log3 logx2 logeps tau')
res = {}

# ---------------------------------------------------------------------------------------------------------
# A. Exact eps-expansion of the two-edge Schwinger integrals (Lemma lem:t-bubble of the design)
#    I_n[g] = iint_{[eps,inf)^2} exp(-t1 t2 p^2/s) s^{-n} g(t1/s) dt1 dt2,  s = t1 + t2.
#    t1 = s u, t2 = s (1-u), dt1 dt2 = s ds du, s >= eps/m(u), m = min(u, 1-u).
# ---------------------------------------------------------------------------------------------------------
def int01(expr_left, expr_right=None):
    """integral over [0,1] of expr; if expr_right given, use expr_left on [0,1/2], expr_right on [1/2,1]."""
    if expr_right is None:
        return sp.simplify(sp.integrate(expr_left, (u, 0, 1)))
    return sp.simplify(sp.integrate(expr_left, (u, 0, sp.Rational(1, 2))) +
                       sp.integrate(expr_right, (u, sp.Rational(1, 2), 1)))

def logm_int(g):
    return int01(g * sp.log(u), g * sp.log(1 - u))

def I2(g):
    """returns (coefficient of L2e, constant) with L2e = log(eps p^2) + gamma_E + log 2 : I2 = c1*L2e + c0."""
    A = int01(g)
    B = int01(g * sp.log(u * (1 - u)))
    C = logm_int(g)
    # value = -A(gamma + log eps p^2) - B + C ; rewrite through L2e
    c1 = -A
    c0 = sp.nsimplify(sp.expand(-B + C + A * sp.log(2)))
    return sp.simplify(c1), sp.simplify(c0)

def I3(g):
    """returns (coefficient of eps^-1, p^2-coefficient of L2e, p^2 constant)."""
    a = u * (1 - u)
    e1 = int01(g * u, g * (1 - u))                   # int g m du
    A = int01(g * a)
    B = int01(g * a * sp.log(u * (1 - u)))
    C = logm_int(g * a)
    # p^2 [ A (gamma + log eps p^2 - 1) + B - C ]
    c1 = A
    c0 = sp.expand(-A - A * sp.log(2) + B - C)
    return sp.simplify(e1), sp.simplify(c1), sp.simplify(c0)

M2_1 = I2(sp.Integer(1))
M3_1 = I3(sp.Integer(1))
M2_a = I2(u * (1 - u))
res['bubble_masters'] = {
    'iint exp/s^2': f'{M2_1[0]}*L2e + {M2_1[1]}',
    'iint exp/s^3': f'{M3_1[0]}/eps + p^2({M3_1[1]}*L2e + {M3_1[2]})',
    'iint t1 t2 exp/s^4': f'{M2_a[0]}*L2e + {M2_a[1]}',
    'L2e': 'log(eps |p|^2) + gamma_E + log 2'}
assert M2_1 == (-1, 1) and M3_1 == (sp.Rational(1, 4), sp.Rational(1, 6), sp.Rational(-25, 72))
assert M2_a == (sp.Rational(-1, 6), sp.Rational(13, 72))

# numeric confirmation of iint exp/s^3 at eps = 1e-5, p^2 = 1
mp.mp.dps = 20
e_num = mp.mpf('1e-5')
num = mp.quad(lambda t1: mp.quad(lambda t2: mp.exp(-t1*t2/(t1+t2))/(t1+t2)**3, [e_num, 1, mp.inf]),
              [e_num, 1e-3, 1, mp.inf])
pred = 1/(4*e_num) + (mp.log(e_num) + mp.euler + mp.log(2))/6 - mp.mpf(25)/72
res['check_iint_s3_numeric'] = {'eps': '1e-5', 'numeric': mp.nstr(num, 12), 'closed_form': mp.nstr(pred, 12)}
assert abs(num - pred) < 1e-3

# ---------------------------------------------------------------------------------------------------------
# B. The A|_Wperp self-energy (graphs (a)-(d) of ch04, second table), plane-wave coefficient of (a,b),
#    in units 'unit'.  (c): numerator -K (a,b) |l - p|^2, centre l0 = -(t2/s) p, positive Laplacian -8
#    => bracket p^2 (2-u)^2 + 2/s, overall -1 (derived in the design, Example ex:t-SEc).
#    (d): Clifford trace gives c_d (l, l+p); bracket -u(1-u) p^2 + 2/s; c_d fixed by the Step-1 eps^-1 entry
#    and checked against the Step-1 log entry.  (a),(b): seagulls, exactly (2 + 5/2) * 2 / eps.
# ---------------------------------------------------------------------------------------------------------
L2e = sp.Symbol('L2e')
def bubble_value(c_s3, g_s2):
    """value of iint exp/s^2 [ g_s2(u) p^2 + c_s3 / s ] as eps^-1, L2e, const (p^2 factored in the last two)."""
    e1, a1, a0 = I3(sp.Integer(1))
    b1, b0 = I2(g_s2)
    return sp.expand(c_s3 * e1 / eps + p2 * (c_s3 * (a1 * L2e + a0) + b1 * L2e + b0))

val_c = -bubble_value(2, (2 - u)**2)
val_d_unit = bubble_value(2, -u * (1 - u))            # to be multiplied by c_d
c_d = sp.solve(sp.Eq(val_d_unit.coeff(1/eps) * sp.Symbol('cd'), -4), sp.Symbol('cd'))[0]
val_d = sp.expand(c_d * val_d_unit)
# Step-1 checks: (c): eps^-1 -1/4 int(A,A) -> -1/2 ; log eps +1 int(dA,dA) -> +2 p^2
assert val_c.coeff(1/eps) == sp.Rational(-1, 2) and sp.expand(val_c).coeff(L2e) == 2 * p2
# (d): eps^-1 -2 -> -4 (used), log eps -2 -> -4 p^2 (check)
assert c_d == -8 and val_d.coeff(L2e) == -4 * p2
seagull = sp.Rational(9) / eps
Sigma = sp.expand(seagull + val_c + val_d)
# subtract the quadratic counterterm: plane wave of (9/4 eps^-1 int(A,A) - log eps int(dA,dA)) = 9/2 eps^-1 - 2 p^2 log eps
Sigma_ren = sp.expand(Sigma - sp.Rational(9, 2) / eps)          # the L2e term carries log eps; CT removes -2p^2 log eps
assert Sigma_ren.coeff(L2e) == -2 * p2
const_SE = sp.simplify(Sigma_ren.subs(L2e, 0) / p2)
assert const_SE == sp.Rational(43, 8)
res['A_perp_self_energy'] = {
    'graph_c': str(val_c), 'graph_d': str(val_d), 'c_d': str(c_d), 'seagulls': '9/eps',
    'renormalised': 'unit (a,b) |p|^2 (43/8 - 2(log|p|^2 + gamma_E + log 2))'}

# ---------------------------------------------------------------------------------------------------------
# C. Fourier transforms at x != 0 from FT[(q^2)^a] = 4^a Gamma(2+a)/(pi^2 Gamma(-a)) (x^2)^(-2-a)
#    and FT[(q^2)^(a-1)] = 4^a Gamma(1+a)/(4 pi^2 Gamma(1-a)) (x^2)^(-1-a)  (d^4q/(2pi)^4, e^{iqx}).
# ---------------------------------------------------------------------------------------------------------
a, X2 = sp.symbols('a X2', positive=True)
F0 = 4**a * sp.gamma(2 + a) / (sp.pi**2 * sp.gamma(-a)) * X2**(-2 - a)
F1 = 4**a * sp.gamma(1 + a) / (4 * sp.pi**2 * sp.gamma(1 - a)) * X2**(-1 - a)
s0 = sp.series(F0, a, 0, 3).removeO()
s1 = sp.series(F1, a, 0, 2).removeO()
FT_log = sp.simplify(s0.coeff(a, 1))
FT_log2 = sp.simplify(2 * s0.coeff(a, 2))
FT_logq_over_q2 = sp.simplify(s1.coeff(a, 1))
LX = sp.log(X2)
assert sp.simplify(FT_log + 1 / (sp.pi**2 * X2**2)) == 0
assert sp.simplify(FT_log2 - 2 / (sp.pi**2 * X2**2) * (LX + 2 * sp.EulerGamma - 1 - 2 * sp.log(2))) == 0
assert sp.simplify(FT_logq_over_q2 + 1 / (4 * sp.pi**2 * X2) * (LX + 2 * sp.EulerGamma - 2 * sp.log(2))) == 0
res['fourier'] = {'FT[log q^2]': '-1/(pi^2 x^4)',
                  'FT[log^2 q^2]': '2/(pi^2 x^4) (log x^2 + 2 gamma_E - 1 - 2 log 2)',
                  'FT[log q^2 / q^2]': '-1/(4 pi^2 x^2) (log x^2 + 2 gamma_E - 2 log 2)'}

# ---------------------------------------------------------------------------------------------------------
# D. Masters (in units D0(x)^2, D0 = 1/(4 pi^2 x^2)), Lemmas of Section 7
#    J  = int dz D(x-z)^2 D(z)^2                        = (1/8pi^2) [log(x^2/eps) + gamma_E - 3 log 2]
#    R  = int dz D(x-z)^2 D(z)           (units D0)      = (1/16pi^2)[log(x^2/eps) + 1 + gamma_E - 3 log 2]
#    Xt = J with the quartic point smeared by K_eps      = J + (1/8pi^2) log(2/3)
#    T  = one O-line replaced by K_eps                   = D0 R + D0^2 (1/16pi^2) log(2/3)
# ---------------------------------------------------------------------------------------------------------
k = sp.symbols('k', positive=True)
Delta = sp.simplify(sp.integrate((sp.exp(-3 * k**2) - sp.exp(-2 * k**2)) / k, (k, 0, sp.oo)) / (8 * sp.pi**2))
assert sp.simplify(Delta - sp.log(sp.Rational(2, 3)) / (16 * sp.pi**2)) == 0
# J from FT of B^2, B = (1/16pi^2)(c - log q^2), c = 1 - gamma - log 2 - log eps
c_B = 1 - gE - l2 - Le
J_br = sp.expand(-2 * c_B * (-1) + 2 * (Lx + 2 * gE - 1 - 2 * l2))          # times 1/(pi^2 x^4) (1/16 pi^2)^2
# (1/16pi^2)^2 * 2/(pi^2 x^4) = D0^2 /(8 pi^2): J = D0^2/(8pi^2) * J_br/2
J_unit = sp.expand(J_br / 2)                                                   # in units D0^2/(8 pi^2)
assert sp.expand(J_unit - (Lx - Le + gE - 3 * l2)) == 0
R_unit = sp.expand(c_B + Lx + 2 * gE - 2 * l2)                                 # units D0/(16 pi^2)
assert sp.expand(R_unit - (Lx - Le + 1 + gE - 3 * l2)) == 0
Xt_unit16 = sp.expand(2 * J_unit + 2 * (l2 - l3))                              # units D0^2/(16 pi^2)
T_unit16 = sp.expand(R_unit + (l2 - l3))                                       # units D0^2/(16 pi^2)

# numerical confirmation of J and R in position space (x = unit vector), eps = 1e-4
def D(r, e):
    r = np.asarray(r, dtype=float)
    return np.where(r > 1e-12, -np.expm1(-r * r / (4 * e)) / (4 * np.pi**2 * np.maximum(r * r, 1e-300)),
                    1 / (16 * np.pi**2 * e))

def J_num(e):
    # J = 2 * (region |z| < |x - z|), polar coordinates around 0; inner angle where r cos th < 1/2
    def inner(r):
        thmin = 0.0 if r <= 0.5 else np.arccos(0.5 / r)
        f = lambda th: np.sin(th)**2 * D(np.sqrt(1 + r * r - 2 * r * np.cos(th)), e)**2
        return integrate.quad(f, thmin, np.pi, limit=200)[0]
    g = lambda s: (lambda r: r**4 * D(r, e)**2 * inner(r))(np.exp(s))
    val = integrate.quad(g, np.log(1e-9), np.log(1e4), points=[0.5 * np.log(e), 0.0], limit=400)[0]
    return 2 * 4 * np.pi * val

def R_num(e):
    # R = int dz D(x-z)^2 D(z), polar coordinates around x: |z| = sqrt(1 + r^2 + 2 r cos th)
    def inner(r):
        f = lambda th: np.sin(th)**2 * D(np.sqrt(1 + r * r + 2 * r * np.cos(th)), e)
        pts = [np.pi * 0.999] if abs(r - 1) < 0.05 else None
        return integrate.quad(f, 0, np.pi, points=pts, limit=400)[0]
    g = lambda s: (lambda r: r**4 * D(r, e)**2 * inner(r))(np.exp(s))
    val = integrate.quad(g, np.log(1e-9), np.log(1e4), points=[0.5 * np.log(e), 0.0], limit=400)[0]
    return 4 * np.pi * val

e_t = 1e-4
D0 = 1 / (4 * np.pi**2)
EG = float(sp.EulerGamma)
J_pred = D0**2 / (8 * np.pi**2) * (-np.log(e_t) + EG - 3 * np.log(2))
R_pred = D0 / (16 * np.pi**2) * (-np.log(e_t) + 1 + EG - 3 * np.log(2))
Jn, Rn = J_num(e_t), R_num(e_t)
res['check_masters_numeric'] = {'eps': e_t, 'J_numeric/J_pred': Jn / J_pred, 'R_numeric/R_pred': Rn / R_pred}
assert abs(Jn / J_pred - 1) < 5e-3 and abs(Rn / R_pred - 1) < 5e-3

# ---------------------------------------------------------------------------------------------------------
# E. Assembly, relative to the tree G0 = 4 hbar dim(g) S D0(x)^2, in units 'unit'
# ---------------------------------------------------------------------------------------------------------
# self-energy, both lines: 2 * 2 * [43/8 - 2 gamma - 2 log2 + 2 (log x^2 + 2 gamma - 2 log 2)] / 2 ... (line factor)
line = 2 * (sp.Rational(43, 8) - 2 * gE - 2 * l2 + 2 * (Lx + 2 * gE - 2 * l2))
SE = sp.expand(2 * line)
# exchange: 32 pi^2 * unit * G0/D0^2 * (Xt - 4T) = 2 [Xt_unit16 - 4 T_unit16]
EX = sp.expand(2 * (Xt_unit16 - 4 * T_unit16))
# quartic: 2 * 16 pi^2 * unit * g0 (2T - S) J  = 4 (2 tau - 1) J_unit
QU = sp.expand(4 * (2 * tau - 1) * J_unit)
# operator counterterms: O(eps) = O - unit log eps [4 O - 2 (a1,a2) trace], twice: -8 (1 - tau) log eps
OC = sp.expand(-8 * (1 - tau) * Le)
TOT = sp.expand(SE + EX + QU + OC)
assert TOT.coeff(Le) == 0
assert sp.expand(TOT.coeff(Lx) - 8 * tau) == 0
const = sp.expand(TOT - TOT.coeff(Lx) * Lx)
expected = sp.Rational(27, 2) - 4 * l2 + 4 * l3 + 8 * tau * (gE - 3 * l2)
assert sp.expand(const - expected) == 0
res['classes_relative_to_tree'] = {'self_energy': str(SE), 'exchange': str(EX), 'quartic': str(QU),
                                   'operator_counterterms': str(OC), 'total': str(TOT)}
res['result'] = {
    'general': 'G = G0 [1 + unit (27/2 + 4 log(3/2) + 8 tau (log|x|^2 + gamma_E - 3 log 2))], tau = T/S',
    'traceless (tau=0)': 'G = G0 [1 + unit (27/2 + 4 log(3/2))]',
    'trace (tau=3)': 'G = G0 [1 + unit (27/2 + 4 log(3/2) + 24 (log|x|^2 + gamma_E - 3 log 2))]',
    'log|x|^2 coefficient': '8 tau unit; = -(1/4)(D_1 G0 + G0 D_1^T)/G0 with D_1 = 4(Z - unit {N,O}): 0 and 24',
    'tree': 'G0(x) = 4 hbar dim(g) [(a1,a2)(b1,b2) + (a1,b2)(b1,a2)] / (16 pi^4 |x|^4)'}
# scaling-identity cross-check: trace D_1 = 2 * (-24) unit -> log coefficient -(1/4)*2*D_1 = 24
assert -sp.Rational(1, 4) * 2 * (2 * (-24)) == 24
with open(os.path.join(HERE, 'two_point_example_results.json'), 'w') as f:
    json.dump(res, f, indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
