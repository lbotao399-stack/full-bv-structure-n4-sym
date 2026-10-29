"""Referee checks of Thm tl-simplex, Prop tl-shape, Lemma tl-special, Lemma tl-k0 (independent numerics).

All integrals are done directly in Schwinger form with mpmath (no use of theory/tl_simplex.py):
  edge momenta l (edge 0) and l - p (edge 1), Gaussian average of the numerator in closed form.
  <|l'|^2> = 2/tau, <|l'|^4> = 6/tau^2, <(l'.c)^2> = |c|^2/(2 tau) for l' centred, d^4l/(2pi)^4 e^{-tau l'^2} = 1/(16 pi^2 tau^2).
Output: check_simplex_indep.json
"""
import json, os
import mpmath as mp
import sympy as sp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))
P2 = mp.mpf(1)            # |p|^2
out = {}


def gauss_avg_l2(tau, c2):      # <|l|^2> with l = l' + c
    return 2 / tau + c2


def gauss_avg_l4(tau, c2):      # <|l|^4>
    return 6 / tau**2 + 6 * c2 / tau + c2**2


# 1. Lemma tl-special example: N = -|l|^2 on the special edge l, P edge l - p. 16 pi^2 w = -int_0^eps dt int_eps^inf dt'
#    <|l|^2> e^{-t t' p^2/tau}/tau^2,  c = t' p/tau.
def special_example(eps):
    def inner(t):
        f = lambda tp: gauss_avg_l2(t + tp, (tp / (t + tp))**2 * P2) * mp.e**(-t * tp * P2 / (t + tp)) / (t + tp)**2
        return mp.quad(f, [eps, 10 * eps, 1, mp.inf])
    return -mp.quad(inner, [0, eps])


rows = []
for eps in [mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')]:
    v = special_example(eps)
    rows.append(dict(eps=str(eps), value=mp.nstr(v, 15), remainder=mp.nstr(v + 1 / (2 * eps) + mp.mpf(3) / 8 * P2, 8)))
out['special_edge_example_theory'] = dict(claim='16pi^2 w = -1/(2 eps) - 3/8 |p|^2 + o(1)', rows=rows)


# 2. Lemma tl-special with a log-producing numerator: N = -|l|^2 (special) * |l|^2.  Prediction from the lemma.
a = sp.symbols('a', positive=True)            # a = alpha_1 (weight of the P edge l - p), alpha_0 = 1 - a (special)
p2 = sp.Integer(1)
c2 = a**2 * p2
Phi = a * (1 - a) * p2
num0, num1, num2 = -c2**2, -6 * c2, sp.Integer(-6)          # (Delta^j N)(c)/(j! 4^j), N = -|l|^4
F2 = num2
F1 = num1 + num2 * (-Phi)
F0 = num0 + num1 * (-Phi) + num2 * Phi**2 / 2
a0 = 1 - a
pred_m2 = sp.Rational(1, 2) * sp.integrate((a**2 - a0**2) * F2, (a, sp.Rational(1, 2), 1))
pred_m1 = sp.integrate((a - a0) * F1, (a, sp.Rational(1, 2), 1))
pred_fin = -sp.integrate(sp.expand(F0) * sp.log(a0 / a), (a, sp.Rational(1, 2), 1))
pred_fin = sp.nsimplify(sp.simplify(pred_fin))
# also: what the same cycle with an ordinary P edge in place of the special one gives for log eps
log_full = -sp.integrate(F0, (a, 0, 1))


def special_l4(eps):
    def inner(t):
        f = lambda tp: gauss_avg_l4(t + tp, (tp / (t + tp))**2 * P2) * mp.e**(-t * tp * P2 / (t + tp)) / (t + tp)**2
        return mp.quad(f, [eps, 10 * eps, 1, mp.inf])
    return -mp.quad(inner, [0, eps])


rows = []
for eps in [mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')]:
    v = special_l4(eps)
    rem = v - mp.mpf(sp.N(pred_m2, 30)) / eps**2 - mp.mpf(sp.N(pred_m1, 30)) / eps - mp.mpf(sp.N(pred_fin, 30))
    rows.append(dict(eps=str(eps), value=mp.nstr(v, 15), remainder_after_lemma_prediction=mp.nstr(rem, 8)))
out['special_edge_l4'] = dict(pred_eps_m2=str(pred_m2), pred_eps_m1=str(pred_m1), pred_finite=str(pred_fin),
                              log_coefficient_of_same_cycle_with_P_edge=str(sp.simplify(log_full)), rows=rows)


# 3. Prop tl-shape example: E=2, N=|l|^4; per-edge vs total-time eps^-2, eps^-1 coefficients (claims 1/4, 1/8 ; 3, 1)
def total_time(eps):
    f = lambda al, tau: tau * gauss_avg_l4(tau, (al**2) * P2) * mp.e**(-tau * al * (1 - al) * P2) / tau**2
    return mp.quad(lambda al: mp.quad(lambda tau: f(al, tau), [eps, 1, mp.inf]), [0, 0.5, 1])


def per_edge(eps):
    f = lambda t1, t2: gauss_avg_l4(t1 + t2, (t2 / (t1 + t2))**2 * P2) * mp.e**(-t1 * t2 * P2 / (t1 + t2)) / (t1 + t2)**2
    return mp.quad(lambda t1: mp.quad(lambda t2: f(t1, t2), [eps, 10 * eps, 1, mp.inf]), [eps, 10 * eps, 1, mp.inf])


rows = []
for eps in [mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')]:
    vt = total_time(eps)
    vp = per_edge(eps)
    rows.append(dict(eps=str(eps), total_time_minus_claim=mp.nstr(vt - 3 / eps**2 - P2 / eps, 10),
                     per_edge_minus_claim=mp.nstr(vp - 1 / (4 * eps**2) - P2 / (8 * eps), 10)))
out['shape_example'] = dict(claim='per-edge 1/4 eps^-2 + |p|^2/8 eps^-1; total-time 3 eps^-2 + |p|^2 eps^-1; log 0',
                            note='remainders must tend to a constant (no log eps): compare successive rows', rows=rows)

# 4. Lemma tl-k0: derivatives of S_eps at 0 vs the pairing formula
x = sp.symbols('x0:4')
eps_s = sp.symbols('eps', positive=True)
r2 = sum(xi**2 for xi in x)
S = sum((-1)**(m + 1) / (16 * sp.pi**2 * sp.factorial(m) * eps_s**m) * (r2 / 4)**(m - 1) for m in range(1, 6))
checks = []
for gamma in [(0, 0, 0, 0), (2, 0, 0, 0), (1, 1, 0, 0), (4, 0, 0, 0), (2, 2, 0, 0), (2, 1, 1, 0), (6, 0, 0, 0), (2, 2, 2, 0)]:
    expr = S
    for i, n in enumerate(gamma):
        expr = sp.diff(expr, x[i], n)
    val = sp.simplify(expr.subs({xi: 0 for xi in x}))
    n = sum(gamma)
    if n % 2:
        pred = 0
    else:
        m = n // 2 + 1
        # number of perfect pairings of the multiset gamma with all pairs equal indices = prod (g_i - 1)!!
        npair = 1
        for g in gamma:
            if g % 2:
                npair = 0
            else:
                npair *= sp.factorial2(g - 1) if g > 0 else 1
        pred = (-1)**(m + 1) / (16 * sp.pi**2 * m * 2**(m - 1) * eps_s**m) * npair
    checks.append(dict(gamma=gamma, value=str(val), formula=str(sp.simplify(pred)), equal=bool(sp.simplify(val - pred) == 0)))
out['k0_formula'] = checks

json.dump(out, open(os.path.join(HERE, 'check_simplex_indep.json'), 'w'), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
