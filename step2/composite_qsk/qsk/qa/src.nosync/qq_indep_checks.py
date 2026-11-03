"""QA: independent analytic checks (sympy, light).
(1) n = 1 defect cycle (Example qs-cubic / Lemma c10-special): special edge t_s in [0, eps] with numerator -|l+q|^2,
    ordinary edge t in [eps, T] with numerator 1, Gaussian loop integral: 16 pi^2 x value = -1/(2 eps) - (3/8)|q|^2
    + O(eps |q|^4) at T = infinity, and the T-dependence is O(eps) at fixed T.
(2) Large-distance coefficient of the composite A|_V--c line along K_b (Prop c8-Tsafe, Thm c8-K): angular average of
    the integrand of int d^4y G(x1-y) (L_v A)_m d_m G(y-x2) at |y| -> infinity, scalar kernels G = 1/(4 pi^2 |y|^2):
    coefficient of log T of the scalar integral = -3/(16 pi^2) b  (R^2 ~ T).
(3) n = 2 constant-numerator closed form of engine B (F8) at n_t = 2: 16 pi^2 [eps^0] = log 2."""
import sympy as sp

s, t, eps, T, q2, u = sp.symbols('s t epsilon T q2 u', positive=True)
tau = s + t
# (1) integrand after the Gaussian integral: -(16 pi^2)^-1 tau^-2 [2/tau + (t/tau)^2 q2] exp(-(s t/tau) q2)
expo = sp.exp(-(s * t / tau) * q2)
integrand = -(2 / tau ** 3 + (t ** 2 / tau ** 4) * q2) * expo
ser = sp.series(integrand, q2, 0, 2).removeO()
c0 = sp.integrate(sp.integrate(ser.coeff(q2, 0), (t, eps, sp.oo)), (s, 0, eps))
c1 = sp.integrate(sp.integrate(sp.simplify(ser.coeff(q2, 1)), (t, eps, sp.oo)), (s, 0, eps))
print('(1) T=oo: q^0 term:', sp.simplify(c0), ' q^2 term:', sp.simplify(c1))
c0T = sp.integrate(sp.integrate(ser.coeff(q2, 0), (t, eps, T)), (s, 0, eps))
c1T = sp.integrate(sp.integrate(sp.simplify(ser.coeff(q2, 1)), (t, eps, T)), (s, 0, eps))
print('    finite T minus T=oo, q^0:', sp.simplify(c0T - c0), ' series in eps:', sp.series(sp.simplify(c0T - c0), eps, 0, 2))
print('    finite T minus T=oo, q^2:', sp.series(sp.simplify(c1T - c1), eps, 0, 2))
# (2) angular average
y = sp.symbols('y0:4', real=True); b = sp.symbols('b0:4', real=True)
Y2 = sum(v ** 2 for v in y); bdy = sum(bi * yi for bi, yi in zip(b, y))
v = [Y2 * b[k] - 2 * bdy * y[k] for k in range(4)]
G = 1 / (4 * sp.pi ** 2 * Y2)
dG = [sp.diff(G, y[m]) for m in range(4)]
# x1 = x2 = 0 at large y: G1 = G(y), dG2 = d/dy G(y - x2) -> dG(y)
res = []
for k in range(4):
    term = sum(v[n] * dG[n] for n in range(4)) * dG[k] + G * sum(sp.diff(v[k], y[m]) * dG[m] for m in range(4))
    res.append(sp.cancel(sp.together(term * Y2 ** 3)))
# angular average on the unit sphere: <y_i y_j> = delta_ij / 4, <1> = 1 (after factoring |y|^-4)
def avg(expr):
    p = sp.Poly(sp.expand(expr), *y)
    tot = 0
    for mon, c in p.terms():
        if sum(mon) == 0:
            tot += c
        elif sum(mon) == 2 and max(mon) == 2:
            tot += c * sp.Rational(1, 4)          # <y_i^2> = 1/4 on the unit sphere of R^4
        elif sum(mon) == 2:
            tot += 0
        else:
            # homogeneous degree 2 only expected after multiplying by |y|^4
            raise ValueError(mon)
    return sp.simplify(tot)
avgs = []
for k in range(4):
    e = sp.expand(res[k])
    # numerator homogeneous of degree 2 in y: unit-sphere average = coefficient of |y|^-4
    avgs.append(sp.simplify(avg(e)))
print('(2) |y|^4 x integrand, unit-sphere average, component k:', avgs)
coef_logR = [sp.simplify(a * 2 * sp.pi ** 2) for a in avgs]       # int d^4y |y|^-4 = 2 pi^2 log R
print('    coefficient of log R:', coef_logR, '; of log T (R^2 ~ T):', [sp.simplify(c / 2) for c in coef_logR])
# (3)
c3 = sp.integrate(sp.integrate(1 / tau ** 2, (t, eps, sp.oo)), (s, 0, eps))
print('(3) n_t = 2 constant numerator, 16 pi^2 x value:', sp.simplify(c3))
