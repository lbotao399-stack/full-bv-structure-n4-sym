"""P8 example, direct numerical evaluation at finite eps (independent of both engines): bubble E=2, numerator |l|^4,
external momentum p (|p|^2 = 1), 16pi^2 I(eps) = int_{t1,t2>=eps} T^-2 [6 T^-2 + 6 c^2 |p|^2 T^-1 + c^4 |p|^4]
exp(-t1 t2 |p|^2/T), T = t1+t2, c = t2/T  (Gaussian loop integral done by hand).
r(eps) = eps dI/deps = -[f-faces];  fit r = -2 a2 eps^-2 - a1 eps^-1 + b; predicted per edge a2 = 1/4, a1 = |p|^2/8, b = 0."""
import json, os
import mpmath as mp
mp.mp.dps = 30
p2 = mp.mpf(1)
def f(t1, t2):
    T = t1 + t2; c = t2 / T
    return T ** -2 * (6 * T ** -2 + 6 * c ** 2 * p2 / T + c ** 4 * p2 ** 2) * mp.e ** (-t1 * t2 * p2 / T)
def r(eps):
    g1 = mp.quad(lambda t: f(eps, t), [eps, 1, 10, mp.inf])
    g2 = mp.quad(lambda t: f(t, eps), [eps, 1, 10, mp.inf])
    return -eps * (g1 + g2)
out = []
for k in (3, 4, 5, 6):
    eps = mp.mpf(10) ** (-k)
    rv = r(eps)
    # remove the predicted singular terms; the remainder should be b + O(eps log eps) with b = 0
    rem = rv + 2 * mp.mpf(1) / 4 / eps ** 2 + (p2 / 8) / eps
    out.append(dict(eps=str(eps), r=mp.nstr(rv, 15), remainder_after_predicted_poles=mp.nstr(rem, 6)))
    print(out[-1])
json.dump(out, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync', 'shape_numeric.json'), 'w'), indent=1)
