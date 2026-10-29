"""Refuter #2 (Part C): exact master integrals over the scaled UV-boundary region.

J(a,b,c,n) = int_0^1 dtau1 int_1^oo dtau2 int_1^oo dtau3  tau1^a tau2^b tau3^c (tau1+tau2+tau3)^(-n)

tau1 = special-edge time / eps in [0,1]; tau2, tau3 = P(eps,infty) times / eps in [1,oo).
Convergence: integrand homogeneous of degree a+b+c-n <= -3 in (tau2,tau3) at large times
(2-dimensional region, degree < -2 suffices); tau1 ranges over a compact interval and the
integrand is bounded there (tau2+tau3 >= 2).

Exact evaluation, by hand-derived iterated formulas (no symbolic integrator):
 (i)  tau3: with S = tau1+tau2, int_1^oo tau3^c (S+tau3)^-n = sum_j C(c,j) (-S)^(c-j) (S+1)^(j-n+1)/(n-j-1)
 (ii) tau2: v = tau1+tau2+1 in [tau1+2, oo): polynomial in v times v^-(n-1); every power v^m, m <= -2
 (iii) tau1: x = tau1+2 in [2,3]: polynomial times x^k; x^-1 gives log(3/2).
Result: r0 + r1*log(3/2), r0, r1 rational.  Cross-check: mpmath quadrature.
Own code: imports no engine/theory/qa module.
"""
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb
import mpmath as mp


def padd(p, q):
    r = dict(p)
    for k, v in q.items():
        r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v != 0}


def pmul(p, q):
    r = {}
    for k1, v1 in p.items():
        for k2, v2 in q.items():
            k = (k1[0] + k2[0], k1[1] + k2[1])
            r[k] = r.get(k, 0) + v1 * v2
    return {k: v for k, v in r.items() if v != 0}


def ppow(p, e):
    r = {(0, 0): Fr(1)}
    for _ in range(e):
        r = pmul(r, p)
    return r

# polynomials in two variables (t1, v) as dict {(i,j): coeff}


@lru_cache(maxsize=None)
def J(a, b, c, n):
    """returns (rational, coefficient of log(3/2))"""
    assert a + b + c - n <= -3, (a, b, c, n)
    # step (i): inner tau3 integral = sum_j C(c,j)(-1)^(c-j)/(n-j-1) * S^(c-j) (S+1)^(j-n+1)
    # S = v - 1, S+1 = v.  => sum_j coef_j (v-1)^(c-j) v^(j-n+1)
    # step (ii): tau2 = v - t1 - 1 ; integrand tau2^b * [..]; poly in (t1, v) times v^(j-n+1)
    tau2 = {(0, 1): Fr(1), (1, 0): Fr(-1), (0, 0): Fr(-1)}  # v - t1 - 1
    vm1 = {(0, 1): Fr(1), (0, 0): Fr(-1)}                  # v - 1
    total = {}  # dict over (t1 power, v power) with v power possibly negative
    tb = pmul(ppow(tau2, b), {(a, 0): Fr(1)})  # tau1^a tau2^b
    for j in range(c + 1):
        coef = Fr(comb(c, j) * (-1) ** (c - j), n - j - 1)
        poly = pmul(tb, ppow(vm1, c - j))
        shift = j - n + 1
        for (i, m), val in poly.items():
            key = (i, m + shift)
            total[key] = total.get(key, 0) + coef * val
    # integrate v from t1+2 to oo: v^m -> -(t1+2)^(m+1)/(m+1), need m <= -2
    # result: sum over (i, m) of val * t1^i * (-(t1+2)^(m+1)/(m+1))
    # step (iii): x = t1 + 2, t1 = x - 2, integrate x in [2,3]
    rat = Fr(0)
    lg = Fr(0)
    for (i, m), val in total.items():
        if val == 0:
            continue
        assert m <= -2, (a, b, c, n, i, m, val)
        pref = -val / (m + 1)
        # t1^i (x)^(m+1) = sum_k C(i,k) x^k (-2)^(i-k) x^(m+1)
        for k in range(i + 1):
            cc = pref * comb(i, k) * Fr(-2) ** (i - k)
            e = k + m + 1
            if e == -1:
                lg += cc
            else:
                rat += cc * (Fr(3) ** (e + 1) - Fr(2) ** (e + 1)) / (e + 1) if e + 1 >= 0 else \
                    cc * (Fr(1, 3 ** (-(e + 1))) - Fr(1, 2 ** (-(e + 1)))) / (e + 1)
    return (rat, lg)


def Jnum(a, b, c, n, dps=30):
    mp.mp.dps = dps
    f = lambda x, y, z: x**a * y**b * z**c / (x + y + z)**n
    return mp.quad(f, [0, 1], [1, mp.inf], [1, mp.inf])


def Jfloat(a, b, c, n, dps=30):
    mp.mp.dps = dps
    r, l = J(a, b, c, n)
    return mp.mpf(r.numerator) / r.denominator + mp.mpf(l.numerator) / l.denominator * mp.log(mp.mpf(3) / 2)


if __name__ == '__main__':
    # regression against sympy-integrated values (computed once, see refute log):
    # (0,0,0,3): log(3/2)/2 ; (0,0,0,4): 1/36 ; (1,0,0,4): log(3/2)/6 - 1/18 ;
    # (0,1,1,5): 17/864 + log(3/2)/24 ; (2,0,0,5): -7/216 + log(3/2)/12 ; (0,0,2,6): 17/2592
    ref = {(0, 0, 0, 3): (Fr(0), Fr(1, 2)), (0, 0, 0, 4): (Fr(1, 36), Fr(0)),
           (1, 0, 0, 4): (Fr(-1, 18), Fr(1, 6)), (0, 1, 1, 5): (Fr(17, 864), Fr(1, 24)),
           (2, 0, 0, 5): (Fr(-7, 216), Fr(1, 12)), (0, 0, 2, 6): (Fr(17, 2592), Fr(0))}
    for k, v in ref.items():
        print(k, J(*k), v, J(*k) == v)
    for k in [(1, 0, 1, 5), (0, 2, 1, 6), (1, 1, 1, 6), (0, 0, 3, 6), (2, 1, 0, 7), (1, 2, 0, 6)]:
        print(k, J(*k), mp.nstr(Jfloat(*k, dps=15), 14), mp.nstr(Jnum(*k, dps=12), 12))
