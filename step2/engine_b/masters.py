"""Engine B, Step 2: exact Schwinger master integrals of the UV-boundary graphs with one special edge.

Special edge: Schwinger time t_1 in [0, eps]; every other edge: t_e in [eps, infty).  After t = eps*u:
  J2(a1, a2; m)     = int_0^1 du1 int_1^oo du2               u1^a1 u2^a2        (u1+u2)^-m
  J3(a1, a2, a3; m) = int_0^1 du1 int_1^oo du2 int_1^oo du3  u1^a1 u2^a2 u3^a3  (u1+u2+u3)^-m
Convergence: J2 needs m > a2 + 1, J3 needs m > a2 + a3 + 2 (large-u power counting; near u = 0 nothing happens since the
region excludes u2, u3 < 1).  Exact values lie in Q + Q log 2 (J2) and Q + Q log(3/2) (J3); returned as
Num = {'1': Fraction, 'log2': Fraction, 'log3': Fraction}.

Method (no special functions): integrate u3 (resp. u2) by the binomial expansion of u^a around (S+u)^-m (rational in S),
then u2 by w = 1+u1+u2 (all powers of w are <= -2, so no logarithm), then u1 by v = 2+u1 (resp. 1+u1): the only
logarithm is int_2^3 dv/v = log(3/2) (resp. int_1^2 dv/v = log 2)."""
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb

def nadd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, Fr(0)) + s * v
    return {k: v for k, v in out.items() if v != 0}

def nscale(a, s):
    return {k: v * s for k, v in a.items() if v * s != 0}

def nzero():
    return {}

def nfloat(a, mp=None):
    import mpmath
    l2 = mpmath.log(2); l3 = mpmath.log(3)
    return sum((mpmath.mpf(v.numerator) / v.denominator) * {'1': 1, 'log2': l2, 'log3': l3}[k] for k, v in a.items())

def nstr(a):
    if not a:
        return '0'
    parts = []
    for k in ('1', 'log2', 'log3'):
        if k in a:
            v = a[k]
            parts.append(str(v) if k == '1' else f'({v})*{"log(2)" if k == "log2" else "log(3)"}')
    return ' + '.join(parts)

# polynomials in one variable: dict power -> Fraction
def pmul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i + j] = out.get(i + j, Fr(0)) + x * y
    return {k: v for k, v in out.items() if v != 0}

def ppow_lin(c0, c1, n):
    """(c0 + c1 X)^n as a polynomial in X."""
    out = {k: Fr(comb(n, k)) * Fr(c0) ** (n - k) * Fr(c1) ** k for k in range(n + 1)}
    return {k: v for k, v in out.items() if v != 0}

def int_v_power(e, lo, hi, logkey):
    """int_lo^hi v^e dv, exact (lo, hi positive integers)."""
    if e == -1:
        return logkey
    return {'1': (Fr(hi) ** (e + 1) - Fr(lo) ** (e + 1)) / (e + 1)}

@lru_cache(maxsize=None)
def J3(a1, a2, a3, m):
    assert m > a2 + a3 + 2, (a1, a2, a3, m)
    res = {}
    # u3 integral: sum_i C(a3,i) (-S)^{a3-i} (S+1)^{i-m+1}/(m-i-1),  S = u1+u2
    for i in range(a3 + 1):
        c_i = Fr(comb(a3, i) * (-1) ** (a3 - i), m - i - 1)
        # integrand in u2: u2^a2 S^{a3-i} (1+S)^{i-m+1}; w = 1+S: u2 = w - 1 - u1, S = w - 1
        # expand (w - 1 - u1)^a2 (w - 1)^{a3-i} in powers of w with coefficients polynomial in u1
        # represent bivariate as dict (pw, pu) -> coeff
        A = {}
        for k in range(a2 + 1):        # (w + (-1 - u1))^a2 = sum_k C(a2,k) w^k (-1-u1)^{a2-k}
            pol = ppow_lin(-1, -1, a2 - k)   # (-1 - u1)^{a2-k}
            for pu, cu in pol.items():
                A[(k, pu)] = A.get((k, pu), Fr(0)) + comb(a2, k) * cu
        B = {}
        n2 = a3 - i
        for k in range(n2 + 1):        # (w - 1)^n2
            B[k] = Fr(comb(n2, k) * (-1) ** (n2 - k))
        AB = {}
        for (pw, pu), c in A.items():
            for k, cb in B.items():
                AB[(pw + k, pu)] = AB.get((pw + k, pu), Fr(0)) + c * cb
        for (pw, pu), c in AB.items():
            if c == 0:
                continue
            e = pw + i - m + 1                     # power of w
            assert e <= -2, (a1, a2, a3, m, e)
            # int_{2+u1}^oo w^e dw = (2+u1)^{e+1} / (-(e+1))
            cc = c_i * c / (-(e + 1))
            # u1 integral: int_0^1 u1^{a1+pu} (2+u1)^{e+1} du1, v = 2+u1
            pol = ppow_lin(-2, 1, a1 + pu)          # (v-2)^{a1+pu}
            for pv, cv in pol.items():
                res = nadd(res, nscale(int_v_power(pv + e + 1, 2, 3, {'log3': Fr(1), 'log2': Fr(-1)}), cc * cv))
    return res

@lru_cache(maxsize=None)
def J2(a1, a2, m):
    assert m > a2 + 1, (a1, a2, m)
    res = {}
    # u2 integral: sum_i C(a2,i) (-u1)^{a2-i} (1+u1)^{i-m+1}/(m-i-1)
    for i in range(a2 + 1):
        c_i = Fr(comb(a2, i) * (-1) ** (a2 - i), m - i - 1)
        # u1 integral: int_0^1 u1^{a1+a2-i} (1+u1)^{i-m+1} du1, v = 1+u1
        pol = ppow_lin(-1, 1, a1 + a2 - i)
        for pv, cv in pol.items():
            res = nadd(res, nscale(int_v_power(pv + i - m + 1, 1, 2, {'log2': Fr(1)}), c_i * cv))
    return res

def check_numeric(nsample=8):
    """second, independent evaluation: direct numerical quadrature (scipy, adaptive) of the defining integrals,
    after u -> 1/x on the infinite ranges; agreement to ~1e-9 relative."""
    import random
    from scipy import integrate
    random.seed(5)
    out = []
    for _ in range(nsample):
        a1, a2, a3 = random.randint(0, 3), random.randint(0, 3), random.randint(0, 3)
        m = a1 + a2 + a3 + random.choice([3, 4])
        def f(x3, x2, u1):
            u2 = 1 / x2; u3 = 1 / x3
            return u1 ** a1 * u2 ** a2 * u3 ** a3 / (u1 + u2 + u3) ** m / (x2 * x2 * x3 * x3)
        num, err = integrate.tplquad(f, 0, 1, 0, 1, 0, 1, epsabs=1e-13, epsrel=1e-11)
        ex = float(nfloat(J3(a1, a2, a3, m)))
        out.append(('J3', (a1, a2, a3, m), nstr(J3(a1, a2, a3, m)), ex, num, abs(num - ex) / abs(ex)))
    for _ in range(nsample):
        a1, a2 = random.randint(0, 4), random.randint(0, 4)
        m = a1 + a2 + random.choice([2, 3])
        def g(x2, u1):
            u2 = 1 / x2
            return u1 ** a1 * u2 ** a2 / (u1 + u2) ** m / (x2 * x2)
        num, err = integrate.dblquad(g, 0, 1, 0, 1, epsabs=1e-14, epsrel=1e-12)
        ex = float(nfloat(J2(a1, a2, m)))
        out.append(('J2', (a1, a2, m), nstr(J2(a1, a2, m)), ex, num, abs(num - ex) / abs(ex)))
    return out

if __name__ == '__main__':
    for r in check_numeric():
        print(r)
