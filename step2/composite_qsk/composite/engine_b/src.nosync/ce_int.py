"""Route B exact integrals for a cycle with one special edge (Schwinger time in [0, eps]) and n_t - 1 ordinary
Schwinger parameters t_f >= c_f eps (c_f = 1, or 2 for the single-segment eps d/d eps term of the square kernel).

After the Gaussian loop integral, the coefficient of eps^{-m} (m >= 0; m = 0 the finite part) of the cycle is
    (16 pi^2)^{-1} int_0^1 dv int_{u_f >= c_f} du  F_m(alpha) tau^{-n_t-m},   t = eps (v, u),  tau = v + sum u,
    alpha = (v, u)/tau,   F_m = sum_j num_j (-Phi)^{j+2-n_t-m}/(j+2-n_t-m)!.
Route B evaluates it with Feynman parameters for the ordinary edges, u_f = U beta_f (beta on the simplex),
du = U^{n_t-2} dU dbeta, u_f >= c_f <=> U >= L(beta) = max_f c_f/beta_f, and the special time v integrated together
with U in closed form.  For a monomial alpha_sigma^a prod_f alpha_f^{n_f}:
    I_m(a,b;L) = int_0^1 dv int_L^oo dU v^a U^b (U+v)^{-(a+b+2+m)},  b = |n| + n_t - 2,
with U = v w, x = w/(1+w), y_L = L/(1+L):
    m = 0 : I_0 = int_{y_L}^1 B(y) dy/(y(1-y)),          B(y) = int_y^1 x^b (1-x)^a dx,
    m >= 1: I_m = L^{-m} int_{y_L}^1 y^{m-1}(1-y)^{-m-1} B_m(y) dy,  B_m(y) = int_y^1 x^b (1-x)^{a+m} dx,
polynomials in y_L plus -B(0) log y_L for m = 0.  The beta-simplex is split by the edge f* attaining L; on that
region L = c_{f*}/s, s = beta_{f*} in [0, 1/Gamma], Gamma = sum_f c_f / c_{f*}, and the other beta's are integrated
exactly (Dirichlet) over the shifted simplex; the last s-integral is done exactly with w = c_{f*} + s.
Values: LogNum = rational + sum_p rational * log p (p prime)."""
from fractions import Fraction as Fr
from functools import lru_cache
from math import comb, factorial
import itertools


def _factor(n):
    n = int(n)
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


class LogNum:
    __slots__ = ('r', 'l')

    def __init__(self, r=0, l=None):
        self.r = Fr(r)
        self.l = {p: Fr(v) for p, v in (l or {}).items() if v != 0}

    @staticmethod
    def log(q):
        """log of a positive rational q."""
        q = Fr(q)
        assert q > 0
        d = {}
        for p, e in _factor(q.numerator).items():
            d[p] = d.get(p, 0) + e
        for p, e in _factor(q.denominator).items():
            d[p] = d.get(p, 0) - e
        return LogNum(0, d)

    def __add__(self, o):
        if not isinstance(o, LogNum):
            o = LogNum(o)
        l = dict(self.l)
        for p, v in o.l.items():
            l[p] = l.get(p, 0) + v
        return LogNum(self.r + o.r, l)

    __radd__ = __add__

    def __neg__(self):
        return LogNum(-self.r, {p: -v for p, v in self.l.items()})

    def __sub__(self, o):
        return self + (-o if isinstance(o, LogNum) else LogNum(-Fr(o)))

    def __mul__(self, c):
        assert not isinstance(c, LogNum)
        c = Fr(c)
        return LogNum(self.r * c, {p: v * c for p, v in self.l.items()})

    __rmul__ = __mul__

    def is_zero(self):
        return self.r == 0 and not self.l

    def __float__(self):
        import math
        return float(self.r) + sum(float(v) * math.log(p) for p, v in self.l.items())

    def items(self):
        """basis -> Fraction: '1', 'log2', 'log3', ..."""
        out = {}
        if self.r != 0:
            out['1'] = self.r
        for p, v in sorted(self.l.items()):
            out[f'log{p}'] = v
        return out

    def __repr__(self):
        s = str(self.r)
        for p, v in sorted(self.l.items()):
            s += f' + ({v})*log({p})'
        return s


# ------------------------------------------------------------------ polynomials in one variable (Fraction lists)
def padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def pmul(p, q):
    out = [Fr(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def ppow(p, n):
    out = [Fr(1)]
    for _ in range(n):
        out = pmul(out, p)
    return out


def pint(p):
    """antiderivative vanishing at 0."""
    return [Fr(0)] + [Fr(c) / (i + 1) for i, c in enumerate(p)]


def peval(p, x):
    acc = Fr(0)
    for c in reversed(p):
        acc = acc * x + c
    return acc


def pdivlin1(p):
    """p(y)/(1-y) for p(1) = 0 (exact)."""
    # p = (1-y) q ; synthetic division by (y - 1): p = (y-1)(-q)
    n = len(p) - 1
    qq = [Fr(0)] * n
    rem = Fr(0)
    for i in range(n, 0, -1):
        rem = p[i] + rem
        qq[i - 1] = rem
    assert p[0] + rem == 0, 'p(1) != 0'
    return [-c for c in qq]


@lru_cache(maxsize=None)
def Ipoly(m, a, b):
    """I_m(a,b;L) = (L^{-m}) * [P(1) - P(y_L)] + logc * (-log y_L), returned as (P, logc) with P a polynomial in y
    (P(1) - P(y_L) the polynomial part) and logc the coefficient of -log y_L (m = 0 only)."""
    # B(y) = int_y^1 x^b (1-x)^(a+m) dx = Bfull(1) - Bfull(y)
    inner = pmul([Fr(0)] * b + [Fr(1)], ppow([Fr(1), Fr(-1)], a + m))
    Bi = pint(inner)
    B1 = peval(Bi, Fr(1))
    Bp = [-c for c in Bi]
    Bp[0] += B1                                   # B(y) = B1 - Bi(y)
    C = Bp
    for _ in range(m + 1):
        C = pdivlin1(C)                           # B(y)/(1-y)^{m+1}
    if m == 0:
        c0 = C[0]
        D = C[1:]                                 # (C - C(0))/y
        P = pint(D)
        return tuple(P), c0
    D = pmul([Fr(0)] * (m - 1) + [Fr(1)], C)      # y^{m-1} C
    P = pint(D)
    return tuple(P), Fr(0)


def Inum(m, a, b, L):
    """numerical I_m(a,b;L) for tests."""
    import math
    P, c0 = Ipoly(m, a, b)
    P = list(P)
    y = L / (1 + L)
    val = float(peval(P, Fr(1))) - float(peval(P, Fr(y)) if isinstance(y, Fr) else sum(float(c) * y ** i for i, c in enumerate(P)))
    val *= L ** (-m)
    if m == 0:
        val += float(c0) * (-math.log(y))
    return val


# ------------------------------------------------------------------ Dirichlet slice
@lru_cache(maxsize=None)
def slice_poly(nexp, gam):
    """int over {beta'_f >= 0, sum beta'_f = R} of prod_f (gam_f s + beta'_f)^{n_f}  (Lebesgue on len-1 coords),
    R = 1 - s*Gamma, Gamma = 1 + sum gam_f, as a polynomial in s (Fraction list).  nexp, gam: tuples over f != f*."""
    d = len(nexp)
    Gam = 1 + sum(gam)
    Rp = [Fr(1), -Fr(Gam)]                       # R(s)
    out = [Fr(0)]
    for ms in itertools.product(*[range(n + 1) for n in nexp]):
        cf = Fr(1)
        spow = 0
        for n, mm, g in zip(nexp, ms, gam):
            cf *= comb(n, mm) * Fr(g) ** (n - mm)
            spow += n - mm
        if cf == 0:
            continue
        M = sum(ms)
        # Dirichlet: R^{M + d - 1} prod m! / (M + d - 1)!
        dc = Fr(1)
        for mm in ms:
            dc *= factorial(mm)
        dc /= factorial(M + d - 1)
        term = pmul([Fr(0)] * spow + [cf * dc], ppow(Rp, M + d - 1))
        out = padd(out, term)
    return tuple(out)


# ------------------------------------------------------------------ s-integrals
def _int_w_powers(coeffs, w0, w1):
    """int_{w0}^{w1} sum_e coeffs[e] w^e dw, e integers (possibly negative)."""
    res = LogNum(0)
    for e, c in coeffs.items():
        if c == 0:
            continue
        if e == -1:
            res = res + LogNum.log(Fr(w1) / Fr(w0)) * c
        else:
            res = res + LogNum((Fr(w1) ** (e + 1) - Fr(w0) ** (e + 1)) / (e + 1) * c)
    return res


@lru_cache(maxsize=None)
def s_integral(spoly, m, a, b, cstar, Gamma):
    """int_0^{1/Gamma} spoly(s) I_m(a,b;cstar/s) ds as LogNum."""
    P, c0 = Ipoly(m, a, b)
    P1 = peval(list(P), Fr(1))
    cs = Fr(cstar)
    w0, w1 = cs, cs + Fr(1) / Fr(Gamma)
    # polynomial part: (s/c)^m [P1 - P(c/w)], s = w - c
    # expand spoly(s) s^m = sum_e g_e (w - c)^e
    sp_ = pmul(list(spoly), [Fr(0)] * m + [Fr(1)])
    wcoef = {}
    for e, g in enumerate(sp_):
        if g == 0:
            continue
        for i in range(e + 1):
            # (w-c)^e = sum_i C(e,i) w^i (-c)^{e-i}
            base = g * comb(e, i) * (-cs) ** (e - i) / cs ** m
            # times [P1 - sum_j P_j c^j w^{-j}]
            wcoef[i] = wcoef.get(i, 0) + base * P1
            for j, pj in enumerate(P):
                if pj:
                    wcoef[i - j] = wcoef.get(i - j, 0) - base * pj * cs ** j
    res = _int_w_powers(wcoef, w0, w1)
    if m == 0 and c0 != 0:
        # c0 * int spoly(s) (-log y_L) ds,  -log y_L = log w - log c
        for e, g in enumerate(spoly):
            if g == 0:
                continue
            for i in range(e + 1):
                cf = g * comb(e, i) * (-cs) ** (e - i) * c0
                # int w^i (log w - log c) dw = [w^{i+1}/(i+1) (log w - 1/(i+1))] - log c [w^{i+1}/(i+1)]
                t1 = LogNum.log(w1) * (w1 ** (i + 1) / (i + 1)) - LogNum.log(w0) * (w0 ** (i + 1) / (i + 1)) \
                    - LogNum((w1 ** (i + 1) - w0 ** (i + 1)) / Fr((i + 1) ** 2))
                t2 = LogNum.log(cs) * ((w1 ** (i + 1) - w0 ** (i + 1)) / (i + 1))
                res = res + (t1 - t2) * cf
    return res


@lru_cache(maxsize=None)
def J(m, a, nvec, cut):
    """coefficient integral of the monomial alpha_sigma^a prod alpha_f^{n_f} (nvec over the ordinary parameters,
    cut = their cutoff multipliers c_f), i.e. int_0^1 dv int_{u_f>=c_f} du alpha^n tau^{-n_t-m}."""
    nt = len(nvec) + 1
    b = sum(nvec) + nt - 2
    res = LogNum(0)
    if len(nvec) == 1:
        # one ordinary parameter: beta = 1, L = c_1, no Feynman simplex
        P, c0 = Ipoly(m, a, b)
        Lc = Fr(cut[0])
        y = Lc / (1 + Lc)
        val = LogNum((peval(list(P), Fr(1)) - peval(list(P), y)) * Lc ** (-m))
        if m == 0 and c0 != 0:
            val = val + LogNum.log(y) * (-c0)
        return val
    for fs in range(len(nvec)):
        cst = Fr(cut[fs])
        others = [f for f in range(len(nvec)) if f != fs]
        gam = tuple(Fr(cut[f]) / cst for f in others)
        Gam = 1 + sum(gam)
        sl = slice_poly(tuple(nvec[f] for f in others), gam)
        spoly = pmul(list(sl), [Fr(0)] * nvec[fs] + [Fr(1)])
        res = res + s_integral(tuple(spoly), m, a, b, cst, Gam)
    return res


# ------------------------------------------------------------------ independent numerical check (mpmath)
def J_numeric(m, a, nvec, cut, dps=20):
    """direct quadrature of int_0^1 dv int_{u_f >= c_f} prod du  v^a prod u_f^{n_f} tau^{-(a+|n|)-n_t-m}
    (alpha-monomial times tau^{-n_t-m}) with u_f = c_f + z_f/(1-z_f)."""
    import mpmath as mp
    mp.mp.dps = dps
    nt = len(nvec) + 1
    deg = a + sum(nvec)

    def f(*xs):
        v = xs[0]
        us = []
        jac = mp.mpf(1)
        for z, c in zip(xs[1:], cut):
            us.append(c + z / (1 - z))
            jac *= 1 / (1 - z) ** 2
        tau = v + sum(us)
        val = v ** a
        for u, n in zip(us, nvec):
            val *= u ** n
        return val * tau ** (-(deg + nt + m)) * jac
    return mp.quad(f, *([[0, 1]] * nt))


if __name__ == '__main__':
    import random
    random.seed(1)
    # I_m against direct 2-d quadrature
    import mpmath as mp
    for (m, a, b, Lv) in [(0, 0, 1, 2.0), (0, 2, 3, 3.0), (1, 1, 2, 2.5), (2, 0, 4, 1.5)]:
        num = mp.quad(lambda v, U: v ** a * U ** b * (U + v) ** (-(a + b + 2 + m)), [0, 1], [Lv, mp.inf])
        print('I', m, a, b, Lv, Inum(m, a, b, Lv), float(num))
    for (m, a, nvec, cut) in [(0, 1, (2,), (2,)), (1, 0, (3,), (2,)), (0, 0, (0, 0), (1, 1)), (0, 1, (2, 0), (1, 1)), (1, 0, (1, 1), (1, 1)),
                              (0, 0, (0, 0, 0), (1, 1, 1)), (0, 2, (1, 0, 1), (1, 1, 1)), (1, 1, (0, 2, 0), (1, 1, 1)),
                              (0, 1, (1, 0), (2, 1)), (1, 0, (0, 1, 1), (2, 1, 1))]:
        ex = J(m, a, nvec, cut)
        nm = J_numeric(m, a, nvec, cut, dps=15)
        print('J', m, a, nvec, cut, ex, float(ex), float(nm), abs(float(ex) - float(nm)))
