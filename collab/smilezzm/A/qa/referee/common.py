"""Common helpers for the referee scripts (pure Python: mpmath).

Two-centre radial reduction in R^4 (|x| = 1):
  int d^4z g(|z|) f(|x-z|) = 2 pi  int int_{|r-s|<=1<=r+s} r s g(r) f(s) sqrt([(r+s)^2-1][1-(r-s)^2]) dr ds
(derived from d^4z = 4 pi r^3 sin^2(theta) dr dtheta and s^2 = 1 + r^2 - 2 r cos(theta)).
"""
import mpmath as mp


def Deps(r, eps):
    """D_eps(r) = (1 - exp(-r^2/4eps)) / (4 pi^2 r^2), with the r -> 0 limit."""
    r = mp.mpf(r)
    if r == 0:
        return 1 / (16 * mp.pi**2 * eps)
    u = r * r / (4 * eps)
    return -mp.expm1(-u) / (4 * mp.pi**2 * r * r)


def D0(r):
    return 1 / (4 * mp.pi**2 * mp.mpf(r)**2)


def Keps(r, eps):
    return mp.e**(-mp.mpf(r)**2 / (4 * eps)) / (4 * mp.pi * eps)**2


def heron(r, s):
    v = ((r + s)**2 - 1) * (1 - (r - s)**2)
    return mp.sqrt(v) if v > 0 else mp.mpf(0)


def kernel(r, s):
    return 2 * mp.pi * r * s * heron(r, s)


def small_breaks(eps, top):
    """breakpoints 0, sqrt(eps)*2^k ..., top for a variable with structure at scale sqrt(eps)."""
    se = mp.sqrt(eps)
    pts = [mp.mpf(0)]
    k = -4
    while se * 2**k < top:
        pts.append(se * mp.mpf(2)**k)
        k += 1
    pts.append(mp.mpf(top))
    return pts


def two_centre(g, f, eps, peak_g=True, peak_f=True, F=None):
    """int d^4z g(|z|) f(|x-z|), |x|=1, split into
       region A: s=|x-z| < 1/2 (then r in [1-s,1+s]),
       region B: s >= 1/2 (outer r in [0,inf), inner s in [max(1/2,|1-r|), 1+r]).
    Returns (value, error_estimate)."""
    half = mp.mpf(1) / 2
    errs = []
    if F is None:
        F = lambda r, s: g(r) * f(s)

    def innerA(s):
        v, e = mp.quad(lambda r: kernel(r, s) * F(r, s), [1 - s, 1, 1 + s], error=True)
        errs.append(e)
        return v

    ptsA = small_breaks(eps, half) if peak_f else [mp.mpf(0), half / 4, half]
    A = mp.mpf(0)
    for a, b in zip(ptsA[:-1], ptsA[1:]):
        A += mp.quad(innerA, [a, b])

    def innerB(r):
        lo = max(half, abs(1 - r))
        hi = 1 + r
        if hi <= lo:
            return mp.mpf(0)
        v, e = mp.quad(lambda s: kernel(r, s) * F(r, s), [lo, hi], error=True)
        errs.append(e)
        return v

    ptsB = small_breaks(eps, half) if peak_g else [mp.mpf(0), half / 4, half]
    ptsB += [mp.mpf(1), mp.mpf(3) / 2, 2, 3, 5, 10, 30, 100, mp.inf]
    B = mp.mpf(0)
    for a, b in zip(ptsB[:-1], ptsB[1:]):
        B += mp.quad(innerB, [a, b])
    return A + B, (max(errs) if errs else 0)


def fit_linear(xs_rows, ys):
    """least squares (exact solve if square) with mpmath matrices."""
    A = mp.matrix(xs_rows)
    y = mp.matrix(ys)
    if A.rows == A.cols:
        return mp.lu_solve(A, y)
    At = A.T
    return mp.lu_solve(At * A, At * y)


def two_centre_F(F, eps, peak_r=True, peak_s=True):
    """int d^4z F(|z|, |x-z|), |x| = 1 (same regions as two_centre, general integrand)."""
    return two_centre(None, None, eps, peak_r, peak_s, F=F)
