"""Exact linear algebra over Q(i) and the plane-wave values of the two-leg basis functionals (translation-invariant
canonical densities, second leg momentum q = -p).  Multilinear convention: coefficient of the odd leg parameters in leg
order, field slots as printed (r2_core / mode-model convention); colour stripped ((Y0,Y1) -> 1).
Clifford data: rho(e_M)|_Sigma = sb[M], pairing (u,w) = u^T w, Gamma(xi,lam)^M = xi^T sb[M] lam (cliff.py)."""
from fractions import Fraction
import numpy as np
import ea_core as E
from r2_core import Q2

SB = [np.asarray(m) for m in E.SGB]       # rho(e_M) on Sigma
SS = [np.asarray(m) for m in E.SG]        # rho(e_M) on Sigma^*
BASIS_NAMES = ['1', 'log2', 'log3', 'log5', 'log7']


def q2(z):
    if isinstance(z, Q2):
        return z
    if isinstance(z, tuple):
        return Q2(z[0], z[1])
    if isinstance(z, (int, Fraction)):
        return Q2(z)
    z = complex(z)
    re, im = Fraction(z.real).limit_denominator(1 << 30), Fraction(z.imag).limit_denominator(1 << 30)
    assert abs(float(re) - z.real) < 1e-9 and abs(float(im) - z.imag) < 1e-9
    return Q2(re, im)


def solve(rows, rhs):
    """exact least-squares-free solve of rows x = rhs (rows: list of lists of Q2-able, rhs: list of Q2-able).
    Returns (x list of Q2 or None if inconsistent, rank, free columns)"""
    m = len(rows); n = len(rows[0])
    A = [[q2(v) for v in r] + [q2(b)] for r, b in zip(rows, rhs)]
    piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, m) if not A[i][c].iszero()), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        inv = Q2(1) / A[r][c]
        A[r] = [v * inv for v in A[r]]
        for i in range(m):
            if i != r and not A[i][c].iszero():
                f = A[i][c]
                A[i] = [vi - f * vr for vi, vr in zip(A[i], A[r])]
        piv.append(c); r += 1
        if r == m:
            break
    for i in range(r, m):
        if not A[i][n].iszero():
            return None, r, None
    x = [Q2(0)] * n
    for i, c in enumerate(piv):
        x[c] = A[i][n]
    free = [c for c in range(n) if c not in piv]
    return x, r, free


def fit(samples, nb):
    """samples: list of (basis values list (length nb), 5-vector of (re,im)); returns ({k: 5-vector of Q2}, rank) or
    raises if inconsistent in any log component"""
    rows = [s[0] for s in samples]
    out = [[Q2(0)] * 5 for _ in range(nb)]
    rank = None
    for comp in range(5):
        rhs = [q2(s[1][comp]) for s in samples]
        if all(v.iszero() for v in rhs):
            _, rank, _ = solve(rows, [0] * len(rows))
            continue
        x, rank, free = solve(rows, rhs)
        if x is None:
            raise ValueError('inconsistent fit (component %s)' % BASIS_NAMES[comp])
        for k in range(nb):
            out[k][comp] = x[k]
    return out, rank


def vstr(v):
    """5-vector of Q2 (or (re,im)) -> string like '-11/3 + 2*log2'"""
    parts = []
    for c, nm in zip(v, BASIS_NAMES):
        c = q2(c)
        if c.iszero():
            continue
        s = str(c)
        parts.append(s if nm == '1' else '(%s)*%s' % (s, nm))
    return ' + '.join(parts) if parts else '0'


def vfrac(v):
    """pure rational real 5-vector -> Fraction (asserts)"""
    c = [q2(x) for x in v]
    assert all(x.im == 0 for x in c) and all(x.iszero() for x in c[1:]), vstr(v)
    return c[0].re


# ------------------------------------------------------------------------------------------------ two-leg bases
def _V(x):
    return np.asarray(x)[:4]


def slash(p):
    return sum(p[m] * SB[m] for m in range(4))


def b_dAVdAV(L):
    (a, p), (b, q) = (L[0]['pol'], L[0]['p0']), (L[1]['pol'], L[1]['p0'])
    a, b = _V(a), _V(b)
    return -2 * ((p @ q) * (a @ b) - (p @ b) * (q @ a))


def b_dsAdsA(L):
    (a, p), (b, q) = (L[0]['pol'], L[0]['p0']), (L[1]['pol'], L[1]['p0'])
    return -2 * (p @ _V(a)) * (q @ _V(b))


def b_dAWdAW(L):
    (a, p), (b, q) = (L[0]['pol'], L[0]['p0']), (L[1]['pol'], L[1]['p0'])
    return -2 * (p @ q) * (np.asarray(a)[4:] @ np.asarray(b)[4:])


def b_ldl(L):
    """int (lam, dslash lam): legs (lam u, p), (lam w, q): i u^T (qslash - pslash) w"""
    (u, p), (w, q) = (L[0]['pol'], L[0]['p0']), (L[1]['pol'], L[1]['p0'])
    return 1j * (u @ ((slash(q) - slash(p)) @ w))


def b_dcAs(L):
    """int (dc, A^*): legs (c, p), (A^*, s, q): i p.s_V"""
    p = L[0]['p0']; s = L[1]['pol']
    return 1j * L[0]['pol'][0] * (p @ _V(s))


def b_AA(L):
    """int (A, A) all ten components: 2 a.b"""
    return 2 * (np.asarray(L[0]['pol']) @ np.asarray(L[1]['pol']))


def b_ddcdA(L):
    """int (d^*dc, d^*A_V): legs (c,p), (A_V, a, q): |p|^2 (-i q.a)"""
    p = L[0]['p0']; a, q = L[1]['pol'], L[1]['p0']
    return -1j * L[0]['pol'][0] * (p @ p) * (q @ _V(a))


# along xi (degree 0, weight 9/2): xi in Sigma (16-vector)
def Gam(xi, u, Ms):
    return np.array([xi @ (SB[M] @ u) for M in Ms])


def b_g1(L, xi):
    """int (Gamma_V, d^*d A_V), legs (A_V a, q), (lam u, p)"""
    a, q = _V(L[0]['pol']), L[0]['p0']; u = L[1]['pol']
    G = Gam(xi, u, range(4))
    return G @ ((q @ q) * a - q * (q @ a))


def b_g2(L, xi):
    a, q = _V(L[0]['pol']), L[0]['p0']; u = L[1]['pol']
    G = Gam(xi, u, range(4))
    return (G @ q) * (q @ a)


def b_g3(L, xi):
    """int (Gamma_perp, d^*d A_perp), legs (A_perp a, q), (lam u, p)"""
    a, q = np.asarray(L[0]['pol'])[4:], L[0]['p0']; u = L[1]['pol']
    G = Gam(xi, u, range(4, 10))
    return (q @ q) * (G @ a)


def b_g4(L, xi):
    """int (d^*dc, (xi, lam^*)), legs (c, p), (lam^* s)"""
    p = L[0]['p0']; s = L[1]['pol']
    return L[0]['pol'][0] * (p @ p) * (xi @ s)


def b_GVA(L, xi):
    """int (Gamma_V, A_V) (pole structure), legs (A_V a), (lam u)"""
    return Gam(xi, L[1]['pol'], range(4)) @ _V(L[0]['pol'])


def b_GPA(L, xi):
    return Gam(xi, L[1]['pol'], range(4, 10)) @ np.asarray(L[0]['pol'])[4:]
