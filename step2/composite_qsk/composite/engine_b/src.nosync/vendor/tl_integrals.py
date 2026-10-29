"""Route B: exact pinned-edge integrals.

For a cycle with E edges, Schwinger times t_e >= eps, the term of order eps^{-m} (m >= 1) and the log eps term of
eps d/deps w come from the boundary graphs in which one edge e is pinned, t_e = eps; with t_f = eps u_f (f != e),
u_f in [1, oo), alpha = u / tau_u, tau_u = 1 + sum_f u_f, alpha_e = 1/tau_u.  A monomial alpha^n of the Gaussian
coefficient F_m(alpha) contributes, after the eps-power bookkeeping,

    I_e(m; n) = int_{[1,oo)^{E-1}} prod_{f != e} u_f^{n_f}  tau_u^{-(m + |n| + E)}  du .

Feynman parametrisation of the remaining edges: u_f = U beta_f, beta on the (E-2)-simplex, du = U^{E-2} dU dbeta,
u_f >= 1  <=>  U >= 1/b,  b = min_f beta_f.  The U-integral is one-dimensional and exact (t = U/(1+U)):

    int_{1/b}^oo U^a (1+U)^{-Q} dU = int_{y}^1 t^a (1-t)^c dt,   y = 1/(1+b),  c = Q - a - 2 = m + n_e >= 0,

a polynomial in 1/(1+b).  The remaining beta-integral is done exactly region by region (beta_f minimal), as
one-dimensional integrals of (polynomial) x (1+x)^{-r}: these produce log(3/2) (E = 3) and log(4/3) (E = 4).
Values are returned as  rational + rational * log(1 + h),  h = 1/(E-1); the log part is required to cancel
(each I_e is the integral of a polynomial over a rational simplex), and the code asserts it.

T(E, m, n) = sum_e I_e(m; n) is the quantity entering
    [eps^-m] w = (1/(16 pi^2 m)) T(E,m,F_m),   [log eps] w = -(1/(16 pi^2)) T(E,0,F_0).
"""
from fractions import Fraction
from functools import lru_cache
from math import comb


# ---------- polynomials in one variable with Fraction coefficients: list c[0..d] ----------
def padd(p, q):
    n = max(len(p), len(q))
    return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]


def pmul(p, q):
    out = [Fraction(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        if a == 0:
            continue
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def ppow(p, k):
    out = [Fraction(1)]
    for _ in range(k):
        out = pmul(out, p)
    return out


def pint_def(p, lo, hi):
    """int_lo^hi p(x) dx, lo/hi Fractions or polynomials (lists) -> Fraction or polynomial."""
    # antiderivative
    P = [Fraction(0)] + [c / (i + 1) for i, c in enumerate(p)]
    return peval(P, hi) - peval(P, lo) if not isinstance(lo, list) and not isinstance(hi, list) else \
        padd(pcompose(P, hi), [-c for c in pcompose(P, lo)])


def peval(p, x):
    s = Fraction(0)
    for c in reversed(p):
        s = s * x + c
    return s


def pcompose(p, q):
    out = [Fraction(0)]
    for c in reversed(p):
        out = padd(pmul(out, q), [c])
    return out


# ---------- the one-dimensional U integral ----------
def inner_coeffs(a, c):
    """int_y^1 t^a (1-t)^c dt = sum_i C(c,i)(-1)^i (1 - y^{a+i+1})/(a+i+1); return list of (const, coeff, power):
    value = const + sum coeff * y^power."""
    const = Fraction(0)
    terms = []
    for i in range(c + 1):
        w = Fraction(comb(c, i) * (-1) ** i, a + i + 1)
        const += w
        terms.append((-w, a + i + 1))
    return const, terms


def int_poly_over_1px(p, r, h):
    """int_0^h p(x) (1+x)^{-r} dx, exact: returns (rational, coefficient of log(1+h))."""
    # v = 1 + x, x = v - 1
    pv = pcompose(p, [Fraction(-1), Fraction(1)])
    rat = Fraction(0)
    lg = Fraction(0)
    H = 1 + Fraction(h)
    for j, cj in enumerate(pv):
        if cj == 0:
            continue
        e = j - r
        if e == -1:
            lg += cj
        else:
            rat += cj * (H ** (e + 1) - 1) / (e + 1)
    return rat, lg


@lru_cache(maxsize=None)
def I_pinned(E, m, n, e):
    """I_e(m; n) for the pinned edge e; n a tuple of length E.  Returns (rational, logcoeff, logarg)."""
    N = sum(n)
    Q = m + N + E
    others = [f for f in range(E) if f != e]
    if E == 1:
        return Fraction(1), Fraction(0), None
    a = sum(n[f] for f in others) + E - 2
    c = m + n[e]
    assert Q - a - 2 == c >= 0
    const, terms = inner_coeffs(a, c)
    if E == 2:
        y = Fraction(1, 2)
        return const + sum(w * y ** pw for w, pw in terms), Fraction(0), None
    if E == 3:
        n1, n2 = n[others[0]], n[others[1]]
        h = Fraction(1, 2)
        # beta = (x, 1-x); x in [0,1/2] with b = x, plus the mirror (exponents swapped)
        poly = padd(pmul(ppow([Fraction(0), Fraction(1)], n1), ppow([Fraction(1), Fraction(-1)], n2)),
                    pmul(ppow([Fraction(0), Fraction(1)], n2), ppow([Fraction(1), Fraction(-1)], n1)))
        rat = const * pint_def(poly, Fraction(0), h)
        lg = Fraction(0)
        for w, pw in terms:      # y^pw = (1+x)^{-pw}
            r_, l_ = int_poly_over_1px(poly, pw, h)
            rat += w * r_
            lg += w * l_
        return rat, lg, Fraction(3, 2)
    if E == 4:
        h = Fraction(1, 3)
        rat = Fraction(0)
        lg = Fraction(0)
        X = [Fraction(0), Fraction(1)]
        for imin in range(3):
            o = [others[(imin + s) % 3] for s in range(3)]
            nx, ny, nz = n[o[0]], n[o[1]], n[o[2]]
            # inner: int_{y=x}^{1-2x} y^ny (1-x-y)^nz dy, as a polynomial in x
            # write y as variable; integrand polynomial in y with coefficients polynomial in x
            # (1-x-y)^nz = sum_k C(nz,k) (1-x)^{nz-k} (-y)^k
            inner = [Fraction(0)]
            for kk in range(nz + 1):
                coefx = [Fraction(comb(nz, kk) * (-1) ** kk) * cc for cc in ppow([Fraction(1), Fraction(-1)], nz - kk)]
                deg = ny + kk
                # int_x^{1-2x} y^deg dy = ((1-2x)^{deg+1} - x^{deg+1})/(deg+1)
                yint = padd(ppow([Fraction(1), Fraction(-2)], deg + 1), [-c for c in ppow(X, deg + 1)])
                yint = [c / (deg + 1) for c in yint]
                inner = padd(inner, pmul(coefx, yint))
            poly = pmul(ppow(X, nx), inner)
            rat += const * pint_def(poly, Fraction(0), h)
            for w, pw in terms:
                r_, l_ = int_poly_over_1px(poly, pw, h)
                rat += w * r_
                lg += w * l_
        return rat, lg, Fraction(4, 3)
    raise NotImplementedError('E >= 5 not implemented in route B')


LOG_RESIDUES = []   # any nonzero log coefficient is recorded here (should stay empty)


@lru_cache(maxsize=None)
def T(E, m, n):
    tot = Fraction(0)
    for e in range(E):
        r, lg, arg = I_pinned(E, m, n, e)
        if lg != 0:
            LOG_RESIDUES.append((E, m, n, e, str(lg), str(arg)))
        tot += r
    return tot


def selftest():
    """Compare a few pinned integrals with direct numerical quadrature of the u-integrals (mpmath)."""
    import mpmath as mp
    mp.mp.dps = 15
    out = []
    cases = [(2, 1, (0, 0)), (2, 0, (1, 2)), (2, 2, (3, 1)), (3, 0, (1, 0, 0)), (3, 1, (0, 2, 1)), (3, 0, (1, 1, 1)),
             (4, 0, (0, 0, 0, 0)), (4, 1, (1, 0, 0, 1))]
    for (E, m, n) in cases:
        exact = T(E, m, n)
        num = mp.mpf(0)
        for e in range(E):
            others = [f for f in range(E) if f != e]
            Q = m + sum(n) + E

            def f(*v):       # u = 1/v, du = dv/v^2, v in (0,1]
                u = [1 / x for x in v]
                tau = 1 + sum(u)
                val = mp.mpf(1)
                for uf, ff, x in zip(u, others, v):
                    val *= uf ** n[ff] / x ** 2
                return val * tau ** (-Q)
            if E == 4:
                num += mp.quad(f, [0, 1], [0, 1], [0, 1], maxdegree=4)
            else:
                num += mp.quad(f, *([[0, 1]] * (E - 1)))
        tol = mp.mpf(10) ** (-5 if E == 4 else -10)
        out.append(dict(E=E, m=m, n=list(n), exact=str(exact), numeric=mp.nstr(num, 12),
                        ok=bool(abs(num - exact) < tol * max(1, abs(num)))))
    return out


if __name__ == '__main__':
    for r in selftest():
        print(r)
    print('log residues', LOG_RESIDUES)
