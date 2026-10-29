"""Exact singular coefficients of one-loop cycles with per-edge Schwinger cutoff t_e >= eps (twoletter_theory.tex,
Theorem thm:tl-simplex).  Pure sympy, exact rational / Gaussian-rational arithmetic.

A cycle with E edges, loop momentum l, edge e carrying e^{-t_e |l - q_e|^2} (t_e in [eps, infinity)), numerator N(l):
   w(eps) = int_{[eps,oo)^E} prod dt_e int d^4l/(2pi)^4 N(l) exp(-sum_e t_e |l - q_e|^2).
With tau = sum t_e, alpha_e = t_e/tau (simplex), c = sum alpha_e q_e, Phi = sum alpha_e |q_e|^2 - |c|^2 and
   (num)_j(alpha) = (Lap_l^j N)(c) / (j! 4^j)      (Lap = sum_mu d^2/dl_mu^2, Gaussian moments of Lemma c4-gauss)
the integrand is sum_j (num)_j / (16 pi^2 tau^{2+j}) e^{-tau Phi}, m_j = j + 2 - E, and (Lemma c4-local)
   16 pi^2 * [coefficient of eps^{-k}] = (1/k) sum_e int_{R_e} alpha_e^k F_k(alpha) d alpha,   k >= 1,
   16 pi^2 * [coefficient of log eps] = - int_simplex F_0(alpha) d alpha,
   F_k = sum_{j: m_j >= k} (num)_j (-Phi)^{m_j - k} / (m_j - k)!,
R_e = {alpha_e = min_f alpha_f} (the cone from the barycentre over the facet alpha_e = 0), d alpha = Lebesgue measure on
{alpha_1..alpha_{E-1} >= 0, sum <= 1}.  The total-time cutoff tau >= eps gives the same log coefficient and
(1/k) int_simplex F_k in place of the eps^{-k} coefficients (scheme dependence)."""
import itertools
from math import factorial
import sympy as sp

L = sp.symbols('l0:4')


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def lap(expr):
    return sum(sp.diff(expr, x, 2) for x in L)


def dirichlet_monomial(exps):
    """int over the standard (E-1)-simplex of prod alpha_e^{n_e}, E = len(exps)"""
    E = len(exps)
    num = 1
    for n in exps:
        num *= factorial(n)
    return sp.Rational(num, factorial(sum(exps) + E - 1))


def simplex_integral(poly_expr, alphas):
    """exact integral over the standard simplex of a polynomial in alphas (all E of them as independent symbols)"""
    P = sp.Poly(sp.expand(poly_expr), *alphas)
    tot = 0
    for mon, cf in P.terms():
        tot += cf * dirichlet_monomial(mon)
    return sp.expand(tot)


def region_integral(poly_expr, alphas, e):
    """integral over R_e = {alpha_e = min} of a polynomial in the alphas: cone from the barycentre over the facet
    alpha_e = 0; R_e has volume 1/E of the simplex, affine map from the standard simplex in mu"""
    E = len(alphas)
    mus = sp.symbols('mu0:%d' % E)
    # vertices: delta_f for f != e, and barycentre
    verts = [[1 if g == f else 0 for g in range(E)] for f in range(E) if f != e] + [[sp.Rational(1, E)] * E]
    sub = {alphas[g]: sum(mus[i] * verts[i][g] for i in range(E)) for g in range(E)}
    expr = sp.expand(sp.sympify(poly_expr).subs(sub, simultaneous=True))
    return sp.expand(sp.Rational(1, E) * simplex_integral(expr, mus))


def cycle_singular(N, qs, kmax=None):
    """N: sympy expression polynomial in L (and external symbols); qs: list of E 4-vectors (edge shifts).
    Returns dict {'log': b, k: a_k} with the common factor 1/(16 pi^2) removed."""
    E = len(qs)
    al = sp.symbols('a0:%d' % E)
    c = [sum(al[e] * qs[e][mu] for e in range(E)) for mu in range(4)]
    Phi = sp.expand(sum(al[e] * dot(qs[e], qs[e]) for e in range(E)) - dot(c, c))
    degl = sp.Poly(sp.expand(N), *L).total_degree() if sp.expand(N) != 0 else 0
    nums = []
    cur = sp.expand(N)
    for j in range(degl // 2 + 1):
        val = cur.subs({L[mu]: c[mu] for mu in range(4)}, simultaneous=True)
        nums.append(sp.expand(val / (factorial(j) * 4 ** j)))
        cur = sp.expand(lap(cur))
    res = {}
    # alpha_E constraint is imposed by the Dirichlet formula (alphas independent symbols, simplex sum = 1)
    F0 = 0
    for j, nj in enumerate(nums):
        m = j + 2 - E
        if m >= 0:
            F0 += nj * (-Phi) ** m / factorial(m)
    res['log'] = sp.expand(-simplex_integral(F0, al))
    mmax = max([j + 2 - E for j in range(len(nums))] + [0])
    if kmax is None:
        kmax = mmax
    for k in range(1, kmax + 1):
        Fk = 0
        for j, nj in enumerate(nums):
            m = j + 2 - E
            if m >= k:
                Fk += nj * (-Phi) ** (m - k) / factorial(m - k)
        if Fk == 0:
            res[k] = sp.Integer(0)
            continue
        tot = 0
        for e in range(E):
            tot += region_integral(sp.expand(al[e] ** k * Fk), al, e)
        res[k] = sp.expand(tot / k)
    return res


def cycle_singular_totaltime(N, qs, kmax=None):
    """same with the total-time cutoff tau >= eps (shape g = 1): eps^{-k} coefficients (1/k) int F_k"""
    E = len(qs)
    al = sp.symbols('a0:%d' % E)
    c = [sum(al[e] * qs[e][mu] for e in range(E)) for mu in range(4)]
    Phi = sp.expand(sum(al[e] * dot(qs[e], qs[e]) for e in range(E)) - dot(c, c))
    degl = sp.Poly(sp.expand(N), *L).total_degree() if sp.expand(N) != 0 else 0
    nums = []
    cur = sp.expand(N)
    for j in range(degl // 2 + 1):
        val = cur.subs({L[mu]: c[mu] for mu in range(4)}, simultaneous=True)
        nums.append(sp.expand(val / (factorial(j) * 4 ** j)))
        cur = sp.expand(lap(cur))
    res = {}
    F0 = sum(nj * (-Phi) ** (j + 2 - E) / factorial(j + 2 - E) for j, nj in enumerate(nums) if j + 2 - E >= 0)
    res['log'] = sp.expand(-simplex_integral(F0, al))
    mmax = max([j + 2 - E for j in range(len(nums))] + [0])
    if kmax is None:
        kmax = mmax
    for k in range(1, kmax + 1):
        Fk = sum(nj * (-Phi) ** (j + 2 - E - k) / factorial(j + 2 - E - k) for j, nj in enumerate(nums)
                 if j + 2 - E >= k)
        res[k] = sp.expand(simplex_integral(Fk, al) / k) if Fk != 0 else sp.Integer(0)
    return res


def special_edge_singular(N, qs, s, kmax=None):
    """cycle with edge s special (t_s in [0,eps], kernel -(K_eps - K_0)-type sign NOT included: the integrand is the
    same e^{-t_s|l-q_s|^2}) and the other edges t >= eps.  Pure negative powers (no log eps):
    16 pi^2 * [eps^{-k}] = (1/k) int_{R_s} [ (min_{f != s} alpha_f)^k - alpha_s^k ] F_k  (Lemma lem:tl-special).
    Implemented for E = 2 (min over the single other edge) and E = 3 (split R_s by the smaller of the two others)."""
    E = len(qs)
    al = sp.symbols('a0:%d' % E)
    c = [sum(al[e] * qs[e][mu] for e in range(E)) for mu in range(4)]
    Phi = sp.expand(sum(al[e] * dot(qs[e], qs[e]) for e in range(E)) - dot(c, c))
    degl = sp.Poly(sp.expand(N), *L).total_degree() if sp.expand(N) != 0 else 0
    nums = []
    cur = sp.expand(N)
    for j in range(degl // 2 + 1):
        val = cur.subs({L[mu]: c[mu] for mu in range(4)}, simultaneous=True)
        nums.append(sp.expand(val / (factorial(j) * 4 ** j)))
        cur = sp.expand(lap(cur))
    mmax = max([j + 2 - E for j in range(len(nums))] + [0])
    if kmax is None:
        kmax = mmax
    res = {'log': sp.Integer(0)}
    others = [f for f in range(E) if f != s]
    for k in range(1, kmax + 1):
        Fk = sum(nj * (-Phi) ** (j + 2 - E - k) / factorial(j + 2 - E - k) for j, nj in enumerate(nums)
                 if j + 2 - E >= k)
        if Fk == 0:
            res[k] = sp.Integer(0)
            continue
        if E == 2:
            f = others[0]
            tot = region_integral(sp.expand((al[f] ** k - al[s] ** k) * Fk), al, s)
        elif E == 3:
            # R_s split into alpha_s <= alpha_f <= alpha_g and alpha_s <= alpha_g <= alpha_f: simplices with
            # vertices (barycentre, midpoint of facet alpha_s=0 ... ) -- implement by the ordered-simplex vertices
            tot = 0
            for (f, g) in [(others[0], others[1]), (others[1], others[0])]:
                # region alpha_s <= alpha_f <= alpha_g : vertices delta_g, (delta_f+delta_g)/2, barycentre
                mus = sp.symbols('nu0:3')
                verts = []
                v1 = [0] * 3; v1[g] = 1
                v2 = [0] * 3; v2[f] = sp.Rational(1, 2); v2[g] = sp.Rational(1, 2)
                v3 = [sp.Rational(1, 3)] * 3
                verts = [v1, v2, v3]
                sub = {al[h]: sum(mus[i] * verts[i][h] for i in range(3)) for h in range(3)}
                integrand = sp.expand(((al[f] ** k) - al[s] ** k) * Fk)
                expr = sp.expand(integrand.subs(sub, simultaneous=True))
                tot += sp.Rational(1, 6) * simplex_integral(expr, mus)   # each ordered region: 1/6 of the simplex
        else:
            raise NotImplementedError
        res[k] = sp.expand(tot / k)
    return res


# ----------------------------------------------------------------------------------------------------------------------
# additions (resumed run): finite part of a special-edge cycle, Bessel generating functions, shape functions
# ----------------------------------------------------------------------------------------------------------------------
def _reduced(N, qs):
    E = len(qs)
    al = sp.symbols('a0:%d' % E)
    c = [sum(al[e] * qs[e][mu] for e in range(E)) for mu in range(4)]
    Phi = sp.expand(sum(al[e] * dot(qs[e], qs[e]) for e in range(E)) - dot(c, c))
    Ne = sp.expand(N)
    degl = sp.Poly(Ne, *L).total_degree() if Ne != 0 else 0
    nums = []
    cur = Ne
    for j in range(degl // 2 + 1):
        nums.append(sp.expand(cur.subs({L[mu]: c[mu] for mu in range(4)}, simultaneous=True)
                              / (factorial(j) * 4 ** j)))
        cur = sp.expand(lap(cur))
    return nums, Phi, al


def F_m(N, qs, m):
    """F_m(alpha) = sum_{j: m_j >= m} num_j (-Phi)^{m_j - m}/(m_j - m)!,  m_j = j + 2 - E"""
    nums, Phi, al = _reduced(N, qs)
    E = len(qs)
    out = 0
    for j, nj in enumerate(nums):
        mj = j + 2 - E
        if mj >= m:
            out += nj * (-Phi) ** (mj - m) / factorial(mj - m)
    return sp.expand(out), al


def special_edge_finite(N, qs, s):
    """eps^0 coefficient (times 16 pi^2) of a cycle whose edge s has Schwinger time in [0, eps] and the other edges in
    [eps, oo):  - int_{R_s} F_0(alpha) log(alpha_s / min_{f != s} alpha_f) d alpha   (Lemma lem:tl-special).
    The special numerator (the factor -|l - q_s|^2 of K_eps - K_0 = -int_0^eps H K_t dt) must be included in N.
    E = 2: exact (sympy); E = 3: mpmath, 30 digits."""
    E = len(qs)
    F0, al = F_m(N, qs, 0)
    others = [f for f in range(E) if f != s]
    if E == 2:
        f = others[0]
        x = sp.Symbol('x', positive=True)
        integrand = F0.subs({al[s]: x, al[f]: 1 - x}, simultaneous=True)
        val = -sp.integrate(sp.expand(integrand) * (sp.log(x) - sp.log(1 - x)), (x, 0, sp.Rational(1, 2)))
        return sp.nsimplify(sp.simplify(val))
    if E == 3:
        import mpmath as mp
        mp.mp.dps = 30
        fF = sp.lambdify(al, F0, 'mpmath')
        tot = mp.mpf(0)
        # R_s = {alpha_s <= alpha_f, alpha_s <= alpha_g}; parametrise alpha_s = x, alpha_f = y, alpha_g = 1 - x - y
        f, g = others

        def integrand(xx, yy):
            a = [None] * 3
            a[s] = xx; a[f] = yy; a[g] = 1 - xx - yy
            if a[g] < xx or yy < xx:
                return mp.mpf(0)
            return fF(*a) * mp.log(xx / min(a[f], a[g]))
        # region: 0 <= x <= 1/3, x <= y <= 1 - 2x
        tot = mp.quad(lambda xx: mp.quad(lambda yy: integrand(xx, yy), [xx, (1 - xx) / 2, 1 - 2 * xx]), [0, mp.mpf(1) / 3])
        return -tot
    raise NotImplementedError


def bessel_series(q, z, order):
    """(|sigma|/(2 sqrt Phi))^q I_q(|sigma| sqrt Phi) = sum_i (|sigma|^2/4)^{i+q} Phi^i / (i!(i+q)!), returned as a
    polynomial in (s2 = |sigma|^2, Phi) truncated at total |sigma|-degree <= order;  z = (s2, Phi)"""
    s2, Ph = z
    out = 0
    i = 0
    while 2 * (i + q) <= order:
        out += (s2 / 4) ** (i + q) * Ph ** i / (factorial(i) * factorial(i + q))
        i += 1
    return sp.expand(out)


def genfun_F(m, qs, sigma, order):
    """Generating function of F_m for the numerator e^{i(sigma,l)} (truncated at |sigma|-degree <= order):
    F_m^gen = (-1)^{m+E} e^{i(sigma,c)} (|sigma|/(2 sqrt Phi))^{m+E-2} I_{m+E-2}(|sigma| sqrt Phi)  (Prop prop:tl-genfun).
    sigma: 4-vector of sympy numbers times a formal scale lam (included by the caller)."""
    E = len(qs)
    al = sp.symbols('a0:%d' % E)
    c = [sum(al[e] * qs[e][mu] for e in range(E)) for mu in range(4)]
    Phi = sp.expand(sum(al[e] * dot(qs[e], qs[e]) for e in range(E)) - dot(c, c))
    q = m + E - 2
    if q < 0:
        return sp.Integer(0), al
    lam = sp.Symbol('lam')
    s2 = dot(sigma, sigma)
    sc = dot(sigma, c)
    expo = sum((sp.I * sc) ** n / factorial(n) for n in range(order + 1))
    bes = bessel_series(q, (s2, Phi), order)
    tot = sp.expand((-1) ** q * expo * bes)
    # truncate in lam
    tot = sp.expand(sum(tot.coeff(lam, n) * lam ** n for n in range(order + 1)))
    return tot, al


def bubble_log_closed(a, b, j):
    """coefficient in the closed form of the scalar bubble log kernel (Prop prop:tl-bubble):
    (a+j)!(b+j)!/(a! b! (j!)^2 (a+b+2j+1)!)"""
    return sp.Rational(factorial(a + j) * factorial(b + j), factorial(a) * factorial(b) * factorial(j) ** 2
                       * factorial(a + b + 2 * j + 1))
