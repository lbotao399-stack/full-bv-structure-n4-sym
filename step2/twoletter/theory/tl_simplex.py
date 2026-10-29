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
