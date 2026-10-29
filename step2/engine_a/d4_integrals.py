"""Step 2, Part D4, engine A: exact eps-expansion of one-loop cycles with all edges P(eps,infty).

For n edges with momenta k_i = k + q_i (q_1 = 0), numerator N(k) (polynomial, exact):
  16 pi^2 int d^4k/(2pi)^4 N(k) prod_i int_eps^infty dt_i e^{-t_i k_i^2}
     = sum_j (4^j j!)^{-1} int_{simplex} du F_j(u) I_{n-3-j}(eps/m(u), q(u)),
  F_j(u) = (Delta^j N)(-K(u)),  K(u) = sum u_i q_i,  q(u) = sum u_i q_i^2 - K(u)^2,  m(u) = min_i u_i,
  I_a(x, q) = int_x^infty T^a e^{-q T} dT      (t_i = T u_i, prod dt_i = T^{n-1} dT du).
Small-x expansion (x = eps/m):
  a >= 0 :  I_a = Gamma(a+1) q^{-a-1} + O(x^{a+1})                      -> finite part only
  a = -n :  I_a = q^{n-1} (-1)^{n-1}/(n-1)! (H_{n-1} - gamma - log q - log x)
                  + sum_{k=0}^{n-2} (-1)^k q^k x^{k-n+1} / (k! (n-1-k)) + O(x log x)
The O(x log x) remainders integrate to O(eps log eps) (the only non-uniformity, m(u) -> 0, is logarithmic).
Two-point cycles: q(u) = V(1-V)|p|^2 with V = sum of u_i over the edges carrying k + (external momentum).
Exact integrals: Dirichlet aggregation in V (Beta functions and digamma differences, rational), and the
piecewise-polynomial integrals of u^beta m(u)^j, u^beta log m(u) over the simplex (n <= 3: sympy, exact; values in
Q + Q log 2 (n = 2), Q + Q log 3 (n = 3)).  Output: dict over the basis
   ('eps^-2'), ('eps^-1'), ('log eps'), ('1'), ('gamma'), ('log 2'), ('log 3'), ('log p^2')   with Q2 coefficients."""
import math, itertools, functools
from fractions import Fraction
import sympy as sp
from r2_core import Q2
from s2_graphs import Poly

BASIS = ['eps^-2', 'eps^-1', 'log eps', '1', 'gamma', 'log 2', 'log 3', 'log p^2']


def H(n):
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


def addto(acc, key, val):
    if isinstance(val, Q2):
        if val.iszero():
            return
    acc[key] = acc[key] + val if key in acc else val


def u_forms(qs, nu):
    """-K(u) as 4 Polys in nu variables (u_2..u_n) and q(u) Poly; qs = offsets (lists/arrays of 4 ints/Fractions)"""
    forms = []
    for mu in range(4):
        forms.append(Poly(nu, {tuple(1 if a == i else 0 for a in range(nu)): Q2(-Fraction(qs[i + 1][mu]))
                               for i in range(nu) if qs[i + 1][mu] != 0}))
    Kf = [f.scale(Q2(-1)) for f in forms]
    q = Poly(nu)
    for i in range(nu):
        qq = sum(Fraction(x) * Fraction(x) for x in qs[i + 1])
        if qq:
            q = q + Poly.var(nu, i, Q2(qq))
    for mu in range(4):
        q = q - Kf[mu] * Kf[mu]
    return forms, Kf, q


def terms(Npoly, qs, degree, dN=None, dqs=None):
    """list of (P(u), b, weight): integrand P(u) I_b(eps/m, q_0(u)); if dN/dqs given, the t-derivative at t=0 of the
    family N_t, q_i(t) = q_i + t dq_i  (dN = d/dt N_t at t=0)"""
    n = len(qs); nu = n - 1
    forms, Kf, q = u_forms(qs, nu)
    lap = [Npoly]
    for j in range(1, degree // 2 + 2):
        lap.append(lap[-1].lap())
    out = []
    if dN is None:
        for j, Lj in enumerate(lap):
            if not Lj.t:
                continue
            Fj = Lj.subs_linear(forms, nu) if nu > 0 else Poly(0, {(): Lj.t.get((0, 0, 0, 0), Q2(0))})
            out.append((Fj, n - 3 - j, Fraction(1, 4 ** j * math.factorial(j))))
        return out
    # derivative terms
    dK = []
    for mu in range(4):
        dK.append(Poly(nu, {tuple(1 if a == i else 0 for a in range(nu)): Q2(Fraction(dqs[i + 1][mu]))
                            for i in range(nu) if dqs[i + 1][mu] != 0}))
    dq = Poly(nu)
    for i in range(nu):
        c = sum(2 * Fraction(qs[i + 1][mu]) * Fraction(dqs[i + 1][mu]) for mu in range(4))
        if c:
            dq = dq + Poly.var(nu, i, Q2(c))
    for mu in range(4):
        dq = dq - (Kf[mu] * dK[mu]).scale(Q2(2))
    dlap = [dN]
    for j in range(1, degree // 2 + 2):
        dlap.append(dlap[-1].lap())
    for j in range(len(lap)):
        w = Fraction(1, 4 ** j * math.factorial(j))
        A = dlap[j].subs_linear(forms, nu) if dlap[j].t else Poly(nu)
        for mu in range(4):
            g = lap[j].deriv(mu)
            if g.t and dK[mu].t:
                A = A - g.subs_linear(forms, nu) * dK[mu]
        if A.t:
            out.append((A, n - 3 - j, w))
        if lap[j].t and dq.t:
            B = lap[j].subs_linear(forms, nu) * dq
            out.append((B.scale(Q2(-1)), n - 3 - j + 1, w))
    return out


# ------------------------------------------------------------------ exact simplex integrals
def digamma_diff(a, b):
    """psi(a) - psi(b) for positive integers a, b"""
    return H(a - 1) - H(b - 1)


def beta_int(A, B):
    """int_0^1 V^A (1-V)^B dV, A, B integers > -1"""
    if A < 0 or B < 0:
        raise ValueError('divergent Beta integral %d %d' % (A, B))
    return Fraction(math.factorial(A) * math.factorial(B), math.factorial(A + B + 1))


def homogenize(P, n):
    """P: Poly in u_2..u_n (u_1 = 1 - sum) -> dict {exponent tuple over u_1..u_n: Q2}, homogeneous of degree max deg P"""
    if n == 1:
        return {(0,): P.t.get((), Q2(0))} if P.t else {}
    D = P.degree()
    out = {}
    for beta, c in P.t.items():
        d = sum(beta)
        # multiply by (u_1 + ... + u_n)^(D - d): multinomial expansion
        for comp in _compositions(D - d, n):
            mc = math.factorial(D - d)
            for x in comp:
                mc //= math.factorial(x)
            e = (comp[0],) + tuple(beta[i] + comp[i + 1] for i in range(n - 1))
            out[e] = out[e] + c * mc if e in out else c * mc
    return {k: v for k, v in out.items() if not v.iszero()}


@functools.lru_cache(maxsize=None)
def _compositions(N, k):
    if k == 1:
        return [(N,)]
    out = []
    for i in range(N + 1):
        for rest in _compositions(N - i, k - 1):
            out.append((i,) + rest)
    return out


def mul_group_powers(mono, c, n, G2full, a, b):
    """multiply the homogeneous monomial by V^a (1-V)^b with V = sum_{G2} u_i, 1-V = sum_{G1} u_i (a, b >= 0)"""
    G1full = [i for i in range(n) if i not in G2full]
    res = {mono: c}
    for grp, pw in ((G2full, a), (G1full, b)):
        for _ in range(pw):
            new = {}
            for e, v in res.items():
                for i in grp:
                    e2 = list(e); e2[i] += 1; e2 = tuple(e2)
                    new[e2] = new[e2] + v if e2 in new else v
            res = new
    return res


def dirichlet_h(beta, G2full, n, a, b, log):
    """int over the simplex (u_1..u_n, sum = 1) of prod u_i^beta_i V^a (1-V)^b {1, log V, log(1-V)}"""
    G1full = [i for i in range(n) if i not in G2full]
    n2 = len(G2full); n1 = len(G1full)
    s2 = sum(beta[i] for i in G2full); s1 = sum(beta[i] for i in G1full)
    c = Fraction(1)
    for i in range(n):
        c *= math.factorial(beta[i])
    c /= math.factorial(n2 - 1 + s2)
    c /= math.factorial(n1 - 1 + s1)
    A = n2 - 1 + s2 + a; B = n1 - 1 + s1 + b
    base = beta_int(A, B)
    if log is None:
        return c * base
    if log == 'V':
        return c * base * digamma_diff(A + 1, A + B + 2)
    if log == '1-V':
        return c * base * digamma_diff(B + 1, A + B + 2)
    raise ValueError(log)


_s, _v, _w = sp.symbols('s v w', nonnegative=True)


@functools.lru_cache(maxsize=None)
def region_int(n, expo, what):
    """exact int over the (n-1)-simplex (u_1..u_n) of prod u_i^expo_i * g(m(u)), g = m^j (what = ('m', j)) or log m
    (what = ('log',)); returns (rational, coefficient of log 2, coefficient of log 3)"""
    g = (lambda x: x ** what[1]) if what[0] == 'm' else (lambda x: sp.log(x))
    if n == 2:
        e1, e2 = expo
        tot = 0
        # min = u2 = s in [0,1/2], u1 = 1 - s ;  min = u1 = s, u2 = 1 - s
        for (U1, U2) in (((1 - _s), _s), (_s, (1 - _s))):
            tot += sp.integrate(sp.expand(U1 ** e1 * U2 ** e2) * g(_s), (_s, 0, sp.Rational(1, 2)))
    elif n == 3:
        tot = 0
        # min = u_i = s; the other two are s + v and 1 - 2s - v, v in [0, 1-3s], s in [0,1/3]
        for imin in range(3):
            others = [j for j in range(3) if j != imin]
            U = [None] * 3
            U[imin] = _s; U[others[0]] = _s + _v; U[others[1]] = 1 - 2 * _s - _v
            mono = 1
            for j in range(3):
                mono *= U[j] ** expo[j]
            inner = sp.integrate(sp.expand(mono), (_v, 0, 1 - 3 * _s))
            tot += sp.integrate(sp.expand(inner) * g(_s), (_s, 0, sp.Rational(1, 3)))
    else:
        raise NotImplementedError(n)
    tot = sp.expand(sp.expand_log(sp.nsimplify(tot), force=True))
    c2 = tot.coeff(sp.log(2)); c3 = tot.coeff(sp.log(3))
    c0 = sp.expand(tot - c2 * sp.log(2) - c3 * sp.log(3))
    assert c0.is_rational and c2.is_rational and c3.is_rational, tot
    return (Fraction(int(c0.p), int(c0.q)), Fraction(int(c2.p), int(c2.q)), Fraction(int(c3.p), int(c3.q)))


def expand_two_point(term_list, n, G2, p2):
    """sum over terms (P, b, w) of w int P I_b(eps/m, V(1-V) p2) -> dict over BASIS (Q2 values); p2 = |p|^2 (int).
    G2: indices (0-based among u_2..u_n) of the edges carrying k + external momentum."""
    acc = {}
    G2full = [i + 1 for i in G2]
    for (P, b, w) in term_list:
        wq = Q2(w)
        if n == 1:
            for beta, c in P.t.items():
                assert b <= -2, ('IR/UV divergent single-edge term', b)
                key = {-2: 'eps^-1', -3: 'eps^-2'}[b]
                addto(acc, key, c * wq * Fraction(1, -b - 1))
            continue
        Ph = homogenize(P, n)
        for beta, c in Ph.items():
            cw = c * wq
            if b >= 0:
                val = dirichlet_h(beta, G2full, n, -b - 1, -b - 1, None) * math.factorial(b) * Fraction(1, p2 ** (b + 1))
                addto(acc, '1', cw * val)
                continue
            N_ = -b
            pref = Fraction((-1) ** (N_ - 1), math.factorial(N_ - 1)) * Fraction(p2) ** (N_ - 1)
            base = dirichlet_h(beta, G2full, n, N_ - 1, N_ - 1, None)
            addto(acc, 'log eps', cw * (-pref * base))
            addto(acc, '1', cw * (pref * base * H(N_ - 1)))
            addto(acc, 'gamma', cw * (-pref * base))
            addto(acc, 'log p^2', cw * (-pref * base))
            lv = dirichlet_h(beta, G2full, n, N_ - 1, N_ - 1, 'V') + dirichlet_h(beta, G2full, n, N_ - 1, N_ - 1, '1-V')
            addto(acc, '1', cw * (-pref * lv))
            for be2, c2 in mul_group_powers(beta, Q2(1), n, G2full, N_ - 1, N_ - 1).items():
                r0, r2, r3 = region_int(n, tuple(be2), ('log',))
                addto(acc, '1', cw * c2 * (pref * r0))
                addto(acc, 'log 2', cw * c2 * (pref * r2))
                addto(acc, 'log 3', cw * c2 * (pref * r3))
            for k in range(0, N_ - 1):
                jpow = N_ - 1 - k
                key = {1: 'eps^-1', 2: 'eps^-2'}[jpow]
                ck = Fraction((-1) ** k, math.factorial(k) * jpow) * Fraction(p2) ** k
                for be2, c2 in mul_group_powers(beta, Q2(1), n, G2full, k, k).items():
                    r0, r2, r3 = region_int(n, tuple(be2), ('m', jpow))
                    assert r2 == 0 and r3 == 0
                    addto(acc, key, cw * c2 * ck * r0)
    return {k: v for k, v in acc.items() if not (isinstance(v, Q2) and v.iszero())}
