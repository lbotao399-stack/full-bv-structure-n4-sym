"""Numerical checks of the simplex formulas of twoletter_theory.tex (Theorem thm:tl-simplex, Proposition
prop:tl-shape, Lemma lem:tl-special) against direct Schwinger integrals at small eps.

For each test cycle (E edges, shifts q_e, numerator N(l)) the Gaussian loop integral is done with the moments of
Lemma c4-gauss (exp(Lap/(4 tau))), and the remaining Schwinger integral is evaluated numerically (mpmath) through the
boundary form of Proposition c4-boundary:
    eps d/deps w = - sum_e eps * w(t_e = eps)   =  - sum_k k a_k eps^{-k} + b + o(1),
i.e. an (E-1)-dimensional integral per face.  A least-squares fit in {eps^{-k}, 1, eps, eps log eps} at several eps
gives (k a_k, b), compared with the exact simplex formulas (tl_simplex.py).  Also:
  * the total-time cutoff tau >= eps: same b, different a_k (scheme dependence),
  * one special edge t_s in [0,eps]: the fit has no log eps, and the eps^{-k} coefficients match lemma lem:tl-special,
  * the Gaussian moments against a direct Monte-Carlo-free radial check.
Writes checks_simplex_results.json."""
import json, os, math
import mpmath as mp
import sympy as sp
import tl_simplex as T

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 20
L = T.L


def reduced_integrand(N, qs):
    """returns (nums_j as sympy in alphas, Phi, alphas) for the Gaussian-integrated integrand
    sum_j nums_j(alpha) tau^{-2-j} exp(-tau Phi(alpha))  (factor 1/16pi^2 dropped)"""
    E = len(qs)
    al = sp.symbols('a0:%d' % E)
    c = [sum(al[e] * qs[e][mu] for e in range(E)) for mu in range(4)]
    Phi = sp.expand(sum(al[e] * T.dot(qs[e], qs[e]) for e in range(E)) - T.dot(c, c))
    deg = sp.Poly(sp.expand(N), *L).total_degree()
    nums = []
    cur = sp.expand(N)
    for j in range(deg // 2 + 1):
        nums.append(sp.expand(cur.subs({L[m]: c[m] for m in range(4)}, simultaneous=True) / (math.factorial(j) * 4 ** j)))
        cur = sp.expand(T.lap(cur))
    return nums, Phi, al


def make_f(N, qs):
    nums, Phi, al = reduced_integrand(N, qs)
    fn = [sp.lambdify(al, nj, 'mpmath') for nj in nums]
    fP = sp.lambdify(al, Phi, 'mpmath')

    def f(ts):
        tau = sum(ts)
        a = [t / tau for t in ts]
        ph = fP(*a)
        return sum(fn[j](*a) * tau ** (-2 - j) for j in range(len(fn))) * mp.e ** (-tau * ph)
    return f


def face_value(f, E, e, eps, special=None, total_time=False):
    """eps * (integral with t_e = eps fixed, other t's in [eps, oo)); with special=s: t_s in [0,eps]"""
    others = [g for g in range(E) if g != e]
    if total_time:
        raise ValueError
    def integrand(*xs):
        ts = [None] * E
        ts[e] = eps
        for g, x in zip(others, xs):
            ts[g] = x
        return f(ts)
    lims = []
    for g in others:
        if special is not None and g == special:
            lims.append([0, eps])
        else:
            lims.append([eps, 10 * eps, 1, mp.inf])
    if len(others) == 1:
        return eps * mp.quad(integrand, lims[0])
    if len(others) == 2:
        return eps * mp.quad(integrand, lims[0], lims[1], maxdegree=7)
    raise NotImplementedError


def eps_derivative(f, E, eps, special=None):
    """eps d/deps of the graph value.  Ordinary: -sum_e eps w(t_e=eps).  Special edge s: the derivative also hits the
    upper limit of t_s: + eps w(t_s = eps) - sum_{e != s} eps w(t_e = eps)."""
    tot = 0
    for e in range(E):
        if special is not None and e == special:
            tot += face_value(f, E, e, eps, special=None)     # t_s = eps (upper limit), others in [eps, oo)
        else:
            tot -= face_value(f, E, e, eps, special=special)
    return tot


def total_time_derivative(N, qs, eps):
    """total-time cutoff tau >= eps: w = int_{tau>=eps} tau^{E-1} dtau d alpha sum_j nums_j tau^{-2-j} e^{-tau Phi};
    eps d/deps w = - eps^E int_simplex sum_j nums_j eps^{-2-j} e^{-eps Phi} d alpha"""
    nums, Phi, al = reduced_integrand(N, qs)
    E = len(qs)
    fn = [sp.lambdify(al, nj, 'mpmath') for nj in nums]
    fP = sp.lambdify(al, Phi, 'mpmath')
    def g(*x):
        a = list(x) + [1 - sum(x)]
        if a[-1] < 0:
            return mp.mpf(0)
        return sum(fn[j](*a) * eps ** (-2 - j) for j in range(len(fn))) * mp.e ** (-eps * fP(*a))
    if E == 2:
        val = mp.quad(lambda x: g(x), [0, 1])
    elif E == 3:
        val = mp.quad(lambda x: mp.quad(lambda y: g(x, y), [0, 1 - x]), [0, 1])
    else:
        raise NotImplementedError
    return -eps ** E * val


def fit(epsl, vals, kmax):
    """vals = eps d/deps w at the eps's; model: -sum_k (k a_k) eps^{-k} + b + d1 eps + d2 eps log eps"""
    rows = []
    for e in epsl:
        rows.append([e ** (-k) for k in range(1, kmax + 1)] + [1, e, e * mp.log(e)])
    A = mp.matrix(rows)
    y = mp.matrix(vals)
    sol = mp.lu_solve(A.T * A, A.T * y)
    kak = {k: -sol[k - 1] for k in range(1, kmax + 1)}
    return kak, sol[kmax]


def numeric_value(x):
    return complex(sp.N(x))


def run_case(name, N, qs, special=None, total_time=False):
    E = len(qs)
    if special is None and not total_time:
        ex = T.cycle_singular(N, qs)
    elif total_time:
        ex = T.cycle_singular_totaltime(N, qs)
    else:
        ex = T.special_edge_singular(N, qs, special)
    kmax = max([k for k in ex if k != 'log'] + [0])
    f = make_f(N, qs)
    epsl = [mp.mpf(10) ** (-x) for x in (2.0, 2.5, 3.0, 3.5, 4.0, 4.5)]
    vals = []
    for e in epsl:
        if total_time:
            vals.append(total_time_derivative(N, qs, e))
        else:
            vals.append(eps_derivative(f, E, e, special=special))
    kak, b = fit(epsl, vals, max(kmax, 1))
    rec = dict(name=name, E=E, exact={str(k): str(v) for k, v in ex.items()},
               numeric_b=str(mp.nstr(b, 10)), exact_b=str(sp.N(ex['log'], 12)),
               numeric_k_a_k={str(k): str(mp.nstr(v, 10)) for k, v in kak.items()},
               exact_k_a_k={str(k): str(sp.N(k * ex.get(k, 0), 12)) for k in kak})
    db = abs(complex(b) - numeric_value(ex['log']))
    dk = max([abs(complex(kak[k]) - numeric_value(k * ex.get(k, 0))) for k in kak] + [0])
    scale = 1 + max([abs(numeric_value(v)) for v in ex.values()] + [0])
    rec['max_rel_dev'] = float(max(db, dk) / scale)
    rec['pass'] = rec['max_rel_dev'] < 1e-4
    print(name, 'pass' if rec['pass'] else 'FAIL', rec['max_rel_dev'], rec['exact'])
    return rec


def main():
    p = [sp.Rational(7, 10), sp.Rational(-3, 10), sp.Rational(1, 5), sp.Rational(11, 10)]
    r = [sp.Rational(-2, 5), sp.Rational(1, 2), sp.Rational(3, 10), sp.Rational(-1, 5)]
    zero = [0, 0, 0, 0]
    pr = [a + b for a, b in zip(p, r)]
    l = L
    res = []
    # E = 2 bubbles
    res.append(run_case('E2: N=1', sp.Integer(1), [zero, p]))
    res.append(run_case('E2: N=l0^2+l1 l2', l[0] ** 2 + l[1] * l[2], [zero, p]))
    res.append(run_case('E2: N=(l.p)^2 + 3 l0^4', T.dot(l, p) ** 2 + 3 * l[0] ** 4, [zero, p]))
    res.append(run_case('E2: N=(l.l)^3 (eps^-3, eps^-2, eps^-1, log)', T.dot(l, l) ** 3, [zero, p]))
    # E = 3 triangles
    res.append(run_case('E3: N=l.l', T.dot(l, l), [zero, p, pr]))
    res.append(run_case('E3: N=(l.l)^2 + l0 (l.r)', T.dot(l, l) ** 2 + l[0] * T.dot(l, r), [zero, p, pr]))
    # total-time cutoff: same log, different eps^-k
    res.append(run_case('E2 total-time cutoff: N=(l.l)^2', T.dot(l, l) ** 2, [zero, p], total_time=True))
    res.append(run_case('E3 total-time cutoff: N=(l.l)^2', T.dot(l, l) ** 2, [zero, p, pr], total_time=True))
    # special edge: no log eps
    res.append(run_case('E2 special edge s=0: N=l.l', T.dot(l, l), [zero, p], special=0))
    res.append(run_case('E2 special edge s=1: N=(l.l)^2+l0^2', T.dot(l, l) ** 2 + l[0] ** 2, [zero, p], special=1))
    res.append(run_case('E3 special edge s=1: N=(l.l)^2', T.dot(l, l) ** 2, [zero, p, pr], special=1))
    # comparison per-edge vs total-time exact coefficients (scheme dependence of eps^-k, universality of log)
    comp = []
    for (nm, N, qs) in [('(l.l)^2 bubble', T.dot(l, l) ** 2, [zero, p]), ('(l.l)^2 triangle', T.dot(l, l) ** 2, [zero, p, pr])]:
        a = T.cycle_singular(N, qs); b = T.cycle_singular_totaltime(N, qs)
        comp.append(dict(case=nm, per_edge={str(k): str(v) for k, v in a.items()},
                         total_time={str(k): str(v) for k, v in b.items()},
                         log_equal=bool(sp.simplify(a['log'] - b['log']) == 0)))
        print('scheme comparison', nm, comp[-1]['log_equal'], a, b)
    # Gaussian moment check: int d^4l/(2pi)^4 l0^4 l1^2 e^{-tau l^2} = (3/(4 tau^2))(1/(2 tau)) /(16 pi^2 tau^2)
    tau = sp.Symbol('tau', positive=True)
    N = l[0] ** 4 * l[1] ** 2
    cur = N; tot = 0
    for j in range(4):
        tot += cur.subs({x: 0 for x in l}) / (math.factorial(j) * (4 * tau) ** j)
        cur = T.lap(cur)
    gauss_ok = sp.simplify(tot - sp.Rational(3, 4) / tau ** 2 * sp.Rational(1, 2) / tau) == 0
    print('Gaussian moments ok:', gauss_ok)
    out = dict(description=__doc__, cases=res, scheme_comparison=comp, gaussian_moments_ok=bool(gauss_ok),
               all_pass=bool(all(r_['pass'] for r_ in res) and gauss_ok and all(c['log_equal'] for c in comp)))
    with open(os.path.join(HERE, 'checks_simplex_results.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print('ALL PASS' if out['all_pass'] else 'SOME FAILED')


if __name__ == '__main__':
    main()
