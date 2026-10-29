"""Checks of the simplex formulas of twoletter_theory.tex (Theorem thm:tl-simplex, Proposition prop:tl-shape,
Lemma lem:tl-special, Proposition prop:tl-genfun, Proposition prop:tl-bubble).

(1) Residual test (replaces the least-squares fit of the interrupted run, which was numerically singular for the
    eps^{-3} case).  For a cycle with E edges, shifts q_e, numerator N(l), the exact formulas give a_k (eps^{-k}) and
    b (log eps).  By Proposition c4-boundary,  eps d/deps w = - sum_e eps w(t_e = eps); the residual
        r(eps) = eps d/deps w + sum_k k a_k eps^{-k} - b
    must be O(eps log eps).  It is evaluated with mpmath (40 digits for E = 2, 25 for E = 3) at eps = 1e-2, 1e-3, 1e-4
    and required to decrease and to be small at the smallest eps.  The same for the total-time cutoff tau >= eps and
    for one special edge t_s in [0, eps].
(2) Special-edge finite part: exact toy (hand computation in the tex) and a numerical eps -> 0 extrapolation-free
    test  w(eps) - sum_k a_k eps^{-k} -> c_0.
(3) Generating functions: F_m for the numerator e^{i(sigma,l)} equals the Bessel closed form, order by order in sigma,
    for E = 1, 2, 3 (exact).
(4) Closed form of the scalar bubble log kernel in the derivative numbers (exact).
(5) Gaussian moments; per-edge vs total-time: equal log, different eps^{-k}.
Writes checks_simplex_results.json."""
import json, os, math
import mpmath as mp
import sympy as sp
import tl_simplex as T

HERE = os.path.dirname(os.path.abspath(__file__))
L = T.L


def reduced_integrand(N, qs):
    nums, Phi, al = T._reduced(N, qs)
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


def face_value(f, E, e, eps, special=None):
    """eps * (integral with t_e = eps, other t's in [eps, oo), or [0, eps] for the special edge)"""
    others = [g for g in range(E) if g != e]

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
    return eps * mp.quad(integrand, lims[0], lims[1], maxdegree=8)


def eps_derivative(f, E, eps, special=None):
    tot = 0
    for e in range(E):
        if special is not None and e == special:
            tot += face_value(f, E, e, eps, special=None)
        else:
            tot -= face_value(f, E, e, eps, special=special)
    return tot


def total_time_derivative(N, qs, eps):
    nums, Phi, al = reduced_integrand(N, qs)
    E = len(qs)
    fn = [sp.lambdify(al, nj, 'mpmath') for nj in nums]
    fP = sp.lambdify(al, Phi, 'mpmath')

    def g(*x):
        a = list(x) + [1 - sum(x)]
        return sum(fn[j](*a) * eps ** (-2 - j) for j in range(len(fn))) * mp.e ** (-eps * fP(*a))
    if E == 2:
        val = mp.quad(lambda x: g(x), [0, 1])
    else:
        val = mp.quad(lambda x: mp.quad(lambda y: g(x, y), [0, 1 - x]), [0, 1])
    return -eps ** E * val


def cnum(x):
    return mp.mpc(complex(sp.N(x, 30)))


def run_residual(name, N, qs, special=None, total_time=False):
    E = len(qs)
    mp.mp.dps = 40 if E <= 2 else 25
    if total_time:
        ex = T.cycle_singular_totaltime(N, qs)
    elif special is not None:
        ex = T.special_edge_singular(N, qs, special)
    else:
        ex = T.cycle_singular(N, qs)
    ks = [k for k in ex if k != 'log']
    f = make_f(N, qs)
    epss = [mp.mpf(10) ** (-2), mp.mpf(10) ** (-3), mp.mpf(10) ** (-4)] if E <= 2 else \
        [mp.mpf(10) ** (-2), mp.mpf('10') ** (-2.5), mp.mpf(10) ** (-3)]
    res = []
    for eps in epss:
        if total_time:
            v = total_time_derivative(N, qs, eps)
        else:
            v = eps_derivative(f, E, eps, special=special)
        pred = -sum(k * cnum(ex[k]) * eps ** (-k) for k in ks) + cnum(ex['log'])
        res.append(abs(v - pred))
    scale = 1 + max([abs(cnum(v)) for v in ex.values()] + [0])
    rel = [float(r / scale) for r in res]
    ok = rel[-1] < (5e-3 if E <= 2 else 3e-2) and rel[-1] < rel[0]
    rec = dict(name=name, E=E, exact={str(k): str(v) for k, v in ex.items()}, residuals_rel=rel,
               eps=[float(e) for e in epss], pass_=bool(ok))
    print(('PASS ' if ok else 'FAIL ') + name, ['%.2e' % r for r in rel], rec['exact'])
    return rec


def run_special_finite_numeric(name, N, qs, s):
    """w(eps) - sum_k a_k eps^{-k} -> c_0, E = 2, direct 2d quadrature"""
    mp.mp.dps = 30
    ex = T.special_edge_singular(N, qs, s)
    c0 = T.special_edge_finite(N, qs, s)
    f = make_f(N, qs)
    E = len(qs)
    other = [g for g in range(E) if g != s][0]
    vals = []
    for eps in [mp.mpf(10) ** (-3), mp.mpf(10) ** (-4), mp.mpf(10) ** (-5)]:
        def integ(ts_, t_):
            tt = [None, None]
            tt[s] = ts_; tt[other] = t_
            return f(tt)
        w = mp.quad(lambda ts_: mp.quad(lambda t_: integ(ts_, t_), [eps, 10 * eps, 1, mp.inf]), [0, eps])
        rem = w - sum(cnum(ex[k]) * eps ** (-k) for k in ex if k != 'log')
        vals.append(rem)
    dev = [float(abs(v - cnum(c0))) for v in vals]
    ok = dev[-1] < 1e-3 * (1 + float(abs(cnum(c0)))) and dev[-1] < dev[0]
    print(('PASS ' if ok else 'FAIL ') + name, 'c0 =', c0, 'deviations', ['%.2e' % d for d in dev])
    return dict(name=name, c0=str(c0), eps_coeffs={str(k): str(v) for k, v in ex.items()}, deviations=dev, pass_=bool(ok))


def check_genfun(E, qs, sigma0, order=6):
    lam = sp.Symbol('lam')
    sigma = [lam * x for x in sigma0]
    Nexp = sum((sp.I * T.dot(sigma, L)) ** n / sp.factorial(n) for n in range(order + 1))
    ok = True
    for m in range(0, 3):
        Fd, al = T.F_m(Nexp, qs, m)
        Fd = sp.expand(sum(sp.expand(Fd).coeff(lam, n) * lam ** n for n in range(order + 1)))
        Fg, al2 = T.genfun_F(m, qs, sigma, order)
        Fg = Fg.subs(dict(zip(al2, al)))
        diff = sp.expand(Fd - Fg)
        if diff != 0:
            ok = False
            print('genfun mismatch E=%d m=%d' % (E, m), sp.factor(diff)[:200] if hasattr(diff, '__len__') else diff)
    print(('PASS ' if ok else 'FAIL ') + 'generating function F_m, E=%d, orders <= %d' % (E, order))
    return ok


def check_bubble_closed(order=6):
    """b^gen(s,u) = - e^{i(u,R)} int_0^1 e^{-i alpha (sigma,R)} I_0(|sigma||R| sqrt(alpha(1-alpha))) d alpha
    = - sum_{a,b,j} (i(s,R))^a (i(u,R))^b |u-s|^{2j} (|R|^2/4)^j (a+j)!(b+j)!/(a!b!(j!)^2(a+b+2j+1)!)  (times 1/16pi^2)
    compared with the direct simplex formula for the bubble with edges l (q_0 = 0) and l + R (q_1 = -R)."""
    lam = sp.Symbol('lam')
    s0 = [sp.Rational(1, 2), sp.Rational(-1, 3), 0, sp.Rational(2, 5)]
    u0 = [sp.Rational(-1, 4), sp.Rational(1, 5), sp.Rational(1, 2), 0]
    R = [sp.Rational(3, 7), sp.Rational(-1, 2), sp.Rational(1, 3), sp.Rational(1, 1)]
    s = [lam * x for x in s0]; u = [lam * x for x in u0]
    sig = [a - b for a, b in zip(u, s)]
    qs = [[0, 0, 0, 0], [-x for x in R]]
    Nexp = sum((sp.I * T.dot(sig, L)) ** n / sp.factorial(n) for n in range(order + 1))
    pref = sum((sp.I * T.dot(u, R)) ** n / sp.factorial(n) for n in range(order + 1))
    direct = T.cycle_singular(Nexp, qs, kmax=0)['log']
    direct = sp.expand(direct * pref)
    direct = sp.expand(sum(direct.coeff(lam, n) * lam ** n for n in range(order + 1)))
    closed = 0
    R2 = T.dot(R, R)
    for a in range(order + 1):
        for b in range(order + 1 - a):
            for j in range((order - a - b) // 2 + 1):
                closed += ((sp.I * T.dot(s, R)) ** a * (sp.I * T.dot(u, R)) ** b * T.dot(sig, sig) ** j * (R2 / 4) ** j
                           * T.bubble_log_closed(a, b, j) / (sp.factorial(a) * sp.factorial(b))) * \
                    (sp.factorial(a) * sp.factorial(b))
    # note: the closed form already contains 1/(a! b!) inside bubble_log_closed
    closed = sp.expand(-closed)
    closed = sp.expand(sum(closed.coeff(lam, n) * lam ** n for n in range(order + 1)))
    ok = sp.expand(direct - closed) == 0
    print(('PASS ' if ok else 'FAIL ') + 'closed form of the scalar bubble log kernel, orders <= %d' % order)
    return ok


def main():
    p = [sp.Rational(7, 10), sp.Rational(-3, 10), sp.Rational(1, 5), sp.Rational(11, 10)]
    r = [sp.Rational(-2, 5), sp.Rational(1, 2), sp.Rational(3, 10), sp.Rational(-1, 5)]
    zero = [0, 0, 0, 0]
    pr = [a + b for a, b in zip(p, r)]
    l = L
    out = dict(description=__doc__)
    res = []
    res.append(run_residual('E2: N=1', sp.Integer(1), [zero, p]))
    res.append(run_residual('E2: N=l0^2+l1 l2', l[0] ** 2 + l[1] * l[2], [zero, p]))
    res.append(run_residual('E2: N=(l.p)^2 + 3 l0^4', T.dot(l, p) ** 2 + 3 * l[0] ** 4, [zero, p]))
    res.append(run_residual('E2: N=(l.l)^3 (eps^-3, eps^-2, eps^-1, log)', T.dot(l, l) ** 3, [zero, p]))
    res.append(run_residual('E3: N=l.l', T.dot(l, l), [zero, p, pr]))
    res.append(run_residual('E3: N=(l.l)^2 + l0 (l.r)', T.dot(l, l) ** 2 + l[0] * T.dot(l, r), [zero, p, pr]))
    res.append(run_residual('E2 total-time cutoff: N=(l.l)^2', T.dot(l, l) ** 2, [zero, p], total_time=True))
    res.append(run_residual('E3 total-time cutoff: N=(l.l)^2', T.dot(l, l) ** 2, [zero, p, pr], total_time=True))
    res.append(run_residual('E2 special edge s=0: N=l.l', T.dot(l, l), [zero, p], special=0))
    res.append(run_residual('E2 special edge s=1: N=(l.l)^2+l0^2', T.dot(l, l) ** 2 + l[0] ** 2, [zero, p], special=1))
    res.append(run_residual('E3 special edge s=1: N=(l.l)^2', T.dot(l, l) ** 2, [zero, p, pr], special=1))
    out['residual_tests'] = res
    # (2) special-edge finite part
    toy = T.special_edge_singular(-T.dot(l, l), [zero, p], 0)
    toy_c0 = T.special_edge_finite(-T.dot(l, l), [zero, p], 0)
    p2 = T.dot(p, p)
    toy_ok = sp.simplify(toy[1] + sp.Rational(1, 2)) == 0 and sp.simplify(toy_c0 + sp.Rational(3, 8) * p2) == 0
    print(('PASS ' if toy_ok else 'FAIL ') + 'special-edge toy: eps^-1 =', toy[1], ' c0 =', toy_c0,
          ' expected -1/2 and', -sp.Rational(3, 8) * p2)
    out['special_toy'] = dict(eps_m1=str(toy[1]), c0=str(toy_c0), expected=['-1/2', str(-sp.Rational(3, 8) * p2)],
                              pass_=bool(toy_ok))
    sf = []
    sf.append(run_special_finite_numeric('special finite E2: N=-(l.l)(1 + l.r)', -T.dot(l, l) * (1 + T.dot(l, r)),
                                         [zero, p], 0))
    sf.append(run_special_finite_numeric('special finite E2: N=-|l-p|^2 (l0^2 + l1 l3)',
                                         -T.dot([a - b for a, b in zip(l, p)], [a - b for a, b in zip(l, p)])
                                         * (l[0] ** 2 + l[1] * l[3]), [zero, p], 1))
    out['special_finite'] = sf
    # (3) generating functions
    sig0 = [sp.Rational(1, 3), sp.Rational(-1, 2), sp.Rational(1, 4), sp.Rational(2, 3)]
    g1 = check_genfun(1, [zero], sig0)
    g2 = check_genfun(2, [zero, p], sig0)
    g3 = check_genfun(3, [zero, p, pr], sig0, order=5)
    out['genfun'] = dict(E1=g1, E2=g2, E3=g3)
    # (4) bubble closed form
    out['bubble_closed_form'] = check_bubble_closed()
    # (5) scheme comparison and Gaussian moments
    comp = []
    for (nm, N, qs) in [('(l.l)^2 bubble', T.dot(l, l) ** 2, [zero, p]),
                        ('(l.l)^2 triangle', T.dot(l, l) ** 2, [zero, p, pr])]:
        a = T.cycle_singular(N, qs); b = T.cycle_singular_totaltime(N, qs)
        comp.append(dict(case=nm, per_edge={str(k): str(v) for k, v in a.items()},
                         total_time={str(k): str(v) for k, v in b.items()},
                         log_equal=bool(sp.simplify(a['log'] - b['log']) == 0),
                         eps_m1_differ=bool(sp.simplify(a.get(1, 0) - b.get(1, 0)) != 0)))
        print('scheme comparison', nm, 'log equal:', comp[-1]['log_equal'], ' eps^-1 differ:', comp[-1]['eps_m1_differ'])
    out['scheme_comparison'] = comp
    tau = sp.Symbol('tau', positive=True)
    N = l[0] ** 4 * l[1] ** 2
    cur = N; tot = 0
    for j in range(4):
        tot += cur.subs({x: 0 for x in l}) / (math.factorial(j) * (4 * tau) ** j)
        cur = T.lap(cur)
    gauss_ok = sp.simplify(tot - sp.Rational(3, 4) / tau ** 2 * sp.Rational(1, 2) / tau) == 0
    out['gaussian_moments_ok'] = bool(gauss_ok)
    allp = (all(r_['pass_'] for r_ in res) and toy_ok and all(s_['pass_'] for s_ in sf) and g1 and g2 and g3
            and out['bubble_closed_form'] and gauss_ok and all(c['log_equal'] and c['eps_m1_differ'] for c in comp))
    out['all_pass'] = bool(allp)
    with open(os.path.join(HERE, 'checks_simplex_results.json'), 'w') as fh:
        json.dump(out, fh, indent=1, default=str)
    print('ALL PASS' if allp else 'SOME FAILED')


if __name__ == '__main__':
    main()
