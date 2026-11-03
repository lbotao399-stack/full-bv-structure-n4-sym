"""Route B numerical eps-extrapolation of single defect-cycle configurations (gates G-num, G-noLog numeric, G-cubic).

A configuration (vertex sequence with S^(1)(x) placed, special edge sigma) of a DJob at rational numerical data
(polarisations, letter components, xi, real Euclidean output momenta) is evaluated at finite eps WITHOUT expansion:
    value(eps) = (16 pi^2)^{-1} int_{alpha_sigma smallest} dalpha sum_j num_j(alpha) int_{tau_lo}^{tau_hi} tau^{n_t-3-j} e^{-tau Phi} dtau,
tau_lo = eps / min_{f != sigma} alpha_f, tau_hi = eps / alpha_sigma (special time in [0, eps], ordinary times >= eps,
T = infinity), the tau integral as an incomplete gamma function and the alpha integral by mpmath tanh-sinh quadrature at
22 digits on the ordered subregions; then fitted at eps = eps_0 2^{-j} to
    a_1/eps + a_0 + c_L log eps + eps (b_1 + b_2 log eps) + b_3 eps^2 + b_4 eps^3
and compared with the exact Feynman-parameter coefficients (ce_int.J) of the same configuration.  c_L must vanish."""
import os, sys, json, time, math, random
from fractions import Fraction as Fr
import mpmath as mp
import sympy as sp
from sympy.polys.domains import QQ_I
import qb_core as B
import qb_job as QJ

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
mp.mp.dps = 22


def Gup(sv, x):
    """upper incomplete gamma Gamma(s, x) for integer s (x > 0, possibly infinite)."""
    if x == mp.inf:
        return mp.mpf(0)
    if sv > 0:
        return mp.gammainc(sv, x)
    return x ** sv * mp.expint(1 - sv, x)


def tau_int(kpow, Phi, tlo, thi):
    sv = kpow + 1
    if Phi == 0 or abs(Phi * tlo) < mp.mpf(10) ** (-mp.mp.dps):
        if thi == mp.inf:
            return mp.mpf(0)          # measure-zero point of the alpha integral
        return mp.log(thi / tlo) if sv == 0 else (thi ** sv - tlo ** sv) / sv
    return Phi ** (-sv) * (Gup(sv, Phi * tlo) - Gup(sv, Phi * thi))


def q2mp(v):
    return mp.mpf(int(v.numerator)) / int(v.denominator)


def numeric_images(job, vals):
    """substitute all generators except q and al by the given rational values (dict name -> Fraction/complex)."""
    R = job.R
    imgs = []
    for nm, g in zip(job.names, R.gens):
        if nm in vals:
            v = vals[nm]
            imgs.append((g, R(QQ_I.from_sympy(sp.nsimplify(v)))))
    return imgs


def config_data(job, seq, sigma, vals):
    types, curs, words = job.config(seq, sigma)
    imgs = numeric_images(job, vals)
    R = job.R
    Num = R.zero
    for w, P in words.items():
        Num += P
    Num = Num.compose(imgs) if imgs else Num
    curs_n = [[c.compose(imgs) if c else R.zero for c in cu] for cu in curs]
    return types, curs_n, Num


def exact_coeffs(job, types, curs, Num, mlist=(0, 1, 2)):
    co = job.integrate(Num, curs, types, mlist)
    out = {}
    for m, d in co.items():
        tot = mp.mpc(0)
        for b, P in d.items():
            assert all(sum(mon) == 0 for mon in P.monoms()), 'non-numeric generator left'
            c = P.coeff(1) if P else 0
            cv = mp.mpc(q2mp(c.x), q2mp(c.y))
            tot += cv * (1 if b == '1' else mp.log(int(b[3:])))
        out[m] = tot
    return out


def numeric_value(job, types, curs, Num, eps):
    """finite-eps value (units hbar/(16 pi^2)), n_t = 2 or 3."""
    R = job.R
    q = job.q
    nt = len(types)
    js = types.index('S')
    order = [js] + [j for j in range(nt) if j != js]
    ks = [q[m] - curs[js][m] for m in range(4)]
    Nt = Num * (-sum((x * x for x in ks), R.zero))
    # Laplacians in q
    laps = [Nt]
    while laps[-1]:
        D = laps[-1]
        laps.append(sum((D.diff(q[m]).diff(q[m]) for m in range(4)), R.zero))
    laps = laps[:-1]
    qi = job.q_idx

    def tonum(P):
        terms = []
        for mon, cf in P.items():
            cv = mp.mpc(q2mp(cf.x), q2mp(cf.y))
            qe = tuple(mon[qi[m]] for m in range(4))
            assert sum(mon) == sum(qe), 'non-numeric generator left'
            terms.append((qe, cv))
        return terms
    lapn = [tonum(D) for D in laps]

    def cnum(c):
        if not c:
            return mp.mpf(0)
        assert all(sum(mon) == 0 for mon in c.monoms())
        v = c.coeff(1)
        return mp.mpc(q2mp(v.x), q2mp(v.y))
    cursn = [[cnum(c) for c in cu] for cu in curs]

    def evalp(terms, x):
        s_ = mp.mpc(0)
        for qe, cv in terms:
            v = cv
            for m in range(4):
                if qe[m]:
                    v *= x[m] ** qe[m]
            s_ += v
        return s_

    def integrand(alpha):        # alpha: list over parameters in 'order'
        cent = [sum(alpha[p] * cursn[order[p]][m] for p in range(nt)) for m in range(4)]
        Phi = sum(alpha[p] * sum(cursn[order[p]][m] ** 2 for m in range(4)) for p in range(nt)) - sum(c * c for c in cent)
        Phi = mp.re(Phi) if abs(mp.im(Phi)) < mp.mpf(10) ** (-30) else Phi
        tlo = eps / min(alpha[1:]); thi = eps / alpha[0]
        tot = mp.mpc(0)
        for j, terms in enumerate(lapn):
            nj = evalp(terms, cent) / (math.factorial(j) * 4 ** j)
            if nj == 0:
                continue
            kpow = nt - 3 - j
            if abs(Phi) < mp.mpf(10) ** (-40):
                if kpow == -1:
                    ti = mp.log(thi / tlo)
                else:
                    ti = (thi ** (kpow + 1) - tlo ** (kpow + 1)) / (kpow + 1)
            else:
                ti = tau_int(kpow, Phi, tlo, thi)
            tot += nj * ti
        return tot

    if nt == 2:
        f = lambda a: integrand([a, 1 - a])
        val = mp.quad(f, [0, mp.mpf(1) / 2])
    elif nt == 3:
        # alpha_s = a, alpha_1 = b, alpha_2 = 1 - a - b; region a <= b, a <= 1 - a - b; split b <= 1-a-b and b >= 1-a-b
        def f1(a, b):
            return integrand([a, b, 1 - a - b])
        v1 = mp.quad(lambda a: mp.quad(lambda b: f1(a, b), [a, (1 - a) / 2]), [0, mp.mpf(1) / 3])
        v2 = mp.quad(lambda a: mp.quad(lambda b: f1(a, b), [(1 - a) / 2, 1 - 2 * a]), [0, mp.mpf(1) / 3])
        val = v1 + v2
    else:
        raise NotImplementedError('n_t > 3')
    return val


def fit(epss, vals):
    """least squares (QR) on 1/x, 1, log x, x, x log x, x^2, x^3 with x = eps/eps_max (well conditioned), converted
    back to a_1/eps + a_0 + c_L log eps + ...: a_1 = c_{1/x} eps_max, c_L = c_{log}, a_0 = c_1 - c_L log eps_max."""
    e0 = max(epss)

    def row(e):
        x = e / e0
        return [1 / x, mp.mpf(1), mp.log(x), x, x * mp.log(x), x ** 2, x ** 3]
    A = mp.matrix([row(e) for e in epss])
    bvec = mp.matrix(vals)
    sol, res = mp.qr_solve(A, bvec)
    cL = sol[2]
    return dict(a_m1=sol[0] * e0, a_0=sol[1] - cL * mp.log(e0), c_log=cL, b1=sol[3] / e0, b2=sol[4] / e0)


def check(job, seq, sigma, vals, eps0=mp.mpf(2) ** -10, npts=9):
    types, curs, Num = config_data(job, seq, sigma, vals)
    if not Num:
        return None
    ex = exact_coeffs(job, types, curs, Num)
    epss = [eps0 / mp.mpf(2) ** j for j in range(npts)]
    vs = [numeric_value(job, types, curs, Num, e) for e in epss]
    fr = fit(epss, vs)
    return dict(types=types, exact_eps0=mp.nstr(ex.get(0, 0), 25), exact_eps_m1=mp.nstr(ex.get(1, 0), 25),
                fit_eps0=mp.nstr(fr['a_0'], 25), fit_eps_m1=mp.nstr(fr['a_m1'], 25), fit_log_eps=mp.nstr(fr['c_log'], 6),
                diff_eps0=mp.nstr(abs(fr['a_0'] - ex.get(0, 0)), 5), diff_eps_m1=mp.nstr(abs(fr['a_m1'] - ex.get(1, 0)), 5),
                eps=[mp.nstr(e, 6) for e in epss])


def rnd_vals(job, seed):
    random.seed(seed)
    vals = {}
    for nm in job.names:
        if nm.startswith(('q', 'al')) and nm[1:].replace('l', '').isdigit():
            continue
        if nm.startswith(('o', 'za_', 'zb_', 'xi')):
            vals[nm] = Fr(random.randint(-3, 3), random.randint(1, 3))
        elif nm.startswith('r'):
            vals[nm] = Fr(random.randint(-4, 4), random.randint(2, 4))
        elif nm.startswith('b'):
            vals[nm] = Fr(random.randint(-2, 2), 1)
    return vals


def run():
    t0 = time.time()
    xi = B.XI[:8] + [0] * 8
    cases = [('cubic example (AW,lamP)->AW, T5 at the special edge', ('AW', 'lamP', ['AW'], 'cubic')),
             ('(AW,AW)->lamM, T1/T3 at the special edge', ('AW', 'AW', ['lamM'], 'quad')),
             ('(lamP,c)->c, n=2 chains', ('lamP', 'c', ['c'], 'all')),
             ('(AV,lamP)->AV, k_I = 0 placements', ('AV', 'lamP', ['AV'], 'all'))]
    out = {}; gates = []
    for name, (sa, sb, outs, terms) in cases:
        job = QJ.DJob(sa, sb, outs, 'xi', xi, terms=terms, tag='num')
        vals = rnd_vals(job, 7)
        rows = []
        for seq in job.sequences():
            for sigma in range(len(seq) + 1):
                if len(seq) + 1 > 3:
                    continue
                if len(seq) + 1 == 3 and sum(1 for r in rows if len(r['types']) == 3) >= 1:
                    continue
                r = check(job, seq, sigma, vals)
                if r is None:
                    continue
                r['seq'] = str(seq); r['sigma'] = sigma
                rows.append(r)
                print(name, seq, sigma, 'eps0 exact', r['exact_eps0'], 'fit', r['fit_eps0'], 'log', r['fit_log_eps'], flush=True)
        out[name] = dict(values={k: str(v) for k, v in vals.items()}, rows=rows)
        def _c(x):
            return abs(complex(str(x).replace(' ', '').replace('(', '').replace(')', '')))
        ok = bool(rows) and all(_c(r['diff_eps0']) < 1e-8 and _c(r['diff_eps_m1']) < 1e-12 and _c(r['fit_log_eps']) < 1e-8 for r in rows)
        gates.append(dict(id=f'G-num {name}', pass_=bool(ok), value=[(r['seq'], r['sigma'], r['diff_eps0'], r['diff_eps_m1'], r['fit_log_eps']) for r in rows],
                          target='numerical eps-extrapolation (22-digit quadrature, 9 points eps = 2^-10 .. 2^-18) = exact eps^0 and eps^-1 coefficients; log eps coefficient 0'))
        json.dump(dict(results=out, gates=gates), open(os.path.join(OUT, 'numeps.json'), 'w'), indent=1, default=str)
    out['seconds'] = round(time.time() - t0, 1)
    return out, gates


if __name__ == '__main__':
    out, gates = run()
    json.dump(dict(results=out, gates=gates), open(os.path.join(OUT, 'numeps.json'), 'w'), indent=1, default=str)
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL')
