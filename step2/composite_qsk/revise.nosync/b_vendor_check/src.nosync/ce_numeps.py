"""Numerical eps-extrapolation (route B check of GC5): the special-edge composite chains with k = 1 vertex evaluated at
finite eps without any expansion, at real Euclidean kinematic points:
    value(eps) = (16 pi^2)^{-1} int_simplex dalpha sum_j num_j(alpha) int_{tau_lo}^{tau_hi} tau^{n_t-3-j} e^{-tau Phi} dtau,
tau_lo = eps max_f(1/alpha_f), tau_hi = eps/alpha_sigma (incomplete gamma, exact), the alpha integral by mpmath
quadrature on R_sigma; fitted to z_2/eps^2 + z_1/eps + z_0 + eps (a + b log eps) and compared with the exact
coefficients of the Feynman-parameter evaluation (ce_int)."""
import os, sys, json, time, math, random
import sympy as sp
import mpmath as mp
from sympy.polys.domains import QQ_I
import ce_job as CJ

I = sp.I


def num_poly(job, P, subs_num):
    """ring element -> complex mpmath number after substituting all generators except q (dict name -> value)."""
    pass


def chain_values(job, seq, c, s, avec, eps_list, dps=20):
    """Pb-mode configuration (seq, c, s) with numeric outputs: exact coefficients and finite-eps values (summed over
    colour words) with the line direction a = avec."""
    R = job.R
    types, curs, words = job.config(seq, c, s, None, 'Pb')
    out = {}
    for w, Num in words.items():
        # substitute a numerically
        imgs = []
        for nm, gen in zip(job.names, R.gens):
            if nm.startswith('a') and nm[1:].isdigit():
                imgs.append((gen, R(QQ_I.from_sympy(sp.nsimplify(avec[int(nm[1:])])))))
        Na = Num.compose(imgs)
        ex = job.integrate(Na, curs, types, (0, 1, 2))
        exact = {m: {b: complex(sp.N(P.as_expr())) if P else 0 for b, P in d.items()} for m, d in ex.items()}
        exact_f = {m: sum(v * (math.log(int(b[3:])) if b != '1' else 1.0) for b, v in d.items()) for m, d in exact.items()}
        # numerical: parameters (special first) as in integrate
        params = []
        for j, t in enumerate(types):
            if t == 'S':
                params.append((j, True))
            elif t == 'K':
                params += [(j, False), (j, False)]
            else:
                params.append((j, False))
        sidx = [i for i, p in enumerate(params) if p[1]][0]
        params = [params[sidx]] + [p for i, p in enumerate(params) if i != sidx]
        nt = len(params)
        assert nt == 3, 'numeric check implemented for k = 1'
        q = job.q
        js = types.index('S')
        ks = [q[m] - curs[js][m] for m in range(4)]
        Nt = Na * (-sum((x * x for x in ks), R.zero)) * R(QQ_I(sp.Rational(1, 4), 0))
        # Laplacians in q
        laps = [Nt]
        while laps[-1]:
            D = laps[-1]
            laps.append(sum((D.diff(q[m]).diff(q[m]) for m in range(4)), R.zero))
        laps = laps[:-1]
        qn = [str(x) for x in q]

        def tonum(P):
            terms = []
            for mon, cf in P.items():
                cv = complex(float(cf.x), float(cf.y))
                qe = tuple(mon[job.q_idx[m]] for m in range(4))
                assert sum(mon) == sum(qe), 'non-numeric generator left'
                terms.append((qe, cv))
            return terms
        lapn = [tonum(D) for D in laps]
        cursn = [[complex(sp.N(cc.as_expr())) if cc else 0j for cc in cu] for cu in curs]

        def evalpoly(terms, x):
            s_ = 0j
            for qe, cv in terms:
                v = cv
                for m in range(4):
                    if qe[m]:
                        v *= x[m] ** qe[m]
                s_ += v
            return s_
        mp.mp.dps = dps

        def value(eps):
            eps = mp.mpf(eps)

            def f(beta, gam):
                al = [beta, gam, 1 - beta - gam]
                cen = [sum(al[p] * cursn[params[p][0]][m] for p in range(nt)) for m in range(4)]
                Phi = sum(al[p] * sum(cursn[params[p][0]][m] ** 2 for m in range(4)) for p in range(nt)) - sum(x * x for x in cen)
                Phi = mp.mpf(Phi.real)
                tlo = eps * max(1 / al[1], 1 / al[2]); thi = eps / al[0]
                tot = mp.mpc(0)
                for j, terms in enumerate(lapn):
                    nj = evalpoly(terms, cen) / (math.factorial(j) * 4 ** j)
                    sexp = nt - 3 - j
                    G = (mp.gammainc(sexp + 1, Phi * tlo) - mp.gammainc(sexp + 1, Phi * thi)) * Phi ** (-(sexp + 1))
                    tot += nj * G
                return tot
            return mp.quad(lambda b_: mp.quad(lambda g_: f(b_, g_), [b_, (1 - b_) / 2, 1 - 2 * b_]), [0, mp.mpf(1) / 3])
        vals = [complex(value(e)) for e in eps_list]
        out[w] = dict(exact=exact_f, values=vals)
    return out


def fit(eps_list, vals, z2, z1):
    """subtract the exact poles, fit z0 + eps (a + b log eps) by least squares (numpy)."""
    import numpy as np
    y = np.array([v - z2 / e ** 2 - z1 / e for v, e in zip(vals, eps_list)])
    A = np.array([[1, e, e * math.log(e)] for e in eps_list], dtype=complex)
    sol, *_ = np.linalg.lstsq(A, y, rcond=None)
    # also a fit with a log eps term to test its absence
    A2 = np.array([[1, math.log(e), e, e * math.log(e)] for e in eps_list], dtype=complex)
    y2 = np.array(vals) - np.array([z2 / e ** 2 + z1 / e for e in eps_list])
    sol2, *_ = np.linalg.lstsq(A2, y2, rcond=None)
    return complex(sol[0]), complex(sol2[1])


def run(npoints=3, seed=7):
    random.seed(seed)
    eps_list = [1e-2, 5e-3, 2e-3, 1e-3, 5e-4]
    report = []
    cases = [('AV', 'AW', ['AW'], {0: 1, 1: 2, 2: -1, 3: 1}, {4: 1, 5: -2, 6: 1, 7: 0, 8: 1, 9: 3}),
             ('AW', 'AV', ['AW'], {4: 2, 5: 1, 6: 0, 7: -1, 8: 1, 9: 1}, {0: 1, 1: 0, 2: 2, 3: -1}),
             ('AV', 'AV', ['AV'], {0: 1, 1: -1, 2: 2, 3: 1}, {0: 2, 1: 1, 2: 0, 3: -1})]
    for (sa, sb, outs, za, zb) in cases:
        for pt in range(npoints):
            p0 = [sp.Rational(random.randint(-9, 9), 7) for _ in range(4)]
            spc = outs[0]
            comps = range(4, 10) if spc == 'AW' else range(4)
            Y = {cc: random.randint(-3, 3) for cc in comps}
            avec = [sp.Rational(random.randint(-5, 5), 3) for _ in range(4)]
            job = CJ.CJob(sa, sb, outs, pol=[Y], za=za, zb=zb)
            # numeric output momentum: substitute r0_m via compose on the config results -> simpler: build job with
            # momenta as generators and substitute before integrating
            R = job.R
            imgs = [(job.g[f'r0_{m}'], R(QQ_I.from_sympy(p0[m]))) for m in range(4)]
            tz = {0: 0, 1: 0, 2: 0}; fs = []
            for seq in job.sequences():
                k = len(seq)
                for c in range(k + 1):
                    for s in range(k + 1):
                        if s == c:
                            continue
                        types, curs, words = job.config(seq, c, s, None, 'Pb')
                        if not words:
                            continue
                        # substitute momenta into numerators and curs
                        curs_n = [[cc.compose(imgs) for cc in cu] for cu in curs]
                        job2 = job
                        for w, Num in words.items():
                            Nn = Num.compose(imgs)
                            res = _single(job2, Nn, curs_n, types, avec, eps_list)
                            for m in (0, 1, 2):
                                tz[m] += res['exact'].get(m, 0)
                            fs.append(res['f'])
            tot_vals = []
            for e in eps_list:
                F = lambda b_, g_, ee=e: sum(ff(b_, g_, ee) for ff in fs)
                tot_vals.append(complex(gl_integral(F, e)))
            z0fit, logc = fit(eps_list, tot_vals, tz[2], tz[1])
            report.append(dict(pair=f'{sa}|{sb}', outs=outs, point=pt, p0=[str(x) for x in p0], a=[str(x) for x in avec],
                               exact_z0=str(tz[0]), exact_z1=str(tz[1]), exact_z2=str(tz[2]), fit_z0=str(z0fit),
                               fit_log_eps_coefficient=str(logc), values=[str(v) for v in tot_vals],
                               rel_err=abs(z0fit - tz[0]) / max(abs(tz[0]), 1e-30)))
            print(report[-1]['pair'], pt, 'exact', tz[0], 'fit', z0fit, 'log coeff', logc, 'poles', tz[1], tz[2], flush=True)
    return report


def _single(job, Nn, curs_n, types, avec, eps_list, dps=15):
    """exact coefficients and a vectorised finite-eps integrand (numpy/scipy, double precision)."""
    import numpy as np
    from scipy import special as SPF
    R = job.R
    imgs = []
    for nm, gen in zip(job.names, R.gens):
        if nm.startswith('a') and nm[1:].isdigit():
            imgs.append((gen, R(QQ_I.from_sympy(sp.nsimplify(avec[int(nm[1:])])))))
    Na = Nn.compose(imgs)
    ex = job.integrate(Na, curs_n, types, (0, 1, 2))
    exact_f = {}
    for m, d in ex.items():
        exact_f[m] = sum(complex(sp.N(P.as_expr())) * (math.log(int(b[3:])) if b != '1' else 1.0) for b, P in d.items())
    params = []
    for j, t in enumerate(types):
        if t == 'S':
            params.append((j, True))
        elif t == 'K':
            params += [(j, False), (j, False)]
        else:
            params.append((j, False))
    sidx = [i for i, p in enumerate(params) if p[1]][0]
    params = [params[sidx]] + [p for i, p in enumerate(params) if i != sidx]
    nt = len(params)
    assert nt == 3
    q = job.q
    js = types.index('S')
    ks = [q[m] - curs_n[js][m] for m in range(4)]
    Nt = Na * (-sum((x * x for x in ks), R.zero)) * R(QQ_I(sp.Rational(1, 4), 0))
    laps = [Nt]
    while laps[-1]:
        D = laps[-1]
        laps.append(sum((D.diff(q[m]).diff(q[m]) for m in range(4)), R.zero))
    laps = laps[:-1]

    def tonum(P):
        exps = []; cfs = []
        for mon, cf in P.items():
            cfs.append(complex(float(cf.x), float(cf.y)))
            exps.append([mon[job.q_idx[m]] for m in range(4)])
        return np.array(exps, dtype=int).reshape(-1, 4), np.array(cfs, dtype=complex)
    lapn = [tonum(D) for D in laps]
    cursn = np.array([[complex(float(sp.re(cc.as_expr())), float(sp.im(cc.as_expr()))) if cc else 0j for cc in cu] for cu in curs_n])

    def upper_gamma(z, x):
        if z > 0:
            return SPF.gammaincc(z, x) * math.gamma(z)
        k = -z
        return x ** (-k) * SPF.expn(k + 1, x)

    def f(beta, gam, eps):
        """arrays beta, gam (same shape) -> integrand array."""
        al = [beta, gam, 1 - beta - gam]
        cur = [cursn[params[p][0]] for p in range(nt)]
        cen = [sum(al[p] * cur[p][m] for p in range(nt)) for m in range(4)]
        Phi = sum(al[p] * np.sum(cur[p] ** 2).real for p in range(nt)) - sum((cen[m] ** 2).real for m in range(4))
        Phi = np.real(Phi)
        tlo = eps * np.maximum(1 / al[1], 1 / al[2]); thi = eps / al[0]
        tot = np.zeros(beta.shape, dtype=complex)
        for j, (E, Cf) in enumerate(lapn):
            nj = np.zeros(beta.shape, dtype=complex)
            for e_, c_ in zip(E, Cf):
                v = c_ * np.ones(beta.shape, dtype=complex)
                for m in range(4):
                    if e_[m]:
                        v = v * cen[m] ** e_[m]
                nj += v
            nj /= (math.factorial(j) * 4 ** j)
            z = nt - 3 - j + 1
            G = (upper_gamma(z, Phi * tlo) - upper_gamma(z, Phi * thi)) * Phi ** (-z)
            tot += nj * G
        return tot
    return dict(exact=exact_f, f=f)


def gl_integral(F, eps, npts=40):
    """int_0^{1/3} dbeta int_beta^{1-2beta} dgam F, beta split geometrically at eps*10^k (the scale of the eps-dependence),
    gam split at (1-beta)/2; tensor Gauss-Legendre on every cell."""
    import numpy as np
    x, w = np.polynomial.legendre.leggauss(npts)
    bps = [0.0]
    b = eps / 10
    while b < 1 / 3:
        bps.append(b); b *= 4
    bps.append(1 / 3)
    tot = 0j
    for b0, b1 in zip(bps[:-1], bps[1:]):
        xb = b0 + (b1 - b0) * (x + 1) / 2; wb = w * (b1 - b0) / 2
        for lo_f, hi_f in ((lambda b: b, lambda b: (1 - b) / 2), (lambda b: (1 - b) / 2, lambda b: 1 - 2 * b)):
            B = np.repeat(xb[:, None], npts, axis=1)
            lo = lo_f(B); hi = hi_f(B)
            G = lo + (hi - lo) * (x[None, :] + 1) / 2
            W = wb[:, None] * w[None, :] / 2 * (hi - lo)
            tot += np.sum(W * F(B, G))
    return tot


if __name__ == '__main__':
    rep = run()
    json.dump(rep, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'out.nosync', 'numeps.json'), 'w'), indent=1)
