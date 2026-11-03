"""G-KT analogue (route B): independence of the eps-expansion coefficients of a defect cycle on the transfer endpoint T
(Thm qs-simplex (iii), (iv)).  Ordinary Schwinger times in [eps, T], the special one in [0, eps]; tau-range
[eps/min_f alpha_f, min(eps/alpha_sigma, T/max_f alpha_f)]; fits as qb_numeps at T|p|^2 in {1, 4, 16} (|p| ~ 1 for the
data used) and T = infinity; the eps^0 and eps^-1 coefficients must agree with the exact T = infinity values."""
import os, json, time, math
import mpmath as mp
import sympy as sp
import qb_core as B
import qb_job as QJ
import qb_numeps as NE

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
mp.mp.dps = 30


def value_T(job, types, curs, Num, eps, T):
    R = job.R; q = job.q
    nt = len(types)
    js = types.index('S')
    order = [js] + [j for j in range(nt) if j != js]
    ks = [q[m] - curs[js][m] for m in range(4)]
    Nt = Num * (-sum((x * x for x in ks), R.zero))
    laps = [Nt]
    while laps[-1]:
        D = laps[-1]
        laps.append(sum((D.diff(q[m]).diff(q[m]) for m in range(4)), R.zero))
    laps = laps[:-1]
    qi = job.q_idx

    def tonum(P):
        return [(tuple(mon[qi[m]] for m in range(4)), mp.mpc(NE.q2mp(cf.x), NE.q2mp(cf.y))) for mon, cf in P.items()]
    lapn = [tonum(D) for D in laps]

    def cnum(c):
        if not c:
            return mp.mpf(0)
        v = c.coeff(1)
        return mp.mpc(NE.q2mp(v.x), NE.q2mp(v.y))
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

    def integrand(alpha):
        if alpha[0] <= 0:
            return mp.mpc(0)
        cent = [sum(alpha[p] * cursn[order[p]][m] for p in range(nt)) for m in range(4)]
        Phi = sum(alpha[p] * sum(cursn[order[p]][m] ** 2 for m in range(4)) for p in range(nt)) - sum(c * c for c in cent)
        Phi = mp.re(Phi)
        tlo = eps / min(alpha[1:])
        thi = eps / alpha[0] if alpha[0] > 0 else mp.inf
        if T is not None:
            thi = min(thi, T / max(alpha[1:]))
        if thi <= tlo:
            return mp.mpc(0)
        tot = mp.mpc(0)
        for j, terms in enumerate(lapn):
            nj = evalp(terms, cent) / (math.factorial(j) * 4 ** j)
            if nj == 0:
                continue
            kp = nt - 3 - j
            ti = NE.tau_int(kp, Phi, tlo, thi)
            tot += nj * ti
        return tot
    if nt == 2:
        pts = [mp.mpf(0)]
        if T is not None:
            astar = eps / (T + eps)          # eps/alpha_sigma = T/(1 - alpha_sigma): kink of the tau range
            if 0 < astar < mp.mpf(1) / 2:
                pts.append(astar)
        pts.append(mp.mpf(1) / 2)
        return mp.quad(lambda a: integrand([a, 1 - a]), pts)
    raise NotImplementedError


def run():
    t0 = time.time()
    xi = B.XI[:8] + [0] * 8
    job = QJ.DJob('AW', 'lamP', ['AW'], 'xi', xi, terms='cubic', tag='numT')
    vals = NE.rnd_vals(job, 7)
    out = {}; ok = True
    for seq in job.sequences():
        for sigma in range(len(seq) + 1):
            types, curs, Num = NE.config_data(job, seq, sigma, vals)
            if not Num:
                continue
            ex = NE.exact_coeffs(job, types, curs, Num)
            row = {}
            for T in (None, 1, 4, 16):
                epss = [mp.mpf(2) ** -10 / mp.mpf(2) ** j for j in range(9)]
                vs = [value_T(job, types, curs, Num, e, (mp.mpf(T) if T else None)) for e in epss]
                fr = NE.fit(epss, vs)
                d0 = abs(fr['a_0'] - ex.get(0, 0)); d1 = abs(fr['a_m1'] - ex.get(1, 0))
                dlin = [mp.nstr((value_T(job, types, curs, Num, e, (mp.mpf(T) if T else None)) - NE.numeric_value(job, types, curs, Num, e)) / e, 10)
                        for e in (mp.mpf(2) ** -12, mp.mpf(2) ** -16, mp.mpf(2) ** -20)] if T else None
                row[str(T)] = dict(T_minus_infinity_over_eps_at_eps_2m12_2m16_2m20=dlin,fit_eps0=mp.nstr(fr['a_0'], 20), fit_eps_m1=mp.nstr(fr['a_m1'], 20), diff_eps0=mp.nstr(d0, 4), diff_eps_m1=mp.nstr(d1, 4),
                                   log_eps=mp.nstr(fr['c_log'], 4), value_at_eps_2m13=mp.nstr(vs[-1], 15))
                ok = ok and d0 < mp.mpf(10) ** -9 and d1 < mp.mpf(10) ** -12
            out[f'{seq} sigma={sigma}'] = dict(exact_eps0=mp.nstr(ex.get(0, 0), 20), exact_eps_m1=mp.nstr(ex.get(1, 0), 20), by_T=row)
            print(seq, sigma, json.dumps(row)[:300], flush=True)
    out['data'] = {k: str(v) for k, v in vals.items()}
    out['seconds'] = round(time.time() - t0, 1)
    p2 = sum(float(vals[f'r0_{m}']) ** 2 for m in range(4))
    out['|p|^2'] = p2
    gates = [dict(id='G-KT (endpoint independence, route B numerics)', pass_=bool(ok), value={k: {T: (r['diff_eps0'], r['diff_eps_m1']) for T, r in v['by_T'].items()} for k, v in out.items() if isinstance(v, dict) and 'by_T' in v},
                  target='eps^0 and eps^-1 coefficients of the T5 cycle at T = 1, 4, 16 (|p|^2 = %.3f) and T = infinity equal the exact T = infinity values (Thm qs-simplex (iii),(iv)); the finite-eps values differ' % p2)]
    return out, gates


if __name__ == '__main__':
    out, gates = run()
    json.dump(dict(results=out, gates=gates), open(os.path.join(OUT, 'numT.json'), 'w'), indent=1)
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL')
