"""Revision (Part Q), gate G-KT along K_b (finding Q14), route B numerics.

Thm qs-simplex (iv): T d/dT of every defect cycle along K_b is O(eps) at fixed finite T, so the eps^-1, eps^0 (and log eps)
coefficients, hence Upsilon_{K_b}, do not depend on T.  Tested here on the defect cycles of the intrinsic special conformal
vertex S^(1)(K_{b,x_O}) (qb_jobinj.register_K_terms: V_{mu nu} multiplying z^mu z^nu, W_mu multiplying z^mu; z = x - x_O)
on the pair (A_V, A_V) with output A_V, every configuration with n_t = 3 edges: S^(1)(K_b) on a line (the A_V--c composite
line on the cycle, special edge at the far end) and at the special edge.  With ordinary Schwinger times in [eps, T] and the
special one in [0, eps], the value is V(eps,T) = int dalpha sum_j n_j(alpha) int_{eps/min alpha_f}^{min(eps/alpha_sigma, T/M)}
tau^{-j} e^{-tau Phi} dtau, M = max_{f != sigma} alpha_f, so that
    T dV/dT = int_{alpha_sigma < eps M/T} dalpha sum_j n_j(alpha) (T/M)^{1-j} e^{-T Phi/M},
an elementary integral over a region of measure O(eps/T) (Gauss-Legendre, double precision, exact region).  The intrinsic
part is the h-derivative at h = 0 of the cycle with momentum h injected at the K_b vertex (z^mu z^nu -> -d_mu d_nu,
z^mu -> -i d_mu): taken by central differences (step 1e-3) of T dV/dT, which is smooth in h at finite T.  Pass: the
derivative combination of T dV/dT vanishes at least like eps (ratio >= ~4 per factor 4 in eps) at T = 1, 4, 16 (|q|^2 ~ 1.4).
Engine B code imported read-only; data rational."""
import os, sys, json, time, math, resource, itertools
from fractions import Fraction as Fr
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_b/src.nosync'))
import numpy as np
import sympy as sp
import qb_core as B
import qb_job as QJ
import qb_jobinj as IJ
import qb_numeps as NE

OUT = os.path.join(QSK, 'revision.nosync/out')
XG, WG = np.polynomial.legendre.leggauss(48)


def gl(lo, hi):
    x = 0.5 * (hi - lo) * XG + 0.5 * (hi + lo)
    return x, 0.5 * (hi - lo) * WG


def prep(job, types, curs, Num):
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
    lapn = []
    for j, P in enumerate(laps):
        terms = [(tuple(mon[qi[m]] for m in range(4)), complex(float(cf.x), float(cf.y)) / (math.factorial(j) * 4 ** j)) for mon, cf in P.items()]
        lapn.append(terms)

    def cnum(c):
        if not c:
            return 0.0
        v = c.coeff(1)
        return complex(float(v.x), float(v.y))
    cursn = np.array([[cnum(c) for c in cu] for cu in curs])           # nt x 4
    return nt, order, lapn, cursn


def tdT(data, eps, T):
    nt, order, lapn, cursn = data
    assert nt == 3
    C = cursn[order]                                                   # rows in 'order' (special first)
    sq = np.sum(C * C, axis=1)

    def f(a, b):
        al = np.stack([a, b, 1 - a - b])                               # 3 x N
        cent = al.T @ C                                                # N x 4
        Phi = (al.T @ sq) - np.sum(cent * cent, axis=1)
        M = np.maximum(b, 1 - a - b)
        tot = np.zeros(a.shape, complex)
        for j, terms in enumerate(lapn):
            nj = np.zeros(a.shape, complex)
            for qe, cv in terms:
                v = np.full(a.shape, cv, complex)
                for m in range(4):
                    if qe[m]:
                        v = v * cent[:, m] ** qe[m]
                nj = nj + v
            tot = tot + nj * (T / M) ** (1 - j)
        live = np.minimum(b, 1 - a - b) > eps * M / T                    # T/M above the lower end eps/min alpha_f
        return tot * np.exp(-T * Phi.real / M) * live

    tot = 0j
    d = 8 * eps / T
    bcuts = [0.0, d, 0.5 - d, 0.5, 0.5 + d, 1.0]
    for b0, b1 in zip(bcuts[:-1], bcuts[1:]):
        bs, wb = gl(b0, b1)
        for b, w in zip(bs, wb):
            ub = min(b, (1 - b) / 2)
            ivs = []
            if 1 - 2 * b > 0:
                ivs.append((0.0, min(1 - 2 * b, eps * (1 - b) / (T + eps), ub)))
            ivs.append((max(0.0, 1 - 2 * b), min(eps * b / T, ub)))
            for lo, hi in ivs:
                if hi <= lo:
                    continue
                aa, wa = gl(lo, hi)
                tot += w * np.sum(wa * f(aa, np.full(aa.shape, b)))
    return tot


def main():
    t0 = time.time()
    b = list(sp.symbols('b0:4'))
    names = IJ.register_K_terms(b, 'gkt')
    picks = [(('V', 0, 1), [((1, 1), 1), ((1, -1), -1), ((-1, 1), -1), ((-1, -1), 1)], (0, 1), 'd^2/dh0 dh1'),
             (('W', 2), [((1,), 1), ((-1,), -1)], (2,), 'd/dh2')]
    dl = 1e-3
    base = None
    rows = []
    ok = True
    for key, stencil, dirs, dname in picks:
        hmax = 2 if key[0] == 'V' else 1
        job = IJ.IJob('AV', 'AV', ['AV'], 'terms', b, hmax=hmax, quad_terms=names[key], Spar=1, tag='gkt%s' % ''.join(map(str, key)))
        if base is None:
            base = NE.rnd_vals(job, 13)
            base.update({'r0_0': Fr(1), 'r0_1': Fr(1, 2), 'r0_2': Fr(-1, 3), 'r0_3': Fr(1, 4),
                         'b0': Fr(1), 'b1': Fr(2), 'b2': Fr(-1), 'b3': Fr(1)})
        for seq in job.sequences():
            nt = len(seq) + 1
            if nt != 3:
                continue
            spos = [i for i, (kd, bl) in enumerate(seq) if kd.startswith('S')][0]
            for sigma in range(nt):
                types0, curs0, words = job.config(seq, sigma)
                Num0 = job.R.zero
                for w, P in words.items():
                    Num0 += P
                if not Num0:
                    continue
                placement = 'line (A_V--c line on the cycle)' if sigma not in (spos, spos + 1) else 'special edge'
                datas = []
                for sh, coef in stencil:
                    v = {k: x for k, x in base.items() if k in job.names}
                    for m in range(4):
                        v['h%d' % m] = Fr(0)
                    for dmu, s in zip(dirs, sh):
                        v['h%d' % dmu] = Fr(s) * Fr(1, 1000)
                    imgs = NE.numeric_images(job, v)
                    Num = Num0.compose(imgs)
                    curs = [[c.compose(imgs) if c else job.R.zero for c in cu] for cu in curs0]
                    datas.append((coef, prep(job, types0, curs, Num)))
                byT = {}
                for T in (1.0, 4.0, 16.0):
                    ds = []
                    for k in (8, 10, 12):
                        e = 2.0 ** -k
                        val = sum(coef * tdT(dat, e, T) for coef, dat in datas) / (2 * dl) ** len(dirs)
                        ds.append(val)
                    ratios = [abs(ds[i] / ds[i + 1]) if abs(ds[i + 1]) > 0 else None for i in range(2)]
                    lin = all(r is not None and r > 3.0 for r in ratios)     # O(eps) or faster (ratio 4: eps, 16: eps^2)
                    tiny = all(abs(x) < 1e-11 for x in ds)
                    okT = bool(lin or tiny)
                    ok = ok and okT
                    byT[str(T)] = dict(T_dV_dT=[complex(x) for x in ds], over_eps=[complex(x / 2.0 ** -k) for x, k in zip(ds, (8, 10, 12))],
                                       ratios_eps_x4=ratios, O_eps=okT)
                rec = dict(vertex=str(key), derivative=dname, seq=str(seq), sigma=sigma, placement=placement, types=types0, by_T=byT,
                           seconds=round(time.time() - t0, 1), rss_MB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, 1))
                rows.append(rec)
                print(key, seq, sigma, placement, {T: (v['ratios_eps_x4'], v['O_eps']) for T, v in byT.items()}, rec['seconds'], 's', flush=True)
                json.dump(dict(rows=rows), open(os.path.join(OUT, 'gkt_Kb.json'), 'w'), indent=1, default=str)
    gate = dict(id='G-KT along K_b (route B numerics): T d/dT of the defect cycles of the intrinsic K_b vertex, at the order of the h-derivative, is O(eps)',
                pass_=bool(ok and rows), n_configs=len(rows), placements=sorted(set(r['placement'] for r in rows)),
                target='T dV/dT -> 0 at least linearly in eps (ratio ~4 or ~16 per factor 4) at T = 1, 4, 16 for every configuration (Thm qs-simplex (iv)): '
                       'the eps^-1, eps^0, log eps coefficients of Upsilon^0_{K_b} are T-independent')
    json.dump(dict(rows=rows, gate=gate, data={k: str(x) for k, x in base.items()}, seconds=round(time.time() - t0, 1)),
              open(os.path.join(OUT, 'gkt_Kb.json'), 'w'), indent=1, default=str)
    print(gate)


if __name__ == '__main__':
    main()
