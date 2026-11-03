"""Task T1: Ob^(1)(Delta) on two legs by the strict route (defect cycles of I + c S^(1)(Delta)), intrinsic part (p = 0)
and position part (p = 1, = 2 Ob^(1)(P_b)), poles; split into line and special-edge placements (cross-check X-scal)."""
import json, os, random, time, sys
from fractions import Fraction as Fr
import numpy as np
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_pos as P, eb_fit as F

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')


def rv(n, lo=-3, hi=3):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]


def gens_delta():
    tdir = G.tmpl_dir()
    tw = G.tmpl_weights(2, 3, 0, 'D')

    def gdir(mu, r):
        e = np.zeros(4, complex)
        e[mu] = 1
        return Y.Gen('Delta', 1, tdir, gvec=e, inj=r)
    return gdir, Y.Gen('Delta', 1, tw)


def make_legs(sector):
    p = rv(4)
    if sector == 'AVAV':
        a1 = np.zeros(10, complex); a2 = np.zeros(10, complex)
        a1[:4] = rv(4, -2, 2); a2[:4] = rv(4, -2, 2)
        sp0, sp1, u0, u1 = 'A', 'A', a1, a2
    elif sector == 'ApAp':
        a1 = np.zeros(10, complex); a2 = np.zeros(10, complex)
        a1[4:] = rv(6, -2, 2); a2[4:] = rv(6, -2, 2)
        sp0, sp1, u0, u1 = 'A', 'A', a1, a2
    elif sector == 'lamlam':
        sp0, sp1 = 'lam', 'lam'
        u0 = np.array(rv(16, -2, 2), complex); u1 = np.array(rv(16, -2, 2), complex)
    elif sector == 'cAsV':
        sp0, sp1 = 'c', 'As'
        u0 = np.array([1], complex); u1 = np.zeros(10, complex); u1[:4] = rv(4, -2, 2)
    else:
        raise ValueError(sector)

    def legs_at(r):
        return [dict(sp=sp0, u=u0, p=p), dict(sp=sp1, u=u1, p=[-p[m] - r[m] for m in range(4)])]
    return legs_at


def clog_get(d, key, m):
    return d.get(key, {}).get(m, E.CLog())


def run(sectors=('AVAV', 'ApAp', 'lamlam', 'cAsV'), seed=7, extra=2, log=print):
    random.seed(seed)
    gdir, g0 = gens_delta()
    B = F.basis_delta2()
    res = {}
    for sec in sectors:
        basis = B[sec]
        labels = [b[0] for b in basis]
        nconf = len(basis) + extra
        vals0, vals1, cands, cands_pole = [], [], [], []
        split = {'line': [], 'special': []}
        betas_nonzero, top_nonzero = 0, 0
        t0 = time.time()
        st = {}
        for ic in range(nconf):
            legs_at = make_legs(sec)
            al, be, top = P.position_split(legs_at, gdir, g0, stats=st)
            for mu, d in be.items():
                for k, dm in d.items():
                    for m, v in dm.items():
                        if not v.is_zero():
                            betas_nonzero += 1
            for k, dm in top.items():
                for m, v in dm.items():
                    if not v.is_zero():
                        top_nonzero += 1
            vals0.append(clog_get(al, 'k', 0))
            vals1.append(clog_get(al, 'k', 1))
            legs0 = legs_at([Fr(0)] * 4)
            cands.append([F.value(t, legs0) for (_, t) in basis])
            cands_pole.append([F.value(t, legs0) for (_, t) in B['AA_pole']])
            if ic == 0:
                for filt in ('line', 'special'):
                    Y.PLACEMENT_FILTER = filt
                    a2, _, _ = P.position_split(legs_at, gdir, g0)
                    Y.PLACEMENT_FILTER = None
                    split[filt] = [clog_get(a2, 'k', 0)]
                    split[filt + '_cand'] = [cands[0]]
        fit0, ok0 = F.fit(vals0, cands, labels)
        if sec in ('AVAV', 'ApAp'):
            fit1, ok1 = F.fit(vals1, cands_pole, ['int(A,A)'])
        else:
            fit1, ok1 = F.fit(vals1, [[(Fr(1), Fr(0))] for _ in vals1], ['zero-check'])
            ok1 = ok1 and all(v.is_zero() for v in vals1)
        # split fits (single configuration: ratio onto the first basis element only when one-dimensional)
        sfit = {}
        for filt in ('line', 'special'):
            v = split[filt][0]
            c = split[filt + '_cand'][0]
            sfit[filt] = {k: (str(sp_div(v, c[0], k))) for k in ('1', 'log2', 'log3')}
        res[sec] = dict(eps0=jsonable(fit0), eps0_consistent=ok0, pole=jsonable(fit1), pole_consistent=ok1,
                        position_part_nonzero_entries=betas_nonzero, topcheck_nonzero=top_nonzero,
                        split_ratio_first_basis=sfit, nconf=nconf, chains=st.get('chains', 0),
                        seconds=round(time.time() - t0, 1))
        log(f'T1 {sec}: {res[sec]}')
    return res


def sp_div(v, c, k):
    import sympy as sp
    re = v.re.items().get(k, Fr(0)); im = v.im.items().get(k, Fr(0))
    z = sp.Rational(str(re)) + sp.I * sp.Rational(str(im))
    cz = sp.Rational(str(c[0])) + sp.I * sp.Rational(str(c[1]))
    return sp.nsimplify(sp.simplify(z / cz)) if cz != 0 else 'nan'


def jsonable(fitres):
    if fitres is None:
        return None
    out = {}
    for k, d in fitres.items():
        out[k] = None if d is None else {lab: str(x) for lab, x in d.items()}
    return out


if __name__ == '__main__':
    r = run()
    with open(os.path.join(OUT, 'T1.json'), 'w') as f:
        json.dump(r, f, indent=1)
