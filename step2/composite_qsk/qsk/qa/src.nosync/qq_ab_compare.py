"""QA (Part Q): entry-by-entry comparison of engine A's TQ2 rows with engine B's code, run here on A's data.
Engine A row data (xi, letters, output polarisations, momenta; A basis) are mapped to engine B's basis with the exact
spinor intertwiner S of revise.nosync/ferm_map/spinor_map.json (u_B = S u_A on Sigma, same vector index M); under S
A's lamP is B's lamM and vice versa.  Letters are covectors (S^{-T}), output polarisations and xi are vectors (S).
S = c S0 with c^2 = C2 = i/16; B is run with S0 and multiplied by C2^{e/2}, e = (n_out_lam - n_let_lam)
- (n_out_lams - n_let_lams) + 1 (xi).  Engine B's code is imported read-only (qsk/engine_b/src.nosync).
Usage: python3 -B qq_ab_compare.py <selector> ; writes qa/out.nosync/ab_compare_<selector>.json incrementally."""
import os, sys, json, time, resource, itertools
import sympy as sp

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
OUT = os.path.join(QSK, 'qa/out.nosync')
sys.path.insert(0, os.path.join(QSK, 'engine_b/src.nosync'))
import qb_job as QJ          # noqa: E402  (engine B, read-only)
from realization import pairing as _pairing   # noqa: E402

SM = json.load(open(os.path.join(ROOT, 'step2/composite_qsk/revise.nosync/ferm_map/spinor_map.json')))
C2 = sp.sympify(SM['c2'])
S = sp.Matrix(16, 16, lambda i, j: sp.nsimplify(sp.simplify(sp.sympify(SM['S'][i][j]) / sp.sqrt(C2))))
SinvT = S.inv().T
Mp = _pairing()
SST0 = (S.T * Mp).inv()
SWAP = {'lamP': 'lamM', 'lamM': 'lamP', 'lamsP': 'lamsM', 'lamsM': 'lamsP'}
COMPS_B = {'AV': list(range(4)), 'AW': list(range(4, 10)), 'AsV': list(range(4)), 'AsW': list(range(4, 10)), 'c': [0], 'cs': [0]}
LOG = {'1': 1, 'log2': sp.log(2), 'log3': sp.log(3), 'log5': sp.log(5), 'log7': sp.log(7)}


def cplx(s):
    z = complex(s)
    return sp.nsimplify(z.real) + sp.I * sp.nsimplify(z.imag)


def to_B(t, z, kind):
    if t in ('lamP', 'lamM', 'lamsP', 'lamsM'):
        zs = sp.Matrix([cplx(x) for x in z])
        if t in ('lamP', 'lamM'):
            zb = (SinvT * zs) if kind == 'letter' else (S * zs)
        else:
            zb = (SST0.inv().T * zs) if kind == 'letter' else (SST0 * zs)
        tb = SWAP[t]
        zb = sp.expand(zb)
        rng_ = range(8, 16) if tb in ('lamM', 'lamsP') else range(0, 8)
        assert all(sp.simplify(zb[i]) == 0 for i in range(16) if i not in rng_), (t, kind)
        return tb, {i: sp.nsimplify(zb[i]) for i in rng_}
    return t, {c: cplx(x) for c, x in zip(COMPS_B[t], z)}


def Aval(s):
    return sp.expand(sp.sympify(s.replace('i', 'I').replace('+-', '-').replace('log2', 'log(2)').replace('log3', 'log(3)')
                                .replace('log5', 'log(5)').replace('log7', 'log(7)')))


def main(sel):
    rowsA = json.load(open(os.path.join(QSK, 'engine_a/out.nosync/tq2_n0.json')))
    resA = json.load(open(os.path.join(QSK, 'engine_a/results.json')))
    xiA = [cplx(x) for x in resA['tasks']['TQ2']['xi']]
    xiB = sp.expand(S * sp.Matrix(xiA))
    assert all(sp.simplify(xiB[i]) == 0 for i in range(8)), 'xi_A not in B components 8..15'
    gvec = [sp.nsimplify(xiB[i]) for i in range(16)]
    pick = []
    for r in rowsA:
        key = '%s|%s:%s' % (r['pair'][0], r['pair'][1], r['leg_set'])
        if sel == 'all' or key in sel.split(','):
            if len(r['ordering']) <= 2:
                pick.append(r)
    out = []
    t0 = time.time()
    fn = os.path.join(OUT, 'ab_compare_%s.json' % ('all' if sel == 'all' else str(abs(hash(sel)) % 10 ** 6)))
    for r in pick:
        ta, tb = r['pair']
        d = r['data']
        tBa, zBa = to_B(ta, d['zeta_a'], 'letter'); tBb, zBb = to_B(tb, d['zeta_b'], 'letter')
        outsB, polB, moms = [], [], []
        for (t, p, m) in d['legs']:
            tt, pp = to_B(t, p, 'out'); outsB.append(tt); polB.append(pp); moms.append(m)
        lamT = ('lamP', 'lamM'); lamsT = ('lamsP', 'lamsM')
        types = [l[0] for l in d['legs']]
        e = (sum(t in lamT for t in types) - sum(t in lamT for t in (ta, tb))
             - sum(t in lamsT for t in types) + sum(t in lamsT for t in (ta, tb)) + 1)
        assert e % 2 == 0, e
        fac = C2 ** (e // 2)
        job = QJ.DJob(tBa, tBb, outsB, 'xi', gvec, pol=polB, za=zBa, zb=zBb, tag='qaab')
        res = job.kernel(mlist=(0,))
        subsd = {}
        for j, m in enumerate(moms):
            for mu in range(4):
                subsd[sp.Symbol(f'r{j}_{mu}')] = int(m[mu])
        bw = {}
        for wd, dm in res['total'].items():
            tot = 0
            for bb, P in dm.get(0, {}).items():
                tot += sp.expand(P.as_expr(*[sp.Symbol(n) for n in job.names]).subs(subsd)) * LOG.get(bb, sp.Symbol('L_' + str(bb)))
            bw[str(list(wd))] = sp.expand(fac * tot)
        a = Aval(r['eps0'])
        match = [w for w, v in bw.items() if sp.expand(v - a) == 0]
        match_neg = [w for w, v in bw.items() if sp.expand(v + a) == 0]
        zeroB = all(v == 0 for v in bw.values())
        rec = {'pair_A': r['pair'], 'leg_set_A': r['leg_set'], 'ordering_A': r['ordering'], 'types_B': [tBa, tBb, outsB],
               'A_eps0': str(a), 'B_words_eps0': {k: str(v) for k, v in bw.items()}, 'match_words': match,
               'match_words_with_minus': match_neg, 'A_zero': a == 0, 'B_all_zero': zeroB,
               'agree': bool(match) or (a == 0 and zeroB),
               'seconds': round(time.time() - t0, 1),
               'rss_MB': round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, 1)}
        print(r['pair'], r['leg_set'], r['ordering'], 'agree:', rec['agree'], 'minus:', bool(match_neg), rec['seconds'], 's', rec['rss_MB'], 'MB', flush=True)
        out.append(rec)
        json.dump(out, open(fn, 'w'), indent=1)
    print('ROWS', len(out), 'AGREE', sum(r['agree'] for r in out), 'AGREE_WITH_MINUS', sum(bool(r['match_words_with_minus']) for r in out))


if __name__ == '__main__':
    QJ._watchdog(1300)
    main(sys.argv[1] if len(sys.argv) > 1 else 'lamP|lamP:AW lamM')
