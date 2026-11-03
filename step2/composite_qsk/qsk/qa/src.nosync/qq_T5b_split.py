"""QA: engine A (read-only) placement split of the (A_V, A_perp) -> (c, lam*) rows at n = 0 and n = 1: which
configurations (S2 quadratic block / S3 = T5 with its output) carry the value.  Engine B places T5 only with an A
output (qb_job.DJob.sequences), so the T5 configurations with the lam* leg as output are absent in route B."""
import os, sys, json
import numpy as np
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_core as QC, qa_util as U
resA = json.load(open(os.path.join(QSK, 'engine_a/results.json')))
xi = np.array([complex(x) for x in resA['tasks']['TQ2']['xi']])
gx = QC.odd_gen(xi)
out = {}
for n, fn in ((0, 'tq2_n0.json'), (1, 'tq2_n1.json')):
    rows = json.load(open(os.path.join(QSK, 'engine_a/out.nosync', fn)))
    for r in rows:
        if r['pair'] != ['AV', 'AW'] or r['leg_set'] != 'c lamsM':
            continue
        d = r['data']
        za = np.array([complex(x) for x in d['zeta_a']]); zb = np.array([complex(x) for x in d['zeta_b']])
        fa = {(1, 0, 0, 0): 1j} if n == 1 else None
        legs = [(t, np.array([complex(x) for x in p]), np.array(m)) for t, p, m in d['legs']]
        tot, parts = QC.pair_value(('AV', za, fa), ('AW', zb, None), legs, gx, what=('fin',), split=True)
        rec = {'total': U.vstr(tot['fin']), 'A_table': r['eps0'], 'parts': {str(k): U.vstr(v) for k, v in parts.items()}}
        print('n=%d' % n, r['ordering'], rec, flush=True)
        out['n%d %s' % (n, r['ordering'])] = rec
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/T5b_split.json'), 'w'), indent=1)
