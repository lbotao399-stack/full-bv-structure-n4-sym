"""QA: emulate route B's omission inside engine A (read-only import): T5 allowed only with an A output
(cubtypes = ('AV','AW')) versus the full set; on the lam*-output rows."""
import os, sys, json
import numpy as np
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_core as QC, qa_util as U
resA = json.load(open(os.path.join(QSK, 'engine_a/results.json')))
xi = np.array([complex(x) for x in resA['tasks']['TQ2']['xi']])
gfull = QC.odd_gen(xi)
gA = QC.odd_gen(xi, name='xi_Aonly'); gA.cubtypes = ('AV', 'AW')
out = []
targets = [(0, ['lamP', 'lamP'], 'AV c lamsM'), (0, ['lamP', 'lamM'], 'AV c lamsP'), (0, ['lamP', 'lamM'], 'AW c lamsM'),
           (0, ['lamM', 'lamM'], 'AW c lamsP'), (1, ['AV', 'AW'], 'c lamsM')]
for n, pair, ls in targets:
    rows = json.load(open(os.path.join(QSK, 'engine_a/out.nosync', 'tq2_n%d.json' % n)))
    for r in rows:
        if r['pair'] != pair or r['leg_set'] != ls:
            continue
        d = r['data']
        za = np.array([complex(x) for x in d['zeta_a']]); zb = np.array([complex(x) for x in d['zeta_b']])
        fa = {(1, 0, 0, 0): 1j} if n == 1 else None
        legs = [(t, np.array([complex(x) for x in p]), np.array(m)) for t, p, m in d['legs']]
        vf = QC.pair_value((pair[0], za, fa), (pair[1], zb, None), legs, gfull)['fin']
        va = QC.pair_value((pair[0], za, fa), (pair[1], zb, None), legs, gA)['fin']
        rec = {'n': n, 'pair': pair, 'leg_set': ls, 'ordering': r['ordering'], 'full': U.vstr(vf), 'T5_A_output_only': U.vstr(va)}
        print(rec, flush=True)
        out.append(rec)
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/T5b_emulate.json'), 'w'), indent=1)
