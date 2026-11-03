"""QA: engine B's code on engine A's twisted data, all beta pairs and all gamma outputs (complement of the first run)."""
import os, sys, json
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'qa/src.nosync'))
import qq_ab_compare as AB
QJ = AB.QJ
d = json.load(open(os.path.join(QSK, 'qa/out.nosync/twistA_data.json')))
xiB = sp.expand(AB.S * sp.Matrix([AB.cplx(x) for x in d['xi']]))
gvec = [sp.nsimplify(xiB[i]) for i in range(16)]
out = {}
for (i, j) in ((0, 1), (0, 2), (1, 2)):
    tBa, za = AB.to_B('lamP', d['beta_z'][i], 'letter'); tBb, zb = AB.to_B('lamP', d['beta_z'][j], 'letter')
    for k in range(3):
        tt0, p0 = AB.to_B('AW', d['gam_out'][k], 'out')
        tt1, p1 = AB.to_B('lamM', d['dc_out'][0], 'out')
        job = QJ.DJob(tBa, tBb, [tt0, tt1], 'xi', gvec, pol=[p0, p1], za=za, zb=zb, tag='qatw2')
        res = job.kernel(mlist=(0,))
        nz = {str(list(w)): {b: str(P.as_expr(*[sp.Symbol(n) for n in job.names])) for b, P in dm.get(0, {}).items()} for w, dm in res['total'].items()}
        print('beta', i, j, 'gamma', k, 'B words:', 'ZERO' if not nz else nz, flush=True)
        out['beta%d beta%d gamma%d' % (i, j, k)] = nz or 'zero'
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/twist_B_on_A2.json'), 'w'), indent=1)
