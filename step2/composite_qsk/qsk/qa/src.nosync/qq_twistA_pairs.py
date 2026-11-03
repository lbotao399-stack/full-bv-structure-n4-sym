"""QA: engine A's own code (read-only) on all beta_I beta_J pairs of A's twisted letters, outputs gamma_k x dc_a,
generic momenta; tests the SU(3) covariance of A's twisted-letter construction."""
import os, sys, json
import numpy as np
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_twist as TW, qa_core as QC, qa_util as U
T = TW.Twist()
gx = QC.odd_gen(T.xi)
rng = np.random.default_rng(9)
print('taus', T.taus, 'z', T.z, 'u', T.u)
res = {}
for (i, j) in ((0, 1), (0, 2), (1, 2)):
    for k in range(3):
        q = [sum(rng.integers(-3, 4) * T.nv[w] for w in T.hol + T.exa) for _ in range(2)]
        legs = [('AW', T.gam_out[k], q[0]), ('lamM', T.dc_out[0] + 2 * T.dc_out[1], q[1])]
        v = QC.pair_value(('lamP', T.beta_z[i], None), ('lamP', T.beta_z[j], None), legs, gx, what=('fin',))
        res['beta%d beta%d -> gamma%d' % (i, j, k)] = U.vstr(v['fin'])
        print(i, j, k, U.vstr(v['fin']), flush=True)
# also: the same pair with the 'exact' S_+ component (u' instead of u) for beta, to see which component carries B's value
up_letters = [T.dual[T.up + t] for t in T.taus]
for (i, j) in ((0, 1),):
    for k in range(3):
        q = [sum(rng.integers(-3, 4) * T.nv[w] for w in T.hol + T.exa) for _ in range(2)]
        legs = [('AW', T.gam_out[k], q[0]), ('lamM', T.dc_out[0] + 2 * T.dc_out[1], q[1])]
        v = QC.pair_value(('lamP', up_letters[i], None), ('lamP', up_letters[j], None), legs, gx, what=('fin',))
        res["u' letters %d %d -> gamma%d" % (i, j, k)] = U.vstr(v['fin'])
        print('uprime', i, j, k, U.vstr(v['fin']), flush=True)
json.dump(res, open(os.path.join(QSK, 'qa/out.nosync/twistA_pairs.json'), 'w'), indent=1)
