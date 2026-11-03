"""QA: dump engine A's twisted-letter data (read-only import of engine A's qa_twist.Twist) to JSON, and evaluate
engine A's own kernel on beta-beta and beta-gamma at generic (not only twisted) output momenta for a cross-check."""
import os, sys, json
import numpy as np

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_twist as TW       # noqa: E402
import qa_core as QC        # noqa: E402
import qa_util as U         # noqa: E402

T = TW.Twist()
c = lambda v: [str(complex(x)) for x in np.asarray(v).ravel()]
d = {'xi': c(T.xi), 'beta_z': [c(v) for v in T.beta_z], 'beta_out': [c(v) for v in T.beta_out],
     'gam_w': [c(v) for v in T.gam_w], 'gam_out': [c(v) for v in T.gam_out], 'dc_z': [c(v) for v in T.dc_z],
     'dc_out': [c(v) for v in T.dc_out], 'nv': {str(k): c(v) for k, v in T.nv.items()}, 'hol': T.hol, 'exa': T.exa,
     'taus': T.taus}
# A's own values: beta_0 beta_1 -> (gamma_k, dc_a) at a few momenta (exact-span and generic)
gx = QC.odd_gen(T.xi)
vals = []
rng = np.random.default_rng(5)
for k in range(3):
    for a in range(2):
        for trial in range(2):
            q = []
            for j in range(2):
                co = rng.integers(-3, 4, 4)
                qq = sum(co[i] * T.nv[w] for i, w in enumerate(T.hol + T.exa))
                q.append(qq)
            legs = [('AW', T.gam_out[k], q[0]), ('lamM', T.dc_out[a], q[1])]
            v = QC.pair_value(('lamP', T.beta_z[0], None), ('lamP', T.beta_z[1], None), legs, gx, what=('fin',))
            vals.append({'k': k, 'a': a, 'q': [c(x) for x in q], 'A_value': U.vstr(v['fin'])})
            print(k, a, trial, U.vstr(v['fin']), flush=True)
d['A_betabeta_values'] = vals
json.dump(d, open(os.path.join(QSK, 'qa/out.nosync/twistA_data.json'), 'w'), indent=1)
