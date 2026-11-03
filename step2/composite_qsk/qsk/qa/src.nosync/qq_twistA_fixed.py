"""QA: engine A's twist entries with beta letters rescaled to Gaussian integers.
Finding: ea_kern.contract (two-letter engine A, inherited) rounds the products of letter components to Gaussian
integers (int(round(.))) without a check; engine A's TQ3 beta letters are dual rows with entries +-1/16, +-i/16,
so every entry with a beta letter evaluated to 0.  Here the letters are 16 x beta (exact Gaussian integers); the
values below are 16 x (per beta letter) the true values.  Also: b-gamma with the gamma output doubled (linearity
check of the output polarisations, which have entries 1/2)."""
import os, sys, json, time
import numpy as np
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_twist as TW, qa_core as QC
T = TW.Twist()
gx = QC.odd_gen(T.xi)
rng = np.random.default_rng(3)
log = []
def L(msg):
    print(msg, flush=True); log.append(msg)
def letters(which, idx=None, scale=16):
    if which == 'beta':
        return [(1, scale * T.beta_z[idx], None, 'lamP')]
    return TW.letters(T, which, idx)
res = {'note': 'beta letters multiplied by 16 (each); divide a value by 16^k for k beta letters', 'rows': {}}
spec = {
    'betabeta(16x16)': (('beta', 0), ('beta', 1), ('gamma', 'dc'), 1),
    'betagamma(16x)': (('beta', 0), ('gamma', 0), ('dc', 'dc'), 1),
    'betab_bdc(16x)': (('beta', 0), ('b', None), ('beta', 'dc'), 1),
    'betab_gg(16x)': (('beta', 0), ('b', None), ('gamma', 'gamma'), 2),
}
t0 = time.time()
for nm, ((wa, ia), (wb, ib), (o1, o2), D) in spec.items():
    La, Lb = letters(wa, ia), letters(wb, ib)
    rows = TW.entry(T, gx, La, Lb, (TW.OUT_TYPE[o1], TW.OUT_TYPE[o2]), [TW.outs(T, o1), TW.outs(T, o2)], D, rng, log=L, label=nm)
    res['rows'][nm] = rows
    L('%s done %.0f s' % (nm, time.time() - t0))
    json.dump(res, open(os.path.join(QSK, 'qa/out.nosync/twistA_fixed.json'), 'w'), indent=1)
# linearity check of b-gamma in the gamma output (entries 1/2): doubled output vector
La, Lb = TW.letters(T, 'b', None), TW.letters(T, 'gamma', 0)
outs_g2 = [('g%d' % a, 2 * v) for a, v in enumerate(T.gam_out)]
rows = TW.entry(T, gx, La, Lb, ('AW', 'lamM'), [outs_g2[:1], TW.outs(T, 'dc')[:1]], 1, rng, log=L, label='bgamma_gout_x2')
res['rows']['bgamma_gout_x2'] = rows
json.dump(res, open(os.path.join(QSK, 'qa/out.nosync/twistA_fixed.json'), 'w'), indent=1)
L('all done %.0f s' % (time.time() - t0))
