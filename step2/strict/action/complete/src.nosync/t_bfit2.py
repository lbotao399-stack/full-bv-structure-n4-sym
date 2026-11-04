import sys, numpy as np
sys.argv = ['x', sys.argv[1], '+']
import xieta as XE, eta1 as H
rng = np.random.default_rng(4)
xi, eta = H.rand_u(XE.XIT, rng), H.rand_u('ls+' if XE.ROUTE == 'B' else 'lamsP', rng)
cD, om, rr = XE.bracket_fit(xi, eta, rng)
cD2, om2, rr2 = XE.bracket_formula(xi, eta)
print('cD', cD, cD2, ' omega antisym', np.allclose(om2, -om2.T), ' omega fit vs formula', np.abs(om - om2).max(), np.abs(om + om2).max(),
      ' r', np.abs(rr - rr2).max(), np.abs(rr + rr2).max(), np.abs(rr - rr2 / 2).max(), np.abs(rr + rr2 / 2).max())
