import sys, numpy as np
sys.argv = ['x', sys.argv[1], '+']
import xieta as XE, eta1 as H
rng = np.random.default_rng(4)
for ech in ('+', '-'):
    ett = ({'+': 'ls+', '-': 'ls-'} if XE.ROUTE == 'B' else {'+': 'lamsP', '-': 'lamsM'})[ech]
    xi, eta = H.rand_u(XE.XIT, rng), H.rand_u(ett, rng)
    r = XE.bracket_fit(xi, eta, rng)
    if r is None:
        print(ech, 'NO FIT'); continue
    cD, om, rr = r
    pair = (xi @ H.MP_ @ eta)
    print(ech, 'c_D', cD, ' (xi,eta)/2 =', pair / 2, ' |omega|', np.abs(om).max(), ' |r|', np.abs(rr).max())
