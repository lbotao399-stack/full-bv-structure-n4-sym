"""Referee check (1): the position-space masters J(x) and R(x) by DIRECT position-space quadrature.

J(x) = int d^4z D_eps(x-z)^2 D_eps(z)^2,   R(x) = int d^4z D_eps(x-z)^2 D_eps(z),  |x| = 1.
Method: exact two-centre reduction to a 2d integral over r=|z|, s=|x-z| (common.two_centre),
mpmath tanh-sinh with breakpoints at the scale sqrt(eps); no Schwinger times, no momentum space.
Claims under review (ip1_theory.tex):
  J = D0^2/(8 pi^2) (log(x^2/eps) + gamma_E - 3 log 2) + o(1)
  R = D0  /(16pi^2) (log(x^2/eps) + 1 + gamma_E - 3 log 2) + o(1)
Usage: python3 xspace_masters.py [J|R] [dps]   (writes xspace_<which>.json next to this file)
"""
import json
import os
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import Deps, D0, two_centre, fit_linear  # noqa: E402

which = sys.argv[1] if len(sys.argv) > 1 else 'J'
dps = int(sys.argv[2]) if len(sys.argv) > 2 else 20
mp.mp.dps = dps
EPS = ['1e-3', '5e-4', '2e-4', '1e-4', '5e-5', '2e-5', '1e-5']
g = mp.euler
L2, L3 = mp.log(2), mp.log(3)
d0 = D0(1)

rows = []
for es in EPS:
    eps = mp.mpf(es)
    t = time.time()
    if which == 'J':
        val, err = two_centre(lambda r: Deps(r, eps)**2, lambda s: Deps(s, eps)**2, eps)
        norm = val * 8 * mp.pi**2 / d0**2           # claimed: -log eps + g - 3 log 2
        claim = -mp.log(eps) + g - 3 * L2
    else:
        val, err = two_centre(lambda r: Deps(r, eps), lambda s: Deps(s, eps)**2, eps)
        norm = val * 16 * mp.pi**2 / d0              # claimed: -log eps + 1 + g - 3 log 2
        claim = -mp.log(eps) + 1 + g - 3 * L2
    diff = norm - claim
    rows.append(dict(eps=es, value=mp.nstr(val, 18), normalised=mp.nstr(norm, 18),
                     claimed_normalised=mp.nstr(claim, 18), diff=mp.nstr(diff, 10),
                     diff_over_epslogeps=mp.nstr(diff / (eps * mp.log(eps)), 8),
                     quad_err=mp.nstr(err, 3), seconds=round(time.time() - t, 1)))
    print(rows[-1], flush=True)

# unbiased fit (no claimed value used): normalised = a log eps + b + c eps log eps + d eps (+ e eps^2 log eps + f eps^2)
fits = {}
for npar in (4, 6):
    X, Y = [], []
    for rw in rows:
        e = mp.mpf(rw['eps'])
        basis = [mp.log(e), 1, e * mp.log(e), e, e * e * mp.log(e), e * e][:npar]
        X.append(basis)
        Y.append(mp.mpf(rw['normalised']))
    sol = fit_linear(X, Y)
    fits['%d-param' % npar] = dict(a_logeps=mp.nstr(sol[0], 14), b_const=mp.nstr(sol[1], 14),
                                    c_epslogeps=mp.nstr(sol[2], 8), d_eps=mp.nstr(sol[3], 8))
b_claim = (g - 3 * L2) if which == 'J' else (1 + g - 3 * L2)
out = dict(which=which, dps=dps, method='direct position-space quadrature, two-centre (r,s) reduction, mpmath',
           rows=rows, fits=fits, a_claimed=-1, b_claimed=mp.nstr(b_claim, 16))
print(json.dumps(fits, indent=1))
print('claimed b =', b_claim)
with open(os.path.join(HERE, 'xspace_%s_dps%d.json' % (which, dps)), 'w') as fh:
    json.dump(out, fh, indent=1)
