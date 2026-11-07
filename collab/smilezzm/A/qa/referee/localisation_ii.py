"""Referee check (2, supplementary): Lemma lem:local (ii) at finite eps.

LHS(ii) = int d^4z1 d^4z2 D(x-z1) D(x-z2) K_eps(z1-z2) D(z1) D(z2),  claimed  -> J(x) + 2 D0(x)^2 c_loc.
With P(z) = D_eps(x-z) D_eps(z):  LHS(ii) - J = <P, (K_eps * P) - P>, symmetric under z -> x - z, hence
  = 2 int_{|x-z2| < |z2|} P(z2) [ (K_eps * P)(z2) - P(z2) ] d^4z2.
In that half space |z2| >= 1/2, and (K_eps * P)(z2) = int K(w) D_eps(x-z2-w) D_eps(z2+w) dw is evaluated by the
harmonic multipole expansion of D_0(z2+w) around z2 (D_eps = D_0 there up to e^{-1/16eps}) and the heat semigroup:
  (K*P)(z2) = D0(z2) D_{2eps}(x-z2) + sum_{l>=1} (2eps)^l [H_l^{(z2)}(grad) D_{2eps}](z2-x),   l = 1, 2 kept.
Usage: python3 localisation_ii.py [dps]
"""
import json
import os
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import Deps, D0, kernel, small_breaks, fit_linear  # noqa: E402
from localisation_helpers import dD, ddD  # noqa: E402

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
mp.mp.dps = dps
pi = mp.pi
cl = mp.log(mp.mpf(2) / 3) / (16 * pi**2)
d0 = D0(1)


def half_space(F, eps):
    """int over {s < r} of d^4z F(r, s), r = |z|, s = |x - z|, |x| = 1."""
    half = mp.mpf(1) / 2

    def inner(s):
        lo = max(s, abs(1 - s))
        return mp.quad(lambda r: kernel(r, s) * F(r, s), [lo, 1 + s])
    pts = small_breaks(eps, half) + [mp.mpf(1), 2, 3, 5, 10, 30, 100, mp.inf]
    tot = mp.mpf(0)
    for a, b in zip(pts[:-1], pts[1:]):
        tot += mp.quad(inner, [a, b])
    return tot


rows = []
for es in ['1e-2', '4e-3', '1e-3', '4e-4', '1e-4', '4e-5', '1e-5']:
    eps = mp.mpf(es)
    t0 = time.time()

    def F(r, s):
        P = Deps(r, eps) * Deps(s, eps)
        KP0 = D0(r) * Deps(s, 2 * eps)
        c1 = (1 - r * r - s * s) / (2 * r * s)            # zhat2 . (x - z2)/|x - z2|
        # l = 1: H_1(w) = -2 (zhat2.w)/(4 pi^2 r^3), evaluated at the argument z2 - x = -(x - z2)
        KP1 = 2 * eps * (2 / (4 * pi**2 * r**3)) * c1 * dD(s, 2 * eps)
        # l = 2: H_2(w) = (4 (zhat2.w)^2 - w^2)/(4 pi^2 r^4)
        KP2 = (2 * eps)**2 * (4 * c1 * c1 - 1) * (ddD(s, 2 * eps) - dD(s, 2 * eps) / s) / (4 * pi**2 * r**4)
        return P * (KP0 + KP1 + KP2 - P)
    val = 2 * half_space(F, eps)
    rows.append(dict(eps=es, LHSii_minus_J_over_D0sq=mp.nstr(val / d0**2, 15), claimed_2c_loc=mp.nstr(2 * cl, 15),
                     difference=mp.nstr(val / d0**2 - 2 * cl, 8),
                     difference_over_epslogeps=mp.nstr((val / d0**2 - 2 * cl) / (eps * mp.log(eps)), 6),
                     seconds=round(time.time() - t0, 1)))
    print(rows[-1], flush=True)
X = [[1, mp.mpf(r['eps']) * mp.log(mp.mpf(r['eps'])), mp.mpf(r['eps'])] for r in rows]
Y = [mp.mpf(r['LHSii_minus_J_over_D0sq']) for r in rows]
sol = fit_linear(X, Y)
out = dict(rows=rows, fit_limit=mp.nstr(sol[0], 14), claimed=mp.nstr(2 * cl, 14),
           note='(K*P) through l<=2 harmonic multipoles; neglected l>=3 terms are O(eps^3 log eps)')
print(out['fit_limit'], out['claimed'])
with open(os.path.join(HERE, 'localisation_ii_dps%d.json' % dps), 'w') as fh:
    json.dump(out, fh, indent=1)
