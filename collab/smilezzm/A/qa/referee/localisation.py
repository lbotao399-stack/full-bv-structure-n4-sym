"""Referee check (2): the localisation constant c_loc and Lemma lem:local (i) at finite eps.

c_loc claimed = log(2/3)/(16 pi^2).
 (a) momentum-space radial integral  int d^4k/(2pi)^4 e^{-2k^2}(e^{-k^2}-1)/k^4
 (b) position-space form             int d^4w D_1(w) [D_2(w) - D_1(w)]   (D_t = D_eps at eps = t)
     (equal to (a) by Parseval; computed independently as a radial x-space integral).

Lemma (i):  LHS(eps) = int d^4z1 K_eps(x-z1) D_eps(z1) T(x,z1),  T(x,z1) = int d^4z2 D_eps(x-z2) D_eps(z1-z2) D_eps(z2),
claimed LHS = D0(x) [R(x) + D0(x) c_loc] + o(1).   Exact reduction used here (referee's own):
  LHS - D0 R = Term1 + Term2,
  Term1 = D0(x) int d^4z2 D_eps(z2) D_eps(x-z2) [D_{2eps}(x-z2) - D_eps(x-z2)]   (heat semigroup K_eps * D_eps = D_{2eps}; exact)
  Term2 = int d^4z1 K_eps(x-z1) [D_eps(z1) - D0(x)] T(x,z1)
        = sum_{l>=1} (2eps)^l int d^4z2 D_eps(z2) D_eps(x-z2) [H_l(grad) D_{2eps}](x-z2)   (+ O(e^{-1/4eps}))
  where D0(x+w) - D0(x) = sum_l H_l(w) (harmonic multipoles; Hobson: H(grad) K = (-1/2eps)^l H K).
  l = 1, 2 are computed (each next order is down by a power of eps); the series is asymptotic in eps.
Usage: python3 localisation.py [dps]
"""
import json
import os
import sys
import time

import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import Deps, D0, two_centre_F, fit_linear  # noqa: E402

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 20
mp.mp.dps = dps
pi = mp.pi
out = {}

# ---- c_loc ----------------------------------------------------------------
cl_claim = mp.log(mp.mpf(2) / 3) / (16 * pi**2)
cl_mom = mp.quad(lambda k: (2 * pi**2) / (2 * pi)**4 * k**3 * mp.e**(-2 * k * k) * mp.expm1(-k * k) / k**4, [0, 1, 3, mp.inf])


def Dt(w, t):
    return Deps(w, t)


cl_pos = mp.quad(lambda w: 2 * pi**2 * w**3 * Dt(w, 1) * (Dt(w, 2) - Dt(w, 1)), [0, 1, 4, 16, mp.inf])
out['c_loc'] = dict(claimed=mp.nstr(cl_claim, 20), momentum_radial=mp.nstr(cl_mom, 20),
                    position_space=mp.nstr(cl_pos, 20),
                    rel_diff_mom=mp.nstr((cl_mom - cl_claim) / cl_claim, 3),
                    rel_diff_pos=mp.nstr((cl_pos - cl_claim) / cl_claim, 3))
print(out['c_loc'], flush=True)


# ---- derivatives of D_t ------------------------------------------------------
def dD(r, t):
    a = 1 / (4 * t)
    E = mp.e**(-a * r * r)
    return (2 * a * E / r - 2 * (-mp.expm1(-a * r * r)) / r**3) / (4 * pi**2)


def ddD(r, t):
    a = 1 / (4 * t)
    E = mp.e**(-a * r * r)
    return (-4 * a * a * E - 6 * a * E / r**2 + 6 * (-mp.expm1(-a * r * r)) / r**4) / (4 * pi**2)


# self-test of the analytic derivatives
_r, _t = mp.mpf('0.37'), mp.mpf('0.05')
assert abs(dD(_r, _t) - mp.diff(lambda r: Deps(r, _t), _r)) < mp.mpf(10)**(-dps + 6)
assert abs(ddD(_r, _t) - mp.diff(lambda r: Deps(r, _t), _r, 2)) < mp.mpf(10)**(-dps + 8)

d0 = D0(1)
rows = []
for es in ['1e-2', '4e-3', '1e-3', '4e-4', '1e-4', '4e-5', '1e-5']:
    eps = mp.mpf(es)
    t0 = time.time()
    # r = |z2|, s = |x - z2|
    T1, e1 = two_centre_F(lambda r, s: Deps(r, eps) * Deps(s, eps) * (Deps(s, 2 * eps) - Deps(s, eps)), eps)
    T1 *= d0

    def F21(r, s):
        c = (1 - r * r + s * s) / (2 * s)          # xhat . (x - z2)/|x - z2|
        return Deps(r, eps) * Deps(s, eps) * (-2 / (4 * pi**2)) * c * dD(s, 2 * eps)
    T21, e21 = two_centre_F(F21, eps)
    T21 *= 2 * eps

    def F22(r, s):
        c = (1 - r * r + s * s) / (2 * s)
        return Deps(r, eps) * Deps(s, eps) * (4 * c * c - 1) * (ddD(s, 2 * eps) - dD(s, 2 * eps) / s) / (4 * pi**2)
    T22, e22 = two_centre_F(F22, eps)
    T22 *= (2 * eps)**2
    target = d0**2 * cl_claim
    tot = T1 + T21 + T22
    rows.append(dict(eps=es, Term1_over_D0sq=mp.nstr(T1 / d0**2, 15), Term2_l1_over_D0sq=mp.nstr(T21 / d0**2, 10),
                     Term2_l2_over_D0sq=mp.nstr(T22 / d0**2, 10),
                     LHS_minus_D0R_over_D0sq=mp.nstr(tot / d0**2, 15), claimed_c_loc=mp.nstr(cl_claim, 15),
                     difference=mp.nstr((tot - target) / d0**2, 8),
                     difference_over_epslogeps=mp.nstr((tot - target) / d0**2 / (eps * mp.log(eps)), 6),
                     quad_err=mp.nstr(max(e1, e21, e22), 3), seconds=round(time.time() - t0, 1)))
    print(rows[-1], flush=True)

# fit (no claimed value used): y = a + b eps log eps + c eps
X = [[1, mp.mpf(r['eps']) * mp.log(mp.mpf(r['eps'])), mp.mpf(r['eps'])] for r in rows]
Y = [mp.mpf(r['LHS_minus_D0R_over_D0sq']) for r in rows]
sol = fit_linear(X, Y)
X1 = [[1, mp.mpf(r['eps']) * mp.log(mp.mpf(r['eps'])), mp.mpf(r['eps'])] for r in rows]
Y1 = [mp.mpf(r['Term1_over_D0sq']) for r in rows]
sol1 = fit_linear(X1, Y1)
out['lemma_i'] = dict(rows=rows, fit_limit_total=mp.nstr(sol[0], 14), fit_limit_Term1=mp.nstr(sol1[0], 14),
                      claimed=mp.nstr(cl_claim, 14),
                      note='Term2 computed through its l=1,2 harmonic multipoles; l>=3 is O(eps^3 log eps)')
print(json.dumps({k: v for k, v in out['lemma_i'].items() if k != 'rows'}, indent=1))
with open(os.path.join(HERE, 'localisation_dps%d.json' % dps), 'w') as fh:
    json.dump(out, fh, indent=1)
