"""Repeat the W-level checks of d1_checks.py (eq:d-D1a, eq:d-D1b, Ward identity) at a second value of hbar
and another random interaction J0, in system II, so that the powers of hbar in the closed forms are tested."""
import json, os, random, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from model import build
from d1_checks import System, wlevel
HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
for hb, seed in ((Fr(7, 3), 101), (Fr(-5, 2), 202)):
    m = build(N=4, gF=Fr(-2))
    sysm = System(m, 'II')
    alg = m['alg']; e, T = m['eps'], m['T']
    P = m['P'](e, T); Ps = m['Pstar'](e, T); Ke = m['Kt'](e); KT = m['Kt'](T); K0 = m['K0']
    SM = sysm.generator(K0, lambda f: sysm.D(m['M'], 0, f), -1)
    Sfree = sysm.generator(K0, lambda f: sysm.D(m['Q'], 1, f), 0)
    Cop = lambda f: alg.add(sysm.dphi(P, alg.mul(SM, f)), (-1, alg.mul(sysm.dphi(P, SM), f)), (-1, alg.mul(SM, sysm.dphi(P, f))))
    Dop = lambda f: alg.add(sysm.dphi(P, Cop(f)), (-1, Cop(sysm.dphi(P, f))))
    out['hbar = %s, seed %d' % (hb, seed)] = wlevel(m, sysm, SM, Sfree, Cop, Dop, P, Ps, Ke, KT, hb, 2, random.Random(seed))
with open(os.path.join(HERE, 'wlevel_second_hbar_results.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
print(json.dumps(out, indent=1))
