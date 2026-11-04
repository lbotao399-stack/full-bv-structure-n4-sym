import sys, time, cProfile, pstats
sys.argv = ['x', 'B', '+']
import xieta as XE, numpy as np, dtree as DT, fitlib as FL
from fractions import Fraction as Fr
raw = XE.run_loops()
key = list(raw)[3]; r = raw[key][0]; t0, t1 = key.split('|')
import eta1 as H
la = H.legs_at_of(t0, t1, H.cvec(r['u0']), H.cvec(r['u1']), [Fr(x) for x in r['p']])
ctx = dict(xi=H.cvec(r['xi']), eta=H.cvec(r['eta']))
c = XE.cands(XE.SECM2[1], None, 0)
print(len(c))
pr = cProfile.Profile(); pr.enable()
t = time.time(); DT.canon(XE.TT.Q(c[0][1](ctx)), la); print('one canon Q', time.time() - t)
pr.disable(); pstats.Stats(pr).sort_stats('cumtime').print_stats(8)
