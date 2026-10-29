import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction
import numpy as np, mpmath as mp
t0=time.time()
import cq_core as Q
print('import', time.time()-t0)
# 1. fin_monomial vs numerics (E=3)
def num_fin(exps, s):
    E=len(exps)
    if E==3:
        others=[g for g in range(3) if g!=s]
        def f(a1, a2):
            a=[0,0,0]
            # coordinates: alpha_0 = 1-a1-a2 etc: use alpha_1=a1, alpha_2=a2
            a=[1-a1-a2, a1, a2]
            if min(a) < 0: return 0
            m=min(a[g] for g in others)
            if a[s] > m: return 0
            return mp.mpf(a[0])**exps[0]*mp.mpf(a[1])**exps[1]*mp.mpf(a[2])**exps[2]*mp.log(m/a[s])
        return mp.quad(lambda x: mp.quad(lambda y: f(x,y), [0, (1-x)/3, (1-x)/2, 1-x]), [0, mp.mpf(1)/3, mp.mpf(1)/2, 1])
for exps, s in [((0,0,0),0), ((1,0,0),0), ((0,2,1),0), ((1,1,0),1), ((2,0,1),2)]:
    r0, r1 = Q.fin_monomial(exps, s)
    ex = float(r0) + float(r1)*mp.log(mp.mpf(3)/2)
    nm = num_fin(exps, s)
    print(exps, s, r0, r1, float(ex), float(nm))
# 2. Lemma tl-special example: k=1, N=-|l|^2, edge 0 special, edge 1 carries l - p
p = np.array([1,2,-1,1])
S = Q.CSimplex([np.zeros(4,int), -p])
N = {}
for mu in range(4):
    e=[0,0,0,0]; e[mu]=2; N[tuple(e)] = complex(-1)
f0, f1 = S.fin(N, 0)
sp_ = S.special(N, 0)
print('fin', f0, f1, 'target -3/8 |p|^2 =', Fraction(-3,8)*int(p@p), 'poles', sp_)
