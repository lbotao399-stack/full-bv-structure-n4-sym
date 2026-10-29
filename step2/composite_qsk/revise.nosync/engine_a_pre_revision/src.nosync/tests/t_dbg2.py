import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_eval as V, ea_ops as EO, ea_core as E, ea_kern as EK
Q.RULE['koszul']='none'
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1]); p3=np.array([2,1,0,-1])
for ta,tb,outs in [('AV','lamP',['AV','lamP']),('AV','lamP',['lamP','AV']),('AW','lamP',['AV','lamP']),('AV','AW',['AV','AW']),('AV','AV',['AV','AV'])]:
    RN={'AV':4,'AW':6,'lamP':16}
    za=rv(RN[ta]); zb=rv(RN[tb]); pols=[rv(RN[t]) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
    Xa=(ta,za,ONE); Xb=(tb,zb,ONE)
    for kind in [('Plog',1),('Ppole',1,1),('P',1)]:
        f=V.pair_on_legs(Xa,Xb,legs,kind); g=V.pair_on_legs(Xb,Xa,legs[::-1],kind)
        print(ta,tb,outs,kind,'a->b',V.vstr(f)[:50],'| b->a',V.vstr(g)[:50])
