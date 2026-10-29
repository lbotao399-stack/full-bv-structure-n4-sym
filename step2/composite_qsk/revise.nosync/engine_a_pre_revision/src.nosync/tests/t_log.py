import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import cq_core as Q, cq_gval as GV, cq_eval as V, ea_ops as EO, ea_core as E
rng=np.random.default_rng(3)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
RN={'AV':4,'AW':6,'lamP':16,'lamM':16,'c':1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
for ta,tb,outs in [('AW','AW',['AW','AW']),('AW','lamP',['AW','lamP']),('AV','lamP',['lamP']),('lamP','lamP',['lamP','lamP']),('AV','c',['c']),('AV','AW',['AV','AW'])]:
    za=rv(RN[ta]); zb=rv(RN[tb]); pols=[rv(RN[t]) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
    a=EO.pair_kernel_on_legs((ta,za,{(1,0,0,0):1j}),(tb,zb,ONE),legs)
    b=GV.logkernel((ta,za,{(1,0,0,0):1j}),(tb,zb,ONE),legs)
    print(ta,tb,outs,E.qstr(a),'|',V.vstr(b), E.qstr(a)==V.vstr(b) or (a==(0,0) and V.viszero(b)))
