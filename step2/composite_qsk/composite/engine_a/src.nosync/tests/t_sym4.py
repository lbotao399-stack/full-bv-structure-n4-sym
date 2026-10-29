import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import cq_core as Q, cq_gval as GV, cq_eval as V, ea_ops as EO
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
rat=lambda x,y: None if (y[0]==0 and y[1]==0) else complex(complex(float(x[0]),float(x[1]))/complex(float(y[0]),float(y[1])))
for ta,tb,outs in [('lamP','lamM',['lamM','lamP']),('lamP','lamM',['lamP','lamM']),('lamP','lamP',['lamP','lamP'])]:
    za=rv(16); zb=rv(16); pols=[rv(16) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
    Xa=(ta,za,ONE); Xb=(tb,zb,ONE)
    f=GV.combine(GV.kernel(Xa,Xb,legs,'D0'),1); g=GV.combine(GV.kernel(Xb,Xa,legs[::-1],'D0'),1)
    l1=EO.pair_kernel_on_legs(Xa,Xb,legs); l2=EO.pair_kernel_on_legs(Xb,Xa,legs[::-1])
    print(ta,tb,outs,'composite',[rat(a,b) for a,b in zip(f,g) if not (b[0]==0 and b[1]==0)],'log',rat(l1,l2))
