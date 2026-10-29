import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_eval as V, ea_ops as EO, ea_core as E
Q.RULE['koszul']=sys.argv[1] if len(sys.argv)>1 else 'none'
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
ROWN={'AV':4,'AW':6,'lamP':16,'lamM':16,'c':1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1]); p3=np.array([2,-1,0,1])
cases=[('AV','lamP',['lamP']), ('AV','AW',['AW']), ('AW','lamP',['AW','lamP']), ('AW','lamP',['lamP','AW']), ('lamP','lamP',['lamP','lamP']), ('lamP','lamM',['lamP','lamM']),('AV','lamP',['AV','lamP']),('AV','c',['c'])]
for ta,tb,outs in cases:
    za=rv(ROWN[ta]); zb=rv(ROWN[tb])
    moms=[p1,p2,p3][:len(outs)]
    pols=[rv(ROWN[t]) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,moms)]
    Xa=(ta,za,ONE); Xb=(tb,zb,ONE)
    f=V.pair_on_legs(Xa,Xb,legs,('P',1)); g=V.pair_on_legs(Xb,Xa,legs[::-1],('P',1))
    l1=EO.pair_kernel_on_legs(Xa,Xb,legs); l2=EO.pair_kernel_on_legs(Xb,Xa,legs[::-1])
    def ratio(x,y):
        x=complex(float(x[0]),float(x[1])); y=complex(float(y[0]),float(y[1]))
        return None if abs(y)<1e-14 else complex(round((x/y).real,9), round((x/y).imag,9))
    print(ta,tb,outs,'composite f/g:', ratio(f[0],g[0]), V.vstr(f)[:60], '| log l1/l2:', ratio(l1,l2), E.qstr(l1)[:40])
