import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import cq_core as Q, cq_eval as V, ea_ops as EO
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
RN={'AV':4,'AW':6,'lamP':16,'lamM':16,'c':1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
rat=lambda x,y: None if (y[0]==0 and y[1]==0) else complex(complex(float(x[0]),float(x[1]))/complex(float(y[0]),float(y[1])))
for ta,tb,outs in [('AV','AW',['AV','AW']),('AV','AV',['AV','AV']),('AV','AV',['AW','AW'])]:
    za=rv(RN[ta]); zb=rv(RN[tb]); pols=[rv(RN[t]) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
    Xa=(ta,za,ONE); Xb=(tb,zb,ONE)
    for kind,args in [('s0',{}),('s',{'v2dir':1,'r':np.array([0,0,0,0])}),('s',{'v2dir':1,'r':np.array([0,1,0,0])})]:
        f=V._raw(Xa,Xb,legs,kind,**args); 
        args2=dict(args)
        g=V._raw(Xb,Xa,legs[::-1],kind,**args2)
        print(ta,tb,outs,kind,args.get('r'),[rat(a,b) for a,b in zip(f,g) if not (b[0]==0 and b[1]==0)])
