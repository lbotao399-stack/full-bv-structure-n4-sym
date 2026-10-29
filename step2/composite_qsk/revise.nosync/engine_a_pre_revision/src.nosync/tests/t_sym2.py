import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import cq_core as Q, cq_eval as V, ea_ops as EO
Q.RULE['spec']='K'; V.VV['kappa']=int(sys.argv[1])
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
RN={'AV':4,'AW':6,'lamP':16,'lamM':16,'c':1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
def u6(i):
    v=np.zeros(6,complex); v[i]=1; return v
# bilinear GC3
tot={'P':V.vzero(),'D0':V.vzero()}
Xa=('AW',u6(0),ONE); Xb=('AW',u6(1),ONE)
for legs in ([('AW',u6(0),p1),('AW',u6(1),p2)],[('AW',u6(1),p2),('AW',u6(0),p1)]):
    tot['P']=V.vadd(tot['P'],V.pair_on_legs(Xa,Xb,legs,('P',0)))
    tot['D0']=V.vadd(tot['D0'],V.pair_on_legs(Xa,Xb,legs,'D0'))
print('bilinear P0',V.vstr(tot['P']),'D0',V.vstr(tot['D0']),' target -1/3*i, -8/3')
for ta,tb,outs in [('AV','lamP',['lamP']),('AV','AW',['AW']),('AW','lamP',['AW','lamP']),('lamP','lamP',['lamP','lamP']),('AV','lamP',['AV','lamP']),('AV','lamP',['lamP','AV']),('AV','AW',['AV','AW']),('AV','AV',['AV','AV']),('AW','AW',['AW','AW'])]:
    za=rv(RN[ta]); zb=rv(RN[tb]); pols=[rv(RN[t]) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
    Xa=(ta,za,ONE); Xb=(tb,zb,ONE)
    f=V.pair_on_legs(Xa,Xb,legs,('P',1)); g=V.pair_on_legs(Xb,Xa,legs[::-1],('P',1))
    l1=EO.pair_kernel_on_legs(Xa,Xb,legs); l2=EO.pair_kernel_on_legs(Xb,Xa,legs[::-1])
    rat=lambda x,y: None if (y[0]==0 and y[1]==0) else complex(complex(float(x[0]),float(x[1]))/complex(float(y[0]),float(y[1])))
    fr=[rat(a,b) for a,b in zip(f,g) if not (b[0]==0 and b[1]==0)]
    print(ta,tb,outs,'composite ratios',[complex(round(z.real,6),round(z.imag,6)) for z in fr],'log ratio',rat(l1,l2))
