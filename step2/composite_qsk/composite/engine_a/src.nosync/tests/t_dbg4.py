import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import cq_core as Q, cq_eval as V
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
RN={'AV':4,'AW':6}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
ta,tb,outs='AV','AW',['AV','AW']
za=rv(4); zb=rv(6); pols=[rv(4),rv(6)]
legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
def P(l,n): return n if l==0 else l-1
def vals(Xa,Xb,legs,kind):
    ta,za,fa=Xa; tb,zb,fb=Xb
    out={}
    for st in V.structures(legs):
        ob=V.blocks_of(st,legs); nb=len(ob); n=nb+1
        for v2pos in range(nb+1):
            for spos in range(nb+2):
                res=Q.config_value(ob,v2pos,spos,kind,ta,tb,want=('fin',))
                if not res: continue
                v=V.from_fin(res['fin'],za,zb)
                if V.viszero(v): continue
                out[(len(st),v2pos,P(spos,n))]=(v,n)
    return out
A=vals((ta,za,ONE),(tb,zb,ONE),legs,'s0'); B=vals((tb,zb,ONE),(ta,za,ONE),legs[::-1],'s0')
for k,(v,n) in sorted(A.items()):
    nst,v2,sp=k; nb=n-1
    k2=(nst,nb-v2,n-sp)
    w=B.get(k2)
    print(k, V.vstr(v)[:50], '| mirror', k2, V.vstr(w[0])[:50] if w else None)
print('B only', [k for k in B if (k[0], (k[0]==1 and 1 or 2)-k[1], 0) is None])
for k,(v,n) in sorted(B.items()): print(' B',k,V.vstr(v)[:50])
