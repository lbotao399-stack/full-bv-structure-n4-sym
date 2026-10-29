import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_eval as V, ea_ops as EO, ea_core as E, ea_kern as EK
Q.RULE['koszul']='none'
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
RN={'AV':4,'AW':6,'lamP':16}
def L(p,n): return p+1 if p<n else 0
def P(l,n): return n if l==0 else l-1
def vals(Xa,Xb,legs):
    ta,za,fa=Xa; tb,zb,fb=Xb
    out={}
    for st in V.structures(legs):
        ob=V.blocks_of(st,legs); nb=len(ob); n=nb+1
        for v2pos in range(nb+1):
            for spos in range(nb+2):
                res=Q.config_value(ob,v2pos,spos,'P',ta,tb,v2dir=1,want=('fin',))
                if not res: continue
                v=V.from_fin(res['fin'],za,zb)
                if V.viszero(v): continue
                out[(tuple(map(tuple,st)),v2pos,P(spos,n))]=v
    return out,n
for ta,tb,outs in [('AV','AW',['AV','AW']),('AV','AV',['AV','AV'])]:
    za=rv(RN[ta]); zb=rv(RN[tb]); pols=[rv(RN[t]) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
    A,n=vals((ta,za,ONE),(tb,zb,ONE),legs); B,_=vals((tb,zb,ONE),(ta,za,ONE),legs[::-1])
    print(ta,tb,outs,len(A),len(B))
    for (st,v2,sp),v in A.items():
        nb=n-1
        # mirror: structure reversed, v2pos -> nb - v2pos, line position p -> n - p
        r=len(legs)
        st2=tuple(sorted(tuple(sorted(r-1-i for i in b)) for b in st))
        key=None
        for k2 in B:
            if tuple(sorted(k2[0]))==st2 and k2[1]==nb-v2 and k2[2]==n-sp: key=k2
        w=B.get(key)
        rat=None
        if w is not None:
            x=complex(float(v[1][0]),float(v[1][1])); y=complex(float(w[1][0]),float(w[1][1]))
            rat=None if abs(y)<1e-12 else round((x/y).real,6)
        print('  st',st,'v2pos',v2,'spos(chainpos)',sp,'ratio log2-part',rat, V.vstr(v)[:40])
    print(' A keys', sorted(A.keys())); print(' B keys', sorted(B.keys()))
    for k,v in sorted(B.items()): print('   B',k,V.vstr(v)[:60])
