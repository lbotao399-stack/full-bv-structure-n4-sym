import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_eval as V, ea_ops as EO, ea_core as E, ea_kern as EK
Q.RULE['koszul']='none'
rng=np.random.default_rng(5)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
ROWN={'AV':4,'AW':6,'lamP':16,'lamM':16,'c':1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
ta,tb,outs='AV','lamP',['AV','lamP']
za=rv(4); zb=rv(16); pols=[rv(4), rv(16)]
legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
def show(Xa,Xb,legs):
    ta,za,fa=Xa; tb,zb,fb=Xb
    tot=V.vzero()
    for st in V.structures(legs):
        ob=V.blocks_of(st,legs); nb=len(ob)
        for v2pos in range(nb+1):
            for spos in range(nb+2):
                cfg=Q.make_config(ob,v2pos,spos,'P',1,None)
                if cfg is None: continue
                res=Q.config_value(ob,v2pos,spos,'P',ta,tb,v2dir=1,want=('fin',))
                if not res: continue
                v=V.from_fin(res['fin'],za,zb)
                print('  st',st,'v2pos',v2pos,'spos',spos,'syms',cfg[4],'lines',cfg[2],V.vstr(v)[:70])
                tot=V.vadd(tot,v)
    print(' total',V.vstr(tot))
show((ta,za,ONE),(tb,zb,ONE),legs)
show((tb,zb,ONE),(ta,za,ONE),legs[::-1])
