import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_gcore as G, cq_eval as V, ea_kern as EK
rng=np.random.default_rng(7)
def rv(n): return np.array([complex(int(x),int(y)) for x,y in zip(rng.integers(-2,3,n), rng.integers(-2,3,n))])
ONE={(0,0,0,0):1}
RN={'AV':4,'AW':6,'lamP':16,'lamM':16,'c':1}
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
nchk=0; bad=0; t0=time.time()
for ta,tb,outs in [('AV','AV',['AV','AV']),('AW','lamP',['AW','lamP']),('AV','lamP',['lamP','AV'])]:
    za=rv(RN[ta]); zb=rv(RN[tb]); pols=[rv(RN[t]) for t in outs]
    legs=[(t,pl,m) for t,pl,m in zip(outs,pols,[p1,p2])]
    fa={(1,0,0,0):1j,(0,1,0,0):-1}; fb={(0,0,1,0):1j}
    for st in V.structures(legs):
        ob=V.blocks_of(st,legs); nb=len(ob)
        for v2pos in range(nb+1):
            for spos in range(nb+2):
                a=Q.config_value(ob,v2pos,spos,'P',ta,tb,fa=fa,fb=fb,v2dir=2,want=('fin','poles'))
                b=G.gconfig(ob,v2pos,spos,'P',ta,tb,fa=fa,fb=fb,v2dir=2,kinds=('fin',('spec',1)))
                if a is None and b is None: continue
                va=V.from_fin(a['fin'],za,zb) if a else V.vzero()
                fb_=b['fin'][0] if b else [None,None]
                vb=V.from_fin((fb_[0],fb_[1],b['E']),za,zb) if b else V.vzero()
                # cq_core has the vv direction sign too
                nchk+=1
                if V.vstr(va)!=V.vstr(vb): bad+=1; print('MISMATCH',ta,tb,outs,v2pos,spos,V.vstr(va),'|',V.vstr(vb))
print('fin checks',nchk,'bad',bad, '%.1fs'%(time.time()-t0))
# t-derivative vs finite differences for D0
Xa=('AV',rv(4),{(0,1,0,0):1j}); Xb=('AW',rv(6),ONE)
legs=[('AV',rv(4),p1),('AW',rv(6),p2)]
t0=time.time(); d_fd=V.pair_on_legs(Xa,Xb,legs,'D0'); t1=time.time()
print('D0 finite-diff', V.vstr(d_fd), '%.1fs'%(t1-t0))
import cq_gval as GV
t0=time.time(); d_t=GV.combine(GV.kernel(Xa,Xb,legs,'D0'),1); t1=time.time()
print('D0 t-deriv   ', V.vstr(d_t), '%.1fs'%(t1-t0), V.vstr(d_t)==V.vstr(d_fd))
# complex momenta: bilinear at null momenta
def u6(i):
    v=np.zeros(6,complex); v[i]=1; return v
q1=np.array([1,1j,0,0]); q2=np.array([0,0,1,-1j])
tot={'P':V.vzero(),'D0':V.vzero()}
Xa=('AW',u6(0),ONE); Xb=('AW',u6(1),ONE)
for lg in ([('AW',u6(0),q1),('AW',u6(1),q2)],[('AW',u6(1),q2),('AW',u6(0),q1)]):
    tot['P']=V.vadd(tot['P'],GV.combine(GV.kernel(Xa,Xb,lg,('P',1)),1))
    tot['D0']=V.vadd(tot['D0'],GV.combine(GV.kernel(Xa,Xb,lg,'D0'),1))
print('bilinear complex: P1', V.vstr(tot['P']), 'target -1/3 i (q1+q2)_1 =', -1/3*1j*(q1+q2)[1], ' D0', V.vstr(tot['D0']))
