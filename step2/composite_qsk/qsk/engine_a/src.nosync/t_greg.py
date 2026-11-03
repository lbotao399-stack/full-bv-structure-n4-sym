import time, numpy as np
from fractions import Fraction
t0=time.time()
import qa_s1 as S1, qa_core as QC
import cq_eval as CV
I6=np.eye(6,dtype=int)
gD0 = QC.even_gen('D0', dv=8)
gDs = [QC.even_gen('Ds%d'%mu, a=2*np.eye(4)[mu]) for mu in range(4)]
q1 = np.array([1,0,2,-1]); q2 = np.array([0,1,-1,1])
def kern(I,J,K,L, which):
    Xa=('AW',I6[I],None); Xb=('AW',I6[J],None)
    tot = CV.vzero()
    for (pk,qa),(pl,qb) in [((K,q1),(L,q2)),((L,q2),(K,q1))]:
        legs=[('AW',I6[pk],qa),('AW',I6[pl],qb)]
        if which=='D0':
            v = QC.pair_value(Xa,Xb,legs,gD0)['fin']
            for mu in range(4):
                f = lambda r, mu=mu: QC.pair_value(Xa,Xb,legs,gDs[mu],r0=r)['fin']
                d = QC.deriv_r(f, mu, 2)
                v = CV.vadd(v, CV.vscale(d,(Fraction(0),Fraction(-1))))
        else:
            mu = which
            v = QC.pair_value(Xa,Xb,legs,gDs[mu])['fin']
        tot = CV.vadd(tot, v)
    return tot
for comp in [(0,1,0,1),(0,0,1,1),(0,0,0,0),(0,1,1,0),(0,1,2,3)]:
    print(comp, 'D0:', CV.vstr(kern(*comp,'D0')), ' Ds0:', CV.vstr(kern(*comp,0)), ' Ds2:', CV.vstr(kern(*comp,2)))
print('time',time.time()-t0, QC.STATS)
