import time, numpy as np
from fractions import Fraction
import sa_core as SA, ea_core as E
t=time.time()
SB=np.array(E.SGB)
gs0 = SA.gen_even('s0', dv=8.0)
gmu = [SA.gen_even('s%d'%m, a=2*np.eye(4)[m]) for m in range(4)]
def delta_cycles(legs):
    G0 = SA.cycle_sum(legs, [gs0], slotnames=['D'])
    alpha = SA.vzero(); beta=[SA.vzero() for _ in range(4)]
    for k,v in G0.items(): alpha = SA.vadd(alpha, v['fin'][0])
    for m in range(4):
        L2 = [dict(l) for l in legs]; L2[-1] = dict(L2[-1]); L2[-1]['p1'] = -np.eye(4,dtype=int)[m]
        Gm = SA.cycle_sum(L2, [gmu[m]], slotnames=['D'], r0=None, r1=np.eye(4,dtype=int)[m])
        for k,v in Gm.items():
            alpha = SA.vadd(alpha, SA.vscale(v['fin'][1], (Fraction(0),Fraction(-1))))
            beta[m] = SA.vadd(beta[m], SA.vscale(v['fin'][0], (Fraction(0),Fraction(-1))))
    return alpha, beta
p = np.array([1,2,-1,1]); q=-p
a = np.array([1,0,2,-1,0,0,0,0,0,0]); b=np.array([0,1,1,2,0,0,0,0,0,0])
al,be = delta_cycles([SA.leg('AV',a,p), SA.leg('AV',b,q)])
print('AVAV alpha', al[0], al[1:3], 'beta', [x[0] for x in be])
dAdA = -2*((p@q)*(a@b)-(p@b[:4])*(q@a[:4])); dsds=-2*(p@a[:4])*(q@b[:4])
print('dAdA',dAdA,'d*d*',dsds)
aw = np.array([0,0,0,0,1,2,0,-1,0,1]); bw=np.array([0,0,0,0,0,1,1,0,2,-1])
al2,be2 = delta_cycles([SA.leg('AW',aw,p), SA.leg('AW',bw,q)])
print('AWAW alpha', al2[0], 'beta', [x[0] for x in be2], 'dAdA', -2*(p@q)*(aw@bw))
print('time', time.time()-t, SA.STATS)
