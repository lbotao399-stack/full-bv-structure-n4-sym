import time, numpy as np
from fractions import Fraction
import sa_core as SA, sa_ob as SO, sa_fit as SF
t=time.time()
rng=np.random.default_rng(5)
def rp(): 
    while True:
        p=rng.integers(-2,3,size=4)
        if p@p: return p
def rpol(t):
    B=SA.E.BASIS[t]; c=rng.integers(-2,3,size=len(B)); return sum(ci*np.asarray(v) for ci,v in zip(c,B)).astype(complex)
spec={'kind':'Delta'}
def samp(t1,t2,n,bfs):
    S=[];P=[]
    for _ in range(n):
        p=rp(); L=[SA.leg(t1,rpol(t1),p), SA.leg(t2,rpol(t2),-p)]
        d=SO.density(spec,L,what=('fin',('pole',1)))
        a=d['fin']['alpha'].get('(Y0,Y1)',SA.vzero()); pa=d[('pole',1)]['alpha'].get('(Y0,Y1)',SA.vzero())
        bnz=[not SA.viszero(d['fin']['beta'][m].get('(Y0,Y1)',SA.vzero())) for m in range(4)]
        S.append(([f(L) for f in bfs],a)); P.append(([f(L) for f in bfs]+[SF.b_AA(L)] if t1[0]=='A' and t2[0]=='A' else [f(L) for f in bfs],pa))
        if any(bnz): print('beta nonzero', t1,t2)
    return S,P
for (t1,t2,bfs) in [('AV','AV',[SF.b_dAVdAV,SF.b_dsAdsA]),('AW','AW',[SF.b_dAWdAW]),('lamP','lamM',[SF.b_ldl]),('c','AsV',[SF.b_dcAs])]:
    S,P=samp(t1,t2,4,bfs)
    c,r=SF.fit(S,len(bfs)); print(t1,t2,[SF.vstr(x) for x in c],'rank',r)
    try:
        cp,r=SF.fit(P,len(P[0][0])); print('  pole',[SF.vstr(x) for x in cp])
    except Exception as e: print(' pole fit fail',e)
print(time.time()-t, SA.STATS)
