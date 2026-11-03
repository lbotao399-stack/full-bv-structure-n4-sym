import time, numpy as np
from fractions import Fraction
import sa_core as SA, sa_ob as SO, sa_fit as SF
t=time.time()
rng=np.random.default_rng(7)
def rp():
    while True:
        p=rng.integers(-2,3,size=4)
        if p@p: return p
def rpol(t):
    B=SA.E.BASIS[t]; c=rng.integers(-2,3,size=len(B)); return sum(ci*np.asarray(v) for ci,v in zip(c,B)).astype(complex)
for chir in ('P','M'):
    xi=SO.xi_generic(chir); spec={'kind':'xi','xi':xi}
    o='M' if chir=='P' else 'P'
    secs=[('AV','lam'+o,[lambda L:SF.b_g1(L,xi),lambda L:SF.b_g2(L,xi)],[lambda L:SF.b_GVA(L,xi)]),
          ('AW','lam'+chir,[lambda L:SF.b_g3(L,xi)],[lambda L:SF.b_GPA(L,xi)]),
          ('c','lams'+chir,[lambda L:SF.b_g4(L,xi)],[lambda L:0*SF.b_g4(L,xi)])]
    for (t1,t2,bfs,pfs) in secs:
        S=[];P=[]
        for _ in range(4):
            p=rp(); L=[SA.leg(t1,rpol(t1),p), SA.leg(t2,rpol(t2),-p)]
            d=SO.density(spec,L,what=('fin',('pole',1)))
            S.append(([f(L) for f in bfs],d['fin']['alpha'].get('(Y0,Y1)',SA.vzero())))
            P.append(([f(L) for f in pfs],d[('pole',1)]['alpha'].get('(Y0,Y1)',SA.vzero())))
        c,r=SF.fit(S,len(bfs)); cp,_=SF.fit(P,len(pfs))
        print(chir,t1,t2,[SF.vstr(x) for x in c],'pole',[SF.vstr(x) for x in cp], 'rank',r)
    # zero sectors
    for (t1,t2) in [('AV','lam'+chir),('AW','lam'+o),('c','lams'+o),('AsV','c')]:
        p=rp(); L=[SA.leg(t1,rpol(t1),p), SA.leg(t2,rpol(t2),-p)]
        d=SO.density(spec,L); print('  zero?',t1,t2,SF.vstr(d['fin']['alpha'].get('(Y0,Y1)',SA.vzero())))
print(time.time()-t, SA.STATS)
