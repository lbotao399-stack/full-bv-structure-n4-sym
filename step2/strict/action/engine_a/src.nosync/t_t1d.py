import numpy as np
import sa_core as SA, sa_ob as SO, sa_fit as SF
rng=np.random.default_rng(5)
spec={'kind':'Delta'}
for _ in range(6):
    p=rng.integers(-2,3,size=4); 
    if not p@p: continue
    s=np.zeros(10,dtype=complex); s[:4]=rng.integers(-2,3,size=4); c=rng.integers(-2,3)
    L=[SA.leg('c',[c],p), SA.leg('AsV',s,-p)]
    d=SO.density(spec,L)
    print(p,c,s[:4],SF.vstr(d['fin']['alpha'].get('(Y0,Y1)',SA.vzero())), SF.b_dcAs(L))
