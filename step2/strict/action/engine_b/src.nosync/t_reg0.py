import time, sys
from fractions import Fraction as Fr
import numpy as np
t0=time.time()
import eb_core as E, eb_cycle as Y
print('import', time.time()-t0, flush=True)
for s in ['A','lam','c','As','lams','cs']:
    t1=time.time(); T=Y.tmpl_I3(s); print('tmpl', s, {k: v.shape for k,v in T.items()}, round(time.time()-t1,1), flush=True)
t1=time.time(); T=Y.tmpl_I4(); print('tmpl I4', {k: v.shape for k,v in T.items()}, round(time.time()-t1,1), flush=True)
E.edge_templates()
p=[Fr(1),Fr(2),Fr(-1),Fr(3)]
a=np.zeros(10,complex); a[:4]=[2,-1,1,1]
legs=[dict(sp='c',u=np.array([1],complex),p=p), dict(sp='A',u=a,p=[-x for x in p])]
st={}
t1=time.time()
res=Y.evaluate(legs, [], mlist=(0,1), cyclic=False, stats=st)
print('time', time.time()-t1, st)
for k,v in res.items():
    for m,x in v.items(): print(k,m,x)
p2=sum(x*x for x in p); pa=sum(p[m]*a[m] for m in range(4))
print('target eps^0: -27/8 i p^2 (p.a) =', -27/8*1j*float(p2)*pa)
