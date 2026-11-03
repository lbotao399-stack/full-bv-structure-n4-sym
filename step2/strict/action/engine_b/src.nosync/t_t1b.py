import time, random
from fractions import Fraction as Fr
import numpy as np
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_pos as P
random.seed(11)
def rv(n,lo=-3,hi=3): return [Fr(random.randint(lo,hi)) for _ in range(n)]
tdir=G.tmpl_dir(); tw=G.tmpl_weights(2,3,0,'D')
def gdir(mu,r):
    e=np.zeros(4,complex); e[mu]=1
    return Y.Gen('D',1,tdir,gvec=e,inj=r)
g0=Y.Gen('D0',1,tw)
def dot(a,b): return sum(complex(a[m])*complex(b[m]) for m in range(4))
p=rv(4); a1=np.zeros(10,complex); a2=np.zeros(10,complex); a1[:4]=rv(4,-2,2); a2[:4]=rv(4,-2,2)
def legs_at(r): return [dict(sp='A',u=a1,p=p),dict(sp='A',u=a2,p=[-p[m]-r[m] for m in range(4)])]
F1=lambda m,n: 1j*(complex(p[m])*a1[n]-complex(p[n])*a1[m])
F2=lambda m,n: 1j*(-complex(p[m])*a2[n]+complex(p[n])*a2[m])
dAdA=2*sum(F1(m,n)*F2(m,n) for m in range(4) for n in range(m+1,4))
print('dAdA value', dAdA)
for filt in ['line','special']:
    Y.PLACEMENT_FILTER=filt
    al,be,top=P.position_split(legs_at,gdir,None)
    a0=Y.evaluate(legs_at([Fr(0)]*4),[g0])
    print(filt,'x-part', {k:{m:complex(x) for m,x in v.items()} for k,v in al.items()}, 'V0 part', {k:{m:complex(x) for m,x in v.items()} for k,v in a0.items()})
    print('   x-part/dAdA', {m:complex(x)/dAdA for m,x in al['k'].items()} if 'k' in al else None, ' V0/dAdA', {m:complex(x)/dAdA for m,x in a0.get('k',{}).items()})
