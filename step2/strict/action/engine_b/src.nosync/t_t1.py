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
for trial in range(2):
    p=rv(4); a1=np.zeros(10,complex); a2=np.zeros(10,complex); a1[:4]=rv(4,-2,2); a2[:4]=rv(4,-2,2)
    def legs_at(r): return [dict(sp='A',u=a1,p=p),dict(sp='A',u=a2,p=[-p[m]-r[m] for m in range(4)])]
    t=time.time(); st={}
    al,be,top=P.position_split(legs_at,gdir,g0,stats=st)
    print('time',round(time.time()-t,1),st)
    for k,v in al.items():
        for m,x in v.items(): print(' alpha',k,m,x, complex(x))
    for mu,d in be.items():
        for k,v in d.items():
            for m,x in v.items(): print(' beta',mu,k,m,complex(x))
    for k,v in top.items():
        for m,x in v.items(): print(' top',k,m,complex(x))
    p2=dot(p,p); a12=dot(a1,a2); pa1=dot(p,a1); pa2=dot(p,a2)
    # (dA,dA) value: 2 * sum_{M<N} F1 F2, F_j,mn = i(p_m a_n - p_n a_m), p2=-p
    dAdA=2*(-(-1)*(p2*a12 - pa1*pa2))  # computed below directly
    F1=lambda m,n: 1j*(complex(p[m])*a1[n]-complex(p[n])*a1[m])
    F2=lambda m,n: 1j*(-complex(p[m])*a2[n]+complex(p[n])*a2[m])
    dAdA=2*sum(F1(m,n)*F2(m,n) for m in range(4) for n in range(m+1,4))
    dsds=2*(-1j*pa1)*(-1j*(-pa2))*1  # d*A = -d.A -> -(i p.a)
    dsds=2*(-(1j*pa1))*(-(1j*(-pa2)))
    print(' target -11/3 (dA,dA) =', -11/3*dAdA, ' (d*A,d*A) value', dsds, ' pole target 9u int(A,A) =', 18*a12)
