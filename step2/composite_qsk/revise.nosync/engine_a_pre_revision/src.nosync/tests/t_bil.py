import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_eval as V, ea_ops as EO
Q.RULE['koszul']=sys.argv[1] if len(sys.argv)>1 else 'front'
t0=time.time()
def u6(i):
    v=np.zeros(6,complex); v[i]=1; return v
p1=np.array([1,0,2,-1]); p2=np.array([0,1,-1,1])
ONE={(0,0,0,0):1}
for (a,b,c,d) in [(0,0,0,0),(0,1,0,1),(0,1,1,0),(0,0,1,1)]:
    Xa=('AW',u6(a),ONE); Xb=('AW',u6(b),ONE)
    res={}
    for kind in [('P',0),('P',2),'D0']:
        tot=V.vzero()
        for legs in ([('AW',u6(c),p1),('AW',u6(d),p2)], [('AW',u6(d),p2),('AW',u6(c),p1)]):
            tot=V.vadd(tot, V.pair_on_legs(Xa,Xb,legs,kind))
        res[str(kind)]=V.vstr(tot)
    Ostr = int(a==c)*int(b==d)+int(a==d)*int(b==c)
    print((a,b,c,d),'O=',Ostr, res, 'i(p1+p2)_0,2 =', (p1+p2)[0], (p1+p2)[2], '%.1fs'%(time.time()-t0), V.STATS)
