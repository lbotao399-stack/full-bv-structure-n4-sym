import time, numpy as np
from fractions import Fraction
import sa_core as SA
t=time.time()
p1 = np.array([1,2,-1,1]); p2=-p1
a = np.array([1,0,2,-1,0,0,0,0,0,0])
legs=[SA.leg('c',[1],p1), SA.leg('AV',a,p2)]
res = SA.cycle_sum(legs, [], what=('fin',('pole',1)))
for k,v in res.items():
    print(k, v['fin'][0], v[('pole',1)][0])
# basis int(d*dc, d*A_V): c odd first leg: (|p1|^2)(-i p2.a)
b = -1j*(p1@p1)*(p2@a)
print('basis', b, 'time', time.time()-t, SA.STATS)
