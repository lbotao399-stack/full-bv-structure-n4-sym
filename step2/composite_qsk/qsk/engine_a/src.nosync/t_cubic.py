import numpy as np
from fractions import Fraction
exec(open('t_gs.py').read().split("def gs(")[0])
for trial in range(3):
    Xa=('AW',letvec('AW'),None); Xb=('lamP',letvec('lamP'),None)
    q=rng.integers(-2,3,4)
    legs=[('AW',legpol('AW'),q)]
    tot,parts=QC.pair_value(Xa,Xb,legs,gx,what=('fin',('pole',1)),split=True)
    f=parts.get(('fin',0,'S3','special')); p=parts.get((('pole',1),0,'S3','special'))
    qq=int(q@q)
    print('q2',qq,'fin',CV.vstr(f),'pole',CV.vstr(p), ' fin/pole =', [ (complex(*map(float,f[0]))/complex(*map(float,p[0]))) ], 'target 3/4 q2 =', 0.75*qq)
    print('   total fin', CV.vstr(tot['fin']), ' total pole', CV.vstr(tot[('pole',1)]))
