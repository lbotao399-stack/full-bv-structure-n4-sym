import numpy as np
exec(open('t_gs.py').read().split("cases=")[0])
QC.RULE['L0fix']=False
for ta,tb,lo in [('lamP','lamP',['AW','lamM']),('AW','lamP',['AW'])]:
    r=gs(ta,tb,lo,npts=2); print('L0fix=False',ta,tb,lo)
    for a,b in r: print('    ab:',a,' | ba:',b)
QC.RULE['L0fix']=True
Xa=('lamP',letvec('lamP'),None); Xb=('AW',letvec('AW'),None)
legs=[('lamM',legpol('lamM'),rng.integers(-2,3,4)) for _ in range(2)]
tot,parts=QC.pair_value(Xa,Xb,legs,gx,split=True)
print('Konishi split', {k:CV.vstr(v) for k,v in parts.items()})
