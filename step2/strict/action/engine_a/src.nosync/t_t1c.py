import numpy as np
import sa_core as SA, sa_ob as SO, sa_fit as SF
spec={'kind':'Delta'}
for p,s in [((1,0,0,0),(1,0,0,0)),((1,0,0,0),(0,1,0,0)),((1,1,0,0),(1,0,0,0)),((1,2,-1,1),(0,1,2,-1))]:
    p=np.array(p); s=np.array(s+(0,)*6,dtype=complex)
    L=[SA.leg('c',[1],p), SA.leg('AsV',s,-p)]
    d=SO.density(spec,L)
    print(p,s[:4],SF.vstr(d['fin']['alpha'].get('(Y0,Y1)',SA.vzero())),[SF.vstr(d['fin']['beta'][m].get('(Y0,Y1)',SA.vzero())) for m in range(4)], SF.b_dcAs(L))
    d0=SA.cycle_sum(L,[SO.gens_of(spec)['0']],slotnames=['x'])
    print('   G0', [ (k,SF.vstr(v['fin'][0])) for k,v in d0.items()])
