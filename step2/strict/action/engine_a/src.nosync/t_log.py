import numpy as np
import sa_core as SA, sa_fit as SF
p=np.array([1,2,-1,1]); aw=np.array([0,0,0,0,1,2,0,-1,0,1],dtype=complex); bw=np.array([0,0,0,0,0,1,1,0,2,-1],dtype=complex)
L=[SA.leg('AW',aw,p),SA.leg('AW',bw,-p)]
r=SA.cycle_sum(L,[],what=('log',('eps',1)))
for k,v in r.items(): print(k, SF.vstr(v['log'][0]), SF.vstr(v[('eps',1)][0]), 'dAdA', SF.b_dAWdAW(L), 'AA', SF.b_AA(L))
