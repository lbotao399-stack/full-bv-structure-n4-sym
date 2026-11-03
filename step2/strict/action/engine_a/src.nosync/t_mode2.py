import numpy as np
import sa_mode as MD, sa_ob as SO, sa_core as SA
p=np.array([1,2,-1,1]); mm=MD.MM(p); xi=SO.xi_generic('P')
dsds=MD.F_named(mm,'dsAdsA'); T1=MD.F_T1(mm,xi)
u=SA.E.BASIS['lamM'][0].astype(complex)+2*SA.E.BASIS['lamM'][3]
L=[('c',np.array([1])),('A',np.array([1,0,2,-1,0,0,0,0,0,0],dtype=complex))]
a=mm.Q(dsds); b=mm.br(mm.Sfree,dsds,0)
print('Q_apply vs {Sfree,.} on dsds:', mm.value(a,L), mm.value(b,L))
a=mm.Q(T1); b=mm.br(mm.Sfree,T1,0); L2=[('A',np.array([1,0,2,-1,0,0,0,0,0,0],dtype=complex)),('l',u)]
print('on T1:', mm.value(a,L2), mm.value(b,L2))
