import time, numpy as np
import sa_mode as MD, sa_ob as SO, sa_core as SA
from r2_core import Q2
t=time.time()
p=np.array([1,2,-1,1]); mm=MD.MM(p)
xi=SO.xi_generic('P')
rng=np.random.default_rng(3)
def rpol(t):
    B=SA.E.BASIS[t]; c=rng.integers(-2,3,size=len(B)); return sum(ci*np.asarray(v) for ci,v in zip(c,B)).astype(complex)
a=rpol('AV'); u=rpol('lamM'); aw=rpol('AW'); up=rpol('lamP'); s=rpol('lamsP')
T1=MD.F_T1(mm,xi); T2=MD.F_T1(mm,xi,range(4,10)); T3=MD.F_T3(mm,xi); T4=MD.F_T3(mm,xi,range(4,10)); T0=MD.F_T0(mm,xi)
g={k:MD.F_g(mm,xi,k) for k in (1,2,3,4)}
AA=MD.F_named(mm,'AA'); dsds=MD.F_named(mm,'dsAdsA'); ddcdA=MD.F_named(mm,'ddcdA')
LV=[('A',a),('l',u)]; LP=[('A',aw),('l',up)]; Lc=[('c',np.array([1])),('ls',s)]
print('{AA,T1}', mm.value(mm.br(AA,T1,0),LV), ' 2 GVA', 2*mm.value(MD.F_GA(mm,xi,range(4)),LV))
print('{AA,T2}', mm.value(mm.br(AA,T2,0),LP), ' 2 GPA', 2*mm.value(MD.F_GA(mm,xi,range(4,10)),LP))
for nm,F in (('T1',T1),('T3',T3)):
    print('Q',nm, mm.value(mm.Q(F),LV), 'g1', mm.value(g[1],LV))
for nm,F in (('T2',T2),('T4',T4)):
    print('Q',nm, mm.value(mm.Q(F),LP), 'g3', mm.value(g[3],LP))
QT0=mm.Q(T0); print('QT0 on V', mm.value(QT0,LV),'g2',mm.value(g[2],LV),' on c', mm.value(QT0,Lc),'g4',mm.value(g[4],Lc))
B=mm.br(dsds,T1,0); print('{dsds,T1} V', mm.value(B,LV),'g2',mm.value(g[2],LV))
B0=mm.br(ddcdA,T1,1); Lcl=[('c',np.array([1])),('l',u)]
print('{ddcdA,T1} on (c,lam)', mm.value(B0,Lcl), ' Q g2', mm.value(mm.Q(g[2]),Lcl), ' Q g4', mm.value(mm.Q(g[4]),Lcl), 'Qg1', mm.value(mm.Q(g[1]),Lcl))
print(time.time()-t)
