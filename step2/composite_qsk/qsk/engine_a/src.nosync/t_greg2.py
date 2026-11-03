import numpy as np, time
from fractions import Fraction
exec(open('t_gs.py').read().split("B = S1.chir_basis")[0])
import ea_ops as EO
t0=time.time()
gD0 = QC.even_gen('D0', dv=8)
gDs = [QC.even_gen('Ds%d'%mu, a=2*np.eye(4)[mu]) for mu in range(4)]
def ups0(Xa,Xb,legs):
    v = QC.pair_value(Xa,Xb,legs,gD0)['fin']
    for mu in range(4):
        f = lambda r, mu=mu: QC.pair_value(Xa,Xb,legs,gDs[mu],r0=r)['fin']
        v = CV.vadd(v, CV.vscale(QC.deriv_r(f, mu, 2),(Fraction(0),Fraction(-1))))
    return v
def ref(Xa,Xb,legs):
    lc = EO.pair_kernel_on_legs((Xa[0],Xa[1],{(0,0,0,0):1}),(Xb[0],Xb[1],{(0,0,0,0):1}),legs,kind='log')
    c0 = CV.pair_on_legs(Xa,Xb,legs,'D0')
    intr = -1 if (QC.Q.PARITY_T[Xa[0]]+QC.Q.PARITY_T[Xb[0]])%2 else 1
    v = CV.vscale(c0,(Fraction(intr),Fraction(0)))
    v[0] = (v[0][0]+4*lc[0], v[0][1]+4*lc[1])
    return v, lc
for ta,tb,lo in [('AW','lamP',['AW','lamP']),('AW','lamP',['lamP','AW']),('lamP','lamP',['lamP','lamP']),('lamP','lamM',['lamP','lamM']),('AW','lamM',['lamM','AW'])]:
    for _ in range(2):
        Xa=(ta,letvec(ta),None); Xb=(tb,letvec(tb),None)
        legs=[(t,legpol(t),rng.integers(-2,3,4)) for t in lo]
        a=ups0(Xa,Xb,legs); b,lc=ref(Xa,Xb,legs)
        print(ta,tb,lo,'Ups0:',CV.vstr(a),'| 4lc+C0:',CV.vstr(b), '| lc', lc)
    print('  t',round(time.time()-t0,1))
