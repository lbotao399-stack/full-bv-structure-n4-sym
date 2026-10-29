import numpy as np, ft, loop
rng=np.random.default_rng(11)
def rv(n): return rng.normal(size=n)+1j*rng.normal(size=n)
def A_leg(p, comps, mono=()):
    f=np.zeros((10,3),complex)
    for M in comps: f[M]=rv(3)
    return ft.Leg('A',f.reshape(30),p,mono)
p1=rng.normal(size=4); p2=rng.normal(size=4); p3=-p1-p2
legs=[ft.Leg('c',rv(3),p1,(1,)),A_leg(p2,range(4,10)),A_leg(p3,range(4,10))]
for j in range(3):
    S=list(legs); S[j]=ft.Q_leg(legs[j])
    tl=loop.ThreeLeg(S,rng)
    print('legset',j+1)
    for it in tl.ol3.items:
        print('  1PI n=%d mono=%s err=%.2e maxc=%.2e'%(it[2],it[0],it[6],np.abs(it[4]).max()))
    for q,Ll,ols in tl.pairs:
        for ol,ch in ols:
            for it in ol.items:
                print('  chi slot=%s n=%d mono=%s err=%.2e maxc=%.2e'%(ch.slot,it[2],it[0],it[6],np.abs(it[4]).max()))
