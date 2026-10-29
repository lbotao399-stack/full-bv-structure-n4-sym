import numpy as np, ft, loop, time
rng=np.random.default_rng(11)
def rv(n): return rng.normal(size=n)+1j*rng.normal(size=n)
def A_leg(p, comps, mono=()):
    f=np.zeros((10,3),complex)
    for M in comps: f[M]=rv(3)
    return ft.Leg('A',f.reshape(30),p,mono)
p1=rng.normal(size=4); p2=rng.normal(size=4); p3=-p1-p2
legs=[ft.Leg('c',rv(3),p1,(1,)),A_leg(p2,range(4,10)),A_leg(p3,range(4,10))]
S=list(legs); S[0]=ft.Q_leg(legs[0])
ol=loop.OneLoop(S,rng)
it=[x for x in ol.items if x[2]==3][0]
mono,w,n,offs,c,cl,err=it
print([o@o for o in offs])
for eps in [1e-3,1e-5,1e-7]:
    vals=[]
    for nt,ny in [(60,(40,60,40)),(80,(60,80,60)),(120,(90,120,90)),(160,(120,160,120))]:
        qd=loop.Quad(nt=nt,ny=ny); t=time.time()
        vals.append(qd.integrate(3,offs,c,cl,eps))
    print(eps, vals[-1], [abs(v-vals[-1])/abs(vals[-1]) for v in vals[:-1]])
