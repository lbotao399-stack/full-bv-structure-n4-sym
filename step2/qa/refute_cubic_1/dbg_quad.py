import numpy as np, ft, loop, mpmath as mp
mp.mp.dps=30
rng=np.random.default_rng(11)
def rv(n): return rng.normal(size=n)+1j*rng.normal(size=n)
def A_leg(p, comps, mono=()):
    f=np.zeros((10,3),complex)
    for M in comps: f[M]=rv(3)
    return ft.Leg('A',f.reshape(30),p,mono)
p1=rng.normal(size=4); p2=rng.normal(size=4); p3=-p1-p2
legs=[ft.Leg('c',rv(3),p1,(1,)),A_leg(p2,range(4,10)),A_leg(p3,range(4,10))]
S=list(legs); S[0]=ft.Q_leg(legs[0])
tl=loop.ThreeLeg(S,rng)
q,Ll,ols=tl.pairs[0]
ol,ch=ols[0]
it=[x for x in ol.items if x[2]==2][0]
mono,w,n,offs,c,cl,err=it
quad=loop.Quad()
def mpval(eps):
    s1=offs[0]; q2=float(s1@s1)
    def f(x):
        u=min(x,1-x)
        if u==0: return 0
        kc=-(float(x)*s1)[None,:]
        Nc=complex(loop.polyval(c,kc)[0]); N2=0.25*complex(loop.polyval(cl,kc)[0])
        D=u*(1-u)*q2; T0=mp.mpf(eps)/u; z=T0*D
        return Nc*mp.e1(z)+N2*(mp.exp(-z)/T0-D*mp.e1(z))
    return complex(mp.quad(f,[0,0.25,0.5,0.75,1]))*loop.INV16PI2
for eps in [2e-3,1e-4,1.5625e-5]:
    a=quad.integrate(2,offs,c,cl,eps); b=mpval(eps)
    print(eps, a, b, abs(a-b)/abs(b))
for nt in (80,):
  for ny in [(60,80,60),(120,160,120)]:
    qd=loop.Quad(nt=nt,ny=ny)
    print(ny, qd.integrate(2,offs,c,cl,1e-4))
