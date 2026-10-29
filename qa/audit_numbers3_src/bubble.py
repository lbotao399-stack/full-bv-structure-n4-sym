import sympy as sp
D=4
p=sp.symbols('p0:4'); q=sp.symbols('q0:4')  # q = p1 loop momentum
a=sp.symbols('a0:4'); b=sp.symbols('b0:4')  # A(p), A(-p) polarisations
def dot(x,y): return sum(x[i]*y[i] for i in range(D))
def e(m): return [1 if i==m else 0 for i in range(D)]
def K(u1,k1,u2,k2,u3,k3):
    # (u3,k2-k1)(u1,u2)+(u1,k3-k2)(u2,u3)+(u2,k1-k3)(u3,u1)
    sub=lambda x,y:[x[i]-y[i] for i in range(D)]
    return dot(u3,sub(k2,k1))*dot(u1,u2)+dot(u1,sub(k3,k2))*dot(u2,u3)+dot(u2,sub(k1,k3))*dot(u3,u1)
mp=[-x for x in p]; mq=[-x for x in q]; qmp=[q[i]-p[i] for i in range(D)]; pmq=[p[i]-q[i] for i in range(D)]
S=0
for m in range(D):
  for n in range(D):
    S+= K(a,p, e(m),mq, e(n),qmp) * K(b,mp, e(m),q, e(n),pmq)
S=sp.expand(S)
p2=dot(p,p); q2=dot(q,q); qmp2=dot(qmp,qmp)
AA=dot(a,b); ipA=dot(p,a)*dot(p,b)
br=dot(q,a)*dot(qmp,b)+dot(qmp,a)*dot(q,b)
N=sp.expand(-(q2+qmp2+4*p2)*AA+2*ipA-5*br)
r=sp.simplify(S/N)
print('ratio K.K / claimed numerator =', r)
# Laplacian in q of N
lap=sum(sp.diff(N,q[i],2) for i in range(D))
print('positive Laplacian of N (-sum d2):', sp.simplify(-lap))
# gaugino numerator check: -16[(q,q-p)AA - (q.a)((q-p).b) - ((q-p).a)(q.b)] Laplacian
G=-16*(dot(q,qmp)*AA - dot(q,a)*dot(qmp,b) - dot(qmp,a)*dot(q,b))
print('gaugino pos Laplacian', sp.simplify(-sum(sp.diff(G,q[i],2) for i in range(D))))
t1,t2=sp.symbols('t1 t2',positive=True); T=t1+t2
cen=[t2/T*x for x in p]
Gc=G.subs(dict(zip(q,cen)))
print('gaugino centre', sp.factor(sp.simplify(Gc)))
Gh=br
print('ghost centre', sp.factor(sp.simplify(Gh.subs(dict(zip(q,cen))))), 'lap', sp.simplify(-sum(sp.diff(Gh,q[i],2) for i in range(D))))
Nc=sp.simplify(N.subs(dict(zip(q,cen))))
print('gauge centre', sp.collect(sp.expand(Nc),[AA,ipA]))
