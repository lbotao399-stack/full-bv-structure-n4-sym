from fractions import Fraction as Fr
from superalg import Alg
A=Alg()
A.gen(0,0); A.gen(2,1); A.gen(3,1)
H,U,G=20,21,22
for g in (H,U,G): A.gen(g,0)
X=lambda i:{((i,1),):1}
def mono(*g):
    e=A.const(1)
    for x in g: e=A.mul(e,X(x))
    return e
# simplest: I = x0^4 (quartic) ; P00 = 1 ; W_1 quadratic should be 1/2 STr(M) with M = P*d d I = 12 x0^2 -> 6 x0^2
I=A.scal(mono(0,0,0,0),Fr(1))
P={(0,0):Fr(1)}
def dP(F):
    out={}
    for (a,b),v in P.items(): out=A.add(out,A.scal(A.dleft(A.dleft(F,a),b),v/2))
    return out
gt=lambda m: dict(m).get(G,0)<=3
expo=A.mul(X(G),A.mul(X(U),I))
E=A.const(1);t=A.const(1)
for m in range(1,4):
    t=A.scal(A.mul(t,expo,gt),Fr(1,m)); E=A.add(E,t)
Z=dict(E);t=dict(E);n=0
while True:
    n+=1; t=A.scal(A.mul(X(H),dP(t),gt),Fr(1,n))
    if not t: break
    Z=A.add(Z,t)
Y=A.add(Z,A.const(1),-1)
print({m:c for m,c in Y.items() if dict(m).get(G,0)==1})
