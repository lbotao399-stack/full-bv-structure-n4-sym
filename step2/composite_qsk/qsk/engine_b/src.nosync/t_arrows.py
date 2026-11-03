import qb_core as B, sympy as sp
import qb_e0 as E
k=E.k
xi=[0]*16; xi[0]=1; xi[5]=sp.Rational(2,3)
quad,cub=B.register_xi_terms(xi, tag='tt')
Xs=E.op_from_functional(quad)
Xf=B.xi_op_fields(k, xi)
OFF=B.OFF
for (t,s) in [('A','lam'),('lam','A')]:
    r=set()
    for i in range(B.V.DIM[t]):
        for j in range(B.V.DIM[s]):
            a=Xs[OFF[t]+i,OFF[s]+j]; b=Xf[OFF[t]+i,OFF[s]+j]
            if a!=0 or b!=0:
                r.add(sp.simplify(a/b) if b!=0 else ('b=0',a))
    print(t,'<-',s, r)
