import qb_core as B, sympy as sp
C=B.C; V=B.V
k=sp.symbols('k0:4'); mk=[-x for x in k]
INT6=C.INT6
pairs=[(a,b) for a in INT6 for b in INT6]
def H(terms, kr, kc):
    bl=V.hessian_blocks(terms, [], list(kr), list(kc), pairs)
    M=sp.zeros(C.NT,C.NT)
    for (a,b),d in bl.items():
        for key,Mx in d.items():
            for i in range(Mx.rows):
                for j in range(Mx.cols):
                    if Mx[i,j]!=0: M[C.OFF[a]+i,C.OFF[b]+j]+=Mx[i,j]
    return M.applyfunc(sp.expand)
Hf=H(['kinA','kinlam','kingh'], mk, k)
q=C.q_matrix(k); K0=C.k0_matrix(); PI=C.PI
cands={'K0 PI H':K0*PI*Hf,'K0 H':K0*Hf,'PI K0 H':PI*K0*Hf,'H K0':Hf*K0,'H PI K0':Hf*PI*K0,'K0 H^T':K0*Hf.T, 'K0 PI H^T':K0*PI*Hf.T}
for nm,X in cands.items():
    X=X.applyfunc(sp.expand)
    for c in (1,-1,2,-2,sp.Rational(1,2),-sp.Rational(1,2)):
        if (X-c*q).applyfunc(sp.expand).is_zero_matrix: print('MATCH',nm,c)
# print a few entries
for (a,b) in [('As','A'),('A','c'),('lams','lam'),('cs','As')]:
    i=C.OFF[a]; j=C.OFF[b]
    print(a,b,'q',q[i,j],' H[i,j]',Hf[i,j],' H[j,i]',Hf[j,i])
