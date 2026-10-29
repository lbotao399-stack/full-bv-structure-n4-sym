"""Check the multilinear vertex evaluation (from the action densities) against the printed vertex
functions eq (c1-AAA), (c1-AAAA) and the table of def c1-vertex."""
import numpy as np, ft
rng=np.random.default_rng(3)
def rv(*s): return rng.normal(size=s)+1j*rng.normal(size=s)
p1,p2=rng.normal(size=4),rng.normal(size=4); p3=-p1-p2
def pad(p): return np.concatenate([p,np.zeros(6)])
A=[rv(10,3) for _ in range(3)]
legs=[ft.Leg('A',A[i].reshape(30),[p1,p2,p3][i],()) for i in range(3)]
mine=ft.multilinear([ft.CUBIC[0]],legs)[()]
br=lambda X,Y: np.einsum('abc,a,b->c',ft.fabc,X,Y)
def iota(p,Ai): return pad(p)@Ai
def pair(X,Y): return X@Y
P=[p1,p2,p3]
def term(i,j,k): # (iota_{p_j - p_i} A_k, [A_i, A_j])
    return pair(iota(P[j]-P[i],A[k]), br(np.einsum('Ma,Mb->ab',A[i],A[j]).sum()*0+0 if False else None, None)) if False else None
def brAA(Ai,Aj): # [A_i,A_j] contracting V10 with metric
    return np.einsum('abc,Ma,Mb->c',ft.fabc,Ai,Aj)
printed=-0.5j*( iota(p2-p1,A[2])@brAA(A[0],A[1]) + iota(p3-p2,A[0])@brAA(A[1],A[2]) + iota(p1-p3,A[1])@brAA(A[2],A[0]) )
print('AAA mine vs printed:', mine, printed, abs(mine-printed)/abs(printed))
A4=[rv(10,3) for _ in range(4)]
legs4=[ft.Leg('A',A4[i].reshape(30),np.zeros(4),()) for i in range(4)]
mine4=ft.multilinear(ft.QUARTIC,legs4)[()]
def wedge(Ai,Aj): # [A_i ^ A_j]_{MN} = [Ai_M,Aj_N]-[Ai_N,Aj_M]
    X=np.einsum('abc,Ma,Nb->MNc',ft.fabc,Ai,Aj); return X-X.transpose(1,0,2)
def hodge(F,G): return 0.5*np.einsum('MNc,MNc->',F,G)
pr4=-0.5*(hodge(wedge(A4[0],A4[1]),wedge(A4[2],A4[3]))+hodge(wedge(A4[0],A4[2]),wedge(A4[1],A4[3]))+hodge(wedge(A4[0],A4[3]),wedge(A4[1],A4[2])))
print('AAAA mine vs printed:', mine4, pr4, abs(mine4-pr4)/abs(pr4))
# lambda lambda A : r1 r2 (lambda(p1), rho([A(p3), lambda(p2)]))
l1,l2=rv(16,3),rv(16,3); A3=rv(10,3)
legs=[ft.Leg('l',l1.reshape(48),p1,(1,)),ft.Leg('l',l2.reshape(48),p2,(2,)),ft.Leg('A',A3.reshape(30),p3,())]
mine=ft.multilinear([ft.CUBIC[1]],legs)[(1,2)]
brAl=np.einsum('abc,Ma,sb->Msc',ft.fabc,A3,l2)   # [A_M, lambda_s]
rho=np.einsum('Mst,Mtc->sc',ft.gam,brAl)
printed=np.einsum('sc,sc->',l1,rho)
print('llA mine vs printed:', mine, printed, abs(mine-printed)/abs(printed))
