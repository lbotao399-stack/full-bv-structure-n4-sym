"""Cheap recomputation: twist-two two-site spectrum from Beisert H'_{12} (B03c 3.14), su(2) and su(1|1) cyclic spectra,
harmonic-number values, invariant combinations."""
import itertools, json
from fractions import Fraction as F
import numpy as np, sympy as sp
def h(n): return sum(F(1,j) for j in range(1,n+1))
out={}
# twist two: basis D^k Z D^{S-k} Z/(k!(S-k)!), H'_{12} matrix M[k,k'] = delta (h(k)+h(S-k)) - (1-delta)/|k-k'|
tw={}
for S in range(0,9):
    n=S+1
    M=sp.zeros(n,n)
    for k in range(n):
        for kp in range(n):
            M[k,kp]= sp.Rational(h(k)+h(S-k)) if k==kp else -sp.Rational(1,abs(k-kp))
    # H21: exchange sites k -> S-k ; P[k,k'] = delta(k', S-k)
    P=sp.zeros(n,n)
    for k in range(n): P[k,S-k]=1
    H=M+P*M*P
    # cyclic (symmetric) states: P v = v
    ev=H.eigenvals()
    # restrict to symmetric subspace
    B=sp.Matrix([[ (1 if (j==k or j==S-k) else 0) for j in range(n)] for k in range(S//2+1)]).T
    Hs=(B.T*B).inv()*B.T*H*B
    tw[S]=sorted([str(-2*e) for e in Hs.eigenvals()])  # -2 units (H12+H21) = map
    tw[S]=sorted(set(str(sp.nsimplify(-2*e)) for e in Hs.eigenvals()))
out['twist_two_map_eigs(-2(H12+H21))']=tw
out['minus8h']={S:str(-8*h(S)) for S in (0,2,4,6,8)}
out['D1=-32h']={S:str(-32*h(S)) for S in range(0,9)}
# su(2) cyclic spectra of sum(1-P), L<=8 small check
def cyc_spec(L,M,graded=False):
    states=[s for s in itertools.product([0,1],repeat=L) if sum(s)==M]
    idx={s:i for i,s in enumerate(states)}
    n=len(states); H=np.zeros((n,n))
    for s in states:
        for l in range(L):
            a,b=s[l],s[(l+1)%L]
            H[idx[s],idx[s]]+=1
            t=list(s); t[l],t[(l+1)%L]=b,a; t=tuple(t)
            sign=1
            if graded:
                if a==1 and b==1: sign=-1
                # moving across the boundary: fermion passing the others in cyclic shift handled below
                if l==L-1 and a!=b:
                    # exchange of last and first: operator ordering crosses the other M-1... fermions
                    nf=M-1 if (a==1 or b==1) else 0
                    sign*= (-1)**(nf) if (a==1 or b==1) else 1
                if l==L-1 and a==1 and b==1:
                    sign = -1*(-1)**(M-2)*(-1)  # psi psi exchange across boundary
            H[idx[t],idx[s]]-=sign
    # cyclic shift operator with grading
    C=np.zeros((n,n))
    for s in states:
        t=s[1:]+s[:1]
        sign=(-1)**((M-1)) if (graded and s[0]==1) else 1
        C[idx[t],idx[s]]+=sign
    # project onto C=1
    w,v=np.linalg.eig(C)
    Vs=v[:,np.isclose(w,1)]
    if Vs.shape[1]==0: return []
    Q,_=np.linalg.qr(Vs)
    Hp=Q.conj().T@H@Q
    return sorted(np.round(np.linalg.eigvals(Hp).real,6).tolist())
out['su2']={f'{L},{M}':cyc_spec(L,M) for L,M in [(4,2),(5,2),(6,2),(6,3),(7,2),(7,3),(8,2)]}
L=F(0)
out['invariants']={'27/4-1/3':str(F(27,4)-F(1,3)),'2*79/12+1/3-27/4':str(2*F(79,12)+F(1,3)-F(27,4)),
  'Konishi -48u': '48 g^2N/(32pi^2)=3g^2N/(2pi^2)', 'D1 traceless at t=-1/3': str(-4*F(-1,3)-F(4,3)),'trace':str(-4*F(-1,3)-F(148,3)),
  'two-eval first':str(F(-4,3)+F(8,3)), 'two-eval second O':str(F(8,3)-F(4,3)),'second S':str(-2+F(5,3))}
json.dump(out,open(__file__.replace('src.nosync/spectra.py','spectra_check.json'),'w'),indent=1)
print(json.dumps(out,indent=1))
