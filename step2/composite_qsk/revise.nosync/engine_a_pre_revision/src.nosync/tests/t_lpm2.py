import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, sympy as sp
import cq_core, cq_gops as GO, cq_gval as GV, cq_families as FM, cq_fsl2 as F2, cq_su11 as S11, ea_core as E
ONE={(0,0,0,0):1}
planes=[(0,1),(2,3),(4,5),(6,7),(8,9)]
Sg=[np.asarray(E.SG[M])@np.asarray(E.SGB[N])/2 for M,N in planes]
# simultaneous eigenvectors in Sigma
states=[]
for signs in itertools.product([1,-1],repeat=5):
    A=np.vstack([Sg[i]-signs[i]*0.5j*np.eye(16) for i in range(5)])
    u,s,vh=np.linalg.svd(A); null=vh[s.size-1:] if s[-1]<1e-9 else None
    if s[-1]<1e-9:
        v=vh[-1].conj(); v=v/v[np.argmax(np.abs(v))]
        chi=np.round((np.asarray(E.G5S)@v)[np.argmax(np.abs(v))]/v[np.argmax(np.abs(v))]).real
        states.append((signs,v,int(chi)))
print('states',len(states))
psiw=None
for sg,v,c in states:
    if np.allclose(v/ v[np.argmax(np.abs(v))], S11.PSI/S11.PSI[np.argmax(np.abs(S11.PSI))]): psiw=sg
psiw=(1,-1,-1,-1,-1); print("psi weights",psiw)
def gauss(v):
    v=v/ v[np.argmax(np.abs(v))]
    for sc in [1,2,4,8]:
        w=v*sc
        if np.allclose(w.real,np.round(w.real)) and np.allclose(w.imag,np.round(w.imag)): return np.round(w.real)+1j*np.round(w.imag)
    raise ValueError
rng=np.random.default_rng(1)
q1,q2=S11.NULLS[2],S11.NULLS[5]
for sg,v,c in states:
    if c!=-1: continue
    lm=gauss(v)
    tot=[psiw[i]+sg[i] for i in range(5)]
    if tot[2:]==[0,0,0]: tag='R-singlet-type'
    else: tag='R-15'
    # need nonzero overlap with on-shell legs
    try:
        u1=S11.onshell_u(q1); 
        Bm=np.column_stack(E.BASIS['lamM']).astype(complex)
        u2=FM._onshell('lamM',q2,rng)
    except Exception as e: print('skip',e); continue
    if abs(np.dot(lm,u2))<1e-9 or abs(np.dot(S11.PSI,u1))<1e-9: print(sg,tag,'zero overlap'); continue
    legs=[('lamP',u1,q1),('lamM',u2,q2)]
    op=[(1,[('lamP',S11.PSI,ONE),('lamM',lm,ONE)])]
    O=F2.symc(GO.op_graded(op,legs))
    lg=GO.planar_graded(op,legs,lambda x,y,w: GV.logkernel(x,y,w)); lcn=F2.symv(lg)-8*O
    a,b=GO.planar_graded(op,legs,lambda x,y,w: GV.kernel(x,y,w,'D0')); c0=F2.symv(GV.combine((a,b),1))
    print(sg, tag, 'spin tot', tot[:2], 'lc^nat/O =', sp.nsimplify(sp.simplify(lcn/O)), ' C0/O =', sp.nsimplify(sp.simplify(c0/O)), '(P7c: -8/3)')
