import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, sympy as sp, warnings; warnings.filterwarnings('ignore'); np.seterr(all='ignore')
import cq_core, cq_gops as GO, cq_gval as GV, cq_families as FM, cq_fsl2 as F2, cq_su11 as S11, ea_core as E
ONE={(0,0,0,0):1}
planes=[(0,1),(2,3),(4,5),(6,7),(8,9)]
Sg=[np.asarray(E.SG[M])@np.asarray(E.SGB[N])/2 for M,N in planes]
def gauss(v):
    v=v/ v[np.argmax(np.abs(v))]
    for sc in [1,2,4,8]:
        w=v*sc
        if np.allclose(w.real,np.round(w.real)) and np.allclose(w.imag,np.round(w.imag)): return np.round(w.real)+1j*np.round(w.imag)
st={1:[],-1:[]}
for signs in itertools.product([1,-1],repeat=5):
    A=np.vstack([Sg[i]-signs[i]*0.5j*np.eye(16) for i in range(5)])
    u,s,vh=np.linalg.svd(A)
    if s[-1]<1e-9:
        v=gauss(vh[-1].conj()); chi=int(np.round(((np.asarray(E.G5S)@v)@v.conj()/(v@v.conj())).real)); st[chi].append((signs,v))
rng=np.random.default_rng(1)
q1,q2=S11.NULLS[2],S11.NULLS[5]
u1=FM._onshell('lamP',q1,rng); u2=FM._onshell('lamM',q2,rng)
legs=[('lamP',u1,q1),('lamM',u2,q2)]
bad=0; n=0; res=[]
for (s1,a),(s2,b) in itertools.product(st[1],st[-1]):
    op=[(1,[('lamP',a,ONE),('lamM',b,ONE)])]
    O=F2.symc(GO.op_graded(op,legs))
    if O==0: continue
    x,y=GO.planar_graded(op,legs,lambda p,q,w: GV.kernel(p,q,w,'D0')); c0=sp.nsimplify(sp.simplify(F2.symv(GV.combine((x,y),1))/O))
    n+=1; bad+= (c0!=sp.Rational(-8,3)); res.append((s1,s2,c0))
print('combos',n,'not -8/3:',bad)
for r in res:
    if r[2]!=sp.Rational(-8,3): print(r)
# random combination
za=FM._rv(16,rng); zb=FM._rv(16,rng)
op=[(1,[('lamP',za,ONE),('lamM',zb,ONE)])]
O=F2.symc(GO.op_graded(op,legs)); x,y=GO.planar_graded(op,legs,lambda p,q,w: GV.kernel(p,q,w,'D0'))
print('random letters C0/O', sp.nsimplify(sp.simplify(F2.symv(GV.combine((x,y),1))/O)))
print('--- invisible combos (O(legs)=0):')
cnt=0
for (s1,a),(s2,b) in itertools.product(st[1],st[-1]):
    op=[(1,[('lamP',a,ONE),('lamM',b,ONE)])]
    O=F2.symc(GO.op_graded(op,legs))
    if O!=0: continue
    x,y=GO.planar_graded(op,legs,lambda p,q,w: GV.kernel(p,q,w,'D0')); c0=sp.nsimplify(F2.symv(GV.combine((x,y),1)))
    if c0!=0: cnt+=1; print(s1,s2,'C0 on legs =',c0, ' overlaps a.u1', np.dot(a,u1), 'b.u2', np.dot(b,u2))
print('nonzero images of invisible operators:',cnt)
print('--- reversed-order operators tr(lam-_b lam+_a) vs -tr(lam+_a lam-_b):')
k=0
for (s1,a),(s2,b) in itertools.product(st[1],st[-1]):
    op=[(1,[('lamP',a,ONE),('lamM',b,ONE)])]; opr=[(1,[('lamM',b,ONE),('lamP',a,ONE)])]
    O=F2.symc(GO.op_graded(op,legs)); Or=F2.symc(GO.op_graded(opr,legs))
    if O==0: continue
    x,y=GO.planar_graded(opr,legs,lambda p,q,w: GV.kernel(p,q,w,'D0')); cr=sp.nsimplify(F2.symv(GV.combine((x,y),1)))
    lr=sp.nsimplify(F2.symv(GO.planar_graded(opr,legs,lambda p,q,w: GV.logkernel(p,q,w))))
    lf=sp.nsimplify(F2.symv(GO.planar_graded(op,legs,lambda p,q,w: GV.logkernel(p,q,w))))
    print(s1,s2,'O',O,'O_rev',Or,' C0(rev)/O_rev', sp.nsimplify(cr/Or), ' log(rev)/log(fwd)', sp.nsimplify(lr/lf) if lf!=0 else None)
    k+=1
    if k>=4: break
print('--- same-chirality operators on mixed legs:')
for ch in (1,-1):
    nz=0
    for (s1,a),(s2,b) in itertools.product(st[ch],st[ch]):
        op=[(1,[('lamP',a,ONE),('lamM',b,ONE)])]
        O=F2.symc(GO.op_graded(op,legs))
        x,y=GO.planar_graded(op,legs,lambda p,q,w: GV.kernel(p,q,w,'D0')); c0=sp.nsimplify(F2.symv(GV.combine((x,y),1)))
        if c0!=0 or O!=0: nz+=1; print(ch,s1,s2,'O',O,'C0',c0) if nz<4 else None
    print('chirality',ch,'nonzero',nz)
