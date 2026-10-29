import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np, sympy as sp
import cq_core, cq_gops as GO, cq_gval as GV, cq_families as FM, cq_fsl2 as F2, cq_su11 as S11
rng=np.random.default_rng(4)
ONE={(0,0,0,0):1}
for it in range(2):
    za=FM._rv(16,rng); zb=FM._rv(16,rng)
    q1,q2=S11.NULLS[2],S11.NULLS[5]
    u1=FM._onshell('lamP',q1,rng); u2=FM._onshell('lamM',q2,rng)
    legs=[('lamP',u1,q1),('lamM',u2,q2)]
    op=[(1,[('lamP',za,ONE),('lamM',zb,ONE)])]
    a,b=GO.planar_graded(op,legs,lambda x,y,w: GV.kernel(x,y,w,'D0'))
    val=F2.symv(GV.combine((a,b),1)); O=F2.symc(GO.op_graded(op,legs))
    print('tr(lam+ lam-) L=2: sum C0 =', sp.nsimplify(val), ' -8/3 O =', sp.nsimplify(-sp.Rational(8,3)*O), ' ratio', sp.nsimplify(val/O) if O!=0 else None)
    # L=3 with a Z spectator: tr(lam+ lam- Z)
    zZ=S11.zZ; eZ=S11.eZ; q3=S11.NULLS[7]
    legs3=[('lamP',u1,q1),('lamM',u2,q2),('AW',eZ,q3)]
    op3=[(1,[('lamP',za,ONE),('lamM',zb,ONE),('AW',zZ,ONE)])]
    a,b=GO.planar_graded(op3,legs3,lambda x,y,w: GV.kernel(x,y,w,'D0'))
    val=F2.symv(GV.combine((a,b),1)); O=F2.symc(GO.op_graded(op3,legs3))
    print('tr(lam+ lam- Z) L=3 (three-field part): sum C0 =', sp.nsimplify(val), ' -2/3 N O = -4 O =', sp.nsimplify(-4*O))
