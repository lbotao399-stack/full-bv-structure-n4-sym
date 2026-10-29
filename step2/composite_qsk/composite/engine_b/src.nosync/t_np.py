import sys, sympy as sp
sys.path.insert(0,'.')
import ce_gluon as GL, ce_su11 as SU
n, nbar = SU.n, SU.nbar
for (A,B) in [('Np','Np'), ('N','N')]:
    job = GL.CJobP('AV','AV',['AV','AV'], pol=[GL.EPS]*2, za=GL.LZ[A], zb=GL.LZ[B], deriv=('null',2,n,nbar), nullw=False, pre={'a':GL.PRE[A],'b':GL.PRE[B]})
    res = job.kernels(mlist=(0,), want=('Pb',), verbose=False)
    tot = {}
    for w, dm in res['Pb'].items():
        for b_, P in dm.get(0,{}).items():
            for mon, cf in P.items():
                key = (mon[job.names.index('ss')], mon[job.names.index('su')])
                tot[key] = tot.get(key,0)+1
    print(A,B,'Pb raw monomial counts by letter orders', tot)
