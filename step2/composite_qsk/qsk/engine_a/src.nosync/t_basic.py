import time, numpy as np
t0=time.time()
import qa_s1 as S1, qa_core as QC
print('spin', S1.check_spin())
xi = S1.chir_basis('lamP')[0]
print('xi', xi)
for w in ('T1','T2','T3','T4'):
    b = S1.quad_block(S1.odd_terms(xi,(w,)))
    nz = np.argwhere(np.abs(b.h0)+sum(np.abs(x) for x in b.Hj)+sum(np.abs(x) for x in b.Hl) > 0)
    print(w, len(nz), sorted(set((QC.Q.E.TYPES[[i for i,t in enumerate(QC.Q.E.TYPES) if QC.Q.OFS[t]<=r][-1]], QC.Q.E.TYPES[[i for i,t in enumerate(QC.Q.E.TYPES) if QC.Q.OFS[t]<=c][-1]]) for r,c in nz)))
print('time', time.time()-t0)
