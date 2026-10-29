import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import json, cq_twist as TW
for N in [1, 2]:
    r = TW.run(N, kinds=('Pg',))
    m=r['Pg']['matrix_by_logbasis']['1']
    print(N, 'Pg rows (b1 a=0..N, b2 a=0..N, q1q2.., q1^2.., q2^2..):', m, 'vv', r['Pg']['vertex_vertex_part_nonzero'], 'res', r['Pg']['fit_residual'], r['Pg']['seconds'])
t0=time.time(); r=TW.run(3, kinds=('D0','Pn')); print(3, r['D0']['matrix_by_logbasis']['1'], r['Pn']['matrix_by_logbasis']['1'], r['D0']['vertex_vertex_part_nonzero'], time.time()-t0)
