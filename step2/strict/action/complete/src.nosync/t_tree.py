import numpy as np, sympy as sp, random
from fractions import Fraction as Fr
import cb_common as B, cb_tree as T
E, C, V, Q = B.E, B.C, B.V, B.Q
random.seed(5)
Sfree = Q.dsum((-0.25, Q.F_dAdA_V), (-0.25, Q.F_dAdA_p), (0.5, Q.F_ldl), (1, Q.F_dcAs))
k = [1, -2, 3, 1]
qn = np.array(C.q_matrix(k).tolist(), dtype=complex)
M = T.Mmat(Sfree, k)
print('M_Sfree == q:', np.allclose(M, qn), np.abs(M - qn).max())
xi = np.zeros(16, complex)
for c in range(8): xi[c] = random.randint(-2, 2)
Tx = Q.xi_terms(xi)
mT1 = Q.map_T1(xi)
MT = T.Mmat(Tx['T1'], k)
# compare on lam legs
worst = 0
for _ in range(5):
    leg = Q.rand_leg('lam'); leg['p'] = k
    a = mT1(leg); v = MT @ T.leg_vec(leg)
    worst = max(worst, np.abs(T.leg_vec(a) - v).max())
print('M_T1 vs map_T1 on lam legs:', worst, 'other species image of M_T1 on A*?')
legA = Q.rand_leg('A'); legA['p'] = k
print('  M_T1 on A leg:', np.abs(MT @ T.leg_vec(legA)).max())
# odd G: N = 2(A,A*) - 4(lam*,lam) - (c*,c); (lam*,lam) = (lam,lam*) (lam* even)
N = Q.dsum((2, Q.F_AAs_V), (2, Q.F_AAs_p), (-4, Q.F_lamslam), (-1, Q.F_csc))
QN = Q.apply_Q(N)
bN = T.bracket(N, 1, Sfree)
bS = T.bracket(Sfree, 0, N)
w1 = w2 = 0
for spcs in (['A', 'A'], ['lam', 'lam'], ['c', 'As'], ['As', 'c'], ['A', 'c']):
    for _ in range(3):
        legs = [Q.rand_leg(s) for s in spcs]
        legs[1]['p'] = [-x for x in legs[0]['p']]
        w1 = max(w1, abs(QN(legs) - bS(legs)))
        w2 = max(w2, abs(QN(legs) - bN(legs)))
print('{Sfree,N} - QN:', w1, '   {N,Sfree} - QN:', w2)
leg = Q.rand_leg('lam'); leg['p'] = k
a = T.leg_vec(mT1(leg)); v = MT @ T.leg_vec(leg)
nz = np.nonzero(np.abs(a) + np.abs(v))[0]
print(nz, a[nz], v[nz])
Pi = np.diag([(-1) ** p for p in C.PARITY]).astype(complex)
for nm, MM in (('-K0 H(p,-p)', MT), ('-K0 H(-p,p)', -T.K0 @ T.hessian(Tx['T1'], [-x for x in k], k)),
               ('-K0 H^T', -T.K0 @ T.hessian(Tx['T1'], k, [-x for x in k]).T), ('Pi', Pi @ MT @ Pi), ('MPi', MT @ Pi), ('PiM', Pi @ MT)):
    print(nm, np.abs(MM @ T.leg_vec(leg) - a).max(), np.abs(MM @ T.leg_vec(leg) + a).max())
