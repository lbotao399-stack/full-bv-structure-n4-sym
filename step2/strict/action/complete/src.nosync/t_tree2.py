import numpy as np, random
import cb_common as B, cb_tree as T
Q, C = B.Q, B.C
random.seed(5)
Sfree = Q.dsum((-0.25, Q.F_dAdA_V), (-0.25, Q.F_dAdA_p), (0.5, Q.F_ldl), (1, Q.F_dcAs))
k = [1, -2, 3, 1]
print('M_Sfree == q:', np.allclose(T.Mmat(Sfree, k), np.array(C.q_matrix(k).tolist(), dtype=complex)))
xi = np.zeros(16, complex)
for c in range(8): xi[c] = random.randint(-2, 2)
Tx = Q.xi_terms(xi); mT1 = Q.map_T1(xi); MT = T.Mmat(Tx['T1'], k)
w = 0
for _ in range(4):
    leg = Q.rand_leg('lam'); leg['p'] = k
    w = max(w, np.abs(T.leg_vec(mT1(leg)) - MT @ T.leg_vec(leg)).max())
print('M_T1 vs map_T1:', w)
N = Q.dsum((2, Q.F_AAs_V), (2, Q.F_AAs_p), (-4, Q.F_lamslam), (-1, Q.F_csc))
QN = Q.apply_Q(N)
for gp in (0, 1):
    bN = T.bracket(N, gp, Sfree)
    ww = 0
    for spcs in (['A', 'A'], ['lam', 'lam'], ['c', 'As'], ['As', 'c'], ['A', 'c'], ['c', 'A']):
        for _ in range(3):
            legs = [Q.rand_leg(s) for s in spcs]; legs[1]['p'] = [-x for x in legs[0]['p']]
            ww = max(ww, abs(QN(legs) - bN(legs)), abs(QN(legs) + bN(legs)) * 0)
    print('gpar', gp, '{N,Sfree} - QN:', ww)
# T1 derivation on a test F vs engine B apply_vf
F = Q.dsum((1, Q.F_dsdc_dsA), (0.3, Q.F_dsA_dsA))
b1 = T.bracket(Tx['T1'], 0, F); b2 = Q.apply_vf(F, mT1)
ww = 0
for spcs in (['c', 'lam'], ['lam', 'c'], ['A', 'lam'], ['lam', 'A']):
    for _ in range(3):
        legs = [Q.rand_leg(s) for s in spcs]; legs[1]['p'] = [-x for x in legs[0]['p']]
        ww = max(ww, abs(b1(legs) - b2(legs)))
print('{T1,F} vs apply_vf:', ww)
