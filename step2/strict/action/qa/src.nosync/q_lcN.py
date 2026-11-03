"""QA: Q N on two legs with engine B's Q (left rule) for N = 2 int(A,A*) - 4 int(lam*,lam) - int(c*,c) (Thm c4-log),
(lam*,lam) = -(lam,lam*); compare with engine A's plain-cycle value lc^(0)|_2 = u[-1, 4, -1] on
[int(dA,dA), int(lam,dslash lam), int(dc,A*)], i.e. lc^(0) = u Lambda N (Prop at-log)."""
import sys, os, random
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', 'engine_b', 'src.nosync')))
import eb_qb as Q
random.seed(11)
# (lam*,lam) with lam* even (Def c1-pairing) equals (lam,lam*): N = 2 N_A - 4 int(lam,lam*) - int(c*,c)
N = Q.dsum((2, Q.F_AAs_V), (2, Q.F_AAs_p), (-4, Q.F_lamslam), (-1, Q.F_csc))
QN = Q.apply_Q(N)
print('convention (lam*,lam) = +(lam,lam*):')
for name, spcs, basis in [('A (V, perp)', ['A', 'A'], [Q.F_dAdA_V, Q.F_dAdA_p]),
                          ('lam', ['lam', 'lam'], [Q.F_ldl]), ('c,A*', ['c', 'As'], [Q.F_dcAs])]:
    print(name, Q.fit_on(QN, basis, [spcs], ntrial=4))
# individual images (sign conventions used in compare.json)
for nm, f, spcs, basis in [('Q int(A_V,A*_V) on (c,A*)', Q.F_AAs_V, ['c', 'As'], [Q.F_dcAs]),
                           ('Q int(A_V,A*_V) on (A,A)', Q.F_AAs_V, ['A', 'A'], [Q.F_dAdA_V, Q.F_dAdA_p]),
                           ('Q int(c*,c) on (c,A*)', Q.F_csc, ['c', 'As'], [Q.F_dcAs]),
                           ('Q int(lam,lam*) on (lam,lam)', Q.F_lamslam, ['lam', 'lam'], [Q.F_ldl])]:
    print(nm, Q.fit_on(Q.apply_Q(f), basis, [spcs], ntrial=4))
