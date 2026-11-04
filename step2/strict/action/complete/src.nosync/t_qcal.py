import numpy as np, sympy as sp
import cb_common as B
E, C, V = B.E, B.C, B.V
def Hnum(terms, kr, kc):
    bl = V.hessian_blocks(terms, [], list(map(sp.Integer, kr)), list(map(sp.Integer, kc)), E.PAIRS)
    H = np.zeros((54, 54), complex)
    for (sJ, sK), d in bl.items():
        for key, M in d.items():
            H[E.OFF[sJ]:E.OFF[sJ] + V.DIM[sJ], E.OFF[sK]:E.OFF[sK] + V.DIM[sK]] += np.array(M.tolist(), dtype=complex)
    return H
k = [1, -2, 3, 1]; mk = [-x for x in k]
q = np.array(C.q_matrix(k).tolist(), dtype=complex)
K0 = np.array(C.k0_matrix().tolist(), dtype=complex)
Pi = np.diag([(-1) ** p for p in C.PARITY]).astype(complex)
for nm, Hm in (('H(k,-k)', Hnum(V.FREE, k, mk)), ('H(-k,k)', Hnum(V.FREE, mk, k))):
    for a, X in (('', np.eye(54)), ('Pi', Pi)):
        for b, Yy in (('', np.eye(54)), ('Pi', Pi)):
            for tr in (False, True):
                Hh = Hm.T if tr else Hm
                for order in ('K0 H', 'H K0'):
                    M = X @ (K0 @ Hh if order == 'K0 H' else Hh @ K0) @ Yy
                    for s in (1, -1, 0.5, -0.5, 2, -2):
                        if np.allclose(s * M, q):
                            print('match', nm, a, b, 'T' if tr else '', order, s)
