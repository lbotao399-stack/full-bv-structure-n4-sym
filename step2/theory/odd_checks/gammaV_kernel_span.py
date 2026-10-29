# Is ker(Gamma|_V : Sym^2 Sigma^* -> V) spanned by symmetric products eta1 eta2 with Gamma(eta1,eta2)|_V = 0 ?
# (odd_components_design.tex, Prop prop:o-Kconsistent). Gamma(eta1,eta2)^M = (rho(e_M) eta1, eta2), rho: Sigma^* -> Sigma = sigbar[M].
import os, sys
import numpy as np
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/superconformal')
import model as Mo
np.seterr(all='ignore')
sigb = [np.array([[complex(*map(float, Mo.SIGB[M][a][b])) for b in range(16)] for a in range(16)]) for M in range(10)]
rng = np.random.default_rng(3)
def GV(e1, e2): return np.array([sigb[M].dot(e1).dot(e2) for M in range(4)])
# symmetric-product coordinates on Sym^2 C^16 (dimension 136)
iu = np.triu_indices(16)
def sym(e1, e2):
    T = np.outer(e1, e2) + np.outer(e2, e1)
    return T[iu]
rows = []
for _ in range(400):
    e1 = rng.normal(size=16) + 1j*rng.normal(size=16)
    A = np.array([sigb[M].T.dot(e1) for M in range(4)])   # Gamma_V(e1, e2) = A e2  (A_{M b} = sum_a e1_a sigb[M]_{a b})
    A = np.array([e1.dot(sigb[M]) for M in range(4)])
    # null space of A (12-dim)
    u, s, vh = np.linalg.svd(A)
    N = vh[4:].conj().T
    e2 = N.dot(rng.normal(size=12) + 1j*rng.normal(size=12))
    assert np.abs(GV(e1, e2)).max() < 1e-9
    rows.append(sym(e1, e2))
R = np.array(rows)
rk = np.linalg.matrix_rank(R, tol=1e-7)
# rank of Gamma_V as a map on Sym^2: build matrix
Gmat = np.array([[ (sigb[M][a, b] + sigb[M][b, a]) if a != b else sigb[M][a, a] for (a, b) in zip(*iu)] for M in range(4)])
print('dim Sym^2 Sigma^* = 136; rank of Gamma_V on Sym^2 =', np.linalg.matrix_rank(Gmat, tol=1e-9))
print('rank of span{eta1 eta2 : Gamma_V(eta1,eta2)=0} =', rk, '(expected 136 - rank Gamma_V)')
