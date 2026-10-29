"""Intertwiner between engine A's and engine B's realisations of Sigma.
A (cliff.py): rho(e_M)|_Sigma = sb[M], pairing u^T w;  B (realization.py): rho(e_M)|_Sigma = R[M], pairing s^T Mp s*.
Find S (u_B = S u_A on Sigma) with  S^T (Mp R_B[M]) S = sb_A[M]  for M = 0..9  (the invariant symmetric forms
(u1, rho(e_M) u2) agree), the vector index M being the same in both engines.  Then on Sigma*: w_B = Sst w_A with
Sst = (S^T Mp)^{-1} (pairing preserved).  Output: S, Sst, the chirality subspaces, written to spinor_map.json."""
import os, sys, json
import sympy as sp
import numpy as np

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'step2', 'twoletter', 'engine_a', 'src.nosync', 'vendor'))
sys.path.insert(0, os.path.join(ROOT, 'step2', 'composite_qsk', 'composite', 'engine_b', 'src.nosync', 'vendor'))
import cliff
import realization as RB

I = sp.I
def tosym(M):
    return sp.Matrix(M.shape[0], M.shape[1], lambda i, j: sp.nsimplify(complex(M[i, j]).real) + I * sp.nsimplify(complex(M[i, j]).imag))
sbA = [tosym(m) for m in cliff.sb]
sA = [tosym(m) for m in cliff.s]
Mp = RB.pairing()
betaB = [Mp * RB.R[M] for M in range(10)]
for M in range(10):
    assert betaB[M].T == betaB[M] and sbA[M].T == sbA[M]
# linear intertwiner: T_M = beta_0^{-1} beta_M (endomorphisms of Sigma); S T^A_M = T^B_M S
TA = [sbA[0].inv() * sbA[M] for M in range(10)]
TB = [betaB[0].inv() * betaB[M] for M in range(10)]
TAn = [np.array(T.evalf(), dtype=complex) for T in TA]; TBn = [np.array(T.evalf(), dtype=complex) for T in TB]
I16 = np.eye(16)
L = np.vstack([np.kron(I16, TAn[M].T) - np.kron(TBn[M], I16) for M in range(1, 10)])   # row-major vec(S)
u, sv, vh = np.linalg.svd(L)
null = vh[np.abs(sv) < 1e-9] if (np.abs(sv) < 1e-9).any() else vh[-1:][np.abs(sv[-1:]) < 1e-9]
print('nullspace dim', int((sv < 1e-9).sum()) + (256 - len(sv)))
v = vh[-1].conj()
Sn = v.reshape(16, 16)
k = np.argmax(np.abs(Sn)); Sn = Sn / Sn.flat[k]
Ssol = sp.Matrix(16, 16, lambda i, j: sp.nsimplify(round(Sn[i, j].real, 10), rational=True) + I * sp.nsimplify(round(Sn[i, j].imag, 10), rational=True))
assert all(Ssol * TA[M] - TB[M] * Ssol == sp.zeros(16, 16) for M in range(1, 10)), 'rationalised S is not exact'
c2 = None
lhs = Ssol.T * betaB[0] * Ssol
for i in range(16):
    for j in range(16):
        if sbA[0][i, j] != 0:
            c2 = sp.nsimplify(sp.simplify(sbA[0][i, j] / lhs[i, j])); break
    if c2 is not None:
        break
Sfin = sp.sqrt(c2) * Ssol
Sfin = sp.simplify(Sfin)
ok = all(sp.simplify(Sfin.T * betaB[M] * Sfin - sbA[M]) == sp.zeros(16, 16) for M in range(10))
print('S^T beta_B S = beta_A for all M:', ok, ' c^2 =', c2)
Sst = sp.simplify((Sfin.T * Mp).inv())
# check rho on Sigma*: Rb_B[M] Sst = Sfin s_A[M]  (rho(e_M): Sigma* -> Sigma)
ok2 = all(sp.simplify(RB.Rb[M] * Sst - Sfin * sA[M]) == sp.zeros(16, 16) for M in range(10))
print('rho on Sigma* intertwined:', ok2)
# chirality: A's G5S = s0 sb1 s2 sb3 on Sigma; B's S_+ (x) W = components 0..7 of Sigma
G5A = sA[0] * sbA[1] * sA[2] * sbA[3]
img = Sfin * (sp.eye(16) + G5A) / 2          # image of A's +1 eigenspace in B coordinates
rowsP = [i for i in range(16) if any(img[i, j] != 0 for j in range(16))]
print("A's lamP subspace maps onto B components", rowsP)
json.dump({'S': [[str(Sfin[i, j]) for j in range(16)] for i in range(16)], 'Sst': [[str(Sst[i, j]) for j in range(16)] for i in range(16)],
           'c2': str(c2), 'forms_intertwined': bool(ok), 'rho_Sigmastar_intertwined': bool(ok2), 'A_lamP_to_B_components': rowsP},
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'spinor_map.json'), 'w'), indent=1)
