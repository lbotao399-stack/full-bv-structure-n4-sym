# Check for odd_components_design.tex, Lemma lem:o-ev (evaluation along a rank-one xi):
#   for xi in S_+ (x) W of rank one: Gamma(xi,xi) = 0 (all ten components),
#   X := {S^(1)(xi), S^(1)(xi)} is a nonzero local functional,
#   X - {S, S^(2)(xi,xi)} = 0 and {S^(1)(xi), S^(2)(xi,xi)} = 0 as local functionals,
# hence the twisted action S + S^(1)(xi) + c S^(2)(xi,xi) solves the CME iff c = -1/2.
# Read-only use of the Step-1 superconformal engine (su(2), exact Gaussian-rational arithmetic).
import os, sys, time, json
import numpy as np
from fractions import Fraction as Fr
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'superconformal'))
import engine as E, model as Mo
t0 = time.time()
# rank-one xi: simultaneous eigenvector of sigbar[M]sig[N], (M,N)=(0,1),(2,3),(4,5),(6,7),(8,9), weights (1,-1,1,1,1),
# in the Omega=+1 eigenspace (S_+ (x) W, ch01 Prop c1-chirality); see find_rank_one.py
np.seterr(all='ignore')
sig = [np.array([[complex(*map(float, Mo.SIG[M][a][b])) for b in range(16)] for a in range(16)]) for M in range(10)]
sigb = [np.array([[complex(*map(float, Mo.SIGB[M][a][b])) for b in range(16)] for a in range(16)]) for M in range(10)]
mm = lambda *A: __import__('functools').reduce(lambda X, Y: X.dot(Y), A)
Om = mm(sigb[0], sig[1], sigb[2], sig[3])
H = [mm(sigb[M], sig[N]) for (M, N) in [(0,1),(2,3),(4,5),(6,7),(8,9)]]
Pj = np.eye(16, dtype=complex)
for s_, h in zip((1, -1, 1, 1, 1), H): Pj = Pj.dot((np.eye(16) - 1j*s_*h)/2)
for a in range(16):
    v = Pj[:, a]
    if np.linalg.norm(v) > 1e-9: break
v = v/v[np.argmax(abs(v))]
Qexact = [(Fr(int(round(z.real))), Fr(int(round(z.imag)))) for z in v]
Qn = np.array([float(r)+1j*float(i) for r, i in Qexact])
a0 = int(a)
res = {'basis_index': a0, 'max_|Gamma(Q,Q)|': float(max(abs(Qn @ sig[M] @ Qn) for M in range(10))),
       'Omega_Q_minus_Q': float(np.linalg.norm(Om @ Qn - Qn)),
       'Gamma_V_rank_of_Gamma(Q,-)_on_Sigma': int(np.linalg.matrix_rank(np.array([sig[M] @ Qn for M in range(4)]), tol=1e-9))}
print(res, flush=True)
V = E.Vars()
S = Mo.build_S(V)
zero = [(Fr(0), Fr(0))]*16
S1 = Mo.S1_odd(V, Qexact, zero)
ks = Mo.killing_spinor(V, Qexact, zero)
S2 = Mo.S2_odd(V, ks, ks)
print('built S, S1, S2: terms', len(S), len(S1), len(S2), round(time.time()-t0, 1), 's', flush=True)
X = E.bracket(V, S1, S1)
Y = E.bracket(V, S, S2)
Z = E.bracket(V, S1, S2)
okX, _ = E.total_derivative_test(V, X)
okXY, _ = E.total_derivative_test(V, E.psub(X, Y))
okZ, _ = E.total_derivative_test(V, Z)
res.update({'X_is_zero_functional': bool(okX), 'X_minus_bracket_S_S2_is_zero': bool(okXY),
            'bracket_S1_S2_is_zero': bool(okZ), 'seconds': round(time.time()-t0, 1),
            'conclusion': 'CME of S+S1+cS2 holds iff (2c+1)X=0, X!=0 => c=-1/2' if (not okX and okXY and okZ) else 'CHECK FAILED'})
print(json.dumps(res, indent=1))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'twist_coeff_check.json'), 'w'), indent=1)
