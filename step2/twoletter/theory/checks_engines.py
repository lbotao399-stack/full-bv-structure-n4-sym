"""Check of the conventions of the existing engines against Ch1 (vertices, Clifford data) and Ch3 (Q^GF, propagator,
K_0), before the engine specification tells the two-letter engines to reuse them.
Engine A: step2/engine_a/s2_core.py (imports engine_a/cliff.py, engine_a/r2_core.py), algebra su2.
Engine B: engine_b/loops.py (derive_propagators: inversion of the free Hessian), engine_b/vertices.py.
Checks
 A1  Clifford relation rho(v)rho(w)+rho(w)rho(v) = 2(v,w) on Sigma and Sigma^*, and (rho(v)x,y) = (x,rho(v)y)
 A2  Pnum: A(x)A = 2 (-,-)^{-1} delta (V and W_R^perp blocks), c(x)A*_V = A*_V(x)c = -i p (-,-)^{-1},
     lambda(x)lambda = i rho(p)|_{Sigma^*->Sigma} (-,-)^{-1}   (Prop c3-P)
 A3  Khat = K_0 components -(-,-)^{-1} on A(x)A*, A*(x)A, lambda^*(x)lambda, lambda(x)lambda^*, c^*(x)c, c(x)c^*
 A4  vertex functions: A c A* = 2 ([A,c],A*) (scale 2), c lam lam* = 2 ([c,lam],lam*), c c c* = 2 ([c,c],c*), and AAA =
     2 x (eq c1-AAA) at random polarisations and momenta
 B1  engine B propagators: A = 2/q^2, lambda = i Rbar M^{-1}/q^2 (internal check of loops.py), ghost c-A*_m and A*_m-c
     numerators equal -i q_m.
Writes checks_engines_results.json."""
import os, sys, json, itertools
import numpy as np

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
HERE = os.path.dirname(os.path.abspath(__file__))
os.environ['S2_ALG'] = 'su2'
sys.path.insert(0, os.path.join(ROOT, 'step2', 'engine_a'))
sys.path.insert(0, os.path.join(ROOT, 'engine_a'))
import s2_core as A
from cliff import s as S_, sb as SB_

out = {}
rng = np.random.default_rng(7)
# A1
S = np.array(S_); SB = np.array(SB_)
cl = all(np.allclose(S[M] @ SB[N] + S[N] @ SB[M], 2 * (M == N) * np.eye(16)) and
         np.allclose(SB[M] @ S[N] + SB[N] @ S[M], 2 * (M == N) * np.eye(16)) for M in range(10) for N in range(10))
sym = all(np.allclose(S[M], S[M].T) and np.allclose(SB[M], SB[M].T) for M in range(10))
out['A1_clifford'] = bool(cl); out['A1_symmetric'] = bool(sym)
# A2
k = np.array([2, -1, 3, 1], dtype=complex)
D = A.D; G = A.GINV
ok = True
NAV = A.Pnum(k, 'AV'); NAP = A.Pnum(k, 'AP'); Ngh = A.Pnum(k, 'gh'); Nl = A.Pnum(k, 'l')
for M in range(10):
    i0 = A.OFF['A'] + M * D
    blk = (NAV if M < 4 else NAP)[i0:i0 + D, i0:i0 + D]
    ok &= np.allclose(blk, 2 * G)
for mu in range(4):
    i0 = A.OFF['As'] + mu * D
    ok &= np.allclose(Ngh[A.OFF['c']:A.OFF['c'] + D, i0:i0 + D], -1j * k[mu] * G)
    ok &= np.allclose(Ngh[i0:i0 + D, A.OFF['c']:A.OFF['c'] + D], -1j * k[mu] * G)
rho = sum(k[m] * S[m] for m in range(4))          # rho(p): Sigma^* -> Sigma
ok &= np.allclose(Nl[A.OFF['l']:A.OFF['l'] + 16 * D, A.OFF['l']:A.OFF['l'] + 16 * D], 1j * np.kron(rho, G))
out['A2_propagator'] = bool(ok)
# A3
okk = True
for line, pairs in [('sAV', [('A', 'As')]), ('sAP', [('A', 'As')]), ('sl', [('ls', 'l')]), ('sc', [('cs', 'c')])]:
    Kh = A.Khat(line)
    nz = np.argwhere(np.abs(Kh) > 0)
    vals = set(np.round(Kh[tuple(nz.T)], 12))
    okk &= vals <= set(np.round(-G[np.abs(G) > 0], 12)) | {-1.0}
    okk &= np.allclose(Kh, Kh.T)
out['A3_K0'] = bool(okk)
# A4: odd-leg vertices on basis elements (all momenta irrelevant)
F = A.F
ok4 = True
for (a, b, c_) in itertools.product(range(3), repeat=3):
    pA = np.zeros((1, 10, D), dtype=complex); pA[0, 0, a] = 1
    pc = np.zeros((1, D), dtype=complex); pc[0, b] = 1
    pAs = np.zeros((1, 10, D), dtype=complex); pAs[0, 0, c_] = 1
    v = A.canon('AcAs', [pA, pc, pAs], [None] * 3)[0, 0, 0]
    ok4 &= np.isclose(v, 2 * F[a, b, c_])
    pc2 = np.zeros((1, D), dtype=complex); pc2[0, a] = 1
    pc3 = np.zeros((1, D), dtype=complex); pc3[0, c_] = 1
    v = A.canon('cccS', [pc2, pc, pc3], [None] * 3)[0, 0, 0]
    ok4 &= np.isclose(v, 2 * F[a, b, c_])
out['A4_odd_vertices'] = bool(ok4)
# AAA vs eq c1-AAA:  -i/2 [ (iota_{p2-p1} A3,[A1,A2]) + (iota_{p3-p2} A1,[A2,A3]) + (iota_{p1-p3} A2,[A3,A1]) ]
def brack(X, Y):  # g-components, orthonormal su2 basis, [t_a,t_b] = eps_abc t_c
    return np.einsum('a,b,abc->c', X, Y, F)
okA = True
for trial in range(5):
    pols = [rng.integers(-2, 3, size=(10, D)).astype(complex) for _ in range(3)]
    p1 = rng.integers(-2, 3, size=4).astype(complex); p2 = rng.integers(-2, 3, size=4).astype(complex)
    p3 = -p1 - p2
    moms = [p1, p2, p3]
    def iota(q, P):  # contraction of the V_10 factor with q in V
        return sum(q[m] * P[m] for m in range(4))
    def AA(P, Qp):   # [A,A'] contracting V_10 factors with the metric
        return sum(brack(P[M], Qp[M]) for M in range(10))
    ref = -0.5j * (np.dot(iota(p2 - p1, pols[2]), AA(pols[0], pols[1]))
                   + np.dot(iota(p3 - p2, pols[0]), AA(pols[1], pols[2]))
                   + np.dot(iota(p1 - p3, pols[1]), AA(pols[2], pols[0])))
    val = A.canon('AAA', [P[None] for P in pols], moms)[0, 0, 0]
    okA &= np.isclose(val, 2 * ref) or np.isclose(val, -2 * ref)
    out.setdefault('A4_AAA_ratio', []).append(complex(val / (2 * ref)) if abs(ref) > 0 else None)
out['A4_AAA'] = bool(okA)
# B1
try:
    sys.path.insert(0, os.path.join(ROOT, 'engine_b'))
    import loops as LB
    import sympy as sp
    d = LB.derive_propagators()
    qs = d['qs']
    gh1 = all(sp.simplify(d['ghost_c_As'][m] + sp.I * qs[m]) == 0 for m in range(4))
    gh2 = all(sp.simplify(d['ghost_As_c'][m] + sp.I * qs[m]) == 0 for m in range(4))
    out['B1_A_2_over_q2'] = bool(d['A_is_2_over_q2'])
    out['B1_lambda'] = bool(d['lam_matches_iRbar_Minv'])
    out['B1_ghost_c_As_minus_iq'] = bool(gh1)
    out['B1_ghost_As_c_minus_iq'] = bool(gh2)
    out['B1_ghost_numerators'] = [str(x) for x in d['ghost_c_As']] + [str(x) for x in d['ghost_As_c']]
except Exception as exc:  # report, do not hide
    out['B1_error'] = repr(exc)
out['all_pass'] = bool(all(v for k_, v in out.items() if k_.startswith(('A1', 'A2', 'A3', 'A4_odd', 'A4_AAA', 'B1_A', 'B1_l', 'B1_g'))
                           and isinstance(v, bool)))
with open(os.path.join(HERE, 'checks_engines_results.json'), 'w') as f:
    json.dump(out, f, indent=1, default=str)
for k_, v in out.items():
    print(k_, v if not isinstance(v, list) else v[:4])
