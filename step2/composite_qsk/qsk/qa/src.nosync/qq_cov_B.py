"""QA: covariance test of engine B's code (read-only) on B's twisted beta-beta data under spin lifts of 90-degree
rotations in the planes (0,1) of V and (4,5) of W_perp.  Lambda = (1 + rho_a rho_b) (scale sqrt 2 x spin lift);
the kernel is multilinear with net Lambda-scale weight 0 here (xi, one spinor output: +2; two spinor letters: -2)."""
import os, sys, json, random
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'qa/src.nosync'))
import qq_mirror_B as MB          # runs the mirror test again on import (cheap) and defines value()
import realization as RZ
R, Rb = RZ.R, RZ.Rb
def lam(a, b):
    return sp.eye(16) + Rb[a] * R[b]
def vrot(a, b):
    O = sp.eye(10)
    # Lambda rho(e_a) Lambda^-1 = -rho(e_b), Lambda rho(e_b) Lambda^-1 = rho(e_a)
    O[a, a] = 0; O[b, b] = 0; O[b, a] = -1; O[a, b] = 1
    return O
xi = sp.zeros(16, 1); xi[3] = 1
zaB = sp.zeros(16, 1); zaB[0] = 1
zbB = sp.zeros(16, 1); zbB[1] = 1
gam = [[RZ.Sig[K][J - 1, 3] for K in range(6)] for J in (1, 2, 3)]
comp = [[RZ.Sig[K][I - 1, J - 1] for K in range(6)] for (I, J) in ((1, 2), (1, 3), (2, 3))]
Linv = sp.Matrix(gam + comp).inv()
Y0 = sp.zeros(10, 1)
for c in range(6):
    Y0[4 + c] = Linv[c, 2]
Y1 = sp.zeros(16, 1); Y1[11] = 1; Y1[15] = 2
random.seed(3)
moms = [sp.Matrix([random.randint(-3, 3) for _ in range(4)]) for _ in range(2)]
def run(L, O):
    xi2 = L * xi; za2 = L.inv().T * zaB; zb2 = L.inv().T * zbB; Y02 = O * Y0; Y12 = L * Y1
    m2 = [(O[:4, :4] * m) for m in moms]
    offl = 0 if any(xi2[i] != 0 for i in range(8)) else 8
    tl = 'lamP' if offl == 0 else 'lamM'; to = 'lamM' if offl == 0 else 'lamP'
    return MB.value(tl, tl, MB.vecd(za2, range(offl, offl + 8)), MB.vecd(zb2, range(offl, offl + 8)), ['AW', to],
                    [{c: Y02[c] for c in range(4, 10)}, MB.vecd(Y12, range(8 - offl, 16 - offl))],
                    [list(m) for m in m2], [xi2[i] for i in range(16)])
# sanity: the vector rotation is the one induced by Lambda
for (a, b) in ((0, 1), (4, 5)):
    L = lam(a, b); O = vrot(a, b)
    ok = all(sp.simplify(L * Rb[M] * R[N] * L.inv() - sum((O[P, M] * O[Q, N] * Rb[P] * R[Q] for P in range(10) for Q in range(10)), sp.zeros(16, 16))) == sp.zeros(16, 16) for M in range(10) for N in range(10) if M != N and (M in (a, b) or N in (a, b)))
    print('plane', (a, b), 'vector rotation consistent:', ok, flush=True)
base = run(sp.eye(16), sp.eye(10))
print('base', base, flush=True)
out = {'base': {k: str(v) for k, v in base.items()}}
for (a, b) in ((0, 1), (4, 5), (1, 3)):
    v = run(lam(a, b), vrot(a, b))
    print('rotated', (a, b), v, flush=True)
    out[str((a, b))] = {k: str(x) for k, x in v.items()}
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/cov_B.json'), 'w'), indent=1)
