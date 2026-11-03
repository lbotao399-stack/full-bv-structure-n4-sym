"""QA: basis-free check of the twisted letters of engines A and B in B's realisation.
(1) Hermitian invariance of the standard inner product on Sigma under real rotations (B realisation).
(2) For each data set (xi, beta letters, gamma outputs, dc outputs), decompose xi = u (x) z and test:
    beta = u^dagger (x) t with t(z) = 0;  dc outputs Y = s (x) l with l proportional to conj(z) (Hermitian dual of z);
    gamma outputs in conj(N), N = ker rho(.) xi."""
import os, sys, json
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'qa/src.nosync'))
import qq_ab_compare as AB
import realization as RZ
R, Rb = RZ.R, RZ.Rb
herm = all(sp.simplify((Rb[a] * R[b]).H + Rb[a] * R[b]) == sp.zeros(16, 16) for a in range(10) for b in range(10) if a != b)
print('standard Hermitian form invariant under rotations (B):', herm)
def mat(v, off):      # 2x4 matrix of a 16-vector restricted to components off..off+7 (index 4a+A)
    return sp.Matrix(2, 4, lambda a, A: v[off + 4 * a + A])
def analyse(name, xi, betas, dcs, gams):
    off = 0 if any(xi[i] != 0 for i in range(8)) else 8
    X = mat(xi, off)
    assert X.rank() == 1
    # u = a nonzero column, z = a nonzero row (rank one: X = u z^T)
    j = [k for k in range(4) if any(X[a, k] != 0 for a in range(2))][0]
    u = X[:, j]
    i = [a for a in range(2) if X[a, j] != 0][0]
    z = X[i, :].T / u[i]
    assert sp.simplify(u * z.T - X) == sp.zeros(2, 4)
    print(name, 'chirality block', off, 'u =', list(u), 'z =', list(z))
    for k, b in enumerate(betas):
        Bm = mat(b, off)
        # beta(u' (x) w) = 0 for u' perp u (Hermitian): u' = (-conj(u1), conj(u0))
        up = sp.Matrix([-sp.conjugate(u[1]), sp.conjugate(u[0])])
        c1 = sp.simplify(up.T * Bm)            # 1x4: must vanish
        c2 = sp.simplify(u.T * Bm * z)         # beta(u (x) z): must vanish
        print('  beta', k, 'kills uperp (x) W:', c1 == sp.zeros(1, 4), ' kills u (x) z:', c2 == sp.zeros(1, 1))
    off2 = 8 - off
    for k, y in enumerate(dcs):
        Y = mat(y, off2)
        # Y = s (x) l, l must annihilate z-perp (Hermitian): l proportional to conj(z)
        zc = sp.Matrix([sp.conjugate(x) for x in z])
        ok = all(sp.simplify(Y[a, A] * zc[B] - Y[a, B] * zc[A]) == 0 for a in range(2) for A in range(4) for B in range(4))
        print('  dc out', k, 'W-part proportional to conj(z):', ok, ' W-part proportional to z:',
              all(sp.simplify(Y[a, A] * z[B] - Y[a, B] * z[A]) == 0 for a in range(2) for A in range(4) for B in range(4)))
    # gamma: N = {w: rho(w) xi = 0}
    Nm = sp.Matrix([[v for v in (sum((w * R[4 + K] * sp.Matrix(xi) for K, w in [(K, 1)]), sp.zeros(16, 1)))] for K in range(6)])
    cols = sp.Matrix.hstack(*[R[4 + K] * sp.Matrix(xi) for K in range(6)])
    N = cols.nullspace()
    Nbar = [sp.Matrix([sp.conjugate(x) for x in n]) for n in N]
    for k, g in enumerate(gams):
        gv = sp.Matrix(g)
        M = sp.Matrix.hstack(*Nbar, gv)
        print('  gamma out', k, 'in conj(N):', M.rank() == len(Nbar), ' in N:', sp.Matrix.hstack(*N, gv).rank() == len(N))
# engine A data in B basis
d = json.load(open(os.path.join(QSK, 'qa/out.nosync/twistA_data.json')))
xiA = list(sp.expand(AB.S * sp.Matrix([AB.cplx(x) for x in d['xi']])))
bA = []
for z_ in d['beta_z']:
    t, dd = AB.to_B('lamP', z_, 'letter'); v = [0] * 16
    for c, x in dd.items(): v[c] = x
    bA.append(v)
dA = []
for y in d['dc_out']:
    t, dd = AB.to_B('lamM', y, 'out'); v = [0] * 16
    for c, x in dd.items(): v[c] = x
    dA.append(v)
gA = [[AB.cplx(x) for x in g] for g in d['gam_out']]
analyse('engine A twist (in B basis)', xiA, bA, dA, gA)
# engine B twist data
xiB = [0] * 16; xiB[3] = 1
bB = []
for I in (1, 2, 3):
    v = [0] * 16; v[I - 1] = 1; bB.append(v)
dB = []
for c in (11, 15):
    v = [0] * 16; v[c] = 1; dB.append(v)
gam = [[RZ.Sig[K][J - 1, 3] for K in range(6)] for J in (1, 2, 3)]
comp = [[RZ.Sig[K][I - 1, J - 1] for K in range(6)] for (I, J) in ((1, 2), (1, 3), (2, 3))]
Linv = sp.Matrix(gam + comp).inv()
gB = [list(Linv[:, k]) for k in range(3)]
analyse('engine B twist', xiB, bB, dB, gB)
