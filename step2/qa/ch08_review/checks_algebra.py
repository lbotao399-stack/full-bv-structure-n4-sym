#!/usr/bin/env python3
"""Referee checks for manuscript/ch08_definitions.tex (algebraic parts).

A  sign of the outer derivation B against ad of diag(1_4,0_4) in gl(4|4) (block table of ch02)
C  Gram-form relation of Cor c8-gram and its covariance under a finite part of the operators
D  the conformal map C of Prop c8-C
E  bound of the derivatives of P(eps,T) by those of P(0,infty) (Def c8-PT)
I  the coefficients a_k = (-1)^{k(k-1)/2}/k! of Lemma c8-ev
J  the log form of Prop c8-2pt (iii)
K  the unit for su(N)
"""
import json, itertools, math, os
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(8)
OUT = {}


def rec(key, ok, **info):
    OUT[key] = dict(passed=bool(ok), **info)
    print(('PASS ' if ok else 'FAIL ') + key, info if len(str(info)) < 400 else '')


# ---------------------------------------------------------------- A: B versus ad diag(1_4,0_4)
blk = {1: [0, 1], 2: [2, 3], 3: [4, 5, 6, 7]}          # SW (17): rows/cols 1,2 conformal 2+2, 3 = R-symmetry 4


def unit(a, b, i=0, j=0):
    X = np.zeros((8, 8)); X[blk[a][i], blk[b][j]] = 1.0; return X


Delta = np.diag([1, 1, -1, -1, 0, 0, 0, 0.])           # diag(1,-1|0) of SW (17)
Z = np.diag([1, 1, 1, 1, 0, 0, 0, 0.])                 # diag(1_4,0_4)
Zp = np.diag([0, 0, 0, 0, 1, 1, 1, 1.])                # diag(0_4,1_4)
table_ch02 = {(2, 1): -2, (1, 2): 2, (3, 1): -1, (2, 3): -1, (3, 2): 1, (1, 3): 1, (1, 1): 0, (2, 2): 0, (3, 3): 0}
ok = True
for (a, b), ev in table_ch02.items():
    X = unit(a, b)
    ok &= np.allclose(Delta @ X - X @ Delta, ev * X)
rec('A1_block_convention_reproduces_ch02_adDelta_table', ok)
names = {(3, 1): 'S+(x)W in Sigma', (2, 3): 'S-(x)W* in Sigma', (3, 2): 'S-(x)W in Sigma*', (1, 3): 'S+(x)W* in Sigma*'}
B_ch08 = {(3, 1): +1, (3, 2): +1, (2, 3): -1, (1, 3): -1}
adZ = {}
adZp = {}
for (a, b) in names:
    X = unit(a, b)
    adZ[names[(a, b)]] = float(np.sum((Z @ X - X @ Z) * X))
    adZp[names[(a, b)]] = float(np.sum((Zp @ X - X @ Zp) * X))
rec('A2_ad_diag(1_4,0_4)_equals_B_of_ch08', all(adZ[names[k]] == B_ch08[k] for k in names),
    ad_diag_1_0=adZ, B_ch08={names[k]: v for k, v in B_ch08.items()})
rec('A3_ad_diag(0_4,1_4)_equals_B_of_ch08', all(adZp[names[k]] == B_ch08[k] for k in names), ad_diag_0_1=adZp)

# ---------------------------------------------------------------- C: Gram relation and covariance
n = 5


def herm_pos(n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); return A @ A.conj().T + n * np.eye(n)


def herm(n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); return A + A.conj().T


def cplx(n):
    return rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))


G0, G1 = herm_pos(n), herm(n)
inv = np.linalg.inv
adj = lambda Y: inv(G0) @ Y.conj().T @ G0
M = inv(G0) @ G1
X0, X1 = cplx(n), cplx(n)
Y0 = adj(X0)
Y1 = adj(X1) + (Y0 @ M - M @ Y0)
h = 1e-4
G = G0 + h * G1; X = X0 + h * X1; Y = Y0 + h * Y1
resid = np.linalg.norm(X.conj().T @ G - G @ Y) / h ** 2
rec('C1_relation_(X^dag)_1=(X_1)^adj+[(X^dag)_0,G^-1G_1]_solves_conj(X)^T G=G Y_to_O(h^2)', resid < 1e3,
    residual_over_h2=float(resid))
# exact covariance: G -> g^dag G g, X -> g^-1 X g, Y -> g^-1 Y g preserves conj(X)^T G = G Y
Psi = cplx(n)
g = np.eye(n) + h * Psi
Gx = G0 + h * G1; Xx = X0 + h * X1; Yx = adj(X0) + h * (adj(X1) + (adj(X0) @ M - M @ adj(X0)))
Gp = g.conj().T @ Gx @ g; Xp = inv(g) @ Xx @ g; Yp = inv(g) @ Yx @ g
rec('C2_transformation_G->g^dag G g,X->g^-1 X g_preserves_relation',
    np.linalg.norm(Xp.conj().T @ Gp - Gp @ Yp) / h ** 2 < 1e3,
    residual_over_h2=float(np.linalg.norm(Xp.conj().T @ Gp - Gp @ Yp) / h ** 2))
# first-order rules
M_new = inv(G0) @ ((Gp - G0) / h)                      # G_tree^-1 G_1 after the change, to O(h)
dM = M_new - M
rule_M = Psi + adj(Psi)
X1_new = (Xp - X0) / h
rule_X_correct = X1 + (X0 @ Psi - Psi @ X0)          # X_1 + [X_0, psi']
rule_X_ch08 = X1 + (Psi @ X0 - X0 @ Psi)             # X_1 + [psi', X_0]   (as printed in ch08)
rec('C3_G^-1G_1_shifts_by_psi+psi^adj', np.linalg.norm(dM - rule_M) < 1e-2, err=float(np.linalg.norm(dM - rule_M)))
rec('C4_X_1_shifts_by_[X_0,psi]_(correct)', np.linalg.norm(X1_new - rule_X_correct) < 1e-2,
    err=float(np.linalg.norm(X1_new - rule_X_correct)))
rec('C5_X_1_shifts_by_[psi,X_0]_(ch08_as_printed)', np.linalg.norm(X1_new - rule_X_ch08) < 1e-2,
    err=float(np.linalg.norm(X1_new - rule_X_ch08)))


# check whether the printed pair of rules keeps the relation (first order, algebraically)
def rel_defect(X1_, Y1_, M_):
    return np.linalg.norm(Y1_ - (adj(X1_) + (Y0 @ M_ - M_ @ Y0)))


Y1_c = Y1 + (Y0 @ Psi - Psi @ Y0)
Y1_p = Y1 + (Psi @ Y0 - Y0 @ Psi)
rec('C6_printed_rules_(X1+[psi,X0],M+psi+psi^adj)_covariant', rel_defect(rule_X_ch08, Y1_p, M + rule_M) < 1e-8,
    defect=float(rel_defect(rule_X_ch08, Y1_p, M + rule_M)))
rec('C7_corrected_rules_(X1+[X0,psi],M+psi+psi^adj)_covariant', rel_defect(rule_X_correct, Y1_c, M + rule_M) < 1e-8,
    defect=float(rel_defect(rule_X_correct, Y1_c, M + rule_M)))
rec('C8_IP-P1_rules_(X1+[Z,X0],M-Z-Z^adj)_covariant',
    rel_defect(X1 + (Psi @ X0 - X0 @ Psi), Y1_p, M - rule_M) < 1e-8,
    defect=float(rel_defect(X1 + (Psi @ X0 - X0 @ Psi), Y1_p, M - rule_M)))
# outer direction: X1 -> X1 - a X0, Y1 -> Y1 + a Y0 keeps the relation iff conj(a) = -a
for a, lab in ((0.7, 'real'), (0.7j, 'imaginary'), (0.3 + 0.5j, 'generic')):
    d = rel_defect(X1 - a * X0, Y1 + a * Y0, M)
    rec('C9_outer_direction_a_' + lab, (d < 1e-8) == (lab == 'imaginary'), defect=float(d))

# ---------------------------------------------------------------- D: the conformal map C
nvec = np.array([1.0, 0, 0, 0])


def vfield(x):
    return nvec * (1 - x @ x) + 2 * (nvec @ x) * x


def flow(x, t, steps=4000):
    dt = t / steps
    for _ in range(steps):
        k1 = vfield(x); k2 = vfield(x + dt / 2 * k1); k3 = vfield(x + dt / 2 * k2); k4 = vfield(x + dt * k3)
        x = x + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return x


def Cmap(x):
    return (2 * x + ((x @ x) - 1 - 2 * (nvec @ x)) * nvec) / ((x + nvec) @ (x + nvec))


def Cinv(x):
    return (2 * x + (1 - (x @ x) - 2 * (nvec @ x)) * nvec) / ((x - nvec) @ (x - nvec))


def rn(x):
    return x - 2 * (nvec @ x) * nvec


errs = []
for _ in range(5):
    x = rng.normal(size=4) * 0.4
    errs.append(np.linalg.norm(flow(x, -np.pi / 4) - Cmap(x)))
rec('D1_C_is_time_-pi/4_flow_of_P_n-K_n', max(errs) < 1e-8, max_err=float(max(errs)))
ok = True
for _ in range(200):
    x = rng.normal(size=4); x[0] = abs(x[0]) + 1e-3
    ok &= np.linalg.norm(Cmap(x)) < 1
rec('D2_C_maps_half_space_into_unit_ball', ok)
errs = []
for _ in range(50):
    y = rng.normal(size=4)
    errs.append(np.linalg.norm(Cmap(rn(Cinv(y))) - y / (y @ y)))
    errs.append(np.linalg.norm(Cmap(Cinv(y)) - y))
rec('D3_C_r_n_C^-1_equals_inversion', max(errs) < 1e-10, max_err=float(max(errs)))

# ---------------------------------------------------------------- E: derivative bounds of P(eps,T)
r = sp.symbols('r', positive=True)
eps_, T_ = sp.Rational(1), sp.Integer(50)
fT = (sp.exp(-r ** 2 / (4 * T_)) - sp.exp(-r ** 2 / (4 * eps_))) / r ** 2
f0 = 1 / r ** 2
ratios = {}
for k in (0, 1, 2, 3):
    a = sp.lambdify(r, sp.diff(fT, r, k), 'numpy'); b = sp.lambdify(r, sp.diff(f0, r, k), 'numpy')
    rr = np.linspace(0.05, 40, 40000)
    ratios[k] = float(np.max(np.abs(a(rr)) / np.abs(b(rr))))
# transverse part of the Hessian, f'(r)/r
a = sp.lambdify(r, sp.diff(fT, r) / r, 'numpy'); b = sp.lambdify(r, sp.diff(f0, r) / r, 'numpy')
rr = np.linspace(0.05, 40, 40000)
ratios['hess_transverse'] = float(np.max(np.abs(a(rr)) / np.abs(b(rr))))
rec('E1_|d^k P(eps,T)|<=|d^k P(0,inf)|_radial_k<=3_and_Hessian', all(v <= 1 + 1e-9 for v in ratios.values()),
    max_ratio=ratios)

# ---------------------------------------------------------------- I: a_k
ok = True
for k in range(8):
    for l in range(8):
        ak = lambda m: sp.Integer(-1) ** (m * (m - 1) // 2) / sp.factorial(m)
        ok &= sp.simplify(ak(k) * ak(l) - ak(k + l) * sp.binomial(k + l, k) * (-1) ** (k * l)) == 0
rec('I1_a_k_a_l=a_{k+l}binom(k+l,k)(-1)^{kl}', ok, a_0_to_4=[str(sp.Integer(-1) ** (m * (m - 1) // 2) / sp.factorial(m)) for m in range(5)])

# ---------------------------------------------------------------- J: log form of Prop c8-2pt (iii)
x, w, ell = sp.symbols('x w ell', positive=True)
C_, D_ = sp.MatrixSymbol('C', 2, 2), sp.MatrixSymbol('D', 2, 2)
Cm = sp.Matrix(2, 2, sp.symbols('c0:4')); Dm = sp.Matrix(2, 2, sp.symbols('d0:4')); Fm = sp.Matrix(2, 2, sp.symbols('f0:4'))
Gm = x ** (-2 * w) * (Cm + Fm - sp.Rational(1, 4) * (Dm * Cm + Cm * Dm.T) * sp.log(x ** 2 / ell ** 2))
# kernel depends on |x-y| = x; x d/dx + y d/dy acts on a function of |x-y| as r d/dr
lhs = (x * sp.diff(Gm, x) + 2 * w * Gm).applyfunc(sp.simplify)
rhs = (-sp.Rational(1, 2) * (Dm * Cm + Cm * Dm.T) * x ** (-2 * w)).applyfunc(sp.simplify)
rec('J1_log_form_solves_scaling_identity', (lhs - rhs).applyfunc(sp.simplify) == sp.zeros(2, 2))

# ---------------------------------------------------------------- K: unit for su(N)
for N in (2, 3, 4):
    # basis of su(N): generalised Gell-Mann, T^a hermitian, tr T^a T^b = 1/2 delta
    Ts = []
    for i in range(N):
        for j in range(i + 1, N):
            E = np.zeros((N, N), complex); E[i, j] = E[j, i] = 0.5; Ts.append(E)
            E = np.zeros((N, N), complex); E[i, j] = -0.5j; E[j, i] = 0.5j; Ts.append(E)
    for k in range(1, N):
        E = np.zeros((N, N), complex)
        for i in range(k): E[i, i] = 1
        E[k, k] = -k
        Ts.append(E / np.sqrt(2 * k * (k + 1)))
    X = sum(rng.normal() * T for T in Ts); Y = sum(rng.normal() * T for T in Ts)
    Xa, Ya = -1j * X, -1j * Y                      # compact form
    def ad(Z):
        cols = []
        for T in Ts:
            Cz = Z @ T - T @ Z
            cols.append([2 * np.trace(Cz @ S).real if np.allclose(Cz, Cz.conj().T) else 0 for S in Ts])
        return None
    # tr_g(ad_X ad_Y) via structure: use the matrix of ad on the basis -i T^a
    B = [-1j * T for T in Ts]
    def coords(Z):
        return np.array([(-2 * np.trace(Z @ b)).real for b in B])   # Z = sum z_a b_a, tr(b_a b_b) = -1/2 delta
    def admat(Z):
        return np.array([coords(Z @ b - b @ Z) for b in B]).T
    killing = np.trace(admat(Xa) @ admat(Ya))
    trXY = np.trace(Xa @ Ya).real
    g, hb = 0.9, 1.0
    pairing = -(4 * hb / g ** 2) * trXY
    unitval = hb * killing / (16 * np.pi ** 2 * pairing)
    rec('K1_unit_su(%d)=-g^2N/(32pi^2)' % N, abs(unitval + g ** 2 * N / (32 * np.pi ** 2)) < 1e-10,
        unit=float(unitval), expected=float(-g ** 2 * N / (32 * np.pi ** 2)),
        konishi_gamma_from_minus24_units=float(-24 * unitval), expected_3g2N_over_4pi2=float(3 * g ** 2 * N / (4 * np.pi ** 2)))

json.dump(OUT, open(os.path.join(HERE, 'checks_algebra.json'), 'w'), indent=1)
print('wrote checks_algebra.json')
