#!/usr/bin/env python3
"""Checks for step2/qa/ip_completeness2.json (completeness critic of the inner-product design).

K1  two-component form of the inversion Theta_I on lambda, lambdabar, lambda*, lambdabar*
    (proposed addition to ch08 prop:c8-C), derived from the index-free statement
    Theta_I lambda(x) = -(x,x)^{-2} rho(x) conj(lambda(I x)),
    Theta_I lambda*(x) = +(x,x)^{-3} rho(x) conj(lambda*(I x)),
    with the formulas of ch08 def:c8-twocomp (eq:c8-pairing) and lem:c8-cc (eq:c8-cc).
    Convention gate: the same code with +rho(n) (resp. -rho(n)) and r_n reproduces the
    printed two-component Theta_n of ch08 def:c8-Theta.  Also Theta_I^2 = 1 on spinors and
    on the vector A_m (weight 1, Jacobian (delta_mk - 2 x_m x_k/(x,x))).
K2  reflection compatibility of the inner direction: the classical dilatation 2(x.d + w)
    on smeared densities commutes with Theta_n, so psi' = c {S^(1)(Delta),-} commutes with
    the conjugation exactly for real c (sign algebra of thm:c8-fields (i), (ii) with the ch01
    symmetry {g,f} = (-1)^{|f||g|}{f,g} gives {S^(1)(Delta),O}^Theta = {S^(1)(Delta),O^Theta}).
K3  K-line bookkeeping (cor:c8-gram, third line) with P_1 = p {S^(1)(P_b),-}, p real:
    -p K + [K, M_tp] = p K + [K, M_tp - 2 p W] in a lowest-weight sl(2) module, i.e. the
    form 'pbar ad + [ad, Gram_tree^{-1} Gram_1]' holds with the radial Gram of that frame,
    (1 - 2 p w) x two-point coefficient, not with the two-point coefficient itself.
K4  grading of Lambda^(1)(x) O on odd O: with Lambda^(1)(x) = (-1)^F ad_{S(x)} (total-parity
    reading of thm:c8-L (iv) for even x) the representation property of cor:c8-homotopy (ii)
    fails on odd O; with Lambda^(1)(x) = ad_{S(x)} (degree reading) it holds (toy model:
    gl(2) acting on C^{1|1} tensor C^2, F the parity operator).
Run: python3 ip_completeness2_checks.py  (seconds, < 100 MB).  Output: ip_completeness2_checks.json
"""
import json, os
import numpy as np
import sympy as sp

rng = np.random.default_rng(20261029)
I2 = np.eye(2, dtype=complex)
s1 = np.array([[0, 1], [1, 0]], dtype=complex)
s2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
s3 = np.array([[1, 0], [0, -1]], dtype=complex)
SIG = [s1, s2, s3, 1j * I2]          # sigma^m_{alpha alphadot}
SIGB = [-s1, -s2, -s3, 1j * I2]      # sigmabar^{m alphadot alpha}

def rand_c(shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)

# spinors: Sigma element (xi[a,i], xib[ad,i]); Sigma* element (eta[a,i], etab[ad,i]); i = 0..3
def rho_on_Sigmastar(v, eta, etab):
    # (rho(v)eta)_a^i = -i v_m sigma^m_{a ad} etab^{ad i};  (rho(v)eta)_{ad i} = i v_m eta^a_i sigma^m_{a ad}
    xi = sum(-1j * v[m] * SIG[m] @ etab for m in range(4))
    xib = sum(1j * v[m] * SIG[m].T @ eta for m in range(4))
    return xi, xib

def rho_on_Sigma(v, xi, xib):
    # (rho(v)xi)^a_i = i v_m xib_{ad i} sigmabar^{m ad a};  (rho(v)xi)^{ad i} = -i v_m sigmabar^{m ad a} xi_a^i
    eta = sum(1j * v[m] * SIGB[m].T @ xib for m in range(4))
    etab = sum(-1j * v[m] * SIGB[m] @ xi for m in range(4))
    return eta, etab

def conj_Sigma(xi, xib):       # Sigma -> Sigma*: eta = conj(xi), etab = -conj(xib)
    return np.conj(xi), -np.conj(xib)

def conj_Sigmastar(eta, etab):  # Sigma* -> Sigma: inverse of the above
    return np.conj(eta), -np.conj(etab)

def close(a, b, tol=1e-11):
    return bool(np.max(np.abs(np.asarray(a) - np.asarray(b))) < tol)

res = {}

# ---------- K1 ----------
k1 = {}
v = rng.normal(size=4)
xi, xib = rand_c((2, 4)), rand_c((2, 4))
eta, etab = rand_c((2, 4)), rand_c((2, 4))
e1, e2 = rho_on_Sigma(v, *rho_on_Sigmastar(v, eta, etab))
k1["clifford_rho(v)^2=(v,v)_on_Sigma*"] = close(e1, (v @ v) * eta) and close(e2, (v @ v) * etab)
x1, x2 = rho_on_Sigmastar(v, *rho_on_Sigma(v, xi, xib))
k1["clifford_rho(v)^2=(v,v)_on_Sigma"] = close(x1, (v @ v) * xi) and close(x2, (v @ v) * xib)
pair = np.sum(xi * np.conj(xi)) - np.sum(xib * (-np.conj(xib)))
k1["pairing_(xi,conj xi)_positive"] = bool(abs(pair.imag) < 1e-12 and pair.real > 0)
a1, a2 = conj_Sigma(*rho_on_Sigmastar(v, eta, etab))
b1, b2 = rho_on_Sigma(v, *conj_Sigmastar(eta, etab))
k1["conj_intertwines_rho_real_v"] = close(a1, b1) and close(a2, b2)

# fields lambda in Sigma (components at a point), lambda* in Sigma*
lam, lamb = rand_c((2, 4)), rand_c((2, 4))
lst, lstb = rand_c((2, 4)), rand_c((2, 4))
n = rng.normal(size=4); n /= np.linalg.norm(n)
# printed ch08 def:c8-Theta (values of the fields at the reflected point are the arguments)
printed_l = sum(1j * n[m] * SIG[m] @ np.conj(lamb) for m in range(4))
printed_lb = sum(1j * n[m] * SIG[m].T @ np.conj(lam) for m in range(4))
printed_s = sum(1j * n[m] * SIGB[m].T @ np.conj(lstb) for m in range(4))
printed_sb = sum(1j * n[m] * SIGB[m] @ np.conj(lst) for m in range(4))
c1, c2 = rho_on_Sigmastar(n, *conj_Sigma(lam, lamb))            # +rho(n) conj(lambda)
d1, d2 = rho_on_Sigma(n, *conj_Sigmastar(lst, lstb)); d1, d2 = -d1, -d2   # -rho(n) conj(lambda*)
k1["convention_gate_Theta_n_lambda_reproduces_ch08"] = close(c1, printed_l) and close(c2, printed_lb)
k1["convention_gate_Theta_n_lambdastar_reproduces_ch08"] = close(d1, printed_s) and close(d2, printed_sb)

x = rng.normal(size=4); xx = x @ x
# index-free Theta_I at x (field values at I(x) are the arguments lam, lamb, lst, lstb)
f1, f2 = rho_on_Sigmastar(x, *conj_Sigma(lam, lamb)); f1, f2 = -f1 / xx**2, -f2 / xx**2
g1, g2 = rho_on_Sigma(x, *conj_Sigmastar(lst, lstb)); g1, g2 = g1 / xx**3, g2 / xx**3
# proposed two-component forms
p_l = sum(-1j * x[m] * SIG[m] @ np.conj(lamb) for m in range(4)) / xx**2
p_lb = sum(-1j * x[m] * SIG[m].T @ np.conj(lam) for m in range(4)) / xx**2
p_s = sum(-1j * x[m] * SIGB[m].T @ np.conj(lstb) for m in range(4)) / xx**3
p_sb = sum(-1j * x[m] * SIGB[m] @ np.conj(lst) for m in range(4)) / xx**3
k1["Theta_I_lambda_two_component_matches"] = close(f1, p_l) and close(f2, p_lb)
k1["Theta_I_lambdastar_two_component_matches"] = close(g1, p_s) and close(g2, p_sb)

# involution: fields as linear functions of the position
L0, L1 = rand_c((2, 4)), rand_c((4, 2, 4))
LB0, LB1 = rand_c((2, 4)), rand_c((4, 2, 4))
def lam_field(y):
    return L0 + np.tensordot(y, L1, axes=1), LB0 + np.tensordot(y, LB1, axes=1)
def theta_I_lam(field):
    def g(y):
        yy = y @ y; Iy = y / yy
        a, b = field(Iy)
        r1, r2 = rho_on_Sigmastar(y, *conj_Sigma(a, b))
        return -r1 / yy**2, -r2 / yy**2
    return g
y = rng.normal(size=4)
tt = theta_I_lam(theta_I_lam(lam_field))(y)
k1["Theta_I_squared_identity_on_lambda"] = close(tt[0], lam_field(y)[0]) and close(tt[1], lam_field(y)[1])
S0, S1 = rand_c((2, 4)), rand_c((4, 2, 4)); SB0, SB1 = rand_c((2, 4)), rand_c((4, 2, 4))
def lst_field(y):
    return S0 + np.tensordot(y, S1, axes=1), SB0 + np.tensordot(y, SB1, axes=1)
def theta_I_lst(field):
    def g(y):
        yy = y @ y; Iy = y / yy
        a, b = field(Iy)
        r1, r2 = rho_on_Sigma(y, *conj_Sigmastar(a, b))
        return r1 / yy**3, r2 / yy**3
    return g
tt = theta_I_lst(theta_I_lst(lst_field))(y)
k1["Theta_I_squared_identity_on_lambdastar"] = close(tt[0], lst_field(y)[0]) and close(tt[1], lst_field(y)[1])
A0, A1 = rand_c(4), rand_c((4, 4))
def A_field(y):
    return A0 + y @ A1
def theta_I_A(field):
    def g(y):
        yy = y @ y; R = np.eye(4) - 2 * np.outer(y, y) / yy
        return R @ np.conj(field(y / yy)) / yy
    return g
k1["Theta_I_squared_identity_on_A_m"] = close(theta_I_A(theta_I_A(A_field))(y), A_field(y))
res["K1_Theta_I_two_component"] = k1

# ---------- K2 ----------
k2 = {}
X = sp.symbols('x1:5', real=True)
w = sp.Rational(5, 2)
coeffs = [sp.Integer(int(a)) + sp.I * sp.Integer(int(b)) for a, b in rng.integers(-5, 6, size=(15, 2))]
monos = [sp.Integer(1)] + list(X) + [X[i] * X[j] for i in range(4) for j in range(i, 4)]
f = sum(c * m for c, m in zip(coeffs, monos))
nn = sp.Matrix([1, 2, -1, 3]) / sp.sqrt(15)
def refl(expr):
    xv = sp.Matrix(X)
    rx = xv - 2 * (nn.dot(xv)) * nn
    return sp.conjugate(expr.subs(dict(zip(X, rx)), simultaneous=True))
def Dil(expr):
    return 2 * (sum(X[i] * sp.diff(expr, X[i]) for i in range(4)) + w * expr)
k2["Theta_n_commutes_with_classical_dilatation"] = sp.simplify(sp.expand(refl(Dil(f)) - Dil(refl(f)))) == 0
cc = sp.symbols('c')
lhs = refl(cc * Dil(f)); rhs = cc * Dil(refl(f))
diff = sp.simplify(sp.expand(lhs - rhs).subs(cc, sp.Symbol('a', real=True) + sp.I * sp.Symbol('b', real=True)))
k2["compatibility_condition"] = "Theta_n(c D f) - c D(Theta_n f) = (conj(c) - c) D(Theta_n f): zero for all f iff c real"
k2["difference_vanishes_for_real_c"] = sp.simplify(diff.subs(sp.Symbol('b', real=True), 0)) == 0
k2["difference_nonzero_for_imaginary_c"] = sp.simplify(diff.subs({sp.Symbol('a', real=True): 0, sp.Symbol('b', real=True): 1})) != 0
# sign algebra: {f,g}^Th = (-1)^{|f|+|g|}{g^Th,f^Th}, S^Th = -S (|S| = 1), {g,f} = (-1)^{|f||g|}{f,g}
signs = {}
for pO in (0, 1):
    s = (-1) ** (1 + pO)          # reflection rule
    s *= -1                       # S^Theta = -S
    s *= (-1) ** (1 * pO)         # {O^Th,S} = (-1)^{|O||S|}{S,O^Th}
    signs["|O|=%d" % pO] = s
k2["sign_of_{S1(Delta),O}^Theta_relative_to_{S1(Delta),O^Theta}"] = signs
res["K2_inner_direction_reality"] = k2

# ---------- K3 ----------
k3 = {}
N = 5
h, p = sp.symbols('h p', real=True)
W = sp.diag(*[h + k for k in range(N)])
Kmat = sp.zeros(N, N)
for k in range(1, N):
    Kmat[k - 1, k] = k * (2 * h + k - 1)
M = sp.Matrix(N, N, lambda i, j: sp.Symbol('m_%d_%d' % (i, j)))
k3["[K,W]=K"] = sp.simplify(Kmat * W - W * Kmat - Kmat) == sp.zeros(N, N)
lhs = -p * Kmat + (Kmat * M - M * Kmat)
Mf = M - 2 * p * W
rhs = p * Kmat + (Kmat * Mf - Mf * Kmat)
k3["-pK+[K,M_twopoint] == pK+[K,M_twopoint-2pW]"] = sp.simplify(lhs - rhs) == sp.zeros(N, N)
k3["-pK+[K,M_twopoint] == pK+[K,M_twopoint] (two-point Gram in the frame formula)"] = sp.simplify(lhs - (p * Kmat + (Kmat * M - M * Kmat))) == sp.zeros(N, N)
res["K3_K_line_frames"] = k3

# ---------- K4 ----------
k4 = {}
# gl(2) acting on V = C^{1|1} (x) C^2 through the second factor; F = parity on the first factor
E = [np.zeros((2, 2)) for _ in range(4)]
E[0][0, 0] = 1; E[1][0, 1] = 1; E[2][1, 0] = 1; E[3][1, 1] = 1
F = np.kron(np.diag([1.0, -1.0]), np.eye(2))
ad = [np.kron(np.eye(2), e) for e in E]       # even, parity-preserving operators
def br(a, b):
    return a @ b - b @ a
ok_deg, ok_tot = True, True
for i in range(4):
    for j in range(4):
        cij = br(E[i], E[j]); Cij = np.kron(np.eye(2), cij)
        # representation with rho = -ad: [rho_i, rho_j] = rho_[i,j]  <=>  [ad_i, ad_j] = -(-ad_[i,j])... use Lambda = -rho
        # degree reading: Lambda = -ad ; required: [Lambda_i, Lambda_j] = -Lambda_[i,j]
        Li, Lj, Lij = -ad[i], -ad[j], -Cij
        ok_deg &= close(br(Li, Lj), -Lij)
        Ti, Tj, Tij = F @ Li, F @ Lj, F @ Lij
        ok_tot &= close(br(Ti, Tj), -Tij)
k4["degree_reading_representation_holds"] = bool(ok_deg)
k4["total_parity_reading_representation_holds"] = bool(ok_tot)
res["K4_grading_of_Lambda1"] = k4

res["summary"] = {
    "K1": all(v for v in k1.values()),
    "K2": bool(k2["Theta_n_commutes_with_classical_dilatation"] and k2["difference_vanishes_for_real_c"]
               and k2["difference_nonzero_for_imaginary_c"] and all(s == 1 for s in signs.values())),
    "K3": bool(k3["[K,W]=K"] and k3["-pK+[K,M_twopoint] == pK+[K,M_twopoint-2pW]"]
               and not k3["-pK+[K,M_twopoint] == pK+[K,M_twopoint] (two-point Gram in the frame formula)"]),
    "K4": bool(k4["degree_reading_representation_holds"] and not k4["total_parity_reading_representation_holds"]),
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ip_completeness2_checks.json")
with open(out, "w") as fh:
    json.dump(res, fh, indent=1, default=str)
print(json.dumps(res["summary"]))
