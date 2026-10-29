#!/usr/bin/env python3
"""Audit ip_audit2 (algebra and signs), part 4.

K1  Reflection compatibility of the kernels (ch08 Thm c8-fields (ii), conjugation.tex Prop ops), position space:
    with (e^a(x))^Theta = sum_b conj(M_ab) e^b(r_n x) (anti-multiplicative conjugation of coordinate functionals),
      sum conj(M_ac M_bd) (-1)^{|c||d|} phi_cd(r_n x, r_n y) = lambda_phi conj(phi_ab(x, y)),
    lambda = +1 for P (components of ch03 Prop c3-P), lambda = -1 for K_t (components of ch01 Def c1-bracket with the
    heat-kernel factor).  Slots: A (10), A* (10), c, c*, lambda (Sigma, 16), lambda* (Sigma^*, 16); colour diagonal.
R1  Rescaling relations (inner_product_design Thm ip-P (ii)) against the representation property with a nonzero D_1:
    on the Konishi pair (K, O_xi) with rho(x) = (1 + r_x u) rho_0(x) for odd x and rho(Delta) = (1 + r_Delta u) rho_0(Delta)
    - D_1, D_1(K) = 2 k u K.
Output: a4_results.json
"""
import json, os, itertools, resource
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(777)
OUT = {}
# ---------------------------------------------------------------- realisation (Def c8-twocomp)
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
sig = [s1, s2, s3, 1j * np.eye(2)]; sigb = [-s1, -s2, -s3, 1j * np.eye(2)]
def ix(block, a, i):
    return {'xi': 0, 'xib': 8, 'eta': 16, 'etab': 24}[block] + 4 * a + i
def rho_V(v):
    M = np.zeros((32, 32), complex)
    vs = sum(v[m] * sig[m] for m in range(4)); vsb = sum(v[m] * sigb[m] for m in range(4))
    for i in range(4):
        for a in range(2):
            for ad in range(2):
                M[ix('xi', a, i), ix('etab', ad, i)] += -1j * vs[a, ad]
                M[ix('xib', ad, i), ix('eta', a, i)] += 1j * vs[a, ad]
                M[ix('eta', a, i), ix('xib', ad, i)] += 1j * vsb[ad, a]
                M[ix('etab', ad, i), ix('xi', a, i)] += -1j * vsb[ad, a]
    return M
Gp = np.diag([1.0] * 8 + [-1.0] * 8)
J = np.zeros((32, 32))           # cc(z) = J conj(z)
for i in range(8):
    J[16 + i, i] = 1; J[24 + i, 8 + i] = -1; J[i, 16 + i] = 1; J[8 + i, 24 + i] = -1

# slots: A (10, even), As (10, odd), c (1, odd), cs (1, even), lam (16, odd), ls (16, even)
slots = [('A', 10, 0), ('As', 10, 1), ('c', 1, 1), ('cs', 1, 0), ('lam', 16, 1), ('ls', 16, 0)]
off = {}; o = 0
for nm, d, p in slots:
    off[nm] = (o, d, p); o += d
NT = o
PAR = np.zeros(NT, int)
for nm, (o_, d, p) in off.items():
    PAR[o_:o_ + d] = p
def Mtheta(n):
    rn = np.eye(4) - 2 * np.outer(n, n); Rn = np.eye(10); Rn[:4, :4] = rn
    RN = rho_V(n) @ J
    M = np.zeros((NT, NT), complex)
    def put(nm, mat):
        o_, d, _ = off[nm]; M[o_:o_ + d, o_:o_ + d] = mat
    put('A', Rn); put('As', -Rn); put('c', np.eye(1)); put('cs', -np.eye(1))
    put('lam', RN[0:16, 0:16]); put('ls', -RN[16:32, 16:32])
    return M
def G(u):
    return 1.0 / (4 * np.pi ** 2 * (u @ u))
def dG(u):
    return -2 * u / (4 * np.pi ** 2 * (u @ u) ** 2)
def heat(u, t=0.37):
    return np.exp(-(u @ u) / (4 * t)) / (4 * np.pi * t) ** 2
def Pker(x, y):
    """components P_ab(x, y) of P(0,infty) (ch03 Prop c3-P), colour factor (-,-)^{-1} = 1."""
    u = x - y; P = np.zeros((NT, NT), complex)
    oA = off['A'][0]; oAs = off['As'][0]; oc = off['c'][0]; ol = off['lam'][0]
    P[oA:oA + 10, oA:oA + 10] = 2 * G(u) * np.eye(10)
    for m in range(4):
        P[oc, oAs + m] = -dG(u)[m]          # c (x) A*|_V : -(,)^{-1} d_x G
        P[oAs + m, oc] = -dG(u)[m]          # A*|_V (x) c : -(,)^{-1} d_x G (d_x on the first argument)
    Mx = sum(dG(u)[m] * rho_V(np.eye(4)[m]) for m in range(4))     # rho(d_x G): Sigma^* -> Sigma
    FD = np.zeros((16, 32))
    for s_ in range(16):
        FD[s_, 16 + s_] = Gp[s_, s_]         # dual basis f_s in Sigma^*
    for s_ in range(16):
        for t_ in range(16):
            P[ol + s_, ol + t_] = (Mx @ FD[t_])[s_]
    return P
def Kker(x, y):
    """components of K_t: each field-antifield slot pair equals -h_t(x-y) times the inverse fibre pairing."""
    h = heat(x - y); K = np.zeros((NT, NT), complex)
    oA = off['A'][0]; oAs = off['As'][0]; oc = off['c'][0]; ocs = off['cs'][0]; ol = off['lam'][0]; ols = off['ls'][0]
    K[oA:oA + 10, oAs:oAs + 10] = -h * np.eye(10); K[oAs:oAs + 10, oA:oA + 10] = -h * np.eye(10)
    K[oc, ocs] = -h; K[ocs, oc] = -h
    K[ol:ol + 16, ols:ols + 16] = -h * Gp; K[ols:ols + 16, ol:ol + 16] = -h * Gp   # canonical element of Sigma (x) Sigma^*
    return K
res = {}
for nm, ker, lam in [('P', Pker, 1.0), ('K_t', Kker, -1.0)]:
    dev = 0.0; devwrong = 0.0
    for _ in range(4):
        n = rng.normal(size=4); n /= np.linalg.norm(n); rn = np.eye(4) - 2 * np.outer(n, n)
        x = rng.normal(size=4); y = rng.normal(size=4)
        M = Mtheta(n); Kr = ker(rn @ x, rn @ y)
        sgn = (-1.0) ** np.outer(PAR, PAR)
        lhs = np.conj(M) @ (sgn * Kr) @ np.conj(M).T
        rhs = lam * np.conj(ker(x, y))
        dev = max(dev, float(np.max(np.abs(lhs - rhs))))
        devwrong = max(devwrong, float(np.max(np.abs(lhs + rhs))))
    res[nm] = {'lambda': lam, 'max_deviation': dev, 'control_opposite_lambda_deviation': devwrong}
# graded symmetry of P as printed (exchange of slots with x <-> y gives the Koszul sign)
x = rng.normal(size=4); y = rng.normal(size=4)
P1 = Pker(x, y); P2 = Pker(y, x)
res['P_graded_symmetric'] = bool(np.allclose(P1, ((-1.0) ** np.outer(PAR, PAR)) * P2.T))
OUT['K1_kernel_reflection'] = res

# ---------------------------------------------------------------- R1 rescalings vs D_1 on the Konishi pair
u, k, rxi, rxid, rD = sp.symbols('u k r_xi r_xidag r_Delta', real=True)
xi2 = sp.symbols('xi2', positive=True)
# classical values (Prop ip-tree): rho0(xi^dag) rho0(xi) K = 2 (xi, conj xi) K ; rho0(Delta) K = -4 K
lhs = (1 + rxid * u) * (1 + rxi * u) * 2 * xi2
rhs = -sp.Rational(1, 2) * xi2 * ((1 + rD * u) * (-4) - 2 * k * u)     # rho([xi^dag,xi]) = -1/2 (xi,conj xi) rho(Delta)
eq1 = sp.expand(lhs - rhs).coeff(u, 1)
sol = sp.solve(sp.Eq(eq1, 0), rxi + rxid) if False else sp.simplify(sp.solve(sp.Eq(eq1, 0), rxid)[0] + rxi)
OUT['R1_rescalings_konishi'] = {
    'order_u_relation_r_xi_plus_r_xidag': str(sol),
    'with_r_Delta_0_and_k_minus_24': str(sol.subs({rD: 0, k: -24})),
    'Thm_ip-P_(ii)_claims': 'r_xi + r_eta = 0 for xi in S_+(x)W, eta in S_+(x)W^* (xi^dag = -conj(xi) is in S_+(x)W^*)',
    'Prop_ip-konishi_T2': 'q + s = -12 for the same rescalings on (K, O_xi)'}
OUT['peak_rss_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
with open(os.path.join(HERE, 'a4_results.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
print(json.dumps(OUT, indent=1, default=str))
