#!/usr/bin/env python3
"""Independent re-derivation (audit ip_audit2, lens: algebra and signs).

Built from scratch from the printed definitions only (no import of project code):
  * ch08 Def c8-twocomp: Euclidean sigma, sigmabar, epsilon, 't Hooft blocks, rho(v), rho(w), pairing, Gamma.
  * ch01 Def c1-spaces (Clifford relation, symmetry of rho, Gamma-pairing), Prop c1-chirality.
  * ch02 Def c2-even, c2-odd, c2-bracket (eq:c2-evenodd, eq:c2-oddodd), Prop c2-special, table of blocks.
  * ch08 Lemma c8-cc, Def c8-dag, Thm c8-dag, Def c8-Theta (dagger_n), Prop c8-C, Def c8-B.
Numerics in double precision; every structure constant is rounded to the lattice (1/8)(Z + iZ) and the
rounding error is recorded, after which the bracket checks are exact integer checks.
Output: a1_results.json next to this file.
"""
import json, os, itertools
import numpy as np
from scipy.linalg import expm

HERE = os.path.dirname(os.path.abspath(__file__))
R = {}
rng = np.random.default_rng(20261028)

# ---------------------------------------------------------------- two-component data (Def c8-twocomp)
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)
sig = [s1, s2, s3, 1j * I2]          # sigma^m_{alpha alphadot}  (row alpha, column alphadot)
sigb = [-s1, -s2, -s3, 1j * I2]      # sigmabar^{m alphadot alpha} (row alphadot, column alpha)
epsU = np.array([[0, 1], [-1, 0]], complex)   # eps^{12} = +1 (both undotted and dotted)
epsL = np.array([[0, -1], [1, 0]], complex)   # eps_{12} = -1
eps3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    eps3[p] = np.linalg.det(np.eye(3)[list(p)])
eps4 = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(4)):
    eps4[p] = np.linalg.det(np.eye(4)[list(p)])
Sig = np.zeros((6, 4, 4), complex)   # Sigma^I_{ij}
for a in range(3):
    for b in range(3):
        for c in range(3):
            Sig[a, b, c] = eps3[a, b, c]
            Sig[a + 3, b, c] = 1j * eps3[a, b, c]
        Sig[a, b, 3] = 1.0 * (a == b); Sig[a, 3, b] = -1.0 * (a == b)
        Sig[a + 3, b, 3] = -1j * (a == b); Sig[a + 3, 3, b] = 1j * (a == b)

def up(wl):            # w^{ij} = 1/2 eps^{ijkl} w_{kl}
    return 0.5 * np.einsum('ijkl,kl->ij', eps4, wl)

# 32-component layout: Sigma = [xi_alpha^i (0..7), xibar_{alphadot i} (8..15)],
#                      Sigma^* = [eta^alpha_i (16..23), etabar^{alphadot i} (24..31)]; index = base + 4*alpha + i
def ix(block, a, i):
    return {'xi': 0, 'xib': 8, 'eta': 16, 'etab': 24}[block] + 4 * a + i

def rho_V(v):
    """rho(v), v in V (components v_m, m = 1..4 -> 0..3), as a 32x32 matrix (Def c8-twocomp)."""
    M = np.zeros((32, 32), complex)
    vs = sum(v[m] * sig[m] for m in range(4)); vsb = sum(v[m] * sigb[m] for m in range(4))
    for i in range(4):
        for a in range(2):
            for ad in range(2):
                # on eta in Sigma^*: (rho(v)eta)_alpha^i = -i v_m sigma^m_{alpha alphadot} etabar^{alphadot i}
                M[ix('xi', a, i), ix('etab', ad, i)] += -1j * vs[a, ad]
                # (rho(v)eta)_{alphadot i} = i v_m eta^alpha_i sigma^m_{alpha alphadot}
                M[ix('xib', ad, i), ix('eta', a, i)] += 1j * vs[a, ad]
                # on xi in Sigma: (rho(v)xi)^alpha_i = i v_m xibar_{alphadot i} sigmabar^{m alphadot alpha}
                M[ix('eta', a, i), ix('xib', ad, i)] += 1j * vsb[ad, a]
                # (rho(v)xi)^{alphadot i} = -i v_m sigmabar^{m alphadot alpha} xi_alpha^i
                M[ix('etab', ad, i), ix('xi', a, i)] += -1j * vsb[ad, a]
    return M

def rho_W(w):
    """rho(w), w in W_R^perp with components w_I (I = 1..6 -> 0..5) (Def c8-twocomp)."""
    M = np.zeros((32, 32), complex)
    wl = np.einsum('I,Iij->ij', w, Sig); wu = up(wl)
    for a in range(2):
        for b in range(2):
            for j in range(4):
                for i in range(4):
                    # (rho(w)xi)^alpha_j = eps^{alpha beta} w_{ji} xi_beta^i
                    M[ix('eta', a, j), ix('xi', b, i)] += epsU[a, b] * wl[j, i]
                    # (rho(w)xi)^{alphadot j} = eps^{alphadot betadot} w^{ji} xibar_{betadot i}
                    M[ix('etab', a, j), ix('xib', b, i)] += epsU[a, b] * wu[j, i]
                    # (rho(w)eta)_alpha^k = -w^{kj} eps_{alpha beta} eta^beta_j   (k -> j, j -> i here)
                    M[ix('xi', a, j), ix('eta', b, i)] += -wu[j, i] * epsL[a, b]
                    # (rho(w)eta)_{betadot k} = -w_{kj} eps_{betadot alphadot} etabar^{alphadot j}
                    M[ix('xib', a, j), ix('etab', b, i)] += -wl[j, i] * epsL[a, b]
    return M

def rho10(u):          # u in V_10 = V (+) W_R^perp, components u[0:4] (V), u[4:10] (W)
    return rho_V(u[:4]) + rho_W(u[4:])

E10 = np.eye(10)
RHO = [rho10(E10[M]) for M in range(10)]
# pairing (xi, eta) = xi_alpha^i eta^alpha_i - xibar_{alphadot i} etabar^{alphadot i}; symmetric extension B
Gp = np.zeros((16, 16)); Gp[:8, :8] = np.eye(8); Gp[8:, 8:] = -np.eye(8)
Bform = np.zeros((32, 32)); Bform[:16, 16:] = Gp; Bform[16:, :16] = Gp.T
def pair(z1, z2):      # symmetric bilinear extension of Sigma (x) Sigma^* -> C
    return z1 @ Bform @ z2
def inS(z):
    return np.allclose(z[16:], 0)
def inSs(z):
    return np.allclose(z[:16], 0)
def gamma(z1, z2):     # (v, Gamma(z1, z2)) = (rho(v) z1, z2) for z1, z2 both in Sigma or both in Sigma^*
    return np.array([pair(RHO[M] @ z1, z2) for M in range(10)])

def rvec(block=None):
    z = rng.normal(size=32) + 1j * rng.normal(size=32)
    if block == 'S': z[16:] = 0
    if block == 'Ss': z[:16] = 0
    return z

# ---------------------------------------------------------------- 1. realisation identities (Prop c8-twocomp)
res = {}
res['clifford'] = all(np.allclose(RHO[M] @ RHO[N] + RHO[N] @ RHO[M], 2 * (M == N) * np.eye(32)) for M in range(10) for N in range(10))
sym = True
for M in range(10):
    for _ in range(3):
        for blk in ('S', 'Ss'):
            z1, z2 = rvec(blk), rvec(blk)
            sym &= np.isclose(pair(RHO[M] @ z1, z2), pair(z1, RHO[M] @ z2))
res['rho_symmetric_for_pairing'] = bool(sym)
chir = RHO[0] @ RHO[1] @ RHO[2] @ RHO[3]
res['rho_e1234_on_SplusW_is_plus1'] = bool(np.allclose(chir[0:8, 0:8], np.eye(8)) and np.allclose(chir[0:8, 8:], 0))
res['rho_e1234_on_SminusWstar_is_minus1'] = bool(np.allclose(chir[8:16, 8:16], -np.eye(8)))
res['rho_e1234_on_SplusWstar_is_plus1'] = bool(np.allclose(chir[16:24, 16:24], np.eye(8)))
res['rho_e1234_on_SminusW_is_minus1'] = bool(np.allclose(chir[24:32, 24:32], -np.eye(8)))
# ESW Assumption 3.1 for 10d (Gamma(Q,Q).Q = 0), checked as rho(Gamma(Q,Q))Q = 0
okA = True
for _ in range(5):
    Q = rvec('S'); G = gamma(Q, Q)
    okA &= np.allclose(sum(G[M] * RHO[M] for M in range(10)) @ Q, 0)
    Q2 = rvec('Ss'); G2 = gamma(Q2, Q2)
    okA &= np.allclose(sum(G2[M] * RHO[M] for M in range(10)) @ Q2, 0)
res['ESW_assumption_3_1_rho_Gamma_QQ_Q_zero'] = bool(okA)
# 't Hooft block identities
res['thooft_sum_2eps'] = bool(np.allclose(np.einsum('Iij,Ikl->ijkl', Sig, Sig), 2 * eps4))
res['thooft_conj_equals_half_eps'] = bool(np.allclose(np.conj(Sig), 0.5 * np.einsum('ijkl,Ikl->Iij', eps4, Sig)))
res['thooft_orthogonality_4delta'] = bool(np.allclose(np.einsum('Iij,Jij->IJ', np.conj(Sig), Sig), 4 * np.eye(6)))

# Gamma component formulas of Def c8-twocomp
def comps(z):
    xi = z[0:8].reshape(2, 4); xib = z[8:16].reshape(2, 4); eta = z[16:24].reshape(2, 4); etab = z[24:32].reshape(2, 4)
    return xi, xib, eta, etab
okG = True
for _ in range(4):
    z1, z2 = rvec('S'), rvec('S')
    G = gamma(z1, z2)
    xi1, xib1, _, _ = comps(z1); xi2, xib2, _, _ = comps(z2)
    GV = np.array([1j * (np.einsum('ai,ab,bi->', xib2, sigb[m], xi1) + np.einsum('ai,ab,bi->', xib1, sigb[m], xi2)) for m in range(4)])
    okG &= np.allclose(G[:4], GV)
    Gkl = 2 * np.einsum('jikl,aj,ab,bi->kl', eps4, xi1, epsU, xi2) - 2 * (np.einsum('ak,ab,bl->kl', xib1, epsU, xib2) - np.einsum('al,ab,bk->kl', xib1, epsU, xib2))
    okG &= np.allclose(np.einsum('I,Ikl->kl', G[4:], Sig), Gkl)
    e1, e2 = rvec('Ss'), rvec('Ss')
    G = gamma(e1, e2)
    _, _, eta1, etab1 = comps(e1); _, _, eta2, etab2 = comps(e2)
    GV = np.array([-1j * (np.einsum('ai,ab,bi->', eta2, sig[m], etab1) + np.einsum('ai,ab,bi->', eta1, sig[m], etab2)) for m in range(4)])
    okG &= np.allclose(G[:4], GV)
    Gu = -2 * np.einsum('klpj,ab,bj,ap->kl', eps4, epsL, eta1, eta2) + 2 * (np.einsum('ba,al,bk->kl', epsL, etab1, etab2) - np.einsum('ba,ak,bl->kl', epsL, etab1, etab2))
    okG &= np.allclose(up(np.einsum('I,Ikl->kl', G[4:], Sig)), Gu)
res['Gamma_component_formulas'] = bool(okG)
# (Gamma,u) = Gamma^m u_m + 1/4 Gamma_ij u^ij
okP = True
for _ in range(3):
    z1, z2 = rvec('S'), rvec('S'); u = rng.normal(size=10)
    G = gamma(z1, z2)
    lhs = G @ u
    rhs = G[:4] @ u[:4] + 0.25 * np.einsum('ij,ij->', np.einsum('I,Iij->ij', G[4:], Sig), up(np.einsum('I,Iij->ij', u[4:], Sig)))
    okP &= np.isclose(lhs, rhs)
res['Gamma_u_pairing_formula'] = bool(okP)
R['1_realisation'] = res

# ---------------------------------------------------------------- 2. spinor conjugation (Lemma c8-cc)
def cc(z):
    """two-component conjugation of eq:c8-cc: Sigma -> Sigma^*, Sigma^* -> Sigma (antilinear)."""
    out = np.zeros(32, complex)
    out[16:24] = np.conj(z[0:8]); out[24:32] = -np.conj(z[8:16])     # xi -> (eta = conj xi, etab = -conj xib)
    out[0:8] = np.conj(z[16:24]); out[8:16] = -np.conj(z[24:32])     # eta -> (xi = conj eta, xib = -conj etab)
    return out
res = {}
res['involutive'] = bool(all(np.allclose(cc(cc(z)), z) for z in [rvec() for _ in range(4)]))
okI = True
for M in range(10):
    for _ in range(2):
        z = rvec(); okI &= np.allclose(cc(RHO[M] @ z), RHO[M] @ cc(z))   # e_M real, so bar(e_M) = e_M
for _ in range(3):
    u = rng.normal(size=10) + 1j * rng.normal(size=10); z = rvec()
    okI &= np.allclose(cc(rho10(u) @ z), rho10(np.conj(u)) @ cc(z))
res['intertwines_rho_vbar'] = bool(okI)
okpos = all((lambda z: np.isclose(pair(z, cc(z)).imag, 0) and pair(z, cc(z)).real > 0)(rvec('S')) for _ in range(5))
okpos2 = all((lambda z: np.isclose(pair(cc(z), z).imag, 0) and pair(cc(z), z).real > 0)(rvec('Ss')) for _ in range(5))
res['positive_on_Sigma_and_SigmaStar'] = bool(okpos and okpos2)
okG = True
for blk in ('S', 'Ss'):
    for _ in range(3):
        z3, z4 = rvec(blk), rvec(blk)
        okG &= np.allclose(gamma(cc(z3), cc(z4)), np.conj(gamma(z3, z4)))
res['Gamma_of_conjugates_is_conj_Gamma'] = bool(okG)
# uniqueness: antilinear maps J o conj with J linear Sigma<->Sigma^*, commuting with all rho(e_M) (e_M real)
# condition: J conj(rho_M) = rho_M J  -> linear equations for J (32x32, off-diagonal blocks)
# restrict J to blocks Sigma->Sigma^* and Sigma^*->Sigma; accumulate the normal matrix (memory-light)
mask = np.zeros((32, 32), bool); mask[16:, :16] = True; mask[:16, 16:] = True
cols = np.where(mask.flatten(order='F'))[0]
Gn = np.zeros((len(cols), len(cols)), complex)
for M in range(10):
    A = np.conj(RHO[M]); Bm = RHO[M]
    # vec(J A - B J) = (A^T kron I - I kron B) vec(J)   (column-major vec)
    Kr = (np.kron(A.T, np.eye(32)) - np.kron(np.eye(32), Bm))[:, cols]
    Gn += Kr.conj().T @ Kr
ev = np.linalg.eigvalsh(Gn)
res['dim_equivariant_antilinear_maps'] = int(np.sum(ev < 1e-9))
R['2_spinor_conjugation'] = res

# ---------------------------------------------------------------- 3. psl(4|4) of Chapter 2
# even V-part: vector field v(x) = a + B x + C(x,x)  (C[k,l,m] symmetric in l,m)
def vf_eval(vf, x):
    a, B, C = vf
    return a + B @ x + np.einsum('klm,l,m->k', C, x, x)
def vf_jac(vf, x):     # J[k,l] = d v^k / d x^l
    a, B, C = vf
    return B + 2 * np.einsum('klm,m->kl', C, x)
def vfP(b):
    return (np.array(b, complex), np.zeros((4, 4), complex), np.zeros((4, 4, 4), complex))
def vfRot(m, n):
    B = np.zeros((4, 4), complex); B[m, n] = 1; B[n, m] = -1     # v = e_m x_n - e_n x_m
    return (np.zeros(4, complex), B, np.zeros((4, 4, 4), complex))
def vfDil():
    return (np.zeros(4, complex), 2 * np.eye(4, dtype=complex), np.zeros((4, 4, 4), complex))
def vfK(b):
    b = np.array(b, complex); C = np.zeros((4, 4, 4), complex)
    for k in range(4):
        for l in range(4):
            for m in range(4):
                C[k, l, m] = b[k] * (l == m) - 0.5 * (b[l] * (k == m) + b[m] * (k == l)) * 2 / 2 * 1.0
    # x^2 b^k - 2 (b.x) x^k : symmetric C[k,l,m] = b^k delta_lm - (b_l delta_km + b_m delta_kl)
    for k in range(4):
        for l in range(4):
            for m in range(4):
                C[k, l, m] = b[k] * (l == m) - (b[l] * (k == m) + b[m] * (k == l))
    return (np.zeros(4, complex), np.zeros((4, 4), complex), C)
# check vfK
_x = rng.normal(size=4); _b = rng.normal(size=4)
assert np.allclose(vf_eval(vfK(_b), _x), (_x @ _x) * _b - 2 * (_b @ _x) * _x)

basis = []      # list of dicts
names = []
for m in range(4):
    basis.append({'t': 'V', 'vf': vfP(np.eye(4)[m])}); names.append('P%d' % (m + 1))
rot_pairs = [(m, n) for m in range(4) for n in range(m + 1, 4)]
for (m, n) in rot_pairs:
    basis.append({'t': 'V', 'vf': vfRot(m, n)}); names.append('M%d%d' % (m + 1, n + 1))
basis.append({'t': 'V', 'vf': vfDil()}); names.append('Delta')
for m in range(4):
    basis.append({'t': 'V', 'vf': vfK(np.eye(4)[m])}); names.append('K%d' % (m + 1))
w_pairs = [(I, J) for I in range(6) for J in range(I + 1, 6)]
for (I, J) in w_pairs:
    Rm = np.zeros((6, 6), complex); Rm[I, J] = 1; Rm[J, I] = -1
    basis.append({'t': 'W', 'R': Rm}); names.append('R%d%d' % (I + 1, J + 1))
for k in range(16):      # Q1 in Sigma (constant spinor)
    z = np.zeros(32, complex); z[k] = 1
    basis.append({'t': 'O', 'Q1': z, 'Q2': np.zeros(32, complex)}); names.append('Sig%d' % k)
for k in range(16):      # Q2 in Sigma^* (x -> rho(x) Q2)
    z = np.zeros(32, complex); z[16 + k] = 1
    basis.append({'t': 'O', 'Q1': np.zeros(32, complex), 'Q2': z}); names.append('SigS%d' % k)
NB = len(basis); assert NB == 62
EVEN = list(range(30)); ODD = list(range(30, 62))
def parity(a): return 0 if a < 30 else 1
def idx_xi(a, i): return 30 + 4 * a + i            # Q^i_alpha
def idx_xib(a, i): return 30 + 8 + 4 * a + i       # Qbar_{alphadot i}
def idx_eta(a, i): return 46 + 4 * a + i           # S^alpha_i
def idx_etab(a, i): return 46 + 8 + 4 * a + i      # Sbar^{alphadot i}

PTS = [rng.normal(size=4) for _ in range(40)]
RV_PTS = [rho_V(x) for x in PTS]
_RVcache = {}
def rhoV_at(x):
    for k, y in enumerate(PTS):
        if y is x:
            return RV_PTS[k]
    return rho_V(x)
# least-squares decomposition matrices
def vf_vec(vf):
    return np.concatenate([vf_eval(vf, x) for x in PTS])
VB = np.stack([vf_vec(basis[a]['vf']) for a in range(15)], axis=1)
VB_pinv = np.linalg.pinv(VB)
def decomp_vf_values(vals):        # vals: concatenated values at PTS
    c = VB_pinv @ vals
    return c, np.max(np.abs(VB @ c - vals))
def eps_of(el, x):                  # conformal Killing spinor value at x (in Sigma)
    return el['Q1'] + rhoV_at(x) @ el['Q2']
_Asp = []
for RV in RV_PTS:
    blk = np.zeros((32, 32), complex)
    blk[:, :16] = np.eye(32)[:, :16]
    blk[:, 16:] = RV[:, 16:]
    _Asp.append(blk)
_Asp = np.vstack(_Asp)
_Asp_pinv = np.linalg.pinv(_Asp)
def decomp_spinor_values(vals_at_pts):
    # unknowns Q1 (Sigma part, 16) and Q2 (Sigma^* part, 16): eps(x) = Q1 + rho(x) Q2
    rhs = np.concatenate(vals_at_pts)
    c = _Asp_pinv @ rhs
    return c, np.max(np.abs(_Asp @ c - rhs))

def elem_to_geom(c):
    """coefficient vector -> geometric data (vf, R, Q1, Q2)"""
    a = sum(c[k] * basis[k]['vf'][0] for k in range(15)); B = sum(c[k] * basis[k]['vf'][1] for k in range(15))
    C = sum(c[k] * basis[k]['vf'][2] for k in range(15))
    Rm = sum(c[k] * basis[k]['R'] for k in range(15, 30))
    Q1 = sum(c[k] * basis[k]['Q1'] for k in range(30, 62)); Q2 = sum(c[k] * basis[k]['Q2'] for k in range(30, 62))
    return (a, B, C), Rm, Q1, Q2

def bracket_basis(p, q):
    """[e_p, e_q] in the realisation of Chapter 2, returned as coefficient vector and fit residual."""
    A, Bq = basis[p], basis[q]
    out = np.zeros(NB, complex); resid = 0.0
    if A['t'] in 'VW' and Bq['t'] in 'VW':
        if A['t'] == 'V' and Bq['t'] == 'V':
            vals = np.concatenate([vf_jac(Bq['vf'], x) @ vf_eval(A['vf'], x) - vf_jac(A['vf'], x) @ vf_eval(Bq['vf'], x) for x in PTS])
            c, resid = decomp_vf_values(vals); out[:15] = c
        elif A['t'] == 'W' and Bq['t'] == 'W':
            Rn = -(A['R'] @ Bq['R'] - Bq['R'] @ A['R'])       # [v_R, v_R'] = v_{R'R - RR'}
            for k, (I, J) in enumerate(w_pairs):
                out[15 + k] = Rn[I, J]
            resid = np.max(np.abs(Rn + Rn.T))
        return out, resid
    if A['t'] == 'O' and Bq['t'] in 'VW':
        o, r = bracket_basis(q, p)
        return -o, r                       # [odd, even] = -[even, odd]
    if A['t'] in 'VW' and Bq['t'] == 'O':
        vals = []
        for x in PTS:
            e = eps_of(Bq, x)
            if A['t'] == 'V':
                v = vf_eval(A['vf'], x); J = vf_jac(A['vf'], x)
                deriv = sum(v[k] * (RHO[k] @ Bq['Q2']) for k in range(4))   # v(eps) = v^k rho(e_k) Q2
                om = J - J.T            # om[k,l] = dv^k/dx^l - dv^l/dx^k ; 2-form coefficient omega_{lk} = J[k,l]-J[l,k]
                q2 = np.zeros((32, 32), complex)
                for l in range(4):
                    for k in range(l + 1, 4):
                        q2 += (J[k, l] - J[l, k]) * (RHO[l] @ RHO[k])
                div = np.trace(J)
                vals.append(deriv + 0.25 * q2 @ e - 0.125 * div * e)
            else:
                Rm = A['R']; q2 = np.zeros((32, 32), complex)
                for I in range(6):
                    for J in range(I + 1, 6):
                        q2 += -2 * Rm[I, J] * (RHO[4 + I] @ RHO[4 + J])
                vals.append(0.25 * q2 @ e)
        c, resid = decomp_spinor_values(vals)
        out[30:46] = c[:16]; out[46:62] = c[16:]
        return out, resid
    # odd-odd
    vals = np.concatenate([gamma(eps_of(A, x), eps_of(Bq, x))[:4] for x in PTS])
    c, resid = decomp_vf_values(vals); out[:15] = c
    Q1, Q2, Q3, Q4 = A['Q1'], A['Q2'], Bq['Q1'], Bq['Q2']
    Rm = np.zeros((6, 6), complex)
    for I in range(6):
        for J in range(6):
            if I != J:
                Mq = RHO[4 + I] @ RHO[4 + J]
                Rm[I, J] = pair(Mq @ Q1, Q4) + pair(Mq @ Q3, Q2)
    for k, (I, J) in enumerate(w_pairs):
        out[15 + k] = Rm[I, J]
    resid = max(resid, np.max(np.abs(Rm + Rm.T)))
    return out, resid

F = np.zeros((NB, NB, NB), complex)       # F[p,q,:] = [e_p, e_q]
maxres = 0.0
for p in range(NB):
    for q in range(NB):
        o, r = bracket_basis(p, q); F[p, q] = o; maxres = max(maxres, r)
Fr = np.round(F * 8) / 8
round_err = float(np.max(np.abs(F - Fr)))
F = Fr
res = {'fit_residual_max': float(maxres), 'rounding_to_eighths_max_error': round_err}
# graded antisymmetry
gas = True
for p in range(NB):
    for q in range(NB):
        sgn = -((-1) ** (parity(p) * parity(q)))
        gas &= np.allclose(F[p, q], sgn * F[q, p])
res['graded_antisymmetry'] = bool(gas)
# super-Jacobi (graded): (-1)^{ac}[x,[y,z]] + cyclic = 0 ; exact on lattice values
def br(c1, c2):
    return np.einsum('p,q,pqr->r', c1, c2, F)
Fi = np.round(F * 8).astype(np.complex128)      # 8F integral
jac_max = 0.0
par = np.array([parity(a) for a in range(NB)])
for a in range(NB):
    # T[b,c,:] = [e_a,[e_b,e_c]] = sum_d F[b,c,d] F[a,d,:]
    T1 = np.einsum('bcd,dr->bcr', F, F[a])
    # [e_b,[e_c,e_a]] and [e_c,[e_a,e_b]]
    T2 = np.einsum('cd,bdr->bcr', F[:, a, :], F)          # [e_c,e_a] = F[c,a,:]; [e_b, .] = F[b,d,:]
    T3 = np.einsum('bd,cdr->bcr', F[a, :, :], F)          # [e_a,e_b] = F[a,b,:]; [e_c, .] = F[c,d,:]
    s_ac = (-1.0) ** (par[a] * par)[None, :]              # (-1)^{|a||c|} over c
    s_ba = (-1.0) ** (par * par[a])[:, None]              # (-1)^{|b||a|} over b
    s_cb = (-1.0) ** np.outer(par, par).T                 # (-1)^{|c||b|}
    J = s_ac[:, :, None] * T1 + s_ba[:, :, None] * T2 + s_cb[:, :, None] * T3
    jac_max = max(jac_max, float(np.max(np.abs(J))))
res['super_jacobi_max_violation'] = jac_max
# brackets of named generators (Prop c2-special), ad_Delta eigenvalues
iD = names.index('Delta')
adD = F[iD].T            # (adD)[r, q] = coefficient r of [Delta, e_q]
evals = {}
for nm, rng_ in [('P', range(0, 4)), ('M', range(4, 10)), ('K', range(11, 15)), ('R', range(15, 30)), ('Sigma', range(30, 46)), ('SigmaStar', range(46, 62))]:
    vals = set()
    ok = True
    for q in rng_:
        col = adD[:, q]; lam = col[q]
        ok &= np.allclose(col - lam * np.eye(NB)[q], 0)
        vals.add(complex(lam))
    evals[nm] = [str(v) for v in vals] + ([] if ok else ['NOT DIAGONAL'])
res['ad_Delta_eigenvalues'] = evals
# [K_b, P_b'] = (b,b') Delta + rotation
okKP = True
for m in range(4):
    for n in range(4):
        c = F[11 + m, n]
        okKP &= np.isclose(c[iD], 1.0 * (m == n))
        okKP &= np.allclose(np.delete(c, list(range(4, 11))), 0) or np.allclose(c[[0, 1, 2, 3] + list(range(11, 62))], 0)
res['K_P_bracket_dilatation_part_is_(b,bprime)'] = bool(okKP)
res['Delta_P_bracket'] = 'ad_Delta on P = ' + evals['P'][0]
# dilatation part of [xi, eta] = 1/2 (xi, eta) Delta
okXE = True
for _ in range(4):
    c1 = np.zeros(NB, complex); c2 = np.zeros(NB, complex)
    c1[30:46] = rng.normal(size=16) + 1j * rng.normal(size=16); c2[46:62] = rng.normal(size=16) + 1j * rng.normal(size=16)
    _, _, Q1, _ = elem_to_geom(c1); _, _, _, Q2 = elem_to_geom(c2)
    okXE &= np.isclose(br(c1, c2)[iD], 0.5 * pair(Q1, Q2))
res['dilatation_part_of_[xi,eta]_is_half_pairing'] = bool(okXE)
R['3_psl44'] = res

# ---------------------------------------------------------------- 4. radial conjugation (Def c8-dag)
def Istar_vf_values(vf, x):      # (I_* u)(x) = (x,x) u(I(x)) - 2 (u(I(x)), x) x
    y = x / (x @ x); u = vf_eval(vf, y)
    return (x @ x) * u - 2 * (u @ x) * x
def dagger_basis(a):
    """X^dagger for a basis element (all basis elements have real geometric data except the odd ones,
    whose conjugation is the two-component map cc)."""
    out = np.zeros(NB, complex); E = basis[a]
    if E['t'] == 'V':
        vf = tuple(np.conj(t) for t in E['vf'])
        vals = np.concatenate([-Istar_vf_values(vf, x) for x in PTS])
        c, r = decomp_vf_values(vals); out[:15] = c
        return out, r
    if E['t'] == 'W':
        Rn = -np.conj(E['R'])
        for k, (I, J) in enumerate(w_pairs):
            out[15 + k] = Rn[I, J]
        return out, 0.0
    Q1n = -cc(E['Q2']); Q2n = -cc(E['Q1'])        # (Q1 + rho(x)Q2)^dagger = -conj Q2 - rho(x) conj Q1
    out[30:46] = Q1n[:16]; out[46:62] = Q2n[16:]
    # check against the defining formula eps^dagger(x) = -rho(x) conj(eps(I(x)))
    r = 0.0
    for x in PTS[:6]:
        lhs = Q1n + rho_V(x) @ Q2n
        rhs = -rho_V(x) @ cc(eps_of(E, x / (x @ x)))
        r = max(r, np.max(np.abs(lhs - rhs)))
    return out, r
T = np.zeros((NB, NB), complex); rmax = 0.0
for a in range(NB):
    col, r = dagger_basis(a); T[:, a] = col; rmax = max(rmax, r)
T = np.round(T * 8) / 8
def dag(c):         # X^dagger for a coefficient vector (antilinear)
    return T @ np.conj(c)
res = {'fit_residual': float(rmax)}
res['involution_TconjT_eq_1'] = bool(np.allclose(T @ np.conj(T), np.eye(NB)))
# anti-automorphism on all ordered basis pairs: [X,Y]^dagger = [Y^dagger, X^dagger]
amax = 0.0
for p in range(NB):
    for q in range(NB):
        lhs = dag(F[p, q]); rhs = br(T[:, q], T[:, p])
        amax = max(amax, float(np.max(np.abs(lhs - rhs))))
res['anti_automorphism_max_violation'] = amax
# control: wrong odd phase c = i (X -> i X^dagger on odd) must fail; also +I_* on even must fail
Tc = T.copy(); Tc[:, ODD] *= 1j
amax_c = max(float(np.max(np.abs(Tc @ np.conj(F[p, q]) - np.einsum('p,q,pqr->r', Tc[:, q], Tc[:, p], F)))) for p in ODD for q in ODD)
res['control_odd_phase_i_fails'] = bool(amax_c > 1e-6)
# named generators
named = {}
ok = True
for a_ in range(2):
    for i in range(4):
        ok &= np.allclose(T[:, idx_xi(a_, i)], -np.eye(NB)[idx_eta(a_, i)])        # (Q^i_alpha)^dag = -S^alpha_i
        ok &= np.allclose(T[:, idx_xib(a_, i)], np.eye(NB)[idx_etab(a_, i)])       # (Qbar)^dag = +Sbar
        ok &= np.allclose(T[:, idx_eta(a_, i)], -np.eye(NB)[idx_xi(a_, i)])        # (S)^dag = -Q
        ok &= np.allclose(T[:, idx_etab(a_, i)], np.eye(NB)[idx_xib(a_, i)])       # (Sbar)^dag = +Qbar
named['Q_dag_eq_minus_S_Qbar_dag_eq_Sbar_S_dag_eq_minus_Q_Sbar_dag_eq_Qbar'] = bool(ok)
okE = True
for m in range(4):
    okE &= np.allclose(T[:, m], -np.eye(NB)[11 + m])      # P_b^dag = -K_b (b real)
    okE &= np.allclose(T[:, 11 + m], -np.eye(NB)[m])      # K_b^dag = -P_b
okE_D = np.allclose(T[:, iD], np.eye(NB)[iD])
okE_rot = all(np.allclose(T[:, k], -np.eye(NB)[k]) for k in range(4, 10)) and all(np.allclose(T[:, k], -np.eye(NB)[k]) for k in range(15, 30))
named['P_dag_eq_minus_K_and_K_dag_eq_minus_P'] = bool(okE)
named['Delta_dag_eq_Delta'] = bool(okE_D)
named['real_rotations_dag_eq_minus_themselves'] = bool(okE_rot)
# complex translation parameter: P_b^dag = -K_{bbar}
bcpx = rng.normal(size=4) + 1j * rng.normal(size=4)
cP = np.zeros(NB, complex); cP[:4] = bcpx
cK = np.zeros(NB, complex); cK[11:15] = np.conj(bcpx)
named['P_b_dag_eq_minus_K_bbar_complex_b'] = bool(np.allclose(dag(cP), -cK))
# two-component form of xi^dagger: (xi^dag)^alpha_i = -conj(xi_alpha^i), (xi^dag)^{alphadot i} = +conj(xibar_{alphadot i})
cx = np.zeros(NB, complex); cx[30:46] = rng.normal(size=16) + 1j * rng.normal(size=16)
d = dag(cx)
named['two_component_xi_dagger'] = bool(np.allclose(d[46:54], -np.conj(cx[30:38])) and np.allclose(d[54:62], np.conj(cx[38:46])) and np.allclose(d[30:46], 0))
ce = np.zeros(NB, complex); ce[46:62] = rng.normal(size=16) + 1j * rng.normal(size=16)
d = dag(ce)
named['two_component_eta_dagger'] = bool(np.allclose(d[30:38], -np.conj(ce[46:54])) and np.allclose(d[38:46], np.conj(ce[54:62])))
# su(4): conjugation of R^i_j  = 1/4 R_IJ (Sigmabar^J Sigma^I)^i_j is the matrix adjoint
res['named'] = named
# dilatation part of [xi, xi^dagger] = -1/2 (xi, conj xi) Delta for all four odd blocks
okd = True; detail = {}
for nm, sl in [('S+W', slice(30, 38)), ('S-W*', slice(38, 46)), ('S-W', slice(54, 62)), ('S+W*', slice(46, 54))]:
    c = np.zeros(NB, complex); c[sl] = rng.normal(size=8) + 1j * rng.normal(size=8)
    g = elem_to_geom(c)
    zz = g[2] if sl.start < 46 else g[3]
    val = br(c, dag(c))[iD]
    norm = pair(zz, cc(zz)) if sl.start < 46 else pair(cc(zz), zz)
    detail[nm] = str(np.round(val / norm, 12))
    okd &= np.isclose(val, -0.5 * norm)
res['Delta_part_of_[X,X_dagger]_over_(X,conjX)'] = detail
res['Delta_part_minus_half_on_all_odd_blocks'] = bool(okd)
# real form: phi(X) = -X^dag (even), -i X^dag (odd) is an antilinear automorphism; even fixed-point Killing signature
Tphi = -T.copy(); Tphi[:, ODD] = -1j * T[:, ODD]
phimax = max(float(np.max(np.abs(Tphi @ np.conj(F[p, q]) - np.einsum('p,q,pqr->r', Tphi[:, p], Tphi[:, q], F)))) for p in range(NB) for q in range(NB))
res['phi_is_automorphism_max_violation'] = phimax
res['phi_involutive'] = bool(np.allclose(Tphi @ np.conj(Tphi), np.eye(NB)))
# real basis of fixed points in the conformal ideal (15) and in so(W) (15): X + phi(X), i(X - phi(X))
def fixed_real_basis(idxs):
    vecs = []
    for k in idxs:
        e = np.eye(NB)[k].astype(complex)
        for v in (e + Tphi @ np.conj(e), 1j * (e - Tphi @ np.conj(e))):
            vecs.append(v)
    Mv = np.array(vecs).T
    # orthonormalise over R: treat as real vectors of length 2 NB
    Rv = np.vstack([Mv.real, Mv.imag])
    U, S_, Vt = np.linalg.svd(Rv, full_matrices=False)
    keep = S_ > 1e-9
    Ur = U[:, keep]
    return Ur[:NB] + 1j * Ur[NB:]
def killing(c1, c2, idxs):
    A1 = np.einsum('p,pqr->rq', c1, F)[np.ix_(idxs, idxs)]
    A2 = np.einsum('p,pqr->rq', c2, F)[np.ix_(idxs, idxs)]
    return np.trace(A1 @ A2)
sigs = {}
for nm, idxs in [('conformal_ideal', list(range(15))), ('soW_ideal', list(range(15, 30)))]:
    Bv = fixed_real_basis(idxs)
    n = Bv.shape[1]
    Kf = np.array([[killing(Bv[:, i], Bv[:, j], idxs) for j in range(n)] for i in range(n)])
    ev = np.linalg.eigvalsh((Kf + Kf.conj().T) / 2)
    sigs[nm] = {'dim': int(n), 'negative': int(np.sum(ev < -1e-9)), 'positive': int(np.sum(ev > 1e-9)),
                'imag_part_max': float(np.max(np.abs(Kf.imag)))}
res['fixed_real_form_even_killing_signature'] = sigs
# Thm c8-dag (iii): antilinear anti-automorphisms equal to dagger on the even part: solve the even-odd condition
# X -> Tn conj(c) with Tn = T on even, Mo on odd: Mo conj(ad_v|odd) = -ad_{v^dag}|odd Mo for all even basis v
Gm = np.zeros((1024, 1024), complex)
for v in EVEN:
    Av = np.conj(F[v][np.ix_(ODD, ODD)].T)              # conj of matrix of ad_v on odd (columns = inputs)
    vd = T[:, v]
    Bv = np.einsum('p,pqr->rq', vd, F)[np.ix_(ODD, ODD)]
    Kr = np.kron(Av.T, np.eye(32)) + np.kron(np.eye(32), Bv)
    Gm += Kr.conj().T @ Kr
ev = np.linalg.eigvalsh(Gm)
res['dim_odd_maps_compatible_with_even_part'] = int(np.sum(ev < 1e-8))
# B (Def c8-B) and the family e^{tB} (+-1)^{|X|} dagger
Bd = np.zeros(NB)
Bd[30:38] = 1; Bd[38:46] = -1; Bd[46:54] = -1; Bd[54:62] = 1
Bm_ = np.diag(Bd)
dermax = max(float(np.max(np.abs((Bm_ @ F[p, q]) - (Bd[p] + Bd[q]) * F[p, q]))) for p in range(NB) for q in range(NB))
res['B_is_derivation_max_violation'] = dermax
res['B_anticommutes_with_dagger'] = bool(np.allclose(T @ Bm_, -Bm_ @ T))
# B not inner: no y with ad_y = B
adall = np.stack([F[p].T.flatten() for p in range(NB)], axis=1)
yv, resid, *_ = np.linalg.lstsq(adall, Bm_.flatten().astype(complex), rcond=None)
res['B_inner_lstsq_residual'] = float(np.linalg.norm(adall @ yv - Bm_.flatten()))
fam_ok = True
for t in (0.3, -0.7):
    for sgn in (1, -1):
        E_ = np.diag(np.exp(t * Bd)) @ np.diag([1.0] * 30 + [sgn] * 32) @ T
        m1 = max(float(np.max(np.abs(E_ @ np.conj(F[p, q]) - np.einsum('p,q,pqr->r', E_[:, q], E_[:, p], F)))) for p in range(NB) for q in ODD)
        fam_ok &= (m1 < 1e-9) and np.allclose(E_ @ np.conj(E_), np.eye(NB))
res['family_exp_tB_pm_dagger_are_involutive_anti_automorphisms'] = bool(fam_ok)
# su(4): action of a real R in so(W_R^perp) on the W index of Q^i_alpha is anti-Hermitian
okah = True
for k in range(15, 30):
    Mw = np.array([[F[k, idx_xi(0, i), idx_xi(0, j)] for i in range(4)] for j in range(4)])   # [R, Q^i_1] = sum_j Mw[j,i] Q^j_1
    okah &= np.allclose(Mw.conj().T, -Mw)
    okah &= np.allclose(np.array([[F[k, idx_xi(1, i), idx_xi(1, j)] for i in range(4)] for j in range(4)]), Mw)
res['real_soW_acts_antihermitian_on_W_index'] = bool(okah)
R['4_radial_conjugation'] = res

# ---------------------------------------------------------------- 5. reflection conjugation (Def c8-Theta) and Prop c8-C
def daggern_matrix(n):
    n = np.array(n, float); n = n / np.linalg.norm(n)
    rn = np.eye(4) - 2 * np.outer(n, n)
    Tn = np.zeros((NB, NB), complex); rmax = 0.0
    for a in range(NB):
        E = basis[a]; out = np.zeros(NB, complex)
        if E['t'] == 'V':
            vf = tuple(np.conj(t) for t in E['vf'])
            vals = np.concatenate([-(rn @ vf_eval(vf, rn @ x)) for x in PTS])
            c, r = decomp_vf_values(vals); out[:15] = c; rmax = max(rmax, r)
        elif E['t'] == 'W':
            Rn = -np.conj(E['R'])
            for k, (I, J) in enumerate(w_pairs):
                out[15 + k] = Rn[I, J]
        else:
            RN = rho_V(n)
            Q1n = RN @ cc(E['Q1']); Q2n = -RN @ cc(E['Q2'])       # (Q1+rho(x)Q2)^{dag_n} = rho(n)conj Q1 - rho(x)rho(n) conj Q2
            out[30:46] = Q1n[:16]; out[46:62] = Q2n[16:]
            for x in PTS[:6]:
                lhs = Q1n + rho_V(x) @ Q2n
                rhs = RN @ cc(eps_of(E, rn @ x))                      # eps^{dag_n}(x) = rho(n) conj eps(r_n x)
                rmax = max(rmax, np.max(np.abs(lhs - rhs)))
        Tn[:, a] = out
    return Tn, rmax, rn
nvec = rng.normal(size=4); nvec /= np.linalg.norm(nvec)
Tn, rmaxn, rn = daggern_matrix(nvec)
res = {'fit_residual': float(rmaxn)}
res['involution'] = bool(np.allclose(Tn @ np.conj(Tn), np.eye(NB)))
amax = max(float(np.max(np.abs(Tn @ np.conj(F[p, q]) - br(Tn[:, q], Tn[:, p])))) for p in range(NB) for q in range(NB))
res['anti_automorphism_max_violation'] = amax
res['Delta_dagn_eq_minus_Delta'] = bool(np.allclose(Tn[:, iD], -np.eye(NB)[iD]))
cPn = np.zeros(NB, complex); cPn[:4] = nvec
res['P_n_dagn_eq_P_n'] = bool(np.allclose(Tn @ np.conj(cPn), cPn))
res['B_anticommutes_with_dagger_n'] = bool(np.allclose(Tn @ Bm_, -Bm_ @ Tn))
# two-component: (xi^{dag_n})_alpha^i = i n_m sigma^m conj(xibar_{alphadot i}), (xi^{dag_n})_{alphadot i} = i n_m conj(xi_alpha^i) sigma^m
cx = np.zeros(NB, complex); cx[30:46] = rng.normal(size=16) + 1j * rng.normal(size=16)
dn = Tn @ np.conj(cx)
xi = cx[30:38].reshape(2, 4); xib = cx[38:46].reshape(2, 4)
ns = sum(nvec[m] * sig[m] for m in range(4))
pred_xi = 1j * np.einsum('ab,bi->ai', ns, np.conj(xib)); pred_xib = 1j * np.einsum('ai,ab->bi', np.conj(xi), ns)
res['two_component_xi_dagn'] = bool(np.allclose(dn[30:38].reshape(2, 4), pred_xi) and np.allclose(dn[38:46].reshape(2, 4), pred_xib))
# [Q1, Q1^{dag_n}] = (Q1, conj Q1) P_n + (translation orthogonal to n)
c = np.zeros(NB, complex); c[30:46] = rng.normal(size=16) + 1j * rng.normal(size=16)
bb = br(c, Tn @ np.conj(c))
g = elem_to_geom(c)
res['[Q,Q^dagn]_P_n_component_over_(Q,conjQ)'] = str(np.round((bb[:4] @ nvec) / pair(g[2], cc(g[2])), 12))
res['[Q,Q^dagn]_is_translation'] = bool(np.allclose(bb[4:], 0))
# Prop c8-C: X^dagger = Ad(C)((Ad(C)^{-1} X)^{dag_n}), Ad(C) = exp(pi/4 ad_{P_n - K_n})
cY = np.zeros(NB, complex); cY[:4] = nvec; cY[11:15] = -nvec
adY = np.einsum('p,pqr->rq', cY, F)
A = expm(np.pi / 4 * adY); Ai = expm(-np.pi / 4 * adY)
# cross-check of expm by eigendecomposition (ad_{P_n-K_n} is diagonalisable with spectrum {0, +-i, +-2i})
wv, Vv = np.linalg.eig(adY)
A_eig = Vv @ np.diag(np.exp(np.pi / 4 * wv)) @ np.linalg.inv(Vv)
res['expm_vs_eig_max_dev'] = float(np.max(np.abs(A - A_eig)))
res['spectrum_of_ad_(P_n-K_n)'] = sorted({str(np.round(z, 9)) for z in wv})
Tcomp = A @ Tn @ np.conj(Ai)
res['Ad(C)_dagn_Ad(C)^-1_equals_dagger_max_dev'] = float(np.max(np.abs(Tcomp - T)))
A2 = expm(-np.pi / 4 * adY); A2i = expm(np.pi / 4 * adY)
Tcomp2 = A2 @ Tn @ np.conj(A2i)
sgnv = np.array([1.0] * 30 + [-1.0] * 32)
res['with_C_pi_over_4_gives_(-1)^|X|_dagger_max_dev'] = float(np.max(np.abs(Tcomp2 - T @ np.diag(sgnv))))
# geometric check of the map C and of C r_n C^{-1} = I
def flow(x, t, steps=4000):
    h = t / steps; y = x.copy()
    f = lambda z: nvec * (1 - z @ z) + 2 * (nvec @ z) * z
    for _ in range(steps):
        k1 = f(y); k2 = f(y + h / 2 * k1); k3 = f(y + h / 2 * k2); k4 = f(y + h * k3)
        y = y + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    return y
Cmap = lambda x: (2 * x + ((x @ x) - 1 - 2 * (nvec @ x)) * nvec) / ((x + nvec) @ (x + nvec))
Cinv = lambda x: (2 * x + (1 - (x @ x) - 2 * (nvec @ x)) * nvec) / ((x - nvec) @ (x - nvec))
xs = rng.normal(size=4) * 0.4
res['C_closed_form_equals_time_minus_pi_over_4_flow'] = float(np.max(np.abs(flow(xs, -np.pi / 4) - Cmap(xs))))
res['C_r_n_Cinv_equals_inversion'] = float(max(np.max(np.abs(Cmap(rn @ Cinv(y)) - y / (y @ y))) for y in [rng.normal(size=4) for _ in range(5)]))
okhs = True
for _ in range(20):
    y = rng.normal(size=4) * 2
    s_ = nvec @ y
    rC = np.linalg.norm(Cmap(y))
    okhs &= (rC < 1) if s_ > 0 else (rC > 1)
    yh = y - s_ * nvec
    okhs &= np.isclose(np.linalg.norm(Cmap(yh)), 1.0)
res['C_maps_half_space_to_ball_hyperplane_to_sphere'] = bool(okhs)
R['5_reflection_and_C'] = res

# ---------------------------------------------------------------- 6. block presentation (Thm ip-P (ii), C9 input)
blocks = {'P': (2, 1), 'K': (1, 2), 'S+W': (3, 1), 'S-W*': (2, 3), 'S+W*': (1, 3), 'S-W': (3, 2)}
def diag_derivation(s):
    d = np.zeros(NB)
    d[0:4] = s[1] - s[0]; d[11:15] = s[0] - s[1]
    d[30:38] = s[2] - s[0]; d[38:46] = s[1] - s[2]; d[46:54] = s[0] - s[2]; d[54:62] = s[2] - s[1]
    return d
res = {}
okder = True
for s in [(1, -1, 0), (0, 0, 1), (0.3, -1.2, 2.5)]:
    d = diag_derivation(s)
    m1 = max(float(np.max(np.abs(d * F[p, q] - (d[p] + d[q]) * F[p, q]))) for p in range(NB) for q in range(NB))
    okder &= m1 < 1e-9
res['block_rescalings_are_derivations'] = bool(okder)
res['s=(1,-1,0)_equals_ad_Delta'] = bool(np.allclose(np.diag(diag_derivation((1, -1, 0))), adD))
res['s=(0,0,1)_equals_B'] = bool(np.allclose(diag_derivation((0, 0, 1)), Bd))
R['6_blocks'] = res

import resource
R['peak_rss_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
with open(os.path.join(HERE, 'a1_results.json'), 'w') as fh:
    json.dump(R, fh, indent=1, default=str)
print(json.dumps(R, indent=1, default=str))
