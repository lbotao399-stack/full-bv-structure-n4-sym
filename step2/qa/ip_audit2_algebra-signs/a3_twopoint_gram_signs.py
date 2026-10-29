#!/usr/bin/env python3
"""Audit ip_audit2 (algebra and signs): two-point functions, Gram relations, Konishi sector numbers.

Independent of the project scripts.  Uses only printed definitions:
  ch01 Def c1-pairing/c1-bracket, ch03 Prop c3-P (propagator: A-A  2(-,-)^{-1} G, lambda-lambda (-,-)^{-1} rho(d_x) G
  read in Sigma (x) Sigma through the pairing, first coordinate functional <-> first tensor factor),
  ch08 Def c8-twocomp, Lemma c8-cc, Prop c8-C (Theta_I with factor -rho(x)/|x| on lambda), Def c8-2pt, Def c8-inner,
  Cor c8-gram; inner_product_design.tex (Prop ip-tree, ip-weight3, Rem ip-radius).
Sections:
  G1  Grassmann bookkeeping of Def c8-2pt: coefficient of t1 t2 in (t1 O1)(t2 O2) for odd O's, and the resulting
      form of the Ward identity for odd generators.
  G2  free radial two-point values: O_zeta (fermionic) and K (bosonic) via the expectation hbar d_P(a^Theta b),
      with s_lambda = -1, at random x in the unit ball; sign and modulus.
  G3  coefficient-column Gram relation: <X a, b> = <a, Y b> <=> conj(X)^T Gram = Gram Y; order expansion; component lines.
  G4  Konishi weight-3 data: two-component rho(1/2[A^A]) xi, B_w, vanishing of the other terms of {S^(1)(xi),O_zeta},
      flavour sum of B_w.
  G5  radial Gram for an inversion sphere of radius R (Remark ip-radius).
Output: a3_results.json
"""
import json, os, itertools, resource
import numpy as np
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(31415)
OUT = {}

# ---------------------------------------------------------------- realisation (Def c8-twocomp)
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
sig = [s1, s2, s3, 1j * np.eye(2)]; sigb = [-s1, -s2, -s3, 1j * np.eye(2)]
epsU = np.array([[0, 1], [-1, 0]], complex); epsL = -epsU
eps3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    eps3[p] = np.linalg.det(np.eye(3)[list(p)])
eps4 = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(4)):
    eps4[p] = np.linalg.det(np.eye(4)[list(p)])
Sig = np.zeros((6, 4, 4), complex)
for a in range(3):
    for b in range(3):
        for c in range(3):
            Sig[a, b, c] = eps3[a, b, c]; Sig[a + 3, b, c] = 1j * eps3[a, b, c]
        Sig[a, b, 3] = 1.0 * (a == b); Sig[a, 3, b] = -1.0 * (a == b)
        Sig[a + 3, b, 3] = -1j * (a == b); Sig[a + 3, 3, b] = 1j * (a == b)
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
def rho_W(w):
    M = np.zeros((32, 32), complex)
    wl = np.einsum('I,Iij->ij', w, Sig); wu = 0.5 * np.einsum('ijkl,kl->ij', eps4, wl)
    for a in range(2):
        for b in range(2):
            for j in range(4):
                for i in range(4):
                    M[ix('eta', a, j), ix('xi', b, i)] += epsU[a, b] * wl[j, i]
                    M[ix('etab', a, j), ix('xib', b, i)] += epsU[a, b] * wu[j, i]
                    M[ix('xi', a, j), ix('eta', b, i)] += -wu[j, i] * epsL[a, b]
                    M[ix('xib', a, j), ix('etab', b, i)] += -wl[j, i] * epsL[a, b]
    return M
RHO = [rho_V(np.eye(4)[m]) for m in range(4)] + [rho_W(np.eye(6)[I]) for I in range(6)]
Gp = np.zeros((16, 16)); Gp[:8, :8] = np.eye(8); Gp[8:, 8:] = -np.eye(8)
Bf = np.zeros((32, 32)); Bf[:16, 16:] = Gp; Bf[16:, :16] = Gp.T
def pair(z1, z2):
    return z1 @ Bf @ z2
def cc(z):
    out = np.zeros_like(z)
    out[16:24] = np.conj(z[0:8]); out[24:32] = -np.conj(z[8:16])
    out[0:8] = np.conj(z[16:24]); out[8:16] = -np.conj(z[24:32])
    return out
def rhox(x):
    return sum(x[m] * RHO[m] for m in range(4))

# ---------------------------------------------------------------- G1 Grassmann bookkeeping
class Grass:
    """exterior algebra on named generators with parities; elements: dict{tuple(names): coeff}"""
    def __init__(self, terms, par):
        self.t = terms; self.par = par
    @staticmethod
    def gen(name, par):
        return Grass({(name,): 1}, par)
    def __mul__(self, o):
        out = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                word = list(k1) + list(k2)
                odd = [w for w in word if self.par[w] == 1]
                if len(set(odd)) < len(odd):
                    continue
                # sort word: order odd generators by name with Koszul sign (even ones commute)
                sign = 1; w = word[:]
                for i in range(len(w)):
                    for j in range(len(w) - 1 - i):
                        if w[j] > w[j + 1]:
                            if self.par[w[j]] == 1 and self.par[w[j + 1]] == 1:
                                sign = -sign
                            w[j], w[j + 1] = w[j + 1], w[j]
                key = tuple(w)
                out[key] = out.get(key, 0) + sign * v1 * v2
        return Grass({k: v for k, v in out.items() if v != 0}, self.par)
par = {'a_t1': 1, 'b_t2': 1, 'c_O1': 1, 'd_O2': 1}
t1, t2, O1, O2 = (Grass.gen(n, par) for n in ['a_t1', 'b_t2', 'c_O1', 'd_O2'])
prod = (t1 * O1) * (t2 * O2)          # (t1 O1)(t2 O2), names sorted t1 < t2 < O1 < O2
OUT['G1_grassmann'] = {
    'coefficient_of_t1t2_times_O1O2_in_(t1O1)(t2O2)_odd_O': prod.t.get(('a_t1', 'b_t2', 'c_O1', 'd_O2'), 0),
    'meaning': 'Def c8-2pt (coefficient of t1 t2) equals (-1)^{|O1||O2|} times the expectation hbar d_P(O1 O2) at tree level',
}
# Ward identity in both normalisations (pure parity algebra): E(L_F G) = -(-1)^{|Y||F|} E(F L_G) (derivation rule);
# <A,B>_t = (-1)^{|A||B|} E(AB).  For |G| = |F| + |Y|:
ward = {}
for Y in (0, 1):
    for Fp in (0, 1):
        Gp_ = (Fp + Y) % 2
        sE = -(-1) ** (Y * Fp)                       # E(L F . G) = sE E(F . L G)
        conv_left = (-1) ** (((Fp + Y) % 2) * Gp_)  # <L F, G>_t = conv_left E(L F G)
        conv_right = (-1) ** (Fp * ((Gp_ + Y) % 2))  # <F, L G>_t = conv_right E(F L G)
        s_t = conv_left * sE * conv_right             # <L F, G>_t = s_t <F, L G>_t
        ward['|Y|=%d,|F|=%d' % (Y, Fp)] = {'expectation_form_sign': sE, 't1t2_form_sign': s_t,
                                           'hypothesis_of_thm_c8_adjoint': -(-1) ** (Y * Fp)}
OUT['G1_ward_signs'] = ward

# ---------------------------------------------------------------- G2 free radial two-point values
def slot_dual(sp_):        # dual basis vector f_s in Sigma^* with (e_s, f_s) = 1 for a Sigma index s
    f = np.zeros(32); f[16 + sp_] = Gp[sp_, sp_]
    return f
FD = np.array([slot_dual(s) for s in range(16)])                 # (16, 32)
def Pll(u):        # lambda-lambda propagator component (s, s') at separation u = x - y, without (-,-)^{-1}
    dG = -2 * u / (4 * np.pi ** 2 * (u @ u) ** 2)
    M = rhox(dG)
    return np.array([[(M @ FD[sp_])[s] for sp_ in range(16)] for s in range(16)])
def GA(u):
    return 1.0 / (4 * np.pi ** 2 * (u @ u))
def O_zeta_coeff(zeta):
    """O_zeta(x) = sum_{a,I} (lambda^a(x), rho(w_I) zeta) A_I^a(x): coefficient T[s, I] (colour diagonal)."""
    return np.array([[(Bf @ (RHO[4 + I] @ zeta))[s] for I in range(6)] for s in range(16)])
def O_zeta_theta_coeff(zeta, x, slam=-1.0):
    """[O_zeta(Theta_I e)(x)]^c as sum T[s,I] lambda_s^a(Ix) A_I^a(Ix): computed by evaluation on unit fields.
    Theta_I: lambda -> slam (x,x)^{-2} rho(x) cc(lambda(Ix)), A_W -> (x,x)^{-1} conj(A_W(Ix)); one odd factor, so ^c is
    complex conjugation of the coefficient."""
    xx = x @ x; T = np.zeros((16, 6), complex)
    for s in range(16):
        lam = np.zeros(32, complex); lam[s] = 1.0
        tl = slam * (rhox(x) @ cc(lam)) / xx ** 2
        for I in range(6):
            Aw = np.zeros(6); Aw[I] = 1.0
            tA = Aw / xx
            val = sum(pair(tl, RHO[4 + J] @ zeta) * tA[J] for J in range(6))
            T[s, I] = np.conj(val)
    return T
def two_point_O(Ta, Tb, u):
    """expectation hbar d_P(O_a(u) O_b(0)) for O = sum T[s,I] lambda_s^a A_I^a, in units of hbar dim g."""
    P = Pll(u)                      # <lambda_s(u) lambda_s'(0)>
    return np.einsum('sI,tJ,st,IJ->', Ta, Tb, P, 2 * GA(u) * np.eye(6))
g2 = {'O_zeta_ratio_to_24(zetabar,zeta)/(16pi^4)': [], 'x_independence_spread': None}
vals = []
for _ in range(6):
    x = rng.normal(size=4); x = 0.8 * rng.random() * x / np.linalg.norm(x)
    zeta = np.zeros(32, complex); zeta[:8] = rng.normal(size=8) + 1j * rng.normal(size=8)
    zp = np.zeros(32, complex); zp[:8] = rng.normal(size=8) + 1j * rng.normal(size=8)
    u = x / (x @ x)
    E = two_point_O(O_zeta_theta_coeff(zeta, x), O_zeta_coeff(zp), u)
    ref = 24 / (16 * np.pi ** 4) * pair(zp, cc(zeta))
    vals.append(E / ref)
g2['O_zeta_ratio_to_24(zetabar,zeta)/(16pi^4)'] = [str(np.round(v, 12)) for v in vals]
g2['all_equal_plus_one'] = bool(np.allclose(vals, 1.0))
# control: s_lambda = +1 gives the opposite sign
x = np.array([0.3, -0.2, 0.1, 0.4]); zeta = np.zeros(32, complex); zeta[:8] = rng.normal(size=8)
E_c = two_point_O(O_zeta_theta_coeff(zeta, x, slam=+1.0), O_zeta_coeff(zeta), x / (x @ x))
g2['control_slambda_plus_ratio'] = str(np.round(E_c / (24 / (16 * np.pi ** 4) * pair(zeta, cc(zeta))), 12))
# bosonic K = sum_{a,I} (A_I^a)^2: K^Theta(x) = (x,x)^{-2} K(Ix);  E = 2 * (2 G)^2 * 6 dim g
x = np.array([0.1, 0.5, -0.3, 0.2]); u = x / (x @ x)
EK = (x @ x) ** -2 * 2 * (2 * GA(u)) ** 2 * 6
g2['K_ratio_to_48/(16pi^4)'] = float(EK / (48 / (16 * np.pi ** 4)))
g2['Def_c8-2pt_value_on_O_zeta'] = '-1 times the expectation value above (G1): negative on the fermionic family'
OUT['G2_free_radial_values'] = g2

# ---------------------------------------------------------------- G3 coefficient-column Gram algebra
def rmat(n, m, cplx=True):
    return sp.Matrix(n, m, lambda i, j: sp.Rational(int(rng.integers(-3, 4)), int(rng.integers(1, 4))) +
                     (sp.I * sp.Rational(int(rng.integers(-3, 4)), int(rng.integers(1, 4))) if cplx else 0))
n = 3
A0 = rmat(n, n); Gt = A0.H * A0 + sp.eye(n)
G1 = rmat(n, n); G1 = (G1 + G1.H) / 2
X0 = rmat(n, n); X1 = rmat(n, n)
h = sp.symbols('h', real=True)
Gram = Gt + h * G1; X = X0 + h * X1
Y = Gram.inv() * X.H * Gram
a_ = rmat(n, 1); b_ = rmat(n, 1)
ip = lambda u_, v_: (u_.H * Gram * v_)[0, 0]          # <a,b> = conj(a)^T Gram b
ok_adj = sp.simplify(ip(X * a_, b_) - ip(a_, Y * b_)) == 0
Y0 = Y.applyfunc(lambda e: sp.limit(e, h, 0)); Y1 = Y.applyfunc(lambda e: sp.diff(e, h).subs(h, 0))
adjm = lambda Mx: Gt.inv() * Mx.H * Gt
Mm = Gt.inv() * G1
ok_ord = sp.simplify(Y0 - adjm(X0)) == sp.zeros(n) and sp.simplify(Y1 - (adjm(X1) + Y0 * Mm - Mm * Y0)) == sp.zeros(n)
# component line for X = Q^i_alpha, X^dag = -S^alpha_i: X = -(C_Q + h q1), Y = matrix of -Lambda(-S) = +(C_S + h s1)
CQ = rmat(n, n); q1 = rmat(n, n)
Gram2 = Gt + h * G1
Xq = -(CQ + h * q1); Yq = Gram2.inv() * Xq.H * Gram2
CS = Yq.applyfunc(lambda e: sp.limit(e, h, 0))
s1_ = Yq.applyfunc(lambda e: sp.diff(e, h).subs(h, 0))
ok_line = sp.simplify(s1_ - (-adjm(q1) + CS * Mm - Mm * CS)) == sp.zeros(n)
OUT['G3_gram_algebra'] = {'adjointness_iff_conjXT_Gram_eq_Gram_Y': bool(ok_adj), 'order_expansion_cor_c8_gram': bool(ok_ord),
                          'S1_line_minus_Q1_adj_plus_commutator': bool(ok_line)}

# ---------------------------------------------------------------- G4 Konishi weight-3 data
g4 = {}
# rho(q(1/2[A^A])) xi, A_W = sum_I A_I^a e_a with su(2) colour; two-component claim -(1/2)[A^{ik},A_{kj}] u_alpha w^j
Acol = rng.normal(size=(6, 3))                               # real field values A_I^a
def br(X_, Y_):
    return np.einsum('abc,b,c->a', eps3, X_, Y_)
u = rng.normal(size=2) + 1j * rng.normal(size=2); w = rng.normal(size=4) + 1j * rng.normal(size=4)
xi = np.zeros(32, complex)
for a in range(2):
    for i in range(4):
        xi[ix('xi', a, i)] = u[a] * w[i]
rF = np.zeros((32, 3), complex)
for I in range(6):
    for J in range(I + 1, 6):
        rF += np.einsum('s,a->sa', RHO[4 + I] @ RHO[4 + J] @ xi, br(Acol[I], Acol[J]))
Al = np.einsum('Ia,Iij->ija', Acol, Sig); Au = 0.5 * np.einsum('ijkl,kla->ija', eps4, Al)
pred = np.zeros((32, 3), complex)
for a in range(2):
    for i in range(4):
        pred[ix('xi', a, i)] = -0.5 * sum(br(Au[i, k], Al[k, j]) * u[a] * w[j] for k in range(4) for j in range(4))
g4['rho_half_bracket_xi_two_component'] = bool(np.allclose(rF, pred))
# B_w = -2 (rho(1/2[A^A]) xi, rho(A_W) zeta) = ([A^{ik}, A_{kj} w^j], A_{il} w^l) for zeta = u' (x) w, u'^alpha u_alpha = 1
up_ = np.array([-np.conj(u[1]), np.conj(u[0])]) / (abs(u[0]) ** 2 + abs(u[1]) ** 2)
upU = epsU @ up_                                             # u'^alpha = eps^{alpha beta} u'_beta
g4['uprime_contracts_to_1'] = str(np.round(upU @ u, 12))
g4['zeta_orthogonal_to_conj_xi'] = bool(np.isclose(np.vdot(np.kron(u, w), np.kron(up_, w)), 0))
zeta = np.zeros(32, complex)
for a in range(2):
    for i in range(4):
        zeta[ix('xi', a, i)] = up_[a] * w[i]
rAz = np.einsum('Ist,t,Ia->sa', np.array(RHO[4:]), zeta, Acol)
Bw_index_free = -2 * np.einsum('sa,st,ta->', rF, Bf, rAz)
Bw_two_comp = sum(np.dot(br(Au[i, k], Al[k, j] * w[j]), Al[i, l] * w[l]) for i in range(4) for k in range(4) for j in range(4) for l in range(4))
g4['B_w_index_free_equals_two_component'] = bool(np.isclose(Bw_index_free, Bw_two_comp))
# other terms of {S^(1)(xi), O_zeta} vanish for xi = u (x) w, zeta = u' (x) w:
# (rho(F_V) xi, rho(A_W) zeta), (rho(F_{V W}) xi, rho(A_W) zeta) and (lambda, rho(Gamma(xi,lambda)|_W) zeta)
okv = True
for m in range(4):
    for n_ in range(m + 1, 4):
        okv &= np.isclose(pair(RHO[m] @ RHO[n_] @ xi, RHO[4] @ zeta), 0) and np.allclose([pair(RHO[m] @ RHO[n_] @ xi, RHO[4 + I] @ zeta) for I in range(6)], 0)
    for I in range(6):
        okv &= np.allclose([pair(RHO[m] @ RHO[4 + I] @ xi, RHO[4 + J] @ zeta) for J in range(6)], 0)
# Gamma(xi, lambda)|_W for arbitrary lambda in Sigma: (rho(Gamma|_W) zeta) paired with any lambda' vanishes
okf = True
for _ in range(3):
    lam = np.zeros(32, complex); lam[:16] = rng.normal(size=16) + 1j * rng.normal(size=16)
    GW = np.array([pair(RHO[4 + I] @ xi, lam) for I in range(6)])
    v = sum(GW[I] * RHO[4 + I] for I in range(6)) @ zeta
    okf &= np.allclose(v[16:24], 0) and np.allclose(v, 0)
g4['other_terms_vanish_V_and_mixed_field_strength'] = bool(okv)
g4['fermion_term_rho(Gamma(xi,lambda)|_W)zeta_vanishes'] = bool(okf)
# flavour sum: c_IJK = conj(Sigma^I_ik) Sigma^J_kj w^j Sigma^K_il w^l (A^{ik} = A_I conj Sigma^I_ik)
assert np.allclose(np.conj(Sig), 0.5 * np.einsum('ijkl,Ikl->Iij', eps4, Sig))
fl = []
for _ in range(4):
    w = rng.normal(size=4) + 1j * rng.normal(size=4)
    cijk = np.einsum('Iik,Jkj,j,Kil,l->IJK', np.conj(Sig), Sig, w, Sig, w)
    ca = sum(np.linalg.det(np.eye(3)[list(p)]) * np.transpose(cijk, p) for p in itertools.permutations(range(3))) / 6
    fl.append(np.sum(abs(ca) ** 2) / np.vdot(w, w).real ** 2)
g4['flavour_sum_over_|w|^4'] = [float(np.round(v, 10)) for v in fl]
# k from (iv) with the positive tree norms (expectation convention): k = |kappa_B|^2 n_B /(unit n_zeta (xi,conj xi))
hb, dimg, lam_ = sp.symbols('hbar dim_g lambda')
w2, u2 = sp.symbols('w2 u2', positive=True)
nB = hb ** 2 * (2 / (4 * sp.pi ** 2)) ** 3 * (-lam_ * dimg) * 6 * 48 * w2 ** 2
nz = 24 * hb * dimg * (1 / u2) * w2 / (16 * sp.pi ** 4)
unit = hb * lam_ / (16 * sp.pi ** 2)
k_E = sp.simplify(sp.Rational(1, 16) * nB / (unit * nz * u2 * w2))
g4['k_expectation_convention'] = str(k_E)
g4['k_if_n_zeta_taken_from_literal_Def_c8-2pt_(sign_flipped)_with_the_same_relations'] = str(sp.simplify(-k_E))
OUT['G4_konishi_weight3'] = g4

# ---------------------------------------------------------------- G5 sphere radius
Dl, R_, c_ = sp.symbols('Delta R c', positive=True)
xs = sp.symbols('x', positive=True)
# scalar of dimension Delta: <conj f(z) f(0)> = c |z|^{-2 Delta}; inversion in radius R: x -> R^2 x/|x|^2, conformal factor R^2/|x|^2
gramR = sp.limit((R_ ** 2 / xs ** 2) ** Dl * c_ * (R_ ** 2 / xs) ** (-2 * Dl), xs, 0)
OUT['G5_radius'] = {'radial_Gram_radius_R': str(sp.simplify(gramR)), 'ratio_to_unit_sphere': str(sp.simplify(gramR / c_))}

OUT['peak_rss_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
with open(os.path.join(HERE, 'a3_results.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=str)
print(json.dumps(OUT, indent=1, default=str))
