#!/usr/bin/env python3
"""Checks for step2/theory/inner_product_design.tex (labels ip-...).

Single process, exact arithmetic (sympy / fractions) except C7 (numpy, double precision, two-component
spinor identities of ch08 Def c8-twocomp).  Writes inner_product_checks.json next to this file.
Code-level names (json/py only): G0 = Gram_tree, G1 = Gram_1, Gh = order-sqrt(unit) part of Gram,
X0, Xh, X1 = matrix of -Lambda^(1)(X) by order, Y* = matrix of -Lambda^(1)(X^dagger), t = sqrt(unit).

C1  graded Gram relation (Thm ip-graded): Y = Gram^{-1} conj(X)^T Gram expanded in t, orders 0, 1/2, 1.
C2  covariance of the relations under a basis change 1 + t^2 psi' (Cor c8-gram, Thm ip-laws).
C3  radial Gram at the origin for a scalar family with a mixing D_1 (Thm ip-origin): log|x|^2 coefficient
    (1/4)(G0 D1 - conj(D1)^T G0); vanishes iff D1 = D1^adj.
C4  descendants (Prop ip-desc): lim d_x d_y [ (x,x)^{-D} |I(x)-y|^{-2D} ] = 2 D delta (level one) and the
    level-two Laplacian pairing 32 D (D-1) in four dimensions.
C5  sl(2) model of the P_1 frames (Thm ip-P): P_1 = r P_0 frame vs P_1 = 0 frame, K line, Gram factor.
C6  Konishi pair (Prop ip-konishi): MC and S line give s+q = kappa/(2 unit), s - conj(q) = m_zeta - m_K,
    tree ratio n_zeta/n_K = -sigma/tau = 2/tau^2; invariance under G_B (imaginary a) and basis changes.
C7  two-component identities (numpy): conj(sigma^m) = -sigmabar^m (index transposed), Gamma(xi,xi')|_V = 0
    and (xi, rho(b) xi') = 0 on S_+ (x) W, (xi, conj xi) > 0, sum_I rho(w_I) rho(x) rho(w_I) = -6 rho(x).
C8  weight 3 (Prop ip-weight3): orthogonality of the Konishi and protected descendants fixes c_K by tree
    data; covariance under psi(B_w) = x Fbar_w; the order-sqrt(unit) relation U = -conj(kappa_B) n_B/n_zeta.
C9  one-loop rescalings of generators compatible with the bracket (Thm ip-P (ii)): on gl(4|4) block
    supermatrices the solutions are r = s_a - s_b (inner by Delta and outer B); r_Delta = r_rot = 0,
    r_K = -r_P, r_Q + r_S = 0.
C10 lemmas of the two-point design used here: Fourier transforms of log|q|^2, log^2|q|^2 from the
    Gamma-function formula; localisation constant log(2/3)/(16 pi^2).
C11 Konishi anomalous dimension from the classical length-changing supercharge: |kappa_B|^2 Gram(B_w,B_w) =
    k unit (xi, conj xi) Gram(O_zeta,O_zeta) with the tree norms; flavour sum 48|w|^4 of B_w with
    c_IJK = conj(Sigma^I_ik) Sigma^J_kj w^j Sigma^K_il w^l and A^{ij} = A_I conj(Sigma^I_ij); predicts k = -24.
C12 two-point convention (revision 2026-10-29): coefficient of t1 t2 vs t2 t1 in (t1 O1)(t2 O2); Ward sign table;
    exact realisation of ch08 Def c8-twocomp: radial value of O_zeta with Theta_I (+24 for the t2t1 form) and the
    separated-point value 24 (rho(x) xi, chi)/(16 pi^4 |x|^6).
C13 radius of the quantisation sphere: the radial Gram scales as R^{-2 Delta}; log R = -c unit for 1 + c unit w.
C14 order-one-half census between B_z (three scalars) and two-fermion weight-3 bilinears: no surviving graph.
C15 frames (route (a), inner direction): two-point coefficients unchanged, radial Gram of the engines' frame
    = (1 - 2 c unit w) x two-point; mixed-frame residuals; the previous IP6 rule fails; Konishi pair residual.
C16 reflection compatibility of universal rescalings: real ones force r_Q = r_P/2 (inner Delta allowed, real B
    excluded); imaginary ones force r_P = 0 (imaginary B).
C17 weight 3: Gram_tree D_1 Hermitian iff the Konishi and protected descendants are orthogonal (T5 = T3 (iii)).
C8 (corrected 2026-10-29): U = -conj(kappa_B) n_B/n_zeta, derived from the matrices and asserted.
"""
import json, os, random, itertools
import sympy as sp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = {}
random.seed(20261028)

def rmat(n, m, cplx=True, rng=3):
    return sp.Matrix(n, m, lambda i, j: sp.Rational(random.randint(-rng, rng), random.randint(1, 3)) +
                     (sp.I * sp.Rational(random.randint(-rng, rng), random.randint(1, 3)) if cplx else 0))

def herm_pd(n):
    A = rmat(n, n)
    return A.H * A + sp.eye(n)

def adj(Y, G0):
    return G0.inv() * Y.H * G0

def is_zero(M):
    return all(sp.simplify(e) == 0 for e in M)

# ------------------------------------------------------------------------------------------------
# C1 graded Gram relation
# ------------------------------------------------------------------------------------------------
def C1():
    t = sp.symbols('t', real=True)
    nL, nL1 = 2, 2          # dimensions of the length-L and length-(L+1) blocks
    n = nL + nL1
    G0 = sp.zeros(n)
    G0[:nL, :nL] = herm_pd(nL); G0[nL:, nL:] = herm_pd(nL1)
    Gh = sp.zeros(n)                       # order-sqrt(unit) entries only between lengths
    B = rmat(nL, nL1); Gh[:nL, nL:] = B; Gh[nL:, :nL] = B.H
    G1 = herm_pd(n) - sp.eye(n)            # arbitrary Hermitian order-unit part
    G1 = (G1 + G1.H) / 2
    G = G0 + t * Gh + t**2 * G1
    X0 = sp.zeros(n); X0[:nL, :nL] = rmat(nL, nL); X0[nL:, nL:] = rmat(nL1, nL1)   # length preserving
    Xh = sp.zeros(n); Xh[nL:, :nL] = rmat(nL1, nL)                                # length L -> L+1
    X1 = rmat(n, n)
    X = X0 + t * Xh + t**2 * X1
    Y = G.inv() * X.H * G
    Yser = [sp.simplify(Y.applyfunc(lambda e: sp.diff(e, t, k).subs(t, 0) / sp.factorial(k))) for k in range(3)]
    M_h = G0.inv() * Gh
    M_1 = G0.inv() * G1
    Y0p = adj(X0, G0)
    Yhp = adj(Xh, G0) + Y0p * M_h - M_h * Y0p
    Y1p = adj(X1, G0) + (Y0p * M_1 - M_1 * Y0p) + (adj(Xh, G0) * M_h - M_h * adj(Xh, G0)) - M_h * (Y0p * M_h - M_h * Y0p)
    res = {'order0': is_zero(Yser[0] - Y0p), 'order_half': is_zero(Yser[1] - Yhp), 'order1': is_zero(Yser[2] - Y1p)}
    # with Gh = 0 the order-1 relation is that of Cor c8-gram
    Gb = G0 + t**2 * G1
    Xb = X0 + t * Xh + t**2 * X1
    Yb = Gb.inv() * Xb.H * Gb
    Yb1 = sp.simplify(Yb.applyfunc(lambda e: sp.diff(e, t, 2).subs(t, 0) / 2))
    Ybh = sp.simplify(Yb.applyfunc(lambda e: sp.diff(e, t, 1).subs(t, 0)))
    res['Gh_zero_order1_is_cor_c8_gram'] = is_zero(Yb1 - (adj(X1, G0) + Y0p * M_1 - M_1 * Y0p))
    res['Gh_zero_order_half_is_tree_adjoint_of_Xh'] = is_zero(Ybh - adj(Xh, G0))
    # Yh is length decreasing when Gh = 0
    res['Yh_length_decreasing'] = is_zero(Ybh[nL:, :]) and is_zero(Ybh[:nL, :nL])
    OUT['C1_graded_relation'] = res
    return all(res.values())

# ------------------------------------------------------------------------------------------------
# C2 covariance under a basis change a -> (1 + hbar psi') a
# ------------------------------------------------------------------------------------------------
def C2():
    n = 3
    G0 = herm_pd(n); G1 = rmat(n, n); G1 = (G1 + G1.H) / 2
    X0 = rmat(n, n); X1 = rmat(n, n); psi = rmat(n, n)
    Y0 = adj(X0, G0); M = G0.inv() * G1
    Y1 = adj(X1, G0) + Y0 * M - M * Y0
    # transformed data (first order): Gram -> conj(1+psi)^T Gram (1+psi), X -> (1+psi)^{-1} X (1+psi)
    X1n = X1 + X0 * psi - psi * X0
    Y1n = Y1 + Y0 * psi - psi * Y0
    G1n = G1 + psi.H * G0 + G0 * psi
    Mn = G0.inv() * G1n
    lhs = Y1n; rhs = adj(X1n, G0) + Y0 * Mn - Mn * Y0
    res = {'relation_preserved': is_zero(lhs - rhs),
           'M_shift_is_psi_plus_psiadj': is_zero(Mn - M - psi - adj(psi, G0))}
    OUT['C2_covariance'] = res
    return all(res.values())

# ------------------------------------------------------------------------------------------------
# C3 radial Gram at the origin with a mixing matrix D_1
# ------------------------------------------------------------------------------------------------
def C3():
    u, L = sp.symbols('u L', real=True)      # L = log|x|^2 (units l = 1)
    d = sp.symbols('d0:4')
    D1 = sp.Matrix(2, 2, [d[0] + sp.I*d[1], d[2], d[3], -d[0]])   # generic complex mixing matrix
    g = sp.symbols('g0:3', real=True)
    H0 = sp.Matrix([[g[0] + 2, g[1] + sp.I*g[2]], [g[1] - sp.I*g[2], 3]])   # Hermitian tree Gram (symbolic)
    f = sp.symbols('f0:4')
    H1 = sp.Matrix(2, 2, f)
    # H(y): bilinear two-point of (conj family at y, family at 0): |y|^{-2w}[H0 + u H1 - (1/4)u(conj(D1)^T H0 + H0 D1) log|y|^2]
    # at y = I(x): log|y|^2 = -L ; conjugation factor (x,x)^{-w - u conj(D1)/2} -> (1 - (u/2) conj(D1) L)
    Hy = H0 + u * H1 + sp.Rational(1, 4) * u * (D1.conjugate().T * H0 + H0 * D1) * L
    Conj = sp.eye(2) - sp.Rational(1, 2) * u * D1.conjugate() * L
    Gram = (Conj.T * Hy).applyfunc(lambda e: sp.expand(e))
    Gram1 = Gram.applyfunc(lambda e: sp.expand(e).coeff(u, 1))
    logcoef = Gram1.applyfunc(lambda e: sp.expand(e).coeff(L, 1))
    pred = sp.Rational(1, 4) * (H0 * D1 - D1.conjugate().T * H0)
    ok1 = is_zero(logcoef - pred)
    # self-adjoint D1 (w.r.t. H0) kills it: build D1 = H0^{-1} A with A Hermitian
    A = sp.Matrix([[1, 2 - sp.I], [2 + sp.I, -3]])
    H0n = H0.subs({g[0]: 1, g[1]: sp.Rational(1, 3), g[2]: sp.Rational(1, 5)})
    D1sa = H0n.inv() * A
    ok2 = is_zero(sp.simplify(H0n * D1sa - D1sa.conjugate().T * H0n))
    OUT['C3_origin_limit'] = {'log_coefficient_formula': ok1, 'self_adjoint_D1_gives_x_independence': ok2}
    return ok1 and ok2

# ------------------------------------------------------------------------------------------------
# C4 descendants: transported conformal two-point kernel
# ------------------------------------------------------------------------------------------------
def C4():
    D = sp.symbols('Delta', positive=True)
    x = sp.symbols('x1:5', real=True); y = sp.symbols('y1:5', real=True)
    xx = sum(xi**2 for xi in x); yy = sum(yi**2 for yi in y); xy = sum(a*b for a, b in zip(x, y))
    Ix = [xi / xx for xi in x]
    dist2 = sp.simplify(sum((a - b)**2 for a, b in zip(Ix, y)))
    ok_dist = sp.simplify(dist2 - (1 - 2*xy + xx*yy) / xx) == 0
    F = (1 - 2*xy + xx*yy)**(-D)          # (x,x)^{-D} |I(x)-y|^{-2D}
    zero = {**{xi: 0 for xi in x}, **{yi: 0 for yi in y}}
    lvl1 = sp.Matrix(4, 4, lambda i, j: sp.simplify(sp.diff(F, x[i], y[j]).subs(zero)))
    ok1 = is_zero(lvl1 - 2*D*sp.eye(4))
    lapx = lambda e: sum(sp.diff(e, xi, 2) for xi in x)
    lapy = lambda e: sum(sp.diff(e, yi, 2) for yi in y)
    # degree-(2,2) part of F suffices: F = 1 - D s + D(D+1)/2 s^2 + ..., s = -2xy + xx yy
    s = -2*xy + xx*yy
    F22 = sp.expand(-D*(xx*yy) + D*(D + 1)/2 * 4*xy**2)
    lvl2 = sp.simplify(lapx(lapy(F22)))
    ok2 = sp.simplify(lvl2 - 32*D*(D - 1)) == 0
    OUT['C4_descendants'] = {'distance_identity': ok_dist, 'level_one_equals_2Delta_delta': ok1,
                             'level_two_laplacian_pairing': str(sp.factor(lvl2)), 'level_two_ok': ok2}
    return ok_dist and ok1 and ok2

# ------------------------------------------------------------------------------------------------
# C5 sl(2) model of P_1 frames
# ------------------------------------------------------------------------------------------------
def C5():
    w, kap, r, u = sp.symbols('w kappa r u', real=True)
    N = 5
    def gram(h):
        return sp.diag(*[sp.factorial(n) * sp.rf(2*h, n) for n in range(N)])
    def Pm():
        M = sp.zeros(N)
        for n in range(N - 1): M[n + 1, n] = 1
        return M
    def Km(h):
        M = sp.zeros(N)
        for n in range(1, N): M[n - 1, n] = n * (2*h + n - 1)
        return M
    h = w + u * kap
    G = gram(h); P = Pm(); K = Km(h)
    ser = lambda M, k: M.applyfunc(lambda e: sp.expand(sp.series(e, u, 0, 2).removeO()).coeff(u, k))
    G0, G1 = ser(G, 0), ser(G, 1)
    # exact adjointness in the model: K = G^{-1} P^T G
    ok_exact = is_zero(sp.simplify(G.inv() * P.T * G - K))
    K0, K1 = ser(K, 0), ser(K, 1)
    Mm = G0.inv() * G1
    ok_Pzero = is_zero(sp.simplify(K1 - (K0 * Mm - Mm * K0)))          # K line with P_1 = 0
    # frame change U = 1 + u r Nop (Nop = level): P -> U^{-1} P U, K -> U^{-1} K U, G -> U^T G U
    Nop = sp.diag(*range(N))
    U = sp.eye(N) + u * r * Nop; Ui = sp.eye(N) - u * r * Nop
    P2 = (Ui * P * U).applyfunc(sp.expand); K2 = (Ui * K * U).applyfunc(sp.expand); G2 = (U.T * G * U).applyfunc(sp.expand)
    P20, P21 = ser(P2, 0), ser(P2, 1); K20, K21 = ser(K2, 0), ser(K2, 1); G20, G21 = ser(G2, 0), ser(G2, 1)
    M2 = G20.inv() * G21
    ok_P1 = is_zero(P21 + r * P20)                       # P_1 = -r P_0 (a rescaling)
    ok_K1 = is_zero(sp.simplify(K21 - (K1 + r * K0)))     # K rescaled oppositely
    ok_Kline = is_zero(sp.simplify(K21 - (adj(P21, G20) + K20 * M2 - M2 * K20)))
    ok_G = is_zero(sp.simplify(G21 - (G1 + 2 * r * Nop * G0)))
    # level one: K_1 on P O equals 2 kappa (the conformal form), level-one Gram ratio kappa/w
    ok_lvl1 = sp.simplify(Mm[1, 1] - Mm[0, 0] - kap / w) == 0
    OUT['C5_sl2_frames'] = {'exact_adjointness': ok_exact, 'K_line_P1_zero': ok_Pzero, 'P1_rescaling': ok_P1,
                            'K1_opposite_rescaling': ok_K1, 'K_line_general': ok_Kline,
                            'Gram_frame_factor_1_plus_2r_level': ok_G, 'level_one_minus_level_zero_kappa_over_w': ok_lvl1}
    return all([ok_exact, ok_Pzero, ok_P1, ok_K1, ok_Kline, ok_G, ok_lvl1])

# ------------------------------------------------------------------------------------------------
# C6 Konishi pair
# ------------------------------------------------------------------------------------------------
def C6():
    # basis (K, O_xi), xi in S_+ (x) W; the unit u is the expansion parameter; m_K, m_zeta, q, s, k are numbers
    u = sp.symbols('u', positive=True)
    nK, nz, mK, mz, tau, sig, k = sp.symbols('n_K n_zeta m_K m_zeta tau sigma k', real=True)
    qr, qi, sr, si, a_im = sp.symbols('q_r q_i s_r s_i a_i', real=True)
    q = qr + sp.I*qi; s = sr + sp.I*si; xi2 = sp.symbols('xi2', positive=True)   # (xi, conj xi) > 0
    # rho(xi) K = -(1 + u q) tau O_xi ;  rho(xi^dagger) O_zeta = +(1 + u s) sigma (zeta, conj xi) K
    X = sp.Matrix([[0, 0], [-(1 + u*q)*tau, 0]])            # matrix of -Lambda^(1)(xi)
    Y = sp.Matrix([[0, (1 + u*s)*sig*xi2], [0, 0]])          # matrix of -Lambda^(1)(xi^dagger) on (K, O_xi)
    G = sp.diag(nK*(1 + u*mK), nz*xi2*(1 + u*mz))            # radial Gram, Gram(O_xi,O_xi) = n_zeta (xi,conj xi)(1+u m_zeta)
    rel = sp.expand((X.H * G - G * Y)[0, 1])                 # exact adjointness conj(X)^T Gram = Gram Y
    r0 = rel.coeff(u, 0); r1 = rel.coeff(u, 1)
    nz_tree = sp.solve(r0, nz)[0]                            # tree adjointness
    # one-loop adjointness divided by the tree relation
    r1n = sp.simplify(sp.expand(r1.subs(nz, nz_tree)) / (sig*nK*xi2))
    # r1n = conj(q) + m_zeta - m_K - s, so r1 = 0 is s - conj(q) = m_zeta - m_K
    adj_ok = sp.simplify(r1n - (sp.conjugate(q) + mz - mK - s)) == 0
    # Maurer-Cartan: rho(xi^dagger) rho(xi) K + rho(xi) rho(xi^dagger) K = rho([xi, xi^dagger]) K,
    # [xi, xi^dagger] = -[xi, conj xi] has dilatation part -(1/2)(xi, conj xi) Delta; rho(Delta) K = -(4 + 2 u k) K
    # (D_1(K) = 2 u k K); rho(xi^dagger) K = 0, so (Y X)_{KK} = -(1/2) xi2 (-(4 + 2 u k)) = (2 + u k) xi2
    mc = sp.expand((Y * X)[0, 0] - (2 + u*k)*xi2)
    mc0 = mc.coeff(u, 0); mc1 = mc.coeff(u, 1)
    tau_sig = sp.simplify(sp.solve(mc0, sig)[0] * tau)
    mc1n = sp.simplify(mc1.subs(sig, -2/tau) / xi2)
    mc_ok = sp.simplify(mc1n - (2*(q + s) - k)) == 0         # q + s = k/2
    ratio = sp.simplify(nz_tree.subs(sig, -2/tau) / nK)
    a = sp.I*a_im
    GB1 = sp.simplify((s - sp.conjugate(q)) - ((s + a) - sp.conjugate(q - a))) == 0
    GB2 = sp.simplify((s + q) - ((s + a) + (q - a))) == 0
    # basis change psi = diag(pK, pz): q -> q + (pK - pz), s -> s + (pz - pK), m -> m + 2 Re p
    pK, pz = sp.symbols('p_K p_z')
    qn = q + (pK - pz); sn = s + (pz - pK)
    mKn = mK + pK + sp.conjugate(pK); mzn = mz + pz + sp.conjugate(pz)
    cov = sp.simplify((sn - sp.conjugate(qn) - (mzn - mKn)) - (s - sp.conjugate(q) - (mz - mK))) == 0
    OUT['C6_konishi_pair'] = {'classical_MC_tau_sigma': str(tau_sig), 'tree_ratio_n_zeta_over_n_K': str(ratio),
                              'one_loop_adjointness': 's - conj(q) = m_zeta - m_K : ' + str(adj_ok),
                              'one_loop_MC': 'q + s = k/2 (D_1(K) = 2 k unit K) : ' + str(mc_ok),
                              'GB_invariance': [GB1, GB2], 'basis_change_covariance': cov,
                              'Konishi_values': 'k = -24 => q + s = -12 ; tree ratio 1/2 for |tau| = 2'}
    return (sp.simplify(tau_sig + 2) == 0) and (sp.simplify(ratio - 2/tau**2) == 0) and adj_ok and mc_ok and GB1 and GB2 and cov

# ------------------------------------------------------------------------------------------------
# C7 two-component identities (Def c8-twocomp)
# ------------------------------------------------------------------------------------------------
def C7():
    s1 = np.array([[0, 1], [1, 0]], dtype=complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], dtype=complex)
    I2 = np.eye(2, dtype=complex)
    sig = [s1, s2, s3, 1j*I2]            # sigma^m_{alpha alphadot}
    sigb = [-s1, -s2, -s3, 1j*I2]        # sigmabar^{m alphadot alpha}
    res = {}
    res['conj_sigma_eq_minus_sigmabar_T'] = all(np.allclose(np.conj(sig[m]), -sigb[m].T) for m in range(4))
    res['clifford'] = all(np.allclose(sig[m] @ sigb[n] + sig[n] @ sigb[m], -2*(m == n)*I2) for m in range(4) for n in range(4))
    rng = np.random.default_rng(7)
    # elements of Sigma: xi[alpha, i], xib[alphadot, i]; Sigma^*: eta[alpha, i] (upper alpha lower i), etab[alphadot, i]
    def Gamma_V(xi1, xib1, xi2, xib2):   # Gamma^m(xi1,xi2)|_V = i(xib2 sigb^m xi1 + xib1 sigb^m xi2)
        return np.array([1j*(np.einsum('ai,ab,bi->', xib2, sigb[m], xi1) + np.einsum('ai,ab,bi->', xib1, sigb[m], xi2)) for m in range(4)])
    xi1 = rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4)); xi2 = rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4))
    z = np.zeros((2, 4), dtype=complex)
    res['Gamma_V_same_chirality_zero'] = np.allclose(Gamma_V(xi1, z, xi2, z), 0)
    # pairing (xi, eta) = xi_a^i eta^a_i - xib_{ad i} etab^{ad i}
    def pair(xi, xib, eta, etab):
        return np.einsum('ai,ai->', xi, eta) - np.einsum('ai,ai->', xib, etab)
    # rho(b) on xi (undotted only): (rho(b)xi)^{ad i} = -i b_m sigb^{m ad a} xi_a^i ; its undotted part 0
    b = rng.normal(size=4) + 1j*rng.normal(size=4)
    rb_dot = -1j*sum(b[m]*np.einsum('da,ai->di', sigb[m], xi2) for m in range(4))
    res['xi_rho_b_xiprime_zero_same_chirality'] = np.isclose(pair(xi1, z, z, rb_dot), 0)
    # conjugation: conj(xi) = (eta = conj(xi), etab = -conj(xib)); (xi, conj xi) = sum |xi|^2 + sum |xib|^2
    xib = rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4))
    val = pair(xi1, xib, np.conj(xi1), -np.conj(xib))
    res['pairing_positive'] = bool(np.isclose(val.imag, 0) and val.real > 0 and np.isclose(val.real, np.sum(abs(xi1)**2) + np.sum(abs(xib)**2)))
    # Clifford with the six transverse directions: build the 32x32 gamma matrices of V_10 abstractly
    # (only the anticommutation matters): sum_I rho(w_I) rho(x) rho(w_I) = -6 rho(x) for x in V
    g = [np.kron(np.kron(s1, np.eye(2)), np.eye(8))]
    # explicit Euclidean Clifford(10): use Jordan-Wigner
    def jw(k, n=5):
        mats = []
        for j in range(n):
            a = [s3]*j + [s1] + [np.eye(2)]*(n - j - 1); bmat = [s3]*j + [s2] + [np.eye(2)]*(n - j - 1)
            A = a[0]; B = bmat[0]
            for mm in a[1:]: A = np.kron(A, mm)
            for mm in bmat[1:]: B = np.kron(B, mm)
            mats += [A, B]
        return mats
    G10 = jw(0)
    ok_cl = all(np.allclose(G10[i] @ G10[j] + G10[j] @ G10[i], 2*(i == j)*np.eye(32)) for i in range(10) for j in range(10))
    xv = rng.normal(size=4)
    rx = sum(xv[m]*G10[m] for m in range(4))
    S6 = sum(G10[4 + I] @ rx @ G10[4 + I] for I in range(6))
    res['sum_rho_w_rho_x_rho_w_eq_minus6_rho_x'] = bool(ok_cl and np.allclose(S6, -6*rx))
    # Euclidean Fierz used for Fbar_w: sum_m sigb^{m ad a} sigb^{m bd b} = -2 eps^{ab} eps^{ad bd}
    eps = np.array([[0, 1], [-1, 0]])
    F = sum(np.einsum('ab,cd->acbd', sigb[m], sigb[m]) for m in range(4))   # indices (ad, bd, a, b)
    res['fierz_sigmabar'] = np.allclose(F, -2*np.einsum('ab,cd->acbd', eps, eps).transpose(0, 2, 1, 3).transpose(0, 2, 1, 3)) or \
        np.allclose(np.einsum('acbd->abcd', F), -2*np.einsum('ac,bd->abcd', eps, eps))
    OUT['C7_two_component'] = {k: bool(v) for k, v in res.items()}
    return all(res.values())

# ------------------------------------------------------------------------------------------------
# C8 weight-3 Konishi sector
# ------------------------------------------------------------------------------------------------
def C8():
    u, hb = sp.symbols('u hbar', positive=True)
    nF, nBp, kB, tau = sp.symbols('n_F n_Bp kappa_B tau', real=True)     # n_B = u n_Bp (relative order unit)
    bF, bB, cK, x, U = sp.symbols('b_F b_B c_K x U')
    # basis (Fbar_w, B_w); Gram_tree = diag(n_F, u n_Bp); minimal scheme: no (B, Fbar) entry at relative order unit
    G = sp.diag(nF, u*nBp)
    KON = sp.Matrix([u*tau*cK, tau*kB])          # Konishi weight-3 descendant: tau(kappa_B B + u c_K Fbar)
    PRO = sp.Matrix([bF, bB])                     # protected descendant b_F Fbar + b_B B
    ip = sp.expand((KON.H * G * PRO)[0, 0])
    lead = ip.coeff(u, 1)
    sol = sp.solve(sp.Eq(lead, 0), sp.conjugate(cK))
    pred = sp.simplify(sp.conjugate(sol[0]))
    # covariance: B_new = B + hbar x Fbar (finite part); in the new basis B = B_new - hbar x Fbar
    # coefficients of the same Konishi operator: B_new coefficient tau kB, Fbar coefficient u tau cK - hbar x tau kB
    # Gram in the new basis: (B_new, B_new) = u n_Bp + O(u^2), (Fbar, B_new) = hbar x n_F
    Gn = sp.Matrix([[nF, hb*x*nF], [hb*sp.conjugate(x)*nF, u*nBp]])
    KONn = sp.Matrix([u*tau*cK - hb*x*tau*kB, tau*kB])
    PROn = sp.Matrix([bF - hb*x*bB, bB])
    ipn = sp.expand((KONn.H * Gn * PROn)[0, 0])
    diff = sp.simplify(sp.expand(ipn - ip).subs(hb, u).coeff(u, 1))
    # order-sqrt(unit) relation, derived from the matrices: X = matrix of -Lambda^(1)(xi) on (O_zeta; Fbar, B),
    # X_{B zeta} = -kappa_B - u c_B', Y_{zeta B} = U; entry (O_zeta, B) of conj(X)^T Gram = Gram Y at leading order
    nz = sp.symbols('n_zeta', positive=True)
    cBp, kBc = sp.symbols('c_Bp kappa_Bc')                   # kappa_Bc: complex classical coefficient
    Xm = sp.zeros(3); Ym = sp.zeros(3)                       # basis (O_zeta, Fbar, B)
    Xm[2, 0] = -kBc - u*cBp; Xm[1, 0] = -u*cK
    U1 = sp.symbols('U1')                                    # U = u U1 (relative order unit in the B column)
    Ym[0, 2] = u*U1
    Gm = sp.diag(nz, nF, u*nBp)
    entry = sp.expand((Xm.H*Gm - Gm*Ym)[0, 2])
    Urel = u*sp.solve(sp.Eq(entry.coeff(u, 1), 0), U1)[0]   # leading (relative order sqrt(unit)) part
    Upred = -sp.conjugate(kBc)*u*nBp/nz                      # Prop ip-weight3 (i)
    ok_U = sp.simplify(Urel - Upred) == 0
    OUT['C8_weight3'] = {'c_K_from_orthogonality': str(pred), 'covariance_under_psi_B_eq_x_Fbar': str(diff),
                         'order_sqrt_unit_relation_U': str(Urel), 'U_equals_minus_conj_kappaB_nB_over_nzeta': ok_U}
    return sp.simplify(pred + kB*sp.conjugate(bB)*nBp/(sp.conjugate(bF)*nF)) == 0 and diff == 0 and ok_U

# ------------------------------------------------------------------------------------------------
# C9 rescalings of generators compatible with the bracket, gl(4|4) block model
# ------------------------------------------------------------------------------------------------
def C9():
    # blocks: rows/columns 1 = (2 dims), 2 = (2 dims), 3 = (4 dims) as in [SW (17)]
    sizes = [2, 2, 4]; off = [0, 2, 4]
    pairs = [(a, b) for a in range(3) for b in range(3)]
    r = {p: sp.Symbol('r_%d%d' % (p[0] + 1, p[1] + 1)) for p in pairs}
    # bracket of block (a,b) with block (b,c) lands in block (a,c) (block rule [SW (19)]); require r_ac = r_ab + r_bc
    eqs = []
    for (a, b) in pairs:
        for (b2, c) in pairs:
            if b2 == b:
                eqs.append(sp.Eq(r[(a, c)], r[(a, b)] + r[(b, c)]))
    sol = sp.solve(eqs, list(r.values()), dict=True)[0]
    free = sorted({str(s) for v in sol.values() for s in sp.sympify(v).free_symbols} | {str(k) for k in r.values() if k not in sol})
    rP = sp.simplify(r[(1, 0)].subs(sol)); rK = sp.simplify(r[(0, 1)].subs(sol))
    rQ = sp.simplify(r[(2, 0)].subs(sol)); rS = sp.simplify(r[(0, 2)].subs(sol))
    rQb = sp.simplify(r[(1, 2)].subs(sol)); rSb = sp.simplify(r[(2, 1)].subs(sol))
    diag_zero = all(sp.simplify(r[(a, a)].subs(sol)) == 0 for a in range(3))
    res = {'diagonal_blocks_and_Delta_unrescaled': diag_zero, 'r_K_plus_r_P': str(sp.simplify(rK + rP)),
           'r_Q_plus_r_S': str(sp.simplify(rQ + rS)), 'r_Qbar_plus_r_Sbar': str(sp.simplify(rQb + rSb)),
           'r_Q_plus_r_Qbar_minus_r_P': str(sp.simplify(rQ + rQb - rP)), 'free_parameters': free}
    # identify the two free directions: inner by Delta = diag(1,-1|0) and diag(0,0|1) (outer B up to sign)
    OUT['C9_rescalings'] = res
    return diag_zero and sp.simplify(rK + rP) == 0 and sp.simplify(rQ + rS) == 0 and sp.simplify(rQ + rQb - rP) == 0

# ------------------------------------------------------------------------------------------------
# C10 lemmas of the two-point design
# ------------------------------------------------------------------------------------------------
def C10():
    a, X = sp.symbols('a X', positive=True)
    F = 4**a * sp.gamma(2 + a) / (sp.pi**2 * sp.gamma(-a)) * X**(-a)      # times |x|^{-4}, X = |x|^2
    ser = sp.series(F, a, 0, 3).removeO()
    c1 = sp.simplify(ser.coeff(a, 1)); c2 = sp.simplify(2 * ser.coeff(a, 2))
    ok1 = sp.simplify(c1 + 1/sp.pi**2) == 0
    ok2 = sp.simplify(c2 - 2/sp.pi**2 * (sp.log(X) + 2*sp.EulerGamma - 1 - 2*sp.log(2))) == 0
    k = sp.symbols('k', positive=True)
    loc = sp.integrate((sp.exp(-3*k**2) - sp.exp(-2*k**2)) / k, (k, 0, sp.oo)) / (8*sp.pi**2)
    ok3 = sp.simplify(loc - sp.log(sp.Rational(2, 3)) / (16*sp.pi**2)) == 0
    OUT['C10_two_point_lemmas'] = {'FT_log': ok1, 'FT_log2': ok2, 'localisation_constant': ok3}
    return ok1 and ok2 and ok3


# ------------------------------------------------------------------------------------------------
# C11 Konishi anomalous dimension from the classical length-changing supercharge (Prop ip-weight3 (iv))
# ------------------------------------------------------------------------------------------------
def C11():
    import itertools
    eps3 = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        eps3[p] = np.linalg.det(np.eye(3)[list(p)])
    Sig = np.zeros((6, 4, 4), dtype=complex)        # 't Hooft blocks of ch08 Def c8-twocomp
    for a in range(3):
        for b in range(3):
            for c in range(3):
                Sig[a, b, c] = eps3[a, b, c]; Sig[a + 3, b, c] = 1j*eps3[a, b, c]
            Sig[a, b, 3] = (a == b); Sig[a, 3, b] = -(a == b)
            Sig[a + 3, b, 3] = -1j*(a == b); Sig[a + 3, 3, b] = 1j*(a == b)
    SigC = np.conj(Sig)                              # complex conjugate entries conj(Sigma^I_ij) (no new symbol)
    eps4 = np.zeros((4, 4, 4, 4))
    for p in itertools.permutations(range(4)):
        eps4[p] = np.linalg.det(np.eye(4)[list(p)])
    # A^{ij} = (1/2) eps^{ijkl} A_kl = A_I conj(Sigma^I_ij)
    ok_blocks = np.allclose(np.einsum('Iij,Ikl->ijkl', Sig, Sig), 2*eps4) and \
        np.allclose(0.5*np.einsum('ijkl,Ikl->Iij', eps4, Sig), SigC)
    rng = np.random.default_rng(11)
    vals = []
    for trial in range(4):
        w = rng.normal(size=4) + 1j*rng.normal(size=4)
        # B_z = ([A^{ik}, A_{kj} z^j], A_{il} z^l) = sum c_IJK ([A_I,A_J],A_K), c_IJK = conj(Sigma^I_ik) Sigma^J_kj z^j Sigma^K_il z^l
        c = np.einsum('Iik,Jkj,j,Kil,l->IJK', SigC, Sig, w, Sig, w)
        ca = np.zeros_like(c)
        for p in itertools.permutations(range(3)):
            ca += np.linalg.det(np.eye(3)[list(p)]) * np.transpose(c, p)
        ca /= 6
        vals.append(float(np.sum(abs(ca)**2) / np.vdot(w, w).real**2))
    flav = vals[0]
    ok_flav = np.allclose(vals, 48.0)
    # colour factor sum_abc ([e_a,e_b],e_c)^2 = -lambda dim g for su(2), lambda = tr(ad ad)/(X,Y)
    f = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        f[p] = np.linalg.det(np.eye(3)[list(p)])        # [e_a,e_b] = eps_abc e_c, (e_a,e_b) = delta_ab
    ad = np.array([[[f[a, b, c] for b in range(3)] for c in range(3)] for a in range(3)])   # (ad_a)_{cb}
    lam = np.trace(ad[0] @ ad[0])                         # tr(ad_X ad_X)/(X,X) with (X,X) = 1
    colour = float(np.sum(f**2)); ok_col = np.isclose(colour, -lam*3)
    # tree norms at unit separation (propagator scalar factor 1/(4 pi^2), A-A factor 2, hbar^(L-1)):
    # n_B = hbar^2 (2/(4pi^2))^3 (-lambda dim g) 6 * 48 |w|^4 ; n_zeta = 24 hbar dim g |u'|^2 |w|^2 / (16 pi^4)
    # unit = hbar lambda/(16 pi^2); kappa_B = 1/4 in modulus [odd design Prop o-konishi]
    hb, dimg = sp.symbols('hbar dim_g', positive=True); lamS = sp.symbols('lambda', negative=True)
    w2, u2, up2 = sp.symbols('w2 u2 up2', positive=True)      # |w|^2, |u|^2, |u'|^2
    nB = hb**2 * (2/(4*sp.pi**2))**3 * (-lamS*dimg) * 6 * sp.Rational(int(round(flav))) * w2**2
    nz = 24*hb*dimg*up2*w2/(16*sp.pi**4)
    unit = hb*lamS/(16*sp.pi**2)
    kapB2 = sp.Rational(1, 16)
    k = sp.simplify(kapB2*nB/(unit*(u2*w2)*nz))
    k_val = sp.simplify(k.subs(up2, 1/u2))                    # u'^a u_a = 1 and u' orthogonal to u: |u||u'| = 1
    OUT['C11_konishi_from_classical'] = {'thooft_blocks': bool(ok_blocks), 'flavour_sum_over_w4': flav,
                                         'colour_su2_equals_minus_lambda_dimg': bool(ok_col),
                                         'n_B_over_unit_n_zeta': str(sp.simplify(nB/(unit*nz))),
                                         'predicted_k (D_1 = 2 k unit)': str(k_val),
                                         'chapter6_value': -24}
    return bool(ok_blocks and ok_flav and ok_col and k_val == -24)

# ------------------------------------------------------------------------------------------------
# C12 two-point convention (coefficient of t2 t1), Ward signs, positivity of the radial form on O_zeta
# ------------------------------------------------------------------------------------------------
class _Ext:
    """exterior algebra over C on named generators with parities; element = dict{tuple(names): coeff}"""
    def __init__(self, terms, par):
        self.t, self.par = terms, par
    @staticmethod
    def gen(name, par):
        return _Ext({(name,): 1}, par)
    def __mul__(self, o):
        out = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                word = list(k1) + list(k2)
                odd = [x for x in word if self.par[x]]
                if len(set(odd)) < len(odd):
                    continue
                sign = 1
                for i in range(len(word)):                 # bubble sort with Koszul signs
                    for j in range(len(word) - 1 - i):
                        if word[j] > word[j + 1]:
                            if self.par[word[j]] and self.par[word[j + 1]]:
                                sign = -sign
                            word[j], word[j + 1] = word[j + 1], word[j]
                out[tuple(word)] = out.get(tuple(word), 0) + sign*v1*v2
        return _Ext({k: v for k, v in out.items() if v != 0}, self.par)
    def coeff(self, word):
        return self.t.get(tuple(word), 0)

def _realisation():
    """Def c8-twocomp: rho on Sigma (+) Sigma^* (32 components), pairing, conjugation of Lemma c8-cc."""
    s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
    sg = [s1, s2, s3, 1j*np.eye(2)]; sgb = [-s1, -s2, -s3, 1j*np.eye(2)]
    eu = np.array([[0, 1], [-1, 0]], complex); el = -eu          # eps^{12} = +1, eps_{12} = -1
    e3 = np.zeros((3, 3, 3)); e4 = np.zeros((4, 4, 4, 4))
    for p in itertools.permutations(range(3)): e3[p] = np.linalg.det(np.eye(3)[list(p)])
    for p in itertools.permutations(range(4)): e4[p] = np.linalg.det(np.eye(4)[list(p)])
    Sg = np.zeros((6, 4, 4), complex)
    for a in range(3):
        for b in range(3):
            for c in range(3):
                Sg[a, b, c] = e3[a, b, c]; Sg[a + 3, b, c] = 1j*e3[a, b, c]
            Sg[a, b, 3] = (a == b); Sg[a, 3, b] = -(a == b); Sg[a + 3, b, 3] = -1j*(a == b); Sg[a + 3, 3, b] = 1j*(a == b)
    # slots: Sigma = [xi_a^i (0..7), xibar_{ad i} (8..15)], Sigma^* = [eta^a_i (16..23), etabar^{ad i} (24..31)]
    X = lambda a, i: 4*a + i; XB = lambda a, i: 8 + 4*a + i; E = lambda a, i: 16 + 4*a + i; EB = lambda a, i: 24 + 4*a + i
    def rhoV(v):
        M = np.zeros((32, 32), complex)
        vs = sum(v[m]*sg[m] for m in range(4)); vsb = sum(v[m]*sgb[m] for m in range(4))
        for i in range(4):
            for a in range(2):
                for ad in range(2):
                    M[X(a, i), EB(ad, i)] += -1j*vs[a, ad]      # (rho(v)eta)_a^i = -i v sigma_{a ad} etabar^{ad i}
                    M[XB(ad, i), E(a, i)] += 1j*vs[a, ad]       # (rho(v)eta)_{ad i} = i v eta^a_i sigma_{a ad}
                    M[E(a, i), XB(ad, i)] += 1j*vsb[ad, a]      # (rho(v)xi)^a_i = i v xibar_{ad i} sigmabar^{ad a}
                    M[EB(ad, i), X(a, i)] += -1j*vsb[ad, a]     # (rho(v)xi)^{ad i} = -i v sigmabar^{ad a} xi_a^i
        return M
    def rhoW(w):
        M = np.zeros((32, 32), complex)
        wl = np.einsum('I,Iij->ij', w, Sg); wu = 0.5*np.einsum('ijkl,kl->ij', e4, wl)
        for a in range(2):
            for b in range(2):
                for j in range(4):
                    for i in range(4):
                        M[E(a, j), X(b, i)] += eu[a, b]*wl[j, i]        # (rho(w)xi)^a_j = eps^{ab} w_ji xi_b^i
                        M[EB(a, j), XB(b, i)] += eu[a, b]*wu[j, i]      # (rho(w)xi)^{ad j} = eps^{ad bd} w^{ji} xibar_{bd i}
                        M[X(a, j), E(b, i)] += -wu[j, i]*el[a, b]       # (rho(w)eta)_a^k = -w^{kj} eps_{ab} eta^b_j
                        M[XB(a, j), EB(b, i)] += -wl[j, i]*el[a, b]     # (rho(w)eta)_{bd k} = -w_kj eps_{bd ad} etabar^{ad j}
        return M
    G16 = np.diag([1.0]*8 + [-1.0]*8)
    B = np.zeros((32, 32)); B[:16, 16:] = G16; B[16:, :16] = G16          # symmetric form (xi,eta) on Sigma (+) Sigma^*
    def cc(z):                                                              # Lemma c8-cc, both directions
        out = np.zeros(32, complex)
        out[16:24] = np.conj(z[0:8]); out[24:32] = -np.conj(z[8:16])
        out[0:8] = np.conj(z[16:24]); out[8:16] = -np.conj(z[24:32])
        return out
    RV = [rhoV(np.eye(4)[m]) for m in range(4)]; RW = [rhoW(np.eye(6)[I]) for I in range(6)]
    return RV, RW, B, cc

def C12():
    res = {}
    # (a) Grassmann coefficients in (t1 O1)(t2 O2) for odd O1, O2 (t_a of the parity of O_a)
    par = {'a_t1': 1, 'b_t2': 1, 'c_O1': 1, 'd_O2': 1}
    t1, t2, O1, O2 = (_Ext.gen(n, par) for n in ['a_t1', 'b_t2', 'c_O1', 'd_O2'])
    prod = (t1*O1)*(t2*O2)
    c_t1t2 = prod.coeff(['a_t1', 'b_t2', 'c_O1', 'd_O2'])            # word t1 t2 O1 O2
    swap = (t2*t1*O1*O2)                                                 # t2 t1 O1 O2 in normal order
    c_t2t1 = prod.coeff(['a_t1', 'b_t2', 'c_O1', 'd_O2'])*swap.coeff(['a_t1', 'b_t2', 'c_O1', 'd_O2'])
    res['coefficient_of_t1t2'] = c_t1t2
    res['coefficient_of_t2t1'] = c_t2t1
    # (b) Ward signs.  E(delta(F G)) = 0 for a derivation delta of total parity y gives E(dF G) = s_E E(F dG),
    #     s_E = -(-1)^{y f}.  The t1t2 form is <A,B>_12 = (-1)^{|A||B|} E(AB), the t2t1 form <A,B>_21 = E(AB).
    table = {}
    for y in (0, 1):
        for f in (0, 1):
            g = (f + y) % 2
            sE = -(-1)**(y*f)
            s12 = (-1)**(((f + y) % 2)*g) * sE * (-1)**(f*((g + y) % 2))
            table['|Y|=%d,|F|=%d' % (y, f)] = {'t2t1_form': sE, 't1t2_form': s12, 'hypothesis_thm_c8_adjoint': -(-1)**(y*f)}
    res['ward_sign_table'] = table
    ok_ward = all(v['t2t1_form'] == v['hypothesis_thm_c8_adjoint'] for v in table.values()) and \
        any(v['t1t2_form'] != v['hypothesis_thm_c8_adjoint'] for v in table.values())
    # (c) free radial value on O_chi = sum_I (lambda, rho(w_I) chi) A_I with Theta_I (spinor factor -rho(x)/|x|)
    RV, RW, B, cc = _realisation()
    ok_cliff = all(np.allclose(Mi @ Mj + Mj @ Mi, 2*(i == j)*np.eye(32)) for i, Mi in enumerate(RV + RW) for j, Mj in enumerate(RV + RW))
    ok_sym = all(np.allclose(B @ M, (B @ M).T) for M in RV + RW)            # (rho(v) z1, z2) = (z1, rho(v) z2)
    rx = lambda x: sum(x[m]*RV[m] for m in range(4))
    S6 = lambda x: sum(RW[I] @ rx(x) @ RW[I] for I in range(6))
    xt = np.array([0.3, -0.7, 0.2, 0.5]); ok_m6 = np.allclose(S6(xt), -6*rx(xt))
    f_dual = [np.eye(32)[16 + s]*(1.0 if s < 8 else -1.0) for s in range(16)]   # lambda_s = (lambda, f_s)
    def E_ll(y):        # <lambda_s(y) lambda_t(0)> = (rho(d_y D(y)) f_s, f_t), D = 1/(4 pi^2 |y|^2)  [conj, prop:RP]
        dD = -2*y/(4*np.pi**2*(y @ y)**2)
        M = rx(dD)
        return np.array([[ (M @ f_dual[s]) @ B @ f_dual[t] for t in range(16)] for s in range(16)])
    def coeffs(chi):    # O_chi = sum_{s,I} T[s,I] lambda_s A_I : T[s,I] = component of rho(w_I)chi dual to lambda_s
        return np.array([[ (np.eye(32)[s] @ B @ (RW[I] @ chi)) for I in range(6)] for s in range(16)])
    def coeffs_theta(chi, x, s_lam=-1.0):   # (O_chi(x))^{Theta_I} = sum T'[s,I] lambda_s(Ix) A_I(Ix), F^Theta(e) = F(Theta e)^c
        xx = x @ x; T = np.zeros((16, 6), complex)
        for s in range(16):
            lam = np.zeros(32, complex); lam[s] = 1.0
            tl = s_lam*(rx(x) @ cc(lam))/xx**2                      # Theta_I lambda (unit field component s)
            for I in range(6):
                T[s, I] = np.conj((tl @ B @ (RW[I] @ chi)) / xx)    # A_I -> (x,x)^{-1} conj A_I(Ix)
        return T
    def expect(Ta, Tb, y):              # E(O_a(y) O_b(0)) / (hbar dim g): lambda contraction times 2 D for A_I A_J
        D = 1/(4*np.pi**2*(y @ y))
        return np.einsum('sI,tJ,st->', Ta, Tb, E_ll(y))*2*D
    np.seterr(all='ignore')          # Accelerate BLAS on this machine emits spurious FP flags in complex matmul
    rng = np.random.default_rng(20261029)
    ratios, ratios_plus, sep = [], [], []
    for trial in range(5):
        x = rng.normal(size=4); x = 0.85*rng.random()*x/np.linalg.norm(x)
        z = np.zeros(32, complex); z[:8] = rng.normal(size=8) + 1j*rng.normal(size=8)      # zeta in S_+ (x) W
        zp = np.zeros(32, complex); zp[:8] = rng.normal(size=8) + 1j*rng.normal(size=8)    # zeta'
        Ix = x/(x @ x)
        val = expect(coeffs_theta(z, x), coeffs(zp), Ix)
        ref = 24/(16*np.pi**4)*(zp @ B @ cc(z))                   # 24 hbar dim g (zeta', conj zeta)/(16 pi^4)
        ratios.append(val/ref)
        ratios_plus.append(expect(coeffs_theta(z, x, +1.0), coeffs(zp), Ix)/ref)
        # separated points: E(O_xi(y) O_chi(0)) = 24 hbar dim g (rho(y) xi, chi)/(16 pi^4 |y|^6), chi in S_- (x) W^*
        chi = np.zeros(32, complex); chi[8:16] = rng.normal(size=8) + 1j*rng.normal(size=8)
        y = rng.normal(size=4)
        sep.append(expect(coeffs(z), coeffs(chi), y)/(24/(16*np.pi**4)*((rx(y) @ z) @ B @ chi)/(y @ y)**3))
    res['realisation'] = {'clifford': bool(ok_cliff), 'rho_symmetric': bool(ok_sym), 'sum_rho_w_rho_x_rho_w_minus6': bool(ok_m6)}
    res['radial_O_zeta_over_24(zeta\',conj zeta)/(16pi^4)'] = [str(np.round(r, 10)) for r in ratios]
    res['radial_with_s_lambda_plus_1'] = [str(np.round(r, 10)) for r in ratios_plus]
    res['separated_E(O_xi(y)O_chi(0))_over_24(rho(y)xi,chi)/(16pi^4|y|^6)'] = [str(np.round(r, 10)) for r in sep]
    ok_pos = np.allclose(ratios, 1.0) and np.allclose(ratios_plus, -1.0) and np.allclose(sep, 1.0)
    res['conclusion'] = ('coefficient of t2 t1 = expectation; coefficient of t1 t2 = (-1)^{|O1||O2|} expectation; '
                         'radial form on O_zeta = +24 hbar dim g (zeta\',conj zeta)/(16 pi^4) in the t2t1 form, -24 in the t1t2 form')
    OUT['C12_two_point_convention'] = res
    return (c_t1t2 == -1) and (c_t2t1 == 1) and ok_ward and ok_cliff and ok_sym and ok_m6 and ok_pos

# ------------------------------------------------------------------------------------------------
# C13 radius of the quantisation sphere (Remark ip-radius)
# ------------------------------------------------------------------------------------------------
def C13():
    r, Rr, D, c = sp.symbols('r R Delta c', positive=True)
    # conjugation with R^2 I: conformal factor (R^2/r^2)^Delta, image at distance R^2/r; two-point c |z|^{-2 Delta}
    val = sp.simplify((Rr**2/r**2)**D * c * (Rr**2/r)**(-2*D))
    ok = sp.simplify(val - c*Rr**(-2*D)) == 0
    u, cc_, w = sp.symbols('u c_in w', real=True)
    # R^{-2w} = 1 - 2 w log R ; basis change 1 + c u w multiplies the Gram by 1 + 2 c u w  =>  log R = -c u
    logR = sp.solve(sp.Eq(-2*w*sp.Symbol('L'), 2*cc_*u*w), sp.Symbol('L'))[0]
    OUT['C13_radius'] = {'radial_Gram_radius_R_over_unit_sphere': str(sp.simplify(val/c)), 'ok': bool(ok),
                         'log_R_for_basis_change_1_plus_c_unit_w': str(logR)}
    return bool(ok) and sp.simplify(logR + cc_*u) == 0

# ------------------------------------------------------------------------------------------------
# C14 order-one-half entries between B_z (three scalars) and the two-fermion weight-3 bilinears
# ------------------------------------------------------------------------------------------------
def C14():
    # half-edge types: 'p' scalar A|_Wperp, 'v' A|_V, 'u' undotted lambda, 'd' dotted lambda; edges: pp, vv, ud
    cubic = {'vvv': ['v', 'v', 'v'], 'vpp': ['v', 'p', 'p'], 'uup': ['u', 'u', 'p'], 'ddp': ['d', 'd', 'p'], 'udv': ['u', 'd', 'v']}
    allowed = {frozenset(['p']), frozenset(['v']), frozenset(['u', 'd'])}
    def matchings(items):
        if not items:
            yield []
            return
        a = items[0]
        for k in range(1, len(items)):
            b = items[k]
            if frozenset([a[1], b[1]]) in allowed and not (a[1] == b[1] and a[1] in 'ud'):
                rest = items[1:k] + items[k + 1:]
                for m in matchings(rest):
                    yield [(a, b)] + m
    report = {}
    total_surviving = 0
    for fermion in ('u', 'd'):                         # F-type (two undotted) or Fbar-type (two dotted)
        for name, legs in cubic.items():
            items = [('B', 'p', 0), ('B', 'p', 1), ('B', 'p', 2), ('F', fermion, 0), ('F', fermion, 1)] + \
                    [('V', t, i) for i, t in enumerate(legs)]
            items = [(it[0] + str(it[2]), it[1], it[0]) for it in items]
            conn = surv = 0
            for m in matchings(items):
                # connectivity of the vertex graph {B, F, V}
                adj = {'B': set(), 'F': set(), 'V': set()}
                selfc = False
                for (a, b) in m:
                    if a[2] == b[2]:
                        selfc = True
                    adj[a[2]].add(b[2]); adj[b[2]].add(a[2])
                seen = {'B'}; stack = ['B']
                while stack:
                    n = stack.pop()
                    for q in adj[n]:
                        if q not in seen:
                            seen.add(q); stack.append(q)
                if len(seen) < 3:
                    continue
                conn += 1
                if not selfc:        # a self-contraction of B (bracket) or of the cubic vertex gives colour 0
                    surv += 1
            report['%s-fermions, vertex %s' % (fermion, name)] = {'connected': conn, 'without_self_contraction': surv}
            total_surviving += surv
    OUT['C14_order_half_census'] = {'classes': report, 'surviving_graphs': total_surviving,
                                    'half_edges': 3 + 2 + 3,
                                    'statement': 'every connected one-cubic-vertex graph has a self-contraction of B_z (colour 0) or does not exist'}
    return total_surviving == 0

# ------------------------------------------------------------------------------------------------
# C15 frames: route (a) leaves the two-point coefficients unchanged; mixed-frame residuals
# ------------------------------------------------------------------------------------------------
def C15():
    w, kap, c, u = sp.symbols('w kappa c u', real=True)
    N = 4
    h = w + u*kap
    P = sp.zeros(N); K = sp.zeros(N)
    for n in range(N - 1): P[n + 1, n] = 1
    for n in range(1, N): K[n - 1, n] = n*(2*h + n - 1)
    GZ = sp.diag(*[sp.factorial(n)*sp.rf(2*h, n) for n in range(N)])       # frame P_1 = 0: radial Gram = two-point coefficient
    ser = lambda M, k: M.applyfunc(lambda e: sp.expand(sp.series(e, u, 0, 2).removeO()).coeff(u, k))
    trunc = lambda M: M.applyfunc(lambda e: sp.expand(sp.series(sp.expand(e), u, 0, 2).removeO()))
    Wt = sp.diag(*[w + n for n in range(N)])                                 # Weyl weight on level n
    U = sp.eye(N) + c*u*Wt; Ui = sp.eye(N) - c*u*Wt                          # route (a) at the origin
    PE, KE = trunc(U*P*Ui), trunc(U*K*Ui)                                     # engines' frame: P_1 = c u P_0
    GE = trunc(Ui.T*GZ*Ui)                                                    # radial Gram in that frame
    two_point = GZ                                                            # unchanged by route (a)
    resid = lambda Pm, Km, G: trunc(Pm.T*G - G*Km)
    ok_E = all(sp.simplify(e) == 0 for e in resid(PE, KE, GE))
    ok_Pscale = all(sp.simplify(e) == 0 for e in trunc(PE - (1 + c*u)*P))
    ok_GE = all(sp.simplify(e) == 0 for e in trunc(GE - (sp.eye(N) - 2*c*u*Wt)*GZ))
    mixed = resid(PE, KE, two_point)                                           # engines' matrices with two-point Grams
    ok_mixed_nonzero = any(sp.simplify(e) != 0 for e in mixed)
    # IP6 as previously written: frame-Z matrices with the two-point Gram multiplied by (1 + 2 c u w)
    G_ip6 = trunc((sp.eye(N) + 2*c*u*Wt)*two_point)
    r_ip6 = resid(P, K, G_ip6)
    K1_res = sp.simplify(trunc(G_ip6.inv()*P.T*G_ip6 - K))
    OUT['C15_frames'] = {'frame_E_relation_with_frame_E_Gram': ok_E, 'P_1_equals_c_unit_P_0': ok_Pscale,
                         'Gram_E_equals_(1-2 c unit w)_two_point': ok_GE,
                         'mixed_frames_residual_nonzero': ok_mixed_nonzero,
                         'IP6_double_rescaling_K_residual_entry_01': str(sp.factor(ser(K1_res, 1)[0, 1])),
                         'K_0_entry_01': str(K[0, 1].subs(u, 0))}
    # Konishi pair (A1 of the audit): frame-Z data satisfying T2 exactly, transported by route (a)
    qr, qi, sr, si, mK, mz, cK, cz = sp.symbols('q_r q_i s_r s_i m_K m_zeta c_K c_zeta', real=True)
    q = qr + sp.I*qi; s = sr + sp.I*si
    sZ = sp.conjugate(q) + mz - mK                                            # T2 radial line in frame Z
    qE = q + sp.Rational(5, 2)*cz - 2*cK; sE = sZ - sp.Rational(5, 2)*cz + 2*cK   # X_E = U X_Z U^{-1}
    res_mixed = sp.simplify(sE - sp.conjugate(qE) - (mz - mK))                 # combined with two-point constants
    res_univ = sp.simplify(res_mixed.subs({cK: c, cz: c}))
    OUT['C15_frames']['T2_residual_mixed_frames'] = str(res_mixed)
    OUT['C15_frames']['T2_residual_universal_c'] = str(res_univ)
    return ok_E and ok_Pscale and ok_GE and ok_mixed_nonzero and sp.simplify(res_mixed - (4*cK - 5*cz)) == 0 \
        and sp.simplify(res_univ + c) == 0

# ------------------------------------------------------------------------------------------------
# C16 reflection compatibility of universal real rescalings (block model of [SW, (17)])
# ------------------------------------------------------------------------------------------------
def C16():
    s1, s2, s3 = sp.symbols('s1 s2 s3', real=True)
    r = lambda a, b: [s1, s2, s3][a - 1] - [s1, s2, s3][b - 1]
    blocks = {'P': (2, 1), 'K': (1, 2), 'Q': (3, 1), 'Qbar': (2, 3), 'S': (1, 3), 'Sbar': (3, 2)}
    # Theta_n exchanges Q <-> Qbar and S <-> Sbar, fixes P_n, K_n; real rescalings: r(X^{dagger_n}) = r(X)
    sol = sp.solve([sp.Eq(r(*blocks['Q']), r(*blocks['Qbar'])), sp.Eq(r(*blocks['S']), r(*blocks['Sbar']))], [s3], dict=True)[0]
    rQ = sp.simplify(r(*blocks['Q']).subs(sol)); rP = sp.simplify(r(*blocks['P']).subs(sol))
    inner = {k: str(r(*v).subs({s1: 1, s2: -1, s3: 0})) for k, v in blocks.items()}     # ad_Delta, Delta = diag(1,-1|0)
    outer = {k: str(r(*v).subs({s1: 0, s2: 0, s3: 1})) for k, v in blocks.items()}      # B = diag(0,0|1)
    cond = lambda v: sp.simplify(([s1, s2, s3][2] - (s1 + s2)/2).subs(v)) == 0   # s3 = (s1+s2)/2
    inner_ok = cond({s1: 1, s2: -1, s3: 0}); outer_ok = cond({s1: 0, s2: 0, s3: 1})
    # imaginary rescalings: conj(r(X)) = r(X^{dagger_n}) with r purely imaginary gives s1 = s2 (r_P = 0), s3 free
    t1, t2, t3 = sp.symbols('t1 t2 t3', real=True)
    ri = lambda a, b: sp.I*([t1, t2, t3][a - 1] - [t1, t2, t3][b - 1])
    soli = sp.solve([sp.Eq(sp.conjugate(ri(3, 1)), ri(2, 3)), sp.Eq(sp.conjugate(ri(1, 3)), ri(3, 2))], [t1], dict=True)[0]
    rPi = sp.simplify(ri(2, 1).subs(soli))
    OUT['C16_reflection_rescalings'] = {'reflection_condition_real': str(sol), 'r_Q_over_r_P': str(sp.simplify(rQ/rP)),
                                        'inner_Delta_pattern': inner, 'outer_B_pattern': outer,
                                        'inner_real_reflection_compatible': bool(inner_ok),
                                        'outer_real_reflection_compatible': bool(outer_ok),
                                        'imaginary_rescalings_condition': str(soli), 'imaginary_r_P': str(rPi)}
    return sp.simplify(rQ - rP/2) == 0 and inner_ok and not outer_ok and rPi == 0

# ------------------------------------------------------------------------------------------------
# C17 weight 3: D_1 = D_1^adj is the orthogonality of the Konishi and protected descendants
# ------------------------------------------------------------------------------------------------
def C17():
    u = sp.symbols('u', positive=True)
    nF, nBp = sp.symbols('n_F n_Bp', positive=True)
    kB, cK, bF, bB = sp.symbols('kappa_B c_K b_F b_B')
    G = sp.diag(nF, u*nBp)                                   # basis (Fbar, B)
    k = sp.Matrix([u*cK, kB]); p = sp.Matrix([bF, bB])     # Konishi and protected descendants
    T = k.row_join(p)
    D1 = T*sp.diag(-48*u, 0)*T.inv()                         # D_1 = -48 unit on k, 0 on p
    A = sp.simplify(G*D1 - D1.H*G)                           # anti-self-adjoint part (Gram D_1 Hermitian iff A = 0)
    orth = sp.expand((k.H*G*p)[0, 0]).coeff(u, 1)            # leading part of <Konishi, protected>
    cK_sol = sp.solve(sp.Eq(orth, 0), sp.conjugate(cK))[0]
    A_on = sp.simplify(A.subs(sp.conjugate(cK), cK_sol).subs(cK, sp.conjugate(cK_sol)))
    lead = A_on.applyfunc(lambda e: sp.limit(sp.simplify(e/u**2), u, 0))
    OUT['C17_T5_vs_T3iii'] = {'orthogonality_solution_conj_c_K': str(cK_sol),
                              'Gram_D1_antiHermitian_part_leading_at_orthogonality': str(lead)}
    return all(sp.simplify(e) == 0 for e in lead)

if __name__ == '__main__':
    summary = {}
    for name, fn in [('C1', C1), ('C2', C2), ('C3', C3), ('C4', C4), ('C5', C5), ('C6', C6), ('C7', C7),
                     ('C8', C8), ('C9', C9), ('C10', C10), ('C11', C11), ('C12', C12), ('C13', C13), ('C14', C14),
                     ('C15', C15), ('C16', C16), ('C17', C17)]:
        try:
            summary[name] = bool(fn())
        except Exception as e:  # record and continue
            summary[name] = 'error: ' + repr(e)
        print(name, summary[name], flush=True)
    OUT['summary'] = summary
    with open(os.path.join(HERE, 'inner_product_checks.json'), 'w') as fh:
        json.dump(OUT, fh, indent=1, default=str)
