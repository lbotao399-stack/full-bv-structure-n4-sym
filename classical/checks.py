"""
checks.py -- machine checks of the classical BV data (writes checks.json).

Sections
  1  explicit model: 10d Clifford blocks, Fierz [ESW Ass. 3.1], su(2)/su(3) structure constants
  2  CME {S,S} = 0 for the full reduced gauge-BV action (su(2), jets), bracket fixed by C07 Sec 10.1
     normalisation; negative controls; free part (S2,S2) (Q^2 = 0)
  3  ESW Lemma 3.4 for the reduced theory with ESW (33): {S, S^(1)(Q)} = 0 (fixes (F,F) normalisation)
  4  free complex in momentum space: Q from S2, Q^GF solved from H = Laplacian, H, [Q^GF,H], (Q^GF)^2,
     adjointness, Im Q cap Im Q^GF
  5  K_t, P(eps,infty): graded symmetry, Q P = K_eps, heat-kernel identities
  6  Hessian inversion (lambda, transverse A, c-A^* exact sector), positivity of A-A covariance, real slice
  7  position space: G_eps, Laplacian identity, derivative kernel
  8  vertex functions (test-field evaluation of I) vs index-free formulas
  9  0-dim check of the expansion convention of W(P,I) (Wick rule and one-loop trace formula)
"""
import sys, os, json, time, random, itertools
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
import numpy as np
import n4_model as nm
import jets as J
import action as Ac

RESULTS = {}
EPS_BRACKET = {'A': 1, 'lam': -1, 'c': -1}   # (phi,phi^*)_eps = eps_phi/kappa : C07 Sec 10.1 even-first


def rec(key, val):
    RESULTS[key] = val
    print(key, '=>', val if not isinstance(val, dict) else json.dumps(val)[:400])


# ======================================================================================
# 1. explicit model
# ======================================================================================
def section1():
    g = nm.gamma_tests()
    sig = [nm.sig_sym(M) for M in range(10)]
    sigb = [nm.sigbar_sym(M) for M in range(10)]
    Om4 = sigb[0] * sig[1] * sigb[2] * sig[3]
    Om6 = sigb[4] * sig[5] * sigb[6] * sig[7] * sigb[8] * sig[9]
    prod = sp.simplify(Om6 * Om4)
    g['Omega6_times_Omega4'] = str(prod[0, 0]) + ' * identity' if (prod - prod[0, 0] * sp.eye(16)).is_zero_matrix else 'not scalar'
    # 4d chirality on Sigma^* = S_-^{10}: Om4s = sig^0 sigbar^1 sig^2 sigbar^3 ; rho(v) for v in V flips, for w in W_perp keeps
    Om4s = sig[0] * sigb[1] * sig[2] * sigb[3]
    g['Om4s_squared_is_1'] = bool((Om4s * Om4s - sp.eye(16)).is_zero_matrix)
    g['rho(e_mu), mu in V: sig^mu Om4 = -Om4s sig^mu (maps S_+(x)W -> S_-(x)W type, flips 4d chirality)'] = all(
        (sig[m] * Om4 + Om4s * sig[m]).is_zero_matrix for m in nm.V4)
    g['rho(e_w), w in Wperp: sig^w Om4 = +Om4s sig^w (keeps 4d chirality)'] = all(
        (sig[m] * Om4 - Om4s * sig[m]).is_zero_matrix for m in nm.WPERP)
    rec('1_gamma', g)
    rec('1_su2', nm.lie_tests('su2'))
    rec('1_su3', nm.lie_tests('su3'))


# ======================================================================================
# 2. CME
# ======================================================================================
def section2():
    out = {}
    V = J.Vars(nm.PARITY)
    d, f = nm.f_int('su2')
    parts = Ac.build_action(V, f, d)
    S = Ac.total(parts)
    fam = Ac.families(d)
    out['n_terms_density_x8'] = {k: len(v) for k, v in parts.items()}
    t = time.time()
    R = J.antibracket_all(V, S, S, nm.ANTI, fam, EPS_BRACKET)
    ok, bad = J.total_derivative_test(V, R)
    out['CME_{S,S}_terms_before_integration'] = len(R)
    out['CME_{S,S}=0_mod_total_derivatives'] = ok
    out['CME_seconds'] = round(time.time() - t, 2)
    # equivalent split: 2 Q I + {I,I}
    I = J.padd(parts['I3'], parts['I4'])
    S2 = parts['S2']
    R2 = J.antibracket_all(V, S2, S2, nm.ANTI, fam, EPS_BRACKET)
    out['free_{S2,S2}=0 (Q^2=0)'] = J.total_derivative_test(V, R2)[0]
    QI = J.padd(J.antibracket_all(V, S2, I, nm.ANTI, fam, EPS_BRACKET), J.antibracket_all(V, I, S2, nm.ANTI, fam, EPS_BRACKET))
    II = J.antibracket_all(V, I, I, nm.ANTI, fam, EPS_BRACKET)
    out['QI_nonzero_density'] = len(QI) > 0
    out['2QI+{I,I}=0'] = J.total_derivative_test(V, J.padd(QI, II))[0]
    # orders separately: cubic order Q I3 = 0, quartic order Q I4 + 1/2{I3,I3} = 0, ...
    # negative controls: other relative signs of the odd pairs
    neg = {}
    for name, eps in [('uniform', {'A': 1, 'lam': 1, 'c': 1}), ('lam_only', {'A': 1, 'lam': -1, 'c': 1}),
                      ('c_only', {'A': 1, 'lam': 1, 'c': -1})]:
        Rn = J.antibracket_all(V, S, S, nm.ANTI, fam, eps)
        neg[name] = J.total_derivative_test(V, Rn)[0]
    out['negative_control_other_pair_signs_CME_holds'] = neg
    # negative control: flip sign of 1/2([c,c],c^*)
    Sbad = dict(S)
    for k in list(Sbad):
        if any(V.key[v][0] == 'cs' for v in k[0] + k[1]):
            Sbad[k] = J.cneg(Sbad[k])
    out['negative_control_flip_ccc*_CME_holds'] = J.total_derivative_test(V, J.antibracket_all(V, Sbad, Sbad, nm.ANTI, fam, EPS_BRACKET))[0]
    rec('2_CME_su2', out)
    return V, parts


# ======================================================================================
# 3. ESW Lemma 3.4 with (33) in the reduced theory
# ======================================================================================
def susy_S1(V, q, fint, dimg, bos_sign=1):
    """8 * S^(1)(Q) = 8*[ -(Gamma(Q,lam),A^*) + 1/2 (rho(F_A) Q, lam^*) ],  Q = q (commuting, Gaussian ints)."""
    out = {}
    fl = list(fint.items())

    def qc(x):
        return x
    # Gamma(Q,lam)^{M,i} = sum_ab sig^M_ab q^b lam^{a,i}
    for M in range(10):
        for (a, b, s) in nm.sig_entries(M):
            coef = J.cmul(s, q[b])
            if coef == (0, 0):
                continue
            for i in range(dimg):
                J.padd_into(out, J.mono(V, (-8 * coef[0], -8 * coef[1]), [V.get('lam', a, i), V.get('As', M, i)]))
    # rho(F) Q = sum_{M<N} F_MN sigbar^M sig^N q ;  F_MN^i = d_M A_N - d_N A_M + f_jki A_M^j A_N^k
    sig = [nm.SIG_RE[M] + 1j * nm.SIG_IM[M] for M in range(10)]
    sigb = [nm.SIGBAR_RE[M] + 1j * nm.SIGBAR_IM[M] for M in range(10)]
    qv = np.array([complex(*x) for x in q])
    for M in range(10):
        for N in range(M + 1, 10):
            vec = sigb[M] @ (sig[N] @ qv)       # exact small Gaussian integers in floats
            for aa in range(16):
                z = vec[aa]
                zr, zi = int(round(z.real)), int(round(z.imag))
                assert abs(z.real - zr) < 1e-9 and abs(z.imag - zi) < 1e-9
                if (zr, zi) == (0, 0):
                    continue
                for i in range(dimg):
                    FM = {}
                    if M in nm.V4:
                        J.padd_into(FM, J.var(V, 'A', N, i, Ac.e_mu(M)))
                    if N in nm.V4:
                        J.padd_into(FM, J.var(V, 'A', M, i, Ac.e_mu(N)), (-1, 0))
                    for (j, k, l), val in fl:
                        if l == i:
                            J.padd_into(FM, J.mono(V, (val, 0), [V.get('A', M, j), V.get('A', N, k)]))
                    term = J.pmul(FM, J.var(V, 'lams', aa, i))
                    J.padd_into(out, term, (4 * zr * bos_sign, 4 * zi * bos_sign))
    return out


def section3():
    out = {}
    V = J.Vars(nm.PARITY)
    d, f = nm.f_int('su2')
    parts = Ac.build_action(V, f, d)
    S = Ac.total(parts)
    fam = Ac.families(d)
    rng = random.Random(7)
    oks = []
    for trial in range(2):
        q = [(rng.randint(-2, 2), rng.randint(-2, 2)) for _ in range(16)]
        S1 = susy_S1(V, q, f, d)
        R = J.antibracket_all(V, S, S1, nm.ANTI, fam, EPS_BRACKET)
        oks.append(J.total_derivative_test(V, R)[0])
    out['{S,S1(Q)}=0 (ESW Lemma 3.4, su2, 2 random Q)'] = all(oks)
    # free abelian normalisation scan: coefficient a of (F,F)_Hodge
    V2 = J.Vars(nm.PARITY)
    p2 = Ac.build_action(V2, {}, 1, include=('S2',))['S2']
    bos = {k: v for k, v in p2.items() if all(V2.key[x][0] == 'A' for x in k[0] + k[1])}
    rest = {k: v for k, v in p2.items() if k not in bos}
    q = [(rng.randint(-2, 2), rng.randint(-2, 2)) for _ in range(16)]
    S1 = susy_S1(V2, q, {}, 1)
    scan = {}
    for amult, label in [(1, 'a=-1/4 (ESW)'), (-1, 'a=+1/4'), (2, 'a=-1/2 (i.e. -1/4 sum_{all M,N})')]:
        Sx = J.padd(rest, bos, (amult, 0))
        R = J.antibracket_all(V2, Sx, S1, nm.ANTI, Ac.families(1), EPS_BRACKET)
        scan[label] = J.total_derivative_test(V2, R)[0]
    out['free_abelian_SUSY_invariance_vs_coefficient_of_(F,F)_Hodge'] = scan
    # with uniform bracket instead
    R = J.antibracket_all(V, S, susy_S1(V, q, f, d), nm.ANTI, fam, {'A': 1, 'lam': 1, 'c': 1})
    out['negative_control_uniform_bracket_Lemma3.4_holds'] = J.total_derivative_test(V, R)[0]
    rec('3_ESW_Lemma3.4', out)


# ======================================================================================
# 4-6. free complex in momentum space
# ======================================================================================
p = sp.symbols('p0:4', real=True)
kap, eps_, t_ = sp.symbols('kappa epsilon t', positive=True)
P2 = sum(x ** 2 for x in p)

COMP = nm.component_list(1)                    # 54 components, single colour
IDX = {(f, a): n for n, (f, a, i) in enumerate(COMP)}
NC = len(COMP)
PAR = [nm.PARITY[f] for (f, a, i) in COMP]
Dg = sp.diag(*[(-1) ** x for x in PAR])


def jet_to_symbol(J4, mom):
    r = sp.Integer(1)
    for mu in range(4):
        r *= (sp.I * mom[mu]) ** J4[mu]
    return r


def linear_to_row(V, P, mom):
    row = [sp.Integer(0)] * NC
    for (e, o), c in P.items():
        vs = e + o
        assert len(vs) == 1
        f, a, i, J4 = V.key[vs[0]]
        row[IDX[(f, a)]] += (sp.Integer(c[0]) + sp.I * c[1]) * jet_to_symbol(J4, mom) / Ac.SCALE
    return row


def free_data():
    V = J.Vars(nm.PARITY)
    S2 = Ac.build_action(V, {}, 1, include=('S2',))['S2']
    Er = J.euler_all(V, S2, 'r')
    Q = sp.zeros(NC, NC)
    for (f, a, i) in COMP:
        row_idx = IDX[(f, a)]
        if f in nm.ANTI.values():
            phi = [k for k, v in nm.ANTI.items() if v == f][0]
            P = Er.get((phi, a, 0), {})
            row = linear_to_row(V, P, p)
            sgn = EPS_BRACKET[phi]
        else:
            P = Er.get((nm.ANTI[f], a, 0), {})
            row = linear_to_row(V, P, p)
            sgn = -EPS_BRACKET[f]
        for n in range(NC):
            Q[row_idx, n] = sp.expand(sgn * row[n])
    # Hessian (momentum symbol) of S2:  S2 = sum n_ab(ik) ;  N_ab(k) = n_ab(k) + (-1)^{|a||b|} n_ba(-k)
    nmat = sp.zeros(NC, NC)
    for (e, o), c in S2.items():
        vs = e + o
        assert len(vs) == 2
        (f1, a1, _, J1), (f2, a2, _, J2) = V.key[vs[0]], V.key[vs[1]]
        val = (sp.Integer(c[0]) + sp.I * c[1]) / Ac.SCALE * kap
        val *= jet_to_symbol(J1, [-x for x in p]) * jet_to_symbol(J2, p)
        nmat[IDX[(f1, a1)], IDX[(f2, a2)]] += val
    nmat_m = nmat.subs({x: -x for x in p}, simultaneous=True)
    N = sp.zeros(NC, NC)
    for a in range(NC):
        for b in range(NC):
            N[a, b] = sp.expand(nmat[a, b] + (-1) ** (PAR[a] * PAR[b]) * nmat_m[b, a])
    return Q, N


def K0_matrix():
    K = sp.zeros(NC, NC)
    for (f, a, i) in COMP:
        if f in nm.ANTI:
            fs = nm.ANTI[f]
            ia, ib = IDX[(f, a)], IDX[(fs, a)]
            e = EPS_BRACKET[f]
            # {phi,phi*} = (-1)^{|phi|} (phi,phi*)_eps = (-1)^{|phi|} e/kappa ; {phi*,phi} = (-1)^{|phi*|}(-e/kappa)
            br_ab = (-1) ** PAR[ia] * e / kap
            br_ba = (-1) ** PAR[ib] * (-e) / kap
            K[ia, ib] = -br_ab
            K[ib, ia] = -br_ba
    return K


def neg_p(M):
    return M.subs({x: -x for x in p}, simultaneous=True)


def QGF_ansatz():
    al, be, ga, de = sp.symbols('alpha beta gamma delta')
    G = sp.zeros(NC, NC)
    for M in range(10):
        G[IDX[('A', M)], IDX[('As', M)]] = be                     # A^* -> A (all ten)
    for mu in range(4):
        G[IDX[('c', 0)], IDX[('A', mu)]] = ga * sp.I * p[mu]       # A -> c : gamma * d_mu A_mu
        G[IDX[('As', mu)], IDX[('cs', 0)]] = de * sp.I * p[mu]     # c^* -> A^* : delta * d_mu c^*
    for mu in range(4):
        sb = nm.sigbar_sym(mu)
        for a in range(16):
            for b in range(16):
                if sb[a, b] != 0:
                    G[IDX[('lam', a)], IDX[('lams', b)]] += al * sb[a, b] * sp.I * p[mu]   # lam^* -> lam
    return G, (al, be, ga, de)


def section4_6():
    out4, out5, out6 = {}, {}, {}
    Q, N = free_data()
    out4['Q^2=0'] = (Q * Q).applyfunc(sp.expand).is_zero_matrix
    # describe Q blocks
    def blk(M, f1, f2, mu_list=None):
        return M[IDX[(f1, 0)]:IDX[(f1, 0)] + nm.FDIM[f1], IDX[(f2, 0)]:IDX[(f2, 0)] + nm.FDIM[f2]]
    out4['Q(c->A)'] = str(blk(Q, 'A', 'c').T)
    out4['Q(A^*->c^*)'] = str(blk(Q, 'cs', 'As'))
    QAA = blk(Q, 'As', 'A')
    out4['Q(A->A^*) = -1/2 (p^2 delta - p p^T on V)'] = sp.simplify(QAA + sp.Rational(1, 2) * (P2 * sp.eye(10) - sp.Matrix(10, 1, lambda m, _: p[m] if m < 4 else 0) * sp.Matrix(1, 10, lambda _, m: p[m] if m < 4 else 0))).is_zero_matrix
    Ql = blk(Q, 'lams', 'lam')
    dsl = sum((nm.sig_sym(mu) * sp.I * p[mu] for mu in range(4)), sp.zeros(16, 16))
    out4['Q(lam->lam^*) = +dslash (sig^mu i p_mu)'] = sp.simplify(Ql - dsl).is_zero_matrix
    out4['Q nonzero blocks only c->A, A->A^*, A^*->c^*, lam->lam^*'] = all(
        Q[a, b] == 0 for a in range(NC) for b in range(NC)
        if (COMP[b][0], COMP[a][0]) not in [('c', 'A'), ('A', 'As'), ('As', 'cs'), ('lam', 'lams')])
    # solve for Q^GF
    G, syms = QGF_ansatz()
    H = (Q * G + G * Q).applyfunc(sp.expand)
    eqs = set()
    for a in range(NC):
        for b in range(NC):
            e = sp.expand(H[a, b] - (P2 if a == b else 0))
            if e != 0:
                for cf in sp.Poly(e, *p).coeffs():
                    eqs.add(cf)
    sol = sp.solve(list(eqs), syms, dict=True)
    out4['QGF_coefficients_solved_from_H=Laplacian'] = {str(k): str(v) for k, v in sol[0].items()} if sol else 'NO SOLUTION'
    assert len(sol) == 1
    G = G.subs(sol[0])
    H = (Q * G + G * Q).applyfunc(sp.expand)
    out4['H=[Q,QGF]=|p|^2*id_on_all_54_components'] = (H - P2 * sp.eye(NC)).applyfunc(sp.expand).is_zero_matrix
    out4['[QGF,H]=0'] = (G * H - H * G).applyfunc(sp.expand).is_zero_matrix
    G2 = (G * G).applyfunc(sp.expand)
    nz = {}
    for a in range(NC):
        for b in range(NC):
            if G2[a, b] != 0:
                nz.setdefault(f'{COMP[b][0]}->{COMP[a][0]}', []).append(f'({COMP[b][1]}->{COMP[a][1]}): {G2[a, b]}')
    out4['(QGF)^2_nonzero_components'] = {k: v[:4] + (['...'] if len(v) > 4 else []) for k, v in nz.items()}
    out4['(QGF)^2!=0'] = not G2.is_zero_matrix
    K0 = K0_matrix()
    out4['K0_entries'] = sorted(set(f'{COMP[a][0]}(x){COMP[b][0]}: {K0[a, b]}' for a in range(NC) for b in range(NC) if K0[a, b] != 0))
    # adjointness: (Q(x)1 + 1(x)Q) K0 = 0  ;  (QGF(x)1 - 1(x)QGF) K0 = 0   [C07 Lemma 6.3.1, Def 4.3.3(1), Sec 4.3]
    def one_x(Xm, K):    # (1 (x) X) K in momentum space: (-1)^{|a|} sum_b X^b_beta(-p) K^{a beta}
        return Dg * K * neg_p(Xm).T
    out4['(Q(x)1+1(x)Q)K0=0'] = (Q * K0 + one_x(Q, K0)).applyfunc(sp.expand).is_zero_matrix
    out4['(QGF(x)1-1(x)QGF)K0=0 (QGF opposite adjointness to Q)'] = (G * K0 - one_x(G, K0)).applyfunc(sp.expand).is_zero_matrix
    out4['(QGF(x)1+1(x)QGF)K0=0'] = (G * K0 + one_x(G, K0)).applyfunc(sp.expand).is_zero_matrix
    # Im Q cap Im QGF != 0 : Q e_c = QGF v
    col = Q[:, IDX[('c', 0)]]
    v = sp.zeros(NC, 1)
    for mu in range(4):
        v[IDX[('As', mu)]] = col[IDX[('A', mu)]] / sol[0][syms[1]]
    out4['ImQ_cap_ImQGF_contains_Qc=-dc (axiom 4.3.3(5) fails)'] = ((G * v - col).applyfunc(sp.expand).is_zero_matrix and not col.is_zero_matrix)
    rec('4_free_complex', out4)

    # ---------------- 5. K_t, P
    Kt = sp.exp(-t_ * P2) * K0
    Lt = G * Kt
    Pm = (sp.exp(-eps_ * P2) / P2) * (G * K0)
    Pm = Pm.applyfunc(sp.expand)
    Pn = neg_p(Pm)
    gs = all(sp.simplify(Pn[b, a] - (-1) ** (PAR[a] * PAR[b]) * Pm[a, b]) == 0 for a in range(NC) for b in range(NC))
    out5['P(eps,inf) graded symmetric: P^ba(-p)=(-1)^{|a||b|}P^ab(p)'] = gs
    QP = (Q * Pm + one_x(Q, Pm)).applyfunc(sp.simplify)
    out5['QP(eps,inf)=K_eps (exact)'] = (QP - sp.exp(-eps_ * P2) * K0).applyfunc(sp.simplify).is_zero_matrix
    out5['(Q(x)1+1(x)Q)K_t=0'] = (Q * Kt + one_x(Q, Kt)).applyfunc(sp.simplify).is_zero_matrix
    out5['dK_t/dt+(Q(x)1+1(x)Q)L_t=0'] = (sp.diff(Kt, t_) + Q * Lt + one_x(Q, Lt)).applyfunc(sp.simplify).is_zero_matrix
    comps = {}
    for a in range(NC):
        for b in range(NC):
            if Pm[a, b] != 0:
                key = f'{COMP[a][0]}-{COMP[b][0]}'
                comps.setdefault(key, set()).add(str(sp.factor(Pm[a, b])))
    out5['P_nonzero_sectors'] = sorted(comps)
    PAA = Pm[IDX[('A', 0)]:IDX[('A', 0)] + 10, IDX[('A', 0)]:IDX[('A', 0)] + 10]
    out5['P_AA = 2/kappa * delta * exp(-eps p^2)/p^2'] = (PAA - 2 / kap * sp.exp(-eps_ * P2) / P2 * sp.eye(10)).applyfunc(sp.simplify).is_zero_matrix
    PcAs = [sp.simplify(Pm[IDX[('c', 0)], IDX[('As', mu)]] * P2 / sp.exp(-eps_ * P2)) for mu in range(10)]
    out5['P_c,A*_mu * p^2 e^{eps p^2}'] = [str(x) for x in PcAs]
    PAsc = [sp.simplify(Pm[IDX[('As', mu)], IDX[('c', 0)]] * P2 / sp.exp(-eps_ * P2)) for mu in range(10)]
    out5['P_A*_mu,c * p^2 e^{eps p^2}'] = [str(x) for x in PAsc]
    Pll = Pm[IDX[('lam', 0)]:IDX[('lam', 0)] + 16, IDX[('lam', 0)]:IDX[('lam', 0)] + 16]
    sbp = sum((nm.sigbar_sym(mu) * p[mu] for mu in range(4)), sp.zeros(16, 16))
    out5['P_lam,lam = (i/kappa) sigbar^mu p_mu exp(-eps p^2)/p^2'] = (Pll - sp.I / kap * sbp * sp.exp(-eps_ * P2) / P2).applyfunc(sp.simplify).is_zero_matrix
    rec('5_kernels_propagator', out5)

    # ---------------- 6. Hessian inversion / positivity / real slice
    P0 = (G * K0 / P2).applyfunc(sp.expand)               # P(0,inf)
    sl = lambda f: slice(IDX[(f, 0)], IDX[(f, 0)] + nm.FDIM[f])
    Nll = N[sl('lam'), sl('lam')]
    out6['P_lamlam(0,inf) = -(N_lamlam)^{-1}'] = (P0[sl('lam'), sl('lam')] + Nll.inv()).applyfunc(sp.simplify).is_zero_matrix
    NAA = N[sl('A'), sl('A')]
    pv = sp.Matrix(10, 1, lambda m, _: p[m] if m < 4 else 0)
    PiT = sp.eye(10) - pv * pv.T / P2
    out6['N_AA = -kappa/2 (p^2 - p p^T)'] = (NAA + kap / 2 * (P2 * sp.eye(10) - pv * pv.T)).applyfunc(sp.simplify).is_zero_matrix
    out6['P_AA(0,inf)*(-N_AA) = transverse projector'] = (P0[sl('A'), sl('A')] * (-NAA) - PiT).applyfunc(sp.simplify).is_zero_matrix
    NGF = (-(P0[sl('A'), sl('A')]).inv() - NAA).applyfunc(sp.simplify)
    out6['Feynman completion: -(N_AA+N_GF)^{-1}=P_AA with N_GF = -kappa/2 p p^T (Hessian of -1/4 (d^*A,d^*A))'] = (NGF + kap / 2 * pv * pv.T).applyfunc(sp.simplify).is_zero_matrix
    # c - A^* sector restricted to Im QGF: A^*_mu = d_mu s  (s odd scalar), c
    Vm = sp.zeros(NC, 2)
    Vm[IDX[('c', 0)], 0] = 1
    for mu in range(4):
        Vm[IDX[('As', mu)], 1] = sp.I * p[mu]
    # B(k) part: S2 restricted: N_L = V(-k)^T N(k) V(k)  (N already graded-symmetrised)
    NL = (neg_p(Vm).T * N * Vm).applyfunc(sp.expand)
    PL = -NL.inv()
    Ppush = (Vm * PL * neg_p(Vm).T).applyfunc(sp.simplify)
    ok_c = all(sp.simplify(Ppush[IDX[('c', 0)], IDX[('As', mu)]] - P0[IDX[('c', 0)], IDX[('As', mu)]]) == 0 and
               sp.simplify(Ppush[IDX[('As', mu)], IDX[('c', 0)]] - P0[IDX[('As', mu)], IDX[('c', 0)]]) == 0 for mu in range(4))
    out6['ghost (c,A^*) propagator = -(N restricted to c + exact A^*)^{-1} pushed forward'] = ok_c
    out6['N_L(c, d s) matrix'] = str(NL)
    # positivity
    ev = (P0[sl('A'), sl('A')] * P2).eigenvals()
    out6['eigenvalues of p^2 P_AA(0,inf)'] = {str(k): int(v) for k, v in ev.items()}
    out6['A-A covariance positive (kappa>0)'] = all(sp.simplify(k).is_positive for k in ev)
    evN = NAA.subs({p[0]: sp.Symbol('q', positive=True), p[1]: 0, p[2]: 0, p[3]: 0}).eigenvals()
    out6['eigenvalues of N_AA at p=(q,0,0,0) [<x,Qx> on real A]'] = {str(k): int(v) for k, v in evN.items()}
    out6['real slice: <A,QA> = 2 S2|_A negative semidefinite, kernel = pure gauge'] = all((k == 0) or sp.simplify(k).is_negative for k in evN)
    rec('6_hessian_positivity', out6)
    return Q, G, K0, Pm, N


# ======================================================================================
# 7. position space
# ======================================================================================
def section7():
    out = {}
    r, u = sp.symbols('r u', positive=True)
    h = (4 * sp.pi * t_) ** -2 * sp.exp(-r ** 2 / (4 * t_))
    Geps = sp.integrate(h, (t_, eps_, sp.oo))
    target = (1 - sp.exp(-r ** 2 / (4 * eps_))) / (4 * sp.pi ** 2 * r ** 2)
    out['int_eps^inf h_t dt = (1-exp(-r^2/4eps))/(4 pi^2 r^2)'] = sp.simplify(Geps - target) == 0
    lap = -(sp.diff(target, r, 2) + 3 / r * sp.diff(target, r))
    out['Delta_+ G_eps = h_eps (positive Laplacian, radial 4d)'] = sp.simplify(lap - h.subs(t_, eps_)) == 0
    dG = sp.simplify(sp.diff(target, r) / r)
    out['(1/r) dG_eps/dr'] = str(dG)
    out['G_eps(0) finite = 1/(16 pi^2 eps)'] = sp.simplify(sp.limit(target, r, 0) - 1 / (16 * sp.pi ** 2 * eps_)) == 0
    # heat kernel normalisation: int h_t d^4x = 1
    out['int h_t d^4x = 1'] = sp.simplify(sp.integrate(h * 2 * sp.pi ** 2 * r ** 3, (r, 0, sp.oo)) - 1) == 0
    rec('7_position_space', out)


# ======================================================================================
# 8. vertex functions
# ======================================================================================
def evaluate_on_test_fields(V, P, legs, nodd):
    """legs: list of dicts {(field,comp,color): value}, momentum, and Grassmann flag.
    Returns coefficient of theta_1...theta_m (odd legs in order) of the multilinear part of P."""
    # Grassmann parameters theta_j as odd variables of a dummy field; momenta enter via jets
    W = J.Vars({'th': 1})
    total = {}
    for (e, o), c in P.items():
        vs = list(e) + list(o)
        if len(vs) != len(legs):
            continue
        # sum over bijections variables -> legs
        for perm in itertools.permutations(range(len(legs))):
            val = sp.Integer(c[0]) + sp.I * c[1]
            ok = True
            thetas = []
            for v, L in zip(vs, perm):
                f, a, i, J4 = V.key[v]
                leg = legs[L]
                x = leg['pol'].get((f, a, i), 0)
                if x == 0:
                    ok = False
                    break
                val *= x * jet_to_symbol(J4, leg['p'])
                if leg['odd']:
                    thetas.append(leg['k'])
            if not ok:
                continue
            # theta ordering: product in the order of appearance in the monomial (evens carry no theta)
            sgn = 1
            arr = list(thetas)
            for x in range(len(arr)):
                for y in range(x + 1, len(arr)):
                    if arr[x] > arr[y]:
                        sgn = -sgn
            total_key = tuple(sorted(arr))
            total[total_key] = total.get(total_key, 0) + sgn * val
    return {k: sp.simplify(v / Ac.SCALE) for k, v in total.items()}


def section8():
    out = {}
    V = J.Vars(nm.PARITY)
    d, f = nm.f_int('su2')
    parts = Ac.build_action(V, f, d)
    I3, I4 = parts['I3'], parts['I4']
    rng = random.Random(11)
    R = lambda: sp.Rational(rng.randint(-5, 5), rng.randint(1, 3))
    F = lambda i, j, k: f.get((i, j, k), 0)

    def br(x, y):   # [x,y] for g-vectors (lists)
        return [sum(F(j, k, i) * x[j] * y[k] for j in range(d) for k in range(d)) for i in range(d)]

    def gp(x, y):   # (x,y) on g with kappa = 1
        return sum(x[i] * y[i] for i in range(d))

    def rand_g():
        return [R() for _ in range(d)]

    def rand_A():
        return [[R() for _ in range(d)] for _ in range(10)]   # a[M][i]

    def rand_p():
        return [R() for _ in range(4)]

    def pol_A(a):
        return {('A', M, i): a[M][i] for M in range(10) for i in range(d)}
    # ---- AAA
    a = [rand_A() for _ in range(3)]
    pp = [rand_p() for _ in range(2)]
    pp.append([-(pp[0][m] + pp[1][m]) for m in range(4)])
    legs = [{'pol': pol_A(a[j]), 'p': pp[j], 'odd': False, 'k': j} for j in range(3)]
    val = evaluate_on_test_fields(V, I3, legs, 0).get((), 0)

    def iota(pv, aa):  # sum_mu p_mu a_mu  in g
        return [sum(pv[m] * aa[m][i] for m in range(4)) for i in range(d)]

    def brV(x, y):     # sum_M [x_M, y_M]
        tot = [0] * d
        for M in range(10):
            b = br(x[M], y[M])
            tot = [tot[i] + b[i] for i in range(d)]
        return tot

    def sub(u, v):
        return [u[m] - v[m] for m in range(4)]
    formula = -sp.I / 2 * (gp(iota(sub(pp[1], pp[0]), a[2]), brV(a[0], a[1])) +
                           gp(iota(sub(pp[2], pp[1]), a[0]), brV(a[1], a[2])) +
                           gp(iota(sub(pp[0], pp[2]), a[1]), brV(a[2], a[0])))
    out['AAA: V=-(i/2) sum_cyc (iota_{p2-p1} a3, [a1,a2])'] = sp.simplify(val - formula) == 0
    # ---- AAAA
    a = [rand_A() for _ in range(4)]
    legs = [{'pol': pol_A(a[j]), 'p': [0, 0, 0, 0], 'odd': False, 'k': j} for j in range(4)]
    val = evaluate_on_test_fields(V, I4, legs, 0).get((), 0)

    def wedge(x, y):   # [x ^ y]_{MN} = [x_M,y_N] - [x_N,y_M]
        return {(M, N): [br(x[M], y[N])[i] - br(x[N], y[M])[i] for i in range(d)] for M in range(10) for N in range(10) if M < N}

    def hodge(X, Y):
        return sum(gp(X[k], Y[k]) for k in X)
    formula = -sp.Rational(1, 2) * (hodge(wedge(a[0], a[1]), wedge(a[2], a[3])) + hodge(wedge(a[0], a[2]), wedge(a[1], a[3])) +
                                    hodge(wedge(a[0], a[3]), wedge(a[1], a[2])))
    out['AAAA: V=-(1/2) sum_pairings ([a1^a2],[a3^a4])_Hodge'] = sp.simplify(val - formula) == 0
    # ---- lam lam A : theta1 theta2 (u1, rho([a,u2]))
    u = [[[R() for _ in range(d)] for _ in range(16)] for _ in range(2)]   # u[j][alpha][i]
    aa = rand_A()
    legs = [{'pol': {('lam', al, i): u[0][al][i] for al in range(16) for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 0},
            {'pol': {('lam', al, i): u[1][al][i] for al in range(16) for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 1},
            {'pol': pol_A(aa), 'p': [0] * 4, 'odd': False, 'k': 2}]
    val = evaluate_on_test_fields(V, I3, legs, 2).get((0, 1), 0)
    sig = [nm.sig_sym(M) for M in range(10)]
    formula = 0
    for M in range(10):
        for al in range(16):
            for be in range(16):
                if sig[M][al, be] != 0:
                    formula += gp(u[0][al], br(aa[M], u[1][be])) * sig[M][al, be]
    out['lam lam A: V=theta1 theta2 (u1, sum_M rho(e^M)[a_M,u2])'] = sp.simplify(val - formula) == 0
    # ---- A c A^* : theta_c theta_* ([a,gamma], alpha)
    gam = rand_g()
    alp = rand_A()
    legs = [{'pol': pol_A(aa), 'p': [0] * 4, 'odd': False, 'k': 0},
            {'pol': {('c', 0, i): gam[i] for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 1},
            {'pol': {('As', M, i): alp[M][i] for M in range(10) for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 2}]
    val = evaluate_on_test_fields(V, I3, legs, 2).get((1, 2), 0)
    formula = sum(gp(br(aa[M], gam), alp[M]) for M in range(10))
    out['A c A*: V=theta_c theta_* ([a,gamma],alpha)'] = sp.simplify(val - formula) == 0
    # ---- c lam lam^*
    uu = [[R() for _ in range(d)] for _ in range(16)]
    vv = [[R() for _ in range(d)] for _ in range(16)]
    legs = [{'pol': {('c', 0, i): gam[i] for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 0},
            {'pol': {('lam', al, i): uu[al][i] for al in range(16) for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 1},
            {'pol': {('lams', al, i): vv[al][i] for al in range(16) for i in range(d)}, 'p': [0] * 4, 'odd': False, 'k': 2}]
    val = evaluate_on_test_fields(V, I3, legs, 2).get((0, 1), 0)
    formula = sum(gp(br(gam, uu[al]), vv[al]) for al in range(16))
    out['c lam lam*: V=theta_c theta_lam ([gamma,u],v)'] = sp.simplify(val - formula) == 0
    # ---- c c c^*
    g1, g2, w = rand_g(), rand_g(), rand_g()
    legs = [{'pol': {('c', 0, i): g1[i] for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 0},
            {'pol': {('c', 0, i): g2[i] for i in range(d)}, 'p': [0] * 4, 'odd': True, 'k': 1},
            {'pol': {('cs', 0, i): w[i] for i in range(d)}, 'p': [0] * 4, 'odd': False, 'k': 2}]
    val = evaluate_on_test_fields(V, I3, legs, 2).get((0, 1), 0)
    out['c c c*: V=theta1 theta2 ([gamma1,gamma2],w)'] = sp.simplify(val - gp(br(g1, g2), w)) == 0
    rec('8_vertex_functions', out)


# ======================================================================================
# 9. 0-dim check of the expansion convention of W(P,I) = hbar log(e^{hbar d_P} e^{I/hbar})
# ======================================================================================
class G0:
    """tiny Grassmann algebra over Q with an hbar power: key (ev, od, h) -> Fraction"""
    def __init__(self, par):
        self.par = par

    def mul(self, A, B):
        R = {}
        for (e1, o1, h1), c1 in A.items():
            for (e2, o2, h2), c2 in B.items():
                s, o = J._merge_odd(o1, o2)
                if s == 0:
                    continue
                k = (tuple(sorted(e1 + e2)), o, h1 + h2)
                R[k] = R.get(k, 0) + s * c1 * c2
        return {k: v for k, v in R.items() if v != 0}

    def add(self, A, B, s=1):
        R = dict(A)
        for k, v in B.items():
            R[k] = R.get(k, 0) + s * v
        return {k: v for k, v in R.items() if v != 0}

    def dl(self, A, v):
        R = {}
        for (e, o, h), c in A.items():
            if self.par[v]:
                if v in o:
                    k = o.index(v)
                    R[(e, o[:k] + o[k + 1:], h)] = R.get((e, o[:k] + o[k + 1:], h), 0) + (c if k % 2 == 0 else -c)
            else:
                m = e.count(v)
                if m:
                    k = e.index(v)
                    key = (e[:k] + e[k + 1:], o, h)
                    R[key] = R.get(key, 0) + m * c
        return {k: v for k, v in R.items() if v != 0}


def _w0dim(par, fl, bgs, P, seed, NORD):
    """exact one-loop part of log(e^{hbar d_P} e^{I/hbar})|_{fluct=0} for random I (quadratic in
    fluctuations x linear in backgrounds), d_P = 1/2 sum P^{ab} d_{e_b} d_{e_a} [C07 Sec 6.4]."""
    G = G0(par)
    rng = random.Random(seed)
    Rr = lambda: Fraction(rng.randint(-3, 3), rng.randint(1, 2))

    def monomial(vs, c):
        A = {((), (), 0): Fraction(1)}
        for x in vs:
            A = G.mul(A, {((x,), (), 0) if not par[x] else ((), (x,), 0): Fraction(1)})
        return {k: c * w for k, w in A.items()}
    I = {}
    for a in fl:
        for b in fl:
            if a > b:
                continue
            for bg in bgs:
                if (par[a] + par[b] + par[bg]) % 2 == 0 and not (a == b and par[a]):
                    I = G.add(I, monomial([a, b, bg], Rr()))

    def dP(F):
        R = {}
        for (a, b), c in P.items():
            if c == 0:
                continue
            T = G.dl(G.dl(F, a), b)
            R = G.add(R, {k: Fraction(1, 2) * c * w for k, w in T.items()})
        return R
    E = {((), (), 0): Fraction(1)}
    term = dict(E)
    for n in range(1, NORD + 1):
        term = G.mul(term, {(k[0], k[1], k[2] - 1): w / n for k, w in I.items()})
        E = G.add(E, term)
    Z = dict(E)
    T = dict(E)
    for k in range(1, 2 * NORD + 1):
        T = dP(T)
        T = {(kk[0], kk[1], kk[2] + 1): w / k for kk, w in T.items()}
        Z = G.add(Z, T)
    Z0 = {k: w for k, w in Z.items() if not any(x in fl for x in k[0] + k[1])}
    Y = G.add(Z0, {((), (), 0): Fraction(1)}, -1)
    logZ = {}
    Yp = {((), (), 0): Fraction(1)}
    for n in range(1, NORD + 1):
        Yp = G.mul(Yp, Y)
        logZ = G.add(logZ, {k: w * Fraction((-1) ** (n + 1), n) for k, w in Yp.items()})
    logZ = {k: w for k, w in logZ.items() if len(k[0]) + len(k[1]) <= NORD}
    hp = sorted(set(k[2] for k in logZ))
    one = {(k[0], k[1], 0): w for k, w in logZ.items() if k[2] == 0}

    def formula(sign_mode, super_):
        # X^a_c = sum_b P^{ab} s(b) (d_{e_b} d_{e_c} I),  d_{e_c} applied first
        X = {}
        for a in fl:
            for c in fl:
                acc = {}
                for b in fl:
                    pab = P.get((a, b), 0)
                    if pab:
                        D = G.dl(G.dl(I, c), b)
                        sg = (-1) ** par[b] if sign_mode == 'D' else 1
                        acc = G.add(acc, {k: pab * sg * w for k, w in D.items()})
                X[(a, c)] = acc
        Xn = {(a, c): ({((), (), 0): Fraction(1)} if a == c else {}) for a in fl for c in fl}
        form = {}
        for n in range(1, NORD + 1):
            new = {}
            for a in fl:
                for c in fl:
                    acc = {}
                    for b in fl:
                        acc = G.add(acc, G.mul(Xn[(a, b)], X[(b, c)]))
                    new[(a, c)] = acc
            Xn = new
            for a in fl:
                sg = (-1) ** par[a] if super_ else 1
                form = G.add(form, {k: w * Fraction(sg, 2 * n) for k, w in Xn[(a, a)].items()})
        form = {k: w for k, w in form.items() if len(k[0]) + len(k[1]) <= NORD}
        return len(G.add(one, form, -1)) == 0
    wick = True
    for a in fl:
        for b in fl:
            if a == b and par[a]:
                continue
            F = monomial([a, b], Fraction(1))
            if dP(F).get(((), (), 0), 0) != P.get((a, b), 0):
                wick = False
    return hp, formula('D', True), formula('1', True), formula('1', False), wick


def section9():
    out = {}
    F_ = Fraction
    configs = [
        ([0, 0, 1, 1, 0, 0, 1, 1], [0, 1, 2, 3], [4, 5, 6, 7],
         {(0, 0): F_(1), (1, 1): F_(-2), (0, 1): F_(1, 2), (1, 0): F_(1, 2), (2, 3): F_(3, 2), (3, 2): F_(-3, 2)}),
        ([0, 1, 1, 0, 1], [0, 1, 2], [3, 4], {(0, 0): F_(1), (1, 2): F_(3), (2, 1): F_(-3)}),
        ([1, 1, 0, 0, 1, 1], [0, 1, 2, 3], [4, 5], {(2, 2): F_(1), (3, 3): F_(2), (0, 1): F_(5), (1, 0): F_(-5)}),
    ]
    main, c1, c2, wk, hps = [], [], [], [], set()
    for (par, fl, bgs, P) in configs:
        for seed in (1, 2, 3):
            hp, a, b, c, w = _w0dim(par, fl, bgs, P, seed, 3)
            hps.update(hp)
            main.append(a); c1.append(b); c2.append(c); wk.append(w)
    out['test spaces: (2|2),(1|2),(2|2) fluctuations with even+odd backgrounds; 3 random I each; order 3'] = True
    out['hbar powers of log Z (I quadratic in fluctuations => one loop only)'] = sorted(hps)
    out['W_1 = hbar sum_n (1/2n) str[(D P I2)^n], (I2)_bc = d_{e_b} d_{e_c} I, D = diag((-1)^{|b|}), str = sum_a (-1)^{|a|}'] = all(main)
    out['(control) str[(P I2)^n] without D'] = all(c1)
    out['(control) tr[(P I2)^n]'] = all(c2)
    out['Wick: d_P(phi^a phi^b) = P^{ab}'] = all(wk)
    rec('9_W_expansion_0dim', out)


if __name__ == '__main__':
    t0 = time.time()
    which = sys.argv[1:] or ['1', '2', '3', '46', '7', '8', '9']
    if '1' in which: section1()
    if '2' in which: section2()
    if '3' in which: section3()
    if '46' in which: section4_6()
    if '7' in which: section7()
    if '8' in which: section8()
    if '9' in which: section9()
    RESULTS['total_seconds'] = round(time.time() - t0, 1)
    fn = os.path.join(HERE, 'checks.json' if not sys.argv[1:] else 'checks_partial.json')
    with open(fn, 'w') as fh:
        json.dump(RESULTS, fh, indent=1, default=str)
    print('wrote', fn)
