#!/usr/bin/env python3
"""Adversarial audit (lens: targets and plan) of step2/theory/inner_product_design.tex and engine_plan_ip.json.

Single process, exact sympy arithmetic (numpy only for the static flavour-sum cross-check), < 0.2 GB, no Wolfram,
no network.  Reads only files of this project folder.  Writes ip_audit2_checks.json next to this file.

Code-level names (json/py only): u = unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)); frame Z = the frame P_1 = 0 in which
Lemma ip-W / Thm ip-origin identify the radial Gram with the two-point coefficient at |y| = ell; frame E = the engines'
frame (translation-invariant arity-one finite counterterm, P_1 = r_P u {S1(P),-}, r_P = -(1/3 + f) on scalar bilinears,
step2/qa/compare_delta.json).

A1  frame bookkeeping of the radial lines T1, T2 (Prop ip-twenty, ip-konishi) and its interplay with the outer
    direction G_B (Cor c8-gram, Thm ip-laws (d)).
A2  the inner direction (Delta) and the outer direction B on the blocks of [SW (17)]: which data separate them.
A3  the sl(2) model of inner_product_checks.py C5: effect of the instruction of IP6 'recompute G1 (multiplied by
    1 + 2 c unit w)' applied to the two-point (frame-Z) Gram.
A4  order/half-edge census of the length-changing Gram entry (B_w, F_n) (Prop ip-orders, tp prop:t-census-w3, IP3).
A5  T5 (D_1 = D_1^adj) on the weight-3 family versus T3 (iii) (orthogonality).
A6  [Ryz] fixed-ratio target versus the IP1 targets.
A7  [BKRS] (24) versus (27): direction of the coupling-normalisation factor (ip_targets.tex Target 11, C7 note).
A8  exact recomputation of the T3 (iv) inputs: flavour sum 48|w|^4 of B_w and k = -24.
A9  static scans: arXiv version stamps of the local copies versus the declared versions; absolute paths in engine
    code (cluster portability); plan items that touch S_- (x) W^* (reflection compatibility of G_B).
"""
import json, os, re, glob, itertools, random
import sympy as sp
from sympy.combinatorics import Permutation

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
OUT = {}

def zero(M):
    return all(sp.simplify(e) == 0 for e in (M if isinstance(M, sp.MatrixBase) else [M]))

# ------------------------------------------------------------------------------------------------------------
# A1 frame bookkeeping of the radial lines
# ------------------------------------------------------------------------------------------------------------
def A1():
    u = sp.symbols('u', real=True)                      # the unit
    nK, x2 = sp.symbols('n_K xi2', positive=True)       # tree norm of K, (xi, conj xi)
    tau = sp.symbols('tau', real=True, nonzero=True)
    mK, mz = sp.symbols('m_K m_zeta', real=True)        # two-point (frame Z) one-loop constants
    qr, qi = sp.symbols('q_r q_i', real=True)
    cK, cz, a = sp.symbols('c_K c_zeta a', real=True)    # per-family translation rescalings r_P, real G_B coefficient
    sig = -2/tau                                         # T0: tau sigma = -2
    nz = -sig*nK/tau                                     # tree adjointness
    q = qr + sp.I*qi
    s = sp.conjugate(q) + (mz - mK)                      # frame-Z data obeying the radial line exactly
    def data(qq, ss):
        X = sp.Matrix([[0, 0], [-(1 + u*qq)*tau, 0]])            # -Lambda^(1)(xi) on (K, O_xi)
        Y = sp.Matrix([[0, (1 + u*ss)*sig*x2], [0, 0]])          # -Lambda^(1)(xi^dagger)
        return X, Y
    G_Z = sp.diag(nK*(1 + u*mK), nz*x2*(1 + u*mz))              # radial Gram = two-point coefficient (frame Z)
    XZ, YZ = data(q, s)
    rel = lambda X, Y, G: (X.H*G - G*Y).applyfunc(lambda e: sp.expand(sp.series(e, u, 0, 2).removeO()))
    ok_Z = zero(rel(XZ, YZ, G_Z))
    # frame E from frame Z: X_Z = U^-1 X_E U, Gram_Z = U^H Gram_E U, U = 1 + u diag(w_K c_K, w_zeta c_zeta)
    U = sp.diag(1 + 2*u*cK, 1 + sp.Rational(5, 2)*u*cz)
    Ui = sp.diag(1 - 2*u*cK, 1 - sp.Rational(5, 2)*u*cz)       # inverse through order u
    XE = (U*XZ*Ui).applyfunc(sp.expand); YE = (U*YZ*Ui).applyfunc(sp.expand)
    G_E = (Ui.H*G_Z*Ui).applyfunc(sp.expand)
    ok_E_consistent = zero(rel(XE, YE, G_E))                    # covariance: frame E with its own Gram
    # read off q_E, s_E and the frame-E Gram constants
    qE = sp.simplify(sp.expand(-XE[1, 0]/tau - 1).coeff(u, 1))
    sE = sp.simplify(sp.expand(YE[0, 1]/(sig*x2) - 1).coeff(u, 1))
    mKE = sp.simplify(sp.expand(G_E[0, 0]/nK - 1).series(u, 0, 2).removeO().coeff(u, 1))
    mzE = sp.simplify(sp.expand(G_E[1, 1]/(nz*x2) - 1).series(u, 0, 2).removeO().coeff(u, 1))
    # (i) radial line evaluated with frame-E generators and the two-point (frame-Z) constants, as IP5 is written
    resid_mixed = sp.simplify(sE - sp.conjugate(qE) - (mz - mK))
    # (ii) same with a real outer-direction coefficient a: q -> q - a, s -> s + a (Prop ip-konishi)
    resid_mixed_GB = sp.simplify((sE + a) - sp.conjugate(qE - a) - (mz - mK))
    a_abs = sp.solve(sp.Eq(resid_mixed_GB, 0), a)[0]
    # (iii) the inner direction (universal c) and G_B act identically on (q, s) of the pair (xi, xi^dagger)
    c = sp.symbols('c', real=True)
    inner_q = sp.simplify(qE.subs({cK: c, cz: c}) - q)          # change of q from frame Z to frame E (universal c)
    inner_s = sp.simplify(sE.subs({cK: c, cz: c}) - s)
    # (iv) tree-level dilatation Ward identity: the inner field redefinition phi ~ S1(Delta) does not change the
    #      two-point function at separated points (x.d_x + y.d_y + 2w) |x-y|^{-2w} = 0
    xs = sp.symbols('x1:5', real=True); ys = sp.symbols('y1:5', real=True); w = sp.symbols('w', positive=True)
    r2 = sum((xi - yi)**2 for xi, yi in zip(xs, ys)); g = r2**(-w)
    ward = sp.simplify(sum(xi*sp.diff(g, xi) for xi in xs) + sum(yi*sp.diff(g, yi) for yi in ys) + 2*w*g)
    OUT['A1_frame_bookkeeping'] = {
        'frame_Z_relation_exact': ok_Z,
        'frame_E_relation_with_frame_E_Gram': ok_E_consistent,
        'q_E_minus_q': str(sp.simplify(qE - q)), 's_E_minus_s': str(sp.simplify(sE - s)),
        'm_K_E_minus_m_K': str(sp.simplify(mKE - mK)), 'm_zeta_E_minus_m_zeta': str(sp.simplify(mzE - mz)),
        'T2_radial_residual_mixed_frames (IP5 as written)': str(resid_mixed),
        'T2_residual_universal_c': str(sp.simplify(resid_mixed.subs({cK: c, cz: c}))),
        'T2_residual_universal_c_engines_value (c = -1/3, f = 0)': str(resid_mixed.subs({cK: sp.Rational(-1, 3), cz: sp.Rational(-1, 3)})),
        'real_GB_coefficient_that_absorbs_it': str(sp.simplify(a_abs)),
        'inner_direction_on_q_s (universal c)': [str(inner_q), str(inner_s)],
        'GB_on_q_s': ['-a', '+a'],
        'same_pattern_inner_vs_GB_on_pair (Delta q = -Delta s for both)': sp.simplify(inner_q + inner_s) == 0 and sp.simplify(inner_q) != 0,
        'two_point_unchanged_by_inner_redefinition (tree dilatation Ward identity)': ward == 0,
        'T1_note': 'identical algebra with (K, O_xi) -> (20prime, 20): residual 4 c_20prime - 5 c_20; for a universal c both T1 and T2 '
                   'are off by -c and both are restored by the same real G_B coefficient a = c/2',
    }
    return ok_Z and ok_E_consistent and sp.simplify(resid_mixed - (4*cK - 5*cz)) == 0 and ward == 0

# ------------------------------------------------------------------------------------------------------------
# A2 inner (Delta) versus outer (B) direction on the blocks of [SW (17)]
# ------------------------------------------------------------------------------------------------------------
def A2():
    s1, s2, s3 = sp.symbols('s1 s2 s3', real=True)
    r = lambda a, b: [s1, s2, s3][a - 1] - [s1, s2, s3][b - 1]     # Thm ip-P (ii): r = s_a - s_b (C9)
    blocks = {'P (2,1)': (2, 1), 'K (1,2)': (1, 2), 'Q xi in S+W (3,1)': (3, 1), 'S eta in S+W* (1,3)': (1, 3),
              'Qbar in S-W* (2,3)': (2, 3), 'Sbar in S-W (3,2)': (3, 2)}
    t, b = sp.symbols('t b', real=True)
    inner = {s1: t, s2: -t, s3: 0}          # Delta = diag(1,-1|0) up to scale
    outer = {s1: 0, s2: 0, s3: b}           # B = diag(0,0|1)
    tab = {k: {'inner': str(sp.simplify(r(*v).subs(inner))), 'outer': str(sp.simplify(r(*v).subs(outer)))} for k, v in blocks.items()}
    # on the pair used by T1, T2 (Q on S+W, S on S+W*) both directions give the pattern (+x, -x)
    same_pair = (sp.simplify(r(3, 1).subs(inner) + r(1, 3).subs(inner)) == 0 and
                 sp.simplify(r(3, 1).subs(outer) + r(1, 3).subs(outer)) == 0)
    # reflection compatibility r(2,3) = conj r(3,1) (Theta_n maps S+W to S-W*) with real s: forces s3 = (s1+s2)/2,
    # i.e. r(3,1) = r(2,1)/2: the real rescaling of xi is fixed by r_P only with data on S-W*
    refl = sp.solve(sp.Eq(r(2, 3), r(3, 1)), s3)
    rxi_given_refl = sp.simplify(r(3, 1).subs(s3, refl[0]))
    OUT['A2_inner_vs_outer'] = {'block_rescalings': tab, 'pair_Q_S_same_pattern_for_both': same_pair,
                                'reflection_condition_solution_s3': str(refl[0]),
                                'r_xi_under_reflection_and_MC': str(rxi_given_refl),
                                'r_P': str(r(2, 1)),
                                'statement': 'inner and outer act on (Q on S+W, S on S+W*) by the same pattern; they are '
                                             'separated only by S-W* / S-W (reflection) or by P, K together with the Gram'}
    return same_pair and sp.simplify(rxi_given_refl - r(2, 1)/2) == 0

# ------------------------------------------------------------------------------------------------------------
# A3 sl(2) model: IP6 instruction applied to the frame-Z (two-point) Gram
# ------------------------------------------------------------------------------------------------------------
def A3():
    w, kap, c, u = sp.symbols('w kappa c u', real=True)
    N = 5
    h = w + u*kap
    G = sp.diag(*[sp.factorial(n)*sp.rf(2*h, n) for n in range(N)])
    P = sp.zeros(N); K = sp.zeros(N)
    for n in range(N - 1): P[n + 1, n] = 1
    for n in range(1, N): K[n - 1, n] = n*(2*h + n - 1)
    ser = lambda M, k: M.applyfunc(lambda e: sp.expand(sp.series(e, u, 0, 2).removeO()).coeff(u, k))
    G0, G1 = ser(G, 0), ser(G, 1); K0, K1 = ser(K, 0), ser(K, 1)
    M = G0.inv()*G1
    ok_Zline = zero((K1 - (K0*M - M*K0)).applyfunc(sp.simplify))
    Nop = sp.diag(*range(N))
    G1_plan = G1 + 2*c*Nop*G0                    # 'G1 multiplied by 1 + 2 c unit (level)' applied to the frame-Z Gram
    Mp = G0.inv()*G1_plan
    resid = (K1 - (K0*Mp - Mp*K0)).applyfunc(sp.simplify)
    ok_resid = zero((resid - 2*c*K0).applyfunc(sp.simplify)) or zero((resid + 2*c*K0).applyfunc(sp.simplify))
    OUT['A3_sl2_IP6_double_frame'] = {'K_line_frame_Z_with_two_point_Gram': ok_Zline,
                                      'K_line_residual_after_IP6_rescaling_of_two_point_G1': str(sp.simplify(resid[0, 1])) + ' (entry (0,1)); equals '
                                      + ('2c K_0' if zero((resid - 2*c*K0).applyfunc(sp.simplify)) else '-2c K_0'),
                                      'nonzero_for_c_ne_0': ok_resid}
    return ok_Zline and ok_resid

# ------------------------------------------------------------------------------------------------------------
# A4 census of the (B_w, F_n) entry
# ------------------------------------------------------------------------------------------------------------
def A4():
    # cubic vertices of the action of ch01 (field content only; d = derivative): from F^2, the Dirac term, Yukawa, ghosts
    cubic = [('AV', 'AV', 'AV'), ('AV', 'AW', 'AW'), ('lam', 'lamb', 'AV'), ('lam', 'lam', 'AW'), ('lamb', 'lamb', 'AW'),
             ('AV', 'c', 'AsV'), ('AW', 'c', 'AsW'), ('c', 'lam', 'lams'), ('c', 'lamb', 'lambs'), ('c', 'c', 'cs')]
    quartic = [('AV',) * 4, ('AV', 'AV', 'AW', 'AW'), ('AW',) * 4]
    partner = {'AV': 'AV', 'AW': 'AW', 'lam': 'lamb', 'lamb': 'lam', 'c': 'AsV', 'AsV': 'c'}
    legs_B = ['AW', 'AW', 'AW']                     # B_w (and its conjugate)
    out = {}
    for name, legs_F in [('(B_w, F_n)', ['lam', 'lam']), ('(Bbar, Fbar_w)', ['lamb', 'lamb'])]:
        rows = []
        for combo in [('c3',), ('c4',), ('c3', 'c3'), ('c3', 'c4'), ('c3', 'c3', 'c3')]:
            val = sum(3 if x == 'c3' else 4 for x in combo)
            half = len(legs_B) + len(legs_F) + val
            if half % 2:
                rows.append({'internal': combo, 'half_edges': half, 'exists': False})
                continue
            E = half // 2; Vp = len(combo); b1 = E - (2 + Vp) + 1
            rows.append({'internal': combo, 'half_edges': half, 'edges': E, 'b1': b1, 'relative_order_unit_power': E - Vp - 2, 'exists': True})
        # one cubic vertex: by Target 18 (i) it carries two legs of B and one of F; vertex content = partners
        need = sorted([partner['AW'], partner['AW'], partner[legs_F[0]]])
        found = [v for v in cubic if sorted(v) == need]
        out[name] = {'classes': rows, 'one_cubic_vertex_needs': need, 'matching_cubic_vertices': found,
                     'one_cubic_class_even_and_at_relative_order_unit': any(r.get('internal') == ('c3',) and r['exists'] and r['relative_order_unit_power'] == 1 and r['half_edges'] % 2 == 0 for r in rows),
                     'vanishes_by_field_content': len(found) == 0}
    OUT['A4_census_B_F_entry'] = out
    return all(v['one_cubic_class_even_and_at_relative_order_unit'] and v['vanishes_by_field_content'] for v in out.values())

# ------------------------------------------------------------------------------------------------------------
# A5 T5 on the weight-3 family versus T3 (iii)
# ------------------------------------------------------------------------------------------------------------
def A5():
    # explicit complex rational data (sympy's series of conjugates of symbolic complex expressions proved unreliable);
    # eigenvalues -48u (Konishi descendant k) and 0 (protected E_w = e); D_1 = D_1^adj <=> G0 D_1 Hermitian
    u = sp.symbols('u', positive=True)
    I_ = sp.I
    trials = [(sp.Rational(1, 4) + I_/7, 2 - I_, 3 + 2*I_, 5, 7), (-sp.Rational(1, 4), 1 + 3*I_, -2 + I_, 3, 11),
              (I_/4, sp.Rational(1, 2) - I_, 4 - 3*I_, 2, 13)]
    res = []
    for kB, bF, bB, nF, nBp in trials:
        def A01(cK):
            G0 = sp.diag(nF, u*nBp)
            k = sp.Matrix([u*cK, kB])
            ell = sp.Matrix([[bB, -bF]])/(u*cK*bB - kB*bF)
            D1 = -48*u*k*ell
            Mx = G0*D1
            A = Mx - Mx.applyfunc(sp.conjugate).T
            return sp.simplify(A[0, 1]), sp.simplify(sp.series(D1[0, 0], u, 0, 2).removeO())
        cK_T3 = -kB*sp.conjugate(bB)*nBp/(sp.conjugate(bF)*nF)       # T3 (iii) with n_B = u n_Bp
        a_T3, dff = A01(cK_T3)
        a_off, _ = A01(cK_T3 + 1)
        res.append({'A01_with_T3iii_cK': str(a_T3), 'A01_with_cK_shifted_by_1 (leading)': str(sp.series(a_off, u, 0, 3).removeO()),
                    'D1_FF_rel_order_unit': str(dff)})
    ok = all(r['A01_with_T3iii_cK'] == '0' and r['D1_FF_rel_order_unit'] == '0' and r['A01_with_cK_shifted_by_1 (leading)'] != '0' for r in res)
    OUT['A5_T5_vs_T3iii'] = {'trials': res, 'equivalent_on_weight3_family': ok,
                             'weights_2_and_5/2': 'D_1 is a real multiple of the identity on each irreducible test family (Example ip-scalars, '
                                                  'Schur): D_1 = D_1^adj holds identically there'}
    return ok

# ------------------------------------------------------------------------------------------------------------
# A6 [Ryz] fixed ratio versus the IP1 targets; A7 [BKRS] (24) vs (27)
# ------------------------------------------------------------------------------------------------------------
def A6_A7():
    # IP1 targets (engine_plan_ip.json): m_20p = 27/2 + 4 log(3/2), log coefficient 0 on the 20'
    m20p = sp.Rational(27, 2) + 4*sp.log(sp.Rational(3, 2))
    ryz = {'Ryz_(33)_structure': 'every O(g^2) entry = Btilde N G_ij, Btilde = -(Y^2/32 pi^2)(log(x^2/eps^2) - 1): an entry with '
                                 'zero log coefficient (G_ij = 0, protected) has zero constant',
           'IP1_target_m_20p': str(m20p), 'numeric': float(m20p), 'IP1_log_coefficient_20p': 0,
           'fixed_ratio_holds_in_project_minimal_scheme': False,
           'after_field_rescaling_4a = -m_20p': 'trace part 24 (log|x|^2 + gamma_E - 3 log 2): ratio gamma_E - 3 log 2 = %.6f (regulator dependent; [Ryz] has -1 - log eps^2)' % float(sp.EulerGamma - 3*sp.log(2))}
    g, N = sp.symbols('g N', positive=True)
    coeff24 = -(3*N/(4*sp.pi**2))*3*(N**2 - 1)/(4*sp.pi**2)**2      # g_std^2 coefficient of ln x^2/x^4 from (24), (21)
    coeff27 = -9*N*(N**2 - 1)/(4*(4*sp.pi**2)**3)                   # printed (27)
    ratio = sp.simplify(coeff24/coeff27)
    # coeff27 g27^2 = coeff24 g_std^2  =>  g27^2 = ratio g_std^2
    OUT['A6_Ryz'] = ryz
    OUT['A7_BKRS_24_vs_27'] = {'coefficient_from_(24)_over_printed_(27)': str(ratio),
                               'implied': 'g^2 in (27) = %s x g_std^2 (BKRS footnote 3: old g of refs [6,7] = 2 x standard g)' % ratio,
                               'ip_targets_annotation': 'consistent with g^2 in (27) = g_std^2/4  (wrong direction)'}
    return ratio == 4

# ------------------------------------------------------------------------------------------------------------
# A8 T3 (iv) inputs recomputed exactly
# ------------------------------------------------------------------------------------------------------------
def A8():
    I_ = sp.I
    eps3 = lambda a, b, c: sp.LeviCivita(a, b, c)
    Sig = [sp.zeros(4) for _ in range(6)]
    for a in range(3):
        for b in range(3):
            for c in range(3):
                Sig[a][b, c] = eps3(a, b, c); Sig[a + 3][b, c] = I_*eps3(a, b, c)
            Sig[a][b, 3] = 1 if a == b else 0; Sig[a][3, b] = -1 if a == b else 0
            Sig[a + 3][b, 3] = -I_ if a == b else 0; Sig[a + 3][3, b] = I_ if a == b else 0
    # check sum_I Sig^I_ij Sig^I_kl = 2 eps_ijkl and A^{ij} = -A_I Sigbar^{Iij}, Sigbar = conj transpose
    ok_blk = True
    for i, j, k, l in itertools.product(range(4), repeat=4):
        if sum(Sig[I][i, j]*Sig[I][k, l] for I in range(6)) != 2*sp.LeviCivita(i, j, k, l):
            ok_blk = False; break
    Sigb = [S.H for S in Sig]
    for I in range(6):
        for i, j in itertools.product(range(4), repeat=2):
            if sp.simplify(sp.Rational(1, 2)*sum(sp.LeviCivita(i, j, k, l)*Sig[I][k, l] for k in range(4) for l in range(4)) + Sigb[I][i, j]) != 0:
                ok_blk = False
    rng = random.Random(28)
    vals = []
    for trial in range(3):
        wv = [rng.randint(-3, 3) + I_*rng.randint(-3, 3) for _ in range(4)]
        if all(x == 0 for x in wv): wv[0] = 1
        wcol = sp.Matrix(wv)
        vJ = [Sig[J]*wcol for J in range(6)]                   # (Sigma^J w)_k
        c = {}
        for I, J, K in itertools.product(range(6), repeat=3):
            # c_IJK = -Sigbar^{I ik} (Sigma^J w)_k (Sigma^K w)_i   (B_w = sum c_IJK ([A_I,A_J],A_K))
            c[(I, J, K)] = sp.expand(-(vJ[K].T*Sigb[I]*vJ[J])[0, 0])
        sgn = {p: Permutation(list(p)).signature() for p in itertools.permutations(range(3))}
        tot = 0
        for trip in itertools.permutations(range(6), 3):
            ca = sum(sgn[p]*c[tuple(trip[x] for x in p)] for p in sgn) / 6
            tot += sp.expand(ca*sp.conjugate(ca))
        w2 = sum(sp.expand(x*sp.conjugate(x)) for x in wv)
        vals.append(sp.nsimplify(sp.simplify(tot/w2**2)))
    flav = vals[0]
    ok_flav = all(v == 48 for v in vals)
    # k from the design's formula with the recomputed flavour sum
    hb, dimg = sp.symbols('hbar dim_g', positive=True); lam = sp.symbols('lambda', negative=True)
    w2s, u2, up2 = sp.symbols('w2 u2 up2', positive=True)
    nB = hb**2*(2/(4*sp.pi**2))**3*(-lam*dimg)*6*flav*w2s**2
    nz = 24*hb*dimg*up2*w2s/(16*sp.pi**4)
    unit = hb*lam/(16*sp.pi**2)
    k = sp.simplify(sp.Rational(1, 16)*nB/(unit*u2*w2s*nz)).subs(up2, 1/u2)
    OUT['A8_T3iv_inputs'] = {'thooft_blocks_ok': ok_blk, 'flavour_sum_over_w4_exact': [str(v) for v in vals], 'k': str(sp.simplify(k))}
    return ok_blk and ok_flav and sp.simplify(k + 24) == 0

# ------------------------------------------------------------------------------------------------------------
# A9 static scans
# ------------------------------------------------------------------------------------------------------------
def A9():
    lit = os.path.join(ROOT, 'step2', 'literature')
    tex = open(os.path.join(lit, 'ip_targets.tex')).read()
    declared = dict(re.findall(r'(?:arXiv:)?((?:hep-th/)?\d{4}\.\d{4,5}|hep-th/\d{7})v(\d+)', tex))
    stamps = {}
    for d in sorted(glob.glob(os.path.join(lit, 'ip_sources', '*'))):
        t = open(os.path.join(d, 'paper.txt'), errors='ignore').read(20000)
        m = re.search(r'arXiv:((?:[a-z-]+/)?[\d.]+)v(\d+)', t)
        if m: stamps[m.group(1)] = m.group(2)
    cmp = {k: {'declared': declared.get(k), 'local': stamps.get(k)} for k in sorted(set(declared) | set(stamps))}
    mismatch = {k: v for k, v in cmp.items() if v['declared'] and v['local'] and v['declared'] != v['local']}
    undeclared_local = [k for k, v in cmp.items() if v['local'] and not v['declared']]
    # absolute paths in engine code (cluster portability)
    abs_hits = {}
    for pat in ['step2/engine_a/*.py', 'step2/engine_b/*.py']:
        for f in glob.glob(os.path.join(ROOT, pat)):
            txt = open(f).read()
            n = len(re.findall(r"/Users/libotao", txt))
            if n: abs_hits[os.path.relpath(f, ROOT)] = n
    # plan items touching S_- (x) W^* (Theta_n image of S_+ (x) W) and the G_B open item
    plan_ip = json.load(open(os.path.join(ROOT, 'step2', 'theory', 'engine_plan_ip.json')))
    plan_o = json.load(open(os.path.join(ROOT, 'step2', 'theory', 'engine_plan.json')))
    sminus = []
    for cpt in plan_o['computations']:
        s = json.dumps(cpt)
        if re.search(r'S_- \(x\) W\^\*|xibar', s): sminus.append(cpt['id'])
    ip_refs = sorted(set(re.findall(r'\bE\d+\b', json.dumps(plan_ip))))
    gb_open = [x for x in plan_o['open_items'] if 'G_B' in x or 'Re a' in x]
    hep0307015_in_targets = 'hep-th/0307015' in tex or '0307015' in tex
    OUT['A9_static'] = {'version_comparison': cmp, 'version_mismatches': mismatch, 'local_copies_not_in_source_keys': undeclared_local,
                        'absolute_path_occurrences_in_engine_code': abs_hits,
                        'odd_plan_items_with_S-W*_generators': sminus, 'E_items_referenced_by_engine_plan_ip': ip_refs,
                        'odd_plan_open_item_on_G_B': gb_open,
                        'hep-th/0307015 (IP7 target) recorded in ip_targets.tex': hep0307015_in_targets}
    return True

if __name__ == '__main__':
    summary = {}
    for name, fn in [('A1', A1), ('A2', A2), ('A3', A3), ('A4', A4), ('A5', A5), ('A6_A7', A6_A7), ('A8', A8), ('A9', A9)]:
        try:
            summary[name] = bool(fn())
        except Exception as ex:
            summary[name] = 'error: ' + repr(ex)
        print(name, summary[name], flush=True)
    OUT['summary'] = summary
    with open(os.path.join(HERE, 'ip_audit2_checks.json'), 'w') as fh:
        json.dump(OUT, fh, indent=1, default=str)
