"""Round-2 revision of step2/theory/engine_plan.json (design odd_components_design.tex, round 2, 2026-10-29).
Loads the plan, applies the changes listed in step2/qa/design_revision_disposition.json (round_2), writes it back.
Deterministic; no inputs other than the plan itself."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, '..', '..', 'theory', 'engine_plan.json')
d = json.load(open(P))

d['date'] = ('2026-10-29 (round-2 revision after the audits odd_rev_audit_algebra, odd_rev_audit_graphs and the K judgement; '
             'aligned with manuscript ch08_definitions.tex; dispositions in step2/qa/design_revision_disposition.json, key round_2)')
d['binding_inputs'] = [
    'DESIGN_STEP2.md (Parts C, D; directives R0-R8; R7 supersedes the IR wording of R5/R6)',
    'DESIGN.md (rule 4 as amended 2026-10-28: transfer endpoint T allowed), DEFINITIONS_v2.md, STORY.md',
    'manuscript ch08_definitions.tex (binding definitions of Step 2: Def c8-PT ... Prop c8-Tsafe, Def c8-K, Thm c8-K, Cor c8-gram; corrections proposed in step2/qa/ch08_corrections_from_design.json)',
    'manuscript ch01-ch07',
    'step2/literature/q1_targets.tex (targets T1-T12, families F0-F5)',
    'step2/notation/esw_to_weinberg_srednicki.tex (dictionary)',
    'step2/theory/cubic_structure.tex (Part C: G_1..G_5, H_1..H_6, Q of the gaugino cubic term)',
    'step2/theory/k_transfer.tex and step2/qa/k_judgement.tex (route (A) for K_1, T-independence, logarithms)',
    'step2/theory/delta_theory.tex (read only: composite remainders, Prop prop:d-P, prop:d-T, prop:d-z)',
    'step2/engine_b/delta_results.json part2 and step2/qa/compare_delta.json (P_1 = -(1/3+f_P) unit on the scalar bilinears)'
]
c = d['constraints']
c['machine'] = ('24 GB shared Mac; at most ONE computational process at a time per agent, each < 1 GB RSS and preferably < 10 min; '
                'never /tmp (work inside step2/); kill what you start; no Wolfram')
c['cluster_portability'] = ('every engine job: pure Python 3 (numpy, sympy, mpmath, fractions) or C99 (libm, GMP/MPFR only); deterministic '
                            '(fixed seeds, exact rationals or mpmath at a fixed stated precision, no wall-clock or thread-order dependence); '
                            'one process per job, single-threaded; inputs and outputs by relative paths under step2/; every job writes one JSON '
                            'with inputs, code version (sha256 of the script), precision, timing and peak RSS; the memory bound of each job is '
                            'stated in its "resources" field and checked with /usr/bin/time -l (macOS) or /usr/bin/time -v (Linux)')
c['rules'] = c['rules'] + ('; order of limits along K_b and along (eta1,eta2), (K_b,x): eps -> 0 at fixed T first, then T -> infinity '
                           '(ch08 Def c8-K; defect graphs at T = infinity and fixed eps are never used); complete composite sums '
                           '(A_V--c lines from (Q^GF)^2 != 0) before any locality or Q-exactness statement (design prop:o-P1, rem:o-Kcomposite)')
c['reduction_rule'] = c['reduction_rule']

cv = d['conventions']
cv['names_ch08'] = ('objects and names of ch08: xi in Sigma, eta in Sigma^*, K_b: x -> (x,x)b - 2(b,x)x, P_b, P_{b\'}, Delta = 2E, '
                    's = d_CE + Q + {I_psl,-}, L = d_CE + Lambda, ad_f = {f,-}, F_psl = hbar F_1 + ..., F^(k)_psl arity parts (degree 0, '
                    'Weyl weight 4), psi scheme for operators (Def c8-renorm), Q^xi_0, Q^xi_1 (Def c8-twisted), N = {S^(1)(R_xi),-}; '
                    'the classical transformation along x is ad_{S^(1)(x)} O = {S^(1)(x),O} (the symbol delta_x is retired)')
cv['weights'] = ('Weyl weights: fields/antifields/derivatives as in ch07; explicit position -1; argument x of a cochain: half its '
                 'ad_Delta eigenvalue (P_b -1, K_b +1, Sigma -1/2, Sigma^* +1/2, rest 0); value weight = cochain weight minus the '
                 'argument weights; finite counterterms have weight 4 (CH08 correction to Def c8-quant)')
cv['sign_rule_products'] = ('products of functional-valued and of operator-valued cochains: the shuffle rule of ch02 Def c2-dce with '
                            'total parities (degree + fermion number mod 2) (CH08 correction to Thm c8-L (iv), Lemma c8-ev); with it the '
                            'product is associative and graded commutative and r_xi is multiplicative (odd_rev_audit_algebra a1, a1b, a1c)')
cv['special_conformal'] = 'K_b: v^m = x^2 b^m - 2 (b.x) x^m, div v = -8 b.x; (1/4) div v = -2 b.x'
cv['evaluation_along_xi'] = cv['evaluation_along_xi'].replace(
    'the coefficient -1/2 for L_1^(2) requires the decalage signs of hypothesis (H1) (thm:o-twist), fixed in E0',
    'the coefficient -1/2 for L_1^(2) follows from ch08 Lemma c8-ev with the total-parity sign rule (former hypothesis (H1) of thm:o-twist)')
cv['outer_derivation_B'] = ('B = ad diag(0_4,1_4) = -ad diag(1_4,0_4) in gl(4|4) (ch08 Def c8-B; the former wording ad diag(1_4,0_4) '
                            'was the opposite sign): +1 on S_+ (x) W (xi_alpha^i) and S_- (x) W (etabar^{alphadot i}), -1 on S_- (x) W^* '
                            '(xibar_{alphadot i}) and S_+ (x) W^* (eta^alpha_i), 0 on the even part; G_B = d/dt (e^{tB})^* S_psl: G_B^(0) = 0, '
                            'G_B^(1)(x) = S^(1)(Bx), G_B^(2)(x1,x2) = S^(2)(Bx1,x2) + S^(2)(x1,Bx2); s-closed, not in C(S - S^(2)) + Im s '
                            '(design prop:o-GB, proof via the translation component G_B^(1)(P_b) = 0)')
cv['R_rotation'] = ('N = {S^(1)(R_xi),-}, r^i_j = 4/3 (w^i n_j - 1/4 delta^i_j); [N, Q^xi_0] = {-S^(1)(xi) + S^(2)(xi,xi), -}, '
                    '-S^(1)(xi) + S^(2)(xi,xi) = image of S^xi under the degree derivation; {S^(1)(R_xi),S^xi} = Q^xi_0 S^(1)(R_xi) '
                    'in Im Q^xi_0; r_xi(G_B) = -Q^xi_0 S^(1)(R_xi); the parameter a of targets F1 equals minus the a of design prop:o-laws (e)')
cv['rescaling_directions'] = {
    'exact_direction': ('phi = kappa S^(1)(Delta) (F -> F + hbar s phi, = finite part kappa ad_{S1(Delta)} on all operators): components along '
                        'P_b, K_b, xi, eta change by -2 kappa hbar, +2 kappa hbar, -kappa hbar, +kappa hbar times the classical transformation; '
                        'D_1, M_1, R_1, L_1^(0) unchanged; (s phi)^(1)(P_b) = +2 kappa S^(1)(P_b), i.e. f_P -> f_P + 2 kappa hbar; reflection '
                        'compatible for real kappa (step2/qa/design_revision_round2/scheme_shift_round2.json; odd_rev_audit_algebra a3 part 4)'),
    'family_finite_part': ('psi\' = (c/2) ad_{S1(Delta)} pi_F on a family F stable under ad_{S1(Delta)}, ad_{S1(P_b)} (pi_F commuting with both): '
                           'components along P_b, K_b, xi, eta of operators in F change by -c hbar, +c hbar, -c/2 hbar, +c/2 hbar times the classical '
                           'transformation, plus (c/2) hbar ad_{S1(Delta)} (1 - pi_F) ad_{S1(x)} O along odd x when ad_{S1(x)} leaves F; D_1, M_1, R_1 '
                           'unchanged; reflection compatible for real c. "c times the Weyl weight of the density" is not a map of local functionals '
                           '(it does not commute with integration by parts in omega) and is not used'),
    'two_component': 'components along P_b, K_b, Q^i_alpha, Qbar_{alphadot i}, S^alpha_i, Sbar^{alphadot i}: -2kappa, +2kappa, -kappa, -kappa, +kappa, +kappa (times hbar and the classical transformation)'
}
cv['translation_component'] = ('P_1(O) = -f_P ad_{S1(P_b)} O - {G_{P_b},O} + R_{O,v} + [(Q+{I,-}),psi]_{P_b}(O) (delta_theory prop:d-P), f_P = coefficient of '
                               'S^(1)(P_b) in F^(1)_psl(P_b); translation-covariant scheme (f_P = 0): P_1 = -(1/3) unit ad_{S1(P_b)} on the scalar '
                               'bilinears (engine B part2, compare_delta: three routes), 0 on the letters; P_1 = 0 is a renormalisation condition '
                               '(design prop:o-P1); two-component on 1/16 int omega (w1^ij A_ij, w2^kl A_kl): -(1/3) unit 1/16 int omega b^m d_m(...)')
cv['transfer_endpoint'] = cv['transfer_endpoint'] + ('; momentum space in scaled Schwinger time: special edge tau in [0,1], P-edges tau in [1, T/eps); '
                                                     'composite line: [d_P,[d_P,S]] = d_{K\'\'}, K\'\' = 2 P* v. P*, fibre factors 2 (,)^{-1} and -(,)^{-1}; '
                                                     'T d/dT of the K_b composite kernel -> 3/(8 pi^2) b (x) (,)^{-1} (= 6 in units hbar/(16 pi^2)) (K judgement)')

# completed checks
cc = d['completed_checks']
cc += [
 {'id': 'KJ', 'what': 'K judgement: T-independence of the eps -> 0 limit of the K-route defect graphs at fixed T (bound T d/dT gamma = O(eps)); composite-line coefficient 3/(8 pi^2) with sign derived; zero-mode trees in four sectors (three realisations); log T of the defect: 0 in two-point action sectors, O(eps) for bilinears, c_eta eps e^{-eps p^2} at arity two; relations (K_b,P_c), (P_c,Delta), (eta1,eta2) with P_1 = p ad',
  'result': 'route (A): K_1 defined by the K-route graphs at finite T, T -> infinity trivial for K_1; flows and O[T]^(1)(K_b) keep log T times Q-closed trees; K_1 == (1/4) div v gamma O - p ad_{S1(K_b)} O (conditional)',
  'files': ['step2/qa/k_judgement.tex', 'step2/qa/k_judgement.json', 'step2/theory/k_transfer.tex', 'step2/qa/k_judgement/indep/', 'step2/qa/k_transfer_numeric/kn/kn_engine.py']},
 {'id': 'AU2', 'what': 'round-2 audit scripts: decalage/product/bracket sign rules (a1, a1b, a1c), BV signs (a2), B block table, G_B non-exactness system, inner-Delta shifts (a3), cubic family with H7, H8 (a4), P_1 scheme bookkeeping (a5), two-component identities (a6), truncation with explicit position (a7), floating-pair census (census_ir), orientation sign of two-ghost-end chains (koszul_orientation), multiplier analysis with P_1 = c ad (scheme_shift_check)',
  'result': 'all rerun in step2/qa/design_revision_round2/rerun/: outputs byte-identical; peak RSS <= 76 MB, <= 10 s each',
  'files': ['step2/qa/odd_rev_audit_algebra/', 'step2/qa/odd_rev_audit_graphs/', 'step2/qa/design_revision_round2/rerun/']},
 {'id': 'R2a', 'what': 'arity-zero cubic-order freedom with H7, H8+, H8-: closed space C (dim 9: 8 even, 1 odd), span V of the five s(N_f) and S (dim 6: 5 even, 1 odd), C = V + Part-C family, dim C/V = 3 (all even), representatives (b1,H1)=(1,1), (b2,H4)=(1/4,1), (b3,H5)=(-1/2,1); Yukawa = -1/2 x first, H5+H7 = 1/2 first + third (mod V)',
  'result': 'exact (sympy rationals), 0.3 s, 60 MB',
  'files': ['step2/qa/design_revision_round2/cubic_quotient.py', 'step2/qa/design_revision_round2/cubic_quotient.json']},
 {'id': 'R2b', 'what': 'sign of the rescaling directions on scalar densities: order reversal reproduces ch08 Cor c8-homotopy (ii); psi\' = kappa ad_{S1(Delta)} shifts the components along P_b, K_b, so(V), Delta by -2 kappa, +2 kappa, 0, 0; ad_{S1(Delta)} commutes with derivatives of the density, c x (density weight) does not',
  'result': 'exact (sympy), 1.7 s, 69 MB',
  'files': ['step2/qa/design_revision_round2/scheme_shift_round2.py', 'step2/qa/design_revision_round2/scheme_shift_round2.json']}
]

comps = {x['id']: x for x in d['computations']}
E0 = comps['E0']
E0['tasks'] = [t for t in E0['tasks'] if not t.startswith('implement r_xi')] + [
    'implement r_xi as the full sum over k with coefficients (-1)^{k(k-1)/2}/k!, and the C^bullet-linear extension of operator maps and the products of W(P,J) with the total-parity sign rule (ch08 Lemma c8-ev, Thm c8-L (iv) as corrected); verify multiplicativity of r_xi on random cochains with fermionic values',
    'canonical form of two-leg outputs in both engines extended from alpha delta + beta^mu d_mu delta to alpha delta + beta^mu d_mu delta + gamma^{mu nu} d_mu d_nu delta (S^(1)(K_b) and S^(2)(eta1,eta2) are quadratic in the position); test on {S^(1)(K_b), int(A,A)} and on {S^(1)(K_b), int (A_W,A_W)} against the classical transformation of design thm:o-KD (two-component v.A_ij = v^n d_n A_ij - 2 (b.x) A_ij)',
    'composite-line kernels P* v. P* for P_b, so(V), Delta, K_b in both engines (finite T), with the fibre factors 2 (,)^{-1} (A (x) A) and -(,)^{-1} (A*_V (x) c)']
E0['resources'] = 'seconds to minutes; < 300 MB; numpy/sympy/fractions; single process'
E1 = comps['E1']
E1['graph_classes']['absent'] = [s.replace('T: delta_xi O_zeta', 'T: ad_{S1(xi)} O_zeta') for s in E1['graph_classes']['absent']] + ['T0: O_zeta gauge invariant']
E1['resources'] = E1.get('resources', '') + '; one process; < 800 MB per graph batch; deterministic'
E2 = comps['E2']
E2['tasks'] = [t.replace('delta_', 'ad_S1') for t in E2['tasks']] + [
    'exact direction kappa S^(1)(Delta) and family finite parts (c/2) ad_{S1(Delta)} pi_F: verify c_K unchanged, c_B\' shifted through kappa_B B_w']
E2['resources'] = 'minutes; < 300 MB'
E3 = comps['E3']
E3['tasks'] = [t for t in E3['tasks'] if not t.startswith('E3b')] + [
    'E3b (corrected, design prop:o-arity0): with the three-leg part of F^(0) in span{H1..H6, H7 = 1/2 int(lam,rho([A_V,lam])), H8+, H8-} (QH7 = -G5, QH8+- = 0) the closed non-equivariant cubic-order freedom C has dim 9 (8 even); the five s(N_f) and S span V (dim 6); dim C/V = 3, representatives (b1,H1) = (1,1), (b2,H4) = (1/4,1), (b3,H5) = (-1/2,1); show that none of the three extends to an s-closed psl(4|4)-cochain (solve {S^(1)(x1), G^(0)} = (-1)^{|x1|} (Q+{I,-}) G^(1)(x1) for x1 = xi, eta, P_b with G^(1) in the arity-one primitives of degree -1), first checking at quartic order which of them extend non-equivariantly (gauge invariance); test hypothesis (H4)',
    'F^(1)(xi) family: contains a S^(1)(xi) (G_B) and kappa S^(1)(xi) from the exact direction s(kappa S^(1)(Delta))']
E3['resources'] = 'tens of minutes; < 1 GB; one process'
E4 = comps['E4']; E4['resources'] = 'minutes; < 500 MB'
E4['tasks'] = [t.replace('G_B and the field redefinitions', 'G_B, the exact direction, the family finite parts and the field redefinitions') for t in E4['tasks']]
E5 = comps['E5']; E5['resources'] = 'minutes; < 500 MB'
E6 = comps['E6']
E6['tasks'] = E6['tasks'] + [
    'letters are gauge variant: include class (T0) (marked vertex (Q+{I,-})O with a vertex S^(1)(x1)) and read off with the corrected design eq:o-L1graphs, subtracting the arity-one trees of hbar L_1(O)^(0) with one S^(1)(x1) (nonlocal; for int omega (X, iota_v A_V) the tree of int omega (X, iota_v dc) joined to S^(1)(xi) through c--A*_V)']
E6['resources'] = 'minutes; < 500 MB'
E7 = comps['E7']; E7['resources'] = 'tens of minutes; < 1 GB'
E8 = comps['E8']
E8['tasks'] = [
    'verify S_1(eta)(O) == Q_1(rho(x)eta)(O) - p ad_{S1(eta)} O + U_eta(O) mod s-exact on the scalar bilinears (design eq:o-SQ), with p the P_1 coefficient of E11 (translation-covariant scheme p = -1/3 unit: S_1 == Q_1(rho(x)eta) + 1/3 unit ad_{S1(eta)}; P_1 = 0 scheme: S_1 == Q_1(rho(x)eta)); U_eta = 0 on the bilinears; two-component (rho(x)eta)_alpha^i = -i x_m sigma^m etabar^i, (rho(x)eta)_{alphadot i} = i x_m eta^alpha_i sigma^m, the same -p along S^alpha_i and Sbar^{alphadot i}',
    'U_eta on the weight-3 Konishi descendants B_w, Fbar_w: coefficients of int omega (lam, rho(A|_W) zeta), zeta in Sigma (two-letter -> one-letter at one loop)',
    'commutant of xi (eta^a_i = u^a m_i with m.w = 0, or eta in S_- (x) W): the R_1 of BK footnote 24; check {Q^xi_0, S_1^xi(eta)} + {S_0^xi(eta), Q^xi_1} = 0 on H(Q^xi_0) (target T10)']
E8['prerequisites'] = ['E1', 'E3 (with eta in place of xi)', 'E11 (P_1 and the scheme)']
E8['resources'] = 'tens of minutes; < 1 GB; one process'
E9 = comps['E9']
E9['tasks'] = ['compute the arity-two component on (xi, eta) with (xi,eta) != 0 for the Konishi primary; verify that its d_CE term carries 1/2 (xi,eta) D_1 + M_1 + R_1 parts with D_1 from Part D; the (xi,eta) P-edge pairs: each orientation alone has degree 0 (S^3 average 2 eta.xi), so combine both orientations and both chirality assignments in ONE integrand before integrating at T = infinity (or sum at finite T, then T -> infinity)']
E9['resources'] = 'tens of minutes; < 1 GB'
E10 = comps['E10']
E10['title'] = 'K_1 (special conformal) at a finite transfer endpoint T: route (A)'
E10['tasks'] = [
    'E10a (done, K judgement Prop prop:kj-defectlog, prop:kj-compare (ii)): log T of the defect along K_b per composite graph c_g (e^{-eps p^2} - 1), c_g = +-9i, +-147/2 i, +-3, +-18 units; summed (3/8 pi^2) hbar Q(zero-mode tree at eps): 0 in every two-point action sector, O(eps) nonzero for the scalar bilinears',
    'E10b (done, K judgement Prop prop:kj-trees, Cor cor:kj-flowlog): zero-mode trees nonzero in (A_W,A*_W), (A_V,A*_V), (c,c*), (lam,lam*) (log T slopes -6i, 46i, 3i, 8i units) and for the bilinears: W^(1)(K_b) and O[T]^(1)(K_b) have no T -> infinity limit; arity two eta-eta: c_eta eps e^{-eps p^2} with c_eta = -48i, 288i, 0, 48i units in the defect',
    'E10c: K_1 by design def:o-K on the letters and on the scalar bilinears (traceless 20\', trace Konishi): eps -> 0 at fixed T of the K-route graphs (classes T, T0, A-E2 with x1 = K_b, edges P(eps,T)), complete composite sums (composite remainder R_{O,K}, kernel K_b.(P*)^2 + [P*,K_b.]P*), read-off by the corrected eq:o-L1graphs (fewest legs; letters: subtract the trees of hbar L_1^(0) with S^(1)(K_b)); exact procedure: for each graph prove that lim_{eps->0} at fixed T equals the eps^0 coefficient of the integrand in scaled Schwinger times with the P-edge upper limit T/eps replaced by infinity (the large-distance part is O(eps), K judgement Thm thm:kj-Tindep), check the convergence of the needed Taylor order at large tau, evaluate exactly (rationals, log 2, log 3, pi^2, Li_2 values) or with >= 40 digits + PSLQ on an explicit basis, confirmed by a second evaluation; numerical cross-check at T|p|^2 = 1, 4, 16 with eps|p|^2 -> 0; targets (design thm:o-KD): K_1 + p ad_{S1(K_b)} == (1/4) div v gamma O = -2 (b.x) gamma O; translation-covariant scheme (p = -1/3 on bilinears, 0 on letters): K_1(20\') == 1/3 unit ad_{S1(K_b)} O, K_1(Konishi) == 48 (b.x) unit O + 1/3 unit ad_{S1(K_b)} O, letters 8, 8, 16, 4 (b.x) unit O; P_1 = 0 scheme: 0, 48 (b.x) unit O; with f_P: p = -(1/3+f_P), -f_P',
    'E10d: relation at (eta1,eta2): K_1(K_{v1}) == [S_1(eta1), ad_{S1(eta2)}] + [ad_{S1(eta1)}, S_1(eta2)] (anticommutators), v1 = Gamma(eta1,eta2)|_V onto V; S_1 from E8 at T = infinity: determines K_1 on cohomology without endpoint',
    'E10e: K = P^dagger cross-check with the inner-product workflow (do not edit its files): K_1 along K_bbar = pbar ad_{S1(K_bbar)} + hbar [ad_{S1(K_bbar)}, Gram_tree^{-1} Gram_1] (K transfer Cor cor:ktr-gram); rescaling part of the commutator -(p + pbar) ad; real p removable by the reflection-compatible finite part (p/2) ad_{S1(Delta)} pi_F',
    'E10f: completion of the arity-one obstruction along K_b in the quadratic sector: add the one-particle-reducible graphs and the counterterm vertices -Q I^CT(K_b), -{S^(1)(K_b), I^CT}_eps to the 1PI special-edge sums (limits 22i, -141/2 i, 44/3, 16 units); check that the eps^{-1} parts 18i b (A_W) and -108i b (A_V) cancel against {S^(1)(K_b), eps^{-1} mass counterterm}; decompose the finite limits into s F^(1)(K_b) + {S^(1)(K_b), F^(0)}',
    'E10g: arity two along (K_b,P_{b\'}): T-independence at finite T and the values L_1^(2)(K_b,P_{b\'}) where needed for the relation at (K_b,P_{b\'}); along (K_b,x) the graphs are defined at finite T only (single floating K log T; chains (K,Delta), (K,M) log, (K,K\') T^{1/2})',
    'arity-one counterterm along K_b: singular parts at finite T (T-independent, design lem:o-Tflow (iv)); no eps^{-1} after the bracket with the eps^{-1} mass counterterm, log eps in the span of (v.A,A*), (v.lam,lam*), (v.c,c*)-type structures']
E10['prerequisites'] = ['E11 (P_1 and the scheme)', 'E8', 'Part D', 'E0 (finite-T propagator, extended canonical form, composite kernels)']
E10['resources'] = 'tens of minutes per sector; < 1 GB; one process per (operator, sector, T) job; deterministic'
E10['notes'] = ('no position cutoff R in any definition (R only as a diagnostic of log coefficients); the limit order is eps -> 0 at fixed T, then T -> infinity '
                '(trivial for K_1); defect graphs at T = infinity and fixed eps are not used; the finite-T propagator (e^{-eps p^2} - e^{-T p^2})/p^2 is used in every K-route job')
E11 = comps['E11']
E11['title'] = 'P_1, M_1, R_1: composite remainders and the scheme (run BEFORE E8 and E10)'
E11['tasks'] = [
    'P_1 on letters, scalar bilinears (traceless and trace) and one three-letter gauge-invariant operator: complete composite sums with the kernel v.(P*)^2, v = b (translation), (Q^GF)^2 A* = 2 d*(A*_V); expected (engine B part2, compare_delta, three routes) -(1/3) unit ad_{S1(P_b)} O on the bilinears and 0 on the letters in the translation-covariant scheme; test the pair conjecture p = -(n-1)/3 unit on operators with n letters (two-body kernel -(1/3)(d_i + d_j))',
    'M_1: the same along so(V) (S^(1)(M) contains (v.A_V, A*_V), composite line v.(P*)^2 with linear v); prediction on scalar families: rescaling part m = 0 on cohomology (design thm:o-KD, conditional)',
    'R_1: verify R_1 = 0 in the so(W_R^perp)-covariant scheme (no floating insertion; design prop:o-P1 (iv))',
    'arity-one counterterms and finite counterterms along P_b, so(V), so(W_R^perp) in the two-leg sectors; record f_P (S^(1)(P_b) coefficient of F^(1)(P_b)); the exact direction s(kappa S^(1)(Delta)) moves f_P by 2 kappa hbar and P_1 by -2 kappa hbar uniformly',
    'fix and record the scheme: (A) translation covariant (F^(1)(P_b) = 0, psi translation covariant); (B) P_1 = 0 by psi\' = (p_F/2) ad_{S1(Delta)} pi_F on each family F with P_1 = p_F ad_{S1(P_b)}; report Q_1, S_1, K_1 in both (shifts on F: -p_F/2, +p_F/2, +p_F times the classical transformation; on the bilinears +1/6, -1/6, -1/3 units relative to (A))',
    'G_B and the coupling direction leave P_1, M_1, R_1 unchanged (prop:o-GB (iii))']
E11['prerequisites'] = ['E0 (composite kernels)', 'Part D4 engines (step2/engine_a, step2/engine_b: reuse, do not edit their outputs)']
E11['resources'] = 'minutes to tens of minutes; < 1 GB; one process per (operator, generator) job'
E12 = comps['E12']; E12['resources'] = 'tens of minutes; < 1 GB'
order = ['E0', 'E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'E7', 'E11', 'E8', 'E9', 'E10', 'E12']
for i, k in enumerate(order):
    comps[k]['order'] = i
d['computations'] = [comps[k] for k in order]

g = d['graph_class_reference']
g['T'] = 'marked vertex ad_{S1(x1)} O (classical transformation), other vertices I: graphs of the operator ad_{S1(x1)} O'
g['T0'] = 'marked vertex (Q+{I,-})O (gauge-variant O), one vertex S^(1)(x1) on a cycle, tree or leg; trees, one loop, counterterm trees'
g['formula'] = ('hbar L_1(O)^(1)(x1)[T]_tree = lim_{eps->0} [ sum of A..E2 graphs + arity-one part along x1 of the tree flow of '
                'd_delta I^CT_psl[eps](s O) ] - (hbar L_1(O)^(0))[T]^(1)(x1) at fixed T (design eq:o-L1graphs; the last term: trees of '
                'hbar L_1^(0) with one S^(1)(x1), nonlocal, same leg count); for gauge-invariant O with L_1^(0) = 0 it reduces to the old '
                'formula; T = infinity where prop:o-Tsafe (i) applies; classes stated at psi = 0 (psi enters by prop:o-laws (b))')
g['placement'] = ('lam*, c*, A*_W never on P(eps,T); floating S^(1): even x with a vector field on V: A_V -> c via (v.A_V, A*_V) (composite line, '
                  'kernel P* v. P*); odd x: lam -> c via -(Gamma(eps(x),lam)|_V, A*_V); so(W_R^perp): never floating; two floating even insertions '
                  'adjacent through an A_V (x) A_V edge (chain c--A*_V c--x1--A_V A_V--x2--A*_V c--c, two ghost ends); two floating odd insertions '
                  'adjacent through lam (x) lam with opposite chiralities (P-edge pair); floating even and odd never adjacent; S^(2) never floating '
                  '(step2/qa/odd_rev_audit_graphs/census_ir.json)')
ld = d['large_distance_decay_of_floating_insertions']
ld['K_b'] = ld.pop('K_{v1}').replace('at finite T: log T with coefficient 3/(8 pi^2) times the cut tree',
            'at finite T: log T with coefficient 3/(8 pi^2) hbar times the zero-mode tree (constant modes A_V = b (x) X_a, c = X^a); per composite defect graph c_g (e^{-eps p^2} - 1), summed 0 in two-point action sectors and O(eps) for bilinears; flows: nonzero trees (K judgement)')
ld['even-even chain (two floating even insertions through A_V (x) A_V)'] = ('degree d1+d2-2 in R^8 per orientation (d = 0,1,1,2 for P, so(V), Delta, K), '
    'relative orientation sign -1 (two ghost ends), leading term cancels: one less after the orientation sum; (P,P), (P,M), (P,Delta), (M,M\'), '
    '(M,Delta) converge after the sum, (Delta,Delta) = 0 by antisymmetry; (K,P) converges as a chain but contains the single floating K log; '
    '(K,Delta), (K,M) log; (K,K\') T^{1/2}; O(eps) coefficients in defect graphs (koszul_orientation.json, census_ir.json)')
ld['P-edge pair Sigma x Sigma^*'] = ('degree 0 per orientation (S^3 average 2 eta.xi), cancels in the orientation sum (degree -1): converges; '
                                      'combine both orientations in one integrand before integrating at T = infinity')
ld['P-edge pair Sigma^* x Sigma^*'] = ('degree 1 per orientation; orientation sum leaves a log with S^3-averaged symbol 6i v1.(f0 g1 - f1 g0), '
                                        'v1 = Gamma(eta1,eta2)|_V; defect: c_eta eps e^{-eps p^2}, c_eta = -48i, 288i, 0, 48i units; O(1) in O[T]^(2) and W^(2), '
                                        'whose Q-image cancels the d_CE log term as eps -> 0; zero pointwise for v1 = 0')
d['open_items'] = [
    'normalisation map to BK and the transport Phi from the collections of ESW Thm 10.9 (design rem:o-transport): needed to compare c_K numerically with BK (3.66) (targets open item 1)',
    'kappa_B sign (E0)',
    'E11: the composite remainders along P_b and so(V) on a three-letter operator (pair conjecture p = -(n-1)/3); M_1 rescaling (predicted 0 on cohomology); choice and record of the scheme before E8, E10',
    'E10c: K_1 on letters and bilinears from the K-route graphs at finite T (the relations predict K_1 == (1/4) div v gamma O - p ad_{S1(K_b)} O); E10f: completion of the arity-one obstruction along K_b (1PR graphs, counterterm vertices)',
    'f_S = -1/3 (S^(1)(Delta)-coefficient of F^(1)(Delta) on the bilinears, engine B normalisation) predicted by the (P,Delta) relation d = 0: to be checked by the dilatation workflow at arity two',
    'hypothesis (H4) / E3b: the three directions of C/V (design prop:o-arity0 (iii)) admit no s-closed psl(4|4)-completion',
    'identification of the one-loop quantisation restricted by r_xi with BK\'s twisted quantisation up to C S^xi + Im Q^xi_0 (design thm:o-twist (iv)); tested on T3-T5, T12',
    'conjugation compatibility: Re a of G_B fixed by the Q^dagger = S test (ch08 Cor c8-gram); the rescaling directions are reflection compatible for real parameters (K transfer prop:ktr-rescale (ii) and K judgement cor:kj-gram state the opposite condition abar = -a; see step2/qa/design_revision_disposition.json round_2, item KJ-3)',
    'gauge-fixing law (c) of prop:o-laws: proof plan only',
    'proposed Chapter 8 changes: step2/qa/ch08_corrections_from_design.json (manuscript not edited)']
d['T_infinity_safe'] = {
    'safe': ['arity zero', 'arity one along Sigma, Sigma^*, Delta, translations, so(V), so(W_R^perp) (with complete composite sums; Q-primitive decompositions of composite contributions made at finite T, delta_theory prop:d-T (ii))',
             'arity two along (xi,xi\'), (xi,eta) (orientation sum in one integrand), (eta1,eta2) with v1 = 0, pairs of translations, so(V), Delta (after the orientation sum)'],
    'constant_in_T': ['counterterms, obstruction and finite counterterm along K_b', 'K_1', 'L_1^(2)(eta1,eta2), L_1^(2)(K_b,P_{b\'})',
                      'd_T O[T]^(1)(K_b) and Q(O[T]^(1)(K_b)) for gauge-invariant O with L_1^(0) = 0'],
    'not_safe_finite_T_only': ['W(P(eps,T),I_psl)^(1)(K_b) and O[T]^(1)(K_b) for the scalar bilinears (3/(8 pi^2) hbar log T x nonzero tree)',
                               'W^(2), O[T]^(2) along (eta1,eta2) with v1 != 0', 'defect graphs at fixed eps along K_b (bilinears) and (eta1,eta2) with v1 != 0: O(eps) log T',
                               'arity two along (K_b, x) for every even x'],
    'rule': 'components in the last two lists are defined at finite T (design def:o-K; ch08 Def c8-K) with eps -> 0 first; K_1 is T-independent (design prop:o-Kclaim, K judgement Thm thm:kj-Tindep), so T -> infinity is trivial for it'}
json.dump(d, open(P, 'w'), indent=1, ensure_ascii=False)
print('ok', [x['id'] for x in d['computations']])
