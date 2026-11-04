"""Reviser: step2/strict/results/action_final.json and operator_final.json (formula-centred, both notations, status per statement)."""
import json, os
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
J = lambda p: json.load(open(os.path.join(ROOT, p)))
ob = J('revise.nosync/onebody_F0.json')
DC = 'derived-from-construction'; C2 = 'computed by two engines + referee'; C1 = 'computed by one engine + referee'
CJ = 'conjectured'; OP = 'open'; XC = 'cross-check only'
U = 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); L = log(3/2)'
WS = 'Weinberg-Srednicki two-component (Def c8-twocomp; sigma^m = (vec sigma, i), sigmabar^m = (-vec sigma, i), eps^12 = eps^{dot1 dot2} = 1); A^{ij} = (1/2) eps^{ijkl} A_kl'

action = {
 'file': 'step2/strict/results/action_final.json', 'date': '2026-11-04', 'units': U,
 'sources': ['step2/strict/action/theory/action_theory.tex (revised, 28 pp, _build/action_theory.pdf)',
             'step2/strict/action/engine_a, engine_b (gates_revised.json)', 'step2/strict/action/complete (results.json, gates.json 73/73)',
             'step2/strict/action/qa (compare.json, findings.json)', 'step2/strict/disposition.json'],
 'conventions': {
   'bracket': {'statement': '{f,g} = Delta0(fg) - Delta0(f)g - (-1)^|f| f Delta0(g) = (-1)^|f| {f,g}_std (Def c1-bracket); a degree-1 derivation D (Q, {I,-}, Lambda, {F,-} with |F| = 0) obeys D{a,b} = -{Da,b} - (-1)^|a| {a,Db}', 'status': DC},
   'component_equations': {'statement': 'Lambda F1^(0) = Ob^(0);  Lambda F1^(1)(x) + {F1^(0),S(x)} = Ob^(1)(x);  Lambda F1^(2)(x,y) + F1^(1)([x,y]) + (-1)^|x| ({F1^(1)(x),S(y)} + {S(x),F1^(1)(y)}) + {F1^(0),S^(2)(x,y)} = Ob^(2)(x,y)', 'status': DC},
   'Wess_Zumino': {'statement': 'Lambda Ob^(1)(x) = -{Ob^(0),S(x)};  Lambda Ob^(2)(x,y) = Ob^(1)([x,y]) - (-1)^|x| {Ob^(1)(x),S(y)} + {S(x),Ob^(1)(y)} - {Ob^(0),S^(2)(x,y)}  (revised 2026-11-04; the earlier printed signs fail on every computed entry)', 'status': DC + '; checked: ' + C2},
   'odd_odd_generating_cochain': {'statement': 'J = I + c S(xi) + c\' S(xi\') - c c\' S^(2)(xi,xi\'), sign forced by the classical master equation: Q S^(2)_2 + S_2(P_{Gamma|V}) - {S_2(xi),S_2(xi\')} = 0', 'status': C2},
   'log_coefficient': {'statement': 'I^CT[eps] = (9/4) u eps^-1 int(A,A) + log(eps) lc;  s lc = 0;  lc^(0) = u Lambda N, N = 2 int(A,A*) - 4 int(lam*,lam) - int(c*,c);  closed combination lc^(1)(x) - {uN,S(x)}', 'status': DC},
 },
 'arity_one': {
  'Ob0': {'value': 'hbar Ob^(0)|_2 = -27/8 u int(d*dc, d*A_V);  hbar Ob^(0)|_3 = u[(11/4 + 19/4 L)(G1+G3) + 27/8 G2 - 9 L G5]',
          'two_component': '-27/8 u int d^4x (d^2 c, d_m A_m)', 'status': C2},
  'F1_0': {'representative': 'F1^(0)|_2 = 27/16 u int(d*A_V, d*A_V) + s-closed arity-zero directions (coupling, beta directions)', 'two_component': '27/16 u int d^4x (d_m A_m, d_n A_n)', 'status': C2},
  'Delta': {
    'Ob1': 'hbar Ob^(1)(Delta)|_2 = u[-11/3 int(dA,dA) + 47/3 int(lam,dslash lam) - 14/3 int(dc,A*)], position part 0',
    'Ob1_two_component': 'u int d^4x { -11/3 [1/2 (F^lin_mn,F^lin_mn) + 1/4 (d_m A^{ij}, d_m A_ij)] + 94/3 i (lambdabar_{alphadot i}, sigmabar^{m alphadot alpha} d_m lambda_alpha^i) - 14/3 (d_m c, A*_m) }',
    'F1_representative': 'u[22/3 int(A,A*) - 47/3 int(lam*,lam) - 8/3 int(c*,c)] + u sum_d theta^Delta_d S(Delta)_d + kappa_V u Q int(iota_x A_V,c*) + kappa_lam u Q int(rho(x)lam*,lam*)',
    'F1_two_component': 'u int d^4x {22/3[(A_m,A*_m) + 1/4 (A_ij,A*^ij)] - 47/3[(lambda_alpha^i,lambda*^alpha_i) - (lambdabar_{alphadot i},lambdabar*^{alphadot i})] - 8/3 (c*,c)} + u sum_d theta^Delta_d S(Delta)_d + kappa-directions',
    'free_directions': ['theta^Delta_perp, theta^Delta_lam, theta^Delta_V (fixed: -1/3 each)', 'kappa_V, kappa_lam (Q-exact, scheme)', 'parity-odd class excluded by reflection'],
    'pole': '(9/4) u {int(A,A),S(Delta)} = -9 u int(A,A)', 'status': C2},
  'P_b': {'Ob1': 'hbar Ob^(1)(P_b)|_2 = 0 (direct cycles)', 'F1_representative': 'u sum_d theta^P_d S(P_b)_d + exact; theta^P uniform (= t_P) by (xi,P_b); t_P = t_Delta = -1/3 in a translation-covariant scheme', 'status': C2},
  'M_omega': {'Ob1': 'hbar Ob^(1)(M_omega)|_2 = 1/3 u Q S(M_omega)_spin, S(M_omega)_spin = int omega_{mu nu}(A_mu,A*_nu) + int((1/4) omega_{mu nu} rho_mu rho_nu lam, lam*); 0 on (A_perp,A_perp)',
              'F1_representative': '1/3 u S(M_omega)_spin + u sum_d theta^M_d S(M_omega)_d; theta^M uniform = t_M = t_Delta = -1/3', 'two_component': 'S(M)_spin = int d^4x {omega_mn (A_m,A*_n) + 1/4 omega_mn [(lambda*^alpha_i (sigma^m sigmabar^n)_alpha^beta lambda_beta^i) - (dotted block with sigmabar^m sigma^n)]}', 'status': C2 + ' (lambda sector: engine B data fitted by the referee replay; engine A in its own realisation, magnitude 1/3 x 1/4)'},
  'R_r': {'Ob1': 'Ob^(1)(R_r) = 0 on all two-leg sectors and three tested three-leg sectors', 'F1_representative': 'u sum_d theta^R_d S(R_r)_d; theta^R uniform = t_R = t_Delta + 1/3 = 0', 'status': C2},
  'xi': {'Ob1': 'hbar Ob^(1)(xi)|_2 = u[(79/24 - L/2)(g1 + g3) + 27/8 g2 + 0 g4] (both chiralities); pole (9/2)[int(Gamma_V,A_V) + int(Gamma_perp,A_perp)] = (9/4){int(A,A),S(xi)}',
         'Q_images': 'Q T1 = -Q T3 = -1/2 g1, Q T2 = -Q T4 = -1/2 g3, Q T0 = -(g2 + g4); {int(d*A_V,d*A_V),T1} = 2 g2; H^0(Q)|_2 = 0; ker Q at degree -1 = span{S^V_2 = T1+T3, S^perp_2 = T2+T4}',
         'F1_representative': 'F1^(1)(xi)|_2 = u(L - 79/12)(T1 + T2) + u theta^xi_V S^V_2(xi) + u theta^xi_perp S^perp_2(xi) (+ beta-shifts {N_f,S_2(xi)}); no T0 term',
         'F1_two_component': 'T1 = -i int (lambdabar_{alphadot i} sigmabar^{m alphadot alpha} xi_alpha^i, A*_m), T2 = -1/2 int eps_{jikl}(xi_alpha^j eps^{alpha beta} lambda_beta^i, A*^{kl}), T3 = -1/4 int ((d_m A_n - d_n A_m) sigma^m sigmabar^n xi, lambda*), T4 = -i/2 int (eps^{alpha beta} xi_beta^j sigma^m_{alpha alphadot} d_m A_ij, lambdabar*^{alphadot i})',
         'free_directions': 'theta^xi_V = theta^xi_perp =: t_xi forced by (xi_pm, xi\'_pm); t_xi+ + t_xi- = t_P + 27/4 + 2L; individual t_xi free (inner and outer directions)',
         'log_coefficient': 'lc^(1)(xi)|_2 = 2u(T1 + T2 - T3 - T4) = {uN,S_2(xi)}; l_xi = 0 on two legs', 'status': C2},
  'eta': {'Ob1': 'hbar Ob^(1)(eta)|_2 = x^mu hbar Ob^(1)(rho(e_mu)eta)|_2 + u[-5/2 int(eta,lam) d*A_V - 5/2 int(eta,rho(dA_V)lam) + (-3/2 + 3L) int(eta,rho(dA_perp)lam) + 5/2 int(rho(dc)eta,lam*)] + Q(.) (canonical density)',
          'F1_representative': 'sum_mu x^mu u(L - 79/12)(T1+T2)(rho(e_mu)eta) - (3L - 3/2) u s0 - 5/2 u u0 + u theta^eta_V S^V(eta) + u theta^eta_perp S^perp(eta); u0 = int(rho(A_V)eta,lam*); S^V(eta) = x.(T1+T3)(rho eta), S^perp(eta) = x.(T2+T4)(rho eta) + s0',
          'cohomology': 'H^0(Q) = 1 at weight 7/2, p <= 1 ((A_perp,lam) sector); residual has no component on it; kernel at degree -1 = span{S^V(eta), S^perp(eta)}',
          'free_directions': 'theta^eta_V = theta^eta_perp =: t_eta (forced by (xi,eta) with theta^xi uniform); t_xi+ + t_eta+ = t_Delta + 27/4 + 2L',
          'two_component': 's0 = int (A^{kj} eps_{alpha beta} eta^beta_j, lambda*^alpha_k), u0 = -i int (A_m eta^alpha_i sigma^m_{alpha alphadot}, lambdabar*^{alphadot i})', 'status': C2},
  'K_b': {'Ob1': 'in the span of canonical data with position degree <= 2; pole sum (9/4){int(A,A),S(K_b)}; WZ1 corrected sign holds',
          'F1_representative': 'u int (b,x){-22/3 (A,A*) + 16 (lam,lam*) + 8/3 (c,c*)} + 2/3 u int{(iota_b A)(iota_x A*) - (iota_x A)(iota_b A*)} - 1/3 u int b_mu x_nu (lam, rho_mu rho_nu lam*) + u sum_d theta^K_d S(K_b)_d + sum_{j=1..4} kappa^K_j G^K_j',
          'two_component': 'u int d^4x { (b.x)[-22/3((A_m,A*_m) + 1/4 (A_ij,A*^ij)) + 16((lambda_alpha^i,lambda*^alpha_i) - (lambdabar,lambdabar*)) + 8/3 (c,c*)] + 2/3 [(b_m A_m)(x_n A*_n) - (x_m A_m)(b_n A*_n)] - 1/3 b_m x_n [rho_m rho_n -> sigma^m sigmabar^n on undotted, sigmabar^m sigma^n on dotted] }',
          'cohomology': 'H^0 = 4 (residual has no component), H^-1 = 7: S(K_b)_perp, S(K_b)_lam, S(K_b)_V and four second-order free-field symmetries G^K_j (not exact, not fixed on two legs)',
          'structure': 'x-linear part = -(b,x) F1part(Delta) + spin term of the rotation part of d_b\' v (16 - 1/3 = 47/3): translation covariant',
          'T_independence': 'structural (Lemma at-special: eps^0 simplex coefficients); no numerical finite-T run', 'status': C2},
 },
 'arity_two_two_legs': {
  '(Delta,P_b)': {'Ob2': 'perp: -2/3 int(d_b A_perp, A*_perp); lam: 2/3 int(d_b lam,lam*) + (-3/2 + log(4/3)) int(rho(b) dslash lam, lam*); V: -2/3 int(d_b A_V,A*_V) + (3/2 - log(4/3)) int(d_mu(b.A_V),A*_mu) + (-5/6 + log(4/3)) int(d_b c,c*); no position part, no pole',
                  'class': 'o_d = -2/3 u for d = perp, lam, V', 'condition': 'theta^Delta_d = o_d/(2u) = -1/3', 'status': C2},
  '(Delta,xi_pm)': {'Ob2': 'hbar Ob^(2)(Delta,xi)|_2 = u[8(T1+T2) - 26/3(T3+T4) + log(4/3) Q(psi1 + psi2)], psi1 = int x^mu (Gamma(xi,lam)^mu, c*), psi2 = int x^mu (A*_mu,(xi,lam*))',
                    'condition': 'solvable iff theta^Delta_perp = theta^Delta_lam = theta^Delta_V = -1/3; theta^xi free', 'F1_2': 'F1^(2)(Delta,xi)|_2 = log(4/3) u (psi1 + psi2) + closed degree -2 directions', 'status': C2 + ' (engine A flag lifted: gate artefact, finding F2)'},
  '(xi_+,xi_-)': {'condition': 'theta^P_d = theta^{xi+}_V + theta^{xi-}_V - 27/4 - 2L (every d); theta^{xi+}_V + theta^{xi-}_V = theta^{xi+}_perp + theta^{xi-}_perp', 'H_minus1': 'V 1, perp 1, lam 3 (incl. two SU(4)-15 free classes)', 'status': C2},
  '(xi_pm,xi\'_pm)': {'condition': 'theta^xi_V = theta^xi_perp', 'status': C2},
  '(xi_pm,P_b)': {'Ob2': '(1/4) log(4/3) x Q-exact functional, no rational part, no pole', 'condition': 'theta^P_V = theta^P_perp = theta^P_lam', 'status': C2},
  '(xi_+,eta_-)': {'condition': 'theta^xi_V - theta^xi_perp = theta^eta_V - theta^eta_perp (parametrisation of arity_one.eta)', 'status': C2},
  '(xi_+,eta_+)': {'condition': 'theta^R_perp = theta^R_lam =: tau; theta^Delta_d = theta^M_d = tau - 1/3; theta^xi_d + theta^eta_d = tau + 77/12 + 2L (d = V, perp)', 'bracket': '[xi,eta] = (xi,eta)/2 Delta + M_omega + R_r (gate G-bracket)', 'status': C2 + ' (basis of (lam,lam*) position sectors one short per chirality: solvability decided on value vectors; pole checked on intrinsic part)'},
  '(K_b,P_b\')': {'status': OP}, '(eta,eta\')': {'status': OP}, '(xi_-,eta_pm)': {'status': OP + ' (expected by chirality exchange)'},
 },
 'one_loop_algebra_relations': {
  '[D,P] = -2P': {'statement': 't_Delta = -1/3 (all sectors)', 'status': C2},
  '[D,Q] = -Q': {'statement': 't_Delta = -1/3 (all sectors), independent determination', 'status': C2},
  '[Q,P] = 0': {'statement': 't_P sector-uniform; Ob2(xi,P_b) exact', 'status': C2},
  '[Q,Q] = P': {'statement': 't_xi+ + t_xi- - t_P = 27/4 + 2 log(3/2); t_xi sector-uniform', 'status': C2},
  '{Q,S} = D + M + R': {'statement': 't_M = t_Delta, t_R = t_Delta + 1/3, t_xi+ + t_eta+ - t_Delta = 27/4 + 2 log(3/2); sector uniformity of t_Delta, t_M, t_R, t_eta', 'status': C2},
  '[K,P] = D + M': {'status': OP}, '{S,S} = K': {'status': OP},
 },
 'invariant_content': {
  'scheme_independent': ['t_Delta = theta = -1/3 (fixes D_1 on operators)', 't_M - t_Delta = 0', 't_R - t_Delta = 1/3',
                         't_xi+ + t_xi- - t_P = 27/4 + 2 log(3/2)', 't_xi+ + t_eta+ - t_Delta = 27/4 + 2 log(3/2)',
                         'sector uniformity of t_P, t_xi, t_eta, t_M, t_R', 'pole structure (9/4) u eps^-1 int(A,A), no poles at arity >= 1', 's lc = 0, l_xi = 0 on two legs'],
  'invariance_group': 'inner direction t_z -> t_z + c lambda_z (lambda = -1, +1, -2, +2 on Sigma, Sigma*, P, K; 0 on Delta, so(V), so(W_perp)); outer direction G_B: t_xi+- -> t_xi+- +- a, t_eta+- -> t_eta+- -+ a',
  'scheme_artefacts': ['individual t_xi+-, t_eta+-, t_P, t_K, theta^K_d, kappa^K_j, kappa_V, kappa_lam', 'particular parts with log(3/2) (L - 79/12, 3L - 3/2, -3/2 + 3L, ...)', 'every representative F1^(k)', 'the split intrinsic/position part (density convention)'],
  'translation_covariant_scheme': 't_P = t_Delta = -1/3 => t_xi+ + t_xi- = 77/12 + 2L, t_eta+ = t_xi-',
 },
 'former_shortcuts': {
  'theta = -1/3': 'now computed from the construction twice ((Delta,P_b), (Delta,xi_pm)); operator-level derivation Thm c9-kappa is a cross-check (X-theta agrees)',
  'scaling identity Ob1(Delta) = 4 lc0 + composite': 'Thm at-second (v), a theorem of the construction; cross-check X-scal holds placement class by placement class (special-edge u[-4,16,-4], line u[1/3,-1/3,-2/3])',
  'K = P^dagger': 'cross-check only, not used, not tested on the action side',
  'twist': 'cross-check only (Rem at-twist)',
  'S from Q by position decomposition': 'F1(eta) computed from its own defect cycles; its position part equals x.F1(rho(e)eta) as a computed gate (G-pos eta)',
  'regulator covariance': 'Ob1(R) = 0 and the P_b/M_omega statements are computed (direct cycles), consistent with Thm at-second',
 },
 'open': ['(K_b,P_b\'), (eta,eta\'), (xi_-,eta_pm) at two legs', 'three- and four-leg classes: S1-cubic, Hypothesis at-H (G-H), GC13 along Delta and P_b and analogues, cubic couplings (Prop at-couplings), the (A_V)^5 sector along P_b',
          'the four closed directions G^K_j and the SU(4)-15 free classes (need three legs)', 'G-log finite-eps extrapolation; G-T numerical finite-T run', 'arity-three formula of s (Thm at-eq2 (iii)) written by the forcing argument only',
          'colour beyond su(2) stripped words', 'Thm at-cubic numerically untested'],
}

opA = J('operator/engine_a/results.json')['stages']
operator = {
 'file': 'step2/strict/results/operator_final.json', 'date': '2026-11-04', 'units': U, 'two_component_conventions': WS,
 'sources': ['step2/strict/operator/theory/operator_theory.tex (revised, 27 pp)', 'step2/strict/operator/engine_a (86/86 gates)', 'step2/strict/operator/engine_b (gates_revised.json: 70 pass, 8 withdrawn, 4 open)', 'step2/strict/operator/qa', 'step2/strict/revise.nosync/onebody_F0.json'],
 'master_formula': {'statement': 'hbar L_1(O)^(1)(x) = sum_{a<b} (+-) Upsilon_x(X_a,X_b) - {hbar F_1^(1)(x), O};  hbar L_1(O)^(0) = sum_{a<b} (+-) Upsilon_0(X_a,X_b) - {hbar F_1^(0), O};  Upsilon_x := [eps^0] D_x (defect cycles), Upsilon_0 := [eps^0] (special-edge cycles)',
                    'revision': 'the finite special vertex F_x and the vacuum-pole term of eq op-master0 are cancelled by the bridges (special edge on a tree adjacent to a divergent cycle): bridge + partner = (e^{-eps H} - 1)(cycle - Sing cycle) -> 0 (QA-F1)',
                    'evidence': 'engine A: pair form with [eps^0] D on 646/646 D>=2 rows, fails on all 68 F_x != 0 single-output rows when F_x is added; engine B same 13 rows with F_x != 0',
                    'status': DC + '; ' + C2},
 'F0_scalar_bilinears': {
  'operator': 'O = int omega (iota_w1 A, iota_w2 A) = (1/16) int omega w1^{ij} w2^{kl} (A_ij, A_kl)',
  'B_functionals_two_component': {'B1 = {S^(1)(xi),O}': '(1/4) int omega [(w1^{kl} A_kl, w2_{ji} xi^{alpha i} lambda_alpha^j) + (1<->2)]',
                                  'B2': '(1/4) w1^{ij} w2_ij int omega (A_ji, xi^{alpha i} lambda_alpha^j)',
                                  'B3': '(1/4) w1^{ij} w2_ij int omega i (lambdabar_{alphadot i}, A_m sigmabar^{m alphadot alpha} xi_alpha^i)',
                                  'B4': '(1/4) w1^{ij} w2_ij int omega (c, xi_alpha^i lambda*^alpha_i)'},
  'Q1': {'representative': 'Q_1(O) = (79/12 - 3L - t_xi) B1 + 4L B2 + 4L B3 + 0 B4 (both chiralities, psi = 0)', 'two_body': '-2L B1 + 4L B2 + 4L B3', 'one_body': '-(L - 79/12 + t_xi) B1',
         'arity_zero': 'hbar L_1(O)^(0) = -(2 + 4L) u (w1,w2) int omega (dc, A_V)', 'c_Q': '79/12 - 3L - t_xi', 'status': C2 + ' (two-body, both engines); one-body ' + DC + ' from action numbers'},
  'S1': {'representative': 'S_1(O) = (79/12 - 3L - t_eta) {S^(1)(eta),O} + x_O^mu (4L B2 + 4L B3)(rho(e_mu) eta); intrinsic part 0 on F0', 'c_S': '79/12 - 3L - t_eta', 'status': DC + ' (one-body from action F1(eta)); two-body position part = Thm op-b, intrinsic zero ' + C2},
  'D1': {'representative': 'D_1(O) = -4 t_Delta - 4/3 (traceless), -4 t_Delta - 148/3 (trace); sum Upsilon^0_Delta = 40/3, -104/3; one-body -44/3 - 4 t_Delta', 'value': '0 (traceless), -48 u = 3 g^2 N/(2 pi^2) (trace) at t_Delta = -1/3', 'status': C2},
  'P1': {'representative': 'P_1(O) = -(1/3 + t_P) u d_b O', 'status': C2},
  'M1_R1': {'representative': 'Upsilon_R = 0, Upsilon^0_M = 0 on F0; R_1 = -{hbar F_1^(1)(R),O} = 0 (t_R = 0)', 'status': C2},
  'K1': {'representative': 'K_1(O) = v(x_O)^mu P_1-part - (b,x_O) D^0_1 + M^0_1 + K^0_1 (Thm op-b); one-body -{F1(K_b)|_2,O} = (44/3) u int omega (b,x)(iota_w1 A, iota_w2 A) - theta^K_perp u {S^(1)(K_b)_perp,O} - sum_j kappa^K_j {G^K_j,O}; intrinsic K^0_1 = 0 on F0',
         'status': DC + ' (one-body); G^K_j action on A_perp letters and t_K open'},
  'invariants': {
    'scheme_independent': 'D_1 = 0 (traceless), -48 u (trace, Konishi 2 gamma_K, gamma_K = 3 g^2 N/(4 pi^2))',
    'free_direction_invariant_at_psi0': {'c_Q(xi+) + c_Q(xi-) - p_1': '27/4 - 8 log(3/2)', 'c_Q(xi+) + c_S(eta+) - D_1(traceless)/4': '27/4 - 8 log(3/2)',
        'note': 'independent of t_Delta, t_P, t_xi, t_eta, inner and outer directions (action relations); NOT invariant under operator finite parts psi\' (coboundaries); the operator relations (R2), (R3) on F0 need the weight-5/2 level components and the arity-two kernels: open'},
    'scheme_artefacts': ['c_Q alone, c_S alone', 'every log(3/2) coefficient (a1, a2, a3, -(2+4L))', 'P_1 = -(1/3 + t_P)', 'L_1^(0)(O) (removed by psi(O) = -(1+2L) u (w1,w2) int omega (A_V,A_V))', 'M_1, R_1 (coboundaries)', 'individual pair kernels'],
  },
 },
 'families_spectra': {
  'F1 su(2) tr(Z^{L-M} X^M)': {'D1': '-8 u sum_l (1 - P_{l,l+1}) = (g^2 N/(4 pi^2)) E_Heisenberg, L <= 8', 'status': C2},
  'F2 su(1|1) tr(Z^{L-M} psi^M)': {'D1': '-8 u sum_l (1 - Pi_{l,l+1}) exact matrix identity on 22 sectors, L <= 6', 'status': C2},
  'F3 twist two tr(D_n^k Z D_n^{S\'-k} Z)': {'D1': 'Spec = {-32 h(j) u = (g^2 N/pi^2) h(j): j even <= S\'}, S\' <= 6: {0}, {0}, {-48,0}, {-48,0}, {-200/3,-48,0}, {-200/3,-48,0}, {-392/5,-200/3,-48,0}', 'status': C2},
  'F4 fermion sl(2)': {'D1': 'Spec = {-32 h(j+1) u: j odd <= S\'}, S\' <= 5', 'status': C2},
  'F5 gaugino pair + Z\'F_11': {'D1': 'S\'=0: [[-32,8],[64,-16]] eigenvalues 0, -48; S\'=1: 0, -48 (3x); S\'=2: 0, -48 (4x), -200/3; S\'=3: 0, -48 (4x), -200/3 (3x)', 'note': 'Z\' = iota_w\' A, w\' = e4 - i e5 (A^{34}); off-diagonal entries depend on the normalisation of Z\'F_11', 'status': C1},
  'F6 gluon sl(2)': {'D1': 'S\'=0,1: -48; S\'=2: -200/3, -48', 'status': C1},
  'F7 F_11 psi mixed': {'D1': 'S\'=0: -48; 1: -48, -48; 2: -200/3, -48, -48; 3: -200/3 (2x), -48 (2x); 4: -392/5, -200/3 (2x), -48 (2x)', 'status': C1},
  'P_1 on F3-F7': {'statement': 'sum Upsilon_{P_n} = -1/3 d_n; P_1 = -(1/3 + t_P) u d_n', 'status': C2 + ' on F3, F4; ' + C1 + ' on F5-F7'},
  '(R7) along P_n': {'statement': 'D_1^(S\'+1) d_n = d_n D_1^(S\') exact matrix identity on F3 (0..6), F4 (1..5), F5 (0..3), F7 (0..4)', 'status': C2 + ' (referee recomputation)'},
  'scheme_independent': 'all spectra above; rational (log 2, log 3 cancel in the pair sums)',
 },
 'F8_tables': {'pair_form': 'Upsilon_Delta = 4 lc + C_Delta, Upsilon_P = C_P, Upsilon_M = C_M, Upsilon_R = 0 pair by pair (Conj op-pair)',
               'evidence': 'engine A 4634/4634 rows D<2 and 646/646 rows D>=2 (n <= 2); engine B every complete row n <= 1; 11 engine-B rows incomplete (no verdict); referee 5280/5280 from engine A raw numbers',
               'status': CJ + ' (proved on (A_perp,A_perp)); ' + C2 + ' on the tables', 'G_noLog': 'structural in both routes', 'G_pole_line_rows': OP},
 'twist': {'class_level_check': 'rank-one xi = u (x) w, w_a in the three-plane of w: Q_1(O) = 0 identically on O = int omega (iota_w1 A, iota_w2 A) (untwisted tr gamma^(I gamma^J)), arity-two (xi,xi) has no two-letter output: n = 2 case of [BK Sec 3.4.1, (3.77)], independent of all free directions, psi and the normalisation of Phi',
           'status_check': DC + ' from computed Q_1 on F0', 'transport_map_Phi': OP + ': BK normalisation of b, beta, gamma, d c and the field map (3.51) in no source; (xi,xi) operator kernels for k = 2 entries (bc, bb) not enumerated',
           'class_comparison_BK_3.7_3.59-3.70': OP, 'note': 'Part Q pair entries carry colour-order-dependent log(3/2) terms; individual entries are not comparable'},
 'relations': {'R1': 'computed in Part Q (both routes) on F0, G-R1', 'R2 [Q,Q]=P': OP + ' two-sided (action-side version computed)', 'R3 {Q,S}=D+M+R': OP + ' two-sided (action-side version computed)', 'R4 twist': OP,
               'R5 [K,P]': OP, 'R6 {S,S}=K': OP, 'R7 [delta(Delta),rho_0(x)]=0': C2 + ' one-sided along P_n (F3-F7), su(2) (F1); along Q, S, K open', 'R8 Whitehead': DC},
 'K_equals_P_dagger': {'status': XC + '; untested (G-KP, G-KPbdry open; hypothesis H3)'},
 'former_shortcuts': {'4 lc / scaling': 'Thm op-ward (i), Thm op-a: theorems; G-D1sum, G-D1pair computed identities', 'theta = -1/3': 'action-side computed input', 'S from Q': 'Thm op-b, S_1 from own defect cycles + action F1(eta) one-body', 'K = P^dagger': XC, 'twist': XC, 'regulator covariance': 'Lemma op-comm, Thm op-ward (ii)'},
 'open': ['G-Cid, G-typeI (locality theorem untested), G-logcov, full G-pole on line rows', 'Q_1, S_1, M_1, K_1, generic-b P_1 on F1-F7 at generic kinematics; route B on F5-F7', 'arity-two operator kernels; two-sided (R2)-(R6)',
          'G-KP, G-KPbdry; G-KT in route A', 'twist transport map Phi and class comparison', 'a priori overall sign of defect cycles (gate-fixed)', 'engine A n=3 tables and F5 at S\'=4 (resumable)', 'manuscript propagation of QA-F1 (Def c11-defect, Thm c11-classes, eq c11-pole, eq c10-finitesp): author action'],
}
os.makedirs(os.path.join(ROOT, 'results'), exist_ok=True)
json.dump(action, open(os.path.join(ROOT, 'results/action_final.json'), 'w'), indent=1)
json.dump(operator, open(os.path.join(ROOT, 'results/operator_final.json'), 'w'), indent=1)
print('written')
