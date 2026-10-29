"""Assemble step2/composite_qsk/results/composite_final.json and qsk_final.json from the engine, QA and revision outputs."""
import json, os, time
W = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk'
C = os.path.join(W, 'composite'); RV = os.path.join(W, 'revise.nosync'); OUT = os.path.join(W, 'results')
os.makedirs(OUT, exist_ok=True)
L = lambda *p: json.load(open(os.path.join(*p)))
A = L(C, 'engine_a', 'results.json'); GA = L(C, 'engine_a', 'gates.json')
AX = L(C, 'engine_a', 'results_ext.json') if os.path.exists(os.path.join(C, 'engine_a', 'results_ext.json')) else {}
GAX = L(C, 'engine_a', 'gates_ext.json') if os.path.exists(os.path.join(C, 'engine_a', 'gates_ext.json')) else {}
B = L(C, 'engine_b', 'results.json'); GB = L(C, 'engine_b', 'gates.json')
QAc = L(C, 'qa', 'compare.json')
TC0 = L(RV, 'tc0', 'tc0.json'); WC = L(RV, 'word_map', 'word_compare.json'); GL3 = L(RV, 'colour_gl3', 'gl3_check.json')
SPM = L(RV, 'ferm_map', 'spinor_map.json')


def eig_B(d):
    return {k: v.get('D1_at_theta_m13_eigenvalues') for k, v in d.items() if k.isdigit() and isinstance(v, dict) and 'D1_at_theta_m13_eigenvalues' in v}


def pass_B(d):
    return {k: (v.get('P7b_pass'), v.get('P7c_pass')) for k, v in d.items() if k.isdigit() and isinstance(v, dict) and 'P7b_pass' in v}


fam_A = A['families']
twA = {N: {'sumC0_Delta': r['D0']['matrix_by_logbasis']['1'], 'P7c': r['P7c_pass'], 'P7b_n': r['P7b_n_pass'], 'D1_eigenvalues': r['D1']['eigenvalues']}
       for N, r in fam_A['F1_twist_two'].items()}
for k, r in AX.items():
    if k.startswith('F1_twist_two_N'):
        N = k[len('F1_twist_two_N'):]
        twA[N] = {'sumC0_Delta': r['D0']['matrix_by_logbasis']['1'], 'P7c': GAX[k]['P7c_sumC0_eq_-2/3(4+N)'], 'P7b_n': GAX[k]['P7b_n'],
                  'D1_eigenvalues': r['D1']['eigenvalues'], 'run': 'ext_families.py'}
fsA = {N: {'sumC0_Delta': r['M0']['matrix'], 'P7b_n': r['Pn']['pass'], 'D1': r['D1']} for N, r in fam_A['F2_fermion_sl2'].items()}
for k, r in AX.items():
    if k.startswith('F2_fermion_sl2_N'):
        N = k[len('F2_fermion_sl2_N'):]
        fsA[N] = {'sumC0_Delta': r['M0']['matrix'], 'P7b_n': r['Pn']['pass'], 'D1': r['D1'], 'run': 'ext_families.py'}
glA = {N: {'sumC0_Delta_letter': r['M0']['letter_matrix'], 'vv': r['M0']['vv_matrix'], 'P7b_n_letter_zero': r['Pn']['letter_zero'],
           'D1': r['D1'].get('kappa=-1')} for N, r in fam_A['F3_gluon'].items()}
for k, r in AX.items():
    if k.startswith('F3_gluon_N'):
        N = k[len('F3_gluon_N'):]
        glA[N] = {'sumC0_Delta_letter': r['M0']['letter_matrix'], 'vv': r['M0']['vv_matrix'], 'P7b_n_letter_zero': r['Pn']['letter_zero'],
                  'D1': r['D1'].get('kappa=-1'), 'run': 'ext_families.py'}

gatesA = [{'id': g['id'], 'pass': g['pass'], 'statement': g.get('statement', g.get('desc', '')), 'value': g.get('value')} for g in GA]
gatesB = [{'id': k, 'pass': v.get('pass_'), 'value': v.get('value'), 'target': v.get('target')} for k, v in GB.items()]

R = {
 'meta': {'title': 'Composite kernels C^(2)_Delta, C^(2)_{P_b}; D_1, P_1 on the test families (final, after revision)', 'date': time.strftime('%Y-%m-%d %H:%M'),
          'theory': 'composite/theory/composite_theory.tex (revised; 0 LaTeX errors, 0 undefined references)',
          'unit': 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); kernel values hbar/(16 pi^2) x colour word; planar factor (-1)^r/2',
          'regulator': 'heat-kernel proper time: cutoff eps on every edge, P(eps,T) = int_eps^T (Q^GF (x) 1) K_t dt, T = infinity (generic momenta); no dimensional regularisation; no IR regulator other than T',
          'routes': {'A': 'Schwinger/simplex exact, vertex form, injected momentum r', 'B': 'momentum space, kernel form (Euler operator), Feynman parameters, exact + numerical eps check',
                     'QA': 'independent momentum-space Taylor expansion (scalar sector) and own harmonic action'}},
 'definitions': {
   'square_kernel': 'X = ((P*)^2 (x) 1) K_0; component c (x) A_V: 2 i p_m e^{-2 eps |p|^2} |p|^{-4} (two Schwinger segments in [eps, infinity))',
   'composite_kernels': 'K_x = (v_x. (x) 1 - 1 (x) v_x.) X;  K_Delta = 2 (v. (x) 1) X + 4 (eps d_eps + T d_T) X;  K_{P_b} = 2 (d_b (x) 1) X',
   'internal_vertex': 'Delta: int((2E+2) A_V, A*_V) = 2 int x^mu (d_mu A_V, A*_V) + 2 int (A_V, A*_V);  P_b: int (d_b A_V, A*_V)',
   'Z_x': 'Z_x(X_a,X_b;eps) = (hbar/2) d_{K_x} d_{K_eps - K_0}(amputated chains X_a -> X_b) = -(hbar/2) d_{K_x}(chains with one bracket-defect vertex)',
   'koszul_sign': 's(c,sigma) = (-1)^{n(<sigma) + n(<c) + [c<sigma]}, n(<j) = number of odd outputs before edge j',
   'C2': 'C^(2)_x(X_a,X_b) = coefficient of eps^0 of Z_x; C^(2),0_Delta its translation-invariant part; (C^(2)_Delta)^b = 2 C^(2)_{P_b}; [C^(2),0_Delta, d_b] = 2 C^(2)_{P_b}',
   'simplex_form': '16 pi^2 [eps^0] = int_{R_sigma} F_0 log(min_{f != sigma} alpha_f / alpha_sigma) d alpha;  16 pi^2 [eps^-m] = (1/m) int_{R_sigma} [(min alpha_f)^m - alpha_sigma^m] F_m;  no log eps',
   'k1_weights': 'W_0(beta) = (1-beta) log((1-beta)/(2 beta)) - (1 - 3 beta), beta in [0, 1/3]; W_m(beta) = (2/m)[(((1-beta)/2)^{m+1} - beta^{m+1})/(m+1) - beta^m ((1-beta)/2 - beta)]',
   'N': 'N(int omega X_1...X_n) = sum_a (2 + |alpha_a|) int omega X_1...X_n',
   'free_directions': 'hbar F_1^(1)(Delta) = 4uF + F^tr + u sum_d theta_d S_d + u Q(kappa_V int(iota_x A_V, c*) + kappa_lambda int(rho(x) lambda*, lambda*)) + cubic; F^tr = u[-2/3 int(A,A*) + 1/3 int(lambda*,lambda) + 4/3 int(c*,c)]; hbar F_1^(1)(P_b) = u sum_d thetaP_d S_d(P_b) + kappaP terms',
   'two_component': {'letters': 'A_V -> A_m, A_perp -> A_ij (6 of SU(4)), lambda_+ -> lambda_alpha^i, lambda_- -> lambdabar_{alphadot i}; antifields A*_m, A*_ij, lambda*^alpha_i, lambdabar*^{alphadot i}, c*',
                     'pairing': '(lambda, lambda*) = (lambda_alpha^i, lambda*^alpha_i) - (lambdabar_{alphadot i}, lambdabar*^{alphadot i}); (u,w) = (1/4) u^ij w_ij on W_R^perp',
                     'propagators': '<A_m A_n> = 2 (-,-)^{-1} delta_mn, <A^*_m c> = -(-,-)^{-1} i p_m, times e^{-eps|p|^2} |p|^{-2}',
                     'internal_vertex': 'int((2E+2)A_m, A*_m) = 2 int x_n (d_n A_m, A*_m) + 2 int (A_m, A*_m)',
                     'S_sectors_P': 'S_perp(P_b) = (1/4) int (b_m d_m A^ij, A*_ij); S_lambda(P_b) = -int[(b.d lambda_alpha^i, lambda*^alpha_i) - (b.d lambdabar_{alphadot i}, lambdabar*^{alphadot i})]; S_V(P_b) = int (b.d A_n, A*_n) - int (b.d c, c*)',
                     'kappa_directions': 'Q int x_m (A_m, c*), Q int (rho(x) lambda*, lambda*) = -2i Q int x_m (lambda*^alpha_i, sigma^m_{alpha alphadot} lambdabar*^{alphadot i})',
                     'Srednicki': 'A_m = -i g A^a_m T^a, u = -g^2 N/(32 pi^2), -u/3 = g^2 N/(96 pi^2), gamma = (g^2 N/(8 pi^2)) E_H, D_1 = 2 gamma'}},
 'closed_forms': {
   'bilinear': {'index_free': 'C^(2)_Delta(A_perp, A_perp)|_{A_perp A_perp} = -(1/3) u v.O - (4/3) u O; C^(2)_{P_b} = -(1/3) u d_b O; C^(2),0_Delta = -(8/3) u O (traceless and trace)',
                'two_component': 'O = (1/16) int omega w1^ij A_ij w2^kl A_kl: C^(2)_Delta = u[-(1/3)((2x.d+2)A_ij) A_kl - (1/3) A_ij (2x.d+2) A_kl - (4/3) A_ij A_kl]',
                'chain': 'A_perp -[K_eps - K_0]- ([A_perp, c], A*_perp) -[K_Delta]- A_V (A_perp)^2 -[P]- A_perp', 'routes': ['A', 'B', 'QA'], 'fitted_constants': 0},
   'universal_two_site_density': {'statement': 'planar adjacent pair, every ordered pair of A_perp, lambda_+, lambda_- (on shell; A_perp also off shell, all SO(6) components)',
                                   'C2_0_Delta': '-(4/3) u X_a X_b (same word), 0 (exchanged word)', 'C2_Pb': '-(1/6) u d_b (X_a X_b) (same word), 0 (exchanged)',
                                   'two_component': '(A_ij, A_kl), (Z, psi), (psi, Z), (psi, psi): -(4/3) u and -(1/6) u b_m d_m on the same word', 'routes': ['A', 'B']},
   'k1_offshell_example': {'pair': '(A_V, A_W) -> A_W', 'values': 'Q + Q log 2 + Q log 3 (kernels_C2.json); eps^-1 pole = route-B symbolic pole (QA)'},
   'off_shell_tables': {'file': 'composite/engine_a/kernels_C2.json (value = letter - vertex_vertex)', 'basis': '{1, log 2, log 3} at k <= 2',
                        'two_route_agreement': {'bosonic_rows': WC['bosonic_summary'], 'fermionic_rows': WC['fermionic_summary'],
                                                'word_map': 'A colour word (outputs in list order) = B word with reversed indices',
                                                'spinor_map': {'S^T (Mp R_B[M]) S = sb_A[M], M = 0..9': SPM['forms_intertwined'], 'rho on Sigma*': SPM['rho_Sigmastar_intertwined'],
                                                               'c^2': SPM['c2'], "A's lambda_+ (rho(e0e1e2e3) = +1) = B components": SPM['A_lamP_to_B_components'],
                                                               'meaning': 'the two routes orient V oppositely: lambda_+ of A = lambda_- of B; parity-even results unaffected'}}}},
 'vertex_vertex_sign': {'value': -1, 'derivation': 'engine_a/notes.txt sec. 2b: (-1)^[c<sigma] of s(c,sigma) carried by S L_0 S for the special edge at X_b, explicit on interior lines',
                        'evidence': {'bosonic_rows_agree_with_+1': WC['bosonic_summary']['agree_with_vv_sign_+1'], 'fermionic_rows_agree_with_+1': WC['fermionic_summary']['agree_with_vv_sign_+1'],
                                     'A*_V c and A*_W c leg sets': 'vanish identically at -1'},
                        'families': 'vertex-vertex parts vanish on shell (gluon) or do not occur: family results independent of the sign'},
 'families': {
   'F0_bilinears_su2': {'A': {'so6_density': fam_A['F0_scalar_density'], 'su2_L<=8': 'Heisenberg spectra x (-8u)'}, 'B': B['F0_so6_pair_AperpAperp_W0'].get('planar_density'),
                        'D1_bilinears': '0, -48 u'},
   'F1_twist_two': {'A': twA, 'A_generic_beta_N<=2': 'sum C_{P_beta} = -(1/3) u d_beta + (q_i^2 terms, free EOM)', 'B_D1_eigenvalues': eig_B(B['F1_Z_sl2']),
                    'B_P7bc_pass': pass_B(B['F1_Z_sl2']), 'QA': 'S\' <= 4, generic complex null legs, generic b: sum C^0 = -2/3(4+S\'), sum C_{P_b} = -1/3 i b.(q1+q2)',
                    'target': 'D_1 = -32 u h(S), S even <= S\''},
   'F2_su11': {'A': fam_A['F2_su11_spectra']['pass'], 'B': {k: v.get('D1_at_theta_m13_eigenvalues') for k, v in B['F2_su11']['cyclic'].items()},
               'target': 'Table tab:cc-su11 x (-8u): (3,2),(3,3): 6; (4,2),(4,3): 4; (5,2),(5,3): 5 +- sqrt5; (5,4),(5,5): 10; (6,2): 2, 6; (6,3): 2, 6, 8, 8; (6,4): 8^3; (6,5): 8'},
   'F2_fermion_sl2': {'A': fsA, 'B_D1_eigenvalues': eig_B(B['F2_psi_sl2']), 'B_P7bc_pass': pass_B(B['F2_psi_sl2']), 'target': 'D_1 = -32 u h(n+1), n odd'},
   'F2_psi3psi4_ZF11': {'B_D1_eigenvalues': eig_B(B['F2_psiZF']), 'B_P7bc_pass': pass_B(B['F2_psiZF']),
                        'sumC0': 'diag(-2/3 u (4+S\') on psi psi\', -2/3 u (5+S\') on Z F11), no cross-sector mixing', 'routes': ['B']},
   'F3_gluon_sl2': {'A': glA, 'B_D1_eigenvalues': eig_B(B['F3_gluon']), 'B_P7bc_pass': pass_B(B['F3_gluon']), 'target': 'D_1 = -32 u h(n+2), n even; tr(F11 F11): -48 u'},
   'F3_F11_psi': {'B_D1_eigenvalues': eig_B(B['F3_mixed']), 'B_P7bc_pass': pass_B(B['F3_mixed']), 'routes': ['B'], 'target': 'E_H = 6; 6^2; 6^2, 25/3; 6^2, 25/3^2; 6^2, 25/3^2, 49/5'}},
 'D1_P1_with_free_directions': {
   'index_free': 'on cohomology of every family: D_1 = 4 lc^natural - u sum_a (theta_{d(a)} + 1/3) v.X_a;  P_1 = -u sum_a (thetaP_{d(a)} + 1/3) d_b X_a;  kappa_V, kappa_lambda s-exact',
   'two_component': 'v.X = (2 x_m d_m + 2 w_X) X on a letter; Xi_Delta = -(1/3 + theta) u (2 x_m d_m + 2 w_O) for equal theta',
   'arity_two': 'position part of D_1 s-exact (Thm cc-arity2) => theta_perp = theta_lambda = theta_V = -1/3 (Prop cc-theta-sector; families F1, F2 fermion, F3 gluon)',
   'result': 'Xi_Delta = 0, D_1 = 4 lc^natural, Spec D_1 = -8 u E_H = (g^2 N/(4 pi^2)) E_H; P_1 = 0 in the translation-covariant scheme (thetaP = theta)',
   'A': A['D1_P1_symbolic'], 'B': B['D1_P1_free_directions']},
 'TC0_quadratic_directions': {'statement': 'H^{-1}(Q) on Spin(4)xSpin(6)-invariant quadratic functionals of weight 4, position degree <= 1: parity even = span{S_perp, S_lambda, S_V}; Q-exact: Q int(iota_x A_V,c*), Q int(rho(x) lambda*,lambda*); parity odd adds int(lambda*, rho(v1v2v3v4) lambda) = int[(lambda_alpha^i,lambda*^alpha_i) + (lambdabar_{alphadot i}, lambdabar*^{alphadot i})]',
                              'dims_parity_even': {'perp': [TC0['perp / parity_even']['rank_structures'], TC0['perp / parity_even']['dim_closed'], TC0['perp / parity_even']['dim_exact']],
                                                   'lambda': [TC0['lambda / parity_even']['rank_structures'], TC0['lambda / parity_even']['dim_closed'], TC0['lambda / parity_even']['dim_exact']],
                                                   'V_and_c': [TC0['V and c / parity_even']['rank_structures'], TC0['V and c / parity_even']['dim_closed'], TC0['V and c / parity_even']['dim_exact']],
                                                   'order': '[structures, Q-closed, Q-exact]'},
                              'dims_all': {k.split(' / ')[0]: [TC0[k]['rank_structures'], TC0[k]['dim_closed'], TC0[k]['dim_exact']] for k in TC0 if k.endswith('/ all')},
                              'known_closed': TC0['known'], 'script': 'revise.nosync/tc0/tc0.py'},
 'colour_gl3': GL3,
 'gates': {'engine_A': {'pass': sum(1 for g in GA if g['pass']), 'total': len(GA), 'list': gatesA},
           'engine_A_ext': GAX,
           'engine_B': {'pass': sum(1 for g in GB.values() if g.get('pass_') is True), 'fail': sum(1 for g in GB.values() if g.get('pass_') is False),
                        'not_run': sum(1 for g in GB.values() if g.get('pass_') is None), 'list': gatesB},
           'two_route_offshell': {'bosonic': WC['bosonic_summary'], 'fermionic': WC['fermionic_summary']},
           'TC0': True, 'GC12_gl3': GL3['all_pass'],
           'vendored_B_modules_regression': {'fermion_sl2_Sprime3_entries_identical': 40, 'density_entries_identical': 109},
           'not_run': ['GC1', 'GC2', 'GC5c', 'GC6', 'GC7', 'GC8 (beyond the bilinear)', 'GC13', 'TC2 (n <= 3 tables)']},
 'status': {
   'proved': ['Lemma cc-kernels', 'Prop cc-master', 'Prop cc-position', 'Lemma cc-canonical', 'Lemma cc-bracketpos', 'Prop cc-census, cc-bound', 'Thm cc-simplex, Prop cc-genfun, cc-k1',
              'Prop cc-arity1 (i),(ii)', 'Thm cc-arity2, Cor cc-theta, Thm cc-P7, Prop cc-theta-sector, Cor cc-noA (under their hypotheses)'],
   'outline': ['Thm cc-local (two-body locality, no one-body term, pole identity, log covariance): consequences checked indirectly; Cor cc-D1local and Prop cc-routes conditional on it'],
   'computed_exact': ['P7 (b), (c) on F0, F1 (S\' <= 6 A / 4 B), F2 su(1|1) L <= 6, F2 fermion sl(2) S\' <= 5, F2 psi3psi4/ZF11 S\' <= 2 (B), F3 gluon S\' <= 2, F3 (F11,psi) S\' <= 4 (B)',
                      'D_1 spectra = -8 u E_H on all these families', 'TC0', 'two-route off-shell agreement (22 bosonic + 26 fermionic rows)'],
   'open': ['GC13: cubic arity-one equation (needs the cubic composite part, composite loops of I with three legs)', 'GC5c, GC6, GC7, GC8: pole identity, log covariance, locality demonstration, two routes beyond the bilinear',
            'TC2 tables with letter derivatives; plain composite cycles and ct^(1)_x untabulated', 'existence of the arity-two quantisation (hypothesis of Thm cc-arity2)',
            '(F11,psi) and psi3psi4/ZF11 families single-route (B)']},
}
json.dump(R, open(os.path.join(OUT, 'composite_final.json'), 'w'), indent=1, default=str)
Q = {'meta': {'title': 'Two-body kernels of Q_1, S_1, K_1 (QSK part)', 'date': time.strftime('%Y-%m-%d %H:%M')},
     'status': 'failed upstream: the QSK theory step, both QSK engines and the QSK referee returned no output (theory Q null, QSK engines [null, null], QSK referee null); no file exists under step2/composite_qsk/qsk/',
     'available_inputs_unchanged': ['manuscript Ch10: Prop c10-arity1 (cycles), Lemma c10-active1 (active letters), Prop c10-bound1 (bound; along xi in Sigma weight w_a + w_b + 1/2), Rem c10-next',
                                    'step2/literature/q1_targets.(tex|json), bgkwwy, gkw (twist targets)', 'step2/theory/odd_components_design.tex, engine_plan.json, k_transfer.(tex|json)'],
     'results': None, 'gates': None,
     'consequence_from_part_C': {'P_1': 'P_1 = 0 on cohomology in the translation-covariant scheme (Cor cc-theta (iii)); otherwise Prop c9-scheme finite part with a = u(theta - theta_P) changes K_1, the Sigma and Sigma* components by a, -a/2, +a/2 times the classical transformation',
                                 'note': 'the QSK kernels must be computed in the same scheme (theta = -1/3, vertex-vertex sign and s(c,sigma) as in Def cc-sp) to be combined with D_1, P_1 of part C'},
     'next': 'rerun the QSK theory step and two engines; nothing to finalise in this revision'}
json.dump(Q, open(os.path.join(OUT, 'qsk_final.json'), 'w'), indent=1)
print('written', os.listdir(OUT))
