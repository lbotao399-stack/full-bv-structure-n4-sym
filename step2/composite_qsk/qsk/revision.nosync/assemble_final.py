"""Revision (Part Q): assemble step2/composite_qsk/results/qsk_final.json from the revised theory, the rerun of engine A,
engine B (first pass + revision reruns), the referee comparison and the revision checks."""
import os, json, time
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
CQ = os.path.join(ROOT, 'step2/composite_qsk')
QSK = os.path.join(CQ, 'qsk')
RV = os.path.join(QSK, 'revision.nosync/out')


def load(p, default=None):
    try:
        return json.load(open(p))
    except Exception as e:
        return default if default is not None else {'unavailable': str(e)}


A = load(os.path.join(QSK, 'engine_a/results.json'))
AG = load(os.path.join(QSK, 'engine_a/gates.json'))
B = load(os.path.join(QSK, 'engine_b/results.json'))
BG = load(os.path.join(QSK, 'engine_b/gates.json'))
QA = load(os.path.join(QSK, 'qa/compare.json'))
QAS = load(os.path.join(QSK, 'qa/out.nosync/ab_strict_summary.json'))
n1 = load(os.path.join(RV, 'ab_compare_n1_lamstar.json'), [])
n1old = load(os.path.join(RV, 'ab_compare_n1_lamstar_oldB.json'), [])
n0 = load(os.path.join(RV, 'ab_compare_n0_lamstar.json'), [])
bl = load(os.path.join(RV, 'b_lamstar_n0.json'), [])
gkt = load(os.path.join(RV, 'gkt_Kb.json'), {})
reg = load(os.path.join(RV, 'tq2_regression_A.json'), {})
DISP = load(os.path.join(QSK, 'disposition.json'))
T = A.get('tasks', {})
L = 'log(3/2)'


def cnt(rows):
    return {'rows': len(rows), 'agree': sum(1 for r in rows if r.get('agree'))}


gates = []
for g in AG.get('gates', []):
    gates.append({'gate': g['gate'], 'route': 'A', 'value': g.get('value'), 'target': g.get('target'), 'pass': g.get('pass'), 'note': g.get('note', '')})
for g in BG.get('gates', []):
    if str(g.get('id', '')).startswith('REV '):
        continue
    gates.append({'gate': g.get('id'), 'route': 'B', 'value': g.get('value'), 'target': g.get('target'),
                  'pass': g.get('pass_or_fail') == 'pass', 'note': g.get('note', '')})
gates.append({'gate': 'G-AB TQ2 n = 0, leg sets <= 2 outputs, all 15 field/ghost pairs, strict colour words (referee)', 'route': 'A vs B',
              'value': QAS.get('ab_compare_all.json'), 'target': 'all rows agree', 'pass': (QAS.get('ab_compare_all.json') or {}).get('strict_agree') == 85})
gates.append({'gate': 'G-AB TQ2 n = 0, rows with a lam* output, route B with the T5 lam*-output placement (revision)', 'route': 'A vs B',
              'value': cnt(n0), 'target': 'all agree', 'pass': bool(n0) and cnt(n0)['agree'] == len(n0)})
gates.append({'gate': 'G-AB TQ2 n = 1 (d_0 on X_a), referee run with the first-pass route B', 'route': 'A vs B',
              'value': {'rows': 85, 'agree': 83, 'disagree': '(A_V, A_perp) -> (c, lam*) both orderings: missing T5 placement in route B'},
              'target': 'all agree', 'pass': 'superseded', 'note': 'first-pass route B without the T5 lam*-output placement; superseded by the next gate'})
gates.append({'gate': 'G-AB TQ2 n = 1, rows with a lam* output, route B with the T5 lam*-output placement (revision)', 'route': 'A vs B',
              'value': {'fixed_B': cnt(n1), 'first_pass_B': cnt(n1old)}, 'target': 'all agree with the fixed route B',
              'pass': bool(n1) and cnt(n1)['agree'] == len(n1)})
gates.append({'gate': 'G-AB intrinsic eta kernels (referee)', 'route': 'A vs B', 'value': '4/4', 'target': 'agree', 'pass': True})
gates.append({'gate': 'G-AB intrinsic K_b kernels (referee)', 'route': 'A vs B', 'value': '2/2 (+ two zero rows)', 'target': 'agree', 'pass': True})
gates.append({'gate': 'B n = 0 lam*-output tables unchanged by the new placement (revision)', 'route': 'B',
              'value': '%d/%d' % (sum(1 for r in bl if r.get('new_equals_old_code')), len(bl)), 'target': 'all',
              'pass': bool(bl) and all(r.get('new_equals_old_code') for r in bl)})
if isinstance(gkt, dict) and gkt.get('gate'):
    g = gkt['gate']
    gates.append({'gate': g['id'], 'route': 'B', 'value': {'configs': g.get('n_configs'), 'placements': g.get('placements')},
                  'target': g.get('target'), 'pass': g.get('pass_')})
gates.append({'gate': 'A TQ2 regression after the exact contraction (revision)', 'route': 'A', 'value': reg.get('summary'),
              'target': 'identical rows', 'pass': 'identical' in str(reg.get('summary'))})
indep = [
    {'gate': 'independent (referee): n = 1 cycle -1/(2 eps) - (3/8)|q|^2, T-dependence O(eps/T)', 'pass': True},
    {'gate': 'independent (referee): K_b composite line T d/dT -> -3 b/(16 pi^2) (scalar), 3/(8 pi^2) with fibre factors (Thm c8-K)', 'pass': True},
    {'gate': 'independent (referee): route B code exactly covariant under the reflection (mirror) and the rotations', 'pass': True},
]
for g in indep:
    g.update(route='independent', value=None, target=None)
gates += indep

tq3 = T.get('TQ3', {})
tq3sum = {k: {kk: vv for kk, vv in v.items()} for k, v in tq3.get('summary', {}).items()}

final = {
 'meta': {'title': 'Part Q: one-loop two-body kernels of Q_1 (Sigma), S_1 (Sigma^*), K_1 (K_b); final after revision',
          'date': time.strftime('%Y-%m-%d %H:%M'),
          'theory': 'step2/composite_qsk/qsk/theory/qsk_theory.tex (_build/qsk_theory.pdf, revised 2026-10-30)',
          'engines': {'A': 'step2/composite_qsk/qsk/engine_a (route A, exact simplex integrals; rerun end to end after the revision)',
                      'B': 'step2/composite_qsk/qsk/engine_b (route B, momentum space, Feynman parameters, eps-extrapolation; T5 placement fixed)'},
          'referee': 'step2/composite_qsk/qsk/qa (compare.json, findings.json)', 'disposition': 'step2/composite_qsk/qsk/disposition.json',
          'revision_checks': 'step2/composite_qsk/qsk/revision.nosync',
          'units': 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); stripped engine values in hbar/(16 pi^2) times colour '
                   'words C_w = (1 (x) ad_{Y_w1}..ad_{Y_wr})(-,-)^{-1}; on a colour-singlet pair with two outputs the stripped value is the coefficient in units u',
          'regulator': 'Schwinger cutoff t >= eps on every edge P(eps,T), special edge K_eps - K_0 with t in [0,eps]; T = infinity along Sigma, '
                       'Sigma^* (Prop c8-Tsafe), finite T along K_b (Def c8-K); no dimensional regularisation in any role; no other IR regulator',
          'numbers': 'exact, Q(i) (x) {1, log2, log3, log5, log7}; L := log(3/2)'},
 'status': 'Kernels computed and cross-checked by two routes: E0, TQ0 (Sigma, Sigma^*), TQ2 (n = 0 all pairs; n = 1..3 route A), '
           'TQ3 twist pair entries at arity one (n <= 1), TQ4 scalar bilinears with the cocycle scheme, G-R1, intrinsic parts along '
           'eta and K_b, G-reg, G-KT along K_b. Open: TQ1 (action side), finite special vertices on route-A rows with D >= 2, '
           'Upsilon_xixi, relations R2-R6, the class-level twist comparison, TQ8.',
 'conventions': {
   'arrows': 'xi.: lambda -> +Gamma(xi,lambda) (slot A), A -> +1/2 rho(dA) xi (slot lambda), A^* -> lambda^*, lambda^* -> A^* (pairing); '
             'Lemma qs-arrows as revised (the first print had the opposite overall sign)',
   'classical_component': '{S^(1)(xi), int omega (iota_w1 A, iota_w2 A)} = +F1 (Cor c7-Q)',
   'sign_rule': {'statement': 'Def qs-defect (ii) without free constant',
                 'route_A': 'on top of the raw chain and the pair sign (-1)^{|X_a|+|X_b|}: (-1)^{n(<sigma) + |S1(x)|(n(<s)+[s<sigma])} (-1)^{[sigma = edge at X_b]}, GLOBAL = +1',
                 'route_B': 'D_x = +(hbar/2) sum s_x STr[...], s_x = (-1)^{n(<sigma) + |S1(x)|(n(<s)+[s<sigma])}, times the reordering sign of the odd outputs',
                 'fixed_by': 'G-reg (Delta, P_b), G-pole, graded symmetry X_a <-> X_b, G-R1'},
   'grassmann_bookkeeping': 'an engine value with ordered outputs is the coefficient of theta_1 theta_2 for chains visiting the outputs in that order; '
                            'the functional coefficient in a fixed order is v(Y1,Y2) - v(Y2,Y1) for two odd outputs, v(Y1,Y2) + v(Y2,Y1) otherwise',
   'orientation': 'dictionary orientation (S_+ = undotted for dx^1 dx^2 dx^3 dx^4); engine A internal labels are mirrored (A lam_+ = dictionary S_- (x) W^*); '
                  'numbers identical by the reflection symmetry (verified in route B code); all reported two-component forms in the dictionary orientation',
   'two_component': 'Def c8-twocomp, Rem qs-twocomp: sigma^m = (sigma^1,sigma^2,sigma^3,i), sigmabar^m = (-sigma^1,-sigma^2,-sigma^3,i), '
                    'eps^{12} = eps^{1dot2dot} = +1; Gamma(xi,lambda)|_V^m = i(lambdabar_{alphadot i} sigmabar^{m alphadot alpha} xi_alpha^i + xibar_{alphadot i} sigmabar^{m alphadot alpha} lambda_alpha^i); '
                    'iota_w A = 1/4 w^{ij} A_ij, (w1,w2) = 1/4 w1^{ij} w2_ij'},
 'proved_theory': [
   'Thm qs-classes: hbar L_1(O)^(1)(x) = sum_{a<b} (Koszul) int omega ... Upsilon_x(X_a,X_b) ... - {hbar F_1^(1)(x), O} for every local O, gauge-variant included; '
   'Upsilon_x = [eps^0] D_x + F_x, F_x = -sum_m ((-1)^m/m!)({I, ct1_{x,m}}_{H^m} + {S^(1)(x), ct_m}_{H^m})',
   'pole identity eq qs-pole (with the vacuum-pole x spectator terms read pair by pair) and arity-one log covariance eq qs-logcov',
   'Cor qs-exist: local representative through arity one exists given the QME through arity one along x',
   'Prop qs-comm (revised): [xi.,Q] = [xi.,H] = 0; C_xi = [xi.,Q^GF]: lambda_- -> c (+i(p,Gamma(xi,lambda)|_V)), A^*|_V -> lambda_+ ((p,A^*|_V) xi), lambda^*_+ -> A|_V, c^* -> lambda^*_- '
   '(xi in S_+ (x) W); eta. = sum x^mu xi_mu. + eta^0. + Y_eta; [eta.,H] = 2 xi_mu. d_mu; [v.,Q^GF] = 4(b,x) Q^GF + C_K, C_K A_V = -4(b,A_V)',
   'Cor qs-edges: V_xi = (C_xi g(H) (x) 1)K_0 in lambda (x) A|_V, c (x) lambda^*; composite line 2(P* xi. P* (x) 1)K_0 = -+2 xi e^{-2 eps p^2}/p^2 (real, negative on c (x) lambda_+); '
   '[Q, P* xi. P*] = E(H) C_xi g(H) != 0',
   'Prop qs-ward (Ward identity of the flow along xi); Prop qs-posdec (position decomposition along eta, K_b); Prop qs-bound (output bounds)',
   'Thm qs-simplex: 16 pi^2 a_0 = int_{alpha_sigma smallest} F_0 log(min_{f!=sigma} alpha_f/alpha_sigma), 16 pi^2 a_m = (1/m) int [(min alpha_f)^m - alpha_sigma^m] F_m, no log eps; '
   'T d_T = O(eps) along K_b',
   'Prop qs-freexi (revised): closed quadratic directions span{S^V_2 = T1+T3, S^perp_2 = T2+T4} along Sigma; along Sigma^* {x.(T1+T3), x.(T2+T4) + s0}; no exact ones',
   'Lemma qs-cocycle (new): with L_1^(0) = 0 on {O, F1} the gauge-variant and BRST coefficients of Q_1 on the trace bilinear vanish; R1 at psi = 0 is a3 + 2 kappa - kappa\' = 0',
   'Prop qs-laws, Thm qs-invariant (renormalisation freedom, invariant content); Prop qs-plane (chiral bilinears); Thm qs-twistcheck (what the twist comparison tests, under D in C S^xi + Im Q^xi_0)'],
 'E0': {'route_A': T.get('E0', {}).get('operators'), 'composite_line_route_A': T.get('E0', {}).get('composite_line'),
        'route_B': {k: B.get('E0', {}).get(k) for k in list(B.get('E0', {}).keys())[:12]} if isinstance(B.get('E0'), dict) else B.get('E0'),
        'agreement': 'both routes: [xi.,Q] = [xi.,H] = 0; C_xi components restricted to V as stated; composite line real'},
 'TQ0': {'Sigma': 'Q-closed quadratic directions span{T1+T3, T2+T4} = span{S^V_2(xi), S^perp_2(xi)}; no exact direction (routes A and B)',
         'Sigma_star': 'span{x.(T1+T3), x.(T2+T4) + s0}: s0 is not an independent closed direction (route B)',
         'K_b': 'not computed'},
 'closed_forms': {
   'simplex': '16 pi^2 [eps^0] = int_{alpha_sigma smallest} F_0 log(min_{f != sigma} alpha_f / alpha_sigma) dalpha; 16 pi^2 [eps^-m] = (1/m) int [(min_{f!=sigma} alpha_f)^m - alpha_sigma^m] F_m dalpha; no log eps',
   'n_t=2': '16 pi^2 a_0 = int_0^{1/2} F_0(alpha) log((1-alpha)/alpha), 16 pi^2 a_m = (1/m) int_0^{1/2} ((1-alpha)^m - alpha^m) F_m: Q + Q log 2',
   'constant_integrands': '16 pi^2 [eps^0] = F_0 log(n_t/(n_t-1))/(n_t-1)! for a constant F_0 (leg sets with D = 0, no pole); 16 pi^2 [eps^-1] = F_1/(n_t^2 (n_t-1)(n_t-1)!) '
                          'for a constant F_1 (D = 2); exact for n_t = 2..6 (route B); e.g. n_t = 3: (1/2) log(3/2) and 1/36',
   'cubic_example': 'k_I = 0 T5 cycle on (A_perp, lambda_+) with output A_perp at the special edge: 16 pi^2 x value = (-1/(2 eps) - (3/8)|q|^2) x factor; eps^0/eps^-1 = (3/4)|q|^2 (both routes, independent check)',
   'D0_legsets': 'chain by chain rational multiples of log(3/2)/2 (n_t = 3) or log(4/3)/6 (n_t = 4) (route B)',
   'n=2_line': 'Prop c11-k1 weights W_0(beta) = (1-beta) log((1-beta)/(2 beta)) - (1 - 3 beta) on beta in [0,1/3]: Q + Q log2 + Q log3'},
 'TQ4_scalar_bilinears': {
   'functionals': {'O': 'int omega (iota_w1 A, iota_w2 A)', 'F1': 'int omega [(iota_w1 A, (rho(w2) xi, lambda)) + (1<->2)] = {S^(1)(xi), O}',
                   'F2': '(w1,w2) int omega (rho(A|_perp) xi, lambda)', 'F3': '(w1,w2) int omega (rho(A|_V) xi, lambda)', 'F4': '(w1,w2) int omega (c, (xi, lambda^*))'},
   'psi_0_two_body': {'a1 (F1)': '-2 L', 'a2 (F2)': '4 L', 'a3 (F3)': '4 L', 'a4 (F4)': '0', 'poles': 'none',
                      'chiralities': 'identical for xi in S_+ (x) W and S_- (x) W^* (both routes)', 'route_A': T.get('TQ4', {}).get('xi_in_S+W'),
                      'route_B': B.get('TQ4', {}).get('decomposition_eps0')},
   'c_Q': {'formula': 'Q_1(O)|_traceless = c_Q {S^(1)(xi), O},  c_Q = -2 L - theta_perp - p_perp  (units u)',
           'two_body': '-2 L', 'one_body': '-theta_perp - p_perp: p_perp from the particular solution of eq qs-arity1xi (TQ1, open), theta_perp the free direction',
           'scheme_family': {'outer direction G_B (a)': 'c_Q -> c_Q - a along S_+ (x) W, + a along S_- (x) W^*',
                             "psi' = (a/2) {S^(1)(Delta), -}": 'c_Q -> c_Q - a/2',
                             'field countings y_f': 'c_Q -> c_Q +- (y_lambda - y_A) (Prop qs-laws (a))',
                             'theta^xi_perp': 'c_Q -> c_Q - theta_perp',
                             'cutoff shape': 'a combination of the above (Prop qs-laws (e))'}},
   'cocycle_scheme': {'L1_at_psi0': 'L_1^(0)(O) = -(2 + 4 L) u (w1,w2) int omega (dc, A_V) (no c d.A_V part); two-field part of L_1^(0)(F1) = 2 u s F3',
                      'psi(O)': '-(1 + 2 L) u (w1,w2) int omega (A_V, A_V)', 'psi(F1)': '-2 u F3 (two-field sector)',
                      'change': '[{S^(1)(xi),-}, psi](O) = (2 kappa - kappa\') F3 = (-2 - 4 L + 2) F3 = -4 L F3',
                      'two_body': {'a1': '-2 L', 'a2': '4 L', 'a3': '0', 'a4': '0'},
                      'statement': 'Q_1(O)|_two-body = -2 L F1 + 4 L F2 on the trace bilinear in a scheme with L_1^(0) = 0 on {O, F1}; the gauge-variant term drops out (Lemma qs-cocycle)',
                      'route_A': T.get('TQ4', {}).get('cocycle_scheme')},
   'two_component': T.get('two_component')},
 'G_R1': {'relation': 's Q_1(O) + L_1^(0)({S^(1)(xi),O}) + {S^(1)(xi), L_1^(0)(O)} = 0 (Cor c8-homotopy (i), psi = 0), two-field sector (c, lambda_-)',
          'terms': 'A = s(a3 F3) = -2 L k, B = L_1^(0)(F1) = -k, C = {S^(1)(xi), L_1^(0)(O)} = (1 + 2 L) k; sum 0, rational and log parts separately',
          'route_A': T.get('G-R1')},
 'TQ3_twist': {
   'setting': 'xi = u (x) z rank one (weight vector), twisted letters basis-free from Cartan weights; Lemma: Q^xi_1 = L_1^(0) + Q_1 - (1/2) L_1^(2)(xi,xi); entries at arity one, psi = 0',
   'normalisation': 'engine weight-vector normalisation (route A); BK normalisation needs the transport map Phi of Thm qs-twistcheck (open)',
   'reduction': 'modulo Im Q^xi_0 only: no antiholomorphic derivative occurs; the derivative on a dc output (holomorphic Dirac component) is dropped, exact up to antifield terms; total derivatives kept',
   'entries': {
     'betagamma [BK (3.66)] (lambda_+, A_perp) -> lambda_- lambda_-': '0 identically for every xi (both routes); c_K = 0 at psi = 0; scheme dependent: c_K -> c_K - kappa_B y',
     'betabeta [BK (3.67)] (lambda_+, lambda_+) -> A_perp lambda_-': 'per unit letters (2/3) d^alphadot gamma^K dc_alphadot + (1/3 - 2 L) (Dirac component, exact); K the complement of I, J (eps_IJK); both colour orders equal; route B: same eps_IJK structure, A code = B code at four momenta',
     'bgamma [BK (3.68)] (A_V with the self-dual u\'u\' derivative, A_perp) -> A_perp lambda_-': 'colour order (gamma, dc): 64/3 - 128 L; order (dc, gamma): 32 + 64 L (coefficient of the derivative on gamma); Dirac components -64/3 - 128 L, -32/3 + 64 L; unequal orders, log(3/2) present',
     'beta b -> beta dc [BK (3.69)]': '-(64/3) c_beta in one colour order, 0 in the other',
     'beta b -> gamma gamma [BK (3.69)]': '(8/3)(c1h1 c2h0 - c1h0 c2h1) on (gamma^J, gamma^K): eps structure, both orders equal',
     'bc [BK (3.64)], b = lambda^* representative u_alpha lambda^{*alpha}_4': '0 (both routes)',
     'gamma gamma (twisted outputs)': '0; G-plane: two-body part on the chiral bilinears exactly 0 (both routes)',
     'bb [BK (3.65)], bc with b = self-dual F': 'not computed (arity two Upsilon_xixi; open)'},
   'route_A_summary': tq3sum,
   'comparison_with_BK': 'BK (3.64)-(3.69) have the form (1/2) f_ACD f_BCE d_alphadot X^D d^alphadot Y^E with equal colour words. At the level of individual '
                         'entries (psi = 0, this regularisation) betagamma and bc(lambda^*) vanish, betabeta and beta b -> gamma gamma have the BK structure with equal words, '
                         'bgamma and beta b -> beta dc have unequal words and (bgamma) log(3/2). Individual entries are scheme dependent (Rem qs-twistfreedom: F1 theta, F3 psi, field rescalings), '
                         'so this is not a disagreement; the class-level test (Thm qs-twistcheck (d): relations [BK, (3.59)-(3.63)] through the tree-level psl(3|3) commutant, '
                         'one-body terms, R_1/S_1 terms) is not done, and Thm qs-twistcheck holds only under the hypothesis D in C S^xi + Im Q^xi_0 (not derived)'},
 'TQ5_eta': {'route_A': T.get('TQ5'), 'route_B': 'intrinsic parts on 5 pairs, agree with route A on 4/4 compared rows; (A_perp, lambda_+) -> A_perp A_perp exceeded the memory cap',
             'S_1_on_bilinears': 'S_1(O) = Q_1 along the Killing spinor rho(x_O) eta (position part); Upsilon^0_eta(A_perp, A_perp) has the single leg set lambda, zero on the colour singlet'},
 'TQ6_Kb': {'route_A': T.get('TQ6'), 'route_B': B.get('TQ5_TQ6'),
            'K_1_on_bilinears': 'K_1(O) = v(x)^mu P_1(O)|_{e_mu} - (b,x) D^0_1(O) + M^0_1(O) (Upsilon^0_K = 0 on the colour-singlet bilinear); in the translation-covariant scheme at theta = -1/3 and with M_1 ~ 0: K_1 ~ -(b,x) D_1',
            'log_T': '(3/8 pi^2) hbar log T x zero-mode tree of the observable at scale T: not part of K_1; O(eps) in the defect graphs',
            'G_KT': gkt.get('gate') if isinstance(gkt, dict) else None},
 'TQ7_regression': T.get('TQ7'),
 'TQ2_tables': {'route_A': {'xi': T.get('TQ2', {}).get('xi'), 'tables': T.get('TQ2', {}).get('tables'), 'orderings': T.get('TQ2', {}).get('orderings'),
                            'finite_special_vertices': 'not included in route A: rows with superficial degree D >= 2 are eps^0 parts of the defect cycles only',
                            'files': 'step2/composite_qsk/qsk/engine_a/out.nosync/tq2_n0.json .. tq2_n3.json', 'regression_after_exact_contraction': reg.get('summary')},
                'route_B': {'n=0': B.get('TQ2', {}).get('summary'), 'finite_special_vertices': 'F_xi = |p|^2 {S^(1)(xi), ct_1(lambda lambda)|_A} on gaugino pairs with a single gaugino output',
                            'files': 'step2/composite_qsk/qsk/engine_b/out.nosync/tq2/ (n = 0); revision recomputation of the lam*-output leg sets in revision.nosync/out/tq2_b_fixed/'},
                'AB_agreement': {'n=0 (referee)': '85/85 strict', 'n=0 lam* rows, fixed B': cnt(n0), 'n=1 (referee, first-pass B)': '83/85',
                                 'n=1 lam* rows, fixed B': cnt(n1), 'n=1 lam* rows, first-pass B': cnt(n1old)}},
 'renormalisation_freedom': {
   'one_body': ['theta^xi_V, theta^xi_perp: -{theta_V S^V_2(xi) + theta_perp S^perp_2(xi), X} = -theta times the classical arrow; the cubic part of eq qs-arity1xi fixes theta_V - theta_perp (TQ1, open), '
                'the common theta^xi is the direction theta S^(1)(xi) = G_B^(1)(xi) / R-rotation',
                'field redefinitions F -> F + hbar s phi: -{(s phi)^(1)(x), -}; field countings multiply each arrow f -> f\' by the counting difference',
                'outer direction G_B: Q_1 -> Q_1 -+ a {S^(1)(xi),-} (S_+ (x) W / S_- (x) W^*), S_1 -> S_1 +- a {S^(1)(eta),-}; K_1, D_1, P_1, M_1, R_1 unchanged',
                "psi' = (a/2){S^(1)(Delta),-}: (a/2) lambda_x {S^(1)(x),-}: -a/2 along Sigma, +a/2 along Sigma^*, +a along K_b"],
   'two_body': ['finite parts psi\' of the operators at the output weight that vanish on O: change -psi\'({S^(1)(x),O}); Konishi entry c_K -> c_K - kappa_B y',
                'cocycle finite part psi with L_1^(0) = 0 on the family: [{S^(1)(xi),-}, psi]; on the trace bilinear a3: 4 L -> 0'],
   'cutoff_shape': 'changes the finite parts of the defect cycles, of the special-edge cycles of L_1^(0) and the finite special vertices; a combination of the above (uniqueness of counterterms up to finite local terms)',
   'twist': 'F0..F5 of q1_targets act as 0 on the cohomology of Q^xi_0; the individual entries are scheme dependent, [Q^xi_1], its kernel, image and rank are not',
   'invariant': ['representation of psl(4|4) on the cohomology of Lambda^(0) up to isomorphism and e^{-a hbar B}', '[Q^xi_1] on H(Q^xi_0)', 'spectrum of D_1 in fixed classical weight',
                 'Cor c8-homotopy and Cor c8-gram relations among the components']},
 'gates': gates,
 'gate_counts': {'route_A': '%s/%s' % (AG.get('n_pass'), AG.get('n_total')), 'route_B': '%s pass, %s fail' % (BG.get('n_pass'), BG.get('n_fail')),
                 'all_listed': '%d pass of %d decided (%d superseded)' % (sum(1 for g in gates if g.get('pass') is True), sum(1 for g in gates if isinstance(g.get('pass'), bool)),
                                                                      sum(1 for g in gates if g.get('pass') == 'superseded')),
                 'failing': [g['gate'] for g in gates if g.get('pass') is False]},
 'open': [
   'TQ1: action side along xi (quadratic Ob^(1)(xi) in (lambda, A_V), (lambda, A_perp), (c, lambda^*), {S^(1)(xi), hbar F_1^(0)}, particular solution p_V, p_perp, cubic part fixing theta_V - theta_perp); needed for c_Q and every one-body kernel of Q_1, S_1, K_1',
   'finite special vertices F_x in route A (rows with D >= 2: single-lambda leg sets of (lambda, lambda), n >= 1 rows, beta b -> A_perp A_perp); route B has them at n = 0',
   'arity-two kernels Upsilon_xixi: bc with b the self-dual F, bb [BK (3.65)]',
   'relations R2, R3 (Konishi relation through D_1), R4 (twisted closure), R5, R6 on the bilinears: need the one-body terms (TQ1) and the cocycle finite part on the descendants beyond the two-field sector',
   'class-level twist comparison: transport map Phi ([ESW, Thm 10.9]) and the hypothesis D in C S^xi + Im Q^xi_0 not derived; relations [BK, (3.59)-(3.63)] not tested',
   'M_1 (rotation component) not evaluated; K_1 ~ -(b,x) D_1 assumes M_1 ~ 0 on the family',
   'TQ8: S_1 = -Q_1^adj + Gram terms (c_S = conj(c_Q) + g_{5/2} - g_2), K = P^dagger: need Gram_1 from Def c8-2pt',
   'TQ0 along K_b; route B TQ2 for n >= 1 (now possible with the fixed T5 placement); the eta intrinsic pair (A_perp, lambda_+) -> A_perp A_perp in route B (memory)',
   'vacuum-pole terms of the finite special vertices (analogue of the G-pole terms) not checked',
   'route A n = 2, 3 tables use one output ordering per leg set'],
 'disposition_summary': [{k: d.get(k) for k in ('id', 'severity', 'decision')} for d in DISP.get('findings', [])],
}
json.dump(final, open(os.path.join(CQ, 'results/qsk_final.json'), 'w'), indent=1, default=str)
print('written', os.path.getsize(os.path.join(CQ, 'results/qsk_final.json')), 'bytes;', final['gate_counts'])
