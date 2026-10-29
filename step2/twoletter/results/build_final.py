"""Assemble step2/twoletter/results/kernels_final.json from the revised theory, the engine outputs, the comparison and
the referee checks (R8, revision of 2026-10-29).  Exact coefficients are kept as the strings of the engines."""
import json, os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
TL = os.path.dirname(HERE)
J = lambda *p: json.load(open(os.path.join(TL, *p)))

A = J('engine_a', 'kernels.json')
GA = J('engine_a', 'gates.json')
GB = J('engine_b', 'gates.json')
TS = J('engine_b', 'two_site.json')
CMP = J('qa', 'compare', 'compare.json')
CMPR = J('results', 'compare_rerun', 'cmp_W0_engineA_D2B.json')['summary']
W2 = J('engine_a', 'out.nosync', 'w2_gate.json')
SU = J('engine_a', 'out.nosync', 'su11_gate.json')
CEN = J('theory', 'census.json')
CENB = J('theory', 'census_tableB.json')
COMP = J('theory', 'census_composite.json')
PRED = J('theory', 'predictions.json')
RT = J('qa', 'refute_theory', 'findings.json')
RN = J('qa', 'refute_numbers', 'findings.json')
DISP = J('qa', 'disposition.json')

U = 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2) for su(N), (X,Y) = -(4 hbar/g^2) tr(XY)'

out = {}
out['meta'] = dict(
    title='R8 two-letter kernels: final agreed results (revision of 2026-10-29)',
    theory='step2/twoletter/theory/twoletter_theory.tex (_build/twoletter_theory.pdf, 22 pages)',
    disposition='step2/twoletter/qa/disposition.json (22 findings, all accepted)',
    unit=U,
    stripped_unit='engine table values are in units hbar/(16 pi^2) for the colour element (1 (x) w)(-,-)^{-1}, '
                  'w = ad_{Y_r}...ad_{Y_1}, Y_1 nearest X_a; a planar adjacent pair of a su(N) single trace gives (-1)^r '
                  '(value/2) u with the pair replaced by Y_1...Y_r',
    conventions='Euclidean R^4, weight e^{S/hbar}, P(eps,T) = int_eps^T (Q^GF (x) 1) K_t dt, per-edge cutoff t_e >= eps, '
                'T = infinity at generic momenta (Prop prop:tl-map (iii), Prop c8-Tsafe), no dimensional regularisation, no IR '
                'regulator besides T; manuscript notation (V, W_R^perp, Sigma, rho, Gamma), letters A|_V, A|_perp, lambda_+-, c '
                'and antifields; two-component form by Def c8-twocomp',
    notation_map_engine_strings='engine A closed-form strings use the engine variables s, u (point splitting), sigma = u - s, '
                                'c(alpha), E = number of edges; theory: y_a, y_b, y_ab = y_b - y_a, l_0(alpha), k+1',
    status_legend=dict(proved='theorem/proposition of the theory with proof', computed_2='computed by both engines, '
                       'equal entry by entry', computed_1='computed by engine A and confirmed by the referee engine',
                       derived='derived in the theory from manuscript inputs (Ch9)', conjecture='stated, untested or '
                       'partially tested', open='not computed'))

out['counterterm_map'] = dict(
    decomposition=dict(
        abstract='d_delta I^CT[eps](M) = hbar Sing_eps d_{P_e(eps,T)} M |_{e=fields} = sum_{a<b} (-1)^{n_ab} int omega '
                 'X_1...K_eps(X_a,X_b)...(X_b omitted)...X_n,  n_ab = |X_b|(|X_{a+1}|+...+|X_{b-1}|);  no one-body terms',
        kernel='K_eps(X_a,X_b) = hbar Sing_eps [(d^alpha_x (x) d^beta_y) P_e(eps,T)^{phi_a phi_b}(x,y)]_{y=x}, '
               'P_e = sum_k P (I''_e P)^k = hbar Sing <X_a(x) X_b(y)>_e (Lemma lem:tl-sign)',
        two_component='same sum over pairs; kernel = hbar Sing of the point-split dressed propagator <X_a(x) X_b(y)>_e '
                      'with <A_m A_n> = 2(-,-)^{-1} delta_mn, <A_ij A_kl> = 4(-,-)^{-1} eps_ijkl, <c A*_m> = '
                      '-(-,-)^{-1} i p_m, <lambda_alpha^i lambdabar_{alphadot j}> = -(-,-)^{-1} p_m sigma^m delta^i_j, '
                      'times (e^{-eps p^2} - e^{-T p^2})/p^2',
        status='proved (Thm thm:tl-decomp)'),
    sign_of_odd_pairs='Hessian-of-M supertrace = (-1)^{|X_a|+|X_b|} K (Lemma lem:tl-sign); engine A derives it, engine B '
                      'equals it entrywise (status proved; computed_2)',
    singular_terms='only eps^{-m} (1 <= m <= w_M/2) and log eps; ct_m of weight w_M - 2m, lc of weight w_M; T-independent; '
                   'scaling (v. - 4 eps d/deps) d_delta I^CT(M) = d_delta I^CT(v.M) - 4 lc(M) (Prop prop:tl-map, proved)',
    finite_special_edge_vertex='-sum_m ((-1)^m/m!) ({I, ct_m(M)}_{H^m} + {ct_m(I), M}_{H^m}), {,}_{H^m} from (H^m (x) 1)K_0 '
                               '(eq tl-finitesp; scheme dependent, enters L_1)')

out['closed_forms'] = dict(
    k0=dict(formula='K^(0)(d^alpha phi_a, d^beta phi_b) = hbar Phat (-1)^{|beta|} (d^{alpha+beta} S_eps)(0), '
                    '(d^gamma S_eps)(0) = (-1)^{m+1}/(16 pi^2 m 2^{m-1} eps^m) sum_{pairings of gamma} prod delta, |gamma| = 2(m-1)',
            generating='16 pi^2 S_eps(y_ab) = sum_m eps^{-m} (-1)^{m+1} (|y_ab|^2/4)^{m-1}/m!',
            content='single power eps^{-w_ab/2}; zero for odd w_ab; no log, no finite part', status='proved; computed_2'),
    simplex=dict(formula='16 pi^2 [eps^-m] = (1/m) int_Delta (min_e alpha_e)^m F_m dalpha, 16 pi^2 [log eps] = -int_Delta F_0 '
                         'dalpha, F_m = sum_{m_j >= m} num_j (-Phi)^{m_j-m}/(m_j-m)!, num_j = (Lap^j N)(l_0)/(j! 4^j), '
                         'm_j = j+1-k, l_0 = -sum alpha_j R_j, Phi = sum alpha_j |R_j|^2 - |l_0|^2', status='proved'),
    shape='log eps independent of the cutoff shape chi; eps^{-m} and finite part depend on it (Prop prop:tl-shape)',
    special_edge='no log eps; [eps^0] = -int_{R_s} F_0 log(alpha_s/min_{f != s} alpha_f) (Lemma lem:tl-special); toy: '
                 '-1/(2 eps) - 3|p|^2/8',
    generating_function='F_m^gen = (-1)^q e^{i(y_ab,l_0)} (|y_ab|/(2 sqrt Phi))^q I_q(|y_ab| sqrt Phi), q = m+k-1; '
                        'sum y_a^alpha y_b^beta/(alpha! beta!) [log]K(d^alpha phi_a, d^beta phi_b) = -(hbar/16pi^2) '
                        'e^{i(y_b,R)} N_V(-i d_{y_ab}) int_Delta F_0^gen; eps^{-m}: (hbar/(16 pi^2 m)) e^{i(y_b,R)} '
                        'N_V(-i d_{y_ab}) int_Delta (min alpha)^m F_m^gen; iR = d on the product of the outputs',
    bubble_closed_form=dict(
        formula='16 pi^2 [log eps] = -sum_{n1,n2,j} (i(y_a,R))^{n1} (i(y_b,R))^{n2} |y_ab|^{2j} (|R|^2/4)^j '
                '(n1+j)!(n2+j)!/(n1! n2! (j!)^2 (n1+n2+2j+1)!)  (constant N_V: quartic bubbles)',
        example='(d_mu phi_a, d_nu phi_b) of a quartic bubble -> (hbar/16pi^2) log eps (-1/6 d_mu d_nu - 1/12 delta_mu nu box)'
                '(Psi_1 Psi_2) x (colour, fibre factor)',
        checks='1350 checks against engine A quartic blocks with up to 3 derivatives (qa/compare); engine A n <= 4',
        status='proved; computed_2'),
    chain_numerator_table='edges: A-A degree 0 (2(-,-)^{-1} eta), c-A*_V degree 1 (-i(l+R_j)), lambda-lambda degree 1 '
                          '(i rho(l+R_j)); vertices (A_V)^3, A_V A_perp^2: one leg momentum; others constant (Rem rem:tl-NV)',
    engine_A_strings=A['closed_forms'])

cen = CEN['output_leg_sets_by_pair_weight']
tableA = {}
for pair, d in cen.items():
    w0 = min(d.keys(), key=lambda x: Fraction(x))
    tableA[pair] = dict(w0=w0, log_leg_sets=d[w0]['log'], eps_leg_sets_upper_bound=d[w0]['eps_max_m'])
out['census'] = dict(
    active_letters='A|_V, A|_perp, lambda_+, lambda_-, c, A*|_V (counterterm); A*|_perp, lambda*_+-, c* spectators of the '
                   'counterterm, active at arity one (Lemma lem:tl-active1)',
    edges='A_V-A_V, A_perp-A_perp, lambda_+-lambda_-, c-A*_V',
    degree='D = w_ab - w_out (Lemma lem:tl-degree)',
    bound='eps^{-m}: weight w_ab - 2m, log: weight w_ab; n_out <= w_out + d_ab <= w_ab + d_ab; k <= n_out (Prop prop:tl-bound)',
    composite_bound='C^(2)_Delta local of weight w_ab (x at weight -1), at most linear in x, n_out <= w_ab + 1 + d_ab; composite '
                    'cycles enter only with pre-Q output weight <= w_ab + 3 (Prop prop:tl-compbound)',
    arity_one_bound='K^(x): weight w_ab - lambda_x/2, degree <= deg_x in x, n_out <= w_ab - lambda_x/2 + deg_x + d_ab '
                    '(Prop prop:tl-bound1)',
    table_A_W0=tableA, table_B_counts=CENB, composite_W0=COMP['pairs'],
    status='proved; enumerations reproduced independently (referee: 108/108 cells; compare C3: 316 kernels, 0 outside)')

out['kernels_W0'] = dict(
    description='exact two-letter kernels at w_ab = w^0_ab for all 21 pair types: every log / eps^{-m} leg set, ordered '
                'word, block structure, letter and output component, output-momentum polynomial (engine A basis, '
                'meta.components); equal to engine B entry by entry (389/389 rows; rerun 2026-10-29: %s/%s)' % (
                    CMPR['n_agree'], CMPR['n_rows']),
    basis=A['meta'], tables=A['pair_weight_W0_tables'],
    identically_zero_leg_sets=A['census_leg_sets_with_identically_zero_kernel_at_W0'],
    invariant_forms_bosonic=A['invariant_forms_bosonic_W0'],
    derivative_tables='W_0+1 (all pairs) and W_0+2 (AV|AV, AV|AW, AW|AW), <= 2 outputs: engine_b/kernels_deriv.json, equal to '
                      'engine A on demand (164/164 words); higher: engine A on demand (ea_kern.kernel), referee 560 + 14144 '
                      'entries, 0 mismatches',
    status='computed_2')

out['planar_densities'] = dict(
    so6_scalar=dict(abstract=A['planar_scalar_two_site']['statement'], lc_natural=A['planar_scalar_two_site']['lc_natural'],
                    two_component=A['planar_scalar_two_site']['two_component'], status='computed_2'),
    su2_sector=dict(statement='Tr(Z^{L-M} X^M): lc^natural = -2u sum_l (1 - P_{l,l+1}); cyclic spectra L <= 8 = Heisenberg '
                              '(predictions.json P6)', engine_A=A['su2_sector'], status='computed_2'),
    su11_sector=dict(statement='{Z, psi}: K(Z,psi) -> u(Z psi) + 2u(psi Z), K(psi,psi)|_{psi psi} = 0, K(Z,Z)|_{ZZ} = 2u; '
                               'lc^natural = -2u(1 - Pi) (graded) on (Z,psi), (Z,Z), (psi,psi)',
                     two_component='Z = A_{I=5} + i A_{I=6}, psi = lambda_alpha^i of so(6) weights (1/2,1/2,1/2)',
                     engine_A=SU, status='computed_1 (engine A gate G8f; referee for (Z,psi))'),
    twist_two=dict(operator=A['twist_two_operator'], engine_A=A['twist_two'], engine_B=TS['twist_two'],
                   statement='eig lc^natural = -8 h(S) u (S even <= N <= 8): 0, -12, -50/3, -98/5, -761/35; matrix level '
                             'lc^natural = -2(H\'_12 + H\'_21)u in the basis D^k Z/k! (qa/compare P6b)',
                   status='computed_2'),
    konishi=dict(lc_natural='-12u', gamma='-24u = 3 g^2 N/(4 pi^2)', D1='-48u (theta_perp = -1/3)', status='computed_2'),
    bilinears=dict(lc='4u O (traceless), -8u O (trace)', lc_natural='0, -12u', status='computed_2'))

D1 = A['D1']
out['D1'] = dict(
    formula='D_1(O) = 4 lc^natural(O) - {F\'_Delta, O} + C_Delta(O) (Thm thm:tl-D1), lc^natural = lc - u{F,-}, '
            'F = 2 int(A,A*) - 4 int(lambda*,lambda) - int(c*,c)',
    F_prime_Delta_quadratic='u[-2/3 int(A,A*) + 1/3 int(lambda*,lambda) + 4/3 int(c*,c)] + u(theta_perp S_perp + theta_lambda '
                            'S_lambda + theta_V S_V) + u Q(kappa_V int(iota_x A_V, c*) + kappa_lambda int(rho(x) lambda*, '
                            'lambda*)) + ...;  S_perp = int(v.A_perp, A*_perp), S_lambda = -int(v.lambda, lambda*), S_V = '
                            'int(v.A_V, A*_V) - int(v.c, c*), S^(1)(Delta) = S_perp + S_lambda + S_V',
    two_body=dict(kernel='4 K^log(X_a,X_b) + C^(2)_Delta(X_a,X_b)',
                  K_log='kernels_W0 (computed_2) and closed forms',
                  C2_Delta_planar_scalar_sector='-u v.(X_a X_b) per adjacent pair of A_perp letters without derivatives, '
                                                'modulo cyclically trivial terms (derived, Prop prop:tl-scalarsector)',
                  C2_Delta_elsewhere='open (task T3); composite leg sets in census.composite_W0',
                  pair_density_planar_scalar_sector_units_u='-8(1 - P + K/2) - (1/2)(1/3 + theta_perp) v.(X_a X_b)',
                  two_component='(A_ij, A_kl) -> u[8 A_kl A_ij - 2 eps_ijkl A^mn A_mn] - 8u A_ij A_kl (from 4 lc^natural) '
                                '- (1/2)(1/3+theta_perp) u ([(2E+2)A_ij] A_kl + A_ij [(2E+2)A_kl]) on a planar adjacent pair'),
    one_body=dict(table_units_u={'A_perp': '-22/3 X - theta_perp v.X; with C^(1)_Delta (local): -8 X + (5/3 - theta_perp) v.X',
                                 'A_V': '-22/3 X - theta_V v.X - kappa_V {Q int(iota_x A_V,c*), X}; without C^(1)_Delta (nonlocal)',
                                 'lambda': '-47/3 X - theta_lambda v.X - kappa_lambda {Q int(rho(x)lambda*,lambda*), X}; without C^(1)_Delta',
                                 'c': '-8/3 X - theta_V v.X - kappa_V {Q int(iota_x A_V,c*), X}; without C^(1)_Delta'},
                  two_component_A_perp='u[-8 + (5/3 - theta_perp)(2E + 2 + 2|alpha|)] d^alpha A_ij',
                  engine_A=D1['one_body_units']),
    planar_scalar_sector=D1['planar_scalar_sector'],
    bilinears=dict(engine_A=D1['bilinears'], statement='D_1 = 4 lc^natural - (1/3 + theta_perp) u {S1(Delta),O}: 0 and -48u '
                                                    'at theta_perp = -1/3 (Thm c9-D1, c9-kappa)', status='derived (Ch9 input)'),
    spectra=D1.get('spectrum_D1_conditional_on_P7_theta=-1/3'),
    conjecture_P7=PRED['P7_composite_part'],
    finite_parts='D_1 unchanged under omega-linear, translation-invariant, weight-preserving finite parts psi\' (Lemma lem:tl-psi)')

P1 = A['P1']
out['P1'] = dict(
    formula='P_1(O) = -{hbar F_1^(1)(P_b), O} + C_{P_b}(O) (Prop prop:tl-P1), no log-kernel term',
    two_body=dict(planar_scalar_sector='-u d_b(X_a X_b) per adjacent pair (symmetric part)', elsewhere='open (task T3)'),
    one_body_planar_scalar='u (5/3 - theta_P) d_b X',
    total_planar_scalar='P_1 = -(1/3 + theta_P) u {S1(P_b),-}',
    finite_part_parameter=dict(psi_prime='(a/2){S1(Delta),-} (Prop c9-scheme)',
                               effect='P_1 -> P_1 - a {S1(P_b),-}; K_b component +a, Sigma -a/2, Sigma* +a/2, D_1, M_1, R_1 unchanged',
                               zero_at='a = -(1/3 + theta_P) u'),
    engine_A=P1, status='derived on the bilinears and the planar scalar sector; open elsewhere')

out['arity_one'] = dict(
    structure='hbar L_1(O)^(1)(x) = sum_pairs K^(x)(X_a,X_b) + sum_letters K^(x)_1(X_a) (Prop prop:tl-arity1, proved)',
    two_body='finite local parts of the cycles through the pair with one defect insertion (special edge K_eps - K_0, S^(1)(x) '
             'on the cycle / composite line, or the arity-one counterterm vertex of O)',
    active_letters='all letters, antifields included (Lemma lem:tl-active1)',
    bound='Prop prop:tl-bound1', lambda_x='P_b -2, K_b +2, Sigma -1, Sigma* +1, Delta/M/R 0',
    scheme='Rem rem:tl-scheme1: F_1^(1)(x) and F_1^(0) one-body; finite parts psi\' via [{S1(x),-},psi\']; cutoff shape in '
           'special-edge finite parts',
    Q1_S1_K1_two_body_kernels='open (task T7)')

out['gates'] = dict(
    engine_A=GA, engine_B=GB, compare=CMP['gates'], compare_rerun_2026_10_29=CMPR,
    W2=dict(summary=W2['summary'], negative_controls=W2['negative_controls'],
            per_case=[dict(pair='|'.join(c['pair']), derivatives=c['derivatives'], W=c['W'],
                           tested=sum(1 for r in c['rows'] if 'pass' in r),
                           passed=sum(1 for r in c['rows'] if r.get('pass')),
                           nontrivial=sum(1 for r in c['rows'] if r.get('nontrivial')),
                           skipped=sum(1 for r in c['rows'] if 'skipped' in r)) for c in W2['cases']],
            rows='engine_a/out.nosync/w2_gate.json'),
    su11=SU,
    referee_numbers=dict(verdict=RN['verdict'], checks=RN['checks']),
    referee_theory=dict(verdict=RT['verdict'], checks=RT['checks']),
    totals=dict(engine_A='%d/%d' % (sum(g['pass'] for g in GA), len(GA)),
                engine_B='%d/%d' % (sum(1 for g in GB.values() if g.get('pass_', g.get('pass'))), len(GB)),
                compare='%d/%d' % (sum(1 for g in CMP['gates'] if g.get('pass_', g.get('pass'))), len(CMP['gates']))))

out['coverage'] = dict(
    two_engines='all 21 pair types at w^0_ab (all leg sets); w^0_ab+1 (all pairs) and w^0_ab+2 (AV|AV, AV|AW, AW|AW) with at '
                'most 2 outputs',
    engine_A_and_referee='w_ab up to w^0_ab+4 on sampled pairs (bosonic W=5,6; fermion pair with 3 derivatives at W=6; '
                         'A|lambda at W=11/2), k = 3 chains, ghost/antifield pairs (3833 entries)',
    W2='engine A on field-letter pairs (741/741, 232 nontrivial); not covered: (A_V,lambda) targets, derivative (lambda,c), '
       'most (A_V,A_V) targets at W_0+1 with three outputs, antifield-letter pairs',
    not_computed='C^(2)_Delta, C_{P_b}, C^(1)_Delta beyond the planar scalar sector; two-body kernels of Q_1, S_1, K_1; '
                 'field-strength families')

out['status'] = dict(
    proved=['Thm thm:tl-decomp', 'Lemma lem:tl-sign', 'Lemma lem:tl-k0', 'Lemma lem:tl-colour', 'Prop prop:tl-census',
            'Prop prop:tl-bound', 'Thm thm:tl-simplex', 'Prop prop:tl-shape', 'Lemma lem:tl-special', 'Prop prop:tl-genfun',
            'Prop prop:tl-bubble', 'Thm thm:tl-brst', 'Cor cor:tl-ward', 'Prop prop:tl-W1', 'Prop prop:tl-reduction',
            'Lemma lem:tl-psi', 'Thm thm:tl-D1 (under H1-H3)', 'Prop prop:tl-compbound (under H3)', 'Prop prop:tl-P1',
            'Prop prop:tl-spec', 'Prop prop:tl-arity1', 'Lemma lem:tl-active1', 'Prop prop:tl-bound1'],
    derived_from_Ch9=['Prop prop:tl-scalarsector (D_1, P_1 on the planar scalar sector)', 'D_1, P_1 on the bilinears'],
    computed_2=['W_0 kernels, 21 pair types', 'W_0+1, W_0+2 kernels (<= 2 outputs)', 'lc^natural spectra: su(2) L <= 8, '
                'twist two N <= 8, Konishi, bilinears', 'W1 gates', 'self-contraction', 'bubble closed form'],
    computed_1=['W2 gate (engine A)', 'su(1|1) density (engine A, referee)', 'higher-W kernels (engine A, referee)'],
    conjecture=['P7 (conj:tl-P7) beyond the planar scalar sector'],
    open=['C^(2)_Delta, C_{P_b} on twist two, fermion and field-strength families (T3)',
          'equality of theta_perp, theta_lambda, theta_V from the cubic arity-one equation',
          'W2 with antifield outputs/letters (T5)', 'Q_1, S_1, K_1 two-body kernels (T7)'])

json.dump(out, open(os.path.join(HERE, 'kernels_final.json'), 'w'), indent=1, default=str)
print('kernels_final.json written', os.path.getsize(os.path.join(HERE, 'kernels_final.json')), 'bytes;', out['gates']['totals'])
