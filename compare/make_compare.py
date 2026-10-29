"""Assemble compare.json from compare_tables.json (engine A vs B vs census vs FeynArts) and third_ghost.json."""
import json
R = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/compare'
T = json.load(open(f'{R}/compare_tables.json'))
G3 = json.load(open(f'{R}/third_ghost.json'))

UNIT = 'hbar tr_g(ad_X ad_Y) / ((X,Y) 16 pi^2)  [= both engines\' hbar c_adj/(16 pi^2)]'
BASIS = {
 'B1': 'int dvol (dA|_V, dA|_V)  [Hodge pairing (x) (-,-)]',
 'B2': 'int dvol (d^* A|_V, d^* A|_V)',
 'B3': 'int dvol (A|_V, A|_V)',
 'B4': 'int dvol (dA|_{W^perp}, dA|_{W^perp})',
 'B5': 'int dvol (A|_{W^perp}, A|_{W^perp})',
 'B6': 'int dvol (lambda, dslash lambda)',
 'B7': 'int dvol (dc, A^*)  (only A^*|_V contributes)'}

final = {
 'unit': UNIT,
 'basis_written_out': BASIS,
 'ICT_quadratic': {k: v['A'] for k, v in T['ICT_quadratic'].items()},
 'ICT_quadratic_zero_sectors': 'A|_V-A|_{W^perp} mixing, (lambda,A), (A^*|_{W^perp}, c): 0',
 'ICT_two_point_per_graph': [{'census': e['census_index'], 'feynarts_id': e['feynarts_id'],
                              'eps_inv': e['eps_inv']['A'], 'log_eps': e['log_eps']['A']} for e in T['two_point']],
 'partial_delta_ICT_scalar_bilinear': {
   'O = int omega (iota_{w1}A, iota_{w2}A); log eps coefficient': {
      'O': '4', '(w1,w2) int omega sum_w (iota_w A, iota_w A)': '-2', '(w1,w2) int omega (A|_V, A|_V)': '0'},
   'eps^-1, field independent [unit hbar/(16 pi^2)]': '2 dim(g) (w1,w2) int omega',
   'eps^-1, field dependent': '0',
   'per graph': [{'census': e['census_index'], 'feynarts_id': e['feynarts_id'], 'eps_inv': e['eps_inv']['A'],
                  'log_eps': e['log_eps']['A']} for e in T['insertion']]},
 'O_eps': {
   'traceless symmetric': 'O(eps) = O - 4 u log(eps) O',
   'trace': 'sum_w O_ww(eps) = (1 + 8 u log eps) sum_w O_ww - 12 hbar dim(g) /(16 pi^2 eps) int omega',
   'single letters (iota_w A, A|_V, lambda, c)': 'O(eps) = O',
   'u': UNIT},
 'dilatation_component': {
   'normalisation used by both engines': 'Delta = -E [SW (3)], on a field of weight w: -(E+w); weights A 1, lambda 3/2, c 0, A^* 3, lambda^* 5/2, c^* 4',
   'form': '[(d_CE+Q) O[infty]](Delta) = -(E + w_O + a u) O ; table gives a',
   'a': {'traceless symmetric bilinear': '0', 'trace bilinear': '-24', 'letter iota_w A': '-4',
         'letter A|_V (transverse)': '-4', 'letter lambda': '-8', 'letter c': '-2 (engine A; = -2 x ghost-vertex weight 1 of the third check)'},
   'routes': 'route (i) and route (ii) give the same a in each engine',
   'conversion to SW (17) normalisation (DESIGN)': 'Delta_(17) = -2 Delta_(3) on bosonic generators [ledger X13]: [(d_CE+Q)O[infty]](Delta_(17)) = 2(E + w_O + a u) O'},
 'konishi': {'u': '-g^2 C_A/(32 pi^2)', 'trace': '-(E + 2 + 3 g^2 C_A/(4 pi^2)) O_tr', 'su(N)': '2 + 3 g^2 N/(4 pi^2) = standard'}}

out = {
 'status': 'complete',
 'inputs_note': 'engine_a/report.md and engine_b/report.md do not exist; compared from results.json, logs and code',
 'census': {
   'two_point_matched': [{k: e[k] for k in ('census_index', 'census_Aut', 'census_omega', 'feynarts_id', 'feynarts_present', 'engine_A_id', 'engine_B_id', 'A_sym', 'B_sym')} for e in T['two_point']],
   'insertion_matched': [{k: e[k] for k in ('census_index', 'feynarts_id', 'feynarts_present', 'engine_A_id', 'engine_B_id')} for e in T['insertion']],
   'extra_engine_entries_identically_zero': T['extra_engine_entries'],
   'feynarts_boundary_counts': T['feynarts_boundary'],
   'mismatch_explanations': [
     'census 7 (lambda loop, external A|_V,A|_{W^perp}): engine A has no per-graph entry (sector total 0 in logs/sector_a_su2.json and _su3.json); engine B entry "AV-Ap:AAA+AAA:all" is a run with all cubic vertices incl. lambda lambda A and no propagator mask, raw value exactly 0 (label "AAA+AAA" and note "no A_V A_V A_perp vertex" do not name the lambda loop)',
     'O2 (tadpole, one external A|_V): neither engine lists it as a graph; both state 0 by colour (bracket contracted with (-,-)^{-1})',
     'omega<0 graph (lambda,lambda with O): census and FeynArts exclude; engine A states degree<0; engine B silent',
     'engine-only zero entries: G^W ghost loop and (A^*,c) with A|_{W^perp} line (no c (x) A^*|_{W^perp} component of P; FeynArts tp_cAsW = 0 graphs); engine A O-triangles with (A|_{W^perp})^3 or A|_V A|_V A|_{W^perp} vertices (no such vertices); engine B quartic AV-Ap (no A|_V (x) A|_{W^perp} component of P)',
     '22 one-particle-reducible insertion graphs (FeynArts): census, engine A and engine B all treat them as finite after the quadratic I^CT (tadpole subgraphs vanish; two-point subgraphs paired with I^CT trees)',
     'FeynArts route (i) = 29 L_eps-pinned graphs = the set evaluated as "route (ii)" by both engines (and used by both for the log coefficients of route (i))',
     'FeynArts route (ii) = 21 1PI K_eps graphs with an S^(1)(Delta) vertex + 4 partial_delta I^CT - K_eps trees: evaluated by neither engine (UNRESOLVED)',
     'symmetry-factor labels: census |Aut| = FeynArts s; engine B prints 1/|Aut|; engine A prints the trace weight 1/(2n) x (orderings of external legs) = 1/2 for every bubble, which times the number of internal label assignments in the trace (1 or 2) equals 1/|Aut|']},
 'coefficients': {
   'two_point_per_graph_all_equal': all(e['eps_inv']['equal'] and e['log_eps']['equal'] for e in T['two_point']),
   'insertion_per_graph_all_equal': all(e['eps_inv']['equal'] and e['log_eps']['equal'] for e in T['insertion']),
   'ICT_quadratic': T['ICT_quadratic'],
   'partial_delta_ICT': T['partial_delta_ICT'], 'partial_delta_ICT_equal': T['partial_delta_ICT_equal'],
   'gamma': T['gamma'], 'gamma_equal': T['gamma_equal'],
   'intermediate_conversions': {'engineA_g_vs_engineB_C0_C1': T['intermediate_conversion'],
                                'engineB_integrals_vs_final': T['engineB_integral_conversion'],
                                'engineA_faces': 'each face value (structure normalisation, structure = 1/2 the basis cross value) is 1/2 of the graph total for every two-point bubble; O triangles 1/3 per face, O bubbles 1/2'},
   'engineA_su2_equals_su3': T['su2_su3_sector_totals_engineA_equal'],
   'disagreements': []},
 'third_independent_check': {
   'script': 'compare/third_ghost.py (exact sympy; own Grassmann algebra compare/grassmann3.py)',
   'formula_test': 'compare/toy3.py: one-loop str formula (right slot first, parity operator between P and d2I) = brute-force hbar log(e^{hbar d_P} e^{I/hbar}) at n=1,2,3 with odd backgrounds; the three other conventions fail',
   'results': G3,
   'agrees_with_engines': {'AV-AV ghost loop': G3['AV-AV ghost loop'] == {'eps_inv': {'B1': '0', 'B2': '0', 'B3': '1/16'}, 'log_eps': {'B1': '1/24', 'B2': '1/8', 'B3': '0'}},
                           'B7': G3['As-c graph (A_V, c-A*)']['log_eps B7'] == '-1'}},
 'sanity': {
   'traceless_gamma_zero': True,
   'coupling_direction': {
     'two_point_alone': 'quadratic part of the coupling direction S - (A^* d_{A^*} + lambda^* d_{lambda^*} + c^* d_{c^*}) S equals {S_free, 1/2 int(A,A^*) - 1/2 int(lambda^*,lambda) - 1/2 int(c^*,c)}: Q-exact; the four two-point log coefficients fix the rescaling weights only up to one free parameter = the coupling-direction coefficient; transversality (B2 = 0) and 10d covariance (B1 = B4 log, B3 = B5 pole) hold',
     'vertex': 'log coefficient of int ([A,c],A^*) at one loop = 1 x tree (third check); rescaling-weight content of this monomial = weight(c) = weight(A|_V) - 1, so coupling-direction coefficient = 0',
     'result': '0'},
   'konishi': final['konishi']},
 'final_tables': final,
 'unresolved': [
   'FeynArts route (ii) set (21 + 4 graphs with K_eps and S^(1)(Delta)) not evaluated by either engine; the UV-boundary identity in its Delta_eps = -partial_{K_eps} form is not yet checked graph by graph',
   'both engines assemble the dilatation coefficients with the assumption that the log-eps part of I^CT is (Q+{I,-}) of 2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c) and that the Delta-component of the action is corrected by 2 hbar times it; checked at quadratic order (both engines) and at the ghost vertex (third check); cubic/quartic vertex log coefficients (AAA, lambda lambda A, A|_V A|_{W^perp}^2, quartics) computed by nobody',
   'dilatation normalisation: engines use SW (3)/(4) Delta = -E, DESIGN asks SW (17); conversion factor -2 (ledger X13) documented, values must be restated before use',
   'ledger notation: engine outputs use symbols not in the ledger (c_adj, B1..B7, u, Psi, R, kappa, gamma, a_k, b_k); compare/report.md writes them out',
   'engine A two-point and O coefficients are float64 least-squares fits with rational reconstruction (limit_denominator 10^4), not exact arithmetic; each value is reproduced exactly by engine B (sympy), three also by the third check',
   'engine B route (ii) numerics cover only the A|_{W^perp} sector; first O-vertex run (route2.log) off by 1.3% at eps p^2 ~ 1, superseded by route2_o.log (momenta /10): 3.7e-6',
   'letter c coefficient reported by engine A only',
   'selection rule for K_mu: engine A "x^2 / two x" vs engine B "-2 a int omega x_mu O + Q-exact" (wording, not a computed coefficient)']}
json.dump(out, open(f'{R}/compare.json', 'w'), indent=1)
print('written', len(json.dumps(out)))
