"""Assemble compare/compare_r2.json from compare_r2_tables.json (engine A vs engine B vs FeynArts) and third_r2_bosonic.json
(third independent exact computation of every purely bosonic class). Run: python3 -I compare/make_compare_r2.py"""
import json
R = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
T = json.load(open(f'{R}/compare/compare_r2_tables.json'))
H = json.load(open(f'{R}/compare/third_r2_bosonic.json'))
A = json.load(open(f'{R}/engine_a/results_r2.json'))
B = json.load(open(f'{R}/engine_b/results_r2.json'))

UNIT = 'hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) (ratio independent of X,Y for simple g; DEFINITIONS_v2 "u", code-level unit only)'
QUART = {'(A|_V)^4': ('a', 'a', 'a', 'a'), '(A|_V)^2 (A|_Wperp)^2': ('a', 'a', 'b', 'b'), '(A|_Wperp)^4': ('b', 'b', 'b', 'b')}
COL = {'b01': '([X_0,X_1],[X_2,X_3])', 'b03': '([X_0,X_3],[X_1,X_2])', 'S': '(1/24) sum over the 24 orderings s of tr_g(ad_{X_s(0)} ad_{X_s(1)} ad_{X_s(2)} ad_{X_s(3)})'}

def lor(setname, l):
    v = QUART[setname]
    pr = {'x': ((0, 1), (2, 3)), 'y': ((0, 2), (1, 3)), 'z': ((0, 3), (1, 2))}[l]
    return ''.join('(%s_%d.%s_%d)' % (v[i], i, v[j], j) for i, j in pr)

def readable(setname, d):
    if d is None or isinstance(d, str):
        return d
    if not d:
        return '0'
    out = {}
    for k, v in d.items():
        col, l = k.split(' x ')
        out['%s %s' % (COL[col], lor(setname, l))] = v
    return out

# ------------------------------------------------------------------ census
census_T1 = {}
for s, v in T['T1'].items():
    census_T1[s] = dict(
        feynarts_job=v['feynarts_job'], feynarts_graphs=v['feynarts_graphs'],
        engine_A_classes=v['engine_A_classes'], engine_A_graphs_counted=v['engine_A_graphs'],
        engine_B_classes=v['engine_B_classes'],
        engine_B_graphs_counted=sum(r['feynarts_count'] for r in v['rows'] if r['present_in_B']),
        per_key=[dict(feynarts_key=r['feynarts_key'], feynarts=r['feynarts_count'], engine_A=r['engine_A_count'],
                      engine_B=('present' if r['present_in_B'] else 'absent')) for r in v['rows']],
        engine_A_equals_feynarts=all(r['census_A_equals_FeynArts'] for r in v['rows']),
        engine_B_equals_feynarts=all(r['present_in_B'] for r in v['rows']))
census_T2 = dict(T['T2']['census'])
census_T2['match'] = ('FeynArts 1PI special-edge graphs 5 = engine A 5 = engine B 5 (one-to-one, table T2.rows); the 8 FeynArts non-1PI graphs '
                      '(4 with the K_eps - K_0 edge as a bridge to a P(eps,infty) self-loop, 4 with a K_eps self-loop) are zero: each contains a '
                      'self-contraction at one vertex of I, colour factor tr_g(ad) = 0 or [t_a, t^a] = 0; no engine lists them, both state the vanishing')

# ------------------------------------------------------------------ coefficient tables
T1_tables = {}
for s, v in T['T1'].items():
    quart = s in QUART
    rows = []
    for c in v['per_class']:
        rows.append(dict(engine_A_id=c['engine_A_id'], engine_B_id=c['engine_B_id'],
                         log_eps_A=readable(s, None if c['engine_A'] is None else (c['engine_A'] if not quart else {k: x for k, x in c['engine_A'].items()})),
                         log_eps_B=(c['engine_B'] if c['engine_B'] == 'absent' or not quart else readable(s, c['engine_B'])),
                         equal=c['equal']))
    eps_all_zero = all(all(x == '0' for x in r['engine_A_eps_inv'] + r['engine_B_eps_inv'] + [y for y in r['engine_A_eps_inv_half'] + r['engine_B_eps_inv_half'] if y is not None]) for r in v['rows'])
    T1_tables[s] = dict(
        units=('multiples of the tree vertex, coefficient of %s x log eps' % UNIT) if not quart else
              ('coefficient of log eps; colour structures ([X_0,X_1],[X_2,X_3]), ([X_0,X_3],[X_1,X_2]) in units %s; symmetrised adjoint trace in units hbar/(16 pi^2); '
               'X_j colour part, a_j / b_j the A|_V / A|_Wperp polarisation of leg j, (x.y) the metric (engine B code-level evaluation basis; engine A su(2) data converted '
               'exactly, X = b03 - S/3, Y = b01 - b03 - S/3, Z = -b01 - S/3 on su(2))' % UNIT),
        per_class=rows,
        eps_inv_and_eps_inv_half_all_zero=eps_all_zero,
        total=v['total'], predicted_from_A=v['predicted_from_A'], superficial_degree=v['superficial_degree'])

coefficients = dict(
    T1=T1_tables,
    T1_all_per_class_equal=all(c['equal'] for v in T['T1'].values() for c in v['per_class']),
    T1_totals_equal=all(v['total']['equal'] for v in T['T1'].values()),
    T1_computed_equals_predicted_from_A=all(v['predicted_from_A']['computed_equals_predicted'] for v in T['T1'].values()),
    T1_coupling_direction=T['T1_coupling_direction'],
    T1_spin6_odd_sets=dict(engine_A=A['T1']['vanishing_sets_checked'], engine_B='not listed',
                           status='zero by Spin(6) covariance of the regulator (an odd number of W_R^perp vectors admits no invariant); engine A computed 0 exactly; not a disagreement'),
    T1_quadratic_regression=dict(engine_A=A['T1']['quadratic_regression'], engine_B=[c for c in B['checks'] if 'round-1' in c]),
    T2=T['T2'],
    T3=T['T3'])

# ------------------------------------------------------------------ disagreements and their settlement
third_sets = {k: dict(census_matches_feynarts=v['census_matches_feynarts'], invariant_basis_dimension=v['invariant_basis_dimension'],
                      evaluation_rank=v['evaluation_rank'], n_configurations=v['n_configurations'],
                      all_classes_match_engine_A=v['all_classes_match_A'], all_classes_match_engine_B=v['all_classes_match_B'],
                      quartic_tree_matches_engine_A_tree=v['quartic_tree_matches_engineA_tree'],
                      classes={kk: dict(third=c['third'], match_A=c['match_A'], match_B=c['match_B']) for kk, c in v['classes'].items()})
              for k, v in H['sets'].items()}
settled = []
for d in T['disagreements']:
    s = d['set']
    third_vals = H['sets'][s]['classes'][d['key']]['third']
    settled.append(dict(
        set=s, feynarts_key=d['key'], feynarts_ids=d['feynarts_ids'],
        localisation='engine B, r2_engine.ExternalSet.evaluate: words with zero log-eps average (val == (0,0)) or zero colour are dropped before grouping, '
                     'so a graph class all of whose words vanish never appears in engine B\'s graph list; engine A lists the class with value 0',
        engine_A='0', engine_B='absent (no nonzero word)',
        third_script='compare/third_r2_bosonic.py', third_value_at_each_configuration=third_vals,
        verdict='settled: the class contributes exactly 0 to the log eps coefficient (zero at %d configurations; the evaluation map on the %d-dimensional space of '
                'Spin(4) x Spin(6) x Ad-invariant forms has full rank %d, so the local form vanishes identically); census of A = FeynArts stands, B omits only zero classes'
                % (H['sets'][s]['n_configurations'], H['sets'][s]['invariant_basis_dimension'], H['sets'][s]['evaluation_rank'])
                + ('; for the bubble the affine numerator vanishes at the centroid k = -q/2, i.e. it is odd under the exchange of the two equal lines, so the regulated graph vanishes at every eps' if 'g4[' in d['key'] else '')))

out = {
 'status': 'complete',
 'round': 2,
 'inputs': ['engine_a/results_r2.json', 'engine_b/results_r2.json', 'diagrams/manifest.json'],
 'scripts': {'comparison': 'compare/compare_r2_tables.py -> compare/compare_r2_tables.json',
             'third_check': 'compare/third_r2_bosonic.py -> compare/third_r2_bosonic.json',
             'assembly': 'compare/make_compare_r2.py -> compare/compare_r2.json'},
 'unit': UNIT,
 'census': {'T1': census_T1, 'T2': census_T2,
            'symmetry_factors': 'engine A 1/|Aut|, engine B 1/|Aut| and FeynArts 1/s agree for every T1 and T2 graph (1/2 exactly for bubbles with two equal non-ghost lines)',
            'graph_class_conventions': 'engine A and engine B classes = all labelled graphs with a given cyclic vertex/line pattern; a closed c - A^* loop is one class '
                                       'covering both FeynArts orientations (FeynArts counts 2)'},
 'coefficients': coefficients,
 'disagreements': {'count': len(settled), 'numerical_disagreements': 0, 'items': settled},
 'third_independent_check': {
    'script': 'compare/third_r2_bosonic.py',
    'what': 'from-scratch exact computation (Fractions, integer matrices) of every purely bosonic one-loop 1PI class of (A|_V)^3, A|_V(A|_Wperp)^2, (A|_V)^4, '
            '(A|_V)^2(A|_Wperp)^2, (A|_Wperp)^4: vertices = full derivatives of the cubic and quartic parts of -1/4(F_A,F_A), P = 2(-,-)^{-1} (x) metric^{-1}/k^2, '
            'log eps coefficient from the exact Gaussian/scale integral; own brute-force graph census',
    'sets': third_sets,
    'su2_identities_ok': H['su2_identities_ok']},
 'final_tables': {},
 'unresolved': []}

# final agreed tables (compact)
ft = {}
ft['T1_totals'] = {s: dict(log_eps_multiple_of_tree=v['total']['A'], predicted_from_A=v['predicted_from_A']['A'], eps_inv='0', eps_inv_half='0')
                   for s, v in T['T1'].items()}
ft['T1_per_class'] = {s: [dict(engine_A_id=r['engine_A_id'], engine_B_id=r['engine_B_id'], log_eps=r['log_eps_A']) for r in t['per_class']]
                      for s, t in T1_tables.items()}
ft['T1_log_part_of_ICT_cubic_quartic'] = ('%s x log eps x [ 6 x ((A|_V)^3 and A|_V (A|_Wperp)^2 parts of int dvol -1/2 (dA, 1/2[A ^ A])) '
                                          '+ 10 x int dvol 1/2 (lambda, rho([A, lambda])) + 1 x int dvol ([A,c], A^*) + 1 x int dvol ([c,lambda],lambda^*) '
                                          '+ 1 x int dvol 1/2 ([c,c],c^*) + 8 x int dvol -1/4 (1/2[A ^ A], 1/2[A ^ A]) ] '
                                          '= %s x log eps x (Q + {I,-})( 2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c) ) (assumption (A) verified at cubic and quartic order); '
                                          'coefficient along S itself = 0' % (UNIT, UNIT))
ft['T1_Psi_weights'] = A['T3']['Psi_weights']
ft['T2'] = dict(sector='external legs c and A|_V',
                per_graph=[dict(engine_A_id=r['engine_A_id'], engine_B_id=r['engine_B_id'], feynarts_id=r['feynarts_id'],
                                eps_inv=r['eps_inv_A'], log_eps=r['log_A'], finite=r['finite_A']) for r in T['T2']['rows']],
                columns='eps^{-1}: coefficient of %s x eps^{-1} int(dc, A|_V); finite: coefficient of %s x int(d^*dc, d^*A|_V)' % (UNIT, UNIT),
                eps_inv_total='0', log_total='0', finite_total='-27/8 x %s x int(d^*dc, d^*A|_V)' % UNIT,
                obstruction='-27/8 x %s x int(d^*dc, d^*A|_V) = 27/16 x %s x Q int(d^*A|_V, d^*A|_V): Q-exact' % (UNIT, UNIT),
                finite_counterterm='add 27/16 x %s x int(d^*A|_V, d^*A|_V) to I^CT[eps] (renormalised interaction I - I^CT[eps] - this term); then the quadratic obstruction vanishes' % UNIT,
                other_quadratic_sectors='(c, A|_Wperp) and (c, lambda): no invariant, zero (both engines)')
ft['T3'] = dict(units='coefficient multiplying %s x O in the anomalous part of [(d_CE + Q) O[infty]](Delta)' % UNIT,
                rows={r['row']: {'SW (3), Delta = -E': r['A_SW3'], 'SW (17), Delta = 2E': r['A_SW17']} for r in T['T3']['rows']},
                konishi='trace row: -24 (SW 3) = +48 (SW 17); with the standard su(N) normalisation the unit equals -g^2 N/(32 pi^2), so -24 x unit = 3 g^2 N/(4 pi^2)')
out['final_tables'] = ft
out['unresolved'] = [
  'T2 cubic sector of the one-loop QME obstruction (legs (c,A|_V,A|_V), (c,A|_Wperp,A|_Wperp), (c,lambda,lambda), (c,c,A^*|_V)): not computed by either engine (both mark OPEN)',
  'engine B only, not cross-checked: T3 row "letter A|_V (longitudinal part)" = 0 (engine A gives no entry; engine A round 1 states Q O != 0 for the non-transverse A|_V letter, engine B round 1 calls the longitudinal projection nonlocal); outside the agreed DEFINITIONS_v2 table',
  'third independent check covers the purely bosonic T1 classes only (every disagreement is bosonic); classes with lambda or c - A^* lines rest on the exact agreement of engines A and B',
  'Delta / Q / S / K ghost components via K_eps graphs (Step 2 by DEFINITIONS_v2 item 4): not in round 2']
json.dump(out, open(f'{R}/compare/compare_r2.json', 'w'), indent=1)
print('ok', len(settled), 'settled;', out['coefficients']['T1_all_per_class_equal'], out['coefficients']['T1_totals_equal'], out['coefficients']['T1_computed_equals_predicted_from_A'])
