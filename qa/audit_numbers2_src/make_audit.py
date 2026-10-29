# Assemble qa/audit_numbers2.json from the fresh outputs and compare with compare/compare_r2.json.
import json
import os
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
UNIT = 'hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2)'

t1 = json.load(open(os.path.join(HERE, 't1_multi_out.json')))
t2 = json.load(open(os.path.join(HERE, 't2_out.json')))
t2b = json.load(open(os.path.join(HERE, 't2_out_v1.json')))
cmp = json.load(open(os.path.join(ROOT, 'compare', 'compare_r2.json')))


def one(block, key):
    v = t1[block][key]
    assert v['configuration_independent'], (block, key)
    return json.loads(v['values'][0])


def F(s):
    return Fraction(s)


# ---------------------------------------------------------------- T1 comparison
claim_AAA = {r['engine_A_id']: r['log_eps'] for r in cmp['final_tables']['T1_per_class']['(A|_V, A|_V, A|_V)']}
map_AAA = {
    'triangle A|_V loop': 'AV,AV,AV | AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~',
    'bubble A|_V loop (quartic)': 'AV,AV,AV | AAAA[A|_V,A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~',
    'triangle A|_Wperp loop': 'AV,AV,AV | AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~',
    'bubble A|_Wperp loop (quartic)': 'AV,AV,AV | AAAA[A|_V,A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~',
    'triangle c-A^* loop': 'AV,AV,AV | AcAs[A|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~',
    'triangle lambda loop': 'AV,AV,AV | llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~',
}
rows_AAA = []
for mine, theirs in map_AAA.items():
    a = one('AVAVAV', mine)
    b = claim_AAA[theirs]
    rows_AAA.append({'class (fresh script)': mine, 'compare_r2 engine_A_id': theirs, 'fresh': a, 'compare_r2': b,
                     'match': F(a) == F(b)})
tot_AAA = one('AVAVAV', 'TOTAL')
claim_tot_AAA = cmp['final_tables']['T1_totals']['(A|_V, A|_V, A|_V)']['log_eps_multiple_of_tree']

claim_Y = {r['engine_A_id']: r['log_eps'] for r in cmp['final_tables']['T1_per_class']['(lambda, lambda, A|_Wperp)']}
map_Y = {
    'external A|_Wperp on boson line (A|_V and A|_Wperp internal)':
        'l,l,AP | AAA[A|_Wperp] ~A|_V-A|_V~ -> llA[lambda] ~lambda-lambda~ -> llA[lambda] ~A|_Wperp-A|_Wperp~',
    'external A|_Wperp on lambda line, internal A|_V':
        'l,l,AP | llA[lambda] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[lambda] ~A|_V-A|_V~',
    'external A|_Wperp on lambda line, internal A|_Wperp':
        'l,l,AP | llA[lambda] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[lambda] ~A|_Wperp-A|_Wperp~',
}
rows_Y = []
for mine, theirs in map_Y.items():
    a = one('lam_lam_AWperp', mine)
    b = claim_Y[theirs]
    rows_Y.append({'class (fresh script)': mine, 'compare_r2 engine_A_id': theirs, 'fresh': a, 'compare_r2': b,
                   'match': F(a) == F(b)})
tot_Y = one('lam_lam_AWperp', 'TOTAL')
claim_tot_Y = cmp['final_tables']['T1_totals']['(lambda, lambda, A|_Wperp)']['log_eps_multiple_of_tree']

# ---------------------------------------------------------------- T2 comparison
S = t2['summary']
mass_tot = S['TOTAL coefficient of unit x eps^-1 int(A|_V,A|_V)']
long_tot = S['TOTAL coefficient of unit x int(d^*A|_V,d^*A|_V) in the renormalised one-loop functional']
log_tot = S['TOTAL longitudinal pol: log eps coefficient (in the unit)']
assert t2b['summary'] == S
PL = t2['per_loop']
KM = 'coefficient of unit x eps^-1 int(A|_V,A|_V)'
KB = 'coefficient of unit x int(d^*A|_V,d^*A|_V) (scheme-free only when the log coefficient is 0)'
KL = 'longitudinal pol: log eps coefficient of the W_1 two-point (in the unit)'


def sector(keys):
    import sympy as sp
    m = sum(sp.sympify(PL[k][KM]) for k in keys)
    b = sp.simplify(sum(sp.sympify(PL[k][KB].replace('gamma_E', 'gE')) for k in keys))
    lg = sum(sp.sympify(PL[k][KL]) for k in keys)
    return m, b, lg


sectors = {
    'A|_Wperp internal lines (cubic bubble + quartic tadpole)': (['A|_Wperp loop (cubic bubble)', 'A|_Wperp loop (quartic tadpole)'],
                                                                 ['QME2-b']),
    'lambda internal lines': (['lambda loop'], ['QME2-d']),
    'A|_V and c - A^* internal lines (gauge bubble + quartic tadpole + ghost bubble)': (
        ['A|_V loop (cubic bubble)', 'A|_V loop (quartic tadpole)', 'c-A^* loop'], ['QME2-a', 'QME2-c', 'QME2-e']),
}
rowsT2 = {r['engine_B_id']: r for r in cmp['final_tables']['T2']['per_graph']}
t2rows = []
for name, (mykeys, theirs) in sectors.items():
    m, b, lg = sector(mykeys)
    ce = sum(Fraction(rowsT2[t]['eps_inv']) for t in theirs)
    cf = sum(Fraction(rowsT2[t]['finite']) for t in theirs)
    # obstruction = Q(local part): eps^-1 column = -2 x (eps^-1 int(A,A) coefficient); finite column = -2 x (int(d^*A,d^*A) coefficient)
    pe = Fraction(str(-2 * m))
    pf = Fraction(str(-2 * b)) if b.is_Rational else None
    t2rows.append({
        'sector': name,
        'fresh: coefficient of unit x eps^-1 int(A|_V,A|_V) in W_1': str(m),
        'fresh: coefficient of unit x int(d^*A|_V,d^*A|_V) in W_1^R': str(b),
        'fresh: longitudinal log eps coefficient': str(lg),
        'fresh, mapped to the (c,A|_V) obstruction via Q c = -dc: eps^-1 int(dc,A|_V) column': str(pe),
        'fresh, mapped: int(d^*dc,d^*A|_V) column': str(pf),
        'compare_r2 special-edge graphs summed': theirs,
        'compare_r2 eps^-1 column (sum)': str(ce),
        'compare_r2 finite column (sum)': str(cf),
        'match': (pe == ce) and (pf == cf),
    })

audit = {
    'role': 'adversarial auditor, lens: numbers, round 2',
    'status': 'complete',
    'date': '2026-10-28',
    'unit': UNIT + ' (written out in full in every reported formula; the code evaluates it for g = su(2), (t_a,t_b) = delta_ab, hbar = 1, as -2/(16 pi^2))',
    'inputs_read': ['DESIGN.md', 'DEFINITIONS_v2.md', 'sources/ledger_SYNTHESIS.md',
                    'classical/classical_spec.md (conventions only: Hodge (F_A,F_A), Q, P(eps,infty) components)',
                    'compare/compare_r2.json (claims under audit)'],
    'independence': 'fresh exact scripts in qa/audit_numbers2_src (core.py, t1.py, t1_multi.py, t2.py, make_audit.py); no project code imported or copied; own 10d chiral blocks (9 real symmetric anticommuting 16x16 Pauli tensor products + one non-real generator in W_R^perp), own su(2) basis, own vertex tensors (symmetrised multilinear coefficients of the action on plane waves), own one-loop formula',
    'method': {
        'one_loop_functional': 'W_1 = (1/2) sum_n (1/n) Tr[(P V)^n] for the A sector and -(1/2) sum_n (1/n) Tr[(P V)^n] for the lambda sector and for the odd pair (c, A^*|_V), with P the Wick contraction values of P(eps,infty) and V the matrix of the background-field quadratic part of I (written 1/2 X^T V X); for two external lambda legs: one-loop quadratic form 1/2 lambda^T <Y G Y> lambda with background-corrected lambda and A lines',
        'propagators_rederived': 'A-A: 2 (-,-)^{-1} e^{-eps p^2}/p^2 (inverse of the Feynman-gauge form of -1/4(dA,dA)); lambda-lambda: i (-,-)^{-1} rho(p) e^{-eps p^2}/p^2 (= (-K)^{-1} for 1/2 lambda^T K lambda); c-A^*|_V and A^*|_V-c: -i (-,-)^{-1} p e^{-eps p^2}/p^2 (pseudo-inverse of (dc,A^*)); all agree with classical_spec section 7',
        'log_eps': 'exact: the log eps coefficient of int d^4k/(2pi)^4 F with every line cut off by e^{-eps k_e^2} is -(1/(16 pi^2)) x (S^3 average of the degree -4 part of the large-k expansion of F); expansion in the external momenta done exactly (polynomial-matrix algebra, dyadic floats converted to Fractions, S^3 moments as exact rationals)',
        'T2_route': 'renormalised self-energy route: exact Schwinger evaluation of every A|_V two-point graph with t_e in [eps,infty) (incomplete-Gamma closed form, exact eps -> 0 expansion, Euler constant and log 2 kept symbolic), longitudinal and transverse polarisations; obstruction = Q of the local longitudinal part, because Q of the transverse part and Q of the (c,A^*) two-point vanish in the (c,A|_V) sector ((dc, d^*dA) integrates to 0)',
        'configurations': 'T1: 20 random integer configurations, entries in [-9,9]; every per-class ratio identical at all 20 (a residual structure of degree <= 4 vanishes at a random point with probability <= 4/19). T2: p along a basis axis with colour t_1, and p = (1/2,1/2,1/2,1/2) with colour t_2; identical results',
    },
    'T1': {
        '(A|_V)^3 vertex, log eps coefficient as a multiple of the tree vertex, in the unit': {
            'per_class': rows_AAA,
            'total_fresh': tot_AAA,
            'total_compare_r2': claim_tot_AAA,
            'match': F(tot_AAA) == F(claim_tot_AAA),
        },
        '(lambda, lambda, A|_Wperp) vertex, log eps coefficient as a multiple of the tree vertex, in the unit': {
            'per_class': rows_Y,
            'total_fresh': tot_Y,
            'total_compare_r2': claim_tot_Y,
            'match': F(tot_Y) == F(claim_tot_Y),
        },
        'regressions (same script, same conventions)': {
            'A|_V two-point log eps coefficient, multiple of tree -1/4 int(dA,dA)': one('AV_two_point', 'TOTAL'),
            'A|_V two-point per internal lines (A|_V, A|_Wperp, lambda, c-A^*)': [one('AV_two_point', k) for k in
                                                                                  ['A|_V loop', 'A|_Wperp loop', 'lambda loop', 'ghost loop']],
            'A|_V two-point longitudinal log eps (A|_V and c-A^* loops cancel)': one('AV_two_point', 'longitudinal log (absolute, must be 0)'),
            'lambda two-point log eps coefficient, multiple of tree 1/2 int(lambda, dslash lambda)': one('lambda_two_point', 'TOTAL'),
            'agreement_with_DEFINITIONS_v2_round1': '4 x (-1/4 int(dA,dA)) = -int(dA|_V,dA|_V) and 8 x (1/2 int(lambda,dslash lambda)) = 4 int(lambda,dslash lambda): agree',
        },
        'consistency_with_s_exact_log_counterterm': 'with weights 2 (A), 4 (lambda): A^2 -> 4, A^3 -> 6, lambda^2 -> 8, lambda^2 A -> 10, as computed; i.e. the log eps parts of the cubic counterterms equal the cubic parts of log eps (Q + {I,-})(2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c)) for these two vertices',
        'verdict': 'CONFIRMED: (A|_V)^3 total 6 and (lambda,lambda,A|_Wperp) total 10 (log eps coefficients as multiples of the tree vertex, in the unit hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2)); every per-class value of compare_r2 reproduced exactly, including the vanishing (A|_V)^2 (A|_Wperp)^2-vertex bubble that engine B omits',
    },
    'T2': {
        'sector': 'external legs c and A|_V',
        'fresh totals': {
            'coefficient of the unit x eps^-1 int(A|_V,A|_V) in W_1': mass_tot,
            'longitudinal log eps coefficient (in the unit)': log_tot,
            'coefficient of the unit x int(d^*A|_V,d^*A|_V) in the renormalised one-loop functional': long_tot,
            'transverse log eps coefficient, multiple of tree -1/4 int(dA,dA)': S['TOTAL transverse pol: log eps coefficient as a multiple of the tree -1/4(dA,dA)'],
            'transverse eps^-1 equals longitudinal eps^-1 (local term proportional to int(A|_V,A|_V))': S['TOTAL transverse pol: eps^-1 coefficient (in the unit)'] == S['TOTAL longitudinal pol: eps^-1 coefficient (in the unit)'],
        },
        'fresh obstruction': 'special-edge part of the eps^-1 column: -2 x 9/4 = -9/2 (unit x eps^-1 int(dc,A|_V)), cancelled by -Q I^CT[eps] (+9/2); log column 0; finite: Q((27/16) x unit x int(d^*A|_V,d^*A|_V)) = -(27/8) x unit x int(d^*dc,d^*A|_V), using Q c = -dc in the A|_V slot (classical_spec section 5)',
        'compare_r2 claim': cmp['final_tables']['T2']['obstruction'],
        'per_sector': t2rows,
        'remark_per_graph': 'Matter sectors match graph by graph (A|_Wperp lines <-> QME2-b: eps^-1 -9/2, finite -5/4; lambda lines <-> QME2-d: eps^-1 2, finite -4/3). In the A|_V / c-A^* sector the fresh route splits into an A|_V-line part (log eps coefficient -1/4, Euler constant and log 2 present) and a c-A^* part (+1/4, opposite constants); only their sum (eps^-1 -2, finite -19/24 in the obstruction columns) is regulator-meaningful, and it equals QME2-a + QME2-c + QME2-e. The individual special-edge values -9/4, 1/8, 1/8 (eps^-1) and -5/8, -1/12, -1/12 (finite) were not recomputed graph by graph.',
        'verdict': 'CONFIRMED: eps^-1 total 0 after the tree counterterm vertex, log total 0, finite total -(27/8) x hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) x int(d^*dc, d^*A|_V) = (27/16) x hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) x Q int(d^*A|_V, d^*A|_V); Q-exact; removed by the finite counterterm (27/16) x hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) x int(d^*A|_V,d^*A|_V)',
    },
    'disagreements': [],
    'caveats': [
        'Convention dependence (not a disagreement): all multiples are in the unit hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) with (F_A,F_A) the Hodge pairing (classical_spec section 4). With -1/4 sum_{M,N}(F_MN,F_MN) instead, every number in this unit halves; ratios such as cubic/quadratic = 3/2 do not change.',
        'su(2) only; adjoint one-loop colour factors of these graphs reduce to f-structure times tr_g(ad ad) for any simple g (the symmetric part of tr(ad ad ad) vanishes for a real representation), so the numbers in this unit are g-independent.',
        'numpy on this machine (Accelerate BLAS) emits spurious floating-point warnings on every matmul, also for trivial integer matrices; exactness checked against Python-integer object arrays (max difference 0) and by all outputs being exact rationals that agree across 20 configurations.',
        'T2 special-edge graphs QME2-a, QME2-c, QME2-e not recomputed individually (sum confirmed); cubic part of the T2 obstruction not computed (also OPEN in compare_r2).',
    ],
    'resources': 'each run < 0.1 GB resident, < 4 s; one process at a time; no Wolfram kernel used; no files written outside the project folder; no processes left running',
    'outputs': ['qa/audit_numbers2_src/t1_multi_out.json', 'qa/audit_numbers2_src/t2_out.json', 'qa/audit_numbers2_src/t2_out_v1.json'],
}
json.dump(audit, open(os.path.join(ROOT, 'qa', 'audit_numbers2.json'), 'w'), indent=1)
print(json.dumps({'T1_AAA_match': audit['T1']['(A|_V)^3 vertex, log eps coefficient as a multiple of the tree vertex, in the unit']['match'],
                  'T1_Y_match': audit['T1']['(lambda, lambda, A|_Wperp) vertex, log eps coefficient as a multiple of the tree vertex, in the unit']['match'],
                  'T1_rows_match': all(r['match'] for r in rows_AAA + rows_Y),
                  'T2_rows': [(r['sector'], r['match']) for r in t2rows]}, indent=1))
