"""Assemble step2/engine_a/cubic_results.json from the logs of s2_verify0d, s2_regress2pt, s2_integrals, s2_Ia, s2_ga_tests,
s2_trees, s2_numeric."""
import json, os, re, glob
import sympy as sp

H = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_a/'
L = lambda f: json.load(open(H + 'logs/' + f))
LOGB = ['1', 'log2', 'log3']
SECT = {'c,AV,AV': '(c, A|_V, A|_V)', 'c,AP,AP': '(c, A|_Wperp, A|_Wperp)', 'c,l,l': '(c, lambda, lambda)',
        'c,c,AsV': '(c, c, A^*|_V)'}
GDEF = {'G1': 'int dvol ([c, A|_V], d^*d A|_V)', 'G2': 'int dvol ([c, A|_V], d d^* A|_V)',
        'G3': 'int dvol ([c, A|_Wperp], d^*d A|_Wperp)', 'G4': 'int dvol ([c, dc], A^*)',
        'G5': 'int dvol ([c, lambda], dslash lambda)'}

runs = {'su2_seed2028': L('s2_Ia_su2_cAVAV-cAPAP-cll-ccAsV_seed2028.json'),
        'su2_seed31': L('s2_Ia_su2_cAVAV-cAPAP-cll-ccAsV_seed31.json'),
        'gl3(sl3)_seed7': L('s2_Ia_gl3_cAVAV-cAPAP-cll-ccAsV_seed7.json')}
trees = {'su2': L('s2_trees_su2.json'), 'gl3(sl3)': L('s2_trees_gl3.json')}
numer = {s: L('s2_numeric_%s_su2.json' % s.replace(',', '')) for s in SECT}
ints = L('s2_integrals.json')
reg = {a: L('s2_regress2pt_%s.json' % a) for a in ('su2', 'gl3')}
gat = {a: L('s2_ga_tests_%s.json' % a) for a in ('su2', 'gl3')}
v0d = open(H + 'logs/s2_verify0d.log').read()


def project(fit, Gin, names):
    """structure-basis dict -> G-basis (exact), Gin = {G: {structure: coeff}}"""
    Gn = sorted(Gin)
    M = sp.Matrix([[sp.sympify(Gin[g].get(nm, '0')) for g in Gn] for nm in names])
    b = sp.Matrix([sp.sympify(fit.get(nm, '0')) for nm in names])
    if all(x == 0 for x in b):
        return {g: '0' for g in Gn}, True
    sol, params = M.gauss_jordan_solve(b)
    sol = sol.subs({p: 0 for p in params})
    ok = all(sp.simplify(x) == 0 for x in (M * sol - b))
    return {g: str(sp.nsimplify(x)) for g, x in zip(Gn, sol)}, ok


out = {'engine': 'A', 'part': 'C (cubic sector of the arity-zero one-loop obstruction)',
       'unit': 'hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))  (code name: unit)',
       'G_basis': GDEF,
       'definition': 'Ob_3 = lim_{eps->0} [I-a + II-a + III(K_0)] (three-leg part of the order-hbar Q W(P(eps,infty), I - I^CT[eps]) '
                     'minus the limit of I-d + II-b, which is the tree of Ob_2 with I_3, and minus I-b + III(K\'), which tends to 0)',
       'sectors': {}}
tot_runtime = 0.0
for s, name in SECT.items():
    base = runs['su2_seed2028'][s]
    sec = {'legs': name, 'structure_basis': base['structures'], 'G_in_structure_basis': base['G_in_structure_basis']}
    names0 = base['structures']['0']
    # ---- I-a classes (agreement across runs)
    cls = []
    others = {k: {c['class']: c for c in r[s]['classes']} for k, r in runs.items()}
    for c in base['classes']:
        agree = all(others[k][c['class']]['coeffs']['eps^0'] == c['coeffs']['eps^0'] for k in others)
        cls.append({'class': c['class'], 'topology': c['topology'], 'vertices': c['vertices'], 'lines': c['lines'],
                    'special_edge': c['special_edge'], 'multiplicity_and_weight': c['multiplicity_and_weight'],
                    'coeffs': {'eps^-1': c['coeffs']['eps^-1'], 'log': c['coeffs']['log'], 'eps^0': c['coeffs']['eps^0']},
                    'eps^0_in_span_G_zero_residual': c['coeffs']['eps^0_in_span_G'],
                    'same_in_su2(2 seeds)_and_sl3': agree})
    sec['graph_classes_I-a'] = cls
    sec['I-a_total'] = {'eps^-1': base['total_eps^-1'], 'log': '0 (every graph separately)', 'eps^0': base['total_eps^0_G'],
                        'eps^0_in_span_G_zero_residual': base['total_eps^0_in_span_G_zero_residual'],
                        'all_class_fits_zero_residual': base['all_fits_zero_residual'],
                        'structure_matrix_full_column_rank': base['structure_matrix_full_column_rank'],
                        'same_total_in_all_runs': all(r[s]['total_eps^0_G'] == base['total_eps^0_G'] for r in runs.values()),
                        'samples_per_run': {k: r[s]['samples'] for k, r in runs.items()}}
    # ---- trees
    T = trees['su2'][s]; Tg = trees['gl3(sl3)'][s]
    fits = T['fits_structure_basis']
    pj = {k: project(fits[k]['structures'], base['G_in_structure_basis'], names0) for k in fits}
    same_gl3 = all(Tg['fits_structure_basis'][k]['structures'] == fits[k]['structures'] for k in fits)
    sec['II-a_and_III(K_0)'] = {
        'II-a = -Q I^CT_3': {'log': pj['IIa_log'][0], 'eps^-1': '0', 'eps^0': '0'},
        'III(K_0) = -{I_3, I^CT_2}_0': {'log': pj['IIIK0_log'][0], 'eps^-1': '0 (identically)' if T['IIIK0_epsinv_identically_zero'] else 'NONZERO', 'eps^0': '0'},
        'sum_log_exactly_zero': T['exact_zero_checks']['IIa+IIIK0_log'],
        'Y_3 = {I_3, Y} with Step-1 multiples': T['exact_zero_checks']['Y3_minus_bracket_I3_Y'],
        'note': 'log eps (Q Y_3 + {I_3, Y_2}) = 0: Q{I_3,Y} = -{I_3,QY}, Q I_3 = 0'}
    sec['I-b_and_III(K_eps-K_0)'] = {
        'I-b eps^0': pj['Ib_eps0'][0], 'III(K\') eps^0': pj['IIIKp_eps0'][0],
        'sum eps^0 exactly zero': T['exact_zero_checks']['Ib+IIIKp_eps0'],
        'eps^-1': '0 (each)', 'log': '0 (each, through eps^0)',
        'expansion': 'I-b + III(K\') = (e^{-eps|p_l|^2}-1) x (renormalised self-energy at p_l) contracted with I_3 by Khat = '
                     '-eps |p_l|^2 x (finite part) + O(eps^2 log eps): eps^0 parts +X and -X, eps log eps parts cancel',
        'self-energy eps^-1 coefficient from graphs (multiple of int(A,A))': trees['su2']['selfenergy_eps^-1_from_graphs']}
    sec['I-d_and_II-b'] = {
        'eps^-1 (nonlocal, 1/|p_l|^2)': 'cancel exactly' if T['exact_zero_checks']['Id+IIb_epsinv'] else 'FAIL',
        'I-d eps^0 local part (-C_Pnum(sb_-1, I_3))': pj['Id_eps0_local'][0],
        'II-b eps^0 local part (+C_Pnum(sb_-1, I_3))': pj['IIb_eps0_local'][0],
        'tree of Ob_2 with I_3 (limit)': pj['tree_Ob2_I3'][0],
        'I-d + II-b - tree at eps^0 exactly zero': T['exact_zero_checks']['Id+IIb_eps0_minus_tree'],
        'note': 'the tree C_{Pnum/|p|^2}(-(27/8) int(d^*dc,d^*A|_V), I_3) is local because Ob_2 carries |p|^2; it is the flow of Ob_2, '
                'not part of Ob_3'}
    sec['tree_checks_same_in_sl3'] = same_gl3
    sec['mode_model_checks'] = {k: v for k, v in T['exact_zero_checks'].items()}
    alpha = base['total_eps^0_G']
    sec['totals'] = {'eps^-1': '0', 'log': '0', 'Ob_3 (eps^0, local)': alpha,
                     'cubic part of lim Q W (Ob_3 + tree(Ob_2,I_3))':
                         {g: str(sp.nsimplify(sp.sympify(alpha[g]) + sp.sympify(pj['tree_Ob2_I3'][0][g]))) for g in alpha}}
    N = numer[s]
    sec['numerical_check_finite_eps'] = {
        'configuration_momenta': N['momenta'], 'predicted_limit': N['predicted_limit'],
        'by_eps': {e: {k: v[k] for k in ('total - predicted limit', 'I-a - alpha.G', 'I-d+II-b - tree', "I-b+III(K')", 'eps*log(eps)', 'seconds')}
                   for e, v in N['eps'].items()},
        'conclusion': 'total - predicted limit -> 0 like eps log eps (no constant, no log eps, no nonlocal remainder)'}
    for k in ('lambda_selfenergy_ratio_polarisation_independent', 'Wperp_selfenergy_ratio_polarisation_independent',
              'AV_selfenergy_two_structure_fit_residual'):
        if k in N:
            sec['numerical_check_finite_eps'][k] = N[k]
    out['sectors'][name] = sec
    tot_runtime += sum(r[s]['seconds'] for r in runs.values()) + T['seconds'] + Tg['seconds'] + sum(v['seconds'] for v in N['eps'].values())

A = {}
for s, name in SECT.items():
    A.update(out['sectors'][name]['totals']['Ob_3 (eps^0, local)'])
out['alpha'] = {'alpha_1': A['G1'], 'alpha_2': A['G2'], 'alpha_3': A['G3'], 'alpha_4': A['G4'], 'alpha_5': A['G5']}
out['alpha_simplified'] = {'alpha_1': '11/4 + (19/4) log(3/2)', 'alpha_2': '27/8', 'alpha_3': '11/4 + (19/4) log(3/2)',
                           'alpha_4': '0', 'alpha_5': '-9 log(3/2)'}
out['Ob_3'] = 'unit * [ (11/4 + (19/4) log(3/2)) (G1 + G3) + (27/8) G2 + 0 G4 - 9 log(3/2) G5 ]'
out['alpha_numeric_20digits'] = {k: str(sp.N(sp.sympify(v), 20)) for k, v in out['alpha'].items()}
out['tree_Ob2_I3_total'] = 'unit * [ (27/8)(G1 - G2) + (27/8) G3 - (27/4) G4 - (27/4) G5 ]'
out['validation'] = {
    '0d_identity (s2_verify0d)': {'runs': [l for l in v0d.splitlines() if l.startswith('{')],
                                  'note': 'mismatch_deg3 = 0: Delta_K W_0|_3 = (1/2)str[KSH3PSH3PSH3] + (1/2)(str[KSH3PSH4]+str[KSH4PSH3]) '
                                          '+ C_P((1/2)str[KSH3PSH3], I_3) (I-a + I-d); the deg-2 d_P((1/2){I,I}_K) comparison differs only '
                                          'by P-tadpoles when d_P I_3 != 0 (absent in the field theory, tr ad = 0) and matches for '
                                          'P-tadpole-free I_3 (seeds 4, 5); the wrong sign of I-a fails (control).'},
    'two-leg regression vs Step 1 (s2_regress2pt)': {a: {'TOTAL': reg[a]['TOTAL'], 'per_graph': {k: v for k, v in reg[a].items() if '|' in k}}
                                                     for a in reg},
    'mode model tests (s2_ga_tests)': gat,
    'exact integrals (s2_integrals)': {k: {'exact [1, log2, log3]': v['exact'], 'independent 50-digit evaluation': v.get('numeric_50digits'),
                                           'abs_diff': v.get('abs_diff'), 'pslq [val,1,log2,log3]': v.get('pslq_[val,1,log2,log3]')}
                                       for k, v in ints.items()}}
out['method_notes'] = [
    'D_0 part of Q W at order hbar = hbar Delta_{K\'} W_0 + {W_0, W_1}_{K\'} (exact identity, s2_notes.txt S1); I-a = 1PI part = '
    'hbar (1/2) str[K\' S H (P S H)^{n-1}], validated as a polynomial identity in a 0-dim super model (s2_verify0d.py) and by the exact '
    'reproduction of the Step-1 two-leg special bubbles (s2_regress2pt.py).',
    'Each oriented graph: numerator N(k) reconstructed exactly as a polynomial (Newton forward differences on the simplex lattice, '
    'checked at extra points), Gaussian loop integral with moments (Delta^j N)(-K(u))/(4^j j! T^j), Schwinger variables scaled by eps '
    '(special tau_1 in [0,1], others in [1,infty)); eps^e coefficient = sum_{j-i=n-2-e} (4^j j! i!)^{-1} int T^{-2-j+i} F_j (-q)^i.',
    'Convergence: integrand <= C T^{-2-j+i}, large-tau region of dimension n-1, T >= n-1 >= 1; converges iff j - i > n - 3, i.e. e < 1. '
    'Remainder O(eps log eps); no log eps term in any I-a graph.',
    'Exact integrals J(n;beta;N) in Q + Q log2 + Q log3 (sympy iterated rational integration), each confirmed by an independent '
    'evaluation in (T,u) variables at 50 digits and PSLQ on {1, log2, log3}.',
    'Projection: each class fitted exactly (sympy over Q(i)) on a complete basis of Spin(4) x Spin(6) x colour-f invariant forms of the '
    'legs (parity-odd and chirality-weighted ones included) with extra samples; full column rank checked; then expressed in G1..G5 with '
    'zero residual.',
    'Colour: every I-a colour factor is a cyclic contraction of three structure constants with (-,-)^{-1}: tr_g(ad_X ad_Y ad_Z) = '
    '(1/2) tr_g(ad_{[X,Y]} ad_Z) = (1/2) [tr_g(ad ad)/(-,-)] ([X,Y],Z) for every simple g (ad_X antisymmetric for (-,-)); '
    'confirmed: identical per-class coefficients for su(2) (orthonormal) and sl(3) = complexified su(3) (embedded in gl(3) with '
    'the unimodular trace form; the centre has zero brackets and decouples), fits with the f-structure only and zero residual.',
    'Tree-level pieces (II-a, III, I-b, I-d, II-b) evaluated exactly in a plane-wave mode model with the BV bracket, Q as the derivation '
    '{S_free,-}, kernels Khat and Pnum (tests: Q^2 = 0, Q Pnum = |p|^2 Khat, derivation = bracket, mode-model vertices = canonical '
    'vertex functions, G_i formulas).',
    'Independent internal check: in all four sectors the full regulated sum (I-a at finite eps without expansion, I-b + III(K\') with '
    'the finite-eps self-energies, I-d + II-b with the finite-eps special bubbles; II-a + III(K_0) = 0 exactly) at eps = 1e-4, 1e-5, '
    '1e-6 and generic momenta approaches alpha.G + tree(Ob_2, I_3) like eps log eps.',
    'No dimensional regularisation; no IR scale; Euclidean; upper Schwinger limit infinity.']
out['runtime_seconds_total'] = round(tot_runtime, 1)
out['files'] = sorted(os.path.basename(f) for f in glob.glob(H + 's2_*.py')) + ['s2_notes.txt']
json.dump(out, open(H + 'cubic_results.json', 'w'), indent=1)
print(json.dumps(out['alpha'], indent=1))
for name, sec in out['sectors'].items():
    print(name, sec['totals'], sec['I-a_total']['same_total_in_all_runs'], sec['tree_checks_same_in_sl3'],
          {k: v['total - predicted limit'] for k, v in sec['numerical_check_finite_eps']['by_eps'].items()})
print('runtime', out['runtime_seconds_total'])
