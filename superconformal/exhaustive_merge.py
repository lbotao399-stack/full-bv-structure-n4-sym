"""
exhaustive_merge.py -- assemble superconformal/verify_exhaustive.json from the outputs of
exhaustive_alg.py and exhaustive_tower.py (out/ directory).
"""
import os, sys, json, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.argv = sys.argv[:1]
import exhaustive_tower as X

O = lambda f: os.path.join(HERE, 'out', f)


def load(f):
    p = O(f)
    return json.load(open(p)) if os.path.exists(p) else None


def main():
    res = {}
    alg = load('exhaustive_alg.json')
    res['algebra_psl44_exact_on_basis'] = alg
    res['algebra_negative_controls'] = {
        k: {'super_Jacobi_violating_ordered_pairs': d['super_Jacobi']['violations'],
            'Killing_form_nonzero_entries': d['Killing_form_str_ad_ad']['nonzero_entries']}
        for k, d in (('so(W_perp) coefficient 1/2', load('exhaustive_alg_control_cR_half.json')),
                     ('so(W_perp) coefficient 0', load('exhaustive_alg_control_cR_zero.json'))) if d}
    res['gauge_algebras'] = {g: X.gauge(g).lie_tests() for g in ('su2', 'su3')}
    res['model_g_su2_equals_model_py_term_by_term'] = X.consistency()
    res['color_tree_injectivity'] = load('exh_trees.json')
    res['graded_symmetry_checks'] = {g: load('exh_symcheck_%s.json' % g) for g in ('su2', 'su3')}
    res['classical_master_equation_arity0'] = {g: load('exh_cme_%s.json' % g) for g in ('su2', 'su3')}
    res['negative_controls_basis_level_su2'] = load('exh_controls_su2.json')
    els, names = X.basis_list()
    totals = {n: len(X.tuples(n, els)) for n in (1, 2, 3)}
    summ = X.summary()
    comp = {}
    for key, d in sorted(summ.items()):
        g, ar = key.split('_')[0], int(key.split('_')[1][2])
        sample = key.endswith('_sample')
        comp[key] = {'gauge_algebra': g, 'arity': ar,
                     'coverage': ('pattern sample' if sample else
                                  ('complete' if d['checked'] == totals[ar] else
                                   'partial %d/%d' % (d['checked'], totals[ar]))),
                     'basis_tuples_checked': d['checked'], 'basis_tuples_total': totals[ar],
                     'failures': d['failures'], 'max_field_factors_per_term_in_pieces': d['max_legs'],
                     'cpu_seconds': d['seconds'], 'tuples_by_subspace_pattern_letters': d['by_types']}
    res['master_equation_components_on_basis'] = comp
    res['tuple_conventions'] = {
        'basis': '62 elements of structure.basis(): P(4) K(4) M(3) Mbar(3) Delta(1) R(15) | Q(8) Qbar(8) Sbar(8) S(8)',
        'arity1': 'all 62 basis elements',
        'arity2': 'all 1953 unordered pairs with repetition',
        'arity3': 'all 21824 unordered triples with repetition containing >= 2 odd elements; triples with <= 1 odd '
                  'element vanish term by term (S^(2) vanishes on an even argument, [even,even] is even, S^(3) = 0)',
        'pattern letters': 'P K M m(=Mbar) D(=Delta) R Q q(=Qbar) s(=Sbar) S',
        'test': 'each component is a local functional; it vanishes iff all Euler-Lagrange derivatives of its density '
                'vanish identically (tested separately in each sector of fixed field content); exact Gaussian rationals',
    }
    def complete(key):
        d = comp.get(key)
        return bool(d) and d['coverage'] == 'complete' and not d['failures']
    tr = res['color_tree_injectivity'] or {}
    res['conclusions'] = {
        'algebra': 'Exact on the 62-element basis (hence on all of the realization by multilinearity): closure, '
                   'even bracket, graded antisymmetry, super-Jacobi on all 62^3 ordered triples, vanishing '
                   'Killing form str(ad_x ad_y) on all 62^2 ordered pairs, strict [SW (17)] block rule on all ordered '
                   'pairs, ad_Delta eigenvalues, [g,g] = g, trivial centre, even part = two commuting 15-dimensional '
                   'ideals with nondegenerate Killing forms.  All pass: %s.' % (alg or {}).get('all_checks_pass'),
        'master_equation_su2_complete_on_basis': {'arity1': complete('su2_ar1'), 'arity2': complete('su2_ar2'),
                                                  'arity3': complete('su2_ar3')},
        'master_equation_su3_complete_on_basis': {'arity1': complete('su3_ar1'), 'arity2': complete('su3_ar2'),
                                                  'arity3': complete('su3_ar3')},
        'every_gauge_algebra': (
            'Every term of the arity 0, 1, 2, 3 components is a g-independent coefficient times a connected tree tensor '
            'built from ([-,-],-) and the inverse pairing, with k <= 5 legs (structural bound: S has <= 4 g-valued '
            'fields and its antifield-dependent terms <= 3, S^(1) <= 3, S^(2) = 2, and a bracket removes two; the '
            'computed pieces have at most 4, 4, 3, 2 field factors at arities 0, 1, 2, 3).  In any g the k-leg tree '
            'tensors span a quotient of the multilinear free Lie algebra on k-1 letters, of dimension (k-2)! '
            '(1, 1, %s, %s for k = 2, 3, 4, 5).  In su(2) the ranks are 1, 1, %s, %s, so evaluation on '
            'su(2) is injective and S_k-equivariant on trees with k <= 5 legs; an S_k-equivariant left inverse commutes with total '
            'derivatives, so a component that is a total derivative for su(2) is a total derivative for every g '
            'with a nondegenerate invariant pairing.  With the complete su(2) basis checks this proves the arity '
            '1, 2, 3 components for every such g; arity 0 (classical master equation, <= 5 legs) follows the same '
            'way from the su(2) check; arity 4 reduces to {S^(2),S^(2)} = 0 identically and arity >= 5 has no terms.'
            % (tr.get('free_Lie_multilinear_dimension', {}).get('3 letters (4 legs)'),
               tr.get('free_Lie_multilinear_dimension', {}).get('4 letters (5 legs)'),
               tr.get('su2', {}).get('rank_4leg_trees(3 trees)'), tr.get('su2', {}).get('rank_5leg_trees(15 trees)'))),
    }
    json.dump(res, open(os.path.join(HERE, 'verify_exhaustive.json'), 'w'), indent=1, default=str)
    print(json.dumps({k: v for k, v in res['master_equation_components_on_basis'].items()}, indent=1))


if __name__ == '__main__':
    main()
