# Audit (odd_rev_audit_graphs), script 1: leg placement, floating insertions, adjacency, and large-distance power counting.
#
# Input data (all from the manuscript / design, nothing fitted):
#   Propagator P(eps,T) components [ch03 prop:c3-P, ch08 prop:c8-twocomp]: A_V x A_V, A_W x A_W, c x Ast_V, Ast_V x c,
#   and lambda x lambda, which in two-component form is <lambda_alpha lambdabar_betadot> only (lam+ x lam-).
#   Quadratic vertex terms of S^(1)(x) [ch02 eq:c2-S1even, eq:c2-S1Q, eq:c2-S1S; design eq:o-S1xi] and S^(2) [ch02 eq:c2-S2].
#   Tails at T = infinity: D ~ |y|^-2, dD ~ |y|^-3 (A-A edge undifferentiated D, A*-c edge dD, lambda-lambda edge dD).
#   Growth of the generator data: vector field v ~ |y|^dv (P: 0, M and Delta: 1, K: 2), dv ~ |y|^(dv-1);
#   Killing spinor ~ |y|^g (Sigma: 0, Sigma^*: 1).
# Output: which insertions float, which floating insertions can be adjacent, the IR degree of single insertions (in R^4)
# and of adjacent pairs (in R^8) per orientation and after the orientation sum, compared with the design's table
# (step2/theory/engine_plan.json, key large_distance_decay_of_floating_insertions, and design Lemma lem:o-legs).
import json, itertools, os

PROP = {('A_V', 'A_V'), ('A_W', 'A_W'), ('c', 'Ast_V'), ('Ast_V', 'c'), ('lam+', 'lam-'), ('lam-', 'lam+')}
NEVER_ON_P = {'lamst', 'cst', 'Ast_W'}
def on_P(leg):
    return any(leg == a for a, b in PROP)
def joinable(l1, l2):
    return (l1, l2) in PROP

# generator types: (name, kind, growth, quadratic terms as lists of legs)
# even: S^(1)(v) = (L_v A + 1/4 div v A_W, A*) - (v.lam, lam*) - (v.c, c*)
EVEN = {'P': 0, 'M(so(V))': 1, 'Delta': 1, 'K': 2}
GEN = {}
for nm, dv in EVEN.items():
    GEN[nm] = {'parity_of_generator': 'even', 'growth': dv,
               'terms': [('A_V', 'Ast_V'), ('A_W', 'Ast_W'), ('lam+', 'lamst'), ('lam-', 'lamst'), ('c', 'cst')]}
GEN['R(so(W))'] = {'parity_of_generator': 'even', 'growth': None,
                   'terms': [('A_W', 'Ast_W'), ('lam+', 'lamst'), ('lam-', 'lamst')]}
# odd: -(Gamma(eps,lam)|_V, A*_V) - (Gamma(eps,lam)|_W, A*_W) + 1/2 (rho(dA) eps, lam*) [+ cubic, + (rho(A_W) eta, lam*) for eta]
# chirality: xi in S+ (x) W couples to lambdabar (lam-); xibar in S- (x) W* to lambda (lam+);
#            eta in S+ (x) W* (undotted) via epsbar(x) = i x eta sigma to lambda (lam+); etabar in S- (x) W via eps(x) to lam-.
ODD = {'xi(S+W)': ('lam-', 0), 'xibar(S-W*)': ('lam+', 0), 'eta(S+W*)': ('lam+', 1), 'etabar(S-W)': ('lam-', 1)}
for nm, (lg, g) in ODD.items():
    GEN[nm] = {'parity_of_generator': 'odd', 'growth': g,
               'terms': [(lg, 'Ast_V'), (lg, 'Ast_W'), ('A_V', 'lamst'), ('A_W', 'lamst')]}
S2_TERMS = [('lamst', 'lamst'), ('A_V', 'cst'), ('A_W', 'cst')]

res = {}
# 1. floating terms (both legs on P)
flo = {}
for nm, d in GEN.items():
    f = [t for t in d['terms'] if on_P(t[0]) and on_P(t[1])]
    flo[nm] = f
res['floating_terms'] = {k: [' x '.join(t) for t in v] for k, v in flo.items()}
res['S2_floating_terms'] = [' x '.join(t) for t in S2_TERMS if on_P(t[0]) and on_P(t[1])]
res['legs_never_on_P_check'] = {l: (not on_P(l)) for l in sorted(NEVER_ON_P)}

# 2. adjacency of two floating insertions: a P edge between one leg of each
adj = {}
names = [n for n in GEN if flo[n]]
for a, b in itertools.combinations_with_replacement(names, 2):
    edges = []
    for ta in flo[a]:
        for tb in flo[b]:
            for la in ta:
                for lb in tb:
                    if joinable(la, lb):
                        rest_a = [x for x in ta if x != la] or [la]
                        rest_b = [x for x in tb if x != lb] or [lb]
                        edges.append({'edge': la + ' x ' + lb, 'chain_ends': (rest_a[0], rest_b[0])})
    if edges:
        adj[a + ' -- ' + b] = edges
res['adjacent_floating_pairs'] = adj

# 3. IR degrees (position space, T = infinity)
def single_exponent(nm):
    d = GEN[nm]
    if d['parity_of_generator'] == 'even':
        dv = d['growth']
        # (v.A, A*): v.d on the A-A edge (|v| |y|^-3) or dv A (|dv| |y|^-2), times the A*-c edge |y|^-3
        return max(dv - 3 - 3, (dv - 1) - 2 - 3)
    g = d['growth']
    return g - 3 - 3          # lambda-lambda edge dD and A*-c edge dD
single = {}
for nm in names:
    e = single_exponent(nm)
    single[nm] = {'integrand_exponent': e, 'R4_degree': 4 + e,
                  'verdict': 'converges' if 4 + e < 0 else ('log' if 4 + e == 0 else 'power %d' % (4 + e))}
res['single_floating_insertion'] = single

def pair_degree(a, b):
    da, db = GEN[a], GEN[b]
    if da['parity_of_generator'] == 'odd' and db['parity_of_generator'] == 'odd':
        # c --A*c-- y1 --lam lam-- y2 --A*c-- c : three dD edges, spinor growth g1+g2
        e = da['growth'] + db['growth'] - 9
    elif da['parity_of_generator'] == 'even' and db['parity_of_generator'] == 'even':
        # c --A*c-- y1 --A_V A_V-- y2 --A*c-- c : two dD edges, one D edge, each v. costs |y|^(dv-1) on the shared A-A edge
        e = (da['growth'] - 1) + (db['growth'] - 1) - 3 - 3 - 2
    else:
        return None
    deg = 8 + e
    return {'integrand_exponent_R8': e, 'degree_per_orientation': deg, 'degree_after_orientation_sum': deg - 1}
pairs = {}
for key in adj:
    a, b = key.split(' -- ')
    pd = pair_degree(a, b)
    ends = adj[key][0]['chain_ends']
    pd['chain_ends'] = ends
    pd['orientation_sum_applies'] = (ends[0] == ends[1])   # both ends c-type (A*_V legs to ghost ends): relabelling partner
    v = pd['degree_after_orientation_sum'] if pd['orientation_sum_applies'] else pd['degree_per_orientation']
    pd['verdict'] = 'converges' if v < 0 else ('log' if v == 0 else 'power %d' % v)
    pairs[key] = pd
res['adjacent_pairs_IR'] = pairs

# 4. comparison with the design
plan = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', 'theory', 'engine_plan.json')))
tab = plan['large_distance_decay_of_floating_insertions']
placement = plan['graph_class_reference']['placement']
even_even = {k: v for k, v in pairs.items() if all(GEN[x]['parity_of_generator'] == 'even' for x in k.split(' -- '))}
res['comparison_with_design'] = {
    'design_placement_rule': placement,
    'design_table_keys': list(tab.keys()),
    'even_even_pairs_found': sorted(even_even.keys()),
    'even_even_pairs_mentioned_in_design': ('even' in placement and 'adjacent' in placement and 'A_V' in placement and 'two floating even' in placement),
    'mixed_even_odd_adjacent_found': [k for k in pairs if len({GEN[x]['parity_of_generator'] for x in k.split(' -- ')}) == 2],
    'same_chirality_odd_pairs_adjacent': [k for k in pairs if k in ('xi(S+W) -- xi(S+W)', 'eta(S+W*) -- eta(S+W*)')],
}
res['statement'] = ('Floating terms: even generators only (v.A_V, A*_V) [A_V -> c composite line], so(W) none, odd generators only '
                    '-(Gamma(eps,lam)|_V, A*_V) [lambda -> c]; S^(2) never floats. Adjacent pairs: odd-odd through '
                    'lambda x lambdabar (opposite chiralities only) AND even-even through A_V x A_V (omitted in design Lemma lem:o-legs '
                    'and in the engine-plan table); even-odd never. Both pair types have two ghost ends, so the orientation sum '
                    '(Koszul sign -1, script koszul_orientation.py) lowers the degree by one.')
print(json.dumps(res, indent=1, default=str))
json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'census_ir.json'), 'w'), indent=1, default=str)
