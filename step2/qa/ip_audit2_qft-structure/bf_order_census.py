"""Audit: graphs of the radial Gram entry between B_w (three A|_Wperp legs, length 3) and Fbar_w (two lambdabar legs,
length 2), Prop ip-orders and engine_plan_ip.json IP3.

For k cubic and m quartic vertices of I (vertex list of the reduced ESW action, ghosts omitted: no antifield or ghost
legs on B_w, Fbar_w, and a ghost loop needs two more cubic vertices), enumerate all perfect matchings of the half-edges
with the propagator components phi-phi, A_V-A_V, lambda-lambdabar (Prop c8-twocomp), keep the connected ones, and
record (i) the half-edge parity, (ii) the order in sqrt(unit) in the tree-normalised basis, E - V' - 1/2
((-,-)^{V'-E} divided by sqrt(n_B n_F) ~ (-,-)^{-1/2}), (iii) whether every graph needs a self-contraction of two legs
of B_w (colour ([e_a,e_b],e_c) delta_bc = 0) or a tadpole at a cubic vertex (colour f_abb = 0).
"""
import json, itertools

PROP = {frozenset(['phi']), frozenset(['AV']), frozenset(['lam', 'lamb'])}
def can_pair(a, b):
    return frozenset([a, b]) in PROP or (a == b and a in ('phi', 'AV'))

CUBIC = {'AV phi phi': ['AV', 'phi', 'phi'], 'lam lam phi': ['lam', 'lam', 'phi'], 'lamb lamb phi': ['lamb', 'lamb', 'phi'],
         'lamb lam AV': ['lamb', 'lam', 'AV'], 'AV AV AV': ['AV', 'AV', 'AV']}
QUART = {'phi^4': ['phi'] * 4, 'phi phi AV AV': ['phi', 'phi', 'AV', 'AV'], 'AV^4': ['AV'] * 4}

def matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        b = items[i]
        if can_pair(a[1], b[1]):
            rest = items[1:i] + items[i + 1:]
            for m in matchings(rest):
                yield [(a, b)] + m

def connected(nv, edges):
    adj = {v: set() for v in range(nv)}
    for (a, b) in edges:
        adj[a[0]].add(b[0]); adj[b[0]].add(a[0])
    seen, st = {0}, [0]
    while st:
        v = st.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w); st.append(w)
    return len(seen) == nv

res = []
for k in range(0, 3):
    for m in range(0, 2):
        for cub in itertools.combinations_with_replacement(sorted(CUBIC), k):
            for qua in itertools.combinations_with_replacement(sorted(QUART), m):
                verts = [('B', ['phi'] * 3), ('Fbar', ['lamb'] * 2)] + [(c, CUBIC[c]) for c in cub] + [(q, QUART[q]) for q in qua]
                legs = [(vi, t, li) for vi, (_, ts) in enumerate(verts) for li, t in enumerate(ts)]
                nh = len(legs)
                if nh % 2:
                    res.append({'k': k, 'm': m, 'vertices': [v[0] for v in verts], 'half_edges': nh, 'parity': 'odd', 'graphs': 0})
                    continue
                E = nh // 2; Vp = len(verts)
                order = E - Vp - 0.5
                n_conn = 0; n_clean = 0; examples = []
                for mt in matchings(legs):
                    if not connected(Vp, mt):
                        continue
                    n_conn += 1
                    selfB = any(a[0] == 0 and b[0] == 0 for a, b in mt)
                    tadpole_cubic = any(a[0] == b[0] and a[0] >= 2 and len(verts[a[0]][1]) == 3 for a, b in mt)
                    if not selfB and not tadpole_cubic:
                        n_clean += 1
                        if len(examples) < 2:
                            examples.append([(verts[a[0]][0] + '.' + a[1], verts[b[0]][0] + '.' + b[1]) for a, b in mt])
                res.append({'k': k, 'm': m, 'vertices': [v[0] for v in verts], 'half_edges': nh, 'parity': 'even',
                            'order_in_sqrt_unit': order, 'connected_graphs': n_conn,
                            'graphs_without_B_selfcontraction_or_cubic_tadpole': n_clean, 'examples': examples})

summary = {}
for r in res:
    if r['parity'] == 'even' and r.get('connected_graphs', 0) > 0:
        key = 'order %.1f' % r['order_in_sqrt_unit']
        summary.setdefault(key, {'connected': 0, 'clean': 0, 'vertex_sets': []})
        summary[key]['connected'] += r['connected_graphs']
        summary[key]['clean'] += r['graphs_without_B_selfcontraction_or_cubic_tadpole']
        summary[key]['vertex_sets'].append(r['vertices'][2:])
out = {'records': res, 'summary_by_order': summary,
       'statement_checked': 'Prop ip-orders proof: "The half-edge count: 3+2+3+3 and 3+2+4 are odd"; IP3: "no graph at '
                            'relative order unit (odd half-edge counts)"; tp prop:t-census-w3: "no graph with an even number '
                            'of half-edges at relative order unit"',
       'finding': 'the order-1/2 entry (one cubic vertex, 3+2+3 = 8 half-edges, even) has connected graphs; all of them '
                  'contain a self-contraction of two legs of B_w (zero by colour), so the entry vanishes by the topology/colour '
                  'argument of tp prop:t-orderg, not by parity; the odd counts concern order 1, where no graph exists at all'}
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
print(json.dumps(summary, indent=1))
