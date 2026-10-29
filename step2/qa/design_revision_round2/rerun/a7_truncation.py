"""Audit script a7: weight/degree census behind design lem:o-actiontrunc and prop:o-trunc.
Fields (degree, Weyl weight): c (+1,0), A (0,1), lambda (0,3/2), A* (-1,3), lambda* (-1,5/2), c* (-2,4);
a derivative adds (0,1); an explicit power of the position x adds (0,-1) (x-dependent cochains,
e.g. S^(1)(Delta), S^(1)(K_b), the inner-Delta coboundary, eta-components).
Every field is in the adjoint of g, so a colour singlet needs at least two field factors.
"""
import json, os, itertools
from fractions import Fraction as Fr

F = {'c': (1, Fr(0)), 'A': (0, Fr(1)), 'lam': (0, Fr(3, 2)), 'A*': (-1, Fr(3)), 'lam*': (-1, Fr(5, 2)), 'c*': (-2, Fr(4))}
names = list(F)


def contents(maxfields=6, maxder=3, maxx=0):
    out = []
    for nf in range(1, maxfields + 1):
        for combo in itertools.combinations_with_replacement(names, nf):
            d = sum(F[f][0] for f in combo); w = sum(F[f][1] for f in combo)
            for nd in range(0, maxder + 1):
                for nx in range(0, maxx + 1):
                    out.append((combo, nd, nx, d, w + nd - nx))
    return out


res = {}
C0 = contents(maxx=0)
C1 = contents(maxx=2)
m = {}
for n in range(1, 6):
    ws = [w for (cb, nd, nx, d, w) in C0 if d == -n]
    m[n] = str(min(ws)) if ws else None
res['m(n)_translation_invariant'] = m
res['design_m(n)'] = {n: str(Fr(2 * n) if n % 2 == 0 else Fr(2 * n) + Fr(1, 2)) for n in range(1, 6)}
res['m(n)_agrees'] = all(m[n] == res['design_m(n)'][n] for n in range(1, 6))
# action side: arity k on Sigma^k: degree -k, weight 4 + k/2
act = {}
for k in range(1, 6):
    hit0 = [(cb, nd, nx) for (cb, nd, nx, d, w) in C0 if d == -k and w == 4 + Fr(k, 2)]
    hit1 = [(cb, nd, nx) for (cb, nd, nx, d, w) in C1 if d == -k and w == 4 + Fr(k, 2)]
    act[k] = {'translation_invariant_contents': [str(h) for h in hit0[:6]], 'n_ti': len(hit0),
              'with_explicit_x_contents': [str(h) for h in hit1 if h[2] > 0][:6], 'n_x': len([h for h in hit1 if h[2] > 0])}
res['action_side_arity_k_on_Sigma'] = act

# prop:o-trunc targets: (w_O, |O|), arity k >= 3: weight w_O + k/2, degree |O| + 1 - k
targets = {'letters A (1,0)': (Fr(1), 0), 'letter lambda (3/2,0)': (Fr(3, 2), 0), 'ghost letter (0,1)': (Fr(0), 1),
           'scalar bilinears (2,0)': (Fr(2), 0), 'Konishi descendant O_zeta (5/2,0)': (Fr(5, 2), 0),
           'bc entry (2,1)': (Fr(2), 1), 'weight-3 B_w, Fbar_w (3,0)': (Fr(3), 0)}
tr = {}
for name, (wO, dO) in targets.items():
    row = {}
    for k in range(3, 7):
        hits = [(cb, nd) for (cb, nd, nx, d, w) in C0 if d == dO + 1 - k and w == wO + Fr(k, 2)]
        row[k] = {'contents': [str(h) for h in hits], 'colour_singlet_capable(>=2 fields)': [str(h) for h in hits if len(h[0]) >= 2]}
    tr[name] = row
res['prop_o_trunc_census'] = tr
# x-dependent primitives relevant to prop:c8-B / prop:o-GB (ii): degree -2, weight 9/2
res['degree-2_weight-9/2_contents_translation_invariant'] = [str((cb, nd)) for (cb, nd, nx, d, w) in C0 if d == -2 and w == Fr(9, 2)]
res['degree-2_weight-9/2_contents_with_x'] = [str((cb, nd, nx)) for (cb, nd, nx, d, w) in C1 if d == -2 and w == Fr(9, 2) and nx > 0][:8]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a7_truncation.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
