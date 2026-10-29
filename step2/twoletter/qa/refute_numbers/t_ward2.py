"""Ward identity W2 (eq tl-brst at two letters) for M = tr(zeta A_I A_J), planar, outputs (A_perp, A_perp, c).
With s A_I = -[A_I, c] (common overall sign), lc(sM) = s lc(M) + {L,M} reduces (derivation in findings) to:
planar K(c, A_perp)|_{(c, A_perp)} = planar K(A_perp, c)|_{(A_perp, c)} = +1 unit (identity); and the A_V outputs of
K(A_perp, A_perp) vanish in total (r = 1: graph (b); r = 2: triangle + bubble)."""
import json, sympy as sp
from ref_core import *
r1 = list(sp.symbols('r1_0:4')); r2 = list(sp.symbols('r2_0:4'))
out = {}
for I in (4, 7):
    for J in (4, 7):
        a = unit(10, I); b = unit(10, J)
        # K(c, A_J): outputs (c, A) in chain order; select output A component K
        vals = {}
        for K in (4, 7):
            r = kernel_chain('c', [1], 'A', b, [('AcAs', 'As', 'A', 'c', [1]), ('AAA', unit(10, K))], [r1, r2])
            vals['K(c,A_%d)->(c,A_%d)' % (J, K)] = str(sp.expand(r.get('log', 0)))
            r = kernel_chain('A', a, 'c', [1], [('AAA', unit(10, K)), ('AcAs', 'A', 'As', 'c', [1])], [r1, r2])
            vals['K(A_%d,c)->(A_%d,c)' % (I, K)] = str(sp.expand(r.get('log', 0)))
        out.update(vals)
# A_V outputs of K(A_perp, A_perp)
yv = list(sp.symbols('y0:4')) + [0] * 6
r = kernel_chain('A', unit(10, 5), 'A', unit(10, 5), [('AAA', yv)], [r1])
out['K(A5,A5)->A_V (r=1)'] = str(sp.expand(r.get('log', 0)))
yv2 = list(sp.symbols('z0:4')) + [0] * 6
tot = 0
for verts, rs in [([('AAA', yv), ('AAA', yv2)], [r1, r2]), ([('AAAA', yv, yv2)], [[r1[m] + r2[m] for m in range(4)]])]:
    tot += kernel_chain('A', unit(10, 5), 'A', unit(10, 6), verts, rs).get('log', 0)
out['K(A5,A6)->A_V A_V total'] = str(sp.expand(tot))
for k, v in out.items():
    print(k, ':', v)
ok = all(out['K(c,A_%d)->(c,A_%d)' % (J, K)] == ('2' if J == K else '0') for J in (4, 7) for K in (4, 7)) and \
     all(out['K(A_%d,c)->(A_%d,c)' % (I, K)] == ('2' if I == K else '0') for I in (4, 7) for K in (4, 7)) and \
     out['K(A5,A5)->A_V (r=1)'] == '0' and out['K(A5,A6)->A_V A_V total'] == '0'
print('W2 instance (A_perp, A_perp) consistent:', ok)
out['ok'] = ok
json.dump(out, open('out_ward2.json', 'w'), indent=1)
