"""Nontrivial W2 instance (two-body, involves Q on A_V outputs): M = tr(zeta A_mu A_I), outputs (c, A_perp) with one
derivative.  With the agreed/verified K(A_V,A_perp)|_{A_perp} = 3i r_mu, K(A_V,A_perp)|_{(A_V,A_perp)} total -1/2,
|_{(A_perp,A_V)} total 3, eq tl-brst at two output fields reduces to (derivation in findings.json):
   tr(zeta A_K c):  identity (checked by the r=2 totals)
   tr(zeta c A_K):  w(q1,q2) + x(q1) delta_IK = (i/2) q1_mu delta_IK
w = K(d_mu c, A_I)|_{(c(q1), A_K(q2))},  x = K(c, A_mu)|_{c(q1)}  (units hbar/16pi^2)."""
import json, sympy as sp
from ref_core import *
q1 = list(sp.symbols('q1_0:4')); q2 = list(sp.symbols('q2_0:4'))
ok = True; rep = {}
for mu in range(4):
    for I, K in [(5, 5), (5, 8)]:
        w = kernel_chain('c', [1], 'A', unit(10, I), [('AcAs', 'As', 'A', 'c', [1]), ('AAA', unit(10, K))], [q1, q2],
                         fa=deriv_factor(tuple(1 if m == mu else 0 for m in range(4)))).get('log', 0)
        x = kernel_chain('c', [1], 'A', unit(10, mu), [('AcAs', 'As', 'A', 'c', [1])], [q1]).get('log', 0)
        lhs = sp.expand(w + x * (1 if I == K else 0))
        rhs = sp.I / 2 * q1[mu] * (1 if I == K else 0)
        good = sp.expand(lhs - rhs) == 0
        ok &= good
        rep['mu=%d I=%d K=%d' % (mu, I, K)] = dict(w=str(sp.expand(w)), x=str(sp.expand(x)), ok=good)
# the two r=2 totals used (from the referee engine, symbolic polarisations)
pol = list(sp.symbols('v0:4')) + [0] * 6
t1 = t2 = 0
for verts, rs in [([('AAA', pol), ('AAA', unit(10, 5))], [q1, q2]), ([('AAAA', pol, unit(10, 5))], [[q1[m] + q2[m] for m in range(4)]])]:
    t1 += kernel_chain('A', unit(10, 2), 'A', unit(10, 5), verts, rs).get('log', 0)
for verts, rs in [([('AAA', unit(10, 5)), ('AAA', pol)], [q1, q2]), ([('AAAA', unit(10, 5), pol)], [[q1[m] + q2[m] for m in range(4)]])]:
    t2 += kernel_chain('A', unit(10, 2), 'A', unit(10, 5), verts, rs).get('log', 0)
rep['K(A_2,A_5)->(A_V,A_5) total'] = str(sp.expand(t1)); rep['K(A_2,A_5)->(A_5,A_V) total'] = str(sp.expand(t2))
ok &= sp.expand(t1 + sp.Rational(1, 2) * pol[2]) == 0 and sp.expand(t2 - 3 * pol[2]) == 0
r = kernel_chain('A', unit(10, 2), 'A', unit(10, 5), [('AAA', unit(10, 5))], [q1]).get('log', 0)
rep['K(A_2,A_5)->A_5'] = str(r); ok &= sp.expand(r - 3 * sp.I * q1[2]) == 0
for k, v in rep.items(): print(k, v)
print('W2 instance (A_V, A_perp) -> (c, A_perp) holds:', ok)
rep['ok'] = ok
json.dump(rep, open('out_ward2b.json', 'w'), indent=1, default=str)
