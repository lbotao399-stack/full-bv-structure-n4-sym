"""Theory Example after Prop prop:tl-bubble: quartic bubble with one derivative on each letter:
K(d_mu X_a, d_nu X_b) = -(W_0 value) * (R_mu R_nu/6 + delta_mu nu |R|^2/12)  (16 pi^2 [log eps])."""
import json, sympy as sp
from ref_core import *
r1 = list(sp.symbols('r1_0:4')); r2 = list(sp.symbols('r2_0:4'))
R = [r1[m] + r2[m] for m in range(4)]
y1 = [0] * 4 + [1, 0, 0, 0, 0, 0]; y2 = [0] * 4 + [0, 1, 0, 0, 0, 0]
ok = True
for (A, B) in [(4, 5), (5, 4), (6, 6)]:
    w0 = kernel_chain('A', unit(10, A), 'A', unit(10, B), [('AAAA', y1, y2)], [R]).get('log', 0)
    for mu in range(4):
        for nu in range(4):
            al = tuple(1 if m == mu else 0 for m in range(4)); be = tuple(1 if m == nu else 0 for m in range(4))
            v = kernel_chain('A', unit(10, A), 'A', unit(10, B), [('AAAA', y1, y2)], [R], fa=deriv_factor(al), fb=deriv_factor(be)).get('log', 0)
            pred = -w0 * (R[mu] * R[nu] / 6 + (1 if mu == nu else 0) * sum(x * x for x in R) / 12)
            ok &= sp.expand(v - pred) == 0
print('bubble example (one derivative on each letter):', ok)
json.dump(dict(ok=ok), open('out_bubble_example.json', 'w'))
