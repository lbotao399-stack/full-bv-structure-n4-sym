"""Ward identity W1 (Prop prop:tl-W1) for the letters lambda and c, full colour (no planar reading), from the referee's
own kernels: for M = int omega (Z, [c, lambda]) eta, lc(M) = (kappa/2) u (Z,[c,lambda]) with kappa the kernel coefficient
(units hbar/16pi^2, colour C_{ad_{Y2} ad_{Y1}}, outputs (c, lambda)); W1 <=> kernel = 2 x identity on the spinor.
For c: M = int omega (Z,[c,c]); W1 <=> K(c,c)|_{cc} = 2.  For A_V: one-field part K(A_mu,c)|_c = i r_mu."""
import json, sympy as sp
from ref_core import *
from ref_cliff import SIG, SIGB
cl = Cliff(SIG, SIGB)
r1 = [1, 2, -1, 3]; r2 = [2, -1, 3, 1]
eta = list(sp.symbols('eta0:16')); x = list(sp.symbols('x0:16'))
res = kernel_chain('c', [1], 'l', eta, [('AcAs', 'As', 'A', 'c', [1]), ('lAl', 'A', 'l', 'l', x)], [r1, r2], cl=cl)
lg = sp.expand(res.get('log', 0))
target = sp.expand(2 * sum(eta[i] * x[i] for i in range(16)))
print('K(c, lambda)|_{c lambda} log:', 'equals 2 eta.x' if sp.expand(lg - target) == 0 else lg, '; other keys', [k for k in res if k != 'log'])
res2 = kernel_chain('c', [1], 'c', [1], [('AcAs', 'As', 'A', 'c', [1]), ('AcAs', 'A', 'As', 'c', [1])], [r1, r2])
print('K(c,c)|_{cc} log:', res2.get('log'))
rs = list(sp.symbols('r0:4'))
out = {}
for mu in range(4):
    r = kernel_chain('A', unit(10, mu), 'c', [1], [('AcAs', 'A', 'As', 'c', [1])], [rs])
    out[mu] = sp.expand(r.get('log', 0))
print('K(A_mu, c)|_c log:', out, ' (expected i r_mu)')
ok = (sp.expand(lg - target) == 0) and res2.get('log') == 2 and all(sp.expand(out[m] - sp.I * rs[m]) == 0 for m in range(4))
json.dump(dict(W1_lambda=str(sp.expand(lg - target) == 0), K_cc=str(res2.get('log')),
               K_Amu_c={str(k): str(v) for k, v in out.items()}, all_ok=ok), open('out_ward1.json', 'w'), indent=1)
print('W1 lambda, c, A_V (one-field part) all consistent:', ok)
