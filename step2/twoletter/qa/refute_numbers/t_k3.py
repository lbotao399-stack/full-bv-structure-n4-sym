"""k = 3 (E = 4) chains with derivative letters, log and eps^{-m}: referee R vs engine A on demand."""
import sys, json
import numpy as np, sympy as sp
from ref_core import *
sys.path.insert(0, '../../engine_a/src.nosync')
import ea_kern as EK
r1 = [1, 2, -1, 3]; r2 = [2, -1, 3, 1]; r3 = [-1, 1, 2, 2]
def ea_fac(alpha): return {tuple(alpha): complex(1j ** sum(alpha))}
rep = []
def run(name, la, lb, rows_a, rows_b, ys, types, alpha, beta, kinds):
    blocks = [('cub', t, np.array(y), np.array(r)) for t, y, r in zip(types, ys, [r1, r2, r3])]
    resA = EK.kernel(blocks, la, lb, fa=ea_fac(alpha), fb=ea_fac(beta))
    nb = nt = nz = 0
    for ia, A in enumerate(rows_a):
        for ib, B in enumerate(rows_b):
            mine = kernel_chain('A', unit(10, A), 'A', unit(10, B), [('AAA', y) for y in ys], [r1, r2, r3],
                                fa=deriv_factor(alpha), fb=deriv_factor(beta))
            for kd in kinds:
                m = mine.get(kd, 0)
                if resA is None or resA.get(kd) is None: t = 0
                else:
                    re, im = resA[kd]
                    t = sp.Rational(re[ia, ib].numerator, re[ia, ib].denominator) + sp.I * sp.Rational(im[ia, ib].numerator, im[ia, ib].denominator)
                nt += 1; nz += (t != 0)
                if sp.simplify(m - t) != 0:
                    nb += 1
                    if nb < 3: print(' mismatch', name, kd, A, B, m, t)
    print('%-50s entries %d nonzero %d mismatches %d' % (name, nt, nz, nb))
    rep.append(dict(case=name, entries=nt, nonzero=nz, mismatches=nb))
yW1 = [0] * 4 + [1, 0, 2, 0, 0, -1]; yW2 = [0] * 4 + [0, 1, 0, 1, 3, 0]; yV = [1, 2, 0, -1] + [0] * 6
run('k=3 Aperp d0d1 | Aperp d3 -> AV Aperp Aperp, W=5', 'AW', 'AW', list(range(4, 10)), list(range(4, 10)),
    [yV, yW1, yW2], ['AV', 'AW', 'AW'], (1, 1, 0, 0), (0, 0, 0, 1), ['log', 1])
run('k=3 AV d2d2 | AV d0d1 -> AV AV AV, W=6', 'AV', 'AV', list(range(4)), list(range(4)),
    [yV, [0, 1, 1, 0] + [0] * 6, [2, 0, 0, 1] + [0] * 6], ['AV', 'AV', 'AV'], (0, 0, 2, 0), (1, 1, 0, 0), ['log', 1, 2])
json.dump(rep, open('out_k3.json', 'w'), indent=1)
