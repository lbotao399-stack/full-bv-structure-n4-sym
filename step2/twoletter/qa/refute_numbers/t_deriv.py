"""Derivative letters (total derivative number >= 3) and eps^{-m} parts: referee engine R vs engine A on demand."""
import sys, json, itertools
import numpy as np, sympy as sp
from fractions import Fraction
from ref_core import *
sys.path.insert(0, '../../engine_a/src.nosync')
import ea_kern as EK

r1 = [1, 2, -1, 3]; r2 = [2, -1, 3, 1]
Vc = list(range(4)); Wc = list(range(4, 10))
comps = {'AV': Vc, 'AW': Wc}
def e10(M): return unit(10, M)
def ea_fac(alpha):
    # (i p)^alpha as {exp: coefficient}
    return {tuple(alpha): complex(1j ** sum(alpha))}
def to_sp(fr_re, fr_im):
    return sp.Rational(fr_re.numerator, fr_re.denominator) + sp.I * sp.Rational(fr_im.numerator, fr_im.denominator)
report = []
def case(name, ta, tb, alpha, beta, blocks_R, blocks_A, kinds_list):
    resA = EK.kernel(blocks_A, ta, tb, fa=ea_fac(alpha), fb=ea_fac(beta))
    nbad = 0; ntot = 0; nz = 0
    for ia, A in enumerate(comps[ta]):
        for ib, B in enumerate(comps[tb]):
            verts, rs = blocks_R
            mine = kernel_chain('A', e10(A), 'A', e10(B), verts, rs, fa=deriv_factor(alpha), fb=deriv_factor(beta))
            for kd in kinds_list:
                m = sp.nsimplify(mine.get(kd, 0))
                if resA is None or kd not in resA or resA[kd] is None:
                    t = 0
                else:
                    re, im = resA[kd]
                    t = to_sp(re[ia, ib], im[ia, ib])
                ntot += 1
                nz += (t != 0)
                if sp.simplify(m - t) != 0:
                    nbad += 1
                    if nbad <= 3: print('  mismatch', name, kd, A, B, m, t)
    print(name, 'entries', ntot, 'nonzero', nz, 'mismatches', nbad)
    report.append(dict(case=name, entries=ntot, nonzero=nz, mismatches=nbad))
y = [0,0,0,0,0,1,0,2,0,-1]; y2 = [0,0,0,0,1,0,1,0,3,0]; yv = [1,0,2,-1,0,0,0,0,0,0]
cases = [
 ('Aperp d0d1 | Aperp d2 -> Aperp Aperp, quartic', 'AW', 'AW', (1,1,0,0), (0,0,1,0),
   ([('AAAA', y, y2)], [[r1[m]+r2[m] for m in range(4)]]), [('qua', (np.array(y), np.array(r1)), (np.array(y2), np.array(r2)))], ['log', 1]),
 ('Aperp d0d1 | Aperp d2 -> Aperp Aperp, gluon exchange', 'AW', 'AW', (1,1,0,0), (0,0,1,0),
   ([('AAA', y), ('AAA', y2)], [r1, r2]), [('cub', 'AW', np.array(y), np.array(r1)), ('cub', 'AW', np.array(y2), np.array(r2))], ['log', 1]),
 ('Aperp d0d0 | Aperp d1d3 -> Aperp Aperp, gluon exchange (W=6)', 'AW', 'AW', (2,0,0,0), (0,1,0,1),
   ([('AAA', y), ('AAA', y2)], [r1, r2]), [('cub', 'AW', np.array(y), np.array(r1)), ('cub', 'AW', np.array(y2), np.array(r2))], ['log', 1, 2]),
 ('Aperp d0d0 | Aperp d1d3 -> Aperp Aperp, quartic (W=6)', 'AW', 'AW', (2,0,0,0), (0,1,0,1),
   ([('AAAA', y, y2)], [[r1[m]+r2[m] for m in range(4)]]), [('qua', (np.array(y), np.array(r1)), (np.array(y2), np.array(r2)))], ['log', 1, 2]),
 ('AV d0d0d1 | Aperp -> Aperp (W=5)', 'AV', 'AW', (2,1,0,0), (0,0,0,0),
   ([('AAA', y)], [r1]), [('cub', 'AW', np.array(y), np.array(r1))], ['log', 1, 2]),
 ('AV d2 | AV d0d3 -> AV (W=5)', 'AV', 'AV', (0,0,1,0), (1,0,0,1),
   ([('AAA', yv)], [r1]), [('cub', 'AV', np.array(yv), np.array(r1))], ['log', 1, 2]),
 ('Aperp d1 | AV d0d2 -> AV AV CC (W=5)', 'AW', 'AV', (0,1,0,0), (1,0,1,0),
   ([('AAA', e10(4)), ('AAA', yv)], [r1, r2]), [('cub', 'AW', np.array(e10(4)), np.array(r1)), ('cub', 'AV', np.array(yv), np.array(r2))], ['log', 1]),
 ('AV d1d1 | AV d3 -> AV AV CCC? no: AV AV quartic (W=5)', 'AV', 'AV', (0,2,0,0), (0,0,0,1),
   ([('AAAA', e10(0), yv)], [[r1[m]+r2[m] for m in range(4)]]), [('qua', (np.array(e10(0)), np.array(r1)), (np.array(yv), np.array(r2)))], ['log', 1]),
]
for c in cases:
    case(*c)
# k = 0 with derivatives (Lemma lem:tl-k0)
for alpha, beta in [((1,1,0,0),(1,1,0,0)), ((2,0,0,0),(0,2,0,0)), ((1,2,0,1),(0,0,2,0)), ((3,0,0,0),(1,0,0,0))]:
    bad = 0
    for A in range(10):
        for B in range(10):
            mine = kernel_k0('A', e10(A), 'A', e10(B), fa=deriv_factor(alpha), fb=deriv_factor(beta))
            # Lemma: 2 delta_AB (-1)^{|beta|} (d^{alpha+beta} S_eps)(0); 16pi^2 S_eps series
            g = [a + b for a, b in zip(alpha, beta)]
            n = sum(g)
            if n % 2:
                lem = {}
            else:
                m = n // 2 + 1
                # number of pairings of multi-index g: product over pairing of deltas
                # sum over pairings of prod delta = prod_mu (g_mu-1)!! if all g_mu even else 0
                from math import prod
                def dfact(k):
                    return 1 if k <= 0 else k * dfact(k - 2)
                pr = prod(dfact(x - 1) for x in g) if all(x % 2 == 0 for x in g) else 0
                val = sp.Rational((-1) ** (m + 1) * pr, m * 2 ** (m - 1)) * (-1) ** sum(beta) * 2 * (1 if A == B else 0)
                lem = {m: val} if val != 0 else {}
            mm = {k: sp.nsimplify(v) for k, v in mine.items()}
            if mm != lem:
                bad += 1
                if bad < 3: print('k0 mismatch', alpha, beta, A, B, mm, lem)
    print('k=0', alpha, beta, 'mismatches vs Lemma lem:tl-k0:', bad)
    report.append(dict(case='k0 %s %s vs lemma' % (alpha, beta), mismatches=bad))
json.dump(report, open('out_deriv.json', 'w'), indent=1)
