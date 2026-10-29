"""Fermion kernels at component level: referee engine R run with ENGINE A's Clifford matrices (rho|_Sigma = sb,
rho|_Sigma* = s; checked here) vs engine A's kernel() on the same chains, all 16 x 16 letter components."""
import sys, json, itertools
import numpy as np, sympy as sp
from ref_core import *
sys.path.insert(0, '../../engine_a/src.nosync'); sys.path.insert(0, '../../engine_a/src.nosync/vendor')
import ea_kern as EK
from cliff import s as Sg, sb as Sgb
cl = Cliff(Sgb, Sg)      # sig = rho|_Sigma : Sigma -> Sigma^*, sigb = rho|_Sigma^*
ok = all((cl.sigb[M] * cl.sig[N] + cl.sigb[N] * cl.sig[M] - (2 if M == N else 0) * sp.eye(16)).is_zero_matrix
         for M in range(10) for N in range(10)) and all(cl.sig[M] == cl.sig[M].T for M in range(10))
print('engine A Clifford matrices: Clifford + symmetric (exact):', ok)
r1 = [1, 2, -1, 3]; r2 = [2, -1, 3, 1]
def e(n, i): return unit(n, i)
def ea_fac(alpha): return {tuple(alpha): complex(1j ** sum(alpha))}
def cmp(name, ta, tb, la, lb, verts_of, blocks_of, outs_list, kinds, alpha=(0, 0, 0, 0), beta=(0, 0, 0, 0), rows_a=None, rows_b=None):
    """outs_list: list of output selector tuples; verts_of(ys) -> (verts, rs) for R; blocks_of(ys) -> engine A blocks"""
    nb = nt = nz = 0
    for ys in outs_list:
        resA = EK.kernel(blocks_of(ys), la, lb, fa=ea_fac(alpha), fb=ea_fac(beta))
        verts, rs = verts_of(ys)
        N, Rs = chain_numerator(ta, sp.Matrix(sp.symbols('za0:%d' % SIZE[ta])).T.tolist()[0], tb,
                                sp.Matrix(sp.symbols('zb0:%d' % SIZE[tb])).T.tolist()[0], verts, rs, cl=cl,
                                fa=deriv_factor(alpha), fb=deriv_factor(beta))
        mine = singular_parts(N, Rs)
        za = sp.symbols('za0:%d' % SIZE[ta]); zb = sp.symbols('zb0:%d' % SIZE[tb])
        for kd in kinds:
            Mk = sp.Poly(mine.get(kd, 0), *(za + zb)) if mine.get(kd, 0) != 0 else None
            for ia, A in enumerate(rows_a):
                for ib, B in enumerate(rows_b):
                    m = Mk.coeff_monomial(za[A] * zb[B]) if Mk is not None else 0
                    if resA is None or resA.get(kd) is None:
                        t = 0
                    else:
                        re, im = resA[kd]
                        t = sp.Rational(re[ia, ib].numerator, re[ia, ib].denominator) + sp.I * sp.Rational(im[ia, ib].numerator, im[ia, ib].denominator)
                    nt += 1; nz += (t != 0)
                    if sp.simplify(m - t) != 0:
                        nb += 1
                        if nb <= 3: print('  mismatch', name, kd, ys, A, B, m, t)
    print('%-58s entries %5d nonzero %4d mismatches %d' % (name, nt, nz, nb))
    return dict(case=name, entries=nt, nonzero=nz, mismatches=nb)
rep = []
L16 = list(range(16)); A10 = list(range(10))
pol16 = lambda i: e(16, i)
# (a) lambda | lambda -> A  (k = 1), log and eps^-1, all 10 output directions
rep.append(cmp('lam|lam -> A (k=1), W=3', 'l', 'l', 'lamP', 'lamM',
               lambda ys: ([('lAl', 'l', 'l', 'A', ys[0])], [r1]),
               lambda ys: [('cub', 'AV' if ys[0].index(1) < 4 else 'AW', np.array(ys[0]), np.array(r1))],
               [(e(10, M),) for M in A10], ['log', 1], rows_a=L16, rows_b=L16))
# (b) with one derivative on letter a: W = 4  (log, eps^-1)
rep.append(cmp('d1 lam|lam -> A (k=1), W=4', 'l', 'l', 'lamP', 'lamM',
               lambda ys: ([('lAl', 'l', 'l', 'A', ys[0])], [r1]),
               lambda ys: [('cub', 'AV' if ys[0].index(1) < 4 else 'AW', np.array(ys[0]), np.array(r1))],
               [(e(10, M),) for M in (0, 3, 5, 9)], ['log', 1], alpha=(0, 1, 0, 0), rows_a=L16, rows_b=L16))
# (c) d0 d2 lam | d3 lam -> A : W = 6 (log, eps^-1, eps^-2) -- total derivative number 3, fermion pair
rep.append(cmp('d0d2 lam|d3 lam -> A (k=1), W=6', 'l', 'l', 'lamP', 'lamM',
               lambda ys: ([('lAl', 'l', 'l', 'A', ys[0])], [r1]),
               lambda ys: [('cub', 'AV' if ys[0].index(1) < 4 else 'AW', np.array(ys[0]), np.array(r1))],
               [(e(10, M),) for M in (1, 6)], ['log', 1, 2], alpha=(1, 0, 1, 0), beta=(0, 0, 0, 1), rows_a=L16, rows_b=L16))
# (d) lambda | lambda -> lambda lambda (gluon exchange), log
rep.append(cmp('lam|lam -> lam lam (gluon exch., k=2), W=3', 'l', 'l', 'lamP', 'lamM',
               lambda ys: ([('lAl', 'l', 'A', 'l', ys[0]), ('lAl', 'A', 'l', 'l', ys[1])], [r1, r2]),
               lambda ys: [('cub', 'lamP', np.array(ys[0]), np.array(r1)), ('cub', 'lamP', np.array(ys[1]), np.array(r2))],
               [(pol16(i), pol16(j)) for i, j in [(0, 3), (5, 5), (7, 12), (15, 2)]], ['log'], rows_a=L16, rows_b=L16))
# (e) lambda | lambda -> A A (k=2), with a derivative: W = 4  (log, eps^-1)
rep.append(cmp('d2 lam|lam -> A A (k=2), W=4', 'l', 'l', 'lamP', 'lamM',
               lambda ys: ([('lAl', 'l', 'l', 'A', ys[0]), ('lAl', 'l', 'l', 'A', ys[1])], [r1, r2]),
               lambda ys: [('cub', 'AV' if ys[0].index(1) < 4 else 'AW', np.array(ys[0]), np.array(r1)),
                           ('cub', 'AV' if ys[1].index(1) < 4 else 'AW', np.array(ys[1]), np.array(r2))],
               [(e(10, M), e(10, N)) for M, N in [(0, 1), (4, 4), (2, 7), (9, 8)]], ['log', 1], alpha=(0, 0, 1, 0), rows_a=L16, rows_b=L16))
# (f) A | lambda -> lambda (k=1, Yukawa/gauge), with derivatives d0 on A, d1d1 on lambda: W = 11/2
rep.append(cmp('d0 A|d1d1 lam -> lam (k=1), W=11/2', 'A', 'l', 'AV', 'lamP',
               lambda ys: ([('lAl', 'A', 'l', 'l', ys[0])], [r1]),
               lambda ys: [('cub', 'lamP', np.array(ys[0]), np.array(r1))],
               [(pol16(i),) for i in (0, 6, 11)], ['log', 1, 2], alpha=(1, 0, 0, 0), beta=(0, 2, 0, 0), rows_a=list(range(4)), rows_b=L16))
rep.append(cmp('A_perp|lam -> lam (k=1, Yukawa), W=5/2', 'A', 'l', 'AW', 'lamP',
               lambda ys: ([('lAl', 'A', 'l', 'l', ys[0])], [r1]),
               lambda ys: [('cub', 'lamP', np.array(ys[0]), np.array(r1))],
               [(pol16(i),) for i in range(16)], ['log'], rows_a=list(range(4, 10)), rows_b=L16))
# (g) c | lambda -> c lambda (k=2, ghost + fermion, odd outputs)
rep.append(cmp('c|lam -> c lam (k=2), W=3/2', 'c', 'l', 'c', 'lamP',
               lambda ys: ([('AcAs', 'As', 'A', 'c', ys[0]), ('lAl', 'A', 'l', 'l', ys[1])], [r1, r2]),
               lambda ys: [('cub', 'c', np.array(ys[0]), np.array(r1)), ('cub', 'lamP', np.array(ys[1]), np.array(r2))],
               [([1], pol16(i)) for i in range(16)], ['log'], rows_a=[0], rows_b=L16))
json.dump(rep, open('out_ferm_engA.json', 'w'), indent=1)
