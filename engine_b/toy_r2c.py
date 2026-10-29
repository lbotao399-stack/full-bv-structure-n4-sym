"""Odd special edge: test STr variants with kernel P + delta K (delta odd, in front of the entry), G-valued Hessians."""
import sympy as sp
import toy_r2 as T0
from grassmann import G
n = T0.n; par = T0.par; P = T0.P; K = T0.K; internal = T0.internal; D = T0.D; I = T0.I; zero = T0.zero
WK = T0.W1_exact(True)
delta = G.odd(99)
def cand(variant):
    Hm = [[D(D(I, internal[j]), internal[i]).subs(zero).drop_odd([0, 1]) for j in range(n)] for i in range(n)]
    M = [[G() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            acc = G()
            for j in range(n):
                if P[i][j] != 0:
                    acc = acc + Hm[j][k] * (P[i][j] * (-1) ** par[j])
                if K[i][j] != 0:
                    s = {'col': (-1) ** par[j], 'row': (-1) ** par[i], 'none': 1, 'neg_col': -(-1) ** par[j]}[variant]
                    if variant == 'deltaright':
                        acc = acc + (G.scalar(K[i][j] * (-1) ** par[j]) * delta) * Hm[j][k]
                    else:
                        acc = acc + delta * Hm[j][k] * (K[i][j] * s)
            M[i][k] = T0.T(acc)
    tot = G(); Mp = M
    for nn in range(1, 4):
        tr = G()
        for i in range(n):
            tr = tr + Mp[i][i] * ((-1) ** par[i])
        tot = tot + tr * (T0.hb / 2 / nn)
        Mp2 = [[G() for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for k in range(n):
                acc = G()
                for j in range(n):
                    acc = acc + Mp[i][j] * M[j][k]
                Mp2[i][k] = T0.T(acc)
        Mp = Mp2
    return T0.T(tot).map(lambda v: sp.expand(v).coeff(T0.hb, 1))
for variant in ('col', 'row', 'none', 'neg_col', 'deltaright'):
    C = cand(variant)
    res = []
    for legs, og, oh in [((0, 1), 2, 0), ((1, 2), 2, 0), ((0, 3), 2, 0), ((0, 1, 2), 3, 0), ((0, 1, 2), 1, 1), ((1, 3), 2, 0), ((2, 3), 2, 0), ((1, 2, 3), 3, 0)]:
        a = T0.coeff_exact(WK, legs, og, oh, delta=True); b = T0.coeff_exact(C, legs, og, oh, delta=True)
        res.append((a, b))
    print(variant, res, flush=True)
