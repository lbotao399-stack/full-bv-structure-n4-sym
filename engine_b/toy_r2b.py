"""Which variant of the one-loop supertrace formula is exact at n = 3 (triangles)?  G-valued Hessians, no stripping."""
import itertools, sympy as sp
import toy_r2 as T0
from grassmann import G
n = T0.n; par = T0.par; P = T0.P; internal = T0.internal; D = T0.D; I = T0.I; zero = T0.zero
W = T0.W1_exact(False)
def cand(S_par, order, strace, Sside):
    Hm = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            Hm[i][j] = (D(D(I, internal[j]), internal[i]) if order == 'ij' else D(D(I, internal[i]), internal[j])).subs(zero).drop_odd([0, 1])
    M = [[G() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            acc = G()
            for j in range(n):
                if P[i][j] != 0:
                    s = ((-1) ** par[j] if Sside == 'mid' else (-1) ** par[k]) if S_par else 1
                    acc = acc + Hm[j][k] * (P[i][j] * s)
            M[i][k] = T0.T(acc)
    tot = G(); Mp = M
    for nn in range(1, 4):
        tr = G()
        for i in range(n):
            tr = tr + Mp[i][i] * ((-1) ** par[i] if strace else 1)
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
for S_par in (True, False):
    for order in ('ij', 'ji'):
        for strace in (True, False):
            for Sside in ('mid', 'right'):
                C = cand(S_par, order, strace, Sside)
                res = []
                for legs, og, oh in [((0, 1, 2), 3, 0), ((1, 2, 3), 3, 0), ((0, 1, 2), 1, 1), ((1, 2), 2, 0), ((0, 1), 2, 0)]:
                    res.append(T0.coeff_exact(W, legs, og, oh) == T0.coeff_exact(C, legs, og, oh))
                print(S_par, order, strace, Sside, res)
