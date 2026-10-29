"""Colour reduction for engine B: every one-loop colour factor of the two-point and O-insertion graphs is reduced to a
rational multiple of c_adj * (c1,c2), with c_adj the scalar by which X -> sum_a [t_a,[t^a,X]] acts (t^a dual w.r.t. (-,-)).
Checked on su(2) and su(3) in non-orthonormal bases with rescaled pairings."""
import itertools, sympy as sp
from realization import su2_data, su3_data

def colour_matrices(L, c):
    n = L['n']; ft = L['ft']
    # (F_c)_{de} = ftil(c, t_d, t_e)
    return sp.Matrix(n, n, lambda d, e: sum(c[a] * ft[a, d, e] for a in range(n)))

def bracket(L, x, y):
    n = L['n']; f = L['f']
    return [sum(x[a] * y[b] * f[a][b][c] for a in range(n) for b in range(n)) for c in range(n)]

def pair(L, x, y):
    return sum(x[a] * L['kap'][a, b] * y[b] for a in range(L['n']) for b in range(L['n']))

def factors(L, c1, c2):
    n = L['n']; kinv = L['kinv']; e = lambda d: [1 if i == d else 0 for i in range(n)]
    F1 = colour_matrices(L, c1); F2 = colour_matrices(L, c2)
    out = {}
    out['bubble'] = (kinv * F1 * kinv * F2).trace()
    # quartic tadpoles: contraction of legs 3,4 with (-,-)^{-1}
    G = {}
    for key in ('K12', 'K13', 'K14'):
        M = sp.zeros(n, n)
        for d in range(n):
            for ee in range(n):
                if key == 'K12':
                    M[d, ee] = pair(L, bracket(L, c1, c2), bracket(L, e(d), e(ee)))
                elif key == 'K13':
                    M[d, ee] = pair(L, bracket(L, c1, e(d)), bracket(L, c2, e(ee)))
                else:
                    M[d, ee] = pair(L, bracket(L, c1, e(ee)), bracket(L, c2, e(d)))
        G[key] = M
        out['tad_' + key] = (kinv * M).trace()
    out['pair'] = pair(L, c1, c2)
    return out

def reduce_all():
    res = {}
    for name, L in (('su2', su2_data()), ('su3', su3_data())):
        n = L['n']
        ratios = {}
        for (i, j) in [(0, 0), (0, 1), (1, 2), (2, 2), (0, n - 1)]:
            c1 = [1 if a == i else 0 for a in range(n)]; c2 = [1 if a == j else 0 for a in range(n)]
            F = factors(L, c1, c2)
            if F['pair'] == 0:
                for k, v in F.items():
                    if k != 'pair':
                        assert sp.simplify(v) == 0, (name, i, j, k, v)
                continue
            for k, v in F.items():
                if k == 'pair':
                    continue
                r = sp.nsimplify(sp.simplify(v / (L['cadj'] * F['pair'])))
                ratios.setdefault(k, set()).add(r)
        res[name] = {k: sorted(v) for k, v in ratios.items()}
        res[name]['dim'] = n
    return res

if __name__ == '__main__':
    print(reduce_all())
