"""Exact check that the Spin(6)-odd external sets have zero one-loop log-eps (and eps^{-1}) coefficients."""
import json, sys
import numpy as np
import r2_core as C
from r2_core import Q2
from r2_graphs import GraphSum
rng = np.random.default_rng(5)
def rpol(t, which=None):
    x = rng.integers(-1, 2, size=C.SHAPE[t]).astype(complex)
    if which == 'V': x[4:] = 0
    if which == 'P': x[:4] = 0
    return x
SETS = {'AV,AV,AP': ([('A', 'V'), ('A', 'V'), ('A', 'P')], 1), 'AP,AP,AP': ([('A', 'P')] * 3, 1),
        'AV,AV,AV,AP': ([('A', 'V')] * 3 + [('A', 'P')], 0), 'AV,AP,AP,AP': ([('A', 'V')] + [('A', 'P')] * 3, 0),
        'c,AV,AsP': ([('c', None), ('A', 'V'), ('As', 'P')], 0), 'c,AP,AsV': ([('c', None), ('A', 'P'), ('As', 'V')], 0)}
out = {}
for nm, (legs, om) in SETS.items():
    vals = []
    for it in range(3):
        ps = [rng.integers(-2, 3, size=4) * 6 for _ in range(len(legs) - 1)]; ps.append(-sum(ps))
        L = [(t, w, rpol(t, w), p) for (t, w), p in zip(legs, ps)]
        res = GraphSum(L).coefficients(om)
        tot = Q2(0)
        for v in res.values():
            tot = tot + v['log']
            assert v['epsinv'].iszero()
        vals.append(str(tot))
    out[nm] = vals
    print(nm, vals, flush=True)
json.dump(out, open('logs/r2_T1_vanishing.json', 'w'), indent=1)
