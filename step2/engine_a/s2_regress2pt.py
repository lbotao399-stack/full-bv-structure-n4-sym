"""Regression of the special-edge engine (s2_graphs) against the Step-1 sector (c, A|_V) of the quadratic obstruction
(engine_a/logs/r2_T2_routeA.json and manuscript ch04 table: eps^{-1} coefficients -9/4, 1/8, -9/2, 2, 1/8 of
eps^{-1} int(dc,A|_V); finite coefficients -5/8, -1/12, -5/4, -4/3, -1/12 of int(d^*dc,d^*A|_V)).
Plane-wave values (legs c: gamma e^{ipx}, A: a e^{-ipx}): int(dc,A|_V) -> i (gamma, iota_p a);
int(d^*dc, d^*A|_V) -> i p^2 (gamma, iota_p a)."""
import json, sys, time
from fractions import Fraction
import numpy as np
import s2_core as C
from r2_core import Q2
import s2_graphs as SG

rng = np.random.default_rng(5)
t0 = time.time()
res = {}
for it in range(2):
    p = rng.integers(-2, 3, size=4)
    if not p.any(): p[0] = 1
    gam = C.traceless(rng.integers(-2, 3, size=C.D).astype(complex))
    a = np.zeros((10, C.D), complex); a[:4] = rng.integers(-2, 3, size=(4, C.D)); a = C.traceless(a)
    Sv = complex(np.einsum('a,ab,m,mb->', gam, C.G, p, a[:4]))
    if abs(Sv) == 0:
        continue
    G = SG.SpecialCycles([('c', None, gam, p), ('A', 'V', a, -p)], degree=3)
    for (seq, Ls, sgn) in G.enumerate():
        v = G.graph_value(seq, Ls, sgn)
        key = (Ls[0], Ls[1])
        for e in (-1, 0):
            b = SG.to_basis(v[e])
            norm = Q2(0, Sv.real) * C.C_ADJ * (1 if e == -1 else int(p @ p))   # i S c_adj (p^2)
            ent = res.setdefault(key, {}).setdefault(e, {})
            ent.setdefault(it, [Q2(0), Q2(0), Q2(0)])
            ent[it] = [x + y / norm for x, y in zip(ent[it], b)]
out = {}
for key, d in res.items():
    o = {}
    for e, per in d.items():
        vals = list(per.values())
        assert all(all((x - y).iszero() for x, y in zip(vals[0], w)) for w in vals), (key, e, vals)
        o['eps^%d' % e] = '%s + (%s) log2 + (%s) log3' % tuple(map(str, vals[0]))
    out['%s|%s' % key] = o
    print(key, o, flush=True)
tot = {e: [Q2(0)] * 3 for e in (-1, 0)}
for key, d in res.items():
    for e in (-1, 0):
        v = list(d[e].values())[0]
        tot[e] = [x + y for x, y in zip(tot[e], v)]
out['TOTAL'] = {'eps^%d' % e: '%s + (%s) log2 + (%s) log3' % tuple(map(str, v)) for e, v in tot.items()}
out['step1_expected'] = {'sAV|AV': (-9 / 4, -5 / 8), 'sAV|gh': (1 / 8, -1 / 12), 'sAP|AP': (-9 / 2, -5 / 4),
                         'sl|l': (2, -4 / 3), 'sc|gh': (1 / 8, -1 / 12), 'total': (-9 / 2, -27 / 8)}
out['seconds'] = round(time.time() - t0, 1)
out['algebra'] = C.ALG
print(json.dumps(out['TOTAL']), out['seconds'], flush=True)
json.dump(out, open('logs/s2_regress2pt_%s.json' % C.ALG, 'w'), indent=1)
