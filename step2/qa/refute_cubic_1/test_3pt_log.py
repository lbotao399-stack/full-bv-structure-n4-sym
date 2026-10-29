"""Independent recomputation of the log(eps) coefficients of the 1PI one-loop three-point functions that
enter the two target sectors, compared with the Step-1 claim (thm c4-log: multiples 6 (AAA), 1 (A c A^*),
1 (c c c^*) of the tree vertices)."""
import json
import numpy as np
import ft
import loop

rng = np.random.default_rng(99)
quad = loop.Quad()
EPS = 2e-4 * 0.5 ** np.arange(9)


def rv(n):
    return rng.normal(size=n) + 1j * rng.normal(size=n)


def leg(slot, comps, p, mono):
    if slot in ('A', 'As'):
        f = np.zeros((10, 3), complex)
        for M in comps:
            f[M] = rv(3)
        return ft.Leg(slot, f.reshape(30), p, mono)
    return ft.Leg(slot, rv(ft.DIM[slot]), p, mono)


def fit(vals):
    B = np.stack([np.ones_like(EPS), np.log(EPS), EPS * np.log(EPS), EPS, EPS ** 2 * np.log(EPS), EPS ** 2], 1)
    return np.linalg.lstsq(B, vals, rcond=None)[0]


p1, p2 = rng.normal(size=4), rng.normal(size=4)
p3 = -p1 - p2
V, Wp = range(4), range(4, 10)
cases = {
    '(A_V,A_W,A_W) vs 6 AAA': ([leg('A', V, p1, ()), leg('A', Wp, p2, ()), leg('A', Wp, p3, ())], 0, 6.0),
    '(A_V,A_V,A_V) vs 6 AAA': ([leg('A', V, p1, ()), leg('A', V, p2, ()), leg('A', V, p3, ())], 0, 6.0),
    '(c,A_W,A*_W) vs 1 AcA*': ([leg('c', None, p1, (1,)), leg('A', Wp, p2, ()), leg('As', Wp, p3, (2,))], 2, 1.0),
    '(c,A_V,A*_V) vs 1 AcA*': ([leg('c', None, p1, (1,)), leg('A', V, p2, ()), leg('As', V, p3, (2,))], 2, 1.0),
    '(c,c,c*) vs 1 ccc*': ([leg('c', None, p1, (1,)), leg('c', None, p2, (2,)), leg('cs', None, p3, ())], 4, 1.0),
    '(l,l,A_V) vs 10 llA': ([leg('l', None, p1, (1,)), leg('l', None, p2, (2,)), leg('A', V, p3, ())], 1, 10.0),
    '(l,l,A_W) vs 10 llA': ([leg('l', None, p1, (1,)), leg('l', None, p2, (2,)), leg('A', Wp, p3, ())], 1, 10.0),
    '(c,l,l*) vs 1 cll*': ([leg('c', None, p1, (1,)), leg('l', None, p2, (2,)), leg('ls', None, p3, ())], 3, 1.0),
}
out = {}
for name, (legs, iv, mult) in cases.items():
    ol = loop.OneLoop(legs, rng)
    vals = {}
    for e in EPS:
        for m, v in ol.value(e, quad).items():
            vals.setdefault(m, []).append(v)
    tree = ft.multilinear([ft.CUBIC[iv]], legs)
    res = {}
    for m, vs in vals.items():
        s = fit(np.array(vs))
        t = tree.get(m, 0)
        res[str(m)] = {'log_coeff/(unit*vertex)': [complex(s[1] / (ft.UNIT * t)).real, complex(s[1] / (ft.UNIT * t)).imag] if t != 0 else None,
                       'abs_log': abs(s[1]), 'abs_vertex': abs(t)}
    out[name] = {'claim': mult, 'fit_err': ol.max_fit_err(), 'result': res}
    print(name, res, flush=True)
json.dump(out, open('test_3pt_log.json', 'w'), indent=1)
