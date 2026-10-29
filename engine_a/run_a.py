import json, sys, time
import numpy as np
from oneloop import Engine
from sector_a import bubble_sing, seagull_sing, Bcross, fit, rat
gname = sys.argv[1] if len(sys.argv) > 1 else 'su2'
eng = Engine(gname); d = eng.d
import check_vertices as cv
cv.d = d; cv.SH = {'c': (d,), 'A': (10, d), 'l': (16, d), 'As': (10, d), 'ls': (16, d), 'cs': (d,)}
rng = np.random.default_rng(11)
def rpol(t, which=None):
    sh = eng.shape[t]
    x = rng.integers(-2, 3, size=sh).astype(complex)
    if t == 'A' and which == 'V':
        x[4:] = 0
    if t == 'A' and which == 'perp':
        x[:4] = 0
    return x
SECT = {
    'AV-AV': (('A', 'V'), ('A', 'V'), ['B1', 'B2', 'B3']),
    'Aperp-Aperp': (('A', 'perp'), ('A', 'perp'), ['B4', 'B5']),
    'AV-Aperp': (('A', 'V'), ('A', 'perp'), ['B1', 'B4']),
    'l-l': (('l', None), ('l', None), ['B6']),
    'As-c': (('As', None), ('c', None), ['B7']),
}
out = {'group': gname, 'c_adj': eng.cadj, 'sectors': {}}
for sname, (t1, t2, Bn) in SECT.items():
    samples = []; parts = []
    for it in range(5):
        p = rng.integers(-3, 4, size=4).astype(float)
        if not p.any():
            p[0] = 1
        b1 = (t1[0], rpol(*t1), p); b2 = (t2[0], rpol(*t2), -p)
        bi, bl, det = bubble_sing(eng, b1, b2, p)
        si, sl = seagull_sing(eng, b1, b2, p)
        crosses = {nm: Bcross(nm, b1, b2, p) for nm in Bn}
        samples.append((bi + si, bl + sl, crosses))
        parts.append(dict(bubble_inv=complex(bi), bubble_log=complex(bl), seagull_inv=complex(si), cubic_check=max(v['cubic_check'] for v in det.values())))
    fr = fit(samples, Bn)
    res = {}
    for key in ('eps_inv', 'log_eps'):
        x, resid, rk = fr[key]
        res[key] = {nm: [complex(xx).real / eng.cadj, complex(xx).imag / eng.cadj] for nm, xx in zip(Bn, x)}
        res[key + '_resid'] = float(resid); res[key + '_rank'] = int(rk)
        res[key + '_rational'] = {nm: str(rat(complex(xx).real / eng.cadj)) for nm, xx in zip(Bn, x)}
    res['max_cubic_check'] = max(pp['cubic_check'] for pp in parts)
    out['sectors'][sname] = res
    print(sname, json.dumps({k: v for k, v in res.items() if 'rational' in k or 'resid' in k}), flush=True)
json.dump(out, open(f'logs/sector_a_{gname}.json', 'w'), indent=1, default=str)
