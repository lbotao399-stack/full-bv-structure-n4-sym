"""Engine B, Step 2: run the I-a graphs of one sector on random plane-wave configurations (resumable), then fit.
Usage: python3 run_ia.py SECTOR NCONF       (raw results: ia_raw_<SECTOR>.json; fits: ia_fit_<SECTOR>.json)"""
import json, os, sys, time
import sympy as sp
import ia_engine as IA
import sectors3 as S

HERE = os.path.dirname(os.path.abspath(__file__))

def run(sector, nconf):
    fn = os.path.join(HERE, f'ia_raw_{sector}.json')
    raw = json.load(open(fn)) if os.path.exists(fn) else {}
    for k in range(nconf):
        if str(k) in raw:
            continue
        t0 = time.time()
        legs = S.config(sector, k)
        res, nt = IA.evaluate(legs)
        raw[str(k)] = dict(terms=nt, seconds=round(time.time() - t0, 1),
                           classes={IA.class_name(cl): {str(n): {b: str(v) for b, v in num.items()} for n, num in d.items()}
                                    for cl, d in res.items()})
        json.dump(raw, open(fn, 'w'), indent=1)
        print(sector, 'config', k, 'terms', nt, 'seconds', raw[str(k)]['seconds'], flush=True)
    return raw

BASIS = ['1', 'log2', 'log3']

def fit(sector, raw):
    confs = sorted(raw.keys(), key=int)
    classes = sorted({cl for k in confs for cl in raw[k]['classes']})
    out = {}
    def fit_values(vals_by_conf, n):
        """vals_by_conf: {conf: {basis: value}} -> fit each basis component on the complete structure list and on G's."""
        names = None; rowsets = {b: [] for b in BASIS}; growsets = {b: [] for b in BASIS}
        for k in confs:
            legs = S.config(sector, int(k))
            st = S.structures(sector, legs, n)
            names = sorted(st)
            gv = {g: S.G_value(g, legs) for g in S.SECTOR_G[sector]}
            for b in BASIS:
                v = sp.sympify(vals_by_conf.get(k, {}).get(b, '0'))
                rowsets[b].append((v, st))
                growsets[b].append((v, gv))
        res = {}
        for b in BASIS:
            if not names:
                allzero = all(r[0] == 0 for r in rowsets[b])
                res[b] = dict(structures={}, G={}, residual_full='0' if allzero else 'nonzero (no structure exists)',
                              residual_G='0' if allzero else 'nonzero')
                continue
            sol, resid = S.solve_exact(rowsets[b], names)
            if n == 0:
                gsol, gres = S.solve_exact(growsets[b], S.SECTOR_G[sector])
            else:
                gsol = {}; gres = [r[0] for r in growsets[b]]
            res[b] = dict(structures={k2: str(v) for k2, v in sol.items() if v != 0},
                          residual_full=str(max([abs(sp.N(x)) for x in resid] + [0])),
                          G={k2: str(v) for k2, v in gsol.items()},
                          residual_G=str(max([abs(sp.N(x)) for x in gres] + [0])))
        return res
    for cl in classes:
        out[cl] = {}
        for n in ('-1', '0'):
            vals = {k: raw[k]['classes'].get(cl, {}).get(n, {}) for k in confs}
            out[cl][n] = fit_values(vals, int(n))
    tot = {}
    for n in ('-1', '0'):
        vals = {}
        for k in confs:
            acc = {}
            for cl in classes:
                for b, v in raw[k]['classes'].get(cl, {}).get(n, {}).items():
                    acc[b] = sp.expand(acc.get(b, 0) + sp.sympify(v))
            vals[k] = {b: str(v) for b, v in acc.items()}
        tot[n] = fit_values(vals, int(n))
    return dict(sector=sector, configs=len(confs), classes=out, total=tot,
                seconds=sum(raw[k]['seconds'] for k in confs))

if __name__ == '__main__':
    sector = sys.argv[1]; nconf = int(sys.argv[2])
    raw = run(sector, nconf)
    res = fit(sector, raw)
    json.dump(res, open(os.path.join(HERE, f'ia_fit_{sector}.json'), 'w'), indent=1)
    print(json.dumps(res['total'], indent=1))
