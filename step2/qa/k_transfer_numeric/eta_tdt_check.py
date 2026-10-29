"""Arity two along (eta1, eta2), defect, sector (A|_W, A^*|_W): T d/dT at finite T, per cycle structure, for one jet
run (a, b) (all four hyper-dual components), Schwinger quadrature as in kv_driver.  Tests T d/dT = O(eps) at fixed T
for every graph of the census (chains, S^(2) graphs, others)."""
import json, time, itertools
import numpy as np
from kt_engine import *
import kv_driver as KD
import eta_driver as ED

if __name__ == '__main__':
    t0 = time.time()
    sec = ED.EtaSector('defect', 'AWAsW')
    cls = [ED.classify(o, L) for (o, L, _) in sec.structs]
    a, b = 1, 2
    ins = sec.eins(a, b)
    combos = [(1e-3, 1.0), (1e-5, 1.0), (1e-3, 100.0), (1e-5, 100.0)]
    out = {'run_(a,b)': [a, b], 'structures': [], 'combos': combos}
    for si, (order, L, pref) in enumerate(sec.structs):
        gp = GraphPoly2(Cycle([ins[o] for o in order], list(L), pref), E4[a], E4[b])
        vals = {}
        for (eps, T) in combos:
            t, w = KD.tdt_points(gp, eps, T, per=6, nsig=4)
            vals['eps=%g,T=%g' % (eps, T)] = gp.integrate(t, w)
        r = {}
        for T in (1.0, 100.0):
            v3 = vals['eps=0.001,T=%g' % T]; v5 = vals['eps=1e-05,T=%g' % T]
            r['T=%g' % T] = {'max|TdT| eps=1e-3': float(np.abs(v3).max()), 'max|TdT| eps=1e-5': float(np.abs(v5).max()),
                             'ratio(1e-3/1e-5)': float(np.abs(v3).max() / max(np.abs(v5).max(), 1e-300))}
        out['structures'].append({'order': list(order), 'lines': list(L), 'class': cls[si], 'TdT': r})
        print('-'.join(order), '-'.join(L), cls[si], {k: round(v['ratio(1e-3/1e-5)'], 2) for k, v in r.items()},
              '%.2e' % r['T=1']['max|TdT| eps=1e-3'], flush=True)
    tot = {}
    out['seconds'] = round(time.time() - t0, 1)
    json.dump(out, open('eta_tdt_check.json', 'w'), indent=1)
    print('seconds', out['seconds'])
