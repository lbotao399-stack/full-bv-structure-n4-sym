"""Analysis of the arity-two (eta1, eta2) runs (kn_eta.py).
Killing combination: graph = - sum_{m1,m2} d_{P1,m1} d_{P2,m2} [delta(Sigma+P1+P2) f_{m1 m2}]
  -> coefficients Y0 (delta), Y1_a (d_a delta), Y2_{m1 m2} (d_m1 d_m2 delta):
  Y0 = sum_{m1,m2} f_{m1m2}[h1 h2],  Y1_a = sum_{m1} f_{m1 a}[h1] + sum_{m2} f_{a m2}[h2],  Y2_{m1m2} = f_{m1m2}[1]
  (overall sign -1 and the factor -1 of the c1c2 coefficient dropped: common to every graph).
Reports per group (total, chain_eta_eta, floating_eta, anchored): T d/dT of Y0, Y1, Y2 at each (eps, T);
log T coefficient B(eps) = T d/dT at T_max and its ratio to exp(-eps p1^2) - 1; values where computed."""
import sys, os, json
import numpy as np
OUT = os.path.dirname(os.path.abspath(__file__))
UNIT = -2.0


def cp(z): return z[0] + 1j * z[1]


def combine(per_m, which, key, sel):
    Y0 = 0j; Y1 = np.zeros(4, complex); Y2 = np.zeros((4, 4), complex)
    have = False
    for mm, ent in per_m.items():
        m1, m2 = int(mm[0]), int(mm[1])
        if key not in ent[which]:
            continue
        for lab, jet in ent[which][key].items():
            if ent['classes'][lab] not in sel:
                continue
            have = True
            j = [cp(z) for z in jet]
            Y0 += j[3]
            Y1[m2] += j[1]
            Y1[m1] += j[2]
            Y2[m1, m2] += j[0]
    return (Y0 / UNIT, Y1 / UNIT, Y2 / UNIT) if have else None


def analyze(sector, tag):
    per_m = {}
    d = None
    for tg in tag.split('+'):
        fn = os.path.join(OUT, 'raw_eta_%s_%s.json' % (sector, tg))
        if not os.path.exists(fn):
            continue
        dd = json.load(open(fn))
        d = dd if d is None else d
        per_m.update(dd['per_m'])
    p1 = np.array(d['p1']); p2 = float(p1 @ p1)
    keysT = sorted({k for e in per_m.values() for k in e['TdT']}, key=lambda k: (float(k.split('|')[0]), float(k.split('|')[1])))
    keysV = sorted({k for e in per_m.values() for k in e['values']})
    allc = {c for e in per_m.values() for c in e['classes'].values()}
    groups = {'total': allc}
    for c in allc: groups[c] = {c}
    out = {'sector': sector, 'm_pairs_included': sorted(per_m.keys()), 'n_m_pairs': len(per_m),
           'complete_killing_sum': len(per_m) == 16, 'Gamma_eta1_eta2_V': d['Gamma_eta1_eta2_V'], 'p1': d['p1'],
           'units': 'unit; Y0, Y1, Y2 as in the module docstring', 'groups': {}}
    for g, sel in groups.items():
        G = {'TdT': {}, 'values': {}}
        for k in keysT:
            r = combine(per_m, 'TdT', k, sel)
            if r is None: continue
            Y0, Y1, Y2 = r
            G['TdT'][k] = {'Y0': [Y0.real, Y0.imag], 'max|Y1|': float(np.max(np.abs(Y1))), 'max|Y2|': float(np.max(np.abs(Y2)))}
        for k in keysV:
            r = combine(per_m, 'values', k, sel)
            if r is None: continue
            Y0, Y1, Y2 = r
            G['values'][k] = {'Y0': [Y0.real, Y0.imag], 'Y1': [[z.real, z.imag] for z in Y1]}
        # log T coefficient at the largest T per eps
        epss = sorted({float(k.split('|')[0]) for k in G['TdT']})
        B = {}
        for e in epss:
            Ts = sorted(float(k.split('|')[1]) for k in G['TdT'] if float(k.split('|')[0]) == e)
            kk = '%g|%g' % (e, Ts[-1])
            y = cp(G['TdT'][kk]['Y0'])
            B['%g' % e] = {'T': Ts[-1], 'B_Y0': [y.real, y.imag], 'B/(exp(-eps p1^2)-1)': [(y / (np.exp(-e * p2) - 1)).real, (y / (np.exp(-e * p2) - 1)).imag],
                           'max_T_|TdT Y0|': float(max(abs(cp(G['TdT']['%g|%g' % (e, T)]['Y0'])) for T in Ts)),
                           'max_T_|TdT Y1|': float(max(G['TdT']['%g|%g' % (e, T)]['max|Y1|'] for T in Ts)),
                           'max_T_|TdT Y2|': float(max(G['TdT']['%g|%g' % (e, T)]['max|Y2|'] for T in Ts))}
        G['logT_and_Tdependence_by_eps'] = B
        out['groups'][g] = G
    json.dump(out, open(os.path.join(OUT, 'ana_eta_%s.json' % sector), 'w'), indent=1)
    return out


if __name__ == '__main__':
    o = analyze(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else 'mA+mB')
    print('m pairs:', o['n_m_pairs'], 'complete:', o['complete_killing_sum'])
    for g, G in o['groups'].items():
        print('==', g)
        for e, v in G['logT_and_Tdependence_by_eps'].items():
            print('  eps', e, 'T', v['T'], 'B(Y0) %.6e%+.6ej' % tuple(v['B_Y0']), ' B/(e^{-eps p^2}-1) %.6f%+.6fj' % tuple(v['B/(exp(-eps p1^2)-1)']),
                  ' max_T|TdT Y0| %.2e Y1 %.2e Y2 %.2e' % (v['max_T_|TdT Y0|'], v['max_T_|TdT Y1|'], v['max_T_|TdT Y2|']))
        for k, v in G['values'].items():
            print('  value', k, 'Y0 %.8f%+.8fj' % tuple(v['Y0']))
