"""Analysis of the raw K_b runs (kn_run_K.py): log T coefficients, eps -> 0 limits at fixed T, T-independence.
Writes ana_<mode>_<sector>.json."""
import sys, os, json
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))
UNIT = -2.0       # unit hbar tr_g(ad ad)/(16 pi^2 (X,Y)) in engine units hbar/(16 pi^2) (su(2))


def cplx(x):
    return x[0] + 1j * x[1]


def load(mode, sector, tag='main'):
    return json.load(open(os.path.join(OUT, 'raw_%s_%s_%s.json' % (mode, sector, tag))))


def grid(d):
    epss = sorted({float(k.split('|')[0]) for k in d['values']})
    Ts = sorted({float(k.split('|')[1]) for k in d['values']})
    return epss, Ts


def key(e, T):
    return '%g|%g' % (e, T)


def collect(d, which, sel, field='X0'):
    """sum over labels with class in sel of the X0 (or X1/X2 flattened) values"""
    out = {}
    for k, v in d[which].items():
        tot = 0
        for lab, x in v.items():
            if d['classes'][lab] in sel:
                if field == 'X0':
                    tot = tot + cplx(x['X0'])
                elif field == 'X1':
                    tot = tot + np.array([cplx(y) for y in x['X1']])
                elif field == 'X2':
                    tot = tot + np.array([[cplx(y) for y in row] for row in x['X2']]).ravel()
        out[k] = tot
    return out


def fitT(Ts, vals):
    """vals(T) = A + B log T + C/T + D log T / T  (T >= 100)"""
    Ts = np.array(Ts); vals = np.array(vals)
    M = np.array([[1, np.log(T), 1 / T, np.log(T) / T] for T in Ts])
    c, *_ = np.linalg.lstsq(M, vals, rcond=None)
    res = vals - M @ c
    return c, float(np.max(np.abs(res))) if len(res) else 0.0


def fiteps(epss, vals):
    """vals(eps) = a/eps + b log eps + c + d eps + e eps log eps (eps <= 1e-2)"""
    epss = np.array(epss); vals = np.array(vals)
    M = np.array([[1 / e, np.log(e), 1, e, e * np.log(e)] for e in epss])
    c, *_ = np.linalg.lstsq(M, vals, rcond=None)
    res = vals - M @ c
    return c, float(np.max(np.abs(res)))


def analyze(mode, sector, tag='main'):
    d = load(mode, sector, tag)
    epss, Ts = grid(d)
    p1 = np.array(d['p1']); p2 = float(p1 @ p1)
    classes = set(d['classes'].values())
    groups = {'total': classes}
    for c in classes:
        groups[c] = {c}
    out = {'mode': mode, 'sector': sector, 'p1': d['p1'], 'b': d['b'], 'units': 'unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))',
           'eps': epss, 'T': Ts, 'groups': {}}
    for gname, sel in groups.items():
        V = collect(d, 'values', sel)
        TD = collect(d, 'TdT', sel)
        g = {}
        # log T coefficient per eps: from T d/dT at the largest T and from a fit in T
        BT = {}
        for e in epss:
            Tl = [T for T in Ts if T >= 100 and key(e, T) in V]
            if len(Tl) < 4:
                continue
            cfit, resid = fitT(Tl, [V[key(e, T)] for T in Tl])
            tdt = TD[key(e, max(Tl))]
            BT['%g' % e] = {'B_from_TdT_at_Tmax_unit': [(tdt / UNIT).real, (tdt / UNIT).imag],
                            'B_from_fit_unit': [(cfit[1] / UNIT).real, (cfit[1] / UNIT).imag],
                            'A_from_fit_unit': [(cfit[0] / UNIT).real, (cfit[0] / UNIT).imag],
                            'fit_residual_unit': resid / abs(UNIT),
                            'B_over_(exp(-eps p1^2)-1)': [(tdt / UNIT / (np.exp(-e * p2) - 1)).real,
                                                         (tdt / UNIT / (np.exp(-e * p2) - 1)).imag],
                            'B_times_exp(eps p1^2)': [(tdt / UNIT * np.exp(e * p2)).real, (tdt / UNIT * np.exp(e * p2)).imag]}
        g['logT_coefficient_X0'] = BT
        # eps -> 0 limit at fixed T (X0)
        lim = {}
        es = [e for e in epss if e <= 1e-2 + 1e-15]
        for T in Ts:
            if not all(key(e, T) in V for e in es) or len(es) < 5:
                continue
            c, resid = fiteps(es, [V[key(e, T)] for e in es])
            lim['%g' % T] = {'eps^-1_unit': [(c[0] / UNIT).real, (c[0] / UNIT).imag],
                             'log_eps_unit': [(c[1] / UNIT).real, (c[1] / UNIT).imag],
                             'finite_unit': [(c[2] / UNIT).real, (c[2] / UNIT).imag],
                             'eps_coeff_unit': [(c[3] / UNIT).real, (c[3] / UNIT).imag],
                             'fit_residual_unit': resid / abs(UNIT)}
        g['eps_to_0_X0_at_fixed_T'] = lim
        if lim:
            fins = np.array([cplx(v['finite_unit']) for v in lim.values()])
            Tk = list(lim.keys())
            ref = cplx(lim[Tk[-1]]['finite_unit'])
            g['finite_X0_spread_over_T_unit'] = float(np.max(np.abs(fins - ref)))
            g['finite_X0_spread_over_T>=1e2_unit'] = float(np.max(np.abs(np.array(
                [cplx(lim[t]['finite_unit']) for t in Tk if float(t) >= 100]) - ref)))
        # T-dependence at fixed eps: X0(eps, T) - X0(eps, Tmax)
        dep = {}
        for e in epss:
            Tl = [T for T in Ts if key(e, T) in V]
            ref = V[key(e, max(Tl))]
            dep['%g' % e] = {('%g' % T): [((V[key(e, T)] - ref) / UNIT).real, ((V[key(e, T)] - ref) / UNIT).imag] for T in Tl}
        g['X0_minus_X0(Tmax)_unit'] = dep
        # X1, X2: T-dependence at the smallest eps and eps -> 0 limits
        for fld in ('X1', 'X2'):
            Vf = collect(d, 'values', sel, fld)
            e0 = min(epss)
            Tl = [T for T in Ts if key(e0, T) in Vf]
            ref = Vf[key(e0, max(Tl))]
            g['%s_max|.-T_max|_at_eps_min_unit' % fld] = {('%g' % T): float(np.max(np.abs((Vf[key(e0, T)] - ref) / UNIT))) for T in Tl}
            limf = {}
            for T in Ts:
                if not all(key(e, T) in Vf for e in es) or len(es) < 5:
                    continue
                M = np.array([[1 / e, np.log(e), 1, e, e * np.log(e)] for e in es])
                Y = np.array([Vf[key(e, T)] for e in es])
                c, *_ = np.linalg.lstsq(M, Y, rcond=None)
                limf['%g' % T] = {'finite_unit': [[x.real, x.imag] for x in (c[2] / UNIT)],
                                  'eps^-1_unit': [[x.real, x.imag] for x in (c[0] / UNIT)],
                                  'log_eps_unit': [[x.real, x.imag] for x in (c[1] / UNIT)]}
            g['%s_eps_to_0' % fld] = limf
            if limf:
                Tk = list(limf.keys())
                fr = np.array([[cplx(x) for x in limf[t]['finite_unit']] for t in Tk])
                g['%s_finite_spread_over_T>=1e2_unit' % fld] = float(np.max(np.abs(
                    fr[[i for i, t in enumerate(Tk) if float(t) >= 100]] - fr[-1])))
        out['groups'][gname] = g
    # per graph at the smallest eps and largest T
    e0 = min(epss); T0 = max(Ts)
    pg = {}
    for lab, x in d['values'][key(e0, T0)].items():
        pg[lab] = {'class': d['classes'][lab], 'X0_unit': [(cplx(x['X0']) / UNIT).real, (cplx(x['X0']) / UNIT).imag],
                   'TdT_X0_unit': [(cplx(d['TdT'][key(e0, T0)][lab]['X0']) / UNIT).real,
                                   (cplx(d['TdT'][key(e0, T0)][lab]['X0']) / UNIT).imag]}
    out['per_graph_at_eps_min_T_max'] = {'eps': e0, 'T': T0, 'graphs': pg}
    # per graph log T coefficients at all eps (TdT at T_max)
    pB = {}
    for lab in d['labels']:
        pB[lab] = {}
        for e in epss:
            Tl = [T for T in Ts if key(e, T) in d['TdT']]
            k_ = key(e, max(Tl))
            if lab in d['TdT'][k_]:
                v = cplx(d['TdT'][k_][lab]['X0']) / UNIT
                pB[lab]['%g' % e] = [v.real, v.imag]
    out['per_graph_TdT_X0_at_Tmax_unit'] = pB
    json.dump(out, open(os.path.join(OUT, 'ana_%s_%s.json' % (mode, sector)), 'w'), indent=1)
    return out


if __name__ == '__main__':
    mode, sector = sys.argv[1], sys.argv[2]
    tag = sys.argv[3] if len(sys.argv) > 3 else 'main'
    o = analyze(mode, sector, tag)
    tot = o['groups']['total']
    print(json.dumps({'logT': tot['logT_coefficient_X0'], 'eps0': tot['eps_to_0_X0_at_fixed_T'],
                      'spread': tot.get('finite_X0_spread_over_T>=1e2_unit')}, indent=1)[:6000])
