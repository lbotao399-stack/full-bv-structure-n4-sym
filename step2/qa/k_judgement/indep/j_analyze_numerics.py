# j_analyze_numerics.py -- judge's own analysis of the stored raw runs of the numerics agent
# (step2/qa/k_transfer_numeric/kn/raw_*.json, read only; a subset was regenerated bit-for-bit by the judge, see
# rerun_numeric/).  Independent fits (different models and point sets from kn_analyze2.py / kn_analyze_eta.py).
#   defect along K_b (arity one, quadratic sector, 1PI): eps -> 0 limit of X0 at every T, T-spread, per-graph T d/dT;
#   flow along K_b: log T coefficient c_f;  arity two along (eta1,eta2): T d/dT Y0 / (eps e^{-eps p^2}).
import json
import numpy as np

KN = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_transfer_numeric/kn/'
OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_judgement/indep/j_analyze_numerics.json'
UNIT = -2.0          # engine units hbar/(16 pi^2) -> unit hbar tr(adad)/(16 pi^2 (X,Y)) for su(2) orthonormal
res = {}


def cp(z):
    return complex(z[0], z[1])


def load(mode, sec, tag='main'):
    return json.load(open(KN + 'raw_%s_%s_%s.json' % (mode, sec, tag)))


def tot(d, which, key, sel=None, fld='X0'):
    s = 0j
    for lab, x in d[which][key].items():
        if sel is not None and d['classes'][lab] not in sel:
            continue
        s += cp(x[fld]) if fld == 'X0' else 0
    return s/UNIT


# ------------------------------------------------------------------ defect along K_b
dres = {}
for sec in ['AW_AW', 'AV_AV', 'l_l', 'c_As']:
    d = load('defect', sec)
    p2 = float(np.dot(d['p1'], d['p1']))
    epss = sorted({float(k.split('|')[0]) for k in d['values']})
    Ts = sorted({float(k.split('|')[1]) for k in d['values']})
    K = lambda e, T: '%g|%g' % (e, T)
    out = {}
    for gname, sel in [('total', None), ('composite', {'D_composite'}), ('anchored', {'B_anchored'})]:
        lims = {}
        for T in Ts:
            es = [e for e in epss if K(e, T) in d['values'] and e <= 3.1e-3]
            Y = np.array([tot(d, 'values', K(e, T), sel) for e in es])
            # model 1 (judge): a + b eps + c eps log eps + d eps^2 through the 4 smallest eps (exact interpolation)
            M = np.array([[1, e, e*np.log(e), e*e] for e in es])
            c1 = np.linalg.solve(M, Y) if len(es) == 4 else np.linalg.lstsq(M, Y, rcond=None)[0]
            # model 2 (judge): with singular terms a + b/eps + c log eps + d eps over all eps <= 3e-3 and 1e-2
            es2 = [e for e in epss if K(e, T) in d['values'] and e <= 1.1e-2]
            Y2 = np.array([tot(d, 'values', K(e, T), sel) for e in es2])
            M2 = np.array([[1, 1/e, np.log(e), e, e*np.log(e)] for e in es2])
            c2 = np.linalg.lstsq(M2, Y2, rcond=None)[0]
            lims['%g' % T] = {'limit_model1': [c1[0].real, c1[0].imag], 'limit_model2': [c2[0].real, c2[0].imag],
                              'eps^-1_model2': [c2[1].real, c2[1].imag], 'log_eps_model2': [c2[2].real, c2[2].imag]}
        L1 = np.array([complex(*v['limit_model1']) for v in lims.values()])
        out[gname] = {'eps_to_0_by_T': lims, 'T_spread_model1': float(np.max(np.abs(L1 - L1[-1])))}
    # per-graph log T structure at the largest T
    Tm = max(Ts)
    per = {}
    for lab, cls in d['classes'].items():
        rows = {}
        for e in epss:
            v = cp(d['TdT'][K(e, Tm)][lab]['X0'])/UNIT
            rows['%g' % e] = {'TdT_X0': [v.real, v.imag],
                              'TdT_X0/(e^{-eps p^2}-1)': [(v/(np.exp(-e*p2) - 1)).real, (v/(np.exp(-e*p2) - 1)).imag],
                              'T*TdT_X0/eps': [(Tm*v/e).real, (Tm*v/e).imag]}
        per[lab] = {'class': cls, 'by_eps': rows}
    # total T d/dT at the largest T, all eps
    tt = {('%g' % e): [tot(d, 'TdT', K(e, Tm)).real, tot(d, 'TdT', K(e, Tm)).imag] for e in epss}
    # T-dependence of the finite value at fixed eps: max_T |X0(eps,T) - X0(eps,Tmax)| for T >= 100, vs eps/T scale
    tdep = {}
    for e in epss:
        ref = tot(d, 'values', K(e, Tm))
        tdep['%g' % e] = float(max(abs(tot(d, 'values', K(e, T)) - ref) for T in Ts if T >= 100 and K(e, T) in d['values']))
    out['per_graph_TdT_at_Tmax'] = per
    out['total_TdT_X0_at_Tmax_by_eps'] = tt
    out['max_T>=100_|X0(eps,T)-X0(eps,Tmax)|'] = tdep
    dres[sec] = out
res['defect_K_arity_one'] = dres

# ------------------------------------------------------------------ flow along K_b
fres = {}
for sec in ['AW_AsW', 'AV_AsV', 'l_ls', 'c_cs']:
    d = load('flow', sec)
    p2 = float(np.dot(d['p1'], d['p1']))
    epss = sorted({float(k.split('|')[0]) for k in d['values']})
    Ts = sorted({float(k.split('|')[1]) for k in d['values']})
    K = lambda e, T: '%g|%g' % (e, T)
    cf = {}
    for e in epss:
        v = tot(d, 'TdT', K(e, max(Ts)))
        cf['%g' % e] = [(v*np.exp(e*p2)).real, (v*np.exp(e*p2)).imag]
    # slope in log T of the eps -> 0 limit (model 1 at each T)
    lim = []
    for T in Ts:
        es = [e for e in epss if e <= 3.1e-3]
        Y = np.array([tot(d, 'values', K(e, T)) for e in es])
        M = np.array([[1, e, e*np.log(e), e*e] for e in es])
        lim.append(np.linalg.solve(M, Y)[0])
    slopes = [((lim[i + 1] - lim[i])/np.log(Ts[i + 1]/Ts[i])) for i in range(len(Ts) - 1)]
    fres[sec] = {'TdT_X0*e^{eps p^2}_at_Tmax': cf, 'eps_to_0_limit_by_T': {('%g' % T): [z.real, z.imag] for T, z in zip(Ts, lim)},
                 'slope_dX0/dlogT_of_limit': [[s.real, s.imag] for s in slopes]}
res['flow_K'] = fres

# ------------------------------------------------------------------ arity two along (eta1, eta2)
eres = {}
for sec in ['AW_AsW', 'AV_AsV', 'l_ls', 'c_cs']:
    per_m = {}
    for tg in ('mA', 'mB'):
        dd = json.load(open(KN + 'raw_eta_%s_%s.json' % (sec, tg)))
        per_m.update(dd['per_m'])
    p1 = np.array(dd['p1']); p2 = float(p1 @ p1)
    keys = sorted({k for e in per_m.values() for k in e['TdT']})

    def Y0(key, sel=None):
        s = 0j
        for mm, ent in per_m.items():
            for lab, jet in ent['TdT'][key].items():
                if sel is not None and ent['classes'][lab] not in sel:
                    continue
                s += cp(jet[3])
        return s/UNIT
    rows = {}
    for key in keys:
        e, T = (float(x) for x in key.split('|'))
        y = Y0(key); yc = Y0(key, {'chain_eta_eta'})
        rows[key] = {'TdT_Y0': [y.real, y.imag], 'TdT_Y0/(eps e^{-eps p^2})': [(y/(e*np.exp(-e*p2))).real, (y/(e*np.exp(-e*p2))).imag],
                     'chain_share': [(yc/y).real, (yc/y).imag] if abs(y) > 0 else None}
    eres[sec] = {'n_m_pairs': len(per_m), 'Gamma_eta1_eta2_V': dd['Gamma_eta1_eta2_V'], 'rows': rows}
res['eta_eta_arity_two'] = eres
json.dump(res, open(OUT, 'w'), indent=1)

# ------------------------------------------------------------------ print a compact summary
for sec, o in dres.items():
    t = o['total']; c = o['composite']; a = o['anchored']
    lt = complex(*t['eps_to_0_by_T']['1e+06']['limit_model1']); lc = complex(*c['eps_to_0_by_T']['1e+06']['limit_model1'])
    la = complex(*a['eps_to_0_by_T']['1e+06']['limit_model1'])
    print('defect', sec, 'limit(T=1e6) total %s = comp %s + anch %s ; T-spread total %.1e comp %.1e anch %.1e' % (
        np.round(lt, 5), np.round(lc, 5), np.round(la, 5), t['T_spread_model1'], c['T_spread_model1'], a['T_spread_model1']))
    print('   model2 eps^-1, log eps at T=1e6:', np.round(complex(*t['eps_to_0_by_T']['1e+06']['eps^-1_model2']), 8),
          np.round(complex(*t['eps_to_0_by_T']['1e+06']['log_eps_model2']), 6))
    for lab, pg in o['per_graph_TdT_at_Tmax'].items():
        if pg['class'] == 'D_composite':
            r = [np.round(complex(*pg['by_eps'][e]['TdT_X0/(e^{-eps p^2}-1)']), 4) for e in ('0.0001', '0.01', '0.1')]
            print('   composite', lab, 'TdT/(e^{-eps p^2}-1) at eps=1e-4,1e-2,1e-1:', r)
    print('   total TdT X0 at T=1e6:', {k: '%.2e' % abs(complex(*v)) for k, v in o['total_TdT_X0_at_Tmax_by_eps'].items()})
for sec, o in fres.items():
    print('flow', sec, 'c_f = TdT X0 e^{eps p^2}:', {k: np.round(complex(*v), 4) for k, v in o['TdT_X0*e^{eps p^2}_at_Tmax'].items() if k in ('0.0001', '0.01', '0.1')},
          'slope of limit:', [np.round(complex(*s), 3) for s in o['slope_dX0/dlogT_of_limit']])
for sec, o in eres.items():
    r = {k: (np.round(complex(*v['TdT_Y0/(eps e^{-eps p^2})']), 3), None if v['chain_share'] is None else np.round(complex(*v['chain_share']), 5))
         for k, v in o['rows'].items() if k.endswith('|10000')}
    print('eta', sec, o['n_m_pairs'], 'Gamma|_V', o['Gamma_eta1_eta2_V'], r)
