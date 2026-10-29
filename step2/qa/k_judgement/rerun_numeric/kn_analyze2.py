"""Robust analysis of raw K_b runs: for each group (total / class) and each field X0, X1 (4), X2 (16):
 (1) T-dependence at fixed eps:  delta_T(eps) = max_T>=Tlo |V(eps,T) - V(eps,Tmax)|  -> scaling with eps;
 (2) log T coefficient B(eps) = T d/dT V at T_max (and its ratio to exp(-eps p1^2)-1 / exp(-eps p1^2));
 (3) eps -> 0 limit at fixed T with two models: S (singular: 1/eps, log eps, 1, eps, eps log eps, eps^2) on all eps,
     R (regular: 1, eps, eps log eps, eps^2) on the 5 smallest eps; residuals reported; T-spread of the limit."""
import sys, os, json
import numpy as np
OUT = os.path.dirname(os.path.abspath(__file__))
UNIT = -2.0
def cp(x): return x[0] + 1j * x[1]
def val(x, fld):
    if fld == 'X0': return np.array([cp(x['X0'])])
    if fld == 'X1': return np.array([cp(y) for y in x['X1']])
    return np.array([cp(y) for row in x['X2'] for y in row])
def analyze(mode, sector, tag='main'):
    d = json.load(open(os.path.join(OUT, 'raw_%s_%s_%s.json' % (mode, sector, tag))))
    epss = sorted({float(k.split('|')[0]) for k in d['values']})
    Ts = sorted({float(k.split('|')[1]) for k in d['values']})
    p1 = np.array(d['p1']); p2 = float(p1 @ p1)
    K = lambda e, T: '%g|%g' % (e, T)
    groups = {'total': set(d['classes'].values())}
    for c in set(d['classes'].values()): groups[c] = {c}
    out = {'mode': mode, 'sector': sector, 'p1': d['p1'], 'b': d['b'], 'eps': epss, 'T': Ts,
           'units': 'unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)); X0, X1_a, X2_ac = coefficients of delta, d_a delta, d_a d_c delta (sign: graph_K = -(...))',
           'groups': {}}
    for gname, sel in groups.items():
        G = {}
        for fld in ('X0', 'X1', 'X2'):
            V = {}; TD = {}
            for k, v in d['values'].items():
                V[k] = sum((val(x, fld) for lab, x in v.items() if d['classes'][lab] in sel), 0) / UNIT
            if fld != 'X2':
                for k, v in d['TdT'].items():
                    TD[k] = sum((val(x, fld) for lab, x in v.items() if d['classes'][lab] in sel), 0) / UNIT
            F = {}
            # (1) T dependence
            F['T_dependence_max_over_T>=10'] = {}
            F['T_dependence_max_over_T>=100'] = {}
            for e in epss:
                Tl = [T for T in Ts if K(e, T) in V]
                ref = V[K(e, max(Tl))]
                F['T_dependence_max_over_T>=10']['%g' % e] = float(max(np.max(np.abs(V[K(e, T)] - ref)) for T in Tl if T >= 10))
                F['T_dependence_max_over_T>=100']['%g' % e] = float(max(np.max(np.abs(V[K(e, T)] - ref)) for T in Tl if T >= 100))
            # (2) log T coefficient
            if fld != 'X2':
                F['logT_coeff_TdT_at_Tmax'] = {}
                for e in epss:
                    Tl = [T for T in Ts if K(e, T) in TD]
                    b_ = TD[K(e, max(Tl))]
                    F['logT_coeff_TdT_at_Tmax']['%g' % e] = {'B': [[z.real, z.imag] for z in b_],
                        'B/(exp(-eps p1^2)-1)': [[z.real, z.imag] for z in b_ / (np.exp(-e * p2) - 1)],
                        'B*exp(eps p1^2)': [[z.real, z.imag] for z in b_ * np.exp(e * p2)]}
            # (3) eps -> 0
            F['eps_to_0'] = {}
            lims = []
            for T in Ts:
                es = [e for e in epss if K(e, T) in V]
                if len(es) < 7: continue
                Y = np.array([V[K(e, T)] for e in es])
                MS = np.array([[1 / e, np.log(e), 1, e, e * np.log(e), e * e] for e in es])
                cS, *_ = np.linalg.lstsq(MS, Y, rcond=None); rS = float(np.max(np.abs(Y - MS @ cS)))
                es5 = sorted(es)[:5]; Y5 = np.array([V[K(e, T)] for e in es5])
                MR = np.array([[1, e, e * np.log(e), e * e] for e in es5])
                cR, *_ = np.linalg.lstsq(MR, Y5, rcond=None); rR = float(np.max(np.abs(Y5 - MR @ cR)))
                F['eps_to_0']['%g' % T] = {'S_eps^-1': [[z.real, z.imag] for z in cS[0]], 'S_log_eps': [[z.real, z.imag] for z in cS[1]],
                                          'S_finite': [[z.real, z.imag] for z in cS[2]], 'S_resid': rS,
                                          'R_finite': [[z.real, z.imag] for z in cR[0]], 'R_eps': [[z.real, z.imag] for z in cR[1]],
                                          'R_resid': rR}
                lims.append((T, cS[2], cR[0]))
            if lims:
                refS = lims[-1][1]; refR = lims[-1][2]
                F['limit_T_spread_T>=100_S'] = float(max(np.max(np.abs(c - refS)) for T, c, _ in lims if T >= 100))
                F['limit_T_spread_T>=100_R'] = float(max(np.max(np.abs(c - refR)) for T, _, c in lims if T >= 100))
                F['limit_T_spread_T>=10_R'] = float(max(np.max(np.abs(c - refR)) for T, _, c in lims if T >= 10))
            G[fld] = F
        out['groups'][gname] = G
    json.dump(out, open(os.path.join(OUT, 'ana2_%s_%s.json' % (mode, sector)), 'w'), indent=1)
    return out
if __name__ == '__main__':
    o = analyze(*sys.argv[1:3])
    for g, G in o['groups'].items():
        for fld, F in G.items():
            lim = F['eps_to_0'].get('%g' % max(o['T']), {})
            print(g, fld, 'Tdep(T>=10) by eps:', {k: '%.1e' % v for k, v in F['T_dependence_max_over_T>=10'].items()})
            if fld == 'X0':
                print('   logT B by eps:', {k: '%.3e' % abs(complex(*v['B'][0])) for k, v in F['logT_coeff_TdT_at_Tmax'].items()})
                print('   B/(e^-eps p^2 -1):', {k: '%.6f' % complex(*v['B/(exp(-eps p1^2)-1)'][0]).imag for k, v in F['logT_coeff_TdT_at_Tmax'].items()})
                print('   B*e^(eps p^2):', {k: '%.6f' % complex(*v['B*exp(eps p1^2)'][0]).imag for k, v in F['logT_coeff_TdT_at_Tmax'].items()})
            if lim:
                print('   limit S:', np.round(np.array([complex(*z) for z in lim['S_finite']]), 6)[:4], 'eps^-1', np.round(np.array([complex(*z) for z in lim['S_eps^-1']]), 6)[:4],
                      'log', np.round(np.array([complex(*z) for z in lim['S_log_eps']]), 6)[:4], 'resS %.1e' % lim['S_resid'])
                print('   limit R:', np.round(np.array([complex(*z) for z in lim['R_finite']]), 6)[:4], 'resR %.1e' % lim['R_resid'],
                      'spread S %.1e R %.1e R(T>=10) %.1e' % (F.get('limit_T_spread_T>=100_S', -1), F.get('limit_T_spread_T>=100_R', -1), F.get('limit_T_spread_T>=10_R', -1)))
