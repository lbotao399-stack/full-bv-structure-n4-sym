"""Assemble results.json for the K-at-finite-T numerics (parts (a), (b), (c) and the validations)."""
import os, json
import numpy as np
import kn_analyze2 as A2
import kn_analyze_eta as AE

OUT = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(OUT)
UNIT = -2.0


def cp(z): return z[0] + 1j * z[1]


def fmt(z, nd=6):
    z = complex(z)
    return [round(z.real, nd), round(z.imag, nd)]


def per_graph_logT(mode, sector, eps_ref=1e-4):
    d = json.load(open(os.path.join(OUT, 'raw_%s_%s_main.json' % (mode, sector))))
    p1 = np.array(d['p1']); p2 = float(p1 @ p1)
    Ts = sorted({float(k.split('|')[1]) for k in d['TdT']})
    out = {}
    for lab in d['labels']:
        rec = {'class': d['classes'][lab], 'logT_coeff_TdT_at_Tmax_unit': {}}
        for e in sorted({float(k.split('|')[0]) for k in d['TdT']}):
            k = '%g|%g' % (e, max(Ts))
            if k in d['TdT'] and lab in d['TdT'][k]:
                v = cp(d['TdT'][k][lab]['X0']) / UNIT
                rec['logT_coeff_TdT_at_Tmax_unit']['%g' % e] = fmt(v, 9)
        k = '%g|%g' % (eps_ref, max(Ts))
        v = cp(d['TdT'][k][lab]['X0']) / UNIT
        if mode == 'defect':
            rec['B_over_(exp(-eps p1^2)-1)_at_eps_1e-4'] = fmt(v / (np.exp(-eps_ref * p2) - 1), 5)
        else:
            rec['B_times_exp(eps p1^2)_at_eps_1e-4'] = fmt(v * np.exp(eps_ref * p2), 5)
        kv = '%g|%g' % (eps_ref, max(Ts))
        rec['X0_at_eps_1e-4_T_max_unit'] = fmt(cp(d['values'][kv][lab]['X0']) / UNIT, 6)
        out[lab] = rec
    return out


def summarize_K(mode, sector):
    o = A2.analyze(mode, sector, 'main')
    res = {'p1': o['p1'], 'b': o['b'], 'eps_grid': o['eps'], 'T_grid': o['T'], 'groups': {}}
    for g, G in o['groups'].items():
        X0 = G['X0']
        rec = {}
        rec['T_dependence_X0_max_over_T>=10_by_eps'] = X0['T_dependence_max_over_T>=10']
        rec['logT_coeff_X0_TdT_at_Tmax_by_eps'] = {e: fmt(complex(*v['B'][0]), 10) for e, v in X0['logT_coeff_TdT_at_Tmax'].items()}
        lim = X0['eps_to_0']
        rec['X0_eps_to_0_limit_by_T_regular_fit'] = {T: fmt(complex(*v['R_finite'][0]), 6) for T, v in lim.items()}
        rec['X0_eps_to_0_fit_residual_regular'] = {T: v['R_resid'] for T, v in lim.items()}
        rec['X0_eps_to_0_singular_fit_eps^-1_and_log_eps_at_Tmax'] = [fmt(complex(*lim[max(lim, key=float)]['S_eps^-1'][0]), 6),
                                                                    fmt(complex(*lim[max(lim, key=float)]['S_log_eps'][0]), 6)]
        rec['X0_limit_T_spread_T>=10'] = X0.get('limit_T_spread_T>=10_R')
        rec['X0_limit_T_spread_T>=100'] = X0.get('limit_T_spread_T>=100_R')
        for fld in ('X1', 'X2'):
            F = G[fld]
            rec['%s_T_dependence_max_over_T>=10_by_eps' % fld] = F['T_dependence_max_over_T>=10']
            limf = F['eps_to_0']
            Tm = max(limf, key=float)
            rec['%s_eps^-1_coefficient_at_Tmax' % fld] = [fmt(complex(*z), 5) for z in limf[Tm]['S_eps^-1']]
            rec['%s_finite_singular_fit_at_Tmax' % fld] = [fmt(complex(*z), 5) for z in limf[Tm]['S_finite']]
            if fld == 'X1':
                rec['X1_logT_coeff_TdT_at_Tmax_by_eps_max'] = {e: max(abs(complex(*z)) for z in v['B']) for e, v in F['logT_coeff_TdT_at_Tmax'].items()}
        res['groups'][g] = rec
    res['per_graph'] = per_graph_logT(mode, sector)
    return res


def main():
    R = {'task': 'independent numerics for K at a finite transfer endpoint T: (a) arity-one defect along K_b, quadratic sector; '
                 '(b) flow W(P(eps,T), I + S^(1)(K_b))^(1); (c) arity-two defect along (eta1, eta2) in Sigma* x Sigma*, Gamma(eta1,eta2)|_V != 0',
         'engine': {'imports': 'engine_a/r2_core.py (Step-1 engine A: field layout, vertex Hessians of I, Pnum, Khat, SDIAG)',
                    'own_code': ['kn/kn_engine.py', 'kn/kn_run_K.py', 'kn/kn_eta.py', 'kn/kn_analyze2.py', 'kn/kn_analyze_eta.py',
                                 'kn/kn_zeromode_kernel.py', 'kn/kn_0d_cycles.py', 'kn/kn_0d_bubble_sign.py', 'kn/kn_test_*.py'],
                    'partial_files_of_earlier_run': 'not used (their 0-dim D-class check had failed; their formula sign for the special edge was not adopted without the own 0-dim check)'},
         'conventions': {
             'K_b': 'v(y) = |y|^2 b - 2 (b.y) y, two-component v^m = x^2 b^m - 2 (b.x) x^m; plane waves S^(1)(e_m e^{iP.y}), y^a -> -i d/dP_a',
             'K_b distribution': 'graph_K = -( X2_ac d_a d_c delta(Sigma) + X1_a d_a delta(Sigma) + X0 delta(Sigma) ), Sigma = total momentum, leg 1 momentum p1 fixed',
             'eta': 'eta1 = e_0 + e_5, eta2 = e_0 - e_5 in the engine-A basis of Sigma^*; eta1 in the +1, eta2 in the -1 eigenspace of rho(e_0 e_1 e_2 e_3) on Sigma^* (two-component: an S^alpha_i-type and an Sbar^{alphadot i}-type special supercharge); Gamma(eta1,eta2)|_V = 2 e_1',
             'eta distribution': 'Y0 = coefficient of delta(Sigma) of -sum_{m1,m2} d_{P1,m1} d_{P2,m2}[...] (overall sign dropped)',
             'special edge': 'K_eps - K_0 = Khat (e^{-eps k^2} - 1): t_0 = eps minus t_0 = 0; special-edge cycle (1/2) sum_sigma str[K S H P S H ...] (0-dim check kn_0d_cycles.json)',
             'propagator': 'P(eps,T)(k) = Pnum(k) (e^{-eps k^2} - e^{-T k^2})/k^2 (ch08 def:c8-PT, prop:c8-twocomp)',
             'units': 'unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -2 x (hbar/(16 pi^2)) for su(2) in an orthonormal basis',
             'logT coefficient': 'T d/dT of the graph at T_max (exact: each P line pinned at t = T), cross-checked by fits in log T',
             'what is summed': '1PI one-loop graphs (cycles) with all orderings and line types; not included: one-particle-reducible graphs (S^(1) on an external leg of a two-point subgraph) and counterterm trees -Q I^CT(K_b), -{S^(1)(K_b), I^CT}_eps; these depend on T only through P(eps,T)(q) at a generic external q, i.e. through e^{-T q^2}'},
         'validations': {}, 'a_defect_K': {}, 'b_flow_K': {}, 'c_eta_eta': {}}
    for fn, key in (('kn_test_step1_cA.json', 'step1_arity0_obstruction_cA'), ('kn_test_step1_AW.json', 'step1_AW_self_energy'),
                    ('kn_test_dil_piint.json', 'dilatation_composite_line_vs_engineA_PartD'), ('kn_0d_cycles.json', 'zero_dim_cycle_formulas'),
                    ('kn_zeromode_kernel.json', 'composite_line_zero_mode_kernel')):
        p = os.path.join(OUT, fn)
        if os.path.exists(p):
            R['validations'][key] = json.load(open(p))
    for s in ('AW_AW', 'AV_AV', 'l_l', 'c_As'):
        if os.path.exists(os.path.join(OUT, 'raw_defect_%s_main.json' % s)):
            R['a_defect_K'][s] = summarize_K('defect', s)
    for s in ('AW_AsW', 'AV_AsV', 'l_ls', 'c_cs'):
        if os.path.exists(os.path.join(OUT, 'raw_flow_%s_main.json' % s)):
            R['b_flow_K'][s] = summarize_K('flow', s)
    for s in ('AW_AsW', 'AV_AsV', 'l_ls', 'c_cs'):
        if os.path.exists(os.path.join(OUT, 'raw_eta_%s_mA.json' % s)):
            o = AE.analyze(s, 'mA+mB')
            p2 = float(np.array(o['p1']) @ np.array(o['p1']))
            rec = {'complete_killing_sum': o['complete_killing_sum'], 'n_m_pairs': o['n_m_pairs'], 'groups': {}}
            for g, G in o['groups'].items():
                Bd = {}
                for e, v in G['logT_and_Tdependence_by_eps'].items():
                    y = complex(*v['B_Y0']); ef = float(e)
                    Bd[e] = {'T_max': v['T'], 'B_Y0_unit': fmt(y, 9), 'B/(eps e^{-eps p1^2})': fmt(y / (ef * np.exp(-ef * p2)), 6),
                             'B/(e^{-eps p1^2}-1)': fmt(y / (np.exp(-ef * p2) - 1), 6),
                             'max_T |T d/dT Y0|': v['max_T_|TdT Y0|'], 'max_T |T d/dT Y1|': v['max_T_|TdT Y1|'], 'max_T |T d/dT Y2|': v['max_T_|TdT Y2|']}
                rec['groups'][g] = {'logT_and_T_dependence': Bd, 'values_Y0': {k: v['Y0'] for k, v in G['values'].items()}}
            R['c_eta_eta'][s] = rec
    # convergence and locality checks
    chk = {}
    try:
        a = json.load(open(os.path.join(OUT, 'raw_defect_AW_AW_main.json'))); b = json.load(open(os.path.join(OUT, 'raw_defect_AW_AW_pd14.json')))
        mx = max(abs(cp(a['values'][k][l]['X0']) - cp(b['values'][k][l]['X0'])) / max(1e-12, abs(cp(b['values'][k][l]['X0'])))
                 for k in b['values'] if k in a['values'] for l in b['values'][k])
        chk['defect_AW_AW_quadrature_pd10_vs_pd14_max_rel'] = mx
        a = json.load(open(os.path.join(OUT, 'raw_flow_AW_AsW_main.json'))); b = json.load(open(os.path.join(OUT, 'raw_flow_AW_AsW_pd14.json')))
        mx = max(abs(cp(a['values'][k][l]['X0']) - cp(b['values'][k][l]['X0'])) / max(1e-12, abs(cp(b['values'][k][l]['X0'])))
                 for k in b['values'] if k in a['values'] for l in b['values'][k])
        chk['flow_AW_AsW_quadrature_pd10_vs_pd14_max_rel'] = mx
    except Exception as ex:
        chk['quadrature_error'] = str(ex)
    for s in ('AW_AW', 'c_As'):
        try:
            res = {}
            for tag in ('main', 'p2x'):
                d = json.load(open(os.path.join(OUT, 'raw_defect_%s_%s.json' % (s, tag))))
                epss = sorted({float(k.split('|')[0]) for k in d['values']})
                T = max(float(k.split('|')[1]) for k in d['values'])
                es = epss[:5]
                Y = np.array([sum(cp(x['X0']) for x in d['values']['%g|%g' % (e, T)].values()) / UNIT for e in es])
                M = np.array([[1, e, e * np.log(e), e * e] for e in es])
                c0 = np.linalg.lstsq(M, Y, rcond=None)[0][0]
                res[tag] = {'p1': d['p1'], 'T': T, 'X0_limit_unit': fmt(c0, 6)}
            res['ratio_X0(2 p1)/X0(p1)'] = fmt(complex(*res['p2x']['X0_limit_unit']) / complex(*res['main']['X0_limit_unit']), 6)
            chk['locality_scaling_%s' % s] = res
        except Exception as ex:
            chk['locality_error_%s' % s] = str(ex)
    R['checks'] = chk
    json.dump(R, open(os.path.join(PARENT, 'results.json'), 'w'), indent=1)
    print('written', os.path.join(PARENT, 'results.json'))


if __name__ == '__main__':
    main()
