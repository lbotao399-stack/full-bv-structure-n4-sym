"""Summary of refute_cubic_1 (independent direct finite-eps evaluation)."""
import json, numpy as np
L32 = np.log(1.5)
claims = {'(c,A_W,A_W)': {'alpha': 11/4 + 19/4*L32, 'tree': 27/8}, '(c,c,A*_V)': {'alpha': 0.0, 'tree': -27/4}}
res = {'method': ('direct evaluation of the order-hbar three-leg part of Q W(P(eps,infty), I - I^CT[eps]) as '
                  'W(e + r Q e) = W(e) + r QW(e) on plane waves (su(2), explicit Clifford module), W from the one-loop '
                  'supertrace formula at the shifted background e + chi(e) (1PI triangles, quartic bubbles, self-energy '
                  'subgraphs joined to I_3 by P(eps,infty)) minus I^CT_2(e + chi) and I^CT_3; eps = 2e-4 * 2^-k, k = 0..8; '
                  'extrapolation by least squares in {1, log eps, eps log eps, eps, eps^2 log eps, eps^2}; '
                  'variant with the finite counterterm F_2 = (27/16) unit int(d^*A_V, d^*A_V) isolates Ob_3'),
       'sectors': {}}
for seed in (11, 23):
    d = json.load(open(f'main_cubic_seed{seed}.json'))
    for ci, c in enumerate(d['configs']):
        for key, name in (('cAWAW', '(c,A_W,A_W)'), ('ccAsV', '(c,c,A*_V)')):
            r = c[key]
            ict3_log = sum(np.array(v['log']) for k, v in r['components_over_unitG'].items() if k.endswith('ICT3'))
            one_pi_log = sum(np.array(v['log']) for k, v in r['components_over_unitG'].items() if k.endswith('1PI'))
            se_inv = sum(np.array(v['eps^-1']) for k, v in r['components_over_unitG'].items() if k.endswith('selfenergy'))
            row = {'seed': seed, 'config': ci,
                   'total_limit/(unit*G)': r['limit_over_unitG'],
                   'with_F2_limit/(unit*G) [= alpha]': r['limit_with_F2_over_unitG'],
                   'tree_part [= total - with_F2]': (np.array(r['limit_over_unitG']) - np.array(r['limit_with_F2_over_unitG'])).tolist(),
                   'log_eps_total': r['log_eps_coeff_over_unitG'],
                   'eps^-1_total(fit incl 1/eps)': r['eps_inv_coeff_over_unitG(fit incl. 1/eps)'],
                   'sum_legs_-QICT3_log [II-a log]': ict3_log.tolist(),
                   'sum_legs_Q1PI_log': one_pi_log.tolist(),
                   'sum_legs_selfenergy_eps^-1': se_inv.tolist(),
                   'tree_Q_I3': r['tree_QI3']}
            res['sectors'].setdefault(name, {'claim_alpha': claims[name]['alpha'], 'claim_tree': claims[name]['tree'],
                                             'claim_total': claims[name]['alpha'] + claims[name]['tree'], 'rows': []})['rows'].append(row)
for name, s in res['sectors'].items():
    tot = np.array([r['total_limit/(unit*G)'][0] for r in s['rows']])
    alp = np.array([r['with_F2_limit/(unit*G) [= alpha]'][0] for r in s['rows']])
    s['max_abs_dev_total'] = float(np.max(np.abs(tot - s['claim_total'])))
    s['max_abs_dev_alpha'] = float(np.max(np.abs(alp - s['claim_alpha'])))
    s['spread_over_configs_total'] = float(tot.max() - tot.min())
    print(name, 'claim total', s['claim_total'], 'found', tot, 'max dev', s['max_abs_dev_total'])
    print(name, 'claim alpha', s['claim_alpha'], 'found', alp, 'max dev', s['max_abs_dev_alpha'])
    print('   II-a log sums', [np.round(r['sum_legs_-QICT3_log [II-a log]'], 8).tolist() for r in s['rows']])
    print('   eps^-1 totals', [np.round(r['eps^-1_total(fit incl 1/eps)'], 6).tolist() for r in s['rows']])
t2 = json.load(open('test_2pt.json'))
res['calibration_two_point'] = 'see test_2pt.log: 9/4, -1 (A_V and A_W), 4 (lambda), -1 (c,A^*), Ob_2 = -27/8 reproduced'
res['calibration_cubic_log_multiples'] = 'see test_3pt_log.log: 6,6,10,10,1,1,1,1 reproduced'
json.dump(res, open('result.json', 'w'), indent=1)
