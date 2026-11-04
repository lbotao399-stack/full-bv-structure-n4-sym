import numpy as np, sympy as sp
import cb_common as B, cb_s2 as S, cb_tree as T, cb_pairfit as PF, fitlib as FL
Q = B.Q
for pair in ('pm', 'pp'):
    tg = PF.deg0_targets(pair)
    mx = {'g13': 0, 'g2': 0, 'ob0': 0}
    for t in tg:
        x1, x2 = t['ctx']['xi1'], t['ctx']['xi2']
        for nm, (c13, c2v) in (('g13', (1.0, 0.0)), ('g2', (0.0, 1.0))):
            f = Q.dsum((1, T.bracket(S.S2xi_dens(x2), 0, PF.Ob1dens(x1, c13, c2v))), (1, T.bracket(S.S2xi_dens(x1), 0, PF.Ob1dens(x2, c13, c2v))))
            mx[nm] = max(mx[nm], abs(f(t['legs'])))
        mx['ob0'] = max(mx['ob0'], abs(T.bracket(S.s2_dens(x1, x2), 0, Q.F_dsdc_dsA)(t['legs'])))
    print(pair, mx)
import json
for pair in ('pm', 'pp', 'mm'):
    r = json.load(open(f'../out.nosync/B_pair_{pair}_fit.json'))
    print(pair, {b: v['ob_closed'] for b, v in r['blocks'].items()})
