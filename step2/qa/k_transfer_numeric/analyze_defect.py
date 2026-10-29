import json, sys
import numpy as np
for s in ['AVAV', 'AWAW', 'll', 'cAs']:
    d = json.load(open('kv_defect_%s_main.json' % s))
    V = d['values']
    eps_l = sorted({v['eps'] for v in V.values()}); T_l = sorted({v['T'] for v in V.values()})
    print('==', s, 'structures', len(d['structures']))
    print('T\\eps ' + ' '.join('%12g' % e for e in eps_l) + '   | TdT kap0 / eps')
    for T in T_l:
        row = []; rat = []
        for e in eps_l:
            v = V['eps=%g,T=%g' % (e, T)]
            k0 = complex(*v['TdT_kap0']); row.append(k0); rat.append(k0 / e)
        print('%6g ' % T + ' '.join('%12.4e' % abs(x) for x in row) + '   | ' + ' '.join('%10.5f%+10.5fi' % (x.real, x.imag) for x in rat[-2:]))
    # lower jets
    mx1 = max(max(abs(complex(*x)) for x in v['TdT_kap1']) / v['eps'] for v in V.values() if v['eps'] <= 1e-4)
    mx2 = max(max(abs(complex(*x)) for r in v['TdT_kap2'] for x in r) / v['eps'] for v in V.values() if v['eps'] <= 1e-4)
    print('max |TdT kap1|/eps, |TdT kap2|/eps (eps<=1e-4):', mx1, mx2)
