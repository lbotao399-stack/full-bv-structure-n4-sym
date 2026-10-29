import json, sys
import numpy as np
for s in ['AVAsV', 'AWAsW', 'lls', 'ccs']:
    d = json.load(open('kv_flow_%s_main.json' % s))
    V = d['values']
    eps_l = sorted({v['eps'] for v in V.values()}); T_l = sorted({v['T'] for v in V.values()})
    print('==', s, 'structures', len(d['structures']))
    print('T\\eps ' + ' '.join('%26g' % e for e in eps_l))
    for T in T_l:
        row = [complex(*V['eps=%g,T=%g' % (e, T)]['TdT_kap0']) for e in eps_l]
        print('%6g ' % T + ' '.join('%12.5e%+12.5ei' % (x.real, x.imag) for x in row))
    for T in T_l:
        k1 = [max(abs(complex(*x)) for x in V['eps=%g,T=%g' % (e, T)]['TdT_kap1']) for e in eps_l]
        k2 = [max(abs(complex(*x)) for r in V['eps=%g,T=%g' % (e, T)]['TdT_kap2'] for x in r) for e in eps_l]
        print('   T=%g  max|TdT kap1| ' % T, ['%.3e' % x for x in k1], ' max|TdT kap2|', ['%.3e' % x for x in k2])
