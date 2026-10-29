import json, sys
import numpy as np
mode = sys.argv[1]; secs = sys.argv[2].split(',')
for s in secs:
    d = json.load(open('kv_%s_%s_main.json' % (mode, s)))
    V = d['values']; st = d['structures']
    eps = min(v['eps'] for v in V.values()); Ts = sorted({v['T'] for v in V.values()})
    print('==', s, 'eps =', eps, ' per-structure TdT kap0 (/eps for defect) at T =', Ts[-3:])
    tot = {T: 0 for T in Ts}
    for si, (o, L, pr) in enumerate(st):
        vals = [complex(*V['eps=%g,T=%g' % (eps, T)]['per_structure_TdT'][str(si)]['kap0']) for T in Ts[-3:]]
        sc = eps if mode == 'defect' else 1.0
        if max(abs(x) for x in vals) / sc > 1e-6:
            print('  %-22s %-24s' % ('-'.join(o), '-'.join(L)), ' '.join('%11.4e%+11.4ei' % ((x / sc).real, (x / sc).imag) for x in vals))
