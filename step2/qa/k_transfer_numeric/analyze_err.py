import json
out = {}
for mode, secs in [('defect', ['AVAV', 'AWAW', 'll', 'cAs']), ('flow', ['AVAsV', 'AWAsW', 'lls', 'ccs'])]:
    for s in secs:
        m = json.load(open('kv_%s_%s_main.json' % (mode, s)))['values']; h = json.load(open('kv_%s_%s_hi.json' % (mode, s)))['values']
        worst = 0; worst1 = 0
        for k, v in h.items():
            if k not in m: continue
            a = complex(*m[k]['TdT_kap0']); b = complex(*v['TdT_kap0'])
            worst = max(worst, abs(a - b) / max(abs(b), 1e-300))
            a1 = [complex(*x) for x in m[k]['TdT_kap1']]; b1 = [complex(*x) for x in v['TdT_kap1']]
            worst1 = max(worst1, max(abs(x - y) for x, y in zip(a1, b1)) / max(max(abs(y) for y in b1), 1e-300))
        out['%s_%s' % (mode, s)] = {'max_rel_dev_TdT_kap0_main_vs_hi': worst, 'max_rel_dev_TdT_kap1': worst1, 'n_compared': len([k for k in h if k in m])}
        print(mode, s, out['%s_%s' % (mode, s)])
json.dump(out, open('error_estimates.json', 'w'), indent=1)
