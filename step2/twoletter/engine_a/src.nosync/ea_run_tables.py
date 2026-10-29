"""T1 at pair weight W_0: every pair type of engine_spec.json, every log / eps^{-m} leg set at W = W_0, every ordered
word and block structure; writes out/kernels_W0.json (and a summary)."""
import json, time, sys, os
from fractions import Fraction
import numpy as np
import ea_core as E, ea_tables as T

SPEC = os.path.join(E.ROOT, 'step2', 'twoletter', 'theory', 'engine_spec.json')


def run(outpath, pairs=None):
    spec = json.load(open(SPEC))
    res = {}
    t00 = time.time()
    for pname, pd in spec['kernels'].items():
        if pairs and pname not in pairs:
            continue
        ta, tb = pd['letters']
        W0 = pd['pair_weight_without_derivatives']
        ls = pd['leg_sets_by_pair_weight'][W0]
        legsets = set(ls['log']) | set(ls['eps_max_m'].keys())
        ent = {'letters': [ta, tb], 'W': W0, 'leg_sets': {}}
        t0 = time.time()
        for legset in sorted(legsets):
            words = T.words_of(legset)
            ent['leg_sets'][legset] = {}
            for w in words:
                for st in T.block_structures(list(w)):
                    tab = T.word_table(ta, tb, w, st)
                    topo = ''.join('Q' if len(b) == 2 else 'C' for b in st) or '0'
                    key = ' '.join(w) + ' | ' + topo if w else '(none) | 0'
                    entry = {}
                    for kind, d in tab.items():
                        kk = 'log' if kind == 'log' else 'eps^-%d' % kind
                        entry[kk] = {'%d,%d,%s' % (ia, ib, ','.join(map(str, oc))): T.poly_str(c, len(w))
                                     for (ia, ib, oc), c in sorted(d.items())}
                    ent['leg_sets'][legset][key] = entry
        ent['seconds'] = round(time.time() - t0, 2)
        res[pname] = ent
        nz = sum(len(v2) for v in ent['leg_sets'].values() for e in v.values() for v2 in e.values())
        print(pname, 'leg sets', len(legsets), 'nonzero entries', nz, 'time', ent['seconds'], flush=True)
    json.dump(res, open(outpath, 'w'))
    print('total', round(time.time() - t00, 1))
    return res


if __name__ == '__main__':
    run(os.path.join(E.HERE, 'out.nosync', 'kernels_W0_raw.json'), sys.argv[1:] or None)
