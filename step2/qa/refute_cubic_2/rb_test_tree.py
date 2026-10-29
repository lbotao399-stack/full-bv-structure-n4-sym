"""Tree-level consistency of my own realisation against Ch.1:
 (i) {S_free, F} = QF for every cubic term of I and every H_j, on random plane waves;
 (ii) Q I_3 = 0 on the four degree-1 sectors;
 (iii) Q I_4 + 1/2 {I_3, I_3} = 0 on (c, A, A, A) and on (c, lam, lam, A).
"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from rb_func import *

rng = np.random.default_rng(20261028)


def rint(shape, lo=-3, hi=4):
    return rng.integers(lo, hi, size=shape).astype(float)


def rand_mom(n):
    ps = [rint(4) for _ in range(n - 1)]
    ps.append(-sum(ps))
    return ps


def leg(slot, p, part=None):
    v = rint(SHAPE[slot])
    if slot in ('A', 'Ast') and part == 'V':
        v[4:] = 0
    if slot in ('A', 'Ast') and part == 'W':
        v[:4] = 0
    return (slot, v, p)


def sector_legs(name):
    if name == 'cAVAV':
        p = rand_mom(3); return [leg('c', p[0]), leg('A', p[1], 'V'), leg('A', p[2], 'V')]
    if name == 'cAWAW':
        p = rand_mom(3); return [leg('c', p[0]), leg('A', p[1], 'W'), leg('A', p[2], 'W')]
    if name == 'cll':
        p = rand_mom(3); return [leg('c', p[0]), leg('lam', p[1]), leg('lam', p[2])]
    if name == 'ccAs':
        p = rand_mom(3); return [leg('c', p[0]), leg('c', p[1]), leg('Ast', p[2], 'V')]
    if name == 'cAAA':
        p = rand_mom(4); return [leg('c', p[0]), leg('A', p[1]), leg('A', p[2]), leg('A', p[3])]
    if name == 'cllA':
        p = rand_mom(4); return [leg('c', p[0]), leg('lam', p[1]), leg('lam', p[2]), leg('A', p[3])]


out = {}
# (i) {S_free, F} = QF; legs: 3-leg sectors chosen so that QF is nonzero generically
tests = {'AAA': (F_AAA, ['cAVAV', 'cAWAW']), 'lamlamA': (F_lamlamA, ['cll']),
         'AcAst': (F_AcAst, ['cAVAV', 'cAWAW', 'ccAs']), 'clamlamst': (F_clamlamst, ['cll']),
         'ccc': (F_ccc, ['ccAs']), 'H1': (F_H1, ['cAVAV']), 'H2': (F_H2, ['cAWAW']),
         'H3': (F_H3, ['cAVAV', 'ccAs']), 'H4': (F_H4, ['cAWAW']), 'H5': (F_H5, ['cll']), 'H6': (F_H6, ['ccAs'])}
res_i = {}
for name, (F, secs) in tests.items():
    worst = 0.0
    nonzero = 0.0
    for sec in secs:
        for _ in range(2):
            L = sector_legs(sec)
            q = Qvalue(F, L)
            b = bracket_value(F_free, 0, 2, F, 3, L)
            worst = max(worst, abs(q - b))
            nonzero = max(nonzero, abs(q))
    res_i[name] = {'max|QF-{S_free,F}|': worst, 'max|QF|': nonzero}
    print(name, res_i[name], flush=True)
out['Sfree_bracket_equals_Q'] = res_i

# (ii) Q I_3 = 0
res_ii = {}
for sec in ['cAVAV', 'cAWAW', 'cll', 'ccAs']:
    vals = [abs(Qvalue(F_I3, sector_legs(sec))) for _ in range(3)]
    res_ii[sec] = max(vals)
print('QI3', res_ii, flush=True)
out['QI3'] = res_ii

# (iii) Q I_4 + 1/2 {I_3,I_3} = 0
res_iii = {}
for sec in ['cAAA', 'cllA']:
    worst = 0.0; scale = 0.0
    for _ in range(2):
        L = sector_legs(sec)
        q = Qvalue(F_I4, L)
        b = bracket_value(F_I3, 0, 3, F_I3, 3, L)
        worst = max(worst, abs(q + 0.5 * b)); scale = max(scale, abs(b))
    res_iii[sec] = {'max|QI4+1/2{I3,I3}|': worst, 'max|{I3,I3}|': scale}
    print(sec, res_iii[sec], flush=True)
out['quartic_CME'] = res_iii
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rb_test_tree.json'), 'w'), indent=1)
