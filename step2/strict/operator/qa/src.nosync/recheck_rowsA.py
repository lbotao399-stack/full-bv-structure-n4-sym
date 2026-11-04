"""Recheck engine A's compact rows from their raw numbers (no engine code):
 pair form along D0: eps0 = 4 lc + line_part; along P, M0: eps0 = line_part; along R: eps0 = 0;
 Ward pole form (non-line placements): [eps^-m] = -4 m ct_m along D0, 0 along P, M0, R;
 the shifted pole identity of eq c11-pole (with the j >= 1 finite-special-vertex shifts) on single-output rows
 with ct_{m+1} != 0: the S-channel shift {S(Delta), ct_{m+1}}_{H} is nonzero there, while [eps^-m] D already
 equals -4 m ct_m, so the shifted identity fails."""
import json, re, collections
from fractions import Fraction as Fr
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/operator'
A = json.load(open(f'{ROOT}/engine_a/out.nosync/tables_compact.json'))

import sympy as sp
def g(s):
    """parse an engine-A value string into an exact sympy number (Gaussian rationals, logs)"""
    if s is None:
        return sp.Integer(0)
    t = str(s).replace(' ', '').replace('+-', '-').replace('*i', '*I')
    t = re.sub(r'log(\d+)', r'log(\1)', t)
    t = re.sub(r'(?<![A-Za-z0-9_)])i(?![A-Za-z(])', 'I', t)
    return sp.sympify(t, rational=True)
def eq(a, b):
    return sp.simplify(sp.expand(a - b)) == 0
res = collections.Counter(); bad = []; shifted = []
for r in A:
    gen = r['generator']
    if gen not in ('D0', 'P', 'M0', 'R'):
        continue
    e0, line = g(r['eps0']), g(r.get('line_part_eps0', '0'))
    if gen == 'D0':
        ok = eq(e0, 4 * g(r.get('lc', '0')) + line)
    elif gen in ('P', 'M0'):
        ok = eq(e0, line)
    else:
        ok = eq(e0, 0) and eq(line, 0)
    res[(gen, 'pair_form', ok)] += 1
    if not ok:
        bad.append(r)
    if 'ct' in r:
        for m, v in r['poles'].items():
            m_ = int(m); p = g(v); c = g(r['ct'].get(m, '0'))
            target = -4 * m_ * c if gen == 'D0' else 0
            res[(gen, 'ward_pole_recomputed', eq(p, target), r['pole_variation_check'].get(m))] += 1
        if gen == 'D0' and len(r['ordering']) == 1:
            nz = [int(m) for m, v in r['ct'].items() if not eq(g(v), 0)]
            if any(m >= 2 for m in nz):
                shifted.append({'pair': r['pair'], 'n': r['n'], 'leg_set': r['leg_set'], 'D': r['D'], 'ct': r['ct'],
                                'poles': r['poles'], 'statement': '[eps^-1] D = -4 ct_1 exactly; eq c11-pole adds '
                                '+{S(Delta), ct_2}_H = 2(w_Y + D - 4)|q|^2 ct_2 != 0 on the left, which no other term '
                                'on this single-output component can cancel'})
out = {'counts': {str(k): v for k, v in res.items()}, 'n_bad': len(bad), 'bad_examples': bad[:5],
       'shifted_pole_identity_counterexamples': shifted}
json.dump(out, open(f'{ROOT}/qa/out.nosync/recheck_rowsA.json', 'w'), indent=1)
print(json.dumps(out['counts'], indent=0)); print('bad', len(bad)); print('shifted', len(shifted))
for s in shifted: print(s['pair'], s['n'], s['leg_set'], s['ct'], s['poles'])
