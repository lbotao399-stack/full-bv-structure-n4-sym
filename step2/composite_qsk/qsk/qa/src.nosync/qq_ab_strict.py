"""QA: strict colour-word matching of the A-vs-B comparison files.  Engine A's ordering (Y_1, ..., Y_r) along the
chain is the colour word ad_{Y_r} ... ad_{Y_1}; engine B's word [w_1, ..., w_r] is ad_{Y_{w_1}} ... ad_{Y_{w_r}} with
B's outputs listed in A's order; hence A's ordering corresponds to B's word [r-1, ..., 0]."""
import os, sys, json
import sympy as sp
QA = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/qsk/qa/out.nosync'
summary = {}
for fn in sys.argv[1:]:
    d = json.load(open(os.path.join(QA, fn)))
    rows = []
    for r in d:
        n = len(r['ordering_A'])
        key = str(list(range(n - 1, -1, -1)))
        a = sp.sympify(r['A_eps0'])
        b = sp.sympify(r['B_words_eps0'].get(key, '0'))
        ok = sp.expand(a - b) == 0
        other = {k: v for k, v in r['B_words_eps0'].items() if k != key}
        rows.append({'pair_A': r['pair_A'], 'leg_set_A': r['leg_set_A'], 'ordering_A': r['ordering_A'],
                     'A_eps0': r['A_eps0'], 'B_same_word': str(b), 'strict_agree': ok, 'B_other_words': other})
    na = sum(x['strict_agree'] for x in rows)
    summary[fn] = {'rows': len(rows), 'strict_agree': na, 'disagree': [x for x in rows if not x['strict_agree']]}
    print(fn, 'rows', len(rows), 'strict agree', na)
    for x in rows:
        if not x['strict_agree']:
            print('  DISAGREE', x)
json.dump(summary, open(os.path.join(QA, 'ab_strict_summary.json'), 'w'), indent=1)
