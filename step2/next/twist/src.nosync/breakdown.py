"""per-arity / per-representative-part breakdown of the raw (code-basis) entries; logs as L = log(3/2), log 2"""
import json, sympy as sp, re
L = sp.Symbol('L'); l2 = sp.Symbol('log2')
def simp(x):
    x = sp.expand(sp.expand_log(sp.sympify(x), force=True))
    x = x.subs({sp.log(3): L + l2, sp.log(2): l2})
    return sp.nsimplify(sp.expand(x))
out = {}
for fn in ('entries.json', 'entries_bb.json'):
    d = json.load(open(fn))
    for nm, rec in d.items():
        agg = {}
        for k, (w, a, b) in rec['parts'].items():
            ar = k.split(':')[0]
            tt = k.split(':')[1]
            part = 'lam*' if re.search(r't4,|,4$|t\d,4', tt) and ('b' in nm or 'bb' in nm) else 'F/plain'
            if nm.startswith('bb'):
                ia, ib = tt[1:].split(',')
                part = ('F' if ia != '4' else 'lam*') + '-' + ('F' if ib != '4' else 'lam*')
            elif nm.startswith('bc') or nm.startswith('bgamma'):
                part = 'F' if not tt.startswith('t4') else 'lam*'
            elif nm.startswith('betab'):
                part = 'F' if not tt.endswith(',4') else 'lam*'
            else:
                part = '-'
            key = ar + ':' + part
            A = sp.sympify(w) * sp.sympify(a); B = sp.sympify(w) * sp.sympify(b)
            o = agg.get(key, (0, 0)); agg[key] = (o[0] + A, o[1] + B)
        out[nm] = {'UX': str(simp(rec['UX'])), 'VX': str(simp(rec['VX'])), 'by_arity_part': {k: [str(simp(v[0])), str(simp(v[1]))] for k, v in sorted(agg.items())}}
json.dump(out, open('breakdown.json', 'w'), indent=1)
for nm, r in out.items():
    print(nm, 'UX', r['UX'], 'VX', r['VX'])
    for k, v in r['by_arity_part'].items(): print('    ', k, v)
