"""Exact word-by-word comparison A vs B at shared off-shell points.  Both runs list the outputs in A's word order; A's colour
word (outputs in list order) corresponds to B's word with the indices reversed, (n-1, ..., 0) (A writes ad_{Y_r}...ad_{Y_1}
with Y_1 first, B lists the word in the order of application).  A value is compared with B's value on that word (0 when B
lists no such word).  Inputs: fermionic rows (ferm_map/ferm_compare.json), bosonic rows (QA b_at_a_*.json vs the regenerated
engine-A tables, value = letter - vertex_vertex)."""
import json, os, glob
import sympy as sp
RV = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/revise.nosync'
QA = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/composite/qa/work.nosync'
TAB = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/composite/engine_a/kernels_C2.json'))
def asym(s):
    return sp.expand(sp.sympify(s.replace('i', 'I').replace('+-', '-').replace('log2', 'log(2)').replace('log3', 'log(3)').replace('log5', 'log(5)').replace('log7', 'log(7)')))
def bword(n):
    return str(tuple(range(n - 1, -1, -1))) if n > 1 else '(0,)'
res = {'fermionic': [], 'bosonic': []}
fc = json.load(open(os.path.join(RV, 'ferm_map', 'ferm_compare.json')))
for r in fc['rows']:
    n = len(r['word_A']); key = bword(n)
    rec = {'pair': r['pair'], 'word_A': r['word_A'], 'B_word': key}
    ok = True; okp = True
    for k in ('Pb', 'D0'):
        if k not in r:
            continue
        B = sp.sympify(r[k]['B_words'].get(key, '0'))
        a = sp.sympify(r[k]['A_value']); ap = sp.expand(sp.sympify(r[k]['A_letter']) + sp.sympify(r[k]['A_vv']))
        rec[k] = {'A': str(a), 'B': str(B), 'equal': sp.expand(a - B) == 0, 'equal_with_vv_sign_+1': sp.expand(ap - B) == 0,
                  'vv_nonzero': sp.sympify(r[k]['A_vv']) != 0}
        ok &= rec[k]['equal']; okp &= rec[k]['equal_with_vv_sign_+1']
    rec['agree'] = bool(ok); rec['agree_vv+1'] = bool(okp)
    res['fermionic'].append(rec)
for f in sorted(glob.glob(os.path.join(QA, 'b_at_a_*.json'))):
    rows = json.load(open(f)); name = os.path.basename(f)[len('b_at_a_'):-5]
    parts = name.split('_'); pair = parts[0] + '|' + parts[1]; ls = ' '.join(parts[2:])
    trows = TAB[pair]['leg_sets'][ls]['rows']
    for r in rows:
        t = trows[r['row']]
        assert t['word'] == r['word_A']
        n = len(t['word']); key = bword(n)
        rec = {'pair': pair, 'leg_set': ls, 'row': r['row'], 'word_A': t['word'], 'B_word': key}
        ok = True; okp = True; vvnz = False
        for ta, tb in (('C2_Pbeta', 'Pb'), ('C2_0_Delta', 'D0')):
            if ta not in t:
                continue
            B = sp.sympify(r['B'].get(tb + ' word ' + key, '0'))
            a = asym(t[ta]['value']); ap = sp.expand(asym(t[ta]['letter']) + asym(t[ta]['vertex_vertex']))
            rec[tb] = {'A': str(a), 'B': str(B), 'equal': sp.expand(a - B) == 0}
            ok &= rec[tb]['equal']; okp &= (sp.expand(ap - B) == 0); vvnz |= asym(t[ta]['vertex_vertex']) != 0
        rec['agree'] = bool(ok); rec['agree_vv+1'] = bool(okp); rec['vv_nonzero'] = bool(vvnz)
        res['bosonic'].append(rec)
for k in ('fermionic', 'bosonic'):
    L = res[k]
    res[k + '_summary'] = {'rows': len(L), 'agree': sum(x['agree'] for x in L), 'agree_with_vv_sign_+1': sum(x['agree_vv+1'] for x in L)}
    print(k, res[k + '_summary'])
    for x in L:
        if not x['agree']:
            print('  DISAGREE', x)
json.dump(res, open('word_compare.json', 'w'), indent=1)
