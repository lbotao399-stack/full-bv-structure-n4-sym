"""Rerun of the QA off-shell comparison against the regenerated engine-A tables (vertex-vertex sign -1):
B values at A's points (QA work.nosync/b_at_a_*.json, same seeded points) vs A's 'value' = letter - vertex_vertex."""
import json, os, glob
import sympy as sp
QA = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/composite/qa/work.nosync'
TAB = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk/composite/engine_a/kernels_C2.json'))
def asym(s):
    return sp.expand(sp.sympify(s.replace('i', 'I').replace('+-', '-').replace('log2', 'log(2)').replace('log3', 'log(3)').replace('log5', 'log(5)').replace('log7', 'log(7)')))
out = {}; nrow = 0; nok = 0
for f in sorted(glob.glob(os.path.join(QA, 'b_at_a_*.json'))):
    rows = json.load(open(f))
    name = os.path.basename(f)[len('b_at_a_'):-5]
    parts = name.split('_')
    pair = parts[0] + '|' + parts[1]; ls = ' '.join(parts[2:])
    trows = TAB[pair]['leg_sets'][ls]['rows']
    recs = []
    for r in rows:
        t = trows[r['row']]
        same_pt = (t['word'] == r['word_A'] and t['C2_Pbeta']['letter'] == r['A_Pbeta']['letter'])
        rec = {'row': r['row'], 'same_point_as_QA': same_pt}
        ok = same_pt
        for ta, tb in (('C2_Pbeta', 'Pb'), ('C2_0_Delta', 'D0')):
            if ta not in t:
                continue
            a = asym(t[ta]['value'])
            bw = {k.split(' word ')[1]: sp.expand(sp.sympify(v)) for k, v in r['B'].items() if k.startswith(tb + ' ')}
            hits = [w for w, v in bw.items() if sp.expand(v - a) == 0] or (['(B zero)'] if (not bw and a == 0) else [])
            rec[tb] = {'A_value': t[ta]['value'], 'B_matching_words': hits}
            ok &= bool(hits)
        rec['agree'] = bool(ok); nrow += 1; nok += bool(ok)
        recs.append(rec)
    out[name] = recs
    print(name, all(x['agree'] for x in recs))
res = {'rows': nrow, 'agree': nok, 'by_leg_set': out}
json.dump(res, open('bos_compare.json', 'w'), indent=1)
print('rows', nrow, 'agree', nok)
