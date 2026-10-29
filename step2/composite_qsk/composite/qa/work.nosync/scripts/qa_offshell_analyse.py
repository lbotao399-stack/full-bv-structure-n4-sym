"""QA: analyse engine B evaluated at engine A's off-shell W0 table points: for each row, which B colour word equals
A's letter + kappa * vertex-vertex part, for kappa in (+1, -1, 0)."""
import json, os, glob
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))


def asym(s):
    if s is None:
        return None
    return sp.expand(sp.sympify(s.replace('i', 'I').replace('+-', '-').replace('log2', 'log(2)').replace('log3', 'log(3)')
                               .replace('log5', 'log(5)').replace('log7', 'log(7)')))


out = {}
for f in sorted(glob.glob(os.path.join(HERE, '..', 'b_at_a_*.json'))):
    rows = json.load(open(f))
    name = os.path.basename(f)[len('b_at_a_'):-5]
    recs = []
    for r in rows:
        rec = {'row': r['row'], 'word_A': r['word_A']}
        for tabA, tabB in (('A_Pbeta', 'Pb'), ('A_D0', 'D0')):
            a = r.get(tabA)
            if a is None:
                continue
            L = asym(a['letter']); Vv = asym(a['vertex_vertex'])
            bw = {k.split(' word ')[1]: sp.expand(sp.sympify(v)) for k, v in r['B'].items() if k.startswith(tabB + ' ')}
            match = {}
            for kap in (1, -1, 0):
                tot = sp.expand(L + kap * Vv)
                hits = [w for w, v in bw.items() if sp.expand(v - tot) == 0]
                if not bw and tot == 0:
                    hits = ['(B: zero)']
                match[str(kap)] = hits
            rec[tabB] = {'A_letter': a['letter'], 'A_vertex_vertex': a['vertex_vertex'], 'B_words': {w: str(v) for w, v in bw.items()},
                         'matching_B_word_for_kappa': match}
        recs.append(rec)
    km1 = all(any(rec[t]['matching_B_word_for_kappa']['-1'] for t in ('Pb', 'D0') if t in rec) and
              all(rec[t]['matching_B_word_for_kappa']['-1'] for t in ('Pb', 'D0') if t in rec) for rec in recs)
    kp1 = all(all(rec[t]['matching_B_word_for_kappa']['1'] for t in ('Pb', 'D0') if t in rec) for rec in recs)
    out[name] = {'rows': recs, 'all_rows_agree_kappa=-1': km1, 'all_rows_agree_kappa=+1': kp1}
    print(name, 'kappa=-1:', km1, 'kappa=+1:', kp1)
json.dump(out, open(os.path.join(HERE, '..', 'offshell_B_at_A.json'), 'w'), indent=1)
