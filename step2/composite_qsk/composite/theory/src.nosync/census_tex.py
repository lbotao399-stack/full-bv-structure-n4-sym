"""LaTeX tables from census_sp.json for composite_theory.tex (field pairs and active pairs)."""
import json, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'census_sp.json')))
NAME = {'AV': 'A|_V', 'AW': 'A|_\\perp', 'lamP': '\\lambda_+', 'lamM': '\\lambda_-', 'c': 'c', 'AsV': 'A^*|_V',
        'AsW': 'A^*|_\\perp', 'lamsP': '\\lambda^*_+', 'lamsM': '\\lambda^*_-', 'cs': 'c^*'}
def ls(s):
    toks = s.split()
    out = []
    i = 0
    while i < len(toks):
        j = i
        while j < len(toks) and toks[j] == toks[i]: j += 1
        n = j - i
        out.append(NAME[toks[i]] + ('^{%d}' % n if n > 1 else ''))
        i = j
    return '\\,'.join(out)
order = ['AV|AV', 'AV|AW', 'AV|lamP', 'AV|lamM', 'AW|AW', 'AW|lamP', 'AW|lamM', 'lamP|lamP', 'lamP|lamM', 'lamM|lamM']
rows = []
for key in order:
    v = d['pairs'][key]
    a, b = key.split('|')
    s0 = v['by_letter_derivatives']['0']
    P = s0['leg_sets_Pb_and_position']; T = s0['leg_sets_translation_invariant']
    ent = []
    for s, info in P.items():
        mark = '' if s in T else '{}^{\\dagger}'
        ent.append('$%s%s$ (%d,%d)' % (ls(s), mark, info['k_min'], info['chains']))
    rows.append('$(%s,%s)$ & $%s$ & %s\\\\' % (NAME[a], NAME[b], v['W0'] if '/' not in v['W0'] else '\\tfrac{%s}{%s}' % tuple(v['W0'].split('/')), ', '.join(ent)))
open(os.path.join(HERE, 'tab_field_pairs.tex'), 'w').write('\n'.join(rows) + '\n')
# counts by n for the 21 active pairs and the 10 field pairs
act = [k for k, v in d['pairs'].items() if v['active_pair']]
lines = []
for key in act:
    v = d['pairs'][key]; a, b = key.split('|')
    c = [str(len(v['by_letter_derivatives'][n]['leg_sets_Pb_and_position'])) + '\\,(' + str(len(v['by_letter_derivatives'][n]['leg_sets_translation_invariant'])) + ')' for n in '0123']
    lines.append('$(%s,%s)$ & %s\\\\' % (NAME[a], NAME[b], ' & '.join(c)))
open(os.path.join(HERE, 'tab_counts.tex'), 'w').write('\n'.join(lines) + '\n')
print('\n'.join(rows)); print(); print('\n'.join(lines)); print(len(act))
