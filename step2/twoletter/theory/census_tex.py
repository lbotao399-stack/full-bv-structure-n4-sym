"""Generates the census tables of twoletter_theory.tex (Sections sec:tl-census, sec:tl-bound) from census.json.
Output: census_tables.txt (LaTeX fragments, pasted into the tex by the author)."""
import json, os
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, 'census.json')))
W = {'c': Fr(0), 'AV': Fr(1), 'AW': Fr(1), 'lamP': Fr(3, 2), 'lamM': Fr(3, 2), 'AsV': Fr(3), 'AsW': Fr(3),
     'lamsP': Fr(5, 2), 'lamsM': Fr(5, 2), 'cs': Fr(4)}
DEG = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1, 'AsW': -1, 'lamsP': -1, 'lamsM': -1, 'cs': -2}
TEX = {'c': 'c', 'AV': 'A|_V', 'AW': 'A|_\\perp', 'lamP': '\\lambda_+', 'lamM': '\\lambda_-', 'AsV': 'A^*|_V',
       'AsW': 'A^*|_\\perp', 'lamsP': '\\lambda^*_+', 'lamsM': '\\lambda^*_-', 'cs': 'c^*'}


def fr(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else '\\tfrac{%d}{%d}' % (x.numerator, x.denominator)


def legset(s):
    if s == '(none)':
        return '\\emptyset'
    parts = s.split()
    from collections import Counter
    c = Counter(parts)
    out = []
    for t in ['c', 'AV', 'AW', 'lamP', 'lamM', 'AsV', 'AsW', 'lamsP', 'lamsM', 'cs']:
        if c[t]:
            out.append(('{' + TEX[t] + '}^{%d}' % c[t]) if c[t] > 1 else TEX[t])
    return '\\,'.join(out)


lines = []
lines.append('% Table A: minimal leg sets (no derivatives on the two letters)')
lines.append('\\begin{tabular}{l|c|c|l|l}')
lines.append('pair & $W_0$ & degree & $\\log\\varepsilon$ leg sets at $W=W_0$ & $\\varepsilon^{-1}$ leg sets at $W=W_0$\\\\\\hline')
for key, v in d['output_leg_sets_by_pair_weight'].items():
    a, b = key.split('|')
    W0 = W[a] + W[b]
    row = v[str(W0)]
    logs = ', '.join('$' + legset(s) + '$' for s in row['log'])
    eps = ', '.join('$' + legset(s) + '$' for s in row['eps_max_m'])
    lines.append('$(%s,%s)$ & $%s$ & $%d$ & %s & %s\\\\' % (TEX[a], TEX[b], fr(W0), DEG[a] + DEG[b], logs, eps or '--'))
lines.append('\\end{tabular}')
lines.append('')
lines.append('% Table B: number of log leg sets and maximal number of output fields, by pair weight W <= 7')
Ws = [Fr(n, 2) for n in range(0, 15)]
lines.append('\\begin{tabular}{l|' + 'c' * 8 + '}')
lines.append('pair & $W_0$ & $W_0+1$ & $W_0+2$ & $W_0+3$ & $W_0+4$ & $W_0+5$ & $W_0+6$ & $W_0+7$\\\\\\hline')
tabB = {}
for key, v in d['output_leg_sets_by_pair_weight'].items():
    a, b = key.split('|')
    W0 = W[a] + W[b]
    cells = []
    for n in range(8):
        Wp = W0 + n
        if Wp > 7 or str(Wp) not in v:
            cells.append('')
            continue
        row = v[str(Wp)]
        nmax = max([len(s.split()) for s in row['log']] + [0])
        cells.append('%d\\,(%d)' % (len(row['log']), nmax))
        tabB.setdefault(key, {})[str(Wp)] = dict(n_log_sets=len(row['log']), max_outputs=nmax,
                                                 bound=str(Wp + DEG[a] + DEG[b]))
    lines.append('$(%s,%s)$ & ' % (TEX[a], TEX[b]) + ' & '.join(cells) + '\\\\')
lines.append('\\end{tabular}')
open(os.path.join(HERE, 'census_tables.txt'), 'w').write('\n'.join(lines) + '\n')
# bound check: max outputs <= W + deg
viol = [(k, w, x) for k, v in tabB.items() for w, x in v.items() if x['max_outputs'] > Fr(x['bound'])]
print('Table B bound violations:', viol)
json.dump(tabB, open(os.path.join(HERE, 'census_tableB.json'), 'w'), indent=1)
