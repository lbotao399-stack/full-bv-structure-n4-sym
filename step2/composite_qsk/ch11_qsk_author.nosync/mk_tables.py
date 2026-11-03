"""Build the two census tables of Section 11.6 from the final theory note (tab:qs-fieldpairs, tab:qs-counts),
rewriting the abbreviations of the note in the manuscript notation."""
import re
from pathlib import Path
T = Path('../qsk/theory/qsk_theory.tex').read_text()
AP = r'A|_{W_{\mathbb R}^\perp}'
def conv(s):
    s = s.replace(r'A|_\perp', 'AP@').replace(r'A^*|_\perp', 'ASP@')
    s = s.replace(r'\lambda^*_+', 'LSP@').replace(r'\lambda^*_-', 'LSM@')
    s = s.replace(r'\lambda_+', 'LP@').replace(r'\lambda_-', 'LM@')
    rep = {'ASP@': r'A^*|_{W_{\mathbb R}^\perp}', 'AP@': AP, 'LSP@': r'\lambda^*|_{S_+\otimes W^*}', 'LSM@': r'\lambda^*|_{S_-\otimes W}',
           'LP@': r'\lambda|_{S_+\otimes W}', 'LM@': r'\lambda|_{S_-\otimes W^*}'}
    # powers: X^{2} -> (X)^2 for the long symbols, A|_V^{2} -> A|_V^2
    for k, v in rep.items():
        s = re.sub(re.escape(k) + r'\^\{(\d)\}', lambda m: '(' + v + ')^' + m.group(1), s)
        s = s.replace(k, v)
    s = re.sub(r'A\|_V\^\{(\d)\}', r'A|_V^\1', s)
    return s
def block(label):
    i = T.index(r'\label{' + label + '}')
    j = T.index(r'\hline', i) + len(r'\hline')
    k = T.index(r'\end{longtable}', j)
    rows = [r.strip() for r in T[j:k].split(r'\\') if r.strip()]
    return rows
fp = block('tab:qs-fieldpairs')
out = [r'\begin{tabular}{l|c|l}', r'pair & output weight at most & leg sets along $\xi\in S_+\otimes W$: fewest vertices of $I$; number of chains\\', r'\hline']
for r in fp:
    cells = [c.strip() for c in r.split('&')]
    pair, wb, legs = cells
    legs = legs.replace(r'$\emptyset$', r'$\emptyset$')
    out.append(conv(pair) + ' & ' + wb.replace('$', '').join(['$', '$']) + ' & ' + conv(legs) + r'\\')
out.append(r'\end{tabular}')
Path('tab_defectpairs.tex').write_text('\n'.join(out) + '\n')
cn = block('tab:qs-counts')
out = [r'\begin{tabular}{l|cccc|cccc|cccc}', r'pair & $\xi$, $n=0$ & $1$ & $2$ & $3$ & $\Upsilon^0_\eta$, $n=0$ & $1$ & $2$ & $3$ & $\Upsilon^0_{K_b}$, $n=0$ & $1$ & $2$ & $3$\\', r'\hline']
for r in cn:
    cells = [c.strip() for c in r.split('&')]
    out.append(conv(cells[0]) + ' & ' + ' & '.join(cells[1:]) + r'\\')
out.append(r'\end{tabular}')
Path('tab_defectcounts.tex').write_text('\n'.join(out) + '\n')
print(len(fp), len(cn))
