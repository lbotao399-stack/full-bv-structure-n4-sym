import json, os
here = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(here, '..', 'sectors.json')))
TEX = {'c': 'c', 'A_V': r'A_V', 'A_p': r'A_\perp', 'lam+': r'\lambda_+', 'lam-': r'\lambda_-', 'A*_V': r'A^*_V',
       'A*_p': r'A^*_\perp', 'lam*+': r'\lambda^*_+', 'lam*-': r'\lambda^*_-', 'c*': r'c^*'}
GTEX = {'Delta': r'\Delta', 'P_b': r'P_b', 'M': r'M_\omega', 'R': r'R_r', 'xi+': r'\xi_+',
        'xi-': r'\xi_-', 'eta+': r'\eta_+', 'eta-': r'\eta_-', 'K_b': r'K_b'}
def cell(secs):
    return ', '.join('$(' + ','.join(TEX[f] for f in s['legs']) + ')_{%d}$' % s['D'] for s in secs)
def rows(block):
    out = []
    for key, v in block.items():
        name = key.split('^')[0]
        ar = '1' if '(1)' in key else '2'
        gl = ','.join(GTEX[g.strip()] for g in key[key.index(')(') + 2:-1].split(','))
        lab = r'$%s^{(%s)}(%s)$' % (r'\Ob' if name == 'Ob' else 'F_1', ar, gl)
        first = True
        for p in range(v['pmax'] + 1):
            for n in (2, 3, 4, 5):
                secs = [s for s in v['sectors'] if s['p'] == p and s['n'] == n]
                if not secs: continue
                out.append(r'%s & %s & %d & %d & %s \\' % (lab if first else '', ('$' + v['weight'] + '$') if first else '', p, n, cell(secs)))
                first = False
        out.append(r'\hline')
    return '\n'.join(out)
open(os.path.join(here, 'tab_arity1.tex'), 'w').write(rows(d['arity1']))
open(os.path.join(here, 'tab_arity2.tex'), 'w').write(rows(d['arity2']))
