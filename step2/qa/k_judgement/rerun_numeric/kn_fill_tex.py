"""fill the (c) table of results.tex from results.json"""
import json, os
import numpy as np
OUT = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(OUT, '..', 'results.json')))
names = {'AW_AsW': '$(A|_W,A^*|_W)$', 'AV_AsV': '$(A|_V,A^*|_V)$', 'l_ls': '$(\\lambda,\\lambda^*)$', 'c_cs': '$(c,c^*)$'}
rows = []
for s in ('AW_AsW', 'AV_AsV', 'l_ls', 'c_cs'):
    if s not in R['c_eta_eta']:
        continue
    g = R['c_eta_eta'][s]['groups']
    tot = g['total']['logT_and_T_dependence']
    c3 = complex(*tot['0.001']['B/(eps e^{-eps p1^2})'])
    c2 = complex(*tot['0.01']['B/(eps e^{-eps p1^2})'])
    B3 = complex(*tot['0.001']['B_Y0_unit'])
    if abs(B3) < 1e-5:
        cet = '$0$ ($|B|<10^{-6}$)'
    else:
        cv = c3
        cet = '$%+.2f%s$' % ((cv.imag if abs(cv.imag) > abs(cv.real) else cv.real), 'i' if abs(cv.imag) > abs(cv.real) else '')
    if 'chain_eta_eta' in g:
        Bc = complex(*g['chain_eta_eta']['logT_and_T_dependence']['0.001']['B_Y0_unit'])
        share = '%.5f' % (abs(Bc) / max(abs(B3), 1e-300)) if abs(B3) > 1e-5 else 'no chain'
    else:
        share = 'no chain graph'
    mt = ' / '.join('%.1e' % tot[e]['max_T |T d/dT Y0|'] for e in ('0.001', '0.01', '0.1'))
    m12 = '%.1e, %.1e' % (tot['0.001']['max_T |T d/dT Y1|'], tot['0.001']['max_T |T d/dT Y2|'])
    rows.append('%s & %s & %s & %s & %s\\\\' % (names[s], cet, share, mt, m12))
tex = open(os.path.join(OUT, '..', 'results.tex')).read()
tex = tex.replace('%CTABLE%', '\n'.join(rows))
open(os.path.join(OUT, '..', 'results.tex'), 'w').write(tex)
print('\n'.join(rows))
