"""Parse the numeric tables of manuscript/ch04, ch05, ch06, ch00 and compare every cell
with compare/compare.json and compare/compare_r2.json (exact Fractions)."""
import json, re, sys
from fractions import Fraction as F
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/'
c1 = json.load(open(ROOT + 'compare/compare.json'))
c2 = json.load(open(ROOT + 'compare/compare_r2.json'))
lines = {ch: open(ROOT + 'manuscript/' + ch).read().split('\n') for ch in
         ['ch00_overview.tex', 'ch04_counterterms.tex', 'ch05_operators.tex', 'ch06_dilatation.tex']}

def num(cell):
    s = cell.strip().strip('$').strip()
    s = s.replace('\\tfrac', '\\frac')
    m = re.fullmatch(r'(-?)\\frac\{?(\d+)\}?\{?(\d+)\}?', s)
    if m:
        return F(int(m.group(2)), int(m.group(3))) * (-1 if m.group(1) else 1)
    m = re.fullmatch(r'(-?)\\frac(\d)(\d)', s)
    if m:
        return F(int(m.group(2)), int(m.group(3))) * (-1 if m.group(1) else 1)
    return F(s)

def row(ch, ln):
    t = lines[ch][ln - 1]
    t = t.replace('\\\\', '').strip()
    return [x.strip() for x in t.split('&')]

out = []
def chk(loc, printed, data, src):
    ok = (F(printed) == F(data))
    out.append({'location': loc, 'printed': str(printed), 'data': str(data), 'source': src, 'match': ok})

# ---- quartic (A|_V)^4: rows 528-536, columns 3..10 classes, 11 total, 12 8xtree
T1 = c2['coefficients']['T1']
def cls(setname, aid):
    for pc in T1[setname]['per_class']:
        if pc['engine_A_id'].split(' | ')[1] == aid:
            return pc['log_eps_A']
    raise KeyError(aid)
TR = '(1/24) sum over the 24 orderings s of tr_g(ad_{X_s(0)} ad_{X_s(1)} ad_{X_s(2)} ad_{X_s(3)}) '
CA = '([X_0,X_1],[X_2,X_3]) '; CB = '([X_0,X_3],[X_1,X_2]) '
def key(colour, lor):
    return {'C12': CA, 'C14': CB, 'tr': TR}[colour] + lor
av4_cols = ['AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~',
            'AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~',
            'llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~',
            'AcAs[A|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~',
            'AAAA[A|_V,A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~',
            'AAAA[A|_V,A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~',
            'AAAA[A|_V,A|_V] ~A|_V-A|_V~ -> AAAA[A|_V,A|_V] ~A|_V-A|_V~',
            'AAAA[A|_V,A|_V] ~A|_Wperp-A|_Wperp~ -> AAAA[A|_V,A|_V] ~A|_Wperp-A|_Wperp~']
rowspec = [('C12', '(a_0.a_1)(a_2.a_3)'), ('C12', '(a_0.a_2)(a_1.a_3)'), ('C12', '(a_0.a_3)(a_1.a_2)'),
           ('C14', '(a_0.a_1)(a_2.a_3)'), ('C14', '(a_0.a_2)(a_1.a_3)'), ('C14', '(a_0.a_3)(a_1.a_2)'),
           ('tr', '(a_0.a_1)(a_2.a_3)'), ('tr', '(a_0.a_2)(a_1.a_3)'), ('tr', '(a_0.a_3)(a_1.a_2)')]
tot_key = {'C12': 'b01', 'C14': 'b03'}; lor_key = {'(a_0.a_1)(a_2.a_3)': 'x', '(a_0.a_2)(a_1.a_3)': 'y', '(a_0.a_3)(a_1.a_2)': 'z'}
tot = T1['(A|_V)^4']['total']['sum_of_per_graph_A_general_basis']
for i, (col, lor) in enumerate(rowspec):
    ln = 528 + i; r = row('ch04_counterterms.tex', ln)
    for j, aid in enumerate(av4_cols):
        d = cls('(A|_V)^4', aid); v = d.get(key(col, lor), '0') if isinstance(d, dict) else d
        chk(f'ch04:{ln} col {j+3}', num(r[2 + j]), F(v), 'compare_r2 T1 (A|_V)^4 per class')
    dv = F(tot.get(f'{tot_key[col]} x {lor_key[lor]}', '0')) if col != 'tr' else F(0)
    chk(f'ch04:{ln} total', num(r[10]), dv, 'compare_r2 T1 (A|_V)^4 total (general basis)')
    chk(f'ch04:{ln} 8xtree', num(r[11]), dv, 'compare_r2 eight_times_tree_general_basis')

# ---- (A|_V)^2 (A|_Wperp)^2 rows 549..560 (12 class rows), columns 5,6,7
mix = [['AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~',
        'llA[A|_V] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~'],
       ['AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_V-A|_V~'],
       ['AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['AAA[A|_Wperp] ~A|_V-A|_V~ -> AAAA[A|_V,A|_V] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['AAAA[A|_V,A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~',
        'AAAA[A|_V,A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['AAAA[A|_V,A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['AAAA[A|_V,A|_V] ~A|_V-A|_V~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_V-A|_V~'],
       ['AAAA[A|_V,A|_V] ~A|_Wperp-A|_Wperp~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_Wperp-A|_Wperp~'],
       ['AAAA[A|_V,A|_Wperp] ~A|_V-A|_V~ -> AAAA[A|_V,A|_Wperp] ~A|_Wperp-A|_Wperp~']]
L = '(a_0.a_1)(b_2.b_3)'
for i, aids in enumerate(mix):
    ln = 549 + i; r = row('ch04_counterterms.tex', ln)
    for j, col in enumerate(['C12', 'C14', 'tr']):
        s = F(0)
        for aid in aids:
            d = cls('(A|_V)^2 (A|_Wperp)^2', aid)
            s += F(d.get(key(col, L), '0')) if isinstance(d, dict) else F(d)
        chk(f'ch04:{ln} col {5+j}', num(r[5 + j]), s, 'compare_r2 T1 (A|_V)^2(A|_Wperp)^2 per class (summed over the listed classes)')
for ln, vals in [(561, ('-4', '-8', '0')), (562, ('-4', '-8', '0'))]:
    r = row('ch04_counterterms.tex', ln)
    for j in range(3):
        chk(f'ch04:{ln} col {5+j}', num(r[5 + j]), F(vals[j]), 'compare_r2 total / eight_times_tree (b01 x x -4, b03 x x -8)')

# ---- (A|_Wperp)^4 rows 587..595
aw4_cols = ['AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~ -> AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~',
            'llA[A|_Wperp] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~',
            'AAA[A|_Wperp] ~A|_V-A|_V~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~',
            'AAAA[A|_Wperp,A|_Wperp] ~A|_Wperp-A|_Wperp~ -> AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~',
            'AAAA[A|_Wperp,A|_Wperp] ~A|_V-A|_V~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_V-A|_V~',
            'AAAA[A|_Wperp,A|_Wperp] ~A|_Wperp-A|_Wperp~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_Wperp-A|_Wperp~']
tot = T1['(A|_Wperp)^4']['total']['sum_of_per_graph_A_general_basis']
for i, (col, lor) in enumerate(rowspec):
    ln = 587 + i; r = row('ch04_counterterms.tex', ln)
    lb = lor.replace('a_', 'b_')
    for j, aid in enumerate(aw4_cols):
        d = cls('(A|_Wperp)^4', aid); v = d.get(key(col, lb), '0') if isinstance(d, dict) else d
        chk(f'ch04:{ln} col {j+3}', num(r[2 + j]), F(v), 'compare_r2 T1 (A|_Wperp)^4 per class')
    dv = F(tot.get(f'{tot_key[col]} x {lor_key[lor]}', '0')) if col != 'tr' else F(0)
    chk(f'ch04:{ln} total', num(r[8]), dv, 'compare_r2 T1 (A|_Wperp)^4 total')
    chk(f'ch04:{ln} 8xtree', num(r[9]), dv, 'compare_r2 eight_times_tree_general_basis')

# ---- cubic tables: (line, set, [engine A ids summed])
cubic = [
 (383, '(A|_V, A|_V, A|_V)', ['AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~']),
 (384, '(A|_V, A|_V, A|_V)', ['AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~']),
 (385, '(A|_V, A|_V, A|_V)', ['llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~']),
 (386, '(A|_V, A|_V, A|_V)', ['AcAs[A|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~']),
 (387, '(A|_V, A|_V, A|_V)', ['AAAA[A|_V,A|_V] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~']),
 (388, '(A|_V, A|_V, A|_V)', ['AAAA[A|_V,A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_V] ~A|_Wperp-A|_Wperp~']),
 (398, '(A|_V, A|_Wperp, A|_Wperp)', ['AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~']),
 (399, '(A|_V, A|_Wperp, A|_Wperp)', ['AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAA[A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~']),
 (400, '(A|_V, A|_Wperp, A|_Wperp)', ['llA[A|_V] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~']),
 (401, '(A|_V, A|_Wperp, A|_Wperp)', ['AAA[A|_V] ~A|_V-A|_V~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_V-A|_V~']),
 (402, '(A|_V, A|_Wperp, A|_Wperp)', ['AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> AAAA[A|_Wperp,A|_Wperp] ~A|_Wperp-A|_Wperp~']),
 (403, '(A|_V, A|_Wperp, A|_Wperp)', ['AAAA[A|_V,A|_Wperp] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~']),
 (413, '(lambda, lambda, A|_V)', ['AAA[A|_V] ~A|_V-A|_V~ -> llA[lambda] ~lambda-lambda~ -> llA[lambda] ~A|_V-A|_V~']),
 (414, '(lambda, lambda, A|_V)', ['AAA[A|_V] ~A|_Wperp-A|_Wperp~ -> llA[lambda] ~lambda-lambda~ -> llA[lambda] ~A|_Wperp-A|_Wperp~']),
 (415, '(lambda, lambda, A|_V)', ['llA[lambda] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[lambda] ~A|_V-A|_V~']),
 (416, '(lambda, lambda, A|_V)', ['llA[lambda] ~lambda-lambda~ -> llA[A|_V] ~lambda-lambda~ -> llA[lambda] ~A|_Wperp-A|_Wperp~']),
 (426, '(lambda, lambda, A|_Wperp)', ['AAA[A|_Wperp] ~A|_V-A|_V~ -> llA[lambda] ~lambda-lambda~ -> llA[lambda] ~A|_Wperp-A|_Wperp~']),
 (427, '(lambda, lambda, A|_Wperp)', ['llA[lambda] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[lambda] ~A|_V-A|_V~']),
 (428, '(lambda, lambda, A|_Wperp)', ['llA[lambda] ~lambda-lambda~ -> llA[A|_Wperp] ~lambda-lambda~ -> llA[lambda] ~A|_Wperp-A|_Wperp~']),
 (438, '(c, A|_V, A^*|_V)', ['AcAs[A^*|_V] ~c-A^*|_V~ -> AcAs[c] ~A|_V-A|_V~ -> AAA[A|_V] ~A|_V-A|_V~']),
 (439, '(c, A|_V, A^*|_V)', ['AcAs[A^*|_V] ~c-A^*|_V~ -> AcAs[A|_V] ~c-A^*|_V~ -> AcAs[c] ~A|_V-A|_V~']),
 (449, '(c, A|_Wperp, A^*|_Wperp)', ['AcAs[A^*|_Wperp] ~c-A^*|_V~ -> AcAs[c] ~A|_V-A|_V~ -> AAA[A|_Wperp] ~A|_Wperp-A|_Wperp~']),
 (459, '(c, lambda, lambda^*)', ['AcAs[c] ~c-A^*|_V~ -> cllS[lambda^*] ~lambda-lambda~ -> llA[lambda] ~A|_V-A|_V~']),
 (469, '(c, c, c^*)', ['AcAs[c] ~c-A^*|_V~ -> cccS[c^*] ~c-A^*|_V~ -> AcAs[c] ~A|_V-A|_V~'])]
for ln, s, aids in cubic:
    r = row('ch04_counterterms.tex', ln)
    chk(f'ch04:{ln} log coefficient', num(r[5]), sum(F(cls(s, a)) for a in aids), f'compare_r2 T1 {s}')
for ln, s in [(389, '(A|_V, A|_V, A|_V)'), (404, '(A|_V, A|_Wperp, A|_Wperp)'), (417, '(lambda, lambda, A|_V)'),
              (429, '(lambda, lambda, A|_Wperp)'), (440, '(c, A|_V, A^*|_V)'), (450, '(c, A|_Wperp, A^*|_Wperp)'),
              (460, '(c, lambda, lambda^*)'), (470, '(c, c, c^*)')]:
    r = row('ch04_counterterms.tex', ln)
    chk(f'ch04:{ln} total', num(r[5]), F(c2['final_tables']['T1_totals'][s]['log_eps_multiple_of_tree']), 'compare_r2 final_tables T1_totals')

# ---- two-point tables (census index -> line, eps_inv basis, log basis)
pg = {g['census']: g for g in c1['final_tables']['ICT_two_point_per_graph']}
def cell_terms(cell):
    """parse '$\\tfrac{19}{24}\\int...(dA|_V,dA|_V)-\\tfrac18\\int...(d^*A|_V,...)$' into {basis: coeff}"""
    s = cell.strip().strip('$')
    if s.strip() == '0': return {}
    terms = re.findall(r'([+-]?)\s*(\\tfrac\{?\d+\}?\{?\d+\}?|\d+)?\s*\\int\\mathrm\{dvol\}\\,\(([^)]*)\)', s)
    res = {}
    for sg, c, b in terms:
        v = num(c) if c else F(1)
        if sg == '-': v = -v
        bas = {'A|_V,A|_V': 'B3', 'dA|_V,dA|_V': 'B1', 'd^*A|_V,d^*A|_V': 'B2', 'A|_{W_{\\mathbb R}^\\perp},A|_{W_{\\mathbb R}^\\perp}': 'B5',
               'dA|_{W_{\\mathbb R}^\\perp},dA|_{W_{\\mathbb R}^\\perp}': 'B4', '\\lambda,{\\not d}\\lambda': 'B6', 'dc,A^*': 'B7'}[b]
        res[bas] = v
    return res
tp = [(136, 3), (137, 4), (138, 1), (139, 2), (140, 6), (141, 5), (151, 9), (152, 10), (153, 8), (154, 11), (156, 7), (165, 13), (166, 14), (168, 12)]
aut = {g['census_index']: g['census_Aut'] for g in c1['census']['two_point_matched']}
for ln, ci in tp:
    r = row('ch04_counterterms.tex', ln)
    chk(f'ch04:{ln} 1/|Aut|', num(r[1]), F(1, aut[ci]), 'compare.json census_Aut')
    e = cell_terms(r[2]); l = cell_terms(r[3])
    de = {k: F(v) for k, v in pg[ci]['eps_inv'].items()}; dl = {k: F(v) for k, v in pg[ci]['log_eps'].items()}
    for k in set(e) | set(de): chk(f'ch04:{ln} eps^-1 {k}', e.get(k, F(0)), de.get(k, F(0)), f'compare.json per graph census {ci}')
    for k in set(l) | set(dl): chk(f'ch04:{ln} log {k}', l.get(k, F(0)), dl.get(k, F(0)), f'compare.json per graph census {ci}')

# ---- QME table lines 759..765
rows2 = {r['engine_B_id']: r for r in c2['coefficients']['T2']['rows']}
for ln, bid in [(759, 'QME2-a'), (760, 'QME2-c'), (761, 'QME2-b'), (762, 'QME2-d'), (763, 'QME2-e'), (764, 'QME2-CT')]:
    r = row('ch04_counterterms.tex', ln)
    chk(f'ch04:{ln} eps^-1', num(r[3]), F(rows2[bid]['eps_inv_A']), f'compare_r2 T2 {bid}')
    chk(f'ch04:{ln} finite', num(r[4]), F(rows2[bid]['finite_A']), f'compare_r2 T2 {bid}')
r = row('ch04_counterterms.tex', 765)
chk('ch04:765 total eps^-1', num(r[3]), F(c2['coefficients']['T2']['totals']['eps_inv']['sum_of_rows']), 'compare_r2 T2 totals')
chk('ch04:765 total finite', num(r[4]), F(c2['coefficients']['T2']['totals']['finite']['sum_of_rows']), 'compare_r2 T2 totals')

# ---- ch05 insertion table lines 252..258 (graphs b..f and total; (a) vacuum checked by text)
og = {g['census']: g for g in c1['final_tables']['partial_delta_ICT_scalar_bilinear']['per graph']}
for ln, ci in [(252, 'O2'), (253, 'O3'), (254, 'O4'), (255, 'O5'), (256, 'O6')]:
    r = row('ch05_operators.tex', ln)
    d = og[ci]['log_eps']
    chk(f'ch05:{ln} O', num(r[2]), F(d.get('O', '0')), f'compare.json {ci}')
    chk(f'ch05:{ln} trace', num(r[3]), F(d.get('trace_mix', '0')), f'compare.json {ci}')
    chk(f'ch05:{ln} A|_V', num(r[4]), F(d.get('AV_mix', '0')), f'compare.json {ci}')
r = row('ch05_operators.tex', 258)
ld = c1['final_tables']['partial_delta_ICT_scalar_bilinear']['O = int omega (iota_{w1}A, iota_{w2}A); log eps coefficient']
chk('ch05:258 O', num(r[2]), F(ld['O']), 'compare.json partial_delta_ICT total')
chk('ch05:258 trace', num(r[3]), F(ld['(w1,w2) int omega sum_w (iota_w A, iota_w A)']), 'compare.json')
chk('ch05:258 A|_V', num(r[4]), F(ld['(w1,w2) int omega (A|_V, A|_V)']), 'compare.json')

# ---- dilatation tables ch06 236..241 and ch00 101..106
T3 = {x['row']: x for x in c2['coefficients']['T3']['rows']}
names = ['traceless symmetric scalar bilinear (iota_w1 A, iota_w2 A)', 'trace scalar bilinear sum_w (iota_w A, iota_w A)',
         'letter iota_w A', 'letter A|_V (transverse part)', 'letter lambda', 'letter c (gauge ghost, computed in round 2)']
wO = [F(2), F(2), F(1), F(1), F(3, 2), F(0)]
logc = [F(4), F(-8), 0, 0, 0, 0]; brk = [F(4), F(4), F(2), F(2), F(4), F(1)]
for i, nm in enumerate(names):
    r = row('ch06_dilatation.tex', 236 + i)
    chk(f'ch06:{236+i} w_O', num(r[1]), wO[i], 'DEFINITIONS_v2 Weyl weights')
    chk(f'ch06:{236+i} log coeff', num(r[2]), logc[i], 'compare.json O_eps / partial_delta_ICT')
    chk(f'ch06:{236+i} bracket', num(r[3]), brk[i], 'compare_r2 T3 psi_weights x field count')
    chk(f'ch06:{236+i} coefficient', num(r[4]), F(T3[nm]['A_SW3']), 'compare_r2 T3 A_SW3 (= B_SW3)')
    chk(f'ch06:{236+i} coefficient = 2(log-bracket)', num(r[4]), 2 * (F(logc[i]) - brk[i]), 'recomputed')
    r0 = row('ch00_overview.tex', 101 + i)
    chk(f'ch00:{101+i} w_O', num(r0[1]), wO[i], 'DEFINITIONS_v2 Weyl weights')
    chk(f'ch00:{101+i} coefficient', num(r0[2]), F(T3[nm]['B_SW3']), 'compare_r2 T3 B_SW3')

bad = [x for x in out if not x['match']]
print(len(out), 'cells checked;', len(bad), 'mismatches')
for b in bad: print(b)
json.dump({'n_checked': len(out), 'n_mismatch': len(bad), 'cells': out},
          open(ROOT + 'qa/audit_numbers3_src/tables_vs_data_out.json', 'w'), indent=1)
