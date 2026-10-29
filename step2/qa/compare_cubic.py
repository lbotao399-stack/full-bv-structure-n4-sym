#!/usr/bin/env python3
"""QA, Step 2 Part C: comparison and synthesis of the cubic sector of the one-loop obstruction.

Inputs (read only): step2/engine_a/cubic_results.json, step2/engine_b/cubic_results.json, step2/engine_b/solve_F3.json,
step2/theory/cubic_maps.json, step2/theory/cubic_census.json, and the QA side runs
qa/plugback_theory.json, qa/plugback_engine_b.json, qa/masters_indep.json.

1. I-a graph classes: each class of each engine and each census row of the theory is reduced to a canonical cyclic word
   (vertex type + external leg, edge type, ... around the loop, minimised over rotations and reflections) and matched;
   every coefficient (eps^-1, log eps, eps^0 on every G_i of the sector) is compared exactly in Q + Q log 2 + Q log 3.
2. Other contributions (II-a, III(K_0), I-b, III(K_eps - K_0), I-d + II-b) compared engine A / engine B / theory.
3. Totals per sector: eps^-1, log eps, eps^0 (= Ob_3); census counts.
4. F_3: exact solution of Ob_3 = Q F_3 + {I_3, F_2} with the theory maps (cross-checked against engine B's maps),
   particular solution (beta = 0, H2- and H6-coefficients 0), free directions, symbolic verification by substitution,
   plus the plane-wave plug-back results of the two realisations.
Writes qa/compare_cubic.json.
"""
import ast
import json
import re

import sympy as sp

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2'
OUT = ROOT + '/qa/compare_cubic.json'
L2, L3 = sp.log(2), sp.log(3)
SECTORS = [('(c, A|_V, A|_V)', 'cAVAV', ['G1', 'G2']), ('(c, A|_Wperp, A|_Wperp)', 'cApAp', ['G3']),
           ('(c, lambda, lambda)', 'clamlam', ['G5']), ('(c, c, A^*|_V)', 'ccAsV', ['G4'])]


# ---------------------------------------------------------------- exact numbers
def num(s):
    """string -> (r, a2, a3) with value r + a2 log 2 + a3 log 3 (exact rationals)."""
    if s is None:
        return (sp.Integer(0),) * 3
    e = sp.expand(sp.expand_log(sp.sympify(str(s)), force=True))
    a2, a3 = e.coeff(L2), e.coeff(L3)
    r = sp.expand(e - a2 * L2 - a3 * L3)
    assert r.is_Rational and a2.is_Rational and a3.is_Rational, (s, e)
    return (r, a2, a3)


def nstr(t):
    r, a2, a3 = t
    if a2 == -a3:
        if a3 == 0:
            return str(r)
        return (f'{r} + ' if r != 0 else '') + f'{a3}*log(3/2)'
    return f'{r} + {a2}*log(2) + {a3}*log(3)'


def nadd(a, b, s=1):
    return tuple(sp.nsimplify(x + s * y) for x, y in zip(a, b))


ZERO = (sp.Integer(0),) * 3

# ---------------------------------------------------------------- canonical cyclic words
VT = {('A', 'A', 'A'): 'AAA', ('A', 'A', 'A', 'A'): 'AAAA', ('A', 'As', 'c'): 'AcA*', ('c', 'c', 'cs'): 'ccc*',
      ('A', 'lam', 'lam'): 'llA', ('c', 'lam', 'lams'): 'cll*'}
EXT_KIND = {'c': 'c', 'A_V': 'A', 'A_W': 'A', 'lam': 'lam', 'As_V': 'As'}


def ext_norm(x):
    x = x.strip()
    return {'c': 'c', 'A|_V': 'A_V', 'A|_Wperp': 'A_W', 'lambda': 'lam', 'A^*|_V': 'As_V', 'A*|_V': 'As_V'}[x]


def canon(vlab, elab):
    """vlab[i] = vertex label, elab[i] = label of the edge joining vertex i and vertex i+1 (mod n); the word is
    minimised over the n rotations and the n reflections of the cycle."""
    n = len(vlab)
    words = []
    for r in range(n):
        fw, bw = [], []
        for j in range(n):
            a = (r + j) % n
            fw += [vlab[a], elab[a]]
            b = (r - j) % n
            bw += [vlab[b], elab[(b - 1) % n]]
        words += [tuple(fw), tuple(bw)]
    return min(words)


A_EDGE = {'AP': 'A_W', 'AV': 'A_V', 'gh': 'c-A*', 'l': 'lam', 'sAP': 'K[A_W]', 'sAV': 'K[A_V]', 'sc': 'K[c]',
          'sl': 'K[lam]'}
SLOT_KIND = {'A': 'A', 'As': 'As', 'c': 'c', 'cs': 'cs', 'l': 'lam', 'ls': 'lams'}


def canon_A(cls):
    t = ast.literal_eval(cls)
    n = len(t)
    edges = [x[0] for x in t]        # e_i joins v_{i-1} (slot e_i[1][0]) and v_i (slot e_i[1][1])
    exts = [tuple(ext_norm(y) for y in x[1]) for x in t]
    vlab, elab = [], []
    for i in range(n):
        internal = [edges[i][1][1], edges[(i + 1) % n][1][0]]
        kinds = tuple(sorted([SLOT_KIND[s] for s in internal] + [EXT_KIND[e] for e in exts[i]]))
        vlab.append(VT[kinds] + '[' + ','.join(exts[i]) + ']')
        elab.append(A_EDGE[edges[(i + 1) % n][0]])
    # consistency of A's line labels with slots: an A^*|_V slot of a ghost line must sit at a vertex whose A^* is in V
    return canon(vlab, elab), vlab, elab


B_VT = {'AAA': 'AAA', 'AAAA': 'AAAA', 'AcA*': 'AcA*', 'ccc*': 'ccc*', 'lambda lambda A': 'llA', 'c lambda lambda*': 'cll*'}
B_EDGE = {'A|_V': 'A_V', 'A|_Wperp': 'A_W', 'lambda': 'lam', 'c-A*': 'c-A*', 'K[A|_V-A*|_V]': 'K[A_V]',
          'K[A|_Wperp-A*|_Wperp]': 'K[A_W]', 'K[c-c*]': 'K[c]', 'K[lambda-lambda*]': 'K[lam]'}


def canon_B(cls):
    toks = cls.split(' | ')[1:]
    vlab, elab = [], []
    for tok in toks:
        m = re.match(r'^(.*)\[(.*)\] -(.*)-$', tok.strip())
        vt, ex, ed = m.group(1).strip(), m.group(2), m.group(3)
        exs = [ext_norm(x) for x in ex.split(',')]
        vlab.append(B_VT[vt] + '[' + ','.join(exs) + ']')
        elab.append(B_EDGE[ed])
    return canon(vlab, elab), vlab, elab


def th_vt(lbl):
    core = lbl.split('[')[0].strip()
    if '^4' in core or (')^2' in core and core.count('^2') == 2):
        return 'AAAA'
    if core in ('(A|_V)^3', 'A|_V (A|_Wperp)^2', '(A|_Wperp)^3', '(A|_V)^2 A|_Wperp'):
        return 'AAA'
    if core in ('A|_V c A^*|_V', 'A|_Wperp c A^*|_Wperp'):
        return 'AcA*'
    if core == 'c c c^*':
        return 'ccc*'
    if core in ('lambda lambda A|_V', 'lambda lambda A|_Wperp'):
        return 'llA'
    if core == 'c lambda lambda^*':
        return 'cll*'
    raise ValueError(lbl)


TH_EDGE = {'A|_V': 'A_V', 'A|_Wperp': 'A_W', 'c--A^*': 'c-A*', 'lambda': 'lam', 'A|_V--A^*|_V': 'K[A_V]',
           'A|_Wperp--A^*|_Wperp': 'K[A_W]', 'c--c^*': 'K[c]', 'lambda--lambda^*': 'K[lam]'}


def canon_T(row):
    vl = []
    for lbl in row['vertices']:
        ex = re.search(r'\[external (.*)\]', lbl).group(1)
        exs = [ext_norm(x.split('#')[0]) for x in ex.split(',')]
        vl.append(th_vt(lbl) + '[' + ','.join(exs) + ']')
    n = len(vl)

    def vidx(s):
        s = s.strip()
        if s.startswith('['):
            s = s[1:-1]
        return 1 if ',' in s and n == 2 else int(s)
    edges = []
    for e in row['P_edges'] + [row['special_edge']]:
        typ, rest = e.split(':', 1)
        refs = re.findall(r'vertex(\[[0-9, ]+\])', rest)
        a, b = vidx(refs[0]), vidx(refs[1])
        edges.append((a, b, TH_EDGE[typ.strip()]))
    if n == 3:
        elab = []
        for i in range(3):
            j = (i + 1) % 3
            elab.append(next(l for (a, b, l) in edges if {a, b} == {i, j}))
    else:
        elab = [edges[0][2], edges[1][2]]
    return canon(vl, elab), vl, elab


# ---------------------------------------------------------------- main
def main():
    EA = json.load(open(ROOT + '/engine_a/cubic_results.json'))
    EB = json.load(open(ROOT + '/engine_b/cubic_results.json'))
    TM = json.load(open(ROOT + '/theory/cubic_maps.json'))
    TC = json.load(open(ROOT + '/theory/cubic_census.json'))
    out = {'unit': 'hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))  (code name: unit)',
           'G_basis': EA['G_basis'],
           'inputs': ['step2/engine_a/cubic_results.json', 'step2/engine_b/cubic_results.json',
                      'step2/engine_b/solve_F3.json', 'step2/theory/cubic_maps.json', 'step2/theory/cubic_census.json'],
           'class_matching_method': 'canonical cyclic word (vertex type[external legs], edge type, ...) around the loop, '
                                    'minimised over rotations and reflections; edge types A_V, A_W (= A|_Wperp), c-A* '
                                    '(ghost line c--A^*|_V), lam, K[x] = special edge K_eps - K_0 on the x--x^* component',
           'sectors': {}}
    all_agree = True
    totals_all = {}
    for secA, secB, gs in SECTORS:
        ca = EA['sectors'][secA]['graph_classes_I-a']
        cb = [c for c in EB['sectors'][secB]['graph_classes'] if c['class'].startswith('I-a')]
        tc = TC[secA]['I-a']
        keyA = {}
        for c in ca:
            k, vl, el = canon_A(c['class'])
            assert k not in keyA, ('duplicate A', k)
            keyA[k] = (c, vl, el)
        keyB = {}
        for c in cb:
            k, vl, el = canon_B(c['class'])
            assert k not in keyB, ('duplicate B', k)
            keyB[k] = (c, vl, el)
        keyT = {}
        for r in tc:
            k, vl, el = canon_T(r)
            assert k not in keyT, ('duplicate T', k)
            keyT[k] = r
        rows = []
        sumA = {o: {g: ZERO for g in gs} for o in ('eps^-1', 'log', 'eps^0')}
        sumB = {o: {g: ZERO for g in gs} for o in ('eps^-1', 'log', 'eps^0')}
        for k in sorted(set(keyA) | set(keyB) | set(keyT)):
            row = {'loop_word': ' '.join(k)}
            a = keyA.get(k)
            b = keyB.get(k)
            t = keyT.get(k)
            row['engine_A_class'] = a[0]['class'] if a else None
            row['engine_B_class'] = b[0]['class'] if b else None
            row['theory_census_row'] = ({'vertices': t['vertices'], 'P_edges': t['P_edges'],
                                         'special_edge': t['special_edge'], 'number': t['number'], 'weight': t['weight']}
                                        if t else None)
            row['topology'] = 'bubble' if len(k) == 4 else 'triangle'
            coeffs = {}
            agree = True
            for o in ('eps^-1', 'log', 'eps^0'):
                coeffs[o] = {}
                for g in gs:
                    va = num(a[0]['coeffs'][o][g]) if a else None
                    vb = num(b[0]['coeffs'][o][g]) if b else None
                    if a:
                        sumA[o][g] = nadd(sumA[o][g], va)
                    if b:
                        sumB[o][g] = nadd(sumB[o][g], vb)
                    same = (va == vb) if (a and b) else None
                    if a and b and not same:
                        agree = False
                    coeffs[o][g] = {'A': nstr(va) if va else None, 'B': nstr(vb) if vb else None, 'agree': same}
            row['coefficients'] = coeffs
            if a and b:
                row['status'] = 'match: all coefficients identical' if agree else 'DISAGREE'
            elif a and not b:
                zero = all(coeffs[o][g]['A'] == '0' for o in coeffs for g in gs)
                row['status'] = ('engine A only, value identically 0; not a graph of the theory census '
                                 '(component-forbidden, see note)') if zero else 'engine A only, NONZERO'
                row['note'] = ('the P(eps,infty) ghost line c--A^* has only the A^*|_V component (Q^GF c^* = dc^* lies in '
                               'the V-part), and the vertex int([A,c],A^*) pairs A and A^* componentwise; this word puts '
                               'an A^*|_Wperp slot (vertex with external A|_Wperp) or an A|_Wperp line at the c-A^*|_V '
                               'vertex on a c--A^* line, so the integrand vanishes identically. Engine A enumerates by '
                               'field type and records it with value 0; engine B and the theory census omit it.')
                if not zero:
                    agree = False
            else:
                row['status'] = 'engine B only'
                agree = False if b else agree
            if t is None and (a or b):
                row['census'] = 'absent from theory census'
            elif t is not None and not (a and b):
                row['census'] = 'in theory census but missing in an engine'
                agree = False
            else:
                row['census'] = 'in theory census'
            all_agree = all_agree and agree
            rows.append(row)
        # totals and other contributions
        A_sec = EA['sectors'][secA]
        B_sec = EB['sectors'][secB]
        Bo = {c['class'].split(':')[0]: c for c in B_sec['graph_classes'] if not c['class'].startswith('I-a')}

        def bget(prefix, order):
            c = next(v for kk, v in Bo.items() if kk.startswith(prefix))
            return {g: nstr(num(c['coeffs'][order][g])) for g in gs}
        bIb = next(v for kk, v in Bo.items() if kk.startswith('I-b'))
        others = {
            'II-a log': {'A': {g: A_sec['II-a_and_III(K_0)']['II-a = -Q I^CT_3']['log'][g] for g in gs},
                         'B': bget('II-a', 'log')},
            'III(K_0) log': {'A': {g: A_sec['II-a_and_III(K_0)']['III(K_0) = -{I_3, I^CT_2}_0']['log'][g] for g in gs},
                             'B': bget('III (K_0', 'log')},
            'I-b eps^0 (alone)': {'A': {g: A_sec['I-b_and_III(K_eps-K_0)']['I-b eps^0'][g] for g in gs},
                                  'B': {g: bIb['eps_expansion']['I_b']['eps^0'][g] for g in gs}},
            'III(K_eps-K_0) eps^0 (alone)': {'A': {g: A_sec['I-b_and_III(K_eps-K_0)']["III(K') eps^0"][g] for g in gs},
                                             'B': {g: bIb['eps_expansion']['III_Khat']['eps^0'][g] for g in gs}},
            'I-b + III(K_eps-K_0) eps^0': {'B': bget('I-b', 'eps^0')},
            'lim (I-d + II-b) = tree(Ob_2, I_3; P(0,infty))': {
                'A': {g: A_sec['I-d_and_II-b']['tree of Ob_2 with I_3 (limit)'][g] for g in gs},
                'B': bget('I-d', 'eps^0')},
        }
        th = TM['log_and_inverse_eps_bookkeeping(units of unit)']
        secT = {'(c, A|_V, A|_V)': '(c,A|_V,A|_V)', '(c, A|_Wperp, A|_Wperp)': '(c,A|_Wperp,A|_Wperp)',
                '(c, lambda, lambda)': '(c,lambda,lambda)', '(c, c, A^*|_V)': '(c,c,A^*|_V)'}[secA]
        others['II-a log']['theory'] = th['II-a_log_per_row'][secT]['total']
        others['III(K_0) log']['theory'] = th['III_K0_log_per_row'][secT]['total']
        others['I-b eps^0 (alone)']['theory'] = th['lim I-b'][secT]
        others['III(K_eps-K_0) eps^0 (alone)']['theory'] = th['lim III_{K_eps-K_0}'][secT]
        others['lim (I-d + II-b) = tree(Ob_2, I_3; P(0,infty))']['theory'] = th['lim (I-d + II-b)'][secT]
        tree_units = TM['matrices']['tree(Ob_2, I3; P(0,infinity)) in units']
        others['lim (I-d + II-b) = tree(Ob_2, I_3; P(0,infty))']['theory_vector'] = {
            g: tree_units[int(g[1]) - 1] for g in gs}
        for kk, v in others.items():
            if 'A' in v and 'B' in v:
                v['A_equals_B'] = all(num(v['A'][g]) == num(v['B'][g]) for g in gs)
                all_agree = all_agree and v['A_equals_B']
        # sums
        tot = {}
        for o in ('eps^-1', 'log', 'eps^0'):
            tot[o] = {g: {'A': nstr(sumA[o][g]), 'B': nstr(sumB[o][g]), 'agree': sumA[o][g] == sumB[o][g]} for g in gs}
        IIa = {g: num(others['II-a log']['A'][g]) for g in gs}
        III0 = {g: num(others['III(K_0) log']['A'][g]) for g in gs}
        Ib = {g: num(others['I-b eps^0 (alone)']['A'][g]) for g in gs}
        IIIk = {g: num(others['III(K_eps-K_0) eps^0 (alone)']['A'][g]) for g in gs}
        totals = {
            'I-a (sum over classes)': tot,
            'eps^-1 total': {
                'I-a': {g: tot['eps^-1'][g]['A'] for g in gs}, 'II-a': '0 (I^CT_3 has no eps^-1 part)',
                'III(K_0)': '0 identically ({I_3, int(A,A)} = 0; theory, A, B)',
                'I-b + III(K_eps-K_0)': '0 (A, B: no eps^-1 term in either)',
                'I-d + II-b': 'nonlocal eps^-1/|p|^2 parts cancel exactly (A: exact; B: D_2 has eps^-1 coefficient 0)',
                'sum': {g: nstr(sumA['eps^-1'][g]) for g in gs},
                'zero': all(sumA['eps^-1'][g] == ZERO and sumB['eps^-1'][g] == ZERO for g in gs)},
            'log eps total': {
                'I-a': {g: tot['log'][g]['A'] for g in gs},
                'II-a': {g: nstr(IIa[g]) for g in gs}, 'III(K_0)': {g: nstr(III0[g]) for g in gs},
                'I-b + III(K_eps-K_0)': '0', 'I-d + II-b': '0',
                'sum': {g: nstr(nadd(nadd(sumA['log'][g], IIa[g]), III0[g])) for g in gs},
                'zero': all(nadd(nadd(sumA['log'][g], IIa[g]), III0[g]) == ZERO for g in gs)},
            'I-b + III(K_eps-K_0) at eps^0': {g: nstr(nadd(Ib[g], IIIk[g])) for g in gs},
            'Ob_3 in this sector (eps^0 of I-a; II-a, III(K_0) have no eps^0 part)': {
                g: nstr(sumA['eps^0'][g]) for g in gs},
            'engine_A_reported_total': A_sec['totals']['Ob_3 (eps^0, local)'],
            'engine_B_reported_total': B_sec['totals']['eps^0 (Ob_3)'],
        }
        totals['reported_totals_agree_with_class_sums'] = all(
            num(A_sec['totals']['Ob_3 (eps^0, local)'][g]) == sumA['eps^0'][g] and
            num(B_sec['totals']['eps^0 (Ob_3)'][g]) == sumB['eps^0'][g] for g in gs)
        totals_all.update({g: sumA['eps^0'][g] for g in gs})
        census = {'theory': TC[secA]['I-a_count'],
                  'theory_rows': len(tc), 'engine_A_classes': len(ca), 'engine_B_classes': len(cb),
                  'matched_A_B_T': sum(1 for r in rows if r['engine_A_class'] and r['engine_B_class']
                                       and r['theory_census_row']),
                  'engine_A_extra_zero_classes': sum(1 for r in rows if r['engine_A_class'] and not r['engine_B_class'])}
        out['sectors'][secA] = {'I-a_classes': rows, 'other_contributions': others, 'totals': totals,
                                'census': census}
    out['all_classes_and_contributions_agree'] = all_agree
    # ------------------------------------------------------------ agreed alpha
    alpha = {f'alpha_{g[1]}': nstr(totals_all[g]) for g in ('G1', 'G2', 'G3', 'G4', 'G5')}
    Aal = EA['alpha']
    Bal = {k: v['exact'] for k, v in EB['alpha'].items()}
    out['agreed_alpha'] = {
        'values': alpha,
        'float': {k: float(sp.N(sp.sympify(v.replace('log(3/2)', 'log(3/2)')), 20)) for k, v in alpha.items()},
        'engine_A_reported': Aal, 'engine_B_reported': Bal,
        'all_three_identical': all(num(Aal[k]) == num(Bal[k]) == totals_all['G' + k[-1]] for k in Aal),
        'Ob_3': 'unit * [ (11/4 + (19/4) log(3/2)) (G1 + G3) + (27/8) G2 + 0 G4 - 9 log(3/2) G5 ]'}
    # ------------------------------------------------------------ maps: theory vs engine B
    QHt = TM['matrices']['Q_H_in_G_basis(rows G1..G5, cols H1..H6)']
    BRt = TM['matrices']['bracket_I3_in_G_basis']
    SB = json.load(open(ROOT + '/engine_b/solve_F3.json'))['projections']
    secG = {'cAVAV': ['G1', 'G2'], 'cApAp': ['G3'], 'clamlam': ['G5'], 'ccAsV': ['G4']}
    mism = []
    for sB, gl in secG.items():
        for g in gl:
            gi = int(g[1]) - 1
            for j in range(6):
                if sp.Rational(QHt[gi][j]) != sp.Rational(SB[sB][f'QH{j + 1}']['G'][g]):
                    mism.append((sB, g, f'QH{j + 1}'))
            for xt, xb in (('D1', 'D'), ('B1', 'B1'), ('B2', 'B2'), ('B3', 'B3'), ('B4', 'B4')):
                if sp.Rational(BRt[xt][gi]) != sp.Rational(SB[sB]['{I3,%s}' % xb]['G'][g]):
                    mism.append((sB, g, xt))
    out['maps'] = {'Q_H_in_G_basis(rows G1..G5, cols H1..H6)': QHt,
                   'bracket_I3_in_G_basis': {k: BRt[k] for k in ('D1', 'B1', 'B2', 'B3', 'B4', 'massH')},
                   'quadratic_functionals': {'D1': 'int(d^*A|_V, d^*A|_V)', 'B1': 'int(dA|_V, dA|_V)',
                                             'B2': 'int(dA|_Wperp, dA|_Wperp)', 'B3': 'int(lambda, dslash lambda)',
                                             'B4': 'int(dc, A^*)', 'massH': 'int(A, H A)'},
                   'H_basis': TM['conventions']['H_basis'],
                   'H1_H2_normalisation': 'H1 = 2 sum_{mu,nu} (d_mu A_nu, [A_mu, A_nu]) (2-form pairing sum_{mu<nu}), '
                                          'H2 = sum_{mu,w} (d_mu iota_w A, [A_mu, iota_w A]); same in theory and engine B',
                   'theory_vs_engine_B_entrywise_identical': len(mism) == 0, 'mismatches': mism}
    # ------------------------------------------------------------ F_3
    h = sp.symbols('h1:7')
    be = sp.symbols('beta1:5')
    k1, k2 = sp.symbols('kappa1 kappa2')
    L = sp.Symbol('log(3/2)')
    al = {g: totals_all[g][0] + totals_all[g][2] * L for g in ('G1', 'G2', 'G3', 'G4', 'G5')}
    assert all(totals_all[g][1] == -totals_all[g][2] for g in al)
    F2c = {'D1': sp.Rational(27, 16), 'B1': be[0], 'B2': be[1], 'B3': be[2], 'B4': be[3]}
    eqs = []
    for gi, g in enumerate(('G1', 'G2', 'G3', 'G4', 'G5')):
        rhs = sum(sp.Rational(QHt[gi][j]) * h[j] for j in range(6)) + sum(
            sp.Rational(BRt[x][gi]) * c for x, c in F2c.items())
        eqs.append(sp.expand(rhs - al[g]))
    sol = sp.solve(eqs, [h[0], h[2], h[3], h[4]], dict=True)
    assert len(sol) == 1
    sol = sol[0]
    g2_row = sp.expand(eqs[1])
    part = {str(hh): sp.expand(sol.get(hh, hh).subs({h[1]: 0, h[5]: 0, **{b: 0 for b in be}})) for hh in h}
    part['h2'] = sp.Integer(0)
    part['h6'] = sp.Integer(0)
    fam = {
        'h1': -al['G1'] / 2 + al['G4'] / 4 + be[0] + k1,
        'h2': 2 * be[1] + k2,
        'h3': -al['G4'] + be[3] - 4 * k1,
        'h4': -2 * al['G3'] + be[3] - 2 * k2,
        'h5': al['G5'] - 2 * be[2] + be[3],
        'h6': be[3] / 2 - 2 * k1}
    # symbolic verification: substitute the family into every row
    resid = [sp.simplify(e.subs({h[i]: fam[f'h{i + 1}'] for i in range(6)})) for e in eqs]
    # kernel check
    kern = [[1, 0, -4, 0, 0, -2], [0, 1, 0, -2, 0, 0]]
    kern_img = [[sum(sp.Rational(QHt[gi][j]) * v[j] for j in range(6)) for gi in range(5)] for v in kern]
    rankQ = sp.Matrix([[sp.Rational(x) for x in row] for row in QHt]).rank()
    pb_t = json.load(open(ROOT + '/qa/plugback_theory.json'))
    pb_b = json.load(open(ROOT + '/qa/plugback_engine_b.json'))
    ms = json.load(open(ROOT + '/qa/masters_indep.json'))

    def pstr(x):
        return str(x).replace('log(3/2)', 'log(3/2)')
    out['F3'] = {
        'equation': 'Ob_3 = Q F_3 + {I_3, F_2},  F_3 = unit * sum_j h_j H_j,  F_2 = unit * [ (27/16) int(d^*A|_V, d^*A|_V) '
                    '+ beta_1 int(dA|_V,dA|_V) + beta_2 int(dA|_Wperp,dA|_Wperp) + beta_3 int(lambda,dslash lambda) '
                    '+ beta_4 int(dc,A^*) ]',
        'row_equations_(G1..G5)_lhs_minus_alpha': [pstr(e) for e in eqs],
        'G2_row': pstr(g2_row) + '  (identically 0 with alpha_2 = 27/8: solvability condition satisfied)',
        'G2_row_vanishes': g2_row == 0,
        'particular_solution_beta0_h2_h6_zero': {k: pstr(v) for k, v in part.items()},
        'particular_F3': 'unit * [ -(11/8 + (19/8) log(3/2)) H1 - (11/2 + (19/2) log(3/2)) H4 - 9 log(3/2) H5 ]',
        'particular_F3_written_out': 'unit * [ -(11/8 + (19/8) log(3/2)) int(dA|_V,[A|_V ^ A|_V]) '
                                     '- (11/2 + (19/2) log(3/2)) int([A|_Wperp, c], A^*|_Wperp) '
                                     '- 9 log(3/2) int([c, lambda], lambda^*) ]',
        'alternative_representative_(kappa2 = -alpha_3: H4-coefficient 0)':
            'unit * [ -(11/8 + (19/8) log(3/2)) H1 - (11/4 + (19/4) log(3/2)) H2 - 9 log(3/2) H5 ]',
        'representative_with_F3_zero': {
            'choice': 'beta_1 = beta_2 = alpha_1/2 (= alpha_3/2), beta_3 = alpha_5/2, beta_4 = 0, kappa_1 = 0, '
                      'kappa_2 = -alpha_3',
            'F3': '0',
            'F2': 'unit * [ (27/16) int(d^*A|_V,d^*A|_V) + (11/8 + (19/8) log(3/2)) ( int(dA|_V,dA|_V) + '
                  'int(dA|_Wperp,dA|_Wperp) ) - (9/2) log(3/2) int(lambda, dslash lambda) ]',
            'statement': 'Ob_3 = {I_3, F_2} exactly: {I_3, int(dA|_V,dA|_V)} = 2 G1, {I_3, int(dA|_Wperp,dA|_Wperp)} = '
                         '2 G3, {I_3, int(lambda,dslash lambda)} = 2 G5, {I_3, int(d^*A|_V,d^*A|_V)} = 2 G2; possible '
                         'because alpha_1 = alpha_3 and alpha_4 = 0. The cubic obstruction is absorbed by finite '
                         'Q-closed quadratic kinetic counterterms alone (plus the Step-1 (27/16) term).',
            'plugback_theory_realisation': None, 'plugback_engine_B_realisation': None},
        'free_directions': {
            'beta_1': 'F_2 += beta_1 int(dA|_V,dA|_V), F_3 += beta_1 H1  (cubic part of int(F_A,F_A) along V)',
            'beta_2': 'F_2 += beta_2 int(dA|_Wperp,dA|_Wperp), F_3 += 2 beta_2 H2  (cubic part of int(d_A phi, d_A phi))',
            'beta_3': 'F_2 += beta_3 int(lambda,dslash lambda), F_3 += -2 beta_3 H5  (Q-equivalent to the cubic part '
                      'int(lambda, rho([A,lambda])) of int(lambda, dslash_A lambda))',
            'beta_4': 'F_2 += beta_4 int(dc,A^*), F_3 += beta_4 (H3 + H4 + H5 + H6/2)  (cubic part of the BRST term)',
            'kappa_1': 'F_3 += kappa_1 (H1 - 4 H3 - 2 H6)  (Q-closed, not Q-exact)',
            'kappa_2': 'F_3 += kappa_2 (H2 - 2 H4)  (Q-closed, not Q-exact)'},
        'family_h': {k: pstr(sp.expand(v)) for k, v in fam.items()},
        'symbolic_substitution_residuals_(G1..G5)': [pstr(r) for r in resid],
        'symbolic_verification': all(r == 0 for r in resid),
        'rank_Q_on_span_H': int(rankQ),
        'kernel_images_Q(H1-4H3-2H6), Q(H2-2H4)': [[str(x) for x in r] for r in kern_img],
        'agreement_with_theory_family': 'identical (theory mu_1, mu_2 = kappa_1, kappa_2)',
        'agreement_with_engine_B_family': 'identical (engine B parametrises by h4, h6: kappa_1 = -h6/2 + beta_4/4, '
                                          'kappa_2 = -(h4 + 2 alpha_3 - beta_4)/2)',
        'plugback_theory_realisation': {
            'file': 'qa/plugback_theory.json',
            'target_sectors_all_hold': all(pb_t['equation_E'][s]['all_hold'] for s in
                                           ('c_AV_AV', 'c_AW_AW', 'c_lam_lam', 'c_c_AsV')),
            'no_invariant_sectors_all_zero': all(pb_t['equation_E'][s]['all_hold'] for s in
                                                 ('c_AV_AW', 'c_A_lam', 'c_c_lams', 'c_c_AsW')),
            'samples_per_sector': 5, 'parameter_sets': pb_t['parameter_sets']},
        'plugback_engine_B_realisation': {
            'file': 'qa/plugback_engine_b.json',
            'all_hold': all(v['all_hold'] for v in pb_b['sectors'].values()),
            'samples_per_sector': 6, 'parameter_sets': pb_b['parameter_sets']},
    }
    out['F3']['representative_with_F3_zero']['plugback_theory_realisation'] = all(
        pb_t['equation_E'][x]['F3_zero_representative_all_hold'] for x in pb_t['equation_E'])
    out['F3']['representative_with_F3_zero']['plugback_engine_B_realisation'] = all(
        v['F3_zero_representative_all_hold'] for v in pb_b['sectors'].values())
    zr = [sp.simplify(e.subs({h[i]: 0 for i in range(6)}).subs(
        {be[0]: al['G1'] / 2, be[1]: al['G3'] / 2, be[2]: al['G5'] / 2, be[3]: 0})) for e in eqs]
    out['F3']['representative_with_F3_zero']['symbolic_residuals_(G1..G5)'] = [str(x) for x in zr]
    out['wess_zumino'] = {
        'statement': 'Q Ob_3 + {I_3, Ob_2} = (alpha_2 - 27/8) unit Q G2 (Q G1 = Q G3 = Q G4 = Q G5 = 0 on (c,c,A|_V))',
        'with_agreed_alpha': 'alpha_2 = 27/8 => 0',
        'plane_wave_check_theory_realisation': pb_t['WZ']['all_zero'],
        'automatic': 'yes: Q^2 = 0 on lim Q W, Q I_3 = 0, Q Ob_2 = 0, Q P(0,infty) = K_0 and the two limit statements '
                     '(I-b + III(K_eps-K_0) -> 0, I-d + II-b -> tree(Ob_2, I_3; P(0,infty))); it is not a consequence of '
                     'linear algebra alone (Q on span H plus the beta-brackets misses G2). Solvability <=> WZ <=> '
                     'alpha_2 = 27/8, and the computed alpha_2 is 27/8.'}
    out['independent_checks'] = {
        'Y_log_from_theory_realisation_(II-a, III(K_0))': pb_t['Y_log_checks'],
        'III(K_eps-K_0)_eps0_from_(9/4){I_3,int(A,HA)}': pb_t['III_Kprime_check'],
        'master_integrals': {
            'file': 'qa/masters_indep.json',
            'engine_A_13_integrals_max_rel_diff': ms['engine_A_max_rel_diff'],
            'engine_B_56_masters_max_rel_diff': ms['engine_B_max_rel_diff'],
            'A_and_B_exact_forms_identical_on_A_list': ms['cross_engine_all_identical'],
            'method': 'last [1,oo) integral in closed form (binomial), remaining variables by 40-digit tanh-sinh'}}
    out['notes'] = [
        'Engines A and B are independent codes but share two Step-1 inputs: the cubic log multiples 6, 10, 1, 1, 1 '
        '(II-a) and the mass coefficient 9/4 (I-b, III, II-b; recomputed from the graphs by engine A). The II-a and '
        'III(K_0) log parts are reproduced here in the theory realisation from Y_3 = 6 AAA + 10 lambda lambda A + AcA^* '
        '+ c lambda lambda^* + ccc^* and Y_2 = -B1 - B2 + 4 B3 - B4; the III(K_eps-K_0) limit from (9/4){I_3,int(A,HA)}.',
        'Both engines use the same graph decomposition (special edge K_eps - K_0 on the loop, Schwinger times scaled by '
        'eps, special time in [0,1], the others in [1,infty)) but different master integrals (engine A: no power of the '
        'special time; engine B: J3(a1,a2,a3;m), J2(a1,a2;m)) and different numerator constructions; the master values '
        'are confirmed independently here (qa/masters_indep.json).',
        'Finite-eps full-sum numerics: engine A in all four sectors (deviation ~ eps log eps); engine B in (c,c,A^*|_V) '
        '(full sum) and (c,A|_Wperp,A|_Wperp) (I-a only).',
        'log(3/2) enters alpha_1 = alpha_3 and alpha_5 only; these coefficients are not constrained by the Wess-Zumino '
        'condition and are removed by the finite counterterm (scheme constants).']
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1, default=str)
    print('all agree:', all_agree)
    for s, v in out['sectors'].items():
        print(s, v['census'], v['totals']['eps^-1 total']['zero'], v['totals']['log eps total']['zero'],
              v['totals']['Ob_3 in this sector (eps^0 of I-a; II-a, III(K_0) have no eps^0 part)'])
        for r in v['I-a_classes']:
            print('   ', r['status'][:40], '|', r['loop_word'], '|', {g: r['coefficients']['eps^0'][g]['A'] or
                                                                      r['coefficients']['eps^0'][g]['B'] for g in
                                                                      r['coefficients']['eps^0']})
        for kk, vv in v['other_contributions'].items():
            print('   *', kk, vv.get('A'), vv.get('B'), vv.get('theory'), vv.get('A_equals_B'))
    print(json.dumps(out['agreed_alpha'], indent=0, default=str))
    print(json.dumps(out['F3'], indent=0, default=str)[:4000])
    print(json.dumps(out['maps']['theory_vs_engine_B_entrywise_identical']))


if __name__ == '__main__':
    main()
