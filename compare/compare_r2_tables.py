"""Round-2 comparison: engine A vs engine B vs FeynArts census (T1, T2, T3).
Writes compare/compare_r2_tables.json. Exact rationals only (sympy parsing of the engines' strings).
Run:  python3 -I compare/compare_r2_tables.py   (from anywhere)"""
import json, re, collections
import sympy as sp

R = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
A = json.load(open(f'{R}/engine_a/results_r2.json'))
B = json.load(open(f'{R}/engine_b/results_r2.json'))
M = json.load(open(f'{R}/diagrams/manifest.json'))

S = lambda x: sp.nsimplify(sp.sympify(str(x).replace('^', '**')))

SETS = [  # (A key, B key, FeynArts job, quartic?)
    ('(A|_V, A|_V, A|_V)', '(A|_V, A|_V, A|_V)', 'c3_AVAVAV', False),
    ('(A|_V, A|_Wperp, A|_Wperp)', '(A|_V, A|_Wperp, A|_Wperp)', 'c3_AVAWAW', False),
    ('(lambda, lambda, A|_V)', '(lambda, lambda, A|_V)', 'c3_lamlamAV', False),
    ('(lambda, lambda, A|_Wperp)', '(lambda, lambda, A|_Wperp)', 'c3_lamlamAW', False),
    ('(c, A|_V, A^*|_V)', '(c, A|_V, A*|_V)', 'c3_cAVAsV', False),
    ('(c, A|_Wperp, A^*|_Wperp)', '(c, A|_Wperp, A*|_Wperp)', 'c3_cAWAsW', False),
    ('(c, lambda, lambda^*)', '(c, lambda, lambda*)', 'c3_clamlams', False),
    ('(c, c, c^*)', '(c, c, c*)', 'c3_cccs', False),
    ('(A|_V)^4', '(A|_V)^4', 'q4_AVAVAVAV', True),
    ('(A|_V)^2 (A|_Wperp)^2', '(A|_V)^2 (A|_Wperp)^2', 'q4_AVAVAWAW', True),
    ('(A|_Wperp)^4', '(A|_Wperp)^4', 'q4_AWAWAWAW', True)]

# ---------------------------------------------------------------- vertex typing (FeynArts model names)
EXT = {'A|_V': 'V', 'A|_Wperp': 'W', 'A^*|_V': 'sV', 'A*|_V': 'sV', 'A^*|_Wperp': 'sW', 'A*|_Wperp': 'sW',
       'lambda': 'lam', 'lambda^*': 'lams', 'lambda*': 'lams', 'c': 'c', 'c^*': 'cs', 'c*': 'cs'}
LINE_A = {'A|_V-A|_V': 'V', 'A|_Wperp-A|_Wperp': 'W', 'lambda-lambda': 'lam', 'c-A^*|_V': 'gh'}
LINE_B = {'A|_V': 'V', 'A|_Wperp': 'W', 'lambda': 'lam', 'c-A*': 'gh'}
FA_LINE = {'V': 'P_AA(V)', 'W': 'P_AA(W)', 'lam': 'P_lamlam', 'gh': 'P_cA*'}

def vtype(name, ext, lines):
    ext = [EXT[e] for e in ext]
    bos = [x for x in ext + lines if x in ('V', 'W')]
    if name in ('AAA', 'AAAA'):
        assert len(bos) == len(ext) + len(lines) == {'AAA': 3, 'AAAA': 4}[name], (name, ext, lines)
        t = ''.join(sorted(bos))
        assert t in ('VVV', 'VWW', 'VVVV', 'VVWW', 'WWWW'), t
        return ('g3[%s]' if name == 'AAA' else 'g4[%s]') % t
    if name in ('llA', 'lambda lambda A'):
        assert len(bos) == 1 and (ext + lines).count('lam') == 2, (name, ext, lines)
        return 'y[%s]' % bos[0]
    if name in ('AcAs', 'AcA*'):
        X = set(x[-1] for x in ext + lines if x in ('V', 'W', 'sV', 'sW'))
        assert len(X) == 1, (name, ext, lines)
        return 'gh[%s]' % X.pop()
    if name in ('cllS', 'c lambda lambda*'):
        return 'clamlams'
    if name in ('cccS', 'ccc*'):
        return 'cccs'
    raise ValueError(name)

def key_from_cycle(verts, lines):
    """verts: list of (name, ext list); lines: list of line kinds; line i joins vertex i and vertex i+1 (cyclic)."""
    n = len(verts)
    vt = []
    for i, (nm, ext) in enumerate(verts):
        vt.append(vtype(nm, ext, [lines[(i - 1) % n], lines[i]]))
    return (tuple(sorted(vt)), tuple(sorted(FA_LINE[l] for l in lines)))

def key_A(g):
    verts = []
    for v in g['vertices']:
        m = re.match(r'^(\w+)\[(.*)\]$', v)
        verts.append((m.group(1), [x.strip() for x in m.group(2).split(',')]))
    return key_from_cycle(verts, [LINE_A[l] for l in g['internal_lines']])

def key_B(g):
    body = g['id'].split(':', 1)[1]
    verts = []; lines = []
    for seg in body.split(' | '):
        m = re.match(r'^\s*(.*?)\((.*)\) -(.*)-\s*$', seg)
        verts.append((m.group(1).strip(), [x.strip() for x in m.group(2).split(',')]))
        lines.append(LINE_B[m.group(3).strip()])
    return key_from_cycle(verts, lines)

def key_FA(fid):
    _, vt, ln = fid.split(':')
    ln = ln.split('#')[0]
    return (tuple(sorted(vt.split(','))), tuple(sorted(ln.split(','))))


def sig_from_cycle(verts, lines):
    n = len(verts)
    seq = []
    for i, (nm, ext) in enumerate(verts):
        seq.append((vtype(nm, ext, [lines[(i - 1) % n], lines[i]]), tuple(sorted(EXT[e] for e in ext))))
    cands = []
    for refl in (False, True):
        if refl:
            vs = seq[::-1]; ls = [lines[(n - 2 - j) % n] for j in range(n)]
        else:
            vs, ls = seq, lines
        for r in range(n):
            cands.append(tuple(zip(vs[r:] + vs[:r], ls[r:] + ls[:r])))
    return min(cands)

def sig_A(g):
    verts = []
    for v in g['vertices']:
        m = re.match(r'^(\w+)\[(.*)\]$', v)
        verts.append((m.group(1), [x.strip() for x in m.group(2).split(',')]))
    return sig_from_cycle(verts, [LINE_A[l] for l in g['internal_lines']])

def sig_B(g):
    body = g['id'].split(':', 1)[1]
    verts = []; lines = []
    for seg in body.split(' | '):
        m = re.match(r'^\s*(.*?)\((.*)\) -(.*)-\s*$', seg)
        verts.append((m.group(1).strip(), [x.strip() for x in m.group(2).split(',')]))
        lines.append(LINE_B[m.group(3).strip()])
    return sig_from_cycle(verts, lines)

# ---------------------------------------------------------------- colour/Lorentz bases (quartic)
# engine A: su(2), X=(12)(34), Y=(13)(24), Z=(14)(23) colour pairings, x,y,z Lorentz; units u (tr_g(ad ad)/(,) = -2)
# engine B: general g, b01 = ([X0,X1],[X2,X3]) (units u), b03 = ([X0,X3],[X1,X2]) (units u),
#           S = (1/24) sum_sigma tr_g(ad ad ad ad) (units hbar/(16 pi^2))
# su(2) specialisation: b01 -> Y - Z, b03 -> X - Y, S -> (2/3)(X+Y+Z), hbar/(16pi^2) = -u/2
#   => X = b03 - S/3, Y = b01 - b03 - S/3, Z = -b01 - S/3 ; inverse: S = -(X+Y+Z), b01 = (X+Y-2Z)/3, b03 = (2X-Y-Z)/3
COLB = {'([X_0,X_1],[X_2,X_3])': 'b01', '([X_0,X_3],[X_1,X_2])': 'b03', '(1/24) sum_sigma tr_g(ad ad ad ad)': 'S'}
def lorB(s):
    s = s.replace(' ', '')
    m = re.findall(r'\((\w)_(\d)\.(\w)_(\d)\)', s)
    if s.startswith('eps('):
        return 'eps'
    idx = tuple(sorted([(int(m[0][1]), int(m[0][3])), (int(m[1][1]), int(m[1][3]))]))
    return {((0, 1), (2, 3)): 'x', ((0, 2), (1, 3)): 'y', ((0, 3), (1, 2)): 'z'}[idx]

def A_to_gen(d):
    """engine A structure dict {'Xx': c, ...} -> {(lor, col): c} in the general basis"""
    out = collections.defaultdict(lambda: sp.Integer(0))
    for lor in 'xyz':
        X, Y, Z = [S(d.get(c + lor, 0)) for c in 'XYZ']
        if X == Y == Z == 0:
            continue
        out[(lor, 'S')] += -(X + Y + Z)
        out[(lor, 'b01')] += (X + Y - 2 * Z) / 3
        out[(lor, 'b03')] += (2 * X - Y - Z) / 3
    return {k: sp.nsimplify(v) for k, v in out.items() if v != 0}

def B_to_gen(d):
    out = collections.defaultdict(lambda: sp.Integer(0))
    for k, v in d.items():
        col, lor = k.split(' x ', 1)
        out[(lorB(lor), COLB[col.strip()])] += S(v)
    return {k: sp.nsimplify(v) for k, v in out.items() if v != 0}

def scal(d, k):
    return {kk: sp.nsimplify(k * v) for kk, v in d.items() if k * v != 0}

def dstr(d):
    return {f'{col} x {lor}': str(v) for (lor, col), v in sorted(d.items())}

COLNAME = {'b01': '([X_0,X_1],[X_2,X_3]) [units u]', 'b03': '([X_0,X_3],[X_1,X_2]) [units u]',
           'S': '(1/24) sum over orderings of tr_g(ad ad ad ad) [units hbar/(16 pi^2)]'}

# ---------------------------------------------------------------- FeynArts census per job
fa = collections.defaultdict(lambda: collections.Counter())
fa_ids = collections.defaultdict(lambda: collections.defaultdict(list))
for d in M['diagrams']:
    if d['feynarts_job'] and (d['feynarts_job'].startswith('c3_') or d['feynarts_job'].startswith('q4_')):
        k = key_FA(d['id'])
        fa[d['feynarts_job']][k] += 1
        fa_ids[d['feynarts_job']][k].append(d['id'])

out = {'T1': {}, 'T2': {}, 'T3': {}}
disagreements = []
for akey, bkey, job, quart in SETS:
    sa = A['T1']['sets'][akey]; sb = B['T1']['external_sets'][bkey]
    # trees (quartic): general basis from engine A's tree_in_basis; check against engine B total/8
    if quart:
        treeG = A_to_gen(sa['tree_in_basis'])
        totB = B_to_gen({k: v for k, v in sb['total_structure_fit'].items() if 'eps(' not in k})
        epsB = {k: v for k, v in sb['total_structure_fit'].items() if 'eps(' in k and S(v) != 0}
        tree_check = (scal(treeG, 8) == totB) and not epsB
    # engine A per key
    A_by = collections.defaultdict(list)
    for g in sa['graphs']:
        A_by[key_A(g)].append(g)
    B_by = collections.defaultdict(list)
    for g in sb['graphs']:
        B_by[key_B(g)].append(g)
    keys = sorted(set(fa[job]) | set(A_by) | set(B_by))
    rows = []
    for k in keys:
        fa_n = fa[job].get(k, 0)
        a_n = 0
        for g in A_by.get(k, []):
            orient = 2 if all(l == 'c-A^*|_V' for l in g['internal_lines']) else 1
            a_n += g['graphs_in_class'] * orient
        # values
        if not quart:
            def a_val(g):
                m = re.match(r'^(.*) x tree vertex$', g['log_eps']); return S(m.group(1))
            def b_val(g):
                m = re.match(r'^(.*) x tree$', g['log_eps']); return S(m.group(1))
            av = sum((a_val(g) for g in A_by.get(k, [])), sp.Integer(0)) if k in A_by else None
            bv = sum((b_val(g) for g in B_by.get(k, [])), sp.Integer(0)) if k in B_by else sp.Integer(0)
            eq = (av == bv) if av is not None else None
            avs, bvs = (str(av) if av is not None else None), str(bv)
        else:
            def a_val(g):
                le = g['log_eps']
                if isinstance(le, str):
                    m = re.match(r'^(.*) x tree vertex$', le); return scal(treeG, S(m.group(1)))
                return A_to_gen(le['structures'])
            def b_val(g):
                le = g['log_eps']
                if isinstance(le, str):
                    m = re.match(r'^(.*) x tree$', le); return scal(treeG, S(m.group(1)))
                return B_to_gen(le)
            av = None
            if k in A_by:
                av = {}
                for g in A_by[k]:
                    for kk, vv in a_val(g).items():
                        av[kk] = sp.nsimplify(av.get(kk, 0) + vv)
                av = {kk: vv for kk, vv in av.items() if vv != 0}
            bv = {}
            for g in B_by.get(k, []):
                for kk, vv in b_val(g).items():
                    bv[kk] = sp.nsimplify(bv.get(kk, 0) + vv)
            bv = {kk: vv for kk, vv in bv.items() if vv != 0}
            eq = (av == bv) if av is not None else None
            avs = dstr(av) if av is not None else None
            bvs = dstr(bv)
        row = dict(feynarts_key='vertices ' + ','.join(k[0]) + ' ; lines ' + ','.join(k[1]),
                   feynarts_ids=fa_ids[job].get(k, []), feynarts_count=fa_n,
                   engine_A_ids=[g['id'] for g in A_by.get(k, [])], engine_A_count=a_n,
                   engine_B_ids=[g['id'] for g in B_by.get(k, [])],
                   engine_A_log_eps=avs, engine_B_log_eps=(bvs if k in B_by else ('absent (no nonzero word)' if not quart else 'absent (no nonzero word)')),
                   engine_B_value_used=(bvs if not quart else bvs),
                   engine_A_eps_inv=[g['eps_inv'] for g in A_by.get(k, [])], engine_A_eps_inv_half=[g.get('eps_inv_half') for g in A_by.get(k, [])],
                   engine_B_eps_inv=[g['eps_inv'] for g in B_by.get(k, [])], engine_B_eps_inv_half=[g.get('eps_inv_half') for g in B_by.get(k, [])],
                   census_A_equals_FeynArts=(a_n == fa_n), present_in_B=(k in B_by), values_equal=eq)
        rows.append(row)
        if a_n != fa_n:
            disagreements.append(dict(set=akey, kind='census A vs FeynArts', key=row['feynarts_key'], A=a_n, FA=fa_n))
        if k not in B_by:
            disagreements.append(dict(set=akey, kind='census: graph class absent from engine B (B lists only classes with a nonzero word)',
                                      key=row['feynarts_key'], feynarts_ids=row['feynarts_ids'], engine_A_value=avs))
        if eq is False:
            disagreements.append(dict(set=akey, kind='per-graph log eps value', key=row['feynarts_key'], A=avs, B=bvs))
        for e in row['engine_A_eps_inv'] + row['engine_B_eps_inv'] + [x for x in row['engine_A_eps_inv_half'] + row['engine_B_eps_inv_half'] if x is not None]:
            if S(e) != 0:
                disagreements.append(dict(set=akey, kind='nonzero eps^-1 or eps^-1/2', key=row['feynarts_key'], value=e))
    # finer check: class by class (cyclic signature up to rotation/reflection)
    sigA = {sig_A(g): g for g in sa['graphs']}; sigB = {sig_B(g): g for g in sb['graphs']}
    assert len(sigA) == len(sa['graphs']) and len(sigB) == len(sb['graphs'])
    per_class = []
    for sg in sorted(set(sigA) | set(sigB), key=str):
        gA = sigA.get(sg); gB = sigB.get(sg)
        if not quart:
            va = S(re.match(r'^(.*) x tree vertex$', gA['log_eps']).group(1)) if gA else None
            vb = S(re.match(r'^(.*) x tree$', gB['log_eps']).group(1)) if gB else sp.Integer(0)
            eqc = (va == vb) if gA else None; vas = str(va) if gA else None; vbs = str(vb)
        else:
            def gen(g, side):
                le = g['log_eps']
                if isinstance(le, str):
                    return scal(treeG, S(le.split(' x ')[0]))
                return A_to_gen(le['structures']) if side == 'A' else B_to_gen(le)
            va = gen(gA, 'A') if gA else None
            vb = gen(gB, 'B') if gB else {}
            eqc = (va == vb) if gA else None; vas = dstr(va) if gA else None; vbs = dstr(vb)
        per_class.append(dict(engine_A_id=gA['id'] if gA else None, engine_B_id=gB['id'] if gB else None,
                              engine_A=vas, engine_B=(vbs if gB else 'absent'), equal=eqc))
        if eqc is False:
            disagreements.append(dict(set=akey, kind='per-class log eps value', A=vas, B=vbs))
        if gA is None:
            disagreements.append(dict(set=akey, kind='class only in engine B', B=gB['id']))
    # totals
    if not quart:
        ta = S(sa['total']['multiple_of_tree'])
        tb = S(re.match(r'^(\S+) x tree', sb['total']).group(1))
        tot = dict(A=str(ta), B=str(tb), equal=(ta == tb))
        sum_rows = sum((S(r['engine_A_log_eps']) for r in rows if r['engine_A_log_eps'] is not None), sp.Integer(0))
        tot['sum_of_per_graph_A'] = str(sum_rows)
        tot['sum_of_per_graph_B'] = str(sum((S(r['engine_B_value_used']) for r in rows), sp.Integer(0)))
    else:
        ta = S(sa['total']['multiple_of_tree'])
        tb = S(re.match(r'^(\S+) x tree', sb['total']).group(1))
        sumA = {}; sumB = {}
        for r in rows:
            for kk, vv in (r['engine_A_log_eps'] or {}).items():
                sumA[kk] = sp.nsimplify(sumA.get(kk, 0) + S(vv))
            for kk, vv in r['engine_B_value_used'].items():
                sumB[kk] = sp.nsimplify(sumB.get(kk, 0) + S(vv))
        sumA = {kk: str(vv) for kk, vv in sumA.items() if vv != 0}; sumB = {kk: str(vv) for kk, vv in sumB.items() if vv != 0}
        tot = dict(A=str(ta), B=str(tb), equal=(ta == tb), sum_of_per_graph_A_general_basis=sumA,
                   sum_of_per_graph_B_general_basis=sumB, sums_equal=(sumA == sumB),
                   eight_times_tree_general_basis=dstr(scal(treeG, 8)), sums_equal_8_tree=(sumA == dstr(scal(treeG, 8))),
                   tree_basis_consistency_A_tree_x8_eq_B_total=tree_check)
    if not tot['equal']:
        disagreements.append(dict(set=akey, kind='total', A=tot['A'], B=tot['B']))
    pa = S(sa['predicted_from_A']); pb = S(re.match(r'^(\S+) x tree', sb['predicted_from_A']).group(1))
    out['T1'][akey] = dict(feynarts_job=job, feynarts_graphs=sum(fa[job].values()),
                           engine_A_classes=len(sa['graphs']), engine_A_graphs=sum(r['engine_A_count'] for r in rows),
                           engine_B_classes=len(sb['graphs']), rows=rows, per_class=per_class, total=tot,
                           predicted_from_A=dict(A=str(pa), B=str(pb), equal=(pa == pb), computed_equals_predicted=(ta == pa and tb == pb)),
                           superficial_degree=dict(A=sa['superficial_degree'], B=sb['superficial_degree']))
    if pa != pb or ta != pa:
        disagreements.append(dict(set=akey, kind='predicted multiple', A=str(pa), B=str(pb)))

# coupling direction / general fit
gfA = A['T1']['test_of_A']['general_fit']; gfB = B['T1']['coupling_direction']['solution_all_sectors']
out['T1_coupling_direction'] = dict(A=gfA[gfA.index('[{'):], B=gfB,
    translation='A weights (w_V,w_P,w_l,w_c) = (alpha_AV, alpha_AWperp, -alpha_lambda, -alpha_c) of B; A beta = B kappa',
    equal=True if ('beta: 0' in gfA and 'wP: 2' in gfA and 'wV: 2' in gfA and 'wc: 1' in gfA and 'wl: 4' in gfA and
                   gfB == "[{alpha_AV: 2, alpha_AWperp: 2, alpha_c: -1, alpha_lambda: -4, kappa: 0}]") else False)

# ---------------------------------------------------------------- T2
T2A = A['T2']; T2B = B['T2']
mapA = {'A|_V - A^*|_V ; P(eps,infty) line A|_V - A|_V': 'QME2-a',
        'A|_Wperp - A^*|_Wperp ; P(eps,infty) line A|_Wperp - A|_Wperp': 'QME2-b',
        'A|_V - A^*|_V ; P(eps,infty) line c - A^*|_V': 'QME2-c',
        'lambda^* - lambda ; P(eps,infty) line lambda - lambda': 'QME2-d',
        'c^* - c ; P(eps,infty) line c - A^*|_V': 'QME2-e'}
FA_T2 = {'QME2-a': 'c,AV:g3[VVV],gh[V]:K-K0_AA*(V),P_AA(V)', 'QME2-b': 'c,AV:g3[VWW],gh[W]:K-K0_AA*(W),P_AA(W)',
         'QME2-c': 'c,AV:gh[V],gh[V]:K-K0_AA*(V),P_cA*', 'QME2-d': 'c,AV:clamlams,y[V]:K-K0_lamlam*,P_lamlam',
         'QME2-e': 'c,AV:cccs,gh[V]:K-K0_cc*,P_cA*'}
fa_qd = [d['id'] for d in M['diagrams'] if d['feynarts_job'] == 'qd_cAV' and d['one_PI']]
fa_qd_red = [d['id'] for d in M['diagrams'] if d['feynarts_job'] in ('qd_cAV', 'qdD_cAV') and not d['one_PI']]
gB = {g['id']: g for g in T2B['graphs']}
num = lambda s: S(re.match(r'^\s*([-+0-9/()]+)', s.replace(' x ', ' ')).group(1).replace('(', '').replace(')', ''))
t2rows = []
for g in T2A['graphs']:
    if 'tree vertex' in g['id']:
        b = gB['QME2-CT']; fid = '(local vertex, no graph)'
        ea = S('9/2'); eb = S('9/2')
        assert '+9/2' in g['eps_inv'] and '+(9/2)' in b['eps_inv']
        fa_ = S(g['finite']); fb_ = S(b['finite'])
        la = S(g['log_eps'].split()[0]); lb = S(b['log_eps'].split()[0])
        bid = 'QME2-CT'
    else:
        tail = g['id'].split('special edge (K_eps - K_0) ')[1]
        bid = mapA[tail]; b = gB[bid]; fid = FA_T2[bid]
        ea = num(g['eps_inv']); eb = num(b['eps_inv'])
        fa_ = num(g['finite']); fb_ = num(b['finite'])
        la = S(g['log_eps']); lb = S(b['log_eps'])
    t2rows.append(dict(engine_A_id=g['id'], engine_B_id=bid, feynarts_id=fid, feynarts_present=(fid in fa_qd or bid == 'QME2-CT'),
                       eps_inv_A=str(ea), eps_inv_B=str(eb), log_A=str(la), log_B=str(lb), finite_A=str(fa_), finite_B=str(fb_),
                       eps_inv_sq_B=b.get('eps_inv_sq', '0'),
                       equal=(ea == eb and fa_ == fb_ and la == lb)))
    if not t2rows[-1]['equal']:
        disagreements.append(dict(kind='T2 per graph', row=t2rows[-1]))
sum_e = sum((S(r['eps_inv_A']) * (1) for r in t2rows), sp.Integer(0))
sum_f = sum((S(r['finite_A']) for r in t2rows), sp.Integer(0))
out['T2'] = dict(
    units='eps^{-1} column: coefficient of u eps^{-1} int(dc, A|_V); finite column: coefficient of u int(d^*dc, d^*A|_V); log column: coefficient of u log eps',
    census=dict(feynarts_1PI=fa_qd, engine_A_special_edge_graphs=len(T2A['graphs']) - 1, engine_B_special_edge_graphs=len(T2B['graphs']) - 1,
                feynarts_1PR_and_K_eps_selfloop=fa_qd_red,
                note_1PR='every listed non-1PI graph contains a self-contraction at one vertex of I; its colour factor is the contraction of a bracket with (-,-)^{-1}, sum_a ([X, t_a], t^a) = 0: both engines state hbar Delta_eps I = 0 by this argument (B: defect_vertex; A: tree-vertex chain); FeynArts manifest: zero'),
    rows=t2rows,
    totals=dict(eps_inv=dict(A=T2A['eps_inv_total'], B=T2B['eps_inv_total'], sum_of_rows=str(sum_e)),
                log=dict(A=T2A['log_total'], B=T2B['log_total']),
                finite=dict(A=T2A['finite_total'], B=T2B['finite_total'], sum_of_rows=str(sum_f)),
                obstruction=dict(A=T2A['obstruction'], B=T2B['obstruction']),
                q_exact=dict(A=T2A['q_exact'], B=T2B['q_exact']),
                finite_counterterm=dict(A=T2A['finite_counterterm'], B=T2B['finite_counterterm']),
                cross_check_route=dict(A=T2A['cross_check_route_b']['result'], B=T2B['cross_check']['agreement']),
                cubic_sector=dict(A=T2A['cubic_sector'], B=T2B['cubic_sector'])),
    totals_equal=(sum_e == 0 and sum_f == S('-27/8')))
if not out['T2']['totals_equal']:
    disagreements.append(dict(kind='T2 totals'))

# ---------------------------------------------------------------- T3
tA = A['T3']['table_units_u_coefficient_multiplying_O']; tB = B['T3']['table']
pairs = [('traceless scalar bilinear (traceless part of int omega (iota_{w1}A, iota_{w2}A))', 'traceless symmetric scalar bilinear (iota_w1 A, iota_w2 A)'),
         ('trace scalar bilinear (sum_w int omega (iota_w A, iota_w A))', 'trace scalar bilinear sum_w (iota_w A, iota_w A)'),
         ('letter iota_w A (w in W_R^perp)', 'letter iota_w A'),
         ('letter A|_V (transverse)', 'letter A|_V (transverse part)'),
         ('letter lambda', 'letter lambda'), ('letter c', 'letter c (gauge ghost, computed in round 2)')]
t3 = []
for ka, kb in pairs:
    a3, a17 = tA[ka]['SW(3)'], tA[ka]['SW(17)']
    b3, b17 = tB[kb]['SW3_Delta_equals_minus_E'], tB[kb]['SW17_Delta_equals_2E']
    t3.append(dict(row=kb, A_SW3=a3, B_SW3=b3, A_SW17=a17, B_SW17=b17, equal=(a3 == b3 and a17 == b17)))
    if not t3[-1]['equal']:
        disagreements.append(dict(kind='T3', row=kb))
out['T3'] = dict(rows=t3, engine_B_only=[{'row': 'letter A|_V (longitudinal part)', 'B': tB['letter A|_V (longitudinal part)']}],
                 psi_weights_A=A['T3']['Psi_weights'], general_rule_B=B['T3']['ghost_letter_computation']['general_rule'],
                 konishi=dict(A=A['T3']['Konishi'], B=B['T3']['konishi_check']))
out['disagreements'] = disagreements
json.dump(out, open(f'{R}/compare/compare_r2_tables.json', 'w'), indent=1, default=str)
print('disagreements:', len(disagreements))
for d in disagreements:
    print(json.dumps(d, default=str)[:600])
for k, v in out['T1'].items():
    print(k, 'FA', v['feynarts_graphs'], 'A', v['engine_A_graphs'], 'Acl', v['engine_A_classes'], 'Bcl', v['engine_B_classes'], 'tot', v['total'].get('equal'), v['total'].get('sums_equal'), v['total'].get('sums_equal_8_tree'), v['total'].get('tree_basis_consistency_A_tree_x8_eq_B_total'))
print('T2 equal rows', all(r['equal'] for r in t2rows), 'totals', out['T2']['totals_equal'])
print('T3 equal', all(r['equal'] for r in t3))
print('coupling', out['T1_coupling_direction']['equal'])
