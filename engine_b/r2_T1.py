"""Engine B round 2, T1: one-loop log eps coefficients of the cubic and quartic vertex functions; test of assumption (A).
Usage: python3 r2_T1.py SETNAME [SETNAME ...]   (writes r2_T1_<SET>.json)"""
import itertools, json, random, sys, time
from fractions import Fraction as Fr
import sympy as sp
import vertices as V
import r2_engine as R

AV = list(range(4)); AP = list(range(4, 10)); LAM = list(range(16)); ONE = [0]
PLUS = list(range(8)); MINUS = list(range(8, 16))      # Sigma = S_+ (x) W (+) S_- (x) W^*  (code basis)
SETS = {
    # name: (legs [(species, comps)], omega, tree term(s), structure family)
    'AV,AV,AV': ([('A', AV), ('A', AV), ('A', AV)], 1, 'AAA', 'AAA_V'),
    'AV,Ap,Ap': ([('A', AV), ('A', AP), ('A', AP)], 1, 'AAA', 'AAA_mixed'),
    'lam,lam,AV': ([('lam', LAM), ('lam', LAM), ('A', AV)], 0, 'lamlamA', 'chiral2'),
    'lam,lam,Ap': ([('lam', LAM), ('lam', LAM), ('A', AP)], 0, 'lamlamA', 'chiral2'),
    'c,AV,As_V': ([('c', ONE), ('A', AV), ('As', AV)], 0, 'AcAs', 'tree'),
    'c,Ap,As_p': ([('c', ONE), ('A', AP), ('As', AP)], 0, 'AcAs', 'tree'),
    'c,lam,lams': ([('c', ONE), ('lam', LAM), ('lams', LAM)], 0, 'clamlams', 'chiral1'),
    'c,c,cs': ([('c', ONE), ('c', ONE), ('cs', ONE)], 0, 'cccs', 'tree'),
    'AV^4': ([('A', AV)] * 4, 0, 'AAAA', 'quartic_V'),
    'AV^2Ap^2': ([('A', AV), ('A', AV), ('A', AP), ('A', AP)], 0, 'AAAA', 'quartic_mixed'),
    'Ap^4': ([('A', AP)] * 4, 0, 'AAAA', 'quartic_p'),
}
PRED_A = {'AAA': 6, 'lamlamA': 10, 'AcAs': 1, 'clamlams': 1, 'cccs': 1, 'AAAA': 8}

def nz():
    return random.choice([-3, -2, -1, 1, 2, 3])

def make_config(setname, seed):
    random.seed(seed)
    legs_spec, omega, _, _ = SETS[setname]
    m = len(legs_spec)
    if omega == 1:
        ps = [[random.randint(-4, 4) for _ in range(4)] for _ in range(m - 1)]
        ps.append([-sum(p[k] for p in ps) for k in range(4)])
    else:
        ps = [[0] * 4 for _ in range(m)]
    legs = []
    for i, (s, comps) in enumerate(legs_spec):
        u = [0] * R.DIM[s]
        for c in comps:
            u[c] = nz()
        legs.append(R.Leg(i, s, comps, u, ps[i]))
    return legs

def sym_leg(l, proj=None):
    u = list(l.u)
    if proj is not None:
        u = [x if k in proj else 0 for k, x in enumerate(u)]
    return dict(sp=l.species, u=[sp.Integer(x) for x in u], p=[sp.Integer(x) for x in l.p])

def tree_vals(term, legs, projs=None):
    sl = [sym_leg(l, None if projs is None else projs[i]) for i, l in enumerate(legs)]
    v = V.vertex(term, sl)
    if term == 'AAAA':
        a = v.get('K12', 0); b = v.get('K13', 0); c = v.get('K14', 0)
        return {'cadj K12': sp.expand(a + b), 'cadj K14': sp.expand(b + c)}
    return {'cadj f': sp.expand(v.get('f', 0))}

def dot(a, b, comps):
    return sum(a[k] * b[k] for k in comps)

def eps4(a, b, c, d):
    tot = 0
    for perm in itertools.permutations(range(4)):
        tot += R.perm_parity(perm) * a[perm[0]] * b[perm[1]] * c[perm[2]] * d[perm[3]]
    return tot

def structures(setname, legs):
    """candidate local structures (complete list for the given external set), as {name: {colour key: value}}"""
    fam = SETS[setname][3]; term = SETS[setname][2]
    u = [l.u for l in legs]; p = [l.p for l in legs]
    S = {}
    if fam == 'AAA_V':
        for k in range(3):
            i, j = [x for x in range(3) if x != k]
            for l in range(2):
                S[f'(a{i}.a{j})(a{k}.p{l})'] = {'cadj f': sp.Integer(dot(u[i], u[j], AV) * dot(u[k], p[l], AV))}
        for l in range(2):
            S[f'eps(a0,a1,a2,p{l})'] = {'cadj f': sp.Integer(eps4(u[0][:4], u[1][:4], u[2][:4], p[l]))}
    elif fam == 'AAA_mixed':
        for l in range(2):
            S[f'(b1.b2)(a0.p{l})'] = {'cadj f': sp.Integer(dot(u[1], u[2], AP) * dot(u[0], p[l], AV))}
        S['(a0.p0)... none'] = None
        del S['(a0.p0)... none']
    elif fam == 'chiral2':
        for c1, P1 in (('+', PLUS), ('-', MINUS)):
            for c2, P2 in (('+', PLUS), ('-', MINUS)):
                S[f'tree[lambda_0 in {c1}, lambda_1 in {c2}]'] = tree_vals(term, legs, [P1, P2, None])
    elif fam == 'chiral1':
        for c1, P1 in (('+', PLUS), ('-', MINUS)):
            S[f'tree[lambda in {c1}]'] = tree_vals(term, legs, [None, P1, None])
    elif fam == 'tree':
        S['tree'] = tree_vals(term, legs)
    elif fam == 'quartic_V':
        for col in ('cadj K12', 'cadj K14', 'str4'):
            for (i, j, k, l) in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
                S[f'{col}:(a{i}.a{j})(a{k}.a{l})'] = {col: sp.Integer(dot(u[i], u[j], AV) * dot(u[k], u[l], AV))}
            S[f'{col}:eps(a0,a1,a2,a3)'] = {col: sp.Integer(eps4(u[0][:4], u[1][:4], u[2][:4], u[3][:4]))}
    elif fam == 'quartic_mixed':
        for col in ('cadj K12', 'cadj K14', 'str4'):
            S[f'{col}:(a0.a1)(b2.b3)'] = {col: sp.Integer(dot(u[0], u[1], AV) * dot(u[2], u[3], AP))}
    elif fam == 'quartic_p':
        for col in ('cadj K12', 'cadj K14', 'str4'):
            for (i, j, k, l) in ((0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)):
                S[f'{col}:(b{i}.b{j})(b{k}.b{l})'] = {col: sp.Integer(dot(u[i], u[j], AP) * dot(u[k], u[l], AP))}
    return S

def to_sp(x):
    return sp.Rational(x[0].numerator, x[0].denominator) + sp.I * sp.Rational(x[1].numerator, x[1].denominator)

def fit(rows):
    """rows: list of (target {colkey: value}, structures {name: {colkey: value}}). Exact solve; returns (solution, residual-free)."""
    names = sorted(rows[0][1].keys())
    cols = sorted({k for t, S in rows for k in list(t.keys()) + [kk for s in S.values() for kk in s.keys()]})
    A = []; b = []
    for t, S in rows:
        for ck in cols:
            A.append([S[nm].get(ck, 0) for nm in names]); b.append(t.get(ck, 0))
    A = sp.Matrix(A); b = sp.Matrix(b)
    try:
        sol, params = A.gauss_jordan_solve(b)
    except ValueError:
        return None, False
    if params.shape[0] > 0:
        sol = sol.subs({p_: 0 for p_ in params})
        undetermined = True
    else:
        undetermined = False
    return {nm: sp.nsimplify(sp.simplify(sol[i])) for i, nm in enumerate(names)}, undetermined

def express_vs_tree(setname, sol):
    """if the fitted structure is a multiple of the tree vertex, return that multiple."""
    fam = SETS[setname][3]
    nzs = {k: v for k, v in sol.items() if v != 0}
    if fam in ('tree',):
        return sol['tree']
    if fam == 'chiral2':
        if setname == 'lam,lam,AV':
            a = sol['tree[lambda_0 in +, lambda_1 in -]']; b = sol['tree[lambda_0 in -, lambda_1 in +]']
            others = [sol['tree[lambda_0 in +, lambda_1 in +]'], sol['tree[lambda_0 in -, lambda_1 in -]']]
        else:
            a = sol['tree[lambda_0 in +, lambda_1 in +]']; b = sol['tree[lambda_0 in -, lambda_1 in -]']
            others = [sol['tree[lambda_0 in +, lambda_1 in -]'], sol['tree[lambda_0 in -, lambda_1 in +]']]
        return a if (a == b) else None
    if fam == 'chiral1':
        a = sol['tree[lambda in +]']; b = sol['tree[lambda in -]']
        return a if a == b else None
    return None

def run_set(setname, nconf=None, verbose=True):
    legs_spec, omega, term, fam = SETS[setname]
    if nconf is None:
        nconf = {'AAA_V': 11, 'AAA_mixed': 4, 'chiral2': 6, 'chiral1': 4, 'tree': 3,
                 'quartic_V': 7, 'quartic_mixed': 3, 'quartic_p': 5}[fam]
    allodd_even = all(R.PAR[s] == 0 for s, _ in legs_spec)
    total_rows = []; tree_rows = []; meta = {}; per_conf = []
    t0 = time.time()
    for ci in range(nconf):
        legs = make_config(setname, 1000 * ci + 7)
        X = R.ExternalSet(setname, legs, omega)
        recs = X.evaluate(use_cyclic=allodd_even)
        S = structures(setname, legs)
        tot = {k: to_sp(v) for k, v in R.assemble(recs, legs).items()}
        total_rows.append((tot, S))
        tree_rows.append((tree_vals(term, legs), S))
        groups = {}
        for r in recs:
            groups.setdefault(R.graph_key(r, legs), []).append(r)
        gv = {}
        for gk, rs in groups.items():
            gv[gk] = {k: to_sp(v) for k, v in R.assemble(rs, legs).items()}
            meta.setdefault(gk, dict(n=rs[0]['n']))
        per_conf.append((gv, S))
        if verbose:
            print(setname, 'config', ci, 'records', len(recs), 'graphs', len(groups), round(time.time() - t0, 1), 's', flush=True)
    per_graph_rows = {gk: [(gv.get(gk, {}), S) for gv, S in per_conf] for gk in meta}
    tree_fit, _ = fit(tree_rows)
    tot_fit, und = fit(total_rows)
    out = dict(set=setname, omega=omega, tree_term=term, n_configs=nconf, tree_fit=str(tree_fit),
               total_fit={k: str(v) for k, v in (tot_fit or {}).items()}, total_fit_undetermined=und)
    # multiple of the tree vertex: fit total = c * tree on all configs (exact check)
    def tree_multiple(rows):
        cs = set()
        for (t, S), (tr, _) in zip(rows, tree_rows):
            keys = set(t) | set(tr)
            for k in keys:
                tv = sp.nsimplify(tr.get(k, 0)); lv = sp.nsimplify(t.get(k, 0))
                if tv == 0:
                    if sp.simplify(lv) != 0:
                        return None
                    continue
                cs.add(sp.nsimplify(sp.simplify(lv / tv)))
        if len(cs) == 1:
            return cs.pop()
        if len(cs) == 0:
            return sp.Integer(0)
        return None
    out['total_multiple_of_tree'] = str(tree_multiple(total_rows))
    out['predicted_from_A'] = PRED_A[term]
    graphs = []
    for gk, rows in sorted(per_graph_rows.items(), key=lambda kv: str(kv[0])):
        f_, u_ = fit(rows)
        mult = tree_multiple(rows)
        verts = [x[0] for x in gk]; lines = [x[2] for x in gk]
        n = meta[gk]['n']
        sym = '1/2' if (n == 2 and lines[0] == lines[1] and 'c-A*' not in lines) else '1'
        graphs.append(dict(id=' | '.join(f'{v}[{",".join(s + ":" + c for s, c in e)}] -{ln}-' for v, e, ln in gk),
                           vertices=verts, internal_lines=lines, topology={2: 'bubble', 3: 'triangle', 4: 'box'}[n],
                           symmetry_factor=sym, log_eps_multiple_of_tree=(str(mult) if mult is not None else 'not proportional'),
                           log_eps_structure={k: str(v) for k, v in (f_ or {}).items() if v != 0},
                           eps_inv='0', eps_inv_half='0'))
    out['graphs'] = graphs
    out['seconds'] = round(time.time() - t0, 1)
    return out

if __name__ == '__main__':
    for nm in sys.argv[1:]:
        res = run_set(nm)
        with open(f'r2_T1_{nm.replace(",", "_").replace("^", "")}.json', 'w') as fh:
            json.dump(res, fh, indent=1)
        print(json.dumps({k: v for k, v in res.items() if k != 'graphs'}, indent=1))
        for g in res['graphs']:
            print('  ', g['id'], '| sym', g['symmetry_factor'], '| x tree:', g['log_eps_multiple_of_tree'], g['log_eps_structure'])
