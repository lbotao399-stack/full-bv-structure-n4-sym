"""T1 (round 2, engine A): cubic and quartic one-loop log-eps coefficients, graph by graph, exact.
Each external set: random integer polarizations/momenta (exact), all 1PI one-loop graphs enumerated by r2_graphs,
log-eps coefficient per graph class fitted EXACTLY (sympy rational Gauss-Jordan over Q(i)) on a complete basis of
Spin(4) x Spin(6) x colour invariant local structures (su(2), orthonormal, eps_abc), with extra samples as an exact
consistency check (zero residual required).  Units: u = hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) (here hbar c_adj/16pi^2,
c_adj = -2)."""
import json, sys, time, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import r2_core as C
from r2_core import Q2
from r2_graphs import GraphSum
from cliff import s, sb

F = C.F.real
with np.errstate(all="ignore"):
    G5 = (s[0] @ sb[1] @ s[2] @ sb[3]).real.round()            # 4d chirality on Sigma (integer entries)
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 2027)

def rpol(t, which=None):
    x = rng.integers(-1, 2, size=C.SHAPE[t]).astype(complex)
    if which == 'V': x[4:] = 0
    if which == 'P': x[:4] = 0
    return x

def toS(z):
    return sp.Rational(Fraction(z.re).numerator, Fraction(z.re).denominator) + sp.I * sp.Rational(Fraction(z.im).numerator, Fraction(z.im).denominator)

def solve_exact(rows, rhs):
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    if all(x == 0 for x in b):
        return [sp.Integer(0)] * M.shape[1], True, M.rank()
    sol, params = M.gauss_jordan_solve(b)
    if params.shape[0] > 0:
        sol = sol.subs({p_: 0 for p_ in params})
    resid = (M * sol - b).applyfunc(sp.expand)
    return [sp.nsimplify(sp.expand(x)) for x in sol], all(r == 0 for r in resid), M.rank()

# ------------------------------------------------------------------ invariant structures (exact, small integers)
def eps4():
    e = np.zeros((4, 4, 4, 4))
    for p in itertools.permutations(range(4)):
        e[p] = np.linalg.det(np.eye(4)[list(p)])
    return e
EPS4 = eps4()

def colour3(x, y, z):
    return np.einsum('abc,...a,...b,...c->...', F, x, y, z)

def st_AV3(L):
    (a1, p1), (a2, p2), (a3, p3) = [(l[2][:4], l[3]) for l in L]
    G = np.einsum('abc,ma,nb,rc->mnr', F, a1, a2, a3)
    out = []
    for p in (p1, p2):
        out += [np.einsum('mmr,r->', G, p), np.einsum('mnm,n->', G, p), np.einsum('mnn,m->', G, p)]
    for p in (p1, p2):
        out.append(np.einsum('mnr,mnrs,s->', G, EPS4, p))
    names = ['(a1.a2)(p1.a3)', '(a1.a3)(p1.a2)', '(a2.a3)(p1.a1)', '(a1.a2)(p2.a3)', '(a1.a3)(p2.a2)', '(a2.a3)(p2.a1)',
             'eps(a1,a2,a3,p1)', 'eps(a1,a2,a3,p2)']
    return out, names

def st_AVAPAP(L):
    a1, p1 = L[0][2][:4], L[0][3]; a2 = L[1][2][4:]; a3 = L[2][2][4:]; p2 = L[1][3]
    G = np.einsum('abc,ma,ib,ic->m', F, a1, a2, a3)
    return [G @ p1, G @ p2], ['(a2.a3)(p1.a1)', '(a2.a3)(p2.a1)']

def st_llA(L, which):
    u1, u2, a = L[0][2], L[1][2], L[2][2]
    Ms = range(4) if which == 'V' else range(4, 10)
    v0 = sum(np.einsum('abc,ia,ij,jb,c->', F, u1, sb[M], u2, a[M]) for M in Ms)
    v5 = sum(np.einsum('abc,ia,ij,jb,c->', F, u1, sb[M] @ G5, u2, a[M]) for M in Ms)
    return [v0, v5], ['(u1,rho(a)[.,u2])', '(u1,rho(a) chirality [.,u2])']

def st_cAAs(L, which):
    g, a, al = L[0][2], L[1][2], L[2][2]
    Ms = range(4) if which == 'V' else range(4, 10)
    return [sum(np.einsum('abc,a,b,c->', F, g, a[M], al[M]) for M in Ms)], ['([gamma,a],alpha)']

def st_cllS(L):
    g, u, v = L[0][2], L[1][2], L[2][2]
    return [np.einsum('abc,a,ib,ic->', F, g, u, v), np.einsum('abc,a,ib,ij,jc->', F, g, u, G5.T, v)], \
           ['([gamma,u],v)', '([gamma,chirality u],v)']

def st_cccS(L):
    return [colour3(L[0][2], L[1][2], L[2][2])], ['([g1,g2],w)']

def colour_pairs4(c1, c2, c3, c4):
    """X=d12 d34, Y=d13 d24, Z=d14 d23 applied to colour vectors"""
    return [(c1 @ c2) * (c3 @ c4), (c1 @ c3) * (c2 @ c4), (c1 @ c4) * (c2 @ c3)]

def st_quartic(L, kinds):
    a = [l[2] for l in L]
    sl = [slice(0, 4) if k == 'V' else slice(4, 10) for k in kinds]
    b = [a[i][sl[i]] for i in range(4)]     # (lorentz, colour)
    out = []; names = []
    LP = {'x': ((0, 1), (2, 3)), 'y': ((0, 2), (1, 3)), 'z': ((0, 3), (1, 2))}
    CP = {'X': ((0, 1), (2, 3)), 'Y': ((0, 2), (1, 3)), 'Z': ((0, 3), (1, 2))}
    for cn, (cp1, cp2) in CP.items():
        for ln, (lp1, lp2) in LP.items():
            ok = all(kinds[i] == kinds[j] for (i, j) in (lp1, lp2))
            if not ok:
                continue
            # contract: lorentz pairs lp, colour pairs cp
            idxL = ['m', 'n', 'r', 's']; idxC = ['a', 'b', 'c', 'd']
            Li = list(idxL); Ci = list(idxC)
            Li[lp1[1]] = Li[lp1[0]]; Li[lp2[1]] = Li[lp2[0]]
            Ci[cp1[1]] = Ci[cp1[0]]; Ci[cp2[1]] = Ci[cp2[0]]
            expr = ','.join(Li[i] + Ci[i] for i in range(4)) + '->'
            out.append(np.einsum(expr, *b)); names.append(cn + ln)
    if all(k == 'V' for k in kinds):
        for cn, (cp1, cp2) in CP.items():
            Ci = ['a', 'b', 'c', 'd']; Ci[cp1[1]] = Ci[cp1[0]]; Ci[cp2[1]] = Ci[cp2[0]]
            expr = 'mnrs,' + ','.join(['m', 'n', 'r', 's'][i] + Ci[i] for i in range(4)) + '->'
            out.append(np.einsum(expr, EPS4, *b)); names.append(cn + '_eps4')
    return out, names

def tree(L, only, scale):
    legs = [(l[0], l[2][None], l[3]) for l in L]
    v = C.V(legs, only=[only])
    return Q2.of(complex(np.asarray(v).reshape(-1)[0])) * Fraction(1, scale)

SETS = {
    'AV,AV,AV':   dict(legs=[('A', 'V'), ('A', 'V'), ('A', 'V')], omega=1, st=st_AV3, tree=('AAA', 2)),
    'AV,AP,AP':   dict(legs=[('A', 'V'), ('A', 'P'), ('A', 'P')], omega=1, st=st_AVAPAP, tree=('AAA', 2)),
    'l,l,AV':     dict(legs=[('l', None), ('l', None), ('A', 'V')], omega=0, st=lambda L: st_llA(L, 'V'), tree=('llA', 2)),
    'l,l,AP':     dict(legs=[('l', None), ('l', None), ('A', 'P')], omega=0, st=lambda L: st_llA(L, 'P'), tree=('llA', 2)),
    'c,AV,AsV':   dict(legs=[('c', None), ('A', 'V'), ('As', 'V')], omega=0, st=lambda L: st_cAAs(L, 'V'), tree=('AcAs', 2)),
    'c,AP,AsP':   dict(legs=[('c', None), ('A', 'P'), ('As', 'P')], omega=0, st=lambda L: st_cAAs(L, 'P'), tree=('AcAs', 2)),
    'c,l,ls':     dict(legs=[('c', None), ('l', None), ('ls', None)], omega=0, st=st_cllS, tree=('cllS', 2)),
    'c,c,cs':     dict(legs=[('c', None), ('c', None), ('cs', None)], omega=0, st=st_cccS, tree=('cccS', 2)),
    'AV,AV,AV,AV': dict(legs=[('A', 'V')] * 4, omega=0, st=lambda L: st_quartic(L, 'VVVV'), tree=('AAAA', 8)),
    'AV,AV,AP,AP': dict(legs=[('A', 'V'), ('A', 'V'), ('A', 'P'), ('A', 'P')], omega=0, st=lambda L: st_quartic(L, 'VVPP'), tree=('AAAA', 8)),
    'AP,AP,AP,AP': dict(legs=[('A', 'P')] * 4, omega=0, st=lambda L: st_quartic(L, 'PPPP'), tree=('AAAA', 8)),
}
# prediction of assumption (A): sum over legs of the weights {Psi,-}: A 2, lambda 4, c 1, A^* -2, lambda^* -4, c^* -1
WEIGHT = {'A': 2, 'l': 4, 'c': 1, 'As': -2, 'ls': -4, 'cs': -1}

def leg_label(t, w):
    return {'A': 'A|_V' if w == 'V' else 'A|_Wperp', 'l': 'lambda', 'c': 'c', 'As': 'A^*|_V' if w == 'V' else 'A^*|_Wperp',
            'ls': 'lambda^*', 'cs': 'c^*'}[t]

def graph_type(key, legs):
    """merge classes related by permutations of identical external legs"""
    parts = []
    for (L, blk) in key:
        parts.append((L, tuple(sorted(leg_label(*legs[j]) for j in blk))))
    n = len(parts); cands = []
    for r in range(n):
        rot = parts[r:] + parts[:r]
        cands.append(tuple(rot))
        Ls = [p[0] for p in rot]; Ss = [p[1] for p in rot]
        cands.append(tuple([(Ls[0], Ss[n - 1])] + [(Ls[n - i], Ss[n - 1 - i]) for i in range(1, n)]))
    return min(cands, key=repr)

LINE_NAME = {'AV': 'A|_V-A|_V', 'AP': 'A|_Wperp-A|_Wperp', 'gh': 'c-A^*|_V', 'l': 'lambda-lambda'}

def vertex_name(ext_types, lin, lout):
    ends = {'AV': ['A'], 'AP': ['A'], 'l': ['l'], 'gh': ['c', 'As']}
    for vt, can in C.CANON.items():
        for e1 in ends[lin]:
            for e2 in ends[lout]:
                if sorted(can) == sorted(list(ext_types) + [e1, e2]):
                    return vt
    return '?'

def run_set(name, nextra=3, verbose=True):
    cfg = SETS[name]
    legs = cfg['legs']; omega = cfg['omega']
    samples = []
    t0 = time.time()
    nst = None
    it = 0
    while True:
        L = []
        if omega == 1:
            ps = [rng.integers(-2, 3, size=4) * 6 for _ in range(len(legs) - 1)]
            ps.append(-sum(ps))
        else:
            ps = [rng.integers(-2, 3, size=4) * 6 for _ in range(len(legs) - 1)]
            ps.append(-sum(ps))
        for (t, w), p in zip(legs, ps):
            L.append((t, w, rpol(t, w), p))
        st, names = cfg['st'](L)
        if nst is None:
            nst = len(st)
        G = GraphSum(L)
        res = G.coefficients(omega)
        tr = tree(L, *cfg['tree'])
        samples.append(dict(st=[Q2.of(x) for x in st], res=res, tree=tr))
        it += 1
        if verbose:
            print(name, 'sample', it, 'classes', len(res), 'time %.1f' % (time.time() - t0), flush=True)
        if it >= nst + nextra:
            break
    # group classes into graph types
    gtypes = {}
    for smp in samples:
        agg = {}
        for key, v in smp['res'].items():
            gt = graph_type(key, legs)
            if gt not in agg:
                agg[gt] = {'log': Q2(0), 'epsinv': Q2(0), 'classes': set()}
            agg[gt]['log'] = agg[gt]['log'] + v['log']; agg[gt]['epsinv'] = agg[gt]['epsinv'] + v['epsinv']
            agg[gt]['classes'].add(key)
        smp['agg'] = agg
        for gt in agg:
            gtypes.setdefault(gt, set()).update(agg[gt]['classes'])
    rows = [[toS(x) for x in smp['st']] for smp in samples]
    out = {'graphs': [], 'structures': names}
    tot = [sp.Integer(0)] * nst
    for gt in sorted(gtypes, key=repr):
        rhs = [toS(smp['agg'][gt]['log'] * Fraction(1, C.C_ADJ)) if gt in smp['agg'] else 0 for smp in samples]
        sol, ok, rk = solve_exact(rows, rhs)
        rhs_e = [toS(smp['agg'][gt]['epsinv']) if gt in smp['agg'] else 0 for smp in samples]
        assert all(x == 0 for x in rhs_e)
        tot = [a + b for a, b in zip(tot, sol)]
        n = len(gt)
        verts = []
        for i in range(n):
            lin = gt[i][0]; lout = gt[(i + 1) % n][0]
            ext_t = [{'A|_V': 'A', 'A|_Wperp': 'A', 'lambda': 'l', 'c': 'c', 'A^*|_V': 'As', 'A^*|_Wperp': 'As',
                      'lambda^*': 'ls', 'c^*': 'cs'}[x] for x in gt[i][1]]
            verts.append(vertex_name(ext_t, lin, lout) + '[' + ','.join(gt[i][1]) + ']')
        ncls = len(gtypes[gt])
        out['graphs'].append(dict(
            id=name + ' | ' + ' -> '.join('%s ~%s~' % (v, LINE_NAME[gt[(i + 1) % n][0]]) for i, v in enumerate(verts)),
            topology={1: 'tadpole', 2: 'bubble', 3: 'triangle', 4: 'box'}[n],
            vertices=verts, internal_lines=[LINE_NAME[gt[(i + 1) % n][0]] for i in range(n)],
            n_graphs_in_class=ncls,
            symmetry_factor='1/|Aut| = 1/2' if (n == 2 and gt[0][0] == gt[1][0] and gt[0][0] != 'gh') else '1/|Aut| = 1',
            log_eps={nm: str(x) for nm, x in zip(names, sol) if x != 0},
            eps_inv='0', exact_fit=ok, rank=rk))
    trs = [toS(smp['tree']) for smp in samples]
    tsol, tok, trk = solve_exact(rows, trs)
    out['tree_in_basis'] = {nm: str(x) for nm, x in zip(names, tsol) if x != 0}
    out['total_log_eps'] = {nm: str(x) for nm, x in zip(names, tot) if x != 0}
    # proportionality to the tree vertex
    ratio = None; prop = True
    for x, t in zip(tot, tsol):
        if t == 0:
            if x != 0: prop = False
            continue
        r = sp.nsimplify(x / t)
        if ratio is None: ratio = r
        elif sp.simplify(r - ratio) != 0: prop = False
    out['total_is_multiple_of_tree'] = prop
    out['computed_multiple'] = str(ratio) if prop else None
    pred = sum(WEIGHT[t] for (t, w) in legs)
    out['predicted_from_A'] = pred
    out['verdict'] = ('agrees with (A)' if prop and ratio == pred else 'DISAGREES with (A)')
    out['samples'] = len(samples)
    out['seconds'] = round(time.time() - t0, 1)
    return out

if __name__ == '__main__':
    which = sys.argv[1].split(';') if len(sys.argv) > 1 and sys.argv[1] != 'all' else list(SETS)
    results = {}
    for nm in which:
        r = run_set(nm)
        results[nm] = r
        print(json.dumps({k: v for k, v in r.items() if k != 'graphs'}), flush=True)
        for g in r['graphs']:
            print('   ', g['id'], g['log_eps'], g['n_graphs_in_class'], flush=True)
        json.dump(results, open('logs/r2_T1_%s.json' % (sys.argv[3] if len(sys.argv) > 3 else 'run'), 'w'), indent=1)
