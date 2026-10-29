"""Step 2, engine A, Part C: graphs of type I-a (special edge K' = K_eps - K_0 on the cycle) in the four cubic sectors
(c, A|_V, A|_V), (c, A|_Wperp, A|_Wperp), (c, lambda, lambda), (c, c, A^*|_V).

Per sample (random integer polarisations, random integer momenta with sum 0): every oriented sequence of
s2_graphs.OrientedCycles is evaluated exactly (eps^{-1} and eps^0 coefficients as exact combinations of the integrals
J(n;b;N) of s2_integrals, hence elements of Q(i) + Q(i) log 2 + Q(i) log 3).  The log eps coefficient of every I-a
graph is 0 identically (expansion theorem in s2_notes.txt).  Classes = oriented cycles up to rotation/reflection, merged
over permutations of identical external legs.  Each class value (and the total) is fitted EXACTLY on a complete basis
of the Spin(4) x Spin(6) x colour-f invariant multilinear forms of the legs (parity-odd ones included), with extra
samples as a zero-residual check; then projected on G1..G5 (plane-wave formulas below, cross-checked in s2_ga.py).

Plane-wave values (legs in the listed order; vertex-function convention = coefficient of the odd parameters in leg
order; T_{mu nu} = ([gamma, a1_mu], a2_nu), s = sum_{w} ([gamma, a1_w], a2_w), B(p) = ([gamma,u1], rho(p) u2)):
  G1 = (p2^2 - p1^2) tr T - p2.T.p2 + p1.T.p1       G2 = p2.T.p2 - p1.T.p1        G3 = (p2^2 - p1^2) s
  G4 = -i p2 . ([gamma0, gamma1], alpha)            G5 = -i B(p0)
Units: values are hbar/(16 pi^2) times the output, divided by tr_g(ad_X ad_Y)/(X,Y) (C_ADJ)."""
import json, sys, time, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import s2_core as C
from r2_core import Q2
import s2_graphs as SG
from cliff import s, sb

SB = np.array(sb)
with np.errstate(all='ignore'):
    CHI = np.rint((s[0] @ sb[1] @ s[2] @ sb[3]).real)
EPS4 = np.zeros((4, 4, 4, 4))
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])


def fcol(x, y, z):
    """([x,y],z) for colour vectors with arbitrary leading axes: returns einsum over the colour axis"""
    return np.einsum('bce,...b,...c,...e->...', C.F, x, y, z)


def Tmat(g, a1, a2):
    return np.einsum('bce,b,mc,ne->mn', C.F, g, a1[:4], a2[:4])


def st_AVAV(L, e):
    g, a1, a2 = L[0][2], L[1][2], L[2][2]; p1, p2 = L[1][3], L[2][3]
    T = Tmat(g, a1, a2); tr = np.trace(T)
    if e == -1:
        return [tr], ['tr T']
    return [tr * (p1 @ p1), tr * (p2 @ p2), tr * (p1 @ p2), p1 @ T @ p1, p1 @ T @ p2, p2 @ T @ p1, p2 @ T @ p2,
            np.einsum('mnrs,mn,r,s->', EPS4, T, p1, p2)], \
           ['trT p1^2', 'trT p2^2', 'trT p1.p2', 'p1.T.p1', 'p1.T.p2', 'p2.T.p1', 'p2.T.p2', 'eps(T,p1,p2)']


def st_APAP(L, e):
    g, a1, a2 = L[0][2], L[1][2], L[2][2]; p1, p2 = L[1][3], L[2][3]
    sv = np.einsum('bce,b,mc,me->', C.F, g, a1[4:], a2[4:])
    if e == -1:
        return [sv], ['s']
    return [sv * (p1 @ p1), sv * (p2 @ p2), sv * (p1 @ p2)], ['s p1^2', 's p2^2', 's p1.p2']


def Bform(g, u1, u2, p, chi=False):
    M = np.einsum('m,mab->ab', np.asarray(p, dtype=complex), SB[:4])
    if chi:
        with np.errstate(all="ignore"):
            M = M @ CHI
        assert np.all(np.isfinite(M))
    return np.einsum('bce,b,Ac,AB,Be->', C.F, g, u1, M, u2)


def st_cll(L, e):
    g, u1, u2 = L[0][2], L[1][2], L[2][2]; p1, p2 = L[1][3], L[2][3]
    if e == -1:
        return [], []
    return [Bform(g, u1, u2, p1), Bform(g, u1, u2, p2), Bform(g, u1, u2, p1, True), Bform(g, u1, u2, p2, True)], \
           ['B(p1)', 'B(p2)', 'Bchi(p1)', 'Bchi(p2)']


def st_ccAs(L, e):
    g0, g1, al = L[0][2], L[1][2], L[2][2]; p0, p1 = L[0][3], L[1][3]
    if e == -1:
        return [], []
    v = np.einsum('bce,b,c,me->m', C.F, g0, g1, al[:4])
    return [p0 @ v, p1 @ v], ['p0.([g0,g1],alpha)', 'p1.([g0,g1],alpha)']


def G_AVAV(L):
    g, a1, a2 = L[0][2], L[1][2], L[2][2]; p1, p2 = L[1][3], L[2][3]
    T = Tmat(g, a1, a2)
    return {'G1': (p2 @ p2 - p1 @ p1) * np.trace(T) - p2 @ T @ p2 + p1 @ T @ p1, 'G2': p2 @ T @ p2 - p1 @ T @ p1}


def G_APAP(L):
    g, a1, a2 = L[0][2], L[1][2], L[2][2]; p1, p2 = L[1][3], L[2][3]
    sv = np.einsum('bce,b,mc,me->', C.F, g, a1[4:], a2[4:])
    return {'G3': (p2 @ p2 - p1 @ p1) * sv}


def G_cll(L):
    g, u1, u2 = L[0][2], L[1][2], L[2][2]; p0 = L[0][3]
    return {'G5': -1j * Bform(g, u1, u2, p0)}


def G_ccAs(L):
    g0, g1, al = L[0][2], L[1][2], L[2][2]; p2 = L[2][3]
    return {'G4': -1j * (p2 @ np.einsum('bce,b,c,me->m', C.F, g0, g1, al[:4]))}


SECTORS = {
    'c,AV,AV': dict(legs=[('c', None), ('A', 'V'), ('A', 'V')], deg=4, st=st_AVAV, G=G_AVAV),
    'c,AP,AP': dict(legs=[('c', None), ('A', 'P'), ('A', 'P')], deg=4, st=st_APAP, G=G_APAP),
    'c,l,l': dict(legs=[('c', None), ('l', None), ('l', None)], deg=3, st=st_cll, G=G_cll),
    'c,c,AsV': dict(legs=[('c', None), ('c', None), ('As', 'V')], deg=3, st=st_ccAs, G=G_ccAs),
}
LEGNAME = {('c', None): 'c', ('A', 'V'): 'A|_V', ('A', 'P'): 'A|_Wperp', ('l', None): 'lambda', ('As', 'V'): 'A^*|_V'}
SLOT = {'A': 'A', 'As': 'A^*', 'l': 'lambda', 'ls': 'lambda^*', 'c': 'c', 'cs': 'c^*'}
LNAME = {'AV': 'A|_V--A|_V', 'AP': 'A|_Wperp--A|_Wperp', 'gh': 'c--A^*|_V', 'l': 'lambda--lambda',
         'sAV': 'K_eps-K_0 (A|_V--A^*|_V)', 'sAP': 'K_eps-K_0 (A|_Wperp--A^*|_Wperp)', 'sl': 'K_eps-K_0 (lambda^*--lambda)',
         'sc': 'K_eps-K_0 (c^*--c)'}


def toS(z):
    return sp.Rational(z.re.numerator, z.re.denominator) + sp.I * sp.Rational(z.im.numerator, z.im.denominator)


RANKS = []


def solve_exact(rows, rhs):
    M = sp.Matrix(rows); b = sp.Matrix(rhs)
    if M.shape[1] == 0:
        return [], all(x == 0 for x in b)
    RANKS.append((M.rank(), M.shape[1]))
    if RANKS[-1][0] < M.shape[1]:
        raise ValueError('structure matrix not of full column rank: %s' % (RANKS[-1],))
    if all(x == 0 for x in b):
        return [sp.Integer(0)] * M.shape[1], True
    sol, params = M.gauss_jordan_solve(b)
    if params.shape[0] > 0:
        sol = sol.subs({p_: 0 for p_ in params})
    resid = (M * sol - b).applyfunc(sp.expand)
    return [sp.nsimplify(sp.expand(x)) for x in sol], all(r == 0 for r in resid)


def vertex_label(ext_labels, slot_in, slot_out):
    return '{' + ','.join(sorted(list(ext_labels) + ['(' + SLOT[slot_in] + ')', '(' + SLOT[slot_out] + ')'])) + '}'


def class_key(seq, Ls, legs):
    lab = [tuple(sorted(LEGNAME[legs[j]] for j in b)) for b in seq]
    return SG.canon_ocycle([(Ls[i], lab[i]) for i in range(len(seq))])


def describe(key):
    n = len(key)
    verts = []; lines = []; special = None
    for i in range(n):
        (Li, oi), lab = key[i]
        (Ln, on), _ = key[(i + 1) % n]
        verts.append(vertex_label(lab, oi[1], on[0]))
    for i in range(n):
        (Li, oi), _ = key[i]
        txt = '%s, slots %s at %s and %s at %s' % (LNAME[Li], SLOT[oi[0]], verts[(i - 1) % n], SLOT[oi[1]], verts[i])
        if Li.startswith('s'):
            special = txt
        else:
            lines.append(txt)
    return {'topology': {2: 'bubble (one quartic vertex)', 3: 'triangle'}[n], 'vertices': verts,
            'lines': lines, 'special_edge': special}


def run_sector(name, rng, nextra=3, verbose=True):
    cfg = SECTORS[name]
    legs = cfg['legs']
    t0 = time.time()
    samples = []
    nst = None
    while True:
        while True:
            ps = [rng.integers(-2, 3, size=4) for _ in range(2)]
            ps = [-(ps[0] + ps[1])] + ps
            if all(p.any() for p in ps):
                break
        L = [(t, w, C.rpol(rng, t, w), p) for (t, w), p in zip(legs, ps)]
        st = {e: cfg['st'](L, e) for e in (-1, 0)}
        if nst is None:
            nst = len(st[0][0])
        Gv = cfg['G'](L)
        G = SG.OrientedCycles(L, cfg['deg'])
        cls = {}
        for (seq, Ls, sgn) in G.enumerate():
            v = G.graph_value(seq, Ls, sgn)
            key = class_key(seq, Ls, legs)
            ent = cls.setdefault(key, {'nseq': 0, -1: [Q2(0)] * 3, 0: [Q2(0)] * 3})
            ent['nseq'] += 1
            for e in (-1, 0):
                b = SG.to_basis(v[e])
                ent[e] = [x + y * Fraction(1, C.C_ADJ) for x, y in zip(ent[e], b)]
        samples.append(dict(L=L, st=st, G=Gv, cls=cls))
        if verbose:
            print(name, C.ALG, 'sample', len(samples), 'classes', len(cls), 'time %.1f' % (time.time() - t0), flush=True)
        if len(samples) >= nst + nextra:
            break
    keys = sorted(set().union(*[set(smp['cls']) for smp in samples]), key=repr)
    names = {e: samples[0]['st'][e][1] for e in (-1, 0)}
    rows = {e: [[toS(Q2.of(x)) for x in smp['st'][e][0]] for smp in samples] for e in (-1, 0)}
    LOGB = ['1', 'log2', 'log3']
    out = {'sector': name, 'algebra': C.ALG, 'structures': {str(e): names[e] for e in (-1, 0)}, 'classes': []}
    tot = {e: [[sp.Integer(0)] * len(names[e]) for _ in LOGB] for e in (-1, 0)}
    allok = True
    for key in keys:
        ent = {'class': repr(key), **describe(key)}
        nseq = {smp['cls'][key]['nseq'] if key in smp['cls'] else 0 for smp in samples}
        ent['n_oriented_sequences'] = sorted(nseq)
        ent['multiplicity_and_weight'] = '%d labelled graphs (sum over legs permutations of identical legs), each with ' \
                                         'weight 1/|Aut| = 1 [(1/2) str(K\'SH(PSH)^{n-1}) x 2 directions]' % (max(nseq) // 2)
        coeffs = {}
        for e in (-1, 0):
            comp = {}
            for ib, lb in enumerate(LOGB):
                rhs = [toS(smp['cls'][key][e][ib]) if key in smp['cls'] else 0 for smp in samples]
                if not names[e]:
                    assert all(x == 0 for x in rhs), (name, key, e, rhs)
                    sol, ok = [], True
                else:
                    sol, ok = solve_exact(rows[e], rhs)
                allok = allok and ok
                tot[e][ib] = [a + b for a, b in zip(tot[e][ib], sol)]
                comp[lb] = {nm: str(x) for nm, x in zip(names[e], sol) if x != 0}
                comp[lb + '_zero_residual'] = ok
            coeffs['eps^-1' if e == -1 else 'eps^0'] = comp
        coeffs['log'] = '0 (identically; no log eps term in any I-a graph)'
        ent['coeffs_in_structure_basis'] = coeffs
        out['classes'].append(ent)
    # G-basis in the structure basis
    Gnames = sorted(samples[0]['G'])
    Gcols = []
    for gn in Gnames:
        sol, ok = solve_exact(rows[0], [toS(Q2.of(smp['G'][gn])) for smp in samples])
        assert ok
        Gcols.append(sol)
    out['G_in_structure_basis'] = {gn: {nm: str(x) for nm, x in zip(names[0], col) if x != 0} for gn, col in zip(Gnames, Gcols)}
    Gm = sp.Matrix(Gcols).T

    def project(vec):
        sol, ok = solve_exact(Gm.tolist(), vec)
        return sol, ok
    for cl in out['classes']:
        pr = {}
        for ib, lb in enumerate(LOGB):
            vec = [sp.Rational(0)] * len(names[0])
            dd = cl['coeffs_in_structure_basis']['eps^0'][lb]
            vec = [sp.sympify(dd.get(nm, '0')) for nm in names[0]]
            sol, ok = project(vec)
            pr[lb] = {gn: str(x) for gn, x in zip(Gnames, sol)}
            pr[lb + '_in_span_G'] = ok
        e1 = cl['coeffs_in_structure_basis']['eps^-1']
        cl['coeffs'] = {'eps^-1': {gn: '0' for gn in Gnames} if all(not e1[lb] for lb in LOGB) else 'NONZERO: ' + json.dumps(e1),
                        'log': {gn: '0' for gn in Gnames},
                        'eps^0': {gn: str(sp.nsimplify(sum(sp.sympify(pr[lb][gn]) * w for lb, w in
                                                         zip(LOGB, [1, sp.log(2), sp.log(3)])))) for gn in Gnames},
                        'eps^0_in_span_G': all(pr[lb + '_in_span_G'] for lb in LOGB)}
    total = {}
    for e in (-1, 0):
        total[str(e)] = {lb: {nm: str(x) for nm, x in zip(names[e], tot[e][ib]) if x != 0} for ib, lb in enumerate(LOGB)}
    out['total_in_structure_basis'] = total
    alpha = {}
    inspan = True
    for ib, lb in enumerate(LOGB):
        sol, ok = project(tot[0][ib])
        inspan = inspan and ok
        for gn, x in zip(Gnames, sol):
            alpha.setdefault(gn, 0)
            alpha[gn] += x * [1, sp.log(2), sp.log(3)][ib]
    out['total_eps^0_G'] = {gn: str(sp.nsimplify(v)) for gn, v in alpha.items()}
    out['total_eps^0_in_span_G_zero_residual'] = inspan
    out['total_eps^-1'] = 'zero' if all(all(x == 0 for x in tot[-1][ib]) for ib in range(3)) else str(tot[-1])
    out['all_fits_zero_residual'] = allok
    out['samples'] = len(samples)
    out['structure_matrix_full_column_rank'] = all(r == c for r, c in RANKS)
    RANKS.clear()
    out['seconds'] = round(time.time() - t0, 1)
    return out


if __name__ == '__main__':
    which = sys.argv[1].split(';') if len(sys.argv) > 1 and sys.argv[1] != 'all' else list(SECTORS)
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 2028
    rng = np.random.default_rng(seed)
    res = {}
    for nm in which:
        r = run_sector(nm, rng)
        res[nm] = r
        print(json.dumps({k: v for k, v in r.items() if k not in ('classes',)}, indent=1), flush=True)
        for cl in r['classes']:
            print('  ', cl['topology'], cl['special_edge'], '|', cl['coeffs']['eps^0'], flush=True)
        json.dump(res, open('logs/s2_Ia_%s_%s_seed%d.json' % (C.ALG, '-'.join(w.replace(',', '') for w in which), seed),
                            'w'), indent=1)
