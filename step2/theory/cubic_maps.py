#!/usr/bin/env python3
"""cubic_maps.py -- theory agent, Step 2 Part C: exact check of the cubic maps on random plane waves.

Computes, in the basis G1..G5 (DESIGN_STEP2.md), the images Q H_j, {I_3, X} for the quadratic functionals X of F_2,
Q G_i, {I_3, Ob_2}, the tree of Ob_2 with I_3 through P(0,infinity), and runs the consistency checks of the
realisation (Clifford relations, cubic vertex of ch01, {S_2,-} = Q, Q I_3 = 0, quartic master equation,
field-count lemma, s-closure of the completions, Q P(0,infinity) = K_0 on trees).
Writes cubic_maps_numeric.json next to this file.  Run: python3 cubic_maps.py
"""
import json
import random
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory')
import cubic_core as C  # noqa: E402
from cubic_core import GQ, ZERO, ONE, Mode  # noqa: E402

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory/cubic_maps_numeric.json'
LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)


SECTORS = {
    'c_AV_AV': [('c', None), ('A', 'V'), ('A', 'V')],
    'c_AW_AW': [('c', None), ('A', 'W'), ('A', 'W')],
    'c_lam_lam': [('c', None), ('lam', None), ('lam', None)],
    'c_c_AsV': [('c', None), ('c', None), ('As', 'V')],
    'deg2_c_c_AV': [('c', None), ('c', None), ('A', 'V')],
    'cAAA4': [('c', None), ('A', 'V'), ('A', 'V'), ('A', 'V')],
    'cAWW4': [('c', None), ('A', 'V'), ('A', 'W'), ('A', 'W')],
    'AAA': [('A', 'all'), ('A', 'all'), ('A', 'all')],
    'c_lam_lams': [('c', None), ('lam', None), ('lams', None)],
    'c_c_cs': [('c', None), ('c', None), ('cs', None)],
    'A_c_As': [('A', 'all'), ('c', None), ('As', 'all')],
    'lam_A_lam': [('lam', None), ('A', 'all'), ('lam', None)],
}


def shash(s):
    return sum((i + 1) * ord(ch) for i, ch in enumerate(s))


def rint(rng, a=-3, b=3, nonzero=False):
    while True:
        x = rng.randint(a, b)
        if x or not nonzero:
            return x


def sample(sector, seed):
    rng = random.Random(seed)
    legs = SECTORS[sector]
    n = len(legs)
    while True:
        ks = [tuple(Fr(rint(rng, -4, 4)) for _ in range(4)) for _ in range(n - 1)]
        ks.append(tuple(-sum((k[i] for k in ks), Fr(0)) for i in range(4)))
        sums = set()
        ok = True
        import itertools
        for r in range(1, n):
            for sub in itertools.combinations(range(n), r):
                s = tuple(sum((ks[j][i] for j in sub), Fr(0)) for i in range(4))
                if C.iszerok(s):
                    ok = False
                sums.add(s)
        if ok:
            break
    modes, gens = [], []
    for j, (s, pol) in enumerate(legs):
        g = C.new_gen(C.SPAR[s], 'ext')
        gens.append(g)
        d = C.SDIM[s]
        if s in ('c', 'cs'):
            u = [ONE]
        elif s in ('A', 'As'):
            idx = C.V4 if pol == 'V' else (C.WP if pol == 'W' else range(10))
            u = [GQ(rint(rng)) if M in idx else ZERO for M in range(10)]
        else:
            u = [GQ(rint(rng)) for _ in range(d)]
        X = [Fr(rint(rng)) for _ in range(C.DIMG)]
        modes.append(Mode(C.gen_elem(g), s, u, X, ks[j]))
    kset = sorted(sums)
    return modes, gens, kset


CACHE = {}


def internal_for(kset):
    key = tuple(kset)
    if key not in CACHE:
        CACHE.clear()
        CACHE[key] = C.make_internal(kset)
    return CACHE[key]


def val_of(elem, gens):
    return C.ext_coeff(elem, gens)


def solve_exact(rows, rhs):
    """least-squares-free exact solve: rows (list of lists GQ) x = rhs; returns x and max residual flag."""
    n = len(rows[0])
    A = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv = []
    r = 0
    for col in range(n):
        p = None
        for i in range(r, len(A)):
            if not A[i][col].iszero():
                p = i
                break
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][col]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and not A[i][col].iszero():
                f = A[i][col]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(col)
        r += 1
    consistent = all(A[i][n].iszero() for i in range(r, len(A)))
    x = [ZERO] * n
    for i, col in enumerate(piv):
        x[col] = A[i][n]
    return x, consistent, len(piv)


def fmt(x):
    if x.im == 0:
        return str(x.re)
    return f'{x.re}+{x.im}*i'


SECTOR_BASIS = {
    'c_AV_AV': ['G1', 'G2'],
    'c_AW_AW': ['G3'],
    'c_lam_lam': ['G5'],
    'c_c_AsV': ['G4'],
}


def fn_K5odd(t):  # ([c, lam], dslash rho(vol4) lam): the parity-odd candidate of the sector (c, lam, lam)
    c, l1, l2 = t
    w = C.rho_vec([C.IU * l2.k[mu] for mu in range(4)], C.vol4(l2.u))
    return c.u[0] * C.spin_pair(l1.u, w) * C.f3(c.X, l1.X, l2.X)


K5odd = [C.term(1, ('c', 'lam', 'lam'), fn_K5odd)]


def fn_w2_AV(t):  # ([c, A_V], A_V): weight 2
    c, a, b = t
    return c.u[0] * C.dot(C.Vpart(a.u), b.u) * C.f3(c.X, a.X, b.X)


def fn_w2_AW(t):
    c, a, b = t
    return c.u[0] * C.dot(C.Wpart(a.u), b.u) * C.f3(c.X, a.X, b.X)


W2 = {'w2_AV': [C.term(1, ('c', 'A', 'A'), fn_w2_AV)], 'w2_AW': [C.term(1, ('c', 'A', 'A'), fn_w2_AW)]}


def decompose(name, func, sectors, nsamp=None, extra_basis=None):
    """func(modes, gens, kset) -> GQ value; decomposed sector by sector in the G basis."""
    res = {}
    for sec in sectors:
        basis = SECTOR_BASIS[sec] + (extra_basis or {}).get(sec, [])
        ns = nsamp or (len(basis) + 3)
        rows, rhs = [], []
        for s in range(ns):
            modes, gens, kset = sample(sec, (1000 * shash(sec)) % 9973 + 17 * s + 3)
            row = []
            for b in basis:
                F = C.G[b] if b in C.G else (K5odd if b == 'K5odd' else W2[b])
                row.append(val_of(C.evaluate(F, modes), gens))
            rows.append(row)
            rhs.append(func(modes, gens, kset))
        x, ok, rank = solve_exact(rows, rhs)
        res[sec] = {b: fmt(v) for b, v in zip(basis, x)}
        res[sec]['_residual_zero'] = ok
        res[sec]['_rank'] = rank
        res[sec]['_samples'] = ns
    log(name, json.dumps(res))
    return res


def main():
    t0 = time.time()
    out = {'checks': {}, 'results': {}}
    # ---------------- C0: realisation
    cl = True
    for M in range(10):
        for N in range(10):
            for v0 in range(16):
                e = [ONE if j == v0 else ZERO for j in range(16)]
                w = [a + b for a, b in zip(C.matvec(C.SIG[M], C.matvec(C.SIGB[N], e)),
                                             C.matvec(C.SIG[N], C.matvec(C.SIGB[M], e)))]
                tgt = [GQ(2) if (M == N and j == v0) else ZERO for j in range(16)]
                if any(not (a - b).iszero() for a, b in zip(w, tgt)):
                    cl = False
    symm = True
    for M in range(10):
        dense = [[ZERO] * 16 for _ in range(16)]
        for a in range(16):
            for (b, c) in C.SIG[M][a]:
                dense[a][b] = c
        if any(not (dense[a][b] - dense[b][a]).iszero() for a in range(16) for b in range(16)):
            symm = False
    out['checks']['clifford_sig_sigbar'] = cl
    out['checks']['rho_symmetric_(rho(v)Q1,Q2)=(Q1,rho(v)Q2)'] = symm
    log('C0 clifford', cl, 'symmetric', symm)

    # ---------------- C1: cubic gauge vertex vs eq:c1-AAA
    okv = True
    for s in range(3):
        modes, gens, kset = sample('AAA', 50 + s)
        val = val_of(C.evaluate(C.I3_terms['gauge'], modes), gens)
        a = modes

        def iota(q, m):
            return sum((q[mu] * m.u[mu] for mu in C.V4), ZERO)

        def brpair(m1, m2, m3):  # (iota_{.}A(p3), [A(p1),A(p2)]) colour: (X3,[X1,X2]) * (u1.u2)
            return C.dot(m1.u, m2.u) * C.f3(m3.X, m1.X, m2.X)
        tot = ZERO
        for (i, j, l) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            q = [a[j].k[mu] - a[i].k[mu] for mu in range(4)]
            tot = tot + iota(q, a[l]) * brpair(a[i], a[j], a[l])
        tot = tot * GQ(0, Fr(-1, 2))
        okv = okv and (val - tot).iszero()
    out['checks']['AAA_vertex_equals_eq_c1_AAA'] = okv
    log('C1 AAA vertex', okv)

    # ---------------- C2: {S_2, F} = Q F ; C3: Q I_3 = 0
    c2 = {}
    c3 = {}
    for sec in ('c_AV_AV', 'c_AW_AW', 'c_lam_lam', 'c_c_AsV'):
        modes, gens, kset = sample(sec, 7)
        cache = internal_for(kset)
        for name, F in list(C.I3_terms.items()) + list(C.H.items()):
            q = val_of(C.applyQ(F, modes), gens)
            b = val_of(C.bracket(C.S2, 0, F, modes, kset, cache), gens)
            c2[f'{sec}:{name}'] = (q - b).iszero()
        q3 = val_of(C.applyQ(C.I3, modes), gens)
        c3[sec] = q3.iszero()
    out['checks']['bracket_S2_equals_Q'] = all(c2.values())
    out['checks']['bracket_S2_equals_Q_detail'] = c2
    out['checks']['Q_I3_zero'] = c3
    log('C2 {S2,F}=QF', all(c2.values()), 'C3 QI3=0', c3)

    # ---------------- C4: field-count lemma
    c4 = {}
    for sec, F, cnt, mult in (('AAA', C.I3_terms['gauge'], 'cntA', 3), ('c_lam_lams', C.I3_terms['clls'], 'cntl', 0),
                              ('c_lam_lams', C.I3_terms['clls'], 'cntc', -1), ('A_c_As', C.I3_terms['AcAs'], 'cntA', 0),
                              ('A_c_As', C.I3_terms['AcAs'], 'cntc', -1), ('c_c_cs', C.I3_terms['cccs'], 'cntc', -1),
                              ('lam_A_lam', C.I3_terms['llA'], 'cntl', -2), ('lam_A_lam', C.I3_terms['llA'], 'cntA', 1)):
        modes, gens, kset = sample(sec, 11)
        cache = internal_for(kset)
        b = val_of(C.bracket(C.Q2[cnt], 1, F, modes, kset, cache), gens)
        v = val_of(C.evaluate(F, modes), gens)
        c4[f'{sec}:{cnt}'] = {'bracket/F': fmt(b / v) if not v.iszero() else ('0/0' if b.iszero() else 'nonzero/0'),
                              'expected': mult}
    out['checks']['field_count_lemma'] = c4
    log('C4 count', c4)
    # ---------------- C5: quartic master equation  Q I_4 + 1/2 {I_3, I_3} = 0 on 4 legs
    c5 = {}
    for sec in ('cAAA4', 'cAWW4'):
        modes, gens, kset = sample(sec, 21)
        cache = internal_for(kset)
        q4 = val_of(C.applyQ(C.I4, modes), gens)
        b33 = val_of(C.bracket(C.I3, 0, C.I3, modes, kset, cache), gens)
        tot = q4 + b33 * Fr(1, 2)
        c5[sec] = {'QI4': fmt(q4), 'half_I3I3': fmt(b33 * Fr(1, 2)), 'sum_zero': tot.iszero()}
    out['checks']['quartic_CME'] = c5
    log('C5 quartic CME', c5)

    # ---------------- results: Q H_j
    res = {}
    for j in range(1, 7):
        name = f'H{j}'
        res[f'Q{name}'] = decompose(f'Q {name}', lambda m, g, k, F=C.H[name]: val_of(C.applyQ(F, m), g),
                                    ['c_AV_AV', 'c_AW_AW', 'c_lam_lam', 'c_c_AsV'])
    # brackets {I_3, X}
    for X in ('D1', 'B1', 'B2', 'B3', 'B4', 'mass', 'massH'):
        def fb(m, g, k, X=X):
            return val_of(C.bracket(C.I3, 0, C.Q2[X], m, k, internal_for(k)), g)
        res[f'{{I3,{X}}}'] = decompose(f'{{I3,{X}}}', fb, ['c_AV_AV', 'c_AW_AW', 'c_lam_lam', 'c_c_AsV'])
    # parity-odd candidate and weight-2 structures
    res['K5odd_vs_G5'] = decompose('K5odd in G5 basis (expect residual nonzero: independent)',
                                   lambda m, g, k: val_of(C.evaluate(K5odd, m), g), ['c_lam_lam'])
    for X in ('B3', 'B4'):
        def fb(m, g, k, X=X):
            return val_of(C.bracket(C.I3, 0, C.Q2[X], m, k, internal_for(k)), g)
        res[f'{{I3,{X}}}_with_odd'] = decompose(f'{{I3,{X}}} with odd', fb, ['c_lam_lam'],
                                                extra_basis={'c_lam_lam': ['K5odd']})
    res['QH5_with_odd'] = decompose('QH5 with odd', lambda m, g, k: val_of(C.applyQ(C.H['H5'], m), g), ['c_lam_lam'],
                                    extra_basis={'c_lam_lam': ['K5odd']})
    w2 = {}
    for sec, nm_ in (('c_AV_AV', 'w2_AV'), ('c_AW_AW', 'w2_AW')):
        vals = []
        for s in range(3):
            modes, gens, kset = sample(sec, 300 + s)
            vals.append(val_of(C.evaluate(W2[nm_], modes), gens).iszero())
        w2[nm_] = all(vals)
    out['checks']['weight2_structures_vanish_identically'] = w2
    log('weight-2 vanish', w2)
    out['results'] = res

    # ---------------- WZ: Q G_i and {I_3, Ob_2} on the degree-2 sector (c, c, A_V)
    wz = {}
    for s in range(3):
        modes, gens, kset = sample('deg2_c_c_AV', 400 + s)
        cache = internal_for(kset)
        qg = {b: val_of(C.applyQ(C.G[b], modes), gens) for b in C.G}
        bo = val_of(C.bracket(C.I3, 0, C.Q2['Ob2'], modes, kset, cache), gens)
        wz[s] = {**{f'Q{b}': fmt(v) for b, v in qg.items()}, '{I3,Ob2int}': fmt(bo),
                 '{I3,Ob2int}-QG2': fmt(bo - qg['G2']), '{I3,Ob2int}+QG2': fmt(bo + qg['G2'])}
    out['wz_degree2'] = wz
    log('WZ', json.dumps(wz))

    # ---------------- trees through P(0, infinity)
    trees = {}
    # check Q tree(D1, I3) = tree(Q D1, I3) + {D1, I3}  (Q P(0,inf) = K_0)
    chk = {}
    for sec in ('c_AV_AV', 'c_c_AsV'):
        modes, gens, kset = sample(sec, 500)
        cache = internal_for(kset)
        dPt = C.dP_terms(kset, cache[1])
        r = C.new_gen(1, 'r')
        qm = [q for q in (C.Qimage(m, r) for m in modes) if q is not None]
        # tree(D1, I3) needs one A^* external absent here; evaluate tree on modes + Q images, then d/dr
        lhs = val_of(C.dleft(C.tree(C.Q2['D1'], 0, C.I3, modes + qm, kset, cache, dPt), r), gens)
        # tree(QD1, I3): QD1 is a degree-1 quadratic functional = -2 (d^*dc, d^*A) = -2 Ob2int
        rhs1 = val_of(C.tree(C.Q2['Ob2'], 1, C.I3, modes, kset, cache, dPt), gens) * GQ(-2)
        rhs2 = val_of(C.bracket(C.Q2['D1'], 0, C.I3, modes, kset, cache), gens)
        chk[sec] = {'Q_tree': fmt(lhs), 'tree_QD1': fmt(rhs1), '{D1,I3}': fmt(rhs2),
                    'identity_holds': (lhs - rhs1 - rhs2).iszero()}
    trees['check_Q_tree_identity'] = chk
    log('tree identity', json.dumps(chk))

    def ftree(m, g, k):
        cache = internal_for(k)
        dPt = C.dP_terms(k, cache[1])
        return val_of(C.tree(C.Q2['Ob2'], 1, C.I3, m, k, cache, dPt), g)
    trees['tree_Ob2int_I3'] = decompose('tree(int(d^*dc,d^*A_V), I3; P(0,inf))', ftree,
                                        ['c_AV_AV', 'c_AW_AW', 'c_lam_lam', 'c_c_AsV'])
    out['trees'] = trees
    out['seconds'] = round(time.time() - t0, 1)
    out['log'] = LOG
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    log('done', out['seconds'], 's')


if __name__ == '__main__':
    main()
