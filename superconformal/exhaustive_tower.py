"""
exhaustive_tower.py -- components of  d_CE S_g + 1/2 {S_g, S_g} = 0  at arities 1, 2, 3 on the full
62-element basis of psl(4|4) (structure.basis()), for a gauge algebra G in {su2, su3} (model_g.py).

The component at (x_1..x_n) is evaluated with tower.Tower (the adopted [ESW Def 1.26] conventions, every shuffle
computed), with S^(1), S^(2) of model_g and S^(k>=3) = 0.  Each component is tested as a local functional:
all Euler-Lagrange derivatives of the density vanish identically (exact Gaussian rationals).
Since every component is multilinear and graded antisymmetric in its arguments, the unordered basis tuples
(with repetition) exhaust it.

Usage:
  python3 exhaustive_tower.py --gauge su2 --arity 2 --shard 0/2
  python3 exhaustive_tower.py --gauge su2 --consistency        (model_g(su2) == model.py, term by term)
  python3 exhaustive_tower.py --summary                         (merge jsonl files)
Output: out/exh_<gauge>_ar<n>_<shard>.jsonl  (one line per basis tuple; resumable).
"""
import os, sys, json, time, itertools, argparse, glob, resource
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
import engine as E
import algebra as Al
import structure as St
import tower as Tw
import model_g as Mg

Al.set_cR(sp.Integer(1))

# ------------------------------------------------------------------------------------------------
# cached BV bracket: Euler derivatives of long-lived functionals (S, S^(1) on basis) computed once
# ------------------------------------------------------------------------------------------------
_EUL = {}
_KEEP = {}
_orig_bracket = E.bracket


def _eul(V, P, side, cache):
    if cache:
        k = (id(P), side)
        r = _EUL.get(k)
        if r is None:
            r = E.euler_all(V, P, side)
            _EUL[k] = r
            _KEEP[id(P)] = P
        return r
    return E.euler_all(V, P, side)


CACHED_IDS = set()


def cached_bracket(V, F, G, dimg=None):
    """identical to engine.bracket, with cached Euler derivatives for registered functionals."""
    s = E.parity(V, F)
    Fr_ = _eul(V, F, 'r', id(F) in CACHED_IDS)
    Gl_ = _eul(V, G, 'l', id(G) in CACHED_IDS)
    out = {}
    for f, fs in E.ANTI.items():
        e = E.EPS_BRACKET[f]
        se = (Fr(e), Fr(0))
        keysF = set(k for k in Fr_ if k[0] in (f, fs))
        keysG = set(k for k in Gl_ if k[0] in (f, fs))
        comps = set((k[1], k[2]) for k in keysF | keysG)
        for (a, col) in comps:
            x1 = Fr_.get((f, a, col)); y1 = Gl_.get((fs, a, col))
            if x1 and y1:
                E.padd_into(out, E.pmul(x1, y1), se)
            x2 = Fr_.get((fs, a, col)); y2 = Gl_.get((f, a, col))
            if x2 and y2:
                E.padd_into(out, E.pmul(x2, y2), (Fr(-e), Fr(0)))
    return E.pscale(out, -1) if s else out


E.bracket = cached_bracket   # tower.Tower calls E.bracket


def key_of(x):
    if isinstance(x, Al.Even):
        return ('e',) + tuple(str(c) for c in x.c)
    return ('o',) + tuple(str(c) for c in list(x.Q) + list(x.eta))


class Ctx:
    def __init__(self, G):
        self.G = G
        self.V = E.Vars()
        self.S = Mg.build_S(self.V, G)
        CACHED_IDS.add(id(self.S))
        self.c1 = {}
        self.c2 = {}
        self.bkeys = set()
        self.cache_s2 = False
        self.T = Tw.Tower(self.V, self.S, {1: self.S1, 2: self.S2})

    def register_basis(self, els):
        self.bkeys = set(key_of(x) for x in els)

    def S1(self, args):
        (x,) = args
        k = key_of(x)
        cache = k in self.bkeys
        r = self.c1.get(k)
        if r is None:
            if isinstance(x, Al.Even):
                vf, R = Al.even_to_engine(self.V, x)
                r = Mg.S1_even(self.V, self.G, vf, R)
            else:
                Q, eta = Al.odd_to_engine(x)
                r = Mg.S1_odd(self.V, self.G, Q, eta)
            if cache:
                self.c1[k] = r
                CACHED_IDS.add(id(r))
        return r

    def S2(self, args):
        x, y = args
        if not (isinstance(x, Al.Odd) and isinstance(y, Al.Odd)):
            return {}
        k = (key_of(x), key_of(y))
        r = self.c2.get(k)
        if r is None:
            q1, e1 = Al.odd_to_engine(x)
            q2, e2 = Al.odd_to_engine(y)
            r = Mg.S2_odd(self.V, self.G, Mg.killing_spinor(self.V, q1, e1), Mg.killing_spinor(self.V, q2, e2))
            if self.cache_s2 and k[0] in self.bkeys and k[1] in self.bkeys:
                self.c2[k] = r
        return r


def legs(V, P):
    """max number of field/antifield factors (all g-valued) in a monomial of P."""
    m = 0
    for (e, o) in P:
        n = sum(1 for v in e if V.kind[v] == 'f') + sum(1 for v in o if V.kind[v] == 'f')
        m = max(m, n)
    return m


def sector_test(V, P):
    """P is a total derivative iff each component of fixed field content (multiset of field names; preserved by
    the total derivatives) is; tested sector by sector to bound memory."""
    sec = {}
    for key, c in P.items():
        e, o = key
        lab = tuple(sorted(V.key[v][0] for v in e + o if V.kind[v] == 'f'))
        sec.setdefault(lab, {})[key] = c
    ok = True
    for lab in list(sec):
        r, bad = E.total_derivative_test(V, sec.pop(lab))
        if not r:
            ok = False
    return ok


def check(ctx, args):
    t = time.time()
    tot, parts = ctx.T.MC(args, split=True)
    ml = max([legs(ctx.V, p) for p in parts.values() if p] + [0])
    nzp = sorted(k for k, v in parts.items() if v)
    del parts
    ok = sector_test(ctx.V, tot)
    return {'holds': ok, 'density_terms': len(tot), 'max_legs_in_pieces': ml, 'nonzero_pieces': nzp,
            'seconds': round(time.time() - t, 2)}


def basis_list():
    B = St.basis()
    order = ['P', 'K', 'M', 'Mbar', 'Delta', 'R', 'Q', 'Qbar', 'Sbar', 'S']
    els, names = [], []
    for n in order:
        for e in B[n]:
            els.append(e)
            names.append(n)
    return els, names


def tuples(arity, els):
    N = len(els)
    par = [Al.parity(x) for x in els]
    if arity == 1:
        return [(i,) for i in range(N)]
    if arity == 2:
        return [(i, j) for i in range(N) for j in range(i, N)]
    if arity == 3:
        # unordered triples with repetition containing at least two odd elements; the others vanish
        # structurally (S^(2) vanishes on an even argument, [even, even] is even, S^(3) = 0)
        out = []
        for (i, j, k) in itertools.combinations_with_replacement(range(N), 3):
            if par[i] + par[j] + par[k] >= 2:
                out.append((i, j, k))
        return out
    raise ValueError


def sample_tuples(arity, els, names, per_pattern, seed=11):
    """per_pattern random basis tuples for every unordered pattern of [SW] subspace names (structurally
    nonzero patterns only)."""
    import random
    rng = random.Random(seed)
    groups = {}
    for tup in tuples(arity, els):
        groups.setdefault(tuple(sorted(names[i] for i in tup)), []).append(tup)
    out = []
    for pat in sorted(groups):
        g = groups[pat]
        out += g if len(g) <= per_pattern else rng.sample(g, per_pattern)
    return out, len(groups)


def gauge(name):
    return Mg.su2() if name == 'su2' else Mg.su3()


def run(gname, arity, shard, nshard, limit=None, per_pattern=0):
    G = gauge(gname)
    ctx = Ctx(G)
    els, names = basis_list()
    ctx.register_basis(els)
    ctx.cache_s2 = (arity == 3 and gname == 'su2' and not per_pattern)
    if per_pattern:
        tl, npat = sample_tuples(arity, els, names, per_pattern)
        print('sample: %d tuples over %d patterns' % (len(tl), npat), flush=True)
        fn = os.path.join(HERE, 'out', 'exhsample_%s_ar%d_%d.jsonl' % (gname, arity, shard))
    else:
        tl = tuples(arity, els)
        fn = os.path.join(HERE, 'out', 'exh_%s_ar%d_%d.jsonl' % (gname, arity, shard))
    done = set()
    if os.path.exists(fn):
        for line in open(fn):
            try:
                done.add(tuple(json.loads(line)['idx']))
            except Exception:
                pass
    f = open(fn, 'a')
    n = 0
    t0 = time.time()
    for ti, tup in enumerate(tl):
        if ti % nshard != shard or tup in done:
            continue
        r = check(ctx, [els[i] for i in tup])
        r['idx'] = list(tup)
        r['types'] = [names[i] for i in tup]
        f.write(json.dumps(r) + '\n')
        f.flush()
        n += 1
        if not r['holds']:
            print('FAIL', r, flush=True)
        if n % 50 == 0:
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20
            print('%s ar%d shard %d: %d done, %.0fs, maxrss %.0f MB' % (gname, arity, shard, n, time.time() - t0, rss), flush=True)
        if limit and n >= limit:
            break
    f.close()
    print('finished', gname, arity, shard, n, round(time.time() - t0, 1), flush=True)


def consistency():
    """model_g with su(2) reproduces model.py (the functionals used in verify.json) exactly."""
    import model as Mo
    import verify as Vf
    G = Mg.su2()
    V = E.Vars()
    out = {'S': E.psub(Mg.build_S(V, G), Mo.build_S(V)) == {}}
    g = Vf.Gen(77)
    for t in ['P', 'M', 'D', 'K', 'R']:
        x = g.make(t)
        vf, R = Al.even_to_engine(V, x)
        out['S1_' + t] = E.psub(Mg.S1_even(V, G, vf, R), Mo.S1_even(V, vf, R, Vf.KL, Vf.WL, wA=Vf.WA)) == {}
    for t in ['Q', 'S']:
        x = g.make(t)
        q, e = Al.odd_to_engine(x)
        out['S1_' + t] = E.psub(Mg.S1_odd(V, G, q, e), Mo.S1_odd(V, q, e, aS=Vf.AS)) == {}
    for a, b in [('Q', 'S'), ('S', 'S')]:
        x, y = g.make(a), g.make(b)
        q1, e1 = Al.odd_to_engine(x)
        q2, e2 = Al.odd_to_engine(y)
        k1, k2 = Mo.killing_spinor(V, q1, e1), Mo.killing_spinor(V, q2, e2)
        out['S2_' + a + b] = E.psub(Mg.S2_odd(V, G, k1, k2), Mo.S2_odd(V, k1, k2)) == {}
    return out


def free_lie_multilinear_dim(n):
    """dimension of the multilinear part of the free Lie algebra on x_1..x_n, computed as the span of all
    left-normed bracketings [[x_s1, x_s2], ...] expanded in the free associative algebra (exact)."""
    def br(P, Q):
        out = {}
        for a, ca in P.items():
            for b, cb in Q.items():
                out[a + b] = out.get(a + b, 0) + ca * cb
                out[b + a] = out.get(b + a, 0) - ca * cb
        return {k: v for k, v in out.items() if v}
    vecs = []
    for perm in itertools.permutations(range(1, n + 1)):
        P = {(perm[0],): 1}
        for x in perm[1:]:
            P = br(P, {(x,): 1})
        vecs.append(P)
    words = sorted(set(w for P in vecs for w in P))
    M = sp.Matrix([[P.get(w, 0) for w in words] for P in vecs])
    return int(M.rank())


def tree_tensors(G, k):
    """all trivalent trees with k = 4 or 5 labelled legs, as tensors in (g^*)^{(x)k} with lowered indices
    f_abc = ([T_a,T_b],T_c) and internal edges contracted with the inverse pairing; returns list of dicts."""
    import numpy as np
    d = G.dim
    fl = np.zeros((d, d, d), dtype=object)
    fl[:] = Fr(0)
    for (j, kk, l), v in G.FL:
        fl[j, kk, l] += v * G.k[l]
    ki = np.array([G.kinv[e] for e in range(d)], dtype=object)
    flk = fl * ki[None, None, :]          # f_ab^p with the last index raised
    trees = []
    if k == 4:
        T = np.einsum('abp,pcd->abcd', flk, fl)
        for perm in [(0, 1, 2, 3), (0, 2, 1, 3), (0, 3, 1, 2)]:
            trees.append(np.transpose(T, np.argsort(perm)))
    else:
        T = np.einsum('abp,pcq,qde->abcde', flk, flk, fl)
        for c in range(5):
            rest = [x for x in range(5) if x != c]
            for (p1, p2) in [((rest[0], rest[1]), (rest[2], rest[3])), ((rest[0], rest[2]), (rest[1], rest[3])),
                             ((rest[0], rest[3]), (rest[1], rest[2]))]:
                perm = (p1[0], p1[1], c, p2[0], p2[1])
                trees.append(np.transpose(T, np.argsort(perm)))
    return [t.reshape(-1) for t in trees]


def gram_rank(vecs):
    n = len(vecs)
    M = sp.Matrix(n, n, lambda i, j: sp.Rational(sum(vecs[i] * vecs[j])))
    return int(M.rank())


def trees():
    """color structures.  Every term of the arity 1-3 components is (g-independent coefficient) x (connected tree
    of the 3-tensor ([-,-],-) and the pairing) with <= 4 legs; the arity-0 component {S,S} has <= 5 legs.
    The span of the k-leg tree tensors in ANY g is a quotient of the multilinear free Lie algebra on k-1 letters
    (a tree is (P(X_1..X_{k-1}), X_k) for a Lie monomial P), of dimension (k-2)!.  If the rank in one g equals that
    dimension, evaluation on that g is injective on trees, and an identity among tree tensors in that g holds in
    every g."""
    res = {'free_Lie_multilinear_dimension': {'1 letter (2 legs)': 1, '2 letters (3 legs)': free_lie_multilinear_dim(2),
                                              '3 letters (4 legs)': free_lie_multilinear_dim(3),
                                              '4 letters (5 legs)': free_lie_multilinear_dim(4)}}
    for gname in ('su2', 'su3'):
        G = gauge(gname)
        r = {}
        r['rank_2leg_trees(pairing)'] = int(any(x != 0 for x in G.k))
        r['rank_3leg_trees(([-,-],-))'] = int(any(v != 0 for _, v in G.FL))
        r['rank_4leg_trees(3 trees)'] = gram_rank(tree_tensors(G, 4))
        r['rank_5leg_trees(15 trees)'] = gram_rank(tree_tensors(G, 5))
        r.update(G.lie_tests())
        res[gname] = r
    return res


def symcheck(gname, nsample=40, seed=5):
    """justifies unordered basis tuples: (i) S^(2)(o,o') = S^(2)(o',o) on all ordered pairs of odd basis elements;
    (ii) the tower component is graded antisymmetric at density level,
         F(..,x_i,x_{i+1},..) = (-1)^{|x_i||x_{i+1}|+1} F(..,x_{i+1},x_i,..),
    on random ordered basis pairs and triples (adjacent transpositions)."""
    import random
    G = gauge(gname)
    ctx = Ctx(G)
    els, names = basis_list()
    ctx.register_basis(els)
    odd = [k for k, x in enumerate(els) if Al.parity(x)]
    bad = 0
    for i in odd:
        for j in odd:
            if j <= i:
                continue
            if E.psub(ctx.S2([els[i], els[j]]), ctx.S2([els[j], els[i]])) != {}:
                bad += 1
    out = {'S2_symmetric_on_all_odd_basis_pairs': {'pairs': len(odd) * (len(odd) - 1) // 2, 'violations': bad}}
    rng = random.Random(seed)
    N = len(els)
    bad2 = 0
    for _ in range(nsample):
        i, j = rng.randrange(N), rng.randrange(N)
        x, y = els[i], els[j]
        sg = (-1) ** (Al.parity(x) * Al.parity(y) + 1)
        if E.psub(ctx.T.MC([x, y]), E.pscale(ctx.T.MC([y, x]), sg)) != {}:
            bad2 += 1
    bad3 = 0
    n3 = 0
    while n3 < nsample:
        i, j, k = rng.randrange(N), rng.randrange(N), rng.randrange(N)
        tri = [els[i], els[j], els[k]]
        if sum(Al.parity(x) for x in tri) < 2:
            continue
        n3 += 1
        p = rng.randrange(2)
        sw = list(tri)
        sw[p], sw[p + 1] = sw[p + 1], sw[p]
        sg = (-1) ** (Al.parity(tri[p]) * Al.parity(tri[p + 1]) + 1)
        if E.psub(ctx.T.MC(tri), E.pscale(ctx.T.MC(sw), sg)) != {}:
            bad3 += 1
    out['graded_antisymmetry_of_component_density'] = {'random_ordered_pairs': nsample, 'violations_arity2': bad2,
                                                       'random_triples_adjacent_swap': nsample, 'violations_arity3': bad3}
    return out


def cme(gname):
    """arity-0 component: (1/2){S,S} = 0 as a local functional (classical master equation)."""
    G = gauge(gname)
    V = E.Vars()
    S = Mg.build_S(V, G)
    t = time.time()
    B = _orig_bracket(V, S, S)
    ok, bad = E.total_derivative_test(V, B)
    return {'holds': ok, 'density_terms': len(B), 'max_legs': legs(V, B), 'S_terms': len(S),
            'seconds': round(time.time() - t, 1)}


class CtxControl(Ctx):
    """negative controls: wrong Weyl weight of lam in S^(1) of even generators, or S^(2) = 0."""
    def __init__(self, G, wl=None, s2_zero=False):
        Ctx.__init__(self, G)
        self.wl_override, self.s2_zero = wl, s2_zero

    def S1(self, args):
        (x,) = args
        if self.wl_override is not None and isinstance(x, Al.Even):
            vf, R = Al.even_to_engine(self.V, x)
            return Mg.S1_even(self.V, self.G, vf, R, wl=self.wl_override)
        return Ctx.S1(self, args)

    def S2(self, args):
        return {} if self.s2_zero else Ctx.S2(self, args)


def controls(gname='su2'):
    G = gauge(gname)
    els, names = basis_list()
    out = {}
    c = CtxControl(G, wl=Fr(1, 4))
    c.register_basis(els)
    fails = [names[i] for i in range(len(els)) if not check(c, [els[i]])['holds']]
    out['arity1, Weyl weight of lam 1/4 instead of 3/8: failing basis elements'] = {n: fails.count(n) for n in sorted(set(fails))}
    odd = [i for i, x in enumerate(els) if Al.parity(x)]
    c = CtxControl(G, s2_zero=True)
    c.register_basis(els)
    nf, nt = 0, 0
    for a in range(len(odd)):
        for b in range(a, len(odd)):
            nt += 1
            if not check(c, [els[odd[a]], els[odd[b]]])['holds']:
                nf += 1
    out['arity2, S^(2) = 0: failing odd-odd basis pairs'] = '%d of %d' % (nf, nt)
    c = CtxControl(G, wl=Fr(1, 4))
    c.register_basis(els)
    iD = names.index('Delta')
    nf, nt = 0, 0
    for a in range(len(odd)):
        for b in range(a, len(odd)):
            nt += 1
            if not check(c, [els[iD], els[odd[a]], els[odd[b]]])['holds']:
                nf += 1
    out['arity3, Weyl weight of lam 1/4 instead of 3/8: failing (Delta, odd, odd) basis triples'] = '%d of %d' % (nf, nt)
    return out


def summary():
    res = {}
    for fn in sorted(glob.glob(os.path.join(HERE, 'out', 'exh_*_ar*_*.jsonl')) +
                     glob.glob(os.path.join(HERE, 'out', 'exhsample_*_ar*_*.jsonl'))):
        base = os.path.basename(fn)
        g, ar = base.split('_')[1], base.split('_')[2]
        key = '%s_%s%s' % (g, ar, '_sample' if base.startswith('exhsample') else '')
        d = res.setdefault(key, {'checked': 0, 'failures': [], 'max_legs': 0, 'seconds': 0.0, 'by_types': {}})
        for line in open(fn):
            r = json.loads(line)
            d['checked'] += 1
            d['seconds'] += r['seconds']
            d['max_legs'] = max(d['max_legs'], r['max_legs_in_pieces'])
            tkey = ''.join(t[0] if t not in ('Mbar', 'Qbar', 'Sbar') else t[0].lower() for t in r['types'])
            d['by_types'][tkey] = d['by_types'].get(tkey, 0) + 1
            if not r['holds']:
                d['failures'].append(r['idx'])
    for d in res.values():
        d['seconds'] = round(d['seconds'], 1)
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--gauge', default='su2')
    ap.add_argument('--arity', type=int, default=1)
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--consistency', action='store_true')
    ap.add_argument('--trees', action='store_true')
    ap.add_argument('--summary', action='store_true')
    ap.add_argument('--count', action='store_true')
    ap.add_argument('--symcheck', action='store_true')
    ap.add_argument('--cme', action='store_true')
    ap.add_argument('--controls', action='store_true')
    ap.add_argument('--per_pattern', type=int, default=0)
    a = ap.parse_args()
    if a.consistency:
        print(json.dumps(consistency()))
    elif a.trees:
        r = trees()
        json.dump(r, open(os.path.join(HERE, 'out', 'exh_trees.json'), 'w'), indent=1)
        print(json.dumps(r))
    elif a.summary:
        print(json.dumps(summary(), indent=1))
    elif a.symcheck:
        r = symcheck(a.gauge)
        json.dump(r, open(os.path.join(HERE, 'out', 'exh_symcheck_%s.json' % a.gauge), 'w'), indent=1)
        print(json.dumps(r))
    elif a.controls:
        r = controls(a.gauge)
        json.dump(r, open(os.path.join(HERE, 'out', 'exh_controls_%s.json' % a.gauge), 'w'), indent=1)
        print(json.dumps(r))
    elif a.cme:
        r = cme(a.gauge)
        json.dump(r, open(os.path.join(HERE, 'out', 'exh_cme_%s.json' % a.gauge), 'w'), indent=1)
        print(json.dumps(r))
    elif a.count:
        els, names = basis_list()
        print({n: len(tuples(n, els)) for n in (1, 2, 3)})
    else:
        sh, ns = [int(x) for x in a.shard.split('/')]
        run(a.gauge, a.arity, sh, ns, a.limit or None, a.per_pattern)
