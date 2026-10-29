"""Audit script a1: the evaluation map r_xi of ch08 Lemma c8-ev / design lem:o-ev, tested on
explicit cochains with the sign conventions of ch02 Def c2-cochains, Def c2-dce, eq:c2-dce.

Model: g = span{e0 (even), e1, e2 (odd)}, [e1,e2] = [e2,e1] = e0 (a Lie superalgebra);
values in the bigraded BV algebra of grass.py (bosonic/fermionic fields, antifields, ghost).
A cochain f of C^bullet-degree D and value-parity offset p has arity-k values of degree D-k and
parity p + k + sum|x_i| (ch02 Def c2-cochains, extended by an independent parity offset).

Tests
 T1 combinatorics: on (k,l)-shuffles of k+l copies of one odd element the sign
    sgn * Koszul * (-1)^{|g|k + sum_{i<=k}|x_sigma(i)|(l+|g|)} equals (-1)^{kl}  (all |g|),
    so {f,g}(xi^n) = sum binom(n,k) (-1)^{kl} {f(xi^k), g(xi^l)}, and a_k a_l = a_{k+l} binom (-1)^{kl}
    for a_k = (-1)^{k(k-1)/2}/k!.
 T2 r_xi intertwines the cochain bracket (Def c2-dce) with the bracket of values, random cochains.
 T3 r_xi intertwines the product of functional-valued cochains defined by the same shuffle rule
    (needed for e^{J/hbar} and W(P,J); not stated in ch08).
 T4 d_CE^2 = 0 (eq:c2-dce) and r_xi(d_CE f) = 0 for xi = e1 ([e1,e1]=0).
 T5 associativity of the shuffle product and its graded commutativity sign.
 T6 twisted action: arity-2 component sign (-1)^{|x1|} of Lemma c2-components re-derived from the
    shuffle formula; coefficient c = -1/2 of S^(2)(xi,xi).
"""
import json, random, sys, os, itertools, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from grass import Alg

names = ['x1', 'X1s', 't1', 'T1s', 'c', 'cs', 'x2', 'X2s']
par = [0, 1, 1, 0, 1, 0, 0, 1]
deg = [0, -1, 0, -1, 1, -2, 0, -1]
pairs = [(0, 1), (2, 3), (4, 5), (6, 7)]
A = Alg(names, par, deg, pairs)
rng = random.Random(7)

gpar = [0, 1, 1]
NB = 3


def gbr(i, j):
    """bracket of basis elements: returns dict {basis index: coefficient}"""
    if {i, j} == {1, 2}:
        return {0: 1}
    return {}


def perm_signs(seq_par, order):
    """order: list giving the original positions in the new order. returns (sgn, koszul)"""
    sgn, kos = 1, 1
    n = len(order)
    for a in range(n):
        for b in range(a + 1, n):
            if order[a] > order[b]:
                sgn = -sgn
                if seq_par[order[a]] * seq_par[order[b]] % 2:
                    kos = -kos
    return sgn, kos


class Cochain:
    def __init__(self, D, p):
        self.D = D  # C^bullet degree
        self.p = p  # parity offset of values
        self.comp = {}  # arity -> {canonical tuple: element}

    def value(self, tup):
        k = len(tup)
        if k not in self.comp:
            return {}
        order = sorted(range(k), key=lambda a: tup[a])
        can = tuple(tup[a] for a in order)
        for a in range(k - 1):
            if can[a] == can[a + 1] and gpar[can[a]] == 0:
                return {}
        val = self.comp[k].get(can, {})
        if not val:
            return {}
        # f(tup) relates to f(can): can = tup permuted by 'order'; f(x_sigma) = sgn*kos f(x)
        sgn, kos = perm_signs([gpar[t] for t in tup], order)
        return A.scale(val, sgn * kos)

    def eval_vec(self, vecs):
        """multilinear evaluation on homogeneous vectors given as dicts {basis: coeff}"""
        k = len(vecs)
        out = {}
        for idx in itertools.product(*[list(v.keys()) for v in vecs]):
            coef = 1
            for v, i in zip(vecs, idx):
                coef *= v[i]
            if coef:
                out = A.add(out, self.value(idx), coef)
        return out


def canonical_tuples(k):
    res = []
    for tup in itertools.combinations_with_replacement(range(NB), k):
        ok = all(not (tup[a] == tup[a + 1] and gpar[tup[a]] == 0) for a in range(k - 1))
        if ok:
            res.append(tup)
    return res


def random_cochain(D, p, arities, density=0.8):
    f = Cochain(D, p)
    for k in arities:
        f.comp[k] = {}
        for tup in canonical_tuples(k):
            if rng.random() > density:
                continue
            vp = (p + k + sum(gpar[t] for t in tup)) % 2
            vd = D - k
            el = A.random_hom(rng, vp, vd, nterms=2, maxdeg=3)
            if el:
                f.comp[k][tup] = el
    return f


def shuffles(n, k):
    for first in itertools.combinations(range(n), k):
        rest = [i for i in range(n) if i not in first]
        yield list(first) + rest


def combine(f, g, op, Dg_for_sign):
    """cochain bracket/product by the shuffle rule of Def c2-dce; op = A.br or A.mul"""
    h = Cochain(None, None)
    h.D = f.D + g.D + (1 if op == 'br' else 0)
    h.p = None
    arities = sorted(set(a + b for a in f.comp for b in g.comp))

    def val(tup):
        n = len(tup)
        out = {}
        for k in f.comp:
            l = n - k
            if l not in g.comp or l < 0:
                continue
            for order in shuffles(n, k):
                xs = [tup[a] for a in order]
                sgn, kos = perm_signs([gpar[t] for t in tup], order)
                e = Dg_for_sign * k + sum(gpar[xs[i]] for i in range(k)) * (l + Dg_for_sign)
                s = sgn * kos * (-1) ** (e % 2)
                fv = f.value(tuple(xs[:k])); gv = g.value(tuple(xs[k:]))
                if not fv or not gv:
                    continue
                r = A.br(fv, gv) if op == 'br' else A.mul(fv, gv)
                out = A.add(out, r, s)
        return out
    h.val = val
    return h


def dCE(f):
    h = Cochain(f.D + 1, None)

    def val(tup):
        n = len(tup)
        out = {}
        if n - 1 not in f.comp:
            return out
        for i in range(n):
            for j in range(i + 1, n):
                xi, xj = tup[i], tup[j]
                e = gpar[xi] * sum(gpar[tup[p]] for p in range(i)) + gpar[xj] * sum(gpar[tup[p]] for p in range(j) if p != i) + (i + 1) + (j + 1) + f.D + 1
                b = gbr(xi, xj)
                rest = [tup[p] for p in range(n) if p not in (i, j)]
                for bi, bc in b.items():
                    out = A.add(out, f.value(tuple([bi] + rest)), bc * (-1) ** (e % 2))
        return out
    h.val = val
    return h


def a(k):
    return Fraction((-1) ** (k * (k - 1) // 2), math.factorial(k))


res = {}

# T1: combinatorics
t1 = {}
allok = True
for n in range(1, 8):
    for k in range(0, n + 1):
        l = n - k
        for Dg in (-2, -1, 0, 1, 2):
            tot = 0
            for order in shuffles(n, k):
                sgn, kos = perm_signs([1] * n, order)
                e = Dg * k + k * (l + Dg)
                tot += sgn * kos * (-1) ** (e % 2)
            if tot != math.comb(n, k) * (-1) ** (k * l):
                allok = False
        if a(k) * a(l) != a(n) * math.comb(n, k) * (-1) ** (k * l):
            allok = False
t1['shuffle_sum_equals_binom_times_(-1)^{kl}_and_a_k_a_l_identity'] = allok
# alternative coefficients fail
alt = all(Fraction(1, math.factorial(k)) * Fraction(1, math.factorial(n - k)) == Fraction(1, math.factorial(n)) * math.comb(n, k) * (-1) ** (k * (n - k)) for n in range(2, 5) for k in range(n + 1))
t1['naive_1/k!_coefficients_intertwine'] = alt
res['T1'] = t1

xi = {1: 1}  # e1, [e1,e1] = 0


def r_xi_component(fv_func, n):
    tup = tuple([1] * n)
    return fv_func(tup)


# T2, T3
t2 = {'bracket': [], 'product': []}
for trial in range(6):
    Df, Dg = rng.choice([0, -1, 1]), rng.choice([0, -1, 1])
    pf, pg = rng.randint(0, 1), rng.randint(0, 1)
    f = random_cochain(Df, pf, [0, 1, 2])
    g = random_cochain(Dg, pg, [0, 1, 2])
    for op in ('br', 'mul'):
        h = combine(f, g, op, g.D)
        ok = True
        for n in range(0, 5):
            lhs = A.scale(h.val(tuple([1] * n)), a(n))
            rhs = {}
            for k in range(0, n + 1):
                l = n - k
                fv = f.value(tuple([1] * k)); gv = g.value(tuple([1] * l))
                if fv and gv:
                    rr = A.br(fv, gv) if op == 'br' else A.mul(fv, gv)
                    rhs = A.add(rhs, rr, a(k) * a(l))
            if not A.is_zero(A.add(lhs, rhs, -1)):
                ok = False
        t2['bracket' if op == 'br' else 'product'].append(ok)
res['T2_T3_r_xi_intertwines'] = t2

# T4: d_CE^2 = 0 and r_xi d_CE = 0
t4 = {'dCE_squared_zero': [], 'r_xi_dCE_zero': []}
for trial in range(4):
    f = random_cochain(rng.choice([0, -1]), rng.randint(0, 1), [0, 1, 2, 3])
    df = dCE(f)
    # wrap df as cochain with comp-like value function for second application
    class Wrap(Cochain):
        pass
    w = Wrap(f.D + 1, None)
    w.comp = {k + 1: None for k in f.comp}
    w.value = df.val
    ddf = dCE(w)
    ok = True
    for n in range(2, 5):
        for tup in itertools.product(range(NB), repeat=n):
            if not A.is_zero(ddf.val(tup)):
                ok = False
                break
    t4['dCE_squared_zero'].append(ok)
    ok2 = all(A.is_zero(df.val(tuple([1] * n))) for n in range(1, 5))
    t4['r_xi_dCE_zero'].append(ok2)
res['T4'] = t4

# T5: associativity and graded commutativity of the shuffle product
t5 = {'associative': [], 'commutativity_sign_found': []}
for trial in range(4):
    fs = [random_cochain(rng.choice([0, -1, 1]), rng.randint(0, 1), [0, 1, 2]) for _ in range(3)]
    f, g, h = fs

    def as_cochain(obj, D, arities):
        w = Cochain(D, None)
        w.comp = {k: None for k in arities}
        w.value = obj.val
        return w
    fg = as_cochain(combine(f, g, 'mul', g.D), f.D + g.D, [a_ + b_ for a_ in f.comp for b_ in g.comp])
    gh = as_cochain(combine(g, h, 'mul', h.D), g.D + h.D, [a_ + b_ for a_ in g.comp for b_ in h.comp])
    L = combine(fg, h, 'mul', h.D)
    R = combine(f, gh, 'mul', gh.D)
    ok = True
    for n in range(0, 4):
        for tup in itertools.product(range(NB), repeat=n):
            if not A.is_zero(A.add(L.val(tup), R.val(tup), -1)):
                ok = False
    t5['associative'].append(ok)
    # commutativity: fg vs gf
    gf = combine(g, f, 'mul', f.D)
    fgc = combine(f, g, 'mul', g.D)
    signs = set()
    for n in range(0, 4):
        for tup in itertools.product(range(NB), repeat=n):
            u, v = fgc.val(tup), gf.val(tup)
            if A.is_zero(u) and A.is_zero(v):
                continue
            if A.is_zero(A.add(u, v, -1)):
                signs.add(+1)
            elif A.is_zero(A.add(u, v, 1)):
                signs.add(-1)
            else:
                signs.add('none')
    t5['commutativity_sign_found'].append({'Df': f.D, 'pf': f.p, 'Dg': g.D, 'pg': g.p, 'signs': sorted(map(str, signs))})
res['T5'] = t5

# T6: arity-two component sign and twisted coefficient (symbolic bookkeeping)
# value parities: S^(1)(x) has parity 1+|x|; symmetry {u,v} = (-1)^{|u||v|}{v,u}
t6 = {}
for x1 in (0, 1):
    for x2 in (0, 1):
        # (1,1) shuffles of (x1,x2) with |f|=|g|=0
        c_id = (-1) ** ((x1 * 1) % 2)  # sigma = id: exponent |x1|(l+|g|) with l=1
        sgn, kos = -1, (-1) ** (x1 * x2)
        c_tr = sgn * kos * (-1) ** (x2 % 2)  # f(x2) g(x1)
        # {S1(x2),S1(x1)} = (-1)^{(1+x1)(1+x2)} {S1(x1),S1(x2)}
        c_tr *= (-1) ** ((1 + x1) * (1 + x2) % 2)
        t6[f'|x1|={x1},|x2|={x2}'] = {'coefficient_of_{S1(x1),S1(x2)}_in_half_bracket': Fraction(c_id + c_tr, 2), 'lemma_c2_components_value': (-1) ** x1}
# twisted: {S + S1 + c S2, same} = (1 + 2c){S1,S1} using {S,S2(xi,xi)} = {S1(xi),S1(xi)} and {S1,S2}=0
t6['twisted_action_coefficient_c_solving_(1+2c)=0'] = str(Fraction(-1, 2))
t6['a_2'] = str(a(2))
res['T6'] = t6

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a1_cochain_decalage.json')
json.dump(res, open(out, 'w'), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
