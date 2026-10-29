"""Audit script a1b: which shuffle sign rule makes the product of functional-valued cochains
associative AND graded commutative (as the product of C^bullet(g) (x) O(E) must be), and does r_xi
stay multiplicative?  Variants of the exponent in Def c2-dce with the bracket replaced by the product:
  V0: |g|k + sum_{i<=k}|x_sigma(i)|(l+|g|)                (Def c2-dce verbatim, |g| = C^bullet degree)
  V1: as V0 with |g| the total parity of g
  V2: V1 + k*l
  V3: |g|k + sum_{i<=k}|x_sigma(i)|(l+|g|) with |g|=parity, plus k*(l+1)
Commutativity is tested in the form f g = (-1)^{P_f P_g} g f with P = total parity of the cochain
(value parity minus arity minus sum of argument parities), cochain by cochain of fixed arities.
"""
import json, random, sys, os, itertools, math
from fractions import Fraction
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a1_cochain_decalage as M

A = M.A
rng = random.Random(11)


def combine_var(f, g, variant):
    def val(tup):
        n = len(tup)
        out = {}
        for k in f.comp:
            l = n - k
            if l not in g.comp or l < 0:
                continue
            for order in M.shuffles(n, k):
                xs = [tup[a] for a in order]
                sgn, kos = M.perm_signs([M.gpar[t] for t in tup], order)
                sx = sum(M.gpar[xs[i]] for i in range(k))
                if variant == 'V0':
                    e = g.D * k + sx * (l + g.D)
                elif variant == 'V1':
                    e = g.p * k + sx * (l + g.p)
                elif variant == 'V2':
                    e = g.p * k + sx * (l + g.p) + k * l
                elif variant == 'V3':
                    e = g.p * k + sx * (l + g.p) + k * (l + 1)
                s = sgn * kos * (-1) ** (e % 2)
                fv = f.value(tuple(xs[:k])); gv = g.value(tuple(xs[k:]))
                if fv and gv:
                    out = A.add(out, A.mul(fv, gv), s)
        return out
    return val


res = {}
for variant in ('V0', 'V1', 'V2', 'V3'):
    comm_ok, assoc_ok, rxi_ok, tot = 0, 0, 0, 0
    for trial in range(10):
        ka, kb = rng.choice([[1], [2], [0, 1], [1, 2]]), rng.choice([[1], [2], [0, 1]])
        Df, Dg = rng.choice([0, -1, 1]), rng.choice([0, -1, 1])
        pf, pg = rng.randint(0, 1), rng.randint(0, 1)
        f = M.random_cochain(Df, pf, ka); g = M.random_cochain(Dg, pg, kb)
        fg = combine_var(f, g, variant); gf = combine_var(g, f, variant)
        tot += 1
        ok = True
        for n in range(0, 4):
            for tup in itertools.product(range(M.NB), repeat=n):
                u, v = fg(tup), gf(tup)
                if not A.is_zero(A.add(u, v, -((-1) ** (pf * pg)))):
                    ok = False
        comm_ok += ok
        # associativity with a third cochain
        h = M.random_cochain(rng.choice([0, -1]), rng.randint(0, 1), rng.choice([[1], [0, 1]]))
        w = M.Cochain(f.D + g.D, (f.p + g.p) % 2)
        w.comp = {a_ + b_: None for a_ in f.comp for b_ in g.comp}
        w.value = fg
        L = combine_var(w, h, variant)
        u2 = M.Cochain(g.D + h.D, (g.p + h.p) % 2)
        u2.comp = {a_ + b_: None for a_ in g.comp for b_ in h.comp}
        u2.value = combine_var(g, h, variant)
        R = combine_var(f, u2, variant)
        ok = True
        for n in range(0, 4):
            for tup in itertools.product(range(M.NB), repeat=n):
                if not A.is_zero(A.add(L(tup), R(tup), -1)):
                    ok = False
        assoc_ok += ok
        # r_xi multiplicativity with a_k
        ok = True
        for n in range(0, 4):
            lhs = A.scale(fg(tuple([1] * n)), M.a(n))
            rhs = {}
            for k in range(0, n + 1):
                fv = f.value(tuple([1] * k)); gv = g.value(tuple([1] * (n - k)))
                if fv and gv:
                    rhs = A.add(rhs, A.mul(fv, gv), M.a(k) * M.a(n - k))
            if not A.is_zero(A.add(lhs, rhs, -1)):
                ok = False
        rxi_ok += ok
    res[variant] = {'graded_commutative_(-1)^{P_f P_g}': f'{comm_ok}/{tot}', 'associative': f'{assoc_ok}/{tot}', 'r_xi_multiplicative_with_a_k': f'{rxi_ok}/{tot}'}

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a1b_product_rule.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
