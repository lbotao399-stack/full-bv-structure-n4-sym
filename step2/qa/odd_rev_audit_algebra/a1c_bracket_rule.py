"""Audit script a1c: the cochain bracket of Def c2-dce with |g| read as (V0) the C^bullet degree or
(V1) the total parity of g, for cochains whose values may be fermionic (value parity offset p != D mod 2).
Tests: graded symmetry, Jacobi (all sign patterns fitted), d_CE derivation of the bracket.
For cochains with p = D mod 2 the two readings coincide; S_psl(4|4) is of this kind, a fermionic
operator such as int omega (zeta, lambda) (degree 0, odd) is not.
"""
import json, random, sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import a1_cochain_decalage as M

A = M.A
rng = random.Random(5)


def cbr(f, g, variant):
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
                G = g.D if variant == 'V0' else g.p
                e = G * k + sx * (l + G)
                s = sgn * kos * (-1) ** (e % 2)
                fv = f.value(tuple(xs[:k])); gv = g.value(tuple(xs[k:]))
                if fv and gv:
                    out = A.add(out, A.br(fv, gv), s)
        return out
    return val


def dCEv(f, variant):
    # eq:c2-dce with |f| read as C^bullet degree (V0) or total parity (V1)
    h = M.Cochain(f.D + 1, None)
    F = f.D if variant == 'V0' else f.p

    def val(tup):
        n = len(tup)
        out = {}
        if n - 1 not in f.comp:
            return out
        for i in range(n):
            for j in range(i + 1, n):
                xi, xj = tup[i], tup[j]
                e = M.gpar[xi] * sum(M.gpar[tup[q]] for q in range(i)) + M.gpar[xj] * sum(M.gpar[tup[q]] for q in range(j) if q != i) + (i + 1) + (j + 1) + F + 1
                b = M.gbr(xi, xj)
                rest = [tup[q] for q in range(n) if q not in (i, j)]
                for bi, bc in b.items():
                    out = A.add(out, f.value(tuple([bi] + rest)), bc * (-1) ** (e % 2))
        return out
    h.val = val
    return h


def wrap(valfun, D, p, arities):
    w = M.Cochain(D, p)
    w.comp = {k: None for k in arities}
    w.value = valfun
    return w


def test_variant(variant, fermionic, trials=8):
    sym_found, jac_ok, der_ok = [], 0, 0
    der_patterns = []
    for t in range(trials):
        def mk():
            D = rng.choice([0, -1, 1])
            p = (D + (rng.randint(0, 1) if fermionic else 0)) % 2
            return M.random_cochain(D, p, rng.choice([[0, 1], [1], [0, 1, 2]]))
        f, g, h = mk(), mk(), mk()
        fg = cbr(f, g, variant); gf = cbr(g, f, variant)
        # symmetry: find sign s with {g,f} = s {f,g}
        signs = set()
        for n in range(0, 4):
            for tup in itertools.product(range(M.NB), repeat=n):
                u, v = fg(tup), gf(tup)
                if A.is_zero(u) and A.is_zero(v):
                    continue
                if A.is_zero(A.add(v, u, -1)):
                    signs.add(1)
                elif A.is_zero(A.add(v, u, 1)):
                    signs.add(-1)
                else:
                    signs.add(0)
        expected = (-1) ** (f.p * g.p)
        sym_found.append(signs <= {expected})
        # Jacobi {f,{g,h}} = alpha {{f,g},h} + beta {g,{f,h}}: with alpha = beta = 1 for odd-shifted? fit
        fgw = wrap(fg, f.D + g.D + 1, (f.p + g.p + 1) % 2, [a_ + b_ for a_ in f.comp for b_ in g.comp])
        ghw = wrap(cbr(g, h, variant), g.D + h.D + 1, (g.p + h.p + 1) % 2, [a_ + b_ for a_ in g.comp for b_ in h.comp])
        fhw = wrap(cbr(f, h, variant), f.D + h.D + 1, (f.p + h.p + 1) % 2, [a_ + b_ for a_ in f.comp for b_ in h.comp])
        L = cbr(f, ghw, variant); T1 = cbr(fgw, h, variant); T2 = cbr(g, fhw, variant)
        good = False
        for al in (1, -1):
            for be in (1, -1):
                ok = True
                for n in range(0, 4):
                    for tup in itertools.product(range(M.NB), repeat=n):
                        if not A.is_zero(A.add(A.add(L(tup), T1(tup), -al), T2(tup), -be)):
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    good = True
        jac_ok += good
        # d_CE derivation: d{f,g} = {df,g} + s {f,dg} for some sign s
        dfg = dCEv(fgw, variant)
        df = wrap(dCEv(f, variant).val, f.D + 1, (f.p + 1) % 2, [k + 1 for k in f.comp])
        dg = wrap(dCEv(g, variant).val, g.D + 1, (g.p + 1) % 2, [k + 1 for k in g.comp])
        R1 = cbr(df, g, variant); R2 = cbr(f, dg, variant)
        good = []
        for r in (1, -1):
            for s in (1, -1):
                ok = True
                for n in range(1, 4):
                    for tup in itertools.product(range(M.NB), repeat=n):
                        if not A.is_zero(A.add(A.add(dfg.val(tup), R1(tup), -r), R2(tup), -s)):
                            ok = False
                            break
                    if not ok:
                        break
                if ok:
                    good.append((r, s, f.D, f.p))
        der_ok += bool(good)
        der_patterns.append(good)
    return {'symmetry_{g,f}=(-1)^{P_fP_g}{f,g}': f'{sum(sym_found)}/{trials}', 'jacobi_some_sign_pattern': f'{jac_ok}/{trials}', 'dCE_derivation_d{f,g}=r{df,g}+s{f,dg}': f'{der_ok}/{trials}', 'patterns_(r,s,Df,pf)': der_patterns}


res = {}
for variant in ('V0', 'V1'):
    for ferm in (False, True):
        res[f'{variant}_fermionic_values={ferm}'] = test_variant(variant, ferm)
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a1c_bracket_rule.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
