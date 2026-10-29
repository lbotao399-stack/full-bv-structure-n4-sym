"""Audit check: with K_infty = 0 (H invertible) the free differential Q on polynomial functionals is acyclic in every
polynomial degree >= 1 (P(0,infty) = Q^GF H^{-1} is a contracting homotopy).  Same toy complex as toy_gf.py.
Computes dim ker Q - dim im Q on the space of homogeneous polynomials of degree 1, 2, 3."""
import json, itertools
import sympy as sp
import toy_gf as T     # reuses the toy (prints its own checks)

alg, N, Qfun = T.alg, T.N, T.Qfun

def monomials(deg):
    out = []
    for combo in itertools.combinations_with_replacement(range(N), deg):
        odd = [i for i in combo if T.par[i] == 1]
        if len(odd) != len(set(odd)):
            continue
        m = alg.one()
        for i in combo:
            m = alg.mul(m, alg.gen(i))
        out.append(m)
    return out

res = {}
for deg in (1, 2, 3):
    mons = monomials(deg)
    keys = [next(iter(m)) for m in mons]
    idx = {k: n for n, k in enumerate(keys)}
    Qmat = sp.zeros(len(keys), len(keys))
    for col, m in enumerate(mons):
        cm = m[keys[col]]
        img = Qfun(m)
        for k, v in img.items():
            Qmat[idx[k], col] = v / cm
    # normalise columns by the monomial's own coefficient sign
    rank = Qmat.rank()
    # cohomology dimension = dim ker - dim im ; Q maps the degree-deg space to itself
    dimker = len(keys) - rank
    res[f'degree {deg}'] = {'dim': len(keys), 'rank Q': int(rank), 'dim H = dim ker - dim im': int(dimker - rank)}
    print(deg, res[f'degree {deg}'], flush=True)
json.dump(res, open('toy_acyclic.json', 'w'), indent=1)
