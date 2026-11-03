"""Engine B, action side: generic two-leg local functionals and the cohomology of Q on them (task S1, two legs).

A candidate density is a full contraction of the vectors {leg momenta P0, P1 (any number), Killing data b (once),
the position index e_mu (once, for the x^mu part)} with the leg polarisations: V-vectors (A_V, A*_V) enter dot
products and the epsilon tensor of V; scalars (c, c*) by their value; W_R^perp vectors (A_perp, A*_perp) by the
SO(6) pairing; spinors by the bilinears (Gamma lam, lam*) (Gamma an even product of rho's, Sigma -> Sigma), (Gamma
lam*_0, lam*_1) (odd product, Sigma* -> Sigma).  Densities are evaluated with eb_qb.dens (vertices.vertex convention:
both slot assignments for equal species), canonical data by eb_qb.canon.  Linear algebra over Q(i) modulo a prime
(ranks) and exactly (sympy) for the final small systems."""
import itertools, random
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import eb_qb as Q

LEVI = np.zeros((4, 4, 4, 4))
for perm in itertools.permutations(range(4)):
    LEVI[perm] = Q._perm_sign(perm)


def matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        rest = items[1:i] + items[i + 1:]
        for m in matchings(rest):
            yield [(a, items[i])] + m


def vec_of(lab, s, ctx):
    """4-vector of a label; s = slot legs, ctx = dict(b=..., mu=...)."""
    if lab == 'P0':
        return Q.pv(s[0])
    if lab == 'P1':
        return Q.pv(s[1])
    if lab == 'b':
        return ctx['b']
    if lab == 'e':
        e = np.zeros(4, complex)
        e[ctx['mu']] = 1
        return e
    if lab == 'a0':
        return s[0]['u'][:4]
    if lab == 'a1':
        return s[1]['u'][:4]
    raise ValueError(lab)


def contraction_value(pairs, eps4, s, ctx):
    v = 1 + 0j
    for (x, y) in pairs:
        v *= np.dot(vec_of(x, s, ctx), vec_of(y, s, ctx))
    if eps4 is not None:
        vs = [vec_of(x, s, ctx) for x in eps4]
        v *= np.einsum('ijkl,i,j,k,l->', LEVI, *vs)
    return v


def scalar_forms(labels):
    """all full contractions (pairs, optionally one epsilon) of a list of labels (as tuples)."""
    out = set()
    n = len(labels)
    if n % 2 == 0:
        for m in matchings(list(range(n))):
            out.add((tuple(sorted(tuple(sorted((labels[i], labels[j]))) for i, j in m)), None))
    if n >= 4 and n % 2 == 0:
        for quad in itertools.combinations(range(n), 4):
            rest = [i for i in range(n) if i not in quad]
            for m in matchings(rest):
                out.add((tuple(sorted(tuple(sorted((labels[i], labels[j]))) for i, j in m)),
                         tuple(labels[i] for i in quad)))
    return sorted(out, key=str)


def momentum_lists(k):
    return [tuple(c) for c in itertools.combinations_with_replacement(('P0', 'P1'), k)]


def make_candidates(slots, kind, nmom, with_e, sel=None):
    """densities (functions of slot legs and ctx) for a two-leg sector.
    kind: 'VV' (two V vectors a0,a1), 'SS' (two scalars), 'PP' (two perp vectors), 'LLs' (lam, lam*), 'LsLs'
    (lam*, lam*), 'SV' (scalar slot 0, V vector slot 1)."""
    out = []
    extra = ['b'] + (['e'] if with_e else [])
    for moms in momentum_lists(nmom):
        base = list(moms) + extra
        if kind in ('VV', 'SS', 'PP', 'SV'):
            labs = base + {'VV': ['a0', 'a1'], 'SS': [], 'PP': [], 'SV': ['a1']}[kind]
            for pairs, eps in scalar_forms(labs):
                def L(s, ctx, pairs=pairs, eps=eps, kind=kind):
                    v = contraction_value(pairs, eps, s, ctx)
                    if kind == 'SS':
                        v *= s[0]['u'][0] * s[1]['u'][0]
                    elif kind == 'PP':
                        v *= np.dot(s[0]['u'][4:], s[1]['u'][4:])
                    elif kind == 'SV':
                        v *= s[0]['u'][0]
                    return v
                out.append((f'{kind}:{moms}:{pairs}:{eps}', L))
        else:
            # spinor bilinears: gamma = ordered product of distinct labels, the rest contracted
            n = len(base)
            parities = (0, 2, 4) if kind == 'LLs' else (1, 3)
            for g in parities:
                for gl in itertools.permutations(range(n), g):
                    if g >= 2 and list(gl) != sorted(gl) and False:
                        continue
                    rest = [i for i in range(n) if i not in gl]
                    for pairs, eps in scalar_forms([base[i] for i in rest]):
                        if eps is not None:
                            continue
                        glabs = [base[i] for i in gl]

                        def L(s, ctx, pairs=pairs, glabs=glabs, kind=kind):
                            v = contraction_value(pairs, None, s, ctx)
                            x = s[0]['u']
                            # build Gamma x: rightmost acts first; alternate R (Sigma->Sigma*) and RB
                            vecs = [vec_of(l, s, ctx) for l in glabs]
                            state = 'S' if kind == 'LLs' else 'Ss'
                            for vv in reversed(vecs):
                                if state == 'S':
                                    x = sum(vv[m] * (Q.R[m] @ x) for m in range(4)); state = 'Ss'
                                else:
                                    x = sum(vv[m] * (Q.RB[m] @ x) for m in range(4)); state = 'S'
                            assert state == 'S'
                            return v * (x @ Q.MP @ s[1]['u'])
                        out.append((f'{kind}:{moms}:G{tuple(glabs)}:{pairs}', L))
    return out


class Fn:
    """a two-leg functional: slots, coefficient, density L(s, ctx) for omega_0 (with_e False) or omega_mu."""

    def __init__(self, label, slots, L, with_e):
        self.label, self.slots, self.L, self.with_e = label, slots, L, with_e

    def as_canon_fn(self, b):
        if self.with_e:
            ws = [Q.dens(1, self.slots, (lambda s, mu=mu: self.L(s, dict(b=b, mu=mu)))) for mu in range(4)]
            return (None, ws)
        return (Q.dens(1, self.slots, lambda s: self.L(s, dict(b=b, mu=None))), None)


def canon_vector(fnpair, configs, Qapply=False):
    """evaluation vector: for each config (legs_at, b) the canonical data (alpha, beta_0..3) as Gaussian rationals."""
    vec = []
    for legs_at, b, fnfac in configs:
        fn = fnfac(b) if callable(fnfac) else fnfac
        pass
    return vec


def evaluate_fns(fns, configs, apply_q=False):
    """fns: list of Fn; configs: list of (legs_at, b).  Returns matrix rows (one per fn) of Gaussian rationals."""
    rows = []
    for fn in fns:
        row = []
        for legs_at, b in configs:
            cf = fn.as_canon_fn(b)
            if apply_q:
                cf = Q.Qfn(cf)
            a, bet = Q.canon(cf, legs_at)
            row += [a] + list(bet)
        rows.append(row)
    return rows


def mod_rows(rows):
    return [[Q.gq_mod(z) for z in r] for r in rows]


def independent(rows):
    """indices of a maximal independent subset (greedy, mod p)."""
    keep = []
    cur = []
    rk = 0
    for i, r in enumerate(mod_rows(rows)):
        t = cur + [r]
        nr = Q.rank_mod(t)
        if nr > rk:
            keep.append(i)
            cur = t
            rk = nr
    return keep
