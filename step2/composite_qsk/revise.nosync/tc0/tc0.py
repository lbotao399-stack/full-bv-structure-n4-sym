"""TC0: Q-cohomology of quadratic local functionals of degree -1, Weyl weight 4 (position at -1), degree <= 1 in the
position, invariant under Spin(4) x Spin(6) (colour: one invariant pairing, suppressed), with Q of Prop c1-Q:
  Q(A_V) = -dc,  Q(A*) = -1/2 d*dA,  Q(c*) = d*(A*|_V),  Q(lambda*) = dslash lambda,  Q(c) = Q(A_perp) = Q(lambda) = 0
(coordinate form of the operator of the table; Q odd derivation, Q(XY) = Q(X)Y + (-1)^{|X|} X Q(Y)).
Total parities |c| = |A*| = |lambda| = 1, |A| = |lambda*| = |c*| = 0.  Spinors: the engine-A Clifford blocks
(cliff.py, vendored read-only from step2/twoletter/engine_a): rho(v)|_Sigma = sb[M] v^M, rho(v)|_Sigma* = s[M] v^M,
pairing u^T w (u in Sigma, w in Sigma*).  A_perp: one representative component (SO(6)-invariant structures are delta_IJ).

A bilinear functional of two plane waves X (momentum p1), Y (momentum p2), int d^4x (k0 + x^mu k_mu)(...) X Y, is zero
iff its graded-symmetrised kernel has the normal form NF = (k_mu(p,-p), k0(p,-p) + i d_{p2}.k(p,-p)) = 0.
Structures (all invariant ones, by Schur): degree -1: (A*_V, D A_V): 1, x.d, x_a d_b, x_b d_a, eps_{ab mu nu} x^mu d^nu;
(A*_perp, D A_perp): 1, x.d; (lambda*, D lambda): 1, G5, x.d, G5 x.d, rho(x)rho(d), G5 rho(x)rho(d); (c*, D c): 1, x.d.
Degree -2: x^mu (A_mu c*), (lambda*, rho(x) lambda*), (lambda*, G5 rho(x) lambda*).
Parity: the reflection of one direction of V and one of W_R^perp (ch09, proof of the (c,lambda,lambda) basis) maps G5 to
-G5 and eps to -eps; parity-even = no G5, no eps."""
import os, sys, json, itertools
import sympy as sp
import numpy as np

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'step2', 'twoletter', 'engine_a', 'src.nosync', 'vendor'))
import cliff  # noqa: E402

I = sp.I
def tosym(M):
    return sp.Matrix(M.shape[0], M.shape[1], lambda i, j: sp.nsimplify(complex(M[i, j]).real) + I * sp.nsimplify(complex(M[i, j]).imag))
s = [tosym(m) for m in cliff.s]
sb = [tosym(m) for m in cliff.sb]
G5 = s[0] * sb[1] * s[2] * sb[3]
assert G5 * G5 == sp.eye(16)

p1 = sp.symbols('a0:4'); p2 = sp.symbols('b0:4'); P = sp.symbols('p0:4')
DIM = {'c': 1, 'AV': 4, 'AW': 1, 'l': 16, 'AsV': 4, 'AsW': 1, 'ls': 16, 'cs': 1}
PAR = {'c': 1, 'AV': 0, 'AW': 0, 'l': 1, 'AsV': 1, 'AsW': 1, 'ls': 0, 'cs': 0}
ORDER = ['c', 'AV', 'AW', 'l', 'AsV', 'AsW', 'ls', 'cs']

# Q in coordinate form: Q(X) = sum_Z M[X][Z](p) Z  (p the common momentum), matrices rows: comps of X, cols: comps of Z
def Qmap(X, p):
    if X == 'AV':
        return {'c': sp.Matrix(4, 1, lambda a, _: -I * p[a])}
    if X == 'AsV':
        pp = sum(q * q for q in p)
        return {'AV': sp.Matrix(4, 4, lambda a, b: -sp.Rational(1, 2) * ((pp if a == b else 0) - p[a] * p[b]))}
    if X == 'AsW':
        pp = sum(q * q for q in p)
        return {'AW': sp.Matrix([[-sp.Rational(1, 2) * pp]])}
    if X == 'cs':
        return {'AsV': sp.Matrix(1, 4, lambda _, b: -I * p[b])}
    if X == 'ls':
        return {'l': I * sum((sb[m] * p[m] for m in range(4)), sp.zeros(16, 16))}
    return {}


class F:
    """bilinear functional: dict (X, Y) -> {'0': K0(p1,p2), mu: Kmu(p1,p2)} (matrices rows X comps, cols Y comps),
    meaning int (k0 + x^mu k_mu) X Y with X at momentum p1, Y at p2."""
    def __init__(self):
        self.t = {}

    def add(self, X, Y, key, M):
        d = self.t.setdefault((X, Y), {})
        d[key] = d.get(key, sp.zeros(DIM[X], DIM[Y])) + M

    def scaled(self, c):
        G = F()
        for k, d in self.t.items():
            for key, M in d.items():
                G.add(k[0], k[1], key, c * M)
        return G

    def __add__(self, o):
        G = F()
        for src in (self, o):
            for k, d in src.t.items():
                for key, M in d.items():
                    G.add(k[0], k[1], key, M)
        return G


def Q(Fn):
    out = F()
    for (X, Y), d in Fn.t.items():
        for key, M in d.items():
            for Z, MX in Qmap(X, p1).items():            # Q(X) Y : X -> MX Z  =>  Z^T MX^T M Y
                out.add(Z, Y, key, MX.T * M)
            for Z, MY in Qmap(Y, p2).items():            # (-1)^{|X|} X Q(Y): X^T M MY Z
                out.add(X, Z, key, (-1) ** PAR[X] * M * MY)
    return out


def canon(Fn):
    """reorder to ORDER (graded swap with p1 <-> p2) and graded-symmetrise equal types"""
    out = F()
    sw = {**{a: b for a, b in zip(p1, p2)}, **{b: a for a, b in zip(p1, p2)}}
    for (X, Y), d in Fn.t.items():
        for key, M in d.items():
            if ORDER.index(X) <= ORDER.index(Y):
                out.add(X, Y, key, M)
            if ORDER.index(X) >= ORDER.index(Y):
                Ms = (-1) ** (PAR[X] * PAR[Y]) * M.T.subs(sw, simultaneous=True)
                if X == Y:
                    out.add(X, Y, key, Ms)               # symmetrised: K + (-1)^{|X|} K^T(p2,p1)
                else:
                    out.add(Y, X, key, Ms)
    return out


def NF(Fn):
    """normal form: list of scalar polynomials in P (the momentum of the first leg), keyed"""
    C = canon(Fn)
    res = {}
    sub = {**{a: P[m] for m, a in enumerate(p1)}, **{b: -P[m] for m, b in enumerate(p2)}}
    for (X, Y), d in C.t.items():
        k0 = d.get('0', sp.zeros(DIM[X], DIM[Y]))
        km = [d.get(m, sp.zeros(DIM[X], DIM[Y])) for m in range(4)]
        n0 = k0.subs(sub, simultaneous=True) + I * sum((km[m].diff(p2[m]) for m in range(4)), sp.zeros(DIM[X], DIM[Y])).subs(sub, simultaneous=True)
        for i in range(DIM[X]):
            for j in range(DIM[Y]):
                res[(X, Y, 'n0', i, j)] = sp.expand(n0[i, j])
                for m in range(4):
                    res[(X, Y, m, i, j)] = sp.expand(km[m].subs(sub, simultaneous=True)[i, j])
    return res


def coeff_rows(nfs):
    """nfs: list of NF dicts -> matrix (rows: (key, monomial), cols: functionals)"""
    keys = set()
    polys = []
    for nf in nfs:
        dd = {}
        for k, v in nf.items():
            if v == 0:
                continue
            pol = sp.Poly(v, *P)
            for mon, c in pol.terms():
                dd[(k, mon)] = c
                keys.add((k, mon))
        polys.append(dd)
    keys = sorted(keys, key=str)
    return sp.Matrix(len(keys), len(nfs), lambda r, c: polys[c].get(keys[r], 0)) if keys else sp.zeros(1, len(nfs))


# ------------------------------------------------------------------------- structures
# position x^mu <-> key mu; a derivative on X (momentum p1) is i p1, on Y is i p2.
def deg1_structures():
    S = {}
    eps = lambda a, b, c, d: sp.LeviCivita(a, b, c, d)
    # (A*_V, D A_V): rows A*_V comps, cols A_V comps; functional int (A*_V, D A_V) = sum_ab A*_a D_ab A_b
    def vec(name, f0=None, fx=None):
        G = F()
        if f0 is not None:
            G.add('AsV', 'AV', '0', f0)
        if fx is not None:
            for m in range(4):
                G.add('AsV', 'AV', m, fx(m))
        S[name] = G
    vec('V:1', f0=sp.eye(4))
    vec('V:x.d', fx=lambda m: I * p2[m] * sp.eye(4))
    vec('V:x_a d_b', fx=lambda m: sp.Matrix(4, 4, lambda a, b: I * p2[b] * (1 if a == m else 0)))
    vec('V:x_b d_a', fx=lambda m: sp.Matrix(4, 4, lambda a, b: I * p2[a] * (1 if b == m else 0)))
    vec('V:eps x d', fx=lambda m: sp.Matrix(4, 4, lambda a, b: sum(eps(a, b, m, nn) * I * p2[nn] for nn in range(4))))
    for nm, f0, fx in (('W:1', sp.eye(1), None), ('W:x.d', None, lambda m: I * p2[m] * sp.eye(1))):
        G = F()
        if f0 is not None:
            G.add('AsW', 'AW', '0', f0)
        if fx is not None:
            for m in range(4):
                G.add('AsW', 'AW', m, fx(m))
        S[nm] = G
    # (lambda*, D lambda) = lambda*^T D lambda  (pairing u^T w, lambda* in Sigma*, D lambda in Sigma)
    rr = lambda m, nn: s[m] * sb[nn]              # rho(e_m) rho(e_n) on Sigma
    lam = {
        'L:1': (sp.eye(16), None), 'L:G5': (G5, None),
        'L:x.d': (None, lambda m: I * p2[m] * sp.eye(16)), 'L:G5 x.d': (None, lambda m: I * p2[m] * G5),
        'L:rho(x)rho(d)': (None, lambda m: sum((rr(m, nn) * I * p2[nn] for nn in range(4)), sp.zeros(16, 16))),
        'L:G5 rho(x)rho(d)': (None, lambda m: G5 * sum((rr(m, nn) * I * p2[nn] for nn in range(4)), sp.zeros(16, 16))),
    }
    for nm, (f0, fx) in lam.items():
        G = F()
        if f0 is not None:
            G.add('ls', 'l', '0', f0)
        if fx is not None:
            for m in range(4):
                G.add('ls', 'l', m, fx(m))
        S[nm] = G
    for nm, f0, fx in (('c:1', sp.eye(1), None), ('c:x.d', None, lambda m: I * p2[m] * sp.eye(1))):
        G = F()
        if f0 is not None:
            G.add('cs', 'c', '0', f0)
        if fx is not None:
            for m in range(4):
                G.add('cs', 'c', m, fx(m))
        S[nm] = G
    return S


def deg2_structures():
    S = {}
    G = F()
    for m in range(4):
        G.add('AV', 'cs', m, sp.Matrix(4, 1, lambda a, _: 1 if a == m else 0))
    S['x^mu (A_mu, c*)'] = G
    for nm, M in (('(lambda*, rho(x) lambda*)', lambda m: s[m]), ('(lambda*, G5 rho(x) lambda*)', lambda m: G5 * s[m])):
        G = F()
        for m in range(4):
            G.add('ls', 'ls', m, M(m))
        S[nm] = G
    return S


PARITY_ODD = {'V:eps x d', 'L:G5', 'L:G5 x.d', 'L:G5 rho(x)rho(d)', '(lambda*, G5 rho(x) lambda*)'}


def analyse(names1, S1, names2, S2):
    N1 = coeff_rows([NF(S1[n]) for n in names1])                    # degree -1 normal forms
    N0 = coeff_rows([NF(Q(S1[n])) for n in names1])                 # their Q-images
    rank1 = N1.rank()
    ker0 = N0.nullspace()
    closed = sp.Matrix.hstack(*[N1 * v for v in ker0]) if ker0 else sp.zeros(N1.rows, 0)
    dim_closed = closed.rank() if ker0 else 0
    # exact: Q of degree -2 structures, normal forms on the degree -1 keys
    E = [NF(Q(S2[n])) for n in names2]
    allk = coeff_rows([NF(S1[n]) for n in names1] + E)
    if not names2:
        return dict(structures=names1, rank_structures=int(rank1), dim_closed=int(dim_closed), dim_exact=0, exact_in_span=True,
                    exact_in_closed=True, dim_cohomology=int(dim_closed), closed_nullspace_basis=[[str(x) for x in v] for v in ker0],
                    exact_in_structure_coordinates=[], degree_minus2=[])
    nE = len(names2)
    exactM = allk[:, len(names1):]
    span1 = allk[:, :len(names1)]
    dim_exact = exactM.rank()
    # exact inside the span of the structures and inside the closed space
    in_span = sp.Matrix.hstack(span1, exactM).rank() == span1.rank()
    closedA = sp.Matrix.hstack(*[allk[:, :len(names1)] * v for v in ker0]) if ker0 else sp.zeros(allk.rows, 0)
    in_closed = sp.Matrix.hstack(closedA, exactM).rank() == (closedA.rank() if ker0 else 0)
    # express closed vectors and exact vectors in the structure basis (reduced)
    closed_vecs = [[str(x) for x in v] for v in ker0]
    ex_coords = []
    for j in range(nE):
        sol = sp.Matrix.hstack(span1, exactM[:, j]).nullspace()
        ex_coords.append([str(x) for x in (sol[0] / -sol[0][-1])[:-1]] if sol else None)
    return dict(structures=names1, rank_structures=int(rank1), dim_closed=int(dim_closed), dim_exact=int(dim_exact),
                exact_in_span=bool(in_span), exact_in_closed=bool(in_closed), dim_cohomology=int(dim_closed - dim_exact),
                closed_nullspace_basis=closed_vecs, exact_in_structure_coordinates=ex_coords, degree_minus2=names2)


def check_known(S1):
    """S_perp, S_lambda, S_V of Lemma c9-SDelta (2E = 2 x.d componentwise) and the kappa directions"""
    out = {}
    Sp = S1['W:x.d'].scaled(2) + S1['W:1'].scaled(2)
    Sl = (S1['L:x.d'].scaled(2) + S1['L:1'].scaled(3)).scaled(-1)
    SV = S1['V:x.d'].scaled(2) + S1['V:1'].scaled(2) + S1['c:x.d'].scaled(-2)
    for nm, G in (('S_perp', Sp), ('S_lambda', Sl), ('S_V', SV), ('S^(1)(Delta)', Sp + Sl + SV)):
        out[nm + ' Q-closed'] = all(v == 0 for v in NF(Q(G)).values())
    for nm, G in (('Q int(iota_x A_V, c*)', None), ('Q int(rho(x) lambda*, lambda*)', None)):
        pass
    return out


if __name__ == '__main__':
    S1 = deg1_structures(); S2 = deg2_structures()
    res = {'known': check_known(S1)}
    print(res['known'], flush=True)
    secs = {'perp': (['W:1', 'W:x.d'], []),
            'lambda': ([n for n in S1 if n.startswith('L:')], [n for n in S2 if 'lambda' in n]),
            'V and c': ([n for n in S1 if n[:2] in ('V:', 'c:')], [n for n in S2 if 'A_mu' in n])}
    for sector, (a1, a2) in secs.items():
        for par in ('parity_even', 'all'):
            b1 = [n for n in a1 if par == 'all' or n not in PARITY_ODD]
            b2 = [n for n in a2 if par == 'all' or n not in PARITY_ODD]
            r = analyse(b1, S1, b2, S2)
            res[f'{sector} / {par}'] = r
            print(sector, par, {k: r[k] for k in ('rank_structures', 'dim_closed', 'dim_exact', 'dim_cohomology', 'exact_in_closed')}, flush=True)
            print('   closed basis', r['closed_nullspace_basis'], flush=True)
            print('   exact coords', r['exact_in_structure_coordinates'], flush=True)
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'tc0.json'), 'w'), indent=1)
