"""QA (PART C referee): independent implementation of the one-loop harmonic action of N. Beisert, hep-th/0307015,
eqs. (4.27)-(4.29), on the oscillator letters (4.9): a^dagger_alpha (alpha=1,2), b^dagger_alphadot (1,2) bosonic,
c^dagger_a (a=1..4) fermionic; central charge C = n_a - n_b + n_c - 2 = 0 per site.
  H_12 |s;A> = sum_{s'} c_{n,n12,n21} delta_{C1,0} delta_{C2,0} |s';A>,
  c_{n,n12,n21} = (-1)^{1+n12 n21} Gamma((n12+n21)/2) Gamma(1+(n-n12-n21)/2) / Gamma(1+n/2),  c_{n,0,0} = h(n/2).
States: unnormalised monomials of creation operators, canonical order per site (a1,a2,b1,b2,c1,c2,c3,c4), site 1 left of
site 2.  Graded swap and graded cyclicity are produced by the fermion reordering signs.  Written without reading
theory/src.nosync/ha_*.py.  Exact arithmetic (fractions + sympy for eigenvalues)."""
import itertools, json, sys, os
from fractions import Fraction as Fr
import sympy as sp

TYPES = ['a1', 'a2', 'b1', 'b2', 'c1', 'c2', 'c3', 'c4']
FERM = [False] * 4 + [True] * 4


def harm(x):
    """h(x) for half-integer or integer x >= 0: sum_{j=1}^{x} 1/j for integer x; for half-integer x use
    h(x) = psi(x+1) + gamma; only differences/combinations with rational results are used; keep sympy exact."""
    return sp.harmonic(sp.Rational(x))


def coef(n, n12, n21):
    if n12 == 0 and n21 == 0:
        return harm(sp.Rational(n, 2))
    return (-1) ** (1 + n12 * n21) * sp.gamma(sp.Rational(n12 + n21, 2)) * sp.gamma(1 + sp.Rational(n - n12 - n21, 2)) / sp.gamma(1 + sp.Rational(n, 2))


def central(site):
    na = site[0] + site[1]; nb = site[2] + site[3]; nc = sum(site[4:])
    return na - nb + nc - 2


def ops_of(state):
    """list of (type index, site) in canonical order for a basis state (site1, site2)"""
    out = []
    for s, site in enumerate(state):
        for t in range(8):
            out += [(t, s)] * site[t]
    return out


def canon(ops):
    """sort (type, site) operators to canonical order; returns (sign, state) or (0, None) if a fermion repeats."""
    # stable sort key (site, type); count fermion inversions
    idx = list(range(len(ops)))
    key = lambda i: (ops[i][1], ops[i][0])
    order = sorted(idx, key=key)
    # sign = parity of the permutation restricted to fermions
    ferm_pos = [i for i in order if FERM[ops[i][0]]]
    inv = 0
    for x in range(len(ferm_pos)):
        for y in range(x + 1, len(ferm_pos)):
            if ferm_pos[x] > ferm_pos[y]:
                inv += 1
    sites = [[0] * 8, [0] * 8]
    for i in order:
        t, s = ops[i]
        sites[s][t] += 1
        if FERM[t] and sites[s][t] > 1:
            return 0, None
    return (-1) ** inv, (tuple(sites[0]), tuple(sites[1]))


def H12(vec):
    """vec: dict state -> coefficient (sympy)"""
    out = {}
    for st, cf in vec.items():
        ops = ops_of(st)
        n = len(ops)
        for snew in itertools.product((0, 1), repeat=n):
            n12 = sum(1 for (t, s), sn in zip(ops, snew) if s == 0 and sn == 1)
            n21 = sum(1 for (t, s), sn in zip(ops, snew) if s == 1 and sn == 0)
            newops = [(t, sn) for (t, s), sn in zip(ops, snew)]
            sg, ns = canon(newops)
            if sg == 0:
                continue
            if central(ns[0]) != 0 or central(ns[1]) != 0:
                continue
            out[ns] = out.get(ns, 0) + cf * sg * coef(n, n12, n21)
    return {k: sp.nsimplify(v) for k, v in out.items() if sp.simplify(v) != 0}


def swap(vec):
    out = {}
    for st, cf in vec.items():
        ops = [(t, 1 - s) for (t, s) in ops_of(st)]
        sg, ns = canon(ops)
        out[ns] = out.get(ns, 0) + cf * sg
    return {k: v for k, v in out.items() if v != 0}


def add(*vs):
    out = {}
    for v in vs:
        for k, c in v.items():
            out[k] = out.get(k, 0) + c
    return {k: sp.nsimplify(c) for k, c in out.items() if sp.simplify(c) != 0}


def scale(v, s):
    return {k: c * s for k, c in v.items()}


def site(**kw):
    s = [0] * 8
    for k, v in kw.items():
        s[TYPES.index(k)] = v
    return tuple(s)


def D(st, k):
    """apply (a1^dag b1^dag)^k to a site"""
    s = list(st); s[0] += k; s[2] += k
    return tuple(s)


Z = site(c3=1, c4=1)
F11 = site(a1=2)
PSI = {i: site(a1=1, **{'c%d' % i: 1}) for i in range(1, 5)}


def two(x, y, cf=1):
    """|x> (x) |y> as a basis state with the canonical sign (x ops then y ops)"""
    ops = [(t, 0) for t in range(8) for _ in range(x[t])] + [(t, 1) for t in range(8) for _ in range(y[t])]
    sg, ns = canon(ops)
    return {ns: cf * sg} if sg else {}


def cyclic_matrix(basis):
    """basis: list of two-site vectors (each Pi-invariant up to normalisation).  Returns matrix M (columns = images)
    of H = H12 + Pi H12 Pi in this basis, exact."""
    imgs = []
    for b in basis:
        h = add(H12(b), swap(H12(swap(b))))
        imgs.append(h)
    keys = sorted(set(k for v in basis + imgs for k in v))
    Bm = sp.Matrix([[v.get(k, 0) for v in basis] for k in keys])
    Im = sp.Matrix([[v.get(k, 0) for v in imgs] for k in keys])
    # solve Bm X = Im exactly
    X = (Bm.T * Bm).inv() * (Bm.T * Im)
    res = sp.simplify(Bm * X - Im)
    assert all(x == 0 for x in res), 'family not closed under H'
    return sp.simplify(X)


def spectrum(M):
    ev = []
    for e, m in M.eigenvals().items():
        ev += [sp.nsimplify(sp.simplify(e))] * m
    return sorted(ev, key=lambda x: float(x))


def family_basis(X, Y, N):
    """single trace tr(D^k X D^{N-k} Y), cyclic (graded): vectors |D^kX, D^{N-k}Y> + Pi(...) for X == Y with k <= N/2,
    or all k for X != Y (cyclicity identifies tr(D^k X D^l Y) with tr(D^l Y D^k X))."""
    out = []
    if X == Y:
        for k in range(N // 2 + 1):
            v = two(D(X, k), D(Y, N - k))
            v = add(v, swap(v))
            if v:
                out.append(v)
    else:
        for k in range(N + 1):
            v = two(D(X, k), D(Y, N - k))
            v = add(v, swap(v))
            if v:
                out.append(v)
    return out


def results():
    R = {}
    # (i) twist two scalars
    R['twist_two_Z'] = {}
    for N in range(0, 7):
        M = cyclic_matrix(family_basis(Z, Z, N))
        R['twist_two_Z'][N] = {'matrix_cols_images_basis_trDkZDN-kZ': [[str(x) for x in row] for row in M.tolist()],
                               'E_H': [str(x) for x in spectrum(M)]}
    # (ii) gluons F11
    R['gluon_F11'] = {}
    for N in range(0, 5):
        M = cyclic_matrix(family_basis(F11, F11, N))
        R['gluon_F11'][N] = {'matrix': [[str(x) for x in row] for row in M.tolist()], 'E_H': [str(x) for x in spectrum(M)]}
    # (iii) one gaugino component psi = a1 c4
    R['fermion_psi'] = {}
    for N in range(0, 6):
        B = family_basis(PSI[4], PSI[4], N)
        if not B:
            R['fermion_psi'][N] = {'dim': 0}
            continue
        M = cyclic_matrix(B)
        R['fermion_psi'][N] = {'matrix': [[str(x) for x in row] for row in M.tolist()], 'E_H': [str(x) for x in spectrum(M)]}
    # (v) (F11, psi)
    R['F11_psi'] = {}
    for N in range(0, 5):
        M = cyclic_matrix(family_basis(F11, PSI[4], N))
        R['F11_psi'][N] = {'matrix': [[str(x) for x in row] for row in M.tolist()], 'E_H': [str(x) for x in spectrum(M)]}
    # (iv) psi3 psi4 with Z F11: basis tr(D^k psi3 D^l psi4) and tr(D^k Z D^l F11), closed check
    R['psi3psi4_ZF11'] = {}
    for N in range(0, 5):
        B = family_basis(PSI[3], PSI[4], N) + family_basis(Z, F11, N)
        M = cyclic_matrix(B)
        R['psi3psi4_ZF11'][N] = {'dim': len(B), 'E_H': [str(x) for x in spectrum(M)]}
    return R


if __name__ == '__main__':
    R = results()
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'qa_harmonic_two_site.json')
    json.dump(R, open(out, 'w'), indent=1)
    for fam, d in R.items():
        for N, v in d.items():
            print(fam, N, v.get('E_H', v))
