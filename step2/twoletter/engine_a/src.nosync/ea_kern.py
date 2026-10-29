"""Engine A: two-letter kernels from colour-stripped chains (see ea_core).  kernel() returns, for a chain with given
output blocks, the singular coefficients of the letter-component matrix K[a, b] (a: components of the slot of letter
a, b: of letter b), in units hbar/(16 pi^2), vertex scales divided out; colour element sum_gamma t^gamma (x) W t_gamma,
W = ad_{Y_k}...ad_{Y_1} (Y_1 nearest to a).  For M = int omega (zeta_a, X_a)(zeta_b, X_b) the pair contributes
zeta_a^T K zeta_b (letter factors for derivatives are applied by kernel_letters)."""
from fractions import Fraction
import numpy as np
import ea_core as E

SLOTROWS = {
    'AV': [E.POS[E.OFS['A'] + M] for M in range(4)],
    'AW': [E.POS[E.OFS['A'] + M] for M in range(4, 10)],
    'A': [E.POS[E.OFS['A'] + M] for M in range(10)],
    'c': [E.POS[E.OFS['c']]],
    'AsV': [E.POS[E.OFS['As'] + M] for M in range(4)],
    'l': [E.POS[E.OFS['l'] + a] for a in range(16)],
}
LSLOT = {'AV': 'AV', 'AW': 'AW', 'c': 'c', 'AsV': 'AsV', 'lamP': 'l', 'lamM': 'l'}


def make_blocks(outs):
    """outs: list of blocks; each block = ('cub', ltype, pol, mom) or ('qua', (pol_first, mom1), (pol_second, mom2))
    where in a quartic block the first output is nearer to letter a.  Returns (affine blocks, momenta, scale)."""
    blks, Ps, scale = [], [], 1
    for b in outs:
        if b[0] == 'cub':
            _, lt, pol, mom = b
            blks.append(E.cubic_block(lt, pol, mom)); Ps.append(np.asarray(mom)); scale *= 2
        else:
            _, (p1, m1), (p2, m2) = b
            # word ad_{second} ad_{first}: quartic_block(polY, polY') = word ad_Y ad_Y'
            blks.append(E.quartic_block(p2, p1)); Ps.append(np.asarray(m1) + np.asarray(m2)); scale *= 8
    return blks, Ps, scale


def chain_poly(outs, ta, tb, line_seq=None):
    blks, Ps, scale = make_blocks(outs)
    T, qs = E.chain_T(blks, Ps, line_seq)
    ra, rb = SLOTROWS[LSLOT.get(ta, ta)], SLOTROWS[LSLOT.get(tb, tb)]
    K = {e: M[np.ix_(ra, rb)] for e, M in T.items()}
    K = {e: M for e, M in K.items() if np.any(M)}
    return K, qs, scale


def apply_letters(K, fa, fb, P0):
    """multiply the polynomial matrix K(k) by letter factors: letter a has incoming momentum k - P0, letter b -k;
    fa, fb: {exp: Gaussian int} polynomials in the letter momentum (i p)^alpha already included by the caller.
    Returns polynomial in k."""
    def shift(f, sgn, P):
        # f(sgn*k + P) as polynomial in k
        out = {(0, 0, 0, 0): 1}
        res = {}
        for e, c in f.items():
            term = {(0, 0, 0, 0): c}
            for mu in range(4):
                for _ in range(e[mu]):
                    lin = {}
                    lin[tuple(np.eye(4, dtype=int)[mu])] = sgn
                    if P[mu]:
                        lin[(0, 0, 0, 0)] = P[mu]
                    new = {}
                    for e1, c1 in term.items():
                        for e2, c2 in lin.items():
                            ee = tuple(x + y for x, y in zip(e1, e2))
                            new[ee] = new.get(ee, 0) + c1 * c2
                    term = new
            for ee, cc in term.items():
                res[ee] = res.get(ee, 0) + cc
        return {e: c for e, c in res.items() if c != 0}
    ga = shift(fa, 1, -np.asarray(P0))
    gb = shift(fb, -1, np.zeros(4, dtype=int))
    out = {}
    for e1, M in K.items():
        for e2, c2 in ga.items():
            for e3, c3 in gb.items():
                ee = tuple(x + y + z for x, y, z in zip(e1, e2, e3))
                out[ee] = out[ee] + M * (c2 * c3) if ee in out else M * (c2 * c3)
    return out


def kernel(outs, ta, tb, kinds=('log', 'eps'), fa=None, fb=None, line_seq=None, mmax=None, intrinsic=True):
    """intrinsic=True: the two-letter kernel K(X_a, X_b) of Def def:tl-dressed.  The str formula with the operator
    Hessian gives the coefficient of the delta-vertex Hessian; d_j d_l (delta M) = (-1)^{|delta|(|j|+|l|)} delta d_j d_l M
    with |delta| = |M| = |X_a| + |X_b| for the pair functional, so the intrinsic kernel is (-1)^{|X_a|+|X_b|} times the
    str value (fixed once here; checked by the Ward identities W1, ea_ward.py)."""
    K, qs, scale = chain_poly(outs, ta, tb, line_seq)
    if intrinsic and (E.PARITY_T[ta] + E.PARITY_T[tb]) % 2 == 1:
        scale = -scale
    if not K:
        return None
    if fa is not None or fb is not None:
        P0 = -sum((np.asarray(b[3]) if b[0] == 'cub' else np.asarray(b[1][1]) + np.asarray(b[2][1]) for b in outs),
                  np.zeros(4, dtype=int))
        K = apply_letters(K, fa or {(0, 0, 0, 0): 1}, fb or {(0, 0, 0, 0): 1}, P0)
    S = E.Simplex(qs)
    res = S.singular(K, mmax=mmax, kinds=kinds)
    out = {}
    for key, v in res.items():
        if v is None:
            continue
        re, im = v
        out[key] = (re / scale, im / scale)
    return out


def contract(val, za, zb):
    """zeta_a^T K zeta_b with exact Gaussian-integer zeta vectors (lists of complex) -> (Fraction, Fraction)"""
    re, im = val
    tr, ti = Fraction(0), Fraction(0)
    for i, a in enumerate(za):
        if a == 0:
            continue
        for j, b in enumerate(zb):
            if b == 0:
                continue
            c = complex(a) * complex(b)
            cr, ci = int(round(c.real)), int(round(c.imag))
            tr += cr * re[i, j] - ci * im[i, j]
            ti += cr * im[i, j] + ci * re[i, j]
    return (tr, ti)
