"""Engine A: explicit kernel tables.  For a pair type (ta, tb) and an ordered output word with a block structure,
the singular coefficients of K(X_a, X_b) as exact polynomials in the output momenta, for every component of the two
letters and of the outputs (component bases: ea_core.BASIS), by exact homogeneous interpolation."""
import itertools, math, json, time
from fractions import Fraction
import numpy as np
import ea_core as E, ea_kern as K

LT = {'AV': 'AV', 'AW': 'AW', 'lamP': 'lamP', 'lamM': 'lamM', 'c': 'c', 'AsV': 'AsV', 'AsW': 'AsW',
      'lamsP': 'lamsP', 'lamsM': 'lamsM', 'cs': 'cs'}


def stripped_pol(t, i):
    """i-th basis polarisation of output type t as a stripped slot vector"""
    b = E.BASIS[t][i]
    sl = E.SLOT[t]
    v = np.zeros(E.SZ[sl], dtype=complex)
    if t in ('AV', 'AW', 'AsV', 'AsW', 'c', 'cs'):
        v[:] = b
    else:
        v[:] = b
    return v


def letter_vecs(t):
    """basis vectors of a letter type inside its kernel slot rows (ea_kern.SLOTROWS)"""
    if t in ('AV', 'AW', 'AsV', 'c'):
        n = len(E.BASIS[t])
        return [np.eye(n, dtype=int)[i] for i in range(n)]
    return [np.asarray(b) for b in E.BASIS[t]]


def block_structures(word):
    """compositions of the word (list of types) into cubic singles and quartic consecutive (A,A) pairs"""
    n = len(word)
    def rec(i):
        if i == n:
            yield []; return
        for tail in rec(i + 1):
            yield [(i,)] + tail
        if i + 1 < n and E.SLOT[word[i]] == 'A' and E.SLOT[word[i + 1]] == 'A':
            for tail in rec(i + 2):
                yield [(i, i + 1)] + tail
    return list(rec(0))


def words_of(legset):
    types = legset.split() if legset not in ('', '(none)') else []
    return sorted(set(itertools.permutations(types)))


def hom_points(nvar, D):
    pts = []
    for c in itertools.combinations_with_replacement(range(nvar), D):
        x = [0] * nvar
        for i in c:
            x[i] += 1
        pts.append(tuple(x))
    return pts


def mono_list(nvar, D):
    return hom_points(nvar, D)


_INV = {}
HOMOG_FAIL = []


def interp_inverse(nvar, D):
    """inverse of the evaluation matrix of degree-D monomials at the simplex-lattice points"""
    key = (nvar, D)
    if key in _INV:
        return _INV[key]
    pts = hom_points(nvar, D); mons = mono_list(nvar, D)
    n = len(pts)
    M = [[Fraction(math.prod(p[i] ** m[i] for i in range(nvar))) for m in mons] for p in pts]
    # Gauss-Jordan
    A = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        piv = A[c][c]
        A[c] = [x / piv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    inv = [row[n:] for row in A]
    _INV[key] = (inv, pts, mons)
    return _INV[key]


def word_table(ta, tb, word, struct, kinds=('log', 'eps'), D_by_kind=None):
    """returns {kind: {(ia, ib, (io...)): {monomial: (re, im)}}} ; kind 'log' or m (eps^{-m})"""
    n = len(word)
    nvar = 4 * n
    W = E.WEYL[ta] + E.WEYL[tb]
    wout = sum((E.WEYL[t] for t in word), Fraction(0))
    Dlog = W - wout
    out = {}
    targets = []
    if Dlog >= 0 and Dlog.denominator == 1:
        targets.append(('log', int(Dlog)))
        m = 1
        while Dlog - 2 * m >= 0:
            targets.append((m, int(Dlog - 2 * m))); m += 1
    if not targets:
        return out
    va, vb = letter_vecs(ta), letter_vecs(tb)
    comps = [range(len(E.BASIS[t])) for t in word]
    Ds = sorted(set(d for _, d in targets))
    for D in Ds:
        inv, pts, mons = interp_inverse(nvar, D) if nvar else ([[Fraction(1)]], [()], [()])
        kinds_here = [k for k, d in targets if d == D]
        for oc in itertools.product(*comps):
            vals = {k: [] for k in kinds_here}
            for p in pts:
                moms = [np.array(p[4 * j:4 * j + 4], dtype=int) for j in range(n)]
                blocks = []
                for b in struct:
                    if len(b) == 1:
                        j = b[0]
                        blocks.append(('cub', word[j], stripped_pol(word[j], oc[j]), moms[j]))
                    else:
                        j, jj = b
                        blocks.append(('qua', (stripped_pol(word[j], oc[j]), moms[j]),
                                       (stripped_pol(word[jj], oc[jj]), moms[jj])))
                mmax = max([k for k in kinds_here if k != 'log'] + [0])
                r = K.kernel(blocks, ta, tb, kinds=('log', 'eps') if mmax else ('log',), mmax=mmax or None) if n else \
                    K.kernel([], ta, tb, kinds=('log', 'eps'), mmax=mmax or None)
                for k in kinds_here:
                    vals[k].append(r.get(k) if r else None)
            # homogeneity / degree check (gate G7): one extra generic point, compared after interpolation
            rngp = np.random.default_rng(len(pts) + 7 * D + sum(oc))
            pchk = tuple(int(x) for x in rngp.integers(-3, 4, nvar))
            moms = [np.array(pchk[4 * j:4 * j + 4], dtype=int) for j in range(n)]
            blocks = []
            for b in struct:
                if len(b) == 1:
                    j = b[0]
                    blocks.append(('cub', word[j], stripped_pol(word[j], oc[j]), moms[j]))
                else:
                    j, jj = b
                    blocks.append(('qua', (stripped_pol(word[j], oc[j]), moms[j]), (stripped_pol(word[jj], oc[jj]), moms[jj])))
            mmax = max([k for k in kinds_here if k != 'log'] + [0])
            rchk = K.kernel(blocks, ta, tb, kinds=('log', 'eps') if mmax else ('log',), mmax=mmax or None)
            for k in kinds_here:
                # letter projections
                for ia, za in enumerate(va):
                    for ib, zb in enumerate(vb):
                        ys = []
                        for v in vals[k]:
                            ys.append(K.contract(v, za, zb) if v is not None else E.Q0)
                        coeffs = {}
                        for mi, mon in enumerate(mons):
                            re = sum(inv[mi][pi] * ys[pi][0] for pi in range(len(pts)))
                            im = sum(inv[mi][pi] * ys[pi][1] for pi in range(len(pts)))
                            if re or im:
                                coeffs[mon] = (re, im)
                        # check at the extra point
                        vchk = K.contract(rchk[k], za, zb) if (rchk and rchk.get(k) is not None) else E.Q0
                        pr, pi_ = Fraction(0), Fraction(0)
                        for mon, (cr, ci) in coeffs.items():
                            w = math.prod(Fraction(pchk[i]) ** mon[i] for i in range(nvar))
                            pr += cr * w; pi_ += ci * w
                        if (pr, pi_) != (vchk[0], vchk[1]):
                            HOMOG_FAIL.append((ta, tb, word, k, ia, ib, oc))
                        if coeffs:
                            out.setdefault(k, {})[(ia, ib, tuple(oc))] = coeffs
    return out


def poly_str(coeffs, n):
    """{monomial over r_{j,mu}: (re,im)} -> string in r1_0..rn_3"""
    terms = []
    for mon, z in sorted(coeffs.items()):
        f = []
        for idx, e in enumerate(mon):
            if e:
                j, mu = divmod(idx, 4)
                f.append('r%d_%d' % (j + 1, mu) + ('^%d' % e if e > 1 else ''))
        terms.append(E.qstr(z) + ('*' + '*'.join(f) if f else ''))
    return ' + '.join(terms) if terms else '0'
