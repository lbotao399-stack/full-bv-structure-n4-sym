"""Step 2, engine A: one-loop graphs with exactly one special edge K' = K_eps - K_0 on the cycle (type I-a).

Formula (validated in 0 dimensions, s2_verify0d.py): the 1PI part of hbar Delta_{K'} W_0 is
    hbar * (1/2) * str[ K' S H (P S H)^{n-1} ]   (H = Hessian of I at the background, S = parity operator),
expanded over ordered sequences of insertions (one external leg on a cubic vertex, or two external A legs on the
quartic vertex) with the special line in position 1.  Line i enters insertion i; k_1 = k (loop momentum),
k_{i+1} = k_i - P_{s_i}; numerator
    N(k) = -k_1^2 str[ Khat_{L1} S h_1(-k_1,k_2) Pnum_{L2}(k_2) S h_2(-k_2,k_3) ... Pnum_{Ln}(k_n) S h_n(-k_n,k_1) ],
since K'(k) = (e^{-eps k^2} - 1) Khat = -int_0^eps dtau k^2 e^{-tau k^2} Khat; P edges int_eps^infty dt e^{-t k^2} Pnum.
Gaussian loop integral: int d^4k/(2pi)^4 N e^{-sum t_i k_i^2} = (16pi^2)^{-1} T^{-2} sum_j (Delta^j N)(-K(u))/(4^j j! T^j)
e^{-T q(u)}, K(u) = sum u_i q_i, q(u) = sum u_i q_i^2 - K(u)^2, q_i = k_i - k, u_i = t_i/T (Delta = sum d^2/dk_mu^2).
Scaling t_i = eps tau_i (tau_1 in [0,1], tau_{i>=2} in [1,infty)) and expanding e^{-eps T q}:
    value = (16pi^2)^{-1} sum_{j,i} eps^{n-2-j+i} /(4^j j! i!) int dtau T^{-2-j+i} F_j(u) (-q(u))^i,
F_j(u) = (Delta^j N)(-K(u)).  The eps^e coefficient (e = -1, 0) collects j - i = n - 2 - e; each integral converges
(large-tau exponent -2-j+i+(n-1) = e-1 < 0), and the remainder is O(eps log eps): no log eps term (s2_notes.txt).
N(k) is reconstructed EXACTLY as a polynomial in k of total degree <= d (d = superficial numerator degree) from its
values on the simplex lattice {alpha in Z_{>=0}^4, |alpha| <= d} (multivariate Newton forward differences), checked at
extra points."""
import itertools, math
from fractions import Fraction
import numpy as np
import s2_core as C
from r2_core import Q2
import s2_integrals as SI

E4 = [np.eye(4, dtype=int)[m] for m in range(4)]

# ---------------------------------------------------------------- polynomials with Q2 coefficients
class Poly:
    __slots__ = ['t', 'nv']
    def __init__(self, nv, t=None):
        self.nv = nv; self.t = t if t is not None else {}
    @staticmethod
    def const(nv, c):
        return Poly(nv, {(0,) * nv: c if isinstance(c, Q2) else Q2(c)})
    @staticmethod
    def var(nv, i, c=1):
        e = [0] * nv; e[i] = 1
        return Poly(nv, {tuple(e): c if isinstance(c, Q2) else Q2(c)})
    def __add__(self, o):
        r = dict(self.t)
        for k, v in o.t.items():
            r[k] = r[k] + v if k in r else v
        return Poly(self.nv, {k: v for k, v in r.items() if not v.iszero()})
    def __sub__(self, o):
        return self + o.scale(Q2(-1))
    def scale(self, c):
        c = c if isinstance(c, Q2) else Q2(c)
        return Poly(self.nv, {k: v * c for k, v in self.t.items() if not (v * c).iszero()})
    def __mul__(self, o):
        r = {}
        for k1, v1 in self.t.items():
            for k2, v2 in o.t.items():
                k = tuple(a + b for a, b in zip(k1, k2))
                r[k] = r[k] + v1 * v2 if k in r else v1 * v2
        return Poly(self.nv, {k: v for k, v in r.items() if not v.iszero()})
    def deriv(self, i):
        r = {}
        for k, v in self.t.items():
            if k[i] == 0:
                continue
            kk = list(k); kk[i] -= 1
            r[tuple(kk)] = v * k[i]
        return Poly(self.nv, r)
    def lap(self):
        out = Poly(self.nv)
        for i in range(self.nv):
            out = out + self.deriv(i).deriv(i)
        return out
    def degree(self):
        return max((sum(k) for k in self.t), default=-1)
    def subs_linear(self, forms, nv_new):
        """substitute x_i = forms[i] (Poly in nv_new variables)"""
        out = Poly(nv_new)
        powcache = {}
        def pw(i, e):
            if (i, e) not in powcache:
                p = Poly.const(nv_new, 1)
                for _ in range(e):
                    p = p * forms[i]
                powcache[(i, e)] = p
            return powcache[(i, e)]
        for k, v in self.t.items():
            term = Poly.const(nv_new, v)
            for i, e in enumerate(k):
                if e:
                    term = term * pw(i, e)
            out = out + term
        return out
    def eval(self, x):
        tot = Q2(0)
        for k, v in self.t.items():
            m = Fraction(1)
            for xi, e in zip(x, k):
                m *= Fraction(xi) ** e
            tot = tot + v * m
        return tot


def binom_poly(m):
    """coefficients of binom(x, m) as a polynomial in x (list of Fractions, index = power)"""
    c = [Fraction(1)]
    for r in range(m):
        new = [Fraction(0)] * (len(c) + 1)
        for i, a in enumerate(c):
            new[i + 1] += a
            new[i] -= a * r
        c = new
    return [a / math.factorial(m) for a in c]


_lat = {}
def lattice(d):
    if d not in _lat:
        _lat[d] = [a for a in itertools.product(range(d + 1), repeat=4) if sum(a) <= d]
    return _lat[d]


def interpolate(f, d):
    """exact polynomial (Poly in 4 vars) of total degree <= d agreeing with f on the simplex lattice"""
    vals = {a: Q2.of(f(np.array(a))) for a in lattice(d)}
    poly = Poly(4)
    bp = [binom_poly(m) for m in range(d + 1)]
    for al in lattice(d):
        c = Q2(0)
        for be in itertools.product(*[range(x + 1) for x in al]):
            w = 1
            for x, y in zip(al, be):
                w *= math.comb(x, y)
            c = c + vals[be] * ((-1) ** (sum(al) - sum(be)) * w)
        if c.iszero():
            continue
        term = Poly.const(4, c)
        for mu in range(4):
            if al[mu]:
                term = term * Poly(4, {tuple(i if nu == mu else 0 for nu in range(4)): Q2(co)
                                       for i, co in enumerate(bp[al[mu]]) if co != 0})
        poly = poly + term
    return poly

# ---------------------------------------------------------------- insertions
class Insertion:
    def __init__(self, bgs):
        z = np.zeros(4, dtype=complex)
        self.h0 = C.hessian(bgs, z, z)
        self.Hj = [C.hessian(bgs, E4[m].astype(complex), z) - self.h0 for m in range(4)]
        self.Hl = [C.hessian(bgs, z, E4[m].astype(complex)) - self.h0 for m in range(4)]
    def __call__(self, kj, kl, rows, cols):
        ix = np.ix_(rows, cols)
        h = self.h0[ix].copy()
        for m in range(4):
            if kj[m]: h = h + kj[m] * self.Hj[m][ix]
            if kl[m]: h = h + kl[m] * self.Hl[m][ix]
        return h


def support(L):
    M = C.linemat(L, np.array([1, 2, 3, 5]))
    return sorted(set(np.nonzero(np.abs(M).sum(axis=1))[0]) | set(np.nonzero(np.abs(M).sum(axis=0))[0]))


SUPP = {L: support(L) for L in C.PLINES + C.SLINES}


def perm_sign(seq):
    s_ = 1; a = list(seq)
    for i in range(len(a)):
        for j in range(len(a) - 1 - i):
            if a[j] > a[j + 1]:
                a[j], a[j + 1] = a[j + 1], a[j]; s_ = -s_
    return s_


def ordered_partitions(legs, types):
    legs = list(legs)
    if not legs:
        yield []
        return
    for i, j in enumerate(legs):
        rest = legs[:i] + legs[i + 1:]
        for tail in ordered_partitions(rest, types):
            yield [(j,)] + tail
        for jj in rest:
            if j < jj and types[j] == 'A' and types[jj] == 'A':
                rest2 = [x for x in rest if x != jj]
                for tail in ordered_partitions(rest2, types):
                    yield [(j, jj)] + tail


class SpecialCycles:
    def __init__(self, ext, degree):
        """ext = list of (type, which, pol, mom); degree = superficial degree of the numerator N(k)"""
        self.ext = ext; self.E = len(ext); self.degree = degree
        self.types = [e[0] for e in ext]
        self.ins = {}
        for j in range(self.E):
            self.ins[(j,)] = Insertion([(ext[j][0], ext[j][2], ext[j][3])])
        for j in range(self.E):
            for jj in range(j + 1, self.E):
                if self.types[j] == 'A' and self.types[jj] == 'A':
                    self.ins[(j, jj)] = Insertion([(ext[j][0], ext[j][2], ext[j][3]),
                                                   (ext[jj][0], ext[jj][2], ext[jj][3])])
        # nonzero pattern: insertion b between line Li (entering) and Lo (leaving)
        self.nz = {}
        kj = np.array([1, -2, 3, 1]); kl = np.array([2, 1, -1, 3])
        for b, ins in self.ins.items():
            for Li in C.PLINES + C.SLINES:
                for Lo in C.PLINES + C.SLINES:
                    h = ins(kj, kl, SUPP[Li], SUPP[Lo])
                    self.nz[(b, Li, Lo)] = bool(np.abs(h).sum() > 0)

    def mom(self, block):
        return sum((np.asarray(self.ext[j][3]) for j in block), np.zeros(4, dtype=int))

    def enumerate(self):
        out = []
        for seq in ordered_partitions(range(self.E), self.types):
            seq = [tuple(sorted(b)) for b in seq]
            n = len(seq)
            if n < 2:
                continue
            odd = [j for b in seq for j in b if C.PARITY[self.types[j]]]
            sgn = perm_sign(odd)
            def rec(i, acc):
                if i == n:
                    if self.nz[(seq[n - 1], acc[n - 1], acc[0])]:
                        yield list(acc)
                    return
                cands = C.SLINES if i == 0 else C.PLINES
                for L in cands:
                    if i > 0 and not self.nz[(seq[i - 1], acc[i - 1], L)]:
                        continue
                    yield from rec(i + 1, acc + [L])
            for Ls in rec(0, []):
                out.append((seq, Ls, sgn))
        return out

    def numerator(self, seq, Ls):
        n = len(seq)
        Ps = [self.mom(b) for b in seq]
        qs = [np.zeros(4, dtype=int)]
        for i in range(1, n):
            qs.append(qs[-1] - Ps[i - 1])
        S = C.SDIAG
        def N(k):
            ks = [k + q for q in qs]
            X = None
            for i in range(n):
                Si = SUPP[Ls[i]]; Sn = SUPP[Ls[(i + 1) % n]]
                Lm = C.linemat(Ls[i], ks[i])[np.ix_(Si, Si)]
                h = self.ins[seq[i]](-ks[i], ks[(i + 1) % n], Si, Sn)
                Mi = C.mm(Lm, S[Si][:, None] * h)
                X = Mi if X is None else C.mm(X, Mi)
            val = np.sum(S[SUPP[Ls[0]]] * np.diag(X))
            return -int(ks[0] @ ks[0]) * val
        return N, qs

    def graph_value(self, seq, Ls, sgn, orders=(-1, 0), check=True):
        """exact coefficients of eps^e (e in orders) of hbar^{-1} 16pi^2 * (graph value), as
        {e: {(n, beta, Npow): Q2}} (before multiplication by the exact integrals)"""
        n = len(seq)
        N, qs = self.numerator(seq, Ls)
        cache = {}
        def Nc(k):
            key = tuple(int(x) for x in k)
            if key not in cache:
                cache[key] = N(np.array(key))
            return cache[key]
        poly = interpolate(Nc, self.degree)
        if check:
            rng = np.random.default_rng(abs(hash((tuple(map(tuple, qs)), tuple(Ls)))) % 2 ** 32)
            for _ in range(2):
                kx = rng.integers(-3, 4, size=4)
                assert (poly.eval(kx) - Q2.of(N(kx))).iszero(), 'interpolation degree too low'
        ncub = sum(1 for b in seq if len(b) == 1); nqua = n - ncub
        pref = Fraction(sgn, 2) * Fraction(1, 2 ** ncub * 8 ** nqua)
        nu = n - 1
        # k = -K(u) = -sum_{i>=2} u_i q_i  (u-variables indexed 0..n-2 for lines 2..n)
        forms = []
        for mu in range(4):
            forms.append(Poly(nu, {tuple(1 if a == i else 0 for a in range(nu)): Q2(-int(qs[i + 1][mu]))
                                   for i in range(nu) if qs[i + 1][mu] != 0}))
        Kf = [f.scale(Q2(-1)) for f in forms]
        q = Poly(nu)
        for i in range(nu):
            q = q + Poly.var(nu, i, Q2(int(qs[i + 1] @ qs[i + 1])))
        for mu in range(4):
            q = q - Kf[mu] * Kf[mu]
        lapN = [poly]
        for j in range(1, self.degree // 2 + 1):
            lapN.append(lapN[-1].lap())
        res = {}
        for e in orders:
            acc = {}
            for j in range(len(lapN)):
                i = j - (n - 2 - e)
                if i < 0:
                    continue
                Fj = lapN[j].subs_linear(forms, nu)
                term = Fj
                mq = q.scale(Q2(-1))
                for _ in range(i):
                    term = term * mq
                c0 = Fraction(1, 4 ** j * math.factorial(j) * math.factorial(i)) * pref
                for beta, v in term.t.items():
                    key = (n, beta, 2 + j - i + sum(beta))
                    acc[key] = acc[key] + v * c0 if key in acc else v * c0
            res[e] = acc
        return res


def to_basis(acc):
    """{(n,beta,N): Q2} -> (Q2 rational part, Q2 coefficient of log 2, Q2 coefficient of log 3)"""
    r0, r2, r3 = Q2(0), Q2(0), Q2(0)
    for (n, beta, Np), v in acc.items():
        c0, c2, c3 = SI.J(n, beta, Np)
        r0 = r0 + v * Fraction(c0.p, c0.q); r2 = r2 + v * Fraction(c2.p, c2.q); r3 = r3 + v * Fraction(c3.p, c3.q)
    return (r0, r2, r3)


def canon_cycle(pairs):
    n = len(pairs); cands = []
    for r in range(n):
        rot = pairs[r:] + pairs[:r]
        cands.append(tuple(rot))
        Ls = [p[0] for p in rot]; Ss = [p[1] for p in rot]
        refl = [(Ls[0], Ss[n - 1])] + [(Ls[n - i], Ss[n - 1 - i]) for i in range(1, n)]
        cands.append(tuple(refl))
    return min(cands, key=lambda x: repr(x))


# ---------------------------------------------------------------- oriented lines (slot assignment of the two ends)
LINE_ORIENT = {'AV': [('A', 'A')], 'AP': [('A', 'A')], 'gh': [('c', 'As'), ('As', 'c')], 'l': [('l', 'l')],
               'sAV': [('A', 'As'), ('As', 'A')], 'sAP': [('A', 'As'), ('As', 'A')], 'sl': [('ls', 'l'), ('l', 'ls')],
               'sc': [('cs', 'c'), ('c', 'cs')]}
WHICH = {'AV': 'V', 'AP': 'P', 'sAV': 'V', 'sAP': 'P'}


def oriented_index(L, t):
    return C.comp_index(t, WHICH.get(L) if t in ('A', 'As') else None)


OLINES_P = [(L, o) for L in C.PLINES for o in LINE_ORIENT[L]]
OLINES_S = [(L, o) for L in C.SLINES for o in LINE_ORIENT[L]]
ROWS = {(L, o): oriented_index(L, o[0]) for (L, o) in OLINES_P + OLINES_S}
COLS = {(L, o): oriented_index(L, o[1]) for (L, o) in OLINES_P + OLINES_S}


class OrientedCycles(SpecialCycles):
    """as SpecialCycles, each line restricted to one slot assignment: rows (first factor) = slot at the previous
    insertion, columns (second factor) = slot at the next insertion; the sum over orientations is the full line."""
    def __init__(self, ext, degree):
        super().__init__(ext, degree)
        self.onz = {}
        kj = np.array([1, -2, 3, 1]); kl = np.array([2, 1, -1, 3])
        for b, ins in self.ins.items():
            for Li in OLINES_P + OLINES_S:
                for Lo in OLINES_P + OLINES_S:
                    h = ins(kj, kl, COLS[Li], ROWS[Lo])
                    self.onz[(b, Li, Lo)] = bool(np.abs(h).sum() > 0)

    def enumerate(self):
        out = []
        for seq in ordered_partitions(range(self.E), self.types):
            seq = [tuple(sorted(b)) for b in seq]
            n = len(seq)
            if n < 2:
                continue
            odd = [j for b in seq for j in b if C.PARITY[self.types[j]]]
            sgn = perm_sign(odd)
            def rec(i, acc):
                if i == n:
                    if self.onz[(seq[n - 1], acc[n - 1], acc[0])]:
                        yield list(acc)
                    return
                cands = OLINES_S if i == 0 else OLINES_P
                for L in cands:
                    if i > 0 and not self.onz[(seq[i - 1], acc[i - 1], L)]:
                        continue
                    yield from rec(i + 1, acc + [L])
            for Ls in rec(0, []):
                out.append((seq, Ls, sgn))
        return out

    def numerator(self, seq, Ls):
        n = len(seq)
        Ps = [self.mom(b) for b in seq]
        qs = [np.zeros(4, dtype=int)]
        for i in range(1, n):
            qs.append(qs[-1] - Ps[i - 1])
        S = C.SDIAG
        def N(k):
            ks = [k + q for q in qs]
            X = None
            for i in range(n):
                Li = Ls[i]; Ln = Ls[(i + 1) % n]
                Lm = C.linemat(Li[0], ks[i])[np.ix_(ROWS[Li], COLS[Li])]
                h = self.ins[seq[i]](-ks[i], ks[(i + 1) % n], COLS[Li], ROWS[Ln])
                Mi = C.mm(Lm, S[COLS[Li]][:, None] * h)
                X = Mi if X is None else C.mm(X, Mi)
            val = np.sum(S[ROWS[Ls[0]]] * np.diag(X))
            return -int(ks[0] @ ks[0]) * val
        return N, qs


def canon_ocycle(pairs):
    """canonical form of a cyclic sequence of (oriented line, block label) under rotation and reflection
    (reflection reverses every line)"""
    n = len(pairs); cands = []
    rev = lambda Lo: (Lo[0], (Lo[1][1], Lo[1][0]))
    for r in range(n):
        rot = pairs[r:] + pairs[:r]
        cands.append(tuple(rot))
        Ls = [p[0] for p in rot]; Ss = [p[1] for p in rot]
        refl = [(rev(Ls[0]), Ss[n - 1])] + [(rev(Ls[n - i]), Ss[n - 1 - i]) for i in range(1, n)]
        cands.append(tuple(refl))
    return min(cands, key=lambda x: repr(x))


# ---------------------------------------------------------------- finite-eps data (for the numerical checks)
def graph_polys(G, seq, Ls, sgn, special=True):
    """exact data of one oriented sequence: n, pref, {j: F_j(u)}, q(u) with u = (u_2..u_n);
    value(eps) = (16pi^2)^{-1} pref sum_j (4^j j!)^{-1} int dt T^{-2-j} F_j(u) e^{-T q(u)} over the Schwinger region"""
    n = len(seq)
    N, qs = G.numerator(seq, Ls)
    cache = {}
    def Nc(k):
        key = tuple(int(x) for x in k)
        if key not in cache:
            cache[key] = N(np.array(key))
        return cache[key]
    poly = interpolate(Nc, G.degree)
    ncub = sum(1 for b in seq if len(b) == 1); nqua = n - ncub
    pref = (Fraction(sgn, 2) if special else Fraction(sgn, 2 * n)) * Fraction(1, 2 ** ncub * 8 ** nqua)
    nu = n - 1
    forms = []
    for mu in range(4):
        forms.append(Poly(nu, {tuple(1 if a == i else 0 for a in range(nu)): Q2(-int(qs[i + 1][mu]))
                               for i in range(nu) if qs[i + 1][mu] != 0}))
    Kf = [f.scale(Q2(-1)) for f in forms]
    q = Poly(nu)
    for i in range(nu):
        q = q + Poly.var(nu, i, Q2(int(qs[i + 1] @ qs[i + 1])))
    for mu in range(4):
        q = q - Kf[mu] * Kf[mu]
    Fj = {}
    lp = poly
    for j in range(G.degree // 2 + 1):
        Fj[j] = lp.subs_linear(forms, nu)
        lp = lp.lap()
    return n, pref, Fj, q


class PCycles(OrientedCycles):
    """ordinary one-loop graphs (all edges P(eps,infty)), weight 1/(2n) str[(P S H)^n]; n = 1 (quartic tadpole)
    included"""
    def enumerate(self):
        out = []
        for seq in ordered_partitions(range(self.E), self.types):
            seq = [tuple(sorted(b)) for b in seq]
            n = len(seq)
            odd = [j for b in seq for j in b if C.PARITY[self.types[j]]]
            sgn = perm_sign(odd)
            def rec(i, acc):
                if i == n:
                    if self.onz[(seq[n - 1], acc[n - 1], acc[0])]:
                        yield list(acc)
                    return
                for L in OLINES_P:
                    if i > 0 and not self.onz[(seq[i - 1], acc[i - 1], L)]:
                        continue
                    yield from rec(i + 1, acc + [L])
            for Ls in rec(0, []):
                out.append((seq, Ls, sgn))
        return out

    def numerator(self, seq, Ls):
        n = len(seq)
        Ps = [self.mom(b) for b in seq]
        qs = [np.zeros(4, dtype=int)]
        for i in range(1, n):
            qs.append(qs[-1] - Ps[i - 1])
        S = C.SDIAG
        def N(k):
            ks = [k + q for q in qs]
            X = None
            for i in range(n):
                Li = Ls[i]; Ln = Ls[(i + 1) % n]
                Lm = C.linemat(Li[0], ks[i])[np.ix_(ROWS[Li], COLS[Li])]
                h = self.ins[seq[i]](-ks[i], ks[(i + 1) % n], COLS[Li], ROWS[Ln])
                Mi = C.mm(Lm, S[COLS[Li]][:, None] * h)
                X = Mi if X is None else C.mm(X, Mi)
            return np.sum(S[ROWS[Ls[0]]] * np.diag(X))
        return N, qs
