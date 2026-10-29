"""Engine A, round 2: exact enumeration and evaluation of one-loop 1PI graphs (method A: momentum space, per-edge
Schwinger parameters t_e in [eps,infty), exact Gaussian loop integral, exact eps-asymptotics).

Graph = cyclic sequence of insertions (one external leg on a cubic vertex, or two external A legs on the quartic
vertex) joined by internal lines of type AV (A|_V - A|_V), AP (A|_{W^perp} - A|_{W^perp}), gh (c - A^*|_V, both
orientations), l (lambda - lambda).  Value = sum over ordered sequences in the class of
  1/(2n) * sign(odd background parameters -> leg order) * str[ N_{L1}(k_1) S h_1 N_{L2}(k_2) S h_2 ... ],
line i carries k_i = k - (P_{s_1}+...+P_{s_{i-1}}), Schwinger weight e^{-t_i k_i^2}, t_i >= eps.
Gaussian integral: int d^4k/(2pi)^4 N e^{-sum t_i k_i^2} = (16 pi^2)^{-1} T^{-2} E_l[N(l - K(u))] e^{-T q(u)},
T = sum t_i, u_i = t_i/T, K(u) = sum u_i q_i (q_i = k_i - k), q(u) = sum u_i q_i^2 - K(u)^2, l Gaussian with covariance
1/(2T); E over the degree-2m part = (4^m m! T^m)^{-1} (Delta^m N)(k=-K(u)).
eps-asymptotics (region t_i >= eps  <=>  T >= eps/min_i u_i):
  log eps coefficient = - sum_j int_simplex du (-q(u))^j/j! (Delta^{n-2+j} N)(-K(u)) / (4^{n-2+j} (n-2+j)!)
  eps^{-1} coefficient = int_simplex du min_i(u_i) (Delta^{n-1} N)/(4^{n-1} (n-1)!)   (degree D = 2n-2 graphs)
Units: hbar/(16 pi^2) (the factor hbar of the one-loop term and (16 pi^2)^{-1} stripped)."""
import itertools, math, json, sys
from fractions import Fraction
import numpy as np
import r2_core as C
from r2_core import Q2

E4 = [np.eye(4, dtype=int)[m] for m in range(4)]
LINE_ENDS = {'AV': ({'A'}, {'A'}), 'AP': ({'A'}, {'A'}), 'gh': ({'c', 'As'}, {'c', 'As'}), 'l': ({'l'}, {'l'})}


class Insertion:
    """affine Hessian h(kj, kl) = h0 + sum_mu kj_mu Hj[mu] + kl_mu Hl[mu]  (exact, linear in momenta), sliced to idx"""
    def __init__(self, bgs, fluct=('c', 'A', 'As', 'l'), idx=None):
        z = np.zeros(4, dtype=complex)
        ix = np.ix_(idx, idx) if idx is not None else (slice(None), slice(None))
        h0 = C.hessian(bgs, z, z, fluct)
        self.h0 = h0[ix]
        self.Hj = [(C.hessian(bgs, E4[m].astype(complex), z, fluct) - h0)[ix] for m in range(4)]
        self.Hl = [(C.hessian(bgs, z, E4[m].astype(complex), fluct) - h0)[ix] for m in range(4)]
    def __call__(self, kj, kl):
        h = self.h0.copy()
        for m in range(4):
            if kj[m]: h = h + kj[m] * self.Hj[m]
            if kl[m]: h = h + kl[m] * self.Hl[m]
        return h


def ordered_partitions(legs, types):
    """ordered sequences of blocks (size 1, or size 2 when both legs are A) covering all legs"""
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
                    # the pair may also be listed with jj first only once (blocks are sets)


def perm_sign(seq):
    s_ = 1; a = list(seq)
    for i in range(len(a)):
        for j in range(len(a) - 1 - i):
            if a[j] > a[j + 1]:
                a[j], a[j + 1] = a[j + 1], a[j]; s_ = -s_
    return s_


def canon_cycle(pairs):
    """canonical representative of a cyclic sequence of (line, block) under rotation and reflection"""
    n = len(pairs); cands = []
    for r in range(n):
        rot = pairs[r:] + pairs[:r]
        cands.append(tuple(rot))
        # reflection: (L1,s_n),(L_n,s_{n-1}),...,(L2,s_1)
        Ls = [p[0] for p in rot]; Ss = [p[1] for p in rot]
        refl = [(Ls[0], Ss[n - 1])] + [(Ls[n - i], Ss[n - 1 - i]) for i in range(1, n)]
        cands.append(tuple(refl))
    return min(cands, key=lambda x: repr(x))


def lap(f, k0):
    tot = 0
    f0 = f(k0)
    for e in E4:
        tot = tot + f(k0 + e) + f(k0 - e) - 2 * f0
    return tot


def lap2(f, k0):
    """sum_{mu,nu} delta_mu^2 delta_nu^2 f : exact Delta^2 for polynomials of degree <= 5"""
    tot = 0
    for m in range(4):
        for n_ in range(4):
            em, en = E4[m], E4[n_]
            if m == n_:
                tot = tot + f(k0 + 2 * em) - 4 * f(k0 + em) + 6 * f(k0) - 4 * f(k0 - em) + f(k0 - 2 * em)
            else:
                for a, ca in ((1, 1), (0, -2), (-1, 1)):
                    for b, cb in ((1, 1), (0, -2), (-1, 1)):
                        tot = tot + ca * cb * f(k0 + a * em + b * en)
    return tot


def cached(f):
    memo = {}
    def g(k):
        key = tuple(int(x) for x in k)
        if key not in memo:
            memo[key] = f(np.array(key))
        return memo[key]
    return g


def lapm(f, k0, m):
    f = cached(f)
    if m == 0: return f(k0)
    if m == 1: return lap(f, k0)
    if m == 2: return lap2(f, k0)
    raise ValueError(m)


class GraphSum:
    def __init__(self, ext, fluct=('c', 'A', 'As', 'l'), lines=('AV', 'AP', 'gh', 'l')):
        """ext = list of (type, which, pol, mom); which in {'V','P',None}"""
        self.ext = ext; self.E = len(ext)
        self.types = [e[0] for e in ext]
        self.lines = lines
        self.fluct = fluct
        idx = set()
        for L in lines:
            Nt = C.Pnum(np.array([1, 2, 3, 5]), L)
            idx |= set(np.nonzero(np.abs(Nt).sum(axis=1))[0]) | set(np.nonzero(np.abs(Nt).sum(axis=0))[0])
        self.idx = sorted(idx)
        self.ix = np.ix_(self.idx, self.idx)
        self.S = C.SDIAG[self.idx]
        self.ins = {}
        for j in range(self.E):
            self.ins[(j,)] = Insertion([(ext[j][0], ext[j][2], ext[j][3])], fluct, self.idx)
        for j in range(self.E):
            for jj in range(j + 1, self.E):
                if self.types[j] == 'A' and self.types[jj] == 'A':
                    self.ins[(j, jj)] = Insertion([(ext[j][0], ext[j][2], ext[j][3]), (ext[jj][0], ext[jj][2], ext[jj][3])], fluct, self.idx)
        self.mask = {}
        for L in lines:
            Nt = C.Pnum(np.array([1, 2, 3, 5]), L)[self.ix]
            self.mask[L] = (np.abs(Nt).sum(axis=1) > 0, np.abs(Nt).sum(axis=0) > 0)   # rows (first factor), cols (second)
        self.nz = {}
        kj = np.array([1, -2, 3, 1]); kl = np.array([2, 1, -1, 3])
        for b, ins in self.ins.items():
            h = ins(kj, kl)
            for Li in lines:
                for Lo in lines:
                    rows = self.mask[Li][1]; cols = self.mask[Lo][0]
                    self.nz[(b, Li, Lo)] = bool(np.abs(h[np.ix_(rows, cols)]).sum() > 0)

    def mom(self, block):
        return sum((np.asarray(self.ext[j][3]) for j in block), np.zeros(4, dtype=int))

    def numerator(self, seq, Ls):
        n = len(seq)
        Ps = [self.mom(b) for b in seq]
        qs = [np.zeros(4, dtype=int)]
        for i in range(1, n):
            qs.append(qs[-1] - Ps[i - 1])
        def N(k):
            ks = [k + q for q in qs]
            X = None
            for i in range(n):
                Pm = C.Pnum(ks[i], Ls[i])[self.ix]
                h = self.ins[seq[i]](-ks[i], ks[(i + 1) % n])
                Mi = C.mm(Pm, self.S[:, None] * h)
                X = Mi if X is None else C.mm(X, Mi)
            return np.sum(self.S * np.diag(X))
        return N, qs

    def enumerate(self):
        out = []
        for seq in ordered_partitions(range(self.E), self.types):
            seq = [tuple(sorted(b)) for b in seq]
            n = len(seq)
            odd = [j for b in seq for j in b if C.PARITY[self.types[j]]]
            sgn = perm_sign(odd)
            # line assignments: line i enters insertion i (j-leg) and leaves insertion i-1 (l-leg)
            def rec(i, acc):
                if i == n:
                    # closing: line 0 leaves insertion n-1
                    if self.nz[(seq[n - 1], acc[n - 1], acc[0])]:
                        yield list(acc)
                    return
                for L in self.lines:
                    if i > 0 and not self.nz[(seq[i - 1], acc[i - 1], L)]:
                        continue
                    yield from rec(i + 1, acc + [L])
            for Ls in rec(0, []):
                out.append((seq, Ls, sgn))
        return out

    def coefficients(self, omega, verbose=False):
        """returns dict class_key -> {'log': Q2, 'epsinv': Q2} in units hbar/(16 pi^2), vertex scale removed"""
        res = {}
        for (seq, Ls, sgn) in self.enumerate():
            n = len(seq)
            N, qs = self.numerator(seq, Ls)
            ncub = sum(1 for b in seq if len(b) == 1); nqua = n - ncub
            scale = Fraction(1, 2 ** ncub * 8 ** nqua)
            pref = Fraction(sgn, 2 * n) * scale
            key = canon_cycle([(Ls[i], seq[i]) for i in range(n)])
            if omega == 0:
                m = n - 2
                val = Q2.of(lapm(N, np.zeros(4, dtype=int), m)) * Fraction(-1, math.factorial(n - 1) * 4 ** m * math.factorial(m))
                ent = {'log': val, 'epsinv': Q2(0)}
            elif omega == 1:
                m = n - 2
                kbar_n = sum((sum(qs[:i + 1], np.zeros(4, dtype=int)) for i in range(n)), np.zeros(4, dtype=int))
                # -K(ubar) = -(1/n) sum q_i ; require divisibility
                kb = -sum(qs, np.zeros(4, dtype=int))
                assert np.all(kb % n == 0), (kb, n)
                kb = kb // n
                val = Q2.of(lapm(N, kb, m)) * Fraction(-1, math.factorial(n - 1) * 4 ** m * math.factorial(m))
                ent = {'log': val, 'epsinv': Q2(0)}
            elif omega == 2:
                ent = self._omega2(N, qs, n)
            else:
                raise ValueError(omega)
            ent = {k: v * pref for k, v in ent.items()}
            if key not in res:
                res[key] = {'log': Q2(0), 'epsinv': Q2(0), 'nseq': 0}
            res[key]['log'] = res[key]['log'] + ent['log']
            res[key]['epsinv'] = res[key]['epsinv'] + ent['epsinv']
            res[key]['nseq'] += 1
        return res

    def _omega2(self, N, qs, n):
        if n == 1:
            return {'epsinv': Q2.of(N(np.zeros(4, dtype=int))), 'log': Q2(0)}
        assert n == 2
        P = -qs[1]                       # q_2 = -P
        assert np.all(P % 2 == 0)
        # -K(u) = u_2 P ; N quadratic in u_2 -> Simpson
        g0 = Q2.of(N(0 * P)); g1 = Q2.of(N(P // 2)); g2 = Q2.of(N(P))
        intg = (g0 + g1 * 4 + g2) * Fraction(1, 6)
        cc = Q2.of(lap(N, 0 * P)) * Fraction(1, 4)     # Delta N / 4 (constant)
        s2 = int(P @ P)
        log = -intg + cc * Fraction(s2, 6)
        epsinv = cc * Fraction(1, 4)
        return {'log': log, 'epsinv': epsinv}
