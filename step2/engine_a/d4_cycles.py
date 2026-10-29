"""Step 2, Part D4, engine A: one-loop cycles with one vertex S_M = S^(1)(Delta) (ex:c2-dil), all edges P(eps,infty).

S_M = int( (2E(A)+2A, A^*) - (2E(lambda)+3 lambda, lambda^*) - (2E(c), c^*) ),  E = x^mu d_mu.
Split  S_M = sum_mu int x^mu s_mu + int s_0,
   s_0  = (2A, A^*) - (3 lambda, lambda^*),     s_mu = (2 d_mu A, A^*) - (2 d_mu lambda, lambda^*) - (2 d_mu c, c^*),
and  int x^mu s_mu = -i d/dr^mu |_{r=0} int e^{i(r,x)} s_mu.  For a graph Gamma with external legs (first leg momentum p
held fixed, the momentum of a designated dependent leg absorbing -r), with values G_0 (vertex s_0) and G_mu(r) (vertex
int e^{i(r,x)} s_mu, momentum r injected), both with (2pi)^4 delta(sum) stripped:
   Gamma[S_M] = alpha(p) delta(Sigma) + beta^mu(p) (d_mu delta)(Sigma),   Sigma = total external momentum,
   alpha = G_0 - i sum_mu dG_mu/dr^mu |_0 ,   beta^mu = -i G_mu(0)          ("canonical form", derivative at fixed p).
Check: S_M itself (legs A at q, A^* at q') has alpha = 2 (a,b), beta^mu = 2 q^mu (a,b).

Cycle formula (validated in 0 dimensions, engine_a/verify_oneloop0d.py, s2_verify0d.py):
   1PI part of W_1 = (1/2) sum over oriented sequences with the S_M insertion first of
   str[ L_1 S h_1 L_2 S h_2 ... L_n S h_n ],  L_i = Pnum_{line i}(k_i) int_eps^infty dt_i e^{-t_i k_i^2},
   k_1 = k, k_{i+1} = k_i - (momentum injected at insertion i), h_i = Hessian at insertion i with fluctuation
   momenta (-k_i, k_{i+1}).  Odd parameters stand to the left in sequence order (background legs: r_j; S_M: xi with
   Hessian of xi S_M = -xi * Hess(S_M), since every monomial of S_M is odd) and are reordered to (xi, r_legs...).
Gaussian integral and Schwinger expansion: see d4_integrals.py."""
import os, sys, itertools, math
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import s2_core as C
from r2_core import Q2
import s2_graphs as SG
from s2_graphs import Poly, interpolate

E4 = [np.eye(4, dtype=int)[m] for m in range(4)]
GG = C.G


def pairing(t1, x, t2, y):
    """fibre pairing of batches x (slot t1) and y (slot t2): A-A^*, lambda-lambda^*, c-c^*; returns (n1, n2) array"""
    if t1 in ('A', 'As'):
        return np.einsum('iMa,kMb,ab->ik', x, y, GG)
    if t1 in ('l', 'ls'):
        return np.einsum('iAa,kAb,ab->ik', x, y, GG)
    return np.einsum('ia,kb,ab->ik', x, y, GG)


class QF:
    """quadratic local functional insertion with injected momentum `mom`:
       kind 's0'  : (2A,A^*) - (3 lambda, lambda^*)
       kind 'smu' : (2 d_mu A, A^*) - (2 d_mu lambda, lambda^*) - (2 d_mu c, c^*)      [x e^{i(r,x)}]
       kind 'O'   : (iota_w1 A, iota_w2 A)                                              [x e^{i(k,x)}]
       kind 'AAs' : (A|_comps, A^*|_comps)  (used for counterterm-type insertions)"""
    def __init__(self, kind, mom, mu=None, w1=None, w2=None, comps=None):
        self.kind = kind; self.mom = np.asarray(mom); self.mu = mu; self.w1 = w1; self.w2 = w2; self.comps = comps
        self.odd = kind in ('s0', 'smu', 'AAs')          # S_M-type functionals are odd
        z = np.zeros(4)
        self.h0 = self.full(z, z)
        self.Hj = [self.full(E4[m].astype(float), z) - self.h0 for m in range(4)]
        self.Hl = [self.full(z, E4[m].astype(float)) - self.h0 for m in range(4)]

    def vfun(self, tj, kj, tl, kl):
        """vertex function (coefficient of odd parameters in leg order) on identity batches; None if absent"""
        pj, pl = C.ident(tj), C.ident(tl)
        k = self.kind
        if k == 'O':
            if tj == 'A' and tl == 'A':
                a = np.einsum('iMa,M->ia', pj, np.eye(10)[self.w1]); b = np.einsum('kNb,N->kb', pl, np.eye(10)[self.w2])
                a2 = np.einsum('iMa,M->ia', pj, np.eye(10)[self.w2]); b2 = np.einsum('kNb,N->kb', pl, np.eye(10)[self.w1])
                return np.einsum('ia,kb,ab->ik', a, b, GG) + np.einsum('ia,kb,ab->ik', a2, b2, GG)
            return None
        pairs = {('A', 'As'), ('l', 'ls'), ('c', 'cs')}
        if (tj, tl) in pairs:
            fld, kf, first = tj, kj, True
        elif (tl, tj) in pairs:
            fld, kf, first = tl, kl, False
        else:
            return None
        if k == 'AAs':
            if fld != 'A':
                return None
            fac = 1.0
            mask = np.zeros(10); mask[list(self.comps)] = 1
        elif k == 's0':
            fac = {'A': 2.0, 'l': -3.0, 'c': 0.0}[fld]
            mask = None
        elif k == 'smu':
            fac = {'A': 2j, 'l': -2j, 'c': -2j}[fld] * kf[self.mu]
            mask = None
        if fac == 0:
            return None
        if fld == 'A' and mask is not None:
            pj = pj * mask[None, :, None] if tj == 'A' else pj
            pl = pl * mask[None, :, None] if tl == 'A' else pl
        if first:
            v = pairing(tj, pj, tl, pl)
        else:
            v = pairing(tl, pl, tj, pj).T
        return fac * v

    def full(self, kj, kl):
        H = np.zeros((C.NTOT, C.NTOT), dtype=complex)
        for tj in C.TYPES:
            for tl in C.TYPES:
                v = self.vfun(tj, kj, tl, kl)
                if v is None:
                    continue
                sg = (-1) ** (C.PARITY[tj] * C.PARITY[tl])
                H[C.OFF[tj]:C.OFF[tj] + C.SIZE[tj], C.OFF[tl]:C.OFF[tl] + C.SIZE[tl]] += sg * v
        return H

    def __call__(self, kj, kl, rows, cols):
        ix = np.ix_(rows, cols)
        h = self.h0[ix].copy()
        for m in range(4):
            if kj[m]: h = h + kj[m] * self.Hj[m][ix]
            if kl[m]: h = h + kl[m] * self.Hl[m][ix]
        return h


class Bg:
    """insertion of the cubic/quartic vertices of I at background legs [(type, pol, mom)] (scaled by 2 / 8 as s2_core)"""
    def __init__(self, bgs):
        self.bgs = bgs
        self.mom = sum((np.asarray(b[2]) for b in bgs), np.zeros(4, dtype=int))
        self.ins = SG.Insertion([(t, pol, np.asarray(m)) for (t, pol, m) in bgs])
        self.nodd = sum(C.PARITY[t] for (t, _, _) in bgs)
        self.odd = self.nodd % 2 == 1

    def __call__(self, kj, kl, rows, cols):
        return self.ins(kj, kl, rows, cols)


SUPP = SG.SUPP


def nz_pattern(ins, Li, Lo):
    kj = np.array([1, -2, 3, 1]); kl = np.array([2, 1, -1, 3])
    h = ins(kj, kl, SUPP[Li], SUPP[Lo])
    return bool(np.abs(h).sum() > 0)


class Cycle:
    """insertions: list; insertion 0 is first in every sequence.  scale: product of vertex scales to divide."""
    def __init__(self, insertions, degree):
        self.ins = insertions; self.n = len(insertions); self.degree = degree

    def sequences(self):
        """yields (order, Ls): order = tuple of insertion indices starting with 0, Ls = line types (line i enters order[i])"""
        n = self.n
        for perm in itertools.permutations(range(1, n)):
            order = (0,) + perm
            insl = [self.ins[i] for i in order]
            def rec(i, acc):
                if i == n:
                    if nz_pattern(insl[n - 1], acc[n - 1], acc[0]):
                        yield list(acc)
                    return
                for L in C.PLINES:
                    if i > 0 and not nz_pattern(insl[i - 1], acc[i - 1], L):
                        continue
                    yield from rec(i + 1, acc + [L])
            for Ls in rec(0, []):
                yield order, Ls

    def offsets(self, order):
        qs = [np.zeros(4, dtype=int)]
        for i in range(1, self.n):
            qs.append(qs[-1] - self.ins[order[i - 1]].mom)
        return qs

    def numerator_fn(self, order, Ls):
        qs = self.offsets(order)
        S = C.SDIAG
        insl = [self.ins[i] for i in order]
        n = self.n
        def N(k):
            ks = [k + q for q in qs]
            X = None
            for i in range(n):
                Si = SUPP[Ls[i]]; Sn = SUPP[Ls[(i + 1) % n]]
                Lm = C.linemat(Ls[i], np.asarray(ks[i], dtype=complex))[np.ix_(Si, Si)]
                h = insl[i](-np.asarray(ks[i], dtype=complex), np.asarray(ks[(i + 1) % n], dtype=complex), Si, Sn)
                Mi = C.mm(Lm, S[Si][:, None] * h)
                X = Mi if X is None else C.mm(X, Mi)
            return np.sum(S[SUPP[Ls[0]]] * np.diag(X))
        return N, qs

    def parity_sign(self, order, canon):
        """sign from moving the odd parameters (in sequence order) to the canonical order `canon` (list of insertion
        indices that carry one odd parameter each)"""
        seq = [i for i in order if self.ins[i].odd]
        assert sorted(seq) == sorted(canon), (seq, canon)
        pos = [canon.index(i) for i in seq]
        return SG.perm_sign(pos)


def poly_from(N, degree, check=True):
    cache = {}
    def Nc(k):
        key = tuple(int(x) for x in k)
        if key not in cache:
            cache[key] = N(np.array(key))
        return cache[key]
    poly = interpolate(Nc, degree)
    if check:
        rng = np.random.default_rng(11)
        for _ in range(2):
            kx = rng.integers(-3, 4, size=4)
            assert (poly.eval(kx) - Q2.of(N(kx))).iszero(), 'interpolation degree too low'
    return poly
