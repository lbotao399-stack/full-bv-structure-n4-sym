"""Mode model (Part-C engine s2_ga, exact Grassmann algebra, manuscript bracket and Q) for the tree-level operations on
quadratic functionals: Q (= {S_free, -}), the bracket {F, G} (C07 formula with Delta_0), plane-wave values.
Copies of s2_ga.py, s2_core.py, grass.py taken from the project's GitHub mirror (the local originals are iCloud-evicted;
sizes equal), imported read-only from gh/.  Colour su(2) orthonormal; legs carry colour e_0 so that the values equal the
colour-stripped values of the cycle engine ((Y0,Y1) -> 1)."""
import os, sys
from fractions import Fraction
import numpy as np
import sa_core as SA            # loads ea_core, r2_core, cliff (vendored) first
GH = os.path.join(SA.SRC, 'gh')
if GH not in sys.path:
    sys.path.insert(0, GH)
import grass                    # noqa: gh copy (must precede s2_ga's own path insertion)
import s2_core as C2            # noqa
import s2_ga as GM              # noqa
from grass import GA            # noqa
from r2_core import Q2          # noqa
import sa_fit as SF             # noqa

SB = [np.asarray(m) for m in SA.E.SGB]
SS = [np.asarray(m) for m in SA.E.SG]
SHP = {'c': 1, 'A': 10, 'As': 10, 'l': 16, 'ls': 16, 'cs': 1}


class MM:
    """two-mode model p, -p (plus optional extra pairs)"""

    def __init__(self, p):
        p = np.asarray(p, dtype=int)
        self.M = GM.Model([p, -p])
        self.pairs = [(0, 1), (1, 0)]
        self.pm = [(0, 1), (1, 0)]
        self.Sfree = self.M.S_free(self.pairs)

    def quad(self, s1, s2, fn):
        """int sum_{i,j} fn(i, j, p1, p2) (phi1_i, phi2_j) with the colour pairing G; s1, s2 slot names
        ('c','A','As','l','ls','cs'); fn returns a complex (Gaussian rational) coefficient, p1, p2 the momenta of the
        two slots (derivatives d -> i p)."""
        acc = {}
        G = C2.G; D = C2.D
        gnz = np.argwhere(G != 0)
        for (m1, m2) in self.pairs:
            p1, p2 = self.M.moms[m1], self.M.moms[m2]
            for i in range(SHP[s1]):
                for j in range(SHP[s2]):
                    c = fn(i, j, p1, p2)
                    if c == 0:
                        continue
                    for (a, b) in gnz:
                        ia = (i, a) if s1 not in ('c', 'cs') else (a,)
                        jb = (j, b) if s2 not in ('c', 'cs') else (b,)
                        self.M.add(acc, SF.q2(c) * GM.q2(G[a, b]), [self.M.y(s1, ia, m1), self.M.y(s2, jb, m2)])
        return self.M.ga(acc)

    def Q(self, F):
        return self.M.Q_apply(F, self.Sfree, self.pm)

    def br(self, F, G_, pf):
        return self.M.bracket0(F, G_, pf, self.pm)

    def value(self, F, legs):
        """legs: [(slot, stripped pol vector)], first at mode 0 (momentum p), second at mode 1 (-p); colour e_0.
        Returns the coefficient of the odd leg parameters in leg order (Q2)."""
        rs = GM.rvars(len(legs))
        L = []
        odd = []
        for k, (s, pol) in enumerate(legs):
            full = np.zeros(C2.SHAPE[s], dtype=complex)
            if s in ('c', 'cs'):
                full[0] = pol[0]
            else:
                full[:, 0] = pol
            r = rs[k] if C2.PARITY[s] else None
            if r is not None:
                odd.append(r)
            L.append((s, full, k, r))
        sub = self.M.subst(F, L, [0, 1])
        return GM.coeff(sub, odd)


# ------------------------------------------------------------------------------------------------ functionals
def d(p, mu):
    return 1j * p[mu]


def F_T1(mm, xi, Ms=range(4)):
    """-int (Gamma(xi,lam)|_Ms, A^*|_Ms)"""
    G = [xi @ SB[M] for M in range(10)]
    return mm.quad('l', 'As', lambda i, j, p1, p2: -G[j][i] if j in Ms else 0)


def F_T3(mm, xi, Ns=range(4)):
    """1/2 int (rho(d A|_Ns) xi, lam^*), rho(F) xi = sum_{M<N} F_MN s[M] sb[N] xi, F_MN = d_M A_N - d_N A_M:
    = 1/2 sum_{M != N} d_M A_N (s[M] sb[N] xi) (M in V)"""
    X = {}
    for M in range(4):
        for N in Ns:
            if M != N:
                X[(M, N)] = SS[M] @ (SB[N] @ xi)
    def fn(i, j, p1, p2):
        return 0.5 * sum(d(p1, M) * X[(M, i)][j] for M in range(4) if (M, i) in X)
    return mm.quad('A', 'ls', fn)


def F_T0(mm, xi):
    """int (d^* A_V, (xi, lam^*)),  d^*A = -d_mu A_mu"""
    return mm.quad('A', 'ls', lambda i, j, p1, p2: (-d(p1, i) * xi[j]) if i < 4 else 0)


def F_g(mm, xi, which):
    G = [xi @ SB[M] for M in range(10)]
    if which == 1:   # int (Gamma_V, d^*d A_V): (d^*dA)_mu -> |p|^2 A_mu - p_mu (p.A)
        return mm.quad('A', 'l', lambda i, j, p1, p2: sum(G[mu][j] * ((p1 @ p1) * (i == mu) - p1[mu] * p1[i])
                                                            for mu in range(4)) if i < 4 else 0)
    if which == 2:   # int (Gamma_V, d d^* A_V): (dd^*A)_mu -> p_mu (p.A)
        return mm.quad('A', 'l', lambda i, j, p1, p2: sum(G[mu][j] * p1[mu] * p1[i] for mu in range(4)) if i < 4 else 0)
    if which == 3:   # int (Gamma_perp, d^*d A_perp)
        return mm.quad('A', 'l', lambda i, j, p1, p2: G[i][j] * (p1 @ p1) if i >= 4 else 0)
    if which == 4:   # int (d^*dc, (xi, lam^*))
        return mm.quad('c', 'ls', lambda i, j, p1, p2: (p1 @ p1) * xi[j])
    raise ValueError(which)


def F_GA(mm, xi, Ms):
    """int (Gamma(xi,lam)|_Ms, A|_Ms)"""
    G = [xi @ SB[M] for M in range(10)]
    return mm.quad('A', 'l', lambda i, j, p1, p2: G[i][j] if i in Ms else 0)


def F_named(mm, name):
    return mm.M.quad(name, mm.pairs)


_M3 = {}


def value3(name, L):
    """value of the cubic functional `name` of s2_ga (G1..G5, AAA, ...) on three cycle-engine legs (colours e_0, e_1, e_2
    so that (Y0,[Y1,Y2]) = 1); coefficient of the odd leg parameters in leg order"""
    import itertools
    moms = tuple(tuple(int(round(x.real)) for x in l['p0']) for l in L)
    if moms not in _M3:
        _M3[moms] = GM.Model([np.array(m) for m in moms])
    M = _M3[moms]
    F = M.cubic(name, list(itertools.permutations(range(3))))
    rs = GM.rvars(3)
    legs, odd = [], []
    for k, l in enumerate(L):
        s = SA.E.SLOT[l['t']]
        full = np.zeros(C2.SHAPE[s], dtype=complex)
        if s in ('c', 'cs'):
            full[k] = l['pol'][0]
        else:
            full[:, k] = l['pol']
        r = rs[k] if C2.PARITY[s] else None
        if r is not None:
            odd.append(r)
        legs.append((s, full, k, r))
    sub = M.subst(F, legs, [0, 1, 2])
    v = GM.coeff(sub, odd)
    return complex(float(v.re), float(v.im)) if v.re.denominator == 1 and v.im.denominator == 1 else (v.re, v.im)
