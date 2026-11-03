"""Step 2, engine A: mode model for tree-level BV operations on cubic and quadratic local functionals.

Functionals are evaluated on finitely many plane-wave modes (momenta p_0, p_1, p_2 with p_0 + p_1 + p_2 = 0, and
their negatives), with (2 pi)^4 delta stripped (as engine_a/free.py, generalised to colour and several modes).  Every
field component at every mode is a variable of the super-commutative algebra of engine_a/grass.py (parity = total
parity of the slot).  Operations, all with the conventions of the manuscript (ch01, ch03, ch04):
  Q F            = {S_free, F}, applied as the odd derivation fixed by Q on each mode variable,
  d_X            = (1/2) sum X_uv d_v d_u  (left derivatives; X a kernel pairing a mode at p with a mode at -p),
  C_X(F,G)       = d_X(FG) - d_X(F) G - F d_X(G),   {F,G}_X = -C_X(F,G),
  K_0            = -sum (x_i (x) xi_i + xi_i (x) x_i) (x) (-,-)^{-1}  (Khat of s2_core), K_t = e^{-t p^2} Khat,
  P(eps,infty)   = e^{-eps p^2}/p^2 Pnum(p)  (first factor at momentum p; s2_core.Pnum).
External legs are substituted by r_j u_j (r_j an odd parameter for an odd leg); the plane-wave value is the coefficient
of the odd parameters in leg order (vertex-function convention)."""
import itertools
from fractions import Fraction
import numpy as np
import s2_core as C
from r2_core import Q2
import sys
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import grass
from grass import GA, PAR, newvar, _mul_keys
from cliff import sb

SB = [np.rint(np.array(m).real).astype(int) + 1j * np.rint(np.array(m).imag).astype(int) for m in sb]
_next = [10 ** 6]


def q2(z):
    if isinstance(z, Q2):
        return z
    if isinstance(z, (int, Fraction)):
        return Q2(z)
    z = complex(z)
    return Q2.of(z)


class Model:
    def __init__(self, moms):
        """moms: list of integer 4-vectors (distinct); mode index = position"""
        self.moms = [np.asarray(m, dtype=int) for m in moms]
        self.var = {}
        for mi in range(len(self.moms)):
            for t in C.TYPES:
                for i in range(C.SIZE[t]):
                    v = _next[0]; _next[0] += 1
                    newvar(v, C.PARITY[t])
                    self.var[(t, i, mi)] = v
        self.inv = {v: k for k, v in self.var.items()}
        self.Qy = {}

    def mode(self, p):
        for i, m in enumerate(self.moms):
            if np.array_equal(m, p):
                return i
        raise KeyError(p)

    def y(self, t, idx, m):
        """variable of component idx (tuple in SHAPE[t]) of field t at mode m"""
        flat = int(np.ravel_multi_index(idx, C.SHAPE[t])) if isinstance(idx, tuple) else idx
        return self.var[(t, flat, m)]

    # ------------------------------------------------------------ monomial accumulation
    @staticmethod
    def add(acc, coef, vs):
        coef = q2(coef)
        if coef.iszero():
            return
        key, sg = (), 1
        for v in vs:
            key, s2 = _mul_keys(key, (v,))
            if key is None:
                return
            sg *= s2
        acc[key] = acc[key] + coef * sg if key in acc else coef * sg

    @staticmethod
    def ga(acc):
        return GA({k: v for k, v in acc.items() if not v.iszero()})

    # ------------------------------------------------------------ local functionals (sum over mode assignments)
    def cubic(self, name, triples):
        """cubic terms of S (I_3) and the basis functionals G1..G5; triples = list of (m1,m2,m3) mode assignments
        to the slots in the printed order"""
        acc = {}
        Fn = np.argwhere(C.F != 0)
        D = C.D
        for (m1, m2, m3) in triples:
            p1, p2, p3 = self.moms[m1], self.moms[m2], self.moms[m3]
            if (p1 + p2 + p3).any():
                continue
            for (a, b, c) in Fn:
                f = C.F[a, b, c]
                if name == 'AAA':            # -1/2 sum_{M,N} (d_M A_N, [A_M, A_N]) : slots (dA_N^c), A_M^a, A_N^b
                    for M in range(4):
                        if not p1[M]:
                            continue
                        for N in range(10):
                            self.add(acc, Q2(0, Fraction(-1, 2) * int(p1[M])) * q2(f),
                                     [self.y('A', (N, c), m1), self.y('A', (M, a), m2), self.y('A', (N, b), m3)])
                elif name == 'llA':          # 1/2 (lambda, rho([A, lambda])) : slots lambda^c, A_M^a, lambda^b
                    for M in range(10):
                        nz = np.argwhere(SB[M] != 0)
                        for (A_, B_) in nz:
                            self.add(acc, q2(SB[M][A_, B_]) * q2(f) * Fraction(1, 2),
                                     [self.y('l', (A_, c), m1), self.y('A', (M, a), m2), self.y('l', (B_, b), m3)])
                elif name == 'AcAs':         # ([A, c], A^*) : A_M^a, c^b, A^*_M^c
                    for M in range(10):
                        self.add(acc, q2(f), [self.y('A', (M, a), m1), self.y('c', (b,), m2), self.y('As', (M, c), m3)])
                elif name == 'cllS':         # ([c, lambda], lambda^*) : c^a, lambda_A^b, lambda^*_A^c
                    for A_ in range(16):
                        self.add(acc, q2(f), [self.y('c', (a,), m1), self.y('l', (A_, b), m2), self.y('ls', (A_, c), m3)])
                elif name == 'cccS':         # 1/2 ([c, c], c^*)
                    self.add(acc, q2(f) * Fraction(1, 2), [self.y('c', (a,), m1), self.y('c', (b,), m2), self.y('cs', (c,), m3)])
                elif name == 'G1':           # ([c, A|_V], d^*d A|_V) ; (d^*dA)_nu = p^2 A_nu - p_nu (p.A)
                    for nu in range(4):
                        self.add(acc, q2(f) * int(p3 @ p3), [self.y('c', (a,), m1), self.y('A', (nu, b), m2), self.y('A', (nu, c), m3)])
                        for mu in range(4):
                            self.add(acc, q2(f) * (-int(p3[nu] * p3[mu])),
                                     [self.y('c', (a,), m1), self.y('A', (nu, b), m2), self.y('A', (mu, c), m3)])
                elif name == 'G2':           # ([c, A|_V], d d^* A|_V) ; (dd^*A)_nu = p_nu (p.A)
                    for nu in range(4):
                        for mu in range(4):
                            self.add(acc, q2(f) * int(p3[nu] * p3[mu]),
                                     [self.y('c', (a,), m1), self.y('A', (nu, b), m2), self.y('A', (mu, c), m3)])
                elif name == 'G3':           # ([c, A|_Wperp], d^*d A|_Wperp)
                    for w in range(4, 10):
                        self.add(acc, q2(f) * int(p3 @ p3), [self.y('c', (a,), m1), self.y('A', (w, b), m2), self.y('A', (w, c), m3)])
                elif name == 'G4':           # ([c, dc], A^*) : c^a, (dc)_mu^b = i p2_mu c^b, A^*_mu^c
                    for mu in range(4):
                        self.add(acc, Q2(0, int(p2[mu])) * q2(f), [self.y('c', (a,), m1), self.y('c', (b,), m2), self.y('As', (mu, c), m3)])
                elif name == 'G5':           # ([c, lambda], dslash lambda) : c^a, lambda_A^b, (i rho(p3) lambda)_A^c
                    for mu in range(4):
                        if not p3[mu]:
                            continue
                        nz = np.argwhere(SB[mu] != 0)
                        for (A_, B_) in nz:
                            self.add(acc, Q2(0, int(p3[mu])) * q2(SB[mu][A_, B_]) * q2(f),
                                     [self.y('c', (a,), m1), self.y('l', (A_, b), m2), self.y('l', (B_, c), m3)])
                else:
                    raise ValueError(name)
        return self.ga(acc)

    def quad(self, name, pairs):
        """quadratic local functionals; pairs = list of (m1, m2) mode assignments to the two printed slots"""
        acc = {}
        G = C.G; D = C.D
        gnz = np.argwhere(G != 0)
        for (m1, m2) in pairs:
            p1, p2 = self.moms[m1], self.moms[m2]
            if (p1 + p2).any():
                continue
            for (a, b) in gnz:
                g = G[a, b]
                if name == 'AA':             # int (A, A), all ten components
                    for M in range(10):
                        self.add(acc, q2(g), [self.y('A', (M, a), m1), self.y('A', (M, b), m2)])
                elif name == 'dAdA':         # int (dA, dA) = 1/2 sum_{MN} F_MN F_MN, F_MN = d_M A_N - d_N A_M
                    for M in range(10):
                        for N in range(10):
                            for (s1, X1, Y1) in ((1, M, N), (-1, N, M)):
                                if X1 >= 4 or not p1[X1]:
                                    continue
                                for (s2, X2, Y2) in ((1, M, N), (-1, N, M)):
                                    if X2 >= 4 or not p2[X2]:
                                        continue
                                    c_ = Fraction(1, 2) * s1 * s2 * (-int(p1[X1]) * int(p2[X2]))   # (i p1)(i p2)
                                    self.add(acc, q2(g) * c_, [self.y('A', (Y1, a), m1), self.y('A', (Y2, b), m2)])
                elif name == 'ldl':          # int (lambda, dslash lambda) = lambda^T rho(e_mu) d_mu lambda
                    for mu in range(4):
                        if not p2[mu]:
                            continue
                        for (A_, B_) in np.argwhere(SB[mu] != 0):
                            self.add(acc, Q2(0, int(p2[mu])) * q2(SB[mu][A_, B_]) * q2(g),
                                     [self.y('l', (A_, a), m1), self.y('l', (B_, b), m2)])
                elif name == 'dcAs':         # int (dc, A^*) (A^*|_V)
                    for mu in range(4):
                        self.add(acc, Q2(0, int(p1[mu])) * q2(g), [self.y('c', (a,), m1), self.y('As', (mu, b), m2)])
                elif name == 'dcAV':         # int (dc, A|_V)
                    for mu in range(4):
                        self.add(acc, Q2(0, int(p1[mu])) * q2(g), [self.y('c', (a,), m1), self.y('A', (mu, b), m2)])
                elif name == 'ddcdA':        # int (d^*dc, d^*A|_V) ; d^*dc -> p^2 c, d^*A -> -i p.A
                    for mu in range(4):
                        self.add(acc, Q2(0, -int(p2[mu]) * int(p1 @ p1)) * q2(g), [self.y('c', (a,), m1), self.y('A', (mu, b), m2)])
                elif name == 'dsAdsA':       # int (d^*A|_V, d^*A|_V), d^*A -> -i p.A
                    for mu in range(4):
                        for nu in range(4):
                            self.add(acc, q2(g) * (-int(p1[mu]) * int(p2[nu])), [self.y('A', (mu, a), m1), self.y('A', (nu, b), m2)])
                elif name == 'AAs':          # int (A, A^*) all ten components
                    for M in range(10):
                        self.add(acc, q2(g), [self.y('A', (M, a), m1), self.y('As', (M, b), m2)])
                elif name == 'lsl':          # int (lambda^*, lambda)
                    for A_ in range(16):
                        self.add(acc, q2(g), [self.y('ls', (A_, a), m1), self.y('l', (A_, b), m2)])
                elif name == 'csc':          # int (c^*, c)
                    self.add(acc, q2(g), [self.y('cs', (a,), m1), self.y('c', (b,), m2)])
                else:
                    raise ValueError(name)
        return self.ga(acc)

    def S_free(self, pairs):
        return self.quad('dAdA', pairs) * Fraction(-1, 4) + self.quad('ldl', pairs) * Fraction(1, 2) + self.quad('dcAs', pairs)

    # ------------------------------------------------------------ kernels and brackets
    def Delta0(self, F, modes_pm):
        """Delta_0 = sum (-,-)^{-1} d_x d_xi over (x at m, xi at -m)"""
        out = GA()
        D = C.D
        for (m, mm_) in modes_pm:
            for (t, ts, n) in (('A', 'As', 10), ('ls', 'l', 16), ('cs', 'c', 1)):
                for k in range(n):
                    for a in range(D):
                        for b in range(D):
                            gi = C.GINV[a, b]
                            if gi == 0:
                                continue
                            ix = (k, a) if t != 'cs' else (a,)
                            iy = (k, b) if ts != 'c' else (b,)
                            d = F.dleft(self.y(ts, iy, mm_))
                            if not d.t:
                                continue
                            out = out + d.dleft(self.y(t, ix, m)) * q2(gi)
        return out

    def bracket0(self, f, g, pf, modes_pm):
        return self.Delta0(f * g, modes_pm) - self.Delta0(f, modes_pm) * g - (f * self.Delta0(g, modes_pm)) * ((-1) ** pf)

    def kernel(self, which, m, mneg, scalar=Q2(1)):
        """kernel dict between mode m (momentum p) and mode mneg (momentum -p): 'K' (Khat), 'P' (Pnum)"""
        p = self.moms[m]
        assert not (p + self.moms[mneg]).any()
        Kd = {}
        if which == 'K':
            M = sum(C.Khat(L) for L in C.SLINES)
            Mneg = M
        else:
            M = sum(C.Pnum(p, L) for L in C.PLINES)
            Mneg = sum(C.Pnum(-p, L) for L in C.PLINES)
        def var_of(gidx, mode):
            for t in C.TYPES:
                if C.OFF[t] <= gidx < C.OFF[t] + C.SIZE[t]:
                    return self.var[(t, gidx - C.OFF[t], mode)]
        for (u, v) in np.argwhere(M != 0):
            Kd[(var_of(u, m), var_of(v, mneg))] = q2(M[u, v]) * scalar
        for (u, v) in np.argwhere(Mneg != 0):
            Kd[(var_of(u, mneg), var_of(v, m))] = q2(Mneg[u, v]) * scalar
        return Kd

    @staticmethod
    def dX(F, X):
        """d_X F = (1/2) sum X_uv d_v d_u F"""
        out = GA()
        firsts = {}
        for (u, v), c in X.items():
            firsts.setdefault(u, []).append((v, c))
        for u, lst in firsts.items():
            Fu = F.dleft(u)
            if not Fu.t:
                continue
            for v, c in lst:
                Fuv = Fu.dleft(v)
                if Fuv.t:
                    out = out + Fuv * (c * Fraction(1, 2))
        return out

    def C(self, X, f, g):
        return self.dX(f * g, X) - self.dX(f, X) * g - f * self.dX(g, X)

    # ------------------------------------------------------------ Q as a derivation
    def Q_var(self, v, Sfree, modes_pm):
        """Q y_v = {S_free, y_v} = Delta_0(S_free y_v) (Delta_0 S_free = 0 = Delta_0 y_v): only the canonical pairs
        containing v contribute"""
        if v not in self.Qy:
            t, i, mo = self.inv[v]
            D = C.D
            pm = dict(modes_pm)
            P = {'A': 'As', 'As': 'A', 'ls': 'l', 'l': 'ls', 'cs': 'c', 'c': 'cs'}
            even = {'A': True, 'As': False, 'ls': True, 'l': False, 'cs': True, 'c': False}
            idx = np.unravel_index(i, C.SHAPE[t])
            a = idx[-1]; rest = idx[:-1]
            partner = P[t]; mo2 = pm[mo]
            Sy = Sfree * GA.var(v)
            out = GA()
            for b in range(D):
                gi = C.GINV[a, b] if even[t] else C.GINV[b, a]
                if gi == 0:
                    continue
                w = self.y(partner, tuple(rest) + (b,), mo2)
                xv, xiv = (v, w) if even[t] else (w, v)
                d = Sy.dleft(xiv)
                if d.t:
                    out = out + d.dleft(xv) * q2(gi)
            self.Qy[v] = out
        return self.Qy[v]

    def Q_var_full(self, v, Sfree, modes_pm):
        return self.bracket0(Sfree, GA.var(v), 0, modes_pm)

    def Q_apply(self, F, Sfree, modes_pm, only_modes=None):
        out = {}
        for key, c in F.t.items():
            nodd = 0
            for i, v in enumerate(key):
                t, _, mo = self.inv.get(v, (None, None, None))
                if t is not None and (only_modes is None or mo in only_modes):
                    qv = self.Q_var(v, Sfree, modes_pm)
                    if qv.t:
                        pre = GA({key[:i]: Q2(1)}); post = GA({key[i + 1:]: Q2(1)})
                        term = pre * qv * post * (c * ((-1) ** nodd))
                        for k2, c2 in term.t.items():
                            out[k2] = out[k2] + c2 if k2 in out else c2
                nodd += PAR[v]
        return GA({k: v for k, v in out.items() if not v.iszero()})

    # ------------------------------------------------------------ external legs
    def subst(self, F, legs, modes):
        """legs: list of (type, pol, mode, rvar or None); variables of these modes are replaced by the leg
        (other field types at that mode -> 0); variables of other modes stay"""
        sub = {}
        for (t, pol, m, r) in legs:
            for tt in C.TYPES:
                for i in range(C.SIZE[tt]):
                    v = self.var[(tt, i, m)]
                    if tt != t:
                        sub[v] = None
                    else:
                        val = complex(np.ravel(pol)[i])
                        if val == 0:
                            sub[v] = None
                        else:
                            sub[v] = GA.var(r, q2(val)) if r is not None else GA.const(q2(val))
        out = GA()
        acc = {}
        for key, c in F.t.items():
            term = GA({(): c})
            dead = False
            rest = []
            for v in key:
                if v in sub:
                    if sub[v] is None:
                        dead = True; break
                    term = term * sub[v]
                else:
                    term = term * GA.var(v)
            if dead:
                continue
            for k2, c2 in term.t.items():
                acc[k2] = acc[k2] + c2 if k2 in acc else c2
        return GA({k: v for k, v in acc.items() if not v.iszero()})


def rvars(n):
    out = []
    for j in range(n):
        v = _next[0]; _next[0] += 1
        newvar(v, 1)
        out.append(v)
    return out


def coeff(F, rs):
    key = tuple(sorted(rs))
    return F.t.get(key, Q2(0))
