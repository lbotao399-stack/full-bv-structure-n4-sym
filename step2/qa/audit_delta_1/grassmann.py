# Independent finite-dimensional BV toy for auditing Part D, D1 (a),(b),(c).
# Written from the manuscript conventions only (ch01 Def c1-pairing, c1-bracket, Prop c1-Q;
# ch03 Def c3-delta, c3-W; [C07 §6.2-6.4] conventions as quoted in the ledger).  Does NOT import delta_toy.py.
#
# Graded-commutative polynomial algebra over Q in:
#   index 0      : xi   (odd parameter = ghost of the dilatation; not a field)
#   index 1..NF  : field coordinates (parity from the slot table)
# plus three even parameters carried as integer exponents: hbar (may be negative), g (>=0, truncated), eps (any int).
# A monomial is a tuple of length NV+3; odd variables have exponent 0/1; canonical order = increasing index.
# Left derivatives throughout:  d/dx_i (x_i m) = m.

from fractions import Fraction as Fr
import itertools, random

class Alg:
    def __init__(self, parities, gmax):
        self.par = list(parities)          # parity of each Grassmann/even variable (index 0 = xi)
        self.NV = len(parities)
        self.gmax = gmax
        self.H, self.G, self.E = self.NV, self.NV + 1, self.NV + 2   # positions of hbar, g, eps exponents

    # ---------- basic arithmetic ----------
    def zero(self):
        return {}

    def const(self, c, h=0, g=0, e=0):
        k = [0] * (self.NV + 3)
        k[self.H], k[self.G], k[self.E] = h, g, e
        return {tuple(k): Fr(c)} if c != 0 else {}

    def var(self, i, c=1):
        k = [0] * (self.NV + 3)
        k[i] = 1
        return {tuple(k): Fr(c)}

    def add(self, *ps):
        r = {}
        for p in ps:
            for k, v in p.items():
                nv = r.get(k, 0) + v
                if nv == 0:
                    r.pop(k, None)
                else:
                    r[k] = nv
        return r

    def scale(self, p, c, h=0, g=0, e=0):
        c = Fr(c)
        if c == 0:
            return {}
        r = {}
        for k, v in p.items():
            kk = list(k)
            kk[self.H] += h
            kk[self.G] += g
            kk[self.E] += e
            if kk[self.G] > self.gmax:
                continue
            r[tuple(kk)] = v * c
        return r

    def sub(self, p, q):
        return self.add(p, self.scale(q, -1))

    def _mulmono(self, a, b):
        par = self.par
        sign = 1
        # sign from moving odd vars of b left past odd vars of a with larger index
        oddb = [i for i in range(self.NV) if b[i] and par[i]]
        if oddb:
            odda = [i for i in range(self.NV) if a[i] and par[i]]
            for j in oddb:
                if a[j]:
                    return None, 0
            cnt = 0
            for j in oddb:
                for i in odda:
                    if i > j:
                        cnt += 1
            if cnt % 2:
                sign = -1
        k = tuple(x + y for x, y in zip(a, b))
        if k[self.G] > self.gmax:
            return None, 0
        return k, sign

    def mul(self, p, q):
        r = {}
        for a, va in p.items():
            for b, vb in q.items():
                k, s = self._mulmono(a, b)
                if k is None:
                    continue
                nv = r.get(k, 0) + s * va * vb
                if nv == 0:
                    r.pop(k, None)
                else:
                    r[k] = nv
        return r

    def parity(self, p):
        ps = set()
        for k in p:
            ps.add(sum(k[i] for i in range(self.NV) if self.par[i]) % 2)
        assert len(ps) <= 1, "inhomogeneous parity"
        return ps.pop() if ps else 0

    # ---------- derivatives ----------
    def d(self, p, i):
        """left derivative d/dx_i"""
        r = {}
        par = self.par
        for k, v in p.items():
            if not k[i]:
                continue
            kk = list(k)
            if par[i]:
                s = sum(k[j] for j in range(i) if par[j]) % 2
                c = -v if s else v
                kk[i] = 0
            else:
                c = v * k[i]
                kk[i] -= 1
            kk = tuple(kk)
            nv = r.get(kk, 0) + c
            if nv == 0:
                r.pop(kk, None)
            else:
                r[kk] = nv
        return r

    def dvec(self, p, vec):
        """partial along a vector: sum_i vec[i] d/dx_i  (vec: dict field-index -> scalar)"""
        return self.add(*[self.scale(self.d(p, i), c) for i, c in vec.items() if c != 0])

    def dker(self, p, phi):
        """partial_phi = 1/2 sum phi^{ij} d_j d_i   [C07 §6.4], phi: dict (i,j)->scalar"""
        terms = []
        for (i, j), c in phi.items():
            if c == 0:
                continue
            terms.append(self.scale(self.d(self.d(p, i), j), Fr(c, 2)))
        return self.add(*terms)

    def vf(self, p, N):
        """derivation induced by the linear operator N (N e_k = sum_j N[j][k] e_j):
           D_N = sum_j (sum_k N[j][k] x_k) d/dx_j, coefficient placed on the left."""
        terms = []
        if isinstance(N, list):
            N = {j: {k: N[j][k] for k in range(len(N[j])) if N[j][k] != 0} for j in range(len(N))}
        for j, row in N.items():
            lin = self.add(*[self.var(k, c) for k, c in row.items() if c != 0])
            if not lin:
                continue
            terms.append(self.mul(lin, self.d(p, j)))
        return self.add(*terms)

    # ---------- exponentials ----------
    def exp_dker(self, p, phi, hpow=1):
        """e^{hbar partial_phi} p (terminates on polynomials)"""
        r = dict(p)
        term = dict(p)
        n = 0
        while term:
            n += 1
            term = self.scale(self.dker(term, phi), Fr(1, n), h=hpow)
            r = self.add(r, term)
        return r

    def exp(self, X):
        """e^X for X of g-order >= 1 or nilpotent (truncated at gmax)"""
        r = self.const(1)
        term = self.const(1)
        n = 0
        while True:
            n += 1
            term = self.scale(self.mul(term, X), Fr(1, n))
            if not term:
                break
            r = self.add(r, term)
            if n > 60:
                raise RuntimeError("exp did not terminate")
        return r

    def log1p(self, Y):
        r = {}
        term = self.const(1)
        n = 0
        while True:
            n += 1
            term = self.mul(term, Y)
            if not term:
                break
            r = self.add(r, self.scale(term, Fr((-1) ** (n + 1), n)))
            if n > 60:
                raise RuntimeError("log did not terminate")
        return r

    def inv1p(self, Y):
        r = self.const(1)
        term = self.const(1)
        n = 0
        while True:
            n += 1
            term = self.scale(self.mul(term, Y), -1)
            if not term:
                break
            r = self.add(r, term)
            if n > 60:
                raise RuntimeError("inverse did not terminate")
        return r

    # ---------- utilities ----------
    def eq(self, p, q):
        return not self.sub(p, q)

    def min_hbar(self, p):
        return min((k[self.H] for k in p), default=None)

    def eps_deriv(self, p):
        """eps d/deps"""
        r = {}
        for k, v in p.items():
            if k[self.E]:
                r[k] = v * k[self.E]
        return r

    def coeff_left(self, p, i):
        """left coefficient of the odd variable i:  p = x_i * a + b  ->  a"""
        return self.d(p, i)
