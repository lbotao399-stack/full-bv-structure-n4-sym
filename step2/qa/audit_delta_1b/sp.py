"""Exact graded-commutative polynomial algebra for the audit of Part D, D1 (a), (b), (c).

Written independently of step2/theory/delta_toy.py (not imported, not copied).

Generators x^i (coordinate functionals on a graded space E with basis e_i), each with a total
parity p_i (Koszul signs use p_i only) and an integer degree as a functional deg_i = -deg(e_i).
A polynomial is a dict {monomial: Fraction}; a monomial is a sorted tuple of generator indices,
odd generators appear at most once.

Operators:
  mul(f, g)          product (Koszul signs from sorting)
  d(i, f)            left derivative d/dx^i (graded derivation of parity p_i)
  D(X, f, conv)      first-order operator induced by a linear map X of E (X[(i, j)] = X^i_j, X e_j = sum_i X^i_j e_i)
                       conv 'naive': x^i -> sum_j X^i_j x^j
                       conv 'canon': x^i -> (-1)^{p_X p_i} sum_j X^i_j x^j  (canonical graded transpose)
  dphi(phi, f, conv) constant-coefficient second-order operator of phi = sum phi^{ij} e_i (x) e_j
                       conv 'A': (1/2) sum phi^{ij} d_i d_j      conv 'B': (1/2) sum phi^{ij} d_j d_i
"""
from fractions import Fraction as Fr
import itertools

class Alg:
    def __init__(self, par, deg=None):
        self.par = list(par)
        self.n = len(par)
        self.deg = list(deg) if deg is not None else [0] * self.n

    # ---------------- basic algebra
    def sort(self, m):
        m = list(m); s = 1
        # insertion sort counting odd-odd transpositions
        for a in range(1, len(m)):
            b = a
            while b > 0 and m[b - 1] > m[b]:
                if self.par[m[b - 1]] and self.par[m[b]]:
                    s = -s
                m[b - 1], m[b] = m[b], m[b - 1]
                b -= 1
        for a in range(1, len(m)):
            if m[a] == m[a - 1] and self.par[m[a]]:
                return None, 0
        return tuple(m), s

    @staticmethod
    def clean(f):
        return {m: c for m, c in f.items() if c != 0}

    def add(self, *terms):
        """add(f1, f2, ...) or add((a1, f1), (a2, f2), ...)"""
        out = {}
        for t in terms:
            if isinstance(t, tuple):
                a, f = t
            else:
                a, f = 1, t
            for m, c in f.items():
                out[m] = out.get(m, 0) + a * c
        return self.clean(out)

    def scale(self, a, f):
        return self.clean({m: a * c for m, c in f.items()})

    def mul(self, f, g):
        out = {}
        for m1, c1 in f.items():
            for m2, c2 in g.items():
                m, s = self.sort(m1 + m2)
                if m is None:
                    continue
                out[m] = out.get(m, 0) + s * c1 * c2
        return self.clean(out)

    def x(self, i):
        return {(i,): Fr(1)}

    def one(self):
        return {(): Fr(1)}

    def parity_of_mono(self, m):
        return sum(self.par[i] for i in m) % 2

    def parity(self, f):
        ps = {self.parity_of_mono(m) for m in f}
        assert len(ps) <= 1, "inhomogeneous parity"
        return ps.pop() if ps else 0

    def degree_of_mono(self, m):
        return sum(self.deg[i] for i in m)

    def degree(self, f):
        ds = {self.degree_of_mono(m) for m in f}
        assert len(ds) <= 1, "inhomogeneous degree"
        return ds.pop() if ds else None

    def d_mono(self, i, m):
        """left derivative of a monomial: returns (factor, monomial) or None"""
        if i not in m:
            return None
        k = m.index(i)
        mult = m.count(i)
        sgn = -1 if (self.par[i] and sum(self.par[j] for j in m[:k]) % 2) else 1
        return sgn * mult, m[:k] + m[k + 1:]

    def d(self, i, f):
        out = {}
        for m, c in f.items():
            r = self.d_mono(i, m)
            if r is None:
                continue
            s, mm = r
            out[mm] = out.get(mm, 0) + s * c
        return self.clean(out)

    # ---------------- operators built from linear data
    def D(self, X, pX, f, conv):
        """first-order operator sum_{i,j} c^i_j x^j d_i with c from X (see module doc)"""
        rows = {}
        for (i, j), v in X.items():
            if v == 0:
                continue
            if conv == 'canon' and pX and self.par[i]:
                v = -v
            rows.setdefault(i, []).append((j, v))
        out = {}
        for m, c in f.items():
            for i in set(m):
                if i not in rows:
                    continue
                s, mm = self.d_mono(i, m)
                for j, v in rows[i]:
                    m2, s2 = self.sort((j,) + mm)
                    if m2 is None:
                        continue
                    out[m2] = out.get(m2, 0) + s * s2 * v * c
        return self.clean(out)

    def dphi(self, phi, f, conv):
        """(1/2) sum phi^{ij} d_i d_j (conv A) or (1/2) sum phi^{ij} d_j d_i (conv B)"""
        # index the kernel by the derivative applied first
        first = {}
        for (i, j), v in phi.items():
            if v == 0:
                continue
            if conv == 'A':
                first.setdefault(j, []).append((i, v))   # d_i d_j : d_j first
            else:
                first.setdefault(i, []).append((j, v))   # d_j d_i : d_i first
        out = {}
        half = Fr(1, 2)
        for m, c in f.items():
            for b in set(m):
                if b not in first:
                    continue
                s1, m1 = self.d_mono(b, m)
                for a, v in first[b]:
                    r = self.d_mono(a, m1)
                    if r is None:
                        continue
                    s2, m2 = r
                    out[m2] = out.get(m2, 0) + half * v * s1 * s2 * c
        return self.clean(out)

    def C(self, phi, pphi, f, g, conv):
        """dphi(fg) - dphi(f) g - (-1)^{|phi||f|} f dphi(g)"""
        pf = self.parity(f)
        s = -1 if (pphi and pf) else 1
        return self.add(self.dphi(phi, self.mul(f, g), conv),
                        (-1, self.mul(self.dphi(phi, f, conv), g)),
                        (-s, self.mul(f, self.dphi(phi, g, conv))))

    def exp_dphi(self, phi, hb, f, conv, nmax=40):
        """e^{hb d_phi} f (finite on polynomials)"""
        out = dict(f); term = dict(f)
        for k in range(1, nmax):
            term = self.scale(Fr(hb) / k, self.dphi(phi, term, conv))
            if not term:
                break
            out = self.add(out, term)
        return out

    def dist(self, f, g):
        h = self.add(f, (-1, g))
        return max([abs(c) for c in h.values()] + [Fr(0)])

    def maxabs(self, f):
        return max([abs(c) for c in f.values()] + [Fr(0)])

    def truncate(self, f, maxdeg):
        return {m: c for m, c in f.items() if len(m) <= maxdeg}


# ---------------- sparse matrices (dict (i,j) -> Fraction)
def mat_mul(A, B):
    Bi = {}
    for (k, j), v in B.items():
        Bi.setdefault(k, []).append((j, v))
    out = {}
    for (i, k), a in A.items():
        for j, b in Bi.get(k, []):
            out[(i, j)] = out.get((i, j), 0) + a * b
    return {k: v for k, v in out.items() if v != 0}

def mat_add(*terms):
    out = {}
    for t in terms:
        a, A = t if isinstance(t, tuple) else (1, t)
        for k, v in A.items():
            out[k] = out.get(k, 0) + a * v
    return {k: v for k, v in out.items() if v != 0}

def mat_scale(a, A):
    return {k: a * v for k, v in A.items() if a * v != 0}

def mat_id(n):
    return {(i, i): Fr(1) for i in range(n)}

def mat_comm(A, B, pA, pB):
    """graded commutator AB - (-1)^{pA pB} BA"""
    s = -1 if (pA and pB) else 1
    return mat_add(mat_mul(A, B), (-s, mat_mul(B, A)))

def mat_maxabs(A):
    return max([abs(v) for v in A.values()] + [Fr(0)])

def mat_inv(A, n):
    """exact inverse of an n x n sparse matrix (Gauss-Jordan)"""
    M = [[A.get((i, j), Fr(0)) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [v / pv for v in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                fct = M[r][c]
                M[r] = [a - fct * b for a, b in zip(M[r], M[c])]
    return {(i, j): M[i][n + j] for i in range(n) for j in range(n) if M[i][n + j] != 0}

def tensor_X1(X, phi):
    """(X (x) 1) phi : components sum_i X^k_i phi^{ij}"""
    return mat_mul(X, phi)

def tensor_1X(X, pX, phi, par):
    """(1 (x) X) phi : components (-1)^{pX p_i} sum_j phi^{ij} X^k_j  at (i,k)"""
    Xi = {}
    for (k, j), v in X.items():
        Xi.setdefault(j, []).append((k, v))
    out = {}
    for (i, j), a in phi.items():
        s = -1 if (pX and par[i]) else 1
        for k, b in Xi.get(j, []):
            out[(i, k)] = out.get((i, k), 0) + s * a * b
    return {k: v for k, v in out.items() if v != 0}

def tensor_act(X, pX, phi, par):
    """(X (x) 1 + 1 (x) X) phi"""
    return mat_add(tensor_X1(X, phi), tensor_1X(X, pX, phi, par))

def graded_swap(phi, par):
    """tau(phi): (tau phi)^{ji} = (-1)^{p_i p_j} phi^{ij}"""
    return {(j, i): (-v if (par[i] and par[j]) else v) for (i, j), v in phi.items()}
