"""Minimal bigraded supercommutative polynomial algebra with an odd BV operator.

Variables carry (parity, degree).  Monomials: (tuple of exponents of even variables,
sorted tuple of odd variable indices).  Coefficients: Fraction.
Delta = sum over (field, antifield) pairs of left derivative d_field d_antifield.
Bracket: {f,g} = Delta(fg) - Delta(f) g - (-1)^{|f|} f Delta(g)   (ch01 def:c1-bracket),
|f| = total parity.
"""
from fractions import Fraction
import itertools, random


class Alg:
    def __init__(self, names, parities, degrees, pairs):
        self.names = list(names)
        self.par = list(parities)
        self.deg = list(degrees)
        self.even = [i for i, p in enumerate(self.par) if p == 0]
        self.odd = [i for i, p in enumerate(self.par) if p == 1]
        self.eidx = {v: k for k, v in enumerate(self.even)}
        self.pairs = list(pairs)  # list of (field_index, antifield_index)

    # element: dict {(evexp_tuple, oddtuple): Fraction}
    def zero(self):
        return {}

    def one(self):
        return {(tuple([0] * len(self.even)), ()): Fraction(1)}

    def var(self, i):
        if self.par[i] == 0:
            e = [0] * len(self.even)
            e[self.eidx[i]] = 1
            return {(tuple(e), ()): Fraction(1)}
        return {(tuple([0] * len(self.even)), (i,)): Fraction(1)}

    @staticmethod
    def add(a, b, cb=1):
        r = dict(a)
        for k, v in b.items():
            r[k] = r.get(k, 0) + cb * v
            if r[k] == 0:
                del r[k]
        return r

    @staticmethod
    def scale(a, c):
        if c == 0:
            return {}
        return {k: v * c for k, v in a.items()}

    @staticmethod
    def merge_odd(o1, o2):
        # sign of concatenation o1+o2 sorted
        if set(o1) & set(o2):
            return None, 0
        seq = list(o1) + list(o2)
        # count inversions
        inv = 0
        for i in range(len(seq)):
            for j in range(i + 1, len(seq)):
                if seq[i] > seq[j]:
                    inv += 1
        return tuple(sorted(seq)), (-1) ** inv

    def mul(self, a, b):
        r = {}
        for (e1, o1), c1 in a.items():
            for (e2, o2), c2 in b.items():
                o, s = self.merge_odd(o1, o2)
                if o is None:
                    continue
                e = tuple(x + y for x, y in zip(e1, e2))
                key = (e, o)
                r[key] = r.get(key, 0) + s * c1 * c2
                if r[key] == 0:
                    del r[key]
        return r

    def dleft(self, a, i):
        """left derivative with respect to variable i"""
        r = {}
        if self.par[i] == 0:
            k = self.eidx[i]
            for (e, o), c in a.items():
                if e[k] == 0:
                    continue
                e2 = list(e)
                e2[k] -= 1
                key = (tuple(e2), o)
                r[key] = r.get(key, 0) + c * e[k]
                if r[key] == 0:
                    del r[key]
        else:
            for (e, o), c in a.items():
                if i not in o:
                    continue
                pos = o.index(i)
                o2 = tuple(x for x in o if x != i)
                key = (e, o2)
                r[key] = r.get(key, 0) + c * (-1) ** pos
                if r[key] == 0:
                    del r[key]
        return r

    def Delta(self, a):
        r = {}
        for (f, fs) in self.pairs:
            r = self.add(r, self.dleft(self.dleft(a, fs), f))
        return r

    def parity(self, a):
        ps = set()
        for (e, o), c in a.items():
            p = sum(self.par[self.even[k]] * x for k, x in enumerate(e)) + len(o)
            ps.add(p % 2)
        if len(ps) > 1:
            raise ValueError('inhomogeneous parity')
        return ps.pop() if ps else 0

    def degree(self, a):
        ds = set()
        for (e, o), c in a.items():
            d = sum(self.deg[self.even[k]] * x for k, x in enumerate(e)) + sum(self.deg[i] for i in o)
            ds.add(d)
        if len(ds) > 1:
            raise ValueError('inhomogeneous degree')
        return ds.pop() if ds else 0

    def br(self, f, g):
        pf = self.parity(f)
        r = self.Delta(self.mul(f, g))
        r = self.add(r, self.mul(self.Delta(f), g), -1)
        r = self.add(r, self.mul(f, self.Delta(g)), -((-1) ** pf))
        return r

    def is_zero(self, a):
        return len(a) == 0

    def random_hom(self, rng, parity, degree, nterms=4, maxdeg=3):
        """random homogeneous element of given parity and degree"""
        out = {}
        tries = 0
        while len(out) < nterms and tries < 5000:
            tries += 1
            e = [rng.randint(0, 1) for _ in self.even]
            o = tuple(sorted(rng.sample(self.odd, rng.randint(0, min(maxdeg, len(self.odd))))))
            if sum(e) + len(o) > maxdeg or sum(e) + len(o) == 0:
                continue
            p = (sum(self.par[self.even[k]] * x for k, x in enumerate(e)) + len(o)) % 2
            d = sum(self.deg[self.even[k]] * x for k, x in enumerate(e)) + sum(self.deg[i] for i in o)
            if p != parity or d != degree:
                continue
            out[(tuple(e), o)] = Fraction(rng.randint(-5, 5) or 1)
        return out
