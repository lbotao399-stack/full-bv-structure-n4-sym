"""Minimal supercommutative polynomial algebra (exact Fractions or complex).

Generators are integers with a parity (0 even, 1 odd). A monomial is a tuple of
(gen, power) sorted by gen; odd generators have power 1. An element is a dict
monomial -> coefficient. Products follow the Koszul rule (odd generators anticommute).
Written from scratch for refute_cubic_1 (no project code imported).
"""
from fractions import Fraction


class Alg:
    def __init__(self):
        self.par = {}

    def gen(self, g, parity):
        self.par[g] = parity
        return {((g, 1),): 1}

    # ---- monomial product with sign
    def mono_mul(self, m1, m2):
        # sign: number of pairs (a in m1 odd, b in m2 odd) with a > b
        odd2 = [g for g, p in m2 if self.par[g]]
        sign = 1
        d = dict(m1)
        for g, p in m1:
            if self.par[g]:
                if g in [x for x, _ in m2]:
                    return None, 0
                cnt = sum(1 for b in odd2 if b < g)
                if cnt % 2:
                    sign = -sign
        for g, p in m2:
            d[g] = d.get(g, 0) + p
        return tuple(sorted(d.items())), sign

    def mul(self, A, B, trunc=None):
        out = {}
        for m1, c1 in A.items():
            for m2, c2 in B.items():
                m, s = self.mono_mul(m1, m2)
                if m is None:
                    continue
                if trunc is not None and not trunc(m):
                    continue
                out[m] = out.get(m, 0) + s * c1 * c2
        return {m: c for m, c in out.items() if c != 0}

    def add(self, A, B, sb=1):
        out = dict(A)
        for m, c in B.items():
            out[m] = out.get(m, 0) + sb * c
        return {m: c for m, c in out.items() if c != 0}

    def scal(self, A, s):
        return {m: s * c for m, c in A.items() if s * c != 0}

    def const(self, c):
        return {(): c} if c != 0 else {}

    # ---- left derivative along generator g
    def dleft(self, A, g):
        out = {}
        for m, c in A.items():
            d = dict(m)
            if g not in d:
                continue
            if self.par[g]:
                cnt = sum(1 for x, _ in m if self.par[x] and x < g)
                s = -1 if cnt % 2 else 1
                del d[g]
                mm = tuple(sorted(d.items()))
                out[mm] = out.get(mm, 0) + s * c
            else:
                p = d[g]
                if p == 1:
                    del d[g]
                else:
                    d[g] = p - 1
                mm = tuple(sorted(d.items()))
                out[mm] = out.get(mm, 0) + p * c
        return {m: c for m, c in out.items() if c != 0}

    # ---- substitution homomorphism: gen -> element (parity preserving)
    def subst(self, A, sub, trunc=None):
        out = {}
        for m, c in A.items():
            term = self.const(c)
            for g, p in m:
                img = sub.get(g, {((g, 1),): 1})
                for _ in range(p):
                    term = self.mul(term, img, trunc)
            out = self.add(out, term)
        return out

    def degree_in(self, m, gens):
        return sum(p for g, p in m if g in gens)
