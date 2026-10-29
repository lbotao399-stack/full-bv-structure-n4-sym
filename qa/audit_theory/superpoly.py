"""Minimal super-polynomial algebra (independent of every engine in this project).
Generators y_0..y_{n-1} with parities par[i] in {0,1}.  A polynomial is a dict
  key = (even_exponents tuple, sorted tuple of odd generator indices) -> sympy coefficient.
Left derivatives, graded products, odd vector fields.  Exact sympy arithmetic."""
import sympy as sp
from itertools import product as iproduct


class Alg:
    def __init__(self, par):
        self.par = list(par)
        self.n = len(par)
        self.ev = [i for i in range(self.n) if par[i] == 0]
        self.od = [i for i in range(self.n) if par[i] == 1]
        self.evpos = {i: k for k, i in enumerate(self.ev)}

    # ---------- construction ----------
    def one(self):
        return {(tuple([0] * len(self.ev)), ()): sp.Integer(1)}

    def gen(self, i):
        if self.par[i] == 0:
            e = [0] * len(self.ev); e[self.evpos[i]] = 1
            return {(tuple(e), ()): sp.Integer(1)}
        return {(tuple([0] * len(self.ev)), (i,)): sp.Integer(1)}

    @staticmethod
    def clean(p):
        out = {}
        for k, v in p.items():
            v = sp.expand(v)
            if v != 0:
                out[k] = v
        return out

    def add(self, *ps):
        out = {}
        for p in ps:
            for k, v in p.items():
                out[k] = out.get(k, 0) + v
        return self.clean(out)

    def scale(self, c, p):
        return self.clean({k: c * v for k, v in p.items()})

    @staticmethod
    def _merge_odd(a, b):
        """sign and merged sorted tuple of a*b (a,b sorted tuples of odd indices)."""
        if set(a) & set(b):
            return 0, None
        lst = list(a) + list(b)
        # bubble sort counting transpositions
        sign = 1
        lst = lst[:]
        for i in range(len(lst)):
            for j in range(len(lst) - 1 - i):
                if lst[j] > lst[j + 1]:
                    lst[j], lst[j + 1] = lst[j + 1], lst[j]
                    sign = -sign
        return sign, tuple(lst)

    def mul(self, p, q):
        out = {}
        for (e1, o1), v1 in p.items():
            for (e2, o2), v2 in q.items():
                s, o = self._merge_odd(o1, o2)
                if s == 0:
                    continue
                e = tuple(x + y for x, y in zip(e1, e2))
                out[(e, o)] = out.get((e, o), 0) + s * v1 * v2
        return self.clean(out)

    def parity(self, p):
        ps = {len(o) % 2 for (e, o) in p}
        if len(ps) > 1:
            raise ValueError('inhomogeneous parity')
        return ps.pop() if ps else 0

    # ---------- left derivative d/dy_i ----------
    def d(self, i, p):
        out = {}
        if self.par[i] == 0:
            k = self.evpos[i]
            for (e, o), v in p.items():
                if e[k] == 0:
                    continue
                e2 = list(e); e2[k] -= 1
                out[(tuple(e2), o)] = out.get((tuple(e2), o), 0) + e[k] * v
        else:
            for (e, o), v in p.items():
                if i not in o:
                    continue
                pos = o.index(i)
                o2 = tuple(x for x in o if x != i)
                out[(e, o2)] = out.get((e, o2), 0) + (-1) ** pos * v
        return self.clean(out)

    # ---------- odd derivation from images of generators ----------
    def vf(self, images, p):
        """Graded derivation V of parity 1 with V(y_i) = images[i]: V = sum_i images[i] * d_i (coefficient on the left)."""
        out = {}
        for i in range(self.n):
            if not images[i]:
                continue
            out = self.add(out, self.mul(images[i], self.d(i, p)))
        return out

    # ---------- order-two contraction partial_phi ----------
    def dphi(self, phi, p):
        """phi = {(i,j): coeff} meaning sum coeff e_i (x) e_j ; partial_phi = 1/2 sum coeff d_j d_i  [C07 sec 6.4]."""
        out = {}
        for (i, j), c in phi.items():
            if c == 0:
                continue
            out = self.add(out, self.scale(sp.Rational(1, 2) * c, self.d(j, self.d(i, p))))
        return out

    def expo(self, hbar, phi, p, maxn=12):
        out = {}
        term = p
        for n in range(maxn + 1):
            if not term:
                break
            out = self.add(out, term)
            term = self.scale(hbar / (n + 1), self.dphi(phi, term))
        return out

    def random_poly(self, rng, deg_max=3, nterms=6, parity=None):
        out = {}
        for _ in range(nterms * 4):
            if len(out) >= nterms:
                break
            mon = self.one()
            k = rng.randint(1, deg_max)
            for _ in range(k):
                mon = self.mul(mon, self.gen(rng.randrange(self.n)))
            if not mon:
                continue
            if parity is not None and self.parity(mon) != parity:
                continue
            c = sp.Rational(rng.randint(-5, 5), rng.randint(1, 3))
            out = self.add(out, self.scale(c, mon))
        return out
