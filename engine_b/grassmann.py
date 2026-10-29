"""Minimal supercommutative polynomial algebra (engine B, sign bookkeeping only).
Elements: dict {tuple(sorted odd generator ids): sympy coefficient (may contain even symbols)}.
Left derivatives.  Used for the 0-dim verification of the one-loop supertrace formula."""
import sympy as sp

def _merge_sign(a, b):
    # sign of reordering the concatenation a+b (each sorted) into sorted order
    inv = 0
    for x in a:
        for y in b:
            if x > y:
                inv += 1
            elif x == y:
                return 0
    return -1 if inv % 2 else 1

class G:
    __slots__ = ("t",)
    def __init__(self, t=None):
        self.t = {} if t is None else {k: v for k, v in t.items() if v != 0}
    @staticmethod
    def scalar(c):
        return G({(): sp.sympify(c)})
    @staticmethod
    def odd(i):
        return G({(i,): sp.Integer(1)})
    def __add__(self, o):
        if not isinstance(o, G):
            o = G.scalar(o)
        r = dict(self.t)
        for k, v in o.t.items():
            r[k] = sp.expand(r.get(k, 0) + v)
        return G(r)
    __radd__ = __add__
    def __neg__(self):
        return G({k: -v for k, v in self.t.items()})
    def __sub__(self, o):
        return self + (-o if isinstance(o, G) else G.scalar(-o))
    def __mul__(self, o):
        if not isinstance(o, G):
            return G({k: sp.expand(v * o) for k, v in self.t.items()})
        r = {}
        for ka, va in self.t.items():
            for kb, vb in o.t.items():
                s = _merge_sign(ka, kb)
                if s == 0:
                    continue
                k = tuple(sorted(ka + kb))
                r[k] = r.get(k, 0) + s * va * vb
        return G({k: sp.expand(v) for k, v in r.items()})
    def __rmul__(self, c):
        return G({k: sp.expand(c * v) for k, v in self.t.items()})
    def d_odd(self, i):
        r = {}
        for k, v in self.t.items():
            if i in k:
                m = k.index(i)
                kk = k[:m] + k[m + 1:]
                r[kk] = r.get(kk, 0) + (-1) ** m * v
        return G(r)
    def d_even(self, x):
        return G({k: sp.expand(sp.diff(v, x)) for k, v in self.t.items()})
    def subs(self, d):
        return G({k: sp.expand(v.subs(d)) for k, v in self.t.items()})
    def drop_odd(self, ids):
        ids = set(ids)
        return G({k: v for k, v in self.t.items() if not (set(k) & ids)})
    def map(self, f):
        return G({k: f(v) for k, v in self.t.items()})
    def iszero(self):
        return all(sp.simplify(v) == 0 for v in self.t.values())
    def __repr__(self):
        return repr(self.t)
