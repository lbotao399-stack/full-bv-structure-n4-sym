"""Minimal exact Grassmann algebra (third, independent implementation; imports nothing from engine_a/engine_b).
Element = dict {sorted tuple of odd generator ids: sympy coefficient}.  Even variables are sympy symbols inside
the coefficients."""
import sympy as sp


def _perm_sign(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


class G:
    __slots__ = ('d',)

    def __init__(self, d=None):
        self.d = {}
        if d:
            for k, v in d.items():
                v = sp.expand(v)
                if v != 0:
                    self.d[k] = v

    @staticmethod
    def c(x):
        return G({(): x}) if x != 0 else G()

    @staticmethod
    def gen(i, coef=1):
        return G({(i,): coef})

    def __add__(self, o):
        if not isinstance(o, G):
            o = G.c(o)
        out = dict(self.d)
        for k, v in o.d.items():
            out[k] = out.get(k, 0) + v
        return G(out)

    __radd__ = __add__

    def __neg__(self):
        return G({k: -v for k, v in self.d.items()})

    def __sub__(self, o):
        if not isinstance(o, G):
            o = G.c(o)
        return self + (-o)

    def __mul__(self, o):
        if not isinstance(o, G):
            return G({k: v * o for k, v in self.d.items()})
        out = {}
        for k1, v1 in self.d.items():
            s1 = set(k1)
            for k2, v2 in o.d.items():
                if s1.intersection(k2):
                    continue
                m = k1 + k2
                key = tuple(sorted(m))
                out[key] = out.get(key, 0) + _perm_sign(m) * v1 * v2
        return G(out)

    def __rmul__(self, o):  # scalar * G
        return G({k: v * o for k, v in self.d.items()})

    def dodd(self, i):
        """left derivative with respect to the odd generator i"""
        out = {}
        for k, v in self.d.items():
            if i in k:
                pos = k.index(i)
                key = k[:pos] + k[pos + 1:]
                out[key] = out.get(key, 0) + (-1) ** pos * v
        return G(out)

    def deven(self, s):
        return G({k: sp.diff(v, s) for k, v in self.d.items()})

    def subs(self, sub):
        return G({k: v.subs(sub) for k, v in self.d.items()})

    def coeff(self, key):
        return self.d.get(tuple(sorted(key)), sp.Integer(0))

    def iszero(self):
        return not self.d
