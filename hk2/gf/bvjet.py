"""
bvjet.py -- a small exact engine for local functionals on R^4 in the conventions of Chapter 1:
jets of component fields with Grassmann parity, polynomial coefficients in the position x,
left Euler-Lagrange derivatives, the bracket of Definition c1-bracket and the zero test of Remark c1-method.

Bracket (derived from {f,g} = Delta_0(fg) - Delta_0(f) g - (-1)^{|f|} f Delta_0(g), Delta_0 = -d_{K_0},
d_phi = 1/2 sum d_{phi''} d_{phi'}, d_e f the directional derivative, f(x + r e) = f(x) + r d_e f):
  {F,G} = -1/2 int sum_{(a,b)} kappa_ab [ (-1)^{|b|(|a|+|F|)} (E_a F)(E_b G) + (-1)^{|a||F|} (E_b F)(E_a G) ],
with K_0 = sum kappa_ab delta(x-y) e_a (x) e_b and E_a the left Euler-Lagrange derivative along the coordinate a.
"""
import sympy as sp
from itertools import product

X = sp.symbols('x1:5', real=True)


class Field:
    def __init__(self, name, parity):
        self.name = name
        self.parity = parity


class Model:
    def __init__(self):
        self.fields = []
        self.index = {}
        self.kappa = {}        # (a,b) -> coefficient of e_a (x) e_b in K_0

    def add(self, name, parity):
        self.index[name] = len(self.fields)
        self.fields.append(Field(name, parity))
        return self.index[name]

    def par(self, a):
        return self.fields[a].parity

    def pair(self, a_name, b_name, kab, kba=None):
        a, b = self.index[a_name], self.index[b_name]
        self.kappa[(a, b)] = kab
        self.kappa[(b, a)] = kab if kba is None else kba

    # ---------------------------------------------------------------- superpolynomials
    # term key: (even, odd); even = tuple(sorted((jet, power))), odd = tuple of jets in canonical order
    # jet = (field_index, (a1,a2,a3,a4))
    def jet(self, name, alpha=(0, 0, 0, 0), coeff=1):
        a = self.index[name]
        j = (a, tuple(alpha))
        if self.par(a) == 0:
            return {(((j, 1),), ()): sp.sympify(coeff)}
        return {((), (j,)): sp.sympify(coeff)}

    def const(self, c):
        return {((), ()): sp.sympify(c)} if c != 0 else {}

    @staticmethod
    def add_into(P, key, c):
        if c == 0:
            return
        v = P.get(key, 0) + c
        if v == 0:
            P.pop(key, None)
        else:
            P[key] = v

    def plus(self, *Ps):
        out = {}
        for P in Ps:
            for k, c in P.items():
                self.add_into(out, k, c)
        return out

    def scale(self, P, c):
        c = sp.sympify(c)
        if c == 0:
            return {}
        return {k: v * c for k, v in P.items()}

    @staticmethod
    def _sort_odd(seq):
        """sort a list of odd jets, return (sign, tuple) or (0, None) if a repetition occurs."""
        s = list(seq)
        sign = 1
        n = len(s)
        for i in range(n):
            for j in range(n - 1 - i):
                if s[j] > s[j + 1]:
                    s[j], s[j + 1] = s[j + 1], s[j]
                    sign = -sign
                elif s[j] == s[j + 1]:
                    return 0, None
        for i in range(n - 1):
            if s[i] == s[i + 1]:
                return 0, None
        return sign, tuple(s)

    @staticmethod
    def _merge_even(e1, e2):
        d = dict(e1)
        for j, pw in e2:
            d[j] = d.get(j, 0) + pw
        return tuple(sorted(d.items()))

    def mul(self, P1, P2):
        out = {}
        for (e1, o1), c1 in P1.items():
            for (e2, o2), c2 in P2.items():
                sign, o = self._sort_odd(o1 + o2)
                if sign == 0:
                    continue
                self.add_into(out, (self._merge_even(e1, e2), o), sign * c1 * c2)
        return out

    def parity_of(self, P):
        ps = set(len(o) % 2 for (e, o) in P.keys())
        if not ps:
            return 0
        assert len(ps) == 1, 'inhomogeneous parity'
        return ps.pop()

    # derivatives
    def dleft(self, P, j):
        """left partial derivative with respect to the jet variable j."""
        a = j[0]
        out = {}
        if self.par(a) == 0:
            for (e, o), c in P.items():
                d = dict(e)
                if j in d:
                    pw = d[j]
                    if pw == 1:
                        del d[j]
                    else:
                        d[j] = pw - 1
                    self.add_into(out, (tuple(sorted(d.items())), o), c * pw)
        else:
            for (e, o), c in P.items():
                if j in o:
                    k = o.index(j)
                    self.add_into(out, (e, o[:k] + o[k + 1:]), c * (-1) ** k)
        return out

    def Dtot(self, P, m):
        """total derivative along x_m (even derivation)."""
        out = {}
        for (e, o), c in P.items():
            dc = sp.diff(c, X[m])
            if dc != 0:
                self.add_into(out, (e, o), dc)
            for (j, pw) in e:
                a, al = j
                al2 = list(al)
                al2[m] += 1
                j2 = (a, tuple(al2))
                d = dict(e)
                if pw == 1:
                    del d[j]
                else:
                    d[j] = pw - 1
                d[j2] = d.get(j2, 0) + 1
                self.add_into(out, (tuple(sorted(d.items())), o), c * pw)
            for k, j in enumerate(o):
                a, al = j
                al2 = list(al)
                al2[m] += 1
                j2 = (a, tuple(al2))
                sign, o2 = self._sort_odd(o[:k] + (j2,) + o[k + 1:])
                if sign == 0:
                    continue
                self.add_into(out, (e, o2), c * sign)
        return out

    def jets_of_field(self, P, a):
        js = set()
        for (e, o) in P.keys():
            for (j, pw) in e:
                if j[0] == a:
                    js.add(j)
            for j in o:
                if j[0] == a:
                    js.add(j)
        return js

    def EL(self, P, a):
        """left Euler-Lagrange derivative along the field a."""
        out = {}
        for j in self.jets_of_field(P, a):
            T = self.dleft(P, j)
            for m in range(4):
                for _ in range(j[1][m]):
                    T = self.scale(self.Dtot(T, m), -1)
            out = self.plus(out, T)
        return out

    def simplify(self, P):
        out = {}
        for k, c in P.items():
            c2 = sp.expand(c)
            if c2 != 0:
                out[k] = c2
        return out

    def bracket(self, F, G):
        pF = self.parity_of(F)
        out = {}
        ELF = {a: self.simplify(self.EL(F, a)) for a in range(len(self.fields))}
        ELG = {a: self.simplify(self.EL(G, a)) for a in range(len(self.fields))}
        for (a, b), kab in self.kappa.items():
            pa, pb = self.par(a), self.par(b)
            if ELF[a] and ELG[b]:
                out = self.plus(out, self.scale(self.mul(ELF[a], ELG[b]),
                                                -sp.Rational(1, 2) * kab * (-1) ** (pb * (pa + pF))))
            if ELF[b] and ELG[a]:
                out = self.plus(out, self.scale(self.mul(ELF[b], ELG[a]),
                                                -sp.Rational(1, 2) * kab * (-1) ** (pa * pF)))
        return self.simplify(out)

    def is_zero_functional(self, P):
        """int P = 0 on compactly supported sections iff all Euler-Lagrange derivatives vanish (Remark c1-method)."""
        P = self.simplify(P)
        for a in range(len(self.fields)):
            if self.simplify(self.EL(P, a)):
                return False
        return True

    def size(self, P):
        return len(P)
