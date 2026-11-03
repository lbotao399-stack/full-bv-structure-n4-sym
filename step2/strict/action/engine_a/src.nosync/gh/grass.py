"""Minimal super-commutative polynomial algebra (engine A).
Variables are integers; PAR[v] in {0,1} is the parity. Monomial key = sorted tuple of variables
(even variables may repeat, odd may not). Coefficients: any ring supporting + and * (sympy / Fraction / complex)."""
PAR = {}

def newvar(v, parity):
    PAR[v] = parity
    return v

def _mul_keys(k1, k2):
    seq = list(k1) + list(k2)
    odd = [x for x in seq if PAR[x]]
    if len(set(odd)) < len(odd):
        return None, 0
    # sign of sorting odd subsequence
    sign = 1
    o = odd[:]
    for i in range(len(o)):
        for j in range(len(o) - 1 - i):
            if o[j] > o[j + 1]:
                o[j], o[j + 1] = o[j + 1], o[j]
                sign = -sign
    return tuple(sorted(seq)), sign

class GA:
    __slots__ = ['t']
    def __init__(self, t=None):
        self.t = t if t is not None else {}
    @staticmethod
    def var(v, c=1):
        return GA({(v,): c})
    @staticmethod
    def const(c):
        return GA({(): c}) if c != 0 else GA()
    def __add__(self, o):
        if not isinstance(o, GA):
            o = GA.const(o)
        r = dict(self.t)
        for k, c in o.t.items():
            r[k] = r.get(k, 0) + c
            if r[k] == 0:
                del r[k]
        return GA(r)
    __radd__ = __add__
    def __neg__(self):
        return GA({k: -c for k, c in self.t.items()})
    def __sub__(self, o):
        return self + (-o)
    def __rsub__(self, o):
        return (-self) + o
    def __mul__(self, o):
        if not isinstance(o, GA):
            if o == 0:
                return GA()
            return GA({k: c * o for k, c in self.t.items()})
        r = {}
        for k1, c1 in self.t.items():
            for k2, c2 in o.t.items():
                k, s = _mul_keys(k1, k2)
                if k is None:
                    continue
                r[k] = r.get(k, 0) + s * c1 * c2
        return GA({k: c for k, c in r.items() if c != 0})
    def __rmul__(self, o):
        # scalar on the left (scalars are even)
        return self.__mul__(o) if not isinstance(o, GA) else o.__mul__(self)
    def dleft(self, v):
        """left derivative d/dv"""
        r = {}
        for k, c in self.t.items():
            if v not in k:
                continue
            idx = k.index(v)
            if PAR[v]:
                nodd = sum(PAR[x] for x in k[:idx])
                s = -1 if nodd % 2 else 1
                nk = k[:idx] + k[idx + 1:]
                r[nk] = r.get(nk, 0) + s * c
            else:
                m = k.count(v)
                nk = k[:idx] + k[idx + 1:]
                r[nk] = r.get(nk, 0) + m * c
        return GA({k: c for k, c in r.items() if c != 0})
    def degree_part(self, d):
        return GA({k: c for k, c in self.t.items() if len(k) == d})
    def coeff(self, key):
        return self.t.get(tuple(sorted(key)), 0)
    def map(self, f):
        return GA({k: f(c) for k, c in self.t.items() if f(c) != 0})

def dP(F, Pmat, varlist):
    """Costello partial_P = 1/2 sum_{a,b} P[a][b] d_b d_a  (left derivatives), so that dP(y_a y_b)=P[a][b]
    when P is graded symmetric."""
    r = GA()
    for a in varlist:
        Fa = F.dleft(a)
        if not Fa.t:
            continue
        for b in varlist:
            pab = Pmat.get((a, b), 0)
            if pab == 0:
                continue
            r = r + Fa.dleft(b) * (pab / 2 if not hasattr(pab, 'free_symbols') else pab / 2)
    return r
