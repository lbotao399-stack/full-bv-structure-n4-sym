"""
jets.py -- exact Grassmann-graded jet-polynomial engine (Gaussian-integer coefficients).

A jet variable is (field, comp, color, J) with J in N^4 (derivatives along R^4 only).
A polynomial is a dict {(ev, od): (re, im)}:
    ev = sorted tuple of even variable ids (with repetition),
    od = strictly increasing tuple of odd variable ids (Grassmann monomial in that order),
    coefficient re + i*im with Python ints (exact).
Monomial value = (product of evens) * od[0] od[1] ... od[-1].
"""
from bisect import bisect_left
from collections import defaultdict

ZERO4 = (0, 0, 0, 0)


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cneg(a):
    return (-a[0], -a[1])


class Vars:
    def __init__(self, parity):
        self.parity_of_field = parity
        self.id = {}
        self.key = []
        self.par = []

    def get(self, field, comp, color, J=ZERO4):
        k = (field, comp, color, J)
        i = self.id.get(k)
        if i is None:
            i = len(self.key)
            self.id[k] = i
            self.key.append(k)
            self.par.append(self.parity_of_field[field])
        return i

    def shift(self, v, mu):
        f, a, c, J = self.key[v]
        J2 = list(J)
        J2[mu] += 1
        return self.get(f, a, c, tuple(J2))


def _merge_odd(o1, o2):
    """product of Grassmann monomials o1*o2 -> (sign, sorted tuple) or (0, None)."""
    if not o1:
        return 1, o2
    if not o2:
        return 1, o1
    s1 = set(o1)
    for y in o2:
        if y in s1:
            return 0, None
    inv = 0
    for y in o2:
        inv += len(o1) - bisect_left(o1, y)
    merged = tuple(sorted(o1 + o2))
    return (-1 if inv % 2 else 1), merged


def padd_into(P, Q, scale=(1, 0)):
    for k, v in Q.items():
        w = cmul(v, scale) if scale != (1, 0) else v
        old = P.get(k)
        if old is None:
            P[k] = w
        else:
            n = cadd(old, w)
            if n == (0, 0):
                del P[k]
            else:
                P[k] = n
    return P


def padd(P, Q, scale=(1, 0)):
    R = dict(P)
    return padd_into(R, Q, scale)


def pscale(P, s):
    return {k: cmul(v, s) for k, v in P.items() if cmul(v, s) != (0, 0)}


def pmul(P, Q):
    R = {}
    for (e1, o1), c1 in P.items():
        for (e2, o2), c2 in Q.items():
            s, o = _merge_odd(o1, o2)
            if s == 0:
                continue
            e = tuple(sorted(e1 + e2))
            c = cmul(c1, c2)
            if s < 0:
                c = cneg(c)
            k = (e, o)
            old = R.get(k)
            if old is None:
                R[k] = c
            else:
                n = cadd(old, c)
                if n == (0, 0):
                    del R[k]
                else:
                    R[k] = n
    return R


def var(V, field, comp, color, J=ZERO4, coeff=(1, 0)):
    v = V.get(field, comp, color, J)
    if V.par[v]:
        return {((), (v,)): coeff}
    return {((v,), ()): coeff}


def const(c):
    return {((), ()): c}


def mono(V, coeff, vlist):
    """monomial coeff * v1 v2 ... vn (product in the given order)."""
    P = const(coeff)
    for v in vlist:
        if V.par[v]:
            P = pmul(P, {((), (v,)): (1, 0)})
        else:
            P = pmul(P, {((v,), ()): (1, 0)})
    return P


def dleft(V, P, v):
    R = {}
    if V.par[v]:
        for (e, o), c in P.items():
            if v in o:
                k = o.index(v)
                o2 = o[:k] + o[k + 1:]
                cc = c if k % 2 == 0 else cneg(c)
                key = (e, o2)
                padd_into(R, {key: cc})
    else:
        for (e, o), c in P.items():
            m = e.count(v)
            if m:
                k = e.index(v)
                e2 = e[:k] + e[k + 1:]
                padd_into(R, {(e2, o): (c[0] * m, c[1] * m)})
    return R


def dright(V, P, v):
    if not V.par[v]:
        return dleft(V, P, v)
    R = {}
    for (e, o), c in P.items():
        if v in o:
            k = o.index(v)
            n = len(o)
            o2 = o[:k] + o[k + 1:]
            cc = c if (n - 1 - k) % 2 == 0 else cneg(c)
            padd_into(R, {(e, o2): cc})
    return R


def Dtot(V, P, mu):
    """total derivative along x^mu (even derivation)."""
    R = {}
    for (e, o), c in P.items():
        # even factors
        done = set()
        for k, v in enumerate(e):
            if v in done:
                continue
            done.add(v)
            m = e.count(v)
            w = V.shift(v, mu)
            # remove one v, add w
            kk = e.index(v)
            e2 = list(e[:kk] + e[kk + 1:])
            e2.append(w)
            e2 = tuple(sorted(e2))
            padd_into(R, {(e2, o): (c[0] * m, c[1] * m)})
        for k, v in enumerate(o):
            w = V.shift(v, mu)
            rest = o[:k] + o[k + 1:]
            if w in rest:
                continue
            pos = bisect_left(rest, w)
            sgn = -1 if abs(pos - k) % 2 else 1
            o2 = rest[:pos] + (w,) + rest[pos:]
            padd_into(R, {(e, o2): c if sgn > 0 else cneg(c)})
    return R


def vars_of(P):
    s = set()
    for (e, o) in P:
        s.update(e)
        s.update(o)
    return s


def euler(V, P, field, comp, color, side='l'):
    """Euler-Lagrange operator  sum_J (-D)^J  d_{side}/d phi_J  applied to density P."""
    out = {}
    present = [v for v in vars_of(P) if V.key[v][0] == field and V.key[v][1] == comp and V.key[v][2] == color]
    for v in present:
        J = V.key[v][3]
        T = dleft(V, P, v) if side == 'l' else dright(V, P, v)
        sign = -1 if sum(J) % 2 else 1
        for mu in range(4):
            for _ in range(J[mu]):
                T = Dtot(V, T, mu)
        padd_into(out, T, (sign, 0))
    return out


def field_families(V, P):
    s = set()
    for v in vars_of(P):
        f, a, c, J = V.key[v]
        s.add((f, a, c))
    return sorted(s)


def is_total_derivative(V, P, fields=None):
    """P is a total divergence iff all Euler operators vanish (polynomial densities, no constant)."""
    bad = []
    for fam in field_families(V, P):
        if fields is not None and fam[0] not in fields:
            continue
        E = euler(V, P, *fam, side='l')
        if E:
            bad.append((fam, len(E)))
    return (len(bad) == 0), bad


def antibracket_density(V, F, G, pairs, families):
    """density of (F,G)_phys = sum_phi [ dr F/dphi * dl G/dphi^* - dr F/dphi^* * dl G/dphi ]
    (Euler derivatives).  families: list of (field, comp, color) for fields phi; antifield from pairs."""
    out = {}
    for (f, a, c) in families:
        fs = pairs[f]
        t1 = pmul(euler(V, F, f, a, c, 'r'), euler(V, G, fs, a, c, 'l'))
        t2 = pmul(euler(V, F, fs, a, c, 'r'), euler(V, G, f, a, c, 'l'))
        padd_into(out, t1)
        padd_into(out, t2, (-1, 0))
    return out


def nterms(P):
    return len(P)


def _index(P):
    idx = defaultdict(dict)
    for key, c in P.items():
        e, o = key
        for v in set(e) | set(o):
            idx[v][key] = c
    return idx


def euler_all(V, P, side='l', only=None):
    """all Euler-Lagrange derivatives of density P: dict (field,comp,color) -> poly.
    only: optional set of field names to restrict to."""
    idx = _index(P)
    res = defaultdict(dict)
    for v, sub in idx.items():
        f, a, c, J = V.key[v]
        if only is not None and f not in only:
            continue
        T = dleft(V, sub, v) if side == 'l' else dright(V, sub, v)
        for mu in range(4):
            for _ in range(J[mu]):
                T = Dtot(V, T, mu)
        padd_into(res[(f, a, c)], T, (-1, 0) if sum(J) % 2 else (1, 0))
    return {k: v for k, v in res.items() if v}


def antibracket_all(V, F, G, pairs, families, eps=None):
    """density of (F,G)_eps = sum_phi eps_phi [ dr F/dphi dl G/dphi^* - dr F/dphi^* dl G/dphi ]
    (Euler derivatives; kappa = 1).  eps: dict field -> +1/-1 (default all +1)."""
    names = set(f for f, _, _ in families) | set(pairs[f] for f, _, _ in families)
    Fr = euler_all(V, F, 'r', names)
    Gl = euler_all(V, G, 'l', names)
    out = {}
    for (f, a, c) in families:
        fs = pairs[f]
        e = 1 if eps is None else eps[f]
        x1 = Fr.get((f, a, c)); y1 = Gl.get((fs, a, c))
        if x1 and y1:
            padd_into(out, pmul(x1, y1), (e, 0))
        x2 = Fr.get((fs, a, c)); y2 = Gl.get((f, a, c))
        if x2 and y2:
            padd_into(out, pmul(x2, y2), (-e, 0))
    return out


def total_derivative_test(V, P):
    E = euler_all(V, P, 'l')
    bad = {k: len(v) for k, v in E.items() if v}
    return len(bad) == 0, bad
