"""
engine.py -- exact Grassmann-graded jet-polynomial engine with explicit coordinates x^mu.

Extension of classical/jets.py (same polynomial layout) by
  * Gaussian-rational coefficients (pairs of fractions.Fraction),
  * coordinate variables ('x', mu, 0, 0000) with total derivative D_nu x^mu = delta,
  * passive variables (ghost coordinates on g[1], kind 'gh') with D_nu = 0,
  * the [C07 Lemma 4.3.2 / Sec 10.1] bracket {F,G} for local functionals.

A polynomial is a dict {(ev, od): (re, im)}:
  ev = sorted tuple of even variable ids (with repetition),
  od = strictly increasing tuple of odd variable ids (Grassmann monomial in that order).
Monomial value = (product of evens) * od[0] od[1] ... od[-1].
"""
from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction as Fr

ZERO4 = (0, 0, 0, 0)
FIELDS = ('c', 'A', 'lam', 'As', 'lams', 'cs')
# total parity [classical_spec Sec 2]
FIELD_PARITY = {'c': 1, 'A': 0, 'lam': 1, 'As': 1, 'lams': 0, 'cs': 0, 'x': 0}
ANTI = {'c': 'cs', 'A': 'As', 'lam': 'lams'}
# C07 even-first canonical pairs (x_i, xi_i) with <x_i, xi_i> = +1 [classical_spec Sec 3]
#   (A, A^*), (lam^*, lam), (c^*, c)
EPS_BRACKET = {'A': 1, 'lam': -1, 'c': -1}

ONE = (Fr(1), Fr(0))
ZEROC = (Fr(0), Fr(0))


def C(re, im=0):
    return (Fr(re), Fr(im))


def cadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def cmul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def cneg(a):
    return (-a[0], -a[1])


def ciszero(a):
    return a[0] == 0 and a[1] == 0


class Vars:
    def __init__(self):
        self.id = {}
        self.key = []
        self.par = []
        self.kind = []          # 'f' field, 'x' coordinate, 'g' ghost
        self.gpar = {}

    def get(self, field, comp, color=0, J=ZERO4):
        k = (field, comp, color, J)
        i = self.id.get(k)
        if i is None:
            i = len(self.key)
            self.id[k] = i
            self.key.append(k)
            if field == 'x':
                self.par.append(0)
                self.kind.append('x')
            elif field == 'gh':
                self.par.append(self.gpar[comp])
                self.kind.append('g')
            else:
                self.par.append(FIELD_PARITY[field])
                self.kind.append('f')
        return i

    def ghost(self, label, parity):
        self.gpar[label] = parity
        return self.get('gh', label)

    def shift(self, v, mu):
        f, a, c, J = self.key[v]
        J2 = list(J)
        J2[mu] += 1
        return self.get(f, a, c, tuple(J2))


def _merge_odd(o1, o2):
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
    return (-1 if inv % 2 else 1), tuple(sorted(o1 + o2))


def padd_into(P, Q, scale=ONE):
    one = (scale == ONE)
    for k, v in Q.items():
        w = v if one else cmul(v, scale)
        old = P.get(k)
        if old is None:
            if not ciszero(w):
                P[k] = w
        else:
            n = (old[0] + w[0], old[1] + w[1])
            if n[0] == 0 and n[1] == 0:
                del P[k]
            else:
                P[k] = n
    return P


def padd(*Ps):
    R = {}
    for P in Ps:
        padd_into(R, P)
    return R


def psub(P, Q):
    R = dict(P)
    return padd_into(R, Q, (Fr(-1), Fr(0)))


def pscale(P, s):
    if not isinstance(s, tuple):
        s = (Fr(s), Fr(0))
    R = {}
    for k, v in P.items():
        w = cmul(v, s)
        if not ciszero(w):
            R[k] = w
    return R


def pmul(P, Q):
    R = {}
    for (e1, o1), c1 in P.items():
        for (e2, o2), c2 in Q.items():
            s, o = _merge_odd(o1, o2)
            if s == 0:
                continue
            e = tuple(sorted(e1 + e2))
            c = (c1[0] * c2[0] - c1[1] * c2[1], c1[0] * c2[1] + c1[1] * c2[0])
            if s < 0:
                c = (-c[0], -c[1])
            k = (e, o)
            old = R.get(k)
            if old is None:
                R[k] = c
            else:
                n = (old[0] + c[0], old[1] + c[1])
                if n[0] == 0 and n[1] == 0:
                    del R[k]
                else:
                    R[k] = n
    return R


def const(c):
    if not isinstance(c, tuple):
        c = (Fr(c), Fr(0))
    if ciszero(c):
        return {}
    return {((), ()): c}


def var(V, field, comp, color=0, J=ZERO4, coeff=ONE):
    v = V.get(field, comp, color, J)
    if V.par[v]:
        return {((), (v,)): coeff}
    return {((v,), ()): coeff}


def vpoly(V, v, coeff=ONE):
    if V.par[v]:
        return {((), (v,)): coeff}
    return {((v,), ()): coeff}


def mono(V, coeff, vlist):
    P = const(coeff)
    for v in vlist:
        P = pmul(P, vpoly(V, v))
    return P


def parity(V, P):
    """total parity of a homogeneous polynomial (asserts homogeneity)."""
    ps = set()
    for (e, o) in P:
        ps.add(len(o) % 2)
    if not ps:
        return 0
    assert len(ps) == 1, 'inhomogeneous parity'
    return ps.pop()


def dleft(V, P, v):
    R = {}
    if V.par[v]:
        for (e, o), c in P.items():
            if v in o:
                k = o.index(v)
                o2 = o[:k] + o[k + 1:]
                padd_into(R, {(e, o2): c if k % 2 == 0 else cneg(c)})
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
            padd_into(R, {(e, o2): c if (n - 1 - k) % 2 == 0 else cneg(c)})
    return R


def Dtot(V, P, mu):
    """total derivative along x^mu: fields shift jets, x^nu -> delta, ghosts constant."""
    R = {}
    for (e, o), c in P.items():
        done = set()
        for v in e:
            if v in done:
                continue
            done.add(v)
            m = e.count(v)
            kd = V.kind[v]
            if kd == 'g':
                continue
            kk = e.index(v)
            if kd == 'x':
                if V.key[v][1] != mu:
                    continue
                e2 = e[:kk] + e[kk + 1:]
                padd_into(R, {(e2, o): (c[0] * m, c[1] * m)})
                continue
            w = V.shift(v, mu)
            e2 = list(e[:kk] + e[kk + 1:])
            e2.append(w)
            padd_into(R, {(tuple(sorted(e2)), o): (c[0] * m, c[1] * m)})
        for k, v in enumerate(o):
            if V.kind[v] != 'f':
                continue
            w = V.shift(v, mu)
            rest = o[:k] + o[k + 1:]
            if w in rest:
                continue
            pos = bisect_left(rest, w)
            sgn = -1 if abs(pos - k) % 2 else 1
            o2 = rest[:pos] + (w,) + rest[pos:]
            padd_into(R, {(e, o2): c if sgn > 0 else cneg(c)})
    return R


def _index(P, V):
    idx = defaultdict(dict)
    for key, c in P.items():
        e, o = key
        for v in set(e) | set(o):
            if V.kind[v] == 'f':
                idx[v][key] = c
    return idx


def euler_all(V, P, side='l', only=None):
    """Euler-Lagrange derivatives sum_J (-D)^J d/dphi_J, over field variables only."""
    idx = _index(P, V)
    res = defaultdict(dict)
    for v, sub in idx.items():
        f, a, col, J = V.key[v]
        if only is not None and f not in only:
            continue
        T = dleft(V, sub, v) if side == 'l' else dright(V, sub, v)
        for mu in range(4):
            for _ in range(J[mu]):
                T = Dtot(V, T, mu)
        padd_into(res[(f, a, col)], T, (Fr(-1), Fr(0)) if sum(J) % 2 else ONE)
    return {k: v for k, v in res.items() if v}


def antibracket(V, F, G, dimg, eps=EPS_BRACKET):
    """sum_phi eps_phi [ dr F/dphi dl G/dphi^* - dr F/dphi^* dl G/dphi ] (Euler derivatives, (-,-) = delta)."""
    Fr_ = euler_all(V, F, 'r')
    Gl_ = euler_all(V, G, 'l')
    out = {}
    for f, fs in ANTI.items():
        e = eps[f]
        se = (Fr(e), Fr(0))
        for (name, a, col) in list(Fr_.keys()):
            pass
        keysF = set(k for k in Fr_ if k[0] in (f, fs))
        keysG = set(k for k in Gl_ if k[0] in (f, fs))
        comps = set((k[1], k[2]) for k in keysF | keysG)
        for (a, col) in comps:
            x1 = Fr_.get((f, a, col)); y1 = Gl_.get((fs, a, col))
            if x1 and y1:
                padd_into(out, pmul(x1, y1), se)
            x2 = Fr_.get((fs, a, col)); y2 = Gl_.get((f, a, col))
            if x2 and y2:
                padd_into(out, pmul(x2, y2), (Fr(-e), Fr(0)))
    return out


def bracket(V, F, G, dimg=None):
    """[C07 Lemma 4.3.2, Sec 10.1] bracket {F,G} = Delta_0(FG) - Delta_0(F)G - (-1)^|F| F Delta_0(G),
    Delta_0 = sum_i d_{x_i} d_{xi_i} with even-first pairs (A,A^*), (lam^*,lam), (c^*,c), <x_i,xi_i> = +1.
    For one canonical pair (x even, xi odd):  {f,g} = (-1)^|f| [ dr f/dx dl g/dxi - dr f/dxi dl g/dx ]."""
    s = parity(V, F)
    R = antibracket(V, F, G, dimg)
    return pscale(R, -1) if s else R


def total_derivative_test(V, P):
    """P (density) is a total divergence iff all field Euler derivatives vanish."""
    E = euler_all(V, P, 'l')
    bad = {k: len(v) for k, v in E.items() if v}
    return len(bad) == 0, bad


def is_zero_functional(V, P):
    return total_derivative_test(V, P)[0]


def nterms(P):
    return len(P)


def evaluate(V, P, values):
    """numeric evaluation of an even-only polynomial: values: var id -> (re,im) Fractions."""
    tot = ZEROC
    for (e, o), c in P.items():
        assert not o
        t = c
        for v in e:
            t = cmul(t, values[v])
        tot = cadd(tot, t)
    return tot
