"""QA (PART C referee), third method for the scalar sector: special-edge composite chains in momentum space, expanded
in the external momenta (lambda) and in one injected momentum (eta), with exact tensor reduction (angular averages on
S^3) and exact radial integrals.  Written from the manuscript conventions (Ch1 vertex functions, Ch3 propagator
components, Lemma lem:cc-kernels kernel form) without engine code.  epsilon = 1 by scaling: the eps^0 coefficient of a
chain of superficial degree d is the lambda^d coefficient (Prop prop:cc-bound, Thm thm:cc-simplex).
Units: hbar/(16 pi^2) times the colour tensor C_w (as the engines).  One overall factor F0 = -1/2 (the hbar/2 of the
closure with the Koszul sign of the odd skeleton V1 = ([A,c],A*) with K_x on c and the special edge on A*) is common to
every chain used here; it is fixed once (see qa_mom_checks.py) and then all other chains are predictions.
Edges (eps = 1): K line X^{c,A_m}(p) = 2i p_m e^{-2p^2}/p^4 (c (x) A|_V, e^{ip(y-z)}), P^{A,A}(p) = 2 e^{-p^2}/p^2,
special edge s(p) = -(e^{-p^2}-1) (components of K_0 = -1 times e^{-p^2}-1).  V1 fibre +1; V2 = -(i/2)(p_P - q_out).a3."""
from fractions import Fraction as Fr
import math, itertools


class C:
    __slots__ = ('r', 'i')

    def __init__(self, r=0, i=0):
        self.r = Fr(r); self.i = Fr(i)

    def __add__(s, o):
        o = cc(o); return C(s.r + o.r, s.i + o.i)
    __radd__ = __add__

    def __sub__(s, o):
        o = cc(o); return C(s.r - o.r, s.i - o.i)

    def __rsub__(s, o):
        return cc(o) - s

    def __neg__(s):
        return C(-s.r, -s.i)

    def __mul__(s, o):
        o = cc(o); return C(s.r * o.r - s.i * o.i, s.r * o.i + s.i * o.r)
    __rmul__ = __mul__

    def __truediv__(s, o):
        o = cc(o); d = o.r * o.r + o.i * o.i
        return C((s.r * o.r + s.i * o.i) / d, (s.i * o.r - s.r * o.i) / d)

    def iszero(s):
        return s.r == 0 and s.i == 0

    def __eq__(s, o):
        o = cc(o); return s.r == o.r and s.i == o.i

    def __repr__(s):
        if s.i == 0:
            return str(s.r)
        return '(%s%+s*I)' % (s.r, s.i) if s.r != 0 else '%s*I' % s.i


def cc(x):
    return x if isinstance(x, C) else C(x, 0)


I = C(0, 1)
ONE = C(1, 0)

# ---------------------------------------------------------------- series
# term key: (c, nexp, lam, eta, mono)  value C ;  meaning  coef * lam^l eta^e * l^mono * (l^2)^(-nexp) * e^{-c l^2}
LMAX = [4]
EMAX = [1]


def S_mul(A, B):
    out = {}
    for (c1, n1, l1, e1, m1), v1 in A.items():
        for (c2, n2, l2, e2, m2), v2 in B.items():
            l = l1 + l2; e = e1 + e2
            if l > LMAX[0] or e > EMAX[0]:
                continue
            k = (c1 + c2, n1 + n2, l, e, (m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2], m1[3] + m2[3]))
            out[k] = out.get(k, C()) + v1 * v2
    return {k: v for k, v in out.items() if not v.iszero()}


def S_add(*As):
    out = {}
    for A in As:
        for k, v in A.items():
            out[k] = out.get(k, C()) + v
    return {k: v for k, v in out.items() if not v.iszero()}


def S_scale(A, s):
    s = cc(s)
    return {k: v * s for k, v in A.items()}


def S_const(x):
    return {(0, 0, 0, 0, (0, 0, 0, 0)): cc(x)}


Z4 = (0, 0, 0, 0)


def unit(mu):
    m = [0, 0, 0, 0]; m[mu] = 1
    return tuple(m)


class Mom:
    """affine momentum sl * l + lam * QL + eta * QE (QL, QE complex 4-vectors), sl in {+1,-1}"""

    def __init__(self, sl, QL, QE=None):
        self.sl = sl
        self.QL = [cc(x) for x in QL]
        self.QE = [cc(x) for x in (QE if QE is not None else [0, 0, 0, 0])]

    def __neg__(self):
        return Mom(-self.sl, [-x for x in self.QL], [-x for x in self.QE])

    def add(self, QL=None, QE=None):
        return Mom(self.sl, [a + cc(b) for a, b in zip(self.QL, QL or [0] * 4)], [a + cc(b) for a, b in zip(self.QE, QE or [0] * 4)])


def lin(v, p):
    """series of v.p (v complex 4-vector, p Mom)"""
    out = {}
    for mu in range(4):
        vm = cc(v[mu])
        if vm.iszero():
            continue
        if p.sl:
            k = (0, 0, 0, 0, unit(mu)); out[k] = out.get(k, C()) + vm * p.sl
        if not p.QL[mu].iszero():
            k = (0, 0, 1, 0, Z4); out[k] = out.get(k, C()) + vm * p.QL[mu]
        if not p.QE[mu].iszero():
            k = (0, 0, 0, 1, Z4); out[k] = out.get(k, C()) + vm * p.QE[mu]
    return {k: v for k, v in out.items() if not v.iszero()}


def comp(p, mu):
    e = [0, 0, 0, 0]; e[mu] = 1
    return lin(e, p)


def dot_lQ(p):
    """2 sl l.Q + Q^2 for (sl l + Q)^2 = l^2 + (2 sl l.Q + Q^2)"""
    out = {}
    for mu in range(4):
        if not p.QL[mu].iszero():
            k = (0, 0, 1, 0, unit(mu)); out[k] = out.get(k, C()) + p.QL[mu] * 2 * p.sl
        if not p.QE[mu].iszero():
            k = (0, 0, 0, 1, unit(mu)); out[k] = out.get(k, C()) + p.QE[mu] * 2 * p.sl
    QQ = {}
    for mu in range(4):
        a, b = p.QL[mu], p.QE[mu]
        QQ[(0, 0, 2, 0, Z4)] = QQ.get((0, 0, 2, 0, Z4), C()) + a * a
        QQ[(0, 0, 1, 1, Z4)] = QQ.get((0, 0, 1, 1, Z4), C()) + a * b * 2
        QQ[(0, 0, 0, 2, Z4)] = QQ.get((0, 0, 0, 2, Z4), C()) + b * b
    return S_add(out, {k: v for k, v in QQ.items() if not v.iszero()})


def S_pow(A, n):
    out = S_const(1)
    for _ in range(n):
        out = S_mul(out, A)
    return out


def expo(p, rate):
    """e^{-rate p^2} as series: e^{-rate l^2} sum_k (-rate u)^k / k!, u = 2 l.Q + Q^2"""
    u = S_scale(dot_lQ(p), -rate)
    out = {}
    term = S_const(1)
    kmax = 2 * (LMAX[0] + EMAX[0]) + 2
    for k in range(kmax + 1):
        out = S_add(out, S_scale(term, Fr(1, math.factorial(k))))
        term = S_mul(term, u)
        if not term:
            break
    return {(c + rate, n, l, e, m): v for (c, n, l, e, m), v in out.items()}


def invsq(p, power):
    """((sl l + Q)^2)^{-power} = (l^2)^{-power} sum_k binom(-power, k) u^k (l^2)^{-k}"""
    u = dot_lQ(p)
    out = {}
    term = S_const(1)
    kmax = 2 * (LMAX[0] + EMAX[0]) + 2
    for k in range(kmax + 1):
        b = Fr(1)
        for j in range(k):
            b *= Fr(-power - j, j + 1)
        out = S_add(out, {(c, n + power + k, l, e, m): v * b for (c, n, l, e, m), v in term.items()})
        term = S_mul(term, u)
        if not term:
            break
    return out


def special(p):
    """s(p) = 1 - e^{-p^2}"""
    return S_add(S_const(1), S_scale(expo(p, 1), -1))


def Kline(p, idx_vec):
    """2i (v.p) e^{-2p^2}/p^4, v the vector contracted with the A_m end"""
    return S_mul(S_scale(lin(idx_vec, p), 2 * I), S_mul(expo(p, 2), invsq(p, 2)))


def Kline_comp(p, m):
    e = [0] * 4; e[m] = 1
    return Kline(p, e)


def Pprop(p):
    return S_mul(S_const(2), S_mul(expo(p, 1), invsq(p, 1)))


# ---------------------------------------------------------------- integration
def dfact(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r


def ang(m):
    if any(x % 2 for x in m):
        return None
    tot = sum(m) // 2
    num = 1
    for x in m:
        num *= dfact(x - 1)
    return Fr(num, 2 ** tot * math.factorial(tot + 1)), tot


def integrate(A, lam, eta):
    """coefficient of lam^l eta^e of the loop integral, units 1/(16 pi^2): returns {'1': C, 'log2': C, 'log3': C, ...}"""
    groups = {}
    for (c, n, l, e, m), v in A.items():
        if l != lam or e != eta:
            continue
        a = ang(m)
        if a is None:
            continue
        w, tot = a
        j = tot - n
        groups.setdefault(j, {})
        groups[j][c] = groups[j].get(c, C()) + v * w
    res = {'1': C()}
    for j, wc in groups.items():
        wc = {c: v for c, v in wc.items() if not v.iszero()}
        if not wc:
            continue
        if j + 1 >= 0:
            for c, v in wc.items():
                assert c > 0, 'UV: zero rate'
                res['1'] = res['1'] + v * Fr(math.factorial(j + 1), c ** (j + 2))
        else:
            mm = -(j + 1) - 1
            for i in range(mm + 1):
                s = C()
                for c, v in wc.items():
                    s = s + v * (c ** i)
                assert s.iszero(), 'IR moment condition fails j=%d i=%d %s' % (j, i, s)
            for c, v in wc.items():
                if c == 1:
                    continue
                coef = v * Fr((-1) ** (mm + 1) * c ** mm, math.factorial(mm))
                for pr, ex in factor_small(c).items():
                    key = 'log%d' % pr
                    res[key] = res.get(key, C()) + coef * ex
    return {k: v for k, v in res.items() if not v.iszero()}


def factor_small(c):
    out = {}
    d = 2
    while c > 1:
        while c % d == 0:
            out[d] = out.get(d, 0) + 1; c //= d
        d += 1
    return out


def vadd(a, b, s=1):
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, C()) + v * s
    return {k: v for k, v in out.items() if not v.iszero()}


def vscale(a, s):
    return {k: v * cc(s) for k, v in a.items()}


F0 = Fr(-1, 2)          # overall factor (hbar/2 with the sign of the odd skeleton), fixed once
