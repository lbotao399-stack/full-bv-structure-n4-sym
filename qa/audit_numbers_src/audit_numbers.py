#!/usr/bin/env python3
"""
Adversarial audit (lens: numbers).  Fresh, self-contained, exact (sympy) recomputation of
  (1) log eps coefficient of the A|_{W^perp} two-point function (and its eps^-1 coefficient),
  (2) log eps coefficient of partial_delta I^CT[eps] for O = (iota_w1 A, iota_w2 A)
      (traceless and trace combinations),
  (3) eps^-1 coefficient of the A|_V two-point function (and its log eps coefficients).
No code from the project is imported or copied.  Explicit bases (su(2) basis, 10d gamma
matrices, Lorentz components) live only here.

Conventions (DESIGN.md / classical_spec.md, stated there, re-derived here):
  weight e^{+S/hbar};  S = int dvol( -1/4 (F_A,F_A) + 1/2 (lambda, dslash_A lambda)
                                     + (d_A c, A^*) + ... )
  (F_A,F_A) = h * sum_{M,N} F_MN^a F_MN^b (X_a,X_b)/kappa ... written below with
  (X,Y) = kappa X^a Y^a, (F,F) = h kappa sum_{M,N} F^a_MN F^a_MN ; Hodge pairing <=> h = 1/2.
  P(eps,infty) on every line = (numerator) * exp(-eps k^2)/k^2  (Feynman gauge, H = k^2 on
  every component).  A-A covariance = delta_MN delta^ab exp(-eps k^2)/(h kappa k^2) > 0.
One-loop background part of W (Gaussian / Pfaffian / determinant expansion):
  bosons  : -(hbar/2) Tr log(1 - P I'')        -> (hbar/2)Tr(P V4) + (hbar/4)Tr(P V3 P V3)
  lambda  : +(hbar/2) Tr log(1 + D^-1 V)       -> -(hbar/4) Tr(D^-1 V D^-1 V)
  ghost   : + hbar    Tr log(1 + G X)          -> -(hbar/2) Tr(G X G X)
  insertion (family I + delta O, O quadratic):
     partial_delta W_1 = (hbar/2) Tr[ P O'' (1 - P I'')^{-1} ]
Singular parts (t_i in (eps, infty), T = sum t_i):
  bubble  : int dt1 dt2 T^{-2}[N0 + N1/(2T)] e^{-x(1-x)T p^2}
            eps^-1 coefficient  int_0^1 (N1/2) min(x,1-x) dx
            log eps coefficient -int_0^1 N0 dx + int_0^1 (N1/2) x(1-x) p^2 dx
  tadpole : int_k e^{-eps k^2}/k^2 = 1/(16 pi^2 eps)
  O-bubble (lines k,k) : log eps coefficient -N0
  triangle (p -> 0)    : log eps coefficient -N1/4
All results below are in units hbar/(16 pi^2).
"""
import json
import sys
from itertools import permutations

import numpy as np
import sympy as sp

R = sp.Rational
h = sp.Symbol('h', positive=True)
kap = sp.Integer(1)
NA = h * kap

ks = sp.symbols('k0:4')
ps = sp.symbols('p0:4')
ls = sp.symbols('l0:4')
x = sp.Symbol('x')
P2 = sum(q**2 for q in ps)

V = list(range(4))
W = list(range(4, 10))
import os
ALG = os.environ.get('ALG', 'su2')


def v10(v4):
    return list(v4) + [0] * 6


def neg(v):
    return [-a for a in v]


def add(u, v):
    return [a + b for a, b in zip(u, v)]


ZERO10 = [0] * 10
K10 = v10(ks)
P10 = v10(ps)

# ---------------------------------------------------------------- Lie algebra: [T^a,T^b] = f^{abc} T^c (real basis)
if ALG == 'su2':
    DIMG = 3
    FNZ = {(0, 1, 2): sp.Integer(1)}
    CAS = -2   # tr(ad ad) = CAS * delta in this basis
elif ALG == 'su3':
    DIMG = 8
    r3 = sp.sqrt(3) / 2
    FNZ = {(0, 1, 2): sp.Integer(1), (0, 3, 6): R(1, 2), (0, 4, 5): -R(1, 2), (1, 3, 5): R(1, 2),
           (1, 4, 6): R(1, 2), (2, 3, 4): R(1, 2), (2, 5, 6): -R(1, 2), (3, 4, 7): r3, (5, 6, 7): r3}
    CAS = -3
else:
    raise SystemExit('ALG must be su2 or su3')
EPS = {(a, b, c): sp.Integer(0) for a in range(DIMG) for b in range(DIMG) for c in range(DIMG)}
for (a, b, c), v in FNZ.items():
    for (i, j, k_), sg in (((a, b, c), 1), ((b, c, a), 1), ((c, a, b), 1), ((b, a, c), -1), ((a, c, b), -1), ((c, b, a), -1)):
        EPS[(i, j, k_)] = sg * v


def trad(a, b):
    # tr_g(ad_{T^a} ad_{T^b}); (ad_{T^a})_{c d} = f^{a d c}
    return sp.nsimplify(sp.expand(sum(EPS[(a, d, c)] * EPS[(b, c, d)] for c in range(DIMG) for d in range(DIMG))))


# Jacobi identity and Killing form
for a in range(DIMG):
    for b in range(DIMG):
        assert trad(a, b) == (CAS if a == b else 0)
for a in range(DIMG):
    for b in range(DIMG):
        for c in range(DIMG):
            for d in range(DIMG):
                jac = sum(EPS[(a, b, e)] * EPS[(e, c, d)] + EPS[(b, c, e)] * EPS[(e, a, d)] + EPS[(c, a, e)] * EPS[(e, b, d)]
                          for e in range(DIMG))
                assert sp.expand(jac) == 0
TRAD00 = trad(0, 0)  # (X,Y) = kappa delta  => unit  hbar tr(ad ad)/((X,Y)16pi^2) = TRAD00 hbar/(16pi^2 kappa)

LABELS = [(a, M) for a in range(DIMG) for M in range(10)]


# ---------------------------------------------------------------- vertex functions from S
# I_3 = -NA int f^{abc} (d_M A^a_N) A^b_M A^c_N ;  I_4 = -(NA/4) int f^{xbc} f^{xde} A^b_M A^c_N A^d_M A^e_N
# multilinear coefficient on plane waves e_{(a,M)} e^{i k.x}, all momenta incoming, (2pi)^4 delta stripped.
def gamma3(legs):
    tot = 0
    for s in permutations(range(3)):
        (a1, M1, k1), (a2, M2, _), (a3, M3, _) = legs[s[0]], legs[s[1]], legs[s[2]]
        if M1 != M3 or M2 >= 4:
            continue
        fv = EPS[(a1, a2, a3)]
        if fv == 0:
            continue
        tot += -NA * fv * sp.I * k1[M2]
    return tot


def gamma4(legs):
    tot = 0
    for s in permutations(range(4)):
        (b, M1, _), (c, N1, _), (d, M2, _), (e, N2, _) = [legs[i] for i in s]
        if M1 != M2 or N1 != N2:
            continue
        ff = sp.expand(sum(EPS[(xx, b, c)] * EPS[(xx, d, e)] for xx in range(DIMG)))
        if ff:
            tot += -NA / 4 * ff
    return tot


# brute-force check of gamma3/gamma4 against direct expansion of -(NA/4) F.F on plane waves
def brute_vertex_check():
    """Expand -(NA/4) sum_{M,N,a} F^a_MN F^a_MN directly on A = sum_j t_j e_{choice_j} e^{i k_j x}
    (d_M -> i k_j[M] on the j-th term) and compare the t_0 t_1 ... multilinear coefficient with
    gamma3 / gamma4, for random and for all-pairs-of-interest index choices."""
    import random
    random.seed(1)
    kk = [[sp.Symbol(f'q{j}_{mu}') for mu in range(4)] + [0] * 6 for j in range(4)]
    t = sp.symbols('t0:4')
    results = []
    for nlegs in (3, 4):
        js = list(range(nlegs))
        choices = [[(random.randrange(DIMG), random.randrange(10)) for _ in js] for _ in range(12)]
        if nlegs == 3:
            choices += [[(0, 0), (1, 1), (2, 0)], [(0, 4), (1, 0), (2, 4)], [(0, 0), (1, 4), (2, 4)]]
        else:
            choices += [[(0, 0), (1, 1), (0, 0), (1, 1)], [(0, 4), (1, 0), (0, 4), (1, 0)],
                        [(0, 4), (1, 5), (1, 4), (0, 5)], [(0, 4), (0, 4), (1, 4), (1, 4)]]
        for choice in choices:
            def A(a, M):
                return sum(t[j] for j in js if choice[j] == (a, M))

            def dA(Md, a, N):
                return sum(t[j] * sp.I * kk[j][Md] for j in js if choice[j] == (a, N))
            Ssum = 0
            for M in range(10):
                for N in range(10):
                    for a in range(DIMG):
                        F = dA(M, a, N) - dA(N, a, M) + sum(EPS[(a, b, c)] * A(b, M) * A(c, N)
                                                            for b in range(DIMG) for c in range(DIMG))
                        if F != 0:
                            Ssum += F * F
            Sexpr = sp.expand(-NA / 4 * Ssum)
            mono = sp.Mul(*[t[j] for j in js])
            coeff = sp.Poly(Sexpr, *t[:nlegs]).coeff_monomial(mono) if Sexpr != 0 else 0
            legs = [(choice[j][0], choice[j][1], kk[j]) for j in js]
            mine = gamma3(legs) if nlegs == 3 else gamma4(legs)
            results.append(sp.simplify(coeff - mine) == 0)
    return all(results), len(results)


# ---------------------------------------------------------------- 10d gamma matrices (exact Gaussian integers)
class GI:
    """exact Gaussian-integer matrix: re + i*im, int64 entries (no BLAS)."""
    def __init__(self, re, im):
        self.re = np.asarray(re, dtype=np.int64)
        self.im = np.asarray(im, dtype=np.int64)

    def __matmul__(self, o):
        return GI(self.re @ o.re - self.im @ o.im, self.re @ o.im + self.im @ o.re)

    def __add__(self, o):
        return GI(self.re + o.re, self.im + o.im)

    def __neg__(self):
        return GI(-self.re, -self.im)

    def eq(self, o):
        return np.array_equal(self.re, o.re) and np.array_equal(self.im, o.im)

    def trace(self):
        return int(np.trace(self.re)), int(np.trace(self.im))


def kronGI(a, b):
    return GI(np.kron(a.re, b.re) - np.kron(a.im, b.im), np.kron(a.re, b.im) + np.kron(a.im, b.re))


def gamma10():
    s1 = GI([[0, 1], [1, 0]], [[0, 0], [0, 0]])
    s2 = GI([[0, 0], [0, 0]], [[0, -1], [1, 0]])
    s3 = GI([[1, 0], [0, -1]], [[0, 0], [0, 0]])
    I2 = GI([[1, 0], [0, 1]], [[0, 0], [0, 0]])

    def kron(lst):
        out = GI([[1]], [[0]])
        for m in lst:
            out = kronGI(out, m)
        return out
    G = []
    for j in range(5):
        G.append(kron([s3] * j + [s1] + [I2] * (4 - j)))
        G.append(kron([s3] * j + [s2] + [I2] * (4 - j)))
    chi = kron([s3] * 5)
    return G, chi


def check_and_trace_identity():
    G, chi = gamma10()
    Id = GI(np.eye(32, dtype=np.int64), np.zeros((32, 32), dtype=np.int64))
    Z = GI(np.zeros((32, 32), dtype=np.int64), np.zeros((32, 32), dtype=np.int64))
    for M in range(10):
        for N in range(10):
            ac = G[M] @ G[N] + G[N] @ G[M]
            want = GI(2 * Id.re, Id.im) if M == N else Z
            assert ac.eq(want)
        assert (chi @ G[M]).eq(-(G[M] @ chi))
    assert (chi @ chi).eq(Id)
    # tr over S_+ (chi=+1, 16 dim) of sigmabar^A sigma^M sigmabar^B sigma^N
    #   = (1/2) tr32((1+chi) G^A G^M G^B G^N)
    bad = 0
    for A_ in V:
        for B_ in V:
            for M in range(10):
                for N in range(10):
                    prod = G[A_] @ G[M] @ G[B_] @ G[N]
                    t1 = prod.trace()
                    t2 = (chi @ prod).trace()
                    tr2 = (t1[0] + t2[0], t1[1] + t2[1])   # = 2 * tr16
                    want = 16 * ((A_ == M) * (B_ == N) - (A_ == B_) * (M == N) + (A_ == N) * (M == B_))
                    if tr2 != (2 * want, 0):
                        bad += 1
    assert bad == 0
    assert (Id.trace()[0] + chi.trace()[0]) == 32  # dim S_+ = 16
    return True


def tr16(A_, M, B_, N):
    return 16 * (int(A_ == M) * int(B_ == N) - int(A_ == B_) * int(M == N) + int(A_ == N) * int(M == B_))


# ---------------------------------------------------------------- Schwinger singular parts
def mint(f):
    return sp.integrate(f * x, (x, 0, R(1, 2))) + sp.integrate(f * (1 - x), (x, R(1, 2), 1))


def bubble_sing(num, sigma):
    """num(k,p) with propagators G(k) G(k + sigma p).  Returns (E, L): singular part
    E/eps + L log eps in units 1/(16 pi^2)."""
    sub = {ks[i]: ls[i] - x * sigma * ps[i] for i in range(4)}
    poly = sp.Poly(sp.expand(sp.sympify(num).subs(sub, simultaneous=True)), *ls)
    N0 = 0
    N1 = 0
    for mon, c in poly.terms():
        d = sum(mon)
        if d == 0:
            N0 += c
        elif d == 2:
            if max(mon) == 2:
                N1 += c
        elif d % 2 == 1:
            continue
        else:
            raise ValueError('unexpected l-degree')
    N1 = sp.expand(N1)
    assert all(sp.diff(N1, q) == 0 for q in ps), 'N1 depends on p'
    E = sp.simplify(mint(N1 / 2))
    L = sp.expand(-sp.integrate(sp.expand(N0), (x, 0, 1)) + sp.integrate(N1 / 2 * x * (1 - x), (x, 0, 1)) * P2)
    return E, L


# ---------------------------------------------------------------- two-point graphs
def isV(r):
    return r[1] < 4


def isW(r):
    return r[1] >= 4


def boson_bubble(al, be, filt):
    num = 0
    for r in LABELS:
        for s in LABELS:
            if not filt(r, s):
                continue
            g1 = gamma3([(al[0], al[1], P10), (r[0], r[1], K10), (s[0], s[1], neg(add(P10, K10)))])
            if g1 == 0:
                continue
            g2 = gamma3([(be[0], be[1], neg(P10)), (s[0], s[1], add(P10, K10)), (r[0], r[1], neg(K10))])
            num += g1 * g2
    num = sp.expand(num / 2 / NA**2)  # (hbar/2) Tr(P V3(al) P V3(be)), P = G/NA
    return bubble_sing(num, +1)


def boson_tadpole(al, be, filt):
    s = 0
    for r in LABELS:
        if not filt(r):
            continue
        s += gamma4([(al[0], al[1], P10), (be[0], be[1], neg(P10)), (r[0], r[1], K10), (r[0], r[1], neg(K10))])
    return sp.simplify(s / 2 / NA), sp.Integer(0)  # (hbar/2) Tr(P V4), int G = 1/(16pi^2 eps)


def fermion_bubble(al, be):
    a, M = al
    b, N = be
    km = [ks[i] - ps[i] for i in range(4)]
    tr = 0
    for A_ in V:
        for B_ in V:
            t = tr16(A_, M, B_, N)
            if t:
                tr += t * ks[A_] * km[B_]
    num = sp.expand(R(1, 2) * tr * trad(a, b))  # -(hbar/2) (1/i)^2 tr16[...] tr_g
    return bubble_sing(num, -1)


def ghost_bubble(al, be):
    a, M = al
    b, N = be
    if M >= 4 or N >= 4:
        return sp.Integer(0), sp.Integer(0)
    num = sp.expand(ks[M] * (ks[N] - ps[N]) * trad(a, b))  # -hbar (i k_M)(i (k-p)_N) tr_g
    return bubble_sing(num, -1)


def two_point(al, be):
    out = {}
    out['VVV bubble (A|_V lines)'] = boson_bubble(al, be, lambda r, s: isV(r) and isV(s))
    out['bubble (A|_W lines)'] = boson_bubble(al, be, lambda r, s: isW(r) and isW(s))
    out['bubble (mixed A|_V,A|_W lines)'] = boson_bubble(al, be, lambda r, s: isV(r) != isV(s))
    out['tadpole (A|_V line)'] = boson_tadpole(al, be, isV)
    out['tadpole (A|_W line)'] = boson_tadpole(al, be, isW)
    out['ghost loop'] = ghost_bubble(al, be)
    out['lambda loop'] = fermion_bubble(al, be)
    E = sp.simplify(sum(v[0] for v in out.values()))
    L = sp.expand(sum(v[1] for v in out.values()))
    return out, (E, L)


# ---------------------------------------------------------------- insertion O = (iota_w1 A, iota_w2 A)
def Opp(r, s, w1, w2):
    if r[0] != s[0]:
        return 0
    return kap * (w1[r[1]] * w2[s[1]] + w2[r[1]] * w1[s[1]])


def insertion(al, be, w1, w2):
    res = {}
    # vacuum: (hbar/2) Tr(P O'')  -> eps^-1
    vac = sum(Opp(r, r, w1, w2) for r in LABELS) / 2 / NA
    res['vacuum eps^-1'] = sp.simplify(vac)
    # O-bubble with quartic vertex: (hbar/2) Tr(P O'' P V4) ; lines k,k ; log eps coeff = -N0
    s4 = 0
    for r in LABELS:
        for s in LABELS:
            o = Opp(r, s, w1, w2)
            if o == 0:
                continue
            s4 += o * gamma4([(s[0], s[1], ZERO10), (r[0], r[1], ZERO10), (al[0], al[1], ZERO10), (be[0], be[1], ZERO10)])
    res['O + quartic: log'] = sp.simplify(-(s4 / 2 / NA**2))
    # triangle: (hbar/2)[Tr(P O'' P V3(al) P V3(be)) + (al<->be)], p -> 0 ; log eps coeff = -N1/4
    num = 0
    for (A1, A2, sg) in ((al, be, 1), (be, al, -1)):
        pp = [sg * q for q in P10]
        pp = [0] * 10  # p -> 0 for the dimension-0 coefficient
        for r1 in LABELS:
            for r3 in LABELS:
                o = Opp(r3, r1, w1, w2)
                if o == 0:
                    continue
                for r2 in LABELS:
                    g1 = gamma3([(A1[0], A1[1], pp), (r1[0], r1[1], K10), (r2[0], r2[1], neg(add(K10, pp)))])
                    if g1 == 0:
                        continue
                    g2 = gamma3([(A2[0], A2[1], neg(pp)), (r2[0], r2[1], add(K10, pp)), (r3[0], r3[1], neg(K10))])
                    num += o * g1 * g2
    num = sp.expand(num / 2 / NA**3)
    poly = sp.Poly(num, *ks)
    N1 = 0
    for mon, c in poly.terms():
        if sum(mon) == 2 and max(mon) == 2:
            N1 += c
        elif sum(mon) == 2:
            pass
        elif sum(mon) != 0 or c != 0:
            raise ValueError('unexpected triangle numerator')
    res['triangle: log'] = sp.simplify(-N1 / 4)
    res['total log'] = sp.simplify(res['O + quartic: log'] + res['triangle: log'])
    return res


# ---------------------------------------------------------------- numerical check of the singular-part rules
def check_extraction():
    import mpmath as mp
    mp.mp.dps = 30
    # F(eps) = int_{eps}^inf int_{eps}^inf dt1 dt2 T^{-2} [N0 + N1/(2T)] e^{-x(1-x) T p2}
    # test N0 = x^2 p2 (log type), N1 = 3 (power type); p2 = 1
    p2 = mp.mpf(1)

    def F(eps):
        def inner(xx):
            m = min(xx, 1 - xx)
            c = eps / m
            a = xx * (1 - xx) * p2
            # int_c^inf T^{-1} e^{-aT} dT = E1(a c) ; int_c^inf T^{-2} e^{-aT} dT = a*Gamma(-1, a c)
            return xx**2 * p2 * mp.e1(a * c) + mp.mpf(3) / 2 * a * mp.gammainc(-1, a * c)
        return mp.quad(inner, [0, 0.5, 1])
    # predicted singular part
    E = R(3, 2) * R(1, 4)          # int (N1/2) min = 3/2 * 1/4
    L = -R(1, 3) + R(3, 2) * R(1, 6)  # -int x^2 + int (3/2) x(1-x)
    vals = []
    for e in (mp.mpf('1e-4'), mp.mpf('1e-5'), mp.mpf('1e-6')):
        vals.append(F(e) - mp.mpf(E.p) / E.q / e - mp.mpf(L.p) / L.q * mp.log(e))
    return [mp.nstr(v, 12) for v in vals], str(E), str(L)


def main():
    out = {}
    ok, n = brute_vertex_check()
    out['vertex_functions_vs_brute_force'] = {'ok': ok, 'trials': n}
    out['gamma_trace_identity_ok'] = check_and_trace_identity()
    out['extraction_check_residual_constant'] = check_extraction()

    unit = TRAD00  # coefficient c in  c * hbar tr(ad ad)/((X,Y)16pi^2) * (structure)  <-> value/(TRAD00 * structure/kappa)

    # ---- A|_V two-point, external colour 0, momenta p, -p
    AV = {}
    Eav = {}
    Lav = {}
    per_graph_AV = {}
    for M in V:
        for N in V:
            g, (E, L) = two_point((0, M), (0, N))
            Eav[(M, N)] = E
            Lav[(M, N)] = L
            if (M, N) == (0, 0) or (M, N) == (0, 1):
                per_graph_AV[(M, N)] = g
    # structures/kappa: B1 -> 2(p^2 d - p p), B2 -> 2 p p, B3 -> 2 d
    cB3 = sp.simplify(Eav[(0, 0)] / (unit * 2))
    # L_MN = unit*[cB1*2(p2 d - pM pN) + cB2*2 pM pN]
    cB1, cB2 = sp.symbols('cB1 cB2')
    eqs = []
    for M in V:
        for N in V:
            target = unit * (cB1 * 2 * (P2 * int(M == N) - ps[M] * ps[N]) + cB2 * 2 * ps[M] * ps[N])
            diff = sp.expand(Lav[(M, N)] - target)
            eqs += sp.Poly(diff, *ps).coeffs()
            # eps^-1 must be proportional to delta
            assert sp.simplify(Eav[(M, N)] - unit * 2 * cB3 * int(M == N)) == 0
    solAV = sp.solve(eqs, [cB1, cB2], dict=True)
    assert len(solAV) == 1
    out['A_V two-point'] = {
        'eps^-1 coefficient of int dvol (A|_V, A|_V)': str(cB3),
        'log eps coefficient of int dvol (dA|_V, dA|_V)': str(solAV[0][cB1]),
        'log eps coefficient of int dvol (d^* A|_V, d^* A|_V)': str(solAV[0][cB2]),
    }
    # per graph (A|_V): eps^-1 -> B3 ; log -> (B1,B2)
    pg = {}
    for name in per_graph_AV[(0, 0)]:
        E00, L00 = per_graph_AV[(0, 0)][name]
        _, L01 = per_graph_AV[(0, 1)][name]
        # L00 = unit*2[cB1 (p2 - p0^2) + cB2 p0^2], L01 = unit*2[-cB1 + cB2] p0 p1
        c1, c2 = sp.symbols('c1 c2')
        e1 = sp.Poly(sp.expand(L00 - unit * 2 * (c1 * (P2 - ps[0]**2) + c2 * ps[0]**2)), *ps).coeffs()
        e2 = sp.Poly(sp.expand(L01 - unit * 2 * (-c1 + c2) * ps[0] * ps[1]), *ps).coeffs()
        sol = sp.solve(e1 + e2, [c1, c2], dict=True)
        sol = sol[0] if sol else {c1: 0, c2: 0}
        pg[name] = {'eps^-1 (A|_V,A|_V)': str(sp.simplify(E00 / (unit * 2))),
                    'log (dA|_V,dA|_V)': str(sol.get(c1, 0)), 'log (d^*A|_V,d^*A|_V)': str(sol.get(c2, 0))}
    out['A_V per graph'] = pg

    # ---- A|_{W^perp} two-point
    g44, (E44, L44) = two_point((0, 4), (0, 4))
    g45, (E45, L45) = two_point((0, 4), (0, 5))
    cB5 = sp.simplify(E44 / (unit * 2))
    cB4 = sp.simplify(sp.expand(L44 / (unit * 2 * P2)))
    assert sp.simplify(E45) == 0 and sp.expand(L45) == 0
    assert sp.expand(L44 - unit * 2 * cB4 * P2) == 0
    out['A_Wperp two-point'] = {
        'log eps coefficient of int dvol (dA|_{W^perp}, dA|_{W^perp})': str(cB4),
        'eps^-1 coefficient of int dvol (A|_{W^perp}, A|_{W^perp})': str(cB5),
        'off-diagonal (iota_w4 A, iota_w5 A) external': 'zero',
    }
    pgw = {}
    for name, (E, L) in g44.items():
        pgw[name] = {'eps^-1 (A|_W,A|_W)': str(sp.simplify(E / (unit * 2))),
                     'log (dA|_W,dA|_W)': str(sp.simplify(sp.expand(L / (unit * 2 * P2))))}
    out['A_Wperp per graph'] = pgw

    # ---- mixed and off-colour checks
    _, (Emx, Lmx) = two_point((0, 0), (0, 4))
    _, (Eoc, Loc) = two_point((0, 0), (1, 0))
    out['checks'] = {
        'A|_V - A|_W mixing zero': bool(sp.simplify(Emx) == 0 and sp.expand(Lmx) == 0),
        'off-colour zero': bool(sp.simplify(Eoc) == 0 and sp.expand(Loc) == 0),
        'log part transverse (d^* structure coefficient 0)': bool(solAV[0][cB2] == 0),
        'log: A|_V coefficient = A|_W coefficient': bool(sp.simplify(solAV[0][cB1] - cB4) == 0),
        'eps^-1: A|_V coefficient = A|_W coefficient': bool(sp.simplify(cB3 - cB5) == 0),
        'two-point coefficients independent of h': bool(not any(sp.sympify(c).has(h) for c in
                                                             (solAV[0][cB1], solAV[0][cB2], cB3, cB4, cB5))),
    }

    # ---- insertion
    def e(i):
        v = [0] * 10
        v[i] = 1
        return v
    resA = insertion((0, 4), (0, 4), e(4), e(4))   # 2 c_O + 2 c_tr
    resB = insertion((0, 5), (0, 5), e(4), e(4))   # 2 c_tr
    resC = insertion((0, 4), (0, 5), e(4), e(5))   # c_O
    resD = insertion((0, 0), (0, 0), e(4), e(4))   # 2 c_V
    resE = insertion((0, 0), (0, 4), e(4), e(4))   # mixed: must vanish
    resF = insertion((0, 6), (0, 7), e(4), e(4))   # must vanish
    cO = sp.simplify(resC['total log'] / unit)
    ctr = sp.simplify(resB['total log'] / (2 * unit))
    cV = sp.simplify(resD['total log'] / (2 * unit))
    assert sp.simplify(resA['total log'] - unit * (2 * cO + 2 * ctr)) == 0
    assert sp.simplify(resE['total log']) == 0 and sp.simplify(resF['total log']) == 0
    cvac = sp.simplify(resA['vacuum eps^-1'] / DIMG)  # per dim g * (w1,w2), units hbar/16pi^2
    trace_coeff = sp.simplify(cO + 6 * ctr)
    hv = R(1, 2)
    out['insertion (general h)'] = {
        'c_O (coefficient of O itself)': str(cO),
        'c_trace_mix (coefficient of (w1,w2) int dvol sum_w (iota_w A, iota_w A))': str(ctr),
        'c_AV_mix (coefficient of (w1,w2) int dvol (A|_V,A|_V))': str(cV),
        'vacuum eps^-1 per dim(g)(w1,w2) [hbar/16pi^2]': str(cvac),
    }
    out['insertion (h=1/2, Hodge)'] = {
        'traceless: log eps coefficient (multiplying O)': str(cO.subs(h, hv)),
        'trace sum_w O_ww: log eps coefficient (multiplying the trace)': str(trace_coeff.subs(h, hv)),
        'trace-mixing coefficient': str(ctr.subs(h, hv)),
        'A|_V mixing coefficient': str(cV.subs(h, hv)),
        'vacuum eps^-1 per dim(g)(w1,w2)': str(cvac.subs(h, hv)),
    }
    out['insertion per graph (h=1/2)'] = {
        'O + quartic (A|_W lines), ext (iota_w4 A,iota_w5 A), w1=e4,w2=e5 -> c_O part': str(sp.simplify(resC['O + quartic: log'] / unit).subs(h, hv)),
        'triangle, ext (iota_w4 A,iota_w5 A), w1=e4,w2=e5 -> c_O part': str(sp.simplify(resC['triangle: log'] / unit).subs(h, hv)),
        'O + quartic, ext (iota_w5 A)^2, w1=w2=e4 -> trace-mix': str(sp.simplify(resB['O + quartic: log'] / (2 * unit)).subs(h, hv)),
        'triangle, ext (iota_w5 A)^2, w1=w2=e4 -> trace-mix': str(sp.simplify(resB['triangle: log'] / (2 * unit)).subs(h, hv)),
        'O + quartic, ext (A|_V)^2 -> A_V mix': str(sp.simplify(resD['O + quartic: log'] / (2 * unit)).subs(h, hv)),
        'triangle, ext (A|_V)^2 -> A_V mix': str(sp.simplify(resD['triangle: log'] / (2 * unit)).subs(h, hv)),
    }
    print(json.dumps(out, indent=1, default=str))
    with open(sys.argv[1] if len(sys.argv) > 1 else 'audit_numbers_results.json', 'w') as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == '__main__':
    main()
