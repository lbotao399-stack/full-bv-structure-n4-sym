"""Refuter #2 (Part C): independent recomputation of the I-a triangle classes of (c, A|_V, A|_V)
whose cycle is  c-vertex ([A,c],A^*) -- (A)^3 -- (A)^3  with only gauge lines:

  class V : AAA[A_V] A_V AAA[A_V] A_V AcA*[c] K[A_V]   (internal lines A|_V)
            claimed eps^0 = (103/144 + 2/3 L) G1 + (41/48 + 1/4 L) G2,  L = log(3/2)
  class W : AAA[A_V] A_W AAA[A_V] A_W AcA*[c] K[A_W]   (internal lines A|_Wperp)
            claimed eps^0 = (11/12 + L) G1 + 5/4 G2

Hand-derived formula (derivation in refute_report.txt):
 defect vertex cross term {I_cAA*, I_AAA}_{K_eps-K_0} = -(derivative of I_AAA along
 A -> A + (e^{-eps H}-1)[A,c]) (the K_0 rule {I_3,X(A)} = -d/dt X(A+t[A,c]) is the one forced by
 the classical master equation Q I_4 + {I_cAA*,I_AAA} = 0, i.e. gauge invariance of -1/4(F,F)).
 Graph (legs c(p0,X0), A(p1,a1,X1) at V1, A(p2,a2,X2) at V2; weight 1):
   Gamma = -r0 hbar int d^4k/(2pi)^4 (e^{-eps|q|^2}-1) * 2e^{-eps k^2}/k^2 * 2e^{-eps l^2}/l^2 * N(k)
   q = k+p0 (special edge), k (line V0-V1), l = k-p1 (line V1-V2),
   N = sum_{e,e'} V3(q, e, [T_a,X0]; p2,a2,X2; -l, e', T_b) * V3(p1,a1,X1; -k, e, T_a; l, e', T_b),
 P(eps,infty)_{AA} = 2 (,)^-1 (x) metric^-1 e^{-eps p^2}/p^2 ; (e^{-eps q^2}-1) = -int_0^eps q^2 e^{-t1 q^2}.
 => Gamma = 4 r0 hbar/(16 pi^2) int_{[0,eps]x[eps,oo)^2} dt T^-2 e^{-E} < q^2 N >_{k0, 1/(2T)}
   T = t1+t2+t3, k0 = (t3 p1 - t1 p0)/T, E = (t1 t2 p0^2 + t2 t3 p1^2 + t1 t3 p2^2)/T.
 Scaling t = eps tau: order-s^j part (p -> s p) carries eps^{j/2-1}; eps^-1: s^0, eps^0: s^2.
 Colour: g = su(2), orthonormal basis, [X,Y] = X x Y, tr(ad_X ad_Y)/(X,Y) = -2, so
 Gamma / unit = 4/(-2) * int(...) per unit r0.
Own code: imports only masters.py of this folder.
"""
import sys, os, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from fractions import Fraction as Fr
from masters import J

u = sp.symbols('u0:4')
K = sp.symbols('K0:4')
t1, t2, t3, w, s = sp.symbols('t1 t2 t3 w s')
DIM = 10


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(X, Y):
    return [X[1] * Y[2] - X[2] * Y[1], X[2] * Y[0] - X[0] * Y[2], X[0] * Y[1] - X[1] * Y[0]]


def f3(X1, X2, X3):
    return dot(X1, cross(X2, X3))


def emb(p4):
    return list(p4) + [0] * (DIM - 4)


def kin3(P1, a1, P2, a2, P3, a3):
    # bracket of eq:c1-AAA without the colour factor and without -i/2
    return (dot(a1, a2) * dot(a3, [P2[i] - P1[i] for i in range(DIM)])
            + dot(a2, a3) * dot(a1, [P3[i] - P2[i] for i in range(DIM)])
            + dot(a3, a1) * dot(a2, [P1[i] - P3[i] for i in range(DIM)]))


def unitvec(i):
    v = [0] * DIM
    v[i] = 1
    return v


TB = [[1 if i == j else 0 for j in range(3)] for i in range(3)]


def graph(p0, p1, a1, X0, X1, p2, a2, X2, internal):
    """eps^-1 and eps^0 coefficients (each in Q + Q log(3/2)), in the unit, per r0."""
    P0 = emb([s * x for x in p0]); P1 = emb([s * x for x in p1]); P2 = emb([s * x for x in p2])
    k = emb([K[i] + u[i] for i in range(4)])
    q = [k[i] + P0[i] for i in range(DIM)]
    l = [k[i] - P1[i] for i in range(DIM)]
    mk = [-x for x in k]; ml = [-x for x in l]
    A1 = emb(a1); A2 = emb(a2)
    # colour factor: (-i/2)^2 * sum_{a,b} f([Ta,X0], X2, Tb) f(X1, Ta, Tb)
    # kinematics depends on (mu,nu) only
    col = {}
    for a in range(3):
        for b in range(3):
            col[(a, b)] = sp.Rational(-1, 4) * f3(cross(TB[a], X0), X2, TB[b]) * f3(X1, TB[a], TB[b])
    C = sum(col.values())
    N = 0
    for mu in internal:
        for nu in internal:
            e = unitvec(mu); ep = unitvec(nu)
            N += kin3(q, e, P2, A2, ml, ep) * kin3(P1, A1, mk, e, l, ep)
    F = sp.expand(dot(q, q) * N * C)
    Fp = sp.Poly(F, *u)
    # Gaussian moments with <u_i u_j> = delta_ij w/2
    avg = 0
    for mon, cf in Fp.terms():
        if any(m % 2 for m in mon):
            continue
        val = 1
        for m in mon:
            # (m-1)!! (w/2)^(m/2)
            df = 1
            for j in range(m - 1, 0, -2):
                df *= j
            val *= df * (w / 2) ** (m // 2)
        avg += cf * val
    # substitute k0
    k0 = [w * (t3 * s * p1[i] - t1 * s * p0[i]) for i in range(4)]
    avg = sp.expand(avg.subs({K[i]: k0[i] for i in range(4)}, simultaneous=True))
    Ep = t1 * t2 * dot(p0, p0) + t2 * t3 * dot(p1, p1) + t1 * t3 * dot(p2, p2)
    integrand = sp.expand(w**2 * avg * (1 - s**2 * w * Ep))
    out = {}
    for order, key in [(0, 'eps^-1'), (2, 'eps^0')]:
        c = integrand.coeff(s, order)
        P = sp.Poly(c, t1, t2, t3, w)
        rat = Fr(0); lg = Fr(0)
        for (a, b, cc, n), cf in P.terms():
            if cf == 0:
                continue
            assert a + b + cc - n == (-4 if order == 0 else -3), (a, b, cc, n)
            r, L = J(a, b, cc, n)
            cf = Fr(int(sp.numer(cf)), int(sp.denom(cf)))
            rat += cf * r; lg += cf * L
        # Gamma/unit = 4/(-2) * integral
        out[key] = (-2 * rat, -2 * lg)
    # odd orders must vanish (no eps^-1/2)
    assert sp.expand(integrand.coeff(s, 1)) == 0
    return out


def G1(p0, p1, a1, X0, X1, p2, a2, X2):
    f = dot(cross(X0, X1), X2)
    return f * ((dot(p2, p2) - dot(p1, p1)) * dot(a1, a2) - dot(a1, p2) * dot(a2, p2) + dot(a1, p1) * dot(a2, p1))


def G2(p0, p1, a1, X0, X1, p2, a2, X2):
    f = dot(cross(X0, X1), X2)
    return f * (dot(a1, p2) * dot(a2, p2) - dot(a1, p1) * dot(a2, p1))


def W2(p0, p1, a1, X0, X1, p2, a2, X2):
    f = dot(cross(X0, X1), X2)
    return f * dot(a1, a2)


def rnd(n, lo=-3, hi=3):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]


def run(internal, nconf=4, seed=11):
    random.seed(seed)
    rows = []
    for _ in range(nconf):
        p0 = rnd(4); p1 = rnd(4); p2 = [-p0[i] - p1[i] for i in range(4)]
        a1 = rnd(4); a2 = rnd(4); X0 = rnd(3); X1 = rnd(3); X2 = rnd(3)
        S = lambda L: [sp.Rational(x.numerator, x.denominator) for x in L]
        args1 = (S(p0), S(p1), S(a1), S(X0), S(X1), S(p2), S(a2), S(X2))
        args2 = (S(p0), S(p2), S(a2), S(X0), S(X2), S(p1), S(a1), S(X1))  # A1 <-> A2
        g1 = graph(*args1, internal=internal)
        g2 = graph(*args2, internal=internal)
        tot = {k: (g1[k][0] + g2[k][0], g1[k][1] + g2[k][1]) for k in g1}
        conf = dict(p0=p0, p1=p1, p2=p2, a1=a1, a2=a2, X0=X0, X1=X1, X2=X2)
        g = [Fr(str(G1(*args1))), Fr(str(G2(*args1))), Fr(str(W2(*args1)))]
        rows.append(dict(single_graph=g1, class_sum=tot, G=g, conf=conf))
    return rows


def solve2(rows, part):
    # solve class eps^0 (part 0 rational, 1 log) = alpha G1 + beta G2 on first two rows, test the rest
    import itertools
    A = [[r['G'][0], r['G'][1]] for r in rows]
    b = [r['class_sum']['eps^0'][part] for r in rows]
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    assert det != 0
    al = (b[0] * A[1][1] - b[1] * A[0][1]) / det
    be = (A[0][0] * b[1] - A[1][0] * b[0]) / det
    resid = [b[i] - al * A[i][0] - be * A[i][1] for i in range(len(rows))]
    return al, be, resid


if __name__ == '__main__':
    import time
    res = {}
    for name, internal in [('class_V_internal_A|_V', list(range(4))), ('class_W_internal_A|_Wperp', list(range(4, 10)))]:
        t = time.time()
        rows = run(internal)
        ar, br, resr = solve2(rows, 0)
        al, bl, resl = solve2(rows, 1)
        # eps^-1: single graph proportional to f (a1,a2); class sum zero
        em1_single = [(str(r['single_graph']['eps^-1'][0]), str(r['single_graph']['eps^-1'][1])) for r in rows]
        em1_single_ratio = [str(r['single_graph']['eps^-1'][0] / r['G'][2]) if r['G'][2] != 0 else None for r in rows]
        em1_class = [(str(r['class_sum']['eps^-1'][0]), str(r['class_sum']['eps^-1'][1])) for r in rows]
        res[name] = dict(
            G1_coefficient=f"{ar} + {al}*log(3/2)", G2_coefficient=f"{br} + {bl}*log(3/2)",
            residuals_rational=[str(x) for x in resr], residuals_log=[str(x) for x in resl],
            eps_inverse_single_graph=em1_single, eps_inverse_single_over_f_a1a2=em1_single_ratio,
            eps_inverse_class_sum=em1_class, seconds=round(time.time() - t, 1))
        print(name, json.dumps(res[name], indent=1), flush=True)
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'triangle_AAA.json'), 'w'), indent=1)
