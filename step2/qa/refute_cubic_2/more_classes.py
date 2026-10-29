"""Refuter #2 (Part C): further I-a classes, own derivation (same rules as ghost_classes.py).

 class (iii)  AcA*[A_V] K[c] ccc*[c] c-A* AcA*[A_V] c-A*   in (c,A_V,A_V)
    D = {1/2([c,c],c^*), ([A,c],A^*)}_{K} = ([A2, kappa[c0,c_int]], A^*_2)   (c -> c + 1/2[c,c], CME-fixed sign)
    V1 = ([A1,c1],A^*_1); contractions (c_int, A^*_1), (c1, A^*_2); Koszul sign -1
    Gamma_iii = - r0 int kappa(p+p0) (-i (a1,p) s(p)) (-i (a2,p-p1) s(p-p1)) C3,
    C3 = sum_{b,g} ([X2,[X0,T_b]],T_g)([X1,T_g],T_b)
 (c, A_W, A_W) classes (W = A|_Wperp polarisations, no coupling to momenta):
    W1 AAA[A_W] A_V AAA[A_W] A_W AcA*[c] K[A_W] : gauge-type, line V0-V1 in W, line V1-V2 in V
    W2 AAA[A_W] A_V AcA*[c] K[A_V] AAA[A_W] A_W : gauge-type, line V0-V1 in V, line V1-V2 in W
    W3 AAA[A_W] A_V AcA*[c] c-A* AcA*[A_W] K[A_W]: class (i)-type with external W polarisations
    claimed: W1 (17/48 + L/2) G3, W2 (-1/48 - L/4) G3, W3 (1/12 + L/2) G3
Own code: imports only masters.py / ghost_classes.py helpers of this folder.
"""
import sys, os, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from fractions import Fraction as Fr
from ghost_classes import (u, t1, t2, t3, w, s, I, cross, f3, integrate_lines, cmul, S, rnd, TB, solve2, G1, G2)

DIM = 10


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def emb(p4):
    return list(p4) + [0] * (DIM - len(p4))


def kin3(P1, a1, P2, a2, P3, a3):
    return (dot(a1, a2) * dot(a3, [P2[i] - P1[i] for i in range(DIM)])
            + dot(a2, a3) * dot(a1, [P3[i] - P2[i] for i in range(DIM)])
            + dot(a3, a1) * dot(a2, [P1[i] - P3[i] for i in range(DIM)]))


def V3(P1, a1, X1, P2, a2, X2, P3, a3, X3):
    return -I / 2 * f3(X1, X2, X3) * kin3(P1, a1, P2, a2, P3, a3)


def unit(i):
    v = [0] * DIM; v[i] = 1; return v


def class_iii_graph(p0, p1, a1, X0, X1, p2, a2, X2):
    sp0 = [s * x for x in p0]; sp1 = [s * x for x in p1]
    p = list(u)
    q = [p[i] + sp0[i] for i in range(4)]
    C3 = 0
    for b in range(3):
        for g in range(3):
            C3 += dot(cross(X2, cross(X0, TB[b])), TB[g]) * dot(cross(X1, TB[g]), TB[b])
    pm = [p[i] - sp1[i] for i in range(4)]
    F = sp.expand(dot(q, q) * dot(a1[:4], p) * dot(a2[:4], pm) * C3)
    res = integrate_lines(F, [(t1, sp0), (t2, [0] * 4), (t3, [-x for x in sp1])], 3, orders=(0, 1, 2))
    # Gamma = -r0 (-1)(-i)(-i) int = -r0 int ... wait: kappa = -int, (-i)^2 = -1: -r0*(-1)*(-1) = -r0
    return {k: cmul((Fr(1, 2), Fr(0)), v) for k, v in res.items()}  # /unit: (-1)/(-2)


def gauge_graph(p0, p1, a1, X0, X1, p2, a2, X2, int_e, int_ep):
    """as triangle_AAA.graph, with the line V0-V1 polarisations in int_e and V1-V2 in int_ep"""
    P0 = emb([s * x for x in p0]); P1 = emb([s * x for x in p1]); P2 = emb([s * x for x in p2])
    K = emb(list(u))
    # here loop momentum k = u (center shift handled by integrate_lines)
    k = K
    q = [k[i] + P0[i] for i in range(DIM)]
    l = [k[i] - P1[i] for i in range(DIM)]
    mk = [-x for x in k]; ml = [-x for x in l]
    A1 = emb(a1); A2 = emb(a2)
    C = 0
    for a in range(3):
        for b in range(3):
            C += sp.Rational(-1, 4) * f3(cross(TB[a], X0), X2, TB[b]) * f3(X1, TB[a], TB[b])
    N = 0
    for mu in int_e:
        for nu in int_ep:
            N += kin3(q, unit(mu), P2, A2, ml, unit(nu)) * kin3(P1, A1, mk, unit(mu), l, unit(nu))
    F = sp.expand(dot(q[:4], q[:4]) * N * C)
    res = integrate_lines(F, [(t1, P0[:4]), (t2, [0] * 4), (t3, [-x for x in P1[:4]])], 3, orders=(0, 1, 2))
    # Gamma = 4 r0 int ; /unit = -2
    return {kk: cmul((Fr(-2), Fr(0)), v) for kk, v in res.items()}


def class_i_10(p0, p1, a1, X0, X1, p2, a2, X2):
    """class (i) with 10-dimensional external polarisations a1, a2 (used with a1,a2 in W)"""
    sp0 = emb([s * x for x in p0]); sp1 = emb([s * x for x in p1]); sp2 = emb([s * x for x in p2])
    p = emb(list(u))
    q = [p[i] + sp2[i] for i in range(DIM)]
    r = [p[i] - sp0[i] for i in range(DIM)]
    mr = [-x for x in r]
    A1 = emb(a1); A2 = emb(a2)
    F = 0
    for b in range(3):
        for d in range(3):
            col = dot(cross(TB[d], X0), TB[b])
            if col == 0:
                continue
            for nu in range(4):
                F += V3(q, A2, cross(X2, TB[b]), sp1, A1, X1, mr, unit(nu), TB[d]) * col * p[nu]
    F = sp.expand(dot(q[:4], q[:4]) * F)
    res = integrate_lines(F, [(t1, sp2[:4]), (t2, [0] * 4), (t3, [-x for x in sp0[:4]])], 3, orders=(0, 1, 2))
    return {k: cmul((Fr(0), Fr(-1)), v) for k, v in res.items()}


def G3(p0, p1, a1, X0, X1, p2, a2, X2):
    f = dot(cross(X0, X1), X2)
    return f * dot(a1, a2) * (dot(p2, p2) - dot(p1, p1))


def run(fn, nconf, Gfun, extW=False, seed=23):
    random.seed(seed)
    rows = []
    while len(rows) < nconf:
        p0 = rnd(4); p1 = rnd(4); p2 = [-p0[i] - p1[i] for i in range(4)]
        if extW:
            a1 = [Fr(0)] * 4 + rnd(6); a2 = [Fr(0)] * 4 + rnd(6)
        else:
            a1 = rnd(4); a2 = rnd(4)
        X0 = rnd(3); X1 = rnd(3); X2 = rnd(3)
        if dot(cross(X0, X1), X2) == 0 or dot(a1, a2) == 0:
            continue
        A1 = (S(p0), S(p1), S(a1), S(X0), S(X1), S(p2), S(a2), S(X2))
        A2 = (S(p0), S(p2), S(a2), S(X0), S(X2), S(p1), S(a1), S(X1))
        r1 = fn(*A1); r2 = fn(*A2)
        tot = {k: [r1[k][i] + r2[k][i] for i in range(4)] for k in r1}
        rows.append(dict(tot=tot, single=r1, G=[Fr(str(g(*A1))) for g in Gfun],
                         fa1a2=Fr(str(dot(cross(A1[3], A1[4]), A1[7]) * dot(A1[2], A1[6])))))
    return rows


def summarise(rows, ng):
    out = {}
    if ng == 2:
        Gm = [r['G'] for r in rows]
        al_r, be_r, res_r = solve2([r['tot'][2][0] for r in rows], Gm)
        al_l, be_l, res_l = solve2([r['tot'][2][1] for r in rows], Gm)
        out['G1_coefficient'] = f"{al_r} + {al_l}*log(3/2)"
        out['G2_coefficient'] = f"{be_r} + {be_l}*log(3/2)"
        out['residuals'] = [str(x) for x in res_r + res_l]
    else:
        g = [r['G'][0] for r in rows]
        cr = rows[0]['tot'][2][0] / g[0]; cl = rows[0]['tot'][2][1] / g[0]
        out['G3_coefficient'] = f"{cr} + {cl}*log(3/2)"
        out['residuals'] = [str(r['tot'][2][0] - cr * r['G'][0]) for r in rows] + \
                           [str(r['tot'][2][1] - cl * r['G'][0]) for r in rows]
    out['imag_eps0'] = [[str(r['tot'][2][2]), str(r['tot'][2][3])] for r in rows]
    out['eps_inverse_class_sum'] = [[str(x) for x in r['tot'][0]] for r in rows]
    out['eps_inverse_single_over_f_a1a2'] = [str(r['single'][0][0] / r['fa1a2']) for r in rows]
    out['eps_inverse_single_log_and_imag'] = [[str(r['single'][0][1]), str(r['single'][0][2])] for r in rows]
    out['order_s1'] = [[str(x) for x in r['tot'][1]] for r in rows]
    return out


if __name__ == '__main__':
    out = {}
    W = list(range(4, 10)); V = list(range(4))
    jobs = [
        ('class_iii_ccc*_(c,A_V,A_V)  claimed (-11/288 - L/24) G1 + (-1/32 - L/8) G2', class_iii_graph, [G1, G2], False),
        ('W1 AAA[A_W] A_V AAA[A_W] A_W AcA*[c] K[A_W]  claimed (17/48 + L/2) G3',
         lambda *A: gauge_graph(*A, int_e=W, int_ep=V), [G3], True),
        ('W2 AAA[A_W] A_V AcA*[c] K[A_V] AAA[A_W] A_W  claimed (-1/48 - L/4) G3',
         lambda *A: gauge_graph(*A, int_e=V, int_ep=W), [G3], True),
        ('W3 AAA[A_W] A_V AcA*[c] c-A* AcA*[A_W] K[A_W]  claimed (1/12 + L/2) G3', class_i_10, [G3], True),
        ('check: class V via gauge_graph  claimed (103/144 + 2/3 L) G1 + (41/48 + L/4) G2',
         lambda *A: gauge_graph(*A, int_e=V, int_ep=V), [G1, G2], False),
    ]
    for name, fn, Gf, extW in jobs:
        rows = run(fn, 5, Gf, extW=extW)
        out[name] = summarise(rows, len(Gf))
        print(name, json.dumps(out[name]), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'more_classes.json'), 'w'), indent=1)
