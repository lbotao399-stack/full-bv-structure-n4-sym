"""Refuter #2 (Part C): the three I-a classes of (c, c, A^*|_V), own derivation.

Legs: c(p0, r0 X0), c(p1, r1 X1), A^*(p2, r2 b X2), b in V.  Value = coefficient of r0 r1 r2.
G4 = int([c,dc],A^*)  ->  -i f (p2,b) r0 r1 r2,  f = ([X0,X1],X2)  (checked by hand).
Rules as in ghost_classes.py.  Koszul bookkeeping (vertex = coordinate sequence A c A^*):
 (c1) engine word AcA*[As_V] A_V AcA*[c] K[A_V] AcA*[c] c-A*   (claimed  L G4)
      D = {V1,V0}_K = -([kappa[A_1,c_a],c_b],A^*_0); V2 = A_2 c_2 A^*_ext
      P: (A_1,A_2), (c_2,A^*_0).  graph(c_a=c1,c_b=c0): sign -1 ; exchanged (c_a=c0,c_b=c1): +1
      value = -2i [ f(a=1,b=0) - f(a=0,b=1) ],  f = int kappa-integrand |p-p_b|^2 (p,b) C,
      exponent t1|p-p_b|^2 + t2|p+p2|^2 + t3|p|^2, C = sum([[T_a,X_a],X_b],T_g)([T_a,T_g],X2)
 (c2) engine word AcA*[As_V] A_V AcA*[c] c-A* ccc*[c] K[c]       (claimed (-1/4 - L/2) G4)
      D = {V1=ccc*,V2}_K = ([A_2, kappa[c_a,c_i]], A^*_ext); V0 = A_0 c_b A^*_0
      P: (A_2,A_0), (c_i,A^*_0).  (a=1,b=0): +1 ; (a=0,b=1): -1
      value = 2i [ f(1,0) - f(0,1) ], f = int |p+p_a|^2 (b,p) C, exponent t1|p+p_a|^2+t2|p-p_b|^2+t3|p|^2,
      C = sum([T_a,[X_a,T_b]],X2)([T_a,X_b],T_b)
 (c3) engine word AcA*[As_V] K[A_V] AcA*[c] A_V AcA*[c] c-A*     (claimed (1/4 - L/2) G4)
      D = {V0,V2}_K = -([kappa[A_0,c_a],c_2],A^*_ext); V1 = A_1 c_b A^*_1
      P: (A_0,A_1), (c_2,A^*_1).  (a=0,b=1): +1 ; (a=1,b=0): -1
      value = 2i [ f(0,1) - f(1,0) ], f = int |p+p2|^2 (b,p) C, exponent t1|p+p2|^2+t2|p-p_b|^2+t3|p|^2,
      C = sum([[T_a,X_a],T_g],X2)([T_a,X_b],T_g)
 (kappa = -int_0^eps |q|^2 e^{-t1 q^2}, P_AA = 2 s, P_cA* = -i p s; the -1 of kappa, the 2 and the -i are
  in the prefactors above.)   /unit: divide by tr(ad ad)/(,) = -2.
"""
import sys, os, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from fractions import Fraction as Fr
from ghost_classes import (u, t1, t2, t3, w, s, I, cross, integrate_lines, cmul, S, rnd, TB)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def f_c1(pa, Xa, pb, Xb, p2, b, X2):
    p = list(u)
    spb = [s * x for x in pb]; sp2 = [s * x for x in p2]
    C = sum(dot(cross(cross(TB[a], Xa), Xb), TB[g]) * dot(cross(TB[a], TB[g]), X2) for a in range(3) for g in range(3))
    q = [p[i] - spb[i] for i in range(4)]
    F = sp.expand(dot(q, q) * dot(p, b) * C)
    return integrate_lines(F, [(t1, [-x for x in spb]), (t2, sp2), (t3, [0] * 4)], 3, orders=(0, 1))


def f_c2(pa, Xa, pb, Xb, p2, b, X2):
    p = list(u)
    spa = [s * x for x in pa]; spb = [s * x for x in pb]
    C = sum(dot(cross(TB[a], cross(Xa, TB[bb])), X2) * dot(cross(TB[a], Xb), TB[bb]) for a in range(3) for bb in range(3))
    q = [p[i] + spa[i] for i in range(4)]
    F = sp.expand(dot(q, q) * dot(b, p) * C)
    return integrate_lines(F, [(t1, spa), (t2, [-x for x in spb]), (t3, [0] * 4)], 3, orders=(0, 1))


def f_c3(pa, Xa, pb, Xb, p2, b, X2):
    p = list(u)
    spb = [s * x for x in pb]; sp2 = [s * x for x in p2]
    C = sum(dot(cross(cross(TB[a], Xa), TB[g]), X2) * dot(cross(TB[a], Xb), TB[g]) for a in range(3) for g in range(3))
    q = [p[i] + sp2[i] for i in range(4)]
    F = sp.expand(dot(q, q) * dot(b, p) * C)
    return integrate_lines(F, [(t1, sp2), (t2, [-x for x in spb]), (t3, [0] * 4)], 3, orders=(0, 1))


def combine(fA, fB, pref):
    # pref * (fA - fB) / (-2)
    z = (pref[0] / -2, pref[1] / -2)
    return {k: cmul(z, [fA[k][i] - fB[k][i] for i in range(4)]) for k in fA}


if __name__ == '__main__':
    random.seed(31)
    out = {}
    claims = {'c1': ('0', '1'), 'c2': ('-1/4', '-1/2'), 'c3': ('1/4', '-1/2')}
    for name in ['c1', 'c2', 'c3']:
        rows = []
        while len(rows) < 4:
            p0 = rnd(4); p1 = rnd(4); p2 = [-p0[i] - p1[i] for i in range(4)]
            b = rnd(4); X0 = rnd(3); X1 = rnd(3); X2 = rnd(3)
            fval = dot(cross(X0, X1), X2)
            if fval == 0 or dot(p2, b) == 0:
                continue
            P0, P1, P2, B, Y0, Y1, Y2 = S(p0), S(p1), S(p2), S(b), S(X0), S(X1), S(X2)
            if name == 'c1':
                v = combine(f_c1(P1, Y1, P0, Y0, P2, B, Y2), f_c1(P0, Y0, P1, Y1, P2, B, Y2), (Fr(0), Fr(-2)))
            elif name == 'c2':
                v = combine(f_c2(P1, Y1, P0, Y0, P2, B, Y2), f_c2(P0, Y0, P1, Y1, P2, B, Y2), (Fr(0), Fr(2)))
            else:
                v = combine(f_c3(P0, Y0, P1, Y1, P2, B, Y2), f_c3(P1, Y1, P0, Y0, P2, B, Y2), (Fr(0), Fr(2)))
            # G4 / (r0 r1 r2) = -i f (p2,b): purely imaginary
            g4_im = -fval * dot(p2, b)
            eps0 = v[1]  # order s^1 = eps^0 (one derivative)
            rows.append(dict(coef_rat=str(eps0[2] / g4_im), coef_log=str(eps0[3] / g4_im),
                             real_part=[str(eps0[0]), str(eps0[1])],
                             eps_inverse=[str(x) for x in v[0]]))
        out[name] = dict(claimed_rat_log=claims[name], samples=rows)
        print(name, json.dumps(out[name]), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ccA_classes.json'), 'w'), indent=1)
