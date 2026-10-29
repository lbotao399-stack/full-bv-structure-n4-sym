"""Refuter #2 (Part C): the five I-a classes of (c, lambda, lambda), own derivation, explicit spinors.

Legs c(p0, r0 X0), lambda(p1, r1 s1 X1), lambda(p2, r2 s2 X2), s_j in Sigma (16 components).
G5 = int([c,lambda], dslash lambda) -> -i f (s1, rho(p0) s2) r0 r1 r2  (checked by hand), f = ([X0,X1],X2).
Realisation: Sigma = Sigma^* = C^16 with the identity pairing, rho(v) = sum v_M sig_M on Sigma (to Sigma^*),
sum v_M sigb_M on Sigma^*; sig_M = sigb_M real symmetric (M<9), sig_9 = i, sigb_9 = -i (checked Clifford).
V = directions 0..3, Wperp = 4..9.
Vertex sequences: A c A^*;  c lambda lambda^*;  lambda_x lambda_y A with fibre B(Q1X,Q2Y,vZ)=(Q1,rho(v)Q2)(X,[Z,Y]).
Koszul results (derivation in refute_report.txt); kappa = -int_0^eps |q|^2 e^{-t1 q^2}:
 (l1) AcA*[c] A_V llA[lam] K[lam] cll*[lam] c-A*    claimed -L G5
      D = B(kappa[c_1,lam_a], lam_b, A_2), V0 = A_0 c0 A^*_0; (c_1,A^*_0), (A_2,A_0); sign (a,b)=(1,2): -1
      value/(r0r1r2) = -2i [f(1,2) - f(2,1)], f = int |p+p_a|^2 (s_a, rho(p) s_b) C ; exps t1|p+p_a|^2,t2|p|^2,t3|p-p0|^2
      C = sum([T_g,X_a],[T_d,X_b])([T_d,X0],T_g)
 (l2/l3) AcA*[c] A llA[lam] lam llA[lam] K[A]  (A = A_V: claimed (-1/2 - L) G5 ; A = A_W: (-3/2 - 3L) G5)
      D = -(A_0 c0)(lam_2i lam_b) B(., ., kappa[A_0,c0]), V1 = lam_1i lam_a A_1; (lam_2i,lam_1i),(A_0,A_1); sign -1
      value = 2i [f(1,2) - f(2,1)], f = int |p+p_b|^2 (s_b^T M(p) s_a) C, M(p) = sum_nu sig(e_nu) sigb(p) sig(e_nu),
      exps t1|p+p_b|^2, t2|p|^2, t3|p-p_a|^2;  C = sum (T_bb,[[T_aa,X0],X_b])(T_bb,[T_aa,X_a])
 (l4/l5) cll*[c] K[lam] llA[lam] A llA[lam] lam  (A_V: claimed (1/2 - L) G5 ; A_W: (3/2 - 3L) G5)
      D = B(kappa[c0,lam_0], lam_b, A_2), V1 = lam_1i lam_a A_1; (lam_0,lam_1i),(A_2,A_1); sign +1
      value = -2i [f(1,2) - f(2,1)], f = int |p+p0|^2 (s_b^T M(p) s_a) C, exps t1|p+p0|^2, t2|p|^2, t3|p-p_a|^2;
      C = sum ([X0,T_aa],[T_d,X_b])(T_aa,[T_d,X_a])
 (P_cA* = -i p s, P_AA = 2 s, P_lamlam = i rho(p)|_{Sigma^*->Sigma} s ; prefactors above include kappa's -1.)
 /unit: divide by -2.
"""
import sys, os, random, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
from ghost_classes import (u, t1, t2, t3, w, s, I, cross, integrate_lines, cmul, S, TB)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def build_sigma():
    Id = sp.eye(2); X = sp.Matrix([[0, 1], [1, 0]]); Z = sp.Matrix([[1, 0], [0, -1]]); E = sp.Matrix([[0, 1], [-1, 0]])
    mats = {'I': Id, 'X': X, 'Z': Z, 'E': E}
    from sympy import kronecker_product as kp
    cands = []
    for w3 in itertools.product('IXZE', repeat=3):
        if w3.count('E') % 2 == 1:
            cands.append(kp(kp(mats[w3[0]], mats[w3[1]]), mats[w3[2]]))
    def anti(A, B):
        return (A * B + B * A) == sp.zeros(8, 8)
    best = []
    def bt(start, cur):
        nonlocal best
        if len(cur) == 7:
            best = list(cur); return True
        for i in range(start, len(cands)):
            if all(anti(cands[i], c) for c in cur):
                cur.append(cands[i])
                if bt(i + 1, cur):
                    return True
                cur.pop()
        return False
    assert bt(0, [])
    Z8 = sp.zeros(8, 8); I8 = sp.eye(8)
    G = [sp.BlockMatrix([[Z8, Ei], [Ei.T, Z8]]).as_explicit() for Ei in best]
    G.append(sp.BlockMatrix([[Z8, I8], [I8, Z8]]).as_explicit())
    G.append(sp.BlockMatrix([[I8, Z8], [Z8, -I8]]).as_explicit())
    sig = G + [I * sp.eye(16)]
    sigb = G + [-I * sp.eye(16)]
    for M in range(10):
        assert sig[M] == sig[M].T and sigb[M] == sigb[M].T
        for N in range(10):
            assert sigb[M] * sig[N] + sigb[N] * sig[M] == (2 if M == N else 0) * sp.eye(16)
    return sig, sigb


SIG, SIGB = build_sigma()


def sigv(v, bar=False):
    L = SIGB if bar else SIG
    out = sp.zeros(16, 16)
    for M, x in enumerate(v):
        if x != 0:
            out += x * L[M]
    return out


def Mmat(p4, dirs):
    """sum_nu sig(e_nu) sigb(p) sig(e_nu), p in V (4 components symbolic)"""
    sb = sigv(list(p4) + [0] * 6, bar=True)
    out = sp.zeros(16, 16)
    for nu in dirs:
        out += SIG[nu] * sb * SIG[nu]
    return out


def bil(sa, Mx, sb):
    return sp.expand((sp.Matrix(sa).T * Mx * sp.Matrix(sb))[0, 0])


def f_l1(pa, Xa, sa, pb, Xb, sb, p0, X0):
    p = list(u); spa = [s * x for x in pa]; sp0 = [s * x for x in p0]
    C = sum(dot(cross(TB[g], Xa), cross(TB[d], Xb)) * dot(cross(TB[d], X0), TB[g]) for g in range(3) for d in range(3))
    q = [p[i] + spa[i] for i in range(4)]
    spin = bil(sa, sigv(p + [0] * 6), sb)
    F = sp.expand(dot(q, q) * spin * C)
    return integrate_lines(F, [(t1, spa), (t2, [0] * 4), (t3, [-x for x in sp0])], 3, orders=(0, 1))


def f_l23(pa, Xa, sa, pb, Xb, sb, p0, X0, dirs):
    p = list(u); spa = [s * x for x in pa]; spb = [s * x for x in pb]
    C = sum(dot(TB[bb], cross(cross(TB[aa], X0), Xb)) * dot(TB[bb], cross(TB[aa], Xa)) for aa in range(3) for bb in range(3))
    q = [p[i] + spb[i] for i in range(4)]
    spin = bil(sb, Mmat(p, dirs), sa)
    F = sp.expand(dot(q, q) * spin * C)
    return integrate_lines(F, [(t1, spb), (t2, [0] * 4), (t3, [-x for x in spa])], 3, orders=(0, 1))


def f_l45(pa, Xa, sa, pb, Xb, sb, p0, X0, dirs):
    p = list(u); spa = [s * x for x in pa]; sp0 = [s * x for x in p0]
    C = sum(dot(cross(X0, TB[aa]), cross(TB[d], Xb)) * dot(TB[aa], cross(TB[d], Xa)) for aa in range(3) for d in range(3))
    q = [p[i] + sp0[i] for i in range(4)]
    spin = bil(sb, Mmat(p, dirs), sa)
    F = sp.expand(dot(q, q) * spin * C)
    return integrate_lines(F, [(t1, sp0), (t2, [0] * 4), (t3, [-x for x in spa])], 3, orders=(0, 1))


def combine(fA, fB, pref):
    z = (pref[0] / -2, pref[1] / -2)
    return {k: cmul(z, [fA[k][i] - fB[k][i] for i in range(4)]) for k in fA}


def rnd(n, lo=-2, hi=2):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]


def to_c(z):
    """sympy Gaussian rational -> (Fraction re, Fraction im)"""
    re_, im_ = sp.re(z), sp.im(z)
    return (Fr(int(sp.numer(re_)), int(sp.denom(re_))), Fr(int(sp.numer(im_)), int(sp.denom(im_))))


def cdiv(a, b):
    # a, b complex as (re, im) Fractions
    den = b[0] * b[0] + b[1] * b[1]
    return ((a[0] * b[0] + a[1] * b[1]) / den, (a[1] * b[0] - a[0] * b[1]) / den)


if __name__ == '__main__':
    random.seed(53)
    out = {}
    V = list(range(4)); W = list(range(4, 10))
    jobs = [('l1', lambda A, B: combine(f_l1(*A), f_l1(*B), (Fr(0), Fr(-2))), ('0', '-1')),
            ('l2_V', lambda A, B: combine(f_l23(*A, V), f_l23(*B, V), (Fr(0), Fr(2))), ('-1/2', '-1')),
            ('l3_W', lambda A, B: combine(f_l23(*A, W), f_l23(*B, W), (Fr(0), Fr(2))), ('-3/2', '-3')),
            ('l4_V', lambda A, B: combine(f_l45(*A, V), f_l45(*B, V), (Fr(0), Fr(-2))), ('1/2', '-1')),
            ('l5_W', lambda A, B: combine(f_l45(*A, W), f_l45(*B, W), (Fr(0), Fr(-2))), ('3/2', '-3'))]
    for name, fn, claim in jobs:
        rows = []
        while len(rows) < 3:
            p0 = rnd(4); p1 = rnd(4); p2 = [-p0[i] - p1[i] for i in range(4)]
            X0 = rnd(3); X1 = rnd(3); X2 = rnd(3)
            s1 = [sp.Integer(random.randint(-2, 2)) + I * random.randint(-2, 2) for _ in range(16)]
            s2 = [sp.Integer(random.randint(-2, 2)) + I * random.randint(-2, 2) for _ in range(16)]
            fv = dot(cross(X0, X1), X2)
            g5 = -I * sp.Rational(fv.numerator, fv.denominator) * bil(s1, sigv(S(p0) + [0] * 6), s2)  # / r0r1r2
            if g5 == 0:
                continue
            P0, P1, P2, Y0, Y1, Y2 = S(p0), S(p1), S(p2), S(X0), S(X1), S(X2)
            A = (P1, Y1, s1, P2, Y2, s2, P0, Y0)
            B = (P2, Y2, s2, P1, Y1, s1, P0, Y0)
            v = fn(A, B)
            e0 = v[1]  # eps^0 (order s^1)
            g = to_c(sp.expand(g5))
            rat = cdiv((e0[0], e0[2]), g); lg = cdiv((e0[1], e0[3]), g)
            rows.append(dict(coef_rat=[str(rat[0]), str(rat[1])], coef_log=[str(lg[0]), str(lg[1])],
                             eps_inverse=[str(x) for x in v[0]]))
        out[name] = dict(claimed_rat_log=claim, samples=rows)
        print(name, json.dumps(out[name]), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'clamlam_classes.json'), 'w'), indent=1)
