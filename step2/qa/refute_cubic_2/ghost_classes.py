"""Refuter #2 (Part C): ghost-line I-a classes of (c, A|_V, A|_V), own derivation.

Rules used (manuscript ch01/ch03 data only):
 * defect cross term {X, Y}_{K_eps-K_0}, X = ([A,c],A^*), Y a functional of A:
   = (d_{A^*}X) kappa (d_A Y) = -( kappa[A,c] , d_A Y ),  kappa = e^{-eps H}-1,
   d_{A^*}X = -[A,c] (left derivative; A^* moved across the odd c). The K_0 version is the
   gauge variation fixed by the classical master equation.
 * edge rule (Prop c3-graphs): contracted coordinate functionals brought adjacent by Koszul
   reordering, replaced by the P(eps,infty) component in that order:
   P_{c (x) A^*|_V}(p) = -(,)^{-1} (x) i p e^{-eps p^2}/p^2   (c at x, incoming momentum p there),
   P_{A (x) A}(p) = 2 (,)^{-1} (x) metric^{-1} e^{-eps p^2}/p^2.
   (checked: (Q (x) 1 + 1 (x) Q) P = K_eps on the A (x) A^* slot, refute_report.txt)
 * odd external leg c(p0) = r0 X0 e^{i p0 x}; value = coefficient of r0 moved to the left.

Koszul signs (derived in refute_report.txt):
 bubble (b) [Step-1 calibration]: D = -([kappa[A,c_i],c0],A^*), self-contraction (c_i,A^*): +1
    Gamma_b = + r0 int kappa(p+pA) (-i (a,p)) s(p) sum_b ([[X,T_b],X0],T_b)
 class (i)  AAA[A_V] A_V AcA*[c] c-A* AcA*[A_V] K[A_V]: D = -( kappa[A2,c_i], d_A I_AAA ),
            V0 = ([A,c0],A^*_0); contractions (c_i, A^*_0), (A,A): +1
    Gamma_i = + r0 int kappa(p+p2) (-i p_mu s(p)) (2 s(p-p0)) sum V3(p+p2,a2,[X2,T_b]; p1,a1,X1; -(p-p0),e_nu,T_d)
              ([T_d,X0],T_b) delta_{mu nu}
 class (ii) AcA*[A_V] K[A_V] AcA*[c] c-A* AcA*[A_V] c-A*: D = -([kappa[A1,c1],c0],A^*_0),
            V2 = ([A2,c2],A^*_2); contractions (c1,A^*_2), (c2,A^*_0): reorder +1, D sign -1
    Gamma_ii = - r0 int kappa(p+p1) (-i (a2,p) s(p)) (-i (a1,p-p2) s(p-p2)) C2
 kappa(q) = -int_0^eps |q|^2 e^{-t1 |q|^2}, s(p) = int_eps^oo e^{-t p^2}.
Own code: imports only masters.py of this folder.
"""
import sys, os, random, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sympy as sp
from fractions import Fraction as Fr
from math import comb
from masters import J

u = sp.symbols('u0:4')
t1, t2, t3, w, s = sp.symbols('t1 t2 t3 w s')
I = sp.I


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(X, Y):
    return [X[1] * Y[2] - X[2] * Y[1], X[2] * Y[0] - X[0] * Y[2], X[0] * Y[1] - X[1] * Y[0]]


def f3(X1, X2, X3):
    return dot(X1, cross(X2, X3))


def kin3(P1, a1, P2, a2, P3, a3):
    return (dot(a1, a2) * dot(a3, [P2[i] - P1[i] for i in range(4)])
            + dot(a2, a3) * dot(a1, [P3[i] - P2[i] for i in range(4)])
            + dot(a3, a1) * dot(a2, [P1[i] - P3[i] for i in range(4)]))


def V3(P1, a1, X1, P2, a2, X2, P3, a3, X3):
    return -I / 2 * f3(X1, X2, X3) * kin3(P1, a1, P2, a2, P3, a3)


TB = [[1 if i == j else 0 for j in range(3)] for i in range(3)]
EV = [[1 if i == j else 0 for j in range(4)] for i in range(4)]


def K2(a, b, n):
    """int_0^1 dtau1 int_1^oo dtau2 tau1^a tau2^b (tau1+tau2)^-n -> (rational, coeff of log 2)"""
    assert a + b - n <= -2
    rat = Fr(0); lg = Fr(0)
    for j in range(b + 1):
        m = j - n
        assert m <= -2
        cf = Fr(comb(b, j)) * Fr(-1) ** (b - j) * Fr(-1, m + 1)
        e1 = a + b - j
        for k in range(e1 + 1):
            c2 = cf * comb(e1, k) * Fr(-1) ** (e1 - k)
            e = k + m + 1
            if e == -1:
                lg += c2
            else:
                rat += c2 * (Fr(2) ** (e + 1) - 1) / (e + 1)
    return rat, lg


def gauss_avg(F, center, T_inv):
    G = sp.expand(F.subs({u[i]: u[i] + center[i] for i in range(4)}, simultaneous=True))
    P = sp.Poly(G, *u)
    out = 0
    for mon, cf in P.terms():
        if any(m % 2 for m in mon):
            continue
        val = 1
        for m in mon:
            df = 1
            for j in range(m - 1, 0, -2):
                df *= j
            val *= df * (T_inv / 2) ** (m // 2)
        out += cf * val
    return sp.expand(out)


def integrate_lines(F, lines, nlines, orders):
    """int dtau T^-2 e^{-E} <F>, Taylor coefficient of s^order; returns complex value as
    dict order -> [re_rat, re_log, im_rat, im_log] (log = log(3/2) for 3 lines, log 2 for 2 lines)."""
    center = [-w * sum(tt * o[i] for tt, o in lines) for i in range(4)]
    avg = gauss_avg(F, center, w)
    E = 0
    for i in range(len(lines)):
        for j in range(i + 1, len(lines)):
            ti, oi = lines[i]; tj, oj = lines[j]
            d = [oi[k] - oj[k] for k in range(4)]
            E += ti * tj * dot(d, d)
    E = sp.expand(w * E)
    integrand = sp.expand(w**2 * avg * (1 - E + E**2 / 2))
    res = {}
    for order in orders:
        c = sp.expand(integrand.coeff(s, order))
        acc = [Fr(0)] * 4
        if c != 0:
            P = sp.Poly(c, t1, t2, t3, w)
            for (a, b, cc, n), cf in P.terms():
                if cf == 0:
                    continue
                if nlines == 3:
                    r, L = J(a, b, cc, n)
                else:
                    assert cc == 0
                    r, L = K2(a, b, n)
                re_, im_ = sp.re(cf), sp.im(cf)
                if re_ != 0:
                    v = Fr(int(sp.numer(re_)), int(sp.denom(re_)))
                    acc[0] += v * r; acc[1] += v * L
                if im_ != 0:
                    v = Fr(int(sp.numer(im_)), int(sp.denom(im_)))
                    acc[2] += v * r; acc[3] += v * L
        res[order] = acc
    return res


def cmul(z, acc):
    """multiply complex number z = (x + i y), x,y Fractions, by acc=[re_rat,re_log,im_rat,im_log]"""
    x, y = z
    return [x * acc[0] - y * acc[2], x * acc[1] - y * acc[3], x * acc[2] + y * acc[0], x * acc[3] + y * acc[1]]


def S(L):
    return [sp.Rational(x.numerator, x.denominator) for x in L]


def rnd(n, lo=-3, hi=3):
    return [Fr(random.randint(lo, hi)) for _ in range(n)]


# ---------------- bubble (b): legs c(p0,X0), A(pA,a,X) -----------------
def bubble_b(p0, a, X0, X):
    pA = [-x for x in p0]
    spA = [s * x for x in pA]
    p = list(u)
    q = [p[i] + spA[i] for i in range(4)]
    Cb = sum(dot(cross(cross(X, TB[b]), X0), TB[b]) for b in range(3))
    # Gamma_b = r0 * (-1)(-i) int |q|^2 (a,p) Cb e^{...}  = i r0 int ...
    F = sp.expand(dot(q, q) * dot(a, p) * Cb)
    res = integrate_lines(F, [(t1, spA), (t2, [0] * 4)], 2, orders=(0, 1, 2, 3))
    # Gamma/unit = (i / (-2)) * int
    return {k: cmul((Fr(0), Fr(-1, 2)), v) for k, v in res.items()}


# ---------------- class (i) ----------------
def class_i_graph(p0, p1, a1, X0, X1, p2, a2, X2):
    sp0 = [s * x for x in p0]; sp1 = [s * x for x in p1]; sp2 = [s * x for x in p2]
    p = list(u)
    q = [p[i] + sp2[i] for i in range(4)]
    r = [p[i] - sp0[i] for i in range(4)]
    mr = [-x for x in r]
    F = 0
    for b in range(3):
        for d in range(3):
            col = dot(cross(TB[d], X0), TB[b])  # ([T_d, X0], T_b)
            if col == 0:
                continue
            for nu in range(4):
                F += V3(q, a2, cross(X2, TB[b]), sp1, a1, X1, mr, EV[nu], TB[d]) * col * p[nu]
    F = sp.expand(dot(q, q) * F)
    res = integrate_lines(F, [(t1, sp2), (t2, [0] * 4), (t3, [-x for x in sp0])], 3, orders=(0, 1, 2))
    # Gamma_i = r0 (-1)(-i)(2) int ... = 2i r0 int ;  /unit: (2i)/(-2) = -i
    return {k: cmul((Fr(0), Fr(-1)), v) for k, v in res.items()}


# ---------------- class (ii) ----------------
def class_ii_graph(p0, p1, a1, X0, X1, p2, a2, X2):
    sp1 = [s * x for x in p1]; sp2 = [s * x for x in p2]
    p = list(u)
    q = [p[i] + sp1[i] for i in range(4)]
    pp = [p[i] - sp2[i] for i in range(4)]
    C2 = 0
    for b in range(3):
        for g in range(3):
            C2 += dot(cross(cross(X1, TB[b]), X0), TB[g]) * dot(cross(X2, TB[g]), TB[b])
    F = sp.expand(dot(q, q) * dot(a1, pp) * dot(a2, p) * C2)
    res = integrate_lines(F, [(t1, sp1), (t2, [0] * 4), (t3, [-x for x in sp2])], 3, orders=(0, 1, 2))
    # Gamma_ii = -r0 * (-1) * (-i)(-i) int = -r0 * (-1)*(-1) int = -r0 int ;  /unit: (-1)/(-2) = 1/2
    return {k: cmul((Fr(1, 2), Fr(0)), v) for k, v in res.items()}


def G1(p0, p1, a1, X0, X1, p2, a2, X2):
    f = dot(cross(X0, X1), X2)
    return f * ((dot(p2, p2) - dot(p1, p1)) * dot(a1, a2) - dot(a1, p2) * dot(a2, p2) + dot(a1, p1) * dot(a2, p1))


def G2(p0, p1, a1, X0, X1, p2, a2, X2):
    f = dot(cross(X0, X1), X2)
    return f * (dot(a1, p2) * dot(a2, p2) - dot(a1, p1) * dot(a2, p1))


def solve2(b, A):
    det = A[0][0] * A[1][1] - A[0][1] * A[1][0]
    al = (b[0] * A[1][1] - b[1] * A[0][1]) / det
    be = (A[0][0] * b[1] - A[1][0] * b[0]) / det
    res = [b[i] - al * A[i][0] - be * A[i][1] for i in range(len(b))]
    return al, be, res


if __name__ == '__main__':
    out = {}
    random.seed(5)
    cal = []
    for _ in range(3):
        p0 = rnd(4); a = rnd(4); X0 = rnd(3); X = rnd(3)
        res = bubble_b(S(p0), S(a), S(X0), S(X))
        pA = [-x for x in p0]
        XX = dot(X0, X)
        # functionals / r0:  int(dc,A_V) = i (p0,a)(X0,X) ;  int(d*dc, d*A_V) = -i |p0|^2 (pA,a)(X0,X)
        g_m1 = dot(p0, a) * XX          # coefficient of i
        g_0 = -dot(p0, p0) * dot(pA, a) * XX  # coefficient of i
        m1 = res[1]; z0 = res[3]
        cal.append(dict(eps_inv_over_int_dcA=str(m1[2] / g_m1) if g_m1 else None,
                        eps_inv_real_part=[str(m1[0]), str(m1[1])], eps_inv_log2=str(m1[3]),
                        eps0_over_int_dsdc_dsA=str(z0[2] / g_0) if g_0 else None,
                        eps0_log2=str(z0[3]), eps0_real=[str(z0[0]), str(z0[1])],
                        even_orders=[[str(x) for x in res[0]], [str(x) for x in res[2]]]))
    out['calibration_bubble_b'] = dict(expected_Step1=dict(eps_inv='1/8', eps0='-1/12'), samples=cal)
    print(json.dumps(out['calibration_bubble_b'], indent=1), flush=True)

    for name, fn in [('class_i', class_i_graph), ('class_ii', class_ii_graph)]:
        random.seed(17)
        rows = []
        while len(rows) < 5:
            p0 = rnd(4); p1 = rnd(4); p2 = [-p0[i] - p1[i] for i in range(4)]
            a1 = rnd(4); a2 = rnd(4); X0 = rnd(3); X1 = rnd(3); X2 = rnd(3)
            if dot(cross(X0, X1), X2) == 0:
                continue
            A1 = (S(p0), S(p1), S(a1), S(X0), S(X1), S(p2), S(a2), S(X2))
            A2 = (S(p0), S(p2), S(a2), S(X0), S(X2), S(p1), S(a1), S(X1))
            r1 = fn(*A1); r2 = fn(*A2)
            tot = {k: [r1[k][i] + r2[k][i] for i in range(4)] for k in r1}
            rows.append(dict(tot=tot, single=r1, G=[Fr(str(G1(*A1))), Fr(str(G2(*A1)))],
                             fa1a2=Fr(str(dot(cross(A1[3], A1[4]), A1[7]) * dot(A1[2], A1[6])))))
        Gm = [r['G'] for r in rows]
        al_r, be_r, res_r = solve2([r['tot'][2][0] for r in rows], Gm)
        al_l, be_l, res_l = solve2([r['tot'][2][1] for r in rows], Gm)
        out[name] = dict(
            G1_coefficient=f"{al_r} + {al_l}*log(3/2)", G2_coefficient=f"{be_r} + {be_l}*log(3/2)",
            residuals=[str(x) for x in res_r + res_l],
            imaginary_parts_eps0=[[str(r['tot'][2][2]), str(r['tot'][2][3])] for r in rows],
            eps_inverse_class_sum=[[str(x) for x in r['tot'][0]] for r in rows],
            eps_inverse_single_over_f_a1a2=[str(r['single'][0][0] / r['fa1a2']) if r['fa1a2'] else None for r in rows],
            eps_inverse_single_log=[str(r['single'][0][1]) for r in rows],
            odd_order_s1=[[str(x) for x in r['tot'][1]] for r in rows])
        print(name, json.dumps(out[name], indent=1), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ghost_classes.json'), 'w'), indent=1)
