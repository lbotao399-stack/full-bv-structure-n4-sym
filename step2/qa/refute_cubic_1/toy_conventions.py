"""Toy (finite-dimensional super vector space) check of the combinatorial conventions used by
refute_cubic_1:

 W(P,J) = hbar log( e^{hbar d_P} e^{J/hbar} ),  d_P = 1/2 sum P^{ab} d_{b_b} d_{b_a}   (left derivatives)

 Claims tested (exact arithmetic):
 (T1) one-loop formula  W_1^{1PI}(x) = sum_n 1/(2n) STr(M^n),  M^a_c = sum_b P^{ab} (-1)^{|b|} d_b d_c J
      at the generic point (coefficients = the coordinates themselves, odd ones included);
 (T2) full order-hbar part of W at cubic order = [W_1^{1PI} - C](x + chi(x)) with
      chi^a = sum_b P^{ab} (d_b I)(x + chi)   (tree field, no parity sign), C an order-hbar quadratic vertex;
 (T3) theta trick: (d_b d_c F)(e) = coefficient of theta_c theta_b in F(e + theta_b b_b + theta_c b_c).
"""
import random
from fractions import Fraction as Fr
from superalg import Alg

random.seed(7)
A = Alg()
NE, NO = 2, 4
coords = list(range(NE + NO))
for i in coords:
    A.gen(i, 0 if i < NE else 1)
H, U, G = 20, 21, 22
for g in (H, U, G):
    A.gen(g, 0)
par = {i: (0 if i < NE else 1) for i in coords}


def X(i):
    return {((i, 1),): 1}


def mono(*gens):
    e = A.const(1)
    for g in gens:
        e = A.mul(e, X(g))
    return e


def rnd():
    return Fr(random.randint(-3, 3))


# random even cubic + quartic interaction
I = {}
ev = [0, 1]
od = [2, 3, 4, 5]
import itertools
for combo in itertools.combinations_with_replacement(ev, 3):
    I = A.add(I, A.scal(mono(*combo), rnd()))
for x in ev:
    for a, b in itertools.combinations(od, 2):
        I = A.add(I, A.scal(mono(x, a, b), rnd()))
for combo in itertools.combinations_with_replacement(ev, 4):
    I = A.add(I, A.scal(mono(*combo), rnd()))
for x, y in itertools.combinations_with_replacement(ev, 2):
    for a, b in itertools.combinations(od, 2):
        I = A.add(I, A.scal(mono(x, y, a, b), rnd()))
for a, b, c, d in itertools.combinations(od, 4):
    I = A.add(I, A.scal(mono(a, b, c, d), rnd()))
# order-hbar quadratic vertex C (counterterm-like)
C = {}
for x, y in itertools.combinations_with_replacement(ev, 2):
    C = A.add(C, A.scal(mono(x, y), rnd()))
for a, b in itertools.combinations(od, 2):
    C = A.add(C, A.scal(mono(a, b), rnd()))

# graded symmetric even propagator
P = {}
for a in ev:
    for b in ev:
        if a <= b:
            v = Fr(random.randint(1, 4))
            P[(a, b)] = v
            P[(b, a)] = v
for a in od:
    for b in od:
        if a < b:
            v = Fr(random.randint(-3, 3))
            P[(a, b)] = v
            P[(b, a)] = -v


def dP(F):
    out = {}
    for (a, b), v in P.items():
        out = A.add(out, A.scal(A.dleft(A.dleft(F, a), b), v / 2))
    return out


def gtrunc(m):
    return dict(m).get(G, 0) <= 3


# e^{J/hbar} with J = I - hbar C, weighted by g:  exp(g (u I - C))
expo = A.add(A.mul(A.const(1), A.scal(A.mul(X(U), I), 1)), C, -1)
expo = A.mul(X(G), expo)
E = A.const(1)
term = A.const(1)
for m in range(1, 4):
    term = A.scal(A.mul(term, expo, gtrunc), Fr(1, m))
    E = A.add(E, term)
# e^{h dP}
Z = dict(E)
term = dict(E)
n = 0
while True:
    n += 1
    term = A.scal(A.mul(X(H), dP(term), gtrunc), Fr(1, n))
    if not term:
        break
    Z = A.add(Z, term)
Y = A.add(Z, A.const(1), -1)
Y2 = A.mul(Y, Y, gtrunc)
Y3 = A.mul(Y2, Y, gtrunc)
L = A.add(A.add(Y, Y2, Fr(-1, 2)), Y3, Fr(1, 3))


def extract(F, hpow, deg):
    """part of hbar*F with hbar power hpow and coordinate degree deg, with g,h,u stripped"""
    out = {}
    for m, c in F.items():
        d = dict(m)
        hp = d.get(H, 0) - d.get(U, 0) + 1
        cd = sum(p for g_, p in m if g_ in coords)
        if hp == hpow and cd == deg:
            mm = tuple((g_, p) for g_, p in m if g_ in coords)
            out[mm] = out.get(mm, 0) + c
    return {m: c for m, c in out.items() if c != 0}


W1_3 = extract(L, 1, 3)
W1_2 = extract(L, 1, 2)
W0_3 = extract(L, 0, 3)
print("tree cubic equals I_3:", W0_3 == {m: c for m, c in I.items() if sum(p for _, p in m) == 3})


def hess(F):
    return {(b, c): A.dleft(A.dleft(F, c), b) for b in coords for c in coords}  # d_b d_c F


def matmul(Ma, Mb, trunc=None):
    out = {}
    for a in coords:
        for c in coords:
            acc = {}
            for b in coords:
                acc = A.add(acc, A.mul(Ma[(a, b)], Mb[(b, c)], trunc))
            out[(a, c)] = acc
    return out


def Mop(Hs):
    M = {}
    for a in coords:
        for c in coords:
            acc = {}
            for b in coords:
                if (a, b) in P:
                    acc = A.add(acc, A.scal(Hs[(b, c)], P[(a, b)] * (-1) ** par[b]))
            M[(a, c)] = acc
    return M


def str_(M):
    acc = {}
    for a in coords:
        acc = A.add(acc, A.scal(M[(a, a)], (-1) ** par[a]))
    return acc


def degpart(F, deg):
    return {m: c for m, c in F.items() if sum(p for g_, p in m if g_ in coords) == deg}


def oneloop(Hs, maxn=3, deg=3):
    M = Mop(Hs)
    Mn = M
    tot = {}
    for n_ in range(1, maxn + 1):
        if n_ > 1:
            Mn = matmul(Mn, M)
        tot = A.add(tot, A.scal(str_(Mn), Fr(1, 2 * n_)))
    return degpart(tot, deg)


I3 = degpart(I, 3)
Hs = hess(I)
W1PI = oneloop(Hs, 3, 3)
# (T1) 1PI part: compare with brute force restricted to ... (brute force has 1PR too)
# tree field chi, solved iteratively: chi^a = sum_b P^{ab} (-1)^{|b|} (d_b I)(x + chi), to cubic order
def degtrunc(F, maxdeg):
    return {m: c for m, c in F.items() if sum(p for g_, p in m if g_ in coords) <= maxdeg}
dI = {b: A.dleft(I, b) for b in coords}
def chi_of(sign_par=True, sgn=1):
    ch = {a: {} for a in coords}
    for it in range(3):
        sub = {a: A.add(X(a), ch[a]) for a in coords}
        new = {}
        for a in coords:
            acc = {}
            for b in coords:
                if (a, b) in P:
                    fac = P[(a, b)] * ((-1) ** par[b] if sign_par else 1) * sgn
                    acc = A.add(acc, A.scal(degtrunc(A.subst(dI[b], sub), 3), fac))
            new[a] = degtrunc(acc, 3)
        ch = new
    return ch
def full_candidate(ch):
    sub = {a: A.add(X(a), ch[a]) for a in coords}
    Hs_sh = {k: degtrunc(A.subst(v, sub), 3) for k, v in Hs.items()}
    out = {}
    for deg in (1, 2, 3):
        out[deg] = A.add(oneloop(Hs_sh, 3, deg), degpart(A.subst(C, sub), deg), -1)
    return out
chi = chi_of(False, 1)   # tree field WITHOUT the parity sign
cand = full_candidate(chi)
print("(T2) W_1 cubic == [1PI - C](x+chi):", len(A.add(W1_3, cand[3], -1)) == 0, " terms in W1_3:", len(W1_3))
print("(T2') W_1 quadratic == [1PI - C](x+chi):", len(A.add(W1_2, cand[2], -1)) == 0)
cand2 = full_candidate(chi_of(True, -1))
print("   (control) chi -> -chi fails:", len(A.add(W1_3, cand2[3], -1)) != 0)
cand3 = full_candidate(chi_of(True, 1))
print("   (control) chi WITH parity sign fails:", len(A.add(W1_3, cand3[3], -1)) != 0)

# (T3) theta trick at a Lambda point
R = [30, 31, 32]
for r in R:
    A.gen(r, 1)
TH = {}
e_pt = {}
# Lambda point: even coords get even elements, odd coords get odd elements of Lambda(R)
e_pt[0] = A.add(A.const(Fr(2)), A.mul(X(30), X(31)))
e_pt[1] = A.add(A.const(Fr(-1)), A.scal(A.mul(X(31), X(32)), 3))
e_pt[2] = A.scal(X(30), 2)
e_pt[3] = A.add(X(31), X(32))
e_pt[4] = A.scal(X(32), -1)
e_pt[5] = A.add(X(30), A.mul(A.mul(X(30), X(31)), X(32)))
ok = True
tb, tc = -2, -1
for b in coords:
    for c in coords:
        A.par[tb] = par[b]
        A.par[tc] = par[c]
        direct = A.subst(Hs[(b, c)], e_pt)
        sub_t = {a: dict(e_pt[a]) for a in coords}
        sub_t[b] = A.add(sub_t[b], {((tb, 1),): 1})
        sub_t[c] = A.add(sub_t[c], {((tc, 1),): 1})
        F = A.subst(I, sub_t)
        # coefficient of theta_c theta_b : canonical order is (tb, tc); theta_c theta_b = (-1)^{|b||c|} theta_b theta_c
        coef = {}
        for m, v in F.items():
            d = dict(m)
            if d.get(tb, 0) == 1 and d.get(tc, 0) == 1:
                mm = tuple((g_, p) for g_, p in m if g_ not in (tb, tc))
                coef[mm] = coef.get(mm, 0) + v
        coef = A.scal(coef, (-1) ** (par[b] * par[c]))
        if A.add(coef, direct, -1):
            ok = False
print("(T3) theta trick for d_b d_c at a Lambda point:", ok)
# first derivative theta trick
ok1 = True
for b in coords:
    A.par[tb] = par[b]
    direct = A.subst(A.dleft(I, b), e_pt)
    sub_t = {a: dict(e_pt[a]) for a in coords}
    sub_t[b] = A.add(sub_t[b], {((tb, 1),): 1})
    F = A.subst(I, sub_t)
    coef = {}
    for m, v in F.items():
        d = dict(m)
        if d.get(tb, 0) == 1:
            mm = tuple((g_, p) for g_, p in m if g_ != tb)
            coef[mm] = coef.get(mm, 0) + v
    if A.add(coef, direct, -1):
        ok1 = False
print("(T3') theta trick for d_b at a Lambda point:", ok1)
