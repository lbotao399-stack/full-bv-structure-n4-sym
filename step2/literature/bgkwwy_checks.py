"""Consistency checks on the one-loop formulas of Appendix C of arXiv:2306.01039v2
(Budzik, Gaiotto, Kulp, Williams, Wu, Yu).  Pure Python, exact rationals, SU(2) with f^{ABC}=eps^{ABC}
(paper eq. (3.20)-(3.22)).  Field components: c_A and b_A, A=1,2,3 (Killing form delta^{AB}, indices
raised/lowered freely, paper below (3.19)).  Derivative orders (d1,d2) w.r.t. z_1, z_2.

Contraction conventions for X_alpha Y^alpha (lower index first):
  K1: X_alpha Y^alpha = X_1 Y_2 - X_2 Y_1   (from footnote 54: psi^alpha chi_alpha = psi_2 chi_1 - psi_1 chi_2 = -psi_alpha chi^alpha)
  K2: X_alpha Y^alpha = X_2 Y_1 - X_1 Y_2   (opposite sign)
Every printed result states the convention.  Output: plain text to stdout.
"""
from fractions import Fraction as Fr
from itertools import product

IDX = (1, 2, 3)


def eps(a, b, c):
    p = (a, b, c)
    if len(set(p)) < 3:
        return 0
    inv = sum(1 for i in range(3) for j in range(i + 1, 3) if p[i] > p[j])
    return -1 if inv % 2 else 1


# ---------- expressions: dict  key=(bosons(sorted tuple), fermions(sorted tuple)) -> Fraction ----------
def mono(bos=(), fer=(), coef=1):
    fer = list(fer)
    if len(set(fer)) < len(fer):
        return {}
    sign = 1
    # bubble sort to canonical order, tracking the sign
    for i in range(len(fer)):
        for j in range(len(fer) - 1 - i):
            if fer[j] > fer[j + 1]:
                fer[j], fer[j + 1] = fer[j + 1], fer[j]
                sign = -sign
    return {(tuple(sorted(bos)), tuple(fer)): Fr(coef) * sign}


def add(*es):
    out = {}
    for e in es:
        for k, v in e.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def scale(e, s):
    return {k: v * Fr(s) for k, v in e.items() if v * Fr(s) != 0}


def mul(e1, e2):
    out = {}
    for (b1, f1), v1 in e1.items():
        for (b2, f2), v2 in e2.items():
            m = mono(b1 + b2, f1 + f2, v1 * v2)
            for k, v in m.items():
                out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v != 0}


def c(A, d1, d2):
    return mono((), (("c", A, d1, d2),))


def b(A, d1, d2):
    return mono((("b", A, d1, d2),), ())


def deriv(e, which):
    """holomorphic derivative d/dz_which (bosonic operator, Leibniz, no signs)"""
    out = {}
    for (bos, fer), v in e.items():
        for i, f in enumerate(bos):
            nb = list(bos)
            nb[i] = (f[0], f[1], f[2] + (which == 1), f[3] + (which == 2))
            for k, w in mono(tuple(nb), fer, v).items():
                out[k] = out.get(k, 0) + w
        for i, f in enumerate(fer):
            nf = list(fer)
            nf[i] = (f[0], f[1], f[2] + (which == 1), f[3] + (which == 2))
            for k, w in mono(bos, tuple(nf), v).items():
                out[k] = out.get(k, 0) + w
    return {k: v for k, v in out.items() if v != 0}


def contr(X, Y, conv):
    """X_alpha Y^alpha for X, Y given as functions (d1,d2)->expression of the undifferentiated slot +1 derivative."""
    s = 1 if conv == "K1" else -1
    return scale(add(mul(X(1, 0), Y(0, 1)), scale(mul(X(0, 1), Y(1, 0)), -1)), s)


def ratio(e1, e2):
    """return r with e1 = r*e2 if it exists, else None"""
    if not e1 and not e2:
        return "both zero"
    keys = set(e1) | set(e2)
    r = None
    for k in keys:
        a, bb = e1.get(k, 0), e2.get(k, 0)
        if bb == 0:
            return None
        q = Fr(a) / Fr(bb)
        if r is None:
            r = q
        elif q != r:
            return None
    return r


def binom(n, k):
    from math import comb
    return comb(n, k)


# ---------- paper formulas ----------
def C3(A, B, conv):
    """(C.3): Q1(b_A b_B) = d_al b_A d^al c^B + d_al b_B d^al c^A - 2 delta_AB d_al b_C d^al c^C"""
    t1 = contr(lambda i, j: b(A, i, j), lambda i, j: c(B, i, j), conv)
    t2 = contr(lambda i, j: b(B, i, j), lambda i, j: c(A, i, j), conv)
    t3 = {}
    if A == B:
        for C in IDX:
            t3 = add(t3, contr(lambda i, j, C=C: b(C, i, j), lambda i, j, C=C: c(C, i, j), conv))
    return add(t1, t2, scale(t3, -2))


def C4(A, B, conv):
    """(C.4): Q1(b_A c^B) = (-1)( d_al c^A d^al c^B - delta^B_A d_al c_D d^al c^D )"""
    t1 = contr(lambda i, j: c(A, i, j), lambda i, j: c(B, i, j), conv)
    t2 = {}
    if A == B:
        for D in IDX:
            t2 = add(t2, contr(lambda i, j, D=D: c(D, i, j), lambda i, j, D=D: c(D, i, j), conv))
    return scale(add(t1, scale(t2, -1)), -1)


def F_(A, B, n1, n2, m1, m2):
    return add(mul(c(A, n1 + 1, n2), c(B, m1, m2 + 1)), scale(mul(c(A, n1, n2 + 1), c(B, m1 + 1, m2)), -1))


def H_(A, B, n1, n2, m1, m2):
    return add(mul(c(A, n1 + 1, n2), b(B, m1, m2 + 1)), scale(mul(c(A, n1, n2 + 1), b(B, m1 + 1, m2)), -1),
               scale(mul(b(A, n1 + 1, n2), c(B, m1, m2 + 1)), -1), mul(b(A, n1, n2 + 1), c(B, m1 + 1, m2)))


def G_(C, D, n1, n2, m1, m2):
    out = {}
    for k in range(m1 + 1):
        for l in range(m2 + 1):
            w = Fr(binom(m1, k) * binom(m2, l), k + l + 1)
            term = add(mul(c(C, m1 - k + 1, m2 - l), c(D, n1 + k, n2 + l + 1)),
                       scale(mul(c(C, m1 - k, m2 - l + 1), c(D, n1 + k + 1, n2 + l)), -1))
            out = add(out, scale(term, w))
    return out


def K_(C, D, n1, n2, m1, m2):
    out = {}
    for k in range(m1 + 1):
        for l in range(m2 + 1):
            w = Fr(binom(m1, k) * binom(m2, l), k + l + 1)
            term = add(mul(c(C, m1 - k + 1, m2 - l), b(D, n1 + k, n2 + l + 1)),
                       scale(mul(c(C, m1 - k, m2 - l + 1), b(D, n1 + k + 1, n2 + l)), -1),
                       scale(mul(b(C, m1 - k + 1, m2 - l), c(D, n1 + k, n2 + l + 1)), -1),
                       mul(b(C, m1 - k, m2 - l + 1), c(D, n1 + k + 1, n2 + l)))
            out = add(out, scale(term, w))
    return out


def C6(A, B, n1, n2, m1, m2, C1, C2):
    """(C.6): Q1(d^{n} b_A d^{m} c_B) = -C2 F(c_A,c_B) + C1 2 eps_ADE eps_BCE/(m1+m2+2) G(c_C,c_D)"""
    out = scale(F_(A, B, n1, n2, m1, m2), -C2)
    for C, D, E in product(IDX, IDX, IDX):
        w = eps(A, D, E) * eps(B, C, E)
        if w:
            out = add(out, scale(G_(C, D, n1, n2, m1, m2), Fr(2 * w, m1 + m2 + 2) * C1))
    return out


def C7(A, B, n1, n2, m1, m2, C1, C2):
    """(C.7): Q1(d^{n} b_A d^{m} b_B) = -C2 H(b_A,b_B,c_A,c_B) + C1 2 eps_ADE eps_BCE/(m1+m2+2) K(b_C,b_D,c_C,c_D)"""
    out = scale(H_(A, B, n1, n2, m1, m2), -C2)
    for C, D, E in product(IDX, IDX, IDX):
        w = eps(A, D, E) * eps(B, C, E)
        if w:
            out = add(out, scale(K_(C, D, n1, n2, m1, m2), Fr(2 * w, m1 + m2 + 2) * C1))
    return out


def show(e):
    if not e:
        return "0"
    parts = []
    for (bos, fer), v in sorted(e.items()):
        name = " ".join(f"{f[0]}{f[1]}[{f[2]},{f[3]}]" for f in bos + fer)
        parts.append(f"{v}*{name}")
    return " + ".join(parts)


if __name__ == "__main__":
    print("== 1. tensor structure of (C.3), (C.4) (both contraction conventions) ==")
    for conv in ("K1", "K2"):
        ok3 = ok4 = True
        for A, B in product(IDX, IDX):
            # claimed: (C.3) = -(eps_ACE eps_BDE + eps_BCE eps_ADE) d b^C . d c^D
            rhs3 = {}
            rhs4 = {}
            for C, D, E in product(IDX, IDX, IDX):
                w = eps(A, C, E) * eps(B, D, E) + eps(B, C, E) * eps(A, D, E)
                if w:
                    rhs3 = add(rhs3, scale(contr(lambda i, j, C=C: b(C, i, j), lambda i, j, D=D: c(D, i, j), conv), -w))
                w4 = eps(A, C, E) * eps(B, D, E)
                if w4:
                    rhs4 = add(rhs4, scale(contr(lambda i, j, C=C: c(C, i, j), lambda i, j, D=D: c(D, i, j), conv), w4))
            ok3 &= ratio(C3(A, B, conv), rhs3) in (1, "both zero")
            ok4 &= ratio(C4(A, B, conv), rhs4) in (1, "both zero")
        print(f"[{conv}] (C.3) == -(f_ACE f_BDE + f_BCE f_ADE) d_al b^C d^al c^D : {ok3}")
        print(f"[{conv}] (C.4) == +f_ACE f_BDE d_al c^C d^al c^D             : {ok4}")

    print("\n== 2. (C.5) ansatz coefficient C1 read off from (C.4) at A=B=1 ==")
    for conv in ("K1", "K2"):
        e = C4(1, 1, conv)
        print(f"[{conv}] Q1(b_1 c_1) from (C.4) = {show(e)}")

    print("\n== 3. (C.6)/(C.7) at n=m=0 with (C1,C2)=(-2,0) versus (C.4)/(C.3) ==")
    for conv in ("K1", "K2"):
        r6 = {ratio(C6(A, B, 0, 0, 0, 0, -2, 0), C4(A, B, conv)) for A, B in product(IDX, IDX)}
        r7 = {ratio(C7(A, B, 0, 0, 0, 0, -2, 0), C3(A, B, conv)) for A, B in product(IDX, IDX)}
        print(f"[{conv}] ratio (C.6)/(C.4) over all A,B: {r6};  ratio (C.7)/(C.3): {r7}")
    print("(C.6) at n=m=0, A=B=1, C1 part (C1=1,C2=0):", show(C6(1, 1, 0, 0, 0, 0, 1, 0)))
    print("(C.6) at n=m=0, A=B=1, C2 part (C1=0,C2=1):", show(C6(1, 1, 0, 0, 0, 0, 0, 1)))
    # is the C2-structure of (C.6) a contribution compatible with (C.4)?  (C.4) has no c_1 c_1 term:
    print("(C.4)[K1] has c1c1 component:", any(k[1][0][1] == 1 and k[1][1][1] == 1 for k in C4(1, 1, "K1")))

    print("\n== 4. Q1 Tr b^2 (Tr t^A t^B = -1/2 delta^AB, (3.20); b^2 = b_A b^A = -2 Tr b^2, (3.23)) ==")
    for conv in ("K1", "K2"):
        q = {}
        for A in IDX:
            q = add(q, C3(A, A, conv))
        q = scale(q, Fr(-1, 2))  # Q1 Tr b^2 = -1/2 Q1(b_A b_A)
        # d_al Tr(b d^al c) = Tr(d_al b d^al c) = -1/2 sum_C d_al b_C d^al c_C  (d_al d^al c = 0)
        dS = {}
        for C in IDX:
            dS = add(dS, contr(lambda i, j, C=C: b(C, i, j), lambda i, j, C=C: c(C, i, j), conv))
        dS = scale(dS, Fr(-1, 2))
        print(f"[{conv}] from (C.3): Q1 Tr b^2 = r * d_al Tr(b d^al c) with r = {ratio(q, dS)}")
        q7 = {}
        for A in IDX:
            q7 = add(q7, C7(A, A, 0, 0, 0, 0, -2, 0))
        q7 = scale(q7, Fr(-1, 2))
        print(f"[{conv}] from (C.7), C1=-2, C2=0: Q1 Tr b^2 = r * d_al Tr(b d^al c) with r = {ratio(q7, dS)}")

    print("\n== 5. Q1 on S_al = Tr b d_al c = -1/2 b_A d_al c_A via (C.6), C1=-2, C2=0 ==")
    for conv in ("K1", "K2"):
        # target: d_al Tr(d_be c d^be c) = -1/2 d_al sum_C (d_be c_C d^be c_C)
        cc = {}
        for C in IDX:
            cc = add(cc, contr(lambda i, j, C=C: c(C, i, j), lambda i, j, C=C: c(C, i, j), conv))
        trcc = scale(cc, Fr(-1, 2))
        for al, (m1, m2) in ((1, (1, 0)), (2, (0, 1))):
            q = {}
            for A in IDX:
                q = add(q, C6(A, A, 0, 0, m1, m2, -2, 0))
            q = scale(q, Fr(-1, 2))
            print(f"[{conv}] Q1 S_{al} = r * d_{al} Tr(d_be c d^be c), r = {ratio(q, deriv(trcc, al))}")
        # same with the C2 structure alone
        for al, (m1, m2) in ((1, (1, 0)), (2, (0, 1))):
            q = {}
            for A in IDX:
                q = add(q, C6(A, A, 0, 0, m1, m2, 0, 1))
            q = scale(q, Fr(-1, 2))
            print(f"[{conv}] C2-structure alone (C1=0,C2=1): Q1 S_{al} = r * d_{al} Tr(d_be c d^be c), r = {ratio(q, deriv(trcc, al))}")

    print("\n== 6. (C.6) normalised to reproduce (C.4) at n=m=0: C1=+1 (K1) or C1=-1 (K2), C2=0 ==")
    for conv, C1n in (("K1", 1), ("K2", -1)):
        ok = all(ratio(C6(A, B, 0, 0, 0, 0, C1n, 0), C4(A, B, conv)) in (1, "both zero") for A, B in product(IDX, IDX))
        ok7 = all(ratio(C7(A, B, 0, 0, 0, 0, C1n, 0), C3(A, B, conv)) in (1, "both zero") for A, B in product(IDX, IDX))
        cc = {}
        for C in IDX:
            cc = add(cc, contr(lambda i, j, C=C: c(C, i, j), lambda i, j, C=C: c(C, i, j), conv))
        trcc = scale(cc, Fr(-1, 2))
        rs = []
        for al, (m1, m2) in ((1, (1, 0)), (2, (0, 1))):
            q = {}
            for A in IDX:
                q = add(q, C6(A, A, 0, 0, m1, m2, C1n, 0))
            rs.append(ratio(scale(q, Fr(-1, 2)), deriv(trcc, al)))
        print(f"[{conv}] C1={C1n}: (C.6)=(C.4) at n=m=0: {ok}; (C.7)=(C.3) at n=m=0: {ok7}; "
              f"Q1 S_al = r d_al Tr(d_be c d^be c), r(al=1,2) = {rs}")
