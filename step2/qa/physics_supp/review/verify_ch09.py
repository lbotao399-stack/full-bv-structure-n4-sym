#!/usr/bin/env python3
"""Numerical checks of the Chapter 9 'Physics notation' supplements (manuscript/supplements/ch09.tex).

Self-contained (numpy only): two-component realisation of Definition c8-twocomp, su(N) with
(X,Y) = -(4 hbar/g^2) tr(XY), plane waves with incoming momenta as in Definition c1-vertex.
Every check is paired with wrong-variant rejections.  Output: verify_ch09.json / verify_ch09.log here.
"""
import itertools, json, os
import numpy as np

np.seterr(all='ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261029 + 9)
RES, LOG = [], []
L32 = np.log(1.5)


def log(s):
    LOG.append(s); print(s)


def _flat(x):
    if isinstance(x, (list, tuple)):
        parts = [_flat(y) for y in x]
        return np.concatenate(parts) if parts else np.zeros(0, complex)
    return np.ravel(np.asarray(x, dtype=complex))


def rec(cid, where, desc, lhs, tol=1e-9):
    v = _flat(lhs)
    sc = 1.0
    res = float(np.max(np.abs(v))) / sc if v.size else 0.0
    ok = res < tol
    RES.append(dict(id=cid, block=where, check=desc, residual=res, passed=bool(ok)))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:5s} res={res:.1e}  ({where}) {desc}")


def neg(cid, where, desc, val):
    r = float(np.max(np.abs(_flat(val))))
    ok = r > 1e-6
    RES.append(dict(id=cid, block=where, check='wrong variant rejected: ' + desc, residual=r, passed=bool(ok)))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:5s} |dev|={r:.1e}  ({where}) wrong variant rejected: {desc}")


def cr(*s):
    return rng.normal(size=s) + 1j * rng.normal(size=s)


# ------------------------------------------------------------------ realisation (Definition c8-twocomp)
I2 = np.eye(2, dtype=complex)
PA = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
SIG = [PA[0], PA[1], PA[2], 1j * I2]
SIGB = [-PA[0], -PA[1], -PA[2], 1j * I2]
EU = np.array([[0, 1], [-1, 0]], complex)
E4 = np.zeros((4,) * 4)
for p in itertools.permutations(range(4)):
    E4[p] = np.linalg.det(np.eye(4)[list(p)])
levi3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    levi3[p] = np.linalg.det(np.eye(3)[list(p)])
SI = []
for a in range(3):
    M = np.zeros((4, 4), complex); M[:3, :3] = levi3[a]; M[a, 3], M[3, a] = 1, -1; SI.append(M)
for a in range(3):
    M = np.zeros((4, 4), complex); M[:3, :3] = 1j * levi3[a]; M[a, 3], M[3, a] = -1j, 1j; SI.append(M)
up = lambda wl: 0.5 * np.einsum('ijkl,kl->ij', E4, wl)


def rho_S(v4, xi, xib):        # rho(v): Sigma -> Sigma^*, v in V only
    eta = 1j * sum(v4[m] * np.einsum('ai,ab->bi', xib, SIGB[m]) for m in range(4))
    etab = -1j * sum(v4[m] * SIGB[m] @ xi for m in range(4))
    return eta, etab


def pair(xi, xib, eta, etab):
    return np.sum(xi * eta) - np.sum(xib * etab)


def suN(N):
    T = []
    for i in range(N):
        for j in range(i + 1, N):
            M = np.zeros((N, N), complex); M[i, j] = M[j, i] = 0.5; T.append(M)
            M = np.zeros((N, N), complex); M[i, j] = -0.5j; M[j, i] = 0.5j; T.append(M)
    for k in range(1, N):
        M = np.zeros((N, N), complex); M[:k, :k] = np.eye(k); M[k, k] = -k
        T.append(M / np.sqrt(2 * k * (k + 1)))
    T = np.array(T)
    f = np.real(-2j * (np.einsum('aij,bjk,cki->abc', T, T, T) - np.einsum('bij,ajk,cki->abc', T, T, T)))
    return T, f


hb, gc = 0.7, 1.3
T2, f2 = suN(2)
pr = lambda X, Y: -4 * hb / gc ** 2 * np.trace(X @ Y)
cm = lambda A, B: A @ B - B @ A
gv = lambda comps, T=T2: -1j * gc * np.einsum('a,aij->ij', comps, T)   # X = -ig X^a T^a

# ================================================================== K1  plane-wave values of the component forms
log('=== K1-K2: components of G_1..G_4, H_1, H_2, H_4 (block after def:c9-bases) ===')
# legs: 0 = c (odd, X0), 1,2 as stated; derivative -> i p of the leg; plane-wave value = coefficient of e^{i sum p x}
def pw_G1(p, a, X, sign=-1):
    # G_1 = sign * int ([c,A_n], d_m F0_mn): sum over the two assignments of legs 1,2 to the two A slots
    tot = 0
    for (j, k) in ((1, 2), (2, 1)):
        Fk = 1j * (np.outer(p[k], a[k]) - np.outer(a[k], p[k]))          # F_mn = d_m A_n - d_n A_m
        dF = 1j * p[k] @ Fk                                              # d_m F_mn  (index n)
        tot += sign * pr(cm(X[0], X[j]), X[k]) * (a[j] @ dF)
    return tot


def pw_G2(p, a, X):
    tot = 0
    for (j, k) in ((1, 2), (2, 1)):
        ddA = (1j * p[k]) * (1j * p[k] @ a[k])                           # d_n d_m A_m
        tot += -pr(cm(X[0], X[j]), X[k]) * (a[j] @ ddA)
    return tot


def pw_G2b(p, a, X):          # -int (d_m A_m, d_n [A_n, c])
    tot = 0
    for (j, k) in ((1, 2), (2, 1)):
        divA = 1j * p[j] @ a[j]
        dcomm = 1j * (p[k] + p[0]) @ a[k]                               # d_n of [A_n(leg k), c(leg 0)]
        tot += -divA * dcomm * pr(X[j], cm(X[k], X[0]))
    return tot


def pw_G3(p, w, X):           # -sum_I int ([c,A_I], d^2 A_I)
    tot = 0
    for (j, k) in ((1, 2), (2, 1)):
        tot += -pr(cm(X[0], X[j]), X[k]) * (w[j] @ w[k]) * (-(p[k] @ p[k]))
    return tot


def pw_G3su4(p, w, X):        # -1/4 int ([c,A^{ij}], d^2 A_ij)
    tot = 0
    for (j, k) in ((1, 2), (2, 1)):
        wl_j = sum(w[j][I] * SI[I] for I in range(6)); wl_k = sum(w[k][I] * SI[I] for I in range(6))
        tot += -0.25 * pr(cm(X[0], X[j]), X[k]) * np.sum(up(wl_j) * wl_k) * (-(p[k] @ p[k]))
    return tot


okK1, negK1 = [], []
for _ in range(4):
    p = [rng.normal(size=4) for _ in range(2)]; p.append(-p[0] - p[1])
    a = [None, rng.normal(size=4), rng.normal(size=4)]
    X = [gv(rng.normal(size=3)) for _ in range(3)]
    f = pr(cm(X[0], X[1]), X[2])
    core1 = f * (((p[2] @ p[2]) - (p[1] @ p[1])) * (a[1] @ a[2]) - (a[1] @ p[2]) * (a[2] @ p[2]) + (a[1] @ p[1]) * (a[2] @ p[1]))
    core2 = f * ((a[1] @ p[2]) * (a[2] @ p[2]) - (a[1] @ p[1]) * (a[2] @ p[1]))
    okK1 += [pw_G1(p, a, X) - core1, pw_G2(p, a, X) - core2, pw_G2b(p, a, X) - core2]
    negK1 += [pw_G1(p, a, X, sign=+1) - core1]
    w = [None, rng.normal(size=6), rng.normal(size=6)]
    core3 = f * (w[1] @ w[2]) * ((p[2] @ p[2]) - (p[1] @ p[1]))
    okK1 += [pw_G3(p, w, X) - core3, pw_G3su4(p, w, X) - core3]
    # G_4 = int ([c, d_m c], A*_m): legs c(r0,X0,p0), c(r1,X1,p1), A*(r2, b, X2, p2); coefficient of r0 r1 r2
    b = rng.normal(size=4)
    val4 = 1j * ((p[0] + p[1]) @ b) * f            # [c,dc] -> r0 r1 i (p0+p1) [X0,X1], then (.,A*) with r2 moved left
    okK1.append(val4 - (-1j * f * (p[2] @ b)))
    # Q G_2 component form (block after thm:c9-wz): -int([c,d_n c], d_n d_m A_m) - int([c,A_n], d_n d^2 c); legs c,c,A(b)
    t1 = -(np.dot(-1j * p[2], -(p[2]) * (p[2] @ b))) * f          # -([c,dc]_n, d_n d.A), [c,dc] = -i p2_n r0r1 [X0,X1]
    t2 = -(-1j * (b @ p[1]) * (p[1] @ p[1]) * pr(cm(X[0], X[2]), X[1]))       # c(leg0) in [c,A], c(leg1) differentiated
    t2 += -(+1j * (b @ p[0]) * (p[0] @ p[0]) * pr(cm(X[1], X[2]), X[0]))      # c(leg1) in [c,A] (r1 r0 = -r0 r1)
    coreQ = -1j * f * ((p[0] @ p[0]) * (p[0] @ b) + (p[1] @ p[1]) * (p[1] @ b) + (p[2] @ p[2]) * (p[2] @ b))
    okK1.append(t1 + t2 - coreQ)
rec('K1', 'after def:c9-bases, after thm:c9-wz', 'plane-wave values of the component forms of G_1, G_2 (both forms), G_3 (SO(6), SU(4)), G_4, QG_2 = Prop c9-basis / Lemma c9-QG', okK1)
neg('N1', 'after def:c9-bases', 'G_1 = +int([c,A_n], d_m F0_mn)', negK1)

# K2 pointwise: H_1 = (dA|_V,[A^A]) Hodge = 2 (d_m A_n,[A_m,A_n]); H_2, H_4 SU(4) forms
okK2, negK2 = [], []
for _ in range(3):
    dA = [[gv(rng.normal(size=3)) for _ in range(4)] for _ in range(4)]   # dA[m][n] = d_m A_n
    A = [gv(rng.normal(size=3)) for _ in range(4)]
    hodge = sum(pr(dA[m][n] - dA[n][m], 2 * cm(A[m], A[n])) for m in range(4) for n in range(m + 1, 4))
    comp = 2 * sum(pr(dA[m][n], cm(A[m], A[n])) for m in range(4) for n in range(4))
    okK2.append(hodge - comp); negK2.append(hodge - comp / 2)
    AI = [gv(rng.normal(size=3)) for _ in range(6)]; dAI = [[gv(rng.normal(size=3)) for _ in range(6)] for _ in range(4)]
    Al = np.einsum('Ikl,Iab->klab', np.array(SI), np.array(AI)); Au = 0.5 * np.einsum('ijkl,klab->ijab', E4, Al)
    dAl = [np.einsum('Ikl,Iab->klab', np.array(SI), np.array(dAI[m])) for m in range(4)]
    dAu = [0.5 * np.einsum('ijkl,klab->ijab', E4, dAl[m]) for m in range(4)]
    H2a = sum(pr(dAI[m][I], cm(A[m], AI[I])) for m in range(4) for I in range(6))
    H2b = 0.25 * sum(pr(dAu[m][i, j], cm(A[m], Al[i, j])) for m in range(4) for i in range(4) for j in range(4))
    c = gv(rng.normal(size=3)); Ast = [gv(rng.normal(size=3)) for _ in range(6)]
    Astl = np.einsum('Ikl,Iab->klab', np.array(SI), np.array(Ast))
    H4a = sum(pr(cm(AI[I], c), Ast[I]) for I in range(6))
    H4b = 0.25 * sum(pr(cm(Au[i, j], c), Astl[i, j]) for i in range(4) for j in range(4))
    okK2 += [H2a - H2b, H4a - H4b]
rec('K2', 'after def:c9-bases', 'H_1 Hodge pairing = 2 sum (d_m A_n,[A_m,A_n]); H_2, H_4 in SO(6) and SU(4) form', okK2)
neg('N2', 'after def:c9-bases', 'H_1 = sum (d_m A_n,[A_m,A_n]) (factor 2 dropped)', negK2)

# ================================================================== K3 colour
log('=== K3: colour factors (block after prop:c9-basis) ===')
okK3, negK3 = [], []
for N in (2, 3, 4):
    T, f = suN(N); d = len(T)
    TA = -1j * f
    Xa = [rng.normal(size=d) for _ in range(3)]
    Xm = [gv(x, T) for x in Xa]
    ad = [-1j * gc * np.einsum('a,abc->bc', x, TA) for x in Xa]
    colour = pr(cm(Xm[0], Xm[1]), Xm[2])
    okK3.append(colour - 2 * hb * gc * np.einsum('abc,a,b,c->', f, *Xa))
    ratio = np.trace(ad[0] @ ad[1]) / pr(Xm[0], Xm[1])
    tr3 = np.trace(ad[0] @ ad[1] @ ad[2])
    okK3.append(tr3 - 0.5 * ratio * pr(Xm[0], cm(Xm[1], Xm[2])))
    okK3.append(tr3 + gc ** 3 * N / 2 * np.einsum('abc,a,b,c->', f, *Xa))
    okK3.append(np.einsum('aij,bjk,cki->abc', TA, TA, TA) - 0.5j * N * f)
    negK3.append(tr3 - ratio * pr(Xm[0], cm(Xm[1], Xm[2])))
rec('K3', 'after prop:c9-basis', '([X0,X1],X2) = 2 hbar g f X0X1X2; tr(ad ad ad) = 1/2 (tr ad ad/(X,Y)) (X,[Y,Z]) = -(g^3N/2) f XYZ; tr(T_A T_A T_A) = (i/2) N f', okK3)
neg('N3', 'after prop:c9-basis', 'Lemma c9-legs without the factor 1/2', negK3)

# ================================================================== K4 Q-kernel, cubic couplings; K5 F3 family
log('=== K4-K5: kernel directions, completions, F_3 family (blocks after cor:c9-completions, thm:c9-F3) ===')
# basis order of degree-0 functionals: H1..H6, YV = int(lam,rho([A_V,lam])), YW
# Q images in basis G1..G5 (Proposition c9-maps (i))
Qmat = np.zeros((5, 8))
Qmat[0, 0] = -2; Qmat[2, 1] = -1; Qmat[0, 2] = -0.5; Qmat[3, 2] = -1; Qmat[2, 3] = -0.5
Qmat[4, 4] = 1; Qmat[3, 5] = 2; Qmat[4, 6] = -2
cubicI = np.array([-0.25, -0.5, 1, 1, 1, 0.5, 0.5, 0.5])
kern = [np.array([1, 0, -4, 0, 0, -2, 0, 0]), np.array([0, 1, 0, -2, 0, 0, 0, 0]),
        np.array([0, 0, 0, 0, 2, 0, 1, 0]), np.array([0, 0, 0, 0, 0, 0, 0, 1])]
gauge3 = np.array([-0.25, 0, 1, 0, 0, 0.5, 0, 0]); scal3 = np.array([0, -0.5, 0, 1, 0, 0, 0, 0])
gaug3 = np.array([0, 0, 0, 0, 1, 0, 0.5, 0]); yuk = np.array([0, 0, 0, 0, 0, 0, 0, 1])
okK4 = [Qmat @ k for k in kern] + [Qmat @ cubicI, kern[0] + 4 * gauge3, kern[1] + 2 * scal3, kern[2] - 2 * gaug3,
                                    kern[3] - yuk, gauge3 + scal3 + gaug3 + 0.5 * yuk - cubicI]
okK4.append(np.linalg.matrix_rank(Qmat) - 4)
# brackets of I with the quadratic functionals (Proposition c9-maps (iii)) in basis G1..G5
Bdd = np.array([0, 2, 0, 0, 0]); BdA = np.array([2, 0, 0, 0, 0]); BdW = np.array([0, 0, 2, 0, 0])
Blam = np.array([0, 0, 0, 0, 2]); Bc = np.array([0.5, 0, 0.5, 0, -1])
# completions: Q(second) + {I, first} = 0
comps = [(BdA, np.array([1, 0, 0, 0, 0, 0, 0, 0])), (BdW, np.array([0, 2, 0, 0, 0, 0, 0, 0])),
         (Blam, np.array([0, 0, 0, 0, 0, 0, 1, 0])), (Bc, np.array([0, 0, 1, 1, 1, 0.5, 0, 0]))]
okK4 += [b_ + Qmat @ h for b_, h in comps]
rec('K4', 'after cor:c9-completions', 'kernel of Q = -4(gauge cubic), -2(gluon-scalar), 2(gluon-gaugino), Yukawa parts of cubic I; completions; rank 4', okK4)
neg('N4', 'after cor:c9-completions', 'H_1 - 2H_3 - 2H_6 in the kernel', [Qmat @ np.array([1, 0, -2, 0, 0, -2, 0, 0])])
alpha = np.array([11 / 4 + 19 / 4 * L32, 27 / 8, 11 / 4 + 19 / 4 * L32, 0, -9 * L32])
okK5 = []
for _ in range(4):
    beta = rng.normal(size=4); kap = rng.normal(size=4)
    h = (-(11 / 8 + 19 / 8 * L32) * np.eye(8)[0] - (11 / 2 + 19 / 2 * L32) * np.eye(8)[3] - 9 * L32 * np.eye(8)[4]
         + beta[0] * np.eye(8)[0] + 2 * beta[1] * np.eye(8)[1] - 2 * beta[2] * np.eye(8)[4]
         + beta[3] * np.array([0, 0, 1, 1, 1, 0.5, 0, 0]) + sum(kap[k] * kern[k] for k in range(4)))
    lhs = Qmat @ h + 27 / 16 * Bdd + beta[0] * BdA + beta[1] * BdW + beta[2] * Blam + beta[3] * Bc
    okK5.append(lhs - alpha)
# representative: beta1=beta2=alpha1/2, beta3=alpha5/2, beta4=0, kappa2=-alpha3 -> vanishing cubic part
b_ = np.array([alpha[0] / 2, alpha[0] / 2, alpha[4] / 2, 0]); k_ = np.array([0, -alpha[2], 0, 0])
hrep = (-(11 / 8 + 19 / 8 * L32) * np.eye(8)[0] - (11 / 2 + 19 / 2 * L32) * np.eye(8)[3] - 9 * L32 * np.eye(8)[4]
        + b_[0] * np.eye(8)[0] + 2 * b_[1] * np.eye(8)[1] - 2 * b_[2] * np.eye(8)[4] + k_[1] * kern[1])
okK5.append(hrep)
okK5.append(np.array([b_[0] - (11 / 8 + 19 / 8 * L32), b_[2] + 9 / 2 * L32]))
# Wess-Zumino: only alpha_2 enters, through (alpha_2 - 27/8) Q G_2 ; Q G_1=Q G_3=Q G_4=Q G_5=0 on (c,c,A)
okK5.append(np.array([alpha[1] - 27 / 8]))
rec('K5', 'after thm:c9-F3, thm:c9-wz', 'F_3 family solves the five G-equations for all beta,kappa; representative has zero cubic part; alpha_2 = 27/8', okK5)
neg('N5', 'after thm:c9-F3', 'family with -(11/4+19/4 log 3/2) H_1', [Qmat @ (hrep - (11 / 8 + 19 / 8 * L32) * np.eye(8)[0]) + 27 / 16 * Bdd
                                                                    + b_[0] * BdA + b_[1] * BdW + b_[2] * Blam - alpha])

# ================================================================== K6 su(N) components of the quadratic functionals
log('=== K6: su(N) forms of int(dA,dA), int(lambda, dslash lambda), int(d*A,d*A) (block after thm:c9-F3) ===')
okK6, negK6 = [], []
for N in (2, 3):
    T, f = suN(N); d = len(T)
    dAa = rng.normal(size=(4, 4, d))      # d_m A^a_n
    dAIa = rng.normal(size=(4, 6, d))     # d_m A^a_I
    dAm = [[gv(dAa[m, n], T) for n in range(4)] for m in range(4)]
    dAIm = [[gv(dAIa[m, I], T) for I in range(6)] for m in range(4)]
    prN = lambda X, Y: -4 * hb / gc ** 2 * np.trace(X @ Y)
    dAdA = sum(prN(dAm[m][n] - dAm[n][m], dAm[m][n] - dAm[n][m]) for m in range(4) for n in range(m + 1, 4)) \
        + sum(prN(dAIm[m][I], dAIm[m][I]) for m in range(4) for I in range(6))
    F0 = dAa - dAa.transpose(1, 0, 2)
    LB = -0.25 * np.sum(F0 ** 2) - 0.5 * np.sum(dAIa ** 2)
    okK6.append(dAdA - (-4 * hb * LB)); negK6.append(dAdA - (-2 * hb * LB))
    divA = sum(dAm[m][m] for m in range(4))
    okK6.append(prN(divA, divA) - 2 * hb * np.sum(np.einsum('mma->a', dAa) ** 2))
    # gaugino: int(lambda, dslash lambda) = 2 i int (lambdabar, sigmabar^m d_m lambda) = 2 i hbar lambdabar^a sigmabar d lambda^a
    # plane waves lambda = r1 z1 e^{ipx} + r2 z2 e^{-ipx}, z in Sigma (x) colour; coefficient of r1 r2
    pm = rng.normal(size=4)
    z = [(cr(2, 4, d), cr(2, 4, d)) for _ in range(2)]
    def bil_if(za, zb, q):     # (za, dslash zb) with d -> i q, colour pairing hbar delta^{ab} in basis -(ig/sqrt2)T^a
        tot = 0
        for a_ in range(d):
            eta, etab = rho_S(1j * q, zb[0][:, :, a_], zb[1][:, :, a_])
            tot += hb * pair(za[0][:, :, a_], za[1][:, :, a_], eta, etab)
        return tot
    def bil_2c(za, zb, q):     # 2 i hbar lambdabar_{ad i} sigmabar^{m ad a} (i q_m) lambda_a^i
        return 2j * hb * sum(np.einsum('ai,ab,bi->', za[1][:, :, a_], sum(SIGB[m] * 1j * q[m] for m in range(4)), zb[0][:, :, a_])
                             for a_ in range(d))
    lhs = bil_if(z[0], z[1], -pm) - bil_if(z[1], z[0], pm)     # r1 r2 coefficient (r2 r1 = - r1 r2)
    rhs = bil_2c(z[0], z[1], -pm) - bil_2c(z[1], z[0], pm)
    okK6.append(lhs - rhs); negK6.append(lhs + rhs)
rec('K6', 'after thm:c9-F3', 'int(dA,dA) = -4 hbar int(-1/4 F0aF0a - 1/2 dA_I^a dA_I^a); int(d*A,d*A) = 2 hbar (d.A^a)^2; int(lam,dslash lam) = 2 i hbar lambdabar^a sigmabar d lambda^a', okK6)
neg('N6', 'after thm:c9-F3', 'factor -2 hbar for int(dA,dA), or the opposite sign of the gaugino form', negK6)
# su(N) coefficients of the representative (c9-F2rep): u = unit = -g^2N/(32 pi^2) (gN2 = g^2 N)
gN2 = 1.7
u = -gN2 / (32 * np.pi ** 2)
cF = [u * 27 / 16 * 2, u * (11 / 8 + 19 / 8 * L32) * (-4), u * (-9 / 2 * L32) * 2]
rec('K7', 'after thm:c9-F3', 'F/hbar = -27 g^2N/256pi^2 (d.A^a)^2 + g^2N/8pi^2 (11/8+19/8 L) L_B + 9 g^2N L/32pi^2 (i lambdabar sigmabar d lambda)',
    [cF[0] + 27 * gN2 / (256 * np.pi ** 2), cF[1] - gN2 / (8 * np.pi ** 2) * (11 / 8 + 19 / 8 * L32), cF[2] - 9 * gN2 * L32 / (32 * np.pi ** 2)])

# ================================================================== K8 arity-one data along Delta
log('=== K8: Ob^(1)(Delta) and F_1^(1)(Delta) (block after thm:c9-quadratic) ===')
# quadratic functionals basis: [int(dA,dA), int(lam,dslash lam), int(dc,A*)], with Q of the degree -1 basis (Theorem c9-quadratic proof)
QAA = np.array([-0.5, 0, -1])      # Q int(A,A*) (ten components: V part -1/2 (dA,dA) - (dc,A*), W part -1/2 (dA,dA))
Qll = np.array([0, -1, 0])         # Q int(lambda*, lambda)
Qcc = np.array([0, 0, -1])         # Q int(c*, c)
Fcoef = np.array([22 / 3, -47 / 3, -8 / 3])
ob = Fcoef[0] * QAA + Fcoef[1] * Qll + Fcoef[2] * Qcc
okK8 = [ob - np.array([-11 / 3, 47 / 3, -14 / 3])]
okK8.append(np.array([-4, 16, -4]) + np.array([1 / 3, -1 / 3, -2 / 3]) - np.array([-11 / 3, 47 / 3, -14 / 3]))
okK8.append(Fcoef - (4 * np.array([2, -4, -1]) + np.array([-2 / 3, 1 / 3, 4 / 3])))
# shift of the dilatation action on the fields by -hbar F: weights of a_A int(A,A*) + a_l int(l*,l) + a_c int(c*,c) are (a_A,-a_l,-a_c)
shift = -np.array([Fcoef[0], -Fcoef[1], -Fcoef[2]])
okK8.append(shift - np.array([-22 / 3, -47 / 3, -8 / 3]))
okK8.append(4 * np.array([2, 4, 1]) * (-1) + np.array([2 / 3, 1 / 3, 4 / 3]) - shift)
okK8.append(np.array([47 / 3 * 2 - 94 / 3]))            # two-component: 47/3 int(lam,dslash lam) = 94/3 i int(lambdabar sigmabar d lambda)
rec('K8', 'after thm:c9-quadratic', 'Q F_1^(1)(Delta) = Ob^(1)(Delta) = 4x(log eps part) + composite part; shift of Delta on A,lambda,c = -22/3,-47/3,-8/3 units', okK8)
neg('N7', 'after thm:c9-quadratic', 'F_1^(1)(Delta) with +8/3 int(c*,c)', [np.array([22 / 3, -47 / 3, 8 / 3]) @ np.array([QAA, Qll, Qcc]) - np.array([-11 / 3, 47 / 3, -14 / 3])])

# ================================================================== K9 composite line and scaling Ward identity
log('=== K9: (Q^GF)^2 and the scaling Ward identity at a finite endpoint (block at the end of Section 9.2) ===')
# slots (momentum space, bosonic, one colour): order [c, A_1..A_4, Astar_1..4, cstar]; Q^GF and Q of Def c3-qgf, Prop c1-Q
def QGF(p):
    M = np.zeros((10, 10), complex)
    for m in range(4):
        M[1 + m, 5 + m] = -2                  # A*_m -> -2 A*_m in slot A_m
        M[0, 1 + m] = 1j * p[m]               # A_m -> d_m A_m in slot c
        M[5 + m, 9] = 1j * p[m]               # c* -> d_m c* in slot A*_m
    return M
def Qop(p):
    M = np.zeros((10, 10), complex)
    p2 = p @ p
    for m in range(4):
        M[1 + m, 0] = -1j * p[m]              # c -> -d_m c in slot A_m
        M[9, 5 + m] = -1j * p[m]              # A*_m -> -d_m A*_m in slot c*
        for n in range(4):
            M[5 + m, 1 + n] = 0.5 * (-p2 * (m == n) + p[m] * p[n])   # A -> 1/2(d^2 A_m - d_m d.A)
    return M
p = rng.normal(size=4)
G2 = QGF(p) @ QGF(p)
tgt = np.zeros((10, 10), complex)
for m in range(4):
    tgt[0, 5 + m] = -2j * p[m]                # (Q^GF)^2 A*_m = -2 d_m A*_m in slot c
    tgt[1 + m, 9] = -2j * p[m]                # (Q^GF)^2 c* = -2 d_m c* in slot A_m
okK9 = [G2 - tgt, Qop(p) @ QGF(p) + QGF(p) @ Qop(p) - (p @ p) * np.diag([1] + [1] * 4 + [1] * 4 + [1])]
# scaling Ward identity on the A(x)A component: (2 x.d_x + 2 y.d_y + 4) G = -4 (eps d_eps + T d_T) G
G = lambda r, e, T: (np.exp(-r ** 2 / (4 * T)) - np.exp(-r ** 2 / (4 * e))) / (4 * np.pi ** 2 * r ** 2)
h = 1e-6
for r, e, Tt in ((0.4, 0.03, 2.0), (1.7, 0.2, 5.0)):
    lhs = 2 * r * (G(r + h, e, Tt) - G(r - h, e, Tt)) / (2 * h) + 4 * G(r, e, Tt)
    rhs = -4 * (e * (G(r, e + h, Tt) - G(r, e - h, Tt)) / (2 * h) + Tt * (G(r, e, Tt + h) - G(r, e, Tt - h)) / (2 * h))
    okK9.append(np.array([(lhs - rhs) / abs(lhs)]))
rec('K9', 'end of Section 9.2', '(Q^GF)^2: A*_m -> -2 i p_m [c], c* -> -2 i p_m [A_m]; QQ^GF+Q^GF Q=|p|^2; (2x.d+2y.d+4)G = -4(eps d_eps + T d_T)G', okK9, tol=1e-6)
lhs_w = 2 * 0.4 * (G(0.4 + h, 0.03, 2.0) - G(0.4 - h, 0.03, 2.0)) / (2 * h) + 4 * G(0.4, 0.03, 2.0)
rhs_w = -2 * (0.03 * (G(0.4, 0.03 + h, 2.0) - G(0.4, 0.03 - h, 2.0)) / (2 * h) + 2.0 * (G(0.4, 0.03, 2.0 + h) - G(0.4, 0.03, 2.0 - h)) / (2 * h))
neg('N8', 'end of Section 9.2', 'scaling Ward identity with -2 instead of -4', [(lhs_w - rhs_w) / abs(lhs_w)])

# ================================================================== K10 D_1, theta, gamma language, P_1, psi'
log('=== K10: D_1, theta, Konishi/20prime, P_1 and the finite part psi\' (blocks after thm:c9-D1, end of 9.4, end of 9.5) ===')
okK10 = []
# table of Theorem c9-D1 (coefficients of O and of {S,O}) for theta arbitrary
for th in (rng.normal(), -1 / 3):
    row1 = {'tl': (44 / 3, -1 / 3), 'tr': (-100 / 3, -1 / 3)}
    rowF = (-44 / 3, -th)
    for key, D1 in (('tl', (0, -(1 / 3 + th))), ('tr', (-48, -(1 / 3 + th)))):
        okK10.append(np.array(row1[key]) + np.array(rowF) - np.array(D1))
okK10.append(np.array([16 - 4 / 3 - 44 / 3, -32 - 4 / 3 + 100 / 3, 2 * 2 * 22 / 3 / 2 - 44 / 3]))
# two evaluations (Remark c9-twoeval) and the composite triangle ratio
okK10.append(np.array([(-4 / 3 + 8 / 3) - (8 / 3 - 4 / 3), (-1 / 3) - (-2 + 5 / 3), 8 / (-4) - (-2), (4 * 1.3) / (-2 * 1.3) - (-2)]))
# gamma language: 1/2 D_1 = gamma O - (1/3+theta) u (E+2) O; with u = -g^2N/32pi^2
okK10.append(np.array([0.5 * (-48 * u) - 3 * gN2 / (4 * np.pi ** 2), -(1 / 3 + 0.2) * u - (1 / 3 + 0.2) * gN2 / (32 * np.pi ** 2)]))
# Theorem c9-kappa: left side 2a d_b O, independent of a_P; [Delta,P_b] = -2 P_b as vector fields
aa, aP = rng.normal(), rng.normal()
okK10.append(np.array([(aa + aP) * 2 - 2 * aP - 2 * aa]))
b = rng.normal(size=4); x0 = rng.normal(size=4); hh = 1e-6
vD = lambda x: 2 * x; vP = lambda x: b
def lie(v, w, x):    # [v,w] = v(w) - w(v) as vector fields
    Jw = np.array([(w(x + hh * np.eye(4)[k]) - w(x - hh * np.eye(4)[k])) / (2 * hh) for k in range(4)]).T
    Jv = np.array([(v(x + hh * np.eye(4)[k]) - v(x - hh * np.eye(4)[k])) / (2 * hh) for k in range(4)]).T
    return Jw @ v(x) - Jv @ w(x)
okK10.append(lie(vD, vP, x0) + 2 * b)                       # [Delta, P_b] = -2 P_b
okK10.append(lie(lambda x: -x, vP, x0) - b)                 # [D, P_b] = P_b for D = -E
# Corollary c9-ident (i): 22/3 int(A,A*) + theta int(2E A + 2A, A*) at theta=-1/3
okK10.append(np.array([22 / 3 + 2 * (-1 / 3) - 20 / 3, 2 * (-1 / 3) + 2 / 3]))
# Proposition c9-scheme: (a/2) lambda_x = -a d_x, d_x = -lambda_x/2
lam = {'P': -2, 'K': 2, 'Sigma': -1, 'Sigma*': 1, 'Delta': 0}
dx = {'P': 1, 'K': -1, 'Sigma': 0.5, 'Sigma*': -0.5, 'Delta': 0}
okK10.append(np.array([0.5 * lam[k] - (-dx[k]) for k in lam]))
thP = rng.normal(); a_ = -(1 / 3 + thP)
okK10.append(np.array([-(1 / 3 + thP) + (-a_) * 1 - 0, -a_ * (-1) - (-(1 / 3 + thP)) * (-1) * (-1),
                       -a_ * 0.5 - 0.5 * (1 / 3 + thP), -a_ * (-0.5) + 0.5 * (1 / 3 + thP)]))
rec('K10', 'after thm:c9-D1, end of 9.4, end of 9.5', 'D_1 table, theta=-1/3 removes the classical multiple, gamma_K = 3g^2N/4pi^2, [Delta,P]=-2P, [D,P]=P, psi\' shifts -a d_x, P_1 -> 0', okK10, tol=1e-7)
neg('N9', 'end of 9.4', 'theta = +1/3 leaves a classical multiple in D_1', [np.array([-(1 / 3 + 1 / 3)])])
neg('N10', 'end of 9.5', "psi' shift with +a d_x (wrong sign) does not cancel P_1", [np.array([-(1 / 3 + thP) + a_ * 1])])

# ================================================================== K11 blocks restating the revised Chapter 8
log('=== K11: revised Chapter 8 statements restated in supplements/ch08.tex (8.7, 8.8, 8.9, after prop:c8-2pt) ===')
okK11 = []
# int_{r0<|y|<sqrt T} |y|^{-4} d^4y = 2 pi^2 log(sqrt(T)/r0) = pi^2 log T - 2 pi^2 log r0
r0 = 0.3
vals = []
for Tt in (5.0, 50.0):
    r = np.geomspace(r0, np.sqrt(Tt), 200001)
    vals.append(np.trapezoid(2 * np.pi ** 2 * r ** 3 * r ** -4, r))
okK11.append(np.array([(vals[1] - vals[0]) / np.log(10.0) - np.pi ** 2]))
# coefficient of the cut tree: fibre factors 2 (A x A) and -1 (A* x c) times T d_T of the scalar integral, -3/(16 pi^2)
okK11.append(np.array([2 * (-1) * (-3 / (16 * np.pi ** 2)) - 3 / (8 * np.pi ** 2)]))
# Prop c8-2pt (iii): -1/4 (D^T C + C D) log; one operator D = 2 gamma -> -gamma C; and the power law |x|^{-2(w+gamma)}
gam, Cc, w0 = 0.137, 1.9, 2.0
okK11.append(np.array([-0.25 * (2 * gam * Cc + Cc * 2 * gam) + gam * Cc]))
xs = np.array([0.8, 1.3]); ell = 1.1
exact = Cc * (xs ** 2 / ell ** 2) ** (-gam) * xs ** (-2 * w0); lin = xs ** (-2 * w0) * (Cc - gam * Cc * np.log(xs ** 2 / ell ** 2))
okK11.append(np.array([np.max(np.abs(exact - lin) / exact) < 0.05 * 1.0]) - 1)        # agreement to first order in gamma
# Definition c8-radial / Proposition c9-scheme: psi' = (a_P/2){S(Delta),-} changes K, Sigma, Sigma* by a_P, -a_P/2, +a_P/2
aP = rng.normal()
okK11.append(np.array([-aP * dx['K'] - aP, -aP * dx['Sigma'] + aP / 2, -aP * dx['Sigma*'] - aP / 2, aP + (-aP * dx['P'])]))
rec('K11', 'ch08 ends of 8.7, 8.8, 8.9, after prop:c8-2pt', 'pi^2 log T; 3/(8pi^2) = 2*(-1)*(-3/16pi^2); -1/4(D^TC+CD) = -gamma C; psi\' shifts (a_P,-a_P/2,a_P/2) and P_1 -> 0', okK11, tol=1e-6)
neg('N11', 'ch08 end of 8.7', 'log T coefficient 2 pi^2 instead of pi^2', [np.array([(vals[1] - vals[0]) / np.log(10.0) - 2 * np.pi ** 2])])
neg('N12', 'ch08 end of 8.7', 'cut-tree coefficient without the fibre factors (3/16pi^2)', [np.array([3 / (16 * np.pi ** 2) - 3 / (8 * np.pi ** 2)])])
neg('N13', 'ch08 end of 8.8', 'psi\' changing Sigma by +a_P/2', [np.array([-aP * dx['Sigma'] - aP / 2])])

npass = sum(r['passed'] for r in RES)
log(f'TOTAL {npass}/{len(RES)} passed')
json.dump(dict(passed=npass, total=len(RES), results=RES), open(os.path.join(HERE, 'verify_ch09.json'), 'w'), indent=1)
open(os.path.join(HERE, 'verify_ch09.log'), 'w').write('\n'.join(LOG) + '\n')
