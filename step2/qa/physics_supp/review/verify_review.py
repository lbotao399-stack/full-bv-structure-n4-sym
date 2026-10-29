#!/usr/bin/env python3
"""Independent re-verification of formulas of the 'Physics notation' supplements (review pass).

Self-contained: the two-component realisation is built literally from Definition c8-twocomp of the
manuscript (sigma_E, sigmabar_E, epsilon, 't Hooft blocks Sigma^I, rho on Sigma and Sigma^*, pairing),
the conjugation from Lemma c8-cc, the radial and reflection conjugations from Definition c8-dag,
Theorem c8-dag and Definition c8-Theta; the su(N) data from (X,Y) = -(4 hbar/g^2) tr(XY).
numpy only, small memory.  Output: verify_review.json, verify_review.log (same directory).
"""
import itertools, json, math, os
import numpy as np

np.seterr(all='ignore')   # Apple Accelerate emits spurious matmul warnings; every check is an explicit residual

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261029)
RES, LOG = [], []
TOL = 1e-9


def log(s):
    LOG.append(s)
    print(s)


def _flat(x):
    if isinstance(x, (list, tuple)):
        parts = [_flat(y) for y in x]
        return np.concatenate(parts) if parts else np.zeros(0, complex)
    return np.ravel(np.asarray(x, dtype=complex))


def rec(cid, where, desc, lhs, rhs=0.0, tol=TOL):
    lhs = _flat(lhs)
    rhs = np.ravel(np.asarray(rhs, dtype=complex)) if np.ndim(rhs) else np.full(lhs.shape, rhs, dtype=complex)
    scale = max(1.0, float(np.max(np.abs(lhs))) if lhs.size else 1.0)
    res = float(np.max(np.abs(lhs - rhs))) / scale if lhs.size else 0.0
    ok = res < tol
    RES.append(dict(id=cid, supplement=where, check=desc, residual=res, passed=bool(ok)))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:6s} res={res:.1e}  ({where}) {desc}")
    return ok


def cr(*s):
    return rng.normal(size=s) + 1j * rng.normal(size=s)


# ---------------------------------------------------------------- Definition c8-twocomp
I2 = np.eye(2, dtype=complex)
PA = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
SIG = [PA[0], PA[1], PA[2], 1j * I2]            # sigma^m_{alpha alphadot}
SIGB = [-PA[0], -PA[1], -PA[2], 1j * I2]        # sigmabar^{m alphadot alpha}
EU = np.array([[0, 1], [-1, 0]], complex)       # eps^{12}=+1
ED = np.array([[0, -1], [1, 0]], complex)       # eps_{12}=-1
E4 = np.zeros((4,) * 4)
for p in itertools.permutations(range(4)):
    E4[p] = np.linalg.det(np.eye(4)[list(p)])
levi3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    levi3[p] = np.linalg.det(np.eye(3)[list(p)])
SI = []
for a in range(3):
    M = np.zeros((4, 4), complex)
    M[:3, :3] = levi3[a]
    M[a, 3], M[3, a] = 1, -1
    SI.append(M)
for a in range(3):
    M = np.zeros((4, 4), complex)
    M[:3, :3] = 1j * levi3[a]
    M[a, 3], M[3, a] = -1j, 1j
    SI.append(M)
SIB = [S.conj().T for S in SI]                   # Sigmabar^I = (Sigma^I)^dagger, index (i,j) = Sigmabar^{I ij}


def up(wl):            # w^{ij} = 1/2 eps^{ijkl} w_{kl}
    return 0.5 * np.einsum('ijkl,kl->ij', E4, wl)


def wlow(w6):
    return sum(w6[I] * SI[I] for I in range(6))


# spinors: xi (2x4) = xi_alpha^i, xib (2x4) = xib_{alphadot i}; eta = eta^alpha_i, etab = etab^{alphadot i}
def rho_S(v10, xi, xib):
    """rho(v): Sigma -> Sigma^*, v10 = (v_1..v_4, w_1..w_6)."""
    v, w6 = v10[:4], v10[4:]
    wl = wlow(w6)
    wu = up(wl)
    eta = 1j * sum(v[m] * np.einsum('ai,ab->bi', xib, SIGB[m]) for m in range(4))
    etab = -1j * sum(v[m] * SIGB[m] @ xi for m in range(4))
    eta = eta + np.einsum('ab,ji,bi->aj', EU, wl, xi)
    etab = etab + np.einsum('ab,ji,bi->aj', EU, wu, xib)
    return eta, etab


def rho_Ss(v10, eta, etab):
    """rho(v): Sigma^* -> Sigma."""
    v, w6 = v10[:4], v10[4:]
    wl = wlow(w6)
    wu = up(wl)
    xi = -1j * sum(v[m] * SIG[m] @ etab for m in range(4))
    xib = 1j * sum(v[m] * np.einsum('ai,ab->bi', eta, SIG[m]) for m in range(4))
    xi = xi - np.einsum('kj,ab,bj->ak', wu, ED, eta)
    xib = xib - np.einsum('kj,ba,aj->bk', wl, ED, etab)
    return xi, xib


def pair(xi, xib, eta, etab):
    return np.sum(xi * eta) - np.sum(xib * etab)


def vecS(xi, xib):
    return np.concatenate([xi.ravel(), xib.ravel()])


def unS(x):
    return x[:8].reshape(2, 4), x[8:].reshape(2, 4)


def mat(f):
    cols = []
    for k in range(16):
        e = np.zeros(16, complex)
        e[k] = 1
        a, b = unS(e)
        cols.append(vecS(*f(a, b)))
    return np.array(cols).T


EM = [np.eye(10)[M] for M in range(10)]
RS = [mat(lambda a, b, M=M: rho_S(EM[M], a, b)) for M in range(10)]     # Sigma -> Sigma^*
RSs = [mat(lambda a, b, M=M: rho_Ss(EM[M], a, b)) for M in range(10)]   # Sigma^* -> Sigma
PAIR = np.array([[pair(*unS(np.eye(16)[r]), *unS(np.eye(16)[c])) for c in range(16)] for r in range(16)])  # (xi,eta)

log('=== A. Two-component realisation (key, ch01, ch08) ===')
rec('A1', 'key/ch08 def:c8-twocomp', 'Clifford rho(e_M)rho(e_N)+rho(e_N)rho(e_M)=2 delta_MN on Sigma and Sigma^*',
    [RSs[M] @ RS[N] + RSs[N] @ RS[M] - 2 * (M == N) * np.eye(16) for M in range(10) for N in range(10)] +
    [RS[M] @ RSs[N] + RS[N] @ RSs[M] - 2 * (M == N) * np.eye(16) for M in range(10) for N in range(10)])
rec('A2', 'ch08 def:c8-twocomp', '(rho(v)zeta1,zeta2)=(zeta1,rho(v)zeta2) with (xi,eta)=xi eta - xib etab',
    [(PAIR @ RS[M]).T - PAIR @ RS[M] for M in range(10)] + [(RSs[M].T @ PAIR).T - RSs[M].T @ PAIR for M in range(10)])
rec('A3', 'ch01 def:c1-spaces', 'u^{ik}w_{kj}+w^{ik}u_{kj}=-2(u,w) delta and w^{ij}=-w_I Sigmabar^{I ij}',
    [(lambda u, w: up(wlow(u)) @ wlow(w) + up(wlow(w)) @ wlow(u) + 2 * (u @ w) * np.eye(4))(cr(6), cr(6))] +
    [up(SI[I]) + SIB[I] for I in range(6)])
# chirality: rho(e1 e2 e3 e4)
P4_on_S = RSs[0] @ RS[1] @ RSs[2] @ RS[3]
P4_on_Ss = RS[0] @ RSs[1] @ RS[2] @ RSs[3]
rec('A4', 'key, ch01 end 1.1', 'rho(e1e2e3e4)=+1 on xi, -1 on xib (Sigma); +1 on eta, -1 on etab (Sigma^*)',
    [P4_on_S - np.diag([1] * 8 + [-1] * 8), P4_on_Ss - np.diag([1] * 8 + [-1] * 8)])
# Dirac column (xi_alpha^i, etab^{alphadot i}) at fixed i: rho(e_m) = [[0,-i sigma^m],[-i sigmabar^m,0]]
okD = []
for m in range(4):
    blk = np.zeros((4, 4), complex)
    blk[:2, 2:] = -1j * SIG[m]
    blk[2:, :2] = -1j * SIGB[m]
    for i in range(4):
        x = cr(2)
        xi = np.zeros((2, 4), complex); xi[:, i] = x
        eb = rho_S(EM[m], xi, np.zeros((2, 4)))[1][:, i]
        y = cr(2)
        etab = np.zeros((2, 4), complex); etab[:, i] = y
        xo = rho_Ss(EM[m], np.zeros((2, 4)), etab)[0][:, i]
        okD.append(blk @ np.concatenate([x, y]) - np.concatenate([xo, eb]))
D = [np.zeros((4, 4), complex) for _ in range(4)]
for m in range(4):
    D[m][:2, 2:] = -1j * SIG[m]; D[m][2:, :2] = -1j * SIGB[m]
okD.append(D[0] @ D[1] @ D[2] @ D[3] - np.diag([1, 1, -1, -1]))
rec('A5', 'ch01 end 1.1', 'Dirac-column form of rho(e_m) and rho(e1)..rho(e4)=diag(1,1,-1,-1)', okD)

# Gamma pairing from (rho(v)Q1,Q2)=(v,Gamma(Q1,Q2))
def Gam_S(x1, x2):
    return np.array([pair(*x2, *rho_S(EM[M], *x1)) for M in range(10)])


def Gam_Ss(e1, e2):
    return np.array([pair(*rho_Ss(EM[M], *e1), *e2) for M in range(10)])


x1, x2 = (cr(2, 4), cr(2, 4)), (cr(2, 4), cr(2, 4))
e1, e2 = (cr(2, 4), cr(2, 4)), (cr(2, 4), cr(2, 4))
g = Gam_S(x1, x2)
G1 = np.array([1j * (np.einsum('ai,ab,bi->', x2[1], SIGB[m], x1[0]) + np.einsum('ai,ab,bi->', x1[1], SIGB[m], x2[0]))
               for m in range(4)])
Gkl = sum(g[4 + I] * SI[I] for I in range(6))
G2 = (2 * np.einsum('jikl,aj,ab,bi->kl', E4, x1[0], EU, x2[0])
      - 2 * (np.einsum('ak,ab,bl->kl', x1[1], EU, x2[1]) - np.einsum('al,ab,bk->kl', x1[1], EU, x2[1])))
gs = Gam_Ss(e1, e2)
G3 = np.array([-1j * (np.einsum('ai,ab,bi->', e2[0], SIG[m], e1[1]) + np.einsum('ai,ab,bi->', e1[0], SIG[m], e2[1]))
               for m in range(4)])
Gup = up(sum(gs[4 + I] * SI[I] for I in range(6)))
G4 = (-2 * np.einsum('klpj,ab,bj,ap->kl', E4, ED, e1[0], e2[0])
      + 2 * (np.einsum('ba,al,bk->kl', ED, e1[1], e2[1]) - np.einsum('ba,ak,bl->kl', ED, e1[1], e2[1])))
rec('A6', 'ch01 def:c1-spaces, ch08 def:c8-twocomp', 'Gamma-pairing components G1,G2 (Sigma) and G3,G4 (Sigma^*)',
    np.concatenate([g[:4] - G1, (Gkl - G2).ravel(), gs[:4] - G3, (Gup - G4).ravel()]))

# rho(q(e_m ^ e_n)) on xi and xib (ch01)
okq = []
for m in range(4):
    for n in range(4):
        if m == n:
            continue
        xi, xib = cr(2, 4), cr(2, 4)
        a, b = rho_Ss(EM[m], *rho_S(EM[n], xi, xib))
        okq.append(a + 0.5 * (SIG[m] @ SIGB[n] - SIG[n] @ SIGB[m]) @ xi)
        okq.append(b - 0.5 * np.einsum('bi,bc->ci', xib, SIGB[m] @ SIG[n] - SIGB[n] @ SIG[m]))
rec('A7', 'ch01 def:c1-spaces, ch02 def:c2-bracket', 'rho(q(e_m^e_n)) = -1/2(s^m sb^n - s^n sb^m) on xi, +1/2 xib(sb^m s^n - sb^n s^m)', okq)

# 't Hooft blocks as SO(6) gamma matrices (ch08 def:c8-twocomp) and Gamma_I formulas
gam = []
for I in range(6):
    G = np.zeros((8, 8), complex)
    G[:4, 4:] = SIB[I]; G[4:, :4] = SI[I]
    gam.append(G)
rec('A8', 'ch08 def:c8-twocomp', 'gamma^I=[[0,Sigmabar^I],[Sigma^I,0]]: Clifford 2 delta and Hermitian',
    [gam[I] @ gam[J] + gam[J] @ gam[I] - 2 * (I == J) * np.eye(8) for I in range(6) for J in range(6)] +
    [gam[I] - gam[I].conj().T for I in range(6)])
okG = []
for I in range(6):
    xi, xib, eta, etab = cr(2, 4), cr(2, 4), cr(2, 4), cr(2, 4)
    o1, o2 = rho_S(EM[4 + I], xi, xib)
    okG += [o1 - np.einsum('ab,ji,bi->aj', EU, SI[I], xi), o2 + np.einsum('ab,ji,bi->aj', EU, SIB[I], xib)]
    o3, o4 = rho_Ss(EM[4 + I], eta, etab)
    okG += [o3 - np.einsum('kj,ab,bj->ak', SIB[I], ED, eta), o4 + np.einsum('kj,ba,aj->bk', SI[I], ED, etab)]
rec('A9', 'ch08 def:c8-twocomp', 'Gamma_I=rho(e_I) component formulas with Sigma^I, Sigmabar^I', okG)

log('=== B. Conjugations (ch08 lem:c8-cc, thm:c8-dag, def:c8-Theta) ===')
def cc_S(xi, xib):            # Sigma -> Sigma^*
    return np.conj(xi), -np.conj(xib)


def cc_Ss(eta, etab):         # Sigma^* -> Sigma (inverse, involutive)
    return np.conj(eta), -np.conj(etab)


okc = []
for M in range(10):
    xi, xib = cr(2, 4), cr(2, 4)
    lhs = cc_Ss(*rho_S(EM[M], xi, xib))
    rhs = rho_Ss(EM[M], *cc_S(xi, xib))
    okc += [lhs[0] - rhs[0], lhs[1] - rhs[1]]
xi, xib = cr(2, 4), cr(2, 4)
okc.append(pair(xi, xib, *cc_S(xi, xib)) - (np.sum(abs(xi) ** 2) + np.sum(abs(xib) ** 2)))
z3, z4 = (cr(2, 4), cr(2, 4)), (cr(2, 4), cr(2, 4))
okc.append(Gam_Ss(cc_S(*z3), cc_S(*z4)) - np.conj(Gam_S(z3, z4)))
rec('B1', 'ch08 lem:c8-cc', 'conj(rho(e_M)zeta)=rho(e_M)conj(zeta); (xi,conj xi)=sum|xi|^2+sum|xib|^2; Gamma(conj,conj)=conj Gamma', okc)


def unitQ(a, i):
    x = np.zeros((2, 4), complex); x[a, i] = 1; return x, np.zeros((2, 4), complex)


def unitQb(ad, i):
    x = np.zeros((2, 4), complex); x[ad, i] = 1; return np.zeros((2, 4), complex), x


# reflection n=e_4: xi^{dagger_n} = rho(n) conj(xi)   (Definition c8-Theta, eta=0)
okr, okr2 = [], []
for a in range(2):
    for i in range(4):
        out = rho_Ss(EM[3], *cc_S(*unitQ(a, i)))                    # element of Sigma
        tgt = unitQb(a, i)                                          # unit Qbar_{alphadot=alpha, i}
        okr += [out[0] + tgt[0], out[1] + tgt[1]]                   # = - Qbar
        out = rho_Ss(EM[3], *cc_S(*unitQb(a, i)))
        tgt = unitQ(a, i)
        okr += [out[0] + tgt[0], out[1] + tgt[1]]                   # = - Q
for _ in range(3):
    n = rng.normal(size=4); n /= np.linalg.norm(n)
    nv = np.concatenate([n, np.zeros(6)])
    xi, xib = cr(2, 4), cr(2, 4)
    out = rho_Ss(nv, *cc_S(xi, xib))
    f1 = 1j * sum(n[m] * SIG[m] @ np.conj(xib) for m in range(4))
    f2 = 1j * sum(n[m] * np.einsum('ai,ab->bi', np.conj(xi), SIG[m]) for m in range(4))
    okr2 += [out[0] - f1, out[1] - f2]
rec('B2', 'ch08 def:c8-Theta', 'two-component (xi^{dagger_n}) = (i n.sigma conj(xib), i conj(xi) n.sigma) from rho(n)conj(xi)', okr2)
rec('B3', 'ch08 def:c8-Theta (suspect)', '(Q^i_alpha)^{dagger e4} = -Qbar_{alphadot i}|_{alphadot=alpha}, (Qbar)^{dagger e4} = -Q|_{alpha=alphadot}', okr)

# radial: xi^dagger = -conj(xi)
okd = []
for a in range(2):
    for i in range(4):
        e = cc_S(*unitQ(a, i)); okd += [-e[0] + unitQ(a, i)[0], -e[1]]         # -S^alpha_i  (unit S: eta[a,i]=1)
        e = cc_S(*unitQb(a, i)); okd += [-e[0], -e[1] - unitQb(a, i)[1]]       # +Sbar^{alphadot i}
rec('B4', 'ch08 thm:c8-dag', '(Q^i_alpha)^dagger = -S^alpha_i, (Qbar_{alphadot i})^dagger = +Sbar^{alphadot i}', okd)


# dilatation part of [xi, rho(x) eta]: x -> Gamma(xi, rho(x)eta)|_V = (xi,eta) E + rotation
def linfield(xs, es):
    Mx = np.zeros((4, 4), complex)
    for nn in range(4):
        rx = rho_Ss(EM[nn], *es)
        Mx[:, nn] = Gam_S(xs, rx)[:4]
    return Mx


okdil = []
for _ in range(3):
    xs, es = (cr(2, 4), cr(2, 4)), (cr(2, 4), cr(2, 4))
    Mx = linfield(xs, es)
    c0 = np.trace(Mx) / 4
    okdil += [c0 - pair(*xs, *es), (Mx - c0 * np.eye(4)) + (Mx - c0 * np.eye(4)).T]
rec('B5', 'ch02 end 2.1, ch08 end 8.1', 'Gamma(xi,rho(x)eta)|_V = (xi,eta) x + antisymmetric: dilatation part 1/2 (xi,eta) Delta', okdil)
# [Q, Q^dagger] contains -1/2 Delta; [Q^i_a, S^b_j] dil part 1/2 delta delta Delta; [Qb, Sb] -1/2 delta delta Delta
okqd = []
for a in range(2):
    for i in range(4):
        xs = unitQ(a, i); es = tuple(-np.asarray(t) for t in cc_S(*xs))     # (Q)^dagger as element of Sigma^*
        okqd.append(np.trace(linfield(xs, es)) / 4 - (-1))                  # coefficient of E = -1: -E = -1/2 Delta
        xs = unitQb(a, i); es = tuple(-np.asarray(t) for t in cc_S(*xs))
        okqd.append(np.trace(linfield(xs, es)) / 4 - (-1))
        for b in range(2):
            for j in range(4):
                S_unit = (np.eye(2)[b][:, None] * np.eye(4)[j][None, :]).astype(complex), np.zeros((2, 4), complex)
                okqd.append(np.trace(linfield(unitQ(a, i), S_unit)) / 4 - (a == b) * (i == j))
                Sb_unit = np.zeros((2, 4), complex), (np.eye(2)[b][:, None] * np.eye(4)[j][None, :]).astype(complex)
                okqd.append(np.trace(linfield(unitQb(a, i), Sb_unit)) / 4 + (a == b) * (i == j))
rec('B6', 'ch08 end 8.1, thm:c8-dag', '[Q,Q^dagger],[Qb,Qb^dagger] contain -1/2 Delta; dil. parts of [Q,S], [Qb,Sb] = +-1/2 delta delta Delta '
    '(coefficient of E checked: (xi,eta))', okqd)
# [Q^i_a, Qbar_{ad j}] = i delta^i_j sigmabar^{m ad a} P_m ; [Q,Q]=0 on V
okqq = []
for a in range(2):
    for i in range(4):
        for ad in range(2):
            for j in range(4):
                gq = Gam_S(unitQ(a, i), unitQb(ad, j))[:4]
                okqq.append(gq - np.array([1j * (i == j) * SIGB[m][ad, a] for m in range(4)]))
        for b in range(2):
            for j in range(4):
                okqq.append(Gam_S(unitQ(a, i), unitQ(b, j))[:4])
rec('B7', 'ch08 end 8.1, ch02 end 2.1', '[Q^i_alpha,Qbar_{alphadot j}] = i delta^i_j sigmabar^{m alphadot alpha} P_m, [Q,Q]=0', okqq)

log('=== C. Superconformal data (ch02, ch07, ch08) ===')
okx = []
for _ in range(3):
    x = rng.normal(size=4)
    xv = np.concatenate([x, np.zeros(6)])
    e1, e2 = (cr(2, 4), cr(2, 4)), (cr(2, 4), cr(2, 4))
    lhs = Gam_S(rho_Ss(xv, *e1), rho_Ss(xv, *e2))
    G0 = Gam_Ss(e1, e2)
    v1 = G0[:4]
    okx.append(lhs[:4] + (x @ x) * v1 - 2 * (v1 @ x) * x)          # V part = -(x^2 v1 - 2(v1.x)x)
    okx.append(lhs[4:] + (x @ x) * G0[4:])                          # W part = -x^2 Gamma^I
    A_ij = wlow(cr(6)); A_m = cr(4)
    AI = np.array([0.25 * np.einsum('ji,ij->', SIB[I], A_ij) for I in range(6)])
    okx.append(lhs[:4] @ A_m + lhs[4:] @ AI - (-(x @ x * v1 - 2 * (v1 @ x) * x) @ A_m - 0.25 * (x @ x) * np.sum(up(sum(G0[4 + I] * SI[I] for I in range(6))) * A_ij)))
rec('C1', 'ch02 end 2.4, end 2.1', 'Gamma(rho(x)eta1,rho(x)eta2): gauge parameter -(x^2v1-2(v1.x)x)A - 1/4 x^2 Gamma^{kl}A_kl', okx)
# pure spinor
okp = []
for _ in range(4):
    u, z = cr(2), cr(4)
    xs = (np.outer(u, z), np.zeros((2, 4), complex))
    okp.append(Gam_S(xs, xs))
rec('C2', 'ch08 lem:c8-ev', 'xi=u(x)z in S_+(x)W is pure: Gamma^M(xi,xi)=0 for all ten M', okp)
# R-symmetry
okR = []
for _ in range(3):
    R = rng.normal(size=(6, 6)); R = R - R.T
    r = 0.25 * sum(R[I, J] * SIB[J] @ SI[I] for I in range(6) for J in range(6))
    okR += [np.trace(r), r + r.conj().T]
    xi, xib = cr(2, 4), cr(2, 4)
    tot = [np.zeros((2, 4), complex), np.zeros((2, 4), complex)]
    for I in range(6):
        for J in range(6):
            a, b = rho_Ss(EM[4 + J], *rho_S(EM[4 + I], xi, xib))
            tot[0] += 0.25 * R[I, J] * a; tot[1] += 0.25 * R[I, J] * b
    okR += [tot[0] - xi @ r.T, tot[1] + xib @ r]                       # r^i_j xi^j ; -xib_j r^j_i
    A6 = cr(6)
    Al = wlow(A6)
    LvA = wlow(np.array([sum(R[K, I] * A6[K] for K in range(6)) for I in range(6)]))
    okR.append(LvA - (-(r.T @ Al) - Al @ r))                           # -r^k_i A_kj - r^k_j A_ik
    okR.append(up(LvA) - (r @ up(Al) + up(Al) @ r.T))                  # r^i_k A^kj + r^j_k A^ik
rec('C3', 'ch02 def:c2-even, def:c2-S1, ch07 def:c7-vdot', 'r=1/4 R_IJ Sb^J S^I traceless anti-Hermitian; = 1/4 rho(d(v,-)) on xi, xib; v.A_ij, v.A^ij', okR)
# Killing spinor Dirac derivative
eta, etab = cr(2, 4), cr(2, 4)
dEps = [-1j * SIG[m] @ etab for m in range(4)]                         # d_m eps_alpha^i
dEpb = [1j * np.einsum('ai,ab->bi', eta, SIG[m]) for m in range(4)]    # d_m epsbar_{alphadot i}
rec('C4', 'ch07 prop:c7-select', 'i d_m epsbar sigmabar^m = 4 eta, -i sigmabar^m d_m eps = 4 etab',
    [1j * sum(np.einsum('ai,ab->bi', dEpb[m], SIGB[m]) for m in range(4)) - 4 * eta,
     -1j * sum(SIGB[m] @ dEps[m] for m in range(4)) - 4 * etab])
# conformal Killing data of K_b
b = rng.normal(size=4); x = rng.normal(size=4)
dv = 2 * np.outer(x, b) - 2 * np.outer(b, x) - 2 * (b @ x) * np.eye(4)   # d_m v_n
X = dv - dv.T
rec('C5', 'ch02 def:c2-even, ch07 prop:c7-dil/cor:c7-K', 'K_b: X_mn=4(x_m b_n - x_n b_m), div v=-8 b.x, d_m v_n+d_n v_m=1/2 div v delta',
    [X - 4 * (np.outer(x, b) - np.outer(b, x)), np.trace(dv) + 8 * b @ x, dv + dv.T - 0.5 * np.trace(dv) * np.eye(4)])

# R-symmetry part of [xi, rho(x)eta]: (e_I, rot e_J) = (rho(q(e_I^e_J))xi, eta); r = 1/4 R_IJ Sb^J S^I vs (B2)
okB2 = []
for _ in range(3):
    xs, es = (cr(2, 4), cr(2, 4)), (cr(2, 4), cr(2, 4))
    R = np.zeros((6, 6), complex)
    for I in range(6):
        for J in range(6):
            if I != J:
                R[I, J] = pair(*rho_Ss(EM[4 + I], *rho_S(EM[4 + J], *xs)), *es)
    r = 0.25 * sum(R[I, J] * SIB[J] @ SI[I] for I in range(6) for J in range(6))
    rB2 = 2 * (np.einsum('ak,ai->ki', xs[0], es[0]) + np.einsum('ai,ak->ki', xs[1], es[1]))
    rB2 = rB2 - np.trace(rB2) / 4 * np.eye(4)
    okB2 += [R + R.T, r - rB2]
rec('C6', 'ch02 def:c2-bracket', 'R-symmetry rotation of [xi,rho(x)eta] in SU(4) form: r^k_i = 2(xi^k eta_i + xib_i etab^k) - trace', okB2)

log('=== D. su(N) normalisation, potentials, Yukawa, bilinears ===')
def suN(N):
    T = []
    for i in range(N):
        for j in range(i + 1, N):
            M = np.zeros((N, N), complex); M[i, j] = M[j, i] = 0.5; T.append(M)
            M = np.zeros((N, N), complex); M[i, j] = -0.5j; M[j, i] = 0.5j; T.append(M)
    for k in range(1, N):
        M = np.zeros((N, N), complex)
        M[:k, :k] = np.eye(k); M[k, k] = -k
        T.append(M / np.sqrt(2 * k * (k + 1)))
    T = np.array(T)
    f = np.real(-2j * np.einsum('aij,bjk,cki->abc', T, T, T) + 2j * np.einsum('bij,ajk,cki->abc', T, T, T))
    return T, f


hb, gc = 0.7, 1.3
okN = []
for N in (2, 3, 4):
    T, f = suN(N)
    d = len(T)
    pr = lambda X, Y: -4 * hb / gc ** 2 * np.trace(X @ Y)
    Xa, Ya = rng.normal(size=d), rng.normal(size=d)
    X = -1j * gc * np.einsum('a,aij->ij', Xa, T); Y = -1j * gc * np.einsum('a,aij->ij', Ya, T)
    TA = -1j * f                                                       # (T^a_A)^{bc} = -i f^{abc}
    adX = -1j * gc * np.einsum('a,abc->bc', Xa, TA); adY = -1j * gc * np.einsum('a,abc->bc', Ya, TA)
    # check adX is the matrix of [X,-] in basis -ig T^a
    Zb = rng.normal(size=d); Z = -1j * gc * np.einsum('a,aij->ij', Zb, T)
    comm = X @ Z - Z @ X
    comp = np.array([np.real(2 * np.trace(comm @ T[a]) / (-1j * gc)) for a in range(d)])
    okN.append(comp - adX @ Zb)
    unit = hb * np.trace(adX @ adY) / (16 * np.pi ** 2 * pr(X, Y))
    okN += [pr(X, Y) - 2 * hb * Xa @ Ya, np.trace(adX @ adY) + gc ** 2 * N * Xa @ Ya,
            np.trace(adX @ adY) - 2 * N * np.trace(X @ Y), unit + gc ** 2 * N / (32 * np.pi ** 2),
            -24 * unit - 3 * gc ** 2 * N / (4 * np.pi ** 2), -4 * unit - gc ** 2 * N / (8 * np.pi ** 2),
            -8 * unit - gc ** 2 * N / (4 * np.pi ** 2), -2 * unit - gc ** 2 * N / (16 * np.pi ** 2),
            unit * 9 / 4 * 2 * hb / hb + 9 * gc ** 2 * N / (64 * np.pi ** 2),
            27 / 16 * unit * 2 * hb + 27 * gc ** 2 * N / (256 * np.pi ** 2) * hb]
    # quartic colour factor and scalar potential
    X4 = [rng.normal(size=d) for _ in range(4)]
    Xm = [-1j * gc * np.einsum('a,aij->ij', x4, T) for x4 in X4]
    cm = lambda A, B: A @ B - B @ A
    okN.append(pr(cm(Xm[0], Xm[1]), cm(Xm[2], Xm[3])) - 2 * hb * gc ** 2 * np.einsum('abe,cde,a,b,c,d->', f, f, *X4))
    Aa = rng.normal(size=(6, d))
    AI = [-1j * gc * np.einsum('a,aij->ij', Aa[I], T) for I in range(6)]
    V1 = -1 / 8 * sum(pr(cm(AI[I], AI[J]), cm(AI[I], AI[J])) for I in range(6) for J in range(6)) / hb
    V2 = -gc ** 2 / 4 * np.einsum('abe,cde,Ia,Jb,Ic,Jd->', f, f, Aa, Aa, Aa, Aa)
    Alow = np.einsum('Ikl,Iab->klab', np.array(SI), np.array(AI))       # A_{kl} matrix-valued
    Aup = 0.5 * np.einsum('ijkl,klab->ijab', E4, Alow)
    V3 = -1 / 128 * sum(pr(cm(Aup[i, j], Aup[k, l]), cm(Alow[i, j], Alow[k, l]))
                       for i in range(4) for j in range(4) for k in range(4) for l in range(4)) / hb
    okN += [V1 - V2, V1 - V3]
    # (A_I,A_J) = 4 hbar tr(phi phi) = 2 hbar phi^a phi^a with A_I = -ig phi_I
    phi = [np.einsum('a,aij->ij', Aa[I], T) for I in range(6)]
    okN.append(pr(AI[0], AI[1]) - 4 * hb * np.trace(phi[0] @ phi[1]))
    okN.append(pr(AI[0], AI[1]) - 2 * hb * Aa[0] @ Aa[1])
rec('D1', 'key, ch04 4.1/thm:c4-quadratic/4.3.2, ch05, ch06, ch01 def:c1-action',
    '(X,Y)=2hbar X^aY^a, tr ad ad=-g^2N X.Y=2N tr XY, unit=-g^2N/32pi^2, Konishi 3g^2N/4pi^2, letters, mass, 27/16, quartic colour, potential (SO(6),SU(4),f f forms)', okN)


# Yukawa su(N) form via Grassmann coefficient matrices (ch01 def:c1-action)
def antisym(C):
    return C - C.T


okY = []
for N in (2, 3):
    T, f = suN(N); d = len(T)
    Ab = rng.normal(size=(6, d))                                       # A^b_I
    Alow_c = np.einsum('Ikl,Ib->klb', np.array(SI), Ab)               # A^b_{kl}
    Aup_c = 0.5 * np.einsum('ijkl,klb->ijb', E4, Alow_c)              # A^{b kl}
    trTTT = np.einsum('aij,bjk,cki->abc', T, T, T)
    comm3 = trTTT - np.einsum('aij,cjk,bki->abc', T, T, T)            # tr(T^a [T^b,T^c])
    pref = -4 * hb / gc ** 2 * (-1j * gc / np.sqrt(2)) * (-1j * gc) * (-1j * gc / np.sqrt(2)) / hb
    # LHS1: 1/2 (lambda_alpha^j, eps^{alpha beta} [A_ji, lambda_beta^i]) / hbar ; index x=(a,alpha,j)
    # ordering: tr(T^d [T^c, T^e]) with lambda^d (first), A^c, lambda^e (second)
    L1 = 0.5 * pref * np.einsum('ab,jic,dce->dajebi', EU, Alow_c, comm3)
    L1 = L1.reshape(d * 2 * 4, d * 2 * 4)
    # RHS1: -g/2 f^{abc} A^b_{ji} eps^{alpha beta} theta_{a beta j} theta_{c alpha i}
    R1 = np.zeros((d, 2, 4, d, 2, 4), complex)
    for a in range(d):
        for c in range(d):
            for j in range(4):
                for i in range(4):
                    coef = -gc / 2 * np.einsum('b,b->', f[a, :, c], Alow_c[j, i, :])
                    for al in range(2):
                        for be in range(2):
                            R1[a, be, j, c, al, i] += coef * EU[al, be]
    R1 = R1.reshape(d * 8, d * 8)
    okY.append(antisym(L1) - antisym(R1))
    # LHS2: -1/2 (lambdabar_{ad j}, eps^{ad bd} [A^{ji}, lambdabar_{bd i}]) / hbar
    L2 = -0.5 * pref * np.einsum('ab,jic,dce->dajebi', EU, Aup_c, comm3).reshape(d * 8, d * 8)
    R2 = np.zeros((d, 2, 4, d, 2, 4), complex)
    for a in range(d):
        for c in range(d):
            for j in range(4):
                for i in range(4):
                    coef = -gc / 2 * np.einsum('b,b->', f[a, :, c], Aup_c[j, i, :])
                    for ad in range(2):
                        for bd in range(2):
                            R2[a, ad, j, c, bd, i] += coef * EU[ad, bd]
    R2 = R2.reshape(d * 8, d * 8)
    okY.append(antisym(L2) - antisym(R2))
rec('D2', 'ch01 def:c1-action', 'Yukawa in su(N) components: -g/2 f A_ji (lam^j lam^i) - g/2 f A^ji (lamb_j lamb_i), Srednicki contractions', okY)

# 20' and Konishi forms (ch05 def:c5-bilinear, ch06 rem:c6-konishi)
ok5 = []
T, f = suN(3)
ph = [np.einsum('a,aij->ij', rng.normal(size=8), T) for _ in range(6)]
w1, w2 = rng.normal(size=6), rng.normal(size=6)
phl = np.einsum('Ikl,Iab->klab', np.array(SI), np.array(ph))
phu = 0.5 * np.einsum('ijkl,klab->ijab', E4, phl)
w1u, w2u = up(wlow(w1)), up(wlow(w2))
lhs = 4 * hb * sum(w1[I] * w2[J] * np.trace(ph[I] @ ph[J] - (I == J) / 6 * sum(ph[K] @ ph[K] for K in range(6)))
                   for I in range(6) for J in range(6))
trmm = sum(np.trace(phu[m, n] @ phl[m, n]) for m in range(4) for n in range(4))
rhs = hb / 4 * sum(w1u[i, j] * w2u[k, l] * (np.trace(phl[i, j] @ phl[k, l]) - E4[i, j, k, l] / 12 * trmm)
                   for i in range(4) for j in range(4) for k in range(4) for l in range(4))
ok5 += [lhs - rhs, 4 * hb * sum(np.trace(ph[I] @ ph[I]) for I in range(6)) - hb * trmm,
        0.25 * np.sum(w1u * wlow(w2)) - w1 @ w2]
# one-loop dilatation operator on the 21 bilinears: tr(phi_I phi_J) -> g^2N/8pi^2 delta_IJ tr(phi_K phi_K)
idx = [(I, J) for I in range(6) for J in range(I, 6)]
Mdil = np.zeros((21, 21))
for c, (I, J) in enumerate(idx):
    if I == J:
        for r_, (K, L) in enumerate(idx):
            if K == L:
                Mdil[r_, c] = 1.0
ev = np.sort(np.linalg.eigvals(Mdil).real)
ok5.append(ev - np.array([0.0] * 20 + [6.0]))
rec('D3', "ch05 def:c5-bilinear, ch06 rem:c6-konishi", "20' SO(6)=SU(4) form, Konishi 4hbar tr phi_I phi_I = hbar tr phi^ij phi_ij, mixing eigenvalues 0 (x20), 6 (Konishi)", ok5)
# propagator SO(6) form
P = np.einsum('jiI,lkJ,ijkl->IJ', np.array(SIB).transpose(1, 2, 0), np.array(SIB).transpose(1, 2, 0), 4 * E4) / 16
rec('D4', 'ch03 prop:c3-P, ch08 prop:c8-twocomp', '<A_ij A_kl> = 4 eps_ijkl (-,-)^{-1}  <=>  <A_I A_J> = 2 delta_IJ (-,-)^{-1}', P - 2 * np.eye(6))

log('=== E. Kernels, logarithms, Feynman-gauge term ===')
# (a) proper-time integral of the heat kernel
okE = []
for r_, eps in ((0.7, 0.05), (1.3, 0.4)):
    u = np.linspace(1e-12, r_ ** 2 / (4 * eps), 200001)              # t = r^2/(4u)
    integrand = (16 * u ** 2 / (16 * np.pi ** 2 * r_ ** 4)) * np.exp(-u) * r_ ** 2 / (4 * u ** 2)
    val = np.trapezoid(integrand, u)
    okE.append(val - (1 - np.exp(-r_ ** 2 / (4 * eps))) / (4 * np.pi ** 2 * r_ ** 2))
rec('E1', 'ch03 def:c3-K/prop:c3-P, ch00 0.1', 'int_eps^inf (4 pi t)^-2 e^{-r^2/4t} dt = (1-e^{-r^2/4eps})/(4pi^2 r^2)', okE, tol=1e-7)
# (b) -Laplacian of G_eps = K_eps (radial, 4d), and G_eps(0)=1/(16 pi^2 eps)
eps = 0.3
G = lambda r: (1 - np.exp(-r ** 2 / (4 * eps))) / (4 * np.pi ** 2 * r ** 2)
rs = np.array([0.2, 0.6, 1.1, 2.0]); h = 1e-4
lap = (G(rs + h) - 2 * G(rs) + G(rs - h)) / h ** 2 + 3 / rs * (G(rs + h) - G(rs - h)) / (2 * h)
rec('E2', 'ch03 end 3.2, prop:c3-P', '-d^2 (1-e^{-r^2/4eps})/(4pi^2r^2) = e^{-r^2/4eps}/(4 pi eps)^2 ; value 1/(16pi^2 eps) at r=0',
    np.concatenate([-lap - np.exp(-rs ** 2 / (4 * eps)) / (4 * np.pi * eps) ** 2, [G(1e-4) - 1 / (16 * np.pi ** 2 * eps)]]), tol=1e-5)
# (c) position-space bubble: d/dlog eps of int_{r<R} d^4x G_eps^2 -> -1/(16 pi^2)
def bub(eps, R=1.0):
    r = np.linspace(1e-9, R, 400001)
    return np.trapezoid(2 * np.pi ** 2 * r ** 3 * ((1 - np.exp(-r ** 2 / (4 * eps))) / (4 * np.pi ** 2 * r ** 2)) ** 2, r)
e1_, e2_ = 1e-5, 2e-5
slope = (bub(e2_) - bub(e1_)) / (np.log(e2_) - np.log(e1_))
rec('E3', 'ch05 end 5.2', 'position-space bubble (G_eps)^2 contains -log(eps)/(16 pi^2) delta(x-y)', slope + 1 / (16 * np.pi ** 2), tol=1e-4)
# (d) momentum-space bubble via Schwinger parameters: eps d/deps at small eps -> -1/(16pi^2)
def mom_bub_deriv(eps, p2=1.0):
    # eps d/deps of int_{[eps,inf)^2} dt1 dt2 e^{-t1 t2 p2/(t1+t2)}/(16 pi^2 (t1+t2)^2) = -2 eps int_eps^inf dt ...
    t = np.geomspace(eps, 1e4, 400001)
    integ = np.exp(-eps * t * p2 / (eps + t)) / (16 * np.pi ** 2 * (eps + t) ** 2)
    return -2 * eps * np.trapezoid(integ, t)
rec('E4', 'ch05 end 5.2, Lemma c5-log', 'momentum bubble e^{-eps p^2}/p^2 x e^{-eps(p+p1)^2}/(p+p1)^2 has log-coefficient -1/(16pi^2)',
    mom_bub_deriv(1e-7) + 1 / (16 * np.pi ** 2), tol=1e-5)
# (e) Prop c4-K pointwise: 2 dA.d(L_vA) - d_k(v_k (dA)^2) = 8 (d.A)(b.A), v = x^2 b - 2(b.x)x
a0 = rng.normal(size=4); B1 = rng.normal(size=(4, 4)); C2 = rng.normal(size=(4, 4, 4))
C2 = C2 + C2.transpose(0, 2, 1)
bb = rng.normal(size=4)
def Afield(x):
    return a0 + B1 @ x + 0.5 * np.einsum('nkl,k,l->n', C2, x, x)
def dA(x):                    # dA[n,k] = d_k A_n
    return B1 + np.einsum('nkl,l->nk', C2, x)
ddA = C2                      # d_k d_l A_n
def vK(x):
    return (x @ x) * bb - 2 * (bb @ x) * x
def dvK(x):                   # dv[n,k] = d_k v_n
    return 2 * np.outer(bb, x) - 2 * np.outer(x, bb) - 2 * (bb @ x) * np.eye(4)
okK = []
for _ in range(5):
    x = rng.normal(size=4)
    divA = np.trace(dA(x))
    grad_divA = np.einsum('nnk->k', ddA)
    # (L_v A)_n = v^k d_k A_n + A_k d_n v^k ; divergence d_n (L_vA)_n
    # exact second derivatives of v_n = x^2 b_n - 2 (b.x) x_n : d_k d_l v_n = 2 b_n delta_kl - 2 b_k delta_nl - 2 b_l delta_nk
    d2v = 2 * np.einsum('n,kl->nkl', bb, np.eye(4)) - 2 * np.einsum('k,nl->nkl', bb, np.eye(4)) - 2 * np.einsum('l,nk->nkl', bb, np.eye(4))
    divLvA = (np.einsum('kn,nk->', dvK(x), dA(x)) + np.einsum('k,nnk->', vK(x), ddA)
              + np.einsum('kn,kn->', dA(x), dvK(x)) + np.einsum('k,knn->', Afield(x), d2v))
    # d_k (v_k (divA)^2) = (div v)(divA)^2 + 2 v_k divA d_k divA
    tot = 2 * divA * divLvA - (np.trace(dvK(x)) * divA ** 2 + 2 * divA * vK(x) @ grad_divA)
    okK.append(tot - 8 * divA * (bb @ Afield(x)))
rec('E5', 'ch04 end 4.4 (Prop c4-K)', '2(d.A)d.(L_vA) = d_k(v_k (d.A)^2) + 8 (d.A)(b.A) for the special conformal v', okK)
# (f) Feynman-gauge Wick map: d^mu A_mu = d_m A_m for A_0 = i A_4, d_0 = i d_4 (signature -+++)
A4, d4 = cr(1)[0], cr(1)[0]
rec('E6', 'ch03 end 3.1', 'Wick: -d_0 A_0 + ... = d_4 A_4 with A_0=iA_4, d_0=i d_4', -(1j * d4) * (1j * A4) - d4 * A4)

log('=== N. Negative controls (must be detected as nonzero) ===')
def neg(cid, desc, val):
    r = float(np.max(np.abs(_flat(val))))
    ok = r > 1e-3
    RES.append(dict(id=cid, supplement='negative control', check=desc, residual=r, passed=bool(ok)))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:6s} |dev|={r:.2e}  (negative control) {desc}")
neg('N1', 'reflection with + sign, (Q)^{dagger e4} = +Qbar, is rejected',
    [rho_Ss(EM[3], *cc_S(*unitQ(0, 0)))[1] - unitQb(0, 0)[1]])
neg('N2', 'Yukawa with +g/2 instead of -g/2 is rejected', [antisym(L1) + antisym(R1)])
neg('N3', 'scalar potential with -1/64 instead of -1/128 in SU(4) form is rejected', [V1 - 2 * V3])
neg('N4', '[Q,Q^dagger] with +1/2 Delta is rejected', [np.trace(linfield(unitQ(0, 0), tuple(-np.asarray(t) for t in cc_S(*unitQ(0, 0))))) / 4 - 1])
neg('N5', "anomalous dimension from twice the log-coefficient of Z=1-4u log eps (instead of the counterterm) is nonzero for the 20'",
    [2 * (-4) - 4 - 4])

npass = sum(r['passed'] for r in RES)
log(f'TOTAL {npass}/{len(RES)} passed')
json.dump(dict(passed=npass, total=len(RES), results=RES), open(os.path.join(HERE, 'verify_review.json'), 'w'), indent=1)
open(os.path.join(HERE, 'verify_review.log'), 'w').write('\n'.join(LOG) + '\n')
