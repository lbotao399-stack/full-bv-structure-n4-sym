#!/usr/bin/env python3
"""Numerical checks of the two-component formulas of the physics-notation supplements of Chapters 0-2
that are not verbatim in the dictionary or in Chapter 8.

Realisation: the two-component realisation D of the dictionary (copied from step2/notation/check_dictionary.py,
functions rho_S, rho_Ss, pair_D; Euclidean sigma_E, sigmabar_E; 't Hooft blocks of engine B).
The index-free side of every check is evaluated with the 16x16 matrices of D (Clifford action, pairing, Gamma),
the two-component side literally as typeset in the supplements.  numpy only, < 50 MB.
Output: verify_supp.log, verify_supp.json in this directory.
"""
import os, json, itertools
import numpy as np

np.seterr(all='ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261029)
TOL = 1e-9
RES, LOG = [], []


def log(s):
    LOG.append(s); print(s)


def flat(x):
    if isinstance(x, (list, tuple)):
        return np.concatenate([flat(y) for y in x]) if x else np.zeros(0, complex)
    return np.ravel(np.asarray(x, dtype=complex))


def record(cid, desc, lhs, rhs=0.0):
    lhs = flat(lhs); rhs = flat(rhs)
    scale = max(1.0, float(np.max(np.abs(lhs))) if lhs.size else 1.0, float(np.max(np.abs(rhs))) if rhs.size else 1.0)
    res = float(np.max(np.abs(lhs - rhs))) / scale if lhs.size else 0.0
    ok = res < TOL
    RES.append(dict(id=cid, check=desc, residual=res, passed=ok))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:6s} res={res:.2e}  {desc}")
    return ok


def rnd(*shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


# ---------------------------------------------------------------- realisation D (dictionary, = check_dictionary.py)
I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EPSU = np.array([[0, 1], [-1, 0]], complex)     # eps^{12} = +1
EPSD = np.array([[0, -1], [1, 0]], complex)     # eps_{12} = -1
SIGE = [PAULI[0], PAULI[1], PAULI[2], 1j * I2]          # sigma^m_{alpha alphadot}
SIGBE = [-PAULI[0], -PAULI[1], -PAULI[2], 1j * I2]      # sigmabar^{m alphadot alpha}


def thooft(selfdual):
    out = []
    for a in range(3):
        M = np.zeros((4, 4), complex)
        for m in range(3):
            for n in range(3):
                M[m, n] = np.linalg.det(np.eye(3)[[a, m, n]]) if len({a, m, n}) == 3 else 0
        s = 1 if selfdual else -1
        M[a, 3] = s; M[3, a] = -s
        out.append(M)
    return out


_eta = thooft(True); _etab = thooft(False)
SIGI = [_eta[0], _eta[1], _eta[2], 1j * _etab[0], 1j * _etab[1], 1j * _etab[2]]       # Sigma^I_{ij}
SIGBI = [-_eta[0], -_eta[1], -_eta[2], 1j * _etab[0], 1j * _etab[1], 1j * _etab[2]]  # Sigmabar^{I ij}
EPS4 = np.zeros((4,) * 4)
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])


def split(x):
    return x[:8].reshape((2, 4) + x.shape[1:]), x[8:].reshape((2, 4) + x.shape[1:])


def join(a, b):
    return np.concatenate([a.reshape((8,) + a.shape[2:]), b.reshape((8,) + b.shape[2:])])


def rho_S(M, x):          # Sigma -> Sigma^*  (R1),(R3)
    q, qb = split(x)
    if M < 4:
        s = 1j * np.einsum('di...,da->ai...', qb, SIGBE[M])
        sb = -1j * np.einsum('da,ai...->di...', SIGBE[M], q)
    else:
        I = M - 4
        s = np.einsum('ab,ji,bi...->aj...', EPSU, SIGI[I], q)
        sb = -np.einsum('ji,de,ei...->dj...', SIGBI[I], EPSU, qb)
    return join(s, sb)


def rho_Ss(M, y):         # Sigma^* -> Sigma  (R2),(R4)
    s, sb = split(y)
    if M < 4:
        q = -1j * np.einsum('ad,di...->ai...', SIGE[M], sb)
        qb = 1j * np.einsum('ai...,ad->di...', s, SIGE[M])
    else:
        I = M - 4
        q = np.einsum('kj,ab,bj...->ak...', SIGBI[I], EPSD, s)
        qb = -np.einsum('kj,ed,dj...->ek...', SIGI[I], EPSD, sb)
    return join(q, qb)


E16 = np.eye(16)
RD = [np.array([rho_S(M, E16[:, k]) for k in range(16)]).T for M in range(10)]     # Sigma -> Sigma^*
RbD = [np.array([rho_Ss(M, E16[:, k]) for k in range(16)]).T for M in range(10)]   # Sigma^* -> Sigma


def pair_D(x, y):
    q, qb = split(x); s, sb = split(y)
    return np.sum(q * s) - np.sum(qb * sb)


PD = np.array([[pair_D(E16[:, a], E16[:, b]) for b in range(16)] for a in range(16)])

log('=== 0. realisation D (sanity: dictionary R1, R2, R4) ===')
record('D1', 'Clifford on Sigma and Sigma^*',
       [RbD[M] @ RD[N] + RbD[N] @ RD[M] - 2 * (M == N) * E16 for M in range(10) for N in range(10)] +
       [RD[M] @ RbD[N] + RD[N] @ RbD[M] - 2 * (M == N) * E16 for M in range(10) for N in range(10)])
record('D2', '(rho(v)Q1,Q2) = (Q1,rho(v)Q2)', [PD @ RD[M] - (PD @ RD[M]).T for M in range(10)])
pp = np.zeros(16); pp[:8] = 1
record('D3', 'rho(e1e2e3e4) = +1 on S_+ (x) W, -1 on S_- (x) W^*', RbD[0] @ RD[1] @ RbD[2] @ RD[3] - np.diag(2 * pp - 1))


# ---------------------------------------------------------------- colour: su(N) matrices, (X,Y) = -(4 hbar/g^2) tr(XY)
def su_gens(N):
    T = []
    for a in range(N):
        for b in range(a + 1, N):
            M = np.zeros((N, N), complex); M[a, b] = M[b, a] = 0.5; T.append(M)
            M = np.zeros((N, N), complex); M[a, b] = -0.5j; M[b, a] = 0.5j; T.append(M)
    for k in range(1, N):
        M = np.zeros((N, N), complex)
        for a in range(k):
            M[a, a] = 1
        M[k, k] = -k
        T.append(M / np.sqrt(2 * k * (k + 1)))
    return T


HBAR = 0.731; GC = 1.37


def gpair(X, Y):          # (X,Y) on g, matrices with arbitrary leading axes contracted
    return -(4 * HBAR / GC ** 2) * np.sum(X * np.swapaxes(Y, -1, -2))


def comm(X, Y):
    return X @ Y - Y @ X


def colour_setup(N):
    T = su_gens(N)
    n = len(T)
    f = np.array([[[(-2j * np.trace(comm(T[a], T[b]) @ T[c])) for c in range(n)] for b in range(n)] for a in range(n)])
    return T, n, f.real if np.max(np.abs(f.imag)) < 1e-12 else f


def to_mat(comp, T, coef):    # X = coef * comp^a T^a, comp has trailing colour axis
    return coef * np.einsum('...a,aij->...ij', comp, np.array(T))


def bis_lo(u):    # u_ij = u_I Sig^I_ij  (u_I on axis 0)
    return np.einsum('I...,Iij->ij...', u, np.array(SIGI))


def bis_up(u):    # u^ij = -u_I Sigb^{I ij} = 1/2 eps^{ijkl} u_kl
    return -np.einsum('I...,Iij->ij...', u, np.array(SIGBI))


# ---------------------------------------------------------------- 1. Clifford data (ch01 after def:c1-spaces, end 1.1)
log('=== 1. rho o q on V, Dirac operator, chirality, canonical element ===')
for m in range(4):
    for n in range(4):
        rq = RbD[m] @ RD[n] - (m == n) * E16           # rho(q(e_m ^ e_n)) on Sigma
        lhs = []; rhs = []
        for k in range(3):
            x = rnd(16); q, qb = split(x)
            y = rq @ x; yq, yqb = split(y)
            Lp = -0.5 * (SIGE[m] @ SIGBE[n] - SIGE[n] @ SIGBE[m])          # (sigma^m sigmabar^n - ...)_alpha^beta
            Lm = 0.5 * (SIGBE[m] @ SIGE[n] - SIGBE[n] @ SIGE[m])           # (sigmabar^m sigma^n - ...)^{betadot}_{alphadot}
            lhs += [yq, yqb]
            rhs += [np.einsum('ab,bi->ai', Lp, q), np.einsum('bi,ba->ai', qb, Lm)]
        record('C1', f'rho(q(e_{m+1}^e_{n+1})): xi -> -1/2(s^m sb^n - s^n sb^m) xi, xibar -> 1/2 xibar(sb^m s^n - sb^n s^m)', lhs, rhs)

# Dirac operator in components and dslash(rho(x)eta) = 4 eta with (K1)
eta = rnd(16); s, sb = split(eta)
# eps(x) = rho(x) eta: eps_a^i = -i x_m sig^m_{a ad} sb^{ad i}, epsbar_{ad i} = i x_m s^a_i sig^m_{a ad}
xr = rnd(4)
eK = -1j * sum(xr[m] * np.einsum('ad,di->ai', SIGE[m], sb) for m in range(4))
ebK = 1j * sum(xr[m] * np.einsum('ai,ad->di', s, SIGE[m]) for m in range(4))
record('C2a', '(K1) = rho(x) eta in D', join(eK, ebK), sum(xr[m] * RbD[m] for m in range(4)) @ eta)
# derivatives of the linear functions: d_m eps_a^i = -i sig^m sb, d_m epsbar = i s sig^m
dE = [-1j * np.einsum('ad,di->ai', SIGE[m], sb) for m in range(4)]
dEb = [1j * np.einsum('ai,ad->di', s, SIGE[m]) for m in range(4)]
# dslash xi: (dslash xi)^a_i = i d_m xibar_{ad i} sbar^{m ad a}, (dslash xi)^{ad i} = -i sbar^{m ad a} d_m xi_a^i
dsl_up = sum(1j * np.einsum('di,da->ai', dEb[m], SIGBE[m]) for m in range(4))
dsl_dn = sum(-1j * np.einsum('da,ai->di', SIGBE[m], dE[m]) for m in range(4))
record('C2b', 'dslash rho(x)eta = 4 eta in two-component form', join(dsl_up, dsl_dn), 4 * eta)
dsl_IF = sum(RD[m] @ join(dE[m], dEb[m]) for m in range(4))
record('C2c', 'dslash components = rho(e_m) d_m (index-free)', join(dsl_up, dsl_dn), dsl_IF)
# dslash on Sigma^*: (dslash eta)_a^i = -i sig^m_{a ad} d_m etabar^{ad i}, (dslash eta)_{ad i} = i d_m eta^a_i sig^m_{a ad}
dS = [rnd(2, 4) for m in range(4)]; dSb = [rnd(2, 4) for m in range(4)]
dslS_lo = sum(-1j * np.einsum('ad,di->ai', SIGE[m], dSb[m]) for m in range(4))
dslS_dn = sum(1j * np.einsum('ai,ad->di', dS[m], SIGE[m]) for m in range(4))
record('C2d', 'dslash|_{Sigma^*} components = rho(e_m) d_m', join(dslS_lo, dslS_dn), sum(RbD[m] @ join(dS[m], dSb[m]) for m in range(4)))

# chirality on the Dirac column (xi_alpha^i, etabar^{alphadot i}) at fixed i, and on all four blocks
Z32 = np.zeros((32, 32), complex)
GAM = []
for m in range(4):
    G = Z32.copy(); G[16:, :16] = RD[m]; G[:16, 16:] = RbD[m]; GAM.append(G)   # on (Sigma, Sigma^*)
vol = GAM[0] @ GAM[1] @ GAM[2] @ GAM[3]
diag_exp = np.concatenate([np.ones(8), -np.ones(8), np.ones(8), -np.ones(8)])   # xi, xibar, eta, etabar
record('C3a', 'rho(e1)rho(e2)rho(e3)rho(e4) = +1 on xi_a^i, eta^a_i; -1 on xibar_{ad i}, etabar^{ad i}', vol - np.diag(diag_exp))
g5S = np.diag([-1, -1, 1, 1]).astype(complex); g5W = -g5S
for i in range(4):
    idx = [0 * 4 + i, 1 * 4 + i, 16 + 8 + 0 * 4 + i, 16 + 8 + 1 * 4 + i]       # (xi_1^i, xi_2^i, etabar^{1 i}, etabar^{2 i})
    col = vol[np.ix_(idx, idx)]
    record('C3b', f'Dirac column i={i+1}: product = -gamma5_S = gamma5_W = diag(1,1,-1,-1)', [col - g5W, g5W + g5S, col - np.diag([1, 1, -1, -1])])
# canonical element of Sigma (x) Sigma^*: C = sum_a e_a (x) e^a, (e_a, e^b) = delta
Cmat = np.linalg.inv(PD).T          # C[a, c]: coefficient of e_a (x) f_c
Cexp = np.diag(np.concatenate([np.ones(8), -np.ones(8)]))
record('C4', 'canonical element: delta_a^b delta^i_j on undotted, -delta delta on dotted blocks', Cmat - Cexp)
# the canonical element reproduces the pairing: sum_a (x, e^a) e_a = x
xt = rnd(16)     # x = sum_{a,c} C[a,c] (x, f_c) e_a  for x in Sigma
record('C4b', 'canonical element reconstructs Sigma: x = sum C[a,c] (x,f_c) e_a', Cmat @ (xt @ PD), xt)


# ---------------------------------------------------------------- 2. action pieces (after def:c1-action, end 1.4, end 1.6)
log('=== 2. action: field strength components, SO(6) forms, su(N) potential and Yukawa, vertices ===')
for N in (2, 3):
    T, nG, f = colour_setup(N)
    # random complexified fields at a point: A_M (10 matrices), dA[m][M] = d_m A_M
    A = [to_mat(rnd(nG), T, 1.0) for _ in range(10)]
    dA = [[to_mat(rnd(nG), T, 1.0) for _ in range(10)] for _ in range(4)]
    Am = A[:4]; AI = A[4:]
    F = [[None] * 10 for _ in range(10)]
    for M in range(10):
        for Nn in range(10):
            dMN = dA[M][Nn] if M < 4 else 0 * A[0]
            dNM = dA[Nn][M] if Nn < 4 else 0 * A[0]
            F[M][Nn] = dMN - dNM + comm(A[M], A[Nn])
    hodge = sum(gpair(F[M][Nn], F[M][Nn]) for M in range(10) for Nn in range(M + 1, 10))
    Fmn = [[F[m][n] for n in range(4)] for m in range(4)]
    DAI = [[dA[m][4 + I] + comm(Am[m], AI[I]) for I in range(6)] for m in range(4)]
    so6 = (0.5 * sum(gpair(Fmn[m][n], Fmn[m][n]) for m in range(4) for n in range(4))
           + sum(gpair(DAI[m][I], DAI[m][I]) for m in range(4) for I in range(6))
           + 0.5 * sum(gpair(comm(AI[I], AI[J]), comm(AI[I], AI[J])) for I in range(6) for J in range(6)))
    record('A1', f'su({N}): (F_A,F_A) Hodge = 1/2(F_mn,F_mn) + (D_mA_I,D_mA_I) + 1/2([A_I,A_J],[A_I,A_J]); F_mI = D_mA_I, F_IJ = [A_I,A_J]', so6, hodge)
    # dictionary SU(4) form
    AIarr = np.array(AI)
    lo = bis_lo(AIarr); up = bis_up(AIarr)
    Dlo = [bis_lo(np.array(DAI[m])) for m in range(4)]; Dup = [bis_up(np.array(DAI[m])) for m in range(4)]
    su4 = (0.5 * sum(gpair(Fmn[m][n], Fmn[m][n]) for m in range(4) for n in range(4))
           + 0.25 * sum(gpair(Dup[m][i, j], Dlo[m][i, j]) for m in range(4) for i in range(4) for j in range(4))
           + (1 / 32) * sum(gpair(comm(up[i, j], up[k, l]), comm(lo[i, j], lo[k, l])) for i in range(4) for j in range(4) for k in range(4) for l in range(4)))
    record('A2', f'su({N}): SO(6) form = dictionary SU(4) form of (F,F)', so6, su4)
    # su(N) potential: -1/8 sum([A_I,A_J],[A_I,A_J]) / hbar = -(g^2/4) f^{abe} f^{cde} A^a_I A^b_J A^c_I A^d_J
    Ac = rnd(6, nG)                                    # physical components A^a_I
    AIs = [to_mat(Ac[I], T, -1j * GC) for I in range(6)]
    pot_IF = -0.125 * sum(gpair(comm(AIs[I], AIs[J]), comm(AIs[I], AIs[J])) for I in range(6) for J in range(6)) / HBAR
    pot_2C = -(GC ** 2 / 4) * np.einsum('abe,cde,Ia,Jb,Ic,Jd->', f, f, Ac, Ac, Ac, Ac)
    record('A3', f'su({N}): potential -1/8([A_I,A_J],[A_I,A_J])/hbar = -(g^2/4) f f A A A A', pot_2C, pot_IF)
    # su(N) Yukawa, Grassmann by polarisation: lambda = th1 u1 + th2 u2, physical components u^a
    u1 = rnd(16, nG); u2 = rnd(16, nG)
    U1 = to_mat(u1, T, -1j * GC / np.sqrt(2)); U2 = to_mat(u2, T, -1j * GC / np.sqrt(2))

    def yuk_IF(X, Y):     # 1/2 (X, rho(e_I)[A_I, Y]) summed over I, X,Y in Sigma (x) g (matrices)
        out = 0
        for I in range(6):
            CY = np.array([comm(AIs[I], Y[a]) for a in range(16)])
            RY = np.einsum('xy,yij->xij', RD[4 + I], CY)        # in Sigma^* (x) g
            q, qb = split(X); s_, sb_ = split(RY)
            out += 0.5 * (gpair(q, s_) - gpair(qb, sb_))
        return out
    vIF = (yuk_IF(U1, U2) - yuk_IF(U2, U1)) / HBAR
    Alo = bis_lo(Ac)          # A^a_{ij} = A^a_I Sig^I_ij  (ij, a)
    Aup = bis_up(Ac)          # A^{a ij}

    def yuk_2C(x1, x2):       # -(g/2) f^{abc} A^b_{ji} (lambda^{a j} lambda^{c i}) - (g/2) f^{abc} A^{b ji} (lambdabar^a_j lambdabar^c_i)
        q1, qb1 = split(x1); q2, qb2 = split(x2)          # q[alpha, i, a], qb[alphadot, i, a]
        # (lambda^{a j} lambda^{c i}) = lambda^{a alpha j} lambda^c_alpha^i = eps^{alpha beta} lambda^a_beta^j lambda^c_alpha^i
        t1 = -(GC / 2) * np.einsum('abc,jib,AB,Bja,Aic->', f, Alo, EPSU, q1, q2)
        # (lambdabar^a_j lambdabar^c_i) = lambdabar^a_{ad j} lambdabar^{c ad}_i = lambdabar^a_{ad j} eps^{ad bd} lambdabar^c_{bd i}
        t2 = -(GC / 2) * np.einsum('abc,jib,AB,Aja,Bic->', f, Aup, EPSU, qb1, qb2)
        return t1 + t2
    v2C = yuk_2C(u1, u2) - yuk_2C(u2, u1)
    record('A4', f'su({N}): Yukawa (1/2)(lambda,rho(A_W)lambda)/hbar = -(g/2)fA_ji(l^j l^i) - (g/2)fA^ji(lb_j lb_i) (Grassmann)', v2C, vIF)
    # cubic and quartic parts of I
    dAf = [[(dA[M][Nn] if M < 4 else 0 * A[0]) - (dA[Nn][M] if Nn < 4 else 0 * A[0]) for Nn in range(10)] for M in range(10)]
    cub_IF = -0.5 * sum(gpair(dAf[M][Nn], comm(A[M], A[Nn])) for M in range(10) for Nn in range(M + 1, 10))
    cub_2C = (-0.5 * sum(gpair(dA[m][n], comm(Am[m], Am[n])) for m in range(4) for n in range(4))
              - 0.5 * sum(gpair(dA[m][4 + I], comm(Am[m], AI[I])) for m in range(4) for I in range(6)))
    record('A5', f'su({N}): -1/2(dA,1/2[A^A]) = -1/2(d_mA_n,[A_m,A_n]) - 1/2(d_mA_I,[A_m,A_I])', cub_2C, cub_IF)
    qua_IF = -0.25 * sum(gpair(comm(A[M], A[Nn]), comm(A[M], A[Nn])) for M in range(10) for Nn in range(M + 1, 10))
    qua_2C = (-0.125 * sum(gpair(comm(Am[m], Am[n]), comm(Am[m], Am[n])) for m in range(4) for n in range(4))
              - 0.25 * sum(gpair(comm(Am[m], AI[I]), comm(Am[m], AI[I])) for m in range(4) for I in range(6))
              - 0.125 * sum(gpair(comm(AI[I], AI[J]), comm(AI[I], AI[J])) for I in range(6) for J in range(6)))
    record('A6', f'su({N}): quartic -1/4(1/2[A^A],1/2[A^A]) = -1/8 AmAn - 1/4 AmAI - 1/8 AIAJ', qua_2C, qua_IF)
    # lambda lambda A vertex (commuting fibre elements): (l1, rho([A,l2])) summed over M
    L1 = np.array([to_mat(rnd(nG), T, 1.0) for _ in range(16)]); L2 = np.array([to_mat(rnd(nG), T, 1.0) for _ in range(16)])
    vtx_IF = 0
    for M in range(10):
        CY = np.array([comm(A[M], L2[a]) for a in range(16)])
        RY = np.einsum('xy,yij->xij', RD[M], CY)
        q, qb = split(L1); s_, sb_ = split(RY)
        vtx_IF += gpair(q, s_) - gpair(qb, sb_)
    q1, qb1 = split(L1); q2, qb2 = split(L2)       # q[alpha, i, N, N]
    vtx_2C = 0
    for m in range(4):
        Cqb2 = np.array([[comm(Am[m], qb2[d, i]) for i in range(4)] for d in range(2)])
        Cq2 = np.array([[comm(Am[m], q2[a, i]) for i in range(4)] for a in range(2)])
        vtx_2C += 1j * sum(SIGBE[m][d, a] * gpair(q1[a, i], Cqb2[d, i]) for a in range(2) for d in range(2) for i in range(4))
        vtx_2C += 1j * sum(SIGBE[m][d, a] * gpair(qb1[d, i], Cq2[a, i]) for a in range(2) for d in range(2) for i in range(4))
    Alo_m = bis_lo(np.array(AI)); Aup_m = bis_up(np.array(AI))
    vtx_2C += sum(EPSU[a, b] * gpair(q1[a, j], comm(Alo_m[j, i], q2[b, i])) for a in range(2) for b in range(2) for i in range(4) for j in range(4))
    vtx_2C -= sum(EPSU[a, b] * gpair(qb1[a, j], comm(Aup_m[j, i], qb2[b, i])) for a in range(2) for b in range(2) for i in range(4) for j in range(4))
    record('A7', f'su({N}): lambda-lambda-A vertex (l1,rho([A,l2])) in two-component form', vtx_2C, vtx_IF)
    # rho(F_A) xi and the supersymmetry transformations xi.A = -Gamma(xi,lambda), xi.lambda = -1/2 rho(F_A) xi
    xi = rnd(16)
    rF = np.zeros((16, N, N), complex)
    for M in range(10):
        for Nn in range(M + 1, 10):
            rF += np.einsum('x,ij->xij', RbD[M] @ RD[Nn] @ xi, F[M][Nn])
    xq, xqb = split(xi)
    Fm = np.array([[Fmn[m][n] for n in range(4)] for m in range(4)])
    low = bis_lo(np.array(AI)); upp = bis_up(np.array(AI))
    DAlo = [bis_lo(np.array(DAI[m])) for m in range(4)]; DAup = [bis_up(np.array(DAI[m])) for m in range(4)]
    # (rho(F_A)xi)_a^i (dictionary F1) and xi.lambda = -1/2 of it (supplement, end of 2.3)
    xl_up = np.zeros((2, 4, N, N), complex); xl_dn = np.zeros((2, 4, N, N), complex)
    for a in range(2):
        for i in range(4):
            t = 0
            for m in range(4):
                for n in range(4):
                    SS = SIGE[m] @ SIGBE[n]
                    t = t + 0.25 * Fm[m, n] * sum(SS[a, b] * xq[b, i] for b in range(2))
            t = t + 0.25 * sum(comm(upp[i, k], low[k, j]) * xq[a, j] for k in range(4) for j in range(4))
            t = t + 0.5j * sum(SIGE[m][a, d] * EPSU[d, e] * DAup[m][i, j] * xqb[e, j] for m in range(4) for d in range(2) for e in range(2) for j in range(4))
            xl_up[a, i] = t
    for d in range(2):
        for i in range(4):
            t = 0
            for m in range(4):
                for n in range(4):
                    BS = SIGBE[m] @ SIGE[n]
                    t = t - 0.25 * Fm[m, n] * sum(xqb[e, i] * BS[e, d] for e in range(2))
            t = t + 0.25 * sum(comm(low[i, k], upp[k, j]) * xqb[d, j] for k in range(4) for j in range(4))
            t = t - 0.5j * sum(EPSU[a, b] * SIGE[m][a, d] * DAlo[m][i, j] * xq[b, j] for m in range(4) for a in range(2) for b in range(2) for j in range(4))
            xl_dn[d, i] = t
    record('A8', f'su({N}): xi.lambda components = -1/2 rho(F_A) xi', join(xl_up, xl_dn), -0.5 * rF)
    # xi.A = -Gamma(xi, lambda): Gamma^M(Q1,Q2) = (Q1, rho(e_M) Q2), lambda commuting here (single odd slot)
    lam = rnd(16)
    GamIF = np.array([pair_D(xi, RD[M] @ lam) for M in range(10)])
    lq, lqb = split(lam)
    xA_m = np.array([-1j * (np.einsum('di,da,ai->', lqb, SIGBE[m], xq) + np.einsum('di,da,ai->', xqb, SIGBE[m], lq)) for m in range(4)])
    xA_kl = (-2 * np.einsum('jikl,aj,ab,bi->kl', EPS4, xq, EPSU, lq)
             + 2 * (np.einsum('ak,ab,bl->kl', xqb, EPSU, lqb) - np.einsum('al,ab,bk->kl', xqb, EPSU, lqb)))
    record('A9', 'xi.A_m = -i(lb sb^m xi + xib sb^m l), xi.A_kl = -Gamma_kl(xi,lambda)', [xA_m, xA_kl], [-GamIF[:4], -bis_lo(GamIF[4:])])

# ---------------------------------------------------------------- 3. pairings, propagator and bilinear normalisations (ch00, ch01)
log('=== 3. su(N) normalisation of the propagator and of the scalar bilinears ===')
for N in (2, 3):
    T, nG, f = colour_setup(N)
    for coef, want, name in ((-1j * GC, 2 * HBAR, 'A_M=-igA^aT^a'), (-1j * GC / np.sqrt(2), HBAR, 'lambda=-(ig/sqrt2)lambda^aT^a')):
        E = [coef * T[a] for a in range(nG)]
        G = np.array([[gpair(E[a], E[b]) for b in range(nG)] for a in range(nG)])
        record('P1', f'su({N}): (e_a,e_b) = {want/HBAR:g} hbar delta_ab for {name}', G, want * np.eye(nG))
    E = [-1j * GC * T[a] for a in range(nG)]
    G = np.array([[gpair(E[a], E[b]) for b in range(nG)] for a in range(nG)]); Ginv = np.linalg.inv(G)
    record('P2', f'su({N}): hbar * 2(-,-)^-1 = delta^ab (A components); hbar*(-,-)^-1 = delta^ab (lambda components)',
           [HBAR * 2 * Ginv, HBAR * np.linalg.inv(G / 2)], [np.eye(nG), np.eye(nG)])
# <A_I A_J> = 2(-,-)^{-1} delta_IJ reproduces the dictionary entries for A_ij
cI = 2.0
lolo = sum(cI * np.einsum('ij,kl->ijkl', SIGI[I], SIGI[I]) for I in range(6))
loup = sum(cI * np.einsum('ij,kl->ijkl', SIGI[I], -SIGBI[I]) for I in range(6))
upup = sum(cI * np.einsum('ij,kl->ijkl', -SIGBI[I], -SIGBI[I]) for I in range(6))
DD = np.einsum('ik,jl->ijkl', np.eye(4), np.eye(4)) - np.einsum('il,jk->ijkl', np.eye(4), np.eye(4))
record('P3', '<A_I A_J> = 2(-,-)^-1 delta_IJ  <=>  <A_ij A_kl> = 4 eps_ijkl, <A_ij A^kl> = 4(dd-dd), <A^ij A^kl> = 4 eps^ijkl',
       [lolo, loup, upup], [4 * EPS4, 4 * DD, 4 * EPS4])
# scalar bilinears: (i_w1 A, i_w2 A) = w1_I w2_J (A_I,A_J) = 2 hbar w1_I w2_J A^a_I A^a_J = 1/16 (w1^ij A_ij, w2^kl A_kl)
for N in (2, 3):
    T, nG, f = colour_setup(N)
    w1 = rnd(6); w2 = rnd(6); Ac = rnd(6, nG)
    AIm = [to_mat(Ac[I], T, -1j * GC) for I in range(6)]
    iw1 = sum(w1[I] * AIm[I] for I in range(6)); iw2 = sum(w2[I] * AIm[I] for I in range(6))
    lhs = gpair(iw1, iw2)
    Alo = bis_lo(np.array(AIm)); w1up = bis_up(w1); w2up = bis_up(w2)
    su4 = (1 / 16) * gpair(np.einsum('ij,ijab->ab', w1up, Alo), np.einsum('kl,klab->ab', w2up, Alo))
    comp = 2 * HBAR * np.einsum('I,J,Ia,Ja->', w1, w2, Ac, Ac)
    record('P4', f'su({N}): (i_w1A,i_w2A) = 2hbar w1_I w2_J A^a_I A^a_J = 1/16(w1^ij A_ij, w2^kl A_kl); (w1,w2) = w1_I w2_I',
           [lhs, lhs, 0.25 * np.einsum('ij,ij->', w1up, bis_lo(w2))], [comp, su4, np.sum(w1 * w2)])
    # Konishi: sum_w (i_w A, i_w A) = (A_I,A_I) = 1/4 (A^ij, A_ij) = 2 hbar A^a_I A^a_I
    Aup = bis_up(np.array(AIm))
    record('P5', f'su({N}): (A_I,A_I) = 1/4(A^ij,A_ij) = 2hbar A^a_I A^a_I = -(4hbar/g^2) tr(A_I A_I)',
           [sum(gpair(AIm[I], AIm[I]) for I in range(6))] * 3,
           [0.25 * sum(gpair(Aup[i, j], Alo[i, j]) for i in range(4) for j in range(4)), 2 * HBAR * np.sum(Ac * Ac),
            -(4 * HBAR / GC ** 2) * sum(np.trace(AIm[I] @ AIm[I]) for I in range(6))])

# ---------------------------------------------------------------- 4. psl(4|4) items (ch02)
log('=== 4. psl(4|4): rotations, S S gauge parameter, special superconformal scalar term, closure ===')
# rotation v^m = omega_mn x_n: 1/4 rho(d(v,-)) = 1/8 X_mn rho(q(e_m^e_n)), X = -2 omega
om = rnd(4, 4); om = om - om.T
X = -2 * om
rot = sum(0.125 * X[m, n] * (RbD[m] @ RD[n] - (m == n) * E16) for m in range(4) for n in range(4))
x0 = rnd(16); q, qb = split(x0); y = rot @ x0; yq, yqb = split(y)
Rp = 0.25 * sum(om[m, n] * SIGE[m] @ SIGBE[n] for m in range(4) for n in range(4))
Rm = -0.25 * sum(om[m, n] * SIGBE[m] @ SIGE[n] for m in range(4) for n in range(4))
record('K1', 'rotation omega: 1/4 rho(d(v,-)) = 1/4 omega_mn sigma^m sigmabar^n on xi, -1/4 omega_mn xibar sigmabar^m sigma^n',
       [yq, yqb], [np.einsum('ab,bi->ai', Rp, q), np.einsum('bi,ba->ai', qb, Rm)])
# Gamma(rho(x)eta1, rho(x)eta2): V part -(x^2 v1 - 2(v1.x)x), W part -x^2 Gamma(eta1,eta2)|_W, with (G3),(G4)
e1 = rnd(16); e2 = rnd(16); xr = rnd(4)
rx = sum(xr[m] * RbD[m] for m in range(4))
GxIF = np.array([pair_D(rx @ e1, RD[M] @ (rx @ e2)) for M in range(10)])
s1, sb1 = split(e1); s2, sb2 = split(e2)
v1 = np.array([-1j * (np.einsum('ai,ad,di->', s2, SIGE[m], sb1) + np.einsum('ai,ad,di->', s1, SIGE[m], sb2)) for m in range(4)])
G4up = (-2 * np.einsum('klpj,ab,bj,ap->kl', EPS4, EPSD, s1, s2)
        + 2 * (np.einsum('ba,al,bk->kl', EPSD, sb1, sb2) - np.einsum('ba,ak,bl->kl', EPSD, sb1, sb2)))
x2 = np.sum(xr * xr)
Vpart = -(x2 * v1 - 2 * np.sum(v1 * xr) * xr)
# check (G3) v1 = Gamma(eta1,eta2)|_V with D: Gamma on Sigma^* = (rho(e_M) eta1, eta2)
GsIF = np.array([pair_D(RbD[M] @ e1, e2) for M in range(10)])
record('S1', '(G3): Gamma(eta1,eta2)|_V^m = -i(eta2 s^m etab1 + eta1 s^m etab2)', v1, GsIF[:4])
record('S2', '(G4): Gamma(eta1,eta2)^{kl} = -Gamma^I Sigmabar^{I kl}', G4up, bis_up(GsIF[4:]))
record('S3', 'Gamma(rho(x)eta1,rho(x)eta2): V part -(x^2 v1 - 2(v1.x)x), W part -x^2 Gamma(eta1,eta2)|_W',
       [GxIF[:4], GxIF[4:]], [Vpart, -x2 * GsIF[4:]])
# gauge parameter: Gamma^m A_m + 1/4 Gamma_ij A^ij = Gamma^m A_m + Gamma^I A_I ; for S S: -(x^2v1-2(v1.x)x)^m A_m - 1/4 x^2 Gamma^kl(eta1,eta2) A_kl
Aw = rnd(10)
lhs = GxIF[:4] @ Aw[:4] + 0.25 * np.einsum('ij,ij->', bis_lo(GxIF[4:]), bis_up(Aw[4:]))
rhs = Vpart @ Aw[:4] - 0.25 * x2 * np.einsum('kl,kl->', G4up, bis_lo(Aw[4:]))
record('S4', 'S S gauge parameter in components; 1/4 Gamma_ij A^ij = Gamma^I A_I', [lhs, 0.25 * np.einsum('ij,ij->', bis_lo(GxIF[4:]), bis_up(Aw[4:]))],
       [rhs, GxIF[4:] @ Aw[4:]])
# special superconformal scalar term: x.lambda contains rho(A_W) eta; -(rho(A_W)eta, lambda^*) = (F2) extra terms
for N in (2, 3):
    T, nG, f = colour_setup(N)
    AIm = np.array([to_mat(rnd(nG), T, 1.0) for _ in range(6)])
    lst = np.array([to_mat(rnd(nG), T, 1.0) for _ in range(16)])     # lambda^* in Sigma^* (x) g
    et = rnd(16)
    rAe = sum(np.einsum('x,ij->xij', RbD[4 + I] @ et, AIm[I]) for I in range(6))   # rho(A_W) eta in Sigma (x) g
    q, qb = split(rAe); s_, sb_ = split(lst)
    IF = -(gpair(q, s_) - gpair(qb, sb_))
    es, esb = split(et); Alo = bis_lo(AIm); Aup = bis_up(AIm)
    F2 = sum(gpair(Aup[k, j] * EPSD[a, b] * es[b, j], s_[a, k]) for k in range(4) for j in range(4) for a in range(2) for b in range(2))
    F2 -= sum(gpair(Alo[k, j] * EPSD[bd, ad] * esb[ad, j], sb_[bd, k]) for k in range(4) for j in range(4) for ad in range(2) for bd in range(2))
    # component form of eta.lambda containing -A^{kj} eps_{ab} eta^b_j and -A_{kj} eps_{bd ad} etabar^{ad j}
    eq = np.array([[sum(-Aup[k, j] * EPSD[a, b] * es[b, j] for j in range(4) for b in range(2)) for k in range(4)] for a in range(2)])
    eqb = np.array([[sum(-Alo[k, j] * EPSD[bd, ad] * esb[ad, j] for j in range(4) for ad in range(2)) for k in range(4)] for bd in range(2)])
    record('S5', f'su({N}): -(rho(A_W)eta,lambda^*) = (F2) scalar terms; rho(A_W)eta components', [F2, join(eq, eqb)], [IF, rAe])
# closure rewriting, abelian, constant v: -v^n F_nm + v^I d_m A_I = -v^n d_n A_m + d_m(v^n A_n + v^I A_I)
v = rnd(10); dA = rnd(4, 10)        # dA[m, M] = d_m A_M
Fab = np.array([[dA[n, m] - dA[m, n] for m in range(4)] for n in range(4)])     # F_nm
lhs = np.array([-sum(v[n] * Fab[n, m] for n in range(4)) + sum(v[4 + I] * dA[m, 4 + I] for I in range(6)) for m in range(4)])
rhs = np.array([-sum(v[n] * dA[n, m] for n in range(4)) + (sum(v[n] * dA[m, n] for n in range(4)) + sum(v[4 + I] * dA[m, 4 + I] for I in range(6))) for m in range(4)])
record('S6', 'abelian closure: -v^n F_nm + v^I d_m A_I = -v^n d_n A_m + d_m(v^n A_n + v^I A_I)', lhs, rhs)

# ---------------------------------------------------------------- 5. negative controls: sign-flipped entries must fail
log('=== 5. negative controls ===')
NEG = []


def neg(cid, desc, lhs, rhs):
    lhs = flat(lhs); rhs = flat(rhs)
    scale = max(1.0, float(np.max(np.abs(lhs))), float(np.max(np.abs(rhs))))
    res = float(np.max(np.abs(lhs - rhs))) / scale
    ok = res > 1e-3
    RES.append(dict(id=cid, check='NEGATIVE: ' + desc, residual=res, passed=ok))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:6s} res={res:.2e}  NEGATIVE (must differ): {desc}")


rq = RbD[0] @ RD[1]
x = rnd(16); q, qb = split(x); yq, yqb = split(rq @ x)
neg('N1', 'rho(q(e1^e2)) on xi with +1/2(s s - s s)', yq, np.einsum('ab,bi->ai', 0.5 * (SIGE[0] @ SIGBE[1] - SIGE[1] @ SIGBE[0]), q))
neg('N2', 'canonical element with +delta on the dotted block', Cmat, np.eye(16))
neg('N3', 'Dirac column product = gamma5_S', vol[np.ix_([0, 4, 24, 28], [0, 4, 24, 28])], g5S)
T, nG, f = colour_setup(2)
Ac = rnd(6, nG); AIs = [to_mat(Ac[I], T, -1j * GC) for I in range(6)]
neg('N4', 'potential with +(g^2/4)', -0.125 * sum(gpair(comm(AIs[I], AIs[J]), comm(AIs[I], AIs[J])) for I in range(6) for J in range(6)) / HBAR,
    (GC ** 2 / 4) * np.einsum('abe,cde,Ia,Jb,Ic,Jd->', f, f, Ac, Ac, Ac, Ac))
neg('N5', 'propagator hbar*<A^a A^b> = 2 delta (missing 1/2 of the pairing)', HBAR * 2 * np.linalg.inv(np.array([[gpair(-1j * GC * T[a], -1j * GC * T[b]) for b in range(nG)] for a in range(nG)])), 2 * np.eye(nG))
neg('N6', 'S S gauge parameter with +x^2 on the W part', GxIF[4:], x2 * GsIF[4:])

# ---------------------------------------------------------------- output
npass = sum(r['passed'] for r in RES)
log(f'RESULT: {npass}/{len(RES)} passed; max residual of the positive checks '
    f'{max(r["residual"] for r in RES if not r["check"].startswith("NEGATIVE")):.2e}')
with open(os.path.join(HERE, 'verify_supp.log'), 'w') as fh:
    fh.write('\n'.join(LOG) + '\n')
with open(os.path.join(HERE, 'verify_supp.json'), 'w') as fh:
    json.dump(RES, fh, indent=1)

