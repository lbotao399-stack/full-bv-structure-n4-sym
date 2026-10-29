#!/usr/bin/env python3
"""Numerical checks of the two-component / su(N) / SO(6) formulas written in the physics-notation
supplements of Chapters 6, 7, 8 that are not verbatim entries of the dictionary or of Chapter 8.

Realisation: the two-component realisation D of step2/notation/check_dictionary.py (copied literally:
SIGE, SIGBE, EPSU, EPSD, Sigma^I = engine B Sig, Sigmabar^I = engine B Sigb, rho_S, rho_Ss, pair_D, G2C),
which check_dictionary.py has verified against engine A / classical model and engine B.
numpy only (+ read-only import of engine_b/realization.py for the 't Hooft blocks), < 100 MB, seconds.
Writes check_supp_ch06_08.log and check_supp_ch06_08.json next to this file.
"""
import os, sys, json, itertools
from fractions import Fraction as Fr
import numpy as np

np.seterr(all='ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'engine_b'))
import realization as RB  # read-only

rng = np.random.default_rng(20261029)
TOL = 1e-10
RESULTS, LOG = [], []


def log(*a):
    s = ' '.join(str(x) for x in a)
    LOG.append(s)
    print(s)


def _flat(x):
    if isinstance(x, (list, tuple)):
        parts = [_flat(y) for y in x]
        return np.concatenate(parts) if parts else np.zeros(0, complex)
    return np.ravel(np.asarray(x, dtype=complex))


def record(cid, desc, lhs, rhs=0.0, tol=TOL):
    lhs = _flat(lhs); rhs = _flat(rhs)
    if rhs.size == 1 and lhs.size > 1:
        rhs = np.full(lhs.shape, rhs[0])
    scale = max(1.0, float(np.max(np.abs(lhs))) if lhs.size else 1.0, float(np.max(np.abs(rhs))) if rhs.size else 1.0)
    res = float(np.max(np.abs(lhs - rhs))) / scale if lhs.size else 0.0
    ok = res < tol
    RESULTS.append(dict(id=cid, check=desc, residual=res, passed=bool(ok)))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:10s} res={res:.2e}  {desc}")
    return ok


def rnd(*shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


# ---------------------------------------------------------------- realisation D (literal copy)
I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EPSU = np.array([[0, 1], [-1, 0]], complex)     # eps^{12} = +1
EPSD = np.array([[0, -1], [1, 0]], complex)     # eps_{12} = -1
SIGE = [PAULI[0], PAULI[1], PAULI[2], 1j * I2]          # sigma^m_{alpha alphadot}, m = 1..4 -> 0..3
SIGBE = [-PAULI[0], -PAULI[1], -PAULI[2], 1j * I2]      # sigmabar^{m alphadot alpha}
SIGI = [np.array(S.tolist(), dtype=complex) for S in RB.Sig]      # Sigma^I_{ij}
SIGBI = [np.array(S.tolist(), dtype=complex) for S in RB.Sigb]    # Sigmabar^{I ij}
EPS4 = np.zeros((4,) * 4)
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])


def split(x):
    return x[:8].reshape((2, 4) + x.shape[1:]), x[8:].reshape((2, 4) + x.shape[1:])


def join(a, b):
    return np.concatenate([a.reshape((8,) + a.shape[2:]), b.reshape((8,) + b.shape[2:])])


def rho_S(M, x):
    q, qb = split(x)
    if M < 4:
        s = 1j * np.einsum('di...,da->ai...', qb, SIGBE[M])
        sb = -1j * np.einsum('da,ai...->di...', SIGBE[M], q)
    else:
        I = M - 4
        s = np.einsum('ab,ji,bi...->aj...', EPSU, SIGI[I], q)
        sb = -np.einsum('ji,de,ei...->dj...', SIGBI[I], EPSU, qb)
    return join(s, sb)


def rho_Ss(M, y):
    s, sb = split(y)
    if M < 4:
        q = -1j * np.einsum('ad,di...->ai...', SIGE[M], sb)
        qb = 1j * np.einsum('ai...,ad->di...', s, SIGE[M])
    else:
        I = M - 4
        q = np.einsum('kj,ab,bj...->ak...', SIGBI[I], EPSD, s)
        qb = -np.einsum('kj,ed,dj...->ek...', SIGI[I], EPSD, sb)
    return join(q, qb)


def pair_D(x, y):          # (Q,S) = q_a^i s^a_i - qb_{d i} sb^{d i},  Q in Sigma, S in Sigma^*
    q, qb = split(x); s, sb = split(y)
    return np.sum(q * s) - np.sum(qb * sb)


def G2C(Q1, Q2):           # dictionary (G1),(G2) in SO(6) components
    q1, qb1 = split(Q1); q2, qb2 = split(Q2)
    out = []
    for m in range(4):
        out.append(1j * (np.einsum('di,da,ai->', qb2, SIGBE[m], q1) + np.einsum('di,da,ai->', qb1, SIGBE[m], q2)))
    for I in range(6):
        out.append(np.einsum('ji,aj,ab,bi->', SIGI[I], q1, EPSU, q2) + np.einsum('ji,aj,ab,bi->', SIGBI[I], qb1, EPSU, qb2))
    return np.array(out)


def Gam_D(Q1, Q2):
    return np.array([pair_D(Q1, rho_S(M, Q2)) for M in range(10)])


def rho_v_S(v, x):   # rho(v) on Sigma, v in V_10 (10 components)
    return sum(v[M] * rho_S(M, x) for M in range(10))


def rho_v_Ss(v, y):
    return sum(v[M] * rho_Ss(M, y) for M in range(10))


def lo(w6):          # w_{ij} = w_I Sigma^I_{ij}
    return sum(w6[I] * SIGI[I] for I in range(6))


def up(w6):          # w^{ij} = -w_I Sigmabar^{I ij} = (1/2) eps^{ijkl} w_{kl}
    return -sum(w6[I] * SIGBI[I] for I in range(6))


E16 = np.eye(16)


def unitQ(alpha, i):      # unit element Q^i_alpha: xi_beta^j = delta_{beta alpha} delta^j_i
    q = np.zeros((2, 4), complex); q[alpha, i] = 1
    return join(q, np.zeros((2, 4), complex))


def unitQb(ad, i):        # unit element Qbar_{alphadot i}
    qb = np.zeros((2, 4), complex); qb[ad, i] = 1
    return join(np.zeros((2, 4), complex), qb)


def unitS(alpha, i):      # unit element S^alpha_i: eta^beta_j = delta delta
    s = np.zeros((2, 4), complex); s[alpha, i] = 1
    return join(s, np.zeros((2, 4), complex))


def unitSb(ad, i):        # unit element Sbar^{alphadot i}
    sb = np.zeros((2, 4), complex); sb[ad, i] = 1
    return join(np.zeros((2, 4), complex), sb)


# ===================================================================== A. su(N) normalisation
log('=== A. su(N) normalisation (Remark c6-konishi) ===')


def gellmann(N):
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


for N in (2, 3):
    T = gellmann(N); dim = len(T)
    g = 0.7 + 0.3 * N; hb = 1.3
    record(f'A0-N{N}', f'su({N}): tr(T^a T^b) = delta/2', [np.trace(T[a] @ T[b]) - 0.5 * (a == b) for a in range(dim) for b in range(dim)])
    pair = lambda X, Y: -4 * hb / g ** 2 * np.trace(X @ Y)
    xa, ya = rng.normal(size=dim), rng.normal(size=dim)
    X = sum(-1j * g * xa[a] * T[a] for a in range(dim)); Y = sum(-1j * g * ya[a] * T[a] for a in range(dim))
    record(f'A1-N{N}', '(X,Y) = 2 hbar X^a Y^a for X = -i g X^a T^a', pair(X, Y), 2 * hb * xa @ ya)
    L1 = sum(-1j * g / np.sqrt(2) * xa[a] * T[a] for a in range(dim)); L2 = sum(-1j * g / np.sqrt(2) * ya[a] * T[a] for a in range(dim))
    record(f'A2-N{N}', '(lambda,lambda\') = hbar lambda^a lambda\'^a for lambda = -(ig/sqrt2) lambda^a T^a', pair(L1, L2), hb * xa @ ya)
    # scalars: A_I = -i g A_I^a T^a, phi_I = A_I^a T^a
    Aa = rng.normal(size=(6, dim))
    AI = [sum(-1j * g * Aa[I, a] * T[a] for a in range(dim)) for I in range(6)]
    phi = [sum(Aa[I, a] * T[a] for a in range(dim)) for I in range(6)]
    record(f'A3-N{N}', '(A_I,A_J) = 4 hbar tr(phi_I phi_J)', [pair(AI[I], AI[J]) - 4 * hb * np.trace(phi[I] @ phi[J]) for I in range(6) for J in range(6)])
    # SU(4) form: A_{ij} = A_I Sigma^I_{ij} (matrices in colour), A^{ij} = -A_I Sigmabar^{I ij}
    Alo = [[sum(AI[I] * SIGI[I][i, j] for I in range(6)) for j in range(4)] for i in range(4)]
    Aup = [[sum(-AI[I] * SIGBI[I][i, j] for I in range(6)) for j in range(4)] for i in range(4)]
    plo = [[sum(phi[I] * SIGI[I][i, j] for I in range(6)) for j in range(4)] for i in range(4)]
    pup = [[0.5 * sum(EPS4[i, j, k, l] * plo[k][l] for k in range(4) for l in range(4)) for j in range(4)] for i in range(4)]
    kon_su4 = 0.25 * sum(pair(Aup[i][j], Alo[i][j]) for i in range(4) for j in range(4))
    record(f'A4-N{N}', 'Konishi: (1/4)(A^{ij},A_{ij}) = sum_I (A_I,A_I) = 4 hbar tr(phi_I phi_I) = hbar tr(phi^{ij} phi_{ij}),  phi^{ij}=(1/2)eps phi_kl',
           [kon_su4 - sum(pair(AI[I], AI[I]) for I in range(6)),
            kon_su4 - 4 * hb * sum(np.trace(phi[I] @ phi[I]) for I in range(6)),
            kon_su4 - hb * sum(np.trace(pup[i][j] @ plo[i][j]) for i in range(4) for j in range(4))])
    # 20' in SO(6) and SU(4) form
    tIJ = np.array([[np.trace(phi[I] @ phi[J]) for J in range(6)] for I in range(6)])
    t20 = tIJ - np.eye(6) * np.trace(tIJ) / 6
    T4 = np.array([[[[np.trace(plo[i][j] @ plo[k][l]) for l in range(4)] for k in range(4)] for j in range(4)] for i in range(4)])
    trup = sum(np.trace(pup[m][n] @ plo[m][n]) for m in range(4) for n in range(4))
    T20_su4 = T4 - EPS4 * trup / 12
    T20_so6 = np.einsum('Iij,Jkl,IJ->ijkl', np.array(SIGI), np.array(SIGI), t20)
    record(f'A5-N{N}', "20': tr(phi_ij phi_kl) - (1/12)eps_ijkl tr(phi^mn phi_mn) = Sigma^I_ij Sigma^J_kl [tr(phi_I phi_J) - (1/6)delta tr(phi phi)], eps-trace 0",
           [T20_su4 - T20_so6, np.einsum('ijkl,ijkl->', EPS4, T20_su4)])
    # unit: hbar tr_g(ad ad)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2)
    def ad(Z):
        M = np.zeros((dim, dim), complex)
        for b in range(dim):
            c = Z @ T[b] - T[b] @ Z
            for a in range(dim):
                M[a, b] = 2 * np.trace(T[a] @ c)       # component along T^a
        return M
    unit = hb * np.trace(ad(X) @ ad(Y)) / (16 * np.pi ** 2 * pair(X, Y))
    record(f'A6-N{N}', 'unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2)', unit, -g ** 2 * N / (32 * np.pi ** 2))
    # inverse pairing: hbar (-,-)^{-1} in basis t_a = -igT^a is (1/2) delta, in basis -ig/sqrt2 T^a it is delta
    G = np.array([[pair(-1j * g * T[a], -1j * g * T[b]) for b in range(dim)] for a in range(dim)])
    Gp = np.array([[pair(-1j * g / np.sqrt(2) * T[a], -1j * g / np.sqrt(2) * T[b]) for b in range(dim)] for a in range(dim)])
    record(f'A7-N{N}', 'hbar*2(-,-)^{-1} -> delta^{ab} (A components);  hbar*(-,-)^{-1} -> delta^{ab} (lambda components)',
           [hb * 2 * np.linalg.inv(G) - np.eye(dim), hb * np.linalg.inv(Gp) - np.eye(dim)])
    # quadratic log counterterm in su(N) components, relative to the kinetic terms (ch06 eq:c6-fieldren, quadratic part)
    Fa = rng.normal(size=dim); dAa = rng.normal(size=(6, dim)); lb = rng.normal(size=dim); dl = rng.normal(size=dim)
    F = sum(-1j * g * Fa[a] * T[a] for a in range(dim))
    record(f'A8-N{N}', '-(1/2)(F0,F0) = 4 hbar (-(1/4)F^a F^a);  -(1/4)(dA^ij,dA_ij) = -(dA_I,dA_I) = 4 hbar (-(1/2) dA^a_I dA^a_I);  8i(lb, d l) = 8 hbar i lb^a dl^a',
           [-0.5 * pair(F, F) - 4 * hb * (-0.25 * Fa @ Fa),
            -sum(pair(sum(-1j * g * dAa[I, a] * T[a] for a in range(dim)), sum(-1j * g * dAa[I, a] * T[a] for a in range(dim))) for I in range(6))
            - 4 * hb * (-0.5 * np.sum(dAa * dAa)),
            8j * pair(sum(-1j * g / np.sqrt(2) * lb[a] * T[a] for a in range(dim)), sum(-1j * g / np.sqrt(2) * dl[a] * T[a] for a in range(dim)))
            - 8 * hb * 1j * lb @ dl])
    # finite counterterm (27/16) unit (dA,dA) = -(27 g^2 N/(256 pi^2)) hbar (d.A^a)^2
    dA = rng.normal(size=dim); D = sum(-1j * g * dA[a] * T[a] for a in range(dim))
    record(f'A9-N{N}', '(27/16) unit int(d_mA_m,d_nA_n) = -(27 g^2 N/(256 pi^2)) hbar (d_m A^a_m)^2',
           27 / 16 * unit * pair(D, D), -27 * g ** 2 * N / (256 * np.pi ** 2) * hb * dA @ dA)

# iota_w A = (1/4) w^{ij} A_{ij} = w_I A_I ;  (u,w) = (1/4) u^{ij} w_{ij} = u_I w_I
w6, A6 = rng.normal(size=6), rng.normal(size=6)
record('A10', 'iota_w A = (1/4) w^{ij} A_{ij} = w_I A_I', 0.25 * np.einsum('ij,ij->', up(w6), lo(A6)), w6 @ A6)
# propagator of A_I from <A_ij A_kl> = 4(-,-)^{-1} eps_ijkl, A_I = (1/4) Sigmabar^{I ji} A_ij
record('A11', '<A_I A_J> = (1/16) Sigb^{I ji} Sigb^{J lk} 4 eps_ijkl (-,-)^{-1} = 2 delta_IJ (-,-)^{-1}',
       [0.25 * np.einsum('ji,lk,ijkl->', SIGBI[I], SIGBI[J], EPS4) - 2 * (I == J) for I in range(6) for J in range(6)])
record('A12', 'u_I = (1/4) Sigmabar^{I ji} u_ij (inverse used in A11)', [0.25 * np.einsum('ji,ij->', SIGBI[I], lo(w6)) - w6[I] for I in range(6)])

# ===================================================================== B. numbers (exact)
log('=== B. exact numbers ===')
unitF = Fr(-1, 32)        # unit in units of g^2 N / pi^2
gam = {'iota_w A': -4 * unitF, 'A|_V': -4 * unitF, 'lambda': -8 * unitF, 'c': -2 * unitF, 'Konishi': -24 * unitF, "20'": 0 * unitF}
record('B1', 'letters: -4,-8,-2 units = g^2N/(8pi^2), g^2N/(4pi^2), g^2N/(16pi^2); Konishi -24 units = 3g^2N/(4pi^2)',
       [float(gam['iota_w A'] - Fr(1, 8)), float(gam['A|_V'] - Fr(1, 8)), float(gam['lambda'] - Fr(1, 4)), float(gam['c'] - Fr(1, 16)),
        float(gam['Konishi'] - Fr(3, 4)), float(gam["20'"])])
# weights of the field-counting functional and the kinetic factors 4,4,8,-1 (dictionary CT-quadratic vs action)
wt = dict(A=2, lam=4, c=1, Astar=-2, lamstar=-4, cstar=-1)
kin = {'gauge': (Fr(-1, 2), Fr(-1, 8)), 'scalar': (Fr(-1, 4), Fr(-1, 16)), 'gaugino': (Fr(8), Fr(1)), 'ghost': (Fr(-1), Fr(1))}
record('B2', 'log counterterm / kinetic term = 4,4,8,-1 = weight sums 2+2, 2+2, 4+4, 1-2',
       [float(kin['gauge'][0] / kin['gauge'][1] - (wt['A'] + wt['A'])), float(kin['scalar'][0] / kin['scalar'][1] - (wt['A'] + wt['A'])),
        float(kin['gaugino'][0] / kin['gaugino'][1] - (wt['lam'] + wt['lam'])), float(kin['ghost'][0] / kin['ghost'][1] - (wt['c'] + wt['Astar']))])
record('B3', 'letter coefficient = -2 x weight (2,2,4,1 -> -4,-4,-8,-2)', [-2 * 2 + 4, -2 * 2 + 4, -2 * 4 + 8, -2 * 1 + 2])
# general column of ex:c6-bilinear: O_{w1w2} -> -4 unit (w1,w2) Trace, O = 4hbar w1I w2J tr phi_I phi_J, Trace = 4hbar tr phi_K phi_K
pairs = [(I, J) for I in range(6) for J in range(I, 6)]
Gm = np.zeros((21, 21))
for a, (I, J) in enumerate(pairs):
    if I == J:
        for b, (K, L) in enumerate(pairs):
            if K == L:
                Gm[b, a] = float(-4 * unitF)          # tr(phi_I phi_I) -> (1/8) sum_K tr(phi_K phi_K)
ev = np.sort(np.linalg.eigvals(Gm).real)
record('B4', "tr(phi_I phi_J) -> (g^2N/8pi^2) delta_IJ tr(phi_K phi_K): eigenvalues 0 (x20) and 6/8 = 3/4 (Konishi)",
       [ev[:20], ev[20] - 0.75])
record('B5', 'D_1 = 2 gamma: (1/2) D_1(Konishi) = -24 units;  two-point: -(1/4)(D^T C + C D) = -gamma C for D = 2 gamma',
       [float(Fr(-48, 2) - (-24)), (lambda gm, C: -0.25 * (2 * gm * C + C * 2 * gm) + gm * C)(0.37, 1.9)])
gm_, ell, r2 = 0.013, 1.7, 3.1
record('B6', '(x^2/l^2)^(-gamma) = 1 - gamma log(x^2/l^2) + O(gamma^2)', ((r2 / ell ** 2) ** (-gm_) - (1 - gm_ * np.log(r2 / ell ** 2))) / gm_ ** 2 < 1.0, True)
record('B7', 'dims: Sym^2 6 = 1+20, 4x4bar = 1+15, wedge^2 6 = 15', [21 - 1 - 20, 16 - 1 - 15, 15 - 15])

# ===================================================================== C. SO(6) / SU(4)
log('=== C. SO(6), SU(4) ===')
gam6 = [np.block([[np.zeros((4, 4)), SIGBI[I]], [SIGI[I], np.zeros((4, 4))]]) for I in range(6)]
record('C1', 'gamma^I = [[0, Sigb^I],[Sig^I, 0]] on W + W^*: {gamma^I, gamma^J} = 2 delta^IJ; gamma^I Hermitian',
       [gam6[I] @ gam6[J] + gam6[J] @ gam6[I] - 2 * (I == J) * np.eye(8) for I in range(6) for J in range(6)] +
       [gam6[I] - gam6[I].conj().T for I in range(6)])
R6 = rng.normal(size=(6, 6)); R6 = R6 - R6.T
r = 0.25 * sum(R6[I, J] * SIGBI[J] @ SIGI[I] for I in range(6) for J in range(6))
record('C2', 'r^i_j = (1/4) R_IJ (Sigb^J Sig^I)^i_j traceless and anti-Hermitian for real antisymmetric R', [np.trace(r), r + r.conj().T])

# ===================================================================== D. ten-dimensional gamma blocks
log('=== D. ten-dimensional gamma matrices ===')
for _ in range(2):
    X = rnd(16); Y = rnd(16)
    q, qb = split(X); s, sb = split(Y)
    out = []
    for I in range(6):
        s_ = np.einsum('ab,ji,bi->aj', EPSU, SIGI[I], q)              # (Gamma_I zeta)^alpha_j = eps^{ab} Sig^I_ji xi_b^i
        sb_ = -np.einsum('ab,ji,bi->aj', EPSU, SIGBI[I], qb)          # (Gamma_I zeta)^{adot j} = -eps^{adot bdot} Sigb^{I ji} xibar_{bdot i}
        out.append(join(s_, sb_) - rho_S(4 + I, X))
        q_ = np.einsum('kj,ab,bj->ak', SIGBI[I], EPSD, s)             # (Gamma_I chi)_alpha^k = Sigb^{I kj} eps_ab eta^b_j
        qb_ = -np.einsum('kj,ba,aj->bk', SIGI[I], EPSD, sb)           # (Gamma_I chi)_{bdot k} = -Sig^I_kj eps_{bdot adot} etabar^{adot j}
        out.append(join(q_, qb_) - rho_Ss(4 + I, Y))
    record('D1', 'Gamma_I blocks in Sigma/Sigmabar form == rho(e_I) of the realisation (R3),(R4)', out)

# ===================================================================== E. Killing spinor, brackets
log('=== E. conformal Killing spinor, odd brackets ===')


def killing(eta, x):       # rho(x) eta, dictionary (K1) with xi = 0
    return sum(x[m] * rho_Ss(m, eta) for m in range(4))


eta = rnd(16)
s, sb = split(eta)
dq = np.array([split(killing(eta, np.eye(4)[m]))[0] for m in range(4)])   # d_m eps_alpha^i
dqb = np.array([split(killing(eta, np.eye(4)[m]))[1] for m in range(4)])  # d_m epsbar_{adot i}
record('E1', 'd_m eps_alpha^i = -i sigma^m_{alpha adot} etabar^{adot i},  d_m epsbar_{adot i} = i eta^alpha_i sigma^m_{alpha adot}',
       [dq[m] - (-1j * SIGE[m] @ sb) for m in range(4)] + [dqb[m] - 1j * np.einsum('ai,ad->di', s, SIGE[m]) for m in range(4)])
dslash = sum(rho_S(m, killing(eta, np.eye(4)[m])) for m in range(4))
lhs_up = 1j * sum(np.einsum('di,da->ai', dqb[m], SIGBE[m]) for m in range(4))     # i d_m epsbar_{adot i} sbar^{m adot alpha}
lhs_dn = -1j * sum(np.einsum('da,ai->di', SIGBE[m], dq[m]) for m in range(4))     # -i sbar^{m adot alpha} d_m eps_alpha^i
record('E2', 'dslash(xi + rho(x)eta) = 4 eta;  components i d epsbar sbar = 4 eta^alpha_i, -i sbar d eps = 4 etabar^{adot i}',
       [dslash - 4 * eta, lhs_up - 4 * s, lhs_dn - 4 * sb])
# odd-odd brackets on unit elements: V-part = Gamma|_V
outQQb, outQQ, outQbQb = [], [], []
for al in range(2):
    for i in range(4):
        for ad in range(2):
            for j in range(4):
                G = G2C(unitQ(al, i), unitQb(ad, j))[:4]
                outQQb.append(G - np.array([1j * (i == j) * SIGBE[m][ad, al] for m in range(4)]))
                outQQb.append(Gam_D(unitQ(al, i), unitQb(ad, j))[:4] - G)
        for be in range(2):
            for j in range(4):
                outQQ.append(G2C(unitQ(al, i), unitQ(be, j))[:4])
                outQbQb.append(G2C(unitQb(al, i), unitQb(be, j))[:4])
record('E3', '[Q^i_alpha, Qbar_{adot j}] = i delta^i_j sbar^{m adot alpha} P_m;  [Q,Q]|_V = [Qbar,Qbar]|_V = 0', outQQb + outQQ + outQbQb)


def dil_part(xi, eta):     # coefficient of Delta = 2x in the V-part x -> Gamma(xi, rho(x)eta)|_V : (1/8) div
    return sum(G2C(xi, killing(eta, np.eye(4)[m]))[m] for m in range(4)) / 8


outD = []
for al in range(2):
    for i in range(4):
        outD.append(dil_part(unitQ(al, i), unitS(al, i)) - 0.5)            # [Q^i_alpha, S^alpha_i]: +1/2 Delta
        outD.append(dil_part(unitQb(al, i), unitSb(al, i)) + 0.5)          # [Qbar_{adot i}, Sbar^{adot i}]: -1/2 Delta
        outD.append(dil_part(unitQ(al, i), -unitS(al, i)) + 0.5)           # [Q, Q^dagger], Q^dagger = -S: -1/2 Delta = D
        outD.append(dil_part(unitQb(al, i), unitSb(al, i)) + 0.5)          # [Qbar, Qbar^dagger], Qbar^dagger = +Sbar: -1/2 Delta = D
        outD.append(dil_part(unitQ(al, i), unitS(al, i)) - 0.5 * pair_D(unitQ(al, i), unitS(al, i)))
for al in range(2):
    for i in range(4):
        for be in range(2):
            for j in range(4):
                outD.append(dil_part(unitQ(al, i), unitS(be, j)) - 0.5 * (al == be) * (i == j))      # (1/2) delta_alpha^beta delta^i_j
                outD.append(dil_part(unitQb(al, i), unitSb(be, j)) + 0.5 * (al == be) * (i == j))    # -(1/2) delta delta
                outD.append(dil_part(unitQ(al, i), unitSb(be, j)))                                   # opposite chirality: no Delta
                outD.append(dil_part(unitQb(al, i), unitS(be, j)))
record('E4', 'Delta-part of [xi,eta] = (1/2)(xi,eta)Delta on unit elements; [Q,Q^dag], [Qbar,Qbar^dag] (no sum) contain -(1/2)Delta = D', outD)
# pure spinor u(x)z and a rank-two element of S_+ (x) W
u, z = rnd(2), rnd(4)
xi_pure = join(np.einsum('a,i->ai', u, z), np.zeros((2, 4), complex))
xi_r2 = join(np.einsum('a,i->ai', rnd(2), rnd(4)) + np.einsum('a,i->ai', rnd(2), rnd(4)), np.zeros((2, 4), complex))
record('E5', 'xi = u(x)z: Gamma(xi,xi) = 0 in all ten components (pure spinor); generic xi in S_+(x)W: Gamma|_V = 0',
       [G2C(xi_pure, xi_pure), G2C(xi_r2, xi_r2)[:4]])
record('E5b', 'generic rank-two xi in S_+(x)W has Gamma(xi,xi)|_{W perp} != 0 (pure-spinor condition is stronger)',
       float(np.linalg.norm(G2C(xi_r2, xi_r2)[4:])) > 1e-3, True)
xi_Q4 = unitQ(1, 3)
record('E6', 'xi_alpha^i = delta_{alpha 2} delta^i_4 (single component Q^4_2) is pure', G2C(xi_Q4, xi_Q4))

# ===================================================================== F. structures of Corollary c7-Q, eta-term of S^(1)(eta)
log('=== F. two-component structures ===')
for _ in range(2):
    q = rnd(2, 4); qb = rnd(2, 4); lam = rnd(16); lstar = rnd(16)
    lq, lqb = split(lam); ls, lsb = split(lstar)
    xiQ = join(q, np.zeros((2, 4), complex)); xiQb = join(np.zeros((2, 4), complex), qb)
    w6 = rng.normal(size=6); v4 = rng.normal(size=4)
    wv = np.concatenate([np.zeros(4), w6]); vv = np.concatenate([v4, np.zeros(6)])
    # (chi, zeta) for chi in Sigma^*, zeta in Sigma := (zeta, chi) of (P1)
    out = [pair_D(lam, rho_v_S(wv, xiQ)) - np.einsum('ji,aj,ab,bi->', lo(w6), q, EPSU, lq),
           pair_D(lam, rho_v_S(vv, xiQ)) - 1j * sum(v4[m] * np.einsum('di,da,ai->', lqb, SIGBE[m], q) for m in range(4)),
           pair_D(xiQ, lstar) - np.einsum('ai,ai->', q, ls),
           pair_D(lam, rho_v_S(wv, xiQb)) - (-np.einsum('kl,ak,ab,bl->', up(w6), qb, EPSU, lqb)),
           pair_D(lam, rho_v_S(vv, xiQb)) - 1j * sum(v4[m] * np.einsum('di,da,ai->', qb, SIGBE[m], lq) for m in range(4)),
           pair_D(xiQb, lstar) - (-np.einsum('di,di->', qb, lsb)),
           pair_D(lam, rho_v_S(wv, xiQ)) - G2C(xiQ, lam)[4:] @ w6,            # = iota_w Gamma(xi, lambda)
           pair_D(lam, rho_v_S(vv, xiQb)) - G2C(xiQb, lam)[:4] @ v4]
    record('F1', "Cor c7-Q structures for Q^i_alpha and Qbar_{adot i}: (rho(w)xi,lambda), (rho(v)xi,lambda), (xi,lambda^*)", out)
    # eta-term: -(rho(w)eta, lambda^*) = (w^{kj} eps_ab eta^b_j) lambda^{*a}_k - (w_kj eps_{bdot adot} etabar^{adot j}) lambdabar^{*bdot k}
    et = rnd(16); s, sb = split(et)
    rw = rho_v_Ss(wv, et); rq, rqb = split(rw)
    f2 = np.einsum('kj,ab,bj,ak->', up(w6), EPSD, s, ls) - np.einsum('kj,ba,aj,bk->', lo(w6), EPSD, sb, lsb)
    record('F2', 'F2 extra term of S^(1)(eta) = -(rho(A_perp)eta, lambda^*); (rho(w)eta)_alpha^k = -w^{kj} eps_ab eta^b_j, (rho(w)eta)_{bdot k} = -w_kj eps etabar',
           [-pair_D(rw, lstar) - f2, rq - (-np.einsum('kj,ab,bj->ak', up(w6), EPSD, s)), rqb - (-np.einsum('kj,ba,aj->bk', lo(w6), EPSD, sb))])

# ===================================================================== G. conjugations
log('=== G. conjugations ===')


def cc(xi):                # Lemma c8-cc: Sigma -> Sigma^*: eta = conj(xi), etabar = -conj(xibar)
    q, qb = split(xi)
    return join(q.conj(), -qb.conj())


n4 = np.zeros(10); n4[3] = 1.0          # n = e_4
out = []
for al in range(2):
    for i in range(4):
        out.append(rho_v_Ss(n4, cc(unitQ(al, i))) - (-unitQb(al, i)))      # (Q^i_alpha)^{dagger_e4} = -Qbar_{adot=alpha, i}
        out.append(rho_v_Ss(n4, cc(unitQb(al, i))) - (-unitQ(al, i)))      # (Qbar_{adot i})^{dagger_e4} = -Q^i_{alpha=adot}
xi = rnd(16); q, qb = split(xi)
nv = rng.normal(size=4); nv /= np.linalg.norm(nv); nn = np.concatenate([nv, np.zeros(6)])
d = rho_v_Ss(nn, cc(xi)); dq_, dqb_ = split(d)
out.append(dq_ - 1j * sum(nv[m] * SIGE[m] @ qb.conj() for m in range(4)))
out.append(dqb_ - 1j * sum(nv[m] * np.einsum('ai,ad->di', q.conj(), SIGE[m]) for m in range(4)))
record('G1', 'xi^{dagger_n} = rho(n) conj(xi): two-component of Def c8-Theta; n = e_4: (Q^i_alpha)^dag = -Qbar_{alpha i}, (Qbar_{adot i})^dag = -Q^i_{adot}', out)
record('G2', '(xi, conj xi) = sum|xi|^2 + sum|xibar|^2 > 0', pair_D(xi, cc(xi)), np.sum(np.abs(xi) ** 2))

# ===================================================================== H. vector fields
log('=== H. conformal Killing fields, physics dilatation D = -Delta/2 ===')


def vf(kind, par, x):
    if kind == 'P':
        return par.copy(), np.zeros((4, 4))
    if kind == 'E':
        return x.copy(), np.eye(4)
    if kind == 'Delta':
        return 2 * x, 2 * np.eye(4)
    if kind == 'K':
        b = par
        v = np.dot(x, x) * b - 2 * np.dot(b, x) * x
        J = 2 * np.outer(x, b) - 2 * np.outer(b, x) - 2 * np.dot(b, x) * np.eye(4)   # J[m,n] = d_m v^n
        return v, J


def bracket(k1, p1, k2, p2, x):          # [u,v]^m = u^n d_n v^m - v^n d_n u^m
    u, Ju = vf(k1, p1, x); v, Jv = vf(k2, p2, x)
    return u @ Jv - v @ Ju


x = rng.normal(size=4); b = rng.normal(size=4); b2 = rng.normal(size=4)
negE = lambda k2, p2: -bracket('E', None, k2, p2, x)
record('H1', '[-E, P_b] = P_b,  [-E, K_b] = -K_b  (D = -E = -Delta/2: dimensions +1, -1)',
       [negE('P', b) - vf('P', b, x)[0], negE('K', b) + vf('K', b, x)[0]])
record('H2', '[K_{v1}, P_{v2}] = (v1,v2) Delta + (2(v1,x)v2 - 2(v2,x)v1)  (Remark c7-KDelta)',
       bracket('K', b, 'P', b2, x) - (np.dot(b, b2) * 2 * x + 2 * np.dot(b, x) * b2 - 2 * np.dot(b2, x) * b))
v, J = vf('K', b, x)
record('H3', 'K_b: div v = -8 b.x, X_mn = d_m v_n - d_n v_m = 4(x_m b_n - x_n b_m), d_m v_n + d_n v_m = (1/2) div v delta_mn',
       [np.trace(J) + 8 * np.dot(b, x), (J - J.T) - 4 * (np.outer(x, b) - np.outer(b, x)), (J + J.T) - 0.5 * np.trace(J) * np.eye(4)])
# scalar of weight w: v^n d_n + (w/4) div v  ==  b^n (x^2 d_n - 2 x_n x.d) - 2 w (b.x)
grad = rng.normal(size=4); wW = 1.0
record('H4', 'K_b on a scalar of weight w: v.d + (w/4)div v = b^n(x^2 d_n - 2x_n x.d) - 2w(b.x)',
       (v @ grad + wW / 4 * np.trace(J)) - (np.dot(x, x) * b @ grad - 2 * np.dot(b, x) * (x @ grad) - 2 * wW * np.dot(b, x)))
# odd: [Delta, eps] = 2x.d eps - (1/8) div(2x) eps = 2x.d eps - eps   ((B3) with X = 0, r = 0)
xi0 = rnd(16); eta0 = rnd(16)
eps_at = lambda y: xi0 + killing(eta0, y)
d_eps = np.array([killing(eta0, np.eye(4)[m]) for m in range(4)])
brD = 2 * sum(x[m] * d_eps[m] for m in range(4)) - eps_at(x)
record('H5', '[Delta, xi] = -xi, [Delta, rho(x)eta] = +rho(x)eta  =>  [D,Q] = Q/2, [D,S] = -S/2', brD - (-xi0 + killing(eta0, x)))

# ===================================================================== I. Q^GF squared, L_eps
log('=== I. Q^GF, L_eps ===')
p = rng.normal(size=4); a = rnd(4); xx = rng.normal(size=4)
ph = np.exp(1j * p @ xx)
A_slot = -2 * a * ph                      # A^*_m -> -2 A^*_m in the slot of A_m
c_slot = np.sum(1j * p * A_slot)          # A_m -> d_m A_m in the slot of c
record('I1', '(Q^GF)^2 A^*_m = -2 d_m A^*_m (slot c) = 2 d^* A^*  (plane wave)', c_slot, -2 * np.sum(1j * p * a) * ph)
epsv, p2 = 0.37, 2.3
h = 1e-6
fP = lambda e: np.exp(-e * p2) / p2
num = 2 * epsv * (fP(epsv + h) - fP(epsv - h)) / (2 * h)
sdds = (fP((1 + h) ** 2 * epsv) - fP((1 - h) ** 2 * epsv)) / (2 * h)
record('I2', 's d/ds|_{s=1} e^{-s^2 eps p^2}/p^2 = 2 eps d/d eps (e^{-eps p^2}/p^2) = -2 eps e^{-eps p^2}',
       [num + 2 * epsv * np.exp(-epsv * p2), sdds + 2 * epsv * np.exp(-epsv * p2)], tol=1e-7)

# ===================================================================== N. negative controls (must differ)
log('=== N. negative controls ===')


def neg(cid, desc, diff):
    res = float(np.max(np.abs(_flat(diff))))
    ok = res > 1e-3
    RESULTS.append(dict(id=cid, check='NEGATIVE: ' + desc, residual=0.0 if ok else 1.0, passed=bool(ok)))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:10s} |diff|={res:.2e}  NEGATIVE CONTROL: {desc}")


neg('N1', '(Q^i_alpha)^{dagger_e4} = +Qbar (wrong sign)', rho_v_Ss(n4, cc(unitQ(0, 0))) - unitQb(0, 0))
neg('N2', '[Q^i_alpha,Qbar_{adot j}] = -i delta sbar P (wrong sign)', G2C(unitQ(0, 1), unitQb(1, 1))[:4] + np.array([1j * SIGBE[m][1, 0] for m in range(4)]))
neg('N3', 'Delta-part of [Q,Q^dag] = +Delta/2 (wrong sign)', dil_part(unitQ(0, 0), -unitS(0, 0)) - 0.5)
neg('N4', '(rho(w)xibar,lambda) with +w^{kl} (wrong sign)',
    pair_D(lam, rho_v_S(wv, xiQb)) - np.einsum('kl,ak,ab,bl->', up(w6), qb, EPSU, lqb))
neg('N5', 'd_m eps = +i sigma etabar (wrong sign)', dq[0] - 1j * SIGE[0] @ sb)
neg('N6', 'Konishi = 2 hbar tr(phi phi) (wrong factor)', kon_su4 - 2 * hb * sum(np.trace(phi[I] @ phi[I]) for I in range(6)))
neg('N7', "20' with 1/6 eps (wrong trace subtraction)", T4 - EPS4 * trup / 6 - T20_so6)

# ===================================================================== summary
npass = sum(r['passed'] for r in RESULTS)
log(f'=== SUMMARY: {npass}/{len(RESULTS)} checks passed; max residual {max(r["residual"] for r in RESULTS):.2e} ===')
json.dump(dict(results=RESULTS), open(os.path.join(HERE, 'check_supp_ch06_08.json'), 'w'), indent=1)
open(os.path.join(HERE, 'check_supp_ch06_08.log'), 'w').write('\n'.join(LOG) + '\n')
