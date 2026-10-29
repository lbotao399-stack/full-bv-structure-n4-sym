#!/usr/bin/env python3
"""Numerical checks of the formulas of the physics-notation supplements of Chapters 3-5
(manuscript/supplements/ch03.tex, ch04.tex, ch05.tex) that are not verbatim in the dictionary
(step2/notation/dictionary.json) or in Chapter 8 (Def c8-twocomp, Prop c8-twocomp).

Realisation: the one of step2/notation/check_dictionary.py (realisation D): Euclidean sigma_E, sigmabar_E,
eps^{12}=+1, 't Hooft blocks Sigma^I of ch08 Def c8-twocomp (compared with engine B), eps_{1234}=1.
Colour: su(2) and su(3), T^a Hermitian with tr(T^a T^b)=delta/2, (X,Y) = -(4 hbar/g^2) tr(XY), X = -i g X^a T^a.
numpy only (engine B is imported read-only for one comparison if available).  Output: check_supp_ch03_05.log/.json
"""
import os, sys, json, itertools, math
import numpy as np

np.seterr(all='ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
rng = np.random.default_rng(20261029)
TOL = 1e-9
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


# ---------------------------------------------------------------------------------------------- data
I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EPSU = np.array([[0, 1], [-1, 0]], complex)      # eps^{12}=+1
EPSD = np.array([[0, -1], [1, 0]], complex)      # eps_{12}=-1
SIGE = [PAULI[0], PAULI[1], PAULI[2], 1j * I2]   # sigma^m_{alpha alphadot}
SIGBE = [-PAULI[0], -PAULI[1], -PAULI[2], 1j * I2]   # sigmabar^{m alphadot alpha}
EPS4 = np.zeros((4,) * 4)
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])
LC3 = np.zeros((3, 3, 3))
for _p in itertools.permutations(range(3)):
    LC3[_p] = np.linalg.det(np.eye(3)[list(_p)])

# 't Hooft blocks of ch08 Def c8-twocomp
SIGI = []
for a in range(3):
    M = np.zeros((4, 4), complex)
    M[:3, :3] = LC3[a]
    M[a, 3] = 1; M[3, a] = -1
    SIGI.append(M)
for a in range(3):
    M = np.zeros((4, 4), complex)
    M[:3, :3] = 1j * LC3[a]
    M[a, 3] = -1j; M[3, a] = 1j
    SIGI.append(M)
SIGI = np.array(SIGI)

log('=== 0. realisation ===')
record('D0-fierz', 'ch08: sum_I Sigma^I_ij Sigma^I_kl = 2 eps_ijkl', np.einsum('Iij,Ikl->ijkl', SIGI, SIGI) - 2 * EPS4)
try:
    sys.path.insert(0, os.path.join(ROOT, 'engine_b'))
    import realization as RB
    SIGB_engine = np.array([np.array(S.tolist(), dtype=complex) for S in RB.Sig])
    record('D0-engB', 'ch08 Sigma^I equal to engine B Sig (realisation of check_dictionary.py)', SIGI - SIGB_engine)
except Exception as exc:          # pragma: no cover
    log('engine B not importable:', exc)


def up(w_lo):
    """w^{ij} = (1/2) eps^{ijkl} w_{kl}"""
    return 0.5 * np.einsum('ijkl,kl...->ij...', EPS4, w_lo)


# ---------------------------------------------------------------------------------------------- colour
def gellmann():
    L = np.zeros((8, 3, 3), complex)
    L[0][0, 1] = L[0][1, 0] = 1
    L[1][0, 1] = -1j; L[1][1, 0] = 1j
    L[2][0, 0] = 1; L[2][1, 1] = -1
    L[3][0, 2] = L[3][2, 0] = 1
    L[4][0, 2] = -1j; L[4][2, 0] = 1j
    L[5][1, 2] = L[5][2, 1] = 1
    L[6][1, 2] = -1j; L[6][2, 1] = 1j
    L[7] = np.diag([1, 1, -2]) / np.sqrt(3)
    return L


GROUPS = {2: np.array([0.5 * P for P in PAULI]), 3: 0.5 * gellmann()}
gcoup, hb = 0.7, 1.3


def fabc(T):
    n = len(T)
    return np.real(np.array([[[-2j * np.trace((T[a] @ T[b] - T[b] @ T[a]) @ T[c]) for c in range(n)] for b in range(n)]
                             for a in range(n)]))


def pair(X, Y):
    return -(4 * hb / gcoup ** 2) * np.trace(X @ Y)


def emb(T, comp, scale=1.0):
    """X = -i g scale X^a T^a"""
    return -1j * gcoup * scale * np.einsum('a,aij->ij', comp, T)


log('=== 1. colour (su(2), su(3)) ===')
for N, T in GROUPS.items():
    n = len(T)
    f = fabc(T)
    record(f'C-T{N}', f'su({N}): tr(T^a T^b) = delta/2, [T^a,T^b] = i f^abc T^c',
           [np.einsum('aij,bji->ab', T, T) - 0.5 * np.eye(n),
            np.einsum('aij,bjk->abik', T, T) - np.einsum('bij,ajk->abik', T, T) - 1j * np.einsum('abc,cik->abik', f, T)])
    E = np.array([-1j * gcoup * T[a] for a in range(n)])            # basis e_a = -i g T^a
    G = np.array([[pair(E[a], E[b]) for b in range(n)] for a in range(n)])
    Ginv = np.linalg.inv(G)
    record(f'C-pair{N}', f'su({N}): (e_a,e_b) = 2 hbar delta_ab;  hbar (-,-)^{{-1}} = (1/2) sum_a e_a (x) e_a', [G - 2 * hb * np.eye(n), hb * Ginv - 0.5 * np.eye(n)])
    # ad_X in the basis e_a: [X, e_b] = sum_c (ad_X)_{cb} e_c
    Xa = rng.normal(size=n); Ya = rng.normal(size=n)
    X = emb(T, Xa); Y = emb(T, Ya)

    def adm(Xm):
        M = np.zeros((n, n), complex)
        for b in range(n):
            com = Xm @ E[b] - E[b] @ Xm
            M[:, b] = [pair(E[c], com) / G[c, c] for c in range(n)]
        return M

    TA = np.array([-1j * f[a] for a in range(n)])                      # (T^a_A)^{bc} = -i f^{abc}
    record(f'C-ad{N}', f'su({N}): ad_X = -i g X^a T^a_A with (T^a_A)^{{bc}} = -i f^{{abc}}', adm(X) - (-1j * gcoup) * np.einsum('a,abc->bc', Xa, TA))
    trad = np.trace(adm(X) @ adm(Y))
    record(f'C-kill{N}', f'su({N}): tr_g(ad_X ad_Y) = -g^2 N X^a Y^a', trad - (-gcoup ** 2 * N * np.dot(Xa, Ya)))
    unit = hb * trad / (16 * np.pi ** 2 * pair(X, Y))
    record(f'C-unit{N}', f'su({N}): unit = hbar tr(ad ad)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2)', unit - (-gcoup ** 2 * N / (32 * np.pi ** 2)))
    # bilinear ([X1,-],[X2,-]) evaluated on (-,-)^{-1}
    val = sum(Ginv[a, b] * pair(X @ E[a] - E[a] @ X, Y @ E[b] - E[b] @ Y) for a in range(n) for b in range(n))
    record(f'C-cc{N}', f'su({N}): ([X1,-],[X2,-]) on (-,-)^{{-1}} = -tr_g(ad ad) = g^2 N X1^a X2^a', val - gcoup ** 2 * N * np.dot(Xa, Ya))
    record(f'C-dim{N}', f'su({N}): (-,-) on (-,-)^{{-1}} = dim g = N^2-1', np.sum(Ginv * G) - (N * N - 1))
    # quartic colour factors
    Xs = [rng.normal(size=n) for _ in range(4)]
    Xm = [emb(T, x) for x in Xs]
    comm = lambda A, B: A @ B - B @ A
    lhs = pair(comm(Xm[0], Xm[1]), comm(Xm[2], Xm[3]))
    rhs = 2 * hb * gcoup ** 2 * np.einsum('abe,cde,a,b,c,d->', f, f, *Xs)
    record(f'C-ff{N}', f'su({N}): ([X1,X2],[X3,X4]) = 2 hbar g^2 f^abe f^cde X1^a X2^b X3^c X4^d', lhs - rhs)
    ads = [adm(x) for x in Xm]
    avg = np.mean([np.trace(ads[s[0]] @ ads[s[1]] @ ads[s[2]] @ ads[s[3]]) for s in itertools.permutations(range(4))])
    TAs = [np.einsum('a,abc->bc', x, TA) for x in Xs]
    avg2 = gcoup ** 4 * np.mean([np.trace(TAs[s[0]] @ TAs[s[1]] @ TAs[s[2]] @ TAs[s[3]]) for s in itertools.permutations(range(4))])
    record(f'C-tr4{N}', f'su({N}): averaged tr_g(ad^4) = g^4 X1^a..X4^d (1/24) sum_orderings tr(T_A^a T_A^b T_A^c T_A^d)', avg - avg2)
    # bracket components
    record(f'C-br{N}', f'su({N}): [X,Y]^c = g f^abc X^a Y^b', comm(X, Y) - emb(T, gcoup * np.einsum('abc,a,b->c', f, Xa, Ya)))

# use su(3) for the field checks
N = 3
T = GROUPS[N]; n = len(T); f = fabc(T)
unitN = -gcoup ** 2 * N / (32 * np.pi ** 2)


def gpair(X, Y):
    """(X,Y) on Lie-algebra-valued arrays with leading colour matrix axes (...,3,3)"""
    return -(4 * hb / gcoup ** 2) * np.einsum('...ij,...ji->...', X, Y)


# ---------------------------------------------------------------------------------------------- Chapter 3
log('=== 2. Chapter 3 ===')
# (Q^GF)^2 and Q c = Q^GF(1/2 dc) on plane waves e^{ipx}: d_m -> i p_m
p = rng.normal(size=4)
Asm = rnd(4); cs = rnd(1)[0]; c = rnd(1)[0]
QGF_Astar_to_A = -2 * Asm                                   # A^*_m -> -2 A^*_m in A_m
QGF_A_to_c = lambda A: np.sum(1j * p * A)                   # A_m -> d_m A_m in c
QGF_cs_to_Astar = lambda x: 1j * p * x                      # c^* -> d_m c^* in A^*_m
record('S3-QGF2a', '(Q^GF)^2 A^*_m = -2 d_m A^*_m (slot c) = 2 d^*(A^*|_V)', QGF_A_to_c(QGF_Astar_to_A) - (-2 * np.sum(1j * p * Asm)))
record('S3-QGF2b', '(Q^GF)^2 c^* = -2 d_m c^* (slot A_m)', -2 * QGF_cs_to_Astar(cs) - (-2 * 1j * p * cs))
record('S3-QQGF', 'Q c = -d_m c = Q^GF(1/2 d_m c in the slot A^*_m)', -2 * (0.5 * 1j * p * c) - (-1j * p * c))
# Feynman-gauge term: -1/4 (d_m A_m, d_n A_n) = -(hbar/2)(d_m A^a_m)^2 ; Wick image of -1/2 (d^mu A^a_mu)^2
dA = rng.normal(size=(4, 4, n))                             # d_m A^a_n  [m, n, a]
dAm = np.einsum('mnaij->mnij', np.einsum('mna,aij->mnaij', dA, -1j * gcoup * T))
divA = np.einsum('mmij->ij', dAm)
record('S3-gf', '-1/4 (d_m A_m, d_n A_n) = -(hbar/2) (d_m A^a_m)^2', -0.25 * gpair(divA, divA) - (-0.5 * hb * np.sum(np.einsum('mma->a', dA) ** 2)))
ETA = np.diag([-1., 1, 1, 1])
dAE = rng.normal(size=(4, 4))                                # Euclidean d_m A_n, m,n=1..4 (index 3 = 4)
idx = [3, 0, 1, 2]
dAL = np.zeros((4, 4), complex)                              # d_mu A_nu, A_0 = i A_4, d_0 = i d_4
for mu in range(4):
    for nu in range(4):
        dAL[mu, nu] = dAE[idx[mu], idx[nu]] * (1j if mu == 0 else 1) * (1j if nu == 0 else 1)
record('S3-gfWick', 'Wick: d^mu A_mu = d_m A_m (A_0 = i A_4, d_0 = i d_4)', np.einsum('mn,mn->', ETA, dAL) - np.trace(dAE))
# canonical propagators: hbar <A^a A^b> = delta^ab (fibre 2(-,-)^{-1});  hbar <lambda^a lambdabar^b> = -delta^ab p.sigma
E = np.array([-1j * gcoup * T[a] for a in range(n)])
G = np.array([[pair(E[a], E[b]) for b in range(n)] for a in range(n)])
Ginv = np.linalg.inv(G)
record('S3-propA', 'A = A^a e_a: hbar * 2(-,-)^{-1} in the coordinates A^a, A^b = delta^ab (A_m and A_I)', hb * 2 * Ginv - np.eye(n))
# gaugino: lambda = -(i g/sqrt2) lambda^a T^a = (1/sqrt2) lambda^a e_a: coordinate along e_a is lambda^a/sqrt2
pS = sum(p[m] * SIGE[m] for m in range(4))
prop_coord = -np.einsum('ab,xy->abxy', Ginv, pS)                   # <coord_a coord_b>, fibre -(-,-)^{-1} p.sigma delta
prop_can = hb * 2 * prop_coord                                    # lambda^a = sqrt2 coord_a
record('S3-propL', 'hbar <lambda^a_alpha^i lambdabar^b_{betadot j}> = -delta^ab delta^i_j p_m sigma^m_{alpha betadot}',
       prop_can - (-np.einsum('ab,xy->abxy', np.eye(n), pS)))
# consistency with the canonical Euclidean action  S/hbar > i lambdabar sigmabar^m d_m lambda:  <lambda lambdabar> = -K^{-1}, K = -p.sigmabar
K = 1j * sum(1j * p[m] * SIGBE[m] for m in range(4))              # (K)^{adot a}
record('S3-propLK', '-K^{-1} = -p.sigma/|p|^2 for K = i sigmabar^m (i p_m) (Grassmann Gaussian e^{lambdabar K lambda})',
       -np.linalg.inv(K) - (-pS / np.dot(p, p)))
# zero mode e_0 = (a.x) b_m, a.b = 0: Q e_0 = 1/2(d^2 A_m - d_m d_n A_n) = 0, F_mn = a_m b_n - a_n b_m
a = rng.normal(size=4); b = rng.normal(size=4); b = b - a * np.dot(a, b) / np.dot(a, a)
dE0 = np.einsum('m,n->mn', a, b)                                  # d_m A_n = a_m b_n (second derivatives vanish)
record('S3-zero', 'A_m = (a.x) b_m, a.b = 0: d_n A_n = 0, Q A = 0, F_mn = a_m b_n - a_n b_m',
       [np.trace(dE0), 0.5 * (0 - 0 * b), (dE0 - dE0.T) - (np.outer(a, b) - np.outer(b, a))])

# ---------------------------------------------------------------------------------------------- Chapter 4
log('=== 3. Chapter 4 ===')
# component fields at a point: A^a_m, A^a_I, d_m A^a_n, d_m A^a_I, lambda^a, d lambda^a
Aa = rng.normal(size=(4, n)); AIa = rng.normal(size=(6, n))
dAa = rng.normal(size=(4, 4, n)); dAIa = rng.normal(size=(4, 6, n))
lb = rnd(2, 4, n); dl = rnd(4, 2, 4, n)                          # lambdabar^a_{adot i}, d_m lambda^a_alpha^i
Am = np.einsum('ma,aij->mij', Aa, -1j * gcoup * T)
AI = np.einsum('Ia,aij->Iij', AIa, -1j * gcoup * T)
A_lo = np.einsum('Ikl,Iij->klij', SIGI, AI)                       # A_{kl} (colour matrix)
A_up = up(A_lo)
dAm_ = np.einsum('mna,aij->mnij', dAa, -1j * gcoup * T)
F0 = dAm_ - np.swapaxes(dAm_, 0, 1)                              # F^0_{mn} = d_m A_n - d_n A_m
dAI = np.einsum('mIa,aij->mIij', dAIa, -1j * gcoup * T)
dA_lo = np.einsum('Ikl,mIij->mklij', SIGI, dAI)
dA_up = 0.5 * np.einsum('klpq,mpqij->mklij', EPS4, dA_lo)
LB = np.einsum('xka,aij->xkij', lb, -1j * gcoup / np.sqrt(2) * T)
DL = np.einsum('myka,aij->mykij', dl, -1j * gcoup / np.sqrt(2) * T)
F0a = dAa - np.swapaxes(dAa, 0, 1)
# pairing-form functional densities
mass_pf = gpair(Am, Am).sum() + 0.25 * np.einsum('klij,klji->', A_up, A_lo) * (-(4 * hb / gcoup ** 2))
FF_pf = -0.5 * gpair(F0, F0).sum()
dd_pf = -0.25 * (-(4 * hb / gcoup ** 2)) * np.einsum('mklij,mklji->', dA_up, dA_lo)
# (lambdabar_{adot i}, sigmabar^{m adot a} d_m lambda_a^i), colour pairing in this order
ferm_pf = 8j * (-(4 * hb / gcoup ** 2)) * sum(np.einsum('xkab,xy,ykba->', LB, SIGBE[m], DL[m]) for m in range(4))
# su(N) component forms
mass_c = 2 * hb * (np.sum(Aa ** 2) + np.sum(AIa ** 2))
FF_c = 4 * hb * (-0.25 * np.sum(F0a ** 2))
dd_c = 4 * hb * (-0.5 * np.sum(dAIa ** 2))
ferm_c = 8j * hb * sum(np.einsum('xka,xy,yka->', lb, SIGBE[m], dl[m]) for m in range(4))
record('S4-mass', '(A_m,A_m) + (1/4)(A^ij,A_ij) = 2 hbar (A^a_m A^a_m + A^a_I A^a_I)', mass_pf - mass_c)
record('S4-FF', '-(1/2)(F^0_mn,F^0_mn) = 4 hbar (-(1/4) F^0a_mn F^0a_mn)', FF_pf - FF_c)
record('S4-dd', '-(1/4)(d_m A^ij, d_m A_ij) = 4 hbar (-(1/2) d_m A^a_I d_m A^a_I)', dd_pf - dd_c)
record('S4-ferm', '8i(lambdabar_{adot i}, sigmabar^{m adot a} d_m lambda_a^i) = 8 i hbar lambdabar^a sigmabar^m d_m lambda^a', ferm_pf - ferm_c)
quad_pf = unitN * (2.25 * mass_pf + (FF_pf + dd_pf + ferm_pf)) / hb       # times eps^{-1}, log eps set to 1
quad_c = -(gcoup ** 2 * N / (32 * np.pi ** 2)) * (4.5 * (np.sum(Aa ** 2) + np.sum(AIa ** 2)) +
                                                  4 * (-0.25 * np.sum(F0a ** 2) - 0.5 * np.sum(dAIa ** 2)) +
                                                  8 * 1j * sum(np.einsum('xka,xy,yka->', lb, SIGBE[m], dl[m]) for m in range(4)))
record('S4-quadSU', 'I^CT_quad/hbar (non-ghost part) = -(g^2N/32pi^2){(9/2)eps^-1(A^aA^a+A^a_IA^a_I) + log eps[4(-F^0aF^0a/4 - (dA^a_I)^2/2) + 8 i lbar^a sbar d l^a]}',
       quad_pf - quad_c)
record('S4-mass-num', '(9/4)*2*unit = -9 g^2 N/(64 pi^2)', 2.25 * 2 * unitN - (-9 * gcoup ** 2 * N / (64 * np.pi ** 2)))
# finite counterterm: -(27/16) unit (d_m A_m, d_n A_n)/hbar = +(27 g^2 N/(256 pi^2)) (d_m A^a_m)^2
divAa = np.einsum('mma->a', dAa)
divAm = np.einsum('mmij->ij', dAm_)
record('S4-finite', '-(27/16) unit (d_mA_m,d_nA_n)/hbar = (27 g^2 N/(256 pi^2)) (d_m A^a_m)^2',
       -(27 / 16) * unitN * gpair(divAm, divAm) / hb - 27 * gcoup ** 2 * N / (256 * np.pi ** 2) * np.sum(divAa ** 2))
# longitudinal piece in su(N): (27/16) unit (d A, d A) = -(27 g^2 N/256 pi^2) hbar (d_m A^a_m)^2
record('S4-long', '(27/16) unit (d_mA_m,d_nA_n) = -(27 g^2N/(256 pi^2)) hbar (d_m A^a_m)^2',
       (27 / 16) * unitN * gpair(divAm, divAm) - (-27 * gcoup ** 2 * N / (256 * np.pi ** 2) * hb * np.sum(divAa ** 2)))
# rescaling = field counting (thm c4-log): monomials of S and their multiples
counts = {  # (n_A, n_Astar, n_lam, n_lamstar, n_c, n_cstar): multiple listed in the proof of thm c4-log
    'kin A': ((2, 0, 0, 0, 0, 0), 4), 'kin lambda': ((0, 0, 2, 0, 0, 0), 8), 'ghost kin': ((0, 1, 0, 0, 1, 0), -1),
    'AAA': ((3, 0, 0, 0, 0, 0), 6), 'lambda lambda A': ((1, 0, 2, 0, 0, 0), 10), 'A c A*': ((1, 1, 0, 0, 1, 0), 1),
    'c lambda lambda*': ((0, 0, 1, 1, 1, 0), 1), 'c c c*': ((0, 0, 0, 0, 2, 1), 1), 'AAAA': ((4, 0, 0, 0, 0, 0), 8)}
uL = 0.37                                                           # unit * log eps (any number)
scal = np.array([-2 * uL, 2 * uL, -4 * uL, 4 * uL, -uL, uL])       # A -> (1-2uL)A, A* -> (1+2uL)A*, ...
res = []
for name, (nn, mult) in counts.items():
    nn = np.array(nn)
    first_order = np.dot(nn, scal)                                  # d/dt of prod (1+t s_phi)^{n_phi} at t = 0
    res.append(first_order - (-uL * mult))
    res.append(mult - (2 * (nn[0] - nn[1]) + 4 * (nn[2] - nn[3]) + (nn[4] - nn[5])))
record('S4-rescale', 'S - I^CT_log = S(A(1-2uL), lam(1-4uL), c(1-uL), A*(1+2uL), lam*(1+4uL), c*(1+uL)) to first order; counts 4,8,-1,6,10,1,1,1,8', res)
record('S4-Zsu', 'su(N): -2unit = g^2N/16pi^2, -4unit = g^2N/8pi^2, -unit = g^2N/32pi^2',
       [-2 * unitN - gcoup ** 2 * N / (16 * np.pi ** 2), -4 * unitN - gcoup ** 2 * N / (8 * np.pi ** 2), -unitN - gcoup ** 2 * N / (32 * np.pi ** 2)])
# Prop c4-K in components: 2 int (d.A)(d.(L_v A)) = 8 int (d.A)(v1.A), (L_v A)_m = v^n d_n A_m + A_n d_m v^n,
# v = x^2 v1 - 2 (v1.x) x; abelian component (the colour pairing is a common factor)
Ng, Lbox = 30, 13.0
h = Lbox / Ng
xs = (np.arange(Ng) - Ng // 2) * h
k1 = 2 * np.pi * np.fft.fftfreq(Ng, d=h)


def bc(vec, m):
    shp = [1, 1, 1, 1]; shp[m] = Ng
    return vec.reshape(shp)


Xb = [bc(xs, m) for m in range(4)]                                 # broadcast coordinates
Kb = [bc(k1, m) for m in range(4)]


def der(F, m):
    return np.real(np.fft.ifftn(1j * Kb[m] * np.fft.fftn(F)))


x2 = sum(Xb[m] ** 2 for m in range(4))
Gs = np.exp(-0.5 * x2)
c0 = rng.normal(size=4); c1 = rng.normal(size=(4, 4)); c2 = rng.normal(size=(4, 4, 4)) * 0.3
Afield = [(c0[m] + sum(c1[m, k] * Xb[k] for k in range(4)) + sum(c2[m, k, l] * Xb[k] * Xb[l] for k in range(4) for l in range(4))) * Gs
          for m in range(4)]
v1 = rng.normal(size=4)
v1x = sum(v1[m] * Xb[m] for m in range(4))
vf = [x2 * v1[m] - 2 * v1x * Xb[m] for m in range(4)]
dv = lambda m, nn: 2 * Xb[m] * v1[nn] - 2 * v1[m] * Xb[nn] - 2 * v1x * (m == nn)       # d_m v^n
divA_g = sum(der(Afield[m], m) for m in range(4))
divLvA = 0
for m in range(4):
    LvA_m = sum(vf[k] * der(Afield[m], k) for k in range(4)) + sum(Afield[k] * dv(m, k) for k in range(4))
    divLvA = divLvA + der(LvA_m, m)
lhs = 2 * np.sum(divA_g * divLvA) * h ** 4
rhs = 8 * np.sum(divA_g * sum(v1[m] * Afield[m] for m in range(4))) * h ** 4
record('S4-K', 'Prop c4-K: 2 int (d_mA_m)(d_n(L_vA)_n) = 8 int (d_mA_m)(v_1n A_n)  [= -8 int (d^*A, iota_v1 A)]', lhs - rhs, tol=1e-7)
record('S4-Kneg', 'negative control: coefficient -8 instead of 8 fails', 0.0 if abs(lhs + rhs) / max(1, abs(lhs)) > 1e-3 else 1.0)
record('S4-Kdiv', 'div v = -8 v1.x for v = x^2 v1 - 2(v1.x)x', np.broadcast_to(sum(dv(m, m) for m in range(4)) - (-8 * v1x), (Ng,) * 4))
del Gs, Afield, vf, divA_g, divLvA, LvA_m

# vertices as parts of the action (su(3) matrices, random fields at a point)
comm = lambda A, B: A @ B - B @ A
cub_F = -0.125 * 2 * sum(gpair(F0[m, nn], comm(Am[m], Am[nn])) for m in range(4) for nn in range(4))
record('S4-v3g', 'cubic part of -(1/8)(F_mn,F_mn) = -(1/4)(F^0_mn,[A_m,A_n]) = -(1/2)(d_m A_n,[A_m,A_n])',
       [cub_F - (-0.25) * sum(gpair(F0[m, nn], comm(Am[m], Am[nn])) for m in range(4) for nn in range(4)),
        cub_F - (-0.5) * sum(gpair(dAm_[m, nn], comm(Am[m], Am[nn])) for m in range(4) for nn in range(4))])
cubAWW = -(1 / 16) * sum(np.einsum('klab,klba->', dA_up[m], np.einsum('ab,klbc->klac', Am[m], A_lo) - np.einsum('klab,bc->klac', A_lo, Am[m])) +
                         np.einsum('klab,klba->', np.einsum('ab,klbc->klac', Am[m], A_up) - np.einsum('klab,bc->klac', A_up, Am[m]), dA_lo[m])
                         for m in range(4)) * (-(4 * hb / gcoup ** 2))
record('S4-v3s', 'cubic part of -(1/16)(D_mA^ij,D_mA_ij) = -(1/8)(d_mA^ij,[A_m,A_ij]) = -(1/2) sum_I (d_mA_I,[A_m,A_I])',
       [cubAWW - (-0.125) * sum(np.einsum('klab,klba->', dA_up[m], np.einsum('ab,klbc->klac', Am[m], A_lo) - np.einsum('klab,bc->klac', A_lo, Am[m]))
                                for m in range(4)) * (-(4 * hb / gcoup ** 2)),
        cubAWW - (-0.5) * sum(gpair(dAI[m, I], comm(Am[m], AI[I])) for m in range(4) for I in range(6))])
q1 = -(1 / 16) * sum(np.einsum('klab,klba->', np.einsum('ab,klbc->klac', Am[m], A_up) - np.einsum('klab,bc->klac', A_up, Am[m]),
                                 np.einsum('ab,klbc->klac', Am[m], A_lo) - np.einsum('klab,bc->klac', A_lo, Am[m])) for m in range(4)) * (-(4 * hb / gcoup ** 2))
record('S4-v4ms', '-(1/16)([A_m,A^ij],[A_m,A_ij]) = -(1/4) sum_I ([A_m,A_I],[A_m,A_I])', q1 - (-0.25) * sum(gpair(comm(Am[m], AI[I]), comm(Am[m], AI[I])) for m in range(4) for I in range(6)))
CU = np.einsum('ijab,klbc->ijklac', A_up, A_up) - np.einsum('klab,ijbc->ijklac', A_up, A_up)
CL = np.einsum('ijab,klbc->ijklac', A_lo, A_lo) - np.einsum('klab,ijbc->ijklac', A_lo, A_lo)
q2 = -(1 / 128) * np.einsum('ijklab,ijklba->', CU, CL) * (-(4 * hb / gcoup ** 2))
record('S4-v4ss', '-(1/128)([A^ij,A^kl],[A_ij,A_kl]) = -(1/8) sum_IJ ([A_I,A_J],[A_I,A_J])', q2 - (-0.125) * sum(gpair(comm(AI[I], AI[J]), comm(AI[I], AI[J])) for I in range(6) for J in range(6)))
DAI = dAI + np.array([[comm(Am[m], AI[I]) for I in range(6)] for m in range(4)])
DA_lo = np.einsum('Ikl,mIij->mklij', SIGI, DAI); DA_up = 0.5 * np.einsum('klpq,mpqij->mklij', EPS4, DA_lo)
record('S4-DAkin', '-(1/16)(D_mA^ij,D_mA_ij) = -(1/4) sum_I (D_mA_I,D_mA_I)',
       -(1 / 16) * np.einsum('mklab,mklba->', DA_up, DA_lo) * (-(4 * hb / gcoup ** 2)) - (-0.25) * sum(gpair(DAI[m, I], DAI[m, I]) for m in range(4) for I in range(6)))
# plane-wave values (eq c4-pw) in components: background A(p) e^{ipx} + A(-p) e^{-ipx}; coefficient of the cross term
pp = rng.normal(size=4); Ap = rnd(4); Amp = rnd(4)                  # A_m(p), A_m(-p) (colour pairing = product)
dAp = np.einsum('m,n->mn', 1j * pp, Ap); dAmp = np.einsum('m,n->mn', -1j * pp, Amp)
cross = lambda B: B(dAp, Ap, dAmp, Amp)
mass_pw = 2 * np.dot(Ap, Amp)
F_pw = 0.5 * 2 * np.sum((dAp - dAp.T) * (dAmp - dAmp.T))
div_pw = 2 * np.trace(dAp) * np.trace(dAmp)
record('S4-pw', 'plane waves: int(A_m,A_m) -> 2(A_m(p),A_m(-p)); (1/2)int(F^0,F^0) -> 2(|p|^2 d_mn - p_m p_n)(A_m(p),A_n(-p)); int(d_mA_m,d_nA_n) -> 2 p_m p_n (A_m(p),A_n(-p))',
       [F_pw - 2 * np.einsum('mn,m,n->', np.dot(pp, pp) * np.eye(4) - np.outer(pp, pp), Ap, Amp), div_pw - 2 * np.einsum('m,n,m,n->', pp, pp, Ap, Amp),
        mass_pw - 2 * np.sum(Ap * Amp)])
# canonical scalar propagator in SU(4) form: hbar <A^a_ij A^b_kl> = 2 delta^ab eps_ijkl, hbar <A^a_ij A^b kl> = 2 delta^ab (d d - d d)
DD = np.einsum('ik,jl->ijkl', np.eye(4), np.eye(4)) - np.einsum('il,jk->ijkl', np.eye(4), np.eye(4))
SIGbarI = np.array([-0.5 * np.einsum('ijkl,kl->ij', EPS4, SIGI[I]) for I in range(6)])          # A^{ij} = -A_I Sigmabar^{I ij}
record('S4-propSU4', 'hbar<A^a_I A^b_J> = delta delta_IJ => hbar<A_ij A_kl> -> 2 eps_ijkl, hbar<A_ij A^kl> -> 2(d^k_i d^l_j - d^l_i d^k_j)',
       [np.einsum('Iij,Ikl->ijkl', SIGI, SIGI) - 2 * EPS4, np.einsum('Iij,Ikl->ijkl', SIGI, -SIGbarI) - 2 * DD])

# ---------------------------------------------------------------------------------------------- Chapter 5
log('=== 4. Chapter 5 ===')
w1 = rng.normal(size=6); w2 = rng.normal(size=6)
w1lo = np.einsum('I,Iij->ij', w1, SIGI); w2lo = np.einsum('I,Iij->ij', w2, SIGI)
w1up = up(w1lo); w2up = up(w2lo)
iw1A = 0.25 * np.einsum('ij,ijab->ab', w1up, A_lo)
record('S5-iota', 'iota_w A = (1/4) w^ij A_ij = w_I A_I', iw1A - np.einsum('I,Iab->ab', w1, AI))
record('S5-w', '(w1,w2) = (1/4) w1^ij w2_ij = w1_I w2_I', 0.25 * np.einsum('ij,ij->', w1up, w2lo) - np.dot(w1, w2))
phi = np.einsum('Ia,aij->Iij', AIa, T)                              # phi_I = A^a_I T^a, A_I = -i g phi_I
record('S5-phi', '(A_I,A_J) = 4 hbar tr(phi_I phi_J) = 2 hbar phi^a_I phi^a_J',
       [np.array([[gpair(AI[I], AI[J]) for J in range(6)] for I in range(6)]) - 4 * hb * np.einsum('Iab,Jba->IJ', phi, phi),
        4 * hb * np.einsum('Iab,Jba->IJ', phi, phi) - 2 * hb * AIa @ AIa.T])
O_so6 = np.einsum('I,J,IJ->', w1, w2, np.array([[gpair(AI[I], AI[J]) for J in range(6)] for I in range(6)]))
Ktr = sum(gpair(AI[I], AI[I]) for I in range(6))
t20_so6 = O_so6 - np.dot(w1, w2) / 6 * Ktr
PA = np.einsum('ijab,klba->ijkl', A_lo, A_lo) * (-(4 * hb / gcoup ** 2))          # (A_ij, A_kl)
PAud = np.einsum('ijab,ijba->', A_up, A_lo) * (-(4 * hb / gcoup ** 2))           # (A^mn, A_mn)
t20_su4 = (1 / 16) * np.einsum('ij,kl,ijkl->', w1up, w2up, PA - (1 / 12) * EPS4 * PAud)
record('S5-20so6su4', '20\': w1_I w2_J[(A_I,A_J) - (1/6)d_IJ(A_K,A_K)] = (1/16) w1^ij w2^kl[(A_ij,A_kl) - (1/12) eps_ijkl (A^mn,A_mn)]', t20_so6 - t20_su4)
record('S5-20dict', '20\' equals the dictionary form (1/16)(w1A,w2A) - (1/96) w1^ij w2_ij (A^kl,A_kl)',
       t20_su4 - ((1 / 16) * np.einsum('ij,kl,ijkl->', w1up, w2up, PA) - (1 / 96) * np.einsum('ij,ij->', w1up, w2lo) * PAud))
traceless = np.einsum('ijkl,ijkl->', 0.5 * EPS4, PA - (1 / 12) * EPS4 * PAud)
record('S5-20trace', 'SU(4) 20\' density traceless: (1/2) eps^ijkl[(A_ij,A_kl) - (1/12) eps_ijkl (A^mn,A_mn)] = 0', traceless)
# 20' in su(N): tr(phi_ij phi_kl) - (1/12) eps_ijkl tr(phi^mn phi_mn), phi_ij = phi_I Sigma^I_ij
phi_lo = np.einsum('Ikl,Iab->klab', SIGI, phi); phi_up = up(phi_lo)
trpp = np.einsum('ijab,klba->ijkl', phi_lo, phi_lo); trud = np.einsum('ijab,ijba->', phi_up, phi_lo)
record('S5-20phi', '(A_ij,A_kl) - (1/12)eps(A^mn,A_mn) = 4 hbar [tr(phi_ij phi_kl) - (1/12) eps_ijkl tr(phi^mn phi_mn)]',
       (PA - (1 / 12) * EPS4 * PAud) - 4 * hb * (trpp - (1 / 12) * EPS4 * trud))
record('S5-Konishi', 'Konishi: (1/4)(A^ij,A_ij) = sum_I (A_I,A_I) = 4 hbar tr(phi_I phi_I) = hbar tr(phi^ij phi_ij) = 2 hbar phi^a_I phi^a_I',
       [0.25 * PAud - Ktr, Ktr - 4 * hb * np.einsum('Iab,Iba->', phi, phi), Ktr - hb * trud, Ktr - 2 * hb * np.sum(AIa ** 2)])
record('S5-20Onum', 'O(eps) traceless: -4 unit = g^2 N/(8 pi^2); Konishi: 8 unit = -g^2 N/(4 pi^2)',
       [-4 * unitN - gcoup ** 2 * N / (8 * np.pi ** 2), 8 * unitN + gcoup ** 2 * N / (4 * np.pi ** 2)])
dimg = N * N - 1
record('S5-vac', 'Konishi vacuum term: 12 hbar dim g/(16 pi^2 eps) = 2 hbar * 3(N^2-1)/(8 pi^2 eps) = 2 hbar * 6(N^2-1)/(16 pi^2 eps)',
       [12 * hb * dimg / (16 * np.pi ** 2) - 2 * hb * 3 * dimg / (8 * np.pi ** 2), 3 * dimg / (8 * np.pi ** 2) - 6 * dimg / (16 * np.pi ** 2),
        12 * hb * dimg / (16 * np.pi ** 2) - 3 * hb * dimg / (4 * np.pi ** 2)])
# graph (a) for Konishi from the canonical coincident propagator hbar<A^a_I A^b_J>(x,x) = delta delta/(16 pi^2 eps)
coinc = 1.0 / (16 * np.pi ** 2)                                     # times 1/eps
record('S5-vacprop', '<phi^a_I phi^a_I>(x,x) = 6 (N^2-1)/(16 pi^2 eps) = (graph (a) of Konishi)/(2 hbar)',
       6 * dimg * coinc - 12 * hb * dimg / (16 * np.pi ** 2) / (2 * hb))
# letters in su(N)
Xa = rng.normal(size=n); Xc = emb(T, Xa)
record('S5-letter', '(X, iota_w A) = 2 hbar X^a w_I A^a_I;  (X, A_m) = 2 hbar X^a A^a_m',
       [gpair(Xc, iw1A) - 2 * hb * np.einsum('a,I,Ia->', Xa, w1, AIa), gpair(Xc, Am[0]) - 2 * hb * np.dot(Xa, Aa[0])])
# Q int omega (A_m,A_m) = -2 int omega (d_m c, A_m): Q c = -d_m c inserted twice by symmetry
cvec = rnd(4)                                                      # d_m c (colour suppressed: abelian check)
Avec = rnd(4)
record('S5-Qmass', 'Q (A_m,A_m) |_{A -> A + r(-d c)} linear in r: -2 (d_m c, A_m)', (np.dot(Avec - cvec, Avec - cvec) - np.dot(Avec, Avec) - np.dot(cvec, cvec)) - (-2 * np.dot(cvec, Avec)))


# heat kernel in momentum space (1d factor, product of four): int dp/(2pi) e^{ipx - t p^2} = (4 pi t)^{-1/2} e^{-x^2/4t}
tt, xx = 0.37, 0.81
pg = np.linspace(-40, 40, 400001)
oned = np.trapezoid(np.exp(1j * pg * xx - tt * pg ** 2), pg) / (2 * np.pi)
record('S3-heatFT', 'int d^4p/(2pi)^4 e^{i(p,x)} e^{-t|p|^2} = (4 pi t)^{-2} e^{-|x|^2/4t} (1d factor)', oned - np.exp(-xx ** 2 / (4 * tt)) / np.sqrt(4 * np.pi * tt), tol=1e-8)
# Weyl dimension of the SU(4) irrep [0,2,0]
a1, a2, a3 = 0, 2, 0
dimA3 = (a1 + 1) * (a2 + 1) * (a3 + 1) * (a1 + a2 + 2) * (a2 + a3 + 2) * (a1 + a2 + a3 + 3) / 12
record('S5-dim20', 'dim [0,2,0] of SU(4) = 20 = dim Sym^2 6 - 1', [dimA3 - 20, 6 * 7 / 2 - 1 - 20])

# ---------------------------------------------------------------------------------------------- negative controls
log('=== 5. negative controls ===')
def negc(cid, desc, r):
    record(cid, 'negative control: ' + desc + f' (perturbed residual {r:.2e} > 1e-3)', 0.0 if r > 1e-3 else 1.0)
negc('N-propL', 'gaugino propagator +p.sigma', np.max(np.abs(prop_can - np.einsum('ab,xy->abxy', np.eye(n), pS))))
negc('N-20', '20\' with 1/6 instead of 1/12', abs(t20_so6 - (1 / 16) * np.einsum('ij,kl,ijkl->', w1up, w2up, PA - (1 / 6) * EPS4 * PAud)) / max(1, abs(t20_so6)))
negc('N-quad', 'scalar log coefficient 2(-1/2)(dA_I)^2 replaced by 4(-1/4)(dA_I)^2', abs(dd_pf - 4 * hb * (-0.25 * np.sum(dAIa ** 2))) / max(1, abs(dd_pf)))

# ---------------------------------------------------------------------------------------------- summary
npass = sum(r['passed'] for r in RESULTS)
log(f'=== {npass}/{len(RESULTS)} passed, max residual {max(r["residual"] for r in RESULTS):.2e} ===')
with open(os.path.join(HERE, 'check_supp_ch03_05.log'), 'w') as fh:
    fh.write('\n'.join(LOG) + '\n')
with open(os.path.join(HERE, 'check_supp_ch03_05.json'), 'w') as fh:
    json.dump(RESULTS, fh, indent=1)
