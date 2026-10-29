#!/usr/bin/env python3
"""check_dictionary.py -- P4 checks of the dictionary  ESW index-free notation -> Weinberg-Srednicki two-component
notation  (step2/notation/esw_to_weinberg_srednicki.tex, dictionary.json).

Every translated formula of the dictionary is implemented here twice:
  (IF) index-free, evaluated with the explicit realisation used by the Step-1 engines
       (classical/n4_model.py == engine_a/cliff.py, 16x16 blocks sig[M], sigbar[M], pairing s^T s*),
  (2C) two-component, implemented literally as typeset (sigma_E, sigmabar_E, epsilon, Sigma^I, Sigmabar^I),
and compared on random data: pointwise jets (value, first and second derivatives) for pointwise identities,
plane waves with total momentum 0 for integrated identities, Grassmann fields by polarisation
(theta_1 theta_2 coefficient = B(u1,u2) - B(u2,u1)).  The two-component realisation D itself is compared
with engine_b/realization.py (exact basis change) and with classical/engine A (intertwiner, computed).

Light: numpy only, < 200 MB, < 1 min.  Output: check_dictionary.log, check_results.json, dictionary.json.
"""
import os, sys, json, itertools, math
import numpy as np

np.seterr(all='ignore')          # Apple Accelerate emits spurious matmul warnings; every check is an explicit residual
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
for sub in ('engine_b', 'classical', 'engine_a'):
    sys.path.insert(0, os.path.join(ROOT, sub))
import realization as RB         # engine B (read-only)
import n4_model as NM            # classical realisation (read-only)
import cliff as CA               # engine A (read-only)

rng = np.random.default_rng(20261028)
TOL = 1e-9
RESULTS = []
LOG = []


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
    lhs = _flat(lhs)
    rhs = _flat(rhs)
    scale = max(1.0, float(np.max(np.abs(lhs))) if lhs.size else 1.0)
    res = float(np.max(np.abs(lhs - rhs))) / scale if lhs.size else 0.0
    ok = res < tol
    RESULTS.append(dict(id=cid, check=desc, residual=res, passed=bool(ok)))
    log(f"[{'PASS' if ok else 'FAIL'}] {cid:8s} res={res:.2e}  {desc}")
    return ok


def rnd(*shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


# ======================================================================================================
# 0. Srednicki two-component data [STANDARD S-SPN-2, S-SPN-3, S-SPN-5] and the Wick map
# ======================================================================================================
I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EPSU = np.array([[0, 1], [-1, 0]], complex)     # eps^{12} = +1      (35.1)
EPSD = np.array([[0, -1], [1, 0]], complex)     # eps_{12} = -1      (35.1)
SIG_L = [I2] + PAULI                            # sigma^mu_{a adot} = (1, sigma)          (34.28)/(35.2)
SIGB_L = [I2] + [-p for p in PAULI]             # sigmabar^{mu adot a} = (1, -sigma)      (35.3)
ETA = np.diag([-1., 1, 1, 1])                   # g_{mu nu} = diag(-1,1,1,1)             (2.4)
# Euclidean sigma matrices: sigma^mu v_mu (Lorentzian) = sigma_E^m v_m (Euclidean) with v^0 = -i v^4, v^k = v^k
SIGE = [PAULI[0], PAULI[1], PAULI[2], 1j * I2]          # sigma_E^m_{alpha alphadot},     m = 1..4 -> 0..3
SIGBE = [-PAULI[0], -PAULI[1], -PAULI[2], 1j * I2]      # sigmabar_E^{m alphadot alpha}
# the six 't Hooft blocks of engine B: Sigma^I_{ij} : W -> W^*, Sigmabar^{I ij} : W^* -> W
SIGI = [np.array(S.tolist(), dtype=complex) for S in RB.Sig]
SIGBI = [np.array(S.tolist(), dtype=complex) for S in RB.Sigb]
EPS4 = np.zeros((4,) * 4)
for _p in itertools.permutations(range(4)):
    EPS4[_p] = np.linalg.det(np.eye(4)[list(_p)])
DD = np.einsum('ik,jl->ijkl', np.eye(4), np.eye(4)) - np.einsum('il,jk->ijkl', np.eye(4), np.eye(4))

log('=== 0. Srednicki data and Wick map ===')
record('W1', 'sigmabar^{mu adot a} = eps^{ab} eps^{adot bdot} sigma^mu_{b bdot}  [S-SPN-5 (35.3)]',
       [SIGB_L[m] - EPSU @ SIG_L[m].T @ EPSU.T for m in range(4)])
for _ in range(2):
    vE = rng.normal(size=4)
    vL_up = np.array([-1j * vE[3], vE[0], vE[1], vE[2]])        # v^0 = -i v^4
    vL_dn = ETA @ vL_up
    record('W2', 'Wick: sigma^mu v_mu = sigma_E^m v_m,  sigmabar^mu v_mu = sigmabar_E^m v_m  (v^0 = -i v^4)',
           [sum(SIG_L[m] * vL_dn[m] for m in range(4)) - sum(SIGE[m] * vE[m] for m in range(4)),
            sum(SIGB_L[m] * vL_dn[m] for m in range(4)) - sum(SIGBE[m] * vE[m] for m in range(4))])
record('W3', 'Euclidean Clifford sigma_E^m sigmabar_E^n + sigma_E^n sigmabar_E^m = -2 delta^{mn}  (Wick image of S-MET-5)',
       [SIGE[m] @ SIGBE[n] + SIGE[n] @ SIGBE[m] + 2 * (m == n) * I2 for m in range(4) for n in range(4)] +
       [SIGBE[m] @ SIGE[n] + SIGBE[n] @ SIGE[m] + 2 * (m == n) * I2 for m in range(4) for n in range(4)])
record('W4', 'sigmabar_E^m = eps sigma_E^{mT} eps^T  and  sigmabar_E^m = -(sigma_E^m)^dagger',
       [SIGBE[m] - EPSU @ SIGE[m].T @ EPSU.T for m in range(4)] + [SIGBE[m] + SIGE[m].conj().T for m in range(4)])
record('W5', "Sigma^I antisymmetric, Sigmabar^I = (Sigma^I)^dagger, (1/2)eps^{ijkl}Sigma^I_{kl} = -Sigmabar^{I ij}",
       [SIGI[I] + SIGI[I].T for I in range(6)] + [SIGBI[I] - SIGI[I].conj().T for I in range(6)] +
       [0.5 * np.einsum('ijkl,kl->ij', EPS4, SIGI[I]) + SIGBI[I] for I in range(6)])
record('W6', 'Sigmabar^I Sigma^J + Sigmabar^J Sigma^I = 2 delta^{IJ},  tr(Sigmabar^I Sigma^J) = 4 delta^{IJ}',
       [SIGBI[I] @ SIGI[J] + SIGBI[J] @ SIGI[I] - 2 * (I == J) * np.eye(4) for I in range(6) for J in range(6)] +
       [np.trace(SIGBI[I] @ SIGI[J]) - 4 * (I == J) for I in range(6) for J in range(6)])
record('W7', "Fierz: sum_I Sig^I_ij Sig^I_kl = 2 eps_ijkl, sum_I Sigb^{I ij} Sigb^{I kl} = 2 eps^ijkl, "
       "sum_I Sigb^{I ij} Sig^I_kl = -2(d^i_k d^j_l - d^i_l d^j_k)",
       [sum(np.einsum('ij,kl->ijkl', SIGI[I], SIGI[I]) for I in range(6)) - 2 * EPS4,
        sum(np.einsum('ij,kl->ijkl', SIGBI[I], SIGBI[I]) for I in range(6)) - 2 * EPS4,
        sum(np.einsum('ij,kl->ijkl', SIGBI[I], SIGI[I]) for I in range(6)) + 2 * DD])


# ======================================================================================================
# 1. The two-component realisation D (dictionary Theorem 1)
#    Sigma  vector = (q_alpha^i [2x4], qb_{alphadot i} [2x4]);   Sigma^* vector = (s^alpha_i [2x4], sb^{alphadot i} [2x4])
#    V_10 index M = 0..3 <-> m = 1..4 (m=4: sigma_E^4 = i), M = 4..9 <-> I = 1..6
# ======================================================================================================
def split(x):
    return x[:8].reshape((2, 4) + x.shape[1:]), x[8:].reshape((2, 4) + x.shape[1:])


def join(a, b):
    return np.concatenate([a.reshape((8,) + a.shape[2:]), b.reshape((8,) + b.shape[2:])])


def rho_S(M, x):
    """rho(e_M) : Sigma -> Sigma^*   (dictionary formulas (R1),(R3))"""
    q, qb = split(x)
    if M < 4:
        s = 1j * np.einsum('di...,da->ai...', qb, SIGBE[M])          # s^a_i  = +i qb_{d i} sbar^{m d a}
        sb = -1j * np.einsum('da,ai...->di...', SIGBE[M], q)          # sb^{d i} = -i sbar^{m d a} q_a^i
    else:
        I = M - 4
        s = np.einsum('ab,ji,bi...->aj...', EPSU, SIGI[I], q)         # s^a_j = eps^{ab} Sig^I_{ji} q_b^i
        sb = -np.einsum('ji,de,ei...->dj...', SIGBI[I], EPSU, qb)     # sb^{d j} = -Sigb^{I ji} eps^{de} qb_{e i}
    return join(s, sb)


def rho_Ss(M, y):
    """rho(e_M) : Sigma^* -> Sigma   (dictionary formulas (R2),(R4))"""
    s, sb = split(y)
    if M < 4:
        q = -1j * np.einsum('ad,di...->ai...', SIGE[M], sb)           # q_a^i = -i sig^m_{a d} sb^{d i}
        qb = 1j * np.einsum('ai...,ad->di...', s, SIGE[M])            # qb_{d i} = +i s^a_i sig^m_{a d}
    else:
        I = M - 4
        q = np.einsum('kj,ab,bj...->ak...', SIGBI[I], EPSD, s)        # q_a^k = Sigb^{I kj} eps_{ab} s^b_j
        qb = -np.einsum('kj,ed,dj...->ek...', SIGI[I], EPSD, sb)      # qb_{e k} = -Sig^I_{kj} eps_{ed} sb^{d j}
    return join(q, qb)


def pair_D(x, y):
    """(Q,S) = q_a^i s^a_i - qb_{d i} sb^{d i}   (dictionary (P1)); trailing axes summed."""
    q, qb = split(x)
    s, sb = split(y)
    return np.sum(q * s) - np.sum(qb * sb)


E16 = np.eye(16)
RD = [np.array([rho_S(M, E16[:, k]) for k in range(16)]).T for M in range(10)]
RbD = [np.array([rho_Ss(M, E16[:, k]) for k in range(16)]).T for M in range(10)]
PD = np.array([[pair_D(E16[:, a], E16[:, b]) for b in range(16)] for a in range(16)])

log('=== 1. realisation D satisfies the manuscript identities ===')
record('R1', 'Clifford rho(v)rho(w)+rho(w)rho(v) = 2(v,w) on Sigma and Sigma^*  [ch01 Def c1-spaces]',
       [RbD[M] @ RD[N] + RbD[N] @ RD[M] - 2 * (M == N) * E16 for M in range(10) for N in range(10)] +
       [RD[M] @ RbD[N] + RD[N] @ RbD[M] - 2 * (M == N) * E16 for M in range(10) for N in range(10)])
record('R2', '(rho(v)Q1,Q2) = (Q1,rho(v)Q2): P_D rho symmetric  [ch01 Def c1-spaces, ESW 2.1]',
       [PD @ RD[M] - (PD @ RD[M]).T for M in range(10)] + [RbD[M].T @ PD - (RbD[M].T @ PD).T for M in range(10)])


def Gam_D(Q1, Q2):
    """Gamma(Q1,Q2)^M = (Q1, rho(e_M) Q2), Q1,Q2 in Sigma [ESW (28)]"""
    return np.array([pair_D(Q1, rho_S(M, Q2)) for M in range(10)])


Qs = [rnd(16) for _ in range(3)]
record('R3', 'Gamma symmetric; ESW Assumption 3.1 rho(G(Q1,Q2))Q3 + cyclic = 0',
       [Gam_D(Qs[0], Qs[1]) - Gam_D(Qs[1], Qs[0])] +
       [sum(sum(Gam_D(Qs[a], Qs[b])[M] * RD[M] for M in range(10)) @ Qs[c] for a, b, c in ((0, 1, 2), (1, 2, 0), (2, 0, 1)))])
vol4 = RbD[0] @ RD[1] @ RbD[2] @ RD[3]
proj_plus = np.zeros(16); proj_plus[:8] = 1
record('R4', 'chirality: rho(e1 e2 e3 e4) = +1 on S_+ (x) W, -1 on S_- (x) W^*  [Prop c1-chirality, orientation e1..e4]',
       vol4 - np.diag(2 * proj_plus - 1))
record('R5', 'rho(e_m) exchanges, rho(e_I) preserves 4d chirality (block structure of D)',
       [RD[M][:8, :8] for M in range(4)] + [RD[M][8:, 8:] for M in range(4)] + [RD[M][:8, 8:] for M in range(4, 10)] +
       [RD[M][8:, :8] for M in range(4, 10)])
# Euclidean Dirac blocks = Weinberg gamma (Wick) = -i Srednicki gamma (Wick)  [Weinberg db D1]
for m in range(4):
    gE = np.block([[np.zeros((2, 2)), -1j * SIGE[m]], [-1j * SIGBE[m], np.zeros((2, 2))]])
gS = [np.block([[np.zeros((2, 2)), SIG_L[mu]], [SIGB_L[mu], np.zeros((2, 2))]]) for mu in range(4)]   # (36.7)
vE = rng.normal(size=4)
vL_dn = ETA @ np.array([-1j * vE[3], vE[0], vE[1], vE[2]])
gE_v = np.block([[np.zeros((2, 2)), -1j * sum(SIGE[m] * vE[m] for m in range(4))],
                 [-1j * sum(SIGBE[m] * vE[m] for m in range(4)), np.zeros((2, 2))]])
gW_v = sum(-1j * gS[mu] * vL_dn[mu] for mu in range(4))            # gamma_W^mu = -i gamma_S^mu  (Weinberg db D1)
record('R6', 'Dirac column (psi_alpha, chibar^alphadot): rho(v) = -i gamma_S^mu v_mu = gamma_W^mu v_mu (Wick)', gE_v - gW_v)
g5S = np.diag([-1, -1, 1, 1])                                       # (36.44)
gE = [np.block([[np.zeros((2, 2)), -1j * SIGE[m]], [-1j * SIGBE[m], np.zeros((2, 2))]]) for m in range(4)]
record('R7', 'gamma_E^1 gamma_E^2 gamma_E^3 gamma_E^4 = -gamma5_S = +1 on undotted  [S-SPN-7]',
       gE[0] @ gE[1] @ gE[2] @ gE[3] + g5S)

# ------------------------------------------------------------------------------------------------------
log('=== 2. D == engine B (exact basis change) ===')
RBe = [np.array(RB.R[M].tolist(), dtype=complex) for M in range(10)]
RbBe = [np.array(RB.Rb[M].tolist(), dtype=complex) for M in range(10)]
# pairing of engine B: Mp R_M symmetric and Rb_M^T Mp symmetric (solve, as engine_b.realization.pairing)
rows = []
for M in range(10):
    for a in range(16):
        for b in range(a + 1, 16):
            r1 = np.zeros((16, 16), complex); r2 = np.zeros((16, 16), complex)
            r1[a, :] += RBe[M][:, b]; r1[b, :] -= RBe[M][:, a]
            r2[:, b] += RbBe[M][:, a]; r2[:, a] -= RbBe[M][:, b]
            rows.append(r1.ravel()); rows.append(r2.ravel())
_u, _s, _vh = np.linalg.svd(np.array(rows))
MpB = _vh[-1].conj().reshape(16, 16)
_nz = np.argwhere(np.abs(MpB) > 1e-9)
MpB = MpB / MpB[tuple(_nz[0])]                      # engine B normalisation: first nonzero entry = +1
record('B0', 'engine B pairing unique up to scale (nullity 1)', np.sum(_s < 1e-8) - 1)
# basis change: engine B Sigma = (q_a^i [a=0,1, A=i], qb^{adot}_i upper dotted) ; Sigma^* = (sb^{adot i}, s_{a i} lower)
BS = np.zeros((16, 16), complex); BSs = np.zeros((16, 16), complex)
for a in range(2):
    for i in range(4):
        BS[a * 4 + i, a * 4 + i] = 1                                   # q_a^i
        for b in range(2):
            BS[8 + a * 4 + i, 8 + b * 4 + i] = EPSD[a, b]              # qb_{a i} = eps_{ab} qb^{b}_i
            BSs[a * 4 + i, 8 + b * 4 + i] = EPSU[a, b]                 # s^a_i = eps^{ab} s_{b i}
        BSs[8 + a * 4 + i, a * 4 + i] = 1                              # sb^{a i}
D2B = [1, 2, 3, 0, 4, 5, 6, 7, 8, 9]                                   # D index (m=1..4, I) -> engine B index
record('B1', 'rho_D(e_M) = B* rho_engineB(e_M) B^{-1}  (M: m=1,2,3,4 -> engine B 1,2,3,0)',
       [RD[M] - BSs @ RBe[D2B[M]] @ np.linalg.inv(BS) for M in range(10)] +
       [RbD[M] - BS @ RbBe[D2B[M]] @ np.linalg.inv(BSs) for M in range(10)])
record('B2', 'pairing_D = pairing of engine B (normalised by engine B: first nonzero entry +1), factor 1',
       PD - np.linalg.inv(BS).T @ MpB @ np.linalg.inv(BSs))

# ------------------------------------------------------------------------------------------------------
log('=== 3. intertwiner D -> classical realisation = engine A ===')
sA = [NM.SIG_RE[M] + 1j * NM.SIG_IM[M] for M in range(10)]
sbA = [NM.SIGBAR_RE[M] + 1j * NM.SIGBAR_IM[M] for M in range(10)]
record('A0', 'classical/n4_model.py realisation == engine_a/cliff.py realisation',
       [sA[M] - CA.s[M] for M in range(10)] + [sbA[M] - CA.sb[M] for M in range(10)])
rows = []
for M in range(10):
    rows.append(np.hstack([np.kron(sA[M], E16), -np.kron(E16, RD[M].T)]))
    rows.append(np.hstack([-np.kron(E16, RbD[M].T), np.kron(sbA[M], E16)]))
_u, _s, _vh = np.linalg.svd(np.vstack(rows))
record('A1', 'intertwiner (U: Sigma_D -> Sigma_A, U*: Sigma*_D -> Sigma*_A) unique up to scale, vector map M -> M',
       np.sum(_s < 1e-8) - 1)
X = _vh[-1].conj()
U = X[:256].reshape(16, 16); Us = X[256:].reshape(16, 16)
mu = (U.T @ Us)                       # classical pairing (a,b)_A = a^T b
lam_ = mu[0, 0] / PD[0, 0] if abs(PD[0, 0]) > 1e-12 else None
kk = np.argwhere(np.abs(PD) > 1e-9)[0]
ratio = mu[tuple(kk)] / PD[tuple(kk)]
U = U / np.sqrt(ratio); Us = Us / np.sqrt(ratio)
record('A2', 'rho_A(e_M) U = U* rho_D(e_M), rho_A U* = U rho_D, (U a, U* b)_A = (a,b)_D',
       [sA[M] @ U - Us @ RD[M] for M in range(10)] + [sbA[M] @ Us - U @ RbD[M] for M in range(10)] +
       [U.T @ Us - PD])
volA = sbA[0] @ sA[1] @ sbA[2] @ sA[3]
record('A3', 'U maps S_+ (x) W onto the +1 eigenspace of rho_A(e_0 e_1 e_2 e_3) (orientation e_0..e_3 of engine A)',
       volA @ U - U @ np.diag(2 * proj_plus - 1))
Uinv = np.linalg.inv(U); Usinv = np.linalg.inv(Us)


# ======================================================================================================
# helpers: colour (su(2), [X_a,X_b] = eps_abc X_c, (X_a,X_b) = kappa delta_ab), jets, plane waves
# ======================================================================================================
KAPPA = 0.83
NC = 3


def br(X, Y):
    X, Y = np.broadcast_arrays(X, Y)
    return np.cross(X, Y)


def gp(X, Y):
    return KAPPA * np.sum(X * Y)


def jet(shape, kind='rand', p=None, amp=None):
    """value, gradient d[m], hessian dd[m,n] at the evaluation point"""
    if kind == 'rand':
        v = rnd(*shape); d = rnd(4, *shape); dd = rnd(4, 4, *shape); dd = 0.5 * (dd + np.swapaxes(dd, 0, 1))
    else:                    # plane wave amp * exp(i p.x) at x = 0
        v = amp; d = np.array([1j * p[m] * amp for m in range(4)]); dd = np.array([[-p[m] * p[n] * amp for n in range(4)] for m in range(4)])
    return dict(v=v, d=d, dd=dd)


def zero_jet(shape):
    return dict(v=np.zeros(shape, complex), d=np.zeros((4,) + shape, complex), dd=np.zeros((4, 4) + shape, complex))


def mapS(Mat, j):
    """apply a 16x16 matrix on the spinor axis of a spinor-valued jet"""
    return dict(v=np.einsum('xy,y...->x...', Mat, j['v']), d=np.einsum('xy,my...->mx...', Mat, j['d']),
                dd=np.einsum('xy,mny...->mnx...', Mat, j['dd']))


def jadd(a, b, ca=1.0, cb=1.0):
    return {k: ca * a[k] + cb * b[k] for k in a}


def waves(p_list, shapes):
    return [jet(sh, kind='wave', p=p, amp=rnd(*sh)) for p, sh in zip(p_list, shapes)]


def momenta(n):
    ps = [rng.normal(size=4) for _ in range(n - 1)]
    ps.append(-sum(ps))
    return ps


# two-component views of a D-basis spinor jet
def qparts(j):
    q, qb = split(j['v'])
    dq = np.array([split(j['d'][m])[0] for m in range(4)]); dqb = np.array([split(j['d'][m])[1] for m in range(4)])
    return q, qb, dq, dqb


# bispinor of a W_perp vector (components u_I on axis 0):  u_ij = sum_I u_I Sig^I_ij ; u^ij = -sum_I u_I Sigb^{I ij}
def bis_lo(u):
    return np.einsum('I...,Iij->ij...', u, np.array(SIGI))


def bis_up(u):
    return -np.einsum('I...,Iij->ij...', u, np.array(SIGBI))


log('=== 4. field dictionary: scalars as an antisymmetric SU(4) bispinor ===')
u = rnd(6); w = rnd(6)
record('S1', 'A^{ij} = (1/2) eps^{ijkl} A_{kl}', bis_up(u) - 0.5 * np.einsum('ijkl,kl->ij', EPS4, bis_lo(u)))
record('S2', 'sum_I u_I w_I = (1/4) u^{ij} w_{ij};  u_I = (1/4) Sigb^{I ji} u_{ij}... inverse map',
       [np.sum(u * w) - 0.25 * np.einsum('ij,ij->', bis_up(u), bis_lo(w)),
        u - 0.25 * np.einsum('Iji,ij->I', np.array(SIGBI), bis_lo(u))])
record('S3', 'Clifford in bispinor form: u^{ik} w_{kj} + w^{ik} u_{kj} = -2 (u,w) delta^i_j',
       bis_up(u) @ bis_lo(w) + bis_up(w) @ bis_lo(u) + 2 * np.sum(u * w) * np.eye(4))
ur = rng.normal(size=6)
record('S4', 'real slice (u_I real): complex conjugate of u_{ij} equals u^{ij}', bis_lo(ur).conj() - bis_up(ur))
X1 = rnd(6, NC); X2 = rnd(6, NC); X3 = rnd(6, NC); X4 = rnd(6, NC)
lhs = 0.5 * sum(gp(br(X1[I], X2[J]), br(X3[I], X4[J])) for I in range(6) for J in range(6))
rhs = (1 / 32) * KAPPA * np.sum(br(bis_up(X1)[:, :, None, None, :], bis_up(X2)[None, None, :, :, :]) *
                                br(bis_lo(X3)[:, :, None, None, :], bis_lo(X4)[None, None, :, :, :]))
record('S5', '(1/2) sum_IJ ([A_I,A_J],[A_I,A_J]) = (1/32) ([A^{ij},A^{kl}],[A_{ij},A_{kl}])  (polarised)', lhs - rhs)


# ======================================================================================================
# 5. Gamma pairing in two-component form
# ======================================================================================================
log('=== 5. Gamma ===')


def G2C(Q1, Q2):
    """(G1) Gamma^m = i(qb2 sbar^m q1 + qb1 sbar^m q2);  (G2) Gamma^I = Sig^I_ji q1_a^j eps^ab q2_b^i + Sigb^{I ji} qb1_{a j} eps^{ab} qb2_{b i}"""
    q1, qb1 = split(Q1); q2, qb2 = split(Q2)
    out = []
    for m in range(4):
        out.append(1j * (np.einsum('di,da,ai->', qb2, SIGBE[m], q1) + np.einsum('di,da,ai->', qb1, SIGBE[m], q2)))
    for I in range(6):
        out.append(np.einsum('ji,aj,ab,bi->', SIGI[I], q1, EPSU, q2) + np.einsum('ji,aj,ab,bi->', SIGBI[I], qb1, EPSU, qb2))
    return np.array(out)


def G2Cs(S1, S2):
    """(G3) Sigma^* x Sigma^*: Gamma^m = -i(s2 sig^m sb1 + s1 sig^m sb2);
            Gamma^I = Sigb^{I kj} eps_ab s1^b_j s2^a_k + Sig^I_kj eps_{ba} sb1^{a j} sb2^{b k}"""
    s1, sb1 = split(S1); s2, sb2 = split(S2)
    out = []
    for m in range(4):
        out.append(-1j * (np.einsum('ai,ad,di->', s2, SIGE[m], sb1) + np.einsum('ai,ad,di->', s1, SIGE[m], sb2)))
    for I in range(6):
        out.append(np.einsum('kj,ab,bj,ak->', SIGBI[I], EPSD, s1, s2) + np.einsum('kj,ba,aj,bk->', SIGI[I], EPSD, sb1, sb2))
    return np.array(out)


def GammaA(Q1, Q2):        # index-free, classical realisation, Q in Sigma_A
    return np.array([Q1 @ (sA[M] @ Q2) for M in range(10)])


def GammaAs(S1, S2):       # Sigma^* x Sigma^*: (rho(e_M)S1, S2)
    return np.array([(sbA[M] @ S1) @ S2 for M in range(10)])


for _ in range(3):
    Q1 = rnd(16); Q2 = rnd(16); S1 = rnd(16); S2 = rnd(16)
    record('G1', 'Gamma(Q1,Q2) two-component (G1),(G2) == index-free (classical realisation)', G2C(Q1, Q2) - GammaA(U @ Q1, U @ Q2))
    record('G3', 'Gamma(S1,S2) on Sigma^* two-component (G3) == index-free', G2Cs(S1, S2) - GammaAs(Us @ S1, Us @ S2))


def Gbis_lo(Q1, Q2):
    """(G4) Gamma_{kl} = 2 eps_{jikl} q1_a^j eps^{ab} q2_b^i - 2(qb1_{a k} eps^{ab} qb2_{b l} - (k<->l))"""
    q1, qb1 = split(Q1); q2, qb2 = split(Q2)
    t = 2 * np.einsum('jikl,aj,ab,bi->kl', EPS4, q1, EPSU, q2)
    x = np.einsum('ak,ab,bl->kl', qb1, EPSU, qb2)
    return t - 2 * (x - x.T)


def Gbis_up_s(S1, S2):
    """(G5) Gamma^{kl}(S1,S2) = -2 eps^{klpj} eps_ab s1^b_j s2^a_p + 2(eps_{ba} sb1^{a l} sb2^{b k} - (k<->l))"""
    s1, sb1 = split(S1); s2, sb2 = split(S2)
    t = -2 * np.einsum('klpj,ab,bj,ap->kl', EPS4, EPSD, s1, s2)
    x = np.einsum('ba,al,bk->kl', EPSD, sb1, sb2)
    return t + 2 * (x - x.T)


Q1 = rnd(16); Q2 = rnd(16); S1 = rnd(16); S2 = rnd(16)
record('G4', 'bispinor Gamma_{kl}(Q1,Q2) = sum_I Gamma^I Sig^I_kl (G4)', Gbis_lo(Q1, Q2) - bis_lo(G2C(Q1, Q2)[4:]))
record('G5', 'bispinor Gamma^{kl}(S1,S2) = -sum_I Gamma^I Sigb^{I kl} (G5)', Gbis_up_s(S1, S2) - bis_up(G2Cs(S1, S2)[4:]))
# square-zero elements: rank-one q_a^i = u_a w^i, qb = 0
uu = rnd(2); ww = rnd(4)
Qr = join(np.einsum('a,i->ai', uu, ww), np.zeros((2, 4)))
record('G6', 'xi_alpha^i = u_alpha w^i, xibar = 0  =>  Gamma(Q,Q) = 0 (square-zero element of Sigma)', G2C(Qr, Qr))


# ======================================================================================================
# 6. The action
# ======================================================================================================
log('=== 6. action ===')


def F_IF(Aj):
    """F_{MN} (10,10,NC) from a jet of A in A-basis (= D basis on vector indices)"""
    A = Aj['v']; d = Aj['d']
    F = np.zeros((10, 10, NC), complex)
    for M in range(10):
        for N in range(10):
            dMN = d[M][N] if M < 4 else 0
            dNM = d[N][M] if N < 4 else 0
            F[M, N] = dMN - dNM + br(A[M], A[N])
    return F


def hodge(F, G):
    return sum(gp(F[M, N], G[M, N]) for M in range(10) for N in range(M + 1, 10))


def kin_gauge_2C(Aj):
    """-1/4(F,F) = -1/8 (F_mn,F_mn) - 1/16 (D_m A^{ij}, D_m A_{ij}) - 1/128 ([A^{ij},A^{kl}],[A_{ij},A_{kl}])"""
    A = Aj['v']; d = Aj['d']
    Am = A[:4]; AI = A[4:]
    Fmn = np.array([[d[m][n] - d[n][m] + br(Am[m], Am[n]) for n in range(4)] for m in range(4)])
    up = bis_up(AI); lo = bis_lo(AI)
    Dup = np.array([bis_up(d[m][4:]) + br(Am[m], up) for m in range(4)])
    Dlo = np.array([bis_lo(d[m][4:]) + br(Am[m], lo) for m in range(4)])
    quart = KAPPA * np.sum(br(up[:, :, None, None, :], up[None, None, :, :, :]) * br(lo[:, :, None, None, :], lo[None, None, :, :, :]))
    return -gp(Fmn, Fmn) / 8 - gp(Dup, Dlo) / 16 - quart / 128


for _ in range(2):
    Aj = jet((10, NC))
    record('L1', 'gauge+scalar kinetic and quartic terms (pointwise)', kin_gauge_2C(Aj) - (-0.25 * hodge(F_IF(Aj), F_IF(Aj))))


def dslashA_IF(lj, Aj):
    """d-slash_A lambda in A-basis: rho(e_m) d_m lambda + rho(e_M)[A_M,lambda]"""
    out = sum(sA[m] @ lj['d'][m] for m in range(4))
    for M in range(10):
        out = out + sA[M] @ br(Aj['v'][M][None, :], lj['v'])
    return out


def ferm_IF(l1, l2, Aj):
    return 0.5 * KAPPA * np.sum(l1['v'] * dslashA_IF(l2, Aj))


def ferm_2C(l1, l2, Aj):
    """(1/2)(lambda, dslash_A lambda) pointwise =
       (i/2)[(lam_a^i, sbar^{m d a} D_m lambar_{d i}) + (lambar_{d i}, sbar^{m d a} D_m lam_a^i)]
       + (1/2)(lam_a^j, eps^{ab}[A_{ji}, lam_b^i]) - (1/2)(lambar_{d j}, eps^{de}[A^{ji}, lambar_{e i}])"""
    q1, qb1, _, _ = qparts(l1)
    q2, qb2, dq2, dqb2 = qparts(l2)
    Am = Aj['v'][:4]; AI = Aj['v'][4:]
    Dq2 = np.array([dq2[m] + br(Am[m], q2) for m in range(4)])
    Dqb2 = np.array([dqb2[m] + br(Am[m], qb2) for m in range(4)])
    kin = 0.5j * KAPPA * (sum(np.einsum('aic,da,dic->', q1, SIGBE[m], Dqb2[m]) for m in range(4)) +
                          sum(np.einsum('dic,da,aic->', qb1, SIGBE[m], Dq2[m]) for m in range(4)))
    lo = bis_lo(AI); up = bis_up(AI)
    yuk1 = 0.5 * KAPPA * np.einsum('ajc,ab,jibc->', q1, EPSU, br(lo[:, :, None, :], np.transpose(q2, (1, 0, 2))[None, :, :, :]))
    yuk2 = -0.5 * KAPPA * np.einsum('djc,de,jiec->', qb1, EPSU, br(up[:, :, None, :], np.transpose(qb2, (1, 0, 2))[None, :, :, :]))
    return kin + yuk1 + yuk2


for _ in range(2):
    Aj = jet((10, NC)); l1 = jet((16, NC)); l2 = jet((16, NC))
    vIF = ferm_IF(mapS(U, l1), mapS(U, l2), Aj) - ferm_IF(mapS(U, l2), mapS(U, l1), Aj)
    v2C = ferm_2C(l1, l2, Aj) - ferm_2C(l2, l1, Aj)
    record('L2', '(1/2)(lambda, dslash_A lambda) pointwise, Grassmann (theta1 theta2 coefficient)', v2C - vIF)
# integrated free form:  int (1/2)(lambda,dslash lambda) = int i lambar_{d i} sbar^{m d a} d_m lambda_a^i
for _ in range(2):
    p1, p2 = momenta(2)
    l1, l2 = waves([p1, p2], [(16, NC), (16, NC)])
    Z = zero_jet((10, NC))
    vIF = ferm_IF(mapS(U, l1), mapS(U, l2), Z) - ferm_IF(mapS(U, l2), mapS(U, l1), Z)

    def kin_int(a, b):
        _, qba, _, _ = qparts(a); _, _, dqb, _ = qparts(b)
        _, _, dq, _ = qparts(b)
        return 1j * KAPPA * sum(np.einsum('dic,da,aic->', qba, SIGBE[m], dq[m]) for m in range(4))
    record('L3', 'int (1/2)(lambda,dslash lambda) = int i lambar_i sbar^m d_m lambda^i  (plane waves, Grassmann)',
           (kin_int(l1, l2) - kin_int(l2, l1)) - vIF)


def ghost_IF(cj, Aj, Asj):
    out = sum(gp(cj['d'][m] + br(Aj['v'][m], cj['v']), Asj['v'][m]) for m in range(4))
    out += sum(gp(br(Aj['v'][M], cj['v']), Asj['v'][M]) for M in range(4, 10))
    return out


def ghost_2C(cj, Aj, Asj):
    """(d_A c, A^*) = (D_m c, A^*_m) + (1/4)([A^{ij}, c], A^*_{ij})"""
    out = sum(gp(cj['d'][m] + br(Aj['v'][m], cj['v']), Asj['v'][m]) for m in range(4))
    out += 0.25 * gp(br(bis_up(Aj['v'][4:]), cj['v']), bis_lo(Asj['v'][4:]))
    return out


cj = jet((NC,)); Aj = jet((10, NC)); Asj = jet((10, NC))
record('L4', '(d_A c, A^*) two-component', ghost_2C(cj, Aj, Asj) - ghost_IF(cj, Aj, Asj))
lj = jet((16, NC)); lsj = jet((16, NC))
vIF = gp(br(cj['v'], (U @ lj['v'])), Us @ lsj['v'])
q, qb = split(lj['v']); s, sb = split(lsj['v'])
v2C = gp(br(cj['v'], q), s) - gp(br(cj['v'], qb), sb)
record('L5', '([c,lambda],lambda^*) = ([c,lambda_a^i],lambda^{*a}_i) - ([c,lambdabar_{d i}],lambdabar^{*d i})', v2C - vIF)


# ======================================================================================================
# 7. Q, Q^GF, H in two-component form (operators on slots)
# ======================================================================================================
log('=== 7. Q, Q^GF, H ===')
for _ in range(2):
    p = rng.normal(size=4)
    lam = rnd(16)
    # Q lambda = dslash lambda in slot of lambda^*: two-component
    q, qb = split(lam)
    s_2c = 1j * sum(1j * p[m] * np.einsum('di,da->ai', qb, SIGBE[m]) for m in range(4))      # lambda^{*a}_i slot
    sb_2c = -1j * sum(1j * p[m] * np.einsum('da,ai->di', SIGBE[m], q) for m in range(4))      # lambdabar^{*d i} slot
    QlamIF = sum(1j * p[m] * sA[m] for m in range(4)) @ (U @ lam)
    record('Q1', 'Q on lambda: lambda^{*a}_i <- i d_m lambar_{d i} sbar^{m d a}; lambar^{*d i} <- -i sbar^{m d a} d_m lambda_a^i',
           Us @ join(s_2c, sb_2c) - QlamIF)
    lams = rnd(16)
    s, sb = split(lams)
    q_2c = 1j * sum(1j * p[m] * np.einsum('ad,di->ai', SIGE[m], sb) for m in range(4))        # lambda_a^i slot
    qb_2c = -1j * sum(1j * p[m] * np.einsum('ai,ad->di', s, SIGE[m]) for m in range(4))       # lambar_{d i} slot
    QgfIF = -sum(1j * p[m] * sbA[m] for m in range(4)) @ (Us @ lams)
    record('Q2', 'Q^GF on lambda^*: lambda_a^i <- i sig^m d_m lambdabar^{*}; lambar_{d i} <- -i d_m lambda^{*a}_i sig^m_{a d}',
           U @ join(q_2c, qb_2c) - QgfIF)
    # H = Q^GF Q on lambda = |p|^2
    q, qb = split(lam)
    s1_ = 1j * sum(1j * p[m] * np.einsum('di,da->ai', qb, SIGBE[m]) for m in range(4))
    sb1_ = -1j * sum(1j * p[m] * np.einsum('da,ai->di', SIGBE[m], q) for m in range(4))
    q2_ = 1j * sum(1j * p[m] * np.einsum('ad,di->ai', SIGE[m], sb1_) for m in range(4))
    qb2_ = -1j * sum(1j * p[m] * np.einsum('ai,ad->di', s1_, SIGE[m]) for m in range(4))
    record('Q3', 'H = Q^GF Q = -d^2 on lambda in two-component form (momentum: |p|^2)', join(q2_, qb2_) - np.dot(p, p) * lam)
    # scalars: Q A_I = (1/2) d^2 A_I  ->  A^*_{ij} slot: (1/2) d^2 A_{ij};  Q^GF A^*_{ij} = -2 A^*_{ij}
    AI = rnd(6)
    record('Q4', 'Q on scalars in bispinor form: A^*_{ij} <- -(1/2)|p|^2 A_{ij}; H = Q^GF Q = |p|^2',
           [bis_lo(-0.5 * np.dot(p, p) * AI) - (-0.5 * np.dot(p, p)) * bis_lo(AI),
            -2 * (-0.5 * np.dot(p, p)) * bis_lo(AI) - np.dot(p, p) * bis_lo(AI)])


# ======================================================================================================
# 8. propagators
# ======================================================================================================
log('=== 8. propagators ===')
for _ in range(2):
    p = rng.normal(size=4)
    T_D = 1j * sum(p[m] * RbD[m] for m in range(4))          # i rho(p) : Sigma^* -> Sigma
    t_D = T_D @ np.linalg.inv(PD)                              # t = sum_b T(f_b) (x) e_b, (e_a,f_b) = delta_ab
    t_2c = np.zeros((16, 16), complex)
    ps = sum(p[m] * SIGE[m] for m in range(4))
    for a in range(2):
        for i in range(4):
            for b in range(2):
                t_2c[a * 4 + i, 8 + b * 4 + i] = -ps[a, b]     # (lambda_a^i, lambar_{b j}) : -p.sig_{a b} delta^i_j
                t_2c[8 + a * 4 + i, b * 4 + i] = -ps[b, a]     # (lambar_{a i}, lambda_b^j) : -p.sig_{b a} delta_i^j
    record('P1', 'lambda (x) lambda kernel: (lambda_a^i, lambar_{b j}) -> -p_m sig^m_{a b} delta^i_j, '
           '(lambar_{a i}, lambda_b^j) -> -p_m sig^m_{b a} delta_i^j', t_D - t_2c)
    T_A = 1j * sum(p[m] * sbA[m] for m in range(4))
    record('P2', 'same kernel in the classical realisation: t_A = U t_D U^T', T_A - U @ t_D @ U.T)
# A (x) A in bispinor form: 2 delta_IJ -> (A_{ij},A_{kl}): 4 eps_{ijkl}; (A_{ij},A^{kl}): 4(d^k_i d^l_j - d^l_i d^k_j)
tII = 2 * np.eye(6)
record('P3', 'A (x) A on W_perp in bispinor form: <A_ij A_kl> = 4 eps_ijkl, <A_ij A^kl> = 4(d_i^k d_j^l - d_i^l d_j^k)',
       [np.einsum('IJ,Iij,Jkl->ijkl', tII, np.array(SIGI), np.array(SIGI)) - 4 * EPS4,
        np.einsum('IJ,Iij,Jkl->ijkl', tII, np.array(SIGI), -np.array(SIGBI)) - 4 * np.einsum('ijkl->klij', DD)])


# ======================================================================================================
# 9. psl(4|4)
# ======================================================================================================
log('=== 9. psl(4|4) ===')


def vfield(kind, par, x):
    """vector field v^m(x) on V, Jacobian J[m,n] = d_m v^n, divergence, so(W_perp) matrix R (v(y)^I = R_IJ y_J)"""
    R6 = np.zeros((6, 6))
    if kind == 'P':
        v = par.copy(); J = np.zeros((4, 4))
    elif kind == 'M':
        w = par; v = w @ x; J = w.T.copy()                        # v^m = w_mn x_n, d_m v^n = w_nm
    elif kind == 'D':
        v = 2 * x; J = 2 * np.eye(4)
    elif kind == 'K':
        b = par; v = np.dot(x, x) * b - 2 * np.dot(b, x) * x
        J = 2 * np.outer(x, b) - 2 * np.outer(b, x) - 2 * np.dot(b, x) * np.eye(4)
    elif kind == 'R':
        v = np.zeros(4); J = np.zeros((4, 4)); R6 = par
    return v, J, np.trace(J), R6


def even_IF(kind, par, x, Aj, lj, cj):
    """eq:c2-S1even: v.A = L_v A + 1/4 div v A|perp, v.lambda = v(lambda) + 1/4 rho(d(v,-)) lambda + 3/8 div v lambda, v.c = v(c)"""
    v, J, div, R6 = vfield(kind, par, x)
    A = Aj['v']; d = Aj['d']
    vA = np.zeros((10, NC), complex)
    for m in range(4):
        vA[m] = sum(v[n] * d[n][m] for n in range(4)) + sum(A[n] * J[m, n] for n in range(4))
    for I in range(6):
        vA[4 + I] = sum(v[n] * d[n][4 + I] for n in range(4)) + sum(R6[J_, I] * A[4 + J_] for J_ in range(6)) + 0.25 * div * A[4 + I]
    Om = np.zeros((10, 10))
    for m in range(4):
        for n in range(4):
            Om[m, n] = J[m, n] - J[n, m]
    for I in range(6):
        for J_ in range(6):
            Om[4 + J_, 4 + I] = R6[I, J_] - R6[J_, I]             # d(v,-) = sum R_IJ dy_J ^ dy_I
    spin = sum(Om[M, N] * sbA[M] @ sA[N] for M in range(10) for N in range(M + 1, 10))
    vl = sum(v[n] * lj['d'][n] for n in range(4)) + 0.25 * spin @ lj['v'] + 0.375 * div * lj['v']
    vc = sum(v[n] * cj['d'][n] for n in range(4))
    return vA, vl, vc


def rmat(R6):
    """so(W_perp) -> sl(4):  r^k_i = (1/4) R_IJ (Sigb^J Sig^I)^k_i"""
    return 0.25 * sum(R6[I, J] * SIGBI[J] @ SIGI[I] for I in range(6) for J in range(6))


def even_2C(kind, par, x, Aj, lj, cj):
    """two-component table (E1)-(E4):
       v.A_m = v^n d_n A_m + A_n d_m v^n
       v.A_{ij} = v^n d_n A_{ij} + (1/4) div v A_{ij} - r^k_i A_{kj} - r^k_j A_{ik}
       v.lambda_a^i = v^n d_n lambda_a^i - (1/8) X_{mn} (sig^m sbar^n)_a^b lambda_b^i + (3/8) div v lambda_a^i + r^i_j lambda_a^j
       v.lambar_{d i} = v^n d_n lambar_{d i} + (1/8) X_{mn} lambar_{e i} (sbar^m sig^n)^e_d + (3/8) div v lambar_{d i} - lambar_{d j} r^j_i
       X_{mn} = d_m v_n - d_n v_m,  v.c = v^n d_n c"""
    v, J, div, R6 = vfield(kind, par, x)
    r = rmat(R6)
    A = Aj['v']; d = Aj['d']
    vAm = np.array([sum(v[n] * d[n][m] for n in range(4)) + sum(A[n] * J[m, n] for n in range(4)) for m in range(4)])
    lo = bis_lo(A[4:]); dlo = np.array([bis_lo(d[n][4:]) for n in range(4)])
    vlo = sum(v[n] * dlo[n] for n in range(4)) + 0.25 * div * lo - np.einsum('ki,kjc->ijc', r, lo) - np.einsum('kj,ikc->ijc', r, lo)
    q, qb, dq, dqb = qparts(lj)
    Xmn = J - J.T
    sq = sum(Xmn[m, n] * SIGE[m] @ SIGBE[n] for m in range(4) for n in range(4))
    sqb = sum(Xmn[m, n] * SIGBE[m] @ SIGE[n] for m in range(4) for n in range(4))
    vq = sum(v[n] * dq[n] for n in range(4)) - 0.125 * np.einsum('ab,bic->aic', sq, q) + 0.375 * div * q + np.einsum('ij,ajc->aic', r, q)
    vqb = sum(v[n] * dqb[n] for n in range(4)) + 0.125 * np.einsum('eic,ed->dic', qb, sqb) + 0.375 * div * qb - np.einsum('djc,ji->dic', qb, r)
    vc = sum(v[n] * cj['d'][n] for n in range(4))
    return vAm, vlo, join(vq, vqb), vc


R6r = rng.normal(size=(6, 6)); R6r = R6r - R6r.T
wr = rng.normal(size=(4, 4)); wr = wr - wr.T
for kind, par in (('P', rng.normal(size=4)), ('M', wr), ('D', None), ('K', rng.normal(size=4)), ('R', R6r)):
    x = rng.normal(size=4)
    Aj = jet((10, NC)); lj = jet((16, NC)); cj = jet((NC,))
    vA, vl, vc = even_IF(kind, par, x, Aj, mapS(U, lj), cj)
    vAm, vlo, vlam, vc2 = even_2C(kind, par, x, Aj, lj, cj)
    v_, J_, div_, R6_ = vfield(kind, par, x); r_ = rmat(R6_)
    up_A = bis_up(Aj['v'][4:]); dup_A = np.array([bis_up(Aj['d'][n][4:]) for n in range(4)])
    vup = sum(v_[n] * dup_A[n] for n in range(4)) + 0.25 * div_ * up_A + np.einsum('ik,kjc->ijc', r_, up_A) + np.einsum('jk,ikc->ijc', r_, up_A)
    record('E-' + kind, f'S^(1)(v) on fields, generator {kind}: two-component table (incl. v.A^ij) vs eq:c2-S1even',
           [vAm - vA[:4], vlo - bis_lo(vA[4:]), vup - bis_up(vA[4:]), U @ vlam - vl, vc2 - vc])


# odd generators ---------------------------------------------------------------------------------------
def rhoF_IF(F, Q):          # Q in Sigma_A (no colour), F (10,10,NC) -> (16,NC)
    out = np.zeros((16, NC), complex)
    for M in range(10):
        for N in range(M + 1, 10):
            out += np.einsum('x,c->xc', sbA[M] @ sA[N] @ Q, F[M, N])
    return out


def S1Q_IF(Q, Aj, lj, Asj, lsj):
    """eq:c2-S1Q density: -(Gamma(Q,lambda),A^*) + 1/2 (rho(F_A)Q, lambda^*)   (Q in Sigma_A)"""
    G = np.array([np.einsum('x,xc->c', Q, sA[M] @ lj['v']) for M in range(10)])
    return -gp(G, Asj['v']) + 0.5 * gp(rhoF_IF(F_IF(Aj), Q), lsj['v'])


def rhoF_2C(Aj, xi):
    """(F1) (rho(F_A) xi)_a^i = -1/2 F_mn (sig^m sbar^n xi)_a^i - 1/2 [A^{ik},A_{kj}] xi_a^j - i sig^m_{ad} eps^{de} D_m A^{ij} xibar_{e j}
            (rho(F_A) xi)_{d i} = +1/2 F_mn (xibar sbar^m sig^n)_{d i} - 1/2 [A_{ik},A^{kj}] xibar_{d j} + i eps^{ab} sig^m_{a d} D_m A_{ij} xi_b^j"""
    A = Aj['v']; d = Aj['d']
    Am = A[:4]
    Fmn = np.array([[d[m][n] - d[n][m] + br(Am[m], Am[n]) for n in range(4)] for m in range(4)])
    lo = bis_lo(A[4:]); up = bis_up(A[4:])
    Dlo = np.array([bis_lo(d[m][4:]) + br(Am[m], lo) for m in range(4)])
    Dup = np.array([bis_up(d[m][4:]) + br(Am[m], up) for m in range(4)])
    xq, xqb = split(xi)
    cuu = br(up[:, :, None, None, :], lo[None, None, :, :, :])          # [A^{ik}, A_{lj}] -> index (i,k,l,j)
    comm_ul = np.einsum('ikkjc->ijc', cuu)                              # [A^{ik},A_{kj}]
    cll = br(lo[:, :, None, None, :], up[None, None, :, :, :])
    comm_lu = np.einsum('ikkjc->ijc', cll)                              # [A_{ik},A^{kj}]
    oq = (-0.5 * sum(np.einsum('c,ab,bi->aic', Fmn[m, n], SIGE[m] @ SIGBE[n], xq) for m in range(4) for n in range(4))
          - 0.5 * np.einsum('ijc,aj->aic', comm_ul, xq)
          - 1j * sum(np.einsum('ad,de,ijc,ej->aic', SIGE[m], EPSU, Dup[m], xqb) for m in range(4)))
    oqb = (0.5 * sum(np.einsum('c,ei,ed->dic', Fmn[m, n], xqb, SIGBE[m] @ SIGE[n]) for m in range(4) for n in range(4))
           - 0.5 * np.einsum('ijc,dj->dic', comm_lu, xqb)
           + 1j * sum(np.einsum('ab,ad,ijc,bj->dic', EPSU, SIGE[m], Dlo[m], xq) for m in range(4)))
    return oq, oqb


def Gam2C_field(xi, lam):
    """Gamma(xi, lambda) with lambda colour valued: Gamma^m and bispinor Gamma_{kl}  (G1),(G4)"""
    xq, xqb = split(xi); lq, lqb = split(lam)
    Gm = np.array([1j * (np.einsum('dic,da,ai->c', lqb, SIGBE[m], xq) + np.einsum('di,da,aic->c', xqb, SIGBE[m], lq)) for m in range(4)])
    t = 2 * np.einsum('jikl,aj,ab,bic->klc', EPS4, xq, EPSU, lq)
    xx = np.einsum('ak,ab,blc->klc', xqb, EPSU, lqb)
    Glo = t - 2 * (xx - np.transpose(xx, (1, 0, 2)))
    return Gm, Glo


def S1Q_2C(xi, Aj, lj, Asj, lsj):
    """S^(1)(Q_1) density = -(Gamma^m(xi,lambda), A^*_m) - (1/4)(Gamma^{ij}(xi,lambda), A^*_{ij})... written with Gamma_{ij} and A^{*ij}:
       -(Gamma^m, A^*_m) - (1/4)(Gamma_{ij}, A^{*ij}) + 1/2[((rho F xi)_a^i, lambda^{*a}_i) - ((rho F xi)_{d i}, lambdabar^{*d i})]"""
    Gm, Glo = Gam2C_field(xi, lj['v'])
    t1 = -gp(Gm, Asj['v'][:4]) - 0.25 * gp(Glo, bis_up(Asj['v'][4:]))
    oq, oqb = rhoF_2C(Aj, xi)
    s, sb = split(lsj['v'])
    t2 = 0.5 * (gp(oq, s) - gp(oqb, sb))
    return t1 + t2


for _ in range(2):
    xi = rnd(16)
    Aj = jet((10, NC)); lj = jet((16, NC)); Asj = jet((10, NC)); lsj = jet((16, NC))
    record('O1', 'S^(1)(Q_1) density, Q_1 = xi in Sigma: two-component (F1) vs eq:c2-S1Q',
           S1Q_2C(xi, Aj, lj, Asj, lsj) - S1Q_IF(U @ xi, Aj, mapS(U, lj), Asj, mapS(Us, lsj)))


def killing_2C(eta, x):
    """(K1) epsilon(x) = rho(x) eta:  eps_a^i = -i x_m sig^m_{ad} etabar^{d i},  epsbar_{d i} = +i x_m eta^a_i sig^m_{a d}"""
    s, sb = split(eta)
    xs = sum(x[m] * SIGE[m] for m in range(4))
    return join(-1j * np.einsum('ad,di->ai', xs, sb), 1j * np.einsum('ai,ad->di', s, xs))


def S1S_IF(eta_A, x, Aj, lj, Asj, lsj):
    """eq:c2-S1S"""
    ex = sum(x[m] * sbA[m] for m in range(4)) @ eta_A
    base = S1Q_IF(ex, Aj, lj, Asj, lsj)
    rhoAperp = np.zeros((16, NC), complex)
    for M in range(4, 10):
        rhoAperp += np.einsum('x,c->xc', sbA[M] @ eta_A, Aj['v'][M])
    return base - gp(rhoAperp, lsj['v'])


def S1S_2C(eta, x, Aj, lj, Asj, lsj):
    """(F2) S^(1)(Q_2) = (F1 with xi -> epsilon(x)) + (A^{kj} eps_ab eta^b_j, lambda^{*a}_k) - (A_{kj} eps_{ed} etabar^{d j}, lambdabar^{*e k})"""
    ex = killing_2C(eta, x)
    base = S1Q_2C(ex, Aj, lj, Asj, lsj)
    s, sb = split(eta)
    ls, lsb = split(lsj['v'])
    lo = bis_lo(Aj['v'][4:]); up = bis_up(Aj['v'][4:])
    t3 = KAPPA * (np.einsum('kjc,ab,bj,akc->', up, EPSD, s, ls) - np.einsum('kjc,ed,dj,ekc->', lo, EPSD, sb, lsb))
    return base + t3


for _ in range(2):
    eta = rnd(16); x = rng.normal(size=4)
    Aj = jet((10, NC)); lj = jet((16, NC)); Asj = jet((10, NC)); lsj = jet((16, NC))
    record('O2', 'Killing spinor rho(x)eta two-component (K1) vs index-free', U @ killing_2C(eta, x) - sum(x[m] * sbA[m] for m in range(4)) @ (Us @ eta))
    record('O3', 'S^(1)(Q_2) density, Q_2 = eta in Sigma^*: two-component (F2) vs eq:c2-S1S',
           S1S_2C(eta, x, Aj, lj, Asj, lsj) - S1S_IF(Us @ eta, x, Aj, mapS(U, lj), Asj, mapS(Us, lsj)))


def S2_IF(Q1A, Q2A, Aj, lsj, csj):
    """eq:c2-S2 density for Q1,Q2 in Sigma_A (Killing spinors already evaluated)"""
    G12 = GammaA(Q1A, Q2A)
    ls = lsj['v']
    Gss = np.array([np.einsum('xc,xd->cd', sbA[M] @ ls, ls) for M in range(10)])      # Gamma(lambda^*,lambda^*) in g (x) g
    t1 = 0.25 * KAPPA * np.einsum('M,Mcc->', G12, Gss)
    p1 = np.einsum('x,xc->c', Q1A, ls); p2 = np.einsum('x,xc->c', Q2A, ls)
    t2 = -0.5 * gp(p1, p2)
    t3 = -gp(np.einsum('M,Mc->c', G12, Aj['v']), csj['v'])
    return t1 + t2 + t3


def S2_2C(Q1, Q2, Aj, lsj, csj):
    """(F3) S^(2) = 1/4[Gamma^m(Q1,Q2)(-2i)(lambda^{*a}_i, sig^m_{ad} lambdabar^{*d i}) + 1/4 Gamma_{kl}(Q1,Q2) (Gamma^{kl}(lambda^*,lambda^*))]
                    - 1/2 ((Q1,lambda^*),(Q2,lambda^*)) - (Gamma^m A_m + 1/4 Gamma_{ij} A^{ij}, c^*)"""
    G = G2C(Q1, Q2); Glo = Gbis_lo(Q1, Q2)
    ls, lsb = split(lsj['v'])
    vm = np.array([-2j * KAPPA * np.einsum('aic,ad,dic->', ls, SIGE[m], lsb) for m in range(4)])
    tt = -2 * np.einsum('klpj,ab,bjc,apc->kl', EPS4, EPSD, ls, ls)
    xx = np.einsum('ba,alc,bkc->kl', EPSD, lsb, lsb)
    Gup_ss = KAPPA * (tt + 2 * (xx - xx.T))
    t1 = 0.25 * (np.sum(G[:4] * vm) + 0.25 * np.sum(Glo * Gup_ss))
    q1, qb1 = split(Q1); q2, qb2 = split(Q2)
    p1 = np.einsum('ai,aic->c', q1, ls) - np.einsum('di,dic->c', qb1, lsb)
    p2 = np.einsum('ai,aic->c', q2, ls) - np.einsum('di,dic->c', qb2, lsb)
    t2 = -0.5 * gp(p1, p2)
    gauge = np.einsum('m,mc->c', G[:4], Aj['v'][:4]) + 0.25 * np.einsum('ij,ijc->c', Glo, bis_up(Aj['v'][4:]))
    t3 = -gp(gauge, csj['v'])
    return t1 + t2 + t3


for case in ('QQ', 'QS', 'SS'):
    x = rng.normal(size=4)
    a = rnd(16); b = rnd(16)
    Q1 = a if case[0] == 'Q' else killing_2C(a, x)
    Q2 = b if case[1] == 'Q' else killing_2C(b, x)
    Aj = jet((10, NC)); lsj = jet((16, NC)); csj = jet((NC,))
    record('O4-' + case, f'S^(2) density ({case}): two-component (F3) vs eq:c2-S2', S2_2C(Q1, Q2, Aj, lsj, csj) - S2_IF(U @ Q1, U @ Q2, Aj, mapS(Us, lsj), csj))

# brackets -----------------------------------------------------------------------------------------------
x = rng.normal(size=4)
xi1, xi3 = rnd(16), rnd(16); eta2, eta4 = rnd(16), rnd(16)
e1 = xi1 + killing_2C(eta2, x); e2 = xi3 + killing_2C(eta4, x)
record('K1', 'odd-odd bracket, V part: v^m(x) = i(epsbar1 sbar^m eps2 + epsbar2 sbar^m eps1)  [eq:c2-oddodd]',
       G2C(e1, e2)[:4] - GammaA(U @ e1, U @ e2)[:4])


def rot_IF(Q1, Q2s, Q3, Q4s):
    RR = np.zeros((6, 6), complex)
    for I in range(6):
        for J in range(6):
            if I != J:
                RR[I, J] = (sbA[4 + I] @ sA[4 + J] @ (U @ Q1)) @ (Us @ Q4s) + (sbA[4 + I] @ sA[4 + J] @ (U @ Q3)) @ (Us @ Q2s)
    return RR


def rot_2C(Q1, Q2s, Q3, Q4s):
    """(K2) r^k_i = 2[(xi1_a^k eta4^a_i + xibar1_{d i} etabar4^{d k}) + (xi3_a^k eta2^a_i + xibar3_{d i} etabar2^{d k})] - trace"""
    def one(Q, S):
        q, qb = split(Q); s, sb = split(S)
        X = np.einsum('ak,ai->ki', q, s) + np.einsum('di,dk->ki', qb, sb)
        return 2 * (X - 0.25 * np.trace(X) * np.eye(4))
    return one(Q1, Q4s) + one(Q3, Q2s)


record('K2', 'odd-odd bracket, rotation of W_perp as sl(4) matrix r^k_i = (1/4)R_IJ(Sigb^J Sig^I)^k_i  (K2)',
       rot_2C(xi1, eta2, xi3, eta4) - rmat(rot_IF(xi1, eta2, xi3, eta4)))
record('K3', 'dilatation part of [xi,eta]: (1/2)(xi,eta) Delta, Delta=2E: trace part of v^m = (xi_a^i eta^a_i - xibar_{d i} etabar^{d i}) x^m',
       0.25 * sum(np.array([G2C(xi1, killing_2C(eta2, y))[:4] for y in np.eye(4)])[m, m] for m in range(4)) - pair_D(xi1, eta2))
# even-odd bracket on the Killing spinor (same table as v.lambda with 3/8 -> -1/8): test at a point via jets of eps(x)


def eps_jet(xi, eta, x):
    v = xi + killing_2C(eta, x)
    d = np.array([killing_2C(eta, np.eye(4)[m]) for m in range(4)])
    return dict(v=v, d=d, dd=np.zeros((4, 4, 16), complex))


for kind, par in (('M', wr), ('D', None), ('K', rng.normal(size=4)), ('R', R6r)):
    x = rng.normal(size=4)
    ej = eps_jet(xi1, eta2, x)
    v, J, div, R6 = vfield(kind, par, x)
    Om = np.zeros((10, 10))
    for m in range(4):
        for n in range(4):
            Om[m, n] = J[m, n] - J[n, m]
    for I in range(6):
        for J_ in range(6):
            Om[4 + J_, 4 + I] = R6[I, J_] - R6[J_, I]
    spin = sum(Om[M, N] * sbA[M] @ sA[N] for M in range(10) for N in range(M + 1, 10))
    ejA = mapS(U, ej)
    IFv = sum(v[n] * ejA['d'][n] for n in range(4)) + 0.25 * spin @ ejA['v'] - 0.125 * div * ejA['v']
    r = rmat(R6); Xmn = J - J.T
    q, qb = split(ej['v']); dq = np.array([split(ej['d'][m])[0] for m in range(4)]); dqb = np.array([split(ej['d'][m])[1] for m in range(4)])
    sq = sum(Xmn[m, n] * SIGE[m] @ SIGBE[n] for m in range(4) for n in range(4))
    sqb = sum(Xmn[m, n] * SIGBE[m] @ SIGE[n] for m in range(4) for n in range(4))
    vq = sum(v[n] * dq[n] for n in range(4)) - 0.125 * sq @ q - 0.125 * div * q + np.einsum('ij,aj->ai', r, q)
    vqb = sum(v[n] * dqb[n] for n in range(4)) + 0.125 * np.einsum('ei,ed->di', qb, sqb) - 0.125 * div * qb - qb @ r
    record('K4-' + kind, f'even-odd bracket [v, Q_1+rho(x)Q_2] two-component (K3), generator {kind}  [eq:c2-evenodd]',
           U @ join(vq, vqb) - IFv)


# ======================================================================================================
# 10. Step-1 results and Step-2 bases
# ======================================================================================================
log('=== 10. Step-1 results and Step-2 bases ===')


def integ(fun, legs_spec, nperm_fields, odd_fields=()):
    """sum over assignments of the legs of each repeated field to its occurrences (Grassmann sign for odd fields)."""
    names = list(legs_spec.keys())
    perms = {n: list(itertools.permutations(range(len(legs_spec[n])))) for n in names}
    total = 0
    for combo in itertools.product(*[perms[n] for n in names]):
        sign = 1
        assign = {}
        for n, pm in zip(names, combo):
            assign[n] = [legs_spec[n][k] for k in pm]
            if n in odd_fields:
                sgn = 1
                pl = list(pm)
                for i_ in range(len(pl)):
                    for j_ in range(i_ + 1, len(pl)):
                        if pl[i_] > pl[j_]:
                            sgn = -sgn
                sign *= sgn
        total += sign * fun(assign)
    return total


def setup(fields):
    """fields: list of (name, shape, count) -> plane-wave legs with total momentum 0"""
    n = sum(c for _, _, c in fields)
    ps = momenta(n)
    out = {}; k = 0
    for name, shape, cnt in fields:
        out[name] = []
        for _ in range(cnt):
            out[name].append(jet(shape, kind='wave', p=ps[k], amp=rnd(*shape))); k += 1
    return out


def lap(j):
    return sum(j['dd'][m][m] for m in range(4))


# quadratic counterterm (thm c4-quadratic), each functional
L = setup([('A', (10, NC), 2)])
mass_IF = integ(lambda a: gp(a['A'][0]['v'], a['A'][1]['v']), L, ['A'])
mass_2C = integ(lambda a: gp(a['A'][0]['v'][:4], a['A'][1]['v'][:4]) + 0.25 * gp(bis_up(a['A'][0]['v'][4:]), bis_lo(a['A'][1]['v'][4:])), L, ['A'])
record('T1', 'int (A,A) = int (A_m,A_m) + (1/4) int (A^{ij},A_{ij})', mass_2C - mass_IF)


def dAdA_IF(a):
    F0 = F_IF(dict(v=np.zeros((10, NC)), d=a['A'][0]['d'])); F1 = F_IF(dict(v=np.zeros((10, NC)), d=a['A'][1]['d']))
    return hodge(F0, F1)


def dAdA_2C(a):
    d0 = a['A'][0]['d']; d1 = a['A'][1]['d']
    f0 = np.array([[d0[m][n] - d0[n][m] for n in range(4)] for m in range(4)])
    f1 = np.array([[d1[m][n] - d1[n][m] for n in range(4)] for m in range(4)])
    return 0.5 * gp(f0, f1) + 0.25 * gp(np.array([bis_up(d0[m][4:]) for m in range(4)]), np.array([bis_lo(d1[m][4:]) for m in range(4)]))


record('T2', 'int (dA,dA) = (1/2) int (d_m A_n - d_n A_m)^2 + (1/4) int (d_m A^{ij}, d_m A_{ij})', integ(dAdA_2C, L, ['A']) - integ(dAdA_IF, L, ['A']))
dstar_IF = integ(lambda a: gp(-sum(a['A'][0]['d'][m][m] for m in range(4)), -sum(a['A'][1]['d'][m][m] for m in range(4))), L, ['A'])
dstar_2C = integ(lambda a: gp(sum(a['A'][0]['d'][m][m] for m in range(4)), sum(a['A'][1]['d'][n][n] for n in range(4))), L, ['A'])
record('T3', 'int (d^*A|_V, d^*A|_V) = int (d_m A_m, d_n A_n)  (d^* = -d_m on one-forms; finite counterterm 27/16)', dstar_2C - dstar_IF)
L = setup([('c', (NC,), 1), ('As', (10, NC), 1)])
record('T4', 'int (dc, A^*) = int (d_m c, A^*_m)', integ(lambda a: gp(a['c'][0]['d'], a['As'][0]['v'][:4]), L, []) -
       integ(lambda a: ghost_IF(a['c'][0], zero_jet((10, NC)), a['As'][0]), L, []))
L = setup([('c', (NC,), 1), ('A', (10, NC), 1)])
obsIF = integ(lambda a: gp(-lap(a['c'][0]), -sum(a['A'][0]['d'][m][m] for m in range(4))), L, [])
record('T5', 'obstruction -(27/8) int (d^*d c, d^*A|_V) = -(27/8) int (d^2 c, d_m A_m)  (d^*d c = -d^2 c)',
       -27 / 8 * integ(lambda a: gp(lap(a['c'][0]), sum(a['A'][0]['d'][m][m] for m in range(4))), L, []) - (-27 / 8) * obsIF)
# field-counting functional
fj = dict(A=jet((10, NC)), As=jet((10, NC)), l=jet((16, NC)), ls=jet((16, NC)), c=jet((NC,)), cs=jet((NC,)))
fc_IF = 2 * gp(fj['A']['v'], fj['As']['v']) - 4 * gp(U @ fj['l']['v'], Us @ fj['ls']['v']) - gp(fj['cs']['v'], fj['c']['v'])
q, qb = split(fj['l']['v']); s, sb = split(fj['ls']['v'])
fc_2C = (2 * (gp(fj['A']['v'][:4], fj['As']['v'][:4]) + 0.25 * gp(bis_up(fj['A']['v'][4:]), bis_lo(fj['As']['v'][4:])))
         - 4 * (gp(q, s) - gp(qb, sb)) - gp(fj['cs']['v'], fj['c']['v']))
record('T6', 'field-counting functional 2(A,A^*) - 4(lambda^*,lambda) - (c^*,c) in two-component form', fc_2C - fc_IF)
# scalar bilinears
w1 = rnd(6); w2 = rnd(6); AI = rnd(6, NC)
O_IF = gp(np.einsum('I,Ic->c', w1, AI), np.einsum('I,Ic->c', w2, AI))
O_2C = (1 / 16) * gp(np.einsum('ij,ijc->c', bis_up(w1), bis_lo(AI)), np.einsum('ij,ijc->c', bis_up(w2), bis_lo(AI)))
record('T7', '(iota_w1 A, iota_w2 A) = (1/16)(w1^{ij} A_{ij}, w2^{kl} A_{kl});  (w1,w2) = (1/4) w1^{ij} w2_{ij}',
       [O_2C - O_IF, np.sum(w1 * w2) - 0.25 * np.einsum('ij,ij->', bis_up(w1), bis_lo(w2))])
K_IF = gp(AI, AI); K_2C = 0.25 * gp(bis_up(AI), bis_lo(AI))
record('T8', 'Konishi (trace) (A|perp,A|perp) = (1/4)(A^{ij},A_{ij})', K_2C - K_IF)
tl_IF = O_IF - np.sum(w1 * w2) / 6 * K_IF
tl_2C = O_2C - (0.25 * np.einsum('ij,ij->', bis_up(w1), bis_lo(w2))) / 6 * K_2C
# dICT for the bilinear (thm c5-Oeps), log part, translated
ww = np.einsum('ij,ij->', bis_up(w1), bis_lo(w2))
lit = 0.25 * gp(np.einsum('ij,ijc->c', bis_up(w1), bis_lo(AI)), np.einsum('ij,ijc->c', bis_up(w2), bis_lo(AI))) - 0.125 * ww * gp(bis_up(AI), bis_lo(AI))
record('T9', 'log part of d_delta I^CT (thm c5-Oeps): 4O - 2(w1,w2)(A|perp,A|perp) = (1/4)(w1^ij A_ij, w2^kl A_kl) - (1/8) w1^ij w2_ij (A^kl,A_kl)',
       lit - (4 * O_IF - 2 * np.sum(w1 * w2) * K_IF))
tl_lit = (1 / 16) * gp(np.einsum('ij,ijc->c', bis_up(w1), bis_lo(AI)), np.einsum('ij,ijc->c', bis_up(w2), bis_lo(AI))) - (1 / 96) * ww * gp(bis_up(AI), bis_lo(AI))
tr_lit = sum((1 / 16) * gp(np.einsum('ij,ijc->c', bis_up(np.eye(6)[I]), bis_lo(AI)), np.einsum('ij,ijc->c', bis_up(np.eye(6)[I]), bis_lo(AI)))
             - (1 / 96) * np.einsum('ij,ij->', bis_up(np.eye(6)[I]), bis_lo(np.eye(6)[I])) * gp(bis_up(AI), bis_lo(AI)) for I in range(6))
record('T10', '20prime: (1/16)(w1^ij A_ij, w2^kl A_kl) - (1/96) w1^ij w2_ij (A^kl,A_kl) = O - (1/6)(w1,w2) Konishi, and its trace over w1=w2=e_I vanishes',
       [tl_lit - tl_IF, tr_lit])
# -------- Step-2 bases G1..G5, H1..H6 (DESIGN_STEP2.md)
L = setup([('c', (NC,), 1), ('A', (10, NC), 2)])


def ddA_V(j):         # (d^*d A)_n = -d^2 A_n + d_n d_m A_m
    return np.array([-lap(j)[n] + sum(j['dd'][n][m][m] for m in range(4)) for n in range(4)])


def ddstarA_V(j):     # (d d^* A)_n = -d_n d_m A_m
    return np.array([-sum(j['dd'][n][m][m] for m in range(4)) for n in range(4)])


G1_IF = integ(lambda a: gp(br(a['c'][0]['v'], a['A'][0]['v'][:4]), ddA_V(a['A'][1])), L, ['A'])
G1_2C = integ(lambda a: gp(br(a['c'][0]['v'], a['A'][0]['v'][:4]),
                           np.array([-sum(a['A'][1]['dd'][m][m][n] - a['A'][1]['dd'][m][n][m] for m in range(4)) for n in range(4)])), L, ['A'])
record('H-G1', 'G1 = int([c,A_n], -d_m F_mn),  F_mn = d_m A_n - d_n A_m', G1_2C - G1_IF)
G2_IF = integ(lambda a: gp(br(a['c'][0]['v'], a['A'][0]['v'][:4]), ddstarA_V(a['A'][1])), L, ['A'])
G2_2C = integ(lambda a: -gp(br(a['c'][0]['v'], a['A'][0]['v'][:4]), np.array([sum(a['A'][1]['dd'][n][m][m] for m in range(4)) for n in range(4)])), L, ['A'])
record('H-G2', 'G2 = -int([c,A_n], d_n d_m A_m)', G2_2C - G2_IF)
G3_IF = integ(lambda a: gp(br(a['c'][0]['v'], a['A'][0]['v'][4:]), -lap(a['A'][1])[4:]), L, ['A'])
G3_2C = integ(lambda a: -0.25 * gp(br(a['c'][0]['v'], bis_up(a['A'][0]['v'][4:])), bis_lo(lap(a['A'][1])[4:])), L, ['A'])
record('H-G3', 'G3 = -(1/4) int([c,A^{ij}], d^2 A_{ij})', G3_2C - G3_IF)
L = setup([('c', (NC,), 2), ('As', (10, NC), 1)])
G4_IF = integ(lambda a: gp(br(a['c'][0]['v'], a['c'][1]['d']), a['As'][0]['v'][:4]), L, ['c'], odd_fields=['c'])
G4_2C = integ(lambda a: sum(gp(br(a['c'][0]['v'], a['c'][1]['d'][m]), a['As'][0]['v'][m]) for m in range(4)), L, ['c'], odd_fields=['c'])
record('H-G4', 'G4 = int([c, d_m c], A^*_m)', G4_2C - G4_IF)
L = setup([('c', (NC,), 1), ('l', (16, NC), 2)])
G5_IF = integ(lambda a: KAPPA * np.sum(br(a['c'][0]['v'], U @ a['l'][0]['v']) * sum(sA[m] @ (U @ a['l'][1]['d'][m]) for m in range(4))),
              L, ['l'], odd_fields=['l'])


def G5_2C_f(a):
    q0, qb0 = split(a['l'][0]['v'])
    dq1 = np.array([split(a['l'][1]['d'][m])[0] for m in range(4)]); dqb1 = np.array([split(a['l'][1]['d'][m])[1] for m in range(4)])
    c = a['c'][0]['v']
    return 1j * KAPPA * sum(np.einsum('aic,da,dic->', br(c, q0), SIGBE[m], dqb1[m]) + np.einsum('dic,da,aic->', br(c, qb0), SIGBE[m], dq1[m])
                            for m in range(4))


record('H-G5', 'G5 = i int{([c,lambda_a^i], sbar^{m d a} d_m lambar_{d i}) + ([c,lambar_{d i}], sbar^{m d a} d_m lambda_a^i)}',
       integ(G5_2C_f, L, ['l'], odd_fields=['l']) - G5_IF)
L = setup([('A', (10, NC), 3)])
H1_IF = integ(lambda a: hodge(F_IF(dict(v=np.zeros((10, NC)), d=a['A'][0]['d'])),
                              np.array([[br(a['A'][1]['v'][M], a['A'][2]['v'][N]) - br(a['A'][1]['v'][N], a['A'][2]['v'][M])
                                         if (M < 4 and N < 4) else np.zeros(NC) for N in range(10)] for M in range(10)])), L, ['A'])
H1_2C = integ(lambda a: 2 * sum(gp(a['A'][0]['d'][m][n], br(a['A'][1]['v'][m], a['A'][2]['v'][n])) for m in range(4) for n in range(4)), L, ['A'])
record('H-H1', 'H1 = int(dA|_V,[A|_V ^ A|_V]) = 2 int(d_m A_n, [A_m, A_n])', H1_2C - H1_IF)
H2_IF = integ(lambda a: sum(gp(a['A'][0]['d'][m][4 + I], br(a['A'][1]['v'][m], a['A'][2]['v'][4 + I])) for m in range(4) for I in range(6)), L, ['A'])
H2_2C = integ(lambda a: 0.25 * sum(gp(bis_up(a['A'][0]['d'][m][4:]), br(a['A'][1]['v'][m], bis_lo(a['A'][2]['v'][4:]))) for m in range(4)), L, ['A'])
record('H-H2', 'H2 (componentwise reading sum_{m,I}(d_m A_I,[A_m,A_I])) = (1/4) int(d_m A^{ij},[A_m, A_{ij}])', H2_2C - H2_IF)
hj = dict(A=jet((10, NC)), As=jet((10, NC)), c=jet((NC,)), l=jet((16, NC)), ls=jet((16, NC)), cs=jet((NC,)))
H3 = gp(br(hj['A']['v'][:4], hj['c']['v']), hj['As']['v'][:4])
H4_IF = gp(br(hj['A']['v'][4:], hj['c']['v']), hj['As']['v'][4:])
H4_2C = 0.25 * gp(br(bis_up(hj['A']['v'][4:]), hj['c']['v']), bis_lo(hj['As']['v'][4:]))
record('H-H4', 'H4 = int([A|perp,c],A^*|perp) = (1/4) int([A^{ij},c],A^*_{ij})', H4_2C - H4_IF)
q, qb = split(hj['l']['v']); s, sb = split(hj['ls']['v'])
H5_IF = gp(br(hj['c']['v'], U @ hj['l']['v']), Us @ hj['ls']['v'])
H5_2C = gp(br(hj['c']['v'], q), s) - gp(br(hj['c']['v'], qb), sb)
record('H-H5', 'H5 = int([c,lambda],lambda^*) = ([c,lambda_a^i],lambda^{*a}_i) - ([c,lambar_{d i}],lambar^{*d i})', H5_2C - H5_IF)


# ======================================================================================================
# 11. round trips
# ======================================================================================================
log('=== 11. round trips ===')
lam = rnd(16, NC)
record('RT1', 'spinor round trip: two-component -> classical realisation (U) -> two-component', Uinv @ (U @ lam) - lam)
p = rng.normal(size=4)
M_A = sum(1j * p[m] * sA[m] for m in range(4))               # bilinear form (l1, dslash l2) in the classical basis: l1^T M_A l2
M_D = U.T @ M_A @ U                                            # pulled back to D
M_2C = np.zeros((16, 16), complex)                             # literal two-component form: q1 (i p sbar qb2) - qb1 (-i sbar p q2)
for k in range(16):
    e = E16[:, k]; q2, qb2 = split(e)
    s_ = 1j * np.einsum('di,da->ai', qb2, sum(1j * p[m] * SIGBE[m] for m in range(4)))
    sb_ = -1j * np.einsum('da,ai->di', sum(1j * p[m] * SIGBE[m] for m in range(4)), q2)
    M_2C[:, k] = PD @ join(s_, sb_)
record('RT2', 'formula round trip (lambda, dslash lambda): index-free(classical) -> two-component -> index-free',
       [M_D - M_2C, Uinv.T @ M_2C @ Uinv - M_A])


# ======================================================================================================
# 12. component checks
# ======================================================================================================
log('=== 12. component checks ===')
ok_pos = True
for _ in range(5):
    p = rng.normal(size=4); AI = rng.normal(size=(6, NC))
    val = -(1 / 16) * KAPPA * np.dot(p, p) * np.sum(bis_up(AI) * bis_lo(AI))           # -(1/16)(d A^{ij}, d A_{ij}) on a real wave
    ok_pos &= (val.real < 0 and abs(val.imag) < 1e-12)
record('C1', 'real slice: scalar kinetic -(1/16)(d_m A^{ij}, d_m A_{ij}) = -(1/4)(dA_perp,dA_perp) < 0  [rem c1-real]', 0 if ok_pos else 1)
# SUSY closure on fields (abelian, free):  {Q_a.,Q_b.} A = -iota_v F,  Q.A = -Gamma(Q,lambda), Q.lambda = -1/2 rho(dA) Q
for _ in range(2):
    xa = rnd(16); xb = rnd(16)
    Aj = jet((10, NC)); Aj['v'] = np.zeros((10, NC), complex)          # value 0: brackets drop out (abelian, free)
    lam_a = -0.5 * join(*rhoF_2C(Aj, xa)); lam_b = -0.5 * join(*rhoF_2C(Aj, xb))
    Gm_b, Glo_b = Gam2C_field(xb, lam_a)
    Gm_a, Glo_a = Gam2C_field(xa, lam_b)
    anti_m = -(Gm_b + Gm_a); anti_lo = -(Glo_b + Glo_a)
    v = G2C(xa, xb)
    exp_m = np.array([-sum(v[n] * (Aj['d'][n][m] - Aj['d'][m][n]) for n in range(4)) +
                      sum(v[4 + I] * Aj['d'][m][4 + I] for I in range(6)) for m in range(4)])
    exp_lo = -sum(v[n] * bis_lo(Aj['d'][n][4:]) for n in range(4))
    record('C2', 'SUSY closure: {Q_a.,Q_b.}A_m = -v^n F_nm + v^I d_m A_I,  {Q_a.,Q_b.}A_ij = -v^n d_n A_ij,  v = Gamma(Q_a,Q_b)',
           [anti_m - exp_m, anti_lo - exp_lo])
# ESW Prop 2.1 in D: Gamma(Q1,rho(X)Q2) + Gamma(Q2,rho(X)Q1) = -2 iota_{Gamma(Q1,Q2)} X
Xf = rng.normal(size=(10, 10)); Xf = Xf - Xf.T
rhoX = sum(Xf[M, N] * RbD[M] @ RD[N] for M in range(10) for N in range(M + 1, 10))
Q1 = rnd(16); Q2 = rnd(16)
record('C3', '[ESW Prop 2.1] Gamma(Q1,rho(X)Q2)+Gamma(Q2,rho(X)Q1) = -2 iota_{Gamma(Q1,Q2)} X  (two-component Gamma)',
       G2C(Q1, rhoX @ Q2) + G2C(Q2, rhoX @ Q1) + 2 * np.einsum('M,MN->N', G2C(Q1, Q2), Xf))

# ------------------------------------------------------------------------------------------------------
log('=== 13. Srednicki-suppressed forms, Wick image of the Lagrangian, canonical normalisation ===')


def up_(x):      # psi^alpha = eps^{alpha beta} psi_beta  (S-SPN-3, (34.24)); first axis spinor
    return np.einsum('ab,b...->a...', EPSU, x)


xi1 = rnd(16); xi2 = rnd(16); eta = rnd(16)
q1, qb1 = split(xi1); q2, qb2 = split(xi2); s, sb = split(eta)
s_lo = np.einsum('ab,b...->a...', EPSD, s)                 # eta_{alpha i} = eps_{alpha beta} eta^beta_i
record('G7', '(xi,eta) = -(xi^i eta_i) - (xibar_i etabar^i)  [Srednicki contractions S-SPN-4]',
       pair_D(xi1, eta) - (-np.einsum('ai,ai->', up_(q1), s_lo) - np.einsum('di,di->', qb1, sb)))
qb2_up = up_(qb2)                                         # xibar^{alphadot}_i = eps^{alphadot betadot} xibar_{betadot i}
supp = (-2 * np.einsum('jikl,aj,ai->kl', EPS4, up_(q1), q2)
        - 2 * (np.einsum('dk,dl->kl', qb1, qb2_up) - np.einsum('dl,dk->kl', qb1, qb2_up)))
record('G8', 'Gamma_kl = -2 eps_jikl (xi1^j xi2^i) - 2[(xibar1_k xibar2_l) - (xibar1_l xibar2_k)]  (suppressed form)',
       Gbis_lo(xi1, xi2) - supp)
# Yukawa in suppressed form (pointwise, Grassmann polarised)
for _ in range(1):
    Aj = jet((10, NC)); l1 = jet((16, NC)); l2 = jet((16, NC))
    lo = bis_lo(Aj['v'][4:]); up = bis_up(Aj['v'][4:])

    def yuk_supp(a, b):
        qa, qba = split(a['v']); qb_, qbb = split(b['v'])
        t1 = -0.5 * KAPPA * np.einsum('ajc,jiac->', up_(qa), br(lo[:, :, None, :], np.transpose(qb_, (1, 0, 2))[None, :, :, :]))
        t2 = -0.5 * KAPPA * np.einsum('djc,jidc->', qba, br(up[:, :, None, :], np.transpose(up_(qbb), (1, 0, 2))[None, :, :, :]))
        return t1 + t2
    Z = dict(v=Aj['v'].copy(), d=np.zeros_like(Aj['d']), dd=np.zeros_like(Aj['dd']))
    Z['v'][:4] = 0
    l1z = dict(v=l1['v'], d=np.zeros_like(l1['d']), dd=l1['dd']); l2z = dict(v=l2['v'], d=np.zeros_like(l2['d']), dd=l2['dd'])
    record('L6', 'Yukawa = -(1/2)(lambda^j [A_ji, lambda^i]) - (1/2)(lambdabar_j [A^ji, lambdabar_i])  (suppressed, Grassmann)',
           (yuk_supp(l1, l2) - yuk_supp(l2, l1)) - (ferm_2C(l1z, l2z, Z) - ferm_2C(l2z, l1z, Z)))
# integrated covariant kinetic term with one gauge leg: int (1/2)(lambda, rho(e_m)[A_m,lambda]) = int i(lambar, sbar^m [A_m, lambda])
for _ in range(2):
    ps = momenta(3)
    l1, l2 = waves(ps[:2], [(16, NC), (16, NC)])
    Aw = jet((10, NC), kind='wave', p=ps[2], amp=rnd(10, NC)); Aw['v'][4:] = 0
    Ad = dict(v=Aw['v'], d=np.zeros_like(Aw['d']), dd=np.zeros_like(Aw['dd']))

    def cub_IF(a, b):
        return 0.5 * KAPPA * np.sum((U @ a['v']) * sum(sA[m] @ br(Ad['v'][m][None, :], U @ b['v']) for m in range(4)))

    def cub_2C(a, b):
        _, qba = split(a['v']); qb_, _ = split(b['v'])
        return 1j * KAPPA * sum(np.einsum('dic,da,aic->', qba, SIGBE[m], br(Ad['v'][m][None, None, :], qb_)) for m in range(4))
    record('L7', 'int (1/2)(lambda,dslash_A lambda)|_{A_m} = int i(lambar_{d i}, sbar^{m d a}[A_m,lambda_a^i])  (Grassmann)',
           (cub_2C(l1, l2) - cub_2C(l2, l1)) - (cub_IF(l1, l2) - cub_IF(l2, l1)))
# Wick image of -1/4 F_{mu nu}F^{mu nu}:  A_0 = i A_4, d_0 = i d_4  =>  F_{0k} = i F_{4k}
FE = rng.normal(size=(4, 4)); FE = FE - FE.T                       # Euclidean F_mn, index 0..3 = m 1..4
FL = np.zeros((4, 4), complex)
idx = {0: 3, 1: 0, 2: 1, 3: 2}                                     # Lorentzian mu -> Euclidean slot (mu=0 -> m=4)
for mu in range(4):
    for nu in range(4):
        FL[mu, nu] = FE[idx[mu], idx[nu]] * ((1j if mu == 0 else 1) * (1j if nu == 0 else 1))
FLup = ETA @ FL @ ETA
record('W8', 'Wick: F_{mu nu}F^{mu nu} = F_mn F_mn (A_0 = i A_4, d_0 = i d_4); e^{iS_L} = e^{S_E} with dt = -i d tau',
       np.sum(FL * FLup) - np.sum(FE * FE))
# canonical normalisation for su(2) in the defining representation (rem c6-konishi)
T = [0.5 * P for P in PAULI]
gcoup = 0.7; hb = 1.3


def pair_su(X, Y):
    return -(4 * hb / gcoup ** 2) * np.trace(X @ Y)


Xa = rng.normal(size=3); Ya = rng.normal(size=3)
Xm = -1j * gcoup * sum(Xa[a] * T[a] for a in range(3)); Ym = -1j * gcoup * sum(Ya[a] * T[a] for a in range(3))
br_m = Xm @ Ym - Ym @ Xm
fS = np.cross(Xa, Ya)                                               # f^{abc} = eps^{abc} for T = sigma/2 (S-GRP-1/2)
record('C4', 'canonical map X = -i g X^a T^a: (X,Y) = 2 hbar X^a Y^a and [X,Y] = -i g (g eps^{abc} X^b Y^c) T^a  [S-GRP-1,2,5,6]',
       [pair_su(Xm, Ym) - 2 * hb * np.dot(Xa, Ya), br_m - (-1j * gcoup) * gcoup * sum(fS[a] * T[a] for a in range(3))])
lamA = rng.normal(size=3); lamB = rng.normal(size=3)
Lm = -1j * gcoup / np.sqrt(2) * sum(lamA[a] * T[a] for a in range(3)); Lm2 = -1j * gcoup / np.sqrt(2) * sum(lamB[a] * T[a] for a in range(3))
record('C5', 'gaugino lambda = -(i g/sqrt 2) lambda^a T^a: (lambda,lambda\') = hbar lambda^a lambda\'^a  =>  (1/2)(lambda,dslash lambda) = i hbar lambdabar^a sbar^m d_m lambda^a',
       pair_su(Lm, Lm2) - hb * np.dot(lamA, lamB))
record('C6', 'Konishi: -24 * (-g^2 N/(32 pi^2)) = 3 g^2 N/(4 pi^2)  [rem c6-konishi]', -24 * (-1 / 32) - 3 / 4)

# ------------------------------------------------------------------------------------------------------
log('=== 14. negative controls (perturbed dictionary entries must fail) ===')


def neg(cid, desc, res):
    record(cid, 'negative control: ' + desc + f'  (residual of perturbed entry {res:.2e} > 1e-3)', 0.0 if res > 1e-3 else 1.0)


PDalt = PD.copy(); PDalt[8:, 8:] *= -1
neg('N1', 'pairing with +qbar.sbar breaks (rho(v)Q1,Q2) = (Q1,rho(v)Q2)',
    max(np.max(np.abs(PDalt @ RD[M] - (PDalt @ RD[M]).T)) for M in range(10)))
p1, p2 = momenta(2)
l1, l2 = waves([p1, p2], [(16, NC), (16, NC)])
Z = zero_jet((10, NC))
vIF = ferm_IF(mapS(U, l1), mapS(U, l2), Z) - ferm_IF(mapS(U, l2), mapS(U, l1), Z)


def kin_wrong(a, b):
    _, qba = split(a['v']); dq = np.array([split(b['d'][m])[0] for m in range(4)])
    return -1j * KAPPA * sum(np.einsum('dic,da,aic->', qba, SIGBE[m], dq[m]) for m in range(4))


neg('N2', 'kinetic term -i lambar sbar d lambda', abs((kin_wrong(l1, l2) - kin_wrong(l2, l1)) - vIF))
xi = rnd(16); x = rng.normal(size=4); eta = rnd(16)
neg('N3', 'Killing spinor with +i x sig etabar', np.max(np.abs(U @ (-killing_2C(eta, x) + 2 * join(split(killing_2C(eta, x))[0] * 0, split(killing_2C(eta, x))[1])) -
                                                       sum(x[m] * sbA[m] for m in range(4)) @ (Us @ eta))))
neg('N4', 'R-rotation coefficient 1 instead of 2', np.max(np.abs(0.5 * rot_2C(xi1, eta2, xi3, eta4) - rmat(rot_IF(xi1, eta2, xi3, eta4)))))
p = rng.normal(size=4)
neg('N5', 'lambda propagator +p.sigma', np.max(np.abs(1j * sum(p[m] * RbD[m] for m in range(4)) @ np.linalg.inv(PD) +
                                                        (1j * sum(p[m] * RbD[m] for m in range(4)) @ np.linalg.inv(PD)))))
Aj = jet((10, NC)); lj = jet((16, NC)); Asj = jet((10, NC)); lsj = jet((16, NC))
neg('N6', 'S^(1)(Q_1) with +(Gamma,A^*)', abs((S1Q_2C(xi, Aj, lj, Asj, lsj) + 2 * gp(Gam2C_field(xi, lj['v'])[0], Asj['v'][:4])) -
                                             S1Q_IF(U @ xi, Aj, mapS(U, lj), Asj, mapS(Us, lsj))))


# ======================================================================================================
# dictionary.json: entries (as typeset in esw_to_weinberg_srednicki.tex) with the residuals of their checks
# ======================================================================================================
ENTRIES = [
 # id, category, esw, ws, source, standard, tag, checks
 ("SPN-eps", "base", "symplectic pairings on S_+, S_- [ESW §10]", r"\varepsilon^{12}=\varepsilon^{\dot1\dot2}=+1,\ \varepsilon_{12}=\varepsilon_{\dot1\dot2}=-1;\ \psi^\alpha=\varepsilon^{\alpha\beta}\psi_\beta,\ \psi_\alpha=\varepsilon_{\alpha\beta}\psi^\beta;\ \chi\psi=\chi^\alpha\psi_\alpha,\ \bar\chi\bar\psi=\bar\chi_{\dot\alpha}\bar\psi^{\dot\alpha}", "Srednicki (35.1),(34.24),(34.25),(35.23)", "S-SPN-2,S-SPN-3,S-SPN-4", "EVIDENCED", ["W1","G7"]),
 ("SPN-sigmaE", "Wick", "Clifford action of V on S_+ (+) S_- (ch01 Def c1-spaces)", r"v^0=-iv^4\ (x^4=it,\ \partial_0=i\partial_4);\ \sigma_E^m=(\vec\sigma,i),\ \bar\sigma_E^m=(-\vec\sigma,i);\ \sigma^\mu v_\mu=\sigma_E^mv_m,\ \bar\sigma^\mu v_\mu=\bar\sigma_E^mv_m;\ \sigma_E^m\bar\sigma_E^n+\sigma_E^n\bar\sigma_E^m=-2\delta^{mn};\ \bar\sigma_E=-\sigma_E^\dagger", "Srednicki (34.28),(35.2),(35.3),(36.8)", "S-SPN-5,S-MET-1,S-MET-5", "DERIVED", ["W2","W3","W4"]),
 ("SU4-blocks", "base", "W, W^*, W_R^perp, rho(w) (ESW §10.3, W = S_+^6)", r"W\ni w^i,\ W^*\ni w_i;\ \Sigma^I_{ij}=(\eta^a,i\bar\eta^a),\ \bar\Sigma^I=(\Sigma^I)^\dagger;\ \tfrac12\varepsilon^{ijkl}\Sigma^I_{kl}=-\bar\Sigma^{Iij};\ \bar\Sigma^I\Sigma^J+\bar\Sigma^J\Sigma^I=2\delta^{IJ};\ \sum_I\Sigma^I_{ij}\Sigma^I_{kl}=2\varepsilon_{ijkl},\ \sum_I\bar\Sigma^{Iij}\bar\Sigma^{Ikl}=2\varepsilon^{ijkl},\ \sum_I\bar\Sigma^{Iij}\Sigma^I_{kl}=-2(\delta^i_k\delta^j_l-\delta^i_l\delta^j_k)", "engine_b/realization.py thooft, Sig, Sigb", "", "DERIVED", ["W5","W6","W7"]),
 ("rho", "realisation", "rho(v), rho(w) on Sigma, Sigma^* (ch01 Def c1-spaces, Prop c1-chirality)", r"(\rho(v)\xi)^\alpha{}_i=iv_m\bar\xi_{\dot\alpha i}\bar\sigma^{m\dot\alpha\alpha},\ (\rho(v)\xi)^{\dot\alpha i}=-iv_m\bar\sigma^{m\dot\alpha\alpha}\xi_\alpha{}^i;\ (\rho(v)\eta)_\alpha{}^i=-iv_m\sigma^m_{\alpha\dot\alpha}\bar\eta^{\dot\alpha i},\ (\rho(v)\eta)_{\dot\alpha i}=iv_m\eta^\alpha{}_i\sigma^m_{\alpha\dot\alpha};\ (\rho(w)\xi)^\alpha{}_j=w_{ji}\xi^{\alpha i},\ (\rho(w)\xi)^{\dot\alpha j}=w^{ji}\bar\xi^{\dot\alpha}{}_i;\ (\rho(w)\eta)_\alpha{}^k=-w^{kj}\eta_{\alpha j},\ (\rho(w)\eta)_{\dot\beta k}=-w_{kj}\bar\eta_{\dot\beta}{}^j", "ch01 Def c1-spaces; engine_b/realization.py; classical/n4_model.py", "S-SPN-3,S-SPN-5,S-SPN-8", "DERIVED", ["R1","R5","B1","A1","A2"]),
 ("pairing", "realisation", "(-,-): Sigma (x) Sigma^* -> C (ch01 Def c1-spaces, [ESW (28)])", r"(\xi,\eta)=\xi_\alpha{}^i\eta^\alpha{}_i-\bar\xi_{\dot\alpha i}\bar\eta^{\dot\alpha i}=-(\xi^i\eta_i)-(\bar\xi_i\bar\eta^i)", "ch01 Def c1-spaces; engine B pairing (factor 1)", "S-SPN-4", "DERIVED", ["R2","B2","A2","G7","N1"]),
 ("chirality", "realisation", "S_+ (x) W = +1 eigenspace of rho(v1 v2 v3 v4) (ch01 Prop c1-chirality)", r"\text{orientation }dx^1\wedge dx^2\wedge dx^3\wedge dx^4;\ S_+=\text{undotted}=(2,1);\ \rho(v)|_{(\psi_\alpha,\bar\chi^{\dot\alpha})}=-i\gamma_S^\mu v_\mu=\gamma_W^\mu v_\mu;\ \gamma_E^1\gamma_E^2\gamma_E^3\gamma_E^4=-\gamma_{5S}", "ch01 Prop c1-chirality; Srednicki (36.7),(36.44); Weinberg D1,D2", "S-SPN-7,S-SPN-8", "DERIVED", ["R4","R6","R7","A3"]),
 ("Gamma", "realisation", "Gamma: Sym^2 Sigma -> V_10, Sym^2 Sigma^* -> V_10 ([ESW (28)])", r"\Gamma^m(\xi_1,\xi_2)=i(\bar\xi_2\bar\sigma^m\xi_1+\bar\xi_1\bar\sigma^m\xi_2);\ \Gamma_{kl}=-2\varepsilon_{jikl}(\xi_1{}^j\xi_2{}^i)-2(\bar\xi_{1k}\bar\xi_{2l})+2(\bar\xi_{1l}\bar\xi_{2k});\ \Gamma^m(\eta_1,\eta_2)=-i(\eta_2\sigma^m\bar\eta_1+\eta_1\sigma^m\bar\eta_2);\ \Gamma^{kl}(\eta_1,\eta_2)=-2\varepsilon^{klpj}\varepsilon_{\alpha\beta}\eta_1{}^\beta{}_j\eta_2{}^\alpha{}_p+2(\varepsilon_{\dot\beta\dot\alpha}\bar\eta_1^{\dot\alpha l}\bar\eta_2^{\dot\beta k}-(k\leftrightarrow l))", "ch01 Def c1-spaces, [ESW (28), Assumption 3.1, Prop 2.1]", "S-CONJ-2 index structure", "DERIVED", ["R3","G1","G3","G4","G5","G8","C3"]),
 ("square-zero", "realisation", "Gamma(Q,Q)=0, Q in Sigma", r"\xi_\alpha{}^i=u_\alpha w^i,\ \bar\xi=0\Rightarrow\Gamma(\xi,\xi)=0", "[ESW Def 2.16]", "", "DERIVED", ["G6"]),
 ("fields", "fields", "c, A|_V, iota_w A, lambda, A^*, lambda^*, c^* (ch01 Def c1-fields)", r"c;\ A_m;\ A_{ij}=A_I\Sigma^I_{ij},\ A^{ij}=-A_I\bar\Sigma^{Iij}=\tfrac12\varepsilon^{ijkl}A_{kl},\ \iota_wA=\tfrac14w^{ij}A_{ij};\ \lambda_\alpha{}^i,\ \bar\lambda_{\dot\alpha i};\ A^*_m,\ A^*_{ij},\ A^{*ij};\ \lambda^{*\alpha}{}_i,\ \bar\lambda^{*\dot\alpha i};\ c^*;\ \text{real slice: c.c.}(A_{ij})=A^{ij}", "ch01 Def c1-fields, rem c1-real", "S-SPN-1,S-CONJ-1,S-CONJ-7", "DERIVED", ["S1","S2","S3","S4"]),
 ("parameters", "fields", "Q_1 in Sigma, Q_2 in Sigma^* (ch02 Def c2-odd)", r"Q_1\leftrightarrow(\xi_\alpha{}^i\ [Q^i_\alpha,(3,1)],\ \bar\xi_{\dot\alpha i}\ [\bar Q_{\dot\alpha i},(2,3)]);\ Q_2\leftrightarrow(\eta^\alpha{}_i\ [S^\alpha{}_i,(1,3)],\ \bar\eta^{\dot\alpha i}\ [\bar S^{\dot\alpha i},(3,2)])", "ch02 table of named subspaces; [SW (17)]", "", "DERIVED", ["K1","K2"]),
 ("pairings-fields", "fields", "(A,A^*), (lambda,lambda^*), (c^*,c)", r"(A,A^*)=(A_m,A^*_m)+\tfrac14(A^{ij},A^*_{ij});\ (\lambda,\lambda^*)=(\lambda_\alpha{}^i,\lambda^{*\alpha}{}_i)-(\bar\lambda_{\dot\alpha i},\bar\lambda^{*\dot\alpha i})", "ch01 Def c1-pairing", "", "DERIVED", ["T6","S2"]),
 ("action", "action", "S (ch01 eq:c1-action)", r"\int\{-\tfrac18(F_{mn},F_{mn})-\tfrac1{16}(D_mA^{ij},D_mA_{ij})-\tfrac1{128}([A^{ij},A^{kl}],[A_{ij},A_{kl}])+\tfrac i2(\lambda_\alpha{}^i,\bar\sigma^{m\dot\alpha\alpha}D_m\bar\lambda_{\dot\alpha i})+\tfrac i2(\bar\lambda_{\dot\alpha i},\bar\sigma^{m\dot\alpha\alpha}D_m\lambda_\alpha{}^i)+\tfrac12(\lambda_\alpha{}^j,\varepsilon^{\alpha\beta}[A_{ji},\lambda_\beta{}^i])-\tfrac12(\bar\lambda_{\dot\alpha j},\varepsilon^{\dot\alpha\dot\beta}[A^{ji},\bar\lambda_{\dot\beta i}])+(D_mc,A^*_m)+\tfrac14([A^{ij},c],A^*_{ij})+([c,\lambda_\alpha{}^i],\lambda^{*\alpha}{}_i)-([c,\bar\lambda_{\dot\alpha i}],\bar\lambda^{*\dot\alpha i})+\tfrac12([c,c],c^*)\}", "ch01 eq:c1-action", "S-ACT-9,S-ACT-12", "DERIVED", ["L1","L2","L4","L5","S5"]),
 ("action-integrated", "action", "int 1/2 (lambda, dslash_A lambda)", r"\int\{i(\bar\lambda_i\bar\sigma^mD_m\lambda^i)-\tfrac12(\lambda^j[A_{ji},\lambda^i])-\tfrac12(\bar\lambda_j[A^{ji},\bar\lambda_i])\}", "ch01 eq:c1-action", "S-ACT-12 (95.77)", "DERIVED", ["L3","L6","L7","N2"]),
 ("wick-canonical", "action", "e^{S/hbar}, su(N), (X,Y)=-(4 hbar/g^2)tr(XY) (ch06 rem c6-konishi)", r"A_m=-igA^a_mT^a,\ A_I=-igA^a_IT^a,\ \lambda=-\tfrac{ig}{\sqrt2}\lambda^aT^a;\ (X,Y)=2\hbar X^aY^a;\ S/\hbar\supset\int\{-\tfrac14F^aF^a-\tfrac12(D_mA_I)^2+i\bar\lambda^a\bar\sigma^mD_m\lambda^a\}=\int\mathcal L_{\rm Sred}|_{x^0=-ix^4}", "ch06 rem c6-konishi", "S-ACT-7,S-ACT-9,S-ACT-12,S-GRP-2,S-GRP-5,S-GRP-6; Weinberg D9", "DERIVED", ["W8","C4","C5"]),
 ("Q", "BV", "Q (ch01 Prop c1-Q)", r"c\mapsto-\partial_mc\,[A_m];\ A_m\mapsto\tfrac12(\partial^2A_m-\partial_m\partial_nA_n)\,[A^*_m];\ A_{ij}\mapsto\tfrac12\partial^2A_{ij}\,[A^*_{ij}];\ A^*_m\mapsto-\partial_mA^*_m\,[c^*];\ \lambda_\alpha{}^i\mapsto-i\bar\sigma^{m\dot\alpha\alpha}\partial_m\lambda_\alpha{}^i\,[\bar\lambda^{*\dot\alpha i}];\ \bar\lambda_{\dot\alpha i}\mapsto i\partial_m\bar\lambda_{\dot\alpha i}\bar\sigma^{m\dot\alpha\alpha}\,[\lambda^{*\alpha}{}_i]", "ch01 Prop c1-Q", "", "DERIVED", ["Q1","Q4"]),
 ("QGF", "BV", "Q^GF (ch03 Def c3-qgf)", r"A^*\mapsto-2A^*\,[A];\ A_m\mapsto\partial_mA_m\,[c];\ c^*\mapsto\partial_mc^*\,[A^*_m];\ \lambda^{*\alpha}{}_i\mapsto-i\partial_m\lambda^{*\alpha}{}_i\sigma^m_{\alpha\dot\alpha}\,[\bar\lambda_{\dot\alpha i}];\ \bar\lambda^{*\dot\alpha i}\mapsto i\sigma^m_{\alpha\dot\alpha}\partial_m\bar\lambda^{*\dot\alpha i}\,[\lambda_\alpha{}^i];\ H=-\partial^2", "ch03 Def c3-qgf, Prop c3-H", "", "DERIVED", ["Q2","Q3","Q4"]),
 ("propagator", "BV", "P(eps,infty) (ch03 Prop c3-P)", r"\times e^{-\varepsilon|p|^2}/|p|^2:\ \langle A_mA_n\rangle=2(-,-)^{-1}\delta_{mn},\ \langle A_{ij}A_{kl}\rangle=4(-,-)^{-1}\varepsilon_{ijkl},\ \langle A_{ij}A^{kl}\rangle=4(-,-)^{-1}(\delta_i^k\delta_j^l-\delta_i^l\delta_j^k),\ \langle A^{ij}A^{kl}\rangle=4(-,-)^{-1}\varepsilon^{ijkl},\ \langle cA^*_m\rangle=\langle A^*_mc\rangle=-(-,-)^{-1}ip_m,\ \langle\lambda_\alpha{}^i\bar\lambda_{\dot\beta j}\rangle=-(-,-)^{-1}p_m\sigma^m_{\alpha\dot\beta}\delta^i_j,\ \langle\bar\lambda_{\dot\alpha i}\lambda_\beta{}^j\rangle=-(-,-)^{-1}p_m\sigma^m_{\beta\dot\alpha}\delta_i^j", "ch03 Prop c3-P", "", "DERIVED", ["P1","P2","P3","N5"]),
 ("generators-even", "psl44", "even part (ch02 Def c2-even)", r"(v^m,X_{mn},\mathrm{div}\,v,r):\ P:(b^m,0,0,0);\ M:(\omega_{mn}x_n,-2\omega_{mn},0,0);\ \Delta:(2x^m,0,8,0);\ K:(x^2b^m-2(b\cdot x)x^m,4(x_mb_n-x_nb_m),-8b\cdot x,0);\ R:(0,0,0,\tfrac14R_{IJ}\bar\Sigma^J\Sigma^I)", "ch02 Def c2-even, rem c2-dilatation", "", "DERIVED", ["E-P","E-M","E-D","E-K","E-R"]),
 ("killing", "psl44", "x -> Q_1 + rho(x) Q_2 (ch02 Def c2-odd)", r"\epsilon_\alpha{}^i=\xi_\alpha{}^i-ix_m\sigma^m_{\alpha\dot\alpha}\bar\eta^{\dot\alpha i},\ \bar\epsilon_{\dot\alpha i}=\bar\xi_{\dot\alpha i}+ix_m\eta^\alpha{}_i\sigma^m_{\alpha\dot\alpha}", "ch02 Def c2-odd", "", "DERIVED", ["O2","N3"]),
 ("brackets", "psl44", "brackets (ch02 eq:c2-evenodd, eq:c2-oddodd, Prop c2-special)", r"[\epsilon_1,\epsilon_2]=i(\bar\epsilon_1\bar\sigma^m\epsilon_2+\bar\epsilon_2\bar\sigma^m\epsilon_1)\partial_m+r,\ r^k{}_i=2(\xi_{1\alpha}{}^k\eta_2{}^\alpha{}_i+\bar\xi_{1\dot\alpha i}\bar\eta_2^{\dot\alpha k})+2(1\leftrightarrow2)-\mathrm{tr};\ [v,\epsilon]_\alpha{}^i=v\cdot\partial\epsilon-\tfrac18X_{mn}\sigma^m\bar\sigma^n\epsilon-\tfrac18\mathrm{div}v\,\epsilon+r\epsilon,\ [v,\epsilon]_{\dot\alpha i}=v\cdot\partial\bar\epsilon+\tfrac18X_{mn}\bar\epsilon\bar\sigma^m\sigma^n-\tfrac18\mathrm{div}v\,\bar\epsilon-\bar\epsilon r", "ch02 eq:c2-evenodd, eq:c2-oddodd, Prop c2-special", "", "DERIVED", ["K1","K2","K3","K4-M","K4-D","K4-K","K4-R","N4"]),
 ("S1-even", "psl44", "S^(1)(v) (ch02 eq:c2-S1even, ch07 Def c7-vdot)", r"\int\{(v.A_m,A^*_m)+\tfrac14(v.A_{ij},A^{*ij})-(v.\lambda_\alpha{}^i,\lambda^{*\alpha}{}_i)+(v.\bar\lambda_{\dot\alpha i},\bar\lambda^{*\dot\alpha i})-(v.c,c^*)\};\ v.A_m=v^n\partial_nA_m+A_n\partial_mv^n,\ v.A_{ij}=v\cdot\partial A_{ij}+\tfrac14\mathrm{div}v\,A_{ij}-r^k{}_iA_{kj}-r^k{}_jA_{ik},\ v.\lambda_\alpha{}^i=v\cdot\partial\lambda-\tfrac18X_{mn}\sigma^m\bar\sigma^n\lambda+\tfrac38\mathrm{div}v\,\lambda+r\lambda,\ v.\bar\lambda_{\dot\alpha i}=v\cdot\partial\bar\lambda+\tfrac18X_{mn}\bar\lambda\bar\sigma^m\sigma^n+\tfrac38\mathrm{div}v\,\bar\lambda-\bar\lambda r,\ v.c=v\cdot\partial c", "ch02 eq:c2-S1even, ex:c2-dil; ch07 Def c7-vdot", "", "DERIVED", ["E-P","E-M","E-D","E-K","E-R"]),
 ("S1-Q", "psl44", "S^(1)(Q_1) (ch02 eq:c2-S1Q)", r"\int\{-(\Gamma^m(\xi,\lambda),A^*_m)-\tfrac14(\Gamma_{ij}(\xi,\lambda),A^{*ij})+\tfrac12((\rho(F_A)\xi)_\alpha{}^i,\lambda^{*\alpha}{}_i)-\tfrac12((\rho(F_A)\xi)_{\dot\alpha i},\bar\lambda^{*\dot\alpha i})\};\ (\rho(F_A)\xi)_\alpha{}^i=-\tfrac12F_{mn}\sigma^m\bar\sigma^n\xi-\tfrac12[A^{ik},A_{kj}]\xi_\alpha{}^j-i\sigma^m_{\alpha\dot\alpha}\varepsilon^{\dot\alpha\dot\beta}D_mA^{ij}\bar\xi_{\dot\beta j},\ (\rho(F_A)\xi)_{\dot\alpha i}=\tfrac12F_{mn}\bar\xi\bar\sigma^m\sigma^n-\tfrac12[A_{ik},A^{kj}]\bar\xi_{\dot\alpha j}+i\varepsilon^{\alpha\beta}\sigma^m_{\alpha\dot\alpha}D_mA_{ij}\xi_\beta{}^j", "ch02 eq:c2-S1Q", "", "DERIVED", ["O1","N6","C2"]),
 ("S1-S", "psl44", "S^(1)(Q_2) (ch02 eq:c2-S1S)", r"S^{(1)}(\xi)|_{\xi\to\rho(x)\eta}+\int\{(A^{kj}\varepsilon_{\alpha\beta}\eta^\beta{}_j,\lambda^{*\alpha}{}_k)-(A_{kj}\varepsilon_{\dot\beta\dot\alpha}\bar\eta^{\dot\alpha j},\bar\lambda^{*\dot\beta k})\}", "ch02 eq:c2-S1S", "", "DERIVED", ["O3"]),
 ("S2", "psl44", "S^(2)(Q_1,Q_2) (ch02 eq:c2-S2)", r"\int\{\tfrac14[-2i\Gamma^m(\epsilon_1,\epsilon_2)(\lambda^{*\alpha}{}_i,\sigma^m_{\alpha\dot\alpha}\bar\lambda^{*\dot\alpha i})+\tfrac14\Gamma_{kl}(\epsilon_1,\epsilon_2)\Gamma^{kl}(\lambda^*,\lambda^*)]-\tfrac12((\epsilon_1,\lambda^*),(\epsilon_2,\lambda^*))-(\Gamma^m A_m+\tfrac14\Gamma_{ij}A^{ij},c^*)\}", "ch02 eq:c2-S2", "", "DERIVED", ["O4-QQ","O4-QS","O4-SS"]),
 ("CT-quadratic", "step1", "I^CT quadratic (ch04 thm:c4-quadratic)", r"\mathrm{unit}\{\tfrac94\varepsilon^{-1}\int[(A_m,A_m)+\tfrac14(A^{ij},A_{ij})]+\log\varepsilon[-\tfrac12\int(F^0_{mn},F^0_{mn})-\tfrac14\int(\partial_mA^{ij},\partial_mA_{ij})+8i\int(\bar\lambda_{\dot\alpha i},\bar\sigma^{m\dot\alpha\alpha}\partial_m\lambda_\alpha{}^i)-\int(\partial_mc,A^*_m)]\}", "ch04 thm:c4-quadratic", "", "DERIVED", ["T1","T2","T4","L3"]),
 ("CT-log", "step1", "log part (ch04 thm:c4-log); field-counting functional (ch06)", r"2\int[(A_m,A^*_m)+\tfrac14(A^{ij},A^*_{ij})]-4\int[(\lambda_\alpha{}^i,\lambda^{*\alpha}{}_i)-(\bar\lambda_{\dot\alpha i},\bar\lambda^{*\dot\alpha i})]-\int(c^*,c)", "ch04 thm:c4-log; ch06 lem:c6-count", "", "DERIVED", ["T6"]),
 ("finite-CT", "step1", "obstruction and finite counterterm (ch04 eq:c4-obstruction, eq:c4-finite)", r"-\tfrac{27}8\mathrm{unit}\int(\partial^2c,\partial_mA_m)=\tfrac{27}{16}\mathrm{unit}\,Q\int(\partial_mA_m,\partial_nA_n);\ I\mapsto I-\tfrac{27}{16}\mathrm{unit}\int(\partial_mA_m,\partial_nA_n);\ \text{mass variation }-\tfrac92\mathrm{unit}\,\varepsilon^{-1}\int(\partial_mc,A_m)", "ch04 thm:c4-obstruction, eq:c4-finite, rem:c4-epsinv", "", "DERIVED", ["T3","T5"]),
 ("bilinears", "step1", "scalar bilinears, 20', Konishi (ch05 def:c5-bilinear, thm:c5-Oeps)", r"\int\omega(\iota_{w_1}A,\iota_{w_2}A)=\tfrac1{16}\int\omega(w_1^{ij}A_{ij},w_2^{kl}A_{kl});\ (w_1,w_2)=\tfrac14w_1^{ij}w_{2ij};\ \mathrm{Konishi}=\tfrac14\int\omega(A^{ij},A_{ij});\ \mathbf{20'}=\tfrac1{16}\int\omega(w_1A,w_2A)-\tfrac1{96}w_1^{ij}w_{2ij}\int\omega(A^{kl},A_{kl});\ \partial_\delta I^{CT}=\mathrm{unit}\log\varepsilon[\tfrac14\int\omega(w_1A,w_2A)-\tfrac18w_1^{ij}w_{2ij}\int\omega(A^{kl},A_{kl})]+\hbar\dim\mathfrak g\,w_1^{ij}w_{2ij}/(32\pi^2\varepsilon)\int\omega", "ch05 def:c5-bilinear, thm:c5-Oeps", "", "DERIVED", ["T7","T8","T9","T10"]),
 ("dilatation", "step1", "dilatation coefficients (ch06 thm:c6-dil, cor:c6-dilcomp, rem:c6-konishi)", r"\Delta=2E;\ \mathbf{20'}:0;\ \mathrm{Konishi}:-24;\ (X,A_{ij}):-4;\ (X,A_m):-4;\ (\zeta,\lambda)=(\lambda_\alpha{}^i,\zeta^\alpha{}_i)-(\bar\lambda_{\dot\alpha i},\bar\zeta^{\dot\alpha i}):-8;\ (X,c):-2\ (\mathrm{units});\ \mathfrak{su}(N):\ -24\,\mathrm{unit}=3g^2N/(4\pi^2)", "ch06 thm:c6-dil, cor:c6-dilcomp, rem:c6-konishi", "", "EVIDENCED", ["C6"]),
 ("G1-G5", "step2", "G1..G5 (DESIGN_STEP2.md Part C)", r"G_1=-\int([c,A_n],\partial_mF^0_{mn}),\ G_2=-\int([c,A_n],\partial_n\partial_mA_m),\ G_3=-\tfrac14\int([c,A^{ij}],\partial^2A_{ij}),\ G_4=\int([c,\partial_mc],A^*_m),\ G_5=i\int\{([c,\lambda_\alpha{}^i],\bar\sigma^{m\dot\alpha\alpha}\partial_m\bar\lambda_{\dot\alpha i})+([c,\bar\lambda_{\dot\alpha i}],\bar\sigma^{m\dot\alpha\alpha}\partial_m\lambda_\alpha{}^i)\}", "DESIGN_STEP2.md Part C", "", "DERIVED", ["H-G1","H-G2","H-G3","H-G4","H-G5"]),
 ("H1-H6", "step2", "H1..H6 (DESIGN_STEP2.md Part C)", r"H_1=2\int(\partial_mA_n,[A_m,A_n]),\ H_2=\tfrac14\int(\partial_mA^{ij},[A_m,A_{ij}])\ (\text{componentwise reading; two-form reading}=2H_2),\ H_3=\int([A_m,c],A^*_m),\ H_4=\tfrac14\int([A^{ij},c],A^*_{ij}),\ H_5=\int([c,\lambda_\alpha{}^i],\lambda^{*\alpha}{}_i)-\int([c,\bar\lambda_{\dot\alpha i}],\bar\lambda^{*\dot\alpha i}),\ H_6=\int([c,c],c^*)", "DESIGN_STEP2.md Part C", "", "DERIVED", ["H-H1","H-H2","H-H4","H-H5"]),
 ("beta-directions", "step2", "Q-closed quadratic finite counterterms (DESIGN_STEP2.md)", r"\beta_1\int(dA|_V,dA|_V)=\tfrac{\beta_1}2\int(F^0_{mn},F^0_{mn}),\ \beta_2\int(dA|_\perp,dA|_\perp)=\tfrac{\beta_2}4\int(\partial_mA^{ij},\partial_mA_{ij}),\ \beta_3\int(\lambda,\not d\lambda)=2i\beta_3\int(\bar\lambda_{\dot\alpha i},\bar\sigma^{m\dot\alpha\alpha}\partial_m\lambda_\alpha{}^i),\ \beta_4\int(\partial_mc,A^*_m)", "DESIGN_STEP2.md Part C", "", "DERIVED", ["T2","L3","T4"]),
 ("closure", "component", "{S^(1)(Q_1),S^(1)(Q_2)} (ch02 rem physical reading)", r"\{\xi_a.,\xi_b.\}A_m=-v^nF_{nm}+v^I\partial_mA_I,\ \{\xi_a.,\xi_b.\}A_{ij}=-v^n\partial_nA_{ij},\ v=\Gamma(\xi_a,\xi_b)", "[ESW Prop 2.1]", "", "DERIVED", ["C2","C3"]),
 ("weinberg", "Weinberg", "Weinberg column", r"\rho(v)=\gamma_W^\mu v_\mu|_{\rm Wick};\ \gamma_{5W}=-\gamma_{5S};\ e_{ab}=\varepsilon^{ab};\ \mathcal Q_a(W)=Q^{\dagger\dot a}(S);\ t_\alpha=gT^a", "notation-db/weinberg-qtf-1995.md D1,D2,D7,D9", "", "DERIVED", ["R6","R7"]),
]

_res = {}
for r_ in RESULTS:
    _res.setdefault(r_['id'], []).append(r_)
_dict = []
for (eid, cat, esw, ws, srcref, std, tag, chk) in ENTRIES:
    rs = [x for c_ in chk for x in _res.get(c_, [])]
    _dict.append(dict(id=eid, category=cat, esw_index_free=esw, weinberg_srednicki=ws, source=srcref, standard_entries=std,
                      tag=tag, checks=chk, max_residual=max((x['residual'] for x in rs), default=None),
                      all_passed=all(x['passed'] for x in rs) and len(rs) > 0))
missing = [e['id'] for e in _dict if not e['all_passed']]
log(f'dictionary.json: {len(_dict)} entries, entries without passing checks: {missing}')
json.dump(dict(title='ESW index-free -> Weinberg-Srednicki two-component dictionary (N=4 Costello-BV, Step 2)',
               date='2026-10-28', conventions=dict(
                   wick='v^0 = -i v^4, x^4 = i t, partial_0 = i partial_4; e^{i S_L} = e^{S/hbar}',
                   sigma_E='sigma_E^m = (sigma^1, sigma^2, sigma^3, i 1), sigmabar_E^m = (-sigma^1, -sigma^2, -sigma^3, i 1)',
                   orientation='dx^1 dx^2 dx^3 dx^4; S_+ = undotted = Srednicki (2,1)',
                   su4='W upper index i, W^* lower index i; Sigma^I_ij = engine B Sig, Sigmabar^I = engine B Sigb = (Sigma^I)^dagger',
                   engine_maps='engine B vector index (m=1,2,3,4) -> (1,2,3,0), I -> 3+I; engine A/classical: m = M+1, I = M-3',
                   not_used='STANDARD.md SIG/REP/CONF blocks'),
               entries=_dict), open(os.path.join(HERE, 'dictionary.json'), 'w'), indent=1)

# ======================================================================================================
npass = sum(r['passed'] for r in RESULTS)
log(f'=== SUMMARY: {npass}/{len(RESULTS)} checks passed; max residual {max(r["residual"] for r in RESULTS):.2e} ===')
json.dump(dict(results=RESULTS, intertwiner_scale_fixed_by='(U a, U* b)_A = (a,b)_D',
               vector_map={'D index 0..3 (m=1..4)': 'classical/engine A M = 0..3; engine B M = 1,2,3,0',
                           'D index 4..9 (I=1..6)': 'classical/engine A and engine B M = 4..9'}),
          open(os.path.join(HERE, 'check_results.json'), 'w'), indent=1)
open(os.path.join(HERE, 'check_dictionary.log'), 'w').write('\n'.join(LOG) + '\n')
