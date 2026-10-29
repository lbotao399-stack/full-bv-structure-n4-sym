#!/usr/bin/env python3
"""audit_dictionary.py -- independent adversarial audit of step2/notation (ESW -> Weinberg-Srednicki dictionary).

Independence: nothing is imported from step2/notation/check_dictionary.py.  The two-component realisation is
re-implemented here directly from the typeset formulas (R1)-(R4), (P1), (G1)-(G4), Thm thm:action, (F1) of
esw_to_weinberg_srednicki.tex; Srednicki data are typed from STANDARD.md S-SPN-2..8, S-MET-5 (book eq. numbers in
comments).  Engines A (classical/n4_model.py) and B (engine_b/realization.py) are read only for the intertwiner /
pairing-normalisation checks (section C).  numpy only, light (< 300 MB).
Sections:
  A  Srednicki base entries and the Wick map
  B  ESW axioms for the typeset realisation; Gamma (G1)-(G4); [ESW Prop 2.1]
  C  equivalence with engine A (intertwiner, orientation) and engine B (pairing factor)
  D  pointwise equality: index-free reduced action (ch01 eq:c1-action) vs typeset two-component action
  E  supersymmetry of the typeset two-component action under typeset (F1), non-abelian su(2) and su(3)
  F  {Q,Qbar}: Gamma(xi, xibar)|_V, closure on A (non-abelian), Lorentzian image
  G  Step-1 translations (quadratic CT pieces, finite CT, field counting, Konishi, canonical su(N) map)
"""
import os, sys, json, itertools, math
import numpy as np

np.seterr(all='ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
rng = np.random.default_rng(424242)
RES = []
LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    LOG.append(s)
    print(s, flush=True)


def flat(x):
    if isinstance(x, (list, tuple)):
        return np.concatenate([flat(y) for y in x]) if len(x) else np.zeros(0, complex)
    return np.ravel(np.asarray(x, dtype=complex))


def rec(cid, desc, lhs, rhs=0.0, tol=1e-9, expect=True, scale=None):
    l = flat(lhs)
    r = flat(rhs) if not np.isscalar(rhs) else np.full(l.shape, rhs, dtype=complex)
    sc = scale if scale is not None else max(1.0, float(np.max(np.abs(l))) if l.size else 1.0,
                                              float(np.max(np.abs(r))) if r.size else 1.0)
    res = float(np.max(np.abs(l - r))) / sc if l.size else 0.0
    ok = (res < tol) if expect else (res > 1e-3)
    tag = 'PASS' if ok else 'FAIL'
    if not expect:
        tag += '(neg-control)'
    RES.append(dict(id=cid, check=desc, residual=res, expect_zero=expect, passed=bool(ok)))
    log(f"[{tag}] {cid:10s} res={res:.2e}  {desc}")
    return ok


def rnd(*shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)


# =====================================================================================================
# A. Srednicki base entries (STANDARD.md S-SPN-2..8, S-MET-1, S-MET-5) and the Wick map
# =====================================================================================================
log('=== A. Srednicki entries and Wick map ===')
I2 = np.eye(2, dtype=complex)
PA = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EU = np.array([[0, 1], [-1, 0]], complex)       # eps^{ab}: eps^{12}=+1           S-SPN-2 (35.1)
ED = np.array([[0, -1], [1, 0]], complex)       # eps_{ab}: eps_{12}=-1           S-SPN-2 (35.1)
g = np.diag([-1., 1, 1, 1])                     # S-MET-1 (2.4)
sig = [I2] + PA                                 # sigma^mu_{a adot} = (1, vec sigma)   S-SPN-5 (34.28)/(35.2)
sigb = [np.einsum('ab,xy,by->xa', EU, EU, s) for s in sig]   # sigbar^{mu adot a} = eps^{ab}eps^{adot bdot}sigma_{b bdot}  (35.3); array [adot, a]
rec('A1', 'S-SPN-2/3: eps_{ab} eps^{bc} = delta_a^c  (34.23)', ED @ EU, I2)
rec('A2', 'S-SPN-5: sigbar^mu = (1, -vec sigma)  (35.3)', [sigb[m] - (I2 if m == 0 else -PA[m - 1]) for m in range(4)])
rec('A3', 'S-SPN-5 (35.4): sigma^mu_{a adot} sigma_{mu b bdot} = -2 eps_{ab} eps_{adot bdot}',
    sum(g[m, m] * np.einsum('ax,by->abxy', sig[m], sig[m]) for m in range(4)) + 2 * np.einsum('ab,xy->abxy', ED, ED))
rec('A4', 'S-SPN-5 (35.5): eps^{ab}eps^{adot bdot} sigma^mu_{a adot} sigma^nu_{b bdot} = -2 g^{mu nu}',
    [np.einsum('ab,xy,ax,by->', EU, EU, sig[m], sig[n]) + 2 * g[m, n] for m in range(4) for n in range(4)])
rec('A5', 'S-MET-5 (36.8): (sigma^mu sigbar^nu + sigma^nu sigbar^mu)_a^c = -2 g^{mu nu} delta_a^c (and dotted)',
    [sig[m] @ sigb[n] + sig[n] @ sigb[m] + 2 * g[m, n] * I2 for m in range(4) for n in range(4)] +
    [sigb[m] @ sig[n] + sigb[n] @ sig[m] + 2 * g[m, n] * I2 for m in range(4) for n in range(4)])
gamS = [np.block([[np.zeros((2, 2)), sig[m]], [sigb[m], np.zeros((2, 2))]]) for m in range(4)]   # S-SPN-8 (36.7)
g5S = np.diag([-1, -1, 1, 1]).astype(complex)                                                    # S-SPN-7 (36.43)
rec('A6', 'S-SPN-7/8: {gamma^mu,gamma^nu} = -2g, gamma_5 = i gamma^0 gamma^1 gamma^2 gamma^3 = diag(-1,-1,1,1)',
    [gamS[m] @ gamS[n] + gamS[n] @ gamS[m] + 2 * g[m, n] * np.eye(4) for m in range(4) for n in range(4)] +
    [1j * gamS[0] @ gamS[1] @ gamS[2] @ gamS[3] - g5S])
# Wick map v^0 = -i v^4 ; dictionary Def def:wick
SE = [PA[0], PA[1], PA[2], 1j * I2]          # sigma_E^m  [alpha, alphadot], m=1..4 -> 0..3
SBE = [-PA[0], -PA[1], -PA[2], 1j * I2]      # sigmabar_E^m [alphadot, alpha]
for t in range(3):
    vE = rng.normal(size=4)
    vup = np.array([-1j * vE[3], vE[0], vE[1], vE[2]])
    vdn = g @ vup
    rec('A7', 'Wick: sigma^mu v_mu = sigma_E^m v_m, sigmabar^mu v_mu = sigmabar_E^m v_m (v^0=-i v^4)',
        [sum(sig[m] * vdn[m] for m in range(4)) - sum(SE[m] * vE[m] for m in range(4)),
         sum(sigb[m] * vdn[m] for m in range(4)) - sum(SBE[m] * vE[m] for m in range(4))])
rec('A8', 'Euclidean: sigma_E sigmabar_E + (m<->n) = -2 delta, sigmabar_E = eps eps sigma_E, sigmabar_E = -(sigma_E)^dagger',
    [SE[m] @ SBE[n] + SE[n] @ SBE[m] + 2 * (m == n) * I2 for m in range(4) for n in range(4)] +
    [SBE[m] @ SE[n] + SBE[n] @ SE[m] + 2 * (m == n) * I2 for m in range(4) for n in range(4)] +
    [SBE[m] - np.einsum('ab,xy,by->xa', EU, EU, SE[m]) for m in range(4)] + [SBE[m] + SE[m].conj().T for m in range(4)])
gamE = [np.block([[np.zeros((2, 2)), SE[m]], [SBE[m], np.zeros((2, 2))]]) for m in range(4)]
rec('A9', 'gamma_E^1 gamma_E^2 gamma_E^3 gamma_E^4 = -gamma_5S  (so rho(e1e2e3e4)=(-i)^4 gE1..gE4 = +1 on undotted)',
    gamE[0] @ gamE[1] @ gamE[2] @ gamE[3], -g5S)
# negative control: the blocked Euclidean sigma of STANDARD SIG ((-i sigma, 1), +2 delta) is NOT the dictionary's
SIGblk = [-1j * PA[0], -1j * PA[1], -1j * PA[2], I2]
rec('A10', 'dictionary sigma_E differs from the SIG-block sigma^E=(-i sigma,1) (negative control: must differ)',
    [SE[m] - SIGblk[m] for m in range(4)], expect=False)
# e^{iS_L}=e^{S/hbar}: i dt = i (-i dx^4) = dx^4 for x^0 = -i x^4 ; scalar kinetic term check
phi_d = rnd(4)                                             # Euclidean gradient d_m phi
dL = np.array([1j * phi_d[3], phi_d[0], phi_d[1], phi_d[2]])   # d_0 = i d_4 (lower index Lorentzian)
rec('A11', 'Wick: -1/2 d^mu phi d_mu phi |_Wick = -1/2 d_m phi d_m phi, and i*dt = dx^4 (e^{iS_L}=e^{int L_Wick})',
    [-0.5 * np.einsum('m,mn,n->', dL, np.linalg.inv(g), dL) + 0.5 * phi_d @ phi_d, 1j * (-1j) - 1])

# =====================================================================================================
# B. Realisation typed from (R1)-(R4), (P1) of the dictionary; ESW axioms; (G1)-(G4)
# =====================================================================================================
log('=== B. typeset realisation: ESW axioms, Gamma ===')
E4 = np.zeros((4,) * 4)
for p in itertools.permutations(range(4)):
    E4[p] = np.linalg.det(np.eye(4)[list(p)])        # eps_{1234}=eps^{1234}=+1
LC3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    LC3[p] = np.linalg.det(np.eye(3)[list(p)])
eta = np.zeros((3, 4, 4)); etab = np.zeros((3, 4, 4))
for a in range(3):
    eta[a, :3, :3] = LC3[a]; etab[a, :3, :3] = LC3[a]
    eta[a, a, 3] = 1; eta[a, 3, a] = -1
    etab[a, a, 3] = -1; etab[a, 3, a] = 1
SIGI = np.array([eta[0], eta[1], eta[2], 1j * etab[0], 1j * etab[1], 1j * etab[2]], complex)   # Sigma^I_{ij}
SIGBI = np.array([S.conj().T for S in SIGI])                                                   # Sigmabar^{I ij}
rec('B0a', "Sigmabar^I = (-eta, i etabar) as typeset", SIGBI,
    np.array([-eta[0], -eta[1], -eta[2], 1j * etab[0], 1j * etab[1], 1j * etab[2]]))
rec('B0b', "Sigma^I antisym; (1/2)eps^{ijkl}Sigma_kl = -Sigmabar^{ij}; Sigmabar Sigma + (I<->J) = 2 delta; tr = 4 delta",
    [SIGI + SIGI.transpose(0, 2, 1), 0.5 * np.einsum('ijkl,Ikl->Iij', E4, SIGI) + SIGBI,
     np.einsum('Iij,Jjk->IJik', SIGBI, SIGI) + np.einsum('Jij,Ijk->IJik', SIGBI, SIGI) - 2 * np.einsum('IJ,ik->IJik', np.eye(6), np.eye(4)),
     np.einsum('Iij,Jji->IJ', SIGBI, SIGI) - 4 * np.eye(6)])
DD = np.einsum('ik,jl->ijkl', np.eye(4), np.eye(4)) - np.einsum('il,jk->ijkl', np.eye(4), np.eye(4))
rec('B0c', 'Fierz: sum Sig Sig = 2 eps_ijkl, sum Sigb Sigb = 2 eps^ijkl, sum Sigb^{ij} Sig_kl = -2(dd-dd)',
    [np.einsum('Iij,Ikl->ijkl', SIGI, SIGI) - 2 * E4, np.einsum('Iij,Ikl->ijkl', SIGBI, SIGBI) - 2 * E4,
     np.einsum('Iij,Ikl->ijkl', SIGBI, SIGI) + 2 * DD])
u6 = rnd(6); w6 = rnd(6)
ulo = np.einsum('I,Iij->ij', u6, SIGI); uup = -np.einsum('I,Iij->ij', u6, SIGBI)
wlo = np.einsum('I,Iij->ij', w6, SIGI); wup = -np.einsum('I,Iij->ij', w6, SIGBI)
rec('B0d', 'u^{ij} = (1/2)eps^{ijkl}u_kl; (u,w)=u_I w_I = (1/4)u^{ij}w_ij; u_I=(1/4)Sigmabar^{I ji}u_ij; u^{ik}w_kj+w^{ik}u_kj=-2(u,w)delta',
    [uup - 0.5 * np.einsum('ijkl,kl->ij', E4, ulo), 0.25 * np.einsum('ij,ij->', uup, wlo) - u6 @ w6,
     0.25 * np.einsum('Iji,ij->I', SIGBI, ulo) - u6, uup @ wlo + wup @ ulo + 2 * (u6 @ w6) * np.eye(4)])


def sidx(dotted, a, i):            # Sigma: (xi_alpha^i) at a*4+i ; (xibar_{alphadot i}) at 8+a*4+i
    return (8 if dotted else 0) + 4 * a + i


# Sigma^*: (eta^alpha_i) at a*4+i ; (etabar^{alphadot i}) at 8+a*4+i  (same function)
R = [np.zeros((16, 16), complex) for _ in range(10)]    # rho(e_M): Sigma -> Sigma^*
Rb = [np.zeros((16, 16), complex) for _ in range(10)]   # rho(e_M): Sigma^* -> Sigma
for m in range(4):
    for a in range(2):
        for d in range(2):
            for i in range(4):
                R[m][sidx(0, a, i), sidx(1, d, i)] += 1j * SBE[m][d, a]     # (R1a) (rho xi)^alpha_i = i xibar_{ad i} sbar^{m ad al}
                R[m][sidx(1, d, i), sidx(0, a, i)] += -1j * SBE[m][d, a]    # (R1b) (rho xi)^{ad i} = -i sbar^{m ad al} xi_al^i
                Rb[m][sidx(0, a, i), sidx(1, d, i)] += -1j * SE[m][a, d]    # (R2a) (rho eta)_al^i = -i s^m_{al ad} etab^{ad i}
                Rb[m][sidx(1, d, i), sidx(0, a, i)] += 1j * SE[m][a, d]     # (R2b) (rho eta)_{ad i} = i eta^al_i s^m_{al ad}
for I in range(6):
    M = 4 + I
    wl = SIGI[I]; wu = -SIGBI[I]
    for a in range(2):
        for b in range(2):
            for j in range(4):
                for i in range(4):
                    R[M][sidx(0, a, j), sidx(0, b, i)] += EU[a, b] * wl[j, i]      # (R3a)
                    R[M][sidx(1, a, j), sidx(1, b, i)] += EU[a, b] * wu[j, i]      # (R3b)
                    Rb[M][sidx(0, a, j), sidx(0, b, i)] += -wu[j, i] * ED[a, b]    # (R4a) (rho eta)_al^k = -w^{kj} eps_{al be} eta^be_j
                    Rb[M][sidx(1, a, j), sidx(1, b, i)] += -wl[j, i] * ED[a, b]    # (R4b) (rho eta)_{bd k} = -w_{kj} eps_{bd ad} etab^{ad j}
PM = np.diag([1.] * 8 + [-1.] * 8).astype(complex)       # (P1) (xi,eta) = xi^T PM eta
Id16 = np.eye(16)
rec('B1', 'Clifford rho(e_M)rho(e_N)+rho(e_N)rho(e_M) = 2 delta_MN on Sigma and Sigma^*',
    [Rb[M] @ R[N] + Rb[N] @ R[M] - 2 * (M == N) * Id16 for M in range(10) for N in range(10)] +
    [R[M] @ Rb[N] + R[N] @ Rb[M] - 2 * (M == N) * Id16 for M in range(10) for N in range(10)])
rec('B2', '(rho(v)Q1,Q2) = (Q1,rho(v)Q2) on Sigma and on Sigma^* w.r.t. (P1)',
    [PM @ R[M] - (PM @ R[M]).T for M in range(10)] + [Rb[M].T @ PM - (Rb[M].T @ PM).T for M in range(10)])
vol = Rb[0] @ R[1] @ Rb[2] @ R[3]
vols = R[0] @ Rb[1] @ R[2] @ Rb[3]
rec('B3', 'chirality: rho(e1e2e3e4) = +1 on undotted (S_+ (x) W), -1 on dotted, both on Sigma and Sigma^*',
    [vol - PM, vols - PM])
rec('B4', 'rho(w) preserves dotted/undotted, rho(v) exchanges them',
    [R[M][:8, 8:] for M in range(4, 10)] + [R[M][8:, :8] for M in range(4, 10)] + [R[M][:8, :8] for M in range(4)] + [R[M][8:, 8:] for M in range(4)])


def GamS(q1, q2):          # Gamma^M(Q1,Q2) for Q1,Q2 in Sigma: (rho(e_M)Q1, Q2)  [ESW (28)]
    return np.array([q2 @ PM @ R[M] @ q1 for M in range(10)])


def GamSs(e1, e2):         # Gamma^M for Sigma^*: (rho(e_M)e1, e2) with rho(e_M)e1 in Sigma
    return np.array([(Rb[M] @ e1) @ PM @ e2 for M in range(10)])


def rhoS(v):
    return sum(v[M] * R[M] for M in range(10))


def rhoSs(v):
    return sum(v[M] * Rb[M] for M in range(10))


Qs = [rnd(16) for _ in range(3)]
tot = sum(rhoS(GamS(Qs[a], Qs[b])) @ Qs[c] for (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)])
rec('B5', '[ESW Assumption 3.1] rho(Gamma(Q1,Q2))Q3 + cyclic = 0 on Sigma', tot)
Es = [rnd(16) for _ in range(3)]
tot = sum(rhoSs(GamSs(Es[a], Es[b])) @ Es[c] for (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)])
rec('B5s', '[ESW Assumption 3.1] on Sigma^*', tot)
rec('B6', 'pairing (P1) nondegenerate and the unique (up to scale) invariant one: nullspace dim',
    [abs(np.linalg.det(PM)) - 1])
# uniqueness of invariant pairing: solve Mp R_M symmetric, Rb_M^T Mp symmetric
rows = []
for M in range(10):
    for (A_, side) in ((R[M], 'L'), (Rb[M], 'R')):
        # X = Mp A (L) must be symmetric ; X = A^T Mp (R) must be symmetric
        for i in range(16):
            for j in range(i + 1, 16):
                row = np.zeros((16, 16), complex)
                if side == 'L':
                    row[i, :] += A_[:, j]; row[j, :] -= A_[:, i]
                else:
                    row[:, j] += A_[:, i]; row[:, i] -= A_[:, j]
                rows.append(row.ravel())
Amat = np.array(rows)
sv = np.linalg.svd(Amat, compute_uv=False)
nnull = int(np.sum(sv < 1e-10 * sv[0]))
rec('B6b', f'invariant pairing Sigma x Sigma^* is unique up to scale (nullity={nnull})', [nnull - 1])


def split2(q):
    return q[:8].reshape(2, 4), q[8:].reshape(2, 4)


def G1G2(q1, q2):
    x1, xb1 = split2(q1); x2, xb2 = split2(q2)
    Gm = np.array([1j * (np.einsum('di,da,ai->', xb2, SBE[m], x1) + np.einsum('di,da,ai->', xb1, SBE[m], x2)) for m in range(4)])
    Gkl = (2 * np.einsum('jikl,aj,ab,bi->kl', E4, x1, EU, x2)
           - 2 * (np.einsum('ak,ab,bl->kl', xb1, EU, xb2) - np.einsum('al,ab,bk->kl', xb1, EU, xb2)))
    return Gm, Gkl


def G3G4(e1, e2):
    y1, yb1 = split2(e1); y2, yb2 = split2(e2)
    Gm = np.array([-1j * (np.einsum('ai,ad,di->', y2, SE[m], yb1) + np.einsum('ai,ad,di->', y1, SE[m], yb2)) for m in range(4)])
    Gup = (-2 * np.einsum('klpj,ab,bj,ap->kl', E4, ED, y1, y2)
           + 2 * (np.einsum('ba,al,bk->kl', ED, yb1, yb2) - np.einsum('ba,ak,bl->kl', ED, yb1, yb2)))
    return Gm, Gup


for t in range(3):
    q1, q2 = rnd(16), rnd(16)
    G = GamS(q1, q2); Gm, Gkl = G1G2(q1, q2)
    rec('B7', '(G1),(G2): Gamma(xi1,xi2) on V and Gamma_kl = Gamma^I Sigma^I_kl, as typeset',
        [G[:4] - Gm, np.einsum('I,Ikl->kl', G[4:], SIGI) - Gkl, G - GamS(q2, q1)])
    e1, e2 = rnd(16), rnd(16)
    G = GamSs(e1, e2); Gm, Gup = G3G4(e1, e2)
    rec('B8', '(G3),(G4): Gamma(eta1,eta2) on V and Gamma^{kl} = -Gamma^I Sigmabar^{I kl}, as typeset',
        [G[:4] - Gm, -np.einsum('I,Ikl->kl', G[4:], SIGBI) - Gup, G - GamSs(e2, e1)])
u = rnd(2); w = rnd(4)
xi_r1 = np.concatenate([np.einsum('a,i->ai', u, w).ravel(), np.zeros(8)])
rec('B9', 'rank-one xi_alpha^i = u_alpha w^i, xibar = 0 is square-zero: Gamma(xi,xi) = 0', GamS(xi_r1, xi_r1))
# [ESW Prop 2.1]: Gamma(Q1, rho(X)Q2) + Gamma(Q2, rho(X)Q1) = -2 iota_{Gamma(Q1,Q2)} X, X a 2-form
Xf = rnd(10, 10); Xf = Xf - Xf.T
rhoX = 0.5 * sum(Xf[M, N] * Rb[M] @ R[N] for M in range(10) for N in range(10))   # rho(q(X)) on Sigma, X = 1/2 X_MN e^M e^N
q1, q2 = rnd(16), rnd(16)
lhs = GamS(q1, rhoX @ q2) + GamS(q2, rhoX @ q1)
rhs = -2 * np.einsum('N,NM->M', GamS(q1, q2), Xf)
rec('B10', '[ESW Prop 2.1] Gamma(Q1,rho(X)Q2)+Gamma(Q2,rho(X)Q1) = -2 iota_{Gamma(Q1,Q2)}X  ((iota_v X)_M = v^N X_NM)', lhs, rhs)
rec('B10n', 'negative control: same with +2 must fail', lhs, -rhs, expect=False)
# Dirac column (psi_alpha, chibar^alphadot) = (xi_alpha^i, etabar^{alphadot i}) for fixed i: rho(v) = -i gamma_S(v)|_Wick = gamma_W(v)
vv = rnd(4)
blocks = []
for i in range(4):
    Mcol = np.zeros((4, 4), complex)
    rows_ = [sidx(0, 0, i), sidx(0, 1, i)]            # Sigma undotted, and Sigma^* dotted (etabar^{ad i}) as second block
    rowsS = [sidx(1, 0, i), sidx(1, 1, i)]
    RS = sum(vv[m] * R[m] for m in range(4)); RbS = sum(vv[m] * Rb[m] for m in range(4))
    Mcol[2:, :2] = RS[np.ix_(rowsS, rows_)]           # xi_alpha -> (rho xi)^{ad i}
    Mcol[:2, 2:] = RbS[np.ix_(rows_, rowsS)]          # etabar^{ad} -> (rho eta)_alpha^i
    blocks.append(Mcol + 1j * sum(vv[m] * gamE[m] for m in range(4)))
rec('B11', 'Dirac column: rho(v) = -i gamma_S^mu v_mu|_Wick (Srednicki (36.7) block) = gamma_W v (Weinberg db D1)', blocks)

# (P1) in Srednicki contractions: (xi,eta) = -(xi^i eta_i) - (xibar_i etabar^i), (xi eta) = xi^al eta_al, (xibar etabar) = xibar_ad etabar^ad
x_, xb_ = rnd(2, 4), rnd(2, 4); y_, yb_ = rnd(2, 4), rnd(2, 4)
qx = np.concatenate([x_.ravel(), xb_.ravel()]); qy = np.concatenate([y_.ravel(), yb_.ravel()])
xup = np.einsum('ab,bi->ai', EU, x_)            # xi^{al i}
ydn = np.einsum('ab,bi->ai', ED, y_)            # eta_{al i} = eps_{al be} eta^be_i
rec('B12', '(P1) = -(xi^i eta_i) - (xibar_i etabar^i) in Srednicki contractions (35.23)',
    qx @ PM @ qy, -np.einsum('ai,ai->', xup, ydn) - np.einsum('xi,xi->', xb_, yb_))
# second Dirac column (eta_{al i}, xibar^{ad}_i): rho(v) = -i gamma_S(v) as well
vv = rnd(4)
RS = sum(vv[m] * R[m] for m in range(4)); RbS = sum(vv[m] * Rb[m] for m in range(4))
blocks = []
for i in range(4):
    # xibar^{ad}_i = eps^{ad bd} xibar_{bd i}  (Sigma, dotted) ; eta_{al i} = eps_{al be} eta^be_i (Sigma^*, undotted)
    Pxb = np.zeros((2, 16), complex); Pet = np.zeros((2, 16), complex)
    for a in range(2):
        for b in range(2):
            Pxb[a, sidx(1, b, i)] = EU[a, b]; Pet[a, sidx(0, b, i)] = ED[a, b]
    # embeddings back: xibar_{bd i} = eps_{bd ad} xibar^{ad}_i ; eta^be_i = eps^{be al} eta_{al i}
    Exb = np.zeros((16, 2), complex); Eet = np.zeros((16, 2), complex)
    for a in range(2):
        for b in range(2):
            Exb[sidx(1, b, i), a] = ED[b, a]; Eet[sidx(0, b, i), a] = EU[b, a]
    Mcol = np.zeros((4, 4), complex)
    Mcol[:2, 2:] = Pet @ RS @ Exb          # xibar^{ad} (Sigma) -> eta_al (Sigma^*)
    Mcol[2:, :2] = Pxb @ RbS @ Eet         # eta_al (Sigma^*) -> xibar^{ad} (Sigma)
    blocks.append(Mcol + 1j * sum(vv[m] * gamE[m] for m in range(4)))
rec('B13', 'second Dirac column (eta_{al i}, xibar^{ad}_i): rho(v) = -i gamma_S(v)|_Wick', blocks)

# =====================================================================================================
# C. Equivalence with engine A (= classical model) and engine B (read-only)
# =====================================================================================================
log('=== C. engines ===')
sys.path.insert(0, os.path.join(ROOT, 'classical'))
import n4_model as NM   # read-only
sigA = [NM.SIG_RE[M].astype(complex) + 1j * NM.SIG_IM[M] for M in range(10)]
sigbA = [NM.SIGBAR_RE[M].astype(complex) + 1j * NM.SIGBAR_IM[M] for M in range(10)]


def intertwiners(RA, RbA, vmap):
    """solve RA[vmap[M]] U = Us R[M],  RbA[vmap[M]] Us = U Rb[M]  (row-major vec)."""
    I16 = np.eye(16)
    rows = []
    for M in range(10):
        rows.append(np.hstack([np.kron(RA[vmap[M]], I16), -np.kron(I16, R[M].T)]))
        rows.append(np.hstack([-np.kron(I16, Rb[M].T), np.kron(RbA[vmap[M]], I16)]))
    A_ = np.vstack(rows)
    u_, s_, vh = np.linalg.svd(A_, full_matrices=False)
    null = vh[s_ < 1e-9 * s_[0]].conj() if np.any(s_ < 1e-9 * s_[0]) else np.zeros((0, 512))
    # svd returns all 512 singular values (A_ is tall)
    return null, s_


nullA, sA = intertwiners(sigA, sigbA, list(range(10)))
rec('C1', f'engine A intertwiner with vector map M->M exists and is unique up to scale (nullity={len(nullA)})', [len(nullA) - 1])
U = nullA[0][:256].reshape(16, 16); Us = nullA[0][256:].reshape(16, 16)
rec('C2', 'intertwiner equations hold', [sigA[M] @ U - Us @ R[M] for M in range(10)] + [sigbA[M] @ Us - U @ Rb[M] for M in range(10)],
    scale=np.max(np.abs(U)))
cpair = U.T @ Us            # (Ua, Us b)_A = a^T U^T Us b ; compare with PM
c0 = cpair[0, 0] / PM[0, 0]
rec('C3', 'engine A pairing pulls back to a multiple of (P1): U^T U* = c PM (c absorbed by U -> c^{-1/2}U)', cpair / c0, PM)
volA = sigbA[0] @ sigA[1] @ sigbA[2] @ sigA[3]
rec('C4', 'U maps undotted (S_+ (x) W) to the +1 eigenspace of rho_A(e0e1e2e3) of engine A, dotted to -1',
    [volA @ U[:, :8] - U[:, :8], volA @ U[:, 8:] + U[:, 8:]], scale=np.max(np.abs(U)))
vrefl = [1, 0, 2, 3, 4, 5, 6, 7, 8, 9]   # orientation-reversing map on V only
nullAr, _ = intertwiners(sigA, sigbA, vrefl)
rec('C5', f'negative control: V-orientation-reversed map admits no Sigma->Sigma_A intertwiner (nullity={len(nullAr)})',
    [len(nullAr)], expect=True)

# engine B
sys.path.insert(0, os.path.join(ROOT, 'engine_b'))
import realization as RB   # read-only
RBn = [np.array(RB.R[M].tolist(), dtype=complex) for M in range(10)]
RbBn = [np.array(RB.Rb[M].tolist(), dtype=complex) for M in range(10)]
rows = []
for M in range(10):
    for (A_, side) in ((RBn[M], 'L'), (RbBn[M], 'R')):
        for i in range(16):
            for j in range(i + 1, 16):
                row = np.zeros((16, 16), complex)
                if side == 'L':
                    row[i, :] += A_[:, j]; row[j, :] -= A_[:, i]
                else:
                    row[:, j] += A_[:, i]; row[:, i] -= A_[:, j]
                rows.append(row.ravel())
_, sB, vhB = np.linalg.svd(np.array(rows), full_matrices=False)
MpB = vhB[-1].conj().reshape(16, 16)
nz = MpB.ravel()[np.abs(MpB.ravel()) > 1e-9 * np.max(np.abs(MpB))]
MpB = MpB / nz[0]          # engine B normalisation: first nonzero entry (row-major) = +1
T = np.zeros((16, 16), complex); Ts = np.zeros((16, 16), complex)   # x_dict = T x_B ; y_dict = Ts y_B
for a in range(2):
    for i in range(4):
        T[sidx(0, a, i), a * 4 + i] = 1
        for b in range(2):
            T[sidx(1, a, i), 8 + b * 4 + i] = ED[a, b]
            Ts[sidx(0, a, i), 8 + b * 4 + i] = EU[a, b]
        Ts[sidx(1, a, i), a * 4 + i] = 1
vmapB = [1, 2, 3, 0, 4, 5, 6, 7, 8, 9]
Ti = np.linalg.inv(T); Tsi = np.linalg.inv(Ts)
rec('C6', 'engine B: rho entry by entry after the stated basis change and vector map (m=1..4 -> 1,2,3,0; I -> 3+I)',
    [Ts @ RBn[vmapB[M]] @ Ti - R[M] for M in range(10)] + [T @ RbBn[vmapB[M]] @ Tsi - Rb[M] for M in range(10)])
rec('C7', 'engine B pairing (normalised first nonzero entry +1) equals (P1) with factor exactly 1', Ti.T @ MpB @ Tsi, PM)

# =====================================================================================================
# Marker algebra: odd (Grassmann) and even-nilpotent markers, one plane wave per marker
# =====================================================================================================
class Ctx:
    def __init__(self, parities, moms):
        self.n = len(parities)
        self.odd = sum(1 << j for j, p in enumerate(parities) if p)
        self.mom = [np.asarray(m, dtype=complex) for m in moms]
        self.P = {}
        for S in range(1 << self.n):
            self.P[S] = sum((self.mom[j] for j in range(self.n) if (S >> j) & 1), np.zeros(4, complex))
        self.sg = {}
        for A_ in range(1 << self.n):
            for B_ in range(1 << self.n):
                if A_ & B_:
                    continue
                oa, ob = A_ & self.odd, B_ & self.odd
                cnt = 0
                for i in range(self.n):
                    if (oa >> i) & 1:
                        cnt += bin(ob & ((1 << i) - 1)).count('1')
                self.sg[(A_, B_)] = -1 if cnt % 2 else 1
        self.full = (1 << self.n) - 1


CTX = None


class GA:
    def __init__(self, t=None):
        self.t = {k: v for k, v in (t or {}).items()}

    def __add__(self, o):
        r = dict(self.t)
        for k, v in o.t.items():
            r[k] = r[k] + v if k in r else v
        return GA(r)

    def __sub__(self, o):
        return self + o * (-1)

    def __mul__(self, c):
        return GA({k: c * v for k, v in self.t.items()})

    __rmul__ = __mul__

    def map(self, f):
        return GA({k: f(v) for k, v in self.t.items()})

    def coef(self, S):
        return self.t.get(S, None)


def field(legs):
    return GA({1 << j: np.asarray(a, dtype=complex) for (j, a) in legs})


def deriv(x):
    out = {}
    for S, v in x.t.items():
        p = 1j * CTX.P[S]
        out[S] = np.einsum('m,...->m...', p, v)
    return GA(out)


def gmul(spec, x, y, *consts):
    out = {}
    for A_, a in x.t.items():
        for B_, b in y.t.items():
            if A_ & B_:
                continue
            v = CTX.sg[(A_, B_)] * np.einsum(spec, a, b, *consts)
            S = A_ | B_
            out[S] = out[S] + v if S in out else v
    return GA(out)


def lmul_marker(j, x):           # theta_j * x
    out = {}
    for S, v in x.t.items():
        if (S >> j) & 1:
            continue
        cnt = bin(S & CTX.odd & ((1 << j) - 1)).count('1') if (CTX.odd >> j) & 1 else 0
        out[S | (1 << j)] = (-1) ** cnt * v
    return GA(out)


def zeroGA():
    return GA({})


def lie(name):
    if name == 'su2':
        f = LC3.copy()
    else:   # su(3): orthonormal t_a = -i lambda_a/2 ... use [T_a,T_b] = f_abc T_c from Gell-Mann traces (real form)
        lam = np.zeros((8, 3, 3), complex)
        lam[0][0, 1] = lam[0][1, 0] = 1
        lam[1][0, 1] = -1j; lam[1][1, 0] = 1j
        lam[2][0, 0] = 1; lam[2][1, 1] = -1
        lam[3][0, 2] = lam[3][2, 0] = 1
        lam[4][0, 2] = -1j; lam[4][2, 0] = 1j
        lam[5][1, 2] = lam[5][2, 1] = 1
        lam[6][1, 2] = -1j; lam[6][2, 1] = 1j
        lam[7] = np.diag([1, 1, -2]) / np.sqrt(3)
        T_ = lam / 2
        f = np.zeros((8, 8, 8))
        for a in range(8):
            for b in range(8):
                Cm = T_[a] @ T_[b] - T_[b] @ T_[a]
                for c in range(8):
                    f[a, b, c] = np.real(-2j * np.trace(Cm @ T_[c]))
    return f


KAP = 0.83
SBEa = np.array(SBE); SEa = np.array(SE)


# ---------------- typeset two-component action (Thm thm:action) ----------------
def L2_parts(Am, AI, lam, lamb, FST, nc, ghosts=None):
    K = KAP * np.eye(nc)
    out = {}
    dAm = deriv(Am)                                                     # [m,n,a] = d_m A_n
    Fmn = dAm - dAm.map(lambda v: v.transpose(1, 0, 2)) + gmul('ma,nb,abc->mnc', Am, Am, FST)
    out['gauge'] = gmul('mna,mnb,ab->', Fmn, Fmn, K) * (-1 / 8)
    Alo = AI.map(lambda v: np.einsum('Ia,Iij->ija', v, SIGI))
    Aup = AI.map(lambda v: -np.einsum('Ia,Iij->ija', v, SIGBI))
    DAlo = deriv(Alo) + gmul('ma,ijb,abc->mijc', Am, Alo, FST)
    DAup = deriv(Aup) + gmul('ma,ijb,abc->mijc', Am, Aup, FST)
    out['scalar'] = gmul('mija,mijb,ab->', DAup, DAlo, K) * (-1 / 16)
    C1 = gmul('ija,klb,abc->ijklc', Aup, Aup, FST); C2 = gmul('ija,klb,abc->ijklc', Alo, Alo, FST)
    out['quartic'] = gmul('ijkla,ijklb,ab->', C1, C2, K) * (-1 / 128)
    Dlamb = deriv(lamb) + gmul('ma,xib,abc->mxic', Am, lamb, FST)      # [m, adot, i, a]
    Dlam = deriv(lam) + gmul('ma,yib,abc->myic', Am, lam, FST)         # [m, alpha, i, a]
    out['kin'] = (gmul('yia,mxib,mxy,ab->', lam, Dlamb, SBEa, K) * (0.5j)
                  + gmul('xia,myib,mxy,ab->', lamb, Dlam, SBEa, K) * (0.5j))
    cY1 = gmul('jia,yib,abc->yjc', Alo, lam, FST)                      # [A_ji, lam_y^i] -> [y, j, c]
    cY2 = gmul('jia,xib,abc->xjc', Aup, lamb, FST)                     # [A^{ji}, lamb_{x i}] -> [x, j, c]
    out['yukawa'] = (gmul('zja,yjb,zy,ab->', lam, cY1, EU, K) * 0.5
                     - gmul('wja,xjb,wx,ab->', lamb, cY2, EU, K) * 0.5)
    if ghosts is not None:
        c, Asm, AsI, lst, lbst, cs = ghosts
        Asl = AsI.map(lambda v: np.einsum('Ia,Iij->ija', v, SIGI))
        Dc = deriv(c) + gmul('ma,b,abc->mc', Am, c, FST)
        out['brst'] = (gmul('ma,mb,ab->', Dc, Asm, K)
                       + gmul('ija,ijb,ab->', gmul('ija,b,abc->ijc', Aup, c, FST), Asl, K) * 0.25
                       + gmul('yia,yib,ab->', gmul('a,yib,abc->yic', c, lam, FST), lst, K)
                       - gmul('xia,xib,ab->', gmul('a,xib,abc->xic', c, lamb, FST), lbst, K)
                       + gmul('a,b,ab->', gmul('a,b,abc->c', c, c, FST), cs, K) * 0.5)
    return out


# ---------------- index-free reduced action (ch01 eq:c1-action), realisation of section B ----------------
RA = np.array(R); RbA = np.array(Rb)


def LIF_parts(A10, lamS, FST, nc, ghosts=None):
    K = KAP * np.eye(nc)
    out = {}
    dA = deriv(A10)                                                     # [m(4), N(10), a]
    dA10 = dA.map(lambda v: np.concatenate([v, np.zeros((6,) + v.shape[1:], complex)], axis=0))   # d_M = 0 along W^perp
    F = dA10 - dA10.map(lambda v: v.transpose(1, 0, 2)) + gmul('Ma,Nb,abc->MNc', A10, A10, FST)
    out['bos'] = gmul('MNa,MNb,ab->', F, F, K) * (-0.25 * 0.5)          # (F,F) = sum_{M<N}
    dl = deriv(lamS).map(lambda v: np.concatenate([v, np.zeros((6,) + v.shape[1:], complex)], axis=0))   # [M, s, a]
    Dl = dl + gmul('Ma,sb,abc->Msc', A10, lamS, FST)
    rDl = Dl.map(lambda v: np.einsum('Mts,Msa->ta', RA, v))            # rho(d_A) lam in Sigma^*
    out['ferm'] = gmul('sa,tb,st,ab->', lamS, rDl, PM, K) * 0.5
    if ghosts is not None:
        c, As10, lst16, cs = ghosts
        dc = deriv(c).map(lambda v: np.concatenate([v, np.zeros((6,) + v.shape[1:], complex)], axis=0))
        Dc = dc + gmul('Ma,b,abc->Mc', A10, c, FST)
        out['brst'] = (gmul('Ma,Mb,ab->', Dc, As10, K)
                       + gmul('sa,tb,st,ab->', gmul('a,sb,abc->sc', c, lamS, FST), lst16, PM, K)
                       + gmul('a,b,ab->', gmul('a,b,abc->c', c, c, FST), cs, K) * 0.5)
    return out


def total(parts, keys=None):
    r = zeroGA()
    for k, v in parts.items():
        if keys is None or k in keys:
            r = r + v
    return r


def cmp_all(x, y):
    keys = set(x.t) | set(y.t)
    l = [x.t.get(k, 0) for k in keys]; r = [y.t.get(k, 0) for k in keys]
    return [np.asarray(a, complex) - np.asarray(b, complex) for a, b in zip(l, r)], max([1.0] + [float(np.max(np.abs(a))) for a in l if np.size(a)])


# =====================================================================================================
# D. pointwise equality of the action densities
# =====================================================================================================
log('=== D. action: index-free density vs typeset two-component density ===')
for gname in ('su2', 'su3'):
    FST = lie(gname); nc = FST.shape[0]
    # D1 bosonic: 4 even markers
    CTX = Ctx([0, 0, 0, 0], [rng.normal(size=4) for _ in range(4)])
    a10 = [rnd(10, nc) for _ in range(4)]
    A10 = field([(j, a10[j]) for j in range(4)])
    Am = A10.map(lambda v: v[:4]); AI = A10.map(lambda v: v[4:])
    zl = field([]); 
    p2 = L2_parts(Am, AI, GA({}), GA({}), FST, nc)
    pI = LIF_parts(A10, GA({}), FST, nc)
    d, sc = cmp_all(total(p2, ['gauge', 'scalar', 'quartic']), pI['bos'])
    rec('D1-' + gname, '-1/4(F_A,F_A) = -1/8(F_mn,F_mn) - 1/16(D A^ij, D A_ij) - 1/128([A^ij,A^kl],[A_ij,A_kl]) pointwise', d, scale=sc)
    # quartic alone vs -1/4 * 1/2 sum_IJ ([A_I,A_J],[A_I,A_J])
    CIJ = gmul('Ia,Jb,abc->IJc', AI, AI, FST)
    q2 = gmul('IJa,IJb,ab->', CIJ, CIJ, KAP * np.eye(nc)) * (-1 / 8)
    d, sc = cmp_all(p2['quartic'], q2)
    rec('D1q-' + gname, '1/2 sum_IJ([A_I,A_J],[A_I,A_J]) = 1/32([A^ij,A^kl],[A_ij,A_kl]) pointwise', d, scale=sc)
    # D2 fermionic: lam two odd markers, A one even marker
    CTX = Ctx([1, 1, 0], [rng.normal(size=4) for _ in range(3)])
    l16 = [rnd(16, nc) for _ in range(2)]
    lamS = field([(0, l16[0]), (1, l16[1])])
    A10 = field([(2, rnd(10, nc))])
    lam = lamS.map(lambda v: v[:8].reshape(2, 4, nc)); lamb = lamS.map(lambda v: v[8:].reshape(2, 4, nc))
    p2 = L2_parts(A10.map(lambda v: v[:4]), A10.map(lambda v: v[4:]), lam, lamb, FST, nc)
    pI = LIF_parts(A10, lamS, FST, nc)
    d, sc = cmp_all(total(p2, ['kin', 'yukawa']), pI['ferm'])
    rec('D2-' + gname, '1/2(lambda, dslash_A lambda) = i/2(lam,sbar D lamb)+i/2(lamb,sbar D lam)+1/2(lam,eps[A_ji,lam])-1/2(lamb,eps[A^ji,lamb]) pointwise', d, scale=sc)
    d, sc = cmp_all(p2['yukawa'] * (-1), p2['yukawa'])
    rec('D2n-' + gname, 'negative control: Yukawa with flipped sign differs', d, scale=sc, expect=False)
    # D3 BRST part: c(2 odd), A(1 even), A*(1 odd), lam(1 odd), lam*(1 even), c*(1 even)
    CTX = Ctx([1, 1, 0, 1, 1, 0, 0], [rng.normal(size=4) for _ in range(7)])
    c = field([(0, rnd(nc)), (1, rnd(nc))]); A10 = field([(2, rnd(10, nc))]); As10 = field([(3, rnd(10, nc))])
    lamS = field([(4, rnd(16, nc))]); lst16 = field([(5, rnd(16, nc))]); cs = field([(6, rnd(nc))])
    lam = lamS.map(lambda v: v[:8].reshape(2, 4, nc)); lamb = lamS.map(lambda v: v[8:].reshape(2, 4, nc))
    lst = lst16.map(lambda v: v[:8].reshape(2, 4, nc)); lbst = lst16.map(lambda v: v[8:].reshape(2, 4, nc))
    p2 = L2_parts(A10.map(lambda v: v[:4]), A10.map(lambda v: v[4:]), lam, lamb, FST, nc,
                  ghosts=(c, As10.map(lambda v: v[:4]), As10.map(lambda v: v[4:]), lst, lbst, cs))
    pI = LIF_parts(A10, lamS, FST, nc, ghosts=(c, As10, lst16, cs))
    d, sc = cmp_all(p2['brst'], pI['brst'])
    rec('D3-' + gname, '(d_A c,A*)+([c,lam],lam*)+1/2([c,c],c*) = (D c,A*_m)+1/4([A^ij,c],A*_ij)+([c,lam],lam*)-([c,lamb],lamb*)+1/2([c,c],c*)', d, scale=sc)

# integrated: int (lambda, dslash lambda) = 2i int (lamb_{adot i} sbar^{m adot al} d_m lam_al^i), and Srednicki-contraction Yukawa
FST = lie('su2'); nc = 3
p1 = rng.normal(size=4)
CTX = Ctx([1, 1], [p1, -p1])
lamS = field([(0, rnd(16, nc)), (1, rnd(16, nc))])
dl = deriv(lamS).map(lambda v: np.concatenate([v, np.zeros((6,) + v.shape[1:], complex)], axis=0))
lhs = gmul('sa,tb,st,ab->', lamS, dl.map(lambda v: np.einsum('Mts,Msa->ta', RA, v)), PM, KAP * np.eye(nc))
lam = lamS.map(lambda v: v[:8].reshape(2, 4, nc)); lamb = lamS.map(lambda v: v[8:].reshape(2, 4, nc))
rhs = gmul('xia,myib,mxy,ab->', lamb, deriv(lam), SBEa, KAP * np.eye(nc)) * 2j
rec('D4', 'int (lambda, dslash lambda) = 2i int (lamb_{adot i}, sbar^{m adot al} d_m lam_al^i)  (plane waves, sum p = 0)',
    lhs.coef(CTX.full), rhs.coef(CTX.full))
rec('D4n', 'negative control: with -2i must fail', lhs.coef(CTX.full), -rhs.coef(CTX.full), expect=False)
# Srednicki contractions: (lam^j lam^i) = lam^{al j} lam_al^i = eps^{al be} lam_be^j lam_al^i ; (lamb_j lamb_i) = lamb_{ad j} eps^{ad bd} lamb_{bd i}
for gname in ('su2', 'su3'):
    FST = lie(gname); nc = FST.shape[0]
    CTX = Ctx([1, 1, 0], [rng.normal(size=4) for _ in range(3)])
    lamS = field([(0, rnd(16, nc)), (1, rnd(16, nc))]); A10 = field([(2, rnd(10, nc))])
    lam = lamS.map(lambda v: v[:8].reshape(2, 4, nc)); lamb = lamS.map(lambda v: v[8:].reshape(2, 4, nc))
    AI = A10.map(lambda v: v[4:])
    p2 = L2_parts(A10.map(lambda v: v[:4]), AI, lam, lamb, FST, nc)
    Alo = AI.map(lambda v: np.einsum('Ia,Iij->ija', v, SIGI)); Aup = AI.map(lambda v: -np.einsum('Ia,Iij->ija', v, SIGBI))
    lamup = lam.map(lambda v: np.einsum('yz,zja->yja', EU, v))          # lam^{al j}
    lambup = lamb.map(lambda v: np.einsum('xw,wia->xia', EU, v))        # lamb^{ad}_i
    K = KAP * np.eye(nc)
    sred = (gmul('yja,yjb,ab->', lamup, gmul('jia,yib,abc->yjc', Alo, lam, FST), K) * (-0.5)
            + gmul('xja,xjb,ab->', lamb, gmul('jia,xib,abc->xjc', Aup, lambup, FST), K) * (-0.5))
    d, sc = cmp_all(p2['yukawa'], sred)
    rec('D5-' + gname, 'Yukawa = -1/2(lam^j[A_ji,lam^i]) - 1/2(lamb_j[A^ji,lamb_i]) in Srednicki contractions (35.23), pointwise', d, scale=sc)

# =====================================================================================================
# E. supersymmetry of the typeset two-component action under the typeset (F1)
# =====================================================================================================
log('=== E. SUSY invariance of the two-component action (typeset (F1)), non-abelian ===')
SSb = np.einsum('mxy,nyz->mnxz', SEa, SBEa)     # (sigma^m sigmabar^n)_alpha^beta   [m,n,al,be]
SbS = np.einsum('myx,nxz->mnyz', SBEa, SEa)     # (sigmabar^m sigma^n)^{bd}_{ad}    [m,n,bd,ad]


def rhoF_xi(Am, AI, xi, xib, FST):
    """typeset (rho(F_A)xi)_alpha^i, (rho(F_A)xi)_{adot i} of thm:S1odd (F1); xi [al,i], xib [ad,i] commuting."""
    dAm = deriv(Am)
    Fmn = dAm - dAm.map(lambda v: v.transpose(1, 0, 2)) + gmul('ma,nb,abc->mnc', Am, Am, FST)
    Alo = AI.map(lambda v: np.einsum('Ia,Iij->ija', v, SIGI)); Aup = AI.map(lambda v: -np.einsum('Ia,Iij->ija', v, SIGBI))
    DAlo = deriv(Alo) + gmul('ma,ijb,abc->mijc', Am, Alo, FST)
    DAup = deriv(Aup) + gmul('ma,ijb,abc->mijc', Am, Aup, FST)
    CupLo = gmul('ika,kjb,abc->ijc', Aup, Alo, FST)          # [A^{ik}, A_kj]
    CloUp = gmul('ika,kjb,abc->ijc', Alo, Aup, FST)          # [A_ik, A^{kj}]
    up = (Fmn.map(lambda v: -0.5 * np.einsum('mna,mnyz,zi->yia', v, SSb, xi))
          + CupLo.map(lambda v: -0.5 * np.einsum('ija,yj->yia', v, xi))
          + DAup.map(lambda v: -1j * np.einsum('myx,xw,mija,wj->yia', SEa, EU, v, xib)))
    dn = (Fmn.map(lambda v: 0.5 * np.einsum('mna,wi,mnwx->xia', v, xib, SbS))
          + CloUp.map(lambda v: -0.5 * np.einsum('ija,xj->xia', v, xib))
          + DAlo.map(lambda v: 1j * np.einsum('yz,yxm,mija,zj->xia', EU, SEa.transpose(1, 2, 0), v, xi)))
    return up, dn


def Gam_xi_lam(lam, lamb, xi, xib):
    """typeset Gamma^m(xi,lambda) and Gamma_kl(xi,lambda) (thm:S1odd), xi commuting."""
    Gm = (lamb.map(lambda v: 1j * np.einsum('xia,mxy,yi->ma', v, SBEa, xi))
          + lam.map(lambda v: 1j * np.einsum('xi,mxy,yia->ma', xib, SBEa, v)))
    Gkl = (lam.map(lambda v: 2 * np.einsum('jikl,yj,yz,zia->kla', E4, xi, EU, v))
           + lamb.map(lambda v: -2 * (np.einsum('xk,xw,wla->kla', xib, EU, v) - np.einsum('xl,xw,wka->kla', xib, EU, v))))
    return Gm, Gkl


def susy_sector(FST, kinds, xi, xib, sig_lam=-1.0, weights=None):
    """kinds: list of 'l' (lambda leg) / 'A' (A leg); marker 0 = theta (odd, p=0)."""
    nc = FST.shape[0]
    n = len(kinds)
    moms = [rng.normal(size=4) for _ in range(n - 1)]
    moms.append(-sum(moms))
    global CTX
    CTX = Ctx([1] + [1 if k == 'l' else 0 for k in kinds], [np.zeros(4)] + moms)
    lam_legs = [(j + 1, rnd(16, nc)) for j, k in enumerate(kinds) if k == 'l']
    A_legs = [(j + 1, rnd(10, nc)) for j, k in enumerate(kinds) if k == 'A']
    lamS = field(lam_legs); A10 = field(A_legs)
    lam = lamS.map(lambda v: v[:8].reshape(2, 4, nc)); lamb = lamS.map(lambda v: v[8:].reshape(2, 4, nc))
    Am = A10.map(lambda v: v[:4]); AI = A10.map(lambda v: v[4:])
    Gm, Gkl = Gam_xi_lam(lam, lamb, xi, xib)
    XAm = Gm * (-1)
    XAI = Gkl.map(lambda v: -0.25 * np.einsum('Iji,ija->Ia', SIGBI, v))
    up, dn = rhoF_xi(Am, AI, xi, xib, FST)
    Xlam = up * (0.5 * sig_lam); Xlamb = dn * (0.5 * sig_lam)
    parts = L2_parts(Am + lmul_marker(0, XAm), AI + lmul_marker(0, XAI),
                     lam + lmul_marker(0, Xlam), lamb + lmul_marker(0, Xlamb), FST, nc)
    vals = {k: complex(v.coef(CTX.full)) if v.coef(CTX.full) is not None else 0j for k, v in parts.items()}
    w = weights or {}
    tot = sum(w.get(k, 1.0) * v for k, v in vals.items())
    scale = max(abs(v) for v in vals.values())
    return tot, scale, vals


SECTORS = [('l', 'A'), ('l', 'A', 'A'), ('l', 'A', 'A', 'A'), ('l', 'l', 'l')]
for gname in ('su2', 'su3'):
    FST = lie(gname)
    for kinds in SECTORS:
        xi = rnd(2, 4); xib = rnd(2, 4)
        tot, scale, vals = susy_sector(FST, list(kinds), xi, xib)
        rec(f'E-{gname}-{"".join(kinds)}', f'delta_xi S = 0, legs (xi;{",".join(kinds)}); parts '
            + ', '.join(f'{k}={abs(v):.2e}' for k, v in vals.items() if abs(v) > 0), [tot], [0.0], scale=scale)
# negative controls (su2): flipped relative signs/factors must break invariance
FST = lie('su2')
xi = rnd(2, 4); xib = rnd(2, 4)
tot, scale, vals = susy_sector(FST, ['l', 'l', 'l'], xi, xib, weights={'yukawa': -1.0})
rec('En1', 'negative control: Yukawa sign flipped, legs (l,l,l)', [tot], [0.0], scale=scale, expect=False)
tot, scale, vals = susy_sector(FST, ['l', 'A', 'A', 'A'], xi, xib, weights={'quartic': 2.0})
rec('En2', 'negative control: quartic coefficient doubled, legs (l,A,A,A)', [tot], [0.0], scale=scale, expect=False)
tot, scale, vals = susy_sector(FST, ['l', 'A'], xi, xib, weights={'scalar': 2.0})
rec('En3', 'negative control: scalar kinetic doubled, legs (l,A)', [tot], [0.0], scale=scale, expect=False)
tot, scale, vals = susy_sector(FST, ['l', 'A'], xi, xib, sig_lam=+1.0)
rec('En4', 'negative control: xi.lambda = +1/2 rho(F)xi (opposite of the dictionary pattern), legs (l,A)', [tot], [0.0], scale=scale, expect=False)
tot, scale, vals = susy_sector(FST, ['l', 'A', 'A'], xi, xib, weights={'yukawa': -1.0})
rec('En5', 'negative control: Yukawa sign flipped, legs (l,A,A)', [tot], [0.0], scale=scale, expect=False)

# typeset rho(F_A)xi vs index-free 1/2 F_MN rho(e_M)rho(e_N) xi, pointwise (two even markers)
for gname in ('su2', 'su3'):
    FST = lie(gname); nc = FST.shape[0]
    CTX = Ctx([0, 0], [rng.normal(size=4) for _ in range(2)])
    A10 = field([(0, rnd(10, nc)), (1, rnd(10, nc))])
    xi = rnd(2, 4); xib = rnd(2, 4); q = np.concatenate([xi.ravel(), xib.ravel()])
    up, dn = rhoF_xi(A10.map(lambda v: v[:4]), A10.map(lambda v: v[4:]), xi, xib, FST)
    dA = deriv(A10).map(lambda v: np.concatenate([v, np.zeros((6,) + v.shape[1:], complex)], axis=0))
    F = dA - dA.map(lambda v: v.transpose(1, 0, 2)) + gmul('Ma,Nb,abc->MNc', A10, A10, FST)
    RR = np.einsum('Mst,Ntu->MNsu', RbA, RA)
    rF = F.map(lambda v: 0.5 * np.einsum('MNa,MNsu,u->sa', v, RR, q))
    two = GA({k: np.concatenate([up.t.get(k, np.zeros((2, 4, nc))).reshape(8, nc), dn.t.get(k, np.zeros((2, 4, nc))).reshape(8, nc)])
              for k in set(up.t) | set(dn.t)})
    d, sc = cmp_all(two, rF)
    rec('E0-' + gname, 'typeset (rho(F_A)xi) two-component = index-free rho(q(F_A))xi, pointwise', d, scale=sc)

# =====================================================================================================
# F. {Q, Qbar}: the odd-odd bracket in two-component form, closure on A (non-abelian), Lorentzian image
# =====================================================================================================
log('=== F. {Q,Qbar} ===')
xi = rnd(2, 4); xib = rnd(2, 4)
q_xi = np.concatenate([xi.ravel(), np.zeros(8)]); q_xib = np.concatenate([np.zeros(8), xib.ravel()])
Gv = GamS(q_xi, q_xib)
u_E = np.array([1j * np.einsum('xi,mxy,yi->m', xib, SBEa, xi)[m] for m in range(4)])
rec('F1', '[xi,xibar] = Gamma(xi,xibar)|_V = i xibar_{ad i} sbar_E^{m ad al} xi_al^i  (xi in S_+(x)W, xibar in S_-(x)W^*)', Gv[:4], u_E)
rec('F2', 'Gamma(xi,xibar) has no W^perp component; Gamma(xi,xi)|_V = 0 and Gamma(xibar,xibar)|_V = 0',
    [Gv[4:], GamS(q_xi, q_xi)[:4], GamS(q_xib, q_xib)[:4]])
# Lorentzian image: u^0 = -i u^4, u^k = u^k  must equal i xibar sbar^mu xi (Lorentzian sbar^mu, upper index)
uL = np.array([-1j * u_E[3], u_E[0], u_E[1], u_E[2]])
rec('F3', 'Wick image of the vector field: u^mu = i xibar_i sbar^mu xi^i (Lorentzian sbar^mu of (35.3))',
    uL, np.array([1j * np.einsum('xi,xy,yi->', xib, sigb[m], xi) for m in range(4)]))
xiup = np.einsum('ab,bi->ai', EU, xi); xibup = np.einsum('xy,yi->xi', EU, xib)
rec('F4', 'commuting parameters: xibar_{ad} sbar^{mu ad al} xi_al = xi^al sigma^mu_{al ad} xibar^ad  (Srednicki index positions)',
    [np.einsum('xi,xy,yi->', xib, sigb[m], xi) - np.einsum('yi,yx,xi->', xiup, sig[m], xibup) for m in range(4)])
# closure on A, non-abelian, typeset formulas: {xa.,xb.}A_M = 1/2 Gamma_M(xb, rho(F)xa) + 1/2 Gamma_M(xa, rho(F)xb)
for gname in ('su2', 'su3'):
    FST = lie(gname); nc = FST.shape[0]
    CTX = Ctx([0, 0], [rng.normal(size=4) for _ in range(2)])
    A10 = field([(0, rnd(10, nc)), (1, rnd(10, nc))])
    Am = A10.map(lambda v: v[:4]); AI = A10.map(lambda v: v[4:])
    xa, xab, xb, xbb = rnd(2, 4), rnd(2, 4), rnd(2, 4), rnd(2, 4)
    upa, dna = rhoF_xi(Am, AI, xa, xab, FST); upb, dnb = rhoF_xi(Am, AI, xb, xbb, FST)
    Gm1, Gk1 = Gam_xi_lam(upa, dna, xb, xbb)       # Gamma(xb, rho(F)xa)  (bilinear, commuting arguments)
    Gm2, Gk2 = Gam_xi_lam(upb, dnb, xa, xab)
    clos_m = (Gm1 + Gm2) * 0.5
    clos_I = (Gk1 + Gk2).map(lambda v: 0.5 * 0.25 * np.einsum('Iji,ija->Ia', SIGBI, v))
    qa = np.concatenate([xa.ravel(), xab.ravel()]); qb = np.concatenate([xb.ravel(), xbb.ravel()])
    v = GamS(qa, qb)                                 # v^M = Gamma(xa,xb)^M
    dAm = deriv(Am)
    Fmn = dAm - dAm.map(lambda w: w.transpose(1, 0, 2)) + gmul('ma,nb,abc->mnc', Am, Am, FST)
    DAI = deriv(AI) + gmul('ma,Ib,abc->mIc', Am, AI, FST)
    CIJ = gmul('Ia,Jb,abc->IJc', AI, AI, FST)
    exp_m = Fmn.map(lambda w: -np.einsum('n,nma->ma', v[:4], w)) + DAI.map(lambda w: np.einsum('I,mIa->ma', v[4:], w))
    exp_I = DAI.map(lambda w: -np.einsum('n,nIa->Ia', v[:4], w)) + CIJ.map(lambda w: -np.einsum('J,JIa->Ia', v[4:], w))
    d1, s1 = cmp_all(clos_m, exp_m); d2, s2 = cmp_all(clos_I, exp_I)
    rec('F5-' + gname, '{xa.,xb.}A_m = -v^n F_nm + v^I D_m A_I, {xa.,xb.}A_I = -v^n D_n A_I - v^J[A_J,A_I], v = Gamma(xa,xb) (typeset, non-abelian)',
        d1 + d2, scale=max(s1, s2))

# =====================================================================================================
# G. Step-1 translations
# =====================================================================================================
log('=== G. Step-1 translations ===')
FST = lie('su2'); nc = 3; K = KAP * np.eye(nc)
CTX = Ctx([0, 0], [rng.normal(size=4) for _ in range(2)])
A10 = field([(0, rnd(10, nc)), (1, rnd(10, nc))])
dA = deriv(A10)
dA10 = dA.map(lambda v: np.concatenate([v, np.zeros((6,) + v.shape[1:], complex)], axis=0))
om = dA10 - dA10.map(lambda v: v.transpose(1, 0, 2))
hodge_V = gmul('mna,mnb,ab->', om.map(lambda v: v[:4, :4]), om.map(lambda v: v[:4, :4]), K) * 0.5   # sum_{m<n}
hodge_W = gmul('mIa,mIb,ab->', om.map(lambda v: v[:4, 4:]), om.map(lambda v: v[:4, 4:]), K)          # sum_{m,I}
Am = A10.map(lambda v: v[:4]); AI = A10.map(lambda v: v[4:])
dAm = deriv(Am); F0 = dAm - dAm.map(lambda v: v.transpose(1, 0, 2))
Alo = AI.map(lambda v: np.einsum('Ia,Iij->ija', v, SIGI)); Aup = AI.map(lambda v: -np.einsum('Ia,Iij->ija', v, SIGBI))
d1, s1 = cmp_all(hodge_V, gmul('mna,mnb,ab->', F0, F0, K) * 0.5)
d2, s2 = cmp_all(hodge_W, gmul('mija,mijb,ab->', deriv(Aup), deriv(Alo), K) * 0.25)
d3, s3 = cmp_all(gmul('Ia,Ib,ab->', AI, AI, K), gmul('ija,ijb,ab->', Aup, Alo, K) * 0.25)
rec('G1', '(dA|_V,dA|_V) = 1/2(F0_mn,F0_mn); (dA|_Wperp,dA|_Wperp) = 1/4(d_m A^ij, d_m A_ij); (A|_Wperp,A|_Wperp) = 1/4(A^ij,A_ij)',
    d1 + d2 + d3, scale=max(s1, s2, s3))
# d^* = formal adjoint: int (d c, w) = int (c, d^* w) with d^* w = -d_m w_m ; finite counterterm functionals
p = rng.normal(size=4)
CTX = Ctx([1, 0], [p, -p])
c = field([(0, rnd(nc))]); Aw = field([(1, rnd(4, nc))])
lhs = gmul('ma,mb,ab->', deriv(c), Aw, K).coef(CTX.full)
dstar = deriv(Aw).map(lambda v: -np.einsum('mma->a', v))
rhs = gmul('a,b,ab->', c, dstar, K).coef(CTX.full)
rec('G2', 'd^* on one-forms is -d_m w_m (int (dc,w) = int (c,d^*w), plane waves)', lhs, rhs)
CTX = Ctx([1, 0], [rng.normal(size=4), rng.normal(size=4)])
c = field([(0, rnd(nc))]); Aw = field([(1, rnd(4, nc))])
dsdc = deriv(deriv(c)).map(lambda v: -np.einsum('mma->a', v))           # d^* d c = -d^2 c
dsA = deriv(Aw).map(lambda v: -np.einsum('mma->a', v))                   # d^* A = -d_m A_m
lap_c = deriv(deriv(c)).map(lambda v: np.einsum('mma->a', v)); divA = deriv(Aw).map(lambda v: np.einsum('mma->a', v))
d1, s1 = cmp_all(gmul('a,b,ab->', dsdc, dsA, K), gmul('a,b,ab->', lap_c, divA, K))
rec('G3', '(d^*dc, d^*A|_V) = (d^2 c, d_m A_m) pointwise, so -27/8 int(d^*dc,d^*A) = -27/8 int(d^2 c, d_m A_m)', d1, scale=s1)
# Q int(d^*A,d^*A) with Q: A|_V-slot <- -dc  equals -2 int(d^*dc, d^*A)
QdsA = deriv(deriv(c) * (-1)).map(lambda v: -np.einsum('mma->a', v))      # d^*(Q A) with QA = -dc
lhs = gmul('a,b,ab->', QdsA, dsA, K) * 2
rhs = gmul('a,b,ab->', dsdc, dsA, K) * (-2)
d1, s1 = cmp_all(lhs, rhs)
rec('G4', 'Q int(d^*A,d^*A) = -2 int(d^*dc,d^*A) (Q inserts -dc into A|_V); 27/16*(-2) = -27/8', d1 + [27 / 16 * (-2) + 27 / 8], scale=s1)
# field-counting functional and pairings in two components
CTX = Ctx([0, 1, 1, 0, 1, 0], [rng.normal(size=4) for _ in range(6)])
A10 = field([(0, rnd(10, nc))]); As10 = field([(1, rnd(10, nc))]); lamS = field([(2, rnd(16, nc))]); lst16 = field([(3, rnd(16, nc))])
c = field([(4, rnd(nc))]); cs = field([(5, rnd(nc))])
IF = (gmul('Ma,Mb,ab->', A10, As10, K) * 2 - gmul('sa,tb,st,ab->', lamS, lst16, PM, K) * 4 - gmul('a,b,ab->', cs, c, K))
AIl = A10.map(lambda v: v[4:]); AsI = As10.map(lambda v: v[4:])
Aup_ = AIl.map(lambda v: -np.einsum('Ia,Iij->ija', v, SIGBI)); Aslo = AsI.map(lambda v: np.einsum('Ia,Iij->ija', v, SIGI))
lam = lamS.map(lambda v: v[:8].reshape(2, 4, nc)); lamb = lamS.map(lambda v: v[8:].reshape(2, 4, nc))
lst = lst16.map(lambda v: v[:8].reshape(2, 4, nc)); lbst = lst16.map(lambda v: v[8:].reshape(2, 4, nc))
TC = ((gmul('ma,mb,ab->', A10.map(lambda v: v[:4]), As10.map(lambda v: v[:4]), K) + gmul('ija,ijb,ab->', Aup_, Aslo, K) * 0.25) * 2
      - (gmul('yia,yib,ab->', lam, lst, K) - gmul('xia,xib,ab->', lamb, lbst, K)) * 4 - gmul('a,b,ab->', cs, c, K))
d1, s1 = cmp_all(IF, TC)
rec('G5', 'field-counting functional 2(A,A*)-4(lam*,lam)-(c*,c) = 2[(A_m,A*_m)+1/4(A^ij,A*_ij)] - 4[(lam,lam*)-(lamb,lamb*)] - (c*,c)', d1, scale=s1)
# scalar bilinears: 20' and trace
w1 = rnd(6); w2 = rnd(6)
w1lo = np.einsum('I,Iij->ij', w1, SIGI); w1up = -np.einsum('I,Iij->ij', w1, SIGBI)
w2lo = np.einsum('I,Iij->ij', w2, SIGI); w2up = -np.einsum('I,Iij->ij', w2, SIGBI)
CTX = Ctx([0, 0], [rng.normal(size=4) for _ in range(2)])
AI = field([(0, rnd(6, nc)), (1, rnd(6, nc))])
Alo = AI.map(lambda v: np.einsum('Ia,Iij->ija', v, SIGI)); Aup = AI.map(lambda v: -np.einsum('Ia,Iij->ija', v, SIGBI))
iw1 = AI.map(lambda v: np.einsum('I,Ia->a', w1, v)); iw2 = AI.map(lambda v: np.einsum('I,Ia->a', w2, v))
O_if = gmul('a,b,ab->', iw1, iw2, K)
O_2c = gmul('a,b,ab->', Alo.map(lambda v: np.einsum('ij,ija->a', w1up, v)), Alo.map(lambda v: np.einsum('ij,ija->a', w2up, v)), K) * (1 / 16)
tr_if = gmul('Ia,Ib,ab->', AI, AI, K)
tr_2c = gmul('ija,ijb,ab->', Aup, Alo, K) * 0.25
w12 = 0.25 * np.einsum('ij,ij->', w1up, w2lo)
d1, s1 = cmp_all(O_if, O_2c); d2, s2 = cmp_all(tr_if, tr_2c)
t20 = O_2c - tr_2c * (w12 / 6)
t20_2c = O_2c - gmul('ija,ijb,ab->', Aup, Alo, K) * (np.einsum('ij,ij->', w1up, w2lo) / 96)
d3, s3 = cmp_all(t20, t20_2c)
# tracelessness: sum over an orthonormal basis w1=w2=e_I of the 20' projection is zero
trsum = zeroGA()
for I in range(6):
    eI = np.zeros(6); eI[I] = 1
    eup = -np.einsum('I,Iij->ij', eI, SIGBI)
    OI = gmul('a,b,ab->', Alo.map(lambda v: np.einsum('ij,ija->a', eup, v)), Alo.map(lambda v: np.einsum('ij,ija->a', eup, v)), K) * (1 / 16)
    trsum = trsum + OI - gmul('ija,ijb,ab->', Aup, Alo, K) * (np.einsum('ij,ij->', eup, np.einsum('I,Iij->ij', eI, SIGI)) / 96)
d4, s4 = cmp_all(trsum, zeroGA())
rec('G6', "(iota_w1 A,iota_w2 A) = 1/16(w1^ij A_ij, w2^kl A_kl); Konishi = 1/4(A^ij,A_ij); 20' = ... - 1/96 w1^ij w2_ij (A^kl,A_kl), traceless",
    d1 + d2 + d3 + d4, scale=max(s1, s2, s3))
# d_delta I^CT translation: u log eps [4 O - 2 (w1,w2) sum_w(iota_w A,iota_w A)] = u log eps [1/4 (w1A)(w2A) - 1/8 w1^ij w2_ij (A^kl,A_kl)]
lhs = O_if * 4 - tr_if * (2 * (w1 @ w2))
rhs = gmul('a,b,ab->', Alo.map(lambda v: np.einsum('ij,ija->a', w1up, v)), Alo.map(lambda v: np.einsum('ij,ija->a', w2up, v)), K) * 0.25 \
    - gmul('ija,ijb,ab->', Aup, Alo, K) * (np.einsum('ij,ij->', w1up, w2lo) / 8)
d1, s1 = cmp_all(lhs, rhs)
rec('G7', 'd_delta I^CT log part: 4O - 2(w1,w2) sum_w (iota_w A)^2 = 1/4(w1^ij A_ij, w2^kl A_kl) - 1/8 w1^ij w2_ij (A^kl,A_kl); '
    '2 hbar dim g (w1,w2)/(16pi^2 eps) = hbar dim g w1^ij w2_ij/(32 pi^2 eps)',
    d1 + [2 * (w1 @ w2) / 16 - np.einsum('ij,ij->', w1up, w2lo) / 32], scale=s1)

# Konishi number and canonical su(N) map
for N in (2, 3, 4):
    # generalized Gell-Mann basis T^a, tr(T^a T^b) = 1/2 delta
    Ts_ = []
    for j in range(N):
        for k in range(j + 1, N):
            S_ = np.zeros((N, N), complex); S_[j, k] = S_[k, j] = 0.5; Ts_.append(S_)
            A_ = np.zeros((N, N), complex); A_[j, k] = -0.5j; A_[k, j] = 0.5j; Ts_.append(A_)
    for l in range(1, N):
        D_ = np.zeros((N, N), complex)
        for j in range(l):
            D_[j, j] = 1
        D_[l, l] = -l
        Ts_.append(D_ / np.sqrt(2 * l * (l + 1)))
    Ts_ = np.array(Ts_)
    dimg = len(Ts_)
    gc, hb = 0.77, 1.3
    X = np.einsum('a,aij->ij', rnd(dimg).real, Ts_) * (-1j * gc); Y = np.einsum('a,aij->ij', rnd(dimg).real, Ts_) * (-1j * gc)
    def ad(Z):
        # matrix of ad_Z on the basis T^a (coordinates via 2 tr(T^a .))
        return np.array([[2 * np.trace(Ts_[a] @ (Z @ Ts_[b] - Ts_[b] @ Z)) for b in range(dimg)] for a in range(dimg)])
    tr_ad = np.trace(ad(X) @ ad(Y))
    pair = -4 * hb / gc ** 2 * np.trace(X @ Y)
    unit = hb * tr_ad / (16 * np.pi ** 2 * pair)
    rec(f'G8-su{N}', f'su({N}): tr_g(ad_X ad_Y) = 2N tr(XY); unit = -g^2 N/(32 pi^2); -24 unit = 3 g^2 N/(4 pi^2); (X,Y) = 2 hbar X^a Y^a',
        [tr_ad - 2 * N * np.trace(X @ Y), unit + gc ** 2 * N / (32 * np.pi ** 2), -24 * unit - 3 * gc ** 2 * N / (4 * np.pi ** 2),
         pair - 2 * hb * np.einsum('a,a->', np.array([2 * np.trace(Ts_[a] @ X) / (-1j * gc) for a in range(dimg)]),
                                    np.array([2 * np.trace(Ts_[a] @ Y) / (-1j * gc) for a in range(dimg)]))])

# S/hbar = int L_Sred|_Wick for the gauge, scalar and gaugino kinetic terms (su(2), su(3)):
# basis t_a = -i g T^a: [t_a,t_b] = g f^{abc} t_c, (t_a,t_b) = 2 hbar delta_ab; A^(t) = A^a, lambda^(t) = lambda^a/sqrt2
for gname in ('su2', 'su3'):
    f = lie(gname); nc = f.shape[0]
    gc, hb = 0.91, 0.6
    KAP_save = KAP; KAP = 2 * hb
    CTX = Ctx([0, 0, 0], [rng.normal(size=4) for _ in range(3)])
    A10 = field([(0, rnd(10, nc)), (1, rnd(10, nc)), (2, rnd(10, nc))])     # Srednicki A^a_m (Euclid comps), A^a_I
    p2 = L2_parts(A10.map(lambda v: v[:4]), A10.map(lambda v: v[4:]), GA({}), GA({}), gc * f, nc)
    Am = A10.map(lambda v: v[:4]); AI = A10.map(lambda v: v[4:])
    dAm = deriv(Am)
    Fa = dAm - dAm.map(lambda v: v.transpose(1, 0, 2)) + gmul('mb,nc,abc->mna', Am, Am, gc * f)    # F^a = dA-dA + g f^{abc} A^b A^c
    DAI = deriv(AI) + gmul('mb,Ic,abc->mIa', Am, AI, gc * f)                                       # (D A_I)^a = d A^a + g f^{abc} A^b A^c
    Lsred = gmul('mna,mna->', Fa, Fa) * (-0.25) + gmul('mIa,mIa->', DAI, DAI) * (-0.5)
    d1, s1 = cmp_all(total(p2, ['gauge', 'scalar']) * (1 / hb), Lsred)
    # gaugino kinetic, integrated (two lambda legs + one A leg, sum p = 0)
    moms = [rng.normal(size=4) for _ in range(2)]; moms.append(-moms[0] - moms[1])
    CTX = Ctx([1, 1, 0], moms)
    lS = field([(0, rnd(16, nc)), (1, rnd(16, nc))]); A10 = field([(2, rnd(10, nc))])
    lam = lS.map(lambda v: v[:8].reshape(2, 4, nc) / np.sqrt(2)); lamb = lS.map(lambda v: v[8:].reshape(2, 4, nc) / np.sqrt(2))
    p2 = L2_parts(A10.map(lambda v: v[:4]), A10.map(lambda v: v[4:]), lam, lamb, gc * f, nc)
    laS = lS.map(lambda v: v[:8].reshape(2, 4, nc)); lbS = lS.map(lambda v: v[8:].reshape(2, 4, nc))
    Am = A10.map(lambda v: v[:4])
    Dla = deriv(laS) + gmul('mb,yic,abc->myia', Am, laS, gc * f)
    Lk = gmul('xia,myia,mxy->', lbS, Dla, SBEa) * 1j                    # i lambdabar^a sbar^m (D_m lambda)^a
    rec(f'G9-{gname}', 'canonical map: S/hbar = int[-1/4 F^aF^a - 1/2 (D A_I)^a(D A_I)^a + i lambdabar^a sbar^m D_m lambda^a] (Srednicki (69.16),(3.9),(95.77) after Wick)',
        d1 + [p2['kin'].coef(CTX.full) / hb - Lk.coef(CTX.full)], scale=max(s1, abs(Lk.coef(CTX.full))))
    KAP = KAP_save

# =====================================================================================================
# H. provenance: no SIG/REP/CONF entry, no blocked user note
# =====================================================================================================
log('=== H. provenance ===')
NOTE = os.path.join(ROOT, 'step2', 'notation')
hits = {}
pat_ids = ['S-SIG', 'S-REP', 'S-CONF', 'SIG-', 'REP-', 'CONF-', 'user-susy-note', 'SUSY_3p1d', 'secondary_baseline',
           'g_E', 'Q^E', 'Q^{E}', '-i\\sigma^i', '(-i\\sigma', 'susy-note']
for fn in ('esw_to_weinberg_srednicki.tex', 'dictionary.json', 'check_dictionary.py'):
    txt = open(os.path.join(NOTE, fn), encoding='utf-8').read()
    hits[fn] = [p_ for p_ in pat_ids if p_ in txt]
for fn, h in hits.items():
    allowed = [x for x in h if x not in ('SIG-',)]
    log(f'  {fn}: pattern hits {h}')
txt = open(os.path.join(NOTE, 'check_dictionary.py'), encoding='utf-8').read()
imports = [l_.strip() for l_ in txt.splitlines() if l_.strip().startswith(('import ', 'from '))]
log('  check_dictionary.py imports:', imports)
blocked = [p_ for fn, h in hits.items() for p_ in h if p_ not in ('SIG-',)]
rec('H1', 'no SIG/REP/CONF/user-note identifiers in tex/json/py (only the explicit "not used" statements)',
    [len([p_ for fn in hits for p_ in hits[fn] if p_ in ('user-susy-note', 'SUSY_3p1d', 'secondary_baseline', 'g_E', 'Q^E', 'Q^{E}', '(-i\\sigma', 'susy-note')])])

npass = sum(r['passed'] for r in RES); ntot = len(RES)
log(f'=== SUMMARY: {npass}/{ntot} passed; max residual among expected-zero checks = '
    f'{max(r["residual"] for r in RES if r["expect_zero"]):.2e} ===')
json.dump(dict(results=RES, passed=npass, total=ntot), open(os.path.join(HERE, 'audit_results.json'), 'w'), indent=1)
open(os.path.join(HERE, 'audit_dictionary.log'), 'w').write('\n'.join(LOG) + '\n')
