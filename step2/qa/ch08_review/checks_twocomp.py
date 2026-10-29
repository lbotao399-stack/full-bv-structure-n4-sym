#!/usr/bin/env python3
"""Two-component checks of ch08 Def c8-twocomp, Prop c8-twocomp, Lemma c8-cc, Thm c8-dag (i), Def c8-Theta,
Def c8-B / Prop c8-B (block vanishing), Lemma c8-ev (rank one).  The realisation is the dictionary one
((R1)-(R4), (P1)), rebuilt here from the formulas; 't Hooft blocks from engine B (read-only import)."""
import os, sys, json, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, os.path.join(ROOT, 'engine_b'))
import realization as RB  # noqa: E402

np.seterr(all='ignore')
rng = np.random.default_rng(11)
OUT = {}


def rec(key, ok, **info):
    OUT[key] = dict(passed=bool(ok), **info)
    print(('PASS ' if ok else 'FAIL ') + key, info if info else '')


I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EPSU = np.array([[0, 1], [-1, 0]], complex)
EPSD = np.array([[0, -1], [1, 0]], complex)
SIGE = [PAULI[0], PAULI[1], PAULI[2], 1j * I2]           # sigma^m_{alpha alphadot}   (ch08 Def c8-twocomp)
SIGBE = [-PAULI[0], -PAULI[1], -PAULI[2], 1j * I2]       # sigmabar^{m alphadot alpha}
SIGI = [np.array(S.tolist(), dtype=complex) for S in RB.Sig]
SIGBI = [np.array(S.tolist(), dtype=complex) for S in RB.Sigb]


def split(x):
    return x[:8].reshape(2, 4), x[8:].reshape(2, 4)


def join(a, b):
    return np.concatenate([a.reshape(8), b.reshape(8)])


def rho_S(M, x):          # Sigma -> Sigma^*   (R1), (R3)
    q, qb = split(x)
    if M < 4:
        s = 1j * np.einsum('di,da->ai', qb, SIGBE[M]); sb = -1j * np.einsum('da,ai->di', SIGBE[M], q)
    else:
        I = M - 4
        s = np.einsum('ab,ji,bi->aj', EPSU, SIGI[I], q); sb = -np.einsum('ji,de,ei->dj', SIGBI[I], EPSU, qb)
    return join(s, sb)


def rho_Ss(M, y):         # Sigma^* -> Sigma   (R2), (R4)
    s, sb = split(y)
    if M < 4:
        q = -1j * np.einsum('ad,di->ai', SIGE[M], sb); qb = 1j * np.einsum('ai,ad->di', s, SIGE[M])
    else:
        I = M - 4
        q = np.einsum('kj,ab,bj->ak', SIGBI[I], EPSD, s); qb = -np.einsum('kj,ed,dj->ek', SIGI[I], EPSD, sb)
    return join(q, qb)


def pair(x, y):           # (xi, eta) = xi_alpha^i eta^alpha_i - xibar_{alphadot i} etabar^{alphadot i}   (P1)
    q, qb = split(x); s, sb = split(y)
    return np.sum(q * s) - np.sum(qb * sb)


E16 = np.eye(16)
R = [np.array([rho_S(M, E16[:, k]) for k in range(16)]).T for M in range(10)]     # Sigma -> Sigma^*
Rb = [np.array([rho_Ss(M, E16[:, k]) for k in range(16)]).T for M in range(10)]   # Sigma^* -> Sigma
P = np.array([[pair(E16[:, a], E16[:, b]) for b in range(16)] for a in range(16)])

rec('T1_Clifford_all_ten_directions',
    max(np.abs(Rb[M] @ R[N] + Rb[N] @ R[M] - 2 * (M == N) * E16).max() for M in range(10) for N in range(10)) < 1e-12
    and max(np.abs(R[M] @ Rb[N] + R[N] @ Rb[M] - 2 * (M == N) * E16).max() for M in range(10) for N in range(10)) < 1e-12)
rec('T2_symmetry_(rho(v)Q1,Q2)=(Q1,rho(v)Q2)', max(np.abs(P @ R[M] - (P @ R[M]).T).max() for M in range(10)) < 1e-12
    and max(np.abs(Rb[M].T @ P - (Rb[M].T @ P).T).max() for M in range(10)) < 1e-12)
vol = Rb[0] @ R[1] @ Rb[2] @ R[3]
rec('T3_rho(e1e2e3e4)=+1_on_S+W_(undotted)', np.allclose(vol[:8, :8], np.eye(8)) and np.allclose(vol[8:, 8:], -np.eye(8)))


# (R1) is fixed by (R2) and the symmetry of the pairing: (rho(v)xi, eta) = (xi, rho(v)eta)
# here (rho(v)xi, eta) means pair(xi', ...) with xi' in Sigma^*: use (a,b) for a in Sigma^*, b in Sigma := pair(b,a)
def R1_from_R2(M):
    # solve (rho(v) e_a, e_b)_{Sigma^* x Sigma} = (e_a, rho(v) e_b): pair(e_b, X e_a) = pair(e_a, Rb e_b)
    A = np.zeros((16, 16), complex)
    # unknown X (16x16): pair(e_b, X e_a) = sum_c P[b,c] X[c,a]  =>  P X = (P Rb)^T
    return np.linalg.solve(P, (P @ Rb[M]).T)


rec('T4_(R1)_is_determined_by_(R2)_through_the_Clifford_relation_rho(e_m)^2=1',
    max(np.abs(np.linalg.inv(Rb[M]) - R[M]).max() for M in range(4)) < 1e-12)
# but rho on W_perp (R3),(R4) is extra data, not fixed by (R2): exhibit a second Clifford module structure agreeing on V
# flip the sign of rho on W_perp: Clifford relation still holds, pairing symmetric, chirality same -> different Gamma|_W
Rf = [R[M] if M < 4 else -R[M] for M in range(10)]; Rbf = [Rb[M] if M < 4 else -Rb[M] for M in range(10)]
ok = max(np.abs(Rbf[M] @ Rf[N] + Rbf[N] @ Rf[M] - 2 * (M == N) * E16).max() for M in range(10) for N in range(10)) < 1e-12
rec('T5_V-data_of_Def_c8-twocomp_do_not_fix_rho_on_W_perp_(sign-flipped_W_perp_part_is_another_realisation)', ok,
    note='rho(e_I) and Gamma|_W are not given in Def c8-twocomp; they enter S^(1)(xi), S^(2) and ESW Assumption 3.1')


# Gamma
def Gam(Q1, Q2):          # Gamma^M = (rho(e_M) Q1, Q2) for Q1,Q2 in Sigma: pair(Q2, rho_S Q1)
    return np.array([pair(Q2, rho_S(M, Q1)) for M in range(10)])


def Gams(S1, S2):         # Sigma^* x Sigma^*: (rho(e_M) S1, S2) = pair(rho_Ss S1, S2)
    return np.array([pair(rho_Ss(M, S1), S2) for M in range(10)])


def rnd(n):
    return rng.normal(size=n) + 1j * rng.normal(size=n)


x1, x2 = rnd(16), rnd(16)
q1, qb1 = split(x1); q2, qb2 = split(x2)
G_print = np.array([1j * (np.einsum('di,da,ai->', qb2, SIGBE[m], q1) + np.einsum('di,da,ai->', qb1, SIGBE[m], q2)) for m in range(4)])
rec('T6_Gamma(xi1,xi2)|_V_as_printed', np.allclose(Gam(x1, x2)[:4], G_print))
y1, y2 = rnd(16), rnd(16)
s1, sb1 = split(y1); s2, sb2 = split(y2)
Gs_print = np.array([-1j * (np.einsum('ai,ad,di->', s2, SIGE[m], sb1) + np.einsum('ai,ad,di->', s1, SIGE[m], sb2)) for m in range(4)])
rec('T7_Gamma(eta1,eta2)|_V_as_printed', np.allclose(Gams(y1, y2)[:4], Gs_print))
# Killing spinor (c8-killing)
xv = rng.normal(size=4)
eps_idx = sum(xv[m] * rho_Ss(m, y1) for m in range(4)) + x1
e, eb = split(eps_idx)
e_print = q1 - 1j * np.einsum('m,mad,di->ai', xv, np.array(SIGE), sb1)
eb_print = qb1 + 1j * np.einsum('m,ai,mad->di', xv, s1, np.array(SIGE))
rec('T8_Killing_spinor_(c8-killing)', np.allclose(e, e_print) and np.allclose(eb, eb_print))


# ESW Assumption 3.1
def assum(Q1, Q2, Q3):
    out = 0
    for a, b, c in ((Q1, Q2, Q3), (Q2, Q3, Q1), (Q3, Q1, Q2)):
        g = Gam(a, b)
        out = out + sum(g[M] * rho_S(M, c) for M in range(10))
    return out


rec('T9_ESW_Assumption_3.1', np.abs(assum(rnd(16), rnd(16), rnd(16))).max() < 1e-10)
# with the W_perp sign flipped, Assumption 3.1 still holds? (it should: it's the same Clifford module up to e_I -> -e_I)


# ---------------- Lemma c8-cc conjugation
def conjS(x):     # Sigma -> Sigma^*: eta = conj(xi), etabar = -conj(xibar)
    q, qb = split(x); return join(q.conj(), -qb.conj())


def conjSs(y):    # Sigma^* -> Sigma
    s, sb = split(y); return join(s.conj(), -sb.conj())


ok = True
for M in range(10):
    for _ in range(3):
        x = rnd(16); y = rnd(16)
        ok &= np.allclose(conjS(rho_Ss(M, y)) if False else conjSs(rho_S(M, x)) * 0 + conjSs(rho_S(M, x)), rho_Ss(M, conjS(x)))
        ok &= np.allclose(conjS(rho_Ss(M, y)), rho_S(M, conjSs(y)))
rec('T10_conj(rho(e_M)Q)=rho(e_M)conj(Q)_all_ten_real_e_M', ok)
x = rnd(16)
rec('T11_(xi,conj xi)=sum|xi|^2>0', np.isclose(pair(x, conjS(x)), np.sum(np.abs(x) ** 2)))
rec('T12_involutive', np.allclose(conjSs(conjS(x)), x))
# Gamma(conj Q3, conj Q4) = conj Gamma(Q3,Q4) for Q3,Q4 in Sigma (conj lands in Sigma^*)
a, b = rnd(16), rnd(16)
rec('T13_Gamma(conj Q3,conj Q4)=conj Gamma(Q3,Q4)', np.allclose(Gams(conjS(a), conjS(b)), Gam(a, b).conj()))
# uniqueness gap: c_a = a*conj on Sigma, (1/a)*conj on Sigma^*: involutive and positive for a>0, Clifford-compatible only at a=1
aa = 2.0
cS = lambda v: aa * conjS(v); cSs = lambda v: conjSs(v) / aa
inv_ok = np.allclose(cSs(cS(x)), x); pos_ok = pair(x, cS(x)).real > 0
cl_ok = np.allclose(cSs(rho_S(0, x)), rho_Ss(0, cS(x)))
rec('T14_uniqueness:_involutivity+positivity_alone_allow_c_a_(a=2);_Clifford_compatibility_excludes_it',
    inv_ok and pos_ok and not cl_ok, involutive=bool(inv_ok), positive=bool(pos_ok), clifford_compatible=bool(cl_ok))


# ---------------- Thm c8-dag (i) two-component
def dagS(x):      # xi^dag = -conj(xi) in Sigma^*
    return -conjS(x)


q, qb = split(x)
s, sb = split(dagS(x))
rec('T15_(xi^dag)^alpha_i=-conj(xi_alpha^i),(xi^dag)^{alphadot i}=+conj(xibar)', np.allclose(s, -q.conj()) and np.allclose(sb, qb.conj()))
# unit elements
xi_unit = np.zeros(16, complex); xi_unit[0] = 1   # Q^1_1 (alpha=1,i=1)
s, sb = split(dagS(xi_unit))
xb_unit = np.zeros(16, complex); xb_unit[8] = 1   # Qbar_{1 1}
s2, sb2 = split(dagS(xb_unit))
rec('T16_(Q^i_alpha)^dag=-S^alpha_i,(Qbar)^dag=+Sbar', np.isclose(s[0, 0], -1) and np.isclose(sb2[0, 0], 1))
# ---------------- Def c8-Theta: xi^{dag_n} = rho(n) conj(xi) two-component
nv = rng.normal(size=4); nv /= np.linalg.norm(nv)
xin = sum(nv[m] * rho_Ss(m, conjS(x)) for m in range(4))
qn, qbn = split(xin)
qn_print = 1j * np.einsum('m,mad,di->ai', nv, np.array(SIGE), qb.conj())
qbn_print = 1j * np.einsum('m,ai,mad->di', nv, q.conj(), np.array(SIGE))
rec('T17_(xi^{dag_n})_two-component_as_printed', np.allclose(qn, qn_print) and np.allclose(qbn, qbn_print))


# ---------------- B: brackets inside one B-eigenspace vanish (Def c2-bracket realised with the dictionary data)
def odd_bracket(Q1, Q2, Q3, Q4, xpt):
    """V-part at x of [Q1+rho(x)Q2, Q3+rho(x)Q4] and the W_perp rotation matrix (w,Rw')."""
    rx = lambda y: sum(xpt[m] * rho_Ss(m, y) for m in range(4))
    e1 = Q1 + rx(Q2); e2 = Q3 + rx(Q4)
    V = Gam(e1, e2)[:4]
    Rw = np.zeros((6, 6), complex)
    for I, J in itertools.product(range(6), repeat=2):
        if I == J:
            continue
        # rho(q(w ^ w')) Q = rho(w) rho(w') Q for orthogonal w,w'  (Sigma -> Sigma)
        qww = lambda Q: rho_Ss(4 + I, rho_S(4 + J, Q))
        Rw[I, J] = pair(qww(Q1), Q4) + pair(qww(Q3), Q2)
    return V, Rw


def proj(v, which):
    a, b = split(v)
    if which == 'undotted':
        return join(a, 0 * b)
    return join(0 * a, b)


zero = np.zeros(16, complex)
xi_p = proj(rnd(16), 'undotted')      # S+ (x) W in Sigma           B = +1
xi_m = proj(rnd(16), 'dotted')        # S- (x) W^* in Sigma         B = -1
eta_m = proj(rnd(16), 'undotted')     # S+ (x) W^* in Sigma^*       B = -1
eta_p = proj(rnd(16), 'dotted')       # S- (x) W in Sigma^*         B = +1
xpt = rng.normal(size=4)
maxin = 0.0
for (A1, A2), (B1, B2) in [((xi_p, zero), (proj(rnd(16), 'undotted'), zero)), ((xi_p, zero), (zero, eta_p)),
                           ((zero, eta_p), (zero, proj(rnd(16), 'dotted'))), ((xi_m, zero), (proj(rnd(16), 'dotted'), zero)),
                           ((xi_m, zero), (zero, eta_m)), ((zero, eta_m), (zero, proj(rnd(16), 'undotted')))]:
    V, Rw = odd_bracket(A1, A2, B1, B2, xpt)
    maxin = max(maxin, np.abs(V).max(), np.abs(Rw).max())
V, Rw = odd_bracket(xi_p, zero, zero, eta_m, xpt)
cross = max(np.abs(V).max(), np.abs(Rw).max())
rec('T18_brackets_inside_each_B-eigenspace_vanish_(B_is_a_derivation_on_odd-odd)', maxin < 1e-10 and cross > 1e-3,
    max_inside=float(maxin), cross_example=float(cross))
# ---------------- Lemma c8-ev: rank one
u = rnd(2); w = rnd(4)
xi1 = join(np.outer(u, w), np.zeros((2, 4)))
rec('T19_rank_one_Gamma(xi,xi)=0_all_ten', np.abs(Gam(xi1, xi1)).max() < 1e-12)
xi2 = proj(rnd(16), 'undotted')
g2 = Gam(xi2, xi2)
rec('T20_generic_S+W:_Gamma|_V=0_but_Gamma|_W_nonzero', np.abs(g2[:4]).max() < 1e-12 and np.abs(g2[4:]).max() > 1e-3,
    W_part=float(np.abs(g2[4:]).max()))
# dilatation part of [xi, xi^dag] = (1/2)(xi, -conj xi) Delta: sign
rec('T21_(xi,xi^dag)=-sum|xi|^2<0', pair(x, dagS(x)).real < 0)
json.dump(OUT, open(os.path.join(HERE, 'checks_twocomp.json'), 'w'), indent=1)
print('wrote checks_twocomp.json')
