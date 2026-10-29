#!/usr/bin/env python3
"""Audit ip_audit2 (algebra and signs): the field realisation of the conjugations.

Independent implementation (no project code imported) of
  * the BV fields of ch01 Def c1-fields in the two-component realisation of ch08 Def c8-twocomp
    (V_10 index M = 0..9, spinors as 32-vectors Sigma (+) Sigma^*), gauge algebra su(2) with real basis,
    [e_a, e_b] = eps_abc e_c, (e_a, e_b) = delta_ab (positive on the compact form);
  * the densities of S (ch01 eq:c1-action), S^(1) (ch02 eq:c2-S1even, eq:c2-S1Q, eq:c2-S1S), S^(2) (eq:c2-S2);
  * the reflection Theta_n of ch08 Def c8-Theta and the radial reflection Theta_I of Prop c8-C
    (factors as listed in step2/theory/conjugation.tex Def Theta), with F^Theta(e) = F(Theta e)^c,
    c = complex conjugation of coefficients and reversal of the order of odd factors.
Checked (pointwise in the densities, ch08 Thm c8-fields (i), Prop c8-C):
  S^Theta = S,  S^(1)(X)^Theta = (-1)^|X| S^(1)(X^dag),  S^(2)(X,Y)^Theta = S^(2)(Y^dag, X^dag),
with the controls: opposite odd sign in dagger_n, opposite sign of s_lambda.
Odd field values are Grassmann-linear (coefficient vectors over NG generators); all densities are at most
bilinear in odd fields, so a Grassmann element is (degree 0 scalar, degree 1 vector, degree 2 antisymmetric
matrix).  Output: a2_results.json.
"""
import json, os, itertools, resource
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(424242)
NG = 5          # number of odd generators

# ------------------------------------------------------------------ realisation (Def c8-twocomp)
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
sig = [s1, s2, s3, 1j * np.eye(2)]; sigb = [-s1, -s2, -s3, 1j * np.eye(2)]
epsU = np.array([[0, 1], [-1, 0]], complex); epsL = -epsU
eps3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    eps3[p] = np.linalg.det(np.eye(3)[list(p)])
eps4 = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(4)):
    eps4[p] = np.linalg.det(np.eye(4)[list(p)])
Sig = np.zeros((6, 4, 4), complex)
for a in range(3):
    for b in range(3):
        for c in range(3):
            Sig[a, b, c] = eps3[a, b, c]; Sig[a + 3, b, c] = 1j * eps3[a, b, c]
        Sig[a, b, 3] = 1.0 * (a == b); Sig[a, 3, b] = -1.0 * (a == b)
        Sig[a + 3, b, 3] = -1j * (a == b); Sig[a + 3, 3, b] = 1j * (a == b)
def ix(block, a, i):
    return {'xi': 0, 'xib': 8, 'eta': 16, 'etab': 24}[block] + 4 * a + i
def rho_V(v):
    M = np.zeros((32, 32), complex)
    vs = sum(v[m] * sig[m] for m in range(4)); vsb = sum(v[m] * sigb[m] for m in range(4))
    for i in range(4):
        for a in range(2):
            for ad in range(2):
                M[ix('xi', a, i), ix('etab', ad, i)] += -1j * vs[a, ad]
                M[ix('xib', ad, i), ix('eta', a, i)] += 1j * vs[a, ad]
                M[ix('eta', a, i), ix('xib', ad, i)] += 1j * vsb[ad, a]
                M[ix('etab', ad, i), ix('xi', a, i)] += -1j * vsb[ad, a]
    return M
def rho_W(w):
    M = np.zeros((32, 32), complex)
    wl = np.einsum('I,Iij->ij', w, Sig); wu = 0.5 * np.einsum('ijkl,kl->ij', eps4, wl)
    for a in range(2):
        for b in range(2):
            for j in range(4):
                for i in range(4):
                    M[ix('eta', a, j), ix('xi', b, i)] += epsU[a, b] * wl[j, i]
                    M[ix('etab', a, j), ix('xib', b, i)] += epsU[a, b] * wu[j, i]
                    M[ix('xi', a, j), ix('eta', b, i)] += -wu[j, i] * epsL[a, b]
                    M[ix('xib', a, j), ix('etab', b, i)] += -wl[j, i] * epsL[a, b]
    return M
RHO = [rho_V(np.eye(4)[m]) for m in range(4)] + [rho_W(np.eye(6)[I]) for I in range(6)]
Gp = np.zeros((16, 16)); Gp[:8, :8] = np.eye(8); Gp[8:, 8:] = -np.eye(8)
Bf = np.zeros((32, 32)); Bf[:16, 16:] = Gp; Bf[16:, :16] = Gp.T
def cc(z):          # Lemma c8-cc, antilinear, acts on the first axis (spinor) of arrays
    out = np.zeros_like(z)
    out[16:24] = np.conj(z[0:8]); out[24:32] = -np.conj(z[8:16])
    out[0:8] = np.conj(z[16:24]); out[8:16] = -np.conj(z[24:32])
    return out
def rhox(x):
    return sum(x[m] * RHO[m] for m in range(4))

# ------------------------------------------------------------------ Grassmann bookkeeping
class G:
    def __init__(self, s=0j, l=None, q=None):
        self.s = complex(s); self.l = np.zeros(NG, complex) if l is None else l; self.q = np.zeros((NG, NG), complex) if q is None else q
    def __add__(self, o):
        return G(self.s + o.s, self.l + o.l, self.q + o.q)
    def __sub__(self, o):
        return G(self.s - o.s, self.l - o.l, self.q - o.q)
    def scale(self, a):
        return G(a * self.s, a * self.l, a * self.q)
    def c(self):     # conjugate coefficients, reverse order of odd factors (degree 2 -> sign -1)
        return G(np.conj(self.s), np.conj(self.l), -np.conj(self.q))
    def norm(self):
        return max(abs(self.s), np.max(np.abs(self.l)), np.max(np.abs(self.q)))
def bil(U, Wmat, V):
    """sum_{alpha,beta} U_alpha W_{alpha beta} V_beta for degree-1 U (shape (n,NG)) and V (shape (m,NG)),
    U first: returns the degree-2 element with antisymmetric matrix."""
    Mq = U.T @ Wmat @ V
    return G(0, None, Mq - Mq.T)

# ------------------------------------------------------------------ densities
def Fstrength(A, dA):
    F = np.zeros((10, 10, 3), complex)
    for M in range(10):
        for N in range(10):
            d = (dA[M, N] if M < 4 else 0) - (dA[N, M] if N < 4 else 0)
            F[M, N] = d + np.einsum('abc,b,c->a', eps3, A[M], A[N])
    return F
def S_terms(f):
    A, dA, lam, dlam, c, dc, As, ls, cs = f['A'], f['dA'], f['lam'], f['dlam'], f['c'], f['dc'], f['As'], f['ls'], f['cs']
    F = Fstrength(A, dA)
    t = {}
    t['YM'] = G(-0.25 * sum(np.sum(F[M, N] ** 2) for M in range(10) for N in range(M + 1, 10)))
    # scalar kinetic part separately (inside YM): (d_A A_W, d_A A_W)
    t['YM_scalar_kin'] = G(-0.25 * sum(np.sum(F[m, 4 + I] ** 2) for m in range(4) for I in range(6)))
    # 1/2 (lambda, dslash_A lambda): Dl[s,a,j]
    Dl = np.zeros((32, 3, NG), complex)
    for M in range(10):
        inner = (dlam[M] if M < 4 else 0) + np.einsum('abc,b,scj->saj', eps3, A[M], lam)
        Dl += np.einsum('st,taj->saj', RHO[M], inner)
    U = lam.reshape(32 * 3, NG); V = Dl.reshape(32 * 3, NG)
    Wm = np.kron(Bf, np.eye(3))
    t['Dirac'] = bil(U, Wm, V).scale(0.5)
    # (d_A c, A^*)
    Dc = np.zeros((10, 3, NG), complex)
    for M in range(10):
        Dc[M] = (dc[M] if M < 4 else 0) + np.einsum('abc,b,cj->aj', eps3, A[M], c)
    t['ghostA'] = bil(Dc.reshape(30, NG), np.eye(30), As.reshape(30, NG))
    # ([c, lambda], lambda^*):  sum eps_abc c^b (lambda^c, lambda*^a)
    Wcl = np.einsum('abc,st,ta->bsc', eps3, Bf, ls).reshape(3, 32 * 3)
    t['ghostLam'] = bil(c, Wcl, lam.reshape(32 * 3, NG))
    # 1/2 ([c,c], c^*)
    Wcc = np.einsum('abc,a->bc', eps3, cs)
    t['ghostGhost'] = bil(c, Wcc, c).scale(0.5)
    return t
def S1_even(f, v, J, Rw):
    """S^(1) of the even element (vector field v on V with Jacobian J at the point, plus Rw in so(W_R^perp))."""
    A, dA, lam, dlam, c, dc, As, ls, cs = f['A'], f['dA'], f['lam'], f['dlam'], f['c'], f['dc'], f['As'], f['ls'], f['cs']
    div = np.trace(J)
    LvA = np.zeros((10, 3), complex)
    for m in range(4):
        LvA[m] = sum(v[k] * dA[k, m] for k in range(4)) + sum(A[k] * J[k, m] for k in range(4))
    for I in range(6):
        LvA[4 + I] = sum(v[k] * dA[k, 4 + I] for k in range(4)) + sum(Rw[J_, I] * A[4 + J_] for J_ in range(6)) + 0.25 * div * A[4 + I]
    E1 = G(0, np.einsum('Ma,Maj->j', LvA, As))
    q2 = np.zeros((32, 32), complex)
    for l in range(4):
        for k in range(l + 1, 4):
            q2 += (J[k, l] - J[l, k]) * (RHO[l] @ RHO[k])
    for I in range(6):
        for J_ in range(I + 1, 6):
            q2 += -2 * Rw[I, J_] * (RHO[4 + I] @ RHO[4 + J_])
    vl = sum(v[k] * dlam[k] for k in range(4)) + 0.25 * np.einsum('st,taj->saj', q2, lam) + 0.375 * div * lam
    E2 = G(0, -np.einsum('saj,st,ta->j', vl, Bf, ls))
    vc = sum(v[k] * dc[k] for k in range(4))
    E3 = G(0, -np.einsum('aj,a->j', vc, cs))
    return E1 + E2 + E3
def S1_odd(f, eps, Q2):
    A, dA, lam, As, ls = f['A'], f['dA'], f['lam'], f['As'], f['ls']
    F = Fstrength(A, dA)
    # -(Gamma(eps, lambda), A^*): Gamma^M(eps, lambda^a) = (rho_M eps, lambda^a)
    Gm = np.array([np.einsum('s,st,taj->aj', RHO[M] @ eps, Bf, lam) for M in range(10)])   # (M, a, j)
    T1 = bil(Gm.reshape(30, NG), -np.eye(30), As.reshape(30, NG))
    rF = np.zeros((3, 32), complex)
    for M in range(10):
        for N in range(M + 1, 10):
            rF += np.einsum('a,s->as', F[M, N], (RHO[M] @ RHO[N]) @ eps)
    T2 = G(0.5 * np.einsum('as,st,ta->', rF, Bf, ls))
    rAQ = np.einsum('Ia,Ist,t->as', A[4:], np.array(RHO[4:]), Q2)
    T3 = G(-np.einsum('as,st,ta->', rAQ, Bf, ls))
    return T1 + T2 + T3
def S2_dens(f, e1, e2):
    A, ls, cs = f['A'], f['ls'], f['cs']
    Gam = np.array([(RHO[M] @ e1) @ Bf @ e2 for M in range(10)])
    GL = np.array([np.einsum('sa,st,ta->', np.einsum('st,ta->sa', RHO[M], ls), Bf, ls) for M in range(10)])
    U1 = 0.25 * Gam @ GL
    p1 = np.einsum('s,st,ta->a', e1, Bf, ls); p2 = np.einsum('s,st,ta->a', e2, Bf, ls)
    U2 = -0.5 * p1 @ p2
    U3 = -np.einsum('M,Ma,a->', Gam, A, cs)
    return G(U1 + U2 + U3)

# ------------------------------------------------------------------ random field jets at a point
def rfield():
    cr = lambda *sh: rng.normal(size=sh) + 1j * rng.normal(size=sh)
    lam = np.zeros((32, 3, NG), complex); lam[:16] = cr(16, 3, NG)
    dlam = np.zeros((4, 32, 3, NG), complex); dlam[:, :16] = cr(4, 16, 3, NG)
    ls = np.zeros((32, 3), complex); ls[16:] = cr(16, 3)
    return {'A': cr(10, 3), 'dA': cr(4, 10, 3), 'lam': lam, 'dlam': dlam, 'c': cr(3, NG), 'dc': cr(4, 3, NG),
            'As': cr(10, 3, NG), 'ls': ls, 'cs': cr(3)}

# ------------------------------------------------------------------ Theta_n on jets: fields at y -> fields at x0 = r_n y
def theta_n(f, n, slam=1.0):
    rn = np.eye(4) - 2 * np.outer(n, n); Rn = np.eye(10); Rn[:4, :4] = rn
    rnn = rho_V(n)
    g = {}
    g['A'] = Rn @ np.conj(f['A'])
    g['dA'] = np.einsum('kj,MN,kNa->jMa', rn, Rn, np.conj(f['dA']))
    g['lam'] = slam * np.einsum('st,taj->saj', rnn, cc(f['lam']))
    g['dlam'] = slam * np.einsum('kj,st,ktag->jsag', rn, rnn, np.array([cc(f['dlam'][k]) for k in range(4)]))
    g['c'] = np.conj(f['c']); g['dc'] = np.einsum('kj,kag->jag', rn, np.conj(f['dc']))
    g['As'] = -np.einsum('MN,Naj->Maj', Rn, np.conj(f['As']))
    g['ls'] = -slam * rnn @ cc(f['ls'])
    g['cs'] = -np.conj(f['cs'])
    return g
# conformal Killing field data
def vf_K(b, x):  # K_b: v = x^2 b - 2 (b.x) x, Jacobian
    v = (x @ x) * b - 2 * (b @ x) * x
    J = 2 * np.outer(b, x) - 2 * np.outer(x, b) - 2 * (b @ x) * np.eye(4)
    return v, J
def even_data(par, x):
    """par = (a, Om, d, b, Rw): v = a + Om x + d 2x + K_b; Om antisymmetric."""
    a, Om, d, b, Rw = par
    vK, JK = vf_K(b, x)
    v = a + Om @ x + 2 * d * x + vK
    J = Om + 2 * d * np.eye(4) + JK
    return v, J, Rw
def even_dagn(par, n):
    """v^{dag_n} = -(r_n)_* conj(v), R^{dag_n} = -conj(R): returns parameters of the image."""
    a, Om, d, b, Rw = par
    rn = np.eye(4) - 2 * np.outer(n, n)
    # (r_n)_* of a + Om x + 2 d x + K_b : r_n a + (r_n Om r_n) x + 2 d x + K_{r_n b}
    return (-(rn @ np.conj(a)), -(rn @ np.conj(Om) @ rn), -np.conj(d), -(rn @ np.conj(b)), -np.conj(Rw))
def even_dag(par):
    """v^dag = -I_* conj(v): P_a -> -K_conj(a), K_b -> -P_conj(b), Delta -> Delta, rotations -> -conj."""
    a, Om, d, b, Rw = par
    return (-np.conj(b), -np.conj(Om), np.conj(d), -np.conj(a), -np.conj(Rw))
def rand_even():
    cr = lambda *sh: rng.normal(size=sh) + 1j * rng.normal(size=sh)
    Om = cr(4, 4); Om = Om - Om.T
    Rw = cr(6, 6); Rw = Rw - Rw.T
    return (cr(4), Om, complex(rng.normal() + 1j * rng.normal()), cr(4), Rw)
def rand_odd():
    Q1 = np.zeros(32, complex); Q1[:16] = rng.normal(size=16) + 1j * rng.normal(size=16)
    Q2 = np.zeros(32, complex); Q2[16:] = rng.normal(size=16) + 1j * rng.normal(size=16)
    return Q1, Q2
def odd_dagn(Q, n, sgn=1.0):
    Q1, Q2 = Q; rnn = rho_V(n)
    return (sgn * rnn @ cc(Q1), -sgn * rnn @ cc(Q2))
def odd_dag(Q):
    Q1, Q2 = Q
    return (-cc(Q2), -cc(Q1))
def eps_at(Q, x):
    return Q[0] + rhox(x) @ Q[1]

OUT = {}
# ------------------------------------------------------------------ 1. Theta_n
def check_theta_n(slam=1.0, odd_sign=1.0, trials=3):
    res = {'S_terms': {}, 'S1_even': 0.0, 'S1_odd': 0.0, 'S2': 0.0}
    for _ in range(trials):
        n = rng.normal(size=4); n /= np.linalg.norm(n)
        rn = np.eye(4) - 2 * np.outer(n, n)
        y = rng.normal(size=4); x0 = rn @ y
        f = rfield(); g = theta_n(f, n, slam)
        tg = S_terms(g); tf = S_terms(f)
        for k in tf:
            res['S_terms'][k] = max(res['S_terms'].get(k, 0.0), (tg[k].c() - tf[k]).norm())
        par = rand_even()
        v, J, Rw = even_data(par, x0)
        lhs = S1_even(g, v, J, Rw).c()
        pd = even_dagn(par, n); v2, J2, Rw2 = even_data(pd, y)
        rhs = S1_even(f, v2, J2, Rw2)
        res['S1_even'] = max(res['S1_even'], (lhs - rhs).norm())
        Q = rand_odd()
        lhs = S1_odd(g, eps_at(Q, x0), Q[1]).c()
        Qd = odd_dagn(Q, n, odd_sign)
        rhs = S1_odd(f, eps_at(Qd, y), Qd[1]).scale(-1.0)          # (-1)^{|X|} = -1
        res['S1_odd'] = max(res['S1_odd'], (lhs - rhs).norm())
        Qa, Qb = rand_odd(), rand_odd()
        lhs = S2_dens(g, eps_at(Qa, x0), eps_at(Qb, x0)).c()
        Qad, Qbd = odd_dagn(Qa, n, odd_sign), odd_dagn(Qb, n, odd_sign)
        rhs = S2_dens(f, eps_at(Qbd, y), eps_at(Qad, y))
        res['S2'] = max(res['S2'], (lhs - rhs).norm())
    return res
OUT['Theta_n'] = check_theta_n()
OUT['Theta_n_control_slambda_minus'] = check_theta_n(slam=-1.0)
OUT['Theta_n_control_opposite_odd_sign_in_dagger_n'] = check_theta_n(odd_sign=-1.0)

# ------------------------------------------------------------------ 2. Theta_I (radial), jets by exact chain rule
def theta_I(f, x0, slam=-1.0):
    """fields at y = I(x0) given by f (value + gradient); returns the jets of Theta_I e at x0 (Def Theta of
    conjugation.tex = Prop c8-C).  Derivatives by the product rule; d/dx_j of e(I(x)) = sum_k (dI)_{kj} (d_k e)."""
    xx = x0 @ x0; xh = x0 / np.sqrt(xx)
    dI = (np.eye(4) - 2 * np.outer(x0, x0) / xx) / xx                  # dI_x (symmetric), (dI)_{kj} = d I^k / d x^j
    def d_dI(j):   # derivative of the matrix dI_x along x_j
        e = np.eye(4)[j]
        return (-2 * (np.outer(e, x0) + np.outer(x0, e)) / xx + 4 * x0[j] * np.outer(x0, x0) / xx ** 2) / xx - 2 * x0[j] * (np.eye(4) - 2 * np.outer(x0, x0) / xx) / xx ** 2
    rx = rhox(x0)
    g = {}
    # A|_V: pullback of conj(A|_V) along I: (I^* a)_m(x) = sum_k (dI_x)_{km} a_k(I x)
    AV = np.conj(f['A'][:4]); dAV = np.conj(f['dA'][:, :4])            # dAV[k, n, a] = d_k a_n at I(x0)
    AW = np.conj(f['A'][4:]); dAW = np.conj(f['dA'][:, 4:])
    A = np.zeros((10, 3), complex); dA = np.zeros((4, 10, 3), complex)
    A[:4] = np.einsum('km,ka->ma', dI, AV)
    A[4:] = AW / xx
    for j in range(4):
        # d/dx_j [ (dI)_{km} a_k(I x) ] = (d_j dI)_{km} a_k + (dI)_{km} sum_l (dI)_{lj} (d_l a_k)
        A_j_V = np.einsum('km,ka->ma', d_dI(j), AV) + np.einsum('km,l,lka->ma', dI, dI[:, j], dAV)
        dA[j, :4] = A_j_V
        dA[j, 4:] = (-2 * x0[j] / xx ** 2) * AW + np.einsum('l,lIa->Ia', dI[:, j], dAW) / xx
    g['A'] = A; g['dA'] = dA
    # lambda: slam (x,x)^{-2} rho(x) cc(lambda(Ix))
    L = cc(f['lam']); dL = np.array([cc(f['dlam'][k]) for k in range(4)])
    g['lam'] = slam * np.einsum('st,taj->saj', rx, L) / xx ** 2
    dl = np.zeros((4, 32, 3, NG), complex)
    for j in range(4):
        dl[j] = slam * (np.einsum('st,taj->saj', RHO[j], L) / xx ** 2 - 4 * x0[j] / xx ** 3 * np.einsum('st,taj->saj', rx, L)
                        + np.einsum('st,taj->saj', rx, np.einsum('l,ltaj->taj', dI[:, j], dL)) / xx ** 2)
    g['dlam'] = dl
    g['c'] = np.conj(f['c']); g['dc'] = np.einsum('lj,lag->jag', dI, np.conj(f['dc']))
    Rx = np.eye(10); Rx[:4, :4] = np.eye(4) - 2 * np.outer(xh, xh)
    g['As'] = -np.einsum('MN,Naj->Maj', Rx, np.conj(f['As'])) / xx ** 3
    g['ls'] = -slam * (rx @ cc(f['ls'])) / xx ** 3        # slam = -1 gives +rho(x)/(x,x)^3 as in Prop c8-C
    g['cs'] = -np.conj(f['cs']) / xx ** 4
    return g
def check_theta_I(slam=-1.0, odd_sign=1.0, trials=3):
    res = {'S_terms': {}, 'S1_even': 0.0, 'S1_odd': 0.0, 'S2': 0.0}
    for _ in range(trials):
        y = rng.normal(size=4); x0 = y / (y @ y)
        jac = (x0 @ x0) ** 4          # dx = (y,y)^{-4} dy, i.e. density at x0 times (x0,x0)^4 ... see below
        # int dx f(x) = int dy |det dI_y| f(I y) and |det dI_y| = (y,y)^{-4} = (x0,x0)^4
        f = rfield(); g = theta_I(f, x0, slam)
        tg = S_terms(g); tf = S_terms(f)
        for k in tf:
            res['S_terms'][k] = max(res['S_terms'].get(k, 0.0), (tg[k].c().scale(jac) - tf[k]).norm())
        par = rand_even()
        v, J, Rw = even_data(par, x0)
        lhs = S1_even(g, v, J, Rw).c().scale(jac)
        pd = even_dag(par); v2, J2, Rw2 = even_data(pd, y)
        rhs = S1_even(f, v2, J2, Rw2)
        res['S1_even'] = max(res['S1_even'], (lhs - rhs).norm())
        Q = rand_odd()
        lhs = S1_odd(g, eps_at(Q, x0), Q[1]).c().scale(jac)
        Qd = odd_dag(Q); Qd = (odd_sign * Qd[0], odd_sign * Qd[1])
        rhs = S1_odd(f, eps_at(Qd, y), Qd[1]).scale(-1.0)
        res['S1_odd'] = max(res['S1_odd'], (lhs - rhs).norm())
        Qa, Qb = rand_odd(), rand_odd()
        lhs = S2_dens(g, eps_at(Qa, x0), eps_at(Qb, x0)).c().scale(jac)
        Qad, Qbd = odd_dag(Qa), odd_dag(Qb)
        Qad = (odd_sign * Qad[0], odd_sign * Qad[1]); Qbd = (odd_sign * Qbd[0], odd_sign * Qbd[1])
        rhs = S2_dens(f, eps_at(Qbd, y), eps_at(Qad, y))
        res['S2'] = max(res['S2'], (lhs - rhs).norm())
    return res
OUT['Theta_I'] = check_theta_I()
# total-derivative form of the scalar kinetic deviation under Theta_I (conjugation.tex Thm fields, ch08 Prop c8-C):
# [s(Theta_I e)(x)]^c = (x,x)^{-4} s(e)(Ix) + TD(x),  TD = -1/4 sum_I d_mu( Om d_mu Om (phi_I, phi_I) ),
# Om = (x,x)^{-1}, phi_I(x) = A_I(I x) (after ^c), Laplacian(Om) = 0 on R^4 minus 0.
td_dev = 0.0
for _ in range(4):
    y = rng.normal(size=4); x0 = y / (y @ y); xx = x0 @ x0
    f = rfield(); g = theta_I(f, x0)
    dev = (S_terms(g)['YM_scalar_kin'].c().scale(xx ** 4) - S_terms(f)['YM_scalar_kin']).s / xx ** 4
    dI = (np.eye(4) - 2 * np.outer(x0, x0) / xx) / xx
    phi = f['A'][4:]                                   # A_I^a at I(x0)
    dphi = np.einsum('lm,lIa->mIa', dI, f['dA'][:, 4:])  # d_mu phi_I^a
    Om = 1 / xx; dOm = -2 * x0 / xx ** 2
    TD = -0.25 * ((dOm @ dOm) * np.sum(phi * phi) + 2 * Om * np.einsum('m,Ia,mIa->', dOm, phi, dphi))
    td_dev = max(td_dev, abs(dev - TD))
OUT['Theta_I_scalar_kinetic_deviation_minus_stated_total_derivative'] = td_dev
OUT['Theta_I_control_slambda_plus'] = check_theta_I(slam=1.0)
OUT['Theta_I_control_opposite_odd_sign'] = check_theta_I(odd_sign=-1.0)
OUT['note'] = ('values are max deviations (0 up to rounding = identity holds pointwise); S_terms are the '
               'separate terms of the action; for Theta_I the scalar kinetic term YM_scalar_kin (inside YM) is '
               'expected to change by a total derivative only')
OUT['peak_rss_MB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6
with open(os.path.join(HERE, 'a2_results.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=float)
print(json.dumps(OUT, indent=1, default=float))
