"""
model.py -- explicit components of the reduced N=4 gauge-BV theory and of the superconformal
Noether functionals S^(1).  Explicit bases live only here (and in verify.py).

Clifford data, Lie algebra and field layout are imported (read-only) from classical/n4_model.py
[classical_spec Secs 1-4]:
  V_10 = R^10, e_M, M = 0..9;  V = {0,1,2,3} (R^4);  W_R^perp = {4..9};
  rho(e_M): Sigma -> Sigma^*  = sig[M];   rho(e_M): Sigma^* -> Sigma = sigbar[M];
  (s, s*) = s^T s*;  Gamma(Q1,Q2)^M = (rho(e_M)Q1, Q2) = Q2^T sig[M] Q1   [ESW (28)];
  g = su(2), (T_i,T_j) = delta_ij (kappa = 1), [T_i,T_j] = f_ijk T_k.
"""
import os, sys
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
CLASSICAL = os.path.join(os.path.dirname(HERE), 'classical')
sys.path.insert(0, CLASSICAL)
import n4_model as nm          # noqa: E402  (read-only import of the classical realization)
sys.path.pop(0)
import engine as E              # noqa: E402

V4 = (0, 1, 2, 3)
WP = (4, 5, 6, 7, 8, 9)
NS = 16

SIG = [[[(Fr(int(nm.SIG_RE[M][a, b])), Fr(int(nm.SIG_IM[M][a, b]))) for b in range(NS)] for a in range(NS)] for M in range(10)]
SIGB = [[[(Fr(int(nm.SIGBAR_RE[M][a, b])), Fr(int(nm.SIGBAR_IM[M][a, b]))) for b in range(NS)] for a in range(NS)] for M in range(10)]
SIG_NZ = [[(a, b, SIG[M][a][b]) for a in range(NS) for b in range(NS) if not E.ciszero(SIG[M][a][b])] for M in range(10)]
SIGB_NZ = [[(a, b, SIGB[M][a][b]) for a in range(NS) for b in range(NS) if not E.ciszero(SIGB[M][a][b])] for M in range(10)]


def matvec(Mt, v):
    """Mt: 16x16 Gaussian-Fraction matrix, v: list of 16 Gaussian Fractions."""
    out = []
    for a in range(NS):
        s = E.ZEROC
        row = Mt[a]
        for b in range(NS):
            if not E.ciszero(row[b]) and not E.ciszero(v[b]):
                s = E.cadd(s, E.cmul(row[b], v[b]))
        out.append(s)
    return out


def pair(u, v):
    """(u, v) : Sigma (x) Sigma^* -> C, u^T v."""
    s = E.ZEROC
    for a in range(NS):
        s = E.cadd(s, E.cmul(u[a], v[a]))
    return s


def rho_vec_on_Sigma(vec, s):
    """rho(v): Sigma -> Sigma^*, v = sum_M vec[M] e_M."""
    out = [E.ZEROC] * NS
    for M in range(10):
        if E.ciszero(vec[M]):
            continue
        t = matvec(SIG[M], s)
        out = [E.cadd(out[a], E.cmul(vec[M], t[a])) for a in range(NS)]
    return out


def rho_vec_on_Sigmastar(vec, s):
    out = [E.ZEROC] * NS
    for M in range(10):
        if E.ciszero(vec[M]):
            continue
        t = matvec(SIGB[M], s)
        out = [E.cadd(out[a], E.cmul(vec[M], t[a])) for a in range(NS)]
    return out


def Gamma_num(Q1, Q2, star=False):
    """Gamma(Q1,Q2)^M = (rho(e_M) Q1, Q2) [ESW (28)]; star=True for Q1,Q2 in Sigma^*."""
    out = []
    for M in range(10):
        t = matvec(SIGB[M] if star else SIG[M], Q1)
        out.append(pair(t, Q2))
    return out


# ------------------------------------------------------------------------------------------
# Lie algebra
# ------------------------------------------------------------------------------------------
DIMG, FINT = nm.f_int('su2')
FL = [((j, k, l), Fr(v)) for (j, k, l), v in FINT.items()]


# ------------------------------------------------------------------------------------------
# field helpers (polynomials)
# ------------------------------------------------------------------------------------------
def e_mu(mu):
    J = [0, 0, 0, 0]
    J[mu] = 1
    return tuple(J)


class Fields:
    def __init__(self, V, dimg=DIMG):
        self.V = V
        self.dimg = dimg

    def A(self, M, i, J=E.ZERO4):
        return E.var(self.V, 'A', M, i, J)

    def lam(self, a, i, J=E.ZERO4):
        return E.var(self.V, 'lam', a, i, J)

    def c(self, i, J=E.ZERO4):
        return E.var(self.V, 'c', 0, i, J)

    def As(self, M, i):
        return E.var(self.V, 'As', M, i)

    def lams(self, a, i):
        return E.var(self.V, 'lams', a, i)

    def cs(self, i):
        return E.var(self.V, 'cs', 0, i)

    def x(self, mu):
        return E.var(self.V, 'x', mu, 0)

    def d(self, P, M):
        """derivative along e_M (zero for M in W_perp)."""
        if M in V4:
            return E.Dtot(self.V, P, M)
        return {}

    def bracket_g(self, X, Y, i):
        """[X,Y]^i = f_jki X^j Y^k ; X, Y: functions j -> poly."""
        out = {}
        for (j, k, l), val in FL:
            if l == i:
                E.padd_into(out, E.pmul(X(j), Y(k)), (val, Fr(0)))
        return out

    def F(self, M, N, i):
        out = {}
        if M in V4:
            E.padd_into(out, self.A(N, i, e_mu(M)))
        if N in V4:
            E.padd_into(out, self.A(M, i, e_mu(N)), (Fr(-1), Fr(0)))
        E.padd_into(out, self.bracket_g(lambda j: self.A(M, j), lambda k: self.A(N, k), i))
        return out

    def Dlam(self, M, a, i):
        out = {}
        if M in V4:
            E.padd_into(out, self.lam(a, i, e_mu(M)))
        E.padd_into(out, self.bracket_g(lambda j: self.A(M, j), lambda k: self.lam(a, k), i))
        return out


# ------------------------------------------------------------------------------------------
# the action [ESW (31)] reduced [classical_spec Sec 4]
# ------------------------------------------------------------------------------------------
def build_S(V, dimg=DIMG):
    Fd = Fields(V, dimg)
    S = {}
    # -1/4 (F,F)_Hodge = -1/8 sum_{M != N} F_MN F_MN
    for M in range(10):
        for N in range(10):
            if M == N:
                continue
            for i in range(dimg):
                f = Fd.F(M, N, i)
                E.padd_into(S, E.pmul(f, f), (Fr(-1, 8), Fr(0)))
    # +1/2 (lam, dslash_A lam) = 1/2 sum_M lam^a sig^M_ab (D_M lam)^b
    for M in range(10):
        for (a, b, s) in SIG_NZ[M]:
            for i in range(dimg):
                E.padd_into(S, E.pmul(Fd.lam(a, i), Fd.Dlam(M, b, i)), E.cmul(s, (Fr(1, 2), Fr(0))))
    # (d_A c, A^*)
    for M in range(10):
        for i in range(dimg):
            Dc = {}
            if M in V4:
                E.padd_into(Dc, Fd.c(i, e_mu(M)))
            E.padd_into(Dc, Fd.bracket_g(lambda j: Fd.A(M, j), lambda k: Fd.c(k), i))
            E.padd_into(S, E.pmul(Dc, Fd.As(M, i)))
    # ([c,lam], lam^*)
    for a in range(NS):
        for i in range(dimg):
            cl = Fd.bracket_g(lambda j: Fd.c(j), lambda k: Fd.lam(a, k), i)
            E.padd_into(S, E.pmul(cl, Fd.lams(a, i)))
    # 1/2 ([c,c], c^*)
    for i in range(dimg):
        cc = Fd.bracket_g(lambda j: Fd.c(j), lambda k: Fd.c(k), i)
        E.padd_into(S, E.pmul(cc, Fd.cs(i)), (Fr(1, 2), Fr(0)))
    return S


# ------------------------------------------------------------------------------------------
# Killing spinors and odd generators
# ------------------------------------------------------------------------------------------
def killing_spinor(V, Q, eta):
    """eps(x) = Q + rho(x) eta, as 16 polynomials in x; Q in Sigma, eta in Sigma^*."""
    Fd = Fields(V)
    out = []
    rows = [matvec(SIGB[mu], eta) for mu in V4]
    for a in range(NS):
        p = E.const(Q[a])
        for mu in V4:
            if not E.ciszero(rows[mu][a]):
                E.padd_into(p, Fd.x(mu), rows[mu][a])
        out.append(p)
    return out


def S1_odd(V, Q, eta, aS=Fr(0), dimg=DIMG, parts=False):
    """S^(1)(Q + rho(x) eta) = int dvol ( -(Gamma(eps,lam),A^*) + 1/2 (rho(F_A) eps, lam^*)
                                           + aS (rho(A|_Wperp) eta, lam^*) )."""
    Fd = Fields(V, dimg)
    eps = killing_spinor(V, Q, eta)
    T1, T2, T3 = {}, {}, {}
    # -(Gamma(eps,lam), A^*) = - sum_M sum_ab lam^a sig^M_ab eps^b A^*_M
    for M in range(10):
        for (a, b, s) in SIG_NZ[M]:
            if not eps[b]:
                continue
            for i in range(dimg):
                t = E.pmul(E.pmul(eps[b], Fd.lam(a, i)), Fd.As(M, i))
                E.padd_into(T1, t, E.cneg(s))
    # 1/2 (rho(F) eps, lam^*) = 1/2 sum_{M<N} F_MN (sigbar^M sig^N eps)^a lam^*_a
    for M in range(10):
        for N in range(M + 1, 10):
            # vec^a = sum_{b,c} sigbar[M]_{ab} sig[N]_{bc} eps^c   (polys in x)
            vec = [dict() for _ in range(NS)]
            for (b, cc, s2) in SIG_NZ[N]:
                if not eps[cc]:
                    continue
                for a in range(NS):
                    s1 = SIGB[M][a][b]
                    if E.ciszero(s1):
                        continue
                    E.padd_into(vec[a], eps[cc], E.cmul(s1, s2))
            for i in range(dimg):
                f = Fd.F(M, N, i)
                if not f:
                    continue
                for a in range(NS):
                    if vec[a]:
                        E.padd_into(T2, E.pmul(E.pmul(vec[a], f), Fd.lams(a, i)), (Fr(1, 2), Fr(0)))
    # aS (rho(A|_Wperp) eta, lam^*) = aS sum_w A_w (sigbar^w eta)^a lam^*_a
    if any(not E.ciszero(e) for e in eta):
        for w in WP:
            t = matvec(SIGB[w], eta)
            for a in range(NS):
                if E.ciszero(t[a]):
                    continue
                for i in range(dimg):
                    E.padd_into(T3, E.pmul(Fd.A(w, i), Fd.lams(a, i)), t[a])
    if parts:
        return T1, T2, T3
    out = E.padd(T1, T2)
    E.padd_into(out, T3, (aS, Fr(0)) if not isinstance(aS, tuple) else aS)
    return out


# ------------------------------------------------------------------------------------------
# even generators: (vector field on V with polynomial components, R in so(W_perp))
# ------------------------------------------------------------------------------------------
def xpoly(V, coeffs):
    """coeffs: dict {exponent tuple (4): Gaussian Fraction} -> polynomial in x."""
    Fd = Fields(V)
    P = {}
    for ex, c in coeffs.items():
        t = E.const(c)
        for mu in V4:
            for _ in range(ex[mu]):
                t = E.pmul(t, Fd.x(mu))
        E.padd_into(P, t)
    return P


def even_data(V, vf, R):
    """vf: list of 4 x-polys (components v^mu), R: 10x10 dict {(i,j): c} with v^i = R^i_j y^j on W_perp.
    returns dict with v, dv[mu][nu] = d_mu v^nu, div, Omega[(M,N)] = d_M v_N - d_N v_M."""
    dv = [[E.Dtot(V, vf[nu], mu) for nu in V4] for mu in V4]
    div = {}
    for mu in V4:
        E.padd_into(div, dv[mu][mu])
    Om = {}
    for M in V4:
        for N in V4:
            if M < N:
                o = E.psub(dv[M][N], dv[N][M])
                if o:
                    Om[(M, N)] = o
    for (i, j), cval in R.items():
        # d_j v^i = R^i_j  ->  Omega_{ji} += R^i_j, Omega_{ij} -= R^i_j
        if j < i:
            Om.setdefault((j, i), {})
            E.padd_into(Om[(j, i)], E.const(cval))
        elif i < j:
            Om.setdefault((i, j), {})
            E.padd_into(Om[(i, j)], E.const(cval), (Fr(-1), Fr(0)))
    Om = {k: v for k, v in Om.items() if v}
    return {'v': vf, 'dv': dv, 'div': div, 'Om': Om, 'R': R}


def S1_even(V, vf, R, kl, wl, wA=Fr(1, 4), dimg=DIMG, parts=False):
    """S^(1)(v) = int dvol ( (L_v A + wA (div v) A|_Wperp, A^*) - (v.lam, lam^*) - (v(c)) c^* ),
       v.lam = v(lam) + kl rho(d(v,-)) lam + wl (div v) lam."""
    Fd = Fields(V, dimg)
    D = even_data(V, vf, R)
    TA, TL, TC = {}, {}, {}
    for i in range(dimg):
        # delta A_M
        for M in range(10):
            dA = {}
            for nu in V4:
                if vf[nu]:
                    E.padd_into(dA, E.pmul(vf[nu], Fd.A(M, i, e_mu(nu))))
            if M in V4:
                for nu in V4:
                    if D['dv'][M][nu]:
                        E.padd_into(dA, E.pmul(Fd.A(nu, i), D['dv'][M][nu]))
            else:
                for (a_, b_), cval in R.items():
                    # (L_v A)_w = A_j R^j_w
                    if b_ == M:
                        E.padd_into(dA, Fd.A(a_, i), cval)
                if D['div']:
                    E.padd_into(dA, E.pmul(D['div'], Fd.A(M, i)), (wA, Fr(0)))
            if dA:
                E.padd_into(TA, E.pmul(dA, Fd.As(M, i)))
        # v.lam
        for a in range(NS):
            dl = {}
            for nu in V4:
                if vf[nu]:
                    E.padd_into(dl, E.pmul(vf[nu], Fd.lam(a, i, e_mu(nu))))
            for (M, N), om in D['Om'].items():
                # (sigbar^M sig^N lam)^a
                for (b, cc, s2) in SIG_NZ[N]:
                    s1 = SIGB[M][a][b]
                    if E.ciszero(s1):
                        continue
                    E.padd_into(dl, E.pmul(om, Fd.lam(cc, i)), E.cmul((kl, Fr(0)), E.cmul(s1, s2)))
            if D['div']:
                E.padd_into(dl, E.pmul(D['div'], Fd.lam(a, i)), (wl, Fr(0)))
            if dl:
                E.padd_into(TL, E.pmul(dl, Fd.lams(a, i)), (Fr(-1), Fr(0)))
        dc = {}
        for nu in V4:
            if vf[nu]:
                E.padd_into(dc, E.pmul(vf[nu], Fd.c(i, e_mu(nu))))
        if dc:
            E.padd_into(TC, E.pmul(dc, Fd.cs(i)), (Fr(-1), Fr(0)))
    if parts:
        return TA, TL, TC
    return E.padd(TA, TL, TC)


# ------------------------------------------------------------------------------------------
# quadratic-antifield functional [ESW (34)] evaluated pointwise on Killing spinors
# ------------------------------------------------------------------------------------------
def Gamma_poly(e1, e2, M):
    """Gamma(e1,e2)^M = sum_ab e2^a sig^M_ab e1^b  (e1, e2: lists of 16 x-polys in Sigma)."""
    out = {}
    for (a, b, s) in SIG_NZ[M]:
        if e1[b] and e2[a]:
            E.padd_into(out, E.pmul(e2[a], e1[b]), s)
    return out


def S2_odd(V, e1, e2, dimg=DIMG, wperp=True, parts=False):
    """S^(2)(eps1,eps2) = int dvol ( 1/4 (Gamma(eps1,eps2), Gamma(lam^*,lam^*)) - 1/2 (eps1,lam^*)(eps2,lam^*)
                                      - iota_{Gamma(eps1,eps2)} A c^* )   [ESW (34), pointwise in x].
    wperp=False drops the W_perp components of Gamma(eps1,eps2) in the first and last terms."""
    Fd = Fields(V, dimg)
    Ms = range(10) if wperp else V4
    G = {M: Gamma_poly(e1, e2, M) for M in range(10)}
    T1, T2, T3 = {}, {}, {}
    for M in Ms:
        if not G[M]:
            continue
        for i in range(dimg):
            gl = {}
            for (a, b, s) in SIGB_NZ[M]:
                E.padd_into(gl, E.pmul(Fd.lams(a, i), Fd.lams(b, i)), s)
            if gl:
                E.padd_into(T1, E.pmul(G[M], gl), (Fr(1, 4), Fr(0)))
            E.padd_into(T3, E.pmul(E.pmul(G[M], Fd.A(M, i)), Fd.cs(i)), (Fr(-1), Fr(0)))
    for i in range(dimg):
        p1, p2 = {}, {}
        for a in range(NS):
            if e1[a]:
                E.padd_into(p1, E.pmul(e1[a], Fd.lams(a, i)))
            if e2[a]:
                E.padd_into(p2, E.pmul(e2[a], Fd.lams(a, i)))
        E.padd_into(T2, E.pmul(p1, p2), (Fr(-1, 2), Fr(0)))
    if parts:
        return T1, T2, T3
    return E.padd(T1, T2, T3)
