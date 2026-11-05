"""
models.py -- local functionals of the abelian theory (one colour component, pairing 1) in the conventions of
Chapters 1-2, with the reduced non-minimal pair (cbar, cbar^*), in three frames:
  'min' : minimal fields only (former Chapter 3),
  'C'   : Costello frame, S' = S - 1/4 int (cbar^*, cbar^*), equivariant extension on (cbar, cbar^*),
  'Psi' : gauge-fermion frame, every functional composed with Phi: A^* -> A^* - d cbar, cbar^* -> cbar^* + d^*A|_V.
Explicit components live only here.
"""
import sympy as sp
import bvjet
import realization as R

X = bvjet.X
Rm_blocks = R.rho_blocks()
Cp = R.Cpair()


def build_model(with_nm=True):
    M = bvjet.Model()
    M.add('c', 1)
    for m in range(4):
        M.add('A%d' % m, 0)
    for I in range(6):
        M.add('P%d' % I, 0)
    for a in range(16):
        M.add('L%d' % a, 1)
    for m in range(4):
        M.add('As%d' % m, 1)
    for I in range(6):
        M.add('Ps%d' % I, 1)
    for a in range(16):
        M.add('Ls%d' % a, 0)
    M.add('cs', 0)
    if with_nm:
        M.add('cb', 1)
        M.add('cbs', 0)
    # K_0 = -delta x inverse fibre pairing on both orders
    for m in range(4):
        M.pair('A%d' % m, 'As%d' % m, -1)
    for I in range(6):
        M.pair('P%d' % I, 'Ps%d' % I, -1)
    for a in range(16):
        M.pair('L%d' % a, 'Ls%d' % a, -Cp[a, a])
    M.pair('c', 'cs', -1)
    if with_nm:
        M.pair('cb', 'cbs', -1)
    return M


E4 = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0), (0, 0, 0, 1)]
Z4 = (0, 0, 0, 0)


class Fields:
    """Field vectors (lists of superpolynomials) in a given frame."""

    def __init__(self, M, frame):
        self.M = M
        self.frame = frame
        j = M.jet
        self.c = j('c')
        self.A = [j('A%d' % m) for m in range(4)] + [j('P%d' % I) for I in range(6)]   # ten components
        self.L = [j('L%d' % a) for a in range(16)]
        self.Ls = [j('Ls%d' % a) for a in range(16)]
        self.cs = j('cs')
        As = [j('As%d' % m) for m in range(4)] + [j('Ps%d' % I) for I in range(6)]
        if frame in ('C', 'Psi'):
            self.cb = j('cb')
            cbs = j('cbs')
        if frame == 'Psi':
            # A^*|_V -> A^*|_V - d cbar ;  cbar^* -> cbar^* + d^* A|_V = cbar^* - sum_m d_m A_m
            for m in range(4):
                As[m] = M.plus(As[m], j('cb', E4[m], -1))
            for m in range(4):
                cbs = M.plus(cbs, j('A%d' % m, E4[m], -1))
        self.As = As
        if frame in ('C', 'Psi'):
            self.cbs = cbs

    # helpers
    def d(self, P, m):
        return self.M.Dtot(P, m)

    def dvec(self, vec, m):
        return [self.M.Dtot(P, m) for P in vec]


def lin(M, Mat, vec):
    out = []
    for i in range(Mat.shape[0]):
        acc = {}
        for k in range(Mat.shape[1]):
            if Mat[i, k] != 0:
                acc = M.plus(acc, M.scale(vec[k], Mat[i, k]))
        out.append(acc)
    return out


def dot(M, u, w):
    acc = {}
    for a, b in zip(u, w):
        if a and b:
            acc = M.plus(acc, M.mul(a, b))
    return acc


def scal(M, P, c):
    return M.scale(P, c)


def vplus(M, *vs):
    return [M.plus(*xs) for xs in zip(*vs)]


def rho2(Mi, Nj):
    """rho(e_M) rho(e_N) on Sigma (16x16): Sigma -> Sigma^* -> Sigma."""
    return Rm_blocks[Mi][1] * Rm_blocks[Nj][0]


def rho2_star(Mi, Nj):
    """rho(e_M) rho(e_N) on Sigma^*: Sigma^* -> Sigma -> Sigma^*."""
    return Rm_blocks[Mi][0] * Rm_blocks[Nj][1]


# ------------------------------------------------------------------- the free action
def S0(F):
    """1/2 <e, Qe> of (c1-free) in the frame of F, plus -1/4 (cbar^*, cbar^*) in frames C, Psi."""
    M = F.M
    out = {}
    # -1/4 (dA,dA): (dA,dA) = sum_{M<N} F_MN^2 with F_mn = d_m A_n - d_n A_m, F_mI = d_m A_I
    for m in range(4):
        for n in range(m + 1, 4):
            Fmn = M.plus(F.d(F.A[n], m), M.scale(F.d(F.A[m], n), -1))
            out = M.plus(out, M.scale(M.mul(Fmn, Fmn), -sp.Rational(1, 4)))
        for I in range(6):
            FmI = F.d(F.A[4 + I], m)
            out = M.plus(out, M.scale(M.mul(FmI, FmI), -sp.Rational(1, 4)))
    # 1/2 (lambda, dslash lambda) = 1/2 lambda^T C (sum_m rho_m d_m lambda) ; rho_m : Sigma -> Sigma^*
    for m in range(4):
        dl = lin(M, Rm_blocks[m][0], F.dvec(F.L, m))
        out = M.plus(out, M.scale(dot(M, F.L, lin(M, Cp, dl)), sp.Rational(1, 2)))
    # (dc, A^*) = sum_m d_m c A^*_m
    for m in range(4):
        out = M.plus(out, M.mul(F.d(F.c, m), F.As[m]))
    if F.frame in ('C', 'Psi'):
        out = M.plus(out, M.scale(M.mul(F.cbs, F.cbs), -sp.Rational(1, 4)))
    return M.simplify(out)


# ------------------------------------------------------------------- conformal Killing data
def div(v):
    return sp.expand(sum(sp.diff(v[m], X[m]) for m in range(4)))


def dv_form(v):
    """components (d(v,-))_{mn} = d_m v_n - d_n v_m on V."""
    return [[sp.expand(sp.diff(v[n], X[m]) - sp.diff(v[m], X[n])) for n in range(4)] for m in range(4)]


def rho_q2form_V(w):
    """rho(q(omega)) on Sigma for omega_{mn} on V: sum_{m<n} omega_mn rho_m rho_n."""
    out = sp.zeros(16, 16)
    for m in range(4):
        for n in range(m + 1, 4):
            if w[m][n] != 0:
                out += w[m][n] * rho2(m, n)
    return out


def rho_q2form_V_star(w):
    out = sp.zeros(16, 16)
    for m in range(4):
        for n in range(m + 1, 4):
            if w[m][n] != 0:
                out += w[m][n] * rho2_star(m, n)
    return out


def rho_q2form_W(r):
    """for an element r of so(W_R^perp): (1/4) rho(d(v,-)) with omega_{IJ} = -2 r_{IJ}, on Sigma."""
    out = sp.zeros(16, 16)
    for I in range(6):
        for J in range(I + 1, 6):
            if r[I, J] != 0:
                out += (-2 * r[I, J]) * rho2(4 + I, 4 + J)
    return out


def rho_q2form_W_star(r):
    out = sp.zeros(16, 16)
    for I in range(6):
        for J in range(I + 1, 6):
            if r[I, J] != 0:
                out += (-2 * r[I, J]) * rho2_star(4 + I, 4 + J)
    return out


def apply_vec_field(F, v, P):
    """v(P) = sum_m v_m D_m P."""
    M = F.M
    out = {}
    for m in range(4):
        if v[m] != 0:
            out = M.plus(out, M.mul(M.const(v[m]), F.d(P, m)))
    return out


def S1_even(F, v=None, r=None, cbar_weight=2, cbar_transforms=True):
    """S^(1) along a conformal Killing field v of V (components polynomial in x), or along r in so(W_R^perp):
    (v.A, A^*) - (v.lambda, lambda^*) - (v.c, c^*)  [+ -(v.cbar, cbar^*) in frames C, Psi]."""
    M = F.M
    out = {}
    if v is not None:
        dvv = div(v)
        # v.A on A|_V : (L_v A)_m = v_n d_n A_m + A_n d_m v_n ; on A_I : v(A_I) + 1/4 div v A_I
        vA = []
        for m in range(4):
            t = apply_vec_field(F, v, F.A[m])
            for n in range(4):
                c_ = sp.diff(v[n], X[m])
                if c_ != 0:
                    t = M.plus(t, M.mul(M.const(c_), F.A[n]))
            vA.append(t)
        for I in range(6):
            t = apply_vec_field(F, v, F.A[4 + I])
            if dvv != 0:
                t = M.plus(t, M.mul(M.const(sp.Rational(1, 4) * dvv), F.A[4 + I]))
            vA.append(t)
        # v.lambda = v(lambda) + 1/4 rho(d(v,-)) lambda + 3/8 div v lambda
        Sp = rho_q2form_V(dv_form(v))
        vL = []
        for a in range(16):
            t = apply_vec_field(F, v, F.L[a])
            for b in range(16):
                c_ = sp.Rational(1, 4) * Sp[a, b] + (sp.Rational(3, 8) * dvv if a == b else 0)
                c_ = sp.expand(c_)
                if c_ != 0:
                    t = M.plus(t, M.mul(M.const(c_), F.L[b]))
            vL.append(t)
        vc = apply_vec_field(F, v, F.c)
    else:
        vA = [{} for _ in range(4)]
        for I in range(6):
            t = {}
            for J in range(6):
                if r[J, I] != 0:
                    t = M.plus(t, M.scale(F.A[4 + J], r[J, I]))
            vA.append(t)
        Sp = sp.Rational(1, 4) * rho_q2form_W(r)
        vL = lin(M, Sp, F.L)
        vc = {}
    out = M.plus(out, dot(M, vA, F.As))
    out = M.plus(out, M.scale(dot(M, vL, lin(M, Cp, F.Ls)), -1))
    out = M.plus(out, M.scale(M.mul(vc, F.cs), -1))
    if F.frame in ('C', 'Psi') and v is not None and cbar_transforms:
        vcb = apply_vec_field(F, v, F.cb)
        if div(v) != 0:
            vcb = M.plus(vcb, M.mul(M.const(sp.Rational(cbar_weight, 4) * div(v)), F.cb))
        out = M.plus(out, M.scale(M.mul(vcb, F.cbs), -1))
    return M.simplify(out)


# ------------------------------------------------------------------- odd generators
def killing(xi, eta):
    """epsilon(x) = xi + rho(x) eta (16-vector of polynomials in x), for xi in Sigma, eta in Sigma^*."""
    eps = sp.Matrix(xi)
    for m in range(4):
        eps += X[m] * (Rm_blocks[m][1] * sp.Matrix(eta))
    return eps.applyfunc(sp.expand)


def Gamma_vec(e1, e2):
    """Gamma(e1,e2)^M = (rho(e_M) e1, e2), e1,e2 in Sigma (functions of x)."""
    out = []
    for Mi in range(10):
        out.append(sp.expand((e2.T * Cp * (Rm_blocks[Mi][0] * e1))[0, 0]))
    return out


def Gamma_lam(F, eps):
    """the ten components of Gamma(eps, lambda) (linear in lambda), eps a Killing spinor (polynomial in x)."""
    M = F.M
    out = []
    for Mi in range(10):
        # Gamma^M(eps,lam) = (rho_M eps, lam) = lam^T C (rho_M eps)
        coeffs = (Cp * (Rm_blocks[Mi][0] * eps)).applyfunc(sp.expand)
        t = {}
        for a in range(16):
            if coeffs[a] != 0:
                t = M.plus(t, M.mul(M.const(coeffs[a]), F.L[a]))
        out.append(t)
    return out


def rhoF_eps(F, eps):
    """rho(F_A) eps for the abelian F_A = dA (Sigma-valued, linear in A): sum_{M<N} F_MN rho_M rho_N eps."""
    M = F.M
    out = [{} for _ in range(16)]
    pairs = []
    for m in range(4):
        for n in range(m + 1, 4):
            Fmn = M.plus(F.d(F.A[n], m), M.scale(F.d(F.A[m], n), -1))
            pairs.append((m, n, Fmn))
        for I in range(6):
            pairs.append((m, 4 + I, F.d(F.A[4 + I], m)))
    for (a, b, Fab) in pairs:
        vecc = (rho2(a, b) * eps).applyfunc(sp.expand)
        for k in range(16):
            if vecc[k] != 0:
                out[k] = M.plus(out[k], M.mul(M.const(vecc[k]), Fab))
    return out


def S1_odd(F, xi, eta):
    """S^(1) along the Killing spinor xi + rho(x) eta:
    -(Gamma(eps,lambda), A^*) + 1/2 (rho(F_A) eps, lambda^*) - (rho(A|_W) eta, lambda^*)."""
    M = F.M
    eps = killing(xi, eta)
    out = {}
    G = Gamma_lam(F, eps)
    out = M.plus(out, M.scale(dot(M, G, F.As), -1))
    rf = rhoF_eps(F, eps)
    out = M.plus(out, M.scale(dot(M, rf, lin(M, Cp, F.Ls)), sp.Rational(1, 2)))
    # -(rho(A|_W) eta, lambda^*) : rho(A|_W) : Sigma^* -> Sigma, sum_I A_I rho_I eta
    etav = sp.Matrix(eta)
    t = [{} for _ in range(16)]
    for I in range(6):
        vecc = (Rm_blocks[4 + I][1] * etav).applyfunc(sp.expand)
        for k in range(16):
            if vecc[k] != 0:
                t[k] = M.plus(t[k], M.mul(M.const(vecc[k]), F.A[4 + I]))
    out = M.plus(out, M.scale(dot(M, t, lin(M, Cp, F.Ls)), -1))
    return M.simplify(out)


def S2_odd(F, xi1, eta1, xi2, eta2, gamma2=None):
    """S^(2) of (c2-S2) on two Killing spinors, plus gamma2 (cbar, u(cbar)) with u = Gamma(eps1,eps2)|_V
    in frames C, Psi."""
    M = F.M
    e1 = killing(xi1, eta1)
    e2 = killing(xi2, eta2)
    G = Gamma_vec(e1, e2)
    out = {}
    # 1/4 (Gamma(e1,e2), Gamma(lambda^*, lambda^*)) ; Gamma(l*,l*)^M = (rho_M l*, l*) = (rho_M l*)^T C l*
    for Mi in range(10):
        if G[Mi] == 0:
            continue
        Rs = Rm_blocks[Mi][1]
        rl = lin(M, Rs, F.Ls)                     # rho_M lambda^* in Sigma
        gl = dot(M, rl, lin(M, Cp, F.Ls))         # (rho_M l*)^T C l*
        out = M.plus(out, M.mul(M.const(sp.Rational(1, 4) * G[Mi]), gl))
    # -1/2 (e1, lambda^*)(e2, lambda^*)
    a1 = {}
    a2 = {}
    v1 = (Cp * e1).applyfunc(sp.expand)
    v2 = (Cp * e2).applyfunc(sp.expand)
    for k in range(16):
        if v1[k] != 0:
            a1 = M.plus(a1, M.mul(M.const(v1[k]), F.Ls[k]))
        if v2[k] != 0:
            a2 = M.plus(a2, M.mul(M.const(v2[k]), F.Ls[k]))
    out = M.plus(out, M.scale(M.mul(a1, a2), -sp.Rational(1, 2)))
    # -(iota_Gamma A, c^*)
    t = {}
    for Mi in range(10):
        if G[Mi] != 0:
            t = M.plus(t, M.mul(M.const(G[Mi]), F.A[Mi]))
    out = M.plus(out, M.scale(M.mul(t, F.cs), -1))
    if gamma2 is not None and F.frame in ('C', 'Psi'):
        u = G[:4]
        out = M.plus(out, M.scale(M.mul(F.cb, apply_vec_field(F, u, F.cb)), gamma2))
    return M.simplify(out)


# ------------------------------------------------------------------- psl(4|4) brackets (Definition c2-bracket)
def bracket_even_even(v, w):
    return [sp.expand(sum(v[n] * sp.diff(w[m], X[n]) - w[n] * sp.diff(v[m], X[n]) for n in range(4)))
            for m in range(4)]


def bracket_even_odd(v, xi, eta, r=None):
    """[v, xi + rho(x) eta] (c2-evenodd) as a pair (xi', eta'); r : element of so(W) instead of v."""
    eps = killing(xi, eta)
    if r is None:
        dv = div(v)
        out = sp.zeros(16, 1)
        for m in range(4):
            out += v[m] * eps.diff(X[m])
        out += sp.Rational(1, 4) * rho_q2form_V(dv_form(v)) * eps - sp.Rational(1, 8) * dv * eps
    else:
        out = sp.Rational(1, 4) * rho_q2form_W(r) * eps
    out = out.applyfunc(sp.expand)
    xi_n = out.subs({t: 0 for t in X})
    lin_part = (out - xi_n).applyfunc(sp.expand)
    # recover eta' from rho(x) eta' : d/dx_m of rho(x) eta' = rho_m eta'; rho_m rho_m = 1 so eta' = rho_m^{-1} d_m
    # use m = 0: rho_0 : Sigma^* -> Sigma, inverse rho_0 : Sigma -> Sigma^*
    d0 = lin_part.diff(X[0])
    eta_n = (Rm_blocks[0][0] * d0).applyfunc(sp.expand)
    chk = (killing(list(xi_n), list(eta_n)) - out).applyfunc(sp.expand)
    assert all(c == 0 for c in chk), 'bracket is not a conformal Killing spinor'
    return list(xi_n), list(eta_n)


def bracket_odd_odd(xi1, eta1, xi2, eta2):
    """[e1, e2] = (x -> Gamma(e1(x), e2(x))|_V) + rotation of W (c2-oddodd)."""
    e1 = killing(xi1, eta1)
    e2 = killing(xi2, eta2)
    G = Gamma_vec(e1, e2)
    v = G[:4]
    Q1, Q2, Q3, Q4 = sp.Matrix(xi1), sp.Matrix(eta1), sp.Matrix(xi2), sp.Matrix(eta2)
    rot = sp.zeros(6, 6)
    for I in range(6):
        for J in range(6):
            if I == J:
                continue
            q = rho2(4 + I, 4 + J)          # rho(q(w_I ^ w_J)) on Sigma
            rot[I, J] = sp.expand((Q4.T * Cp.T * (q * Q1))[0, 0] + (Q2.T * Cp.T * (q * Q3))[0, 0])
    return v, rot


# ------------------------------------------------------------------- operators on coordinates
def apply_momentum_matrix_row(M, Qmat, a, p):
    """the jet expression sum_b Q_ab(p -> -i d) phi^b (coordinates of the model have the freecx layout)."""
    out = {}
    for b in range(Qmat.shape[1]):
        q = sp.expand(Qmat[a, b])
        if q == 0:
            continue
        poly = sp.Poly(q, *p)
        for mon, cf in poly.terms():
            coeff = cf * (-sp.I) ** sum(mon)
            out = M.plus(out, M.jet(M.fields[b].name, mon, coeff))
    return out


def calibrate_Q(M, F, Qmat, p, f=None, sign_rule=False):
    """{S0, l_a} = l_a o Q for l_a = int f phi^a, every coordinate a (l_a o Q: the coordinates of l_a replaced
    by the components of Q applied to the coordinates).  sign_rule=True tests (-1)^{|a|} l_a o Q instead."""
    if f is None:
        f = X[0] * X[1] ** 2 + 3 * X[2] - X[3] ** 3 + X[0] * X[3] + 2
    S = S0(F)
    bad = []
    for a in range(Qmat.shape[0]):
        la = M.mul(M.const(f), M.jet(M.fields[a].name))
        lhs = M.bracket(S, la)
        rhs = M.scale(M.mul(M.const(f), apply_momentum_matrix_row(M, Qmat, a, p)), (-1) ** M.par(a) if sign_rule else 1)
        if not M.is_zero_functional(M.plus(lhs, M.scale(rhs, -1))):
            bad.append(M.fields[a].name)
    return bad
