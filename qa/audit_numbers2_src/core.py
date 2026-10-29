# Fresh independent audit core (round 2, numbers lens).  No project code is imported.
# Conventions taken only from classical/classical_spec.md (the definition of the theory):
#   S = int dvol( -1/4 (F_A,F_A)_Hodge + 1/2 (lambda, dslash_A lambda) + (d_A c, A^*) + ... )
#   (F,F)_Hodge = 1/2 sum_{M,N} (F_MN, F_MN);  F_MN = d_M A_N - d_N A_M + [A_M, A_N]  (d_M = 0 for M > 4)
#   propagators = Wick contraction values of P(eps,infty) (spec section 7), re-derived here as
#   covariances / pseudo-inverses of the free quadratic forms (see derive_checks()).
# Explicit bases live only in this code: g = su(2), [t_a,t_b] = eps_abc t_c, (t_a,t_b) = delta_ab.
# Then tr_g(ad_X ad_Y)/(X,Y) = -2, so the unit hbar tr_g(ad ad)/((X,Y) 16 pi^2) = -2/(16 pi^2) (hbar = 1).
import itertools
from fractions import Fraction
import numpy as np

# ---------------------------------------------------------------- Lie algebra su(2)
NC = 3
EPS3 = np.zeros((3, 3, 3))
for a, b, c in itertools.permutations(range(3)):
    EPS3[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])

# ---------------------------------------------------------------- 10d chiral gamma blocks
I2 = np.eye(2)
PX = np.array([[0., 1.], [1., 0.]])
PZ = np.array([[1., 0.], [0., -1.]])
PE = np.array([[0., 1.], [-1., 0.]])  # real antisymmetric, square -1


def _kron(*ms):
    out = np.array([[1.]])
    for m in ms:
        out = np.kron(out, m)
    return out


def _spin9_real_symmetric():
    cands = []
    for combo in itertools.product(range(4), repeat=4):
        if all(c == 0 for c in combo):
            continue
        if sum(1 for c in combo if c == 3) % 2:
            continue
        mats = [I2, PX, PZ, PE]
        cands.append(_kron(*[mats[c] for c in combo]))

    def anti(a, b):
        return np.allclose(a @ b + b @ a, 0)

    best = []

    def bt(chosen, start):
        nonlocal best
        if len(chosen) == 9:
            best = list(chosen)
            return True
        for i in range(start, len(cands)):
            if all(anti(cands[i], cands[j]) for j in chosen):
                chosen.append(i)
                if bt(chosen, i + 1):
                    return True
                chosen.pop()
        return False

    bt([], 0)
    return [cands[i] for i in best]


_G9 = _spin9_real_symmetric()
# sigmabar^M : Sigma -> Sigma^*, symmetric for the identity pairing Sigma x Sigma^* -> C
SB = [g.astype(complex) for g in _G9] + [-1j * np.eye(16)]
SU = [g.astype(complex) for g in _G9] + [1j * np.eye(16)]  # sigma^M : Sigma^* -> Sigma
ND = 10      # V_10 directions; M = 0..3 are R^4 (real generators), M = 4..9 are W_R^perp
NV = 4


def check_clifford():
    ok = True
    for M in range(10):
        ok &= np.allclose(SB[M], SB[M].T)
        for N in range(10):
            ok &= np.allclose(SU[M] @ SB[N] + SU[N] @ SB[M], 2 * (M == N) * np.eye(16))
            ok &= np.allclose(SB[M] @ SU[N] + SB[N] @ SU[M], 2 * (M == N) * np.eye(16))
    return bool(ok)


# ---------------------------------------------------------------- index maps
# boson fluctuation index i = M*3 + a (M in 0..9, a in 0..2): 30 components
NB = ND * NC
# fermion index (b, alpha) -> b*16 + alpha : 48 components
NF = NC * 16
# ghost sector X = (c^b [3], A^*_{mu}^b [12]) ; A^*_mu^b -> 3 + mu*3 + b
NG = NC + NV * NC


def bidx(M, a):
    return M * NC + a


# ---------------------------------------------------------------- bosonic vertex tensors
def T3_tensor(p1):
    """Trilinear tensor of T3(A1,A2,A3) = -1/2 sum_{MN} (i p1_M A1_N, [A2_M, A3_N]); slot-1 momentum p1."""
    T = np.zeros((NB, NB, NB), dtype=complex)
    for M in range(NV):
        if p1[M] == 0:
            continue
        for N in range(ND):
            for a in range(3):
                for b in range(3):
                    for c in range(3):
                        e = EPS3[a, b, c]
                        if e == 0:
                            continue
                        # (A1_N)^c [A2_M]^a [A3_N]^b eps_abc
                        T[bidx(N, c), bidx(M, a), bidx(N, b)] += -0.5 * 1j * p1[M] * e
    return T


def gamma3(moms):
    """Fully symmetrised cubic vertex tensor Gamma3[i1,i2,i3] for momenta moms (all incoming)."""
    G = np.zeros((NB, NB, NB), dtype=complex)
    for perm in itertools.permutations(range(3)):
        # slot s holds leg perm[s]; T3 uses the momentum of the leg in slot 0
        T = T3_tensor(moms[perm[0]])
        # T indices are (slot0, slot1, slot2) = legs (perm[0], perm[1], perm[2]); reorder to legs (0,1,2)
        inv = np.argsort(perm)
        G += np.transpose(T, axes=inv)
    return G


_T4 = None


def T4_tensor():
    global _T4
    if _T4 is not None:
        return _T4
    T = np.zeros((NB, NB, NB, NB))
    # -1/8 sum_{MN} ([A1_M, A2_N], [A3_M, A4_N]) = -1/8 eps_abe eps_cde A1_M^a A2_N^b A3_M^c A4_N^d
    for M in range(ND):
        for N in range(ND):
            for a, b, e in itertools.permutations(range(3)):
                for c, d in itertools.permutations([x for x in range(3) if x != e], 2):
                    T[bidx(M, a), bidx(N, b), bidx(M, c), bidx(N, d)] += -0.125 * EPS3[a, b, e] * EPS3[c, d, e]
    _T4 = T
    return T


_G4 = None


def gamma4():
    global _G4
    if _G4 is not None:
        return _G4
    T = T4_tensor()
    G = np.zeros_like(T)
    for perm in itertools.permutations(range(4)):
        G += np.transpose(T, axes=np.argsort(perm))
    _G4 = G
    return G


def tree_quadratic_A(a1, a2, p):
    """Multilinear coefficient of -1/4 (dA,dA)_Hodge on a1 e^{ipx} + a2 e^{-ipx}  (dA_{MN} only for M,N in R^4 or one index)."""
    # -1/8 sum_{MN} (f_MN, f_MN), f_MN = d_M A_N - d_N A_M ;  bilinear coefficient = -1/4 sum_{MN}(f1_MN, f2_MN)
    tot = 0
    a1 = a1.reshape(ND, NC)
    a2 = a2.reshape(ND, NC)
    for M in range(ND):
        for N in range(ND):
            f1 = (1j * p[M] if M < NV else 0) * a1[N] - (1j * p[N] if N < NV else 0) * a1[M]
            f2 = (-1j * p[M] if M < NV else 0) * a2[N] - (-1j * p[N] if N < NV else 0) * a2[M]
            tot += -0.25 * np.dot(f1, f2)
    return tot


# ---------------------------------------------------------------- fermion and ghost matrices
def Y_of(a):
    """Matrix Y with 1/2 lambda^T Y lambda = 1/2 (lambda, rho([A,lambda])), A = a (shape (10,3))."""
    a = np.asarray(a).reshape(ND, NC)
    Y = np.zeros((NF, NF), dtype=complex)
    for M in range(ND):
        for aa in range(3):
            if a[M, aa] == 0:
                continue
            for c in range(3):
                for b in range(3):
                    e = EPS3[aa, c, b]
                    if e == 0:
                        continue
                    Y[b * 16:(b + 1) * 16, c * 16:(c + 1) * 16] += a[M, aa] * e * SB[M]
    return Y


def YA(M, a):
    v = np.zeros((ND, NC))
    v[M, a] = 1
    return Y_of(v)


def ghost_V(a):
    """Antisymmetric matrix of 1/2 X^T V X = ([A|_V, c], A^*|_V) = A_mu^a c^b A^*_mu^e eps_abe."""
    a = np.asarray(a).reshape(ND, NC)
    V = np.zeros((NG, NG), dtype=complex)
    for mu in range(NV):
        for aa in range(3):
            if a[mu, aa] == 0:
                continue
            for b in range(3):
                for e in range(3):
                    m = a[mu, aa] * EPS3[aa, b, e]
                    V[b, 3 + mu * 3 + e] += m
                    V[3 + mu * 3 + e, b] -= m
    return V


# propagator numerators (hat G), G(k) = hatG(k)/k^2, regulator exp(-eps k^2)
def hatG_boson_const():
    return 2.0 * np.eye(NB, dtype=complex)


def hatG_fermion_lin():
    """list over mu of matrices: hatG_F(k) = sum_mu k_mu * out[mu]  (= i k.sigmabar)."""
    return [np.kron(np.eye(NC), 1j * SB[mu]) for mu in range(NV)]


def hatG_ghost_lin():
    out = []
    for mu in range(NV):
        g = np.zeros((NG, NG), dtype=complex)
        for b in range(3):
            g[b, 3 + mu * 3 + b] = -1j
            g[3 + mu * 3 + b, b] = -1j
        out.append(g)
    return out


# ---------------------------------------------------------------- polynomial-matrix algebra
# A PM is a dict: key (tau, m, mono) -> ndarray ; value = sum tau^t (k^2)^(-m) k^mono * matrix
def unitvec(mu):
    v = [0, 0, 0, 0]
    v[mu] = 1
    return tuple(v)


def addmono(x, y):
    return tuple(i + j for i, j in zip(x, y))


Z4 = (0, 0, 0, 0)


def pm_mul(A, B, taumax):
    out = {}
    for (ta, ma, xa), Ma in A.items():
        for (tb, mb, xb), Mb in B.items():
            t = ta + tb
            if t > taumax:
                continue
            key = (t, ma + mb, addmono(xa, xb))
            prod = Ma @ Mb
            if key in out:
                out[key] = out[key] + prod
            else:
                out[key] = prod
    return out


def pm_scale(A, s):
    return {k: v * s for k, v in A.items()}


def pm_add(*As):
    out = {}
    for A in As:
        for k, v in A.items():
            out[k] = out[k] + v if k in out else v.copy()
    return out


def inv_sq_shift(a, taumax, n):
    """Scalar PM (times identity of size n) for 1/(k + tau a)^2 expanded to tau^taumax.
    1/(k+ta)^2 = (1/k^2) sum_j (-1)^j (2 tau k.a + tau^2 a^2)^j / (k^2)^j."""
    I = np.eye(n, dtype=complex)
    # X = 2 tau k.a + tau^2 a^2 as scalar PM with m = 0
    X = {}
    for mu in range(4):
        if a[mu] != 0:
            X[(1, 0, unitvec(mu))] = 2 * a[mu] * I
    a2 = sum(x * x for x in a)
    if a2 != 0:
        X[(2, 0, Z4)] = a2 * I
    total = {(0, 1, Z4): I.copy()}
    term = {(0, 1, Z4): I.copy()}
    for j in range(1, taumax + 1):
        # term_j = term_{j-1} * (-X) / k^2
        nx = {(t, m + 1, mono): -M for (t, m, mono), M in X.items()}
        term = pm_mul(term, nx, taumax)
        if not term:
            break
        total = pm_add(total, term)
    return total


def prop_pm(kind, a, taumax):
    """G(k + tau a) as PM."""
    if kind == 'B':
        n = NB
        num = {(0, 0, Z4): hatG_boson_const()}
    else:
        lin = hatG_fermion_lin() if kind == 'F' else hatG_ghost_lin()
        n = lin[0].shape[0]
        num = {}
        for mu in range(4):
            num[(0, 0, unitvec(mu))] = lin[mu]
        const = sum(a[mu] * lin[mu] for mu in range(4))
        if np.any(const != 0):
            num[(1, 0, Z4)] = const
    den = inv_sq_shift(a, taumax, n)
    return pm_mul(num, den, taumax)


def const_pm(M, tau=0):
    return {(tau, 0, Z4): M}


def dfact2(n):  # (n-1)!! for even n, as integer: product of odd numbers < n
    r = 1
    for k in range(n - 1, 0, -2):
        r *= k
    return r


def ang_avg(mono):
    """< k^mono > over the unit 3-sphere in R^4 (exact Fraction)."""
    if any(e % 2 for e in mono):
        return Fraction(0)
    num = 1
    for e in mono:
        num *= dfact2(e)
    n = sum(mono) // 2
    den = 2 ** n
    for k in range(2, n + 2):
        den *= k
    return Fraction(num, den)


def cfrac(z):
    return (Fraction(float(z.real)), Fraction(float(z.imag)))


def log_avg(P, tau_target, trace_proj=None):
    """Angular average of the homogeneous degree -4 part (k^4 F_{-4}) at order tau_target, traced.
    Returns exact complex as (Fraction re, Fraction im)."""
    re = Fraction(0)
    im = Fraction(0)
    for (t, m, mono), M in P.items():
        if t != tau_target:
            continue
        if sum(mono) - 2 * m != -4:
            continue
        w = ang_avg(mono)
        if w == 0:
            continue
        tr = np.trace(M if trace_proj is None else trace_proj @ M)
        r, i = cfrac(tr)
        re += w * r
        im += w * i
    return (re, im)


def cdiv(a, b):
    ar, ai = a
    br, bi = b
    d = br * br + bi * bi
    return ((ar * br + ai * bi) / d, (ai * br - ar * bi) / d)


def cadd(*zs):
    return (sum((z[0] for z in zs), Fraction(0)), sum((z[1] for z in zs), Fraction(0)))


def cscale(z, s):
    return (z[0] * s, z[1] * s)


def cstr(z):
    if z[1] == 0:
        return str(z[0])
    return f"{z[0]} + {z[1]} i"
