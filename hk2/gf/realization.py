"""
realization.py -- explicit realisation of V_10, Sigma, Sigma^*, rho, the pairing and Gamma,
built from the two-component formulas of the manuscript (Definition c8-twocomp).
Explicit bases live only in code; reported formulas are index-free.

Index maps (code):
  Sigma   : xi_alpha^i  -> 4*alpha + i        (0..7),   xibar_{alphadot i} -> 8 + 4*alphadot + i  (8..15)
  Sigma^* : eta^alpha_i -> 4*alpha + i        (0..7),   etabar^{alphadot i} -> 8 + 4*alphadot + i (8..15)
  V_10    : e_m, m = 0..3 (V);  e_I, I = 4..9 (W_R^perp, the 't Hooft label I-3 = 1..6)
"""
import sympy as sp
from functools import lru_cache

I_ = sp.I
s1 = sp.Matrix([[0, 1], [1, 0]])
s2 = sp.Matrix([[0, -I_], [I_, 0]])
s3 = sp.Matrix([[1, 0], [0, -1]])
one2 = sp.eye(2)
# Euclidean sigma matrices of Definition c8-twocomp
SIG = [s1, s2, s3, I_ * one2]                 # sigma^m_{alpha alphadot}
SIGB = [-s1, -s2, -s3, I_ * one2]             # sigmabar^{m alphadot alpha}  (row alphadot, column alpha)
EPS_UP = sp.Matrix([[0, 1], [-1, 0]])         # epsilon^{12} = +1
EPS_DN = sp.Matrix([[0, -1], [1, 0]])         # epsilon_{12} = -1


def _levi4(i, j, k, l):
    return sp.LeviCivita(i, j, k, l)


@lru_cache(None)
def thooft():
    """Sigma^I_{ij}, I = 0..5 (labels 1..6), i,j = 0..3 (labels 1..4)."""
    out = []
    for a in range(3):
        M = sp.zeros(4, 4)
        for b in range(3):
            for c in range(3):
                M[b, c] = sp.LeviCivita(a, b, c)
            M[b, 3] = 1 if a == b else 0
            M[3, b] = -1 if a == b else 0
        out.append(M)
    for a in range(3):
        M = sp.zeros(4, 4)
        for b in range(3):
            for c in range(3):
                M[b, c] = I_ * sp.LeviCivita(a, b, c)
            M[b, 3] = -I_ if a == b else 0
            M[3, b] = I_ if a == b else 0
        out.append(M)
    return tuple(out)


def xi(al, i):
    return 4 * al + i


def xib(ad, i):
    return 8 + 4 * ad + i


@lru_cache(None)
def rho_blocks():
    """For each M = 0..9: (R, Rs) with R : Sigma -> Sigma^* (16x16) and Rs : Sigma^* -> Sigma (16x16)."""
    TH = thooft()
    blocks = []
    for M in range(10):
        R = sp.zeros(16, 16)    # (rho(e_M) xi) in Sigma^*
        Rs = sp.zeros(16, 16)   # (rho(e_M) eta) in Sigma
        if M < 4:
            m = M
            # (rho(v)xi)^alpha_i = i v_m xibar_{ad i} sigmabar^{m ad alpha}
            # (rho(v)xi)^{ad i} = -i v_m sigmabar^{m ad alpha} xi_alpha^i
            for i in range(4):
                for al in range(2):
                    for ad in range(2):
                        R[xi(al, i), xib(ad, i)] += I_ * SIGB[m][ad, al]
                        R[xib(ad, i), xi(al, i)] += -I_ * SIGB[m][ad, al]
            # (rho(v)eta)_alpha^i = -i v_m sigma^m_{alpha ad} etabar^{ad i}
            # (rho(v)eta)_{ad i} = i v_m eta^alpha_i sigma^m_{alpha ad}
            for i in range(4):
                for al in range(2):
                    for ad in range(2):
                        Rs[xi(al, i), xib(ad, i)] += -I_ * SIG[m][al, ad]
                        Rs[xib(ad, i), xi(al, i)] += I_ * SIG[m][al, ad]
        else:
            Iw = M - 4
            wdn = TH[Iw]                       # w_{ij}
            wup = sp.zeros(4, 4)               # w^{ij} = 1/2 eps^{ijkl} w_{kl}
            for i in range(4):
                for j in range(4):
                    wup[i, j] = sp.Rational(1, 2) * sum(_levi4(i, j, k, l) * wdn[k, l]
                                                        for k in range(4) for l in range(4))
            # (rho(w)xi)^alpha_j = eps^{alpha beta} w_{ji} xi_beta^i
            # (rho(w)xi)^{ad j} = eps^{ad bd} w^{ji} xibar_{bd i}
            for j in range(4):
                for i in range(4):
                    for al in range(2):
                        for be in range(2):
                            R[xi(al, j), xi(be, i)] += EPS_UP[al, be] * wdn[j, i]
                            R[xib(al, j), xib(be, i)] += EPS_UP[al, be] * wup[j, i]
            # (rho(w)eta)_alpha^k = - w^{kj} eps_{alpha beta} eta^beta_j
            # (rho(w)eta)_{bd k} = - w_{kj} eps_{bd ad} etabar^{ad j}
            for k in range(4):
                for j in range(4):
                    for al in range(2):
                        for be in range(2):
                            Rs[xi(al, k), xi(be, j)] += -wup[k, j] * EPS_DN[al, be]
                            Rs[xib(al, k), xib(be, j)] += -wdn[k, j] * EPS_DN[al, be]
        blocks.append((R, Rs))
    return tuple(blocks)


@lru_cache(None)
def Cpair():
    """(xi, eta) = xi^T C eta, C = diag(+1 x 8, -1 x 8)."""
    return sp.diag(*([1] * 8 + [-1] * 8))


def rho_full(M):
    """32x32 matrix of rho(e_M) on Sigma (+) Sigma^* (first 16 = Sigma)."""
    R, Rs = rho_blocks()[M]
    F = sp.zeros(32, 32)
    F[16:32, 0:16] = R
    F[0:16, 16:32] = Rs
    return F


def rho_vec(vV, vW=None):
    """(R, Rs) for v = sum_m vV[m] e_m + sum_I vW[I] e_I."""
    R = sp.zeros(16, 16)
    Rs = sp.zeros(16, 16)
    blocks = rho_blocks()
    for m in range(4):
        if vV[m] != 0:
            R += vV[m] * blocks[m][0]
            Rs += vV[m] * blocks[m][1]
    if vW is not None:
        for Iw in range(6):
            if vW[Iw] != 0:
                R += vW[Iw] * blocks[4 + Iw][0]
                Rs += vW[Iw] * blocks[4 + Iw][1]
    return R, Rs


def Gamma(z1, z2):
    """Gamma(z1,z2)^M = (rho(e_M) z1, z2) for z1, z2 both in Sigma (16-vectors) or both in Sigma^*.
    Returned as a list of 10 components (V first, then W_R^perp)."""
    C = Cpair()
    blocks = rho_blocks()
    out = []
    for M in range(10):
        R, Rs = blocks[M]
        out.append(sp.expand((z2.T * C * (R * z1))[0, 0]))
    return out


def Gamma_star(z1, z2):
    C = Cpair()
    blocks = rho_blocks()
    out = []
    for M in range(10):
        R, Rs = blocks[M]
        out.append(sp.expand(((Rs * z1).T * C * z2)[0, 0]))
    return out


def verify_realization():
    res = {}
    ok = True
    F = [rho_full(M) for M in range(10)]
    for M in range(10):
        for N in range(M, 10):
            X = F[M] * F[N] + F[N] * F[M] - (2 if M == N else 0) * sp.eye(32)
            if any(sp.simplify(x) != 0 for x in X):
                ok = False
    res['clifford_V10'] = ok
    C = Cpair()
    ok = True
    for M in range(10):
        R, Rs = rho_blocks()[M]
        if any(sp.simplify(x) != 0 for x in (C * R - (C * R).T)):
            ok = False
        if any(sp.simplify(x) != 0 for x in (Rs.T * C - C * Rs)):
            ok = False
    res['rho_symmetric_for_pairing'] = ok
    chi = F[0] * F[1] * F[2] * F[3]
    plus = all(sp.simplify(chi[xi(al, i), xi(al, i)] - 1) == 0 for al in range(2) for i in range(4))
    res['chirality_rho(e1e2e3e4)=+1_on_S+W'] = bool(plus)
    TH = thooft()
    ok = True
    for i in range(4):
        for j in range(4):
            for k in range(4):
                for l in range(4):
                    s = sum(TH[Iw][i, j] * TH[Iw][k, l] for Iw in range(6))
                    if sp.simplify(s - 2 * _levi4(i, j, k, l)) != 0:
                        ok = False
    res['thooft_sum_SigmaSigma=2eps'] = ok
    return res


if __name__ == '__main__':
    print(verify_realization())
