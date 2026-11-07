"""Engine B -- the graph classes T, SE, B, QU, OC (exact) and the numerical confirmations, plus the EX
integrand.  All values are relative to the tree 4 hbar dim(g) S12 D_0(x)^2, at |x|^2 = 1 = ell^2, in units
(unit = -g^2 N/(32 pi^2) = -lambda_hat/2).  Derivations: notes.txt sections S (SE), Q (QU), E (EX).
"""
import math
import mpmath as mp
import numpy as np
import sympy as sp

EG = sp.EulerGamma
L2 = sp.log(2)
L3 = sp.log(3)
eps = sp.Symbol('epsilon', positive=True)


# ------------------------------------------------------------------ SE: exact pieces
def se_exact():
    """Self-energy line (both lines, with the counterterm vertex -I^CT and the seagulls).
    Sigma = -9 lh/eps delta + g^2 C_A [2 Laplacian(D_eps^2) + 9 D_eps K_eps],  Sing Sigma = lh[-9/(2eps) + 2 p^2 log eps].
    Line:  L = L_i + L_ii,
      L_i  = -32 pi^2 lh [D_eps^2 * D_2eps](x) - 2 lh log eps D_2eps(x),
      L_ii = 9 lh int_eps^inf dt (t+eps)^-2 int int_{[eps,inf)^2} dt1 dt4 [K_{t1+t4+c(t)}(x) - K_{t1+t4}(x)].
    [D_eps^2 * D_2eps](x) = D_0(x)/(16 pi^2) * I_J(eps),
      I_J = int_0^{2/eps} dw rho(w) w^-2 (1 - exp(-w/(4(1+2 eps w)))),  rho(w) = min(w, 2/eps - w)  (x^2 = 1).
    """
    Y = sp.Symbol('Y', positive=True)
    y = sp.Symbol('y', positive=True)
    w = sp.Symbol('w', positive=True)
    t = sp.Symbol('t', positive=True)
    # int_0^Y (1 - e^-y)/y dy = log Y + EulerGamma + E1(Y)   (sympy: Ein)
    ein = sp.limit(sp.integrate((1 - sp.exp(-y)) / y, (y, 0, Y)) - sp.log(Y), Y, sp.oo)
    region2 = sp.simplify(sp.integrate((2 / eps - w) / w ** 2, (w, 1 / eps, 2 / eps)))
    # I_J = -log(4 eps) + ein + region2 + o(1)
    IJ_const = sp.simplify(-sp.log(4) + ein + region2)
    kappa = sp.simplify(IJ_const)  # I_J = -log eps + kappa
    # L_i / D_0 = -2 lh (I_J) + 2 lh log eps = -2 lh kappa  (+ -2 lh log x^2 at general x)
    # L_ii / D_0 = -9 lh * eps * int_eps^inf t/(t+eps)^3 dt
    m = sp.simplify(sp.integrate(eps * t / (t + eps) ** 3, (t, eps, sp.oo)))
    Lii = -9 * m
    # relative, both lines:  2 L / D_0 in units of lh;  unit = -lh/2
    rel_lh = 2 * (-2 * kappa + Lii)
    c_S = sp.nsimplify(sp.simplify(rel_lh * (-2)))
    d_S = sp.Integer(2 * (-2) * (-2))  # 2 lines * (-2 lh) * log x^2 -> units * (-2)
    return dict(c_S=sp.simplify(c_S), d_S=d_S, kappa=kappa, Lii_over_lhD0=Lii, ein=ein, region2=region2,
                log_eps_without_ct=sp.Integer(-8), log_eps_ct=sp.Integer(8))


def se_numeric(eps_list):
    """Finite-eps evaluation of I_J(eps) and of L_ii/(lh D_0) (mpmath, 30 digits); returns the deviations from
    the asymptotic forms -log eps + kappa and -27/8."""
    mp.mp.dps = 30
    kappa = mp.euler + 1 - 3 * mp.log(2)
    out = []
    for e in eps_list:
        e = mp.mpf(e)
        Lam = 1 / e

        def fJ(w):
            r = w if w <= Lam else 2 * Lam - w
            return r / w ** 2 * (-mp.expm1(-w / (4 * (1 + 2 * e * w))))
        pts = [mp.mpf(0)] + [mp.mpf(2) ** k for k in range(-6, 200) if mp.mpf(2) ** k < Lam] + [Lam, 1.5 * Lam, 2 * Lam]
        IJ = mp.quad(fJ, pts)

        # L_ii: Phi(c) = -int_{2e}^{2e+c} (T-2e) K_T dT - c int_{2e+c}^inf K_T dT,  K_T(1) = e^{-1/4T}/(16 pi^2 T^2)
        def K(Tv):
            return mp.exp(-1 / (4 * Tv)) / (16 * mp.pi ** 2 * Tv ** 2)

        def Phi(c):
            a = -mp.quad(lambda Tv: (Tv - 2 * e) * K(Tv), [2 * e, 2 * e + c])
            b = -c * (-mp.expm1(-1 / (4 * (2 * e + c)))) / (4 * mp.pi ** 2)
            return a + b
        Lii = 9 * mp.quad(lambda tv: Phi(tv * e / (tv + e)) / (tv + e) ** 2, [e, 10 * e, 1000 * e, mp.inf])
        D0 = 1 / (4 * mp.pi ** 2)
        out.append(dict(eps=float(e), IJ=float(IJ), IJ_dev=float(IJ - (-mp.log(e) + kappa)),
                        Lii=float(Lii / D0), Lii_dev=float(Lii / D0 + mp.mpf(27) / 8)))
    return out


# ------------------------------------------------------------------ QU: exact
def qu_exact():
    """QU = g^2 C_A dim g (S12 - 2 T12) int d^4z D_eps(x-z)^2 D_eps(z)^2  (census weight).
    Q/D_0^2 = I2/(256 pi^2),  I2 = int int dw1 dw2 rho(w1) rho(w2) (w1+w2)^-2 exp(-w1 w2/(4(w1+w2))) (x^2 = 1),
    I2 = 32 int_0^{1/2} du G(u/(4 eps))/(u(1-u)) + 2*16 int_{1/eps}^{2/eps} (2/eps - w)/w^2 dw + o(1),
    G(y) = 1 - (1+y) e^-y,  int_0^Y G(y)/y dy = log Y + EulerGamma - 1 + o(1)."""
    Y = sp.Symbol('Y', positive=True)
    y = sp.Symbol('y', positive=True)
    w = sp.Symbol('w', positive=True)
    u = sp.Symbol('u', positive=True)
    G = 1 - (1 + y) * sp.exp(-y)
    CG = sp.limit(sp.integrate(G / y, (y, 0, Y)) - sp.log(Y), Y, sp.oo)
    regA = 32 * (sp.log(1 / (8 * eps)) + CG + sp.integrate(1 / (1 - u), (u, 0, sp.Rational(1, 2))))
    regB = 2 * 16 * sp.integrate((2 / eps - w) / w ** 2, (w, 1 / eps, 2 / eps))
    I2 = sp.expand(sp.simplify(sp.expand_log(regA + regB, force=True)))
    # relative = g^2 C_A (S12 - 2 T12) I2/(256 pi^2) = lh (S12 - 2T12) I2/16 ;  units: * (-2)
    rel_units = sp.expand(-2 * I2 / 16)  # coefficient of (S12 - 2 T12)
    a = rel_units.coeff(sp.log(eps))
    c = sp.simplify(rel_units - a * sp.log(eps))
    return dict(I2=I2, CG=CG, log_eps_S=a, log_eps_T=-2 * a, c_S=c, c_T=sp.simplify(-2 * c),
                d_S=-a, d_T=2 * a)


def qu_numeric(eps_list):
    """I2(eps) at finite eps from the exact one-dimensional representation (inner W-integral in closed form),
    mpmath 30 digits; deviation from 32(-log eps + EulerGamma - 3 log 2)."""
    mp.mp.dps = 30
    out = []
    for e in eps_list:
        e = mp.mpf(e)
        Lam = 1 / e

        def P1(a, A, B):
            f = lambda Wv: -(Wv / a + 1 / a ** 2) * mp.exp(-a * Wv)
            return f(B) - f(A)

        def P0(a, A, B):
            return (mp.exp(-a * A) - mp.exp(-a * B)) / a

        def Pm1(a, A, B):
            return mp.e1(a * A) - mp.e1(a * B)

        def inner(u):
            # u in (0, 1/2]; w1 = W u <= w2 = W (1-u)
            a = u * (1 - u) / 4
            v = 1 - u
            s = u * v * P1(a, 0, Lam / v)
            B2 = min(2 * Lam / v, Lam / u)
            s += u * (2 * Lam * P0(a, Lam / v, B2) - v * P1(a, Lam / v, B2))
            if Lam / u < 2 * Lam / v:
                A3, B3 = Lam / u, 2 * Lam / v
                s += 4 * Lam ** 2 * Pm1(a, A3, B3) - 2 * Lam * P0(a, A3, B3) + u * v * P1(a, A3, B3)
            return s
        pts = [mp.mpf(0)] + [mp.mpf(2) ** k * e for k in range(-4, 400) if mp.mpf(2) ** k * e < 0.5] + [mp.mpf(1) / 3, mp.mpf(1) / 2]
        I2 = 2 * mp.quad(inner, pts)
        ex = 32 * (-mp.log(e) + mp.euler - 3 * mp.log(2))
        out.append(dict(eps=float(e), I2=float(I2), dev=float(I2 - ex)))
    return out


# ------------------------------------------------------------------ OC, B, T
def oc_exact(g2):
    """Trees with -d_delta I^CT[eps](O_a), a = 1, 2; d_delta I^CT(O_ab) = unit log eps [alpha O_ab + beta (a,b) K]
    (alpha = 4, beta = -2 from G2).  tree(O_ab, O_cd) ~ S12, tree(K, O_cd) ~ 2 (c,d)."""
    alpha, beta = g2['total']['alpha'], g2['total']['beta']
    per_end_S = -alpha
    per_end_T = -beta * 2
    return dict(log_eps_S=2 * per_end_S, log_eps_T=2 * per_end_T, c_S=0, c_T=0)


def b_exact():
    """-F_0 = -beta_2 unit int(dA_perp, dA_perp) = 4 beta_2 unit (-1/4 int (dA,dA)): on each regulated line
    D_eps -> D_eps - 4 beta_2 unit D_2eps."""
    b2 = sp.Symbol('beta_2')
    return dict(c_S=-8 * b2, c_T=0, d_S=0, d_T=0, log_eps_S=0, log_eps_T=0)
