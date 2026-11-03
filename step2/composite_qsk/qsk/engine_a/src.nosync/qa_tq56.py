"""Engine A, part Q, tasks TQ5 (special supercharges) and TQ6 (special conformal K_b): intrinsic kernels.

Along eta in Sigma^*: S^(1)(eta) = sum_mu int x^mu S^(1)(xi_mu) + int s0, xi_mu = rho(e_mu) eta = s[mu] eta (Prop
prop:qs-posdec (i)).  Position part along b: the xi-kernel at rho(b) eta (graph by graph identical by construction).
Intrinsic part  Upsilon^0_eta = sum_mu (-i d/dr_mu) [cycles with S^(1)(xi_mu), momentum r injected at the S vertex]_{r=0}
                               + [cycles with s0].
Along K_b (v(z) = z^2 b - 2 (b.z) z, z = x - x_O): S^(1)(K_b) = S^(1)(P_{v(x_O)}) - (b,x_O) S^(1)(Delta_{x_O})
   + S^(1)(M_{omega_{x_O}}) + S^(1)(K_{b,x_O}) (eq qs-Kdec); the intrinsic vertex S^(1)(K_{b,x_O}) is
   sum_rho z^rho G_rho + sum_{rho<=sigma} z^rho z^sigma G_{rho sigma},
   G_rho: v. = 1/2 om^(rho) Sigma + (w/4) dv^(rho),  om^(rho)_{mu nu} = 2 (d_{rho mu} b_nu - d_{rho nu} b_mu), dv^(rho) = -8 b_rho;
   G_{rho sigma}: v. = a^(rho sigma) . d,  a^(rho rho) = b - 2 b_rho e_rho,  a^(rho sigma) = -2 (b_rho e_sigma + b_sigma e_rho);
   z^rho -> -i d/dr_rho, z^rho z^sigma -> -d^2/dr_rho dr_sigma at r = 0 (Thm thm:qs-simplex (iii)).
The eps -> 0 coefficients are T-independent (Thm qs-simplex (iv)): route A evaluates them directly (no endpoint)."""
from fractions import Fraction
import numpy as np
import qa_s1 as S1, qa_core as QC, qa_util as U

SS = S1.SS
MI = (Fraction(0), Fraction(-1))
M1 = (Fraction(-1), Fraction(0))


def eta_gens(eta):
    return [QC.odd_gen(SS[mu] @ np.asarray(eta, complex), name='xi_mu%d' % mu) for mu in range(4)], QC.s0_gen(eta)


def eta_intrinsic(Xa, Xb, legs, eta, D, gens=None):
    gxs, g0 = gens if gens is not None else eta_gens(eta)
    tot = QC.pair_value(Xa, Xb, legs, g0)['fin']
    for mu in range(4):
        f = lambda r, mu=mu: QC.pair_value(Xa, Xb, legs, gxs[mu], r0=r)['fin']
        tot = U.vadd(tot, U.vscale(QC.deriv_r(f, mu, D + 2), MI))
    return tot


def K_gens(b):
    b = np.asarray(b, dtype=float)
    G1 = []
    for rho in range(4):
        om = np.zeros((4, 4))
        for mu in range(4):
            for nu in range(4):
                om[mu, nu] = 2 * ((rho == mu) * b[nu] - (rho == nu) * b[mu])
        G1.append(QC.even_gen('K1_%d' % rho, om=om, dv=-8 * b[rho]))
    G2 = {}
    for rho in range(4):
        for sig in range(rho, 4):
            a = np.zeros(4)
            if rho == sig:
                a = b.copy(); a[rho] -= 2 * b[rho]
            else:
                a[sig] -= 2 * b[rho]; a[rho] -= 2 * b[sig]
            if np.any(a):
                G2[(rho, sig)] = QC.even_gen('K2_%d%d' % (rho, sig), a=a)
    return G1, G2


def K_intrinsic(Xa, Xb, legs, b, D, gens=None):
    G1, G2 = gens if gens is not None else K_gens(b)
    tot = U.vzero()
    for rho in range(4):
        f = lambda r, rho=rho: QC.pair_value(Xa, Xb, legs, G1[rho], r0=r)['fin']
        tot = U.vadd(tot, U.vscale(QC.deriv_r(f, rho, D + 2), MI))
    for (rho, sig), g in G2.items():
        f = lambda r, g=g: QC.pair_value(Xa, Xb, legs, g, r0=r)['fin']
        tot = U.vadd(tot, U.vscale(QC.deriv_r2(f, rho, sig, D + 3), M1))
    return tot


def D_intrinsic(Xa, Xb, legs, D):
    """Upsilon^0_Delta by the route of the special edges (regression, TQ7)"""
    g0 = QC.even_gen('D0', dv=8)
    tot = QC.pair_value(Xa, Xb, legs, g0)['fin']
    for mu in range(4):
        g = QC.even_gen('Ds%d' % mu, a=2 * np.eye(4)[mu])
        f = lambda r, g=g: QC.pair_value(Xa, Xb, legs, g, r0=r)['fin']
        tot = U.vadd(tot, U.vscale(QC.deriv_r(f, mu, D + 2), MI))
    return tot
