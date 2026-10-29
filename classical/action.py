"""
action.py -- component expansion of the reduced gauge-BV action [ESW (31)] on R^4.

Densities are returned multiplied by SCALE = 8 so that all coefficients are Gaussian integers.
kappa = 1 in the jet densities (the pairing (T_i,T_j) = kappa delta_ij scales S by kappa and the
bracket by 1/kappa; {S,S} = 0 is homogeneous, so kappa = 1 is no loss for the CME).

Literal reading of (31) (pairings written with their left entry first):
  -1/4 (F_A,F_A),  (F,F) = sum_{M<N} (F_MN, F_MN)   [Hodge pairing: (F,F) dvol = F ^ *F]
  +1/2 (lam, dslash_A lam) = 1/2 sum_M lam^a sig^M_ab (D_M lam)^b
  +(d_A c, A^*)            = sum_M (D_M c)^i A^*_M^i
  +([c,lam], lam^*)        = sum_a [c,lam]^{a,i} lam^*_a^i
  +1/2 ([c,c], c^*)        = 1/2 [c,c]^i c^*^i
with [X,Y]^i = f_jki X^j Y^k, D_M = d_M + [A_M, -], d_M = 0 for M >= 4.
"""
import itertools
from jets import Vars, var, mono, pmul, padd_into, pscale, const, ZERO4
import n4_model as nm

SCALE = 8


def e_mu(mu):
    J = [0, 0, 0, 0]
    J[mu] = 1
    return tuple(J)


def build_action(V, fint, dimg, include=('S2', 'I3', 'I4')):
    """returns dict name -> density (scaled by SCALE)."""
    f = fint
    out = {'S2': {}, 'I3': {}, 'I4': {}}
    A = lambda M, i, J=ZERO4: V.get('A', M, i, J)
    lam = lambda a, i, J=ZERO4: V.get('lam', a, i, J)
    c = lambda i, J=ZERO4: V.get('c', 0, i, J)
    As = lambda M, i: V.get('As', M, i)
    lams = lambda a, i: V.get('lams', a, i)
    cs = lambda i: V.get('cs', 0, i)
    V4 = nm.V4
    fl = list(f.items())

    # ---------------- bosonic: -1/4 (F,F)_Hodge = -1/8 sum_{MN} F_MN F_MN ; x8 -> - sum F F
    def fMN(M, N, i):
        P = {}
        if M in V4:
            padd_into(P, var(V, 'A', N, i, e_mu(M)))
        if N in V4:
            padd_into(P, var(V, 'A', M, i, e_mu(N)), (-1, 0))
        return P

    def bMN(M, N, i):
        P = {}
        for (j, k, l), val in fl:
            if l == i:
                padd_into(P, mono(V, (val, 0), [A(M, j), A(N, k)]))
        return P

    for M in range(10):
        for N in range(10):
            if M == N:
                continue
            for i in range(dimg):
                fm = fMN(M, N, i)
                bm = bMN(M, N, i)
                padd_into(out['S2'], pmul(fm, fm), (-1, 0))
                padd_into(out['I3'], pmul(fm, bm), (-2, 0))
                padd_into(out['I4'], pmul(bm, bm), (-1, 0))

    # ---------------- fermion: 1/2 (lam, dslash_A lam); x8 -> 4 lam^a sig^M_ab (D_M lam)^b
    for M in range(10):
        ent = nm.sig_entries(M)
        for (a, b, s) in ent:
            coef = (4 * s[0], 4 * s[1])
            for i in range(dimg):
                if M in V4:
                    padd_into(out['S2'], mono(V, coef, [lam(a, i), lam(b, i, e_mu(M))]))
                for (j, k, l), val in fl:
                    if l == i:
                        padd_into(out['I3'], mono(V, (coef[0] * val, coef[1] * val), [lam(a, i), A(M, j), lam(b, k)]))

    # ---------------- ghost: (d_A c, A^*); x8
    for M in range(10):
        for i in range(dimg):
            if M in V4:
                padd_into(out['S2'], mono(V, (8, 0), [c(i, e_mu(M)), As(M, i)]))
            for (j, k, l), val in fl:
                if l == i:
                    padd_into(out['I3'], mono(V, (8 * val, 0), [A(M, j), c(k), As(M, i)]))
    # ([c,lam], lam^*); x8
    for a in range(16):
        for (j, k, l), val in fl:
            padd_into(out['I3'], mono(V, (8 * val, 0), [c(j), lam(a, k), lams(a, l)]))
    # 1/2 ([c,c], c^*); x8 -> 4
    for (j, k, l), val in fl:
        padd_into(out['I3'], mono(V, (4 * val, 0), [c(j), c(k), cs(l)]))

    return {k: v for k, v in out.items() if k in include}


def total(parts):
    S = {}
    for P in parts.values():
        padd_into(S, P)
    return S


def families(dimg):
    fam = []
    for (name, dim) in (('c', 1), ('A', 10), ('lam', 16)):
        for a in range(dim):
            for i in range(dimg):
                fam.append((name, a, i))
    return fam
