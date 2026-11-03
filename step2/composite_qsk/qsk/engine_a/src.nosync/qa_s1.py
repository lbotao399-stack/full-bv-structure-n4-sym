"""Engine A, part Q (route A): vertex blocks of S^(1)(x) in the conventions of the two-letter engine A.

Conventions (inherited, read-only): step2/twoletter/engine_a/src.nosync/ea_core.py and its vendored Step-1 core r2_core.py:
plane-wave legs, all momenta incoming, vertex function = coefficient of the odd parameters in leg order, Hessian
H_jl = d_j d_l (left derivatives, l first) with the background odd parameters to the left, vertex functions scaled by 2
(divided out by the caller), colour su(2) orthonormal, [t_a, t_b] = eps_abc t_c.
Clifford data (cliff.py): rho(e_M)|_Sigma = sb[M] (Sigma -> Sigma^*), rho(e_M)|_Sigma^* = s[M] (Sigma^* -> Sigma),
pairing (u, w) = u^T w (u in Sigma, w in Sigma^*), Gamma(Q1, Q2)^M = Q1^T sb[M] Q2 (symmetric).

The functionals (theory qsk_theory.tex Def def:qs-terms, Table tab:qs-terms; eq c2-S1Q; Def c7-vdot):
  S^(1)(xi) = int dvol ( -(Gamma(xi,lam), A^*) + 1/2 (rho(F_A) xi, lam^*) ),
     F_MN = d_M A_N - d_N A_M + [A_M, A_N],  rho(F) xi = sum_{M<N} F_MN s[M] sb[N] xi,
     T1 = -(Gamma(xi,lam)|_V, A^*|_V), T2 = -(Gamma(xi,lam)|_perp, A^*|_perp), T3 = 1/2 (rho(d A|_V) xi, lam^*),
     T4 = 1/2 (rho(d A|_perp) xi, lam^*), T5 = 1/2 (rho(1/2 [A ^ A]) xi, lam^*);
  s0(eta) = -(rho(A|_perp) eta, lam^*),  eta in Sigma^*,  rho(A|_perp) eta = sum_{M>=4} A_M s[M] eta;
  even x with vector field v:  S^(1)(v) = int ( (v.A, A^*) - (v.lam, lam^*) - (v.c, c^*) ),
     v. = v^mu d_mu + 1/2 omega_{mu nu} Sigma^{mu nu} + (w/4) div v  on a slot of Weyl weight w.
Every block is an affine matrix function h(kj, kl; mom) = h0 + sum mom Hb + kj Hj + kl Hl on the 54-component
colour-stripped space of ea_core (rows: j-leg, columns: l-leg); quadratic blocks carry the identity colour, the cubic
term T5 the colour form ([Z1, Z2], Z3) stripped as the cubic vertices of I (Y = t0 at (in t1, out t2))."""
import os, sys, itertools, hashlib, warnings
import numpy as np
np.seterr(all='ignore')

SRC = os.path.dirname(os.path.abspath(__file__))
HERE = os.path.dirname(SRC)                                              # qsk/engine_a
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))       # project root
TLSRC = os.path.join(ROOT, 'step2', 'twoletter', 'engine_a', 'src.nosync')
CQSRC = os.path.join(ROOT, 'step2', 'composite_qsk', 'composite', 'engine_a', 'src.nosync')
for _p in (TLSRC, CQSRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ea_core as E          # noqa: E402
import cq_core as Q          # noqa: E402

C = E.C                       # vendored r2_core
SB = np.array(E.SGB)          # rho(e_M) on Sigma
SS = np.array(E.SG)           # rho(e_M) on Sigma^*
F3 = C.F                      # structure constants eps_abc
D = C.D
NS = E.NS
OFS = E.OFS
V4 = [0, 1, 2, 3]
WP = [4, 5, 6, 7, 8, 9]
Z4 = np.zeros(4, dtype=complex)
E4 = [np.eye(4)[m].astype(complex) for m in range(4)]


def provenance():
    out = Q.provenance()
    out['qa_s1.py'] = hashlib.sha256(open(os.path.abspath(__file__), 'rb').read()).hexdigest()
    return out


def mom10(k):
    v = np.zeros(10, dtype=complex); v[:4] = k[:4]; return v


# ------------------------------------------------------------------------------------------------ spin generators
def spin_sigma(mu, nu):
    """Sigma^{mu nu} on Sigma and on Sigma^* fixed by the equivariance of rho:
    S*^{mu nu} rho(v) - rho(v) S^{mu nu} = rho(Sigma^{mu nu} v),  (Sigma^{mu nu} v)_m = d_{m mu} v_nu - d_{m nu} v_mu"""
    Sg = 0.25 * (SS[mu] @ SB[nu] - SS[nu] @ SB[mu])
    Sgs = 0.25 * (SB[mu] @ SS[nu] - SB[nu] @ SS[mu])
    return Sg, Sgs


def check_spin():
    """returns the sign c with S^{mu nu}_Sigma = c/4 (s[mu] sb[nu] - s[nu] sb[mu]) making rho equivariant"""
    ok = {}
    for c in (1, -1):
        good = True
        for mu in range(4):
            for nu in range(4):
                if mu == nu:
                    continue
                Sg, Sgs = spin_sigma(mu, nu)
                Sg, Sgs = c * Sg, c * Sgs
                for M in range(10):
                    v = np.zeros(10); v[M] = 1
                    Sv = np.zeros(10)
                    if M == mu:
                        Sv[mu] += 0; Sv[mu] += 0
                    # (Sigma v)_m = delta_{m mu} v_nu - delta_{m nu} v_mu
                    Sv = np.zeros(10); Sv[mu] += v[nu]; Sv[nu] -= v[mu]
                    lhs = Sgs @ SB[M] - SB[M] @ Sg
                    rhs = sum(Sv[N] * SB[N] for N in range(10))
                    if not np.allclose(lhs, rhs):
                        good = False
        ok[c] = good
    return ok


SPIN_C = None


def spin(mu, nu):
    global SPIN_C
    if SPIN_C is None:
        ok = check_spin()
        assert sum(ok.values()) == 1, ok
        SPIN_C = 1 if ok[1] else -1
    Sg, Sgs = spin_sigma(mu, nu)
    return SPIN_C * Sg, SPIN_C * Sgs


# ------------------------------------------------------------------------------------------------ canonical functions
def _X(xi):
    """X[M, N] = s[M] sb[N] xi (16-vectors in Sigma), M != N"""
    X = np.zeros((10, 10, 16), dtype=complex)
    for M in range(10):
        for N in range(10):
            if M != N:
                X[M, N] = SS[M] @ (SB[N] @ xi)
    return X


def odd_terms(xi, which=('T1', 'T2', 'T3', 'T4', 'T5')):
    """scaled (x2) canonical vertex functions of the terms of S^(1)(xi), xi in Sigma (16-vector)"""
    xi = np.asarray(xi, dtype=complex)
    X = _X(xi)
    G = np.einsum('A,MAB->MB', xi, SB)          # Gamma(xi, lam)^M = xi^T sb[M] lam  ->  G[M] . lam
    T = {}
    def t12(Ms):
        GM = G[Ms]
        def fn(pols, moms):
            l, As = pols
            return -2 * np.einsum('MB,iBe,jMe->ij', GM, l, As[:, Ms, :], optimize=True)
        return fn
    if 'T1' in which:
        T['T1'] = (['l', 'As'], t12(V4))
    if 'T2' in which:
        T['T2'] = (['l', 'As'], t12(WP))
    def t34(Ns):
        # i sum_{M in V, N in Ns, M != N} k_M A_N X[M,N] . lam^*
        def fn(pols, moms):
            A, ls = pols
            k = mom10(moms[0])
            Xs = np.zeros((10, 16), dtype=complex)        # sum_M k_M X[M, N] for N in Ns
            for N in Ns:
                for M in V4:
                    if M != N:
                        Xs[N] += k[M] * X[M, N]
            return 1j * np.einsum('Na,iNe,jae->ij', Xs, A, ls, optimize=True)
        return fn
    if 'T3' in which:
        T['T3'] = (['A', 'ls'], t34(V4))
    if 'T4' in which:
        T['T4'] = (['A', 'ls'], t34(WP))
    if 'T5' in which:
        def fn5(pols, moms):
            A1, A2, ls = pols
            # sum_{M != N} f_pqe A1_M^p A2_N^q X[M,N] . lam^*^e
            return np.einsum('iMp,jNq,pqe,MNa,kae->ijk', A1, A2, F3, X, ls, optimize=True)
        T['T5'] = (['A', 'A', 'ls'], fn5)
    return T


def s0_term(eta):
    """scaled (x2) -(rho(A|_perp) eta, lam^*),  eta in Sigma^*"""
    eta = np.asarray(eta, dtype=complex)
    Y = np.zeros((10, 16), dtype=complex)
    for M in WP:
        Y[M] = SS[M] @ eta
    def fn(pols, moms):
        A, ls = pols
        return -2 * np.einsum('Ma,iMe,jae->ij', Y, A, ls, optimize=True)
    return {'s0': (['A', 'ls'], fn)}


WEIGHT = {'A': 1.0, 'l': 1.5, 'c': 0.0}


def even_terms(a=None, om=None, dv=0.0, slots=('A', 'l', 'c'), Wperp=None):
    """scaled (x2) canonical functions of int ( (v.A, A^*) - (v.lam, lam^*) - (v.c, c^*) ) for the constant-coefficient
    vector-field data at the vertex: v.f = a^mu d_mu f + 1/2 om_{mu nu} Sigma^{mu nu} f + (w/4) dv f, with a (4-vector),
    om (4x4 antisymmetric), dv (number) the value of v, of its rotation part and of div v at the vertex.
    Wperp: optional 6x6 antisymmetric R (so(W_R^perp) rotation, (v.A)_w = R_{w w'} A_{w'}, on lam the spinor action)."""
    a = np.zeros(4) if a is None else np.asarray(a, dtype=complex)
    om = np.zeros((4, 4)) if om is None else np.asarray(om, dtype=complex)
    # vector action on A_M (10-dim): (v.A)_M = a.d A_M + (om A)_M (M in V) + (w/4) dv A_M (+ R on perp)
    RA = np.zeros((10, 10), dtype=complex)
    for m in range(4):
        for n in range(4):
            RA[m, n] += om[m, n]          # 1/2 om_{mu nu} (Sigma^{mu nu} A)_m = om_{m nu} A_nu
    if Wperp is not None:
        RA[4:, 4:] += np.asarray(Wperp, dtype=complex)
    RA += 0.25 * WEIGHT['A'] * dv * np.eye(10)
    Rl = np.zeros((16, 16), dtype=complex)
    for mu in range(4):
        for nu in range(4):
            if om[mu, nu] != 0:
                Sg, _ = spin(mu, nu)
                Rl += 0.5 * om[mu, nu] * Sg
    if Wperp is not None:
        for i in range(6):
            for j in range(6):
                if Wperp[i][j] != 0:
                    Sg, _ = spin(4 + i, 4 + j)
                    Rl += 0.5 * Wperp[i][j] * Sg
    Rl += 0.25 * WEIGHT['l'] * dv * np.eye(16)
    T = {}
    if 'A' in slots:
        def fA(pols, moms):
            A, As = pols
            k = moms[0]
            der = 1j * complex(np.dot(a, k[:4]))
            return 2 * np.einsum('MN,iNe,jMe->ij', RA + der * np.eye(10), A, As, optimize=True)
        T['vA'] = (['A', 'As'], fA)
    if 'l' in slots:
        def fl(pols, moms):
            l, ls = pols
            k = moms[0]
            der = 1j * complex(np.dot(a, k[:4]))
            return -2 * np.einsum('ab,ibe,jae->ij', Rl + der * np.eye(16), l, ls, optimize=True)
        T['vl'] = (['l', 'ls'], fl)
    if 'c' in slots:
        def fc(pols, moms):
            c, cs = pols
            k = moms[0]
            der = 1j * complex(np.dot(a, k[:4]))
            return -2 * der * np.einsum('ie,je->ij', c, cs, optimize=True)
        T['vc'] = (['c', 'cs'], fc)
    return T


# ------------------------------------------------------------------------------------------------ Hessians
def V_S(legs, terms):
    """as r2_core.V with the canonical functions `terms` {name: (types, fn)}"""
    types = [L[0] for L in legs]
    out = 0
    for nm, (can, fn) in terms.items():
        if sorted(can) != sorted(types):
            continue
        slot = [None] * len(legs); used = set()
        for i, t in enumerate(types):
            for j, ct in enumerate(can):
                if ct == t and j not in used:
                    slot[i] = j; used.add(j); break
        pols = [None] * len(can); moms = [None] * len(can)
        for i, j in enumerate(slot):
            pols[j] = legs[i][1]; moms[j] = legs[i][2]
        val = np.transpose(fn(pols, moms), slot)
        odd_req = [slot[i] for i in range(len(legs)) if C.PARITY[types[i]]]
        sgn = 1
        for x in range(len(odd_req)):
            for y in range(x + 1, len(odd_req)):
                if odd_req[x] > odd_req[y]:
                    sgn = -sgn
        out = out + sgn * val
    return out


def hessian_S(terms, bgs, kj, kl):
    H = np.zeros((C.NTOT, C.NTOT), dtype=complex)
    for tj in C.TYPES:
        for tl in C.TYPES:
            legs = [(tj, C.ident(tj), kj), (tl, C.ident(tl), kl)] + [(t, pol[None], m) for (t, pol, m) in bgs]
            v = V_S(legs, terms)
            if isinstance(v, int):
                continue
            v = v.reshape(C.SIZE[tj], C.SIZE[tl])
            sg = (-1) ** (C.PARITY[tj] * C.PARITY[tl])
            H[C.OFF[tj]:C.OFF[tj] + C.SIZE[tj], C.OFF[tl]:C.OFF[tl] + C.SIZE[tl]] += sg * v
    return H


def _round(M):
    R = np.round(M.real) + 1j * np.round(M.imag)
    assert np.allclose(M, R, atol=1e-9), np.abs(M - R).max()
    return R


def quad_block(terms):
    """two-leg block (identity colour) of quadratic terms: Aff(h0, Hj, Hl) on the 54-space (scaled x2)"""
    def H(kj, kl):
        return _round(E._strip_fast(hessian_S(terms, [], kj, kl), 0, 0))
    h0 = H(Z4, Z4)
    Hj = [H(E4[m], Z4) - h0 for m in range(4)]
    Hl = [H(Z4, E4[m]) - h0 for m in range(4)]
    t = H(np.array([1, -1, 2, 0], complex), np.array([0, 1, 1, -2], complex))
    pred = h0 + sum(c * Hj[m] for m, c in enumerate([1, -1, 2, 0])) + sum(c * Hl[m] for m, c in enumerate([0, 1, 1, -2]))
    assert np.allclose(t, pred)
    # colour diagonal check
    H01 = E._strip_fast(hessian_S(terms, [], Z4, Z4), 0, 1)
    assert np.allclose(H01, 0)
    return Q.Aff(h0, Hj, Hl, 'S2')


_CUB = {}


def cubic_block(terms_key, terms, ltype, pol, mom):
    """block of a cubic term (T5) with output of letter type ltype, stripped polarisation pol, momentum mom; colour
    (ad_Y t_in, t_out) as the cubic vertices of I.  Cached per unit polarisation."""
    slot = E.SLOT[ltype]
    h0 = np.zeros((NS, NS), complex)
    for comp in range(E.SZ[slot]):
        cpol = pol[comp]
        if cpol == 0:
            continue
        key = (terms_key, slot, comp)
        if key not in _CUB:
            v = np.zeros(E.SZ[slot]); v[comp] = 1
            P = E._full_pol(slot, v, 0)
            Hm = _round(E._strip_fast(hessian_S(terms, [(slot, P, Z4)], Z4, Z4), 1, 2))
            # T5 has no derivatives: check momentum independence
            Ht = _round(E._strip_fast(hessian_S(terms, [(slot, P, E4[1])], E4[2], E4[3]), 1, 2))
            assert np.allclose(Hm, Ht)
            _CUB[key] = Hm
        h0 = h0 + cpol * _CUB[key]
    z = [np.zeros((NS, NS), complex) for _ in range(4)]
    return Q.Aff(h0, z, [x.copy() for x in z], 'S3')


# ------------------------------------------------------------------------------------------------ chirality data
G5S = np.asarray(E.G5S)
G5SS = np.asarray(E.G5SS)


def chir_basis(name):
    return [np.asarray(v, dtype=int) for v in E.BASIS[name]]
