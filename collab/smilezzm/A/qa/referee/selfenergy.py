"""Referee check (3): the one-loop scalar self-energy of N=4 SYM with the heat-kernel cutoff, derived independently.

Conventions (Euclidean, weight e^{-S_E}, Feynman gauge, coupling g, su(N), f^{acd} f^{bcd} = N delta^{ab}):
  S_E = int 1/4 F^a_MN F^a_MN + 1/2 (d_m A^a_M)^2 + fermions,  F^a_MN = d_M A_N - d_N A_M + g f^{abc} A^b_M A^c_N,
  d_I = 0 (dimensional reduction), A_I = scalars, I = 5..10.  Propagator delta^{ab} delta_MN e^{-eps k^2}/k^2.
  Sigma defined by  dressed line = Dhat + Dhat Sigma Dhat  (Sigma = sum of 1PI graphs with vertex factors -S_E).

Part A (algebra, from the action itself, numerically with explicit su(2), su(3) structure constants and
        explicit 10d gamma matrices):
   * cubic and quartic vertex functions obtained by polarisation of the cubic/quartic terms of -S_E;
   * seagull contraction factors (4-dim A|_V loop, 6-dim A_perp loop);
   * bubble numerator;  * gaugino trace tr_16(G^I G^K G^J G^L);
   * fermion-loop sign check by supersymmetry: with an identical cutoff on every loop momentum the
     quadratic divergences of the scalar self-energy must cancel (seagulls + bubble + gaugino).
Part B (referee's reduction): with Khat(k)=e^{-eps k^2}=k^2 Dhat(k) the four graphs reduce exactly to
     T = Deps(0) = 1/(16 pi^2 eps),  B(p) = int Dhat(l) Dhat(p-l),  C(p) = int Dhat(l) Khat(p+l)
   bubble  = g^2 N [2 p^2 B + C],   gaugino = g^2 N 4[2C - p^2 B],   seagulls = -9 g^2 N T,
   C(p) = (e^{-eps p^2/2} - e^{-eps p^2})/(16 pi^2 eps^2 p^2) exactly (one Gaussian + elementary integral),
   B(p) evaluated by a ONE-dimensional position-space Hankel transform, 16 pi^2 B(p=1) = 4 int_0^inf (1-e^{-r^2/4eps})^2 J_1(r)/r^2 dr.
Part C: direct numerical evaluation of the bubble and gaugino momentum integrals at many eps (two-centre (|l|,|l+p|)
     reduction, mpmath), and a free fit a/eps + b log eps + c + (eps log eps, eps, eps^2 log eps, eps^2).
Usage: python3 selfenergy.py [dps]
"""
import itertools
import json
import os
import sys
import time

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from common import fit_linear, heron  # noqa: E402

dps = int(sys.argv[1]) if len(sys.argv) > 1 else 22
mp.mp.dps = dps
out = {}
rng = np.random.default_rng(12345)


# ============================== Part A ===================================
def su_structure(N):
    gens = []
    for i in range(N):
        for j in range(i + 1, N):
            m = np.zeros((N, N), complex); m[i, j] = m[j, i] = 0.5; gens.append(m)
            m = np.zeros((N, N), complex); m[i, j] = -0.5j; m[j, i] = 0.5j; gens.append(m)
    for k in range(1, N):
        m = np.zeros((N, N), complex)
        for i in range(k):
            m[i, i] = 1
        m[k, k] = -k
        gens.append(m / np.sqrt(2 * k * (k + 1)))
    d = len(gens)
    f = np.zeros((d, d, d))
    for a, b, c in itertools.product(range(d), repeat=3):
        f[a, b, c] = np.real(-2j * np.trace((gens[a] @ gens[b] - gens[b] @ gens[a]) @ gens[c]))
    return f


def cubic_density(fields, f):
    """-S_E cubic term at one point for plane waves: fields = [(A (dimG x 10 array), k (10-vector, k_I = 0))].
    -S_E,3 = - f^{abc} (d_M A^a_N) A^b_M A^c_N  with d -> i k (coupling g = 1)."""
    (A1, k1), (A2, k2), (A3, k3) = fields
    dA1 = 1j * np.einsum('M,aN->aMN', k1, A1)
    return -np.einsum('abc,aMN,bM,cN->', f, dA1, A2, A3)


def V3(fields, f):
    """vertex function = symmetrised derivative = sum over the 3! assignments of the fields to the slots."""
    return sum(cubic_density([fields[i] for i in perm], f) for perm in itertools.permutations(range(3)))


def quartic_density(As, f):
    """-S_E quartic: -(1/4) f^{abe} f^{cde} A^a_M A^b_N A^c_M A^d_N (g = 1)."""
    A1, A2, A3, A4 = As
    return -0.25 * np.einsum('abe,cde,aM,bN,cM,dN->', f, f, A1, A2, A3, A4)


def V4(As, f):
    return sum(quartic_density([As[i] for i in perm], f) for perm in itertools.permutations(range(4)))


def unit_field(dimG, a, M):
    A = np.zeros((dimG, 10)); A[a, M] = 1.0
    return A


partA = {}
for N in (2, 3):
    f = su_structure(N)
    dimG = f.shape[0]
    cas = np.einsum('acd,bcd->ab', f, f)
    partA['su%d_casimir_is_N' % N] = bool(np.allclose(cas, N * np.eye(dimG)))
    # seagull: (1/2) sum over internal (c,K) of V4(ext a I, ext b J, c K, c K), K in V (0..3) or in W_perp (4..9)
    a, b, I, J = 0, 0, 5, 5
    for name, Ks in (('V', range(4)), ('Wperp', range(4, 10))):
        s = 0.0
        for c in range(dimG):
            for K in Ks:
                s += V4([unit_field(dimG, a, I), unit_field(dimG, b, J), unit_field(dimG, c, K), unit_field(dimG, c, K)], f)
        partA['su%d_seagull_%s_factor_over_minus_N' % (N, name)] = float(np.real(s) / (-N))
    # off-diagonal flavour / colour must vanish
    s = sum(V4([unit_field(dimG, 0, 5), unit_field(dimG, 0, 6), unit_field(dimG, c, K), unit_field(dimG, c, K)], f)
            for c in range(dimG) for K in range(10))
    partA['su%d_seagull_offdiag_flavour' % N] = float(abs(s))
    # bubble numerator from the cubic vertex functions; internal line X carries l (v1 -> v2), Y carries p-l
    errs = []
    for trial in range(3):
        p = np.zeros(10); p[:4] = rng.normal(size=4)
        l = np.zeros(10); l[:4] = rng.normal(size=4)
        tot = 0.0
        for cx in range(dimG):
            for Mx in range(10):
                for cy in range(dimG):
                    for My in range(10):
                        X, Y = unit_field(dimG, cx, Mx), unit_field(dimG, cy, My)
                        v1 = V3([(unit_field(dimG, 0, 5), p), (X, -l), (Y, -(p - l))], f)
                        if v1 == 0:
                            continue
                        v2 = V3([(unit_field(dimG, 0, 5), -p), (X, l), (Y, p - l)], f)
                        tot += v1 * v2
        tot *= 0.5
        expect = 0.5 * N * (np.dot(2 * p - l, 2 * p - l) + np.dot(p + l, p + l))  # gluon on l, or gluon on p-l
        errs.append(abs(tot - expect))
    partA['su%d_bubble_numerator_maxerr' % N] = float(max(errs))
    # V4 normalisation: V4(A,A,A,A) = 24 * quartic(A)
    Ar = rng.normal(size=(dimG, 10))
    partA['su%d_V4_polarisation_check' % N] = float(abs(V4([Ar] * 4, f) - 24 * quartic_density([Ar] * 4, f)))

# 10d Euclidean gamma matrices (hermitian, {G^M,G^N} = 2 delta), chirality projection to 16 dims
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
i2 = np.eye(2)


def kron(*ms):
    r = np.array([[1]], complex)
    for m in ms:
        r = np.kron(r, m)
    return r


G = []
for k in range(5):
    pre = [s3] * k
    post = [i2] * (4 - k)
    G.append(kron(*(pre + [s1] + post)))
    G.append(kron(*(pre + [s2] + post)))
G = np.array(G)
clif = max(abs(G[m] @ G[n] + G[n] @ G[m] - 2 * (m == n) * np.eye(32)).max() for m in range(10) for n in range(10))
G11 = np.eye(32, dtype=complex)
for m in range(10):
    G11 = G11 @ G[m]
G11 = G11 * (1j)**(-5) if abs((G11 @ G11 - np.eye(32))).max() > 1e-9 else G11
if abs((G11 @ G11 - np.eye(32))).max() > 1e-9:
    G11 = G11 * 1j
Pp = (np.eye(32) + G11) / 2
tr16 = lambda M: np.trace(Pp @ M)
maxdev = 0.0
for I, J in itertools.product(range(4, 10), repeat=2):
    for K, L in itertools.product(range(4), repeat=2):
        val = tr16(G[I] @ G[K] @ G[J] @ G[L])
        maxdev = max(maxdev, abs(val - (-16) * (I == J) * (K == L)))
partA['clifford_err'] = float(clif)
partA['chiral_trace_tr16(GI GK GJ GL)+16dIJdKL_maxdev'] = float(maxdev)
partA['tr16(1)'] = float(np.real(np.trace(Pp)))

# gaugino numerator: Sigma_F^{bb'}_{IJ} = (1/2) * (-1)_loop * sum colour * tr16[(-f Gamma^I) S(l) (-f Gamma^J) S(l-p)],
# S(k) = (i Gamma.k)^{-1} = -i Gamma.k / k^2 (numerators only here); colour: vertex psibar^a phi^b psi^c : -f^{abc} Gamma^I
f3 = su_structure(3); N = 3
p = np.zeros(10); p[:4] = rng.normal(size=4)
l = np.zeros(10); l[:4] = rng.normal(size=4)
Sl = -1j * np.einsum('m,mij->ij', l, G)
Slp = -1j * np.einsum('m,mij->ij', l - p, G)
colour = np.einsum('abc,cda->bd', f3, f3)   # f^{abc} f^{c b' a}
spin = tr16(G[5] @ Sl @ G[5] @ Slp)
SigmaF = 0.5 * (-1) * colour[0, 0] * spin   # the two (-1)'s of the vertices give +1
partA['gaugino_numerator_over_N_l.(l-p)'] = float(np.real(SigmaF / (N * np.dot(l, l - p))))   # expect +8
# SUSY check: quadratic divergences with an identical cutoff on |l|: coefficients of int 1/l^2
seag = 0.5 * (-partA['su3_seagull_V_factor_over_minus_N'] - partA['su3_seagull_Wperp_factor_over_minus_N'])
bub = 1.0      # (2p+l)^2/(l^2 (p+l)^2) -> 1/l^2 ; both assignments already counted (sum = 2 * 1/2)
ferm = partA['gaugino_numerator_over_N_l.(l-p)']
partA['quadratic_divergence_identical_cutoff_seagulls_bubble_gaugino'] = [seag, bub, ferm, seag + bub + ferm]
out['partA'] = partA
print(json.dumps(partA, indent=1), flush=True)

# ============================== Part B ===================================
pi = mp.pi
g_E, L2 = mp.euler, mp.log(2)


def b_hankel(eps):
    """16 pi^2 B(p=1) = 4 int_0^inf (1 - e^{-r^2/4eps})^2 J1(r)/r^2 dr  (position space, Hankel)."""
    se = mp.sqrt(eps)
    pts = [0] + [se * 2**k for k in range(-3, 7)]
    head = mp.quad(lambda r: (-mp.expm1(-r * r / (4 * eps)))**2 * mp.besselj(1, r) / r**2, pts)
    a = pts[-1]
    # beyond a, (1-e^{-r^2/4eps})^2 = 1 - O(e^{-16}) ; keep the exact factor anyway up to a few zeros
    mid_end = mp.besseljzero(1, int(a / pi) + 3)
    mid = mp.quad(lambda r: (-mp.expm1(-r * r / (4 * eps)))**2 * mp.besselj(1, r) / r**2, [a, mid_end])
    tail = mp.quadosc(lambda r: mp.besselj(1, r) / r**2, [mid_end, mp.inf], zeros=lambda n: mp.besseljzero(1, n + int(a / pi) + 3))
    return 4 * (head + mid + tail)


rowsB = []
for es in ['1e-2', '3e-3', '1e-3', '3e-4', '1e-4', '3e-5', '1e-5', '1e-6']:
    eps = mp.mpf(es)
    bh = b_hankel(eps)
    claim = 1 - g_E - L2 - mp.log(eps)
    rowsB.append(dict(eps=es, b_hankel=mp.nstr(bh, 16), b_claimed_M002=mp.nstr(claim, 16),
                      diff=mp.nstr(bh - claim, 8), diff_over_epslogeps=mp.nstr((bh - claim) / (eps * mp.log(eps)), 6)))
    print(rowsB[-1], flush=True)
X = [[mp.log(mp.mpf(r['eps'])), 1, mp.mpf(r['eps']) * mp.log(mp.mpf(r['eps'])), mp.mpf(r['eps'])] for r in rowsB]
Y = [mp.mpf(r['b_hankel']) for r in rowsB]
sol = fit_linear(X, Y)
out['B_hankel'] = dict(rows=rowsB, fit_log=mp.nstr(sol[0], 14), fit_const=mp.nstr(sol[1], 14),
                       claimed_const=mp.nstr(1 - g_E - L2, 14))
print(out['B_hankel']['fit_log'], out['B_hankel']['fit_const'], out['B_hankel']['claimed_const'], flush=True)


# ============================== Part C ===================================
def mom_two_centre(F, eps):
    """int d^4l/(2pi)^4 F(|l|, |l+p|), |p| = 1: two-centre reduction, outer r=|l|, inner s=|l+p| in [|1-r|, 1+r]."""
    ise = 1 / mp.sqrt(eps)
    pts = [mp.mpf(0), mp.mpf(1) / 4, mp.mpf(1) / 2, mp.mpf(1), mp.mpf(3) / 2, 2, 4]
    k = 3
    while 2**k < 12 * ise:
        pts.append(mp.mpf(2)**k); k += 1
    pts += [12 * ise, mp.inf]

    def inner(r):
        lo, hi = abs(1 - r), 1 + r
        return mp.quad(lambda s: 2 * pi * r * s * heron(r, s) * F(r, s), [lo, hi])
    tot = mp.mpf(0)
    for a_, b_ in zip(pts[:-1], pts[1:]):
        tot += mp.quad(inner, [a_, b_])
    return tot / (2 * pi)**4


def C_exact(eps, p2=1):
    y = eps * p2
    return (mp.e**(-y / 2) - mp.e**(-y)) / (16 * pi**2 * eps**2 * p2)


units = 16 * pi**2          # results below in units of g^2 N / (16 pi^2)
rowsC = []
for es in ['2e-2', '1e-2', '5e-3', '2e-3', '1e-3', '5e-4', '2e-4', '1e-4']:
    eps = mp.mpf(es)
    t0 = time.time()
    # r = |l| (gluon, A|_V), s = |l+p| (scalar);  numerator (2p+l)^2 = 2 + 2 s^2 - r^2 ;  l.(l+p) = (r^2+s^2-1)/2
    bub = mom_two_centre(lambda r, s: (2 + 2 * s * s - r * r) * mp.e**(-eps * (r * r + s * s)) / (r * r * s * s), eps) * units
    gau = mom_two_centre(lambda r, s: 8 * (r * r + s * s - 1) / 2 * mp.e**(-eps * (r * r + s * s)) / (r * r * s * s), eps) * units
    Cn = mom_two_centre(lambda r, s: mp.e**(-eps * (r * r + s * s)) / (r * r), eps)
    rowsC.append(dict(eps=es, bubble=mp.nstr(bub, 20), gaugino=mp.nstr(gau, 20),
                      C_numeric=mp.nstr(Cn, 16), C_exact=mp.nstr(C_exact(eps), 16), seconds=round(time.time() - t0, 1)))
    print(rowsC[-1], flush=True)

fitsC = {}
for name in ('bubble', 'gaugino'):
    X, Y = [], []
    for r in rowsC:
        e = mp.mpf(r['eps'])
        X.append([1 / e, mp.log(e), 1, e * mp.log(e), e, e * e * mp.log(e), e * e])
        Y.append(mp.mpf(r[name]))
    sol = fit_linear(X, Y)
    fitsC[name] = dict(inv_eps=mp.nstr(sol[0], 12), log_eps=mp.nstr(sol[1], 12), const=mp.nstr(sol[2], 12))
# predictions at p = 1 (log p^2 = 0):
pred = {
    'bubble_const_referee(13/8 - 2g - 2log2)': mp.nstr(mp.mpf(13) / 8 - 2 * g_E - 2 * L2, 12),
    'bubble_const_note(167/72 - 2g - 2log2)': mp.nstr(mp.mpf(167) / 72 - 2 * g_E - 2 * L2, 12),
    'gaugino_const_referee(-7 + 4g + 4log2)': mp.nstr(-7 + 4 * g_E + 4 * L2, 12),
    'gaugino_const_note(-13/9 + 4g + 4log2)': mp.nstr(-mp.mpf(13) / 9 + 4 * g_E + 4 * L2, 12),
}
# total renormalised Sigma at p=1 from the direct numerics (seagulls exact: -9/eps)
tot_const = mp.mpf(fitsC['bubble']['const']) + mp.mpf(fitsC['gaugino']['const'])
out['direct_momentum'] = dict(rows=rowsC, fits=fitsC, predictions=pred,
                              total_inv_eps=mp.nstr(-9 + mp.mpf(fitsC['bubble']['inv_eps']) + mp.mpf(fitsC['gaugino']['inv_eps']), 12),
                              total_log_eps=mp.nstr(mp.mpf(fitsC['bubble']['log_eps']) + mp.mpf(fitsC['gaugino']['log_eps']), 12),
                              Sigma_ren_const_minus_2g_2log2=mp.nstr(tot_const - 2 * g_E - 2 * L2, 12),
                              claimed=mp.nstr(-mp.mpf(43) / 8, 12))
print(json.dumps({k: v for k, v in out['direct_momentum'].items() if k != 'rows'}, indent=1))
with open(os.path.join(HERE, 'selfenergy_dps%d.json' % dps), 'w') as fh:
    json.dump(out, fh, indent=1, default=str)
