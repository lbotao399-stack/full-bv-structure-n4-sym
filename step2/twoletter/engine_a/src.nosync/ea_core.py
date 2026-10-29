"""Engine A (R8, route A): exact colour-stripped chains through two letters, Schwinger-parameter Gaussian integral and
exact simplex / boundary-face integrals (twoletter_theory.tex, Thm thm:tl-simplex).

Conventions (inherited from the Step-1 engine A, engine_a/r2_core.py, verified A1-A4 in checks_engines_results.json):
  plane-wave legs, all momenta incoming; vertex function = coefficient of the odd parameters in leg order; Hessian
  H_jl = d_j d_l I (left derivatives, l first) at the background, background odd parameters to the left; cubic vertex
  functions scaled by 2, quartic by 8 (divided out exactly); one loop W_1 = sum_n 1/(2n) str[(P S H)^n],
  S = parity sign, P(eps,infty)(k) = int_eps^infty dt e^{-t k^2} Phat(k), Phat: A-A 2, c-A*_V and A*_V-c -i k_mu,
  lambda-lambda i rho(k) (each times (-,-)^{-1}).
Colour stripping.  Every cubic vertex is (kinematic) x ([Z_1,Z_2],Z_3); with output Y, cycle slots in (row, j-leg) and
  out (column, l-leg) it is h3[in,out] (ad_Y t_in, t_out).  The quartic vertex is h4[in,out] (ad_Y ad_Y' t_in, t_out)
  + h4'[in,out] (ad_Y' ad_Y t_in, t_out).  Both are read off exactly from the su(2) Hessians of r2_core
  ([t_a,t_b] = eps_abc t_c, orthonormal): cubic Y = t0 -> entry (in t1, out t2); quartic Y = t0, Y' = t1 ->
  word ad_Y ad_Y' at (in t0, out t1), word ad_Y' ad_Y at (in t1, out t0).
Chain and colour word.  With the operator Hessian at position 0 and blocks 1..k,
  str[P_0 S H_M P_1 S h_1 ... P_k S h_k] = sum_{r,c} H_M[r,c] T[c,r],   T = P_1 S h_1 ... P_k S h_k S P_0 S;
  the column (l-leg) letter of H_M is letter a, the row (j-leg) letter is letter b, the blocks are met in the order
  1..k from a to b, and the colour element of the chain is sum_gamma t^gamma (x) W(t_gamma),
  W = ad_{Y_k} ... ad_{Y_1} (Y_1 adjacent to a), first factor in the colour slot of a (twoletter_theory Lemma
  lem:tl-colour (i)).  A quartic block with outputs (o_j, o_{j+1}) in this order contributes ad_{o_{j+1}} ad_{o_j}.
Momenta.  Line i (i = 0..k) carries k_i = k + q_i, q_0 = 0, q_i = q_{i-1} - P_{i-1}, P_i the momentum entering block i
  (P_0 = momentum of omega = minus the sum of the outputs); line i enters block i as its j-leg (incoming -k_i) and
  leaves block i-1 as its l-leg (incoming k_i).  Letter b carries incoming momentum -k, letter a carries k - P_0.
All arithmetic exact: matrices are Gaussian integers in complex128 with the 2^52 bound asserted on every product;
rational weights are Fractions; results are Gaussian rationals (re, im) of Fractions."""
import os, sys, itertools, math, pickle
from fractions import Fraction
import numpy as np

SRC = os.path.dirname(os.path.abspath(__file__))                 # engine_a/src.nosync (excluded from iCloud sync)
HERE = os.path.dirname(SRC)                                       # engine_a
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))      # project root
# Verbatim copies of the Step-1 engine A modules engine_a/r2_core.py and engine_a/cliff.py (vendor/): the originals live
# on an iCloud-synced Desktop and were evicted (dataless) during this run; run_all compares the copies with the originals
# byte by byte whenever the originals are readable.
sys.path.insert(0, os.path.join(SRC, 'vendor'))
import r2_core as C                     # su(2) Step-1 engine A core (vertex functions, Hessians)
from cliff import s as SG, sb as SGB    # rho(v)|_Sigma* = s[M] v^M, rho(v)|_Sigma = sb[M] v^M

assert C.D == 3
TYPES = ['c', 'A', 'As', 'l', 'ls', 'cs']
SZ = {'c': 1, 'A': 10, 'As': 10, 'l': 16, 'ls': 16, 'cs': 1}
PAR = {'c': 1, 'A': 0, 'As': 1, 'l': 1, 'ls': 0, 'cs': 0}
OFS = {}
_o = 0
for _t in TYPES:
    OFS[_t] = _o; _o += SZ[_t]
NS = _o                                  # 54 colour-stripped components
SDIAG = np.concatenate([np.full(SZ[t], (-1) ** PAR[t]) for t in TYPES]).astype(complex)
# active components (those met by a nonzero component of P): c, A (10), A*|_V (4), lambda (16)
ACT = list(range(OFS['c'], OFS['c'] + 1)) + list(range(OFS['A'], OFS['A'] + 10)) + \
      list(range(OFS['As'], OFS['As'] + 4)) + list(range(OFS['l'], OFS['l'] + 16))
NA = len(ACT)
POS = {g: i for i, g in enumerate(ACT)}
SACT = SDIAG[ACT]


def full_index(t, comp, colour):
    """index into r2_core's colour-full field space"""
    return C.OFF[t] + comp * C.D + colour


def mm(A, B):
    return C.mm(A, B)


def gq(z):
    """exact Gaussian-integer complex -> (Fraction, Fraction)"""
    r, i = C.to_int(z)
    return (Fraction(r), Fraction(i))

# ------------------------------------------------------------------------------------------------ chirality bases
with np.errstate(all='ignore'):
    G5S = np.rint((SG[0] @ SGB[1] @ SG[2] @ SGB[3]).real).astype(int)     # rho(v1v2v3v4) on Sigma
    G5SS = np.rint((SGB[0] @ SG[1] @ SGB[2] @ SG[3]).real).astype(int)    # on Sigma^*
assert np.array_equal(G5S @ G5S, np.eye(16, dtype=int)) and np.array_equal(G5SS @ G5SS, np.eye(16, dtype=int))


def _int_basis(P2):
    """8 independent integer columns of 2*projector"""
    cols = []
    for j in range(16):
        v = P2[:, j]
        if not v.any():
            continue
        M = np.array(cols + [v], dtype=float)
        if np.linalg.matrix_rank(M) == len(cols) + 1:
            cols.append(v)
        if len(cols) == 8:
            break
    assert len(cols) == 8
    return [np.array(c, dtype=int) for c in cols]


BASIS = {
    'lamP': _int_basis(np.eye(16, dtype=int) + G5S),      # S_+ (x) W   (+1 on Sigma)
    'lamM': _int_basis(np.eye(16, dtype=int) - G5S),      # S_- (x) W^* (-1 on Sigma)
    'lamsP': _int_basis(np.eye(16, dtype=int) + G5SS),    # S_+ (x) W^* (+1 on Sigma^*)
    'lamsM': _int_basis(np.eye(16, dtype=int) - G5SS),    # S_- (x) W   (-1 on Sigma^*)
}
BASIS['AV'] = [np.eye(10, dtype=int)[M] for M in range(4)]
BASIS['AW'] = [np.eye(10, dtype=int)[M] for M in range(4, 10)]
BASIS['AsV'] = [np.eye(10, dtype=int)[M] for M in range(4)]
BASIS['AsW'] = [np.eye(10, dtype=int)[M] for M in range(4, 10)]
BASIS['c'] = [np.array([1])]
BASIS['cs'] = [np.array([1])]
SLOT = {'AV': 'A', 'AW': 'A', 'AsV': 'As', 'AsW': 'As', 'lamP': 'l', 'lamM': 'l', 'lamsP': 'ls', 'lamsM': 'ls',
        'c': 'c', 'cs': 'cs'}
WEYL = {'c': Fraction(0), 'AV': Fraction(1), 'AW': Fraction(1), 'lamP': Fraction(3, 2), 'lamM': Fraction(3, 2),
        'AsV': Fraction(3), 'AsW': Fraction(3), 'lamsP': Fraction(5, 2), 'lamsM': Fraction(5, 2), 'cs': Fraction(4)}
DEG = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1, 'AsW': -1, 'lamsP': -1, 'lamsM': -1, 'cs': -2}
PARITY_T = {k: PAR[v] for k, v in SLOT.items()}


def embed(t, v):
    """stripped component vector (length SZ[slot]) of a letter type from its basis coefficients / raw vector"""
    sl = SLOT[t]
    out = np.zeros(SZ[sl], dtype=complex)
    out[:] = v
    return out

# ------------------------------------------------------------------------------------------------ vertex blocks
E4 = [np.eye(4, dtype=int)[m] for m in range(4)]
Z4 = np.zeros(4, dtype=complex)


def _strip(Hfull, cin, cout):
    out = np.zeros((NS, NS), dtype=complex)
    for t1 in TYPES:
        for t2 in TYPES:
            for a in range(SZ[t1]):
                ia = full_index(t1, a, cin)
                row = Hfull[ia]
                for b in range(SZ[t2]):
                    v = row[full_index(t2, b, cout)]
                    if v != 0:
                        out[OFS[t1] + a, OFS[t2] + b] = v
    return out


def _strip_fast(Hfull, cin, cout):
    idx_in = np.array([full_index(t, a, cin) for t in TYPES for a in range(SZ[t])])
    idx_out = np.array([full_index(t, a, cout) for t in TYPES for a in range(SZ[t])])
    return Hfull[np.ix_(idx_in, idx_out)].copy()


def _full_pol(slot, v, colour):
    shp = C.SHAPE[slot]
    x = np.zeros(shp, dtype=complex)
    if slot in ('c', 'cs'):
        x[colour] = v[0]
    else:
        x[:, colour] = v
    return x


_CACHE_FILE = os.path.join(HERE, 'out.nosync', 'blocks_cache.pkl')
_BLK = None


def _build_blocks():
    """affine data for unit polarisations: cubic[(slot, comp)] = (h0, Hb[4], Hj[4], Hl[4]) with
    h(kj, kl; mom) = h0 + sum mom_mu Hb_mu + kj_mu Hj_mu + kl_mu Hl_mu ; quartic[(M,N)] = (w1, w2) constant
    (word ad_Y ad_Y' and ad_Y' ad_Y for Y = e_M, Y' = e_N)."""
    cub = {}
    for slot in TYPES:
        for comp in range(SZ[slot]):
            v = np.zeros(SZ[slot]); v[comp] = 1
            pol = _full_pol(slot, v, 0)
            def H(mom, kj, kl):
                return _strip_fast(C.hessian([(slot, pol, mom)], kj, kl), 1, 2)
            h0 = H(Z4, Z4, Z4)
            Hb = [H(E4[m].astype(complex), Z4, Z4) - h0 for m in range(4)]
            Hj = [H(Z4, E4[m].astype(complex), Z4) - h0 for m in range(4)]
            Hl = [H(Z4, Z4, E4[m].astype(complex)) - h0 for m in range(4)]
            # affinity check
            test = H(np.array([1, -1, 2, 0], complex), np.array([0, 1, 1, -2], complex), np.array([2, 0, -1, 1], complex))
            pred = h0 + sum(c * Hb[m] for m, c in enumerate([1, -1, 2, 0])) + sum(c * Hj[m] for m, c in enumerate([0, 1, 1, -2])) \
                + sum(c * Hl[m] for m, c in enumerate([2, 0, -1, 1]))
            assert np.array_equal(test, pred), (slot, comp)
            if np.abs(h0).sum() + sum(np.abs(x).sum() for x in Hb + Hj + Hl) == 0:
                continue
            cub[(slot, comp)] = (h0, Hb, Hj, Hl)
    qua = {}
    for M in range(10):
        for N in range(10):
            vY = np.zeros(10); vY[M] = 1; vY2 = np.zeros(10); vY2[N] = 1
            bgs = [('A', _full_pol('A', vY, 0), Z4), ('A', _full_pol('A', vY2, 1), Z4)]
            Hf = C.hessian(bgs, np.array([1, 2, -1, 3], complex), np.array([0, -2, 1, 1], complex))
            Hf0 = C.hessian(bgs, Z4, Z4)
            assert np.array_equal(Hf, Hf0)
            w1 = _strip_fast(Hf0, 0, 1)
            w2 = _strip_fast(Hf0, 1, 0)
            qua[(M, N)] = (w1, w2)
    return {'cub': cub, 'qua': qua}


def blocks():
    global _BLK
    if _BLK is None:
        if os.path.exists(_CACHE_FILE):
            with open(_CACHE_FILE, 'rb') as f:
                _BLK = pickle.load(f)
        else:
            _BLK = _build_blocks()
            with open(_CACHE_FILE, 'wb') as f:
                pickle.dump(_BLK, f)
    return _BLK


class Affine:
    """h(kj, kl) = h0 + sum kj_mu Hj_mu + kl_mu Hl_mu, restricted to the active components"""
    __slots__ = ['h0', 'Hj', 'Hl']

    def __init__(self, h0, Hj, Hl):
        self.h0, self.Hj, self.Hl = h0, Hj, Hl


def _restrict(X):
    return X[np.ix_(ACT, ACT)]


def cubic_block(ltype, pol, mom):
    """colour-stripped Hessian block of a cubic vertex with one output of letter type ltype, stripped polarisation
    `pol` (vector of the slot, Gaussian-integer entries), incoming momentum mom (integer 4-vector)"""
    B = blocks()['cub']
    slot = SLOT[ltype]
    h0 = np.zeros((NA, NA), complex); Hj = [np.zeros((NA, NA), complex) for _ in range(4)]
    Hl = [np.zeros((NA, NA), complex) for _ in range(4)]
    for comp in range(SZ[slot]):
        c = pol[comp]
        if c == 0 or (slot, comp) not in B:
            continue
        b0, Hb, bj, bl = B[(slot, comp)]
        h0 = h0 + c * _restrict(b0 + sum(mom[m] * Hb[m] for m in range(4)))
        for m in range(4):
            Hj[m] = Hj[m] + c * _restrict(bj[m]); Hl[m] = Hl[m] + c * _restrict(bl[m])
    return Affine(h0, Hj, Hl)


def quartic_block(polY, polY2):
    """word ad_Y ad_Y' (Y' met first, i.e. nearer to letter a): both outputs A, pol 10-vectors"""
    Q = blocks()['qua']
    h0 = np.zeros((NA, NA), complex)
    for M in range(10):
        if polY[M] == 0:
            continue
        for N in range(10):
            if polY2[N] == 0:
                continue
            h0 = h0 + polY[M] * polY2[N] * _restrict(Q[(M, N)][0])
    zero = [np.zeros((NA, NA), complex) for _ in range(4)]
    return Affine(h0, zero, zero)

# ------------------------------------------------------------------------------------------------ propagator
def _line_mats():
    """Phat(k) = P0 + sum_mu k_mu P1[mu] on the active components (colour stripped), per line type"""
    out = {}
    P0 = np.zeros((NS, NS), complex)
    for M in range(4):
        P0[OFS['A'] + M, OFS['A'] + M] = 2
    out['AV'] = (P0[np.ix_(ACT, ACT)], [np.zeros((NA, NA), complex)] * 4)
    P0 = np.zeros((NS, NS), complex)
    for M in range(4, 10):
        P0[OFS['A'] + M, OFS['A'] + M] = 2
    out['AP'] = (P0[np.ix_(ACT, ACT)], [np.zeros((NA, NA), complex)] * 4)
    P1 = []
    for mu in range(4):
        X = np.zeros((NS, NS), complex)
        X[OFS['c'], OFS['As'] + mu] = -1j
        X[OFS['As'] + mu, OFS['c']] = -1j
        P1.append(X[np.ix_(ACT, ACT)])
    out['gh'] = (np.zeros((NA, NA), complex), P1)
    P1 = []
    for mu in range(4):
        X = np.zeros((NS, NS), complex)
        X[OFS['l']:OFS['l'] + 16, OFS['l']:OFS['l'] + 16] = 1j * np.asarray(SG[mu])
        P1.append(X[np.ix_(ACT, ACT)])
    out['l'] = (np.zeros((NA, NA), complex), P1)
    return out


LINES = _line_mats()
LINE_NAMES = ['AV', 'AP', 'gh', 'l']

# ------------------------------------------------------------------------------------------------ polynomial matrices
def monos(deg, n=4):
    out = []
    for d in range(deg + 1):
        for c in itertools.combinations_with_replacement(range(n), d):
            e = [0] * n
            for i in c:
                e[i] += 1
            out.append(tuple(e))
    return out


def padd(A, B):
    out = dict(A)
    for e, M in B.items():
        out[e] = out[e] + M if e in out else M
    return out


def pmul(A, B):
    out = {}
    for e1, M1 in A.items():
        for e2, M2 in B.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            P = mm(M1, M2)
            out[e] = out[e] + P if e in out else P
    return {e: M for e, M in out.items() if np.any(M)}


def affine_poly(const, lin):
    """const + sum_mu k_mu lin[mu] as polynomial matrix in k"""
    P = {}
    if np.any(const):
        P[(0, 0, 0, 0)] = const
    for mu in range(4):
        if np.any(lin[mu]):
            P[tuple(E4[mu])] = lin[mu]
    return P


def line_poly(q, lines=None):
    """Phat(k + q) restricted to the given line types (default all)"""
    lines = LINE_NAMES if lines is None else lines
    const = np.zeros((NA, NA), complex); lin = [np.zeros((NA, NA), complex) for _ in range(4)]
    for L in lines:
        P0, P1 = LINES[L]
        const = const + P0 + sum(q[mu] * P1[mu] for mu in range(4))
        for mu in range(4):
            lin[mu] = lin[mu] + P1[mu]
    return affine_poly(const, lin)


def block_poly(blk, qin, qout):
    """S h(kj = -(k+qin), kl = k+qout) as polynomial matrix in k (S applied from the left)"""
    const = blk.h0 - sum(qin[mu] * blk.Hj[mu] for mu in range(4)) + sum(qout[mu] * blk.Hl[mu] for mu in range(4))
    lin = [blk.Hl[mu] - blk.Hj[mu] for mu in range(4)]
    S = SACT[:, None]
    return affine_poly(S * const, [S * x for x in lin])


def chain_T(blks, Ps, line_seq=None):
    """T(k) = P_1 S h_1 ... P_k S h_k S P_0 S as polynomial matrix in the loop momentum k (active components).
    blks: list of Affine blocks (1..k); Ps: list of their incoming momenta (integer 4-vectors).
    Returns (T, qs) with qs[i] the shift of line i (k_i = k + qs[i]), P_0 = -sum(Ps)."""
    kk = len(blks)
    P0 = -sum((np.asarray(p) for p in Ps), np.zeros(4, dtype=int))
    allP = [P0] + [np.asarray(p) for p in Ps]
    qs = [np.zeros(4, dtype=int)]
    for i in range(1, kk + 1):
        qs.append(qs[-1] - allP[i - 1])
    line_seq = line_seq or [None] * (kk + 1)
    I = {(0, 0, 0, 0): np.eye(NA, dtype=complex)}
    T = I
    for i in range(1, kk + 1):
        T = pmul(T, line_poly(qs[i], line_seq[i]))
        qn = qs[i + 1] if i < kk else qs[0]
        T = pmul(T, block_poly(blks[i - 1], qs[i], qn))
    S = SACT
    T = {e: M * S[None, :] for e, M in T.items()}
    T = pmul(T, line_poly(qs[0], line_seq[0]))
    T = {e: M * S[None, :] for e, M in T.items()}
    return T, qs

# ------------------------------------------------------------------------------------------------ exact algebra
def q_add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def q_mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def q_scale(a, r):
    return (a[0] * r, a[1] * r)


Q0 = (Fraction(0), Fraction(0))


def lap_poly(P):
    """Laplacian in k of a polynomial {exp: value} (values any ring supporting * int and +)"""
    out = {}
    for e, v in P.items():
        for mu in range(4):
            if e[mu] >= 2:
                e2 = list(e); e2[mu] -= 2; e2 = tuple(e2)
                w = v * (e[mu] * (e[mu] - 1))
                out[e2] = out[e2] + w if e2 in out else w
    return out


# scalar polynomials in alpha: dict exps(E) -> Fraction
def sp_mul(A, B):
    out = {}
    for e1, a in A.items():
        for e2, b in B.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            out[e] = out.get(e, 0) + a * b
    return {e: v for e, v in out.items() if v != 0}


def sp_pow(A, n, E):
    out = {tuple([0] * E): Fraction(1)}
    for _ in range(n):
        out = sp_mul(out, A)
    return out


def dirichlet(exps):
    E = len(exps)
    num = 1
    for n in exps:
        num *= math.factorial(n)
    return Fraction(num, math.factorial(sum(exps) + E - 1))


_REG = {}


def region_monomial(exps, e):
    """int_{R_e} alpha^exps d alpha (R_e = {alpha_e = min}, cone from the barycentre over the facet alpha_e = 0)"""
    key = (exps, e)
    if key in _REG:
        return _REG[key]
    E = len(exps)
    others = [g for g in range(E) if g != e]
    # alpha_g = mu_g + mu_b / E (g != e), alpha_e = mu_b / E ; Jacobian 1/E
    tot = Fraction(0)
    ranges = [range(exps[g] + 1) for g in others]
    for ks in itertools.product(*ranges):
        coef = Fraction(1)
        nb = exps[e]
        mexps = []
        for g, kg in zip(others, ks):
            coef *= math.comb(exps[g], kg)
            nb += kg
            mexps.append(exps[g] - kg)
        coef *= Fraction(1, E ** nb)
        tot += coef * dirichlet(tuple(mexps + [nb]))
    tot *= Fraction(1, E)
    _REG[key] = tot
    return tot


_MIN2 = {}


def specialregion_monomial(exps, s, kpow):
    """int_{R_s} [ (min_{f!=s} alpha_f)^kpow - alpha_s^kpow ] alpha^exps   (Lemma lem:tl-special), via the ordered
    regions alpha_s <= alpha_f1 <= alpha_f2 <= ... (simplices with rational vertices)"""
    key = (exps, s, kpow)
    if key in _MIN2:
        return _MIN2[key]
    E = len(exps)
    others = [g for g in range(E) if g != s]
    tot = Fraction(0)
    for order in itertools.permutations(others):
        # region alpha_s <= alpha_{o0} <= alpha_{o1} <= ... ; vertices v_j = (1/(E-j)) * indicator of the last E-j
        seq = [s] + list(order)
        verts = []
        for j in range(E):
            v = [Fraction(0)] * E
            for g in seq[j:]:
                v[g] = Fraction(1, E - j)
            verts.append(v)
        # integrand: (alpha_{o0}^kpow - alpha_s^kpow) alpha^exps ; alpha = sum_j mu_j v_j
        for which, sign in ((seq[1], 1), (s, -1)):
            ex = list(exps); ex[which] += kpow
            tot += sign * _ordered_simplex_monomial(tuple(ex), verts)
    _MIN2[key] = tot
    return tot


def _ordered_simplex_monomial(exps, verts):
    """int over the simplex with vertices verts (E points in the alpha-simplex) of alpha^exps (Lebesgue d alpha)"""
    E = len(exps)
    # alpha_g = sum_j mu_j verts[j][g]; expand
    poly = {tuple([0] * E): Fraction(1)}
    for g in range(E):
        lin = {}
        for j in range(E):
            if verts[j][g] != 0:
                e = [0] * E; e[j] = 1
                lin[tuple(e)] = verts[j][g]
        poly = sp_mul(poly, sp_pow(lin, exps[g], E))
    # Jacobian: |det| of the vertex matrix relative to the standard simplex
    Mx = [[verts[j][g] for g in range(E - 1)] + [Fraction(1)] for j in range(E)]
    det = _det(Mx)
    tot = sum(v * dirichlet(e) for e, v in poly.items())
    return tot * abs(det)


def _det(M):
    M = [row[:] for row in M]
    n = len(M); d = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if M[r][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            M[i], M[p] = M[p], M[i]; d = -d
        d *= M[i][i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n):
                M[r][c] -= f * M[i][c]
    return d


class Simplex:
    """exact singular coefficients of a cycle with E edges, edge e weight e^{-t_e |k + qs[e]|^2}, t_e in [eps,oo),
    for a numerator polynomial in k with matrix (or scalar) coefficients.  In units hbar/(16 pi^2) x (stripped)."""

    def __init__(self, qs):
        self.E = len(qs)
        self.qs = [np.asarray(q, dtype=int) for q in qs]
        E = self.E
        # c(alpha) = - sum alpha_e qs[e]  (completing the square of sum t_e |k+q_e|^2)
        self.cpoly = []
        for mu in range(4):
            P = {}
            for e in range(E):
                if self.qs[e][mu] != 0:
                    ex = [0] * E; ex[e] = 1
                    P[tuple(ex)] = Fraction(-int(self.qs[e][mu]))
            self.cpoly.append(P)
        Phi = {}
        for e in range(E):
            ex = [0] * E; ex[e] = 1
            v = int(self.qs[e] @ self.qs[e])
            if v:
                Phi[tuple(ex)] = Fraction(v)
        cc = {}
        for mu in range(4):
            for k_, v in sp_mul(self.cpoly[mu], self.cpoly[mu]).items():
                cc[k_] = cc.get(k_, 0) + v
        for e_, v in cc.items():
            Phi[e_] = Phi.get(e_, 0) - v
        self.Phi = {e_: v for e_, v in Phi.items() if v != 0}
        self._cpow = {}
        self._phipow = {}
        self._w = {}

    def cpow(self, mu, n):
        if (mu, n) not in self._cpow:
            self._cpow[(mu, n)] = sp_pow(self.cpoly[mu], n, self.E)
        return self._cpow[(mu, n)]

    def phipow(self, n):
        if n not in self._phipow:
            self._phipow[n] = sp_pow({e: -v for e, v in self.Phi.items()}, n, self.E)   # (-Phi)^n
        return self._phipow[n]

    def weight(self, kexp, j, kind, m=0):
        """scalar weight of the k-monomial kexp of Delta^j N:  kind 'log': -int_Delta c^kexp (-Phi)^{m_j}/m_j!/(j!4^j);
        kind ('eps', m): (1/m) sum_e int_{R_e} alpha_e^m c^kexp (-Phi)^{m_j-m}/(m_j-m)!/(j!4^j);
        kind ('tot', m): (1/m) int_Delta (...) (total-time cutoff); kind ('spec', s, m): special edge s, eps^{-m}"""
        key = (kexp, j, kind, m)
        if key in self._w:
            return self._w[key]
        E = self.E
        mj = j + 2 - E
        if kind == 'log':
            p = 0 if mj < 0 else mj
            if mj < 0:
                self._w[key] = Fraction(0); return Fraction(0)
        else:
            if mj < m:
                self._w[key] = Fraction(0); return Fraction(0)
            p = mj - m
        poly = self.phipow(p)
        for mu in range(4):
            if kexp[mu]:
                poly = sp_mul(poly, self.cpow(mu, kexp[mu]))
        pref = Fraction(1, math.factorial(p) * math.factorial(j) * 4 ** j)
        if kind == 'log':
            val = -sum(v * dirichlet(e) for e, v in poly.items())
        elif kind == 'eps':
            val = Fraction(0)
            for e_, v in poly.items():
                for e in range(E):
                    ex = list(e_); ex[e] += m
                    val += v * region_monomial(tuple(ex), e)
            val /= m
        elif kind == 'tot':
            val = sum(v * dirichlet(e) for e, v in poly.items()) / m
        elif isinstance(kind, tuple) and kind[0] == 'spec':
            s = kind[1]
            val = Fraction(0)
            for e_, v in poly.items():
                val += v * specialregion_monomial(e_, s, m)
            val /= m
        else:
            raise ValueError(kind)
        val *= pref
        self._w[key] = val
        return val

    def singular(self, N, mmax=None, kinds=('log', 'eps')):
        """N: {kexp: coefficient} with coefficients numpy arrays (exact Gaussian integers) or Gaussian rationals.
        Returns {'log': value, m: value (eps^{-m}), ('tot', m): value} with exact Fraction-valued arrays/pairs."""
        degk = max((sum(e) for e in N), default=0)
        jmax = degk // 2
        laps = [N]
        for j in range(1, jmax + 1):
            laps.append(lap_poly(laps[-1]))
        if mmax is None:
            mmax = max(0, jmax + 2 - self.E)
        out = {}
        targets = []
        if 'log' in kinds:
            targets.append(('log', 0))
        if 'eps' in kinds:
            targets += [('eps', m) for m in range(1, mmax + 1)]
        if 'tot' in kinds:
            targets += [('tot', m) for m in range(1, mmax + 1)]
        for (kind, m) in targets:
            acc = None
            for j, Lj in enumerate(laps):
                for kexp, coef in Lj.items():
                    w = self.weight(kexp, j, kind, m)
                    if w == 0:
                        continue
                    acc = _acc(acc, coef, w)
            key = 'log' if kind == 'log' else (m if kind == 'eps' else ('tot', m))
            out[key] = _finish(acc)
        return out

    def special(self, N, s, mmax=None):
        degk = max((sum(e) for e in N), default=0)
        jmax = degk // 2
        laps = [N]
        for j in range(1, jmax + 1):
            laps.append(lap_poly(laps[-1]))
        if mmax is None:
            mmax = max(0, jmax + 2 - self.E)
        out = {}
        for m in range(1, mmax + 1):
            acc = None
            for j, Lj in enumerate(laps):
                for kexp, coef in Lj.items():
                    w = self.weight(kexp, j, ('spec', s), m)
                    if w:
                        acc = _acc(acc, coef, w)
            out[m] = _finish(acc)
        return out


def _acc(acc, coef, w):
    """accumulate coef * w with coef a Gaussian-integer numpy array (or complex scalar) and w a Fraction: kept as
    {denominator: integer object array} for exactness"""
    if acc is None:
        acc = {}
    d = w.denominator
    arr = np.asarray(coef)
    re = np.vectorize(lambda z: int(round(z.real)), otypes=[object])(arr) if arr.dtype != object else arr
    im = np.vectorize(lambda z: int(round(z.imag)), otypes=[object])(arr)
    assert np.all(np.abs(arr.real - np.rint(arr.real)) < 1e-9) and np.all(np.abs(arr.imag - np.rint(arr.imag)) < 1e-9)
    n = w.numerator
    if d in acc:
        acc[d] = (acc[d][0] + re * n, acc[d][1] + im * n)
    else:
        acc[d] = (re * n, im * n)
    return acc


def _finish(acc):
    """-> (re, im) object arrays of Fractions (or a pair of Fractions for scalars)"""
    if acc is None:
        return None
    re = None; im = None
    for d, (r, i) in acc.items():
        fr = np.vectorize(lambda x: Fraction(x, d), otypes=[object])(r)
        fi = np.vectorize(lambda x: Fraction(x, d), otypes=[object])(i)
        re = fr if re is None else re + fr
        im = fi if im is None else im + fi
    if np.ndim(re) == 0:
        return (Fraction(re[()]) if isinstance(re, np.ndarray) else re, Fraction(im[()]) if isinstance(im, np.ndarray) else im)
    return (re, im)


def qstr(z):
    """Gaussian rational (re, im) -> string"""
    re, im = Fraction(z[0]), Fraction(z[1])
    if im == 0:
        return str(re)
    if re == 0:
        return '%s*i' % im
    return '(%s+%s*i)' % (re, im)
