"""Engine A, round 2: exact one-loop core.
All arithmetic exact: every array entry is a Gaussian integer stored in complex128 (|.| < 2^52 is asserted on every
product) or as a pair of Python-int/Fraction values.  Cubic vertex functions are scaled by 2 and quartic ones by 8 so that
all vertex functions are Gaussian integers; the scale is divided out exactly at the end (Fraction arithmetic).
Conventions (classical/classical_spec.md): plane-wave legs, all momenta incoming; vertex function = coefficient of the
odd parameters in leg order; Hessian H_jl = d_j d_l I (left derivatives, l first) evaluated at the background, background
odd parameters factored to the left of the matrix entry; one loop W_1 = sum_n 1/(2n) str[(P (-1)^|.| H)^n]."""
import itertools
from fractions import Fraction
import numpy as np
from cliff import s, sb

SARR = np.array(s); SBARR = np.array(sb)
TYPES = ['c', 'A', 'As', 'l', 'ls', 'cs']
PARITY = {'c': 1, 'A': 0, 'As': 1, 'l': 1, 'ls': 0, 'cs': 0}
D = 3                                   # su(2), orthonormal basis for (-,-), [t_a,t_b] = eps_abc t_c
F = np.zeros((3, 3, 3))
for (a, b, c) in itertools.permutations(range(3)):
    F[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])
F = F.astype(complex)
C_ADJ = -2                              # tr_g(ad_X ad_Y)/(X,Y) for this basis
SHAPE = {'c': (D,), 'A': (10, D), 'As': (10, D), 'l': (16, D), 'ls': (16, D), 'cs': (D,)}
SIZE = {t: int(np.prod(SHAPE[t])) for t in TYPES}
OFF = {}
_o = 0
for _t in TYPES:
    OFF[_t] = _o; _o += SIZE[_t]
NTOT = _o
SDIAG = np.concatenate([np.full(SIZE[t], (-1) ** PARITY[t]) for t in TYPES]).astype(complex)
LIM = 2.0 ** 52


def mm(A, B):
    # integer-valued complex128 product; exact while every partial sum stays below 2^52 (asserted via the
    # entrywise bound |A|@|B|).  (numpy+Accelerate on macOS raises spurious FP flags in matmul: silenced.)
    with np.errstate(all='ignore'):
        bound = (np.abs(A.real) + np.abs(A.imag)) @ (np.abs(B.real) + np.abs(B.imag))
        if bound.size and not (bound.max() < LIM):
            raise OverflowError('exactness bound exceeded: %g' % bound.max())
        out = A @ B
    if not np.all(np.isfinite(out)):
        raise FloatingPointError('non-finite product')
    return out


def to_int(z):
    zr, zi = complex(z).real, complex(z).imag
    assert abs(zr - round(zr)) < 1e-6 and abs(zi - round(zi)) < 1e-6, z
    return (int(round(zr)), int(round(zi)))


class Q2:
    """exact Gaussian rational"""
    __slots__ = ['re', 'im']
    def __init__(self, re=0, im=0):
        self.re = Fraction(re); self.im = Fraction(im)
    @staticmethod
    def of(z):
        r, i = to_int(z); return Q2(r, i)
    def __add__(s_, o):
        o = o if isinstance(o, Q2) else Q2(o); return Q2(s_.re + o.re, s_.im + o.im)
    __radd__ = __add__
    def __neg__(s_):
        return Q2(-s_.re, -s_.im)
    def __sub__(s_, o):
        return s_ + (-(o if isinstance(o, Q2) else Q2(o)))
    def __mul__(s_, o):
        o = o if isinstance(o, Q2) else Q2(o)
        return Q2(s_.re * o.re - s_.im * o.im, s_.re * o.im + s_.im * o.re)
    __rmul__ = __mul__
    def __truediv__(s_, o):
        o = o if isinstance(o, Q2) else Q2(o)
        n = o.re ** 2 + o.im ** 2
        return Q2((s_.re * o.re + s_.im * o.im) / n, (s_.im * o.re - s_.re * o.im) / n)
    def iszero(s_):
        return s_.re == 0 and s_.im == 0
    def __eq__(s_, o):
        o = o if isinstance(o, Q2) else Q2(o); return s_.re == o.re and s_.im == o.im
    def __repr__(s_):
        if s_.im == 0: return str(s_.re)
        if s_.re == 0: return '%s*i' % s_.im
        return '(%s+%s*i)' % (s_.re, s_.im)


def mom10(k):
    v = np.zeros(10, dtype=complex); v[:4] = k[:4]; return v

# ---------------------------------------------------------------- canonical vertex functions (scaled)
CUBIC = {'AAA': ['A', 'A', 'A'], 'llA': ['l', 'l', 'A'], 'AcAs': ['A', 'c', 'As'], 'cllS': ['c', 'l', 'ls'],
         'cccS': ['c', 'c', 'cs']}
QUARTIC = {'AAAA': ['A', 'A', 'A', 'A']}
CANON = dict(CUBIC); CANON.update(QUARTIC)
SCALE = {'AAA': 2, 'llA': 2, 'AcAs': 2, 'cllS': 2, 'cccS': 2, 'AAAA': 8}


def canon(vt, pols, moms):
    """scaled canonical vertex function; pols[j] has a leading batch axis"""
    f = F
    if vt == 'AAA':                       # 2 * (-1/2) sum_perm (d_M A_N, [A_M, A_N]),  d -> i k
        out = 0
        for perm in itertools.permutations(range(3)):
            a, b, c = perm
            t = 1j * np.einsum('M,iNe,jMp,kNq,pqe->ijk', mom10(moms[a]), pols[a], pols[b], pols[c], f, optimize=True)
            out = out + (-1) * np.transpose(t, np.argsort(perm))
        return out
    if vt == 'AAAA':                      # 8 * (-1/8) sum_perm ([A_M,A_N],[A_M,A_N])
        out = 0
        for perm in itertools.permutations(range(4)):
            a, b, c, e = perm
            t = np.einsum('iMx,jNy,kMz,lNw,xyg,zwg->ijkl', pols[a], pols[b], pols[c], pols[e], f, f, optimize=True)
            out = out + (-1) * np.transpose(t, np.argsort(perm))
        return out
    if vt == 'llA':                       # 2 * (u1, rho([a,u2]))
        return 2 * np.einsum('iAe,MAB,kMp,jBq,pqe->ijk', pols[0], SBARR, pols[2], pols[1], f, optimize=True)
    if vt == 'AcAs':                      # 2 * ([a,gamma],alpha)
        return 2 * np.einsum('iMp,jq,kMe,pqe->ijk', pols[0], pols[1], pols[2], f, optimize=True)
    if vt == 'cllS':                      # 2 * ([gamma,u],v)
        return 2 * np.einsum('ip,jAq,kAe,pqe->ijk', pols[0], pols[1], pols[2], f, optimize=True)
    if vt == 'cccS':                      # 2 * ([gamma1,gamma2],w)
        return 2 * np.einsum('ip,jq,ke,pqe->ijk', pols[0], pols[1], pols[2], f, optimize=True)
    raise ValueError(vt)


def V(legs, only=None):
    """scaled multilinear coefficient (odd parameters in the order of `legs`); legs = [(type, pol_batch, mom)]"""
    types = [L[0] for L in legs]
    out = 0
    for vt, can in CANON.items():
        if only is not None and vt not in only:
            continue
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
        val = np.transpose(canon(vt, pols, moms), slot)
        odd_req = [slot[i] for i in range(len(legs)) if PARITY[types[i]]]
        sgn = 1
        for x in range(len(odd_req)):
            for y in range(x + 1, len(odd_req)):
                if odd_req[x] > odd_req[y]:
                    sgn = -sgn
        out = out + sgn * val
    return out


def ident(t):
    sz = SIZE[t]
    return np.eye(sz, dtype=complex).reshape((sz,) + SHAPE[t])


FLUCT_ALL = ['c', 'A', 'As', 'l', 'ls', 'cs']


def hessian(bgs, kj, kl, fluct=FLUCT_ALL):
    """H_jl = d_j d_l I at the background(s) bgs=[(type,pol,mom)] (scaled), j-leg incoming momentum kj, l-leg kl;
    the odd background parameters stand to the left in bgs order."""
    H = np.zeros((NTOT, NTOT), dtype=complex)
    for tj in fluct:
        for tl in fluct:
            legs = [(tj, ident(tj), kj), (tl, ident(tl), kl)] + [(t, pol[None], m) for (t, pol, m) in bgs]
            v = V(legs)
            if isinstance(v, int):
                continue
            v = v.reshape(SIZE[tj], SIZE[tl])
            sg = (-1) ** (PARITY[tj] * PARITY[tl])
            H[OFF[tj]:OFF[tj] + SIZE[tj], OFF[tl]:OFF[tl] + SIZE[tl]] += sg * v
    return H

# ---------------------------------------------------------------- kernels
def comp_index(t, which=None):
    """component indices of type t; for A/As: which='V' (M<4) or 'P' (M>=4)"""
    if t in ('A', 'As') and which is not None:
        Ms = range(4) if which == 'V' else range(4, 10)
        return [OFF[t] + M * D + a for M in Ms for a in range(D)]
    return list(range(OFF[t], OFF[t] + SIZE[t]))


def Pnum(k, line):
    """numerator of P(eps,infty)(k) restricted to one line type: 'AV','AP','gh','l'"""
    N = np.zeros((NTOT, NTOT), dtype=complex)
    if line in ('AV', 'AP'):
        idx = comp_index('A', line[1])
        N[idx, idx] = 2
    elif line == 'gh':
        for mu in range(4):
            for a in range(D):
                N[OFF['c'] + a, OFF['As'] + mu * D + a] = -1j * k[mu]
                N[OFF['As'] + mu * D + a, OFF['c'] + a] = -1j * k[mu]
    elif line == 'l':
        sl = np.einsum('m,mab->ab', np.asarray(k[:4], dtype=complex), SARR[:4])
        N[OFF['l']:OFF['l'] + SIZE['l'], OFF['l']:OFF['l'] + SIZE['l']] = 1j * np.kron(sl, np.eye(D))
    else:
        raise ValueError(line)
    return N


def Khat(line):
    """momentum-independent factor of K_t(k) = e^{-t k^2} Khat: -sum (x_i (x) xi_i + xi_i (x) x_i) (x) (-,-)^{-1}
    pairs (A,A^*), (lambda^*,lambda), (c^*,c); restricted to one pair type 'sAV','sAP','sl','sc'"""
    N = np.zeros((NTOT, NTOT), dtype=complex)
    if line in ('sAV', 'sAP'):
        i1 = comp_index('A', line[2]); i2 = comp_index('As', line[2])
    elif line == 'sl':
        i1 = comp_index('ls'); i2 = comp_index('l')
    elif line == 'sc':
        i1 = comp_index('cs'); i2 = comp_index('c')
    else:
        raise ValueError(line)
    for x, y in zip(i1, i2):
        N[x, y] = -1; N[y, x] = -1
    return N


PLINES = ['AV', 'AP', 'gh', 'l']
SLINES = ['sAV', 'sAP', 'sl', 'sc']


def strace(X):
    return np.sum(SDIAG * np.diag(X))
