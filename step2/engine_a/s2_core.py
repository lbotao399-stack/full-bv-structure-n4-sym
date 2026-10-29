"""Step 2, engine A: exact one-loop core, generalised from engine_a/r2_core.py to a non-orthonormal colour basis.
Gauge algebra selected by the environment variable S2_ALG:
  'su2' : orthonormal basis of su(2) for (-,-), [t_a,t_b] = eps_abc t_c; tr_g(ad_X ad_Y)/(X,Y) = -2.
  'gl3' : basis E_ij of gl(3), (X,Y) = tr(XY) (unimodular: g^{-1} = g); sl(3) = complexified su(3) is the traceless
          part, the centre has vanishing brackets and decouples from every vertex; for X,Y in sl(3)
          tr_{sl(3)}(ad_X ad_Y)/(X,Y) = 6.
Fields carry upper colour components; F[a,b,c] = ([t_a,t_b],t_c) (lowered, totally antisymmetric); every component of
P(eps,infty) and K_t carries (-,-)^{-1} = GINV.  All arithmetic exact (Gaussian integers in complex128 with the 2^52
bound asserted, as in r2_core)."""
import os, sys, itertools
from fractions import Fraction
import numpy as np
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
from cliff import s, sb
from r2_core import Q2, to_int, mom10

ALG = os.environ.get('S2_ALG', 'su2')
SARR = np.array(s); SBARR = np.array(sb)
TYPES = ['c', 'A', 'As', 'l', 'ls', 'cs']
PARITY = {'c': 1, 'A': 0, 'As': 1, 'l': 1, 'ls': 0, 'cs': 0}
if ALG == 'su2':
    D = 3
    F = np.zeros((3, 3, 3))
    for (a, b, c) in itertools.permutations(range(3)):
        F[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])
    G = np.eye(3)
    C_ADJ = -2
elif ALG == 'gl3':
    D = 9
    def E(i, j):
        m = np.zeros((3, 3)); m[i, j] = 1; return m
    basis = [E(i, j) for i in range(3) for j in range(3)]
    G = np.array([[np.trace(x @ y) for y in basis] for x in basis])
    F = np.array([[[np.trace((x @ y - y @ x) @ z) for z in basis] for y in basis] for x in basis])
    C_ADJ = 6
else:
    raise ValueError(ALG)
F = np.rint(F).astype(complex)
GINV = np.rint(np.linalg.inv(G)).astype(complex)
assert np.allclose(np.linalg.inv(G), GINV.real)
G = np.rint(G).astype(complex)
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
    with np.errstate(all='ignore'):
        bound = (np.abs(A.real) + np.abs(A.imag)) @ (np.abs(B.real) + np.abs(B.imag))
        if bound.size and not (bound.max() < LIM):
            raise OverflowError('exactness bound exceeded: %g' % bound.max())
        out = A @ B
    if not np.all(np.isfinite(out)):
        raise FloatingPointError('non-finite product')
    return out

# ---------------------------------------------------------------- canonical vertex functions (scaled), as r2_core
CUBIC = {'AAA': ['A', 'A', 'A'], 'llA': ['l', 'l', 'A'], 'AcAs': ['A', 'c', 'As'], 'cllS': ['c', 'l', 'ls'],
         'cccS': ['c', 'c', 'cs']}
QUARTIC = {'AAAA': ['A', 'A', 'A', 'A']}
CANON = dict(CUBIC); CANON.update(QUARTIC)
SCALE = {'AAA': 2, 'llA': 2, 'AcAs': 2, 'cllS': 2, 'cccS': 2, 'AAAA': 8}


def canon(vt, pols, moms):
    f = F
    if vt == 'AAA':
        out = 0
        for perm in itertools.permutations(range(3)):
            a, b, c = perm
            t = 1j * np.einsum('M,iNe,jMp,kNq,pqe->ijk', mom10(moms[a]), pols[a], pols[b], pols[c], f, optimize=True)
            out = out + (-1) * np.transpose(t, np.argsort(perm))
        return out
    if vt == 'AAAA':
        out = 0
        for perm in itertools.permutations(range(4)):
            a, b, c, e = perm
            # ([A_M,A_N],[A_M,A_N]) = f_{xy}^g f_{zw}^h g_{gh} ... = F[x,y,g'] GINV[g',h'] F[z,w,h']
            t = np.einsum('iMx,jNy,kMz,lNw,xyg,gh,zwh->ijkl', pols[a], pols[b], pols[c], pols[e], f, GINV, f,
                          optimize=True)
            out = out + (-1) * np.transpose(t, np.argsort(perm))
        return out
    if vt == 'llA':
        return 2 * np.einsum('iAe,MAB,kMp,jBq,pqe->ijk', pols[0], SBARR, pols[2], pols[1], f, optimize=True)
    if vt == 'AcAs':
        return 2 * np.einsum('iMp,jq,kMe,pqe->ijk', pols[0], pols[1], pols[2], f, optimize=True)
    if vt == 'cllS':
        return 2 * np.einsum('ip,jAq,kAe,pqe->ijk', pols[0], pols[1], pols[2], f, optimize=True)
    if vt == 'cccS':
        return 2 * np.einsum('ip,jq,ke,pqe->ijk', pols[0], pols[1], pols[2], f, optimize=True)
    raise ValueError(vt)


def V(legs, only=None):
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

# ---------------------------------------------------------------- kernels with (-,-)^{-1} = GINV
def comp_index(t, which=None):
    if t in ('A', 'As') and which is not None:
        Ms = range(4) if which == 'V' else range(4, 10)
        return [OFF[t] + M * D + a for M in Ms for a in range(D)]
    return list(range(OFF[t], OFF[t] + SIZE[t]))


def Pnum(k, line):
    """numerator of P(eps,infty)(k) (first factor at momentum k) for one line type: 'AV','AP','gh','l'"""
    N = np.zeros((NTOT, NTOT), dtype=complex)
    if line in ('AV', 'AP'):
        Ms = range(4) if line == 'AV' else range(4, 10)
        for M in Ms:
            i0 = OFF['A'] + M * D
            N[i0:i0 + D, i0:i0 + D] = 2 * GINV
    elif line == 'gh':
        for mu in range(4):
            i0 = OFF['As'] + mu * D
            N[OFF['c']:OFF['c'] + D, i0:i0 + D] = -1j * k[mu] * GINV
            N[i0:i0 + D, OFF['c']:OFF['c'] + D] = -1j * k[mu] * GINV
    elif line == 'l':
        sl = np.einsum('m,mab->ab', np.asarray(k[:4], dtype=complex), SARR[:4])
        N[OFF['l']:OFF['l'] + SIZE['l'], OFF['l']:OFF['l'] + SIZE['l']] = 1j * np.kron(sl, GINV)
    else:
        raise ValueError(line)
    return N


def Khat(line):
    """momentum-independent factor of K_t(k) = e^{-t k^2} Khat, -sum(x_i (x) xi_i + xi_i (x) x_i) (x) (-,-)^{-1},
    restricted to one pair type 'sAV','sAP','sl','sc'"""
    N = np.zeros((NTOT, NTOT), dtype=complex)
    if line in ('sAV', 'sAP'):
        Ms = range(4) if line == 'sAV' else range(4, 10)
        blocks = [(OFF['A'] + M * D, OFF['As'] + M * D) for M in Ms]
    elif line == 'sl':
        blocks = [(OFF['ls'] + A * D, OFF['l'] + A * D) for A in range(16)]
    elif line == 'sc':
        blocks = [(OFF['cs'], OFF['c'])]
    else:
        raise ValueError(line)
    for (x0, y0) in blocks:
        N[x0:x0 + D, y0:y0 + D] = -GINV
        N[y0:y0 + D, x0:x0 + D] = -GINV
    return N


PLINES = ['AV', 'AP', 'gh', 'l']
SLINES = ['sAV', 'sAP', 'sl', 'sc']


def linemat(L, k):
    if L in PLINES:
        return Pnum(k, L)
    return Khat(L)


def strace(X):
    return np.sum(SDIAG * np.diag(X))


def traceless(x):
    """project the colour (last) axis onto sl(3) for gl3 (integer-preserving: X^{00} := -X^{11}-X^{22});
    identity for su2"""
    if ALG != 'gl3':
        return x
    x = np.array(x, dtype=complex)
    x[..., 0] = -x[..., 4] - x[..., 8]
    return x


def rpol(rng, t, which=None, lo=-1, hi=2):
    while True:
        x = rng.integers(lo, hi, size=SHAPE[t]).astype(complex)
        if which == 'V': x[4:] = 0
        if which == 'P': x[:4] = 0
        x = traceless(x)
        if np.abs(x).sum() > 0:
            return x
