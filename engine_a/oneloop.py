"""One-loop engine A: W_1 = sum_n 1/(2n) str[(P S d2I)^n] (verified in verify_oneloop0d.py), momentum space,
per-edge Schwinger parameters t_e in [eps, infty), exact Gaussian loop integral, exact eps-asymptotics."""
import itertools
import warnings
import numpy as np
from cliff import s, sb
from vertices import structure_constants, c_adj
warnings.filterwarnings('ignore')
SBARR = np.array(sb); SARR = np.array(s)
TYPES = ['A', 'c', 'As', 'l']
CANON = {'AAA': ['A', 'A', 'A'], 'AAAA': ['A', 'A', 'A', 'A'], 'llA': ['l', 'l', 'A'], 'AcAs': ['A', 'c', 'As']}
PARITY = {'A': 0, 'c': 1, 'As': 1, 'l': 1, 'ls': 0, 'cs': 0}

def mom10(k):
    v = np.zeros(10, dtype=complex); v[:4] = k[:4]; return v

class Engine:
    def __init__(self, gname='su2'):
        self.f = structure_constants(gname).astype(complex)
        self.d = d = self.f.shape[0]
        self.cadj = c_adj(self.f.real)
        self.shape = {'A': (10, d), 'c': (d,), 'As': (10, d), 'l': (16, d)}
        self.size = {t: int(np.prod(sh)) for t, sh in self.shape.items()}
        self.off = {}; o = 0
        for t in TYPES:
            self.off[t] = o; o += self.size[t]
        self.n = o
        self.S = np.concatenate([np.full(self.size[t], (-1) ** PARITY[t]) for t in TYPES]).astype(complex)

    def ident(self, t):
        sz = self.size[t]
        return np.eye(sz, dtype=complex).reshape((sz,) + self.shape[t])

    # ---------- canonical vertex functions with batch axes (pols have a leading batch axis)
    def canon(self, vt, pols, moms):
        f = self.f
        if vt == 'AAA':
            out = 0
            for perm in itertools.permutations(range(3)):
                a, b, c = perm
                t = 1j * np.einsum('M,iNe,jMp,kNq,pqe->ijk', mom10(moms[a]), pols[a], pols[b], pols[c], f, optimize=True)
                inv = np.argsort(perm)
                out = out + (-0.5) * np.transpose(t, inv)
            return out
        if vt == 'AAAA':
            out = 0
            for perm in itertools.permutations(range(4)):
                a, b, c, e = perm
                t = np.einsum('iMx,jNy,kMz,lNw,xyg,zwg->ijkl', pols[a], pols[b], pols[c], pols[e], f, f, optimize=True)
                out = out + (-1 / 8) * np.transpose(t, np.argsort(perm))
            return out
        if vt == 'llA':
            return np.einsum('iAe,MAB,kMp,jBq,pqe->ijk', pols[0], SBARR, pols[2], pols[1], f, optimize=True)
        if vt == 'AcAs':
            return np.einsum('iMp,jq,kMe,pqe->ijk', pols[0], pols[1], pols[2], f, optimize=True)
        raise ValueError(vt)

    def V(self, legs):
        """multilinear coefficient with odd parameters in the order of `legs`; legs = [(type, pol_batch, mom)].
        Returns array with one axis per leg (batch axes)."""
        types = [L[0] for L in legs]
        out = 0
        for vt, can in CANON.items():
            if sorted(can) != sorted(types):
                continue
            # first-fit assignment of requested legs to canonical slots
            slot = [None] * len(legs); used = set()
            for i, t in enumerate(types):
                for j, ct in enumerate(can):
                    if ct == t and j not in used:
                        slot[i] = j; used.add(j); break
            pols = [None] * len(can); moms = [None] * len(can)
            for i, j in enumerate(slot):
                pols[j] = legs[i][1]; moms[j] = legs[i][2]
            val = self.canon(vt, pols, moms)
            # back to requested leg order
            val = np.transpose(val, slot)
            # Koszul sign: odd legs in canonical order vs requested order
            odd_req = [slot[i] for i in range(len(legs)) if PARITY[types[i]]]
            sgn = 1
            for x in range(len(odd_req)):
                for y in range(x + 1, len(odd_req)):
                    if odd_req[x] > odd_req[y]:
                        sgn = -sgn
            out = out + sgn * val
        return out

    def hessian(self, bgs, kj, kl):
        """H_{jl} = d_j d_l I at background(s) bgs=[(type,pol,mom)], j-leg momentum kj, l-leg momentum kl,
        with the background Grassmann parameters factored out to the left in bgs order."""
        H = np.zeros((self.n, self.n), dtype=complex)
        for tj in TYPES:
            for tl in TYPES:
                legs = [(tj, self.ident(tj), kj), (tl, self.ident(tl), kl)] + [(t, pol[None], m) for (t, pol, m) in bgs]
                v = self.V(legs)
                if isinstance(v, int):
                    continue
                v = v.reshape(self.size[tj], self.size[tl])
                # term = y_j y_l r_b... V ;  d_j d_l (y_j y_l R V) = (-1)^{|j||l|} R V
                sg = (-1) ** (PARITY[tj] * PARITY[tl])
                # moving R (background params, product parity pb) left past nothing: R stands after y_j y_l
                H[self.off[tj]:self.off[tj] + self.size[tj], self.off[tl]:self.off[tl] + self.size[tl]] += sg * v
        return H

    def Nprop(self, k):
        d = self.d; N = np.zeros((self.n, self.n), dtype=complex)
        oA, oc, oS, ol = self.off['A'], self.off['c'], self.off['As'], self.off['l']
        N[oA:oA + 10 * d, oA:oA + 10 * d] = 2 * np.eye(10 * d)
        for mu in range(4):
            for a in range(d):
                N[oc + a, oS + mu * d + a] = -1j * k[mu]
                N[oS + mu * d + a, oc + a] = -1j * k[mu]
        sl = np.einsum('m,mab->ab', k[:4], SARR[:4])
        N[ol:ol + 16 * d, ol:ol + 16 * d] = 1j * np.kron(sl, np.eye(d))
        return N

    def str_(self, X):
        return np.sum(self.S * np.diag(X))
