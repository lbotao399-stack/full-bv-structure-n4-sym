"""Colour data and momentum-space vertex functions of I (engine A).
Vertex function = coefficient of r_1 r_2 ... (odd legs in leg order) in the multilinear part of I evaluated on
sum_j r_j u_j e^{i p_j.x} (even legs: nilpotent even parameters), (2 pi)^4 delta stripped.
Closed forms below are checked against a brute-force Grassmann evaluation of [ESW (31)] built directly from
F_A = dA + 1/2[A^A] (Hodge pairing), (lambda, dslash_A lambda), (d_A c, A^*), ([c,lambda],lambda^*), 1/2([c,c],c^*)."""
import itertools
import warnings
import numpy as np
from cliff import s, sb
warnings.filterwarnings('ignore')

def structure_constants(name):
    """orthonormal basis t_a w.r.t. (-,-) (positive definite on the compact form), [t_a,t_b] = f_abc t_c"""
    if name == 'su2':
        f = np.zeros((3, 3, 3))
        for (a, b, c) in itertools.permutations(range(3)):
            f[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])
        return f
    if name == 'su3':
        f = np.zeros((8, 8, 8))
        vals = {(0, 1, 2): 1, (0, 3, 6): .5, (0, 4, 5): -.5, (1, 3, 5): .5, (1, 4, 6): .5, (2, 3, 4): .5,
                (2, 5, 6): -.5, (3, 4, 7): np.sqrt(3) / 2, (5, 6, 7): np.sqrt(3) / 2}
        for (a, b, c), v in vals.items():
            for perm in itertools.permutations(range(3)):
                idx = tuple((a, b, c)[i] for i in perm)
                f[idx] = v * np.linalg.det(np.eye(3)[list(perm)])
        return f
    raise ValueError(name)

def c_adj(f):
    """scalar by which X -> sum_a [t_a,[t^a,X]] acts"""
    d = f.shape[0]
    M = np.einsum('abc,acd->bd', f, f)   # [t_a,[t_a,t_b]] = f_abc f_acd t_d
    val = M[0, 0]
    assert np.allclose(M, val * np.eye(d))
    return val

class Vtx:
    def __init__(self, f):
        self.f = f
        self.d = f.shape[0]
    def br(self, x, y):           # [x,y] for colour vectors (last axis colour)
        return np.einsum('abc,...a,...b->...c', self.f, x, y)
    # ---- cubic
    def AAA(self, a1, p1, a2, p2, a3, p3):
        legs = [(a1, p1), (a2, p2), (a3, p3)]
        tot = 0
        for (x, px), (y, py), (z, pz) in itertools.permutations(legs):
            # -1/2 sum_{MN} (d_M A_N, [A_M, A_N]),  d_M only along V
            dA = np.zeros((10, 10, self.d), dtype=complex)
            for M in range(4):
                dA[M] = 1j * px[M] * x
            br = np.einsum('abc,Ma,Nb->MNc', self.f, y, z)
            tot += -0.5 * np.einsum('MNc,MNc->', dA, br)
        return tot
    def AAAA(self, a1, a2, a3, a4):
        tot = 0
        for (x, y, z, w) in itertools.permutations([a1, a2, a3, a4]):
            b1 = np.einsum('abc,Ma,Nb->MNc', self.f, x, y)
            b2 = np.einsum('abc,Ma,Nb->MNc', self.f, z, w)
            tot += -1 / 8 * np.einsum('MNc,MNc->', b1, b2)
        return tot
    def llA(self, u1, u2, a):       # coefficient of r1 r2
        # (u1, rho([a,u2])) = sum_M u1^T sb_M [a_M, u2]
        b = np.einsum('abc,Ma,ib->Mic', self.f, a, u2)
        return sum(np.einsum('ic,ij,jc->', u1, sb[M], b[M]) for M in range(10))
    def AcAs(self, a, g, al):       # coefficient of r_c r_*
        return np.einsum('abc,Ma,b,Mc->', self.f, a, g, al)
    def cllS(self, g, u, v):        # coefficient of r_c r_lambda
        return np.einsum('abc,a,ib,ic->', self.f, g, u, v)
    def cccS(self, g1, g2, w):      # coefficient of r1 r2
        return np.einsum('abc,a,b,c->', self.f, g1, g2, w)
