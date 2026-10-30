"""Vertex functions of the dimensionally reduced 10d Yang-Mills theory, obtained by polarisation of the
cubic and quartic terms of -S_E (Euclidean, coupling g = 1), with explicit su(N) structure constants.
-S_E,3 = - f^{abc} (d_M A^a_N) A^b_M A^c_N,  -S_E,4 = -(1/4) f^{abe} f^{cde} A^a_M A^b_N A^c_M A^d_N,  d -> i k, k_I = 0.
"""
import itertools

import numpy as np


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
    (A1, k1), (A2, k2), (A3, k3) = fields
    dA1 = 1j * np.einsum('M,aN->aMN', k1, A1)
    return -np.einsum('abc,aMN,bM,cN->', f, dA1, A2, A3)


def V3(fields, f):
    return sum(cubic_density([fields[i] for i in perm], f) for perm in itertools.permutations(range(3)))


def quartic_density(As, f):
    A1, A2, A3, A4 = As
    return -0.25 * np.einsum('abe,cde,aM,bN,cM,dN->', f, f, A1, A2, A3, A4)


def V4(As, f):
    return sum(quartic_density([As[i] for i in perm], f) for perm in itertools.permutations(range(4)))


def unit_field(dimG, a, M, coef=1.0):
    A = np.zeros((dimG, 10), complex); A[a, M] = coef
    return A


def scalar_field(dimG, a, w):
    """colour a, flavour six-vector w placed in the components 4..9."""
    A = np.zeros((dimG, 10), complex); A[a, 4:] = w
    return A


def vec(k4):
    k = np.zeros(10); k[:4] = k4
    return k


# ---- fast tensor forms (validated against the polarisation functions above in fourier_opct.py) ----
def V3_tensor(ks, f):
    """V[x0,x1,x2], x = (colour, M) flattened, for incoming momenta ks = (k0,k1,k2); same definition as V3."""
    dimG = f.shape[0]
    eye = np.eye(10)
    out = 0
    for perm in itertools.permutations(range(3)):
        k = ks[perm[0]]
        T = -1j * np.einsum('abc,M,NP->aNbMcP', f, k, eye).reshape(dimG * 10, dimG * 10, dimG * 10)
        inv = [perm.index(i) for i in range(3)]
        out = out + np.transpose(T, axes=inv)
    return out


def V4_tensor(f, comps):
    """V[x0..x3] restricted to the index components `comps` (list of M), x = (colour, M in comps) flattened."""
    dimG = f.shape[0]
    n = len(comps)
    C1 = np.einsum('abe,cde->abcd', f, f)
    d = np.eye(n)
    T = (-0.25 * np.einsum('abcd,MP,NQ->aMbNcPdQ', C1, d, d)).reshape(dimG * n, dimG * n, dimG * n, dimG * n)
    out = 0
    for perm in itertools.permutations(range(4)):
        inv = [perm.index(i) for i in range(4)]
        out = out + np.transpose(T, axes=inv)
    return out
