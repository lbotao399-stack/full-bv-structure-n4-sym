"""Refuter #2 (resumed run), own realisation, written from Ch.1 / Ch.3 only.

Realisation:
  V_10 = C^10, V = directions 0..3, W_perp = directions 4..9.
  Sigma, Sigma^* each C^16, pairing (xi_s, psi_t) = delta_st.
  rho(v): Sigma -> Sigma^* has matrix A(v) = sum v_M a_M, rho(v): Sigma^* -> Sigma has B(v) = sum v_M b_M,
  a_M, b_M complex symmetric 16x16 with a_M b_N + a_N b_M = 2 delta_MN = b_M a_N + b_N a_M.
  (symmetric a_M, b_M  <=>  (rho(v)Q1,Q2) = (Q1, rho(v)Q2) with C = identity)
  g = su(2): [e_a,e_b] = eps_abc e_c, (e_a,e_b) = delta_ab, tr(ad_X ad_Y)/(X,Y) = -2.

Grassmann bookkeeping: an element is a dict {tuple_of_generators_in_order: complex}; products
concatenate tuples; 'canon' sorts a tuple and returns the Koszul sign.
Floating complex arithmetic is used only with dyadic rationals of small size (exact in binary).
"""
import itertools
import numpy as np

# ---------------------------------------------------------------- Clifford
s0 = np.eye(2)
s1 = np.array([[0., 1.], [1., 0.]])
s3 = np.array([[1., 0.], [0., -1.]])
ep = np.array([[0., 1.], [-1., 0.]])          # i sigma_2 (real antisymmetric)


def _kron(*ms):
    out = np.array([[1.0]])
    for m in ms:
        out = np.kron(out, m)
    return out


def _real_sym_cl9():
    cands = []
    for combo in itertools.product([s0, s1, ep, s3], repeat=4):
        M = _kron(*combo)
        if np.allclose(M, M.T) and np.allclose(M @ M, np.eye(16)) and not np.allclose(M, np.eye(16)):
            cands.append(M)
    # backtracking search for 9 mutually anticommuting
    best = []

    def bt(chosen, start):
        nonlocal best
        if len(chosen) == 9:
            best = list(chosen)
            return True
        for i in range(start, len(cands)):
            M = cands[i]
            if all(np.allclose(M @ X + X @ M, 0) for X in chosen):
                chosen.append(M)
                if bt(chosen, i + 1):
                    return True
                chosen.pop()
        return False

    bt([], 0)
    assert len(best) == 9
    return best


_g9 = _real_sym_cl9()
a_mat = np.zeros((10, 16, 16), dtype=complex)
b_mat = np.zeros((10, 16, 16), dtype=complex)
for M in range(9):
    a_mat[M] = _g9[M]
    b_mat[M] = _g9[M]
a_mat[9] = 1j * np.eye(16)
b_mat[9] = -1j * np.eye(16)

for M in range(10):
    assert np.allclose(a_mat[M], a_mat[M].T) and np.allclose(b_mat[M], b_mat[M].T)
    for N in range(10):
        assert np.allclose(a_mat[M] @ b_mat[N] + a_mat[N] @ b_mat[M], 2 * (M == N) * np.eye(16))
        assert np.allclose(b_mat[M] @ a_mat[N] + b_mat[N] @ a_mat[M], 2 * (M == N) * np.eye(16))

VDIR = [0, 1, 2, 3]
WDIR = [4, 5, 6, 7, 8, 9]


def Amat(v):
    """rho(v): Sigma -> Sigma^*  (v in C^10)."""
    return np.tensordot(np.asarray(v, dtype=complex), a_mat, axes=(0, 0))


def Bmat(v):
    """rho(v): Sigma^* -> Sigma."""
    return np.tensordot(np.asarray(v, dtype=complex), b_mat, axes=(0, 0))


def vol4_on_Sigma():
    """rho(v1 v2 v3 v4) on Sigma = B(e0) A(e1) B(e2) A(e3)."""
    e = np.eye(10)
    return Bmat(e[0]) @ Amat(e[1]) @ Bmat(e[2]) @ Amat(e[3])


# ---------------------------------------------------------------- su(2)
EPS = np.zeros((3, 3, 3))
for (i, j, k), sg in {(0, 1, 2): 1, (1, 2, 0): 1, (2, 0, 1): 1, (0, 2, 1): -1, (2, 1, 0): -1, (1, 0, 2): -1}.items():
    EPS[i, j, k] = sg
KILLING_RATIO = -2.0     # tr(ad_X ad_Y)/(X,Y) for this basis/pairing
assert np.allclose(np.einsum('abc,dcb->ad', EPS, EPS), KILLING_RATIO * np.eye(3))


def gbr(X, Y):
    """[X,Y] for colour vectors X, Y (last axis colour)."""
    return np.einsum('abc,...a,...b->...c', EPS, X, Y)


# ---------------------------------------------------------------- Grassmann
def canon(gens):
    """Sort a tuple of odd generators; return (sign, sorted tuple) or (0, None) if repeated."""
    g = list(gens)
    if len(set(g)) < len(g):
        return 0, None
    sign = 1
    # bubble sort counting transpositions
    n = len(g)
    for i in range(n):
        for j in range(n - 1 - i):
            if g[j] > g[j + 1]:
                g[j], g[j + 1] = g[j + 1], g[j]
                sign = -sign
    return sign, tuple(g)


def gadd(*elems):
    out = {}
    for e in elems:
        for k, v in e.items():
            out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if abs(v) > 1e-13}


def gscale(e, c):
    return {k: v * c for k, v in e.items()}


def gmul(e1, e2):
    out = {}
    for k1, v1 in e1.items():
        for k2, v2 in e2.items():
            sg, key = canon(k1 + k2)
            if sg == 0:
                continue
            out[key] = out.get(key, 0) + sg * v1 * v2
    return out


def gcoef(e, gens):
    """Coefficient of the ordered product of generators 'gens' in e."""
    sg, key = canon(tuple(gens))
    return sg * e.get(key, 0)


# ---------------------------------------------------------------- plane-wave fields
# A field configuration: dict slot -> list of terms; a term is (gens_tuple, vec, p, legs_frozenset)
#   gens_tuple: odd generators in order (Grassmann factor placed to the LEFT of the fibre vector)
#   vec: numpy array (A, Ast: (10,3); c, cst: (3,); lam, lamst: (16,3))
#   p: momentum (4,) numpy array (V components)
#   legs: frozenset of leg labels carried by this term
SLOTS = ['A', 'Ast', 'c', 'cst', 'lam', 'lamst']
SHAPE = {'A': (10, 3), 'Ast': (10, 3), 'c': (3,), 'cst': (3,), 'lam': (16, 3), 'lamst': (16, 3)}
PARITY = {'A': 0, 'Ast': 1, 'c': 1, 'cst': 0, 'lam': 1, 'lamst': 0}


def empty_fields():
    return {s: [] for s in SLOTS}


def p10(p):
    """embed a V-momentum (4,) into C^10."""
    v = np.zeros(10, dtype=complex)
    v[:4] = p
    return v


# derivative operations on a single term (return new vec)
def d_scalar(vec, p):
    """d on a g-valued scalar (shape (3,)) -> one-form (10,3) with V components i p."""
    out = np.zeros((10, 3), dtype=complex)
    out[:4, :] = 1j * np.outer(p, vec)
    return out


def dstar_oneform(vec, p):
    """d^* on a one-form (10,3): -i (p, A|_V) -> (3,)."""
    return -1j * np.einsum('m,ma->a', p.astype(complex), vec[:4, :])


def dstard(vec, p):
    """d^* d componentwise: V part |p|^2 a - p (p,a), W part |p|^2 a."""
    out = np.array(vec, dtype=complex) * np.dot(p, p)
    out[:4, :] -= np.outer(p, np.einsum('m,ma->a', p, vec[:4, :]))
    return out


def ddstar(vec, p):
    """d d^* on a one-form: V part p (p,a); W part 0."""
    out = np.zeros((10, 3), dtype=complex)
    out[:4, :] = np.outer(p, np.einsum('m,ma->a', p, vec[:4, :]))
    return out


def dA2form(vec, p):
    """(dA)_{MN} = i p_M A_N - i p_N A_M, shape (10,10,3)."""
    pp = p10(p)
    return 1j * (np.einsum('m,na->mna', pp, vec) - np.einsum('n,ma->mna', pp, vec))


def dslash(vec, p):
    """Dirac operator on a Sigma-valued field: i A(p) lam  (Sigma^* valued)."""
    return 1j * np.einsum('st,ta->sa', Amat(p10(p)), vec)


def scalar_lap(vec, p):
    return np.dot(p, p) * np.array(vec, dtype=complex)
