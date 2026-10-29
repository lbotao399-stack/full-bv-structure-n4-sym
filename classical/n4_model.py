"""
n4_model.py -- explicit internal realization of the classical BV data of 4d N=4 SYM
obtained from 10d N=(1,0) SYM [ESW (31)] by dimensional reduction R^10 -> R^4
[ESW Def 1.55, Prop 1.56, Prop 2.20, Sec 2.6, Sec 10.3].

Explicit bases live ONLY here (and in checks.py); reported formulas are index-free.

Conventions (all cited in classical_spec.md):
  * V_10 = R^10 Euclidean, orthonormal basis e_M, M = 0..9.
    V (the R^4 directions) = M in {0,1,2,3};  W_R^perp (six transverse) = M in {4..9}.
    Fields are constant along W_R^perp: d_M = 0 for M >= 4.
  * Clifford [ESW Sec 2.1]: v1 v2 + v2 v1 = 2 (v1, v2).
    Sigma = S_+^{10} (16-dim), Sigma^* = S_-^{10} (16-dim), pairing (s, s*) = s^T s* (identity).
    rho(e_M) : Sigma -> Sigma^*  is the 16x16 matrix  sig[M];
    rho(e_M) : Sigma^* -> Sigma  is the 16x16 matrix  sigbar[M];
    sig[M] sigbar[N] + sig[N] sigbar[M] = 2 delta_MN = sigbar[M] sig[N] + sigbar[N] sig[M];
    all sig[M], sigbar[M] symmetric, so (rho(v)Q1, Q2) = (Q1, rho(v)Q2) [ESW Sec 2.1].
    Euclidean signature forces one complex direction: sig[9] = i*1, sigbar[9] = -i*1
    (placed in W_R^perp so that the 4d Dirac operator has real coefficients).
  * Gamma pairing [ESW (28)]: Gamma(Q1,Q2)^M = (rho(e_M) Q1, Q2) = Q2^T sig[M] Q1.
  * Lie algebra g: orthonormal basis T_i of the compact real form, [T_i,T_j] = f_ijk T_k,
    f real and totally antisymmetric; invariant pairing (T_i,T_j) = kappa delta_ij, kappa > 0
    (the coupling sits in (-,-): gap G12 option a).
"""
from fractions import Fraction
import itertools
import numpy as np
import sympy as sp

# ----------------------------------------------------------------------------------
# 10d Euclidean chiral gamma blocks
# ----------------------------------------------------------------------------------
_I2 = np.array([[1, 0], [0, 1]], dtype=np.int64)
_X2 = np.array([[0, 1], [1, 0]], dtype=np.int64)
_Z2 = np.array([[1, 0], [0, -1]], dtype=np.int64)
_E2 = np.array([[0, 1], [-1, 0]], dtype=np.int64)   # i*sigma_2, real antisymmetric
_PAULI = {'I': _I2, 'X': _X2, 'Z': _Z2, 'E': _E2}


def _kron(labels):
    m = np.array([[1]], dtype=np.int64)
    for l in labels:
        m = np.kron(m, _PAULI[l])
    return m


def _anticommute(a, b):
    # tensor products of {I,X,Z,E}: X,Z,E pairwise anticommute, I commutes with all
    n = sum(1 for x, y in zip(a, b) if x != 'I' and y != 'I' and x != y)
    return n % 2 == 1


def _real_symmetric_cl9():
    """Nine mutually anticommuting real symmetric 16x16 matrices squaring to +1
    (Cl(9,0) real representation), found by deterministic backtracking over
    4-fold tensor products of {I,X,Z,E} with an even number of E's."""
    cands = []
    for labels in itertools.product('IXZE', repeat=4):
        if labels == ('I',) * 4:
            continue
        if sum(1 for l in labels if l == 'E') % 2 == 0:
            cands.append(labels)

    best = []

    def bt(start, chosen):
        nonlocal best
        if len(chosen) == 9:
            best = list(chosen)
            return True
        for k in range(start, len(cands)):
            c = cands[k]
            if all(_anticommute(c, d) for d in chosen):
                chosen.append(c)
                if bt(k + 1, chosen):
                    return True
                chosen.pop()
        return False

    bt(0, [])
    assert len(best) == 9
    return [_kron(l) for l in best], best


_GAM9, GAMMA_LABELS = _real_symmetric_cl9()

DIM_V10 = 10
V4 = (0, 1, 2, 3)            # the R^4 directions (V)
WPERP = (4, 5, 6, 7, 8, 9)   # W_R^perp (six transverse directions)
SPIN = 16                    # dim S_+^{10} = dim S_-^{10}

# Gaussian-integer matrices: represent as pair (real int matrix, imag int matrix)
SIG_RE = [g.copy() for g in _GAM9] + [np.zeros((16, 16), dtype=np.int64)]
SIG_IM = [np.zeros((16, 16), dtype=np.int64) for _ in range(9)] + [np.eye(16, dtype=np.int64)]
SIGBAR_RE = [g.copy() for g in _GAM9] + [np.zeros((16, 16), dtype=np.int64)]
SIGBAR_IM = [np.zeros((16, 16), dtype=np.int64) for _ in range(9)] + [-np.eye(16, dtype=np.int64)]


def sig_sym(M):
    """sympy Matrix of rho(e_M): Sigma -> Sigma^*."""
    return sp.Matrix(SIG_RE[M].tolist()) + sp.I * sp.Matrix(SIG_IM[M].tolist())


def sigbar_sym(M):
    """sympy Matrix of rho(e_M): Sigma^* -> Sigma."""
    return sp.Matrix(SIGBAR_RE[M].tolist()) + sp.I * sp.Matrix(SIGBAR_IM[M].tolist())


def sig_entries(M):
    """sparse list of (a, b, (re, im)) for sig[M]_{ab} (Gaussian integers)."""
    out = []
    for a in range(16):
        for b in range(16):
            r, i = int(SIG_RE[M][a, b]), int(SIG_IM[M][a, b])
            if r or i:
                out.append((a, b, (r, i)))
    return out


def sigbar_entries(M):
    out = []
    for a in range(16):
        for b in range(16):
            r, i = int(SIGBAR_RE[M][a, b]), int(SIGBAR_IM[M][a, b])
            if r or i:
                out.append((a, b, (r, i)))
    return out


# ----------------------------------------------------------------------------------
# Lie algebras (structure constants in an orthonormal basis of the compact form)
# ----------------------------------------------------------------------------------
def structure_constants(name):
    """Returns (dim, f) with f a dict {(i,j,k): sympy number}, [T_i,T_j] = f_ijk T_k."""
    if name == 'u1':
        return 1, {}
    if name == 'su2':
        f = {}
        for (i, j, k) in itertools.permutations(range(3)):
            f[(i, j, k)] = sp.Integer(sp.LeviCivita(i, j, k))
        return 3, f
    if name == 'su3':
        h = sp.Rational(1, 2)
        s = sp.sqrt(3) / 2
        base = {(0, 1, 2): sp.Integer(1), (0, 3, 6): h, (0, 4, 5): -h, (1, 3, 5): h,
                (1, 4, 6): h, (2, 3, 4): h, (2, 5, 6): -h, (3, 4, 7): s, (5, 6, 7): s}
        f = {}
        for (i, j, k), v in base.items():
            for perm in itertools.permutations(range(3)):
                idx = tuple((i, j, k)[p] for p in perm)
                sgn = sp.combinatorics.Permutation(list(perm)).signature()
                f[idx] = sgn * v
        return 8, f
    raise ValueError(name)


def f_int(name):
    """Integer structure constants (only su2 / u1); dict (i,j,k)->int."""
    d, f = structure_constants(name)
    out = {}
    for k, v in f.items():
        assert v.is_integer
        out[k] = int(v)
    return d, out


def lie_tests(name):
    d, f = structure_constants(name)
    return lie_tests_from(d, f)


def lie_tests_from(d, f):
    """generic g: f = {(i,j,k): value}, [T_i,T_j] = sum_k f_ijk T_k, orthonormal for (T_i,T_j) = kappa delta_ij."""

    def F(i, j, k):
        return f.get((i, j, k), sp.Integer(0))
    # total antisymmetry (=> invariance of (T_i,T_j) = kappa delta_ij)
    anti = all(sp.simplify(F(i, j, k) + F(j, i, k)) == 0 and sp.simplify(F(i, j, k) - F(j, k, i)) == 0
               for i in range(d) for j in range(d) for k in range(d))
    # Jacobi: f_ijm f_mkl + f_jkm f_mil + f_kim f_mjl = 0
    jac = True
    for i, j, k, l in itertools.product(range(d), repeat=4):
        s = sum(F(i, j, m) * F(m, k, l) + F(j, k, m) * F(m, i, l) + F(k, i, m) * F(m, j, l) for m in range(d))
        if sp.simplify(s) != 0:
            jac = False
            break
    # tr_g(ad_X ad_Y) in this basis: sum_{k,l} f_ikl f_jlk = -C delta_ij
    kil = sp.Matrix(d, d, lambda i, j: sum(F(i, k, l) * F(j, l, k) for k in range(d) for l in range(d)))
    return {'dim': d, 'totally_antisymmetric': bool(anti), 'jacobi': bool(jac),
            'tr_ad_ad_matrix': str(sp.simplify(kil))}


# ----------------------------------------------------------------------------------
# Field content: (name, Z-degree, fermionic Z/2, total parity, fibre dimension, rep)
# [ESW Def 3.2, Rem 1.33, Rem 1.36, (31)]
# ----------------------------------------------------------------------------------
FIELDS = [
    # name   deg  Z2  parity  dim   description
    ('c',     -1,  0,   1,     1,   'Omega^0(R^4; g[1])'),
    ('A',      0,  0,   0,    10,   'Omega^0(R^4; V_10^vee (x) g) = Omega^1(V-part) + W_R^perp-part'),
    ('lam',    0,  1,   1,    16,   'Omega^0(R^4; Pi Sigma_10 (x) g), Sigma_10 = S_+ (x) W + S_- (x) W^*'),
    ('As',     1,  0,   1,    10,   'antifield of A'),
    ('lams',   1,  1,   0,    16,   'antifield of lambda (Sigma_10^* (x) g)'),
    ('cs',     2,  0,   0,     1,   'antifield of c'),
]
PARITY = {f[0]: f[3] for f in FIELDS}
FDIM = {f[0]: f[4] for f in FIELDS}
PAIRS = [('c', 'cs'), ('A', 'As'), ('lam', 'lams')]
ANTI = {'c': 'cs', 'A': 'As', 'lam': 'lams'}


def component_list(dimg=1):
    """Ordered list of (field, comp, color) components of the BV field space."""
    out = []
    for (name, deg, z2, par, dim, _) in FIELDS:
        for a in range(dim):
            for i in range(dimg):
                out.append((name, a, i))
    return out


# ----------------------------------------------------------------------------------
# Self-tests of the explicit realization
# ----------------------------------------------------------------------------------
def gamma_tests():
    res = {}
    sig = [sig_sym(M) for M in range(10)]
    sigb = [sigbar_sym(M) for M in range(10)]
    one = sp.eye(16)
    ok1 = all((sig[M] * sigb[N] + sig[N] * sigb[M] - 2 * (1 if M == N else 0) * one).is_zero_matrix
              for M in range(10) for N in range(10))
    ok2 = all((sigb[M] * sig[N] + sigb[N] * sig[M] - 2 * (1 if M == N else 0) * one).is_zero_matrix
              for M in range(10) for N in range(10))
    sym = all((sig[M] - sig[M].T).is_zero_matrix and (sigb[M] - sigb[M].T).is_zero_matrix for M in range(10))
    res['clifford_sig_sigbar'] = bool(ok1)
    res['clifford_sigbar_sig'] = bool(ok2)
    res['symmetric_blocks'] = bool(sym)
    # ESW Assumption 3.1 (10d N=(1,0)): sum_M sig^M_(ab sig^M_c)d = 0 (sym. over a,b,c)
    S = [np.array(SIG_RE[M], dtype=object) + 1j * 0 for M in range(10)]
    fierz_ok = True
    # use exact integers: sig = re + i im
    def e(M, a, b):
        return complex(int(SIG_RE[M][a, b]), int(SIG_IM[M][a, b]))
    for a in range(16):
        for b in range(a, 16):
            for c in range(b, 16):
                for d in range(16):
                    tot = 0
                    for (x, y, z) in set(itertools.permutations((a, b, c))):
                        tot += sum(e(M, x, y) * e(M, z, d) for M in range(10))
                    if tot != 0:
                        fierz_ok = False
                        break
                if not fierz_ok:
                    break
            if not fierz_ok:
                break
        if not fierz_ok:
            break
    res['fierz_ESW_assumption_3_1'] = bool(fierz_ok)
    # 4d chirality on S_+^{10}: Omega = sigbar^0 sig^1 sigbar^2 sig^3 ; Omega^2 = 1, eigenspaces 8 + 8
    Om = sigb[0] * sig[1] * sigb[2] * sig[3]
    res['Omega4_squared_is_1'] = bool((Om * Om - one).is_zero_matrix)
    ev = Om.eigenvals()
    res['Omega4_eigenvalues'] = {str(k): int(v) for k, v in ev.items()}
    # 6d chirality on S_+^{10}
    Om6 = sigb[4] * sig[5] * sigb[6] * sig[7] * sigb[8] * sig[9]
    pr = sp.simplify(Om6 * Om)
    res['Omega6_times_Omega4'] = (str(pr[0, 0]) + ' * identity') if (pr - pr[0, 0] * one).is_zero_matrix else 'not scalar' 
    res['gamma_labels'] = [''.join(l) for l in GAMMA_LABELS] + ['i*1']
    return res


if __name__ == '__main__':
    import json
    out = {'gamma': gamma_tests(), 'su2': lie_tests('su2'), 'su3': lie_tests('su3')}
    print(json.dumps(out, indent=1))
