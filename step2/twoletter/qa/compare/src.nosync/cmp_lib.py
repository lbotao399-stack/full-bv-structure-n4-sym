"""Comparator library (R8, qa/compare): parse engine A and engine B two-letter kernel tables into one canonical form
in engine A's component basis, so that they can be compared entry by entry.

Basis map (derived from step2/notation/check_dictionary.py sections 1-3, recomputed here, read-only imports):
  vector directions: engine A M (= realisation D index M) <-> engine B index D2B[M], D2B = [1,2,3,0,4,...,9];
  Sigma:   lambda_A = U BS lambda_B ;  Sigma^*: lambda*_A = Us BSs lambda*_B ;
  (U, Us) the intertwiner D -> A with engine A's own convention rho|_Sigma = sb[M] (ea_core, r2_core Yukawa vertex)
  [variant 'dict': rho|_Sigma = s[M] as in check_dictionary.py A1-A3, kept for the convention test].
Letters are linear functionals: A letter ia = BASIS_ia^T lambda_A -> zeta_B = (U BS)^T BASIS_ia.
Outputs are field values: A output BASIS_io -> lambda_B = (U BS)^{-1} BASIS_io.
"""
import os, sys, re, json, itertools
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
QA = os.path.dirname(HERE)                       # qa/compare
TL = os.path.dirname(os.path.dirname(QA))        # step2/twoletter
ROOT = os.path.dirname(os.path.dirname(TL))      # project root
EA = os.path.join(TL, 'engine_a')
EB = os.path.join(TL, 'engine_b')

sys.path.insert(0, os.path.join(EA, 'src.nosync', 'vendor'))
import cliff as CA                                # engine A Clifford blocks (read-only)
_cwd = os.getcwd()
sys.path.insert(0, os.path.join(ROOT, 'engine_b'))
import realization as RB                          # engine B realisation (read-only)

np.seterr(all='ignore')

# ------------------------------------------------------------------ realisation D (check_dictionary.py, section 0-1)
I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EPSU = np.array([[0, 1], [-1, 0]], complex)
EPSD = np.array([[0, -1], [1, 0]], complex)
SIGE = [PAULI[0], PAULI[1], PAULI[2], 1j * I2]
SIGBE = [-PAULI[0], -PAULI[1], -PAULI[2], 1j * I2]
SIGI = [np.array(S.tolist(), dtype=complex) for S in RB.Sig]
SIGBI = [np.array(S.tolist(), dtype=complex) for S in RB.Sigb]


def split(x):
    return x[:8].reshape((2, 4) + x.shape[1:]), x[8:].reshape((2, 4) + x.shape[1:])


def join(a, b):
    return np.concatenate([a.reshape((8,) + a.shape[2:]), b.reshape((8,) + b.shape[2:])])


def rho_S(M, x):
    q, qb = split(x)
    if M < 4:
        s = 1j * np.einsum('di...,da->ai...', qb, SIGBE[M])
        sb = -1j * np.einsum('da,ai...->di...', SIGBE[M], q)
    else:
        I = M - 4
        s = np.einsum('ab,ji,bi...->aj...', EPSU, SIGI[I], q)
        sb = -np.einsum('ji,de,ei...->dj...', SIGBI[I], EPSU, qb)
    return join(s, sb)


def rho_Ss(M, y):
    s, sb = split(y)
    if M < 4:
        q = -1j * np.einsum('ad,di...->ai...', SIGE[M], sb)
        qb = 1j * np.einsum('ai...,ad->di...', s, SIGE[M])
    else:
        I = M - 4
        q = np.einsum('kj,ab,bj...->ak...', SIGBI[I], EPSD, s)
        qb = -np.einsum('kj,ed,dj...->ek...', SIGI[I], EPSD, sb)
    return join(q, qb)


def pair_D(x, y):
    q, qb = split(x)
    s, sb = split(y)
    return np.sum(q * s) - np.sum(qb * sb)


E16 = np.eye(16)
RD = [np.array([rho_S(M, E16[:, k]) for k in range(16)]).T for M in range(10)]
RbD = [np.array([rho_Ss(M, E16[:, k]) for k in range(16)]).T for M in range(10)]
PD = np.array([[pair_D(E16[:, a], E16[:, b]) for b in range(16)] for a in range(16)])

BS = np.zeros((16, 16), complex); BSs = np.zeros((16, 16), complex)
for a in range(2):
    for i in range(4):
        BS[a * 4 + i, a * 4 + i] = 1
        for b in range(2):
            BS[8 + a * 4 + i, 8 + b * 4 + i] = EPSD[a, b]
            BSs[a * 4 + i, 8 + b * 4 + i] = EPSU[a, b]
        BSs[8 + a * 4 + i, a * 4 + i] = 1
D2B = [1, 2, 3, 0, 4, 5, 6, 7, 8, 9]
B2D = [D2B.index(m) for m in range(10)]

# check D == engine B (B1 of check_dictionary) -- recomputed
RBe = [np.array(RB.R[M].tolist(), dtype=complex) for M in range(10)]
RbBe = [np.array(RB.Rb[M].tolist(), dtype=complex) for M in range(10)]
CHECK_B1 = max(np.abs(RD[M] - BSs @ RBe[D2B[M]] @ np.linalg.inv(BS)).max() for M in range(10))


def intertwiner(variant='engineA', flip9=None):
    """(U, Us): Sigma_D -> Sigma_A, Sigma*_D -> Sigma*_A with rho_A|Sigma U = Us rho_D|Sigma, rho_A|Sigma* Us = U rho_D|Sigma*,
    U^T Us = PD.  variant 'engineA': rho_A|Sigma = sb[M] (ea_core / r2_core);  'dict': rho_A|Sigma = s[M] (check_dictionary)."""
    sS = [CA.sb[M] if variant == 'engineA' else CA.s[M] for M in range(10)]       # Sigma -> Sigma*
    sSs = [CA.s[M] if variant == 'engineA' else CA.sb[M] for M in range(10)]      # Sigma* -> Sigma
    if flip9 is None:
        flip9 = (variant == 'engineA')
    sg = [(-1 if (flip9 and M == 9) else 1) for M in range(10)]      # D e_9 <-> -A e_9 when flipped
    RDx = [sg[M] * RD[M] for M in range(10)]; RbDx = [sg[M] * RbD[M] for M in range(10)]
    rows = []
    for M in range(10):
        rows.append(np.hstack([np.kron(sS[M], E16), -np.kron(E16, RDx[M].T)]))
        rows.append(np.hstack([-np.kron(E16, RbDx[M].T), np.kron(sSs[M], E16)]))
    A = np.vstack(rows)
    _u, _s, _vh = np.linalg.svd(A, full_matrices=False)
    null = int(np.sum(_s < 1e-8))
    X = _vh[-1].conj()
    U = X[:256].reshape(16, 16); Us = X[256:].reshape(16, 16)
    mu = U.T @ Us
    kk = np.argwhere(np.abs(PD) > 1e-9)[0]
    ratio = mu[tuple(kk)] / PD[tuple(kk)]
    U = U / np.sqrt(ratio); Us = Us / np.sqrt(ratio)
    res = max([np.abs(sS[M] @ U - Us @ RDx[M]).max() for M in range(10)] +
              [np.abs(sSs[M] @ Us - U @ RbDx[M]).max() for M in range(10)] + [np.abs(U.T @ Us - PD).max()])
    # fix the overall sign: make the largest-modulus entry of U have positive real part (sign cancels in kernels)
    k = np.unravel_index(np.argmax(np.abs(U)), U.shape)
    if U[k].real < 0:
        U, Us = -U, -Us
    return U, Us, null, res, sg


# ------------------------------------------------------------------ engine A bases (ea_core.BASIS, recomputed)
G5S = np.rint((CA.s[0] @ CA.sb[1] @ CA.s[2] @ CA.sb[3]).real).astype(int)
G5SS = np.rint((CA.sb[0] @ CA.s[1] @ CA.sb[2] @ CA.s[3]).real).astype(int)


def _int_basis(P2):
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
    return [np.array(c, dtype=int) for c in cols]


BASIS_A = {
    'lamP': _int_basis(np.eye(16, dtype=int) + G5S),
    'lamM': _int_basis(np.eye(16, dtype=int) - G5S),
    'lamsP': _int_basis(np.eye(16, dtype=int) + G5SS),
    'lamsM': _int_basis(np.eye(16, dtype=int) - G5SS),
}
NA = {'AV': 4, 'AW': 6, 'AsV': 4, 'AsW': 6, 'c': 1, 'cs': 1, 'lamP': 8, 'lamM': 8, 'lamsP': 8, 'lamsM': 8}

# ------------------------------------------------------------------ engine B slots (tl_engine.SLOT / comp_name)
SLOT_B = {'AV': ('A', list(range(4))), 'AW': ('A', list(range(4, 10))),
          'lamP': ('lam', list(range(8))), 'lamM': ('lam', list(range(8, 16))),
          'c': ('c', [0]), 'AsV': ('As', list(range(4))), 'AsW': ('As', list(range(4, 10))),
          'lamsM': ('lams', list(range(8))), 'lamsP': ('lams', list(range(8, 16))), 'cs': ('cs', [0])}
FERMI = {'lamP', 'lamM', 'c', 'AsV', 'AsW'}          # odd letters / outputs (parity of the field)
PARITY = {'AV': 0, 'AW': 0, 'lamP': 1, 'lamM': 1, 'c': 1, 'AsV': 1, 'AsW': 1, 'lamsP': 0, 'lamsM': 0, 'cs': 0}


def comp_name(slot, i):
    sp_, comps = SLOT_B[slot]
    if sp_ in ('A', 'As'):
        nm = 'A' if sp_ == 'A' else 'As'
        return f'{nm}_{i}' if i < 4 else f'{nm}_I{i - 3}'
    if sp_ == 'lam':
        if i < 8:
            return f'lam_a{i // 4}^{i % 4 + 1}'
        return f'lambar_ad{(i - 8) // 4}_{i % 4 + 1}'
    if sp_ == 'lams':
        if i < 8:
            return f'lambars^ad{i // 4}{i % 4 + 1}'
        return f'lams^a{(i - 8) // 4}_{i % 4 + 1}'
    return sp_


NAME2 = {}       # (slot type, name) -> position within the slot comps list
for t, (sp_, comps) in SLOT_B.items():
    for pos, c in enumerate(comps):
        NAME2[(t, comp_name(t, c))] = pos


# ------------------------------------------------------------------ transforms B-slot-coordinates -> A basis
def transforms(variant='engineA', vecmap='D2B'):
    """returns dict (role, type) -> matrix T (n_A x n_Bslot): A_value[iA] = sum_j T[iA, j] B_coeff[j]"""
    U, Us, null, res, sg = intertwiner(variant)
    UB = U @ BS; UsB = Us @ BSs
    vm = D2B if vecmap == 'D2B' else list(range(10))
    # vector sign: A e_M = sg[M] D e_M  (sg[9] = -1 for the engine-A convention)
    T = {}
    for role in ('letter', 'output'):
        for t in ('AV', 'AsV'):
            M = np.zeros((4, 4), complex)
            for i in range(4):
                M[i, SLOT_B[t][1].index(vm[i])] = sg[i]
            T[(role, t)] = M
        for t in ('AW', 'AsW'):
            M = np.zeros((6, 6), complex)
            for i in range(6):
                M[i, SLOT_B[t][1].index(vm[4 + i])] = sg[4 + i]
            T[(role, t)] = M
        for t in ('c', 'cs'):
            T[(role, t)] = np.ones((1, 1), complex)
    for t in ('lamP', 'lamM'):
        comps = SLOT_B[t][1]
        Z = np.array([UB.T @ b for b in BASIS_A[t]])                       # zeta_B rows (8 x 16)
        T[('letter', t)] = Z[:, comps]
        T[('letter_leak', t)] = float(np.abs(np.delete(Z, comps, axis=1)).max())
        L = np.array([np.linalg.solve(UB, b) for b in BASIS_A[t]])         # lambda_B rows (8 x 16)
        T[('output', t)] = L[:, comps]
        T[('output_leak', t)] = float(np.abs(np.delete(L, comps, axis=1)).max())
    for t in ('lamsP', 'lamsM'):
        comps = SLOT_B[t][1]
        L = np.array([np.linalg.solve(UsB, b) for b in BASIS_A[t]])
        T[('output', t)] = L[:, comps]
        T[('output_leak', t)] = float(np.abs(np.delete(L, comps, axis=1)).max())
    T['_meta'] = dict(null=null, res=float(res), vecmap=vecmap, variant=variant, sg=sg)
    T['_U'] = (U, Us)
    return T


# ------------------------------------------------------------------ parsing
_TERM = re.compile(r'^\(([^)]*)\)(.*)$')


def parse_B(s, outs, nlab):
    """engine B string -> {(pos_a, pos_b, pos_o0, ...): {mon: complex}}, mon = tuple of exponents over p{lab}_{m}
    (lab-major, 4 per label, engine B momentum components)."""
    out = {}
    if not s or s.strip() == '0':
        return out
    for term in s.split(' + '):
        m = _TERM.match(term.strip())
        coef = Fraction(m.group(1))
        rest = m.group(2)
        im = False
        ia = ib = None
        opos = [None] * nlab
        mon = [0] * (4 * nlab)
        if rest:
            facs = rest[1:].split('*') if rest.startswith('*') else rest.split('*')
            k = 0
            while k < len(facs):
                f = facs[k]
                e = 1
                if k + 2 < len(facs) and facs[k + 1] == '' :
                    e = int(facs[k + 2]); k += 2
                if f == 'I':
                    assert e == 1
                    im = not im if True else im
                elif f.startswith('Xa['):
                    ia = f[3:-1]
                elif f.startswith('Xb['):
                    ib = f[3:-1]
                elif f.startswith('Y'):
                    lab = int(f[1:f.index('[')]); opos[lab] = f[f.index('[') + 1:-1]
                elif f.startswith('p'):
                    lab, mm = f[1:].split('_'); mon[4 * int(lab) + int(mm)] += e
                else:
                    raise ValueError('token ' + f + ' in ' + term)
                k += 1
        val = complex(coef) * (1j if im else 1)
        key = (ia, ib) + tuple(opos)
        d = out.setdefault(key, {})
        d[tuple(mon)] = d.get(tuple(mon), 0) + val
    return out


def parse_A(s, n):
    """engine A polynomial string in r{j}_{mu} (j = 1..n), 'i' imaginary unit, '^' powers -> {mon: complex}"""
    out = {}
    if s.strip() == '0':
        return out
    for term in s.split(' + '):
        facs = term.strip().split('*')
        c = complex(1)
        mon = [0] * (4 * n)
        for f in facs:
            if f == 'i':
                c *= 1j
            elif f.startswith('r'):
                if '^' in f:
                    f, e = f.split('^'); e = int(e)
                else:
                    e = 1
                j, mu = f[1:].split('_')
                mon[4 * (int(j) - 1) + int(mu)] += e
            else:
                c *= complex(Fraction(f))
        out[tuple(mon)] = out.get(tuple(mon), 0) + c
    return out
