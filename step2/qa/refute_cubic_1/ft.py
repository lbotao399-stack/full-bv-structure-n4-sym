"""Independent field-theory realisation for refute_cubic_1 (written from scratch; imports nothing
from engine_a/, engine_b/, theory/ or any other project code).

Conventions are taken from the manuscript text only:
  ch01: fields, pairing, action (eq c1-action), Q table (prop c1-Q), vertex-function rule (def c1-vertex)
  ch03: Q^GF, propagator components (prop c3-P), one-loop formula (prop c3-oneloop), graph rules
The combinatorial formulas (one-loop operator, tree field, theta trick) were checked on a finite-dimensional
super toy model in toy_conventions.py.

Realisation:  g = su(2), [e_a,e_b] = eps_abc e_c, (X,Y) = X.Y  (so (-,-)^{-1} = delta);
              V_10 = C^10 with M = 0..3 the directions of V (spacetime), M = 4..9 those of W_R^perp;
              Sigma = C^16, Sigma^* = C^16 with the pairing (s, s*) = sum_a s_a s*_a;
              rho(e_M): Sigma -> Sigma^* is gam[M], rho(e_M): Sigma^* -> Sigma is gamt[M].
Grassmann bookkeeping: a leg carries a canonical monomial (sorted tuple of odd generator ids).
Ids: theta_b = -2, theta_c = -1, r (for Q) = 0, r1 = 1, r2 = 2, r3 = 3.
"""
import itertools
import numpy as np

# ---------------------------------------------------------------- Clifford module
_I2 = np.eye(2)
_X = np.array([[0, 1], [1, 0]], float)
_Z = np.array([[1, 0], [0, -1]], float)
_Yr = np.array([[0, 1], [-1, 0]], float)  # real antisymmetric
_PAULI = {'I': _I2, 'X': _X, 'Z': _Z, 'Y': _Yr}


def _string_mat(s):
    m = np.array([[1.0]])
    for ch in s:
        m = np.kron(m, _PAULI[ch])
    return m


def _anticommute(s1, s2):
    n = sum(1 for a, b in zip(s1, s2) if a != 'I' and b != 'I' and a != b)
    return n % 2 == 1


def _find_nine():
    cands = [''.join(t) for t in itertools.product('IXZY', repeat=4)]
    cands = [s for s in cands if s != 'IIII' and s.count('Y') % 2 == 0]

    def bt(chosen, start):
        if len(chosen) == 9:
            return chosen
        for i in range(start, len(cands)):
            s = cands[i]
            if all(_anticommute(s, c) for c in chosen):
                r = bt(chosen + [s], i + 1)
                if r:
                    return r
        return None
    return bt([], 0)


NINE = _find_nine()
gam = np.zeros((10, 16, 16), complex)   # Sigma -> Sigma^*
gamt = np.zeros((10, 16, 16), complex)  # Sigma^* -> Sigma
gam[0] = 1j * np.eye(16)
gamt[0] = -1j * np.eye(16)
for j, s in enumerate(NINE):
    m = _string_mat(s)
    gam[j + 1] = m
    gamt[j + 1] = m


def check_clifford():
    err = 0.0
    for u in range(10):
        for v in range(10):
            a = gamt[u] @ gam[v] + gamt[v] @ gam[u]   # on Sigma
            b = gam[u] @ gamt[v] + gam[v] @ gamt[u]   # on Sigma^*
            t = 2.0 * (u == v) * np.eye(16)
            err = max(err, abs(a - t).max(), abs(b - t).max())
        err = max(err, abs(gam[u] - gam[u].T).max(), abs(gamt[u] - gamt[u].T).max())
    return err


# ---------------------------------------------------------------- colour
fabc = np.zeros((3, 3, 3))
for a, b, c in itertools.permutations(range(3)):
    fabc[a, b, c] = np.linalg.det(np.eye(3)[[a, b, c]])
NC = 3
UNIT = -2.0 / (16 * np.pi ** 2)   # hbar tr(ad ad)/(16 pi^2 (X,Y)) with hbar = 1, tr(ad_X ad_Y) = -2 X.Y

# ---------------------------------------------------------------- slots
DIM = {'A': 30, 'As': 30, 'c': 3, 'cs': 3, 'l': 48, 'ls': 48}
PAR = {'A': 0, 'As': 1, 'c': 1, 'cs': 0, 'l': 1, 'ls': 0}
# flattened index conventions: A, As: M*3 + a ; l, ls: s*3 + a ; c, cs: a


def A_index(M, a):
    return M * 3 + a


def l_index(s, a):
    return s * 3 + a


# fluctuation space (slots met by the propagator): A (30), c (3), A^*|_V (12), lambda (48)
FBLOCKS = [('A', np.arange(30)), ('c', np.arange(3)), ('As', np.array([A_index(M, a) for M in range(4) for a in range(3)])),
           ('l', np.arange(48))]
FOFF = {}
_o = 0
for _s, _idx in FBLOCKS:
    FOFF[_s] = (_o, _idx)
    _o += len(_idx)
NF = _o  # 93
FPAR = np.zeros(NF)
for _s, _idx in FBLOCKS:
    o, idx = FOFF[_s]
    FPAR[o:o + len(idx)] = PAR[_s]
PI = np.diag((-1.0) ** FPAR)


def Phat(q):
    """numerator of P(eps,infty) in momentum space, fluctuation x fluctuation (rows: first slot)."""
    P = np.zeros((NF, NF), complex)
    oA, _ = FOFF['A']
    P[oA:oA + 30, oA:oA + 30] = 2.0 * np.eye(30)
    oc, _ = FOFF['c']
    oS, idxS = FOFF['As']
    for k, ind in enumerate(idxS):
        M, a = divmod(ind, 3)
        # c (x) A*|_V : -(,)^{-1} i q ;  A*|_V (x) c : -(,)^{-1} i q
        P[oc + a, oS + k] = -1j * q[M]
        P[oS + k, oc + a] = -1j * q[M]
    ol, _ = FOFF['l']
    gq = sum(q[mu] * gamt[mu] for mu in range(4))
    P[ol:ol + 48, ol:ol + 48] = np.kron(1j * gq, np.eye(3))
    return P


# ---------------------------------------------------------------- Grassmann monomials
def canon(seq):
    """seq: list of generator ids (odd) in product order -> (sign, sorted tuple) or (0, None)"""
    if len(set(seq)) != len(seq):
        return 0, None
    s = list(seq)
    sign = 1
    for i in range(len(s)):
        for j in range(len(s) - 1 - i):
            if s[j] > s[j + 1]:
                s[j], s[j + 1] = s[j + 1], s[j]
                sign = -sign
    return sign, tuple(s)


# ---------------------------------------------------------------- vertex tensors (densities of I)
def T_AAA(p1):
    """-1/2 sum_{MN} (d_M A_N, [A_M, A_N]) ; positions (A[derivative], A, A); depends on momentum of position 1"""
    T = np.zeros((10, 3, 10, 3, 10, 3), complex)
    for M in range(4):
        if p1[M] == 0:
            continue
        for N in range(10):
            # T[(N,a),(M,b),(N,c)] = -(i/2) p1_M f_{bca}
            T[N, :, M, :, N, :] += -0.5j * p1[M] * np.transpose(fabc, (2, 0, 1))  # [a,b,c] -> f[b,c,a]
    return T.reshape(30, 30, 30)


def _fix_llA():
    T = np.zeros((16, 3, 10, 3, 16, 3), complex)
    for M in range(10):
        x = 0.5 * np.einsum('st,bca->sabtc', gam[M], fabc)  # indices s,a,b,t,c
        T[:, :, M, :, :, :] = x
    return T.reshape(48, 30, 48)


T_llA = _fix_llA()


def _T_AcAs():
    T = np.zeros((10, 3, 3, 10, 3))
    for M in range(10):
        T[M, :, :, M, :] = fabc  # ([A_M, c], A*_M) = f_abc A^a c^b A*^c
    return T.reshape(30, 3, 30)


T_AcAs = _T_AcAs()


def _T_clls():
    T = np.zeros((3, 16, 3, 16, 3))
    for s in range(16):
        T[:, s, :, s, :] = fabc  # ([c, lambda_s], lambda*_s)
    return T.reshape(3, 48, 48)


T_clls = _T_clls()
T_cccs = 0.5 * fabc.copy()


def _T_AAAA():
    ff = np.einsum('abe,cde->abcd', fabc, fabc)
    T = np.zeros((10, 3, 10, 3, 10, 3, 10, 3))
    for M in range(10):
        for N in range(10):
            T[M, :, N, :, M, :, N, :] += -0.125 * ff
    return T.reshape(30, 30, 30, 30)


T_AAAA = _T_AAAA()

# list of interaction terms: (slots, tensor or callable(momenta)->tensor)
CUBIC = [
    (('A', 'A', 'A'), lambda ps: T_AAA(ps[0])),
    (('l', 'A', 'l'), lambda ps: T_llA),
    (('A', 'c', 'As'), lambda ps: T_AcAs),
    (('c', 'l', 'ls'), lambda ps: T_clls),
    (('c', 'c', 'cs'), lambda ps: T_cccs),
]
QUARTIC = [(('A', 'A', 'A', 'A'), lambda ps: T_AAAA)]


# ---------------------------------------------------------------- legs
class Leg:
    def __init__(self, slot, fib, p, mono):
        self.slot = slot
        self.fib = np.asarray(fib, complex).reshape(DIM[slot])
        self.p = np.asarray(p, float)
        self.mono = tuple(mono)

    def __repr__(self):
        return f"Leg({self.slot},p={self.p},m={self.mono})"


def contract(T, vecs):
    out = T
    for v in vecs:
        out = np.tensordot(out, v, axes=([0], [0]))
    return out


def multilinear(terms, legs):
    """Multilinear plane-wave value of sum of density terms on the given legs (each leg once).
    Returns dict canonical monomial -> complex."""
    n = len(legs)
    res = {}
    for slots, Tf in terms:
        if len(slots) != n:
            continue
        for perm in itertools.permutations(range(n)):
            # perm[i] = leg at position i
            if any(legs[perm[i]].slot != slots[i] for i in range(n)):
                continue
            ps = [legs[perm[i]].p for i in range(n)]
            T = Tf(ps)
            val = contract(T, [legs[perm[i]].fib for i in range(n)])
            seq = []
            for i in range(n):
                seq += list(legs[perm[i]].mono)
            sg, m = canon(seq)
            if sg == 0:
                continue
            res[m] = res.get(m, 0) + sg * val
    return res


def mono_parity(m):
    return len(m) % 2


# ---------------------------------------------------------------- Hessian blocks at a background
def hessian(terms, bglegs, kb, kc):
    """H_{beta gamma} = coefficient of theta_c theta_b in I(e + theta_b b_beta + theta_c b_gamma),
    multilinear in the background legs bglegs (1 leg for cubic, 2 for quartic).
    kb, kc: momenta of the fluctuation directions beta, gamma (plane waves e^{i k x}).
    Returns (H as NF x NF complex, canonical monomial of the background)."""
    nb = len(bglegs)
    H = np.zeros((NF, NF), complex)
    mono_out = None
    for slots, Tf in terms:
        n = len(slots)
        if n != nb + 2:
            continue
        # choose positions for background legs (ordered assignment) and for beta, gamma
        for pos in itertools.permutations(range(n)):
            # pos[0..nb-1]: positions of bg legs, pos[nb]: position of beta, pos[nb+1]: position of gamma
            ok = all(slots[pos[i]] == bglegs[i].slot for i in range(nb))
            if not ok:
                continue
            sb, sc = slots[pos[nb]], slots[pos[nb + 1]]
            if sb not in FOFF or sc not in FOFF:
                continue
            ps = [None] * n
            for i in range(nb):
                ps[pos[i]] = bglegs[i].p
            ps[pos[nb]] = kb
            ps[pos[nb + 1]] = kc
            T = Tf(ps)
            # contract background fibres at their positions
            # move axes: bring beta-axis, gamma-axis to the end
            order = [pos[i] for i in range(nb)] + [pos[nb], pos[nb + 1]]
            Tt = np.transpose(T, order)
            mat = contract(Tt, [bglegs[i].fib for i in range(nb)])  # shape (dim sb, dim sc)
            # Koszul: coefficients in position order
            seq = []
            for i in range(n):
                if i == pos[nb]:
                    if PAR[sb]:
                        seq.append(-2)
                elif i == pos[nb + 1]:
                    if PAR[sc]:
                        seq.append(-1)
                else:
                    j = pos.index(i)
                    seq += list(bglegs[j].mono)
            sg, m = canon(seq)
            if sg == 0:
                continue
            # canonical starts with (-2,-1) if both odd: theta_b theta_c X = (-1)^{|b||c|} theta_c theta_b X
            sg *= (-1) ** (PAR[sb] * PAR[sc])
            mb = tuple(x for x in m if x >= 0)
            if mono_out is None:
                mono_out = mb
            assert mono_out == mb
            ob, ib = FOFF[sb]
            oc, ic = FOFF[sc]
            H[ob:ob + len(ib), oc:oc + len(ic)] += sg * mat[np.ix_(ib, ic)]
    return H, mono_out


def gradient(terms, bglegs, kb):
    """g_beta = coefficient of theta_b (in front) in I(e + theta_b b_beta), multilinear in bglegs."""
    nb = len(bglegs)
    g = np.zeros(NF, complex)
    mono_out = None
    for slots, Tf in terms:
        n = len(slots)
        if n != nb + 1:
            continue
        for pos in itertools.permutations(range(n)):
            if not all(slots[pos[i]] == bglegs[i].slot for i in range(nb)):
                continue
            sb = slots[pos[nb]]
            if sb not in FOFF:
                continue
            ps = [None] * n
            for i in range(nb):
                ps[pos[i]] = bglegs[i].p
            ps[pos[nb]] = kb
            T = Tf(ps)
            order = [pos[i] for i in range(nb)] + [pos[nb]]
            vec = contract(np.transpose(T, order), [bglegs[i].fib for i in range(nb)])
            seq = []
            for i in range(n):
                if i == pos[nb]:
                    if PAR[sb]:
                        seq.append(-2)
                else:
                    seq += list(bglegs[pos.index(i)].mono)
            sg, m = canon(seq)
            if sg == 0:
                continue
            mb = tuple(x for x in m if x >= 0)
            if mono_out is None:
                mono_out = mb
            assert mono_out == mb
            ob, ib = FOFF[sb]
            g[ob:ob + len(ib)] += sg * vec[ib]
    return g, mono_out


def tree_field(bglegs, eps):
    """chi^alpha = sum_beta P^{alpha beta}(q) (d_beta I_3)(e_j + e_m), q = p_j + p_m; returns list of Legs (one per slot)."""
    q = bglegs[0].p + bglegs[1].p
    g, m = gradient(CUBIC, bglegs, -q)
    if m is None or not np.any(g):
        return []
    s = np.exp(-eps * q @ q) / (q @ q)
    x = s * (Phat(q) @ g)
    legs = []
    for slot, idx in FBLOCKS:
        o, _ = FOFF[slot]
        part = x[o:o + len(idx)]
        if np.any(abs(part) > 0):
            fib = np.zeros(DIM[slot], complex)
            fib[idx] = part
            legs.append(Leg(slot, fib, q, m))
    return legs


# ---------------------------------------------------------------- Q on legs (prop c1-Q), r prepended
def Q_leg(L, rid=0):
    p = L.p
    p2 = p @ p
    if L.slot == 'c':
        f = np.zeros((10, 3), complex)
        for mu in range(4):
            f[mu] = -1j * p[mu] * L.fib
        newslot = 'A'
        fib = f.reshape(30)
    elif L.slot == 'A':
        A = L.fib.reshape(10, 3)
        f = np.zeros((10, 3), complex)
        pA = sum(p[mu] * A[mu] for mu in range(4))
        for mu in range(4):
            f[mu] = -0.5 * (p2 * A[mu] - p[mu] * pA)
        for w in range(4, 10):
            f[w] = -0.5 * p2 * A[w]
        newslot = 'As'
        fib = f.reshape(30)
    elif L.slot == 'As':
        A = L.fib.reshape(10, 3)
        newslot = 'cs'
        fib = -1j * sum(p[mu] * A[mu] for mu in range(4))
    elif L.slot == 'l':
        lam = L.fib.reshape(16, 3)
        gp = sum(p[mu] * gam[mu] for mu in range(4))
        newslot = 'ls'
        fib = (1j * gp @ lam).reshape(48)
    else:
        return None
    sg, m = canon([rid] + list(L.mono))
    return Leg(newslot, sg * fib, p, m)


# ---------------------------------------------------------------- local quadratic functionals (plane-wave densities)
def T2_AA(ps):          # int (A, A), all ten components
    return np.eye(30)


def T2_AVAV(ps):
    T = np.zeros((10, 3, 10, 3))
    for mu in range(4):
        T[mu, :, mu, :] = np.eye(3)
    return T.reshape(30, 30)


def T2_AWAW(ps):
    T = np.zeros((10, 3, 10, 3))
    for w in range(4, 10):
        T[w, :, w, :] = np.eye(3)
    return T.reshape(30, 30)


def _dAdA(ps, comps):
    p1, p2 = ps
    T = np.zeros((10, 3, 10, 3), complex)
    # (dA1, dA2) Hodge = -( (p1.p2)(A1.A2) - (p2.A1)(p1.A2) )  with i*i = -1
    for M in comps:
        T[M, :, M, :] += -(p1 @ p2) * np.eye(3)
    for M in comps:
        for N in comps:
            pm2 = p2[M] if M < 4 else 0.0
            pn1 = p1[N] if N < 4 else 0.0
            T[M, :, N, :] += pm2 * pn1 * np.eye(3)
    return T.reshape(30, 30)


def T2_dAdA(ps):
    return _dAdA(ps, range(10))


def T2_dAVdAV(ps):
    return _dAdA(ps, range(4))


def T2_dAWdAW(ps):
    return _dAdA(ps, range(4, 10))


def T2_dsAdsA(ps):  # int (d^*A|_V, d^*A|_V), d^* -> -i iota_p
    p1, p2 = ps
    T = np.zeros((10, 3, 10, 3), complex)
    for mu in range(4):
        for nu in range(4):
            T[mu, :, nu, :] = (-1j * p1[mu]) * (-1j * p2[nu]) * np.eye(3)
    return T.reshape(30, 30)


def T2_ldl(ps):  # int (lambda, dslash lambda): positions (l, l), derivative on the second
    p2 = ps[1]
    gp = sum(p2[mu] * gam[mu] for mu in range(4))
    return np.kron(1j * gp, np.eye(3))


def T2_dcAs(ps):  # int (dc, A^*): positions (c, As)
    p1 = ps[0]
    T = np.zeros((3, 10, 3), complex)
    for mu in range(4):
        T[:, mu, :] = 1j * p1[mu] * np.eye(3)
    return T.reshape(3, 30)


def T2_ddc_dsA(ps):  # int (d^*dc, d^*A|_V): positions (c, A)
    p1, p2 = ps
    T = np.zeros((3, 10, 3), complex)
    for mu in range(4):
        T[:, mu, :] = (p1 @ p1) * (-1j * p2[mu]) * np.eye(3)
    return T.reshape(3, 30)


# G basis (cubic)
def T3_G3(ps):  # int ([c, A|_W], d^*d A|_W): positions (c, A, A)
    p3 = ps[2]
    T = np.zeros((3, 10, 3, 10, 3))
    for w in range(4, 10):
        T[:, w, :, w, :] = fabc * (p3 @ p3)
    return T.reshape(3, 30, 30)


def T3_G4(ps):  # int ([c, dc], A^*): positions (c, c, As), A^* paired with the V components
    p2 = ps[1]
    T = np.zeros((3, 3, 10, 3), complex)
    for mu in range(4):
        T[:, :, mu, :] = fabc * (1j * p2[mu])
    return T.reshape(3, 3, 30)


def T3_G1(ps):  # int ([c, A|_V], d^*d A|_V)
    p3 = ps[2]
    T = np.zeros((3, 10, 3, 10, 3), complex)
    for mu in range(4):
        for nu in range(4):
            T[:, mu, :, nu, :] = fabc * ((p3 @ p3) * (mu == nu) - p3[mu] * p3[nu])
    return T.reshape(3, 30, 30)


def T3_G2(ps):  # int ([c, A|_V], d d^* A|_V)
    p3 = ps[2]
    T = np.zeros((3, 10, 3, 10, 3), complex)
    for mu in range(4):
        for nu in range(4):
            T[:, mu, :, nu, :] = fabc * (p3[mu] * p3[nu])
    return T.reshape(3, 30, 30)


# ---------------------------------------------------------------- polynomial numerators of one-loop words
EXPS = [e for e in itertools.product(range(4), repeat=4) if sum(e) <= 3]
EXPS = np.array(EXPS)  # (35, 4)


def poly_basis(K):
    """K: (n,4) -> (n,35) monomials"""
    return np.prod(K[:, None, :] ** EXPS[None, :, :], axis=2)


def lap_coeffs(c):
    """coefficients of sum_mu d^2/dk_mu^2 of the polynomial with coefficients c (over EXPS)"""
    out = np.zeros(len(EXPS), complex)
    idx = {tuple(e): i for i, e in enumerate(EXPS)}
    for i, e in enumerate(EXPS):
        for mu in range(4):
            if e[mu] >= 2:
                e2 = list(e)
                e2[mu] -= 2
                out[idx[tuple(e2)]] += c[i] * e[mu] * (e[mu] - 1)
    return out


def word_numerator(word_ops, kpts):
    """word_ops: list of (callable m(kin) -> NF x NF matrix, momentum shift); returns STr values at kpts"""
    vals = []
    for k in kpts:
        kin = np.array(k, float)
        mat = np.eye(NF, dtype=complex)
        for mfun, shift in word_ops:
            mat = mfun(kin) @ mat
            kin = kin + shift
        vals.append(np.sum(np.diag(mat) * (-1.0) ** FPAR))
    return np.array(vals)


def fit_poly(word_ops, rng, npts=70, maxdeg=3):
    """least-squares reconstruction of the (polynomial) numerator in the monomials of degree <= maxdeg;
    validated at fresh points (a numerator of higher degree would fail the validation)."""
    K = rng.uniform(-1.3, 1.3, size=(npts, 4))
    v = word_numerator(word_ops, K)
    mask = EXPS.sum(1) <= maxdeg
    B = poly_basis(K)[:, mask]
    cm, res, rk, sv = np.linalg.lstsq(B, v, rcond=None)
    c = np.zeros(len(EXPS), complex)
    c[mask] = cm
    K2 = rng.uniform(-1.3, 1.3, size=(12, 4))
    v2 = word_numerator(word_ops, K2)
    err = np.max(np.abs(poly_basis(K2) @ c - v2)) / max(1e-300, np.max(np.abs(v2)) + 1e-300)
    return c, err
