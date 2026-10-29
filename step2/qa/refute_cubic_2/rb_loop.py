"""Own one-loop engine (refuter #2, resumed run).

Graph rule (Ch.3 Prop c3-graphs): vertices = full multilinear vertex functions (Def c1-vertex), edges:
two contracted coordinate functionals brought adjacent by Koszul reordering and replaced by the kernel
component in that order.  Q W (D_0-part) = d/ds W(P + s X, I), s odd, X = K_0 - K_eps
  [from D_0 = 1/2({I,I}_eps - {I,I}) = 1/2(d_X(I^2) - 2 (d_X I) I),  Delta = -d_K]
so the odd X-pair is moved to the FRONT (left action of d_X) before replacement.
Kernel numerators (momentum k entering the vertex of the FIRST slot), Ch.3 Prop c3-P, Def c1-bracket:
  P_AA = 2 delta, P_{c A*} = P_{A* c} = -i k delta (V only), P_{lam lam} = i B(k) delta  [x e^{-t k^2}, t>=eps]
  X_{(x,x*)} = X_{(x*,x)} = -k^2 delta                                                   [x e^{-t k^2}, 0<=t<=eps]
Loop: cycle v_0 -> v_1 -> ... -> v_{n-1} -> v_0, edge j from v_j.next to v_{j+1}.prev, momentum q_j from v_j
to v_{j+1}: q_0 = l, q_j = l + sum_{i=1..j} P_i; first slot = v_j.next, so k = -q_j.
Scaled Schwinger: c_j = order-j Taylor coefficient at eps = 1, X time in [0,1], P times in [1,oo):
  c_j = 1/(16 pi^2) int dtau T^-2 sum_n (-Q0)^n/n! (4T)^-i/i! (nabla^{2i} N)(l0),  i = (deg N - j + 2n)/2,
  l0 = -sum tau_e c_e / T,  Q0 = sum_{e<e'} tau_e tau_e' |c_e - c_e'|^2 / T   (q_e = l + c_e).
Output coefficients are returned WITHOUT the 1/(16 pi^2) and as (rational, coefficient of log(3/2)).
"""
import itertools, warnings
warnings.filterwarnings('ignore')
import numpy as np
from fractions import Fraction as Fr
from math import comb, factorial
from functools import lru_cache
from rb_alg import EPS, a_mat, b_mat, Amat, Bmat, p10, PARITY, SHAPE

DIM = {s: int(np.prod(SHAPE[s])) for s in SHAPE}
ID3 = np.eye(3)


# ---------------------------------------------------------------- vertex tensors (canonical slot order)
def T_AcAst(ks):
    T = np.einsum('MN,abc->MabNc', np.eye(10), EPS)
    return T.reshape(30, 3, 30)


def T_clamlamst(ks):
    T = np.einsum('st,abc->asbtc', np.eye(16), EPS)
    return T.reshape(3, 48, 48)


def T_ccc(ks):
    return EPS.astype(complex)


def T_lamlamA(ks):
    T = np.einsum('Mst,cba->satbMc', a_mat, EPS)
    return T.reshape(48, 48, 30)


def T_AAA(ks):
    k1, k2, k3 = [p10(k) for k in ks]
    d = np.eye(10)
    T = (np.einsum('P,MN->MNP', k2 - k1, d) + np.einsum('M,NP->MNP', k3 - k2, d)
         + np.einsum('N,PM->MNP', k1 - k3, d))
    T = -0.5j * np.einsum('MNP,abc->MaNbPc', T, EPS)
    return T.reshape(30, 30, 30)


def T_AAAA(ks):
    d = np.eye(10)
    C = np.einsum('abe,cde->abcd', EPS, EPS)          # ([e_a,e_b],[e_c,e_d])
    def phi(i, j, k, l):
        # ([A_i ^ A_j],[A_k ^ A_l]) with A_n = e_{M_n} x e_{a_n}
        dd = np.einsum('ik,jl->ijkl', d, d) - np.einsum('il,jk->ijkl', d, d)    # (a_i,a_k)(a_j,a_l)-(a_i,a_l)(a_j,a_k)
        return dd, (i, j, k, l)
    # build T[M1 a1, M2 a2, M3 a3, M4 a4]
    T = np.zeros((10, 3, 10, 3, 10, 3, 10, 3), dtype=complex)
    # pairing (12|34)
    dd = np.einsum('MP,NR->MNPR', d, d) - np.einsum('MR,NP->MNPR', d, d)
    T += np.einsum('MNPR,abcd->MaNbPcRd', dd, C)
    # (13|24): ([A1^A3],[A2^A4]) : (a1,a2)(a3,a4)-(a1,a4)(a3,a2) x ([X1,X3],[X2,X4])
    dd = np.einsum('MN,PR->MNPR', d, d) - np.einsum('MR,PN->MNPR', d, d)
    T += np.einsum('MNPR,acbd->MaNbPcRd', dd, C)
    # (14|23): ([A1^A4],[A2^A3]) : (a1,a2)(a4,a3)-(a1,a3)(a4,a2) x ([X1,X4],[X2,X3])
    dd = np.einsum('MN,RP->MNPR', d, d) - np.einsum('MP,RN->MNPR', d, d)
    T += np.einsum('MNPR,adbc->MaNbPcRd', dd, C)
    return (-0.5 * T).reshape(30, 30, 30, 30)


VERTS = {'AcAst': (['A', 'c', 'Ast'], T_AcAst), 'clamlamst': (['c', 'lam', 'lamst'], T_clamlamst),
         'ccc': (['c', 'c', 'cst'], T_ccc), 'lamlamA': (['lam', 'lam', 'A'], T_lamlamA),
         'AAA': (['A', 'A', 'A'], T_AAA), 'AAAA': (['A', 'A', 'A', 'A'], T_AAAA)}

PSET = {('A', 'A'), ('c', 'Ast'), ('Ast', 'c'), ('lam', 'lam')}
XSET = {('A', 'Ast'), ('Ast', 'A'), ('lam', 'lamst'), ('lamst', 'lam'), ('c', 'cst'), ('cst', 'c')}


def proj_A(part):
    P = np.zeros((10, 10))
    for M in (range(4) if part == 'V' else range(4, 10)):
        P[M, M] = 1
    return np.kron(P, np.eye(3))


def kernel(kind, t1, t2, k, part=None):
    """numerator of the kernel component (t1, t2), k entering the vertex of the first slot."""
    k = np.asarray(k, dtype=float)
    if kind == 'P':
        if (t1, t2) == ('A', 'A'):
            return 2.0 * proj_A(part).astype(complex)
        if (t1, t2) == ('c', 'Ast'):
            m = np.zeros((3, 10, 3), dtype=complex)
            for N in range(4):
                m[:, N, :] = -1j * k[N] * ID3
            return m.reshape(3, 30)
        if (t1, t2) == ('Ast', 'c'):
            m = np.zeros((10, 3, 3), dtype=complex)
            for N in range(4):
                m[N, :, :] = -1j * k[N] * ID3
            return m.reshape(30, 3)
        if (t1, t2) == ('lam', 'lam'):
            return np.kron(1j * Bmat(p10(k)), ID3)
    if kind == 'X':
        k2 = float(k @ k)
        n = DIM[t1]
        if t1 in ('A', 'Ast'):
            return -k2 * proj_A(part).astype(complex)
        return -k2 * np.eye(n, dtype=complex)
    raise ValueError((kind, t1, t2))


def koszul(seq, target, parity):
    odd_seq = [x for x in seq if parity[x]]
    odd_tgt = [x for x in target if parity[x]]
    pos = [odd_seq.index(x) for x in odd_tgt]
    inv = sum(1 for i in range(len(pos)) for j in range(i + 1, len(pos)) if pos[i] > pos[j])
    return -1 if inv % 2 else 1


# ---------------------------------------------------------------- enumeration of cycle graphs
def role_assignments(vtype, ext_types):
    """distinct ways to map roles (ext_1..ext_m, next, prev) to canonical slots of vtype, one per type pattern.
    returns list of (slot_index_for_each_ext, slot_next, slot_prev, type_next, type_prev)."""
    slots, _ = VERTS[vtype]
    n = len(slots)
    m = len(ext_types)
    seen = set()
    out = []
    for perm in itertools.permutations(range(n)):
        exts, nxt, prv = perm[:m], perm[m], perm[m + 1]
        if any(slots[exts[i]] != ext_types[i] for i in range(m)):
            continue
        pattern = (tuple(slots[e] for e in exts), slots[nxt], slots[prv])
        if pattern in seen:
            continue
        seen.add(pattern)
        out.append((exts, nxt, prv, slots[nxt], slots[prv]))
    return out


def enumerate_cycles(vertex_ext_types):
    """vertex_ext_types: list over cycle positions of tuples of external leg slot types.
    Yields dicts describing each graph (one X edge, others P)."""
    choices = []
    for ext_types in vertex_ext_types:
        cj = []
        for vt, (slots, _) in VERTS.items():
            if len(slots) != len(ext_types) + 2:
                continue
            for ra in role_assignments(vt, ext_types):
                cj.append((vt,) + ra)
        choices.append(cj)
    n = len(vertex_ext_types)
    for combo in itertools.product(*choices):
        kinds = []
        ok = True
        for j in range(n):
            t1 = combo[j][4]                 # type of v_j.next
            t2 = combo[(j + 1) % n][5]       # type of v_{j+1}.prev
            if (t1, t2) in PSET:
                kinds.append('P')
            elif (t1, t2) in XSET:
                kinds.append('X')
            else:
                ok = False
                break
        if not ok or kinds.count('X') != 1:
            continue
        yield {'verts': combo, 'kinds': kinds}


def a_parts(graph):
    """all V/W choices for the A-type edges."""
    n = len(graph['verts'])
    types = [(graph['verts'][j][4], graph['verts'][(j + 1) % n][5]) for j in range(n)]
    isA = [t in (('A', 'A'), ('A', 'Ast'), ('Ast', 'A')) for t in types]   # ghost lines (c,A*), (A*,c) are not A-type
    for parts in itertools.product(*[('V', 'W') if a else (None,) for a in isA]):
        yield parts


# ---------------------------------------------------------------- evaluation of N(l)
def graph_N(graph, parts, ext, l):
    """ext: list over cycle positions of lists of (vec_flat, p) for each external leg at that vertex
    (in leg order). Returns (sign * contraction) at loop momentum l (4,)."""
    verts, kinds = graph['verts'], graph['kinds']
    n = len(verts)
    Ptot = [sum((np.asarray(p, dtype=float) for (_, p) in ext[j]), np.zeros(4)) for j in range(n)]
    q = [np.asarray(l, dtype=float)]
    for j in range(1, n):
        q.append(q[-1] + Ptot[j])
    # vertex matrices M_j[next, prev]
    Ms = []
    for j in range(n):
        vt, exts, nxt, prv, tn, tp = verts[j]
        slots, Tf = VERTS[vt]
        moms = [None] * len(slots)
        for i, e in enumerate(exts):
            moms[e] = np.asarray(ext[j][i][1], dtype=float)
        moms[nxt] = -q[j]                    # edge j leaves v_j
        moms[prv] = q[(j - 1) % n]           # edge j-1 enters v_j
        T = Tf(moms)
        # contract externals
        order = list(exts) + [nxt, prv]
        T = np.transpose(T, order)
        for i in range(len(exts)):
            T = np.tensordot(ext[j][i][0], T, axes=(0, 0))
        Ms.append(T)
    # kernels K_j[ v_j.next , v_{j+1}.prev ]
    Ks = []
    for j in range(n):
        t1 = verts[j][4]; t2 = verts[(j + 1) % n][5]
        Ks.append(kernel(kinds[j], t1, t2, -q[j], parts[j]))
    # contraction  sum M0[a,b] K0[a,c] M1[d,c] K1[d,e] M2[f,e] K2[f,b] ...
    # sum M0[a,b] K0[a,c] M1[d,c] K1[d,e] M2[f,e] K2[f,b] = sum(M0 * (K0 M1^T K1 M2^T K2))
    if n == 3:
        val = np.sum(Ms[0] * (Ks[0] @ Ms[1].T @ Ks[1] @ Ms[2].T @ Ks[2]))
    elif n == 2:
        val = np.sum(Ms[0] * (Ks[0] @ Ms[1].T @ Ks[1]))
    else:
        raise ValueError
    # Koszul sign: sequence of coordinate functionals in canonical slot order, vertex by vertex
    seq, par = [], {}
    for j in range(n):
        slots = VERTS[verts[j][0]][0]
        for s in range(len(slots)):
            seq.append((j, s)); par[(j, s)] = PARITY[slots[s]]
    target = []
    jx = kinds.index('X')
    edge_pairs = [((j, verts[j][2]), ((j + 1) % n, verts[(j + 1) % n][3])) for j in range(n)]
    target += list(edge_pairs[jx])
    for j in range(n):
        if j != jx:
            target += list(edge_pairs[j])
    for j in range(n):
        for e in verts[j][1]:
            target.append((j, e))
    sg = koszul(seq, target, par)
    return sg * val


# ---------------------------------------------------------------- exact interpolation in l
def monomials(deg, nvar=4):
    out = []
    for d in range(deg + 1):
        for c in itertools.combinations_with_replacement(range(nvar), d):
            e = [0] * nvar
            for i in c:
                e[i] += 1
            out.append(tuple(e))
    return out


MONO = monomials(4)
PTS = [m for m in MONO]                          # simplex lattice |l|_1 <= 4 (unisolvent for P_4)


def _inv_vandermonde():
    import sympy as sp
    V = sp.Matrix([[sp.Integer(np.prod([p[i] ** m[i] for i in range(4)])) for m in MONO] for p in PTS])
    Vi = V.inv()
    return [[Fr(int(sp.fraction(Vi[i, j])[0]), int(sp.fraction(Vi[i, j])[1])) for j in range(len(PTS))]
            for i in range(len(MONO))]


_VI = None


def interpolate(fn):
    """fn(l) -> complex (exact dyadic). Returns dict monomial -> (Fr re, Fr im)."""
    global _VI
    if _VI is None:
        _VI = _inv_vandermonde()
    vals = [complex(fn(np.array(p, dtype=float))) for p in PTS]
    vr = [Fr(v.real) for v in vals]
    vi = [Fr(v.imag) for v in vals]
    # check dyadic exactness: values must be multiples of 2^-8 with small size
    for v in vals:
        for x in (v.real, v.imag):
            assert abs(x * 256 - round(x * 256)) < 1e-6 and abs(x) < 2 ** 45, x
    poly = {}
    for i, m in enumerate(MONO):
        re = sum((_VI[i][j] * vr[j] for j in range(len(PTS))), Fr(0))
        im = sum((_VI[i][j] * vi[j] for j in range(len(PTS))), Fr(0))
        if re != 0 or im != 0:
            poly[m] = (re, im)
    return poly


def p_eval_check(poly, fn, npts=3, rng=None):
    rng = rng or np.random.default_rng(5)
    worst = 0.0
    for _ in range(npts):
        l = rng.integers(-5, 6, size=4).astype(float)
        v = complex(fn(l))
        pv = sum((float(c[0]) + 1j * float(c[1])) * np.prod([l[i] ** m[i] for i in range(4)])
                 for m, c in poly.items())
        worst = max(worst, abs(v - pv))
    return worst


def lap(poly):
    out = {}
    for m, (re, im) in poly.items():
        for i in range(4):
            if m[i] >= 2:
                e = list(m); e[i] -= 2; e = tuple(e)
                c = m[i] * (m[i] - 1)
                o = out.get(e, (Fr(0), Fr(0)))
                out[e] = (o[0] + c * re, o[1] + c * im)
    return {k: v for k, v in out.items() if v != (0, 0)}


def degree(poly):
    return max((sum(m) for m in poly), default=-1)


# ---------------------------------------------------------------- master integrals
# triangle: J3(a,b,c;m) = int_0^1 dt1 int_1^oo dt2 int_1^oo dt3 t1^a t2^b t3^c / (t1+t2+t3)^m
# bubble:   J2(a,b;m)   = int_0^1 dt1 int_1^oo dt2 t1^a t2^b / (t1+t2)^m
# values as (rational, coefficient of log(3/2)) for J3; J2 gives (rational, coefficient of log 2)
def _poly_mul(p, q):
    out = {}
    for i, a in p.items():
        for j, b in q.items():
            out[i + j] = out.get(i + j, 0) + a * b
    return out


def _poly_pow_binom(c0, c1, n):
    """(c0 + c1 v)^n as dict power->coef."""
    return {k: Fr(comb(n, k)) * (c0 ** (n - k)) * (c1 ** k) for k in range(n + 1)}


@lru_cache(maxsize=None)
def J3(a, b, c, m):
    assert m - b - c >= 3
    # step 1: int_1^oo t3^c (S+t3)^-m = sum_k C(c,k) (-S)^{c-k} (S+1)^{k-m+1}/(m-k-1)
    # step 2: S = t1 + t2, v = S + 1: t2 = v - 1 - t1 ; int_{t1+2}^oo dv (v-1-t1)^b (v-1)^{c-k} v^{k-m+1}
    # represent polynomials in (t1, v) as dict {(i,j): coef} = t1^i v^j
    total_rat = Fr(0)
    total_log = Fr(0)
    for k in range(c + 1):
        pref = Fr(comb(c, k)) * Fr((-1) ** (c - k), m - k - 1)
        # (v - 1 - t1)^b (v - 1)^{c-k}
        P = {(0, 0): Fr(1)}
        for _ in range(b):
            P = {(i + di, j + dj): P.get((i, j), 0) * 0 for (i, j) in P for (di, dj) in [(0, 0)]} if False else P
        # expand explicitly
        P = {(0, 0): Fr(1)}
        for _ in range(b):
            Q = {}
            for (i, j), co in P.items():
                for (di, dj, cc) in [(0, 1, 1), (0, 0, -1), (1, 0, -1)]:
                    key = (i + di, j + dj)
                    Q[key] = Q.get(key, 0) + co * cc
            P = Q
        for _ in range(c - k):
            Q = {}
            for (i, j), co in P.items():
                for (di, dj, cc) in [(0, 1, 1), (0, 0, -1)]:
                    key = (i + di, j + dj)
                    Q[key] = Q.get(key, 0) + co * cc
            P = Q
        # int_{t1+2}^oo v^{j + k - m + 1} dv = (t1+2)^{j+k-m+2} / (m-k-j-2)
        for (i, j), co in P.items():
            if co == 0:
                continue
            e = j + k - m + 2                   # exponent of (t1+2), negative
            assert e <= -1
            c2 = pref * co / (m - k - j - 2)
            # step 3: int_0^1 t1^{a+i} (t1+2)^e dt1 ; w = t1 + 2 in [2,3]: (w-2)^{a+i} w^e
            n = a + i
            for r in range(n + 1):
                cr = Fr(comb(n, r)) * Fr(-2) ** (n - r)       # coefficient of w^r
                ex = r + e
                if ex == -1:
                    total_log += c2 * cr
                else:
                    total_rat += c2 * cr * (Fr(3) ** (ex + 1) - Fr(2) ** (ex + 1)) / (ex + 1) if ex + 1 >= 0 else \
                        c2 * cr * (Fr(1, 3 ** (-(ex + 1))) - Fr(1, 2 ** (-(ex + 1)))) / (ex + 1)
    return total_rat, total_log


@lru_cache(maxsize=None)
def J2(a, b, m):
    """int_0^1 dt1 int_1^oo dt2 t1^a t2^b (t1+t2)^-m, m - b >= 2. Returns (rational, coef of log 2)."""
    assert m - b >= 2
    # inner: t2 = u - t1, u from 1+t1: sum_k C(b,k) (-t1)^{b-k} (1+t1)^{k-m+1}/(m-k-1)
    rat = Fr(0); lg = Fr(0)
    for k in range(b + 1):
        pref = Fr(comb(b, k)) * Fr((-1) ** (b - k), m - k - 1)
        e = k - m + 1                       # exponent of (1+t1), <= -1
        n = a + b - k                       # power of t1
        # int_0^1 t1^n (1+t1)^e ; w = 1 + t1 in [1,2]: (w-1)^n w^e
        for r in range(n + 1):
            cr = Fr(comb(n, r)) * Fr(-1) ** (n - r)
            ex = r + e
            if ex == -1:
                lg += pref * cr
            else:
                val = (Fr(2) ** (ex + 1) - 1) / (ex + 1) if ex + 1 >= 0 else (Fr(1, 2 ** (-(ex + 1))) - 1) / (ex + 1)
                rat += pref * cr * val
    return rat, lg


# ---------------------------------------------------------------- assemble c_j
def poly_subst_l0(poly, w):
    """substitute l = sum_e x_e w_e (w: list of numeric vectors, rational) into poly (dict monomial->(re,im));
    returns dict: tuple of exponents of x_e -> (re, im)."""
    E = len(w)
    out = {}
    wf = [[Fr(x) for x in we] for we in w]
    for m, (re, im) in poly.items():
        # prod_i l_i^{m_i}, l_i = sum_e x_e w_e[i]
        terms = {tuple([0] * E): Fr(1)}
        for i in range(4):
            for _ in range(m[i]):
                new = {}
                for ex, co in terms.items():
                    for e in range(E):
                        if wf[e][i] == 0:
                            continue
                        ex2 = list(ex); ex2[e] += 1; ex2 = tuple(ex2)
                        new[ex2] = new.get(ex2, 0) + co * wf[e][i]
                terms = new
        for ex, co in terms.items():
            o = out.get(ex, (Fr(0), Fr(0)))
            out[ex] = (o[0] + co * re, o[1] + co * im)
    return out


def coefficient(poly, cvecs, jorder, x_index, kind):
    """Order-j Taylor coefficient c_j (without 1/16pi^2) for a cycle with edge offsets cvecs (q_e = l + c_e).
    x_index: index of the X edge (gets tau in [0,1]); kind: 'tri' or 'bub'.
    Returns ((rat_re, log_re), (rat_im, log_im))."""
    E = len(cvecs)
    dN = degree(poly)
    w = [[-Fr(x) for x in c] for c in cvecs]
    # Q0 = sum_{e<e'} tau_e tau_e' |c_e - c_e'|^2 / T
    q0 = {}
    for e in range(E):
        for f in range(e + 1, E):
            d2 = sum((Fr(cvecs[e][i]) - Fr(cvecs[f][i])) ** 2 for i in range(4))
            if d2 != 0:
                ex = [0] * E; ex[e] += 1; ex[f] += 1
                q0[tuple(ex)] = q0.get(tuple(ex), 0) + d2
    # integrand: sum over n: (-Q0)^n/n! * (4T)^-i/i! * (lap^i N)(l0) * T^-2
    # monomial tau^ex / T^power  -> masters
    acc = {}          # (ex tuple, Tpower) -> (re, im)
    lapN = [poly]
    for _ in range(3):
        lapN.append(lap(lapN[-1]))
    n = 0
    while True:
        twoi = dN - jorder + 2 * n
        if twoi > dN:
            break
        if twoi < 0 or twoi % 2:
            n += 1
            continue
        i = twoi // 2
        Li = lapN[i]
        if Li:
            sub = poly_subst_l0(Li, w)          # in x_e = tau_e/T
            # (-Q0)^n
            qn = {tuple([0] * E): Fr(1)}
            for _ in range(n):
                new = {}
                for ex, co in qn.items():
                    for ex2, co2 in q0.items():
                        k = tuple(a + b for a, b in zip(ex, ex2))
                        new[k] = new.get(k, 0) + co * co2
                qn = new
            fac = Fr((-1) ** n, factorial(n)) * Fr(1, 4 ** i * factorial(i))
            for exq, cq in qn.items():
                for exs, (re, im) in sub.items():
                    ex = tuple(a + b for a, b in zip(exq, exs))
                    Tpow = 2 + n + i + sum(exs)      # T^-2 * T^-n (from Q0) * T^-i * T^-(deg x)
                    key = (ex, Tpow)
                    o = acc.get(key, (Fr(0), Fr(0)))
                    acc[key] = (o[0] + fac * cq * re, o[1] + fac * cq * im)
        n += 1
    rr = Fr(0); rl = Fr(0); ir = Fr(0); il = Fr(0)
    for (ex, Tpow), (re, im) in acc.items():
        if re == 0 and im == 0:
            continue
        others = [e for e in range(E) if e != x_index]
        if kind == 'tri':
            mr, ml = J3(ex[x_index], ex[others[0]], ex[others[1]], Tpow)
        else:
            mr, ml = J2(ex[x_index], ex[others[0]], Tpow)
        rr += re * mr; rl += re * ml; ir += im * mr; il += im * ml
    return (rr, rl), (ir, il)
