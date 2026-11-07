"""Engine B -- graph census by explicit Wick enumeration, with explicit su(2), su(3) colour and O(6) flavour.

Canonical normalisation (derived in notes.txt, section F):  S/hbar = -S_E,
  S_E = int 1/4 F^a_mn F^a_mn + 1/2 (D_m phi_I)^a (D_m phi_I)^a + 1/2 (d.A)^2 (Feynman gauge)
        + g^2/4 f^{abe} f^{cde} phi^a_I phi^b_J phi^c_I phi^d_J + fermions,
  D_m phi^a = d_m phi^a + g f^{abc} A^b_m phi^c,  so the vertices of -S_E at order g^2 are
  V3 = -g f^{abc} (d_m phi^a_I) A^b_m phi^c_I,
  Vs = -(g^2/2) f^{abc} f^{ab'c'} A^b_m A^b'_m phi^c_I phi^c'_I            (seagull),
  V4 = -(g^2/4) f^{abe} f^{cde} phi^a_I phi^b_J phi^c_I phi^d_J            (quartic).
  Propagators <phi^a_I phi^b_J> = delta^ab delta_IJ D, <A^a_m A^b_n> = delta^ab delta_mn D (same regulated D).
The operators are O_k = a_k,I b_k,J phi^c_I phi^c_J at x (k = 1) and at the origin (k = 2).
The fermion loop is treated separately (Pfaffian; fermion_trace_check below).
"""
import itertools
import numpy as np
from fractions import Fraction


# ---------------------------------------------------------------- Lie algebra data
def sun_f(N):
    """Structure constants f^{abc} of su(N) for hermitian generators T^a, tr T^a T^b = delta/2,
    [T^a, T^b] = i f^{abc} T^c; real, totally antisymmetric, f^{acd} f^{bcd} = N delta^{ab}."""
    Ts = []
    for j in range(N):
        for k in range(j + 1, N):
            m = np.zeros((N, N), complex)
            m[j, k] = m[k, j] = 0.5
            Ts.append(m)
            m = np.zeros((N, N), complex)
            m[j, k] = -0.5j
            m[k, j] = 0.5j
            Ts.append(m)
    for l in range(1, N):
        m = np.zeros((N, N), complex)
        for j in range(l):
            m[j, j] = 1.0
        m[l, l] = -l
        m *= 1.0 / np.sqrt(2.0 * l * (l + 1))
        Ts.append(m)
    d = len(Ts)
    f = np.zeros((d, d, d))
    for a in range(d):
        for b in range(d):
            comm = Ts[a] @ Ts[b] - Ts[b] @ Ts[a]
            for c in range(d):
                f[a, b, c] = (-2j * np.trace(comm @ Ts[c])).real
    cas = np.einsum('acd,bcd->ab', f, f)
    assert np.allclose(cas, N * np.eye(d))
    return f


# ---------------------------------------------------------------- vertices
def make_vertices(f, kind):
    """Return the list of vertices (name, position, coefficient symbol, fields, colour tensor, flavour pairs).
    Each field: (kind 'phi'/'A', deriv flag).  Flavour pairs: pairs of field slots joined by delta_IJ."""
    d = f.shape[0]
    delta = np.eye(d)
    T4 = np.einsum('pqe,rse->pqrs', f, f)
    Ts = np.einsum('abc,ade->bdce', f, f)  # legs (A^b, A^b', phi^c, phi^c')
    V = {
        'V4': dict(fields=[('phi', 0)] * 4, colour=T4, flav=[(0, 2), (1, 3)], coef=('g2', Fraction(-1, 4))),
        'V3': dict(fields=[('phi', 1), ('A', 0), ('phi', 0)], colour=f, flav=[(0, 2)], coef=('g', Fraction(-1))),
        'Vs': dict(fields=[('A', 0), ('A', 0), ('phi', 0), ('phi', 0)], colour=Ts, flav=[(2, 3)],
                   coef=('g2', Fraction(-1, 2))),
    }
    return V[kind]


def operator_vertex(f):
    d = f.shape[0]
    return dict(fields=[('phi', 0), ('phi', 0)], colour=np.eye(d), flav=None, coef=('1', Fraction(1)))


def perfect_matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        b = items[i]
        rest = items[1:i] + items[i + 1:]
        for m in perfect_matchings(rest):
            yield [(a, b)] + m


def enumerate_graphs(f, vertex_kinds, flav_vectors):
    """Enumerate all Wick contractions of O1(x) O2(o) prod(vertices); return list of dicts with
    edges, derivative markers, colour value, flavour value (for each flavour vector set), Lorentz factor,
    combinatorial weight (1/n! of identical vertices times vertex coefficients)."""
    verts = [('x', operator_vertex(f)), ('o', operator_vertex(f))]
    for i, k in enumerate(vertex_kinds):
        verts.append(('z%d' % (i + 1), make_vertices(f, k)))
    # field list
    fields = []
    for vi, (name, V) in enumerate(verts):
        for si, (kd, dv) in enumerate(V['fields']):
            fields.append((vi, si, kd, dv))
    phis = [i for i, fl in enumerate(fields) if fl[2] == 'phi']
    As = [i for i, fl in enumerate(fields) if fl[2] == 'A']
    import math
    nfact = 1
    for k in set(vertex_kinds):
        nfact *= math.factorial(vertex_kinds.count(k))
    out = []
    for mphi in perfect_matchings(phis):
        for mA in perfect_matchings(As):
            edges = mphi + mA
            # connectivity of the vertex graph
            adj = {i: set() for i in range(len(verts))}
            for a, b in edges:
                va, vb = fields[a][0], fields[b][0]
                adj[va].add(vb)
                adj[vb].add(va)
            seen = {0}
            stack = [0]
            while stack:
                v = stack.pop()
                for w in adj[v]:
                    if w not in seen:
                        seen.add(w)
                        stack.append(w)
            if len(seen) != len(verts):
                continue
            # colour value: einsum with one label per edge
            label = {}
            letters = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
            for ei, (a, b) in enumerate(edges):
                label[a] = letters[ei]
                label[b] = letters[ei]
            ops, subs = [], []
            for vi, (name, V) in enumerate(verts):
                sl = [fi for fi, fl in enumerate(fields) if fl[0] == vi]
                ops.append(V['colour'])
                subs.append(''.join(label[fi] for fi in sl))
            colour = float(np.einsum(','.join(subs) + '->', *ops))
            # flavour value for each vector set
            flvals = []
            for vecs in flav_vectors:
                # union-find over phi slots: vertex flavour deltas and edges
                parent = {i: i for i in phis}

                def find(i):
                    while parent[i] != i:
                        parent[i] = parent[parent[i]]
                        i = parent[i]
                    return i

                def union(i, j):
                    parent[find(i)] = find(j)

                for a, b in mphi:
                    union(a, b)
                for vi, (name, V) in enumerate(verts):
                    if V['flav']:
                        sl = [fi for fi, fl in enumerate(fields) if fl[0] == vi]
                        for (p, q) in V['flav']:
                            union(sl[p], sl[q])
                # each chain contains exactly two operator slots (or none -> closed loop factor 6)
                chains = {}
                for i in phis:
                    chains.setdefault(find(i), []).append(i)
                val = 1.0
                for r, mem in chains.items():
                    opslots = [i for i in mem if fields[i][0] in (0, 1)]
                    if len(opslots) == 0:
                        val *= 6.0
                    elif len(opslots) == 2:
                        u = vecs[(fields[opslots[0]][0], fields[opslots[0]][1])]
                        w = vecs[(fields[opslots[1]][0], fields[opslots[1]][1])]
                        val *= float(np.dot(u, w))
                    else:
                        raise RuntimeError('flavour chain with %d operator slots' % len(opslots))
                flvals.append(val)
            # Lorentz: seagull A-A self loop gives 4; V3-V3 A edge gives a dot product of derivatives
            lor = 1.0
            dots = []
            for a, b in mA:
                va, vb = fields[a][0], fields[b][0]
                if va == vb:
                    lor *= 4.0
                else:
                    dots.append((va, vb))
            # derivative markers: for each V3, the edge (neighbour vertex) carrying its derivative
            derivs = {}
            for a, b in mphi:
                for s, o in ((a, b), (b, a)):
                    if fields[s][3]:
                        derivs[fields[s][0]] = (fields[o][0], o)
            coef = Fraction(1, nfact)
            gpow = 0
            for vi, (name, V) in enumerate(verts):
                cs, cv = V['coef']
                coef *= cv
                gpow += {'1': 0, 'g': 1, 'g2': 2}[cs]
            out.append(dict(edges=[(fields[a][0], fields[b][0], fields[a][2]) for a, b in edges],
                            phi_edges=[(a, b) for a, b in mphi], colour=colour, flav=flvals, lor=lor,
                            dots=dots, derivs=derivs, coef=coef, gpow=gpow,
                            names=[n for n, _ in verts]))
    return out


def flavour_basis(vecs):
    a1, b1, a2, b2 = vecs[(0, 0)], vecs[(0, 1)], vecs[(1, 0)], vecs[(1, 1)]
    S12 = np.dot(a1, a2) * np.dot(b1, b2) + np.dot(a1, b2) * np.dot(b1, a2)
    A12 = np.dot(a1, a2) * np.dot(b1, b2) - np.dot(a1, b2) * np.dot(b1, a2)
    T12 = np.dot(a1, b1) * np.dot(a2, b2)
    return S12, A12, T12


def random_flavour_sets(nsets=6, seed=12345):
    rng = np.random.default_rng(seed)
    sets = []
    for _ in range(nsets):
        sets.append({(0, 0): rng.normal(size=6), (0, 1): rng.normal(size=6),
                     (1, 0): rng.normal(size=6), (1, 1): rng.normal(size=6)})
    return sets


def project_flavour(values, sets):
    """Least-squares coefficients of values on (S12, A12, T12) and the residual."""
    Mtx = np.array([flavour_basis(s) for s in sets])
    coef, res, rk, sv = np.linalg.lstsq(Mtx, np.array(values), rcond=None)
    resid = float(np.max(np.abs(Mtx @ coef - np.array(values))))
    return coef, resid


# ---------------------------------------------------------------- 10d chiral gamma matrices
def gamma10():
    """16x16 chiral blocks gamma^M, gammabar^M (M = 0..9) of a Euclidean Cl(10) representation:
    Gamma^M = [[0, gamma^M], [gammabar^M, 0]], gamma^M gammabar^N + gamma^N gammabar^M = 2 delta^{MN}."""
    s1 = np.array([[0, 1], [1, 0]], complex)
    s2 = np.array([[0, -1j], [1j, 0]], complex)
    s3 = np.array([[1, 0], [0, -1]], complex)
    I2 = np.eye(2)

    def kron(*ms):
        out = np.array([[1.0 + 0j]])
        for m in ms:
            out = np.kron(out, m)
        return out

    G = []
    for k in range(5):
        G.append(kron(*([s3] * k + [s1] + [I2] * (4 - k))))
        G.append(kron(*([s3] * k + [s2] + [I2] * (4 - k))))
    chir = np.eye(32, dtype=complex)
    for M in range(10):
        chir = chir @ G[M]
    chir = chir * (1j) ** 5  # (-i)^5 Gamma^0...Gamma^9 squares to 1
    if not np.allclose(chir @ chir, np.eye(32)):
        chir = chir * -1j
    w, v = np.linalg.eigh((chir + chir.conj().T) / 2)
    P = v[:, w > 0]
    Mn = v[:, w < 0]
    gam = [P.conj().T @ G[M] @ Mn for M in range(10)]   # maps (-) -> (+)
    gbar = [Mn.conj().T @ G[M] @ P for M in range(10)]  # maps (+) -> (-)
    return gam, gbar


def fermion_trace_check():
    """tr_16(gamma^M gammabar^N gamma^P gammabar^Q) = 16 (d^MN d^PQ - d^MP d^NQ + d^MQ d^NP) on all index sets
    used (M,P in R^4, N,Q in W^perp, and all of 0..9 on a sample)."""
    gam, gbar = gamma10()
    err = 0.0
    for M in range(10):
        for N in range(10):
            a = gam[M] @ gbar[N] + gam[N] @ gbar[M]
            err = max(err, np.max(np.abs(a - 2 * (M == N) * np.eye(16))))
    rng = np.random.default_rng(1)
    for _ in range(300):
        M, N, P, Q = rng.integers(0, 10, size=4)
        tr = np.trace(gam[M] @ gbar[N] @ gam[P] @ gbar[Q])
        ex = 16 * (int(M == N) * int(P == Q) - int(M == P) * int(N == Q) + int(M == Q) * int(N == P))
        err = max(err, abs(tr - ex))
    # the combination of the scalar self-energy: sum_{m} tr(gamma^m gammabar^I gamma^m gammabar^J) = -64 delta^IJ
    sc = np.zeros((6, 6))
    for I in range(6):
        for J in range(6):
            sc[I, J] = sum(np.trace(gam[m] @ gbar[4 + I] @ gam[m] @ gbar[4 + J]).real for m in range(4))
    return err, sc


# ---------------------------------------------------------------- classification
def classify(g):
    """Class label of a connected graph with vertices 0 = x, 1 = o, 2.. = internal."""
    from collections import Counter
    E = Counter()
    for a, b, k in g['edges']:
        E[(min(a, b), max(a, b), k)] += 1
    nint = len(g['names']) - 2
    if nint == 1:
        z = 2
        if E[(0, z, 'phi')] == 2 and E[(1, z, 'phi')] == 2:
            return 'QU'
        if E[(0, 1, 'phi')] == 1 and E[(0, z, 'phi')] == 1 and E[(1, z, 'phi')] == 1:
            if E[(z, z, 'phi')] == 1:
                return 'SE_seagull_perp'
            if E[(z, z, 'A')] == 1:
                return 'SE_seagull_V'
        return 'other'
    if nint == 2:
        if E[(2, 3, 'A')] != 1:
            return 'other'
        if all(E[(0, z, 'phi')] == 1 and E[(1, z, 'phi')] == 1 for z in (2, 3)):
            return 'EX'
        if E[(0, 1, 'phi')] == 1 and E[(2, 3, 'phi')] == 1:
            if (E[(0, 2, 'phi')] == 1 and E[(1, 3, 'phi')] == 1) or (E[(0, 3, 'phi')] == 1 and E[(1, 2, 'phi')] == 1):
                return 'SE_bubble'
        return 'other'
    return 'other'


def census_summary(groups=(2, 3)):
    """Per class and gauge algebra: total colour-flavour weight / (C_A dim g) projected on (S12, A12, T12),
    in units of g^2 (the vertex coefficients included; g^1 per V3)."""
    sets = random_flavour_sets()
    summary = {}
    patterns = {}
    for N in groups:
        f = sun_f(N)
        dim = f.shape[0]
        for kinds in (['V4'], ['Vs'], ['V3', 'V3']):
            for g in enumerate_graphs(f, kinds, sets):
                cl = classify(g)
                w = float(g['coef']) * g['colour'] * g['lor'] / (N * dim)
                vals = [w * v for v in g['flav']]
                key = (N, cl, tuple(kinds))
                summary.setdefault(key, {'n': 0, 'vals': np.zeros(len(sets)), 'zero_colour': 0})
                summary[key]['n'] += 1
                summary[key]['vals'] += np.array(vals)
                if abs(g['colour']) < 1e-12:
                    summary[key]['zero_colour'] += 1
                if cl in ('EX', 'SE_bubble') and abs(g['colour']) > 1e-12:
                    coef, res = project_flavour(vals, sets)
                    patterns.setdefault((N, cl), []).append(dict(g=g, coef=coef, res=res))
    out = {}
    for key, v in summary.items():
        coef, res = project_flavour(list(v['vals']), sets)
        out[key] = dict(n=v['n'], zero_colour=v['zero_colour'], S12=coef[0], A12=coef[1], T12=coef[2], resid=res)
    return out, patterns


def ex_integrand_terms(patterns, tsyms):
    """EX numerator from the enumerated patterns: list of (coeff_S12, [(d1, d2)]) for gauss.gaussian with
    internal vertices z1 -> 0, z2 -> 1 and edges (x,z1,t1), (z1,o,t2), (x,z2,t3), (z2,o,t4), (z1,z2,t5)."""
    import sympy as sp
    from gauss import vec_deriv
    t1, t2, t3, t4, t5 = tsyms
    tmap = {(0, 2): t1, (1, 2): t2, (0, 3): t3, (1, 3): t4}
    pos = {0: 'x', 1: 'o', 2: 0, 3: 1}
    acc = {}
    for p in patterns:
        g = p['g']
        cS = p['coef'][0]
        d = []
        for z in (2, 3):
            nb = g['derivs'][z][0]
            d.append((z, nb))
        key = tuple(d)
        acc[key] = acc.get(key, 0.0) + cS
    terms = []
    for key, c in acc.items():
        (z1, n1), (z2, n2) = key
        cr = sp.nsimplify(round(c, 10))
        if cr == 0:
            continue
        va = vec_deriv(pos[z1], pos[n1], tmap[(min(z1, n1), max(z1, n1))])
        vb = vec_deriv(pos[z2], pos[n2], tmap[(min(z2, n2), max(z2, n2))])
        terms.append((cr, [(va, vb)]))
    return terms, acc


def se_bubble_terms(patterns, ssyms):
    """SE bubble line numerator from the enumerated patterns, relabelled so that u (internal 0) is the vertex
    adjacent to x and v (internal 1) the vertex adjacent to o; edges (x,u,s1), (u,v,s2 scalar), (u,v,s3 gluon),
    (v,o,s4).  Returns terms for gauss.gaussian (coefficients per g^2 C_A dim g S12, both lines included)."""
    import sympy as sp
    from gauss import vec_deriv
    s1, s2, s3, s4 = ssyms
    acc = {}
    for p in patterns:
        g = p['g']
        cS = p['coef'][0]
        adjx = [b if a == 0 else a for (a, b, k) in g['edges'] if k == 'phi' and 0 in (a, b) and 1 not in (a, b)]
        u = adjx[0]
        v = 5 - u
        role = {0: 'x', 1: 'o', u: 0, v: 1}
        d = []
        for z in (u, v):
            nb = g['derivs'][z][0]
            d.append((role[z], role[nb]))
        key = tuple(d)
        acc[key] = acc.get(key, 0.0) + cS
    terms = []
    tm = {frozenset(['x', 0]): s1, frozenset([0, 1]): s2, frozenset([1, 'o']): s4}
    for key, c in acc.items():
        cr = sp.nsimplify(round(c, 10))
        if cr == 0:
            continue
        (za, na), (zb, nb) = key
        va = vec_deriv(za, na, tm[frozenset([za, na])])
        vb = vec_deriv(zb, nb, tm[frozenset([zb, nb])])
        terms.append((cr, [(va, vb)]))
    return terms, acc
