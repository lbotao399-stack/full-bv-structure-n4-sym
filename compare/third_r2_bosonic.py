"""Third, independent exact check (round 2): every purely bosonic one-loop 1PI graph of the five bosonic external sets
(A|_V)^3, A|_V (A|_Wperp)^2, (A|_V)^4, (A|_V)^2 (A|_Wperp)^2, (A|_Wperp)^4, built from scratch from the action.

Inputs from the sources only:  S = int dvol( -1/4 (F_A,F_A) ),  F_A = dA + 1/2[A ^ A],  d only along V  [ESW (31), DESIGN];
weight e^{+S/hbar}, W(P,I) = hbar log(e^{hbar d_P} e^{I/hbar})  [C07/CG2];  P(eps,infty) = int_eps^infty dt e^{-t k^2} x 2 (-,-)^{-1} (x) metric^{-1}
on A - A (Feynman gauge: H = Laplacian on every component, positive A-A covariance, DESIGN / DEFINITIONS_v2 item 1).
Code-level choices (not reported notation): g = su(2), orthonormal basis for (-,-), [e_a,e_b] = eps_abc e_c (tr_g(ad_X ad_Y)/(X,Y) = -2);
V_10 = R^4 (+) R^6 with Euclidean metric, components 0..3 along V, 4..9 along W_R^perp.

Method (independent of both engines' codes): Feynman rules = full functional derivatives of the cubic and quartic terms at plane waves;
graph value = (1/|Aut|) int d^4k/(2pi)^4 N(k) prod_e int_eps^infty dt_e e^{-t_e (k+q_e)^2}; Gaussian k-integral and the scale integral
over T = sum t_e give for the log eps coefficient (exact):
    L_G = -(16 pi^2)^{-1} (1/|Aut|) vol(simplex_{n-1}) [ (4^m m!)^{-1} (d_k.d_k)^m N ](k = -(1/n) sum_e q_e),  m = n - 2,
valid because (d_k.d_k)^m N is affine in k for omega in {0,1}; (d_k.d_k)^{n-1} N = 0 (no eps^{-1}); only integer powers of T occur (no eps^{-1/2}).
Comparison: engines report log eps part of I^CT[eps] as (number) x u x (structure), u = hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) = -2 hbar/(16 pi^2) here, so
    sum_G J_G/|Aut_G| = 2 x (engine value),  J_G := vol (4^m m!)^{-1} (d.d)^m N (k_centroid).
Graph census: brute-force enumeration of cyclic vertex sequences (cubic: 1 external leg, quartic: 2), line species V / W_R^perp, vertex
selection rule (even number of W_R^perp legs), modulo rotation and reflection; |Aut| = 2 for a bubble with two equal lines, else 1.
Run: python3 -I compare/third_r2_bosonic.py   -> compare/third_r2_bosonic.json"""
import itertools, json, random, collections
from fractions import Fraction as Fr
import numpy as np

R = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
VC = list(range(4)); WC = list(range(4, 10))

def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0]) + a[2] * (b[0] * c[1] - b[1] * c[0]))

def cross(a, b):
    return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]

def dot(a, b):
    return sum(x * y for x, y in zip(a, b))

class Leg:
    def __init__(self, c, v, p):
        self.c = c; self.v = v; self.p = list(p) + [0] * (10 - len(p))

# S_3 = -1/2 sum_{M,N} (d_M A_N, [A_M, A_N]),  S_4 = -1/8 sum_{M,N} ([A_M,A_N],[A_M,A_N]);  d_M -> i p_M on e^{i p x}
def T3r(x, y, z):          # value = i * T3r
    return Fr(-1, 2) * dot(x.p, y.v) * dot(x.v, z.v) * det3(x.c, y.c, z.c)

def V3r(a, b, c):          # full third derivative / i
    return sum(T3r(*s) for s in itertools.permutations((a, b, c)))

def T4(a, b, c, d):
    return Fr(-1, 8) * dot(a.v, c.v) * dot(b.v, d.v) * dot(cross(a.c, b.c), cross(c.c, d.c))

def V4(a, b, c, d):
    return sum(T4(*s) for s in itertools.permutations((a, b, c, d)))

def basis(comp, p):
    a, M = comp
    c = [0, 0, 0]; c[a] = 1
    v = [0] * 10; v[M] = 1
    return Leg(c, v, p)

def comps(species):
    return [(a, M) for a in range(3) for M in (VC if species == 'V' else WC)]

def neg(p):
    return [-x for x in p]

def add(p, q):
    return [x + y for x, y in zip(p, q)]

UNIT = [[1 if i == m else 0 for i in range(4)] for m in range(4)]

def vertex_poly(ext, sp_in, sp_out, q_in, q_out):
    """vertex matrix as polynomial in k: (i-power, denominator, dict exponent-tuple -> int64 matrix). rows: incoming line
    (momentum k+q_in), cols: outgoing line (incoming momentum -(k+q_out)). Integer matrices = denominator x exact Fraction matrices."""
    rc, cc = comps(sp_in), comps(sp_out)
    def to_int(Mf, den):
        Mi = np.zeros(Mf.shape, dtype=np.int64)
        for i in range(Mf.shape[0]):
            for j in range(Mf.shape[1]):
                x = Mf[i, j] * den
                assert x.denominator == 1
                Mi[i, j] = int(x)
        return Mi
    if len(ext) == 1:
        E = ext[0]
        def mat(pin, pout):
            M = np.empty((len(rc), len(cc)), dtype=object)
            for i, r in enumerate(rc):
                br = basis(r, pin)
                for j, s in enumerate(cc):
                    M[i, j] = Fr(V3r(E, br, basis(s, pout)))
            return M
        M0 = mat(q_in, neg(q_out))
        out = {(0, 0, 0, 0): to_int(M0, 2)}
        for m in range(4):
            Mm = mat(add(q_in, UNIT[m]), neg(add(q_out, UNIT[m]))) - M0
            e = [0, 0, 0, 0]; e[m] = 1
            out[tuple(e)] = to_int(Mm, 2)
        return 1, 2, {e: X for e, X in out.items() if X.any()}
    E1, E2 = ext
    M = np.empty((len(rc), len(cc)), dtype=object)
    for i, r in enumerate(rc):
        br = basis(r, [0] * 4)
        for j, s in enumerate(cc):
            M[i, j] = Fr(V4(E1, E2, br, basis(s, [0] * 4)))
    X = to_int(M, 8)
    return 0, 8, ({(0, 0, 0, 0): X} if X.any() else {})

def polymul(P, Q):
    out = {}
    for e1, A_ in P.items():
        for e2, B_ in Q.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            C = A_.dot(B_)
            assert np.abs(C).max() < 2 ** 50
            out[e] = out[e] + C if e in out else C
    return out

def laplacian(p):
    out = collections.defaultdict(Fr)
    for e, c in p.items():
        for m in range(4):
            if e[m] >= 2:
                f = list(e); f[m] -= 2
                out[tuple(f)] += c * e[m] * (e[m] - 1)
    return {e: c for e, c in out.items() if c != 0}

def evalp(p, k):
    tot = Fr(0)
    for e, c in p.items():
        t = c
        for m in range(4):
            t *= Fr(k[m]) ** e[m]
        tot += t
    return tot

def graph_J(legs, verts, lines):
    """verts: list of tuples of leg indices (cyclic order); lines[i] species of the line from vertex i to vertex i+1.
    returns J as complex Fraction pair (re, im)."""
    n = len(verts)
    P = [[0] * 4 for _ in range(n)]
    for i, vs in enumerate(verts):
        for l in vs:
            P[i] = add(P[i], legs[l].p[:4])
    q = []; acc = [0] * 4
    for i in range(n):
        acc = add(acc, P[i]); q.append(acc)
    assert all(x == 0 for x in q[-1])
    ipow = 0; poly = None; den = 1
    for i in range(n):
        ip, dn, vp = vertex_poly([legs[l] for l in verts[i]], lines[(i - 1) % n], lines[i], q[(i - 1) % n], q[i])
        ipow += ip; den *= dn
        poly = vp if poly is None else polymul(poly, vp)
    N = {e: Fr(2 ** n * int(np.trace(M)), den) for e, M in poly.items()}
    N = {e: c for e, c in N.items() if c != 0}
    m = n - 2
    L = N
    for _ in range(m):
        L = laplacian(L)
    if m >= 0:
        L2 = laplacian(L)
        assert all(sum(e) <= 0 for e in L2) and not L2, 'no eps^-1 expected'
    deg_left = max([sum(e) for e in L] + [0])
    assert deg_left <= 1
    kc = [Fr(-sum(q[e][mu] for e in range(n)), n) for mu in range(4)]
    import math
    val = Fr(1, math.factorial(n - 1)) * Fr(1, 4 ** m * math.factorial(m)) * evalp(L, kc)
    ph = [(1, 0), (0, 1), (-1, 0), (0, -1)][ipow % 4]
    return (val * ph[0], val * ph[1])

# ------------------------------------------------------------------ census by brute force
def vtype(ext_sp, lin, lout):
    f = list(ext_sp) + [lin, lout]
    nW = f.count('W')
    if nW % 2:
        return None
    return ('g3[%s]' if len(ext_sp) == 1 else 'g4[%s]') % ''.join(sorted(f))

def canon(verts, lines):
    n = len(verts); cands = []
    for refl in (False, True):
        if refl:
            vs = verts[::-1]; ls = [lines[(n - 2 - j) % n] for j in range(n)]
        else:
            vs, ls = verts, lines
        for r in range(n):
            cands.append((tuple(tuple(sorted(x)) for x in vs[r:] + vs[:r]), tuple(ls[r:] + ls[:r])))
    return min(cands)

def enumerate_graphs(species):
    nl = len(species); out = {}
    legs = list(range(nl))
    # partitions into ordered cyclic groups of size 1 or 2
    def groupings(rem):
        if not rem:
            yield []
            return
        for size in (1, 2):
            for grp in itertools.combinations(rem, size):
                rest = [x for x in rem if x not in grp]
                for g in groupings(rest):
                    yield [grp] + g
    seen = set()
    for grp in groupings(legs):
        n = len(grp)
        if n < 2:
            continue
        for perm in itertools.permutations(grp):
            for lines in itertools.product('VW', repeat=n):
                vt = []
                ok = True
                for i, g in enumerate(perm):
                    t = vtype([species[l] for l in g], lines[(i - 1) % n], lines[i])
                    if t is None:
                        ok = False; break
                    vt.append(t)
                if not ok:
                    continue
                c = canon(list(perm), list(lines))
                if c in seen:
                    continue
                seen.add(c)
                key = (tuple(sorted(vt)), tuple(sorted('P_AA(%s)' % l for l in lines)))
                aut = 2 if (n == 2 and lines[0] == lines[1]) else 1
                out[c] = dict(key=key, verts=[list(x) for x in c[0]], lines=list(c[1]), aut=aut)
    return out

# ------------------------------------------------------------------ comparison data
CMP = json.load(open(f'{R}/compare/compare_r2_tables.json'))
A = json.load(open(f'{R}/engine_a/results_r2.json'))
B = json.load(open(f'{R}/engine_b/results_r2.json'))
M = json.load(open(f'{R}/diagrams/manifest.json'))

SETS = {
    '(A|_V, A|_V, A|_V)': ('VVV', 'c3_AVAVAV', 1),
    '(A|_V, A|_Wperp, A|_Wperp)': ('VWW', 'c3_AVAWAW', 1),
    '(A|_V)^4': ('VVVV', 'q4_AVAVAVAV', 0),
    '(A|_V)^2 (A|_Wperp)^2': ('VVWW', 'q4_AVAVAWAW', 0),
    '(A|_Wperp)^4': ('WWWW', 'q4_AWAWAWAW', 0)}

def keystr(k):
    return 'vertices ' + ','.join(k[0]) + ' ; lines ' + ','.join(k[1])

def fa_bosonic_counts(job):
    c = collections.Counter()
    for d in M['diagrams']:
        if d['feynarts_job'] == job:
            _, vt, ln = d['id'].split(':')
            ln = ln.split('#')[0]
            if all(x.startswith('g3[') or x.startswith('g4[') for x in vt.split(',')):
                c[keystr((tuple(sorted(vt.split(','))), tuple(sorted(ln.split(',')))))] += 1
    return c

def rand_vec(rng, comps_, dim):
    v = [0] * dim
    for i in comps_:
        v[i] = rng.choice([-3, -2, -1, 1, 2, 3])
    return v

def make_legs(species, omega, seed):
    rng = random.Random(seed)
    n = len(species)
    if omega == 1:
        ps = [[rng.randint(-3, 3) for _ in range(4)] for _ in range(n - 1)]
        ps.append([-sum(p[m] for p in ps) for m in range(4)])
    else:
        ps = [[0] * 4 for _ in range(n)]
    return [Leg(rand_vec(rng, range(3), 3), rand_vec(rng, VC if s == 'V' else WC, 10), ps[i]) for i, s in enumerate(species)]

# engine A quartic structures: X=(12)(34), Y=(13)(24), Z=(14)(23) colour; x,y,z Lorentz (1-indexed legs), units u
def eval_A_struct(d, legs):
    c = [l.c for l in legs]; v = [l.v for l in legs]
    col = {'X': dot(c[0], c[1]) * dot(c[2], c[3]), 'Y': dot(c[0], c[2]) * dot(c[1], c[3]), 'Z': dot(c[0], c[3]) * dot(c[1], c[2])}
    lor = {'x': dot(v[0], v[1]) * dot(v[2], v[3]), 'y': dot(v[0], v[2]) * dot(v[1], v[3]), 'z': dot(v[0], v[3]) * dot(v[1], v[2])}
    return sum(Fr(str(val)) * col[k[0]] * lor[k[1]] for k, val in d.items())

def ad(x):
    return np.array([[sum((1 if (a, b, cc) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)) else -1 if (a, b, cc) in ((0, 2, 1), (2, 1, 0), (1, 0, 2)) else 0) * x[a]
                          for a in range(3)) for b in range(3)] for cc in range(3)], dtype=object)   # (ad_x)_{cb}: [x, e_b] = sum_c eps_{a b c} x_a e_c

# engine B general basis evaluated directly on su(2): b01=([X0,X1],[X2,X3]), b03=([X0,X3],[X1,X2]) (units u),
# S = (1/24) sum_sigma tr(ad ad ad ad) (units hbar/(16pi^2) = u/(-2))
def eval_B_struct(d, legs):
    c = [l.c for l in legs]; v = [l.v for l in legs]
    b01 = dot(cross(c[0], c[1]), cross(c[2], c[3])); b03 = dot(cross(c[0], c[3]), cross(c[1], c[2]))
    S = Fr(0)
    for s in itertools.permutations(range(4)):
        Mx = ad(c[s[0]]).dot(ad(c[s[1]])).dot(ad(c[s[2]])).dot(ad(c[s[3]]))
        S += sum(Mx[i, i] for i in range(3))
    S = S / 24
    col = {'b01': b01, 'b03': b03, 'S': S / (-2)}
    lor = {'x': dot(v[0], v[1]) * dot(v[2], v[3]), 'y': dot(v[0], v[2]) * dot(v[1], v[3]), 'z': dot(v[0], v[3]) * dot(v[1], v[2])}
    tot = Fr(0)
    for k, val in d.items():
        cn, ln = k.split(' x ')
        tot += Fr(str(val)) * col[cn] * lor[ln]
    return tot


COLB = {'([X_0,X_1],[X_2,X_3])': 'b01', '([X_0,X_3],[X_1,X_2])': 'b03', '(1/24) sum_sigma tr_g(ad ad ad ad)': 'S'}
def lorB(s):
    import re
    s = s.replace(' ', '')
    m = re.findall(r'\((\w)_(\d)\.(\w)_(\d)\)', s)
    idx = tuple(sorted([(int(m[0][1]), int(m[0][3])), (int(m[1][1]), int(m[1][3]))]))
    return {((0, 1), (2, 3)): 'x', ((0, 2), (1, 3)): 'y', ((0, 3), (1, 2)): 'z'}[idx]

def rawB(d):
    out = {}
    for k, v in d.items():
        if 'eps(' in k:
            assert Fr(str(v)) == 0
            continue
        cn, ln = k.split(' x ', 1)
        kk = COLB[cn.strip()] + ' x ' + lorB(ln)
        out[kk] = out.get(kk, 0) + Fr(str(v))
    return out

def pred_A(setname, gid, legs):
    sa = A['T1']['sets'][setname]
    g = [x for x in sa['graphs'] if x['id'] == gid][0]
    le = g['log_eps']
    if isinstance(le, str):
        return Fr(le.split(' x ')[0]) * eval_A_struct(sa['tree_in_basis'], legs)
    return eval_A_struct(le['structures'], legs)

def pred_B(setname, gid, legs):
    sb = B['T1']['external_sets'][setname]
    g = [x for x in sb['graphs'] if x['id'] == gid][0]
    le = g['log_eps']
    if isinstance(le, str):
        treeB = {k: v / 8 for k, v in rawB(sb['total_structure_fit']).items()}
        return Fr(le.split(' x ')[0]) * eval_B_struct(treeB, legs)
    return eval_B_struct(rawB(le), legs)


def eps4(a, b, c, d):
    tot = 0
    for s in itertools.permutations(range(4)):
        sg = 1
        for i in range(4):
            for j in range(i + 1, 4):
                if s[i] > s[j]:
                    sg = -sg
        tot += sg * a[s[0]] * b[s[1]] * c[s[2]] * d[s[3]]
    return tot

def inv_values(setname, legs):
    """values of a basis of Spin(4) x Spin(6) x Ad-invariant multilinear forms (linear in momenta for omega=1) at this configuration"""
    c = [l.c for l in legs]; v = [l.v for l in legs]; p = [l.p[:4] for l in legs]; a = [x[:4] for x in v]
    if setname == '(A|_V, A|_V, A|_V)':
        col = det3(c[0], c[1], c[2])
        vals = [col * dot(a[i], a[j]) * dot(a[k], p[l]) for (i, j, k) in ((0, 1, 2), (0, 2, 1), (1, 2, 0)) for l in (0, 1)]
        return vals + [col * eps4(a[0], a[1], a[2], p[l]) for l in (0, 1)]
    if setname == '(A|_V, A|_Wperp, A|_Wperp)':
        col = det3(c[0], c[1], c[2])
        return [col * dot(v[1], v[2]) * dot(a[0], p[l]) for l in (0, 1)]
    cols = [dot(c[0], c[1]) * dot(c[2], c[3]), dot(c[0], c[2]) * dot(c[1], c[3]), dot(c[0], c[3]) * dot(c[1], c[2])]
    lors = [dot(v[0], v[1]) * dot(v[2], v[3]), dot(v[0], v[2]) * dot(v[1], v[3]), dot(v[0], v[3]) * dot(v[1], v[2])]
    if setname == '(A|_V)^4':
        lors = lors + [eps4(a[0], a[1], a[2], a[3])]
    if setname == '(A|_V)^2 (A|_Wperp)^2':
        lors = lors[:1]
    return [x * y for x in cols for y in lors]

def rank_fr(rows):
    Mx = [[Fr(x) for x in r] for r in rows]
    rk = 0; ncol = len(Mx[0])
    for col in range(ncol):
        piv = next((i for i in range(rk, len(Mx)) if Mx[i][col] != 0), None)
        if piv is None:
            continue
        Mx[rk], Mx[piv] = Mx[piv], Mx[rk]
        for i in range(len(Mx)):
            if i != rk and Mx[i][col] != 0:
                f = Mx[i][col] / Mx[rk][col]
                Mx[i] = [x - f * y for x, y in zip(Mx[i], Mx[rk])]
        rk += 1
    return rk, ncol

def check_ad():
    # su(2) identity used by engine A, re-derived here: tr(ad_a ad_b ad_c ad_d) = (a,b)(c,d) + (a,d)(b,c); tr(ad_a ad_b) = -2 (a,b)
    e = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    ok = True
    for a, b, c_, d in itertools.product(range(3), repeat=4):
        Mx = ad(e[a]).dot(ad(e[b])).dot(ad(e[c_])).dot(ad(e[d]))
        ok &= sum(Mx[i, i] for i in range(3)) == (a == b) * (c_ == d) + (a == d) * (b == c_)
    for a, b in itertools.product(range(3), repeat=2):
        Mx = ad(e[a]).dot(ad(e[b]))
        ok &= sum(Mx[i, i] for i in range(3)) == -2 * (a == b)
    return bool(ok)

out = {'method': __doc__, 'su2_identities_ok': check_ad(), 'sets': {}}
DIMS = {'(A|_V, A|_V, A|_V)': 8, '(A|_V, A|_Wperp, A|_Wperp)': 2, '(A|_V)^4': 12, '(A|_V)^2 (A|_Wperp)^2': 3, '(A|_Wperp)^4': 9}
all_ok = True
for setname, (species, job, omega) in SETS.items():
    graphs = enumerate_graphs(species)
    bykey = collections.defaultdict(list)
    for c, g in graphs.items():
        bykey[keystr(g['key'])].append(g)
    fa = fa_bosonic_counts(job)
    census = {k: dict(third=len(bykey.get(k, [])), feynarts=fa.get(k, 0)) for k in sorted(set(bykey) | set(fa))}
    census_ok = all(v['third'] == v['feynarts'] for v in census.values())
    rows_cmp = {r['feynarts_key']: r for r in CMP['T1'][setname]['rows']}
    res = {}
    tree_check = []; inv_rows = []
    for seed in range(101, 101 + DIMS[setname] + 2):
        legs = make_legs(species, omega, seed)
        inv_rows.append(inv_values(setname, legs))
        if omega == 1:
            tr = V3r(*legs)        # tree = i * tr
            assert tr != 0
        for k, gs in sorted(bykey.items()):
            Jre = Fr(0); Jim = Fr(0)
            for g in gs:
                jr, ji = graph_J(legs, g['verts'], g['lines'])
                Jre += jr / g['aut']; Jim += ji / g['aut']
            r = rows_cmp[k]
            ent = res.setdefault(k, dict(third=[], engine_A=r['engine_A_log_eps'], engine_B=(r['engine_B_log_eps']), match_A=True, match_B=True))
            if omega == 1:
                # sum J/|Aut| = 2 kappa x tree,  tree = i tr  ->  kappa = (Jre + i Jim)/(2 i tr) = Jim/(2 tr) - i Jre/(2 tr)
                assert Jre == 0, 'non-real multiple'
                kappa = Jim / (2 * tr)
                ent['third'].append(str(kappa))
                a_val = Fr(str(r['engine_A_log_eps'])); b_val = Fr(str(r['engine_B_value_used']))
                ent['match_A'] &= (kappa == a_val); ent['match_B'] &= (kappa == b_val)
            else:
                assert Jim == 0
                val = Jre / 2                               # = engine value in units u at this configuration
                ent['third'].append(str(val))
                a_pred = sum((pred_A(setname, gid, legs) for gid in r['engine_A_ids']), Fr(0))
                b_pred = sum((pred_B(setname, gid, legs) for gid in r['engine_B_ids']), Fr(0))
                ent.setdefault('engine_A_at_config', []).append(str(a_pred)); ent.setdefault('engine_B_at_config', []).append(str(b_pred))
                ent['match_A'] &= (val == a_pred); ent['match_B'] &= (val == b_pred)
        if omega == 0:
            # tree normalisation check against engine A's tree_in_basis (independent of the loop)
            tv = V4(*legs)
            tree_check.append(tv == eval_A_struct(A['T1']['sets'][setname]['tree_in_basis'], legs))
    rk, ncol = rank_fr(inv_rows)
    out['sets'][setname] = dict(census=census, census_matches_feynarts=census_ok, classes=res, n_configurations=len(inv_rows),
                                invariant_basis_dimension=ncol, evaluation_rank=rk, exact_identity_of_forms=(rk == ncol),
                                quartic_tree_matches_engineA_tree=(all(tree_check) if tree_check else None))
for s in out['sets'].values():
    s['all_classes_match_A'] = all(c['match_A'] for c in s['classes'].values())
    s['all_classes_match_B'] = all(c['match_B'] for c in s['classes'].values())
json.dump(out, open(f'{R}/compare/third_r2_bosonic.json', 'w'), indent=1, default=str)
for k, s in out['sets'].items():
    print(k, 'rank', s['evaluation_rank'], '/', s['invariant_basis_dimension'], 'census ok', s['census_matches_feynarts'], 'match A', s['all_classes_match_A'], 'match B', s['all_classes_match_B'], 'tree4', s['quartic_tree_matches_engineA_tree'])
    for kk, c in s['classes'].items():
        print('   ', kk, c['third'], c.get('engine_A_at_config', c['engine_A']), c.get('engine_B_at_config', c['engine_B']), c['match_A'], c['match_B'])
print('su2 identities', out['su2_identities_ok'])
