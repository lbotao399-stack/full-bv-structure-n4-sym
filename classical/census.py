"""
census.py -- enumeration of one-loop 1PI graph classes of W(P(eps,infty), I) (and of I + delta O).

Leg types (refined by the V / W_R^perp split of V_10):
  AV, AW   : A along V, A along W_R^perp (iota_w A)
  lam      : lambda
  c        : ghost
  AsV, AsW : A^* along V, along W_R^perp
  lams, cs : lambda^*, c^*
Nonzero components of P(eps,infty) (checks.json, section 5): AV-AV, AW-AW, lam-lam, c-AsV (= AsV-c).
Vertices of I (index-free forms in classical_spec.md):
  g3  : -1/2 (dA, 1/2[A^A])        legs {AV,AV,AV}, {AV,AW,AW}
  g4  : -1/4 (1/2[A^A], 1/2[A^A])  legs {AV^4}, {AV,AV,AW,AW}, {AW^4}
  y   : 1/2 (lam, (dslash_A - dslash) lam)   legs {lam,lam,AV}, {lam,lam,AW}
  gh  : ([A,c], A^*)               legs {AV,c,AsV}, {AW,c,AsW}
  cl  : ([c,lam], lam^*)           legs {c,lam,lams}
  cc  : 1/2([c,c], c^*)            legs {c,c,cs}
Operator insertion (Step 1(b)): O = (iota_{w1}A, iota_{w2}A), legs {AW,AW}, colour singlet, dimension 2.
Superficial degree: omega = 4 - sum_ext d (two-point), omega = 2 - sum_ext d (one O),
  d: AV=AW=1, lam=3/2, c=0, AsV=AsW=3, lams=5/2, cs=4 (all vertices of I have dimension 4).
"""
import itertools, json
from fractions import Fraction

VERTS = {
    'g3[VVV]': ('AV', 'AV', 'AV'), 'g3[VWW]': ('AV', 'AW', 'AW'),
    'g4[VVVV]': ('AV', 'AV', 'AV', 'AV'), 'g4[VVWW]': ('AV', 'AV', 'AW', 'AW'), 'g4[WWWW]': ('AW', 'AW', 'AW', 'AW'),
    'y[V]': ('lam', 'lam', 'AV'), 'y[W]': ('lam', 'lam', 'AW'),
    'gh[V]': ('AV', 'c', 'AsV'), 'gh[W]': ('AW', 'c', 'AsW'),
    'cl': ('c', 'lam', 'lams'), 'cc': ('c', 'c', 'cs'),
}
OVERT = {'O': ('AW', 'AW')}
PROP = {('AV', 'AV'): 'P_AA(V)', ('AW', 'AW'): 'P_AA(W)', ('lam', 'lam'): 'P_lamlam',
        ('c', 'AsV'): 'P_cA*', ('AsV', 'c'): 'P_A*c'}
DIM = {'AV': Fraction(1), 'AW': Fraction(1), 'lam': Fraction(3, 2), 'c': Fraction(0), 'AsV': Fraction(3), 'AsW': Fraction(3),
       'lams': Fraction(5, 2), 'cs': Fraction(4)}
PAR = {'AV': 0, 'AW': 0, 'lam': 1, 'c': 1, 'AsV': 1, 'AsW': 1, 'lams': 0, 'cs': 0}


def vertex_choices(vset):
    """(vertex name, in-leg index, out-leg index, external legs tuple)"""
    out = []
    for name, legs in vset.items():
        n = len(legs)
        for i in range(n):
            for j in range(n):
                if i == j:
                    continue
                ext = tuple(sorted(legs[k] for k in range(n) if k not in (i, j)))
                out.append((name, legs[i], legs[j], ext))
    return out


def canon(cyc):
    """canonical form of a cycle [(name, in, out, ext), ...] under rotation and reflection."""
    n = len(cyc)
    forms = []
    for r in range(n):
        rot = cyc[r:] + cyc[:r]
        forms.append(tuple(rot))
        refl = [(v[0], v[2], v[1], v[3]) for v in reversed(rot)]
        forms.append(tuple(refl))
    return min(forms)


def aut(cyc):
    """automorphisms fixing every external leg: identity, and the reflection fixing all vertices
    (possible only for n<=2: n=1 loop flip, n=2 edge swap), if in-leg ~ out-leg at every vertex."""
    n = len(cyc)
    if n <= 2 and all(v[1] == v[2] for v in cyc):
        return 2
    # n >= 3: a reflection fixes at most two vertices; with >= 3 vertices carrying distinguishable
    # roles (external legs or the O insertion) the only automorphism is the identity
    if n >= 3:
        # reflection fixing vertex k maps v_{k+i} -> v_{k-i}; allowed only if mapped vertices carry
        # no external legs and are the same type
        cnt = 1
        for k in range(n):
            ok = True
            for i in range(1, n):
                a, b = cyc[(k + i) % n], cyc[(k - i) % n]
                if (k + i) % n != (k - i) % n and (a[3] or b[3] or a[0] != b[0]):
                    ok = False
            if ok and all(v[1] == v[2] for v in cyc):
                cnt = 2
        return cnt
    return 1


def enumerate_cycles(nmax, nO, ext_total):
    """one-loop cycles with n <= nmax vertices, exactly nO O-vertices, total external legs = ext_total."""
    chI = vertex_choices(VERTS)
    chO = vertex_choices(OVERT)
    found = {}
    for n in range(1, nmax + 1):
        for pattern in itertools.product([0, 1], repeat=n):
            if sum(pattern) != nO:
                continue
            pools = [chO if x else chI for x in pattern]
            for cyc in itertools.product(*pools):
                if sum(len(v[3]) for v in cyc) != ext_total:
                    continue
                ok = True
                for i in range(n):
                    a, b = cyc[i], cyc[(i + 1) % n]
                    if (a[2], b[1]) not in PROP:
                        ok = False
                        break
                if not ok:
                    continue
                key = canon(list(cyc))
                found[key] = key
    return list(found.values())


def describe(cyc, insertion=False):
    n = len(cyc)
    ext = tuple(sorted(sum((list(v[3]) for v in cyc), [])))
    lines = [PROP[(cyc[i][2], cyc[(i + 1) % n][1])] for i in range(n)]
    om = (2 if insertion else 4) - sum(DIM[x] for x in ext)
    return {'external': ext, 'vertices': [v[0] for v in cyc],
            'internal_lines': lines, 'Aut': aut(cyc), 'omega': str(om)}


def notes(d, insertion):
    ext = d['external']
    nb = []
    vs = d['vertices']
    if insertion and len(ext) == 1:
        nb.append('vanishes: colour factor f(X,i,i) of a singlet insertion (contraction of the structure constants with the inverse pairing) = 0')
    if not insertion and set(ext) == {'AV', 'AW'}:
        nb.append('vanishes: Spin(4) x Spin(6) invariance (regulator H = Laplacian (x) id preserves it); the fermion trace has one W_R^perp Clifford factor')
    if 'cl' in vs or 'cc' in vs:
        nb.append('cannot occur: no propagator ends on lambda^* or c^*')
    return nb


def run():
    res = {}
    # (i) two-point
    two = enumerate_cycles(2, 0, 2)
    rows = []
    for cyc in two:
        d = describe(cyc)
        d['notes'] = notes(d, False)
        rows.append(d)
    rows.sort(key=lambda d: (d['external'], d['vertices']))
    res['two_point'] = rows
    # (ii) one O insertion, external legs 0, 1, 2 (cycles up to 3 vertices suffice: O 2-valent + <= 2 cubic)
    ins = []
    for k in (0, 1, 2):
        for cyc in enumerate_cycles(3, 1, k):
            d = describe(cyc, True)
            d['notes'] = notes(d, True)
            ins.append(d)
    ins.sort(key=lambda d: (len(d['external']), d['external'], d['vertices']))
    res['one_O_insertion_all'] = ins
    res['one_O_insertion_omega>=0'] = [d for d in ins if Fraction(d['omega']) >= 0]
    # all two-point external pairs of total ghost degree 0 and even parity: which have graphs
    deg = {'AV': 0, 'AW': 0, 'lam': 0, 'c': 1, 'AsV': -1, 'AsW': -1, 'lams': -1, 'cs': -2}   # coordinate degrees
    allpairs = []
    for a, b in itertools.combinations_with_replacement(sorted(DIM), 2):
        if deg[a] + deg[b] == 0 and (PAR[a] + PAR[b]) % 2 == 0:
            have = [r for r in rows if tuple(sorted((a, b))) == r['external']]
            allpairs.append({'pair': (a, b), 'omega': str(4 - DIM[a] - DIM[b]), 'n_graph_classes': len(have)})
    res['degree0_even_two_point_sectors'] = allpairs
    return res


if __name__ == '__main__':
    r = run()
    print(json.dumps(r, indent=1, default=str))
