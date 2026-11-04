"""Candidate two-leg local functionals bilinear in two constant spinors (xi1, xi2 in Sigma), for the arity-two sectors
of odd pairs, in any realisation (R_M: Sigma -> Sigma^*, RB_M: Sigma^* -> Sigma, MP: pairing (s, s*) = s^T MP s*).
A candidate density is a full contraction (Kronecker deltas over the ten-dimensional index; momenta and A|_V carry
only V components) of
  bosonic sectors:  the bilinear (xi1, Gamma_{M1..Mk} xi2), k = 1, 3, with the vectors p (momentum of the first slot
                     leg, D = 1 or 2 factors), a (A or A* polarisation of slot 0), b (of slot 1);
  spinor sectors:   Fierz products (xi_i, X l1)(xi_j, Y l2), {i,j} = {1,2}, {l1,l2} = the two spinor legs,
                     X, Y gamma strings of the parity required by the pairing, contracted with p and each other.
Densities are passed to eb_qb.dens (vertices.vertex convention: all slot assignments, odd-leg permutation sign), so
the same candidates serve both routes; completeness is checked against the representation count of each sector."""
import itertools
import numpy as np


class Real:
    def __init__(self, R, RB, MP):
        self.R = [np.asarray(x, complex) for x in R]
        self.RB = [np.asarray(x, complex) for x in RB]
        self.MP = np.asarray(MP, complex)

    def string_on(self, idx, vec, start):
        """apply rho(e_{idx[-1]}) first ... rho(e_{idx[0]}) last to vec in space start ('S' Sigma or 'Ss' Sigma*);
        returns (vector, final space)."""
        v, sp_ = vec, start
        for M in reversed(idx):
            if sp_ == 'S':
                v = self.R[M] @ v; sp_ = 'Ss'
            else:
                v = self.RB[M] @ v; sp_ = 'S'
        return v, sp_

    def pair(self, s1, sp1, s2, sp2):
        """(s1, s2) with one in Sigma and the other in Sigma*."""
        if sp1 == 'S' and sp2 == 'Ss':
            return s1 @ self.MP @ s2
        if sp1 == 'Ss' and sp2 == 'S':
            return s2 @ self.MP @ s1
        raise ValueError((sp1, sp2))


def matchings(items):
    if not items:
        yield []
        return
    a = items[0]
    for i in range(1, len(items)):
        rest = items[1:i] + items[i + 1:]
        for m in matchings(rest):
            yield [(a, items[i])] + m


def _vec10(x):
    v = np.zeros(10, complex)
    x = np.asarray(x, complex)
    v[:len(x)] = x
    return v


def index_sum(slots, vectors, fn):
    """sum over 10d indices of the gamma slots: slots = list of labels, each 'g0','g1',... (gamma positions) or vector
    names; pairing = list of pairs over slot positions.  Helper below."""
    raise NotImplementedError


class Cand:
    def __init__(self, label, slots, L):
        self.label, self.slots, self.L = label, slots, L


def _contract_eval(pattern, gam_lists, vecs, evalfn):
    """pattern: matching of the abstract index carriers; gam_lists: list of gamma-position names; vecs: name -> 10-vector.
    evalfn(assign) -> value for an assignment {gamma name: index}; sums the contraction exactly."""
    gnames = [g for gl in gam_lists for g in gl]
    # pairs: (x, y); gamma-gamma pairs summed over M, gamma-vector pairs fix M by weighting with vector components
    total = 0j
    gg = [(x, y) for (x, y) in pattern if x in gnames and y in gnames]
    gv = [(x, y) if x in gnames else (y, x) for (x, y) in pattern if (x in gnames) != (y in gnames)]
    vv = [(x, y) for (x, y) in pattern if x not in gnames and y not in gnames]
    pre = 1 + 0j
    for (x, y) in vv:
        pre *= vecs[x] @ vecs[y]
    if pre == 0:
        return 0j
    ranges = []
    for (g, v) in gv:
        nz = [M for M in range(10) if vecs[v][M] != 0]
        ranges.append(nz)
    for choice in itertools.product(*ranges):
        w = 1 + 0j
        assign = {}
        for (g, v), M in zip(gv, choice):
            assign[g] = M; w *= vecs[v][M]
        for Ms in itertools.product(range(10), repeat=len(gg)):
            a2 = dict(assign)
            for (x, y), M in zip(gg, Ms):
                a2[x] = M; a2[y] = M
            val = evalfn(a2)
            if val != 0:
                total += w * val
    return pre * total


def boson_cands(rl, kind, nmom, maxk=3):
    """kind: 'AAs' (slots A, As), 'cc' (c, cs), 'Acs' (A, cs), 'cAs' (c, As), 'cA' etc.; nmom momenta of slot 0 leg."""
    slots = {'AAs': ('A', 'As'), 'ccs': ('c', 'cs'), 'Acs': ('A', 'cs'), 'AA': ('A', 'A'), 'cAs': ('c', 'As'),
             'Ac': ('A', 'c')}[kind]
    vnames = ['p%d' % i for i in range(nmom)]
    if slots[0] in ('A', 'As'):
        vnames.append('a')
    if slots[1] in ('A', 'As'):
        vnames.append('b')
    out = []
    for k in range(1, maxk + 1, 2):
        g = ['g%d' % i for i in range(k)]
        carriers = g + vnames
        if len(carriers) % 2:
            continue
        seen = set()
        for m in matchings(carriers):
            key = tuple(sorted(tuple(sorted(x)) for x in m))
            if key in seen:
                continue
            seen.add(key)
            # skip gamma-gamma contractions inside the string of adjacent slots (they reduce to shorter strings)
            if any(x.startswith('g') and y.startswith('g') for (x, y) in m):
                continue

            def L(s, ctx, m=m, g=g, slots=slots, nmom=nmom):
                vecs = {}
                p0 = _vec10(np.array([complex(x) for x in s[0]['p']])[:4])
                for i in range(nmom):
                    vecs['p%d' % i] = p0
                if slots[0] in ('A', 'As'):
                    vecs['a'] = _vec10(s[0]['u'])
                if slots[1] in ('A', 'As'):
                    vecs['b'] = _vec10(s[1]['u'])
                xi1, xi2 = ctx['xi1'], ctx['xi2']

                def ev(asg):
                    idx = [asg[x] for x in g]
                    v, spc = rl.string_on(idx, xi2, 'S')
                    return rl.pair(xi1, 'S', v, spc)
                val = _contract_eval(m, [g], vecs, ev)
                for t, leg in zip(slots, s):
                    if t in ('c', 'cs'):
                        val *= leg['u'][0]
                return val
            out.append(Cand('%s|k%d|%s' % (kind, k, key), slots, L))
    return out


def spinor_cands(rl, kind, nmom, maxk=3):
    """kind: 'llS' (slots lam, lams), 'lsls' (lams, lams); Fierz products (xi_i, X l_a)(xi_j, Y l_b)."""
    slots = {'llS': ('lam', 'lams'), 'lsls': ('lams', 'lams')}[kind]
    spaces = {'lam': 'S', 'lams': 'Ss'}
    vnames = ['p%d' % i for i in range(nmom)]
    out = []
    for (i1, i2) in ((1, 2), (2, 1)):
        for k1 in range(0, maxk + 1):
            for k2 in range(0, maxk + 1):
                # parity: (xi, X l) needs X: space(l) -> Sigma^*; from S: odd length, from Ss: even length
                if (k1 % 2) != (1 if spaces[slots[0]] == 'S' else 0):
                    continue
                if (k2 % 2) != (1 if spaces[slots[1]] == 'S' else 0):
                    continue
                g1 = ['x%d' % i for i in range(k1)]; g2 = ['y%d' % i for i in range(k2)]
                carriers = g1 + g2 + vnames
                if len(carriers) % 2 or len(carriers) > 6:
                    continue
                seen = set()
                for m in matchings(carriers):
                    if any((x in g1 and y in g1) or (x in g2 and y in g2) for (x, y) in m):
                        continue
                    key = tuple(sorted(tuple(sorted(x)) for x in m))
                    if key in seen:
                        continue
                    seen.add(key)

                    def L(s, ctx, m=m, g1=g1, g2=g2, i1=i1, i2=i2, nmom=nmom):
                        vecs = {}
                        p0 = _vec10(np.array([complex(x) for x in s[0]['p']])[:4])
                        for i in range(nmom):
                            vecs['p%d' % i] = p0
                        xa, xb = ctx['xi%d' % i1], ctx['xi%d' % i2]
                        la, lb = s[0]['u'], s[1]['u']

                        def ev(asg):
                            va, sa = rl.string_on([asg[x] for x in g1], la, spaces[slots[0]])
                            vb, sb = rl.string_on([asg[x] for x in g2], lb, spaces[slots[1]])
                            return rl.pair(xa, 'S', va, sa) * rl.pair(xb, 'S', vb, sb)
                        return _contract_eval(m, [g1, g2], vecs, ev)
                    out.append(Cand('%s|%d%d|k%d,%d|%s' % (kind, i1, i2, k1, k2, key), slots, L))
    return out


def spinor_cands_b(rl, kind, nmom, maxk=3):
    """second Fierz ordering: (xi1, X xi2)(l1, Y l2), X odd string (Sigma -> Sigma^*), Y with the parity required by the
    leg spaces; contracted with p and each other."""
    slots = {'llS': ('lam', 'lams'), 'lsls': ('lams', 'lams')}[kind]
    spaces = {'lam': 'S', 'lams': 'Ss'}
    vnames = ['p%d' % i for i in range(nmom)]
    out = []
    for k1 in (1, 3):
        for k2 in range(0, maxk + 1):
            # (l1, Y l2): l2 in space s2; Y l2 must land in the dual of space(l1)
            s1, s2 = spaces[slots[0]], spaces[slots[1]]
            need = 'Ss' if s1 == 'S' else 'S'
            land = s2 if k2 % 2 == 0 else ('Ss' if s2 == 'S' else 'S')
            if land != need:
                continue
            g1 = ['x%d' % i for i in range(k1)]; g2 = ['y%d' % i for i in range(k2)]
            carriers = g1 + g2 + vnames
            if len(carriers) % 2 or len(carriers) > 6:
                continue
            seen = set()
            for m in matchings(carriers):
                if any((x in g1 and y in g1) or (x in g2 and y in g2) for (x, y) in m):
                    continue
                key = tuple(sorted(tuple(sorted(x)) for x in m))
                if key in seen:
                    continue
                seen.add(key)

                def L(s, ctx, m=m, g1=g1, g2=g2, nmom=nmom, s1=s1, s2=s2):
                    vecs = {}
                    p0 = _vec10(np.array([complex(x) for x in s[0]['p']])[:4])
                    for i in range(nmom):
                        vecs['p%d' % i] = p0

                    def ev(asg):
                        va, sa = rl.string_on([asg[x] for x in g1], ctx['xi2'], 'S')
                        vb, sb = rl.string_on([asg[x] for x in g2], s[1]['u'], s2)
                        return rl.pair(ctx['xi1'], 'S', va, sa) * rl.pair(s[0]['u'], s1, vb, sb)
                    return _contract_eval(m, [g1, g2], vecs, ev)
                out.append(Cand('%s|b|k%d,%d|%s' % (kind, k1, k2, key), slots, L))
    return out


def lin_cands(rl, slots, nmom, vecs_extra=('b',), maxk=3, ckey='xi', cspace='S'):
    """candidates linear in one constant spinor xi (ctx['xi']) and in the vectors vecs_extra (ctx[name], 4-vectors):
    (xi, X l) with l the spinor leg (lam in Sigma or lam* in Sigma^*) times full contractions of the gamma indices of X,
    the momenta p (nmom factors, momentum of slot 0), the vectors of vecs_extra and the vector legs (A, A*)."""
    spaces = {'lam': 'S', 'lams': 'Ss'}
    sidx = [i for i, t in enumerate(slots) if t in spaces][0]
    oidx = 1 - sidx
    vnames = ['p%d' % i for i in range(nmom)] + list(vecs_extra)
    if slots[oidx] in ('A', 'As'):
        vnames.append('a')
    out = []
    dual = 'Ss' if cspace == 'S' else 'S'
    par_needed = 0 if spaces[slots[sidx]] == dual else 1
    for k in range(0, maxk + 1):
        if k % 2 != par_needed:
            continue
        g = ['g%d' % i for i in range(k)]
        carriers = g + vnames
        if len(carriers) % 2:
            continue
        seen = set()
        for m in matchings(carriers):
            if any(x.startswith('g') and y.startswith('g') for (x, y) in m):
                continue
            key = tuple(sorted(tuple(sorted(x)) for x in m))
            if key in seen:
                continue
            seen.add(key)

            def L(s, ctx, m=m, g=g, nmom=nmom):
                vecs = {}
                p0 = _vec10(np.array([complex(x) for x in s[0]['p']])[:4])
                for i in range(nmom):
                    vecs['p%d' % i] = p0
                for nm in vecs_extra:
                    vecs[nm] = _vec10(ctx[nm])
                if slots[oidx] in ('A', 'As'):
                    vecs['a'] = _vec10(s[oidx]['u'])
                sp_ = spaces[slots[sidx]]

                def ev(asg):
                    v, spc = rl.string_on([asg[x] for x in g], s[sidx]['u'], sp_)
                    return rl.pair(ctx[ckey], cspace, v, spc)
                val = _contract_eval(m, [g], vecs, ev)
                if slots[oidx] in ('c', 'cs'):
                    val *= s[oidx]['u'][0]
                return val
            out.append(Cand('lin%s|%s|k%d|%s' % (''.join(slots), ''.join(vecs_extra), k, key), slots, L))
    return out


def _string_parity(start, end):
    """parity of the gamma-string length mapping space start to space end"""
    return 0 if start == end else 1


def bil2_cands(rl, slots, nmom, vecs_extra=(), sp1=('xi', 'S'), sp2=('eta', 'Ss'), maxk=4, maxc=6):
    """candidates bilinear in two constant spinors of arbitrary spaces (sp1, sp2): bosonic legs: (s1, X s2) times
    contractions with p (nmom), vecs_extra, vector legs; spinor legs (both slots spinors): Fierz products
    (s1, X l_a)(s2, Y l_b) for both assignments and (s1, X s2)(l_0, Y l_1)."""
    spaces = {'lam': 'S', 'lams': 'Ss'}
    dual = {'S': 'Ss', 'Ss': 'S'}
    vn = ['p%d' % i for i in range(nmom)] + list(vecs_extra)
    vec_slots = [i for i, t in enumerate(slots) if t in ('A', 'As')]
    sc_slots = [i for i, t in enumerate(slots) if t in ('c', 'cs')]
    spn_slots = [i for i, t in enumerate(slots) if t in spaces]
    vnames = vn + ['a%d' % i for i in vec_slots]
    out = []

    def setvecs(s, ctx):
        vecs = {}
        p0 = _vec10(np.array([complex(x) for x in s[0]['p']])[:4])
        for i in range(nmom):
            vecs['p%d' % i] = p0
        for nm in vecs_extra:
            vecs[nm] = _vec10(ctx[nm])
        for i in vec_slots:
            vecs['a%d' % i] = _vec10(s[i]['u'])
        return vecs

    def scal(s):
        v = 1
        for i in sc_slots:
            v *= s[i]['u'][0]
        return v

    def gen_matchings(groups, extra):
        carriers = [g for gl in groups for g in gl] + extra
        if len(carriers) % 2 or len(carriers) > maxc + 2:
            return []
        res, seen = [], set()
        for m in matchings(carriers):
            if any(any(x in gl and y in gl for gl in groups) for (x, y) in m):
                continue
            key = tuple(sorted(tuple(sorted(x)) for x in m))
            if key not in seen:
                seen.add(key); res.append(m)
        return res
    # (s1, X s2) * (bosonic or spinor-leg bilinear)
    k1par = _string_parity(sp2[1], dual[sp1[1]])
    for k1 in range(0, maxk + 1):
        if k1 % 2 != k1par:
            continue
        g1 = ['x%d' % i for i in range(k1)]
        if not spn_slots:
            for m in gen_matchings([g1], vnames):
                def L(s, ctx, m=m, g1=g1):
                    vecs = setvecs(s, ctx)
                    def ev(asg):
                        v, spc = rl.string_on([asg[x] for x in g1], ctx[sp2[0]], sp2[1])
                        return rl.pair(ctx[sp1[0]], sp1[1], v, spc)
                    return _contract_eval(m, [g1], vecs, ev) * scal(s)
                out.append(Cand('b2|%s|k%d|%s' % (''.join(slots), k1, m), slots, L))
        elif len(spn_slots) == 2:
            i0, i1 = spn_slots
            s0, s1_ = spaces[slots[i0]], spaces[slots[i1]]
            k2par = _string_parity(s1_, dual[s0])
            for k2 in range(0, 3):
                if k2 % 2 != k2par:
                    continue
                g2 = ['y%d' % i for i in range(k2)]
                for m in gen_matchings([g1, g2], vnames):
                    def L(s, ctx, m=m, g1=g1, g2=g2, i0=i0, i1=i1, s0=s0, s1_=s1_):
                        vecs = setvecs(s, ctx)
                        def ev(asg):
                            va, sa = rl.string_on([asg[x] for x in g1], ctx[sp2[0]], sp2[1])
                            vb, sb = rl.string_on([asg[x] for x in g2], s[i1]['u'], s1_)
                            return rl.pair(ctx[sp1[0]], sp1[1], va, sa) * rl.pair(s[i0]['u'], s0, vb, sb)
                        return _contract_eval(m, [g1, g2], vecs, ev) * scal(s)
                    out.append(Cand('b2|%s|ss|k%d,%d|%s' % (''.join(slots), k1, k2, m), slots, L))
    # Fierz products (s1, X l_a)(s2, Y l_b)
    if len(spn_slots) == 2:
        for (ia, ib) in (tuple(spn_slots), tuple(reversed(spn_slots))):
            la, lb = spaces[slots[ia]], spaces[slots[ib]]
            pa = _string_parity(la, dual[sp1[1]]); pb = _string_parity(lb, dual[sp2[1]])
            for k1 in range(0, 4):
                for k2 in range(0, 4):
                    if k1 % 2 != pa or k2 % 2 != pb:
                        continue
                    g1 = ['x%d' % i for i in range(k1)]; g2 = ['y%d' % i for i in range(k2)]
                    for m in gen_matchings([g1, g2], vnames):
                        def L(s, ctx, m=m, g1=g1, g2=g2, ia=ia, ib=ib, la=la, lb=lb):
                            vecs = setvecs(s, ctx)
                            def ev(asg):
                                va, sa = rl.string_on([asg[x] for x in g1], s[ia]['u'], la)
                                vb, sb = rl.string_on([asg[x] for x in g2], s[ib]['u'], lb)
                                return rl.pair(ctx[sp1[0]], sp1[1], va, sa) * rl.pair(ctx[sp2[0]], sp2[1], vb, sb)
                            return _contract_eval(m, [g1, g2], vecs, ev)
                        out.append(Cand('b2|%s|F%d%d|k%d,%d|%s' % (''.join(slots), ia, ib, k1, k2, m), slots, L))
    return out


LEVI4 = None


def _levi():
    global LEVI4
    if LEVI4 is None:
        L = np.zeros((10, 10, 10, 10))
        for perm in itertools.permutations(range(4)):
            s = 1
            p = list(perm)
            for i in range(4):
                for j in range(i + 1, 4):
                    if p[i] > p[j]:
                        s = -s
            L[perm] = s
        LEVI4 = L
    return LEVI4


def vec_cands(rl, slots, nmom, vecs_extra=('b',), maxk=3):
    """candidates with no constant spinor: bosonic legs -> full contractions (deltas, at most one epsilon of V) of the
    vectors p (nmom), vecs_extra (ctx), vector legs; two spinor legs (lam, lam) or (lam, lam*) -> (l0, X l1) with X a
    gamma string contracted with the vectors."""
    spaces = {'lam': 'S', 'lams': 'Ss'}
    vec_slots = [i for i, t in enumerate(slots) if t in ('A', 'As')]
    sc_slots = [i for i, t in enumerate(slots) if t in ('c', 'cs')]
    spn = [i for i, t in enumerate(slots) if t in spaces]
    vnames = ['p%d' % i for i in range(nmom)] + list(vecs_extra) + ['a%d' % i for i in vec_slots]
    out = []

    def setvecs(s, ctx):
        vecs = {}
        p0 = _vec10(np.array([complex(x) for x in s[0]['p']])[:4])
        for i in range(nmom):
            vecs['p%d' % i] = p0
        for nm in vecs_extra:
            vecs[nm] = _vec10(ctx[nm])
        for i in vec_slots:
            vecs['a%d' % i] = _vec10(s[i]['u'])
        return vecs

    def scal(s):
        v = 1
        for i in sc_slots:
            v *= s[i]['u'][0]
        return v
    if not spn:
        n = len(vnames)
        seen = set()
        if n % 2 == 0:
            for m in matchings(vnames):
                key = tuple(sorted(tuple(sorted(x)) for x in m))
                if key in seen:
                    continue
                seen.add(key)
                def L(s, ctx, m=m):
                    vecs = setvecs(s, ctx)
                    v = 1 + 0j
                    for (x, y) in m:
                        v *= vecs[x] @ vecs[y]
                    return v * scal(s)
                out.append(Cand('v|%s|%s' % (''.join(slots), key), slots, L))
            for quad in itertools.combinations(range(n), 4):
                rest = [vnames[i] for i in range(n) if i not in quad]
                for m in (matchings(rest) if rest else [[]]):
                    key = ('eps', tuple(vnames[i] for i in quad), tuple(sorted(tuple(sorted(x)) for x in m)))
                    if key in seen:
                        continue
                    seen.add(key)
                    def L(s, ctx, m=m, quad=quad):
                        vecs = setvecs(s, ctx)
                        v = 1 + 0j
                        for (x, y) in m:
                            v *= vecs[x] @ vecs[y]
                        vs = [vecs[vnames[i]] for i in quad]
                        return v * np.einsum('ijkl,i,j,k,l->', _levi(), *vs) * scal(s)
                    out.append(Cand('v|%s|%s' % (''.join(slots), key), slots, L))
        return out
    i0, i1 = spn
    s0, s1 = spaces[slots[i0]], spaces[slots[i1]]
    par = _string_parity(s1, 'Ss' if s0 == 'S' else 'S')
    for k in range(0, maxk + 1):
        if k % 2 != par:
            continue
        g = ['g%d' % i for i in range(k)]
        carriers = g + vnames
        if len(carriers) % 2:
            continue
        seen = set()
        for m in matchings(carriers):
            if any(x.startswith('g') and y.startswith('g') for (x, y) in m):
                continue
            key = tuple(sorted(tuple(sorted(x)) for x in m))
            if key in seen:
                continue
            seen.add(key)
            def L(s, ctx, m=m, g=g):
                vecs = setvecs(s, ctx)
                def ev(asg):
                    v, spc = rl.string_on([asg[x] for x in g], s[i1]['u'], s1)
                    return rl.pair(s[i0]['u'], s0, v, spc)
                return _contract_eval(m, [g], vecs, ev)
            out.append(Cand('v|%s|k%d|%s' % (''.join(slots), k, key), slots, L))
    return out
