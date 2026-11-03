"""Engine B, action side: the defect cycles of the cochain I + sum_x c_x S^(1)(x) (+ c_x c_y S^(2)(x,y)).

A configuration is a list of external legs (species, numeric polarisation, exact momentum) and a list of generator
letters.  Every letter is a vertex of the cycle with two cycle slots:
  I3 (one external leg, cubic part of I), I4 (two external A legs, quartic part of I, colour keys K12, K13, K14),
  G  (quadratic part of S^(1)(x), no external leg, optional injected momentum r), G3 (cubic part T5 of S^(1)(xi) with
  one external leg).
The CE coordinate c_x of an even generator is an odd parameter placed in front of all leg parameters (in the order of
the generator list); the coefficient of (c's, odd leg parameters in leg order) is returned, i.e. the arity component in
the normalisation of Def at-names (i) (theory Lemma at-defect; Prop at-classes2 for arity two).
Colour: the endomorphism M of each letter (I3, G3: ad of its leg; I4: K12 -> ad_i ad_j - ad_j ad_i, K13 -> -ad_j ad_i,
K14 -> -ad_i ad_j; G: identity), colour = tr_g(M_n ... M_1), reduced to c_adj (c0,c1), c_adj ([c0,c1],c2), and for four
legs c_adj K12, c_adj K14, str4 (step-1 engine B r2_engine.reduce_traces).  Units u = hbar tr(ad ad)/(16 pi^2 (X,Y))."""
import itertools, math, time
from fractions import Fraction as Fr
import numpy as np
import eb_core as E
import sympy as sp

V = E.V
_TM = {}
PLACEMENT_FILTER = None   # None, 'line' (no generator letter adjacent to the special edge), 'special'


def tmpl_I3(spc):
    key = ('I3', spc)
    if key not in _TM:
        _TM[key] = E.build_template(V.CUBIC, [spc], name=f'I3_{spc}')
    return _TM[key]


def tmpl_I4():
    key = ('I4',)
    if key not in _TM:
        _TM[key] = E.build_template(V.QUARTIC, ['A', 'A'], name='I4_AA')
    return _TM[key]


def perm_parity(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def colour_of_word(Ms):
    acc = {(): Fr(1)}
    for M in Ms:
        new = {}
        for w1, c1 in M.items():
            for w2, c2 in acc.items():
                w = tuple(w1) + tuple(w2)
                new[w] = new.get(w, 0) + c1 * c2
        acc = new
    return {w: c for w, c in acc.items() if c != 0}


def reduce_traces(T, legs):
    pos = {l: i for i, l in enumerate(legs)}
    out = {}

    def add(k, v):
        out[k] = out.get(k, 0) + v
    for w, c in T.items():
        assert sorted(w) == sorted(legs) and len(set(w)) == len(w), w
        s = [pos[x] for x in w]
        if len(s) == 2:
            add('k', c)
        elif len(s) == 3:
            add('f', c * Fr(1, 2) * perm_parity(s))
        elif len(s) == 4:
            i0 = s.index(0)
            r = s[i0:] + s[:i0]
            if r[1] > r[3]:
                r = [r[0], r[3], r[2], r[1]]
            t = tuple(r)
            if t == (0, 1, 2, 3):
                add('str4', c); add('K12', c * Fr(1, 6)); add('K14', -c * Fr(1, 6))
            elif t == (0, 1, 3, 2):
                add('str4', c); add('K12', -c * Fr(1, 3)); add('K14', -c * Fr(1, 6))
            elif t == (0, 2, 1, 3):
                add('str4', c); add('K12', c * Fr(1, 6)); add('K14', c * Fr(1, 3))
            else:
                raise ValueError(t)
        else:
            raise ValueError(len(s))
    return {k: v for k, v in out.items() if v != 0}


def colour_M(kind, ext, key):
    if kind in ('I3', 'G3'):
        assert key == 'f', key
        return {(ext[0],): Fr(1)}
    if kind == 'G':
        assert key == 'k', key
        return {(): Fr(1)}
    i, j = ext
    if key == 'K12':
        return {(i, j): Fr(1), (j, i): Fr(-1)}
    if key == 'K13':
        return {(j, i): Fr(-1)}
    if key == 'K14':
        return {(i, j): Fr(-1)}
    raise ValueError(key)


class Gen:
    """a generator letter source: name, CE parity codd (1 for even generators), quadratic template (dict key->T),
    cubic template per external species (dict spc -> {key: T}) or None, numeric generator vector, injected momentum."""

    def __init__(self, name, codd, tq, gvec=None, t3=None, inj=None):
        self.name, self.codd, self.tq, self.gvec, self.t3 = name, codd, tq, gvec, t3
        self.inj = [Fr(0)] * 4 if inj is None else [Fr(x) for x in inj]


def letters_for(legs, gens):
    """all letter sets: partitions of the legs into I3, I4, G3 letters, with one letter per generator."""
    idx = list(range(len(legs)))
    out = []
    # choose which generators are cubic (G3) and which leg they take
    gchoices = []
    for g in gens:
        opts = [('G', None)]
        if g.t3:
            for i in idx:
                if legs[i]['sp'] in g.t3:
                    opts.append(('G3', i))
        gchoices.append(opts)
    for gc in itertools.product(*gchoices):
        used = [i for (k, i) in gc if k == 'G3']
        if len(set(used)) != len(used):
            continue
        rest = [i for i in idx if i not in used]
        Aidx = [i for i in rest if legs[i]['sp'] == 'A']
        # partitions of rest into I3 singles and I4 pairs of A legs
        def parts(rem):
            if not rem:
                yield []
                return
            i0 = rem[0]
            for p in parts(rem[1:]):
                yield [('I3', (i0,))] + p
            if legs[i0]['sp'] == 'A':
                for j in rem[1:]:
                    if legs[j]['sp'] == 'A':
                        r2 = [x for x in rem[1:] if x != j]
                        for p in parts(r2):
                            yield [('I4', (i0, j))] + p
        for p in parts(rest):
            lets = []
            for gi, (k, i) in enumerate(gc):
                lets.append(('G' if k == 'G' else 'G3', gi, () if i is None else (i,)))
            for (k, ids) in p:
                lets.append((k, None, ids))
            out.append(lets)
    return out


def word_evaluate(word, legs, gens, mlist=(0, 1), stats=None):
    """sum over special positions and colour keys of one word (ordered list of letters).
    Returns {colour basis key: {m: CLog}} without the 1/(2n) and sign factors (applied by the caller per pos)."""
    n = len(word)
    # momenta
    curs = [[Fr(0)] * 4]
    for (kind, gi, ids) in word:
        inj = [Fr(0)] * 4
        for i in ids:
            inj = [inj[m] + Fr(legs[i]['p'][m]) for m in range(4)]
        if gi is not None:
            inj = [inj[m] + gens[gi].inj[m] for m in range(4)]
        curs.append([curs[-1][m] + inj[m] for m in range(4)])
    assert all(x == 0 for x in curs[-1]), curs[-1]
    curs = curs[:n]
    # letter templates and colour keys
    ltemps = []
    for (kind, gi, ids) in word:
        if kind == 'I3':
            ltemps.append(tmpl_I3(legs[ids[0]]['sp']))
        elif kind == 'I4':
            ltemps.append(tmpl_I4())
        elif kind == 'G':
            ltemps.append(gens[gi].tq)
        else:
            ltemps.append(gens[gi].t3[legs[ids[0]]['sp']])
    keysets = [sorted(t.keys()) for t in ltemps]
    res = {}
    for keys in itertools.product(*keysets):
        Ms = [colour_M(word[kk][0], word[kk][2], keys[kk]) for kk in range(n)]
        col = reduce_traces(colour_of_word(Ms), sorted(range(len(legs))))
        if not col:
            continue
        # letter factors
        lf = []
        for kk, (kind, gi, ids) in enumerate(word):
            T = ltemps[kk][keys[kk]]
            us = [legs[i]['u'] for i in ids]
            ps = [[complex(x) for x in legs[i]['p']] for i in ids]
            gvec = gens[gi].gvec if gi is not None else None
            cin = [complex(x) for x in curs[kk]]
            cout = [complex(x) for x in curs[(kk + 1) % n]]
            lf.append(E.letter_factor(T, us, ps, cin, cout, gen=gvec))
        gpos = [kk for kk in range(n) if word[kk][1] is not None]
        for pos in range(n):
            if PLACEMENT_FILTER is not None and gpos:
                adj = any(pos == g or pos == (g + 1) % n for g in gpos)
                if (PLACEMENT_FILTER == 'line') == adj:
                    continue
            facs = []
            for kk in range(n):
                facs.append(E.edge_factor([complex(x) for x in curs[kk]], kk == pos))
                facs.append(lf[kk])
            d = sum(1 for (c, lc) in facs if lc is not None)
            pts = np.array(E.simplex_grid(d), dtype=float)
            prod = None
            for (c, lc) in facs:
                bl = E.eval_blocks(E.to_blocks(c, lc), pts)
                prod = bl if prod is None else E.bmul(prod, bl)
                if not prod:
                    break
            if not prod:
                continue
            vals = E.bstr(prod, pts.shape[0])
            if not np.any(vals != 0):
                continue
            poly = E.interpolate(vals, d)
            if not poly:
                continue
            sc = E.schwinger(poly, curs, pos, mlist)
            if stats is not None:
                stats['chains'] = stats.get('chains', 0) + 1
            for ck, cv in col.items():
                d2 = res.setdefault(pos, {}).setdefault(ck, {})
                for m, v in sc.items():
                    d2[m] = d2.get(m, E.CLog()) + v.scale((cv, 0))
    return res, curs


def evaluate(legs, gens, mlist=(0, 1), cyclic=True, stats=None, verbose=False):
    """total of all defect cycles (class C) of the configuration: {colour key: {m: CLog}}, units u (coefficient of the
    colour structure), the eps^{-m} coefficients.  cyclic=True: words with the first letter fixed (the 1/(2n) sum over
    all n! words equals (1/2) times the sum over the (n-1)! words with a fixed first letter; checked by gate)."""
    # odd parameters: CE coordinates first (generator order), then odd legs in leg order
    oddkeys = [('c', gi) for gi, g in enumerate(gens) if g.codd] + [('l', i) for i, l in enumerate(legs) if V.PAR[l['sp']]]
    rank = {k: t for t, k in enumerate(oddkeys)}
    total = {}
    for lets in letters_for(legs, gens):
        n = len(lets)
        perms = itertools.permutations(range(n)) if not cyclic else ((0,) + p for p in itertools.permutations(range(1, n)))
        for perm in perms:
            word = [lets[i] for i in perm]
            seq = []
            per_letter_odd = []
            for (kind, gi, ids) in word:
                cnt = 0
                if gi is not None and gens[gi].codd:
                    seq.append(rank[('c', gi)]); cnt += 1
                for i in ids:
                    if V.PAR[legs[i]['sp']]:
                        seq.append(rank[('l', i)]); cnt += 1
                per_letter_odd.append(cnt)
            rsign = perm_parity(seq)
            res, curs = word_evaluate(word, legs, gens, mlist, stats)
            for pos, d in res.items():
                ob = sum(per_letter_odd[:pos])
                w = Fr(1, 2 * n) * rsign * (-1) ** ob
                if cyclic:
                    w *= n
                for ck, dm in d.items():
                    for m, v in dm.items():
                        tt = total.setdefault(ck, {})
                        tt[m] = tt.get(m, E.CLog()) + v.scale((w, 0))
            if verbose:
                print('  word', [(k, gi, ids) for (k, gi, ids) in word], flush=True)
    return total
