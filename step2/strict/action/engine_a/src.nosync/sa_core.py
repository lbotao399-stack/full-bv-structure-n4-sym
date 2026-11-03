"""Strict route, action side, engine A: closed defect cycles (Thm at-classes, Prop at-classes2), exact.

A defect cycle: all vertices on the cycle, external legs at generic (integer) momenta, exactly one special edge
K' = K_eps - K_0 (numerator +|l+q|^2 Kn, Schwinger time in [0, eps]), every other edge P(eps, T = infinity).
Formula (Part-C engine, validated in 0 dimensions there): the special-edge part of the one-loop flow is
    (1/2) str[ K' S H (P S H)^(n-1) ],   H = Hessian of the interaction at the background,
expanded over ordered words of insertions with the special line in front.  For arity k the interaction is
I + S(x) c^x (+ S(y) c^y, or S^(2)(x,y) c^y c^x); the ghost c^x is odd for even x and commuting for odd x and stands to
the RIGHT of S(x) (Def at-names (i): {S, S(x) c} = {S, S(x)} c).  Hessian entries carry their odd parameters to the left
(r2_core convention), so a word reads (odd parameters in word order) x (number); the value is brought to the canonical
order theta_1 ... theta_n c^y c^x (legs in list order, then ghosts) by the permutation sign.
Ob^(1)(x) = coefficient of theta_1..theta_n c^x;  Ob^(2)(x,y) = (-1)^{|x|} coefficient of theta_1..theta_n c^y c^x
(the c^y c^x coefficient of 1/2{J,J}_eps is the special part of {S(x),S(y)}_eps for every parity combination).

Colour: colour-stripped blocks of the two-letter engine A (cubic (ad_Y t_in, t_out), quartic words, T5 as a cubic
vertex, quadratic S-blocks colour-diagonal).  A word whose external legs are met in the order Y_0, Y_1, ... carries the
colour factor tr_g(ad_{Y_m-1} ... ad_{Y_0}) (ea_core docstring); words are accumulated per cyclic colour word.

Momenta: line i (i = 0..n-1) carries k + q_i, q_0 = 0, q_{i+1} = q_i - P_i, block i between line i and line i+1 (mod n);
P = P0 + t P1 (t first order: the injected momentum r = t e_mu of a position factor, absorbed by the last leg).
Values: units hbar/(16 pi^2) x colour word, as {t-order: [rational (re,im), log(E/(E-1)) coefficient (re,im)]}."""
import os, sys, itertools, math, time
from fractions import Fraction
import numpy as np

SRC = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(SRC, '..', '..', '..', '..', '..'))
TLSRC = os.path.join(ROOT, 'step2', 'twoletter', 'engine_a', 'src.nosync')
CQSRC = os.path.join(ROOT, 'step2', 'composite_qsk', 'composite', 'engine_a', 'src.nosync')
QSSRC = os.path.join(ROOT, 'step2', 'composite_qsk', 'qsk', 'engine_a', 'src.nosync')
for _p in (QSSRC, CQSRC, TLSRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)
import ea_core as E          # noqa
import cq_core as Q          # noqa
import cq_gcore as G         # noqa
import qa_s1 as S1           # noqa
np.seterr(all='ignore')

NS = E.NS
SD = Q.SD
PARITY_T = E.PARITY_T
SLOT = E.SLOT
Z4 = np.zeros(4, dtype=complex)
Z4i = np.zeros(4, dtype=int)
STATS = {'words': 0, 'nonzero': 0, 'time': 0.0}


class Blk:
    """affine block h = h0 + t h1 + kj.Hj + kl.Hl on the 54-space (scaled), with bookkeeping"""
    __slots__ = ['h0', 'h1', 'Hj', 'Hl', 'scale', 'legs', 'colour', 'ghost', 'name']

    def __init__(self, h0, h1, Hj, Hl, scale, legs, colour, ghost, name):
        self.h0, self.h1, self.Hj, self.Hl = h0, h1, Hj, Hl
        self.scale, self.legs, self.colour, self.ghost, self.name = scale, legs, colour, ghost, name


def _zero4():
    return [np.zeros((NS, NS), complex) for _ in range(4)]


def leg(lt, pol, p0, p1=None):
    return {'t': lt, 'pol': np.asarray(pol, dtype=complex), 'p0': np.asarray(p0, dtype=complex),
            'p1': np.zeros(4, dtype=int) if p1 is None else np.asarray(p1, dtype=int)}


_CUBC = {}


def cubic_I(L, i):
    lt, pol = L['t'], L['pol']
    stp = E.SLOT[lt]
    spol = pol
    b0 = Q.cubic_full(lt, spol, L['p0'])
    h1 = None
    if np.any(L['p1']):
        h1 = Q.cubic_full(lt, spol, L['p1']).h0 - Q.cubic_full(lt, spol, Z4).h0
    return Blk(b0.h0, h1, b0.Hj, b0.Hl, 2, [i], [i], None, 'I1')


def quartic_I(L1, L2, i, j):
    """quartic vertex, colour word: Y_i met first, then Y_j (contributes ad_j ad_i)"""
    b = Q.quartic_full(L2['pol'], L1['pol'])
    return Blk(b.h0, None, b.Hj, b.Hl, 8, [i, j], [i, j], None, 'I2')


class Gen:
    """vertex data of a generator: quad = Aff (scaled x2) of the quadratic part at r (callable r-kind -> Aff),
    cub = callable(ltype, pol) -> Aff for T5; odd = |x| (1 for xi, eta), name."""

    def __init__(self, name, quad, cub=None, odd=0, cubtypes=('AV', 'AW', 'lamsP', 'lamsM')):
        self.name, self.quad, self.cub, self.odd, self.cubtypes = name, quad, cub, odd, cubtypes


def gen_even(name, **kw):
    return Gen(name, S1.quad_block(S1.even_terms(**kw)), None, 0)


def gen_odd(xi, name='xi', which=('T1', 'T2', 'T3', 'T4', 'T5')):
    xi = np.asarray(xi)
    qt = S1.odd_terms(xi, tuple(w for w in which if w != 'T5'))
    quad = S1.quad_block(qt) if qt else None
    cub = None
    if 'T5' in which:
        t5 = S1.odd_terms(xi, ('T5',))
        key = ('T5', name, tuple(np.asarray(xi).tolist()))
        cub = lambda lt, pol: S1.cubic_block(key, t5, lt, pol, Z4)
    return Gen(name, quad, cub, 1)


def gen_s0(eta, name='s0'):
    return Gen(name, S1.quad_block(S1.s0_term(eta)), None, 1)


def gen_from_terms(name, terms, odd):
    return Gen(name, S1.quad_block(terms), None, odd)


def S_block(gen, slotname, r0=None, r1=None):
    a = gen.quad
    blk = Blk(a.h0, None, a.Hj, a.Hl, 2, [], [], slotname if gen.odd == 0 else None, 'S2:' + slotname)
    blk.legs = []
    return blk, (Z4 if r0 is None else np.asarray(r0, complex)), (Z4i if r1 is None else np.asarray(r1, int))


def S3_block(gen, slotname, L, i):
    a = gen.cub(L['t'], L['pol'])
    return Blk(a.h0, None, a.Hj, a.Hl, 2, [i], [i], slotname if gen.odd == 0 else None, 'S3:' + slotname)


# ------------------------------------------------------------------------------------------------ polynomial cycle
def block_poly(b, qin0, qin1, qout0, qout1):
    c0 = b.h0 - sum(qin0[mu] * b.Hj[mu] for mu in range(4)) + sum(qout0[mu] * b.Hl[mu] for mu in range(4))
    c1 = None
    if b.h1 is not None or np.any(qin1) or np.any(qout1):
        c1 = (b.h1 if b.h1 is not None else 0) - sum(qin1[mu] * b.Hj[mu] for mu in range(4)) \
            + sum(qout1[mu] * b.Hl[mu] for mu in range(4))
        if not np.any(c1):
            c1 = None
    lin = [b.Hl[mu] - b.Hj[mu] for mu in range(4)]
    S = SD[:, None]
    return G.affine(S * c0, None if c1 is None else S * c1, [S * x for x in lin])


def cycle_N(blks, P0s, P1s, plain=False):
    """scalar polynomial {(k0..k3,t): complex} of str[K' S h_0 P S h_1 ... P S h_{n-1}] and the shifts (q0, q1)"""
    n = len(blks)
    q0 = [np.zeros(4, dtype=complex)]; q1 = [np.zeros(4, dtype=int)]
    for i in range(n - 1):
        q0.append(q0[-1] - P0s[i]); q1.append(q1[-1] - P1s[i])
    assert np.allclose(q0[-1] - P0s[-1], 0) and not np.any(q1[-1] - P1s[-1]), 'momentum conservation'
    T = None
    for i in range(n):
        Lp = G.line_poly('S' if (i == 0 and not plain) else 'P', q0[i], q1[i])
        T = Lp if T is None else G.pmul(T, Lp)
        if not T:
            return None, (q0, q1)
        j = (i + 1) % n
        T = G.pmul(T, block_poly(blks[i], q0[i], q1[i], q0[j], q1[j]))
        if not T:
            return None, (q0, q1)
    N = {}
    for e, M in T.items():
        v = complex(np.sum(SD * np.diag(M)))
        if v != 0:
            N[e] = np.array(v)
    return (N or None), (q0, q1)


def perm_sign(seq):
    a = list(seq); s_ = 1
    for i in range(len(a)):
        for j in range(len(a) - 1 - i):
            if a[j] > a[j + 1]:
                a[j], a[j + 1] = a[j + 1], a[j]; s_ = -s_
    return s_


def vzero():
    return [(Fraction(0), Fraction(0)) for _ in range(5)]


def vadd(a, b):
    return [(x[0] + y[0], x[1] + y[1]) for x, y in zip(a, b)]


def vscale(a, c):
    if isinstance(c, tuple):
        return [(x[0] * c[0] - x[1] * c[1], x[0] * c[1] + x[1] * c[0]) for x in a]
    return [(x[0] * c, x[1] * c) for x in a]


def viszero(a):
    return all(x[0] == 0 and x[1] == 0 for x in a)


def _rl_to_v(rl, E_, fac):
    """[rational (re,im), log-coefficient (re,im)] -> 5-vector on (1, log2, log3, log5, log7), times Fraction fac"""
    out = vzero()
    r0, r1 = rl
    if r0 is not None:
        out[0] = (Fraction(r0[0]) * fac, Fraction(r0[1]) * fac)
    if r1 is not None:
        c = (Fraction(r1[0]) * fac, Fraction(r1[1]) * fac)
        for bi, nn in enumerate(Q.log_basis(E_)):
            if nn:
                out[1 + bi] = (out[1 + bi][0] + nn * c[0], out[1 + bi][1] + nn * c[1])
    return out


def eval_word(blks, P0s, P1s, sign, what=('fin',)):
    """value of one word: {what: {t-order: 5-vector}} (units hbar/16pi^2, times 1/2, sign, 1/scales).
    what entries: 'fin', ('pole', m) (special-edge cycles), 'log', ('eps', m) (plain cycles, factor 1/(2n))"""
    plain = any(w == 'log' or (isinstance(w, tuple) and w[0] == 'eps') for w in what)
    N, qs = cycle_N(blks, P0s, P1s, plain)
    if N is None:
        return None
    E_ = len(blks)
    scale = 1
    for b in blks:
        scale *= b.scale
    fac = Fraction(sign, 2 * scale * (E_ if plain else 1))
    Sx = G.GSimplex(qs)
    out = {}
    for w in what:
        if w == 'fin':
            r = Sx.evaluate(N, 'fin', s=0)
        elif w == 'log':
            r = Sx.evaluate(N, 'log')
        elif w[0] == 'eps':
            r = Sx.evaluate(N, 'eps', m=w[1])
        else:
            r = Sx.evaluate(N, 'spec', s=0, m=w[1])
        out[w] = {tt: _rl_to_v(r[tt], E_, fac) for tt in (0, 1)}
    return out


# ------------------------------------------------------------------------------------------------ enumeration
def partitions(legs, gens, s2=None):
    """block multisets: every leg in exactly one block; one S block per generator in gens (quadratic, or T5 absorbing
    one leg); s2: optional generator-like object for an S^(2) vertex replacing the S blocks (arity two).
    Yields lists of block descriptors ('I1',i) ('I2',i,j) ('S2',g) ('S3',g,i) ('X2',)"""
    n = len(legs)
    Atypes = ('AV', 'AW')

    def rec(i, rem, used3):
        if i == n:
            yield []
            return
        if i in rem:
            yield from rec(i + 1, rem, used3)
            return
        for tail in rec(i + 1, rem, used3):
            yield [('I1', i)] + tail
        if legs[i]['t'] in Atypes:
            for j in range(i + 1, n):
                if j not in rem and legs[j]['t'] in Atypes:
                    for tail in rec(i + 1, rem | {j}, used3):
                        yield [('I2', i, j)] + tail
        for g, gen in enumerate(gens):
            if gen.cub is not None and g not in used3 and legs[i]['t'] in gen.cubtypes:
                for tail in rec(i + 1, rem, used3 | {g}):
                    yield [('S3', g, i)] + tail

    if s2 is not None:
        for st in rec(0, frozenset(), frozenset(range(len(gens)))):
            yield st + [('X2',)]
        return
    for st in rec(0, frozenset(), frozenset()):
        used = {b[1] for b in st if b[0] == 'S3'}
        extra = [('S2', g) for g in range(len(gens)) if g not in used]
        if any(gens[g].quad is None for g in range(len(gens)) if g not in used):
            continue
        yield st + extra


def cycle_sum(legs, gens, slotnames=None, r0=None, r1=None, s2=None, what=('fin',), split=False, colour='auto'):
    """sum over all words (closed cycles, special line in front, 1/2) of the given legs with one S block per
    generator (gens) or one S^(2) block (s2 = Gen).  r0, r1: injected momentum at the S block of gens[0] (or at the
    S^(2) block).  Returns {colour-word: {what: {t: 5-vector}}} with colour words as tuples of leg indices in cycle
    order (canonicalised by rotation)."""
    t0 = time.time()
    nl = len(legs)
    slotnames = slotnames or ['g%d' % g for g in range(len(gens))]
    odd_legs = [i for i in range(nl) if PARITY_T[legs[i]['t']] == 1]
    ghosts = [slotnames[g] for g in range(len(gens)) if gens[g].odd == 0]
    if s2 is not None:
        ghosts = []                      # S^(2) of two odd generators: commuting ghosts
    # canonical order of odd symbols: legs, then ghosts in reverse list order (c^y c^x)
    canon = [('L', i) for i in odd_legs] + [('G', gname) for gname in reversed(ghosts)]
    rank = {s: k for k, s in enumerate(canon)}
    res = {}
    for st in partitions(legs, gens, s2):
        nb = len(st)
        items = []
        for b in st:
            if b[0] == 'I1':
                items.append([('I1', b[1])])
            elif b[0] == 'I2':
                items.append([('I2', b[1], b[2]), ('I2', b[2], b[1])])
            else:
                items.append([b])
        for perm in itertools.permutations(range(nb)):
            seqs = [items[k] for k in perm]
            for choice in itertools.product(*seqs):
                blks, P0s, P1s, syms, cw = [], [], [], [], []
                for d in choice:
                    if d[0] == 'I1':
                        L = legs[d[1]]
                        blks.append(cubic_I(L, d[1])); P0s.append(L['p0']); P1s.append(L['p1'])
                        if PARITY_T[L['t']]:
                            syms.append(('L', d[1]))
                        cw.append(d[1])
                    elif d[0] == 'I2':
                        L1, L2 = legs[d[1]], legs[d[2]]
                        blks.append(quartic_I(L1, L2, d[1], d[2])); P0s.append(L1['p0'] + L2['p0'])
                        P1s.append(L1['p1'] + L2['p1']); cw += [d[1], d[2]]
                    elif d[0] == 'S2':
                        g = d[1]
                        blk, rr0, rr1 = S_block(gens[g], slotnames[g], r0 if g == 0 else None, r1 if g == 0 else None)
                        blks.append(blk); P0s.append(rr0); P1s.append(rr1)
                        if gens[g].odd == 0:
                            syms.append(('G', slotnames[g]))
                    elif d[0] == 'S3':
                        g, i = d[1], d[2]
                        L = legs[i]
                        blks.append(S3_block(gens[g], slotnames[g], L, i))
                        rr0 = L['p0'] + (r0 if (g == 0 and r0 is not None) else 0)
                        rr1 = L['p1'] + (r1 if (g == 0 and r1 is not None) else 0)
                        P0s.append(rr0); P1s.append(rr1)
                        if gens[g].odd == 0:
                            syms.append(('G', slotnames[g]))
                        cw.append(i)
                    elif d[0] == 'X2':
                        blk, rr0, rr1 = S_block(s2, 'X2', r0, r1)
                        blks.append(blk); P0s.append(rr0); P1s.append(rr1)
                STATS['words'] += 1
                sign = perm_sign([rank[s] for s in syms])
                v = eval_word(blks, P0s, P1s, sign, what)
                if v is None:
                    continue
                STATS['nonzero'] += 1
                key = canon_word(tuple(cw))
                key_s = key if not split else (key, tuple(d[0] for d in choice))
                acc = res.setdefault(key_s, {w: {0: vzero(), 1: vzero()} for w in what})
                for w in what:
                    for tt in (0, 1):
                        acc[w][tt] = vadd(acc[w][tt], v[w][tt])
    STATS['time'] += time.time() - t0
    return res


def canon_word(w):
    if not w:
        return w
    n = len(w)
    rots = [w[k:] + w[:k] for k in range(n)]
    return min(rots)


# ------------------------------------------------------------------------------------------------ colour reduction
def colour_reduce(res, nl):
    """reduce per-word sums to the independent colour structures:
    n=2: tr(ad_1 ad_0) = C (Y0,Y1)                -> {'(Y0,Y1)': v}
    n=3: tr(ad_c ad_b ad_a) (word (a,b,c)) = (C/2) (Y_c,[Y_b,Y_a]) -> {'(Y0,[Y1,Y2])': v}
    n>=4: kept per cyclic word (tr(ad_{w_m-1}...ad_{w_0}))"""
    out = {}
    for key, val in res.items():
        w = key[0] if isinstance(key[0], tuple) else key
        if nl == 2:
            out.setdefault('(Y0,Y1)', []).append(val)
        elif nl == 3:
            # (Y_c,[Y_b,Y_a]) = -(Y_a,[Y_b,Y_c]) ; for w=(0,1,2): (Y2,[Y1,Y0]) = -(Y0,[Y1,Y2]); orientation sign
            s = perm_sign(list(w))     # w cyclic rotation-canonical, starts with 0: (0,1,2) -> +1, (0,2,1) -> -1
            out.setdefault('(Y0,[Y1,Y2])', []).append((Fraction(-s, 2), val))
        else:
            out.setdefault(w, []).append(val)
    return out
