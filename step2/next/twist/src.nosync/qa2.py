"""Arity-two extension of engine A's defect cycles (qa_core): chains with exactly nS = 2 vertices of S^(1)(xi)
(quadratic blocks S2 at gaps, or cubic T5 blocks S3 with one output), one special line, all placements.
For odd xi |S^(1)(xi)| = 0, so the configuration sign is that of qa_core.config_sign (no S-dependent Koszul factor);
the twisted interaction I + S^(1)(xi) - 1/2 S^(2)(xi,xi) (Lemma c8-ev, q = 1) carries S^(1)(xi) with coefficient 1, so
the sum over chains with two S^(1)(xi) vertices is the (xi,xi) part of Q^xi_1 directly (no extra 1/2: the two
vertices sit at distinct chain positions).  S^(2) cannot sit on a cycle (both legs lambda*), cf. lemma S2."""
import env
import itertools
from fractions import Fraction
import numpy as np
import qa_core as QC
import ea_core as E, ea_ops as EO, cq_core as Q, cq_gcore as G

def structures2(legs, gen, nS=2):
    n = len(legs)
    def rec(i, used):
        if i == n:
            yield [], used; return
        lt = legs[i][0]
        for tail, u in rec(i + 1, used):
            yield [('I1', i)] + tail, u
        if gen.cub is not None and used < nS and lt in gen.cubtypes:
            for tail, u in rec(i + 1, used + 1):
                yield [('S3', i)] + tail, u
        if i + 1 < n and E.SLOT[lt] == 'A' and E.SLOT[legs[i + 1][0]] == 'A':
            for tail, u in rec(i + 2, used):
                yield [('I2', i, i + 1)] + tail, u
    for st, used in rec(0, 0):
        k = nS - used
        if k == 0:
            yield st; continue
        if gen.quad is None:
            continue
        # insert k quadratic blocks at gaps 0..len(st) (multiset of gap positions)
        for gaps in itertools.combinations_with_replacement(range(len(st) + 1), k):
            out = []
            gi = 0
            for pos in range(len(st) + 1):
                while gi < len(gaps) and gaps[gi] == pos:
                    out.append(('S2',)); gi += 1
                if pos < len(st):
                    out.append(st[pos])
            yield out

def configs2(legs, gen, nS=2):
    for st in structures2(legs, gen, nS):
        blks, P0s, P1s, scale, odd, spos = QC.build(st, legs, gen)
        n = len(blks)
        for sig in range(n + 1):
            lines = ['P'] * (n + 1); lines[sig] = 'S'
            yield st, blks, P0s, P1s, scale, odd, spos, sig, lines

def pair_value2(Xa, Xb, legs, gen, what=('fin',), nS=2):
    ta, za, fa = Xa; tb, zb, fb = Xb
    tot = {w: QC.vzero() for w in what}
    intr = -1 if (QC.RULE['intr'] and (Q.PARITY_T[ta] + Q.PARITY_T[tb]) % 2 == 1) else 1
    for st, blks, P0s, P1s, scale, odd, spos, sig, lines in configs2(legs, gen, nS):
        K, qs, P0 = G.chain(blks, P0s, P1s, lines, ta, tb)
        if not K: continue
        if fa is not None or fb is not None:
            K = G.apply_letters(K, fa or {(0, 0, 0, 0): 1}, fb or {(0, 0, 0, 0): 1}, P0)
            if not K: continue
        n = len(blks)
        sgn = QC.config_sign(odd, spos, sig, n, gen.par) * intr
        Sx = G.GSimplex(qs)
        E_ = len(qs[0])
        for w in what:
            r = Sx.evaluate(K, 'fin', s=sig) if w == 'fin' else Sx.evaluate(K, 'spec', s=sig, m=w[1])
            rl = r[0]
            rl = [None if v is None else (v[0] * Fraction(sgn, scale), v[1] * Fraction(sgn, scale)) for v in rl]
            tot[w] = QC.vadd(tot[w], QC._from_rl(rl, E_, za, zb))
    return tot
