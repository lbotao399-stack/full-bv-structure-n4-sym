"""Engine A cross-check route: explicit su(2) colour (r2_core Hessians, colour-full), operator insertion
  value = 1/2 sum_{ordered block sequences} sgn_odd * str[P_0 S H_M P_1 S h_1 ... P_k S h_k] / scales,
both orientations of H_M included (no colour stripping).  Used for gate G4 (colour reading, orientation symmetry,
fermion signs) against the colour-stripped kernels, and for the Ward gates W1 of the antifield letters
(value_vertex_operator: the operator {I_T, N} as the vertex T with one leg fixed)."""
import itertools
from fractions import Fraction
import numpy as np
import ea_core as E
C = E.C
E4 = [np.eye(4, dtype=int)[m] for m in range(4)]
Z4 = np.zeros(4, dtype=complex)


class Ins:
    def __init__(self, bgs):
        self.h0 = C.hessian(bgs, Z4, Z4)
        self.Hj = [C.hessian(bgs, E4[m].astype(complex), Z4) - self.h0 for m in range(4)]
        self.Hl = [C.hessian(bgs, Z4, E4[m].astype(complex)) - self.h0 for m in range(4)]


def Pfull_poly(q):
    const = np.zeros((C.NTOT, C.NTOT), complex); lin = [np.zeros((C.NTOT, C.NTOT), complex) for _ in range(4)]
    for L in ['AV', 'AP', 'gh', 'l']:
        P0 = C.Pnum(np.zeros(4), L)
        const = const + P0 + sum(q[m] * (C.Pnum(E4[m], L) - P0) for m in range(4))
        for m in range(4):
            lin[m] = lin[m] + (C.Pnum(E4[m], L) - P0)
    return E.affine_poly(const, lin)


def ins_poly(ins, qin, qout):
    const = ins.h0 - sum(qin[m] * ins.Hj[m] for m in range(4)) + sum(qout[m] * ins.Hl[m] for m in range(4))
    lin = [ins.Hl[m] - ins.Hj[m] for m in range(4)]
    S = C.SDIAG[:, None]
    return E.affine_poly(S * const, [S * x for x in lin])


def HM_bilinear(ta_slot, za, ca, tb_slot, zb, cb):
    """Hessian of M = (za.X_a)^alpha (zb.X_b)^beta mu_{alpha beta}, mu = outer(ca, cb) (colour vectors);
    H[j,l] = (-1)^{p_j p_l} V(j,l): H[b,a] = za zb mu, H[a,b] = (-1)^{p_a p_b} za zb mu."""
    H = np.zeros((C.NTOT, C.NTOT), complex)
    pa, pb = C.PARITY[ta_slot], C.PARITY[tb_slot]
    def idx(t, comp, col):
        return C.OFF[t] + comp * C.D + col if t not in ('c', 'cs') else C.OFF[t] + col
    for i, x in enumerate(za):
        if x == 0: continue
        for j, y in enumerate(zb):
            if y == 0: continue
            for al in range(3):
                for be in range(3):
                    v = x * y * ca[al] * cb[be]
                    if v == 0: continue
                    ia, ib = idx(ta_slot, i, al), idx(tb_slot, j, be)
                    H[ib, ia] += v
                    H[ia, ib] += (-1) ** (pa * pb) * v
    return H


def perm_sign(seq):
    s_ = 1; a = list(seq)
    for i in range(len(a)):
        for j in range(len(a) - 1 - i):
            if a[j] > a[j + 1]:
                a[j], a[j + 1] = a[j + 1], a[j]; s_ = -s_
    return s_


def _seqs(legs, types):
    if not legs:
        yield []; return
    for i, j in enumerate(legs):
        rest = legs[:i] + legs[i + 1:]
        for tail in _seqs(rest, types):
            yield [(j,)] + tail
        for jj in rest:
            if types[j] == 'A' and types[jj] == 'A' and j < jj:
                rest2 = [x for x in rest if x != jj]
                for tail in _seqs(rest2, types):
                    yield [(j, jj)] + tail


def value(HM, outs, kinds=('log',)):
    """outs: list of (slot type, full pol (SHAPE), mom (int 4-vector)); returns dict of Gaussian rationals"""
    n = len(outs)
    types = [o[0] for o in outs]
    tot = {}
    for seq in _seqs(list(range(n)), types):
        odd = [j for b in seq for j in b if C.PARITY[types[j]]]
        sgn = perm_sign(odd)
        Ps = [sum((np.asarray(outs[j][2]) for j in b), np.zeros(4, dtype=int)) for b in seq]
        P0 = -sum(Ps, np.zeros(4, dtype=int))
        allP = [P0] + Ps
        qs = [np.zeros(4, dtype=int)]
        for i in range(1, len(seq) + 1):
            qs.append(qs[-1] - allP[i - 1])
        scale = 1
        X = {(0, 0, 0, 0): np.eye(C.NTOT, dtype=complex)}
        X = E.pmul(X, Pfull_poly(qs[0]))
        X = E.pmul(X, {(0, 0, 0, 0): C.SDIAG[:, None] * HM})
        for i, b in enumerate(seq, start=1):
            X = E.pmul(X, Pfull_poly(qs[i]))
            ins = Ins([(outs[j][0], outs[j][1], outs[j][2]) for j in b])
            qn = qs[i + 1] if i < len(seq) else qs[0]
            X = E.pmul(X, ins_poly(ins, qs[i], qn))
            scale *= 2 if len(b) == 1 else 8
        N = {e: np.sum(C.SDIAG * np.diag(M)) for e, M in X.items()}
        N = {e: v for e, v in N.items() if v != 0}
        if not N:
            continue
        S = E.Simplex(qs)
        r = S.singular(N, kinds=kinds)
        for k, v in r.items():
            if v is None: continue
            w = (Fraction(sgn, 2 * scale) * v[0], Fraction(sgn, 2 * scale) * v[1])
            tot[k] = E.q_add(tot.get(k, E.Q0), w)
    return tot


def hessian_only(bgs, kj, kl, only):
    H = np.zeros((C.NTOT, C.NTOT), dtype=complex)
    for tj in C.TYPES:
        for tl in C.TYPES:
            legs = [(tj, C.ident(tj), kj), (tl, C.ident(tl), kl)] + [(t, pol[None], m) for (t, pol, m) in bgs]
            v = C.V(legs, only=only)
            if isinstance(v, int):
                continue
            v = v.reshape(C.SIZE[tj], C.SIZE[tl])
            sg = (-1) ** (C.PARITY[tj] * C.PARITY[tl])
            H[C.OFF[tj]:C.OFF[tj] + C.SIZE[tj], C.OFF[tl]:C.OFF[tl] + C.SIZE[tl]] += sg * v
    return H


class InsOnly:
    def __init__(self, bgs, only):
        self.h0 = hessian_only(bgs, Z4, Z4, only)
        self.Hj = [hessian_only(bgs, E4[m].astype(complex), Z4, only) - self.h0 for m in range(4)]
        self.Hl = [hessian_only(bgs, Z4, E4[m].astype(complex), only) - self.h0 for m in range(4)]


def value_vertex_operator(vt, zleg, outs, kinds=('log',)):
    """lc of the operator M = (vertex vt with one leg fixed to zleg = (type, full pol, momentum)), i.e. {I_vt, N} up
    to a common constant, on the output legs outs.  Cycles through M once (M with 0 or 1 spectator legs)."""
    n = len(outs)
    types = [o[0] for o in outs]
    arity = len(C.CANON[vt])
    nspec = arity - 3
    tot = {}
    for Sm in itertools.combinations(range(n), nspec):
        rest = [j for j in range(n) if j not in Sm]
        for seq in _seqs(rest, types):
            seq = [tuple(Sm)] + seq
            odd = [j for b in seq for j in b if C.PARITY[types[j]]]
            sgn = perm_sign(odd)
            Ps = [sum((np.asarray(outs[j][2]) for j in b), np.zeros(4, dtype=int)) for b in seq]
            Ps[0] = Ps[0] + np.asarray(zleg[2])
            qs = [np.zeros(4, dtype=int)]
            for i in range(1, len(seq)):
                qs.append(qs[-1] - Ps[i - 1])
            assert np.all(sum(Ps, np.zeros(4, dtype=int)) == 0)
            scale = 2 if arity == 3 else 8
            X = {(0, 0, 0, 0): np.eye(C.NTOT, dtype=complex)}
            for i, b in enumerate(seq):
                X = E.pmul(X, Pfull_poly(qs[i]))
                qn = qs[i + 1] if i + 1 < len(seq) else qs[0]
                if i == 0:
                    ins = InsOnly([zleg] + [(outs[j][0], outs[j][1], outs[j][2]) for j in b], [vt])
                else:
                    ins = Ins([(outs[j][0], outs[j][1], outs[j][2]) for j in b])
                    scale *= 2 if len(b) == 1 else 8
                X = E.pmul(X, ins_poly(ins, qs[i], qn))
            N = {e: np.sum(C.SDIAG * np.diag(M)) for e, M in X.items()}
            N = {e: v for e, v in N.items() if v != 0}
            if not N:
                continue
            S = E.Simplex(qs)
            r = S.singular(N, kinds=kinds)
            for k, v in r.items():
                if v is None: continue
                w = (Fraction(sgn, 2 * scale) * v[0], Fraction(sgn, 2 * scale) * v[1])
                tot[k] = E.q_add(tot.get(k, E.Q0), w)
    return tot


def tree_vertex_operator(vt, zleg, outs):
    """tree value V_vt(zleg, outs) / scale (odd parameters of the legs in leg order, zleg first)"""
    legs = [(zleg[0], zleg[1][None], zleg[2])] + [(o[0], o[1][None], o[2]) for o in outs]
    v = C.V(legs, only=[vt])
    if isinstance(v, int):
        return 0
    return complex(v.reshape(-1)[0]) / (2 if len(C.CANON[vt]) == 3 else 8)
