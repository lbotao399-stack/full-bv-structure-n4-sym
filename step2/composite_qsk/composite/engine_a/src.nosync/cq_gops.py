"""Engine A, part C: graded planar evaluation of single-trace operators with odd letters on odd legs.

term c tr(X_1 ... X_L); a pair (X_p, X_{p+1}) (cyclic; both orders for L = 2) is brought to the front by the graded
cyclic move tr(X_1..X_L) = (-1)^{(sum_{i<p}|X_i|)(sum_{i>=p}|X_i|)} tr(X_p .. X_{p-1}) and replaced by the kernel outputs
Y_1..Y_r with the planar factor (-1)^r/2 (Lemma lem:tl-colour (iii), (iv)); the colour-ordered evaluation on legs l_1..l_n
sums over the cyclic rotations rho of the legs with the sign of the permutation of the odd legs (theta order of the
factors of the functional versus the reference order of the legs).  The kernel value of an output word with legs is the
coefficient of the legs' odd parameters in chain order (engine convention, as ea_su11)."""
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_eval as V, cq_ops as CO

vzero, vadd, vscale = V.vzero, V.vadd, V.vscale


def par_letter(X):
    return Q.PARITY_T[X[0]]


def perm_sign(order, ref):
    """sign of the permutation of the odd legs: order, ref are lists of leg indices (odd ones only)"""
    pos = {x: i for i, x in enumerate(ref)}
    seq = [pos[x] for x in order]
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def odd_legs(legs, idx):
    return [i for i in idx if Q.PARITY_T[legs[i][0]] == 1]


def trace_graded(letters, legs):
    L = len(letters)
    if len(legs) != L:
        return 0
    ref = odd_legs(legs, list(range(L)))
    tot = 0
    for r in range(L):
        idx = [(j + r) % L for j in range(L)]
        v = 1
        for j in range(L):
            v *= CO.letter_on_leg(letters[j], legs[idx[j]])
            if v == 0:
                break
        if v == 0:
            continue
        tot += perm_sign(odd_legs(legs, idx), ref) * v
    return tot


def op_graded(op, legs):
    return sum(complex(c) * trace_graded(lt, legs) for c, lt in op)


def planar_graded(op, legs, kfun, cache=None):
    """kfun(Xa, Xb, outs) -> (letter part, vertex-vertex part) or a single value (then returned as such)"""
    n = len(legs)
    ref = odd_legs(legs, list(range(n)))
    tl, tv = vzero(), vzero()
    single = False
    for coef, letters in op:
        L = len(letters)
        nout = n - (L - 2)
        if nout < 1 or L < 2:
            continue
        pars = [par_letter(X) for X in letters]
        pairs = range(L) if L > 2 else range(2)
        for p in pairs:
            eps = (-1) ** ((sum(pars[:p]) * sum(pars[p:])) % 2)
            W = letters[p:] + letters[:p]
            Xa, Xb = W[0], W[1]
            spect = W[2:]
            for r in range(n):
                idx = [(j + r) % n for j in range(n)]
                rot = [legs[i] for i in idx]
                outs = rot[:nout]
                sv = 1
                for j, X in enumerate(spect):
                    sv *= CO.letter_on_leg(X, rot[nout + j])
                    if sv == 0:
                        break
                if sv == 0:
                    continue
                key = None
                if cache is not None:
                    key = (CO._lkey(Xa), CO._lkey(Xb), tuple(CO._legkey(x) for x in outs))
                if key is not None and key in cache:
                    res = cache[key]
                else:
                    res = kfun(Xa, Xb, outs)
                    if key is not None:
                        cache[key] = res
                sg = eps * perm_sign(odd_legs(legs, idx), ref)
                f = Fraction(sg * (-1) ** nout, 2)
                c = CO.gauss(complex(coef) * sv)
                if isinstance(res, tuple):
                    a, b = res
                    tl = vadd(tl, vscale(vscale(a, (f, Fraction(0))), c))
                    tv = vadd(tv, vscale(vscale(b, (f, Fraction(0))), c))
                else:
                    single = True
                    tl = vadd(tl, vscale(vscale(res, (f, Fraction(0))), c))
    return tl if single else (tl, tv)


def covariant_op(letter_k, k, N, An):
    """tr(D^k X D^{N-k} X) to first order in A_n, D = d + [A_n, -] (X any parity; A_n even):
    letter_k(m) = d^m X, An(m) = d^m A_n"""
    op = [(1, [letter_k(k), letter_k(N - k)])]
    for i in range(k):
        c = __import__('math').comb(k, i + 1)
        op.append((c, [An(i), letter_k(k - 1 - i), letter_k(N - k)]))
        op.append((-c, [letter_k(k - 1 - i), An(i), letter_k(N - k)]))
    for i in range(N - k):
        c = __import__('math').comb(N - k, i + 1)
        op.append((c, [letter_k(k), An(i), letter_k(N - k - 1 - i)]))
        op.append((-c, [letter_k(k), letter_k(N - k - 1 - i), An(i)]))
    return op
