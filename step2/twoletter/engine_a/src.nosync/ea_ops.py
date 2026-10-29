"""Engine A: planar (single-trace) action of the two-letter kernels on local operators, and the general-colour
reading.  For su(N) with (X,Y) = kappa tr(XY): an adjacent pair X_a X_b of tr(... X_a X_b ...) is replaced by
(-1)^r (c_W / 2) units times tr(... Y_1 ... Y_r ...), c_W the stripped kernel value in units hbar/(16 pi^2) for the
ordered outputs Y_1..Y_r (Y_1 nearest to a) (twoletter_theory Lemma lem:tl-colour (ii),(iii)).
Operators: list of (coefficient, [letters]); a letter = (ltype, zeta, deriv), zeta a Gaussian-integer vector on the
kernel rows of its slot (AV: 4, AW: 6, lam: 16, c: 1, AsV: 4), deriv = {exp: complex} the polynomial f(p) of the
letter's incoming momentum p (d_mu <-> i p_mu).  Legs: (ltype, pol (vector on the same rows), momentum)."""
import itertools
from fractions import Fraction
import numpy as np
import ea_core as E, ea_kern as K

ONE = {(0, 0, 0, 0): 1}


def slotvec_to_stripped(t, v):
    """letter-row vector -> stripped slot vector (used as an output polarisation)"""
    sl = E.SLOT[t]
    out = np.zeros(E.SZ[sl], dtype=complex)
    if t in ('AV', 'AsV'):
        out[:4] = v
    elif t in ('AW', 'AsW'):
        out[4:] = v
    else:
        out[:] = v
    return out


def poly_eval(f, q):
    tot = 0
    for e, c in f.items():
        t = c
        for mu in range(4):
            t = t * (int(q[mu]) ** e[mu])
        tot += t
    return complex(tot)


def letter_on_leg(X, leg):
    t, z, f = X
    lt, pol, q = leg
    if t != lt:
        return 0
    return complex(np.dot(z, pol)) * poly_eval(f, q)


def pair_kernel_on_legs(Xa, Xb, legs, kind='log', topo=None):
    """stripped kernel value (units hbar/16pi^2) of the pair (Xa, Xb) with the ordered outputs `legs` (Y_1 = legs[0]),
    summed over block structures (cubic singles, quartic consecutive A pairs); kind 'log' or m (eps^-m)."""
    ta, za, fa = Xa; tb, zb, fb = Xb
    n = len(legs)
    def rec(i):
        if i == n:
            yield []; return
        for tail in rec(i + 1):
            yield [(i,)] + tail
        if i + 1 < n and E.SLOT[legs[i][0]] == 'A' and E.SLOT[legs[i + 1][0]] == 'A':
            for tail in rec(i + 2):
                yield [(i, i + 1)] + tail
    tot = E.Q0
    for st in rec(0):
        if topo is not None and ''.join('Q' if len(b) == 2 else 'C' for b in st) != topo:
            continue
        blocks = []
        for b in st:
            if len(b) == 1:
                lt, pol, q = legs[b[0]]
                blocks.append(('cub', lt, slotvec_to_stripped(lt, pol), np.asarray(q)))
            else:
                (l1, p1, q1), (l2, p2, q2) = legs[b[0]], legs[b[1]]
                blocks.append(('qua', (slotvec_to_stripped(l1, p1), np.asarray(q1)), (slotvec_to_stripped(l2, p2), np.asarray(q2))))
        kinds = ('log',) if kind == 'log' else ('eps',)
        r = K.kernel(blocks, ta, tb, kinds=kinds, fa=fa, fb=fb, mmax=None if kind == 'log' else kind)
        if not r or kind not in r:
            continue
        tot = E.q_add(tot, K.contract(r[kind], za, zb))
    return tot


def trace_on_legs(letters, legs):
    """colour-ordered evaluation of tr(X_1...X_L) (bosonic letters) on legs in cyclic order"""
    L = len(letters)
    if len(legs) != L:
        return 0
    tot = 0
    for r in range(L):
        v = 1
        for j in range(L):
            v *= letter_on_leg(letters[(j + r) % L], legs[j])
            if v == 0:
                break
        tot += v
    return tot


def lc_planar_on_legs(op, legs, kind='log'):
    """colour-ordered evaluation, on legs in cyclic order, of the planar part of the counterterm of op (bosonic
    letters); for a term with L letters every adjacent pair is replaced by r = len(legs) - (L - 2) >= 1 outputs, with
    the planar factor (-1)^r/2; units u."""
    n = len(legs)
    tot = E.Q0
    for coef, letters in op:
        L = len(letters)
        nout = n - (L - 2)
        if nout < 1 or L < 2:
            continue
        pairs = range(L) if L > 2 else range(2)
        for p in pairs:
            Xa, Xb = letters[p], letters[(p + 1) % L]
            spect = [letters[(p + 2 + j) % L] for j in range(L - 2)]
            for r in range(n):
                rot = [legs[(j + r) % n] for j in range(n)]
                outs = rot[:nout]
                sv = 1
                for j, X in enumerate(spect):
                    sv *= letter_on_leg(X, rot[nout + j])
                    if sv == 0:
                        break
                if sv == 0:
                    continue
                kv = pair_kernel_on_legs(Xa, Xb, outs, kind=kind)
                f = Fraction((-1) ** nout, 2)
                c = complex(coef) * sv
                cr, ci = Fraction(c.real).limit_denominator(10 ** 12), Fraction(c.imag).limit_denominator(10 ** 12)
                kv = (kv[0] * f, kv[1] * f)
                tot = E.q_add(tot, E.q_mul(kv, (cr, ci)))
    return tot


NF = {'c': 1, 'AV': 2, 'AW': 2, 'lamP': 4, 'lamM': 4, 'AsV': -2, 'AsW': -2, 'lamsP': -4, 'lamsM': -4, 'cs': -1}
