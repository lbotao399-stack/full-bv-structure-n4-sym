"""Engine A, part C: composite kernels evaluated on legs.

A value is an exact element of (Q + iQ) (x) span{1, log2, log3, log5, log7}: a list of five Gaussian rationals (re, im).
pair_on_legs(Xa, Xb, legs, kind) returns the stripped value (units hbar/(16 pi^2), colour word ad_{Y_r}...ad_{Y_1} with
Y_1 = legs[0] nearest to X_a) of C^(2)_x(X_a, X_b) with the ordered outputs `legs`, summed over the block structures of
I (cubic singles, quartic consecutive A pairs), the positions of the internal vertex and of the special edge.
kind: ('P', b) the translation kernel C^(2)_{P_b}; 'D0' the translation-invariant part C^(2),0_Delta (internal vertex
s_0 at r = 0 plus -i d/dr_mu of the vertex s_mu with injected momentum r); ('Ds', mu) the vertex s_mu at r = 0 (= 2 x
the P_mu graph, position part); poles: kind ('Pp', b, m) etc. are the eps^-m coefficients.
Letters: (ltype, zeta, deriv) with deriv = {exp: complex} the polynomial of the letter's incoming momentum (d <-> i p).
Legs: (ltype, pol (vector on the letter rows), momentum (integer 4-vector))."""
import itertools
from fractions import Fraction
import numpy as np
import cq_core as Q
import ea_core as E
import ea_kern as EK
import ea_ops as EO

NB = 5
ZERO = [(Fraction(0), Fraction(0)) for _ in range(NB)]


def vzero():
    return [(Fraction(0), Fraction(0)) for _ in range(NB)]


def vadd(a, b):
    return [(x[0] + y[0], x[1] + y[1]) for x, y in zip(a, b)]


def vscale(a, z):
    """multiply by a Gaussian rational z = (re, im)"""
    zr, zi = z
    return [(x[0] * zr - x[1] * zi, x[0] * zi + x[1] * zr) for x in a]


def viszero(a):
    return all(x[0] == 0 and x[1] == 0 for x in a)


def vstr(a):
    parts = []
    for nm, (r, i) in zip(Q.BASIS_NAMES, a):
        if r == 0 and i == 0:
            continue
        c = E.qstr((r, i))
        parts.append(c if nm == '1' else '%s*%s' % (c, nm))
    return ' + '.join(parts) if parts else '0'


def from_fin(fin, za, zb):
    """configuration finite part (rational, log(E/(E-1)) part, E) contracted with the letter vectors -> value"""
    f0, f1, E_ = fin
    out = vzero()
    if f0 is not None:
        out[0] = EK.contract(f0, za, zb)
    if f1 is not None:
        c = EK.contract(f1, za, zb)
        for bi, n in enumerate(Q.log_basis(E_)):
            if n:
                out[1 + bi] = (out[1 + bi][0] + n * c[0], out[1 + bi][1] + n * c[1])
    return out


def structures(legs):
    n = len(legs)
    def rec(i):
        if i == n:
            yield []; return
        for tail in rec(i + 1):
            yield [(i,)] + tail
        if i + 1 < n and E.SLOT[legs[i][0]] == 'A' and E.SLOT[legs[i + 1][0]] == 'A':
            for tail in rec(i + 2):
                yield [(i, i + 1)] + tail
    return list(rec(0))


def blocks_of(st, legs):
    out = []
    for b in st:
        if len(b) == 1:
            lt, pol, q = legs[b[0]]
            out.append(('cub', lt, EO.slotvec_to_stripped(lt, pol), np.asarray(q)))
        else:
            (l1, p1, q1), (l2, p2, q2) = legs[b[0]], legs[b[1]]
            out.append(('qua', (EO.slotvec_to_stripped(l1, p1), np.asarray(q1)), (EO.slotvec_to_stripped(l2, p2), np.asarray(q2))))
    return out


STATS = {'configs': 0, 'nonzero': 0}


VV = {'kappa': 1}          # vertex-vertex special edges: kappa = +1 (VV_a), -1 (VV_b), 0 (VV_0)


def _raw(Xa, Xb, legs, v2kind, v2dir=None, r=None, what='fin', m=None, plain=False, split=False):
    """sum over structures and positions of one internal-vertex kind; returns a value (what='fin'),
    an eps^-m value (what='pole'), or for plain cycles {'log': value, m: value}"""
    ta, za, fa = Xa; tb, zb, fb = Xb
    tot = vzero(); tvv = vzero()
    for st in structures(legs):
        outs = blocks_of(st, legs)
        nb = len(outs)
        for v2pos in range(nb + 1):
            if plain:
                res = Q.plain_config_value(outs, v2pos, v2kind, ta, tb, fa=fa, fb=fb, v2dir=v2dir, r=r)
                STATS['configs'] += 1
                if not res:
                    continue
                STATS['nonzero'] += 1
                key = 'log' if what == 'log' else m
                if key in res and res[key] is not None:
                    c = EK.contract(res[key], za, zb)
                    tot[0] = (tot[0][0] + c[0], tot[0][1] + c[1])
                continue
            for spos in range(nb + 2):
                STATS['configs'] += 1
                want = ('fin',) if what == 'fin' else ('poles',)
                res = Q.config_value(outs, v2pos, spos, v2kind, ta, tb, fa=fa, fb=fb, v2dir=v2dir, r=r, want=want)
                if not res:
                    continue
                STATS['nonzero'] += 1
                tgt = tvv if res['vv'] else tot
                if what == 'fin':
                    v = from_fin(res['fin'], za, zb)
                else:
                    v = vzero()
                    if m in res and res[m] is not None:
                        v[0] = EK.contract(res[m], za, zb)
                for i_ in range(NB):
                    tgt[i_] = (tgt[i_][0] + v[i_][0], tgt[i_][1] + v[i_][1])
    if split:
        return tot, tvv
    return vadd(tot, vscale(tvv, (Fraction(VV['kappa']), Fraction(0))))


def _deriv0(f, dmax):
    """exact d/dt at t = 0 of a polynomial of degree <= dmax given by its values f(t) (callable), symmetric points"""
    mm = max(1, (dmax + 1) // 2)
    ts = [t for t in range(-mm, mm + 1) if t != 0]
    # odd part only: f(t) - f(-t) = 2 sum_{n odd} c_n t^n ; solve for c_1
    ts = list(range(1, mm + 1))
    vals = [vadd(f(t), vscale(f(-t), (Fraction(-1), Fraction(0)))) for t in ts]
    # Vandermonde in odd powers: sum_j c_{2j+1} t^{2j+1} = vals/2
    n = len(ts)
    A = [[Fraction(t) ** (2 * j + 1) for j in range(n)] for t in ts]
    # invert A exactly (small), take the first row of A^{-1}
    inv = _inv(A)
    out = vzero()
    for i in range(n):
        out = vadd(out, vscale(vals[i], (inv[0][i] / 2, Fraction(0))))
    return out


def _inv(A):
    n = len(A)
    M = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(A)]
    for i in range(n):
        p = next(r for r in range(i, n) if M[r][i] != 0)
        M[i], M[p] = M[p], M[i]
        piv = M[i][i]
        M[i] = [x / piv for x in M[i]]
        for r in range(n):
            if r != i and M[r][i] != 0:
                f = M[r][i]
                M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [row[n:] for row in M]


def pair_on_legs(Xa, Xb, legs, kind, dmax=6):
    if isinstance(kind, tuple) and kind[0] == 'P':
        return _raw(Xa, Xb, legs, 'P', v2dir=kind[1])
    if isinstance(kind, tuple) and kind[0] == 'Ds':
        return _raw(Xa, Xb, legs, 's', v2dir=kind[1], r=np.zeros(4, dtype=int))
    if kind == 'D0':
        tot = _raw(Xa, Xb, legs, 's0')
        for mu in range(4):
            def f(t, mu=mu):
                r = np.zeros(4, dtype=int); r[mu] = t
                return _raw(Xa, Xb, legs, 's', v2dir=mu, r=r)
            d = _deriv0(f, dmax)
            tot = vadd(tot, vscale(d, (Fraction(0), Fraction(-1))))      # -i d/dr_mu
        return tot
    if isinstance(kind, tuple) and kind[0] == 'Pp':          # poles of Z along P_b: ('Pp', b, m)
        return _raw(Xa, Xb, legs, 'P', v2dir=kind[1], what='pole', m=kind[2])
    if isinstance(kind, tuple) and kind[0] == 'Plog':        # plain composite cycles along P_b, log eps
        return _raw(Xa, Xb, legs, 'P', v2dir=kind[1], what='log', plain=True)
    if isinstance(kind, tuple) and kind[0] == 'Ppole':       # plain composite cycles along P_b, eps^-m
        return _raw(Xa, Xb, legs, 'P', v2dir=kind[1], what='pole', m=kind[2], plain=True)
    raise ValueError(kind)


def letter_on_leg(X, leg):
    return EO.letter_on_leg(X, leg)
