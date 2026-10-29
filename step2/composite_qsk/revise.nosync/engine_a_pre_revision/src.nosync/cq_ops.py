"""Engine A, part C: planar action of the composite kernels on single-trace operators (bosonic letters), as
ea_ops.lc_planar_on_legs: every adjacent pair (both orders for two letters) is replaced by r = len(legs) - (L - 2)
outputs with the planar factor (-1)^r/2 (Lemma lem:tl-colour (iv)); units u.  Returns (letter part, vertex-vertex part)."""
from fractions import Fraction
import numpy as np
import cq_gval as GV, cq_eval as V
import ea_ops as EO

vzero, vadd, vscale = V.vzero, V.vadd, V.vscale


def poly_eval(f, q):
    tot = 0
    for e, c in f.items():
        t = complex(c)
        for mu in range(4):
            t = t * (complex(q[mu]) ** e[mu])
        tot += t
    return complex(tot)


def letter_on_leg(X, leg):
    t, z, f = X
    lt, pol, q = leg
    if t != lt:
        return 0
    return complex(np.dot(z, pol)) * poly_eval(f, q)


def gauss(z):
    z = complex(z)
    r, i = round(z.real), round(z.imag)
    assert abs(z.real - r) < 1e-9 and abs(z.imag - i) < 1e-9, z
    return (Fraction(r), Fraction(i))


def trace_on_legs(letters, legs):
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


def op_on_legs(op, legs):
    return sum(complex(c) * trace_on_legs(lt, legs) for c, lt in op)


def planar_on_legs(op, legs, kind, cache=None):
    n = len(legs)
    tl, tv = vzero(), vzero()
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
                key = None
                if cache is not None:
                    key = (_lkey(Xa), _lkey(Xb), tuple(_legkey(x) for x in outs), str(kind))
                if key is not None and key in cache:
                    a, b = cache[key]
                else:
                    a, b = GV.kernel(Xa, Xb, outs, kind)
                    if key is not None:
                        cache[key] = (a, b)
                f = Fraction((-1) ** nout, 2)
                c = gauss(complex(coef) * sv)
                tl = vadd(tl, vscale(vscale(a, (f, Fraction(0))), c))
                tv = vadd(tv, vscale(vscale(b, (f, Fraction(0))), c))
    return tl, tv


def _lkey(X):
    t, z, f = X
    return (t, tuple(complex(x) for x in z), tuple(sorted((e, complex(c)) for e, c in f.items())))


def _legkey(l):
    t, pol, q = l
    return (t, tuple(complex(x) for x in pol), tuple(complex(x) for x in q))
