"""Engine A, part Q: small exact helpers (random Gaussian-integer data, letters, legs, value arithmetic)."""
from fractions import Fraction
import numpy as np
import qa_s1 as S1
import qa_core as QC
import cq_eval as CV

vzero, vadd, vscale, vstr, viszero = CV.vzero, CV.vadd, CV.vscale, CV.vstr, CV.viszero
G5S, G5SS = S1.G5S, S1.G5SS
I16 = np.eye(16, dtype=int)


class Rng:
    def __init__(self, seed):
        self.r = np.random.default_rng(seed)

    def gv(self, n, lo=-2, hi=3):
        return np.array([complex(int(a), int(b)) for a, b in zip(self.r.integers(lo, hi, n), self.r.integers(lo, hi, n))])

    def iv(self, n, lo=-2, hi=3):
        return self.r.integers(lo, hi, n)

    def vec(self, t):
        if t in ('AW', 'AsW'):
            return self.gv(6)
        if t in ('AV', 'AsV'):
            return self.gv(4)
        if t in ('c', 'cs'):
            return np.array([complex(int(self.r.integers(1, 3)), 0)])
        if t == 'lamP':
            return (I16 + G5S) @ self.gv(16)
        if t == 'lamM':
            return (I16 - G5S) @ self.gv(16)
        if t == 'lamsP':
            return (I16 + G5SS) @ self.gv(16)
        if t == 'lamsM':
            return (I16 - G5SS) @ self.gv(16)
        raise ValueError(t)

    def mom(self):
        while True:
            q = self.iv(4)
            if np.any(q):
                return q


def gdiv(v, z):
    """value / Gaussian rational z"""
    zr, zi = Fraction(z[0]), Fraction(z[1])
    n = zr * zr + zi * zi
    return vscale(v, (zr / n, -zi / n))


def gq(z):
    return QC.gq(z)


def vsub(a, b):
    return vadd(a, vscale(b, (Fraction(-1), Fraction(0))))


def vneg(a):
    return vscale(a, (Fraction(-1), Fraction(0)))


def vjson(v):
    return vstr(v)


def both_orders(Xa, Xb, legs, gen, **kw):
    """sum of the kernel value over the two orderings of two outputs (colour-singlet bilinear reading)"""
    a = QC.pair_value(Xa, Xb, legs, gen, **kw)
    b = QC.pair_value(Xa, Xb, legs[::-1], gen, **kw)
    return {k: vadd(a[k], b[k]) for k in a}
