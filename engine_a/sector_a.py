"""Item (a): quadratic one-loop counterterms (engine A).  All W_1 coefficients are in units 1/(16 pi^2) (hbar stripped)."""
import json, sys
import numpy as np
from fractions import Fraction
from oneloop import Engine, PARITY
import check_vertices as cv
from grass import GA

def poly_parts(F, k0, dirs):
    """F quadratic in k: return F(k0), Laplacian, and a cubic-term check"""
    f0 = F(k0)
    lap = 0; chk = 0
    for e in dirs:
        fp, fm = F(k0 + e), F(k0 - e)
        lap += fp + fm - 2 * f0
        chk = max(chk, abs(F(k0 + 2 * e) - 3 * fp + 3 * f0 - fm))
    return f0, lap, chk

E4 = [np.eye(4)[m] for m in range(4)]

def bubble_sing(eng, b1, b2, p):
    """sum of n=2 bubble terms 1/4[str(X1X2)+str(X2X1)] for backgrounds b1 at p, b2 at -p: returns (eps_inv, log, details)"""
    S = np.diag(eng.S)
    t1, pol1, _ = b1; t2, pol2, _ = b2
    sgn21 = (-1) ** (PARITY[t1] * PARITY[t2])
    def Num12(k):
        h1 = eng.hessian([(t1, pol1, p)], -k, k - p)
        h2 = eng.hessian([(t2, pol2, -p)], p - k, k)
        return eng.str_(eng.Nprop(k) @ S @ h1 @ eng.Nprop(k - p) @ S @ h2)
    def Num21(k):
        h2 = eng.hessian([(t2, pol2, -p)], -k, k + p)
        h1 = eng.hessian([(t1, pol1, p)], -k - p, k)
        return sgn21 * eng.str_(eng.Nprop(k) @ S @ h2 @ eng.Nprop(k + p) @ S @ h1)
    tot_inv = 0; tot_log = 0; det = {}
    for name, F, sg in (('12', Num12, +1), ('21', Num21, -1)):
        # k0 = sg * u2 * p ; g(u2) quadratic -> Simpson-exact integral over [0,1]
        g0 = F(0 * p); g1 = F(sg * 0.5 * p); g2 = F(sg * 1.0 * p)
        intg = (g0 + 4 * g1 + g2) / 6
        _, lap, chk = poly_parts(F, 0 * p, E4)
        c = lap / 4
        s2 = float(p @ p)
        e_inv = 0.25 * c / 4
        e_log = 0.25 * (-intg + s2 * c / 6)
        tot_inv += e_inv; tot_log += e_log
        det[name] = dict(cubic_check=float(chk))
    return tot_inv, tot_log, det

def seagull_sing(eng, b1, b2, p):
    S = np.diag(eng.S)
    k = np.array([0.3, -0.7, 1.1, 0.2])
    h4 = eng.hessian([b1[:2] + (p,), b2[:2] + (-p,)], -k, k)
    val = eng.str_(eng.Nprop(k) @ S @ h4)
    return 0.5 * val, 0.0

# ---------- target functionals (cross coefficient on the two-wave background), brute force from densities
def Bdensity(name, fld, der):
    d = cv.d; tot = GA()
    A = fld['A']; dA = der['A']
    if name == 'B1':
        for m in range(4):
            for n_ in range(m + 1, 4):
                for a in range(d):
                    F = dA[m, n_, a] - dA[n_, m, a]
                    tot = tot + F * F
    elif name == 'B2':
        for a in range(d):
            div = GA()
            for m in range(4):
                div = div - dA[m, m, a]
            tot = tot + div * div
    elif name == 'B3':
        for m in range(4):
            for a in range(d):
                tot = tot + A[m, a] * A[m, a]
    elif name == 'B4':
        for m in range(4):
            for i in range(4, 10):
                for a in range(d):
                    tot = tot + dA[m, i, a] * dA[m, i, a]
    elif name == 'B5':
        for i in range(4, 10):
            for a in range(d):
                tot = tot + A[i, a] * A[i, a]
    elif name == 'B6':
        from cliff import sb
        for m in range(4):
            for al in range(16):
                for be in range(16):
                    if sb[m][al, be] != 0:
                        for a in range(d):
                            tot = tot + fld['l'][al, a] * der['l'][m, be, a] * complex(sb[m][al, be])
    elif name == 'B7':
        for m in range(4):
            for a in range(d):
                tot = tot + der['c'][m, a] * fld['As'][m, a]
    return tot

def Bcross(name, b1, b2, p):
    fld, der, vs = cv.build([(b1[0], b1[1], p), (b2[0], b2[1], -p)])
    return complex(Bdensity(name, fld, der).coeff(vs))

def rat(x, tol=1e-7):
    fr = Fraction(x).limit_denominator(10000)
    return fr if abs(float(fr) - x) < tol else None

def fit(samples, names):
    """samples: list of (w_inv, w_log, {B: cross}) -> least squares coefficients"""
    M = np.array([[smp[2][nm] for nm in names] for smp in samples])
    out = {}
    for idx, key in ((0, 'eps_inv'), (1, 'log_eps')):
        rhs = np.array([smp[idx] for smp in samples])
        x, res, rk, _ = np.linalg.lstsq(M, rhs, rcond=None)
        resid = np.linalg.norm(M @ x - rhs)
        out[key] = (x, resid, rk)
    return out
