"""Brute-force Grassmann evaluation of the [ESW (31)] density on plane waves; compares with vertices.Vtx."""
import numpy as np
from grass import GA, newvar
from cliff import sb
from vertices import Vtx, structure_constants
rng = np.random.default_rng(7)
f = structure_constants('su2'); d = 3
nz = [(a, b, c, f[a, b, c]) for a in range(d) for b in range(d) for c in range(d) if f[a, b, c] != 0]
SH = {'c': (d,), 'A': (10, d), 'l': (16, d), 'As': (10, d), 'ls': (16, d), 'cs': (d,)}
ODD = {'c': 1, 'A': 0, 'l': 1, 'As': 1, 'ls': 0, 'cs': 0}
_cnt = [5000]

def zero(shape):
    z = np.empty(shape, dtype=object)
    for idx in np.ndindex(shape):
        z[idx] = GA()
    return z

def build(legs):
    """legs: list of (name, pol, mom). returns field arrays and derivative arrays, and list of leg vars"""
    fld = {k: zero(v) for k, v in SH.items()}
    der = {k: zero((4,) + v) for k, v in SH.items()}
    vs = []
    for (nm, pol, mom) in legs:
        v = newvar(_cnt[0], ODD[nm]); _cnt[0] += 1; vs.append(v)
        g = GA.var(v)
        for idx in np.ndindex(SH[nm]):
            if pol[idx] != 0:
                fld[nm][idx] = fld[nm][idx] + g * complex(pol[idx])
                for mu in range(4):
                    if mom[mu] != 0:
                        der[nm][(mu,) + idx] = der[nm][(mu,) + idx] + g * complex(1j * mom[mu] * pol[idx])
    return fld, der, vs

def br(x, y):          # colour bracket of GA colour vectors
    out = [GA() for _ in range(d)]
    for a, b, c, v in nz:
        if x[a].t and y[b].t:
            out[c] = out[c] + x[a] * y[b] * v
    return out

def density(fld, der):
    A, l, c = fld['A'], fld['l'], fld['c']
    tot = GA()
    F = {}
    for M in range(10):
        for N in range(M + 1, 10):
            bb = br(A[M], A[N])
            F[(M, N)] = [(der['A'][M, N, a] if M < 4 else GA()) - (der['A'][N, M, a] if N < 4 else GA()) + bb[a] for a in range(d)]
            for a in range(d):
                tot = tot + F[(M, N)][a] * F[(M, N)][a] * (-0.25)
    for M in range(10):
        bl = [br(A[M], l[al]) for al in range(16)]
        for al in range(16):
            for be in range(16):
                if sb[M][al, be] == 0:
                    continue
                for a in range(d):
                    Dl = (der['l'][M, be, a] if M < 4 else GA()) + bl[be][a]
                    tot = tot + l[al, a] * Dl * (0.5 * complex(sb[M][al, be]))
    for M in range(10):
        bc = br(A[M], c)
        for a in range(d):
            Dc = (der['c'][M, a] if M < 4 else GA()) + bc[a]
            tot = tot + Dc * fld['As'][M, a]
    for al in range(16):
        bcl = br(c, l[al])
        for a in range(d):
            tot = tot + bcl[a] * fld['ls'][al, a]
    bcc = br(c, c)
    for a in range(d):
        tot = tot + bcc[a] * fld['cs'][a] * 0.5
    return tot

def brute(legs):
    fld, der, vs = build(legs)
    return density(fld, der).coeff(vs)

def ri(*shape):
    return rng.integers(-2, 3, size=shape).astype(float)

if __name__ == '__main__':
    V = Vtx(f)
    res = {}
    p1, p2 = ri(4), ri(4); p3 = -p1 - p2
    a1, a2, a3, a4 = ri(10, d), ri(10, d), ri(10, d), ri(10, d)
    res['AAA'] = (brute([('A', a1, p1), ('A', a2, p2), ('A', a3, p3)]), V.AAA(a1, p1, a2, p2, a3, p3))
    z = np.zeros(4)
    res['AAAA'] = (brute([('A', a1, z), ('A', a2, z), ('A', a3, z), ('A', a4, z)]), V.AAAA(a1, a2, a3, a4))
    u1, u2 = ri(16, d), ri(16, d)
    res['llA'] = (brute([('l', u1, p1), ('l', u2, p2), ('A', a1, p3)]), V.llA(u1, u2, a1))
    g1, g2, al, v, w = ri(d), ri(d), ri(10, d), ri(16, d), ri(d)
    res['AcAs'] = (brute([('A', a1, p1), ('c', g1, p2), ('As', al, p3)]), V.AcAs(a1, g1, al))
    res['cllS'] = (brute([('c', g1, p1), ('l', u1, p2), ('ls', v, p3)]), V.cllS(g1, u1, v))
    res['cccS'] = (brute([('c', g1, p1), ('c', g2, p2), ('cs', w, p3)]), V.cccS(g1, g2, w))
    for k, (x, yv) in res.items():
        print(k, complex(x), complex(yv), 'OK' if abs(complex(x) - complex(yv)) < 1e-9 else 'MISMATCH')
