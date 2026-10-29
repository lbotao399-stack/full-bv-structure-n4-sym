"""Plug-back (refuter #2, own script): the maps Q H_j, {I_3, B_k}, the F_3 family, the kernel, the
F_3 = 0 representative, the degree-2 sector and the Wess-Zumino identity, on random plane waves in my own
realisation (rb_alg / rb_func).  Exact: dyadic floats, rational reconstruction with zero residual.
"""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
from rb_func import *

rng = np.random.default_rng(777)


def rint(shape, lo=-3, hi=4):
    return rng.integers(lo, hi, size=shape).astype(float)


def rand_mom(n):
    ps = [rint(4) for _ in range(n - 1)]
    ps.append(-sum(ps))
    return ps


def leg(slot, p, part=None):
    v = rint(SHAPE[slot])
    if slot in ('A', 'Ast') and part == 'V':
        v[4:] = 0
    if slot in ('A', 'Ast') and part == 'W':
        v[:4] = 0
    return (slot, v, p)


def sector_legs(name):
    p = rand_mom(3)
    if name == 'cAVAV':
        return [leg('c', p[0]), leg('A', p[1], 'V'), leg('A', p[2], 'V')]
    if name == 'cAWAW':
        return [leg('c', p[0]), leg('A', p[1], 'W'), leg('A', p[2], 'W')]
    if name == 'cll':
        return [leg('c', p[0]), leg('lam', p[1]), leg('lam', p[2])]
    if name == 'ccAs':
        return [leg('c', p[0]), leg('c', p[1]), leg('Ast', p[2], 'V')]
    if name == 'ccAV':       # degree-2 sector
        return [leg('c', p[0]), leg('c', p[1]), leg('A', p[2], 'V')]
    # no-invariant degree-1 sectors
    if name == 'cAVAW':
        return [leg('c', p[0]), leg('A', p[1], 'V'), leg('A', p[2], 'W')]
    if name == 'cAl':
        return [leg('c', p[0]), leg('A', p[1]), leg('lam', p[2])]
    if name == 'cclst':
        return [leg('c', p[0]), leg('c', p[1]), leg('lamst', p[2])]
    if name == 'ccAsW':
        return [leg('c', p[0]), leg('c', p[1]), leg('Ast', p[2], 'W')]


SECT_G = {'cAVAV': ['G1', 'G2'], 'cAWAW': ['G3'], 'cll': ['G5'], 'ccAs': ['G4']}
GF = {'G1': F_G1, 'G2': F_G2, 'G3': F_G3, 'G4': F_G4, 'G5': F_G5}
HF = {'H1': F_H1, 'H2': F_H2, 'H3': F_H3, 'H4': F_H4, 'H5': F_H5, 'H6': F_H6}
BF = {'D1': F_D1, 'B1': F_B1, 'B2': F_B2, 'B3': F_B3, 'B4': F_B4, 'AA': F_AA, 'AH': F_AH}


def to_fr(x):
    f = Fr(x).limit_denominator(10 ** 6)
    assert abs(float(f) - x) < 1e-9, x
    return f


def decompose(vals_by_sample, gvals_by_sample):
    """solve v = sum_i a_i G_i over samples (complex), exact rational coefficients, check residual."""
    A = np.array(gvals_by_sample, dtype=complex)
    b = np.array(vals_by_sample, dtype=complex)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    coef = []
    for s in sol:
        assert abs(s.imag) < 1e-8, s
        coef.append(to_fr(s.real))
    resid = np.max(np.abs(A @ np.array([float(c) for c in coef]) - b)) if len(b) else 0.0
    return coef, float(resid), float(np.max(np.abs(b)))


NS = 5
samples = {sec: [sector_legs(sec) for _ in range(NS)] for sec in SECT_G}
Gvals = {sec: [[value(GF[g], L) for g in SECT_G[sec]] for L in samples[sec]] for sec in SECT_G}
out = {'G_closed_forms': {}, 'QH': {}, 'brackets_I3': {}, 'extra': {}}

# closed forms of Prop s2c-dims
cf = []
for L in samples['cAVAV']:
    (_, X0, p0), (_, a1, p1), (_, a2, p2) = L
    # f = ([X0, X1], X2): legs carry rank-n colour; use the full tensor contraction
    def fval(K):
        return np.einsum('abc,ma,nb->mn', EPS, np.ones((1, 3)) * 0 + X0[None, :], a1) if False else None
    # direct: G1 = sum ([X0, a1_m], a2_n) * K_mn
    T = np.einsum('abc,a,mb,nc->mn', EPS, X0, a1, a2)    # ([X0, a1_m], a2_n)
    p1v, p2v = p1, p2
    G1c = np.trace(T[:4, :4]) * (p2v @ p2v - p1v @ p1v) - p2v @ T[:4, :4] @ p2v + p1v @ T[:4, :4] @ p1v
    G2c = p2v @ T[:4, :4] @ p2v - p1v @ T[:4, :4] @ p1v
    v1 = value(F_G1, L); v2 = value(F_G2, L)
    cf.append([abs(v1 - G1c), abs(v2 - G2c)])
out['G_closed_forms']['cAVAV_max_dev'] = float(np.max(cf))
cf = []
for L in samples['cll']:
    (_, X0, p0), (_, l1, p1), (_, l2, p2) = L
    # G5 = -i f (lam1, rho(p0) lam2) with colour: sum ([X0, e_b], e_c) l1_{sb} (A(p0) l2)_{sc}
    G5c = -1j * np.einsum('abc,a,sb,sc->', EPS, X0, l1, np.einsum('st,tc->sc', Amat(p10(p0)), l2))
    cf.append(abs(value(F_G5, L) - G5c))
out['G_closed_forms']['cll_max_dev'] = float(np.max(cf))
cf = []
for L in samples['ccAs']:
    (_, X0, p0), (_, X1, p1), (_, b, p2) = L
    G4c = -1j * np.einsum('abc,a,b,mc->m', EPS, X0, X1, b)[:4] @ p2
    cf.append(abs(value(F_G4, L) - G4c))
out['G_closed_forms']['ccAs_max_dev'] = float(np.max(cf))
print('closed forms', out['G_closed_forms'], flush=True)

# Q H_j and Q of the lambda-lambda-A terms
for name, F in list(HF.items()) + [('YV', F_YV), ('YW', F_YW)]:
    row = {}
    for sec in SECT_G:
        vals = [Qvalue(F, L) for L in samples[sec]]
        coef, res, mx = decompose(vals, Gvals[sec])
        row[sec] = {'coef': [str(c) for c in coef], 'resid': res, 'max': mx}
    out['QH'][name] = row
    print('Q', name, {s: row[s]['coef'] for s in row}, flush=True)

# {I_3, B}
for name, F in BF.items():
    row = {}
    for sec in SECT_G:
        vals = [bracket_value(F_I3, 0, 3, F, 2, L) for L in samples[sec]]
        coef, res, mx = decompose(vals, Gvals[sec])
        row[sec] = {'coef': [str(c) for c in coef], 'resid': res, 'max': mx}
    out['brackets_I3'][name] = row
    print('{I3,', name, '}', {s: row[s]['coef'] for s in row}, flush=True)

# no-invariant sectors: Q H_j and {I_3,B} vanish there
noinv = {}
for sec in ['cAVAW', 'cAl', 'cclst', 'ccAsW']:
    mx = 0.0
    for _ in range(2):
        L = sector_legs(sec)
        for F in list(HF.values()) + [F_YV, F_YW]:
            mx = max(mx, abs(Qvalue(F, L)))
        for F in BF.values():
            mx = max(mx, abs(bracket_value(F_I3, 0, 3, F, 2, L)))
    noinv[sec] = mx
out['no_invariant_sectors_max'] = noinv
print('no-invariant sectors', noinv, flush=True)

# degree-2 sector (c, c, A_V): Q G_i, {I_3, int(d*dc, d*A_V)}
d2 = {}
L2 = [sector_legs('ccAV') for _ in range(4)]
QG = {g: [Qvalue(GF[g], L) for L in L2] for g in GF}
br_E2 = [bracket_value(F_I3, 0, 3, F_E2, 2, L) for L in L2]
d2['max|QG_i|'] = {g: float(np.max(np.abs(QG[g]))) for g in GF}
d2['max|{I3,E2} - QG2|'] = float(np.max(np.abs(np.array(br_E2) - np.array(QG['G2']))))
d2['max|{I3,E2}|'] = float(np.max(np.abs(br_E2)))
out['degree2'] = d2
print('degree 2', d2, flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rb_plugback_maps.json'), 'w'), indent=1)
