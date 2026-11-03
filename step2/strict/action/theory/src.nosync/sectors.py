"""Leg-set census for the action-level one-loop cochains Ob^(k), F^(k) along the generators of psl(4|4).
Necessary conditions only: ghost degree, Weyl weight (position at -1), position degree <= that of the generators,
derivative count D >= 0, colour (>= 2 adjoint legs), the centre of Spin(4) = SU(2)_L x SU(2)_R (Z2 x Z2) and the
centre Z4 of SU(4)_R.  Pure Python, exact fractions.  Output: sectors.json, sectors_tables.tex (next to this folder)."""
import itertools, json, os
from fractions import Fraction as Fr

# name: (degree, weight, nV, nS+, nS-, Z4, grassmann)
FIELDS = {
 'c':      (1, Fr(0),    0, 0, 0, 0, 1),
 'A_V':    (0, Fr(1),    1, 0, 0, 0, 0),
 'A_p':    (0, Fr(1),    0, 0, 0, 2, 0),
 'lam+':   (0, Fr(3,2),  0, 1, 0, 1, 1),
 'lam-':   (0, Fr(3,2),  0, 0, 1, 3, 1),
 'A*_V':   (-1, Fr(3),   1, 0, 0, 0, 1),
 'A*_p':   (-1, Fr(3),   0, 0, 0, 2, 1),
 'lam*+':  (-1, Fr(5,2), 0, 1, 0, 3, 0),
 'lam*-':  (-1, Fr(5,2), 0, 0, 1, 1, 0),
 'c*':     (-2, Fr(4),   0, 0, 0, 0, 0),
}
ORDER = list(FIELDS)
# generator: (ad_Delta eigenvalue, max position degree, nV of data, nS+, nS-, Z4)
GENS = {
 'Delta': (0, 1, 0, 0, 0, 0),
 'P_b':   (-2, 0, 1, 0, 0, 0),
 'M':     (0, 1, 2, 0, 0, 0),
 'R':     (0, 0, 0, 0, 0, 0),
 'xi+':   (-1, 0, 0, 1, 0, 1),
 'xi-':   (-1, 0, 0, 0, 1, 3),
 'eta+':  (1, 1, 0, 1, 0, 3),
 'eta-':  (1, 1, 0, 0, 1, 1),
 'K_b':   (2, 2, 1, 0, 0, 0),
}
PAIRS = [('Delta','P_b'), ('K_b','P_b'), ('xi+','xi-'), ('xi+','xi+'), ('xi+','eta+'), ('xi+','eta-'),
         ('xi+','P_b'), ('eta+','eta-'), ('eta+','eta+'), ('Delta','xi+')]

def census(gens, degree, maxlegs=6):
    lam = sum(GENS[g][0] for g in gens); pmax = sum(GENS[g][1] for g in gens)
    Wtot = 4 - Fr(lam, 2)
    dV = sum(GENS[g][2] for g in gens); dSp = sum(GENS[g][3] for g in gens)
    dSm = sum(GENS[g][4] for g in gens); dZ = sum(GENS[g][5] for g in gens)
    out = []
    for n in range(2, maxlegs + 1):
        for legs in itertools.combinations_with_replacement(ORDER, n):
            deg = sum(FIELDS[f][0] for f in legs)
            if deg != degree: continue
            wf = sum(FIELDS[f][1] for f in legs)
            for p in range(0, pmax + 1):
                D = Wtot + p - wf
                if D < 0 or D.denominator != 1: continue
                D = int(D)
                nV = sum(FIELDS[f][2] for f in legs) + D + p + dV
                nSp = sum(FIELDS[f][3] for f in legs) + dSp
                nSm = sum(FIELDS[f][4] for f in legs) + dSm
                z = (sum(FIELDS[f][5] for f in legs) + dZ) % 4
                if (nSp + nV) % 2 or (nSm + nV) % 2 or z: continue
                out.append({'legs': list(legs), 'n': n, 'p': p, 'D': D})
    return Wtot, pmax, out

def main():
    res = {'conventions': 'degree: c +1, A,lam 0, A*,lam* -1, c* -2; Weyl weight: c 0, A 1, lam 3/2, A* 3, lam* 5/2, c* 4, derivative 1, '
           'position -1; lam+ = S_+ (x) W, lam- = S_- (x) W*, lam*+ = S_+ (x) W*, lam*- = S_- (x) W; total weight 4 - (sum of ad_Delta '
           'eigenvalues)/2; D = number of derivatives; p = degree in the explicit position', 'arity1': {}, 'arity2': {}}
    for g in GENS:
        for name, deg in (('Ob', 0), ('F', -1)):
            W, pm, out = census([g], deg)
            res['arity1'][f'{name}^(1)({g})'] = {'weight': str(W), 'pmax': pm, 'sectors': out}
    for a, b in PAIRS:
        for name, deg in (('Ob', -1), ('F', -2)):
            W, pm, out = census([a, b], deg)
            res['arity2'][f'{name}^(2)({a},{b})'] = {'weight': str(W), 'pmax': pm, 'sectors': out}
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(res, open(os.path.join(here, '..', 'sectors.json'), 'w'), indent=1)
    for k, v in list(res['arity1'].items()) + list(res['arity2'].items()):
        s = v['sectors']
        print(k, 'weight', v['weight'], 'pmax', v['pmax'], 'count', len(s))
        for e in s:
            print('   p=%d n=%d D=%d  %s' % (e['p'], e['n'], e['D'], ' '.join(e['legs'])))
main()
