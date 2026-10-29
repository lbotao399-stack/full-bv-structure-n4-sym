"""Engine A, part C: W_0 tables of the composite kernels for the ten field pairs (leg sets with at most two outputs of
census_sp.json at the pair weight without derivatives), every ordered output word, exact values at seeded random off-shell
kinematic points (real integer momenta) and Gaussian-integer polarisations: C^(2),0_Delta (if the leg set enters the
translation-invariant part) and C^(2)_{P_beta}, each as letter part and vertex-vertex part, plus the poles z_1 of Z_{P_beta}."""
import itertools, json, os
from fractions import Fraction
import numpy as np
import cq_core as Q, cq_gval as GV, cq_eval as V

ONE = {(0, 0, 0, 0): 1}
RN = {'AV': 4, 'AW': 6, 'lamP': 16, 'lamM': 16, 'c': 1, 'cs': 1, 'AsV': 4, 'AsW': 6, 'lamsP': 16, 'lamsM': 16}
FIELD_PAIRS = ['AV|AV', 'AV|AW', 'AV|lamP', 'AV|lamM', 'AW|AW', 'AW|lamP', 'AW|lamM', 'lamP|lamP', 'lamP|lamM', 'lamM|lamM']


def rv(n, rng):
    return np.array([complex(int(x), int(y)) for x, y in zip(rng.integers(-2, 3, n), rng.integers(-2, 3, n))])


def run(census_path, npts=2, seed=17, beta=(1, 2, -1, 3), maxout=2, log=print):
    cen = json.load(open(census_path))['pairs']
    rng = np.random.default_rng(seed)
    out = {}
    for pname in FIELD_PAIRS:
        ta, tb = pname.split('|')
        ent = cen[pname]['by_letter_derivatives']['0']
        lsP = ent['leg_sets_Pb_and_position']; lsT = ent['leg_sets_translation_invariant']
        out[pname] = {'W0': cen[pname]['W0'], 'leg_sets': {}}
        for ls, info in lsP.items():
            types = ls.split()
            if len(types) > maxout:
                continue
            words = sorted(set(itertools.permutations(types)))
            rows = []
            for w in words:
                for it in range(npts):
                    za, zb = rv(RN[ta], rng), rv(RN[tb], rng)
                    moms = [rng.integers(-3, 4, 4) for _ in w]
                    legs = [(t, rv(RN[t], rng), m) for t, m in zip(w, moms)]
                    Xa, Xb = (ta, za, ONE), (tb, zb, ONE)
                    row = {'word': list(w), 'zeta_a': [str(x) for x in za], 'zeta_b': [str(x) for x in zb],
                           'legs': [{'type': t, 'pol': [str(x) for x in p], 'mom': [int(x) for x in m]} for t, p, m in legs]}
                    a, b = GV.kernel(Xa, Xb, legs, ('Pv', beta))
                    row['C2_Pbeta'] = {'letter': V.vstr(a), 'vertex_vertex': V.vstr(b)}
                    if ls in lsT:
                        a, b = GV.kernel(Xa, Xb, legs, 'D0')
                        row['C2_0_Delta'] = {'letter': V.vstr(a), 'vertex_vertex': V.vstr(b)}
                    rows.append(row)
            nz = any(r['C2_Pbeta']['letter'] != '0' or r['C2_Pbeta']['vertex_vertex'] != '0' for r in rows)
            out[pname]['leg_sets'][ls] = {'census_chains': info['chains'], 'k_min': info['k_min'], 'w_out': info['w_out'],
                                          'translation_invariant_part': ls in lsT, 'rows': rows, 'nonzero': nz}
            log('  table %s %s: %d rows, nonzero %s' % (pname, ls, len(rows), nz))
    return out
