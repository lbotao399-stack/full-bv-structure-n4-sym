"""Engine A, part Q, task TQ2: tables of Upsilon_xi (eps^0 part of the defect cycles, and their poles) for the pairs of
field and ghost letters of Table tab:qs-fieldpairs, every leg set of census_qsk.json at n letter derivatives, every
ordering of the outputs (colour word ad_{Y_r}..ad_{Y_1}), at seeded Gaussian-integer data (letters, polarisations,
momenta recorded).  xi in S_+ (x) W generic (seeded).  Values exact, units hbar/(16 pi^2), stripped colour."""
import itertools, json, os, time
from fractions import Fraction
import numpy as np
import qa_s1 as S1, qa_core as QC, qa_util as U, ea_core as E

CENSUS = os.path.join(QC.S1.ROOT, 'step2', 'composite_qsk', 'qsk', 'theory', 'src.nosync', 'census_qsk.json')
PAIRS = [('AV', 'AV'), ('AV', 'AW'), ('AV', 'lamP'), ('AV', 'lamM'), ('AV', 'c'), ('AW', 'AW'), ('AW', 'lamP'),
         ('AW', 'lamM'), ('AW', 'c'), ('lamP', 'lamP'), ('lamP', 'lamM'), ('lamP', 'c'), ('lamM', 'lamM'), ('lamM', 'c'),
         ('c', 'c')]


def load_census(gen='xi'):
    return json.load(open(CENSUS))['generators'][gen]['pairs']


def cstr(v):
    return [str(complex(x)) for x in np.asarray(v).ravel()]


DERIVS = [{(0, 0, 0, 0): 1}, {(1, 0, 0, 0): 1j}, {(0, 1, 0, 0): 1j}]


def run(xi, n=0, pairs=PAIRS, seed=101, log=print, maxperm=6):
    cen = load_census('xi')
    R = U.Rng(seed + n)
    gx = QC.odd_gen(xi)
    rows = []
    for ta, tb in pairs:
        key = '%s|%s' % (ta, tb)
        ent = cen[key]['by_letter_derivatives'].get(str(n))
        if ent is None:
            continue
        W = Fraction(ent['W'])
        for ls, info in ent['defect_leg_sets'].items():
            types = [] if ls == '(none)' else ls.split()
            wout = Fraction(info['w_out'])
            D = W + Fraction(1, 2) - wout
            mmax = int(D // 2) if D >= 2 else 0
            za, zb = R.vec(ta), R.vec(tb)
            fa = DERIVS[1] if n >= 1 else None
            fb = DERIVS[2] if n >= 2 else None
            if n == 3:
                fa = {(1, 1, 0, 0): -1}
            Xa = (ta, za, fa); Xb = (tb, zb, fb)
            legs0 = [(t, R.vec(t), R.mom()) for t in types]
            perms = sorted(set(itertools.permutations(range(len(types)))))[:maxperm]
            what = ('fin',) + tuple(('pole', m) for m in range(1, mmax + 1))
            for pm in perms:
                legs = [legs0[i] for i in pm]
                QC.STATS['nonzero'] = 0
                v = QC.pair_value(Xa, Xb, legs, gx, what=what)
                rows.append({'pair': [ta, tb], 'n': n, 'W': str(W), 'leg_set': ls, 'ordering': [types[i] for i in pm],
                             'D': str(D), 'eps0': U.vstr(v['fin']),
                             'poles': {str(m): U.vstr(v[('pole', m)]) for m in range(1, mmax + 1)},
                             'nonzero_configs': QC.STATS['nonzero'], 'census_chains': info['chains'],
                             'census_kImin': info['k_I_min'],
                             'data': {'zeta_a': cstr(za), 'zeta_b': cstr(zb), 'deriv_a': str(fa), 'deriv_b': str(fb),
                                      'legs': [[l[0], cstr(l[1]), [int(x) for x in l[2]]] for l in legs]}})
        log('  TQ2 n=%d %s done (%d rows so far)' % (n, key, len(rows)))
    return rows
