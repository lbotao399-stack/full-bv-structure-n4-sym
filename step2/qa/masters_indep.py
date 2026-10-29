#!/usr/bin/env python3
"""QA, Step 2 Part C: independent evaluation of the Schwinger master integrals of the I-a graphs.

  J(a1, a2, a3; m) = int_0^1 du1 int_1^oo du2 int_1^oo du3 u1^a1 u2^a2 u3^a3 (u1+u2+u3)^-m     (three edges)
  J(a1, a2; m)     = int_0^1 du1 int_1^oo du2 u1^a1 u2^a2 (u1+u2)^-m                          (two edges)
(u1 = special edge).  Route used here (different from both engines): the last [1,oo) integral in closed form,
  int_1^oo u^a (s+u)^-m du = sum_k binom(a,k) (-s)^(a-k) (1+s)^(k-m+1)/(m-k-1),
then mpmath tanh-sinh quadrature at 40 digits in the remaining variables.  Compared with
  engine A: step2/engine_a/logs/s2_integrals.json (13 integrals J(n; b; N), a1 = 0),
  engine B: step2/engine_b/masters.py J3/J2 (exact closed forms) on the 16 masters of masters_check.json and on 40
            further random convergent index sets.  Writes qa/masters_indep.json."""
import json
import random
import sys
from fractions import Fraction as Fr
from math import comb

import mpmath as mp

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_b')
import masters as MB  # noqa: E402

mp.mp.dps = 40
OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/masters_indep.json'
L2, L3 = mp.log(2), mp.log(3)


def tail(a, m, s):
    return mp.fsum(comb(a, k) * (-s) ** (a - k) * (1 + s) ** (k - m + 1) / (m - k - 1) for k in range(a + 1))


def J3(a1, a2, a3, m):
    assert m > a2 + a3 + 2
    return mp.quad(lambda u1, u2: u1 ** a1 * u2 ** a2 * tail(a3, m, u1 + u2), [0, 1], [1, 2, 4, mp.inf])


def J2(a1, a2, m):
    assert m > a2 + 1
    return mp.quad(lambda u1: u1 ** a1 * tail(a2, m, u1), [0, 1])


def val3(c):  # (rational, log2, log3) strings or Fractions
    return mp.mpf(Fr(c[0]).numerator) / Fr(c[0]).denominator + (mp.mpf(Fr(c[1]).numerator) / Fr(c[1]).denominator) * L2 \
        + (mp.mpf(Fr(c[2]).numerator) / Fr(c[2]).denominator) * L3


def main():
    out = {'engine_A': {}, 'engine_B_listed': [], 'engine_B_random': []}
    A = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_a/logs/s2_integrals.json'))
    worst = mp.mpf(0)
    for key, rec in A.items():
        n, b, N = key.split('|')
        n, N = int(n), int(N)
        b = [int(x) for x in b.split(',')]
        q = J2(0, b[0], N) if n == 2 else J3(0, b[0], b[1], N)
        ex = val3(rec['exact'])
        d = abs(q - ex)
        worst = max(worst, d / max(abs(ex), mp.mpf('1e-30')))
        out['engine_A'][key] = {'exact_claimed': rec['exact'], 'quadrature_40d': mp.nstr(q, 30), 'abs_diff': mp.nstr(d, 3)}
    out['engine_A_max_rel_diff'] = mp.nstr(worst, 3)
    print('A max rel diff', out['engine_A_max_rel_diff'], flush=True)
    Bl = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_b/masters_check.json'))
    worstB = mp.mpf(0)

    def bnum(num):
        return MB.nfloat(num)

    for rec in Bl:
        args = rec['args']
        if rec['integral'] == 'J3':
            ex = bnum(MB.J3(*args))
            q = J3(*args)
        else:
            ex = bnum(MB.J2(*args))
            q = J2(*args)
        d = abs(q - ex) / abs(ex)
        worstB = max(worstB, d)
        out['engine_B_listed'].append({'integral': rec['integral'], 'args': args, 'exact_B': MB.nstr(
            MB.J3(*args) if rec['integral'] == 'J3' else MB.J2(*args)), 'quadrature_40d': mp.nstr(q, 30),
            'rel_diff': mp.nstr(d, 3)})
    rng = random.Random(2718)
    for _ in range(40):
        if rng.random() < 0.6:
            a1, a2, a3 = rng.randint(0, 3), rng.randint(0, 3), rng.randint(0, 3)
            m = a2 + a3 + 3 + rng.randint(0, 4)
            args, name = [a1, a2, a3, m], 'J3'
            ex = bnum(MB.J3(*args))
            q = J3(*args)
            s = MB.nstr(MB.J3(*args))
        else:
            a1, a2 = rng.randint(0, 4), rng.randint(0, 4)
            m = a2 + 2 + rng.randint(0, 4)
            args, name = [a1, a2, m], 'J2'
            ex = bnum(MB.J2(*args))
            q = J2(*args)
            s = MB.nstr(MB.J2(*args))
        d = abs(q - ex) / abs(ex)
        worstB = max(worstB, d)
        out['engine_B_random'].append({'integral': name, 'args': args, 'exact_B': s, 'quadrature_40d': mp.nstr(q, 30),
                                       'rel_diff': mp.nstr(d, 3)})
    out['engine_B_max_rel_diff'] = mp.nstr(worstB, 3)
    # cross-engine identity: engine A J(3; b2, b3; N) = engine B J3(0, b2, b3, N), J(2; b; N) = J2(0, b, N)
    cross = {}
    for key, rec in A.items():
        n, b, N = key.split('|')
        n, N = int(n), int(N)
        b = [int(x) for x in b.split(',')]
        numB = MB.J2(0, b[0], N) if n == 2 else MB.J3(0, b[0], b[1], N)
        fa = [Fr(x) for x in rec['exact']]
        fb = [numB.get('1', Fr(0)), numB.get('log2', Fr(0)), numB.get('log3', Fr(0))]
        cross[key] = {'A': [str(x) for x in fa], 'B': [str(x) for x in fb], 'identical': fa == fb}
    out['cross_engine_exact_identity'] = cross
    out['cross_engine_all_identical'] = all(v['identical'] for v in cross.values())
    print('B max rel diff', out['engine_B_max_rel_diff'], 'cross identical', out['cross_engine_all_identical'], flush=True)
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)


if __name__ == '__main__':
    main()
