"""Audit #3: fit-free check of the exchange antisymmetry of engine A's Y1a, with engine A's own graph code.

Calls step2/engine_a/d4_y1a.y1a (read-only use; nothing written there) at the logged point
  p1 = (1,2,0,1), k = (0,-1,1,2)   (dependent leg p2 = (-1,-1,-1,-3))
and at the exchanged point p1' = p2 (dependent leg p2' = p1), with the polarisations of d4_y1a_ward.py.
For the kernel of Q(Y1a) (both legs A|_Wperp after Q), the beta part of the vertex value is
  beta_Q(p1) + beta_Q(p1')  with beta_Q(q) = -1/2 |q + k|^2 beta(q),
which must be the local 4k (units unit*T_O) if the nonlocal part is exchange-antisymmetric (check_y1a_symmetry.py).
This uses only raw graph values (no fitted formula).  su(2), mpmath dps 15 as in engine A.
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
EA = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_a'
sys.path.insert(0, EA)
os.environ.setdefault('S2_ALG', 'su2')
import s2_core as C
import d4_y1a as Y


def main():
    w1, w2 = 4, 5
    a = np.zeros((10, C.D)); a[w1] = C.traceless(np.array([1, 2, -1, 0, 1, 0, 2, 0, -1][:C.D], float))
    a[w2] = C.traceless(np.array([0, 1, 1, 2, 0, -1, 1, 1, 0][:C.D], float))
    b = np.zeros((10, C.D)); b[w1] = C.traceless(np.array([2, 0, 1, 1, 0, -1, 0, 1, 1][:C.D], float))
    b[w2] = C.traceless(np.array([1, -1, 2, 0, 1, 1, -1, 0, 2][:C.D], float))
    k = np.array([0, -1, 1, 2]); p1 = np.array([1, 2, 0, 1]); p2 = -k - p1
    t0 = time.time()
    r1 = Y.y1a(p1, k, a, b, w1, w2, dps=15)
    r2 = Y.y1a(p2, k, a, b, w1, w2, dps=15)
    be1 = [complex(x) / r1['T_O'] for x in r1['beta']]
    be2 = [complex(x) / r2['T_O'] for x in r2['beta']]
    n1, n2 = int(p1 @ p1), int(p2 @ p2)
    bQ1 = [-n2 / 2 * x for x in be1]          # first leg p1, second (Q-image) leg p2
    bQ2 = [-n1 / 2 * x for x in be2]          # first leg p2, second leg p1
    vv = [bQ1[m] + bQ2[m] for m in range(4)]
    out = {'p1': p1.tolist(), 'k': k.tolist(), 'T_O(p1)': r1['T_O'], 'T_O(p2)': r2['T_O'],
           'beta(Y1a)/T_O at p1': [str(x) for x in be1], 'beta(Y1a)/T_O at p1 swapped': [str(x) for x in be2],
           'alpha(Y1a)/T_O at p1, swapped': [str(complex(r1['alpha']) / r1['T_O']), str(complex(r2['alpha']) / r2['T_O'])],
           'VV beta_Q / (unit T_O) from raw graphs': [str(x) for x in vv], '4k': (4 * k).tolist(),
           'max |VV - 4k|': max(abs(vv[m] - 4 * k[m]) for m in range(4)),
           'antisym beta_Q / (unit T_O)': [str(bQ1[m] - bQ2[m]) for m in range(4)],
           'time_s': time.time() - t0}
    print(json.dumps(out, indent=1))
    json.dump(out, open(os.path.join(HERE, 'check_y1a_raw_swap.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
