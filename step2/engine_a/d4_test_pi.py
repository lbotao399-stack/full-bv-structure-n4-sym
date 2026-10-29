import os, sys, time, json
import numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s2_core as C
import d4_two_point as TP
from r2_core import Q2
p = np.array([int(x) for x in sys.argv[1].split(',')]) if len(sys.argv) > 1 else np.array([1, 2, 0, 1])
a = np.zeros((10, C.D)); a[4] = C.traceless(np.array([1, 2, -1, 0, 1, 0, 2, 0, -1][:C.D], dtype=float))
b = np.zeros((10, C.D)); b[4] = C.traceless(np.array([2, 0, 1, 1, 0, -1, 0, 1, 1][:C.D], dtype=float))
ab = int(round(float(np.einsum('ma,mb,ab->', a, b, C.G.real))))
t0 = time.time()
alpha, beta, n0 = TP.pi_int('AsP,AP', p, a, b, dependent=(sys.argv[2] if len(sys.argv) > 2 else 'antifield'))
print('Pi_int (A^*_P,A_P): (a,b)=', ab, 'p=', p, 'p^2=', p @ p, 'nseq', n0, 'time %.1f' % (time.time() - t0))
print('alpha per (a,b):', {k: str(v * Q2(Fraction(1, ab))) for k, v in alpha.items()})
for mu in range(4):
    print('beta^%d per (a,b) / p_mu:' % mu, {k: str(v * Q2(Fraction(1, ab * int(p[mu])))) for k, v in beta[mu].items()} if p[mu] else {k: str(v) for k, v in beta[mu].items()})
