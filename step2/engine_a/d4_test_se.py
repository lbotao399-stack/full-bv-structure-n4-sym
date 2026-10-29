import os, sys, time, json
import numpy as np
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s2_core as C
import d4_two_point as TP
from r2_core import Q2
rng = np.random.default_rng(5)
p = np.array([1, 1, 1, -1])
a = np.zeros((10, C.D)); a[4] = [1, 2, -1]
b = np.zeros((10, C.D)); b[4] = [2, 0, 1]
ab = float(np.einsum('ma,mb,ab->', a, b, C.G.real))
t0 = time.time()
se = TP.self_energy([('A', a), ('A', b)], p)
print('A_P self-energy, (a,b) =', ab, 'p^2 =', p @ p, 'time %.1f' % (time.time() - t0))
for k, v in se.items():
    print('  ', k, v, ' per (a,b):', v * Q2(Fraction(1, int(ab))) if ab else '')
