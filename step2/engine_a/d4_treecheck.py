"""D4 engine A: sign of the tree insertion of lambda*int(dA|_P,dA|_P) into O = int(iota_w1 A, iota_w2 A) in the mode
model (C_P with P(0,infty)); Step-1 convention (ex:c6-bilinear): each A|_P leg of O multiplied by 4*lambda."""
import os, sys, json
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s2_core as C
import d4_signs as DS
from r2_core import Q2
p = np.array([1, 2, -1, 1])
M, pm, pairs = DS.build(p)
Pk = DS.Pkernel(M)
PC = range(4, 10)
# int(dA_P,dA_P) restricted: p^2 (A(p),A(-p)) per ordered mode pair
W = DS.quadP(M, 'AA', pairs, PC) * Fraction(int(p @ p))
acc = {}
for (m1, m2) in pairs:
    for a in range(C.D):
        for b in range(C.D):
            if C.G[a, b] == 0: continue
            M.add(acc, DS.q2(C.G[a, b]), [M.y('A', (4, a), m1), M.y('A', (5, b), m2)])
O = M.ga(acc)
t = M.C(Pk, W, O)
r = DS.ratio(t, O)
print(json.dumps({'C_P(int(dA_P,dA_P), O)/O (both legs; per leg = half)': str(r)}))
