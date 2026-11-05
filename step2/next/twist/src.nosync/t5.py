import env, sector as S, twdriver as D
from fractions import Fraction as Fr
T = S.T
def h(k): return tuple(Fr(x, 2) for x in k)
gch = []; bch = []
for i, t in enumerate(T.taus):
    wt = tuple(a + b for a, b in zip(t, T.z))
    gch.append(S.wt5((0, 0), tuple(Fr(x, 2) for x in wt)))     # weight of null vector (entries +-2/2 -> +-1)
    bch.append(tuple(-x for x in h(T.u + t)))
    print('i', i, 'tau', t, 'gamma charge', gch[-1], 'beta charge', bch[-1])
Q = S.QXI
print('q_xi', Q)
def add(*cs): return tuple(sum(x) for x in zip(*cs))
for i in range(3):
    for j in range(3):
        for k in range(j + 1, 3):
            d = tuple(a - b for a, b in zip(add(gch[j], gch[k]), add(bch[i], Q)))
            if all(x == 0 for x in d): print('Q0 beta(%d) ~ [gamma(%d), gamma(%d)]' % (i, j, k))
for i in range(3):
    for j in range(3):
        print('beta(%d)gamma(%d) charge' % (i, j), add(bch[i], gch[j]), 'k =', S.in_lattice(add(bch[i], gch[j], tuple(-x for x in S.wt5((-1,1))))))
