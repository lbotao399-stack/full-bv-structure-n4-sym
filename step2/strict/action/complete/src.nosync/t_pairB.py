import numpy as np, random, time, resource
from fractions import Fraction as Fr
import cb_common as B, cb_s2 as S
E, Y, G = B.E, B.Y, B.G
random.seed(3)
def rspin(rng):
    x = np.zeros(16, complex)
    for c in rng: x[c] = random.randint(-2, 2)
    return x
xi1, xi2 = rspin(B.LAMP), rspin(B.LAMM)
tq, t3 = G.tmpl_xi()
t0 = time.time()
names = S.s2_terms(xi1, xi2, 'test')
tS2 = E.build_template(names, [], extra_odd=0, name=None)
print('S2 template', time.time() - t0, 's', list(tS2.keys()))
p = [Fr(random.randint(-3, 3)) for _ in range(4)]
u0 = np.zeros(16, complex); u1 = np.zeros(16, complex)
for c in B.LAMP: u0[c] = random.randint(-2, 2)
for c in B.LAMS_P: u1[c] = random.randint(-2, 2)
legs = [dict(sp='lam', u=u0, p=p), dict(sp='lams', u=u1, p=[-x for x in p])]
st = {}
t0 = time.time()
raw2 = Y.evaluate(legs, [Y.Gen('x1', 0, tq, gvec=xi1, t3=t3), Y.Gen('x2', 0, tq, gvec=xi2, t3=t3)], mlist=(0, 1), stats=st)
print('two-letter', time.time() - t0, 's', st, {k: {m: complex(v) for m, v in d.items()} for k, d in raw2.items()})
t0 = time.time()
raw1 = Y.evaluate(legs, [Y.Gen('S2', 0, tS2)], mlist=(0, 1), stats=st)
print('S2 letter', time.time() - t0, 's', {k: {m: complex(v) for m, v in d.items()} for k, d in raw1.items()})
print('RSS MB', resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20)
