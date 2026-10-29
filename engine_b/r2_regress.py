"""Regression: the round-2 exact engine reproduces the round-1 omega=1 two-point log coefficients (lam-lam: 4 B6, As-c: -1 B7)."""
import random, time
from fractions import Fraction as Fr
import sympy as sp
import r2_engine as R
import sectors as S

print('24-cell 5-design:', R.check_design())
random.seed(5)
def ri(n, lo=-3, hi=3):
    return [random.randint(lo, hi) for _ in range(n)]
p = [1, -2, 3, 1]
for name, s1, s2, basis in [('lam-lam', ('lam', list(range(16))), ('lam', list(range(16))), 'B6'),
                            ('As-c', ('As', list(range(4))), ('c', [0]), 'B7')]:
    t0 = time.time()
    u1 = [0] * R.DIM[s1[0]]; u2 = [0] * R.DIM[s2[0]]
    for i in s1[1]: u1[i] = random.choice([-3,-2,-1,1,2,3])
    for i in s2[1]: u2[i] = random.choice([-3,-2,-1,1,2,3])
    legs = [R.Leg(0, s1[0], s1[1], u1, p), R.Leg(1, s2[0], s2[1], u2, [-x for x in p])]
    X = R.ExternalSet(name, legs, omega=1)
    recs = X.evaluate()
    tot = R.assemble(recs, legs)
    e1 = dict(sp=s1[0], u=[sp.Integer(x) for x in u1], p=[sp.Integer(x) for x in p])
    e2 = dict(sp=s2[0], u=[sp.Integer(x) for x in u2], p=[sp.Integer(-x) for x in p])
    bv = S.basis_value(basis, e1, e2)
    v = tot.get('cadj k', (0, 0))
    ratio = sp.nsimplify((sp.Rational(v[0].numerator, v[0].denominator) + sp.I * sp.Rational(v[1].numerator, v[1].denominator)) / bv)
    print(name, 'log coefficient / basis', ratio, ' (round 1:', {'B6': 4, 'B7': -1}[basis], ')', round(time.time() - t0, 1), 's', flush=True)
