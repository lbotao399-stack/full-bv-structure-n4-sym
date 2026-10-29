"""Regression: the I-a machinery on two legs (c, A|_V) must reproduce the Step-1 special bubbles (engine_b/r2_T2.py):
classes (a)..(e): eps^-1 coefficients of i(p.a): -9/4, 1/8, -9/2, 2, 1/8; eps^0 coefficients of i(p.a)|p|^2: -5/8, -1/12, -5/4, -4/3, -1/12."""
import sympy as sp, time
import ia_engine as IA
t0 = time.time()
p = [sp.Integer(x) for x in (1, -2, 3, 1)]
a = [sp.Integer(x) for x in (2, 1, 1, 3)]
legs = [dict(sp='c', u=[sp.Integer(1)], p=p), dict(sp='A', u=a + [sp.Integer(0)] * 6, p=[-x for x in p])]
res, nt = IA.evaluate(legs)
pa = sum(x * y for x, y in zip(p, a)); p2 = sum(x * x for x in p)
for cl, d in res.items():
    print(IA.class_name(cl))
    for n, num in d.items():
        print('   eps^%d' % n, {k: sp.nsimplify(v / (sp.I * pa * (p2 if n == 0 else 1))) for k, v in num.items()})
print('terms', nt, 'time', round(time.time() - t0, 1))
import json
json.dump(dict(statement='I-a machinery (special u in [0,1], others [1,oo), masters J2) on the two legs (c at p, A|_V at -p) reproduces Step-1 '
               'special bubbles (a)-(e) of the (c, A|_V) obstruction (computed in Step 1 by a different eps-expansion, engine_b/r2_T2.py)',
               classes={IA.class_name(cl): {('eps^%d' % n): {k: str(sp.nsimplify(v / (sp.I * pa * (p2 if n == 0 else 1)))) for k, v in num.items()} for n, num in d.items()} for cl, d in res.items()},
               expected={'eps^-1 / i(p.a)': '-9/4, 1/8, -9/2, 2, 1/8', 'eps^0 / i(p.a)|p|^2': '-5/8, -1/12, -5/4, -4/3, -1/12'}),
          open('regress_2pt.json', 'w'), indent=1)
