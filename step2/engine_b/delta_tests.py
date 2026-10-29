"""Regression tests of delta_engine against Step-1 / Part C results (no new physics)."""
import json, time, random
import sympy as sp
import delta_engine as D
import masters as MS
import ia_engine as IA

out = {}
t0 = time.time()
# 1. masters
ok = True
for (a1, a2, m) in [(0, 0, 3), (1, 2, 5), (2, 0, 4), (0, 3, 6)]:
    J, r = D.master_special(a1, [a2], m)
    ref = MS.J2(a1, a2, m)
    ok &= (J['1'] == ref.get('1', 0) and J['L'] == ref.get('log2', 0))
for (a1, a2, a3, m) in [(0, 0, 0, 4), (1, 1, 0, 5), (0, 2, 1, 7), (2, 0, 0, 5)]:
    J, r = D.master_special(a1, [a2, a3], m)
    ref = MS.J3(a1, a2, a3, m)
    # J3 uses log3 - log2 = log(3/2)
    ok &= (J['1'] == ref.get('1', 0) and J['L'] == ref.get('log3', 0) and -J['L'] == ref.get('log2', 0))
out['masters_J2_J3'] = ok

# 2. special bubbles (c, A|_V): Step-1 values
p = [sp.Integer(x) for x in (1, -2, 3, 1)]
a = [sp.Integer(x) for x in (2, 1, 1, 3)]
legs = [dict(sp='c', u=[sp.Integer(1)], p=p), dict(sp='A', u=a + [sp.Integer(0)] * 6, p=[-x for x in p])]
part = [('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC)]
recs = D.enumerate_terms(legs, part, 'K', [0, 1])
pa = sum(x * y for x, y in zip(p, a)); p2 = sum(x * x for x in p)
tot = {-1: 0, 0: 0}
for r in recs:
    v = D.special_coeffs(r['Num'], r['curs'], r['s'])
    for n in (-1, 0):
        assert v[n]['L'] == 0
        tot[n] += r['fac'] * v[n]['1']
out['special_bubbles_c_AV'] = {'eps^-1 / i(p.a)': str(sp.nsimplify(tot[-1] / (sp.I * pa))),
                               'eps^0 / i(p.a)|p|^2': str(sp.nsimplify(tot[0] / (sp.I * pa * p2))),
                               'expected': '-9/2 (sum of a-e), -27/8 (sum of a-e)'}

# 3. log eps coefficients of the two-point functions (Theorem c4-quadratic): symbolic momenta
P = list(sp.symbols('P0:4'))
mP = [-x for x in P]
random.seed(3)
u = [0] * 4 + [random.randint(-3, 3) for _ in range(6)]
v = [0] * 4 + [random.randint(-3, 3) for _ in range(6)]
legsA = [dict(sp='A', u=[sp.Integer(x) for x in u], p=P), dict(sp='A', u=[sp.Integer(x) for x in v], p=mP)]
tot = 0
for part in ([('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC)], [('quartic', (0, 1), D.V.QUARTIC)]):
    for r in D.enumerate_terms(legsA, part, 'P', [0, 1]):
        tot += r['fac'] * D.log_coeff(r['Num'], r['curs'], 2, P)
uv = sum(x * y for x, y in zip(u, v)); PP = sum(x * x for x in P)
out['log_coeff_Aperp'] = {'value / (2|p|^2 (u,v))': str(sp.simplify(tot / (2 * PP * uv))), 'expected': '-1'}
# lambda sector
legsL = [dict(sp='lam', u=[sp.Integer(random.randint(-2, 2)) for _ in range(16)], p=P),
         dict(sp='lam', u=[sp.Integer(random.randint(-2, 2)) for _ in range(16)], p=mP)]
tot = 0
for r in D.enumerate_terms(legsL, [('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC)], 'P', [0, 1]):
    tot += r['fac'] * D.log_coeff(r['Num'], r['curs'], 1, P)
ref = D.V.vertex('kinlam', legsL)['k'] * 2      # 4 int(lambda, dslash lambda) = 8 * (1/2 int(...)), kinlam term = 1/2 int
out['log_coeff_lambda'] = {'value / plane-wave value of int(lambda,dslash lambda)': str(sp.simplify(tot / ref)), 'expected': '4'}
out['time_s'] = round(time.time() - t0, 1)
json.dump(out, open('delta_tests.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
