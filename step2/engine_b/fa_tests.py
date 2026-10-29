"""Checks of the tree-level algebra fa.py (field-count lemma, Q from the table, CME, Q-identity of a P-tree)."""
import random, time, json, sys
import sympy as sp
import fa
import vertices as V

I = sp.I; Z = sp.Integer(0)
rng = random.Random(7)
def ri(): return sp.Integer(rng.choice([-3, -2, -1, 1, 2, 3]))
def mom(): return [ri() for _ in range(4)]
def leg(sp_, p, comps=None):
    u = [Z] * V.DIM[sp_]
    for c in (range(V.DIM[sp_]) if comps is None else comps):
        u[c] = ri()
    if sp_ in ('c', 'cs'):
        u = [sp.Integer(1)]
    return dict(sp=sp_, u=u, p=p)
def eq(a, b):
    a = fa.reduce_K13(a); b = fa.reduce_K13(b)
    return all(sp.expand(a.get(k, 0) - b.get(k, 0)) == 0 for k in set(a) | set(b))
out = {}
t0 = time.time()
# 1. field count lemma
def count_test(Zf, weights, F, legs):
    n = {s: sum(1 for l in legs if l['sp'] == s) for s in ('A', 'As', 'lam', 'lams', 'c', 'cs')}
    w = weights(n)
    lhs = fa.bracket(Zf, F)(legs); rhs = {k: w * v for k, v in F(legs).items()}
    return eq(lhs, rhs), w
p0, p1 = mom(), mom(); p2 = [-x - y for x, y in zip(p0, p1)]
cases = [('AcAs', [leg('A', p0), leg('c', p1), leg('As', p2)]), ('clamlams', [leg('c', p0), leg('lam', p1), leg('lams', p2)]),
         ('cccs', [leg('c', p0), leg('c', p1), leg('cs', p2)]), ('AAA', [leg('A', p0), leg('A', p1), leg('A', p2)]),
         ('lamlamA', [leg('lam', p0), leg('A', p1), leg('lam', p2)])]
res = {}
for t, legs in cases:
    F = fa.I3_TERMS[t]
    res[t] = dict(ZA=count_test(fa.ZA, lambda n: n['A'] - n['As'], F, legs),
                  Zlam_raw=count_test(fa.Zlam_raw, lambda n: n['lams'] - n['lam'], F, legs),
                  Zc_raw=count_test(fa.Zc_raw, lambda n: n['cs'] - n['c'], F, legs))
out['field_count'] = {t: {k: str(v) for k, v in d.items()} for t, d in res.items()}
print('field count', out['field_count'], round(time.time() - t0, 1), flush=True)
# 2. Q from the table: Q int(d^*A_V, d^*A_V) = -2 int(d^*dc, d^*A_V);  Q int(dc, A*) on (c, A_V): -1/2 int(dc, d^*dA) = 0
fa.add_term('dsAdsA', 1, ('A', 'A'), lambda s: (I * V.dot4(s[0]['p'], s[0]['u'])) * (I * V.dot4(s[1]['p'], s[1]['u'])), ('k', (0, 1)))
fa.add_term('dsdc_dsA', 1, ('c', 'A'), lambda s: V.dot4(s[0]['p'], s[0]['p']) * s[0]['u'][0] * (-I * V.dot4(s[1]['p'], s[1]['u'])), ('k', (0, 1)))
F = fa.local([(1, 'dsAdsA')], 'int(d*A,d*A)', 0)
G = fa.local([(1, 'dsdc_dsA')], 'int(d*dc,d*A)', 1)
p = mom()
legs = [leg('c', p), leg('A', [-x for x in p], range(4))]
out['Q_table_dsAdsA'] = eq(fa.Q(F)(legs), {k: -2 * v for k, v in G(legs).items()})
print('Q int(d*A,d*A) = -2 int(d*dc,d*A):', out['Q_table_dsAdsA'], flush=True)
# 3. CME: Q I_3 = 0 and Q I_4 + 1/2 {I_3, I_3} = 0
cme3 = {}
for nm, legs in [('c,A,A', [leg('c', p0), leg('A', p1), leg('A', p2)]), ('c,lam,lam', [leg('c', p0), leg('lam', p1), leg('lam', p2)]),
                 ('c,c,As', [leg('c', p0), leg('c', p1), leg('As', p2)]), ('A,A,As', [leg('A', p0), leg('A', p1), leg('As', p2)]),
                 ('c,A,lams', [leg('c', p0), leg('A', p1), leg('lams', p2)])]:
    cme3[nm] = eq(fa.Q(fa.I3)(legs), {})
out['QI3=0'] = cme3
print('Q I_3 = 0:', cme3, round(time.time() - t0, 1), flush=True)
q0, q1, q2 = mom(), mom(), mom(); q3 = [-a - b - c for a, b, c in zip(q0, q1, q2)]
cme4 = {}
for nm, legs in [('c,A,A,A', [leg('c', q0), leg('A', q1), leg('A', q2), leg('A', q3)]),
                 ('c,A,lam,lam', [leg('c', q0), leg('A', q1), leg('lam', q2), leg('lam', q3)])]:
    lhs = fa.addd(fa.Q(fa.I4)(legs), {k: v / 2 for k, v in fa.bracket(fa.I3, fa.I3)(legs).items()})
    cme4[nm] = eq(lhs, {})
out['QI4+1/2{I3,I3}=0'] = cme4
print('Q I_4 + 1/2{I_3,I_3} = 0:', cme4, round(time.time() - t0, 1), flush=True)
# 4. Q-identity of a P-tree: Q B_P(F, I_3) = B_P(Q F, I_3) - B_{K_eps}(F, I_3)   (Q I_3 = 0, F even)
eps = sp.Rational(1, 3)
KP = fa.kernel_P(eps)
Keps = fa.kernel_scaled(fa.kernel_K0(), lambda q2: sp.exp(-eps * q2))
tree = fa.Bfun(KP, F, fa.I3, 'B_P(F,I3)')
legs = [leg('c', p0), leg('A', p1, range(4)), leg('A', p2, range(4))]
lhs = fa.Q(tree)(legs)
rhs = fa.addd(fa.B(KP, fa.Q(F), fa.I3, legs), fa.B(Keps, F, fa.I3, legs), -1)
out['Q B_P(F,I3) = B_P(QF,I3) - B_Keps(F,I3)'] = eq({k: sp.nsimplify(v) for k, v in lhs.items()}, {k: sp.nsimplify(v) for k, v in rhs.items()})
print('P-tree identity:', out['Q B_P(F,I3) = B_P(QF,I3) - B_Keps(F,I3)'], {k: sp.N(v) for k, v in lhs.items()}, {k: sp.N(v) for k, v in rhs.items()}, round(time.time() - t0, 1), flush=True)
json.dump(out, open('fa_tests.json', 'w'), indent=1, default=str)
