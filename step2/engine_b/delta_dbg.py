import sympy as sp
import delta_engine as D, delta_canon as C
u1 = [sp.Integer(x) for x in [0]*4 + [-1, -3, -3, 2, -1, 1]]; u2 = [sp.Integer(x) for x in [0]*4 + [1, 2, -1, -2, 3, -2]]
p1 = [sp.Integer(x) for x in (-1, 1, -3, -2)]
legsB = [dict(sp='A', u=u1, p=list(p1)), dict(sp='A', u=u2, p=[-x for x in p1]), dict(sp='xi', u=[sp.Integer(1)], p=[0]*4)]
for part in ([('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC), ('SM', (2,), D.T_TERMS(0))],):
    recs = D.enumerate_terms(legsB, part, 'K', [0, 1])
    for r in recs:
        v = D.special_coeffs(r['Num'], r['curs'], r['s'])
        print([w[0]+str(w[1]) for w in r['word']], r['edges'], 's=', r['s'], 'fac', r['fac'], {n: v[n] for n in v})
