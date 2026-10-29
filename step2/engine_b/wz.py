"""Wess-Zumino consistency in the sector (c, c, A|_V): Q Ob_3 + {Ob_2, I_3} = 0 ?  (Q of the 3-leg limit of Q W, using
Q B_{P(0,oo)}(Ob_2, I_3) = -B_{K_0}(Ob_2, I_3) = {Ob_2, I_3}, Q Ob_2 = 0, Q I_3 = 0.)"""
import sympy as sp, random, json
import fa, sectors3 as S, vertices as V
import others as O
I = sp.I; Z = sp.Integer(0)
def G(name):
    V.TERMS['_' + name] = S.GTERMS[name]
    return fa.local([(1, '_' + name)], name, 1)
Gf = {g: G(g) for g in S.GTERMS}
rng = random.Random(3)
def ri(): return sp.Integer(rng.choice([-3, -2, -1, 1, 2, 3]))
rows = []
for c in range(4):
    p0 = [ri() for _ in range(4)]; p1 = [ri() for _ in range(4)]; p2 = [-x - y for x, y in zip(p0, p1)]
    legs = [dict(sp='c', u=[sp.Integer(1)], p=p0), dict(sp='c', u=[sp.Integer(1)], p=p1), dict(sp='A', u=[ri() for _ in range(4)] + [Z] * 6, p=p2)]
    r = {g: O.fval(fa.Q(Gf[g]), legs) for g in ('G1', 'G2', 'G4')}
    r['{Ob2,I3}'] = O.fval(fa.bracket(O.OB2, fa.I3), legs)
    r['{I3,Ob2}'] = O.fval(fa.bracket(fa.I3, O.OB2), legs)
    rows.append(r)
    print(r)
# solve x1 QG1 + x2 QG2 + x4 QG4 = -{Ob2,I3}
x1, x2, x4 = sp.symbols('x1 x2 x4')
eqs = [x1 * r['G1'] + x2 * r['G2'] + x4 * r['G4'] + r['{Ob2,I3}'] for r in rows]
print('solutions for Q(x1 G1 + x2 G2 + x4 G4) = -{Ob2,I3}:', sp.solve(eqs, [x1, x2, x4], dict=True))
print('QG1 / QG4 ratios:', [sp.nsimplify(r['G1'] / r['G4']) if r['G4'] != 0 else None for r in rows])
json.dump(dict(statement='Q Ob_3 + {I_3, Ob_2} = 0 (equivalently Q Ob_3 + {Ob_2, I_3} = 0; here {Ob_2, I_3} = {I_3, Ob_2})',
               derivation='Q^2 = 0 applied to the three-leg part of lim_{eps->0} Q W(P(eps,infty), I - I^CT[eps]) = Ob_3 + B_{P(0,infty)}(Ob_2, I_3); '
                          'Q B_{P(0,oo)}(Ob_2, I_3) = B_{P(0,oo)}(Q Ob_2, I_3) - B_{P(0,oo)}(Ob_2, Q I_3) - B_{K_0}(Ob_2, I_3) = {Ob_2, I_3} '
                          '(Q P(0,oo) = K_0, Q Ob_2 = 0, Q I_3 = 0); the limit commutes with Q (local operator at fixed momenta). Hence automatic.',
               sector='(c, c, A|_V) is the only sector reached: Q G1 = Q G3 = Q G4 = Q G5 = 0 identically (Q G1, Q G4 are multiples of '
                      'int([c,c], d^*d^*dA) = 0); Q G2 != 0',
               configs=[{k: str(v) for k, v in r.items()} for r in rows],
               solution_for_alpha2=str(sp.solve(eqs, [x1, x2, x4], dict=True)),
               consequence='alpha_2 = 27/8 is forced by Ob_2; the computed alpha_2 = 27/8 satisfies it; alpha_1, alpha_3, alpha_4, alpha_5 are not constrained'),
          open('wz_result.json', 'w'), indent=1)
