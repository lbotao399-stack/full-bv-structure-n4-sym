"""Exactness of the cubic obstruction: solve Ob_3 = Q F_3 + {I_3, F_2} exactly (sector by sector, on G1..G5).
F_2 = (27/16) int(d^*A|_V, d^*A|_V) + beta_1 int(dA|_V,dA|_V) + beta_2 int(dA|_Wperp,dA|_Wperp) + beta_3 int(lambda,dslash lambda)
      + beta_4 int(dc,A^*);   F_3 = sum_i h_i H_i,
H1 = int(dA|_V, [A|_V ^ A|_V]) (= 2 sum_{mu,nu} (d_mu A_nu, [A_mu,A_nu])), H2 = int(dA|_Wperp, [A|_V, A|_Wperp]) (= sum_{mu,w}
(d_mu iota_w A, [A_mu, iota_w A])), H3 = int([A|_V,c],A^*|_V), H4 = int([A|_Wperp,c],A^*|_Wperp), H5 = int([c,lambda],lambda^*),
H6 = int([c,c],c^*).  2-form pairing (B,C) = sum_{mu<nu} B_{mu nu} C_{mu nu} (as in -1/4 (dA,dA) of the action)."""
import json, sys, sympy as sp
import fa, sectors3 as S, vertices as V, others as O
I = sp.I
dot4 = V.dot4
fa.add_term('H1', 1, ('A', 'A', 'A'), lambda s: 2 * sum(I * s[0]['p'][m] * s[0]['u'][n] * s[1]['u'][m] * s[2]['u'][n] for m in range(4) for n in range(4)), ('f', (0, 1, 2)))
fa.add_term('H2', 1, ('A', 'A', 'A'), lambda s: sum(I * s[0]['p'][m] * s[0]['u'][w] * s[1]['u'][m] * s[2]['u'][w] for m in range(4) for w in range(4, 10)), ('f', (0, 1, 2)))
fa.add_term('H3', 1, ('A', 'c', 'As'), lambda s: sum(s[0]['u'][m] * s[2]['u'][m] for m in range(4)) * s[1]['u'][0], ('f', (0, 1, 2)))
fa.add_term('H4', 1, ('A', 'c', 'As'), lambda s: sum(s[0]['u'][w] * s[2]['u'][w] for w in range(4, 10)) * s[1]['u'][0], ('f', (0, 1, 2)))
fa.add_term('H5', 1, ('c', 'lam', 'lams'), lambda s: s[0]['u'][0] * (sp.Matrix(s[1]['u']).T * V.Mp * sp.Matrix(s[2]['u']))[0], ('f', (0, 1, 2)))
fa.add_term('H6', 1, ('c', 'c', 'cs'), lambda s: s[0]['u'][0] * s[1]['u'][0] * s[2]['u'][0], ('f', (0, 1, 2)))
def dA(x, M, N):
    v = 0
    if M < 4: v += I * x['p'][M] * x['u'][N]
    if N < 4: v -= I * x['p'][N] * x['u'][M]
    return v
fa.add_term('B1', 1, ('A', 'A'), lambda s: sum(dA(s[0], m, n) * dA(s[1], m, n) for m in range(4) for n in range(m + 1, 4)), ('k', (0, 1)))
fa.add_term('B2', 1, ('A', 'A'), lambda s: sum(dA(s[0], m, w) * dA(s[1], m, w) for m in range(4) for w in range(4, 10)), ('k', (0, 1)))
fa.add_term('dsAdsA', 1, ('A', 'A'), lambda s: (I * dot4(s[0]['p'], s[0]['u'])) * (I * dot4(s[1]['p'], s[1]['u'])), ('k', (0, 1)))
H = {h: fa.local([(1, h)], h, 0) for h in ('H1', 'H2', 'H3', 'H4', 'H5', 'H6')}
F2 = {'D': fa.local([(1, 'dsAdsA')], 'int(d*A,d*A)', 0), 'B1': fa.local([(1, 'B1')], 'int(dA_V,dA_V)', 0),
      'B2': fa.local([(1, 'B2')], 'int(dA_perp,dA_perp)', 0), 'B3': fa.local([(2, 'kinlam')], 'int(lam,dslash lam)', 0),
      'B4': fa.local([(1, 'kingh')], 'int(dc,A*)', 0)}
items = {f'Q{h}': fa.Q(F) for h, F in H.items()}
items.update({'{I3,%s}' % k: fa.bracket(fa.I3, F) for k, F in F2.items()})
NC = {'cAVAV': 4, 'cApAp': 3, 'clamlam': 4, 'ccAsV': 3}
proj = {}
for sec in S.SECTORS:
    rows = {k: [] for k in items}
    for c in range(NC[sec]):
        legs = S.config(sec, c)
        gv = {g: S.G_value(g, legs) for g in S.SECTOR_G[sec]}
        for k, F in items.items():
            rows[k].append((O.fval(F, legs), gv))
    proj[sec] = {}
    for k, rr in rows.items():
        sol, res = S.solve_exact(rr, S.SECTOR_G[sec])
        proj[sec][k] = ({g: v for g, v in sol.items()}, max(abs(sp.N(x)) for x in res))
    print(sec, {k: ({g: str(v) for g, v in d.items()}, str(r)) for k, (d, r) in proj[sec].items()}, flush=True)
# alpha from cubic_results (I-a totals)
L32 = sp.log(sp.Rational(3, 2))
alpha = {'G1': sp.Rational(11, 4) + sp.Rational(19, 4) * L32, 'G2': sp.Rational(27, 8), 'G3': sp.Rational(11, 4) + sp.Rational(19, 4) * L32,
         'G4': sp.Integer(0), 'G5': -9 * L32}
hs = sp.symbols('h1:7'); bs = sp.symbols('beta1:5')
F2coef = {'D': sp.Rational(27, 16), 'B1': bs[0], 'B2': bs[1], 'B3': bs[2], 'B4': bs[3]}
eqs = []
for sec in S.SECTORS:
    for g in S.SECTOR_G[sec]:
        rhs = sum(hs[i] * proj[sec][f'QH{i + 1}'][0].get(g, 0) for i in range(6)) + sum(F2coef[k] * proj[sec]['{I3,%s}' % k][0].get(g, 0) for k in F2)
        eqs.append(sp.Eq(alpha[g], rhs))
sol = sp.solve(eqs, list(hs), dict=True)
print('equations:', eqs)
print('solution (h in terms of beta):', sol)
json.dump(dict(projections={sec: {k: dict(G={g: str(v) for g, v in d.items()}, residual=str(r)) for k, (d, r) in proj[sec].items()} for sec in proj},
               equations=[str(e) for e in eqs], solution=[{str(k): str(v) for k, v in s.items()} for s in sol]),
          open('solve_F3.json', 'w'), indent=1)
