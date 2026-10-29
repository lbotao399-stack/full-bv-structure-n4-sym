"""The two-leg defect D_2(eps) in the sector (c, A|_V): special bubbles - Q I^CT_2, as a series in eps.
Each special bubble is a power series in eps (exponent t1 t2 |p|^2/(t1+t2) <= eps |p|^2 on the special region, Lemma
c4-special); here the coefficients of eps^-1, eps^0, eps^1 are computed exactly (masters J2), giving
D_2(eps) = i(p.a) [ 0 * eps^-1 - (27/8) |p|^2 + d_1 eps |p|^4 + O(eps^2) ]  (c at p, A|_V at -p)."""
import json, sympy as sp
import ia_engine as IA, masters as M, loops as L
q = L.q
def coeffs(Num, curs, s, ns=(-1, 0, 1)):
    E = len(curs)
    ks = [q[m] - curs[s][m] for m in range(4)]
    Nt = sp.expand(-Num * sum(x * x for x in ks))
    shift = [sum(IA.BET[e] * curs[e][m] for e in range(E)) for m in range(4)]
    C = L.gaussian_moments(Nt, shift)
    rho = sp.expand(sum(IA.BET[e] * sum(x * x for x in curs[e]) for e in range(E)) - sum(x * x for x in shift))
    order = [s] + [e for e in range(E) if e != s]
    out = {}
    for n in ns:
        acc = {}
        for k, Ck in C.items():
            j = k + 2 - E + n
            if j < 0:
                continue
            expr = sp.expand(Ck * rho ** j) * sp.Rational((-1) ** j, 2 ** k * sp.factorial(j))
            if expr == 0:
                continue
            P = sp.Poly(expr, *IA.BET[:E]); D = P.total_degree(); sb = sum(IA.BET[:E])
            expr = sp.expand(sum(cf * sp.Mul(*[IA.BET[i] ** e for i, e in enumerate(mon)]) * sb ** (D - sum(mon)) for mon, cf in P.terms()))
            for mon, cf in sp.Poly(expr, *IA.BET[:E]).terms():   # homogenised with sum(beta) = 1
                al = [mon[e] for e in order]
                J = M.J2(al[0], al[1], E - n + sum(al))
                for kk, v in J.items():
                    acc[kk] = acc.get(kk, 0) + cf * sp.Rational(v.numerator, v.denominator)
        out[n] = {kk: sp.expand(v) for kk, v in acc.items() if sp.expand(v) != 0}
    return out
p = [sp.Integer(x) for x in (1, -2, 3, 1)]
a = [sp.Integer(x) for x in (2, 1, 1, 3)]
legs = [dict(sp='c', u=[sp.Integer(1)], p=p), dict(sp='A', u=a + [sp.Integer(0)] * 6, p=[-x for x in p])]
terms = []
IA.evaluate(legs, collect=terms)
pa = sum(x * y for x, y in zip(p, a)); p2 = sum(x * x for x in p)
tot = {}
for t in terms:
    c = coeffs(t['Num'], t['curs'], t['s'])
    for n, d in c.items():
        for b, v in d.items():
            tot.setdefault(n, {})[b] = sp.expand(tot.get(n, {}).get(b, 0) + t['fac'] * v)
res = {}
for n, d in tot.items():
    res[str(n)] = {b: str(sp.nsimplify(v / (sp.I * pa * p2 ** (n + 1)))) for b, v in d.items()}
res['-1']['1'] = str(sp.nsimplify(sp.sympify(res['-1'].get('1', '0')) + sp.Rational(9, 2)))   # + (-Q I^CT_2) = +9/2 eps^-1
res['note'] = 'coefficient of i(p.a)|p|^{2n+2} eps^n in D_2 = special bubbles - Q I^CT_2 (c at p, A|_V at -p), unit*(c_0,c_1)'
print(json.dumps(res, indent=1))
json.dump(res, open('d2_series.json', 'w'), indent=1)
