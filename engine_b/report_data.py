"""Engine B: intermediate quantities for the formula chains of report.md.
For each two-point graph word (leg at p first) the Gaussian average over the loop momentum of the supertrace numerator,
   <STr[...]>(t) = C_1(x)/(2T) + C_0(x),  x = t_2/(t_1+t_2), T = t_1+t_2,
is fitted to index-free bilinears of the background fields; the pinned-edge (UV boundary) integrals are then written in
Feynman-parameter form:
   eps dX/deps = -(16pi^2)^{-1} { eps^{-1} [int_{1/2}^1 (1-x) C_1/2 dx + int_0^{1/2} x C_1/2 dx] + int_0^1 [C_0 - x(1-x) p^2 C_1/2] dx } + o(1),
which is checked against the direct u-integrals of loops.boundary_singular."""
import json
import sympy as sp
import vertices as V
import loops as L
import sectors as S
import o_graphs as OG

x = sp.Symbol('x')
p = S.p; q = L.q

def structs(e1, e2):
    u = e1['u']; v = e2['u']
    p2 = sum(pp * pp for pp in p)
    if e1['sp'] == 'A' and e2['sp'] == 'A':
        pu = sum(p[m] * u[m] for m in range(4)); pv = sum(p[m] * v[m] for m in range(4))
        uv = sum(u[M] * v[M] for M in range(10))
        return {'(A(-p),A(p))': uv, 'p^2 (A(-p),A(p))': p2 * uv, '(p.A(-p))(p.A(p))': pu * pv}
    if e1['sp'] == 'lam':
        val = (sp.Matrix(u).T * sum((p[m] * V.MR[m] for m in range(4)), sp.zeros(16, 16)) * sp.Matrix(v))[0]
        return {'(lam(-p), rho(p) lam(p))': sp.expand(val)}
    if e1['sp'] == 'As':
        return {'(p.A*(-p)) c(p)': sp.expand(sum(p[m] * u[m] for m in range(4)) * v[0])}

def fitC(expr, st):
    names = list(st)
    cs = sp.symbols(f'K0:{len(names)}')
    diff = sp.expand(expr - sum(c * st[n] for c, n in zip(cs, names)))
    syms = sorted(diff.free_symbols - set(cs) - {x}, key=lambda s: s.name)
    eqs = sp.Poly(diff, *syms).coeffs() if syms else [diff]
    sol = sp.solve([e for e in eqs if e != 0], cs, dict=True)
    if not sol and all(e == 0 for e in eqs):
        return {}
    return {n: sp.factor(sp.simplify(sol[0].get(c, 0))) for c, n in zip(cs, names) if sp.simplify(sol[0].get(c, 0)) != 0}

def word12(gid):
    for g in S.GRAPHS:
        if g[0] == gid:
            _, sector, kind, s1, s2, mask = g
    e1 = S.leg(s1[0], 'u', p, s1[1]); e2 = S.leg(s2[0], 'v', S.mp, s2[1])
    if kind == 'bubble':
        h1 = L.full_hessian(V.CUBIC, [e1], [-w for w in q], [q[m] - p[m] for m in range(4)], S.CUBW)
        h2 = L.full_hessian(V.CUBIC, [e2], [p[m] - q[m] for m in range(4)], list(q), S.CUBW)
        Num = L.str_product([L.prop_numerator(list(q), mask), L.PiH(h1), L.prop_numerator([q[m] - p[m] for m in range(4)], mask), L.PiH(h2)])
        lines = [[0] * 4, list(p)]
    else:
        h = L.full_hessian(V.QUARTIC, [e1, e2], [-w for w in q], list(q), S.QUAW)
        Num = L.str_product([L.prop_numerator(list(q), mask), L.PiH(h)]); lines = [[0] * 4]
    return e1, e2, kind, Num, lines

def chain(gid):
    e1, e2, kind, Num, lines = word12(gid)
    st = structs(e1, e2)
    if Num == 0:
        return dict(numerator='0')
    if kind == 'tadpole':
        C = L.gaussian_moments(Num, [0] * 4)
        return dict(word='STr[P(q) Pi d2I_4[A(p),A(-p)]]', C0=str(fitC(C.get(0, 0), st)), C1=str(fitC(C.get(1, 0), st)),
                    pinned='eps dX/deps = -(16pi^2)^{-1} C_0/eps')
    betas = sp.symbols('beta0:2')
    shift = [betas[1] * p[m] for m in range(4)]
    C = L.gaussian_moments(Num, shift)
    C0 = sp.expand(sp.sympify(C.get(0, 0)).subs({betas[1]: x, betas[0]: 1 - x}))
    C1 = sp.expand(sp.sympify(C.get(1, 0)).subs({betas[1]: x, betas[0]: 1 - x}))
    p2 = sum(pp * pp for pp in p)
    Am1 = sp.integrate((1 - x) * C1 / 2, (x, sp.Rational(1, 2), 1)) + sp.integrate(x * C1 / 2, (x, 0, sp.Rational(1, 2)))
    A0 = sp.integrate(C0 - x * (1 - x) * p2 * C1 / 2, (x, 0, 1))
    a_chk, b_chk = L.boundary_singular(Num, lines)
    ok = sp.simplify(Am1 / (16 * sp.pi**2) - a_chk) == 0 and sp.simplify(-A0 / (16 * sp.pi**2) - b_chk) == 0
    return dict(word='STr[P(q) Pi d2I_3[leg at p] P(q-p) Pi d2I_3[leg at -p]]',
                C1=str(fitC(C1, st)), C0=str(fitC(C0, st)),
                pole_integral=str(fitC(sp.expand(Am1), st)), log_integral_A0=str(fitC(sp.expand(A0), st)),
                feynman_form_matches_u_integrals=bool(ok))

if __name__ == '__main__':
    out = {}
    for g in S.GRAPHS:
        if g[0].endswith(':all') or g[1] == 'AV-Ap':
            continue
        out[g[0]] = chain(g[0])
        print(g[0], out[g[0]], flush=True)
    json.dump(out, open('report_data.json', 'w'), indent=1)
