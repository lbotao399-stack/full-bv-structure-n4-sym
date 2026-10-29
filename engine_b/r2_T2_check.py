"""Independent check of T2: Q W_1 in the (c, A|_V) sector equals -V^{AA}(leg at +p with polarisation i p, leg at -p with polarisation a),
V^{AA} the regulated one-loop A|_V self-energy (round-1 words, no special edges).  Exact eps -> 0 expansion of each bubble:
  X = (16pi^2)^{-1} int_0^1 dx sum_k 2^{-k} C_k(x) int_{eps/min(x,1-x)}^infty T^{-1-k} e^{-z T} dT,  z = x(1-x)p^2,
  k=0: E_1(w), k=1: z Gamma(-1,w), w = p^2 eps max(x,1-x)  ->
  X = eps^{-1} int (1/2)C_1 min(x,1-x) + int (-C_0 + z C_1/2)(log(p^2 eps) + gamma_E) + int [(-C_0 + z C_1/2) log max(x,1-x) - z C_1/2] + o(1).
Tadpoles: X = C_0/eps exactly."""
import json
import sympy as sp
import vertices as V
import loops as L
import sectors as S

I = sp.I
x = sp.Symbol('x')
p = S.p; q = L.q
a = sp.symbols('a0:4')
def leg_long():   # leg at +p, polarisation i p (A_V)
    return dict(sp='A', u=[I * p[m] for m in range(4)] + [0] * 6, p=list(p))
def leg_a():
    return dict(sp='A', u=list(a) + [0] * 6, p=[-v for v in p])
p2 = sum(v * v for v in p); pa = sum(p[m] * a[m] for m in range(4))
def int01_logmax(expr):
    expr = sp.expand(expr)
    return sp.integrate(expr * sp.log(1 - x), (x, 0, sp.Rational(1, 2))) + sp.integrate(expr * sp.log(x), (x, sp.Rational(1, 2), 1))
def int01_min(expr):
    return sp.integrate(expr * x, (x, 0, sp.Rational(1, 2))) + sp.integrate(expr * (1 - x), (x, sp.Rational(1, 2), 1))
tot = {'pole': 0, 'log': 0, 'finite': 0}
per = {}
masks = {'gauge loop A|_V': {('AV', 'AV')}, 'scalar loop A|_Wperp': {('Ap', 'Ap')}, 'ghost loop': {('c', 'As'), ('As', 'c')}, 'fermion loop': {('lam', 'lam')}}
e1 = leg_long(); e2 = leg_a()
for name, mask in masks.items():
    pole = 0; logc = 0; fin = 0
    for (x1, x2, px, sgn) in [(e1, e2, list(p), 1), (e2, e1, [-v for v in p], 1)]:
        h1 = L.full_hessian(V.CUBIC, [x1], [-v for v in q], [q[m] - px[m] for m in range(4)], S.CUBW)
        h2 = L.full_hessian(V.CUBIC, [x2], [px[m] - q[m] for m in range(4)], list(q), S.CUBW)
        Num = L.str_product([L.prop_numerator(list(q), mask), L.PiH(h1), L.prop_numerator([q[m] - px[m] for m in range(4)], mask), L.PiH(h2)])
        if Num == 0:
            continue
        b = sp.Symbol('beta1')
        C = L.gaussian_moments(Num, [b * px[m] for m in range(4)])
        C0 = sp.expand(sp.sympify(C.get(0, 0)).subs(b, x)); C1 = sp.expand(sp.sympify(C.get(1, 0)).subs(b, x))
        z = x * (1 - x) * p2
        pole += sgn * int01_min(C1 / 2)
        logc += sgn * sp.integrate(-C0 + z * C1 / 2, (x, 0, 1))
        fin += sgn * (int01_logmax(-C0 + z * C1 / 2) - sp.integrate(z * C1 / 2, (x, 0, 1)))
    w = sp.Rational(1, 4)
    per[name] = dict(pole=sp.simplify(w * pole / (I * pa)), log=sp.simplify(w * logc / (I * pa)), finite=sp.simplify(w * fin / (I * pa * p2)))
    tot['pole'] += w * pole; tot['log'] += w * logc; tot['finite'] += w * fin
# quartic tadpoles: pure eps^{-1}
for name, mask in {'quartic tadpole A|_V': {('AV', 'AV')}, 'quartic tadpole A|_Wperp': {('Ap', 'Ap')}}.items():
    h = L.full_hessian(V.QUARTIC, [e1, e2], [-v for v in q], list(q), S.QUAW)
    Num = L.str_product([L.prop_numerator(list(q), mask), L.PiH(h)])
    C = L.gaussian_moments(Num, [0] * 4)
    assert set(C.keys()) <= {0}
    val = sp.Rational(1, 2) * C.get(0, 0)
    per[name] = dict(pole=sp.simplify(val / (I * pa)), log=0, finite=0)
    tot['pole'] += val
out = {k: dict(pole=str(v['pole']), log=str(v['log']), finite=str(v['finite'])) for k, v in per.items()}
out['V_AA_longitudinal_total'] = dict(pole_over_i_pa=str(sp.simplify(tot['pole'] / (I * pa))), log_over_i_pa=str(sp.simplify(tot['log'] / (I * pa))),
                                      finite_over_i_pa_p2=str(sp.simplify(tot['finite'] / (I * pa * p2))))
out['QW1 = -V_AA(ip, a): prediction for the D-graph sum'] = dict(eps_inv=str(sp.simplify(-tot['pole'] / (I * pa))) + ' i(p.a)',
                                                                 finite=str(sp.simplify(-tot['finite'] / (I * pa * p2))) + ' i(p.a)p^2')
print(json.dumps(out, indent=1))
json.dump(out, open('r2_T2_check.json', 'w'), indent=1)
