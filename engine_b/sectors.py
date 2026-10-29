"""Engine B: one-loop two-point graphs and the quadratic counterterm I^CT[eps] (units hbar c_adj/(16 pi^2))."""
import json, sys, time
import sympy as sp
import vertices as V
import loops as L

I = sp.I
p = sp.symbols('p0:4')
mp = [-x for x in p]
q = L.q

def leg(spc, name, mom, comps=None):
    u = []
    for i in range(V.DIM[spc]):
        u.append(sp.Symbol(f'{name}{i}') if (comps is None or i in comps) else sp.Integer(0))
    return dict(sp=spc, u=u, p=list(mom))

AVc = list(range(4)); APc = list(range(4, 10))
CUBW = {'f': 1}                                    # bubble colour factor = c_adj (c1,c2)
QUAW = {'K12': 0, 'K13': -1, 'K14': -1}             # tadpole colour factors / (c_adj (c1,c2))

def bubble(e1, e2, mask):
    """(1/4)(X_12 + (-1)^{|1||2|} X_21): coefficient of r1 r2 hbar c_adj (c1,c2) in W_1; returns (a,b) in units 1/(16pi^2)."""
    p1 = e1['p']; p2 = e2['p']
    res_a = 0; res_b = 0
    for (x, y, px, py, sgn) in [(e1, e2, p1, p2, 1), (e2, e1, p2, p1, (-1) ** (V.PAR[e1['sp']] * V.PAR[e2['sp']]))]:
        h1 = L.full_hessian(V.CUBIC, [x], [-v for v in q], [q[m] - px[m] for m in range(4)], CUBW)
        h2 = L.full_hessian(V.CUBIC, [y], [px[m] - q[m] for m in range(4)], list(q), CUBW)
        N1 = L.prop_numerator(list(q), mask); N2 = L.prop_numerator([q[m] - px[m] for m in range(4)], mask)
        Num = L.str_product([N1, L.PiH(h1), N2, L.PiH(h2)])
        if Num == 0:
            continue
        a, b = L.boundary_singular(Num, [[0] * 4, list(px)])
        res_a += sgn * a; res_b += sgn * b
    f = 16 * sp.pi**2 * sp.Rational(1, 4)
    return sp.expand(res_a * f), sp.expand(res_b * f)

def tadpole(e1, e2, mask):
    h = L.full_hessian(V.QUARTIC, [e1, e2], [-v for v in q], list(q), QUAW)
    N = L.prop_numerator(list(q), mask)
    Num = L.str_product([N, L.PiH(h)])
    if Num == 0:
        return sp.Integer(0), sp.Integer(0)
    a, b = L.boundary_singular(Num, [[0] * 4])
    f = 16 * sp.pi**2 * sp.Rational(1, 2)
    return sp.expand(a * f), sp.expand(b * f)

# ---- counterterm basis (colour (c1,c2)); value = coefficient of r1 r2 on the plane-wave background ----
def _dA(x, M, N):
    v = 0
    if M < 4: v += I * x['p'][M] * x['u'][N]
    if N < 4: v -= I * x['p'][N] * x['u'][M]
    return v
BASIS = {
    'B1': (('A', 'A'), lambda s: sum(_dA(s[0], m, n) * _dA(s[1], m, n) for m in range(4) for n in range(m + 1, 4))),
    'B2': (('A', 'A'), lambda s: sum(I * s[0]['p'][m] * s[0]['u'][m] for m in range(4)) * sum(I * s[1]['p'][m] * s[1]['u'][m] for m in range(4))),
    'B3': (('A', 'A'), lambda s: sum(s[0]['u'][m] * s[1]['u'][m] for m in range(4))),
    'B4': (('A', 'A'), lambda s: sum(_dA(s[0], m, n) * _dA(s[1], m, n) for m in range(4) for n in range(4, 10))),
    'B5': (('A', 'A'), lambda s: sum(s[0]['u'][n] * s[1]['u'][n] for n in range(4, 10))),
    'B6': (('lam', 'lam'), lambda s: sum((sp.Matrix(s[0]['u']).T * V.MR[m] * sp.Matrix(s[1]['u']))[0] * I * s[1]['p'][m] for m in range(4))),
    'B7': (('c', 'As'), lambda s: sum(I * s[0]['p'][m] * s[0]['u'][0] * s[1]['u'][m] for m in range(4))),
}
def basis_value(name, e1, e2):
    slots, Lf = BASIS[name]
    V.TERMS['_b'] = (sp.Integer(1), slots, Lf, ('k', (0, 1)))
    v = V.vertex('_b', [e1, e2]).get('k', 0)
    del V.TERMS['_b']
    return sp.expand(v)

def fit(target, e1, e2, names):
    cs = sp.symbols(f'C0:{len(names)}')
    expr = sp.expand(target - sum(c * basis_value(n, e1, e2) for c, n in zip(cs, names)))
    syms = sorted((expr.free_symbols - set(cs)), key=lambda s: s.name)
    eqs = sp.Poly(expr, *syms).coeffs() if syms else [expr]
    sol = sp.solve(eqs, cs, dict=True)
    if not sol:
        return None
    return {n: sp.nsimplify(sol[0].get(c, c)) for c, n in zip(cs, names)}

GRAPHS = [
    # id, sector, kind, ext species/comps, mask
    ('AV-AV:AAA+AAA:AV,AV', 'AV-AV', 'bubble', ('A', AVc), ('A', AVc), {('AV', 'AV')}),
    ('AV-AV:AAA+AAA:Ap,Ap', 'AV-AV', 'bubble', ('A', AVc), ('A', AVc), {('Ap', 'Ap')}),
    ('AV-AV:AcAs+AcAs:c-As,c-As', 'AV-AV', 'bubble', ('A', AVc), ('A', AVc), {('c', 'As'), ('As', 'c')}),
    ('AV-AV:lamlamA+lamlamA:lam,lam', 'AV-AV', 'bubble', ('A', AVc), ('A', AVc), {('lam', 'lam')}),
    ('AV-AV:AAAA:AV', 'AV-AV', 'tadpole', ('A', AVc), ('A', AVc), {('AV', 'AV')}),
    ('AV-AV:AAAA:Ap', 'AV-AV', 'tadpole', ('A', AVc), ('A', AVc), {('Ap', 'Ap')}),
    ('Ap-Ap:AAA+AAA:Ap,AV', 'Ap-Ap', 'bubble', ('A', APc), ('A', APc), {('AV', 'AV'), ('Ap', 'Ap')}),
    ('Ap-Ap:AcAs+AcAs:c-As,c-As', 'Ap-Ap', 'bubble', ('A', APc), ('A', APc), {('c', 'As'), ('As', 'c')}),
    ('Ap-Ap:lamlamA+lamlamA:lam,lam', 'Ap-Ap', 'bubble', ('A', APc), ('A', APc), {('lam', 'lam')}),
    ('Ap-Ap:AAAA:AV', 'Ap-Ap', 'tadpole', ('A', APc), ('A', APc), {('AV', 'AV')}),
    ('Ap-Ap:AAAA:Ap', 'Ap-Ap', 'tadpole', ('A', APc), ('A', APc), {('Ap', 'Ap')}),
    ('AV-Ap:AAA+AAA:all', 'AV-Ap', 'bubble', ('A', AVc), ('A', APc), None),
    ('AV-Ap:AAAA:all', 'AV-Ap', 'tadpole', ('A', AVc), ('A', APc), None),
    ('lam-lam:lamlamA+lamlamA:lam,AV', 'lam-lam', 'bubble', ('lam', None), ('lam', None), {('lam', 'lam'), ('AV', 'AV')}),
    ('lam-lam:lamlamA+lamlamA:lam,Ap', 'lam-lam', 'bubble', ('lam', None), ('lam', None), {('lam', 'lam'), ('Ap', 'Ap')}),
    ('lam-lam:all', 'lam-lam', 'bubble', ('lam', None), ('lam', None), None),
    ('As-c:AcAs+AcAs:AV,c-As', 'As-c', 'bubble', ('As', AVc), ('c', None), {('AV', 'AV'), ('c', 'As'), ('As', 'c')}),
    ('As-c:all', 'As-c', 'bubble', ('As', None), ('c', None), None),
]

def run(selected=None):
    out = {}
    for (gid, sector, kind, s1, s2, mask) in GRAPHS:
        if selected and gid not in selected:
            continue
        t0 = time.time()
        e1 = leg(s1[0], 'u', p, s1[1]); e2 = leg(s2[0], 'v', mp, s2[1])
        a, b = (bubble if kind == 'bubble' else tadpole)(e1, e2, mask)
        names = {'AV-AV': ['B1', 'B2', 'B3'], 'Ap-Ap': ['B4', 'B5'], 'AV-Ap': ['B3'], 'lam-lam': ['B6'],
                 'As-c': ['B7']}[sector]
        fa = fit(a, e1, e2, names); fb = fit(b, e1, e2, names)
        out[gid] = dict(sector=sector, kind=kind, eps_inv=({k: str(v) for k, v in fa.items()} if fa is not None else 'NONLOCAL:' + str(a)),
                        log_eps=({k: str(v) for k, v in fb.items()} if fb is not None else 'NONLOCAL:' + str(b)),
                        raw_a=str(a)[:400], raw_b=str(b)[:400], seconds=round(time.time() - t0, 1))
        print(gid, out[gid]['eps_inv'], out[gid]['log_eps'], out[gid]['seconds'], flush=True)
    return out

if __name__ == '__main__':
    sel = sys.argv[1:] or None
    res = run(sel)
    tag = 'all' if not sel else 'part'
    with open(f'two_point_{tag}.json', 'w') as fh:
        json.dump(res, fh, indent=1)
