"""Engine B, dilatation route (ii): direct UV-boundary graphs at finite eps and finite external momenta.
For a graph X(eps) = int prod_e dt_e f(t), the action of the dilatation Delta = -(E+w) on the regulated graph is
   Delta.X = sum_e X[edge e -> Delta.P(eps,infty)],   Delta.P(eps,infty) = -2 eps d/deps P(eps,infty) = 2 eps L_eps,
i.e. D(eps) = 2 eps^E sum_e int_{[1,inf)^{E-1}} f(t_e = eps, t_f = eps u_f) du.  Here D(eps) is evaluated numerically with the
exact Gaussian factor e^{-R} (no Taylor expansion), for several eps, and the eps -> 0 behaviour D = 2a/eps - 2b + o(1) is
extracted by a fit; a, b are then compared with the counterterm route.  Rational random momenta/polarisations."""
import json, math, random
import numpy as np
import sympy as sp
from scipy import integrate
import vertices as V
import loops as L
import sectors as S
import o_graphs as OG

random.seed(11)
def rnd(n):
    return [sp.Rational(random.randint(-9, 9), random.randint(1, 5)) for _ in range(n)]

def f_builder(Num, lines, subs):
    E = len(lines)
    betas = sp.symbols(f'beta0:{E}')
    shift = [sum(betas[e] * lines[e][m] for e in range(E)) for m in range(4)]
    C = L.gaussian_moments(Num, shift)
    P2 = lambda v: sum(x * x for x in v)
    RoverT = sp.expand(sum(betas[e] * P2(lines[e]) for e in range(E)) - P2(shift))
    T = sp.Symbol('T', positive=True)
    expr = sum(Ck * (2 * T) ** (-k) for k, Ck in C.items()) * T**-2 / (16 * sp.pi**2)
    expr = expr.subs(subs); RT = RoverT.subs(subs)
    fn = sp.lambdify(list(betas) + [T], sp.expand(expr), 'numpy')
    rn = sp.lambdify(list(betas), RT, 'numpy')
    def f(ts):
        Tt = sum(ts)
        b = [t / Tt for t in ts]
        return complex(fn(*b, Tt)) * math.exp(-Tt * float(rn(*b)))
    return f, E

def D_of_eps(f, E, eps):
    tot = 0
    for e in range(E):
        if E == 1:
            tot += 2 * eps * f([eps])
        elif E == 2:
            g = lambda s: f([eps if i == e else eps / s for i in range(2)]) / s**2
            re = integrate.quad(lambda s: g(s).real, 0, 1, limit=200, epsabs=1e-13, epsrel=1e-11)[0]
            im = integrate.quad(lambda s: g(s).imag, 0, 1, limit=200, epsabs=1e-13, epsrel=1e-11)[0]
            tot += 2 * eps**2 * (re + 1j * im)
        else:
            others = [i for i in range(3) if i != e]
            def g(s1, s2):
                ts = [0, 0, 0]; ts[e] = eps; ts[others[0]] = eps / s1; ts[others[1]] = eps / s2
                return f(ts) / (s1**2 * s2**2)
            re = integrate.dblquad(lambda s2, s1: g(s1, s2).real, 0, 1, 0, 1, epsabs=1e-12, epsrel=1e-9)[0]
            im = integrate.dblquad(lambda s2, s1: g(s1, s2).imag, 0, 1, 0, 1, epsabs=1e-12, epsrel=1e-9)[0]
            tot += 2 * eps**3 * (re + 1j * im)
    return tot

def fit_limit(epss, vals, with_pole):
    # model: D = A/eps + B + eps*(C log eps + Dc)
    rows = []
    for e in epss:
        r = ([1 / e] if with_pole else []) + [1, e * math.log(e), e]
        rows.append(r)
    M = np.array(rows); y = np.array(vals)
    sol = np.linalg.lstsq(M, y, rcond=None)[0]
    return (sol[0], sol[1]) if with_pole else (0.0, sol[0])

def self_energy_check():
    """A_perp two-point (all graphs): D(eps)*16pi^2 -> 2a/eps - 2b, with a = 9/4 B5, b = -1 B4 (units)."""
    k = rnd(4); u = [0] * 4 + rnd(6); v = [0] * 4 + rnd(6)
    e1 = S.leg('A', 'u', S.p, S.APc); e2 = S.leg('A', 'v', S.mp, S.APc)
    subs = {**dict(zip(S.p, k)), **{sp.Symbol(f'u{i}'): u[i] for i in range(4, 10)}, **{sp.Symbol(f'v{i}'): v[i] for i in range(4, 10)}}
    graphs = []
    for (x, y, px, sgn) in [(e1, e2, S.p, 1), (e2, e1, S.mp, 1)]:
        h1 = L.full_hessian(V.CUBIC, [x], [-w for w in L.q], [L.q[m] - px[m] for m in range(4)], S.CUBW)
        h2 = L.full_hessian(V.CUBIC, [y], [px[m] - L.q[m] for m in range(4)], list(L.q), S.CUBW)
        N1 = L.prop_numerator(list(L.q)); N2 = L.prop_numerator([L.q[m] - px[m] for m in range(4)])
        Num = L.str_product([N1, L.PiH(h1), N2, L.PiH(h2)])
        graphs.append((sp.Rational(1, 4) * sgn, Num, [[0] * 4, list(px)]))
    h = L.full_hessian(V.QUARTIC, [e1, e2], [-w for w in L.q], list(L.q), S.QUAW)
    Num = L.str_product([L.prop_numerator(list(L.q)), L.PiH(h)])
    graphs.append((sp.Rational(1, 2), Num, [[0] * 4]))
    epss = [2e-3, 1e-3, 5e-4, 2.5e-4, 1.25e-4]
    vals = []
    fs = [(wt, f_builder(Num, lines, subs)) for (wt, Num, lines) in graphs]
    for eps in epss:
        vals.append(sum(float(wt) * D_of_eps(f, E, eps) for wt, (f, E) in fs) * 16 * math.pi**2)
    A, B = fit_limit(epss, [x.real for x in vals], True)
    k2 = sum(x * x for x in k)
    uv = sum(u[i] * v[i] for i in range(4, 10))
    # predicted: a = (9/4) * B5-value, b = (-1) * B4-value ; B5-value = 2 u.v, B4-value = 2 k^2 u.v
    pa = float(sp.Rational(9, 4) * 2 * uv); pb = float(-1 * 2 * k2 * uv)
    return dict(fit_pole=A, pred_pole=2 * pa, fit_finite=B, pred_finite=-2 * pb, max_imag=max(abs(x.imag) for x in vals))

def o_vertex_check():
    """1PI O-graphs, A_perp legs: D(eps)*16pi^2 -> -2 b with b = 4 O - 2 (w1,w2) sum_w (iota_w A, iota_w A)."""
    p1v = [x / 10 for x in rnd(4)]; p2v = [x / 10 for x in rnd(4)]; u = [0] * 4 + rnd(6); v = [0] * 4 + rnd(6); a1 = rnd(6); a2 = rnd(6)
    subs = {**dict(zip(OG.p1, p1v)), **dict(zip(OG.p2, p2v)),
            **{sp.Symbol(f'u{i}'): u[i] for i in range(4, 10)}, **{sp.Symbol(f'v{i}'): v[i] for i in range(4, 10)},
            **{sp.Symbol(f'w1_{i}'): a1[i] for i in range(6)}, **{sp.Symbol(f'w2_{i}'): a2[i] for i in range(6)}}
    e1 = S.leg('A', 'u', OG.p1, S.APc); e2 = S.leg('A', 'v', OG.p2, S.APc)
    words = [(sp.Rational(1, 2), [('cubic', [e1], list(OG.p1)), ('cubic', [e2], list(OG.p2)), ('O', [], OG.qO)], {('AV', 'AV'), ('Ap', 'Ap')}),
             (sp.Rational(1, 2), [('cubic', [e2], list(OG.p2)), ('cubic', [e1], list(OG.p1)), ('O', [], OG.qO)], {('AV', 'AV'), ('Ap', 'Ap')}),
             (sp.Rational(1, 4), [('quartic', [e1, e2], [OG.p1[m] + OG.p2[m] for m in range(4)]), ('O', [], OG.qO)], {('Ap', 'Ap')}),
             (sp.Rational(1, 4), [('O', [], OG.qO), ('quartic', [e1, e2], [OG.p1[m] + OG.p2[m] for m in range(4)])], {('Ap', 'Ap')})]
    fs = []
    for wt, word, mask in words:
        mats = []; lines = []; cur = [sp.Integer(0)] * 4
        for (kind, ext, inj) in word:
            lines.append(list(cur))
            Nq = L.prop_numerator([L.q[m] - cur[m] for m in range(4)], mask)
            nxt = [cur[m] + inj[m] for m in range(4)]
            hh = OG.H_of(kind, ext, [-(L.q[m] - cur[m]) for m in range(4)], [L.q[m] - nxt[m] for m in range(4)])
            mats += [Nq, L.PiH(hh)]; cur = nxt
        Num = L.str_product(mats)
        fs.append((wt, f_builder(Num, lines, subs)))
    epss = [2e-2, 1e-2, 5e-3, 2.5e-3, 1.25e-3, 6.25e-4]
    vals = [sum(float(wt) * D_of_eps(f, E, eps) for wt, (f, E) in fs) * 16 * math.pi**2 for eps in epss]
    A, B = fit_limit(epss, [x.real for x in vals], False)
    Ov = sum(u[i] * a1[i - 4] for i in range(4, 10)) * sum(v[i] * a2[i - 4] for i in range(4, 10)) + \
         sum(u[i] * a2[i - 4] for i in range(4, 10)) * sum(v[i] * a1[i - 4] for i in range(4, 10))
    tr = sum(a1[i] * a2[i] for i in range(6)) * 2 * sum(u[i] * v[i] for i in range(4, 10))
    pb = float(4 * Ov - 2 * tr)
    return dict(values=[x.real for x in vals], eps=epss, fit_finite=B, pred_finite=-2 * pb, max_imag=max(abs(x.imag) for x in vals))

if __name__ == '__main__':
    import sys
    if 'oonly' in sys.argv:
        r2 = o_vertex_check(); print('O-vertex boundary check', r2, flush=True)
        prev = json.load(open('route2_numeric.json'))
        prev['O_vertex_Aperp'] = r2
        json.dump(prev, open('route2_numeric.json', 'w'), indent=1, default=float); sys.exit()
    r1 = self_energy_check(); print('self-energy boundary check', r1, flush=True)
    r2 = o_vertex_check(); print('O-vertex boundary check', r2, flush=True)
    with open('route2_numeric.json', 'w') as fh:
        json.dump({'self_energy_Aperp': r1, 'O_vertex_Aperp': r2}, fh, indent=1, default=float)
