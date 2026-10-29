"""Engine A, part C: gates (value, target, pass) recorded in GATES."""
import itertools, math, time, json
from fractions import Fraction
import numpy as np
import sympy as sp
import mpmath as mp
import cq_core as Q, cq_gcore as G, cq_gval as GV, cq_eval as V, cq_fsl2 as F2
import ea_core as E, ea_ops as EO

GATES = []
ONE = {(0, 0, 0, 0): 1}


def record(gid, what, value, target, ok, kind='gate'):
    GATES.append({'id': gid, 'what': what, 'value': value, 'target': target, 'pass': bool(ok), 'kind': kind})
    print('[%s] %s: %s' % (gid, 'PASS' if ok else 'FAIL', what), flush=True)
    return ok


def u6(i):
    v = np.zeros(6, complex); v[i] = 1; return v


def rv(n, rng):
    return np.array([complex(int(x), int(y)) for x, y in zip(rng.integers(-2, 3, n), rng.integers(-2, 3, n))])


def g_weights():
    """finite-part simplex weights vs mpmath quadrature (E = 3), Lemma tl-special example"""
    mp.mp.dps = 12
    rows = []; ok = True
    for exps, s in [((0, 0, 0), 0), ((0, 2, 1), 0), ((2, 0, 1), 2)]:
        r0, r1 = Q.fin_monomial(exps, s)
        ex = mp.mpf(r0.numerator) / r0.denominator + mp.mpf(r1.numerator) / r1.denominator * mp.log(mp.mpf(3) / 2)
        others = [g for g in range(3) if g != s]
        def f(a1, a2):
            a = [1 - a1 - a2, a1, a2]
            if min(a) < 0:
                return 0
            m = min(a[g] for g in others)
            if a[s] > m:
                return 0
            return a[0] ** exps[0] * a[1] ** exps[1] * a[2] ** exps[2] * mp.log(m / a[s])
        nm = mp.quad(lambda x: mp.quad(lambda y: f(x, y), [0, (1 - x) / 3, (1 - x) / 2, 1 - x]), [0, mp.mpf(1) / 3, mp.mpf(1) / 2, 1])
        good = abs(ex - nm) < 1e-6
        ok &= good
        rows.append({'exps': exps, 's': s, 'exact': '%s + %s log(3/2)' % (r0, r1), 'numeric': float(nm), 'diff': float(abs(ex - nm))})
    record('GW1', 'finite-part weights int_{R_s} alpha^n log(min_{f!=s} alpha_f/alpha_s) exact vs mpmath (E=3)', rows, '|diff| < 1e-6', ok)
    p = np.array([1, 2, -1, 1])
    S = Q.CSimplex([np.zeros(4, int), -p])
    N = {}
    for mu in range(4):
        e = [0, 0, 0, 0]; e[mu] = 2; N[tuple(e)] = complex(-1)
    f0, f1 = S.fin(N, 0); sp_ = S.special(N, 0)
    record('GW2', 'Lemma tl-special example: k=1, N=-|l|^2: eps^-1 and finite part', {'eps^-1': str(sp_[1][0]), 'finite': str(f0[0]), 'log part': str(f1[0])},
           {'eps^-1': '-1/2', 'finite': str(Fraction(-3, 8) * int(p @ p))}, sp_[1][0] == Fraction(-1, 2) and f0[0] == Fraction(-3, 8) * int(p @ p) and f1[0] == 0)


def g_k1():
    """Prop prop:cc-k1: k = 1 closed form (sympy) vs the simplex engine, scalar numerators"""
    b = sp.Symbol('beta', positive=True)
    W0 = (1 - b) * sp.log((1 - b) / (2 * b)) - (1 - 3 * b)
    W1 = 2 * ((((1 - b) / 2) ** 2 - b ** 2) / 2 - b * ((1 - b) / 2 - b))
    p = np.array([1, -1, 2, 0]); a = np.array([0, 1, 1, -1])
    rows = []; ok = True
    ks = sp.symbols('k0:4')
    for name, Nsym in [('-|k|^2', -sum(k ** 2 for k in ks)),
                       ('-|k|^2 (a.k)^2', -sum(k ** 2 for k in ks) * sum(int(a[i]) * ks[i] for i in range(4)) ** 2),
                       ('-|k|^4', -sum(k ** 2 for k in ks) ** 2)]:
        # engine: E = 3, special index 0 (q = 0), segments 1, 2 (q = -p)
        poly = sp.Poly(sp.expand(Nsym), *ks)
        N = {tuple(m): complex(int(c)) for m, c in zip(poly.monoms(), poly.coeffs())}
        S = Q.CSimplex([np.zeros(4, int), -p, -p])
        f0, f1 = S.fin(N, 0)
        eng = (sp.Rational(f0[0]) if f0 else 0) + (sp.Rational(f1[0]) if f1 else 0) * sp.log(sp.Rational(3, 2))
        pol = S.special(N, 0)
        eng1 = sp.Rational(pol[1][0]) if 1 in pol and pol[1] is not None else 0
        # closed form: centre c = (1-beta) p, Phi = beta (1-beta) |p|^2, F_m = sum_{j>=m+1} (Lap^j N)(c)/(j! 4^j) (-Phi)^{j-1-m}/(j-1-m)!
        c = [(1 - b) * int(p[i]) for i in range(4)]
        Phi = b * (1 - b) * int(p @ p)
        def lap(f):
            return sum(sp.diff(f, k, 2) for k in ks)
        def Fm(m):
            tot = 0; L = sp.expand(Nsym); j = 0
            while L != 0:
                if j >= m + 1:
                    tot += L.subs(dict(zip(ks, c))) / (sp.factorial(j) * 4 ** j) * (-Phi) ** (j - 1 - m) / sp.factorial(j - 1 - m)
                L = sp.expand(lap(L)); j += 1
            return sp.expand(tot)
        cf = sp.nsimplify(sp.integrate(sp.expand(Fm(0) * W0), (b, 0, sp.Rational(1, 3))))
        cf = sp.nsimplify(sp.expand(sp.expand_log(cf, force=True)))
        c1 = sp.nsimplify(sp.integrate(sp.expand(Fm(1) * W1), (b, 0, sp.Rational(1, 3))))
        good = sp.simplify(sp.expand_log(cf - eng, force=True)) == 0 and sp.simplify(c1 - eng1) == 0
        ok &= bool(good)
        rows.append({'N': name, 'engine_finite': str(eng), 'closed_form_finite': str(cf), 'engine_eps^-1': str(eng1), 'closed_form_eps^-1': str(c1)})
    record('GK1', 'Prop cc-k1 (k=1 closed form W_0, W_1 over [0,1/3]) vs the ordered-subsimplex engine', rows, 'equal', ok)


def g_cores(rng):
    """general core (Gaussian momenta, t-derivative) vs int core; D0 by t-derivative vs finite differences;
    pure two-letter log kernel vs ea_kern"""
    p1 = np.array([1, 0, 2, -1]); p2 = np.array([0, 1, -1, 1])
    n = 0; bad = 0
    for ta, tb, outs in [('AV', 'AV', ['AV', 'AV']), ('AW', 'lamP', ['AW', 'lamP']), ('AV', 'lamP', ['lamP', 'AV']), ('AW', 'AW', ['AV', 'AW', 'AW'])]:
        RN = {'AV': 4, 'AW': 6, 'lamP': 16}
        za, zb = rv(RN[ta], rng), rv(RN[tb], rng)
        legs = [(t, rv(RN[t], rng), m) for t, m in zip(outs, [p1, p2, p1 - 2 * p2])]
        fa = {(1, 0, 0, 0): 1j, (0, 1, 0, 0): -1}; fb = {(0, 0, 1, 0): 1j}
        for st in V.structures(legs):
            ob = V.blocks_of(st, legs); nb = len(ob)
            for v2pos in range(nb + 1):
                for spos in range(nb + 2):
                    a = Q.config_value(ob, v2pos, spos, 'P', ta, tb, fa=fa, fb=fb, v2dir=2, want=('fin', 'poles'))
                    b = G.gconfig(ob, v2pos, spos, 'P', ta, tb, fa=fa, fb=fb, v2dir=2, kinds=('fin', ('spec', 1)))
                    if a is None and b is None:
                        continue
                    va = V.from_fin(a['fin'], za, zb) if a else V.vzero()
                    vb = V.from_fin((b['fin'][0][0], b['fin'][0][1], b['E']), za, zb) if b else V.vzero()
                    n += 1
                    if V.vstr(va) != V.vstr(vb):
                        bad += 1
    record('GE1', 'general core (Gaussian momenta, t-expansion) = integer core on every configuration (finite parts)', '%d configs, %d mismatches' % (n, bad), '0 mismatches', bad == 0)
    Xa = ('AV', rv(4, rng), {(0, 1, 0, 0): 1j}); Xb = ('AW', rv(6, rng), ONE)
    legs = [('AV', rv(4, rng), p1), ('AW', rv(6, rng), p2)]
    Q.RULE['koszul'] = 'R2'
    d_fd = V.pair_on_legs(Xa, Xb, legs, 'D0')
    GV.INTR['on'] = False
    d_t = GV.combine(GV.kernel(Xa, Xb, legs, 'D0'), GV.VVSIGN)
    GV.INTR['on'] = True
    record('GE2', 'translation-invariant part along Delta: exact t^1 coefficient = finite differences in r (int core)', V.vstr(d_t), V.vstr(d_fd), V.vstr(d_t) == V.vstr(d_fd))
    rows = []; ok = True
    for ta, tb, outs in [('AW', 'AW', ['AW', 'AW']), ('AW', 'lamP', ['AW', 'lamP']), ('AV', 'lamP', ['lamP']), ('lamP', 'lamP', ['lamP', 'lamP']),
                         ('AV', 'c', ['c']), ('AV', 'AW', ['AV', 'AW'])]:
        RN = {'AV': 4, 'AW': 6, 'lamP': 16, 'c': 1}
        za, zb = rv(RN[ta], rng), rv(RN[tb], rng)
        legs = [(t, rv(RN[t], rng), m) for t, m in zip(outs, [p1, p2])]
        a = EO.pair_kernel_on_legs((ta, za, {(1, 0, 0, 0): 1j}), (tb, zb, ONE), legs)
        b = GV.logkernel((ta, za, {(1, 0, 0, 0): 1j}), (tb, zb, ONE), legs)
        good = E.qstr(a) == V.vstr(b) or (a == (0, 0) and V.viszero(b))
        ok &= good
        rows.append([ta, tb, outs, E.qstr(a), V.vstr(b)])
    record('GE3', 'pure two-letter log kernel of the general core = ea_kern.kernel (two-letter engine A)', rows, 'equal', ok)


def g_bilinear():
    """GC3: C^(2)(A_perp, A_perp)|_{A_perp A_perp} on the bilinears, traceless and trace, real and null momenta"""
    out = {}; ok = True
    for lab, (q1, q2) in [('real', (np.array([1, 0, 2, -1]), np.array([0, 1, -1, 1]))), ('null', (np.array([1, 1j, 0, 0]), np.array([0, 0, 1, -1j])))]:
        for part in ('traceless', 'trace'):
            pairs = [(0, 1)] if part == 'traceless' else [(I, I) for I in range(6)]
            leg_pol = (0, 1) if part == 'traceless' else (0, 0)
            tot = {'P': V.vzero(), 'D0': V.vzero()}
            for (a, b) in pairs:
                Xa = ('AW', u6(a), ONE); Xb = ('AW', u6(b), ONE)
                for lg in ([('AW', u6(leg_pol[0]), q1), ('AW', u6(leg_pol[1]), q2)], [('AW', u6(leg_pol[1]), q2), ('AW', u6(leg_pol[0]), q1)]):
                    tot['P'] = V.vadd(tot['P'], GV.combine(GV.kernel(Xa, Xb, lg, ('P', 1)), GV.VVSIGN))
                    tot['D0'] = V.vadd(tot['D0'], GV.combine(GV.kernel(Xa, Xb, lg, 'D0'), GV.VVSIGN))
            Oleg = sum(int(a == leg_pol[0]) * int(b == leg_pol[1]) + int(a == leg_pol[1]) * int(b == leg_pol[0]) for a, b in pairs)
            tP = sp.Rational(-1, 3) * F2.symc(1j * (q1 + q2)[1]) * Oleg
            vP = F2.symv(tot['P']); v0 = F2.symv(tot['D0'])
            good = sp.simplify(vP - tP) == 0 and sp.simplify(v0 + sp.Rational(8, 3) * Oleg) == 0
            ok &= good
            out['%s_%s' % (lab, part)] = {'C_P1 (units u)': str(vP), 'target -1/3 u d_1 O': str(tP), 'C0_Delta': str(v0), 'target -8/3 u O': str(-sp.Rational(8, 3) * Oleg)}
    record('GC3', 'bilinears: C^(2)_{P_b} = -1/3 u d_b O and C^(2),0_Delta = -8/3 u O (= -1/3 u v.^0 - 4/3 u), traceless and trace, real and null momenta, no free constant',
           out, 'Example ex:cc-bilinear', ok)


def g_mirror(rng):
    """GS: graded symmetry C(X_b,X_a)[reversed] / C(X_a,X_b) = the same ratio as the two-letter log kernel, every kappa"""
    p1 = np.array([1, 0, 2, -1]); p2 = np.array([0, 1, -1, 1])
    RN = {'AV': 4, 'AW': 6, 'lamP': 16, 'lamM': 16, 'c': 1}
    rows = []; ok = True
    rat = lambda x, y: None if (y[0] == 0 and y[1] == 0) else complex(complex(float(x[0]), float(x[1])) / complex(float(y[0]), float(y[1])))
    for ta, tb, outs in [('AV', 'lamP', ['lamP']), ('AV', 'AW', ['AW']), ('AW', 'lamP', ['AW', 'lamP']), ('lamP', 'lamP', ['lamP', 'lamP']),
                         ('AV', 'lamP', ['AV', 'lamP']), ('AV', 'lamP', ['lamP', 'AV']), ('AV', 'AW', ['AV', 'AW']), ('AV', 'AV', ['AV', 'AV']),
                         ('AW', 'AW', ['AW', 'AW']), ('lamP', 'lamM', ['lamM', 'lamP'])]:
        za, zb = rv(RN[ta], rng), rv(RN[tb], rng)
        legs = [(t, rv(RN[t], rng), m) for t, m in zip(outs, [p1, p2])]
        Xa = (ta, za, ONE); Xb = (tb, zb, ONE)
        l1 = EO.pair_kernel_on_legs(Xa, Xb, legs); l2 = EO.pair_kernel_on_legs(Xb, Xa, legs[::-1])
        lr = rat(l1, l2)
        f = GV.kernel(Xa, Xb, legs, ('P', 1)); g = GV.kernel(Xb, Xa, legs[::-1], ('P', 1))
        for kap in (1, -1, 0):
            fk = GV.combine(f, kap); gk = GV.combine(g, kap)
            rs = [rat(a, b) for a, b in zip(fk, gk) if not (b[0] == 0 and b[1] == 0)]
            good = all(abs(r - lr) < 1e-9 for r in rs) if lr is not None else True
            ok &= good
            rows.append({'pair': [ta, tb], 'outputs': outs, 'kappa': kap, 'composite_ratios': [str(complex(round(r.real, 9), round(r.imag, 9))) for r in rs],
                         'log_ratio': None if lr is None else str(complex(round(lr.real, 9), round(lr.imag, 9))), 'pass': bool(good)})
    record('GS', 'graded symmetry of the composite kernels (chain from X_b with reversed outputs) has the log-kernel ratio, kappa = +1, -1, 0', rows, 'equal ratios', ok)


def g_position(rng):
    """GC4: position part (s_mu vertex at r=0) = 2 x P_mu graph; GC4': [C^0, d_b] = C^b = 2 C_{P_b} on letters with derivatives"""
    p1 = np.array([1, 0, 2, -1]); p2 = np.array([0, 1, -1, 1])
    rows = []; ok = True
    rows2 = []; ok2 = True
    for ta, tb, outs, fa in [('AV', 'AW', ['AW'], {(0, 0, 1, 0): 1j}), ('AW', 'AW', ['AW', 'AW'], {(1, 0, 0, 0): 1j}),
                             ('AV', 'lamP', ['lamP'], {(0, 1, 0, 0): 1j}), ('AW', 'lamP', ['AW', 'lamP'], {(0, 0, 0, 1): 1j})]:
        RN = {'AV': 4, 'AW': 6, 'lamP': 16}
        za, zb = rv(RN[ta], rng), rv(RN[tb], rng)
        legs = [(t, rv(RN[t], rng), m) for t, m in zip(outs, [p1, p2])]
        Xa = (ta, za, fa); Xb = (tb, zb, ONE)
        for b in (0, 3):
            s_ = GV.combine(GV.kernel(Xa, Xb, legs, ('Ds', b)), GV.VVSIGN)
            P_ = GV.combine(GV.kernel(Xa, Xb, legs, ('P', b)), GV.VVSIGN)
            good = V.vstr(s_) == V.vstr(V.vscale(P_, (Fraction(2), Fraction(0))))
            ok &= good
            rows.append({'pair': [ta, tb], 'b': b, 'pos_part': V.vstr(s_), '2C_Pb': V.vstr(V.vscale(P_, (Fraction(2), Fraction(0))))})
            # [M^0, d_b] = M^b
            def dlet(X):
                t, z, f = X
                out = {}
                for e, c in f.items():
                    e2 = list(e); e2[b] += 1
                    out[tuple(e2)] = out.get(tuple(e2), 0) + 1j * c
                return (t, z, out)
            c_a = GV.combine(GV.kernel(dlet(Xa), Xb, legs, 'D0'), GV.VVSIGN)
            c_b = GV.combine(GV.kernel(Xa, dlet(Xb), legs, 'D0'), GV.VVSIGN)
            c_0 = GV.combine(GV.kernel(Xa, Xb, legs, 'D0'), GV.VVSIGN)
            qtot = sum(np.asarray(l[2]) for l in legs)
            lhs = V.vadd(V.vadd(c_a, c_b), V.vscale(c_0, F2.CO.gauss(-1j * qtot[b]) if hasattr(F2, 'CO') else (Fraction(0), Fraction(-int(qtot[b])))))
            rhs = V.vscale(P_, (Fraction(2), Fraction(0)))
            good2 = V.vstr(lhs) == V.vstr(rhs)
            ok2 &= good2
            rows2.append({'pair': [ta, tb], 'b': b, '[C0,d_b]': V.vstr(lhs), '2 C_Pb': V.vstr(rhs), 'pass': good2})
    record('GC4', 'position identity (C^(2)_Delta)^b = 2 C^(2)_{P_b} graph by graph (s_mu vertex at r = 0)', rows, 'equal', ok)
    record('GC4c', 'canonical form (Lemma lem:cc-canonical): C^0(d_b X_a, X_b) + C^0(X_a, d_b X_b) - d_b C^0(X_a, X_b) = (C_Delta)^b = 2 C_{P_b}', rows2, 'equal', ok2)
