"""Step 2, Part D4, engine A, task (2): the one-particle-irreducible graph class Y1a of the arity-one part of O[infty]
along Delta, legs (A|_Wperp, A^*|_Wperp): one loop through the vertices
   S_M (composite A|_V--c line: S_M's A leg joined by the A--A component, its A^* leg by the A^*|_V--c component of P),
   O = int omega (iota_w1 A, iota_w2 A)  (omega = e^{i(k,x)} dvol),
   X = cubic vertex with the external A|_Wperp leg (momentum p1, held fixed),
   Y = cubic vertex ([A,c],A^*) with the external A^*|_Wperp leg (momentum p2 = -k - p1 - r, the dependent leg).
All edges P(eps,infty).  The class is UV convergent (superficial degree -2 for alpha, -1 for beta); its eps -> 0 limit
is the eps = 0 integral (only Schwinger exponents a >= 0 occur; asserted).  The Q-image (A^* leg -> -(1/2) d^*d A,
s = +1 of d4_signs) in canonical form (fixed p1, k):
   alpha_Q = -(1/2)|p2|^2 alpha + p2 . beta,     beta_Q = -(1/2)|p2|^2 beta.
Numerics: per oriented sequence the exact polynomial integrand (d4_cycles / d4_integrals) is homogenised, the split of
the composite pair (edges 1, 2: equal momenta at r = 0) is integrated exactly, and the remaining triangle integral
int_{simplex} poly(U, u3, u4) q^{-a-1} is done with mpmath (tanh-sinh, 2-dim) at mp.dps digits."""
import os, sys, json, time, math
from fractions import Fraction
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import s2_core as C
from r2_core import Q2
import d4_cycles as DC
import d4_integrals as DI
import d4_two_point as TP
from s2_graphs import Poly


def reduce_composite(Ph):
    """homogeneous monomials in (u1,u2,u3,u4) -> dict over (E, e3, e4) with U = u1 + u2:
    int_0^U u1^e1 (U-u1)^e2 du1 = e1! e2!/(e1+e2+1)! U^{e1+e2+1}"""
    out = {}
    for (e1, e2, e3, e4), c in Ph.items():
        f = Fraction(math.factorial(e1) * math.factorial(e2), math.factorial(e1 + e2 + 1))
        key = (e1 + e2 + 1, e3, e4)
        out[key] = out[key] + c * f if key in out else c * f
    return out


def tri_integral(polys_by_pow, q3, q4):
    """int over the triangle U + u3 + u4 = 1 (measure dU du3) of sum_c sum c U^E u3^e3 u4^e4 * q^{-c}  (all powers c in
    ONE integrand, so that cancellations between Schwinger exponents near the vertices are kept),
    q = u3 |q3|^2 + u4 |q4|^2 - |u3 q3 + u4 q4|^2"""
    a33 = mp.mpf(int(q3 @ q3)); a44 = mp.mpf(int(q4 @ q4)); a34 = mp.mpf(int(q3 @ q4))
    items = []
    for cpow, polyd in polys_by_pow.items():
        for (E, e3, e4), c in polyd.items():
            items.append((cpow, E, e3, e4, mp.mpf(int(c.re.numerator)) / int(c.re.denominator),
                          mp.mpf(int(c.im.numerator)) / int(c.im.denominator)))
    def f(x, y, part):
        U = x; u3 = (1 - x) * y; u4 = (1 - x) * (1 - y)
        q = u3 * a33 + u4 * a44 - (u3 * u3 * a33 + 2 * u3 * u4 * a34 + u4 * u4 * a44)
        s = mp.mpf(0)
        for (cpow, E, e3, e4, cr, ci) in items:
            s += (cr if part == 0 else ci) * U ** E * u3 ** e3 * u4 ** e4 / q ** cpow
        return s * (1 - x)
    re = mp.quad(lambda x, y: f(x, y, 0), [0, 1], [0, 1])
    im = mp.quad(lambda x, y: f(x, y, 1), [0, 1], [0, 1])
    return mp.mpc(re, im)


def conv_value(term_lists_with_q, dps):
    """sum over (terms, qs, factor) of factor * sum_terms w Gamma(b+1) int P q^{-b-1}; returns mpc"""
    tot = mp.mpc(0)
    for (terms, qs, factor) in term_lists_with_q:
        q3 = np.array([int(x) for x in qs[2]]); q4 = np.array([int(x) for x in qs[3]])
        assert all(int(x) == 0 for x in qs[1]), 'composite edges must carry equal momenta at r = 0'
        bypow = {}
        for (P, b, w) in terms:
            assert b >= 0, ('non-convergent Schwinger exponent in a convergent class', b)
            Ph = DI.homogenize(P, 4)
            red = reduce_composite(Ph)
            d = bypow.setdefault(b + 1, {})
            for k_, v in red.items():
                vv = v * Q2(Fraction(w) * math.factorial(b))
                d[k_] = d[k_] + vv if k_ in d else vv
        bypow = {c: {k_: v for k_, v in d.items() if not v.iszero()} for c, d in bypow.items()}
        bypow = {c: d for c, d in bypow.items() if d}
        if bypow:
            tot += factor * tri_integral(bypow, q3, q4)
    return tot


def y1a(p1, k, a, bpol, w1, w2, degree=4, dps=20, verbose=True):
    """canonical alpha, beta^mu of Y1a (before Q), and of its Q-image, at eps -> 0, in units; also the O-structure"""
    mp.mp.dps = dps
    p1 = np.asarray(p1, dtype=int); k = np.asarray(k, dtype=int)
    X = DC.Bg([('A', a, p1)])
    O = DC.QF('O', k, w1=w1, w2=w2)
    def ins(qf, rmom):
        Y = DC.Bg([('As', bpol, -k - p1 - rmom)])
        return [qf, O, X, Y]
    canon = [0, 3]            # (xi, r of the A^* leg at Y)
    scale = 4
    def collect(insertions, deriv_build=None):
        out = []
        if deriv_build is None:
            cyc = DC.Cycle(insertions, degree)
            for order, Ls in cyc.sequences():
                N, qs = cyc.numerator_fn(order, Ls)
                poly = DC.poly_from(N, degree)
                if not poly.t:
                    continue
                sgn = cyc.parity_sign(order, canon) * (-1)
                out.append((DI.terms(poly, qs, degree), qs, mp.mpf(sgn) / (2 * scale)))
            return out
        cycs = {t: DC.Cycle(deriv_build(t), degree) for t in TP.TVALS}
        c0 = cycs[0]
        for order, Ls in c0.sequences():
            polys = {}; qss = {}
            for t in TP.TVALS:
                N, qs = cycs[t].numerator_fn(order, Ls)
                polys[t] = DC.poly_from(N, degree); qss[t] = qs
            if all(not polys[t].t for t in TP.TVALS):
                continue
            dN = Poly(4)
            for wl, t in zip(TP.LW, TP.TVALS):
                dN = dN + polys[t].scale(Q2(wl))
            qs0 = qss[0]
            dqs = [np.array([Fraction(int(qss[1][i][m]) - int(qs0[i][m])) for m in range(4)]) for i in range(len(qs0))]
            sgn = c0.parity_sign(order, canon) * (-1)
            out.append((DI.terms(polys[0], qs0, degree, dN=dN, dqs=dqs), qs0, mp.mpf(sgn) / (2 * scale)))
        return out
    z = np.zeros(4, dtype=int)
    t0 = time.time()
    G0 = conv_value(collect(ins(DC.QF('s0', z), z)), dps)
    Gmu = []; dG = mp.mpc(0)
    for mu in range(4):
        e = DC.E4[mu]
        Gmu.append(conv_value(collect(ins(DC.QF('smu', z, mu=mu), z)), dps))
        dG += conv_value(collect(None, deriv_build=lambda t, mu=mu, e=e: ins(DC.QF('smu', t * e, mu=mu), t * e)), dps)
    alpha = G0 - 1j * dG
    beta = [-1j * g for g in Gmu]
    p2 = -k - p1
    alphaQ = -mp.mpf(int(p2 @ p2)) / 2 * alpha + sum(int(p2[m]) * beta[m] for m in range(4))
    betaQ = [-mp.mpf(int(p2 @ p2)) / 2 * bm for bm in beta]
    cu = mp.mpf(1) / C.C_ADJ
    res = dict(alpha=alpha * cu, beta=[b_ * cu for b_ in beta], alphaQ=alphaQ * cu, betaQ=[b_ * cu for b_ in betaQ])
    # O structure on (A at p1: a, A at p2: a' = polarisation of the Q-image leg = bpol)
    TO = float(np.einsum('a,b,ab->', a[w1], bpol[w2], C.G.real) + np.einsum('a,b,ab->', a[w2], bpol[w1], C.G.real))
    res['T_O'] = TO
    res['time'] = time.time() - t0
    return res


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--p1', default='1,2,0,1'); ap.add_argument('--k', default='0,-1,1,2')
    ap.add_argument('--w', default='4,5'); ap.add_argument('--dps', type=int, default=15)
    ap.add_argument('--out', default='')
    A = ap.parse_args()
    p1 = np.array([int(x) for x in A.p1.split(',')]); k = np.array([int(x) for x in A.k.split(',')])
    w1, w2 = [int(x) for x in A.w.split(',')]
    a = np.zeros((10, C.D)); a[w1] = C.traceless(np.array([1, 2, -1, 0, 1, 0, 2, 0, -1][:C.D], dtype=float))
    a[w2] = C.traceless(np.array([0, 1, 1, 2, 0, -1, 1, 1, 0][:C.D], dtype=float)) if w2 != w1 else a[w1]
    b = np.zeros((10, C.D)); b[w1] = C.traceless(np.array([2, 0, 1, 1, 0, -1, 0, 1, 1][:C.D], dtype=float))
    b[w2] = C.traceless(np.array([1, -1, 2, 0, 1, 1, -1, 0, 2][:C.D], dtype=float)) if w2 != w1 else b[w1]
    r = y1a(p1, k, a, b, w1, w2, dps=A.dps)
    p2 = -k - p1
    out = {'p1': p1.tolist(), 'k': k.tolist(), 'p2': p2.tolist(), 'w': [w1, w2], 'T_O': r['T_O'],
           'alpha': str(r['alpha']), 'beta': [str(x) for x in r['beta']],
           'alphaQ/T_O': str(r['alphaQ'] / r['T_O']), 'betaQ/T_O': [str(x / r['T_O']) for x in r['betaQ']],
           'time': r['time']}
    print(json.dumps(out, indent=1))
    if A.out:
        json.dump(out, open(A.out, 'w'), indent=1)
