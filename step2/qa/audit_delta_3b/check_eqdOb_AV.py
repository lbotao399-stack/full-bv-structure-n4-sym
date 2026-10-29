"""Audit #3: independent test of delta_theory eq:d-Ob, Ob(Delta) = 4 C_log + Q Pi_Psi, in the output sector
(A|_V, A|_V), with engine A's composite cycle Pi_int (the only input) against engine B's Ob(Delta) (the target).

Inputs (step2/engine_a/delta_results.json, sector A^*_V,A_V; canonical form: field leg A_V at p, polarisation a,
dependent leg A^*_V at p' = -p, polarisation b; Lo = (p.a)(p.b)/|p|^2; L3 = log|p|^2 + gamma + log 3 after removing the
singular part log eps x X_V):
   alpha = (a.b)(13/6 + L3) + 2 Lo,   beta^mu = (10/3) p^mu (a.b) + (L3 - 1/2) a^mu (p.b) - 2 p^mu Lo.
Q on the A^*_V leg: b_nu -> g_{nu rho}(p') b'_rho, g(q) = -1/2(|q|^2 - q q^T) (the measured -1/2 d^*d of engine A);
the canonical-form rule for a function of the dependent momentum: alpha_Q = g alpha - (d_mu g) beta^mu, beta_Q = g beta.
Both legs are then A_V; the functional sees VV = D + D^T with alpha^T(p) = alpha(-p)[a<->b'] - div beta(-p)[a<->b'],
beta^T(p) = beta(-p)[a<->b'].
Calibration: the same procedure on engine A's A_perp entry (alpha = 8/3, beta = 10/3 p) must give VV = (2/3)|p|^2,
i.e. (1/3) int(dA_P,dA_P) (VV of int(dA_P,dA_P) = 2|p|^2), the cross-check quoted in the theory note.
Prediction of eq:d-Ob with engine B's Ob(Delta)|_(A_V,A_V) = -(11/3) int(dA_V,dA_V) + 0 int(d^*A_V,d^*A_V) and
4 C_log|_(A_V,A_V) = -4 int(dA_V,dA_V):   VV[Q Pi_Psi]_(A_V,A_V) = (1/3) x 2(|p|^2 a.b' - p.a p.b'),   beta part 0.
"""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
P = sp.symbols('p0:4', real=True); A = sp.symbols('a0:4'); Bp = sp.symbols('bp0:4'); Bv = sp.symbols('b0:4')
C0 = sp.Symbol('C0')    # gamma + log 3
dot = lambda u, v: sum(x * y for x, y in zip(u, v))


def make_kernel(alpha_fn, beta_fn):
    """returns functions (alpha_Q, beta_Q) of (q, first-leg pol, second-leg pol) after Q on the second (A^*) leg"""
    def g(q, nu, rho):
        return -sp.Rational(1, 2) * ((dot(q, q) if nu == rho else 0) - q[nu] * q[rho])

    def dg(q, mu, nu, rho):
        return -sp.Rational(1, 2) * (2 * q[mu] * (1 if nu == rho else 0) - (1 if mu == nu else 0) * q[rho]
                                     - q[nu] * (1 if mu == rho else 0))

    def kernel(q, pol1, pol2):
        qd = [-x for x in q]
        # alpha, beta linear in the A^* polarisation b: coefficients by differentiation
        al = alpha_fn(q, pol1, Bv); be = [beta_fn(q, pol1, Bv, mu) for mu in range(4)]
        Acoef = [sp.diff(al, Bv[nu]) for nu in range(4)]
        Bcoef = [[sp.diff(be[mu], Bv[nu]) for nu in range(4)] for mu in range(4)]
        aQ = sum(pol2[rho] * (g(qd, nu, rho) * Acoef[nu]
                              - sum(dg(qd, mu, nu, rho) * Bcoef[mu][nu] for mu in range(4)))
                 for nu in range(4) for rho in range(4))
        bQ = [sum(pol2[rho] * g(qd, nu, rho) * Bcoef[mu][nu] for nu in range(4) for rho in range(4)) for mu in range(4)]
        return sp.simplify(aQ), [sp.simplify(x) for x in bQ]
    return kernel


def vv(kernel, p, a, b2):
    Q = sp.symbols('q0:4', real=True)
    aQ, bQ = kernel(list(p), a, b2)
    aS, bS = kernel(list(Q), b2, a)
    div = sum(sp.diff(bS[mu], Q[mu]) for mu in range(4))
    sub = {Q[m]: -p[m] for m in range(4)}
    aT = sp.simplify(aS.subs(sub) - div.subs(sub)); bT = [sp.simplify(x.subs(sub)) for x in bS]
    return sp.simplify(aQ + aT), [sp.simplify(bQ[m] + bT[m]) for m in range(4)]


def main():
    res = {}
    # calibration: A_perp (scalar polarisations: (a,b) -> a.b with the same machinery, no transverse projector)
    def al_perp(q, a, b):
        return sp.Rational(8, 3) * dot(a, b)

    def be_perp(q, a, b, mu):
        return sp.Rational(10, 3) * q[mu] * dot(a, b)

    def kernel_perp(q, pol1, pol2):
        qd = [-x for x in q]
        al = al_perp(q, pol1, pol2); be = [be_perp(q, pol1, pol2, mu) for mu in range(4)]
        g = -sp.Rational(1, 2) * dot(qd, qd); dgm = [-qd[mu] for mu in range(4)]   # d/dq' of -1/2 |q'|^2
        return sp.simplify(g * al - sum(dgm[mu] * be[mu] for mu in range(4))), [sp.simplify(g * x) for x in be]
    va, vb = vv(kernel_perp, P, A, Bp)
    res['calibration A_perp: VV alpha'] = str(sp.factor(va))
    res['calibration A_perp: VV beta'] = [str(x) for x in vb]
    res['calibration target'] = '(2/3)|p|^2 (a.b\'), beta 0'
    # A_V sector
    pp = dot(P, P)

    def L3(q):
        return sp.log(dot(q, q)) + C0

    def al_V(q, a, b):
        return dot(a, b) * (sp.Rational(13, 6) + L3(q)) + 2 * dot(q, a) * dot(q, b) / dot(q, q)

    def be_V(q, a, b, mu):
        return (sp.Rational(10, 3) * q[mu] * dot(a, b) + (L3(q) - sp.Rational(1, 2)) * a[mu] * dot(q, b)
                - 2 * q[mu] * dot(q, a) * dot(q, b) / dot(q, q))
    kV = make_kernel(al_V, be_V)
    va, vb = vv(kV, P, A, Bp)
    va = sp.simplify(sp.expand(va)); vb = [sp.simplify(sp.expand(x)) for x in vb]
    target = sp.Rational(1, 3) * 2 * (pp * dot(A, Bp) - dot(P, A) * dot(P, Bp))
    res['A_V: VV alpha'] = str(sp.factor(va))
    res['A_V: VV beta'] = [str(sp.factor(x)) for x in vb]
    res['A_V: target alpha (eq:d-Ob with engine B)'] = str(sp.factor(target))
    diff = sp.simplify(va - target)
    res['A_V: VV alpha - target'] = str(sp.factor(diff))
    # decompose the result on the local structures |p|^2 a.b', p.a p.b' (coefficients must be constants if local)
    x1, x2 = sp.symbols('x1 x2')
    sol = sp.solve([sp.expand(va - (x1 * 2 * (pp * dot(A, Bp) - dot(P, A) * dot(P, Bp)) + x2 * 2 * dot(P, A) * dot(P, Bp)))
                    .subs({P[0]: 1, P[1]: 2, P[2]: 0, P[3]: 1, A[0]: 1, A[1]: 0, A[2]: 2, A[3]: -1,
                           Bp[0]: 0, Bp[1]: 1, Bp[2]: 1, Bp[3]: 3}),
                    sp.expand(va - (x1 * 2 * (pp * dot(A, Bp) - dot(P, A) * dot(P, Bp)) + x2 * 2 * dot(P, A) * dot(P, Bp)))
                    .subs({P[0]: 2, P[1]: -1, P[2]: 1, P[3]: 3, A[0]: 0, A[1]: 1, A[2]: -1, A[3]: 2,
                           Bp[0]: 1, Bp[1]: 1, Bp[2]: 0, Bp[3]: -2})], [x1, x2], dict=True)
    res['A_V: fit on [int(dA_V,dA_V), int(d^*A_V,d^*A_V)] at two samples (must be momentum independent)'] = \
        {str(k): str(v) for k, v in sol[0].items()} if sol else 'no solution'
    # (c, A^*_V) output: Q on the A_V field leg of the (A^*_V,A_V) entry (a_mu -> sigma1 p_mu gamma, function of the
    # fixed momentum: no delta' correction) plus Q on the dependent c^* leg of the (c^*,c) entry (w -> sigma2 (p'.b'),
    # delta' correction -(d/dp'^mu (p'.b')) beta^mu = -b'.beta).  Engine A: (c^*,c) entry alpha = -8/3, beta = p(-5/6 - L3)
    # (gamma,w) after removing log eps; int(Ec,c^*) <-> beta = +p (gamma,w); int(c^*,c) <-> alpha = s_c (gamma,w).
    # Calibration (engine B): Q int(A_V,A^*_V) and Q int(c^*,c) have the same (c,A^*_V) image -int(dc,A^*),
    # hence sigma1 (p.b') = sigma2 s_c (-p.b'), i.e. sigma2 = -sigma1 s_c.
    s1, sc = sp.symbols('sigma1 s_c')
    L3p = sp.log(pp) + C0
    alV = sp.expand(dot(P, Bp) * (sp.Rational(13, 6) + L3p) + 2 * dot(P, P) * dot(P, Bp) / pp)      # a -> p
    beV = [sp.expand(sp.Rational(10, 3) * P[m] * dot(P, Bp) + (L3p - sp.Rational(1, 2)) * P[m] * dot(P, Bp)
                     - 2 * P[m] * dot(P, P) * dot(P, Bp) / pp) for m in range(4)]
    s2 = -s1 * sc
    al_c, be_c = sp.Rational(-8, 3), [P[m] * (-sp.Rational(5, 6) - L3p) for m in range(4)]
    aQc = s2 * (dot([-x for x in P], Bp) * al_c - sum(Bp[m] * be_c[m] for m in range(4)))
    bQc = [s2 * dot([-x for x in P], Bp) * be_c[m] for m in range(4)]
    tot_a = sp.simplify(s1 * alV + aQc); tot_b = [sp.simplify(s1 * beV[m] + bQc[m]) for m in range(4)]
    image_cc = sp.simplify(s2 * dot([-x for x in P], Bp) * sc)       # Q int(c^*,c) = -int(dc,A^*) in this normalisation
    for scv in (1, -1):
        ta = sp.simplify(tot_a.subs(sc, scv)); tb = [sp.simplify(x.subs(sc, scv)) for x in tot_b]
        res['(c,A^*_V) output, s_c = %d: alpha' % scv] = str(sp.factor(ta))
        res['(c,A^*_V) output, s_c = %d: beta' % scv] = [str(sp.factor(x)) for x in tb]
        res['(c,A^*_V) output, s_c = %d: alpha / image of int(c^*,c) (prediction 2/3)' % scv] = \
            str(sp.simplify(ta / image_cc.subs(sc, scv)))
    res['reading'] = ('(A_V,A_V): VV[Q Pi_int] = (1/3) int(dA_V,dA_V) + 0 int(d^*A_V,d^*A_V) exactly: the log|p|^2 and '
                      '(p.a)(p.b)/|p|^2 terms of Pi_int drop out, eq:d-Ob holds there with engine A input and engine B '
                      'target.  (c,A^*_V): only s_c = +1 gives a local image (logs cancel), and then Q Pi_int = (2/3) Q '
                      'int(c^*,c) = -(2/3) int(dc,A^*), the value eq:d-Ob requires.  Pi_Psi itself is nonlocal in these '
                      'sectors (def:d-PiPsi calls it local).')
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'check_eqdOb_AV.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
