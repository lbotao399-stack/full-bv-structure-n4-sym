"""Audit #3: the scheme that removes P_1 (delta_theory.tex prop:d-P (ii)) and the (Delta,P) component of L o L = 0.

Model: a scalar letter phi of Weyl weight 1 (stand-in for iota_w A), local operators int omega g with omega an arbitrary
test function and g a polynomial in the jets of phi with explicit x-dependent coefficients.  Two such operators are
equal for all omega iff their densities are equal, so all comparisons are pointwise identities of densities.
Classical actions (def:c7-vdot, weight term 1/4 div v on A|_Wperp):
  P_b : dphi = b.d phi ;  Delta (v = 2E): dphi = 2 x.d phi + 2 phi ;
  K_b : v = |x|^2 b - 2(b.x) x,  dphi = v.d phi - 2 (b.x) phi ;  rotation: v = Omega x.
ad_v g = sum_alpha (dg/dJ_alpha) D^alpha(dphi)  (jet calculus; ad acts on the fields, never on explicit x).
Weyl weight of a monomial = sum over jet factors of (1 + order); explicit x does not count (lem:c6-data).
Schemes compared:  psi_w = c x (Weyl weight)  [the text of prop:d-P (ii)]  and  psi_D = (c/2) ad_Delta = c(E + w).
Change of the component along x (prop:c8-laws (b), as written in the proof of prop:d-P): ch_x = ad_x psi g - psi ad_x g.
The text claims: P: -c, Delta, rotations: 0, K: +c times the classical transformation.
Also: order-hbar part of [Lambda(Delta),Lambda(P)] + Lambda([Delta,P]) with Lambda(x) = (1 + hbar kappa_x) ad_x.
"""
import itertools, json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
X = sp.symbols('x0:4'); B = sp.symbols('b0:4'); c, hb = sp.symbols('c hbar')
kD, kP = sp.symbols('kappa_Delta kappa_P')
MAXO = 3
alphas = [a for n in range(MAXO + 1) for a in itertools.product(range(n + 1), repeat=4) if sum(a) == n]
alphas = sorted(set(alphas), key=lambda a: (sum(a), a))
Jsym = {a: sp.Symbol('J_' + ''.join(map(str, a))) for a in alphas}
WEIGHT = {Jsym[a]: 1 + sum(a) for a in alphas}


def shift(a, mu):
    b = list(a); b[mu] += 1; return tuple(b)


def Dtot(expr, mu):
    out = sp.diff(expr, X[mu])
    for a in alphas:
        if sum(a) < MAXO:
            out += Jsym[shift(a, mu)] * sp.diff(expr, Jsym[a])
    return sp.expand(out)


def Dmulti(expr, a):
    for mu in range(4):
        for _ in range(a[mu]):
            expr = Dtot(expr, mu)
    return expr


def ad(v, w, g):
    dphi = sum(v[m] * Jsym[shift((0, 0, 0, 0), m)] for m in range(4)) + w * Jsym[(0, 0, 0, 0)]
    out = 0
    for a in alphas:
        dg = sp.diff(g, Jsym[a])
        if dg != 0:
            out += dg * Dmulti(dphi, a)
    return sp.expand(out)


def weight_grade(g):
    """returns dict weight -> part of g"""
    g = sp.expand(g); parts = {}
    for term in sp.Add.make_args(g):
        wt = 0
        for f, e in term.as_powers_dict().items():
            if f in WEIGHT:
                wt += WEIGHT[f] * e
        parts[wt] = parts.get(wt, 0) + term
    return parts


def psi_w(g):
    return sp.expand(sum(c * wt * part for wt, part in weight_grade(g).items()))


xx = sum(x * x for x in X); bx = sum(b * x for b, x in zip(B, X))
GEN = {
    'P': ([B[m] for m in range(4)], sp.Integer(0)),
    'Delta': ([2 * X[m] for m in range(4)], sp.Integer(2)),
    'K': ([xx * B[m] - 2 * bx * X[m] for m in range(4)], -2 * bx),
    'rot': ([X[1], -X[0], 0, 0], sp.Integer(0)),
}


def adg(name, g):
    v, w = GEN[name]
    return ad(v, w, g)


def psi_D(g):
    return sp.expand(c / 2 * adg('Delta', g))


def is_multiple(a, b):
    """is a = lam * b with lam independent of x, b and the jets?  returns lam or None"""
    a = sp.expand(a); b = sp.expand(b)
    if b == 0:
        return 0 if a == 0 else None
    gens = list(Jsym.values()) + list(X) + list(B)
    Pb = sp.Poly(b, *gens); Pa = sp.Poly(a, *gens)
    m0 = Pb.monoms()[0]
    lam = sp.simplify(Pa.coeff_monomial(m0) / Pb.coeff_monomial(m0))
    return lam if sp.expand(a - lam * b) == 0 else None


def main():
    res = {}
    g0 = Jsym[(0, 0, 0, 0)] ** 2                     # (iota_w1 A, iota_w2 A)-type density, Weyl weight 2
    g1 = Jsym[(0, 0, 0, 0)] * Jsym[(1, 0, 0, 0)]     # a weight-3 density, second test
    for gname, g in (('phi^2', g0), ('phi d_0 phi', g1)):
        r = {}
        # classical brackets
        for x in ('P', 'K', 'rot'):
            com = sp.expand(adg('Delta', adg(x, g)) - adg(x, adg('Delta', g)))
            r['[ad_Delta, ad_%s] / ad_%s' % (x, x)] = str(is_multiple(com, adg(x, g)))
        # scheme changes
        for sname, psi in (('psi_w = c x Weyl weight', psi_w), ('psi_D = (c/2) ad_Delta', psi_D)):
            rr = {}
            for x in ('P', 'Delta', 'K', 'rot'):
                ch = sp.expand(adg(x, psi(g)) - psi(adg(x, g)))
                lam = is_multiple(ch, adg(x, g)) if adg(x, g) != 0 else None
                rr[x] = {'change / (classical transformation)': (str(lam / c) + ' c') if lam is not None else 'NOT a multiple'}
                if lam is None:
                    rr[x]['change (density)'] = str(sp.factor(ch))
                    rr[x]['classical transformation (density)'] = str(sp.factor(adg(x, g)))
            r[sname] = rr
        res[gname] = r
    # weight decomposition of the classical K transformation of phi^2 (ch07 table: output weights w_O+1, w_O, w_O-1)
    res['Weyl weights present in ad_K(phi^2)'] = [int(w) for w in sorted(weight_grade(adg('K', g0)).keys())]
    res['Weyl weights present in ad_P(phi^2)'] = [int(w) for w in sorted(weight_grade(adg('P', g0)).keys())]
    # (Delta, P) component of L o L = 0 at order hbar on the family int omega phi^2:
    # Lambda(x) = (1 + hbar kappa_x) ad_x, [Delta,P] = -2P (ad_Delta eigenvalue of P, ch07 table)
    LD = lambda h: sp.expand((1 + hb * kD) * adg('Delta', h))
    LP = lambda h: sp.expand((1 + hb * kP) * adg('P', h))
    mc = sp.expand(LD(LP(g0)) - LP(LD(g0)) + (-2) * LP(g0))
    res['[Lambda(Delta),Lambda(P)] + Lambda([Delta,P]) on phi^2, order hbar^0'] = str(sp.factor(mc.coeff(hb, 0)))
    res['... order hbar^1 / ad_P(phi^2)'] = str(is_multiple(mc.coeff(hb, 1), adg('P', g0)))
    res['reading'] = ('The order-hbar part of the (Delta,P) component is 2 kappa_Delta ad_P O.  For O in the 20 of '
                      'so(W_R^perp), ad_P O is not Lambda^(0)-exact (no local functional of degree -1 and Weyl weight 3 '
                      'is a colour singlet), so kappa_Delta = 0 on the cohomology for every one-loop quantisation: '
                      'the multiple of {S^(1)(Delta),O} in D_1 must vanish, independently of P_1 (kappa_P drops out).')
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'check_scheme_P1.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
