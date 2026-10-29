# T2 audit, route "renormalised self-energy": the one-loop A|_V two-point function of
# W(P(eps,infty), I) computed exactly (Schwinger parameters on [eps, infty) per edge, exact eps -> 0
# asymptotics), split into transverse / longitudinal parts.  Since Q W^R = obstruction and, in the
# (c, A|_V) sector, Q of the (c, A^*) two-point vanishes ((dc, d^*dA) = 0) and Q of the transverse part
# vanishes, the (c, A|_V) obstruction is Q of the local longitudinal part of the renormalised A|_V
# two-point function.  Fresh script; imports only core.py.
import json
import warnings
import numpy as np
import sympy as sp
from fractions import Fraction
import core as C

np.seterr(all='ignore')
warnings.filterwarnings('ignore')

z, Lz, gam, eps = sp.symbols('z Lz gamma_E eps', positive=True)
xs = sp.symbols('x', positive=True)
LOGEPS = sp.Symbol('log_eps')
LOG2 = sp.log(2)


def gamma_neg_laurent(h):
    """Gamma(-h, z) up to o(1): dict {(r, s): coef} meaning coef * z^(-r) * (log z)^s (gamma_E symbolic).
    Gamma(-h,z) = ((-1)^h/h!) [E1(z) - e^{-z} sum_{j<h} (-1)^j j! z^(-j-1)],  E1(z) = -gamma_E - log z + O(z)."""
    pref = sp.Integer(-1) ** h / sp.factorial(h)
    d = {}
    def add(k, v):
        d[k] = d.get(k, 0) + pref * v
    add((0, 0), -gam)
    add((0, 1), -1)
    for j in range(h):
        for n in range(j + 2):
            r = j + 1 - n
            add((r, 0), -(sp.Integer(-1) ** j) * sp.factorial(j) * (sp.Integer(-1) ** n) / sp.factorial(n))
    return d


def J(j, h):
    """int_0^1 dx (1-x)^j (x(1-x))^h Gamma(-h, eps max(x,1-x)) up to o(1), in terms of eps, log_eps."""
    lap = gamma_neg_laurent(h)
    base = (1 - xs) ** j * (xs * (1 - xs)) ** h

    def xint(g):
        return (sp.integrate(sp.expand(base * g(1 - xs)), (xs, 0, sp.Rational(1, 2)))
                + sp.integrate(sp.expand(base * g(xs)), (xs, sp.Rational(1, 2), 1)))
    res = 0
    for (r, s), coef in lap.items():
        if coef == 0:
            continue
        if s == 0:
            res += coef * eps ** (-r) * xint(lambda m: m ** (-r))
        else:
            assert r == 0
            res += coef * (LOGEPS * xint(lambda m: 1) + xint(lambda m: sp.log(m)))
    return sp.expand(res)


_Jcache = {}


def Jc(j, h):
    if (j, h) not in _Jcache:
        _Jcache[(j, h)] = J(j, h)
    return _Jcache[(j, h)]


def unit(mu):
    v = np.zeros(4)
    v[mu] = 1
    return v


def G3(moms):
    return C.gamma3([np.asarray(m, dtype=float) for m in moms])


def vert_exact(a_bg, p_bg, P):
    """v(q + P) = Gamma3(p_bg; -(q+P+p_bg); q+P) as PM (tau ignored): linear part + constant part."""
    out = {}
    zz = np.zeros(4)
    for mu in range(4):
        e = unit(mu)
        L = np.einsum('i,iwu->wu', a_bg, G3([zz, -e, e]))
        if np.any(L != 0):
            out[(0, 0, C.unitvec(mu))] = L
    Cc = np.einsum('i,iwu->wu', a_bg, G3([p_bg, -(P + p_bg), P]))
    if np.any(Cc != 0):
        out[(0, 0, C.Z4)] = out.get((0, 0, C.Z4), 0) + Cc
    return out


def hatG_exact(kind, a):
    if kind == 'B':
        return {(0, 0, C.Z4): C.hatG_boson_const()}
    lin = C.hatG_fermion_lin() if kind == 'F' else C.hatG_ghost_lin()
    out = {(0, 0, C.unitvec(mu)): lin[mu] for mu in range(4)}
    const = sum(a[mu] * lin[mu] for mu in range(4))
    if np.any(const != 0):
        out[(0, 0, C.Z4)] = const
    return out


def numerator(factors, proj=None):
    P = factors[0]
    for f in factors[1:]:
        P = C.pm_mul(P, f, 0)
    poly = {}
    for (t, m, mono), M in P.items():
        tr = np.trace(M if proj is None else proj @ M)
        r, i = C.cfrac(tr)
        val = sp.Rational(r.numerator, r.denominator) + sp.I * sp.Rational(i.numerator, i.denominator)
        if val != 0:
            poly[mono] = poly.get(mono, 0) + val
    return poly


lsym = sp.symbols('l1:5')
ysym = sp.Symbol('y')


def schwinger_bubble(poly, shift):
    """int_q N(q) e^{-eps q^2} e^{-eps (q+shift)^2} / (q^2 (q+shift)^2), shift^2 = 1, returned as
    16 pi^2 x (expansion in eps^-1, log_eps, constants)."""
    assert abs(np.dot(shift, shift) - 1) < 1e-15
    qexpr = 0
    for mono, c in poly.items():
        t = c
        for mu in range(4):
            sh = Fraction(float(shift[mu]))
            t *= (lsym[mu] - ysym * sp.Rational(sh.numerator, sh.denominator)) ** mono[mu]
        qexpr += t
    P = sp.Poly(sp.expand(qexpr), *lsym, ysym)
    total = 0
    for monom, c in P.terms():
        al = monom[:4]
        j = monom[4]
        if any(a % 2 for a in al):
            continue
        h = sum(al) // 2
        w = 1
        for a in al:
            w *= C.dfact2(a)
        coef = c * w / sp.Integer(2) ** h
        total += coef * Jc(j, h)
    return sp.expand(total)


def tadpole(v12):
    """int_q tr[hatG_B v12] e^{-eps q^2}/q^2 = tr[2 v12]/(16 pi^2 eps); returns 16pi^2 x value."""
    r, i = C.cfrac(np.trace(2 * v12))
    assert i == 0
    return sp.Rational(r.numerator, r.denominator) / eps


def proj(indices, n):
    Pm = np.zeros((n, n), dtype=complex)
    for i in indices:
        Pm[i, i] = 1
    return Pm


VIDX = [C.bidx(M, a) for M in range(C.NV) for a in range(C.NC)]
WIDX = [C.bidx(M, a) for M in range(C.NV, C.ND) for a in range(C.NC)]
PV = proj(VIDX, C.NB)
PW = proj(WIDX, C.NB)


def two_point(a1, a2, p):
    """B(a1,a2;p) * 16 pi^2 per loop type (multilinear coefficient of W_1 on a1 e^{ipx} + a2 e^{-ipx})."""
    p1, p2 = p, -p
    out = {}
    # bosonic bubble: +1/2 Tr(G V1 G V2) = 1/2 int tr[G(q) v1(q+p2) G(q+p2) v2(q)]
    fB = [hatG_exact('B', np.zeros(4)), vert_exact(a1, p1, p2), hatG_exact('B', p2), vert_exact(a2, p2, np.zeros(4))]
    out['A|_V loop (cubic bubble)'] = sp.Rational(1, 2) * schwinger_bubble(numerator(fB, PV), p2)
    out['A|_Wperp loop (cubic bubble)'] = sp.Rational(1, 2) * schwinger_bubble(numerator(fB, PW), p2)
    # quartic tadpoles: +1/2 Tr(G V12)
    v12 = np.einsum('i,j,ijwu->wu', a1, a2, C.gamma4())
    out['A|_V loop (quartic tadpole)'] = sp.Rational(1, 2) * tadpole(PV @ v12)
    out['A|_Wperp loop (quartic tadpole)'] = sp.Rational(1, 2) * tadpole(PW @ v12)
    # gaugino bubble: -1/2 Tr
    fF = [hatG_exact('F', np.zeros(4)), {(0, 0, C.Z4): C.Y_of(a1)}, hatG_exact('F', p2), {(0, 0, C.Z4): C.Y_of(a2)}]
    out['lambda loop'] = -sp.Rational(1, 2) * schwinger_bubble(numerator(fF), p2)
    # ghost bubble: -1/2 Tr
    fG = [hatG_exact('G', np.zeros(4)), {(0, 0, C.Z4): C.ghost_V(a1)}, hatG_exact('G', p2),
          {(0, 0, C.Z4): C.ghost_V(a2)}]
    out['c-A^* loop'] = -sp.Rational(1, 2) * schwinger_bubble(numerator(fG), p2)
    return out


def split(expr):
    expr = sp.expand(expr)
    e_inv = expr.coeff(eps, -1)
    rest = sp.expand(expr - e_inv / eps)
    lg = rest.coeff(LOGEPS, 1)
    fin = sp.expand(rest - lg * LOGEPS)
    return sp.nsimplify(e_inv), sp.expand(lg), sp.expand(fin)


def in_units(expr):
    # unit = -2/(16 pi^2); expr is 16 pi^2 x value
    return sp.expand(expr / (-2))


def run(variant=0):
    if variant == 0:
        p = unit(0)
        tvec = unit(1)
        col = 0
    else:
        p = np.array([0.5, 0.5, 0.5, 0.5])
        tvec = np.array([0.5, -0.5, 0.5, -0.5])
        col = 1
    # longitudinal: a1 = a2 = p (x) t_col  -> (a1,a2) = 1, (p.a1)(p.a2) = 1
    aL = np.zeros((C.ND, C.NC))
    aL[:C.NV, col] = p
    aL = aL.reshape(-1).astype(complex)
    L = two_point(aL, aL, p)
    # transverse: a1 = a2 = tvec (x) t_col, tvec.p = 0, |tvec| = 1 -> (a1,a2) = 1, tree coefficient = -1/2
    aT = np.zeros((C.ND, C.NC))
    aT[:C.NV, col] = tvec
    aT = aT.reshape(-1).astype(complex)
    T = two_point(aT, aT, p)
    tree_T = C.tree_quadratic_A(aT, aT, p)
    tree_L = C.tree_quadratic_A(aL, aL, p)
    rows = {}
    totL = 0
    totT = 0
    for k in L:
        eL, lL, fL = split(in_units(L[k]))
        eT, lT, fT = split(in_units(T[k]))
        totL += in_units(L[k])
        totT += in_units(T[k])
        # B/u = 2 (m/eps + beta) on longitudinal pols ; m from eps^-1 / 2 ; beta from finite / 2
        rows[k] = {
            'longitudinal pol: eps^-1 coefficient of the W_1 two-point (in the unit)': str(eL),
            'longitudinal pol: log eps coefficient of the W_1 two-point (in the unit)': str(sp.simplify(lL)),
            'longitudinal pol: finite part of the W_1 two-point (in the unit)': str(sp.simplify(fL)),
            'coefficient of unit x eps^-1 int(A|_V,A|_V)': str(sp.simplify(eL / 2)),
            'coefficient of unit x int(d^*A|_V,d^*A|_V) (scheme-free only when the log coefficient is 0)': str(sp.simplify(fL / 2)),
            'transverse pol: eps^-1 coefficient of the W_1 two-point (in the unit)': str(eT),
            'transverse pol: log eps coefficient as a multiple of the tree -1/4(dA,dA)': str(sp.simplify(lT / sp.Rational(-1, 2))),
        }
    eL, lL, fL = split(totL)
    eT, lT, fT = split(totT)
    summary = {
        'tree -1/4(dA,dA) on transverse pol': str(tree_T), 'tree on longitudinal pol': str(tree_L),
        'TOTAL longitudinal pol: eps^-1 coefficient (in the unit)': str(eL),
        'TOTAL longitudinal pol: log eps coefficient (in the unit)': str(sp.simplify(lL)),
        'TOTAL longitudinal pol: finite part (in the unit)': str(sp.simplify(fL)),
        'TOTAL coefficient of unit x eps^-1 int(A|_V,A|_V)': str(sp.simplify(eL / 2)),
        'TOTAL coefficient of unit x int(d^*A|_V,d^*A|_V) in the renormalised one-loop functional': str(sp.simplify(fL / 2)),
        'TOTAL transverse pol: eps^-1 coefficient (in the unit)': str(eT),
        'TOTAL transverse pol: log eps coefficient as a multiple of the tree -1/4(dA,dA)': str(sp.simplify(lT / sp.Rational(-1, 2))),
    }
    return rows, summary


if __name__ == '__main__':
    import sys
    v = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    rows, summary = run(v)
    print(json.dumps({'variant': v, 'unit': 'hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2) (hbar = 1; su(2), (t_a,t_b) = delta_ab: -2/(16 pi^2))', 'gamma_E': 'Euler constant', 'per_loop': rows, 'summary': summary}, indent=1))
