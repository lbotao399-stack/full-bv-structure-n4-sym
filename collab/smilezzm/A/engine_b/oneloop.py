"""Engine B -- one-loop singular parts (gates G1, G2) from the canonical Feynman rules, in the Schwinger
representation, exact (sympy).

Conversions (notes.txt, section F):  lambda_hat := g^2 N/(16 pi^2) = -2 unit;  for su(N) in the spec normalisation
int(A_perp, A_perp) = 2 hbar int phi^a_I phi^a_I,  int(dA_perp, dA_perp) = 2 hbar int (d phi)^2,
O_ab (manuscript) = 2 hbar (a.phi^c)(b.phi^c),  K (manuscript) = 2 hbar phi^c.phi^c.
A one-loop quadratic vertex V2 = 1/2 int int phi(u) Sigma(u-v) phi(v) (entering exp(+V2)) with plane-wave
Sigma^(p) = s_m/eps + s_l p^2 log eps + (finite) corresponds to
I^CT  contains  unit [X eps^-1 int(A_perp,A_perp) + Y log eps int(dA_perp,dA_perp)],  X = -s_m/(2 lh), Y = -s_l/(2 lh).
"""
import sympy as sp
from gauss import gaussian, vec_deriv

eps, t, T, u = sp.symbols('epsilon t T u', positive=True)
PI = sp.pi


def D_at_zero():
    """D_eps(0) = int_eps^inf (4 pi t)^-2 dt."""
    return sp.integrate((4 * PI * t) ** -2, (t, eps, sp.oo))


def bubble_log_coefficient():
    """Plane-wave value of D_eps^2:  Bhat(p) = (16 pi^2)^-1 int int_{[eps,inf)^2} (t1+t2)^-2 exp(-p^2 t1 t2/(t1+t2)).
    Boundary form: eps d/deps Bhat -> -2 eps int_eps^inf (eps+t)^-2 dt / (16 pi^2) = coefficient of log eps.
    (The O(p^2 eps) remainder of the pinned integral vanishes as eps -> 0.)"""
    pinned = -2 * eps * sp.integrate((eps + t) ** -2, (t, eps, sp.oo)) / (16 * PI ** 2)
    return sp.simplify(pinned)


def DK_moments():
    """D_eps K_eps(w) = (16 pi^2)^-1 int_eps^inf (t+eps)^-2 K_{c(t)}(w) dt,  c = t eps/(t+eps)  (Gaussian product).
    Plane wave: (16 pi^2)^-1 int (t+eps)^-2 exp(-p^2 c) dt = m0 + p^2 m2 + O(eps p^4):  m0 exact, m2 finite."""
    c = t * eps / (t + eps)
    m0 = sp.integrate((t + eps) ** -2, (t, eps, sp.oo)) / (16 * PI ** 2)
    m2 = -sp.integrate(sp.simplify(c / (t + eps) ** 2), (t, eps, sp.oo)) / (16 * PI ** 2)
    return sp.simplify(m0), sp.simplify(m2)


def check_DK_product():
    """Verify the Gaussian product identity K_t(w) K_eps(w) (4 pi)^0 ... symbolically."""
    w2 = sp.Symbol('w2', positive=True)
    Kt = (4 * PI * t) ** -2 * sp.exp(-w2 / (4 * t))
    Ke = (4 * PI * eps) ** -2 * sp.exp(-w2 / (4 * eps))
    c = t * eps / (t + eps)
    Kc = (4 * PI * c) ** -2 * sp.exp(-w2 / (4 * c))
    return sp.simplify(Kt * Ke - Kc / (16 * PI ** 2 * (t + eps) ** 2))


def g1_table():
    """Gate G1: singular parts of the one-loop A_perp self-energy graphs, in units (X: eps^-1 int(A,A),
    Y: log eps int(dA,dA)).  Colour-flavour weights from the census (both lines -> per line /2)."""
    lh = sp.Symbol('lh', positive=True)  # lambda_hat = g^2 N / (16 pi^2);  g^2 C_A = 16 pi^2 lh
    g2CA = 16 * PI ** 2 * lh
    D0 = D_at_zero()
    Blog = bubble_log_coefficient()  # coefficient of log eps in Bhat(p)
    m0, m2 = DK_moments()
    rows = {}
    # seagull with A|_V loop: census weight -8 (both lines) per g^2 C_A dim g S12, spatial factor D_eps(0)
    # -> Sigma = (w/2) g^2 C_A D_eps(0) delta
    for name, w in (('seagull_AV_loop', -8), ('seagull_Aperp_loop', -10)):
        sm = sp.simplify(sp.Rational(w, 2) * g2CA * D0 * eps)
        rows[name] = dict(X=sp.simplify(-sm / (2 * lh)), Y=sp.Integer(0))
    # bubble A|_V + A_perp: Sigma_b = g^2 C_A [-2 Laplacian(D^2) + D K]  ->  Sigma^ = g^2 C_A [2 p^2 Bhat + (DK)^]
    sm = sp.simplify(g2CA * m0 * eps)
    sl = sp.simplify(g2CA * 2 * Blog)
    rows['bubble_AV_Aperp'] = dict(X=sp.simplify(-sm / (2 * lh)), Y=sp.simplify(-sl / (2 * lh)))
    # gaugino loop: Sigma_F = 8 g^2 C_A (dD)^2 = g^2 C_A [4 Laplacian(D^2) + 8 D K]
    sm = sp.simplify(g2CA * 8 * m0 * eps)
    sl = sp.simplify(g2CA * (-4) * Blog)
    rows['gaugino_loop'] = dict(X=sp.simplify(-sm / (2 * lh)), Y=sp.simplify(-sl / (2 * lh)))
    tot = dict(X=sum(r['X'] for r in rows.values()), Y=sum(r['Y'] for r in rows.values()))
    rows['total'] = tot
    return rows, dict(D0=D0, Blog=Blog, DK_m0=m0, DK_m2=m2)


def triangle_log(edges, terms, k=2):
    """Log-eps coefficient of a superficially log-divergent one-loop vertex graph at zero external momenta:
    F(eps) = int_{[eps,inf)^n} f(t),  eps dF/deps = -sum_e int_{[1,inf)^(n-1)} f|_{t_e=1}  (f homogeneous of
    degree -n); exact (sympy)."""
    pref, nx, n1, phi = gaussian(k, edges, terms)
    f = sp.factor(pref * n1)
    ts = [e[2] for e in edges]
    L = 0
    for te in ts:
        others = [s for s in ts if s != te]
        g = f.subs(te, 1)
        for s in others:
            g = sp.simplify(sp.integrate(g, (s, 1, sp.oo)))
        L += -g
    return sp.simplify(L), f


def g2_table():
    """Gate G2: singular parts of the one-loop graphs with one O = (a.phi^c)(b.phi^c) on the cycle, converted to
    the manuscript basis unit log eps [alpha O + beta (a,b) K] (+ vacuum eps^-1 term)."""
    lh = sp.Symbol('lh', positive=True)
    g2CA = 16 * PI ** 2 * lh
    t1, t2, t3 = sp.symbols('t1 t2 t3', positive=True)
    Blog = bubble_log_coefficient()
    out = {}
    # (a) vacuum: <O> = (a,b) dim g D_eps(0); manuscript value 2 hbar times this
    out['a'] = dict(vacuum_coefficient_of_hbar_dimg_ab_over_16pi2_eps=sp.simplify(2 * D_at_zero() * 16 * PI ** 2 * eps))
    # (b) one external A|_V: colour f^{abc} delta^{bc} = 0 (checked numerically in census)
    out['b'] = dict(alpha=0, beta=0)
    # (e) bubble with V4: V^(e) = -g^2 C_A [(a,b) phi.phi - (a.phi)(b.phi)] D_eps(z)^2  (hand Wick count, notes F.4)
    #     Sing D_eps^2 -> Blog * log eps;   in units: lh * X = -2 unit * X
    ce_ab = -g2CA * Blog  # coefficient of (a,b) phi.phi log eps
    ce_O = g2CA * Blog    # coefficient of (a.phi)(b.phi) log eps
    out['e'] = dict(alpha=sp.simplify(ce_O / (-lh / 2)), beta=sp.simplify(ce_ab / (-lh / 2)))
    # (f) triangle: V^(f) = g^2 C_A int int dD(z1).dD(z2) D(z1-z2) (a.phi(z1))(b.phi(z2))
    Lf, ff = triangle_log([(0, 'o', t1), (1, 'o', t2), (0, 1, t3)],
                          [(1, [(vec_deriv(0, 'o', t1), vec_deriv(1, 'o', t2))])])
    out['f'] = dict(alpha=sp.simplify(g2CA * Lf / (-lh / 2)), beta=0, triangle_log=Lf)
    # (c) seagull with external A|_V: V^(c) = -g^2 C_A (a,b) D_eps^2 A.A  ->  coefficient of (a,b) A.A log eps
    cc = -g2CA * Blog
    # (d) triangle with external A|_V, scalar lines: V^(d) = g^2 C_A (a,b) int int A_m(z1) A_n(z2) [ ... ]_{mn},
    #     isotropic average 1/4 of the trace (zero external momenta).  Edges: O-z1 (t1), O-z2 (t2), z1-z2 (t3).
    d1o = vec_deriv(0, 'o', t1)
    d1z = vec_deriv(0, 1, t3)
    d2o = vec_deriv(1, 'o', t2)
    d2z = vec_deriv(1, 0, t3)
    Ld, fd = triangle_log([(0, 'o', t1), (1, 'o', t2), (0, 1, t3)],
                          [(1, [(d1z, d2z)]), (1, [(d1o, d2o)]), (-1, [(d1o, d2z)]), (-1, [(d1z, d2o)])])
    cd = g2CA * Ld / 4
    # manuscript: (w1,w2) int omega (A_V, A_V) = 2 hbar (a,b) A.A  ->  same conversion as the other columns
    out['c'] = dict(coef_ab_AVAV=sp.simplify(cc / (-lh / 2)))
    out['d'] = dict(coef_ab_AVAV=sp.simplify(cd / (-lh / 2)), triangle_log=Ld)
    out['total'] = dict(alpha=out['e']['alpha'] + out['f']['alpha'], beta=out['e']['beta'],
                        coef_ab_AVAV=out['c']['coef_ab_AVAV'] + out['d']['coef_ab_AVAV'])
    return out


if __name__ == '__main__':
    print(check_DK_product())
    print(g1_table())
    print(g2_table())
