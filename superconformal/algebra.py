"""
algebra.py -- explicit realization of g = psl(4|4) (complexified) by conformal Killing data on V = R^4.

Even part:  conformal Killing vector fields v on V (polynomial, degree <= 2)  (+)  so(W_R^perp),
            the latter as linear vector fields y -> R y on W_R^perp.  Every even element is one
            vector field on V_10 = V (+) W_R^perp; the bracket is the vector-field bracket
            [v,w] = v(w) - w(v)  ([SW Prop 2.2 (3)-(4)] convention: generators = vector fields).
Odd part:   Pi(Sigma (+) Sigma^*):  Q in Sigma (constant spinor), eta in Sigma^* (S-type), realized as the
            conformal Killing spinor eps(x) = Q + rho(x) eta.
Brackets (derived and machine-checked in verify.py):
  [v, eps]   = v(eps) + 1/4 rho(d(v,-)) eps - 1/8 (div v) eps
  [eps1,eps2]= Gamma(eps1(x), eps2(x))|_V  (vector field on V)  +  cR * Rpart(eps1,eps2) in so(W_R^perp)
               (w_i, Rpart w_j) = (rho(e_i)rho(e_j) Q1, eta2) + (rho(e_i)rho(e_j) Q2, eta1), i != j in W_perp.
All arithmetic exact (sympy).
"""
import sympy as sp
import model as Mo

X = sp.symbols('x0:4')
Y = {w: sp.Symbol('y%d' % w) for w in Mo.WP}
COORD = list(X) + [Y[w] for w in Mo.WP]

SIGM = [sp.Matrix(16, 16, lambda a, b: sp.Rational(Mo.SIG[M][a][b][0].numerator, Mo.SIG[M][a][b][0].denominator)
                  + sp.I * sp.Rational(Mo.SIG[M][a][b][1].numerator, Mo.SIG[M][a][b][1].denominator)) for M in range(10)]
SIGBM = [sp.Matrix(16, 16, lambda a, b: sp.Rational(Mo.SIGB[M][a][b][0].numerator, Mo.SIGB[M][a][b][0].denominator)
                   + sp.I * sp.Rational(Mo.SIGB[M][a][b][1].numerator, Mo.SIGB[M][a][b][1].denominator)) for M in range(10)]

CR = sp.Rational(1, 2)   # placeholder; fixed by super-Jacobi in verify.py (set_cR)


def set_cR(v):
    global CR
    CR = v


class Even:
    """vector field on V_10: comps[M] polynomial in x (M in V) or linear in y (M in W_perp)."""
    def __init__(self, comps):
        self.c = [sp.expand(c) for c in comps]

    def __add__(self, o):
        return Even([a + b for a, b in zip(self.c, o.c)])

    def scale(self, s):
        return Even([s * a for a in self.c])

    def is_zero(self):
        return all(sp.expand(a) == 0 for a in self.c)

    def apply(self, f):
        return sp.expand(sum(self.c[M] * sp.diff(f, COORD[M]) for M in range(10)))

    def div4(self):
        return sp.expand(sum(sp.diff(self.c[mu], X[mu]) for mu in range(4)))

    def Omega(self):
        """Omega[M][N] = d_M v_N - d_N v_M (10d, orthonormal)."""
        return [[sp.expand(sp.diff(self.c[N], COORD[M]) - sp.diff(self.c[M], COORD[N])) for N in range(10)] for M in range(10)]

    def vpart(self):
        return [self.c[mu] for mu in range(4)]

    def Rmatrix(self):
        """R^i_j with v^i = R^i_j y^j on W_perp."""
        return {(i, j): sp.expand(sp.diff(self.c[i], Y[j])) for i in Mo.WP for j in Mo.WP if sp.diff(self.c[i], Y[j]) != 0}


class Odd:
    def __init__(self, Q, eta):
        self.Q = sp.Matrix(Q)
        self.eta = sp.Matrix(eta)

    def __add__(self, o):
        return Odd(self.Q + o.Q, self.eta + o.eta)

    def scale(self, s):
        return Odd(s * self.Q, s * self.eta)

    def is_zero(self):
        return all(sp.expand(a) == 0 for a in list(self.Q) + list(self.eta))

    def field(self):
        """eps(x) = Q + rho(x) eta (16-vector of polynomials in x)."""
        e = sp.Matrix(self.Q)
        for mu in range(4):
            e += X[mu] * (SIGBM[mu] * self.eta)
        return e.applyfunc(sp.expand)


def odd_from_field(eps):
    """decompose a Killing-spinor polynomial eps(x) = Q + rho(x) eta; raises if not of this form."""
    eps = eps.applyfunc(sp.expand)
    Q = eps.subs({x: 0 for x in X})
    lin = [eps.applyfunc(lambda f: sp.diff(f, X[mu])).subs({x: 0 for x in X}) for mu in range(4)]
    # rho(e_mu) eta = lin[mu]  ->  eta = sig[mu] lin[mu] (sig sigbar = 1)
    eta = SIGM[0] * lin[0]
    o = Odd(Q, eta)
    res = (o.field() - eps).applyfunc(sp.expand)
    if any(r != 0 for r in res):
        raise ValueError('not a conformal Killing spinor of the form Q + rho(x) eta')
    return o


def ebracket(v, w):
    return Even([sp.expand(v.apply(w.c[M]) - w.apply(v.c[M])) for M in range(10)])


def rho2(Om, s):
    """rho(sum_{M<N} Om_MN e^M ^ e^N) on Sigma-valued s: sum_{M<N} Om_MN sigbar^M sig^N s."""
    out = sp.zeros(16, 1)
    for M in range(10):
        for N in range(M + 1, 10):
            if Om[M][N] != 0:
                out += Om[M][N] * (SIGBM[M] * (SIGM[N] * s))
    return out


def eo_bracket(v, o):
    """[v, eps] = v(eps) + 1/4 rho(d(v,-)) eps - 1/8 (div v) eps."""
    eps = o.field()
    r = eps.applyfunc(lambda f: v.apply(f))
    r += sp.Rational(1, 4) * rho2(v.Omega(), eps)
    r += -sp.Rational(1, 8) * v.div4() * eps
    return odd_from_field(r.applyfunc(sp.expand))


def Gamma_fields(e1, e2):
    """Gamma(e1,e2)^M = (rho(e_M) e1, e2) = e2^T sig^M e1, e1, e2 in Sigma (polynomial entries)."""
    return [sp.expand((e2.T * SIGM[M] * e1)[0, 0]) for M in range(10)]


def Rpart(o1, o2):
    """R^i_j = (rho(e_i) rho(e_j) Q1, eta2) + (rho(e_i) rho(e_j) Q2, eta1), i != j in W_perp."""
    R = {}
    for i in Mo.WP:
        for j in Mo.WP:
            if i == j:
                continue
            val = (o2.eta.T * (SIGBM[i] * (SIGM[j] * o1.Q)))[0, 0] + (o1.eta.T * (SIGBM[i] * (SIGM[j] * o2.Q)))[0, 0]
            R[(i, j)] = sp.expand(val)
    return R


def oo_bracket(o1, o2, cR=None):
    if cR is None:
        cR = CR
    G = Gamma_fields(o1.field(), o2.field())
    comps = [G[mu] for mu in range(4)] + [0] * 6
    R = Rpart(o1, o2)
    for i in Mo.WP:
        comps[i] = sp.expand(cR * sum(R.get((i, j), 0) * Y[j] for j in Mo.WP))
    return Even(comps)


def Gamma_Wperp_fields(o1, o2):
    G = Gamma_fields(o1.field(), o2.field())
    return [G[w] for w in Mo.WP]


def bracket(a, b):
    if isinstance(a, Even) and isinstance(b, Even):
        return ebracket(a, b)
    if isinstance(a, Even) and isinstance(b, Odd):
        return eo_bracket(a, b)
    if isinstance(a, Odd) and isinstance(b, Even):
        return eo_bracket(b, a).scale(-1)
    return oo_bracket(a, b)


def parity(a):
    return 1 if isinstance(a, Odd) else 0


def add(a, b):
    return a + b


def is_zero(a):
    return a.is_zero()


# ---------------- standard generators ----------------
def P(a):
    return Even([a[mu] for mu in range(4)] + [0] * 6)


def Mrot(om):
    """rotation v^nu = om[nu][mu] x^mu (om antisymmetric)."""
    return Even([sum(om[nu][mu] * X[mu] for mu in range(4)) for nu in range(4)] + [0] * 6)


def Dil():
    """Delta of [SW (17)]: ad has weight -2 on translations -> vector field 2E = 2 x^mu d_mu."""
    return Even([2 * X[mu] for mu in range(4)] + [0] * 6)


def Kb(b):
    """K_b = b^mu K_mu, K_mu = |x|^2 d_mu - 2 x_mu E  [SW (3)]."""
    x2 = sum(x ** 2 for x in X)
    bx = sum(b[mu] * X[mu] for mu in range(4))
    return Even([x2 * b[nu] - 2 * bx * X[nu] for nu in range(4)] + [0] * 6)


def Rsym(R):
    """R: dict (i,j)->value with v^i = R^i_j y^j, antisymmetric."""
    comps = [0] * 10
    for i in Mo.WP:
        comps[i] = sum(R.get((i, j), 0) * Y[j] for j in Mo.WP)
    return Even(comps)


def Qgen(Q):
    return Odd(Q, [0] * 16)


def Sgen(eta):
    return Odd([0] * 16, eta)


# ---------------- conversion to engine data ----------------
def to_frac(z):
    from fractions import Fraction as Fr
    z = sp.nsimplify(z)
    re, im = sp.re(z), sp.im(z)
    return (Fr(int(sp.fraction(re)[0]), int(sp.fraction(re)[1])), Fr(int(sp.fraction(im)[0]), int(sp.fraction(im)[1])))


def even_to_engine(V, v):
    """returns (vf list of x-polys for engine, R dict for engine)."""
    import engine as E
    vf = []
    for mu in range(4):
        poly = sp.Poly(v.c[mu], *X) if v.c[mu] != 0 else None
        co = {}
        if poly is not None:
            for mon, cf in poly.terms():
                co[tuple(mon)] = to_frac(cf)
        vf.append(Mo.xpoly(V, co))
    Rd = {k: to_frac(val) for k, val in v.Rmatrix().items()}
    return vf, Rd


def odd_to_engine(o):
    return [to_frac(z) for z in o.Q], [to_frac(z) for z in o.eta]
