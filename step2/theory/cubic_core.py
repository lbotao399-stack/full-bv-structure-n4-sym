"""cubic_core.py -- exact plane-wave machinery for the theory agent (Step 2, Part C).

Functionals are evaluated on finite sums of plane-wave modes e = sum_m kappa_m b_m, where b_m is a numeric section
(slot, fibre vector, colour vector, momentum k: b_m = fibre * colour * e^{i(k,x)}) and kappa_m a generator of a
supercommutative algebra (all generators square-zero; parity = total parity of the slot).  The (2 pi)^4 delta of
momentum conservation is stripped: a local functional int dvol(density) evaluated on modes is the sum over mode tuples
with total momentum 0 of (product of coefficients in density order) * (numeric density value), d -> i k.

Conventions (manuscript ch01, ch03):
  Q F       = left d/dr of F(e + r Q e), Q on numeric sections: c -> -dc (A|_V slot), A -> -1/2 d^*dA (A^* slot),
              A^* -> d^*(A^*|_V) (c^* slot), lambda -> dslash lambda (lambda^* slot)            [Prop c1-Q]
  {f,g}     = Delta_0(fg) - Delta_0(f) g - (-1)^{|f|} f Delta_0(g), Delta_0 = sum d_x d_xi           [Def c1-bracket]
              x = A, lambda^*, c^* (even members), xi = A^*, lambda, c, <x-mode at k, xi-mode at -k> = 1
              => {f,g} = sum [ (d_xi f)(d_x g) + (-1)^{|f|} (d_x f)(d_xi g) ]   (left derivatives)
  K_0       = - sum (b^x_k (x) b^xi_{-k} + b^xi_{-k} (x) b^x_k)
  P(0,inf)  = int_0^inf (Q^GF (x) 1) K_t dt,  Q^GF: A^* -> -2A^*, A -> -d^*(A|_V), c^* -> dc^*, lambda^* -> -dslash lambda^*
  d_P       = 1/2 sum d_{phi''} d_{phi'}   [C07 6.4]
Realisation: V_10 orthonormal e_0..e_9 (V = 0..3, W_R^perp = 4..9); Sigma, Sigma^* 16-dim with rho(e_M) = sig[M]
(Sigma -> Sigma^*), sigbar[M] (Sigma^* -> Sigma) from classical/n4_model.py, pairing (s, s*) = s^T s*;
g = su(2), orthonormal basis T_a, [T_a,T_b] = eps_abc T_c, (T_a,T_b) = delta_ab.
"""
import sys
import itertools
from fractions import Fraction as Fr

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/classical')
import n4_model as nm  # noqa: E402


# ------------------------------------------------------------------ Gaussian rationals
class GQ:
    __slots__ = ('re', 'im')

    def __init__(self, re=0, im=0):
        self.re = re if isinstance(re, Fr) else Fr(re)
        self.im = im if isinstance(im, Fr) else Fr(im)

    def __add__(self, o):
        if not isinstance(o, GQ):
            o = GQ(o)
        return GQ(self.re + o.re, self.im + o.im)
    __radd__ = __add__

    def __sub__(self, o):
        if not isinstance(o, GQ):
            o = GQ(o)
        return GQ(self.re - o.re, self.im - o.im)

    def __rsub__(self, o):
        return GQ(o) - self

    def __neg__(self):
        return GQ(-self.re, -self.im)

    def __mul__(self, o):
        if not isinstance(o, GQ):
            o = Fr(o)
            return GQ(self.re * o, self.im * o)
        return GQ(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
    __rmul__ = __mul__

    def __truediv__(self, o):
        if not isinstance(o, GQ):
            o = GQ(o)
        n = o.re * o.re + o.im * o.im
        return GQ((self.re * o.re + self.im * o.im) / n, (self.im * o.re - self.re * o.im) / n)

    def iszero(self):
        return self.re == 0 and self.im == 0

    def __repr__(self):
        if self.im == 0:
            return str(self.re)
        return f'({self.re}+{self.im}i)'


ZERO = GQ(0)
ONE = GQ(1)
IU = GQ(0, 1)

# ------------------------------------------------------------------ realisation
V4 = (0, 1, 2, 3)
WP = (4, 5, 6, 7, 8, 9)


def _sparse(re, im):
    out = []
    for a in range(16):
        row = []
        for b in range(16):
            r, i = int(re[a, b]), int(im[a, b])
            if r or i:
                row.append((b, GQ(r, i)))
        out.append(row)
    return out


SIG = [_sparse(nm.SIG_RE[M], nm.SIG_IM[M]) for M in range(10)]      # rho(e_M): Sigma -> Sigma^*
SIGB = [_sparse(nm.SIGBAR_RE[M], nm.SIGBAR_IM[M]) for M in range(10)]  # rho(e_M): Sigma^* -> Sigma


def matvec(S, v):
    return [sum((c * v[b] for (b, c) in S[a]), ZERO) for a in range(16)]


def dot(u, v):
    return sum((a * b for a, b in zip(u, v)), ZERO)


def spin_pair(u, v):  # (s, s*) = s^T s*
    return dot(u, v)


def rho_vec(q, v, bar=False):
    """rho(q) v for q in V (4-vector of GQ or Fr), v in Sigma (bar=False) or Sigma^* (bar=True)."""
    S = SIGB if bar else SIG
    out = [ZERO] * 16
    for mu in V4:
        if q[mu] == 0:
            continue
        w = matvec(S[mu], v)
        out = [o + q[mu] * x for o, x in zip(out, w)]
    return out


def vol4(v):
    """rho(v_0 v_1 v_2 v_3) on Sigma: sigbar0 sig1 sigbar2 sig3."""
    w = matvec(SIG[3], v)
    w = matvec(SIGB[2], w)
    w = matvec(SIG[1], w)
    return matvec(SIGB[0], w)


# su(2), orthonormal
def br(X, Y):
    return [X[1] * Y[2] - X[2] * Y[1], X[2] * Y[0] - X[0] * Y[2], X[0] * Y[1] - X[1] * Y[0]]


def cp(X, Y):
    return sum((a * b for a, b in zip(X, Y)), Fr(0))


def f3(X, Y, Z):  # (X,[Y,Z]) = ([X,Y],Z)
    return cp(X, br(Y, Z))


DIMG = 3
SLOTS = ('c', 'A', 'lam', 'As', 'lams', 'cs')
SPAR = {'c': 1, 'A': 0, 'lam': 1, 'As': 1, 'lams': 0, 'cs': 0}
SDIM = {'c': 1, 'A': 10, 'lam': 16, 'As': 10, 'lams': 16, 'cs': 1}
XMEMBER = {'A': 'As', 'lams': 'lam', 'cs': 'c'}          # even member -> odd partner

# ------------------------------------------------------------------ algebra
PAR = {}       # generator id -> parity
KIND = {}      # generator id -> 'ext' | 'r' | 'int'
_next = [0]


def new_gen(parity, kind):
    g = _next[0]
    _next[0] += 1
    PAR[g] = parity
    KIND[g] = kind
    return g


def mono_mul(m1, m2):
    s1 = set(m1)
    for g in m2:
        if g in s1:
            return 0, None
    sgn = 0
    for b in m2:
        if PAR[b]:
            for a in m1:
                if a > b and PAR[a]:
                    sgn ^= 1
    return (-1 if sgn else 1), tuple(sorted(m1 + m2))


def nint(m):
    return sum(1 for g in m if KIND[g] == 'int')


def amul(x, y, maxint=2):
    out = {}
    for m1, c1 in x.items():
        for m2, c2 in y.items():
            s, m = mono_mul(m1, m2)
            if m is None or nint(m) > maxint:
                continue
            v = c1 * c2
            if s < 0:
                v = -v
            out[m] = out.get(m, ZERO) + v
    return {m: c for m, c in out.items() if not c.iszero()}


def aadd(x, y, s=1):
    out = dict(x)
    for m, c in y.items():
        out[m] = out.get(m, ZERO) + (c if s == 1 else -c)
    return {m: c for m, c in out.items() if not c.iszero()}


def ascale(x, s):
    return {m: c * s for m, c in x.items() if not (c * s).iszero()}


def dleft(x, g):
    """left derivative d/dg."""
    out = {}
    for m, c in x.items():
        if g not in m:
            continue
        j = m.index(g)
        sgn = 1
        if PAR[g]:
            if sum(PAR[h] for h in m[:j]) % 2:
                sgn = -1
        mm = m[:j] + m[j + 1:]
        out[mm] = out.get(mm, ZERO) + (c if sgn == 1 else -c)
    return {m: c for m, c in out.items() if not c.iszero()}


def gen_elem(g):
    return {(g,): ONE}


# ------------------------------------------------------------------ modes
class Mode:
    __slots__ = ('coef', 'slot', 'u', 'X', 'k', 'internal')

    def __init__(self, coef, slot, u, X, k, internal=False):
        self.coef, self.slot, self.u, self.X, self.k, self.internal = coef, slot, u, X, tuple(k), internal


def ksum(ms):
    return tuple(sum((m.k[i] for m in ms), Fr(0)) for i in range(4))


def iszerok(k):
    return all(x == 0 for x in k)


def ik(m, mu):  # i k_mu
    return IU * m.k[mu]


def k2(m):
    return sum((m.k[mu] * m.k[mu] for mu in V4), Fr(0))


def kdotV(m, u=None):
    u = m.u if u is None else u
    return sum((m.k[mu] * u[mu] for mu in V4), ZERO)


# ------------------------------------------------------------------ local functionals as terms
# term = (coefficient (Fr), slots tuple, function(list of modes) -> GQ)
def evaluate(F, modes, internal=(), maxint=1):
    """F = list of terms.  modes: external modes (any coefficient elements); internal: internal modes.
    Returns algebra element (monomials with at most maxint internal generators)."""
    out = {}
    byslot = {s: [m for m in modes if m.slot == s] for s in SLOTS}
    ibyslot = {s: [m for m in internal if m.slot == s] for s in SLOTS}
    for coef, slots, fn in F:
        n = len(slots)
        choices = [None] + (list(range(n)) if maxint >= 1 else [])
        for ipos in choices:
            lists = [(ibyslot[slots[j]] if j == ipos else byslot[slots[j]]) for j in range(n)]
            for tup in itertools.product(*lists):
                if not iszerok(ksum(tup)):
                    continue
                c = {(): ONE}
                for m in tup:
                    c = amul(c, m.coef, maxint=maxint)
                    if not c:
                        break
                if not c:
                    continue
                v = fn(tup)
                if v.iszero():
                    continue
                out = aadd(out, ascale(c, v * coef))
    return out


def term(coef, slots, fn):
    return (Fr(coef), tuple(slots), fn)


# ---- densities ----
def dA(m, M, N):
    v = ZERO
    if M in V4:
        v = v + ik(m, M) * m.u[N]
    if N in V4:
        v = v - ik(m, N) * m.u[M]
    return v


def fn_AAA(t):  # sum_{mu,N} (d_mu A_N)(A_mu)(A_N) colour (X1,[X2,X3])
    a, b, c = t
    col = f3(a.X, b.X, c.X)
    if col == 0:
        return ZERO
    return sum((ik(a, mu) * a.u[N] * b.u[mu] * c.u[N] for mu in V4 for N in range(10)), ZERO) * col


def fn_AAAA(t):
    a, b, c, d = t
    col = cp(br(a.X, b.X), br(c.X, d.X))
    if col == 0:
        return ZERO
    s1 = sum((a.u[M] * c.u[M] for M in range(10)), ZERO)
    s2 = sum((b.u[N] * d.u[N] for N in range(10)), ZERO)
    return s1 * s2 * col


def fn_llA(t):  # (lam1, rho([A, lam2]))
    l1, A, l2 = t
    col = f3(l1.X, A.X, l2.X)
    if col == 0:
        return ZERO
    v = ZERO
    for M in range(10):
        if A.u[M].iszero():
            continue
        v = v + A.u[M] * spin_pair(l1.u, matvec(SIG[M], l2.u))
    return v * col


def fn_AcAs(t):
    A, c, As = t
    return dot(A.u, As.u) * c.u[0] * f3(A.X, c.X, As.X)  # ([A,c],A^*) = (A,[c,A^*])


def fn_clls(t):
    c, l, ls = t
    return c.u[0] * spin_pair(l.u, ls.u) * f3(c.X, l.X, ls.X)


def fn_cccs(t):
    a, b, cs = t
    return a.u[0] * b.u[0] * cs.u[0] * f3(a.X, b.X, cs.X)


def fn_kinA(t):  # sum_{M,N} F1_MN F2_MN (X1,X2)
    a, b = t
    return sum((dA(a, M, N) * dA(b, M, N) for M in range(10) for N in range(10)), ZERO) * cp(a.X, b.X)


def fn_kinl(t):
    a, b = t
    return spin_pair(a.u, rho_vec([IU * b.k[mu] for mu in range(4)], b.u)) * cp(a.X, b.X)


def fn_king(t):
    c, As = t
    return sum((ik(c, mu) * As.u[mu] for mu in V4), ZERO) * c.u[0] * cp(c.X, As.X)


I3_terms = {
    'gauge': [term(Fr(-1, 2), ('A', 'A', 'A'), fn_AAA)],
    'llA': [term(Fr(1, 2), ('lam', 'A', 'lam'), fn_llA)],
    'AcAs': [term(1, ('A', 'c', 'As'), fn_AcAs)],
    'clls': [term(1, ('c', 'lam', 'lams'), fn_clls)],
    'cccs': [term(Fr(1, 2), ('c', 'c', 'cs'), fn_cccs)],
}
I3 = sum(I3_terms.values(), [])
I4 = [term(Fr(-1, 8), ('A', 'A', 'A', 'A'), fn_AAAA)]
S2 = [term(Fr(-1, 8), ('A', 'A'), fn_kinA), term(Fr(1, 2), ('lam', 'lam'), fn_kinl), term(1, ('c', 'As'), fn_king)]


def Vpart(u):
    return [u[M] if M in V4 else ZERO for M in range(10)]


def Wpart(u):
    return [u[M] if M in WP else ZERO for M in range(10)]


def dsdA(m):  # (d^*d A|_V)_mu = k^2 u_mu - k_mu (k.u_V)
    kk = k2(m)
    kd = kdotV(m)
    return [m.u[mu] * kk - kd * m.k[mu] for mu in V4]


def ddsA(m):  # (d d^* A|_V)_mu = k_mu (k.u_V)
    kd = kdotV(m)
    return [kd * m.k[mu] for mu in V4]


# ---- target basis G ----
def fn_G1(t):
    c, a1, a2 = t
    w = dsdA(a2)
    return c.u[0] * sum((a1.u[mu] * w[mu] for mu in range(4)), ZERO) * f3(c.X, a1.X, a2.X)


def fn_G2(t):
    c, a1, a2 = t
    w = ddsA(a2)
    return c.u[0] * sum((a1.u[mu] * w[mu] for mu in range(4)), ZERO) * f3(c.X, a1.X, a2.X)


def fn_G3(t):
    c, a1, a2 = t
    return c.u[0] * sum((a1.u[w] * a2.u[w] for w in WP), ZERO) * k2(a2) * f3(c.X, a1.X, a2.X)


def fn_G4(t):  # ([c, dc], A^*)
    c1, c2, As = t
    return c1.u[0] * c2.u[0] * sum((ik(c2, mu) * As.u[mu] for mu in V4), ZERO) * f3(c1.X, c2.X, As.X)


def fn_G5(t):  # ([c, lam], dslash lam)
    c, l1, l2 = t
    return c.u[0] * spin_pair(l1.u, rho_vec([IU * l2.k[mu] for mu in range(4)], l2.u)) * f3(c.X, l1.X, l2.X)


G = {
    'G1': [term(1, ('c', 'A', 'A'), fn_G1)],
    'G2': [term(1, ('c', 'A', 'A'), fn_G2)],
    'G3': [term(1, ('c', 'A', 'A'), fn_G3)],
    'G4': [term(1, ('c', 'c', 'As'), fn_G4)],
    'G5': [term(1, ('c', 'lam', 'lam'), fn_G5)],
}


# ---- H basis ----
def fn_H1(t):  # (dA_V, [A_V ^ A_V]) = 2 sum_{mu nu} (d_mu a_nu, [a_mu, a_nu])
    a, b, c = t
    return 2 * sum((ik(a, mu) * a.u[nu] * b.u[mu] * c.u[nu] for mu in V4 for nu in V4), ZERO) * f3(a.X, b.X, c.X)


def fn_H2(t):  # sum_{mu,w} (d_mu phi_w, [a_mu, phi_w])
    a, b, c = t
    return sum((ik(a, mu) * a.u[w] * b.u[mu] * c.u[w] for mu in V4 for w in WP), ZERO) * f3(a.X, b.X, c.X)


def fn_H3(t):
    A, c, As = t
    return dot(Vpart(A.u), As.u) * c.u[0] * f3(A.X, c.X, As.X)


def fn_H4(t):
    A, c, As = t
    return dot(Wpart(A.u), As.u) * c.u[0] * f3(A.X, c.X, As.X)


H = {
    'H1': [term(1, ('A', 'A', 'A'), fn_H1)],
    'H2': [term(1, ('A', 'A', 'A'), fn_H2)],
    'H3': [term(1, ('A', 'c', 'As'), fn_H3)],
    'H4': [term(1, ('A', 'c', 'As'), fn_H4)],
    'H5': [term(1, ('c', 'lam', 'lams'), fn_clls)],
    'H6': [term(1, ('c', 'c', 'cs'), fn_cccs)],
}


# ---- quadratic functionals ----
def fn_D1(t):  # (d^*A_V, d^*A_V), d^*A = -i k.u_V
    a, b = t
    return (-IU * kdotV(a)) * (-IU * kdotV(b)) * cp(a.X, b.X)


def fn_B1(t):  # (dA_V, dA_V) = 1/2 sum_{mu nu} F F
    a, b = t
    return Fr(1, 2) * sum((dA(a, M, N) * dA(b, M, N) for M in V4 for N in V4), ZERO) * cp(a.X, b.X)


def fn_B2(t):  # (dA_W, dA_W)
    a, b = t
    return sum((ik(a, mu) * a.u[w] * ik(b, mu) * b.u[w] for mu in V4 for w in WP), ZERO) * cp(a.X, b.X)


def fn_mass(t):
    a, b = t
    return dot(a.u, b.u) * cp(a.X, b.X)


def fn_Ob2(t):  # (d^*d c, d^*A_V)
    c, A = t
    return k2(c) * c.u[0] * (-IU * kdotV(A)) * cp(c.X, A.X)


def fn_cntA(t):
    A, As = t
    return dot(A.u, As.u) * cp(A.X, As.X)


def fn_cntl(t):  # <lam^*, lam> = (lam, lam^*), slot order (lams, lam)
    ls, l = t
    return spin_pair(l.u, ls.u) * cp(ls.X, l.X)


def fn_cntc(t):  # <c^*, c>
    cs, c = t
    return cs.u[0] * c.u[0] * cp(cs.X, c.X)


Q2 = {
    'D1': [term(1, ('A', 'A'), fn_D1)],
    'B1': [term(1, ('A', 'A'), fn_B1)],
    'B2': [term(1, ('A', 'A'), fn_B2)],
    'B3': [term(1, ('lam', 'lam'), fn_kinl)],
    'B4': [term(1, ('c', 'As'), fn_king)],
    'mass': [term(1, ('A', 'A'), fn_mass)],
    'Ob2': [term(1, ('c', 'A'), fn_Ob2)],
    'massH': [term(1, ('A', 'A'), lambda t: dot(t[0].u, t[1].u) * k2(t[1]) * cp(t[0].X, t[1].X))],
    'cntA': [term(1, ('A', 'As'), fn_cntA)],
    'cntl': [term(1, ('lams', 'lam'), fn_cntl)],
    'cntc': [term(1, ('cs', 'c'), fn_cntc)],
}


# ------------------------------------------------------------------ Q, bracket, propagator
def Qimage(m, r):
    """mode of Q(m) with coefficient r*kappa, or None."""
    coef = amul(gen_elem(r), m.coef)
    s = m.slot
    if s == 'c':
        u = [-ik(m, mu) * m.u[0] if mu in V4 else ZERO for mu in range(10)]
        return Mode(coef, 'A', u, m.X, m.k, m.internal)
    if s == 'A':
        w = dsdA(m)
        kk = k2(m)
        u = [Fr(-1, 2) * w[M] if M in V4 else Fr(-1, 2) * kk * m.u[M] for M in range(10)]
        return Mode(coef, 'As', u, m.X, m.k, m.internal)
    if s == 'As':
        return Mode(coef, 'cs', [-IU * kdotV(m)], m.X, m.k, m.internal)
    if s == 'lam':
        return Mode(coef, 'lams', rho_vec([IU * m.k[mu] for mu in range(4)], m.u), m.X, m.k, m.internal)
    return None


def applyQ(F, modes, internal=(), maxint=1):
    r = new_gen(1, 'r')
    extra = [q for q in (Qimage(m, r) for m in modes) if q is not None]
    val = evaluate(F, list(modes) + extra, internal, maxint)
    return dleft(val, r)


def make_internal(kset):
    """internal basis modes for every slot, fibre basis vector, colour basis vector and momentum in kset."""
    out = []
    idx = {}
    for s in SLOTS:
        for i in range(SDIM[s]):
            for a in range(DIMG):
                for k in kset:
                    g = new_gen(SPAR[s], 'int')
                    u = [ONE if j == i else ZERO for j in range(SDIM[s])]
                    X = [Fr(1) if b == a else Fr(0) for b in range(DIMG)]
                    m = Mode(gen_elem(g), s, u, X, k, internal=True)
                    out.append(m)
                    idx[(s, i, a, tuple(k))] = g
    return out, idx


def neg(k):
    return tuple(-x for x in k)


def derivs(val):
    """map internal generator -> left derivative (an element without internal generators); single pass over the
    monomials with exactly one internal generator."""
    out = {}
    for m, c in val.items():
        ints = [g for g in m if KIND[g] == 'int']
        if len(ints) != 1:
            continue
        g = ints[0]
        j = m.index(g)
        sgn = -1 if (PAR[g] and sum(PAR[h] for h in m[:j]) % 2) else 1
        mm = m[:j] + m[j + 1:]
        d = out.setdefault(g, {})
        d[mm] = d.get(mm, ZERO) + (c if sgn == 1 else -c)
    return {g: {m: c for m, c in d.items() if not c.iszero()} for g, d in out.items()}


def bracket(f, pf, g, modes, kset, internal_cache):
    """{f, g} on external modes; pf = parity of f."""
    internal, idx = internal_cache
    df = derivs(evaluate(f, modes, internal, 1))
    dg = derivs(evaluate(g, modes, internal, 1))
    out = {}
    for x, xi in XMEMBER.items():
        for i in range(SDIM[x]):
            for a in range(DIMG):
                for k in kset:
                    gx = idx[(x, i, a, tuple(k))]
                    gxi = idx[(xi, i, a, neg(k))]
                    if gxi in df and gx in dg:
                        out = aadd(out, amul(df[gxi], dg[gx]))
                    if gx in df and gxi in dg:
                        t = amul(df[gx], dg[gxi])
                        out = aadd(out, t if pf == 0 else ascale(t, -1))
    return out


def QGF_numeric(slot, i, k):
    """Q^GF applied to the basis vector e_i of slot at momentum k: returns (target slot, numeric vector) or None."""
    if slot == 'As':
        return 'A', [GQ(-2) if j == i else ZERO for j in range(10)]
    if slot == 'A':
        if i in V4:
            return 'c', [IU * k[i]]          # -d^*(A_V) = i k.u
        return None
    if slot == 'cs':
        return 'As', [IU * k[j] if j in V4 else ZERO for j in range(10)]
    if slot == 'lams':
        e = [ONE if j == i else ZERO for j in range(16)]
        w = rho_vec([IU * k[mu] for mu in range(4)], e, bar=True)
        return 'lam', [-x for x in w]
    return None


def dP_terms(kset, idx):
    """list of (coef, gen_a, gen_b): d_P = sum coef * d_{gen_a} d_{gen_b} (gen_b applied first), P = P(0,inf)."""
    out = []
    for x, xi in XMEMBER.items():
        for i in range(SDIM[x]):
            for a in range(DIMG):
                for k in kset:
                    kk = sum((k[mu] * k[mu] for mu in V4), Fr(0))
                    # K_0 = - sum (b^x_k (x) b^xi_{-k} + b^xi_{-k} (x) b^x_k);  P = -(1/k^2)[Q^GF b^x_k (x) b^xi_{-k} + Q^GF b^xi_{-k} (x) b^x_k]
                    for (s1, k1, s2, k2_) in ((x, k, xi, neg(k)), (xi, neg(k), x, k)):
                        img = QGF_numeric(s1, i, k1)
                        if img is None:
                            continue
                        ts, vec = img
                        g2 = idx[(s2, i, a, tuple(k2_))]
                        for j, cj in enumerate(vec):
                            if cj.iszero():
                                continue
                            g1 = idx[(ts, j, a, tuple(k1))]
                            # d_P = 1/2 sum d_{phi''} d_{phi'}, phi' = -(1/k^2) cj b_{g1}, phi'' = b_{g2}
                            out.append((cj * Fr(-1, 2) / kk, g2, g1))
    return out


def tree(f, pf, g, modes, kset, internal_cache, dPt):
    """cross part of d_P(f g), P = P(0,infinity): sum coef * d_a d_b (f g) with one derivative on each factor,
    d_b(fg) = (d_b f) g + (-1)^{|b||f|} f d_b g, hence the cross part of d_a d_b (fg) is
    (-1)^{|a|(|f|+|b|)} (d_b f)(d_a g) + (-1)^{|b||f|} (d_a f)(d_b g)."""
    internal, idx = internal_cache
    df = derivs(evaluate(f, modes, internal, 1))
    dg = derivs(evaluate(g, modes, internal, 1))
    out = {}
    for coef, ga, gb in dPt:
        if gb in df and ga in dg:
            s1 = -1 if (PAR[ga] * (pf + PAR[gb])) % 2 else 1
            out = aadd(out, ascale(amul(df[gb], dg[ga]), coef * s1))
        if ga in df and gb in dg:
            s2 = -1 if (PAR[gb] * pf) % 2 else 1
            out = aadd(out, ascale(amul(df[ga], dg[gb]), coef * s2))
    return out


def dP_apply(val, dPt):
    out = {}
    for coef, ga, gb in dPt:
        t = dleft(dleft(val, gb), ga)
        if t:
            out = aadd(out, ascale(t, coef))
    return out


def ext_coeff(val, gens):
    m = tuple(sorted(gens))
    return val.get(m, ZERO)
