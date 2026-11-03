"""Engine A, part Q, task TQ3: the pair entries of the holomorphic twist along a rank-one supercharge xi = u (x) z.

Weight basis (basis-free construction).  Cartan generators Sigma^{01}, Sigma^{23} (so(4)) and Sigma^{45}, Sigma^{67},
Sigma^{89} (so(6)) act on Sigma with eigenvalues +-i/2; the 16 joint eigenvectors are the weight vectors (exact, Gaussian
integers up to normalisation).  S_+ (x) W = {s1 = -s2, t1 t2 t3 = -1}, S_- (x) W^* = {s1 = s2, t1 t2 t3 = +1}.
  xi = weight vector (u; z) with u = (+,-), z = (+,+,-): Gamma(xi, xi) = 0 (weight argument; checked);
  u' = (-,+) (the other S_+ weight);  vectors of V (x) C: null n_w with so(4) weights (+-1,0), (0,+-1); the directions
  u (x) S_-  (weights (1,0), (0,-1)) are Q^xi_0-exact derivatives, u' (x) S_- (weights (0,1), (-1,0)) are d_alphadot.
Twisted letters (Rem c11-twist, the dictionary of Rem qs-twistpairs, index-free):
  gamma : iota_w A with rho(w) xi = 0 (w in W_R^perp (x) C, 3-dim);
  beta  : coefficient of lam_+ along the weight vectors (u; tau), tau != z (3);
  dc    : coefficient of lam_- along (sbar; -z), sbar = (+,+), (-,-) (2) [lambdabar_{alphadot i} z^i];
  b     : u'u' component of the self-dual field strength, F(d_1, d_2), d_1, d_2 the two d_alphadot null vectors;
          or the antifield representative xi . lam^*_+ (= (xi, lam^*) = u_alpha lam^{*alpha}_4).
Outputs are evaluated on the weight vectors (dual to the letters); output momenta are expanded in the four null
coordinates c_w(q) = (n_{-w}.q)/(n_w.n_{-w}); c_w multiplies the derivative along n_{-w}, so the holomorphic
coordinates (labels h0 = derivative along d_1 = n_(0,1), h1 = along d_2 = n_(-1,0)) are c_(0,-1), c_(1,0), the exact
ones (labels e0 = along n_(1,0), e1 = along n_(0,-1)) are c_(-1,0), c_(0,1)."""
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import qa_s1 as S1, qa_core as QC, qa_util as U

SB = S1.SB
CART = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]


def _sym(M):
    return sp.Matrix(M.shape[0], M.shape[1], lambda i, j: sp.nsimplify(round(M[i, j].real, 8)) + sp.I * sp.nsimplify(round(M[i, j].imag, 8)))


def spinor_weights():
    T = [_sym(S1.spin(a, b)[0]) for a, b in CART]
    out = {}
    for signs in itertools.product([1, -1], repeat=5):
        M = sp.Matrix.vstack(*[T[k] - signs[k] * sp.I / 2 * sp.eye(16) for k in range(5)])
        ns = M.nullspace()
        if not ns:
            continue
        v = ns[0]
        piv = [x for x in v if x != 0][0]
        v = v / piv
        den = sp.ilcm(*[sp.fraction(sp.nsimplify(sp.re(x)))[1] for x in v] + [sp.fraction(sp.nsimplify(sp.im(x)))[1] for x in v])
        v = v * den
        out[signs] = np.array([complex(sp.re(x), sp.im(x)) for x in v])
    assert len(out) == 16
    return out


def vector_null(k, s):
    """null vector of weight s (= +-1) in the Cartan plane k: e_a + i s e_b  (eigenvalue i s of Sigma^{ab})"""
    a, b = CART[k]
    v = np.zeros(10, complex); v[a] = 1; v[b] = 1j * s
    # check: (Sigma v)_m = delta_{m a} v_b - delta_{m b} v_a
    Sv = np.zeros(10, complex); Sv[a] += v[b]; Sv[b] -= v[a]
    assert np.allclose(Sv, 1j * s * v)
    return v


class Twist:
    def __init__(self):
        W = spinor_weights()
        self.W = W
        self.u, self.up = (1, -1), (-1, 1)
        self.z = (1, 1, -1)
        self.xi = W[self.u + self.z]
        for M in range(10):
            assert abs(self.xi @ SB[M] @ self.xi) < 1e-12
        Wp = [t for t in itertools.product([1, -1], repeat=3) if t[0] * t[1] * t[2] == -1]
        self.taus = [t for t in Wp if t != self.z]                      # W weights != z
        # basis of Sigma: weight vectors; dual functionals (rows of the inverse)
        keys = sorted(W)
        Bm = np.column_stack([W[k] for k in keys])
        Binv = np.linalg.inv(Bm)
        self.dual = {k: Binv[i] for i, k in enumerate(keys)}
        # gamma letters: w in W_perp (x) C with rho(w) xi = 0: null vectors in the perp planes with weight tau + z
        self.gam_w, self.gam_out = [], []
        for t in self.taus:
            wt = tuple(a + b for a, b in zip(t, self.z))               # entries in {0, +-2}/2 -> one nonzero +-1 (x2/2)
            k = [i for i, x in enumerate(wt) if x != 0]
            assert len(k) == 1
            s = wt[k[0]] // 2
            w10 = vector_null(2 + k[0], s)
            w6 = w10[4:]
            assert np.allclose(sum(w6[i] * (SB[4 + i] @ self.xi) for i in range(6)), 0)
            self.gam_w.append(w6)
            wbar = vector_null(2 + k[0], -s)[4:]
            self.gam_out.append(wbar / (w6 @ wbar))                   # w.out = 1: the output leg of the letter iota_w A
        # beta letters (u; tau), tau != z: functional = dual row; output legs = the weight vectors
        self.beta_z = [self.dual[self.u + t] for t in self.taus]
        self.beta_out = [W[self.u + t] for t in self.taus]
        # dc letters: lam_- weight (sbar; -z), sbar in {(1,1), (-1,-1)}
        mz = tuple(-x for x in self.z)
        self.dc_keys = [(1, 1) + mz, (-1, -1) + mz]
        self.dc_z = [self.dual[k] for k in self.dc_keys]
        self.dc_out = [W[k] for k in self.dc_keys]
        # vectors: so(4) null directions; holomorphic (u') d_1 = weight (0,1), d_2 = (-1,0); exact (1,0), (0,-1)
        self.nv = {}
        for k in (0, 1):
            for s in (1, -1):
                w = [0, 0]; w[k] = s
                self.nv[tuple(w)] = vector_null(k, s)[:4]
        self.hol = [(0, 1), (-1, 0)]
        self.exa = [(1, 0), (0, -1)]
        # b letter (self-dual u'u' field strength) = F(d1, d2)
        self.d1, self.d2 = self.nv[(0, 1)], self.nv[(-1, 0)]

    def coords(self, q):
        """null coordinates c_w(q) with q = sum_w c_w n_w"""
        out = {}
        for w, n in self.nv.items():
            nb = self.nv[(-w[0], -w[1])]
            out[w] = complex(np.dot(q, nb) / np.dot(n, nb))
        return out


def b_letters(T):
    """the b letter F(d1,d2) = sum_nu (sum_mu c_{mu nu} d_mu) A_nu, c_{mu nu} = d1^mu d2^nu - d2^mu d1^nu, grouped by nu:
    list of (1, zeta = e_nu, deriv = sum_mu c_{mu nu} (i p_mu))"""
    terms = []
    for nu in range(4):
        f = {}
        for mu in range(4):
            c = T.d1[mu] * T.d2[nu] - T.d2[mu] * T.d1[nu]
            if abs(c) < 1e-12:
                continue
            e = [0, 0, 0, 0]; e[mu] = 1
            f[tuple(e)] = 1j * c
        if f:
            terms.append((1, np.eye(4)[nu].astype(complex), f))
    return terms


def multi_value(La, Lb, legs, gen, what=('fin',)):
    """letters given as lists of (coef, zeta, deriv)"""
    tot = {w: U.vzero() for w in what}
    for ca, za, fa, ta in La:
        for cb, zb, fb, tb in Lb:
            v = QC.pair_value((ta, za, fa), (tb, zb, fb), legs, gen, what=what)
            c = complex(ca) * complex(cb)
            for w in what:
                tot[w] = U.vadd(tot[w], U.vscale(v[w], U.gq(c)) if abs(c - round(c.real) - 1j * round(c.imag)) < 1e-9
                                else vscale_c(v[w], c))
    return tot


def vscale_c(v, c):
    fr = Fraction(c.real).limit_denominator(10 ** 6), Fraction(c.imag).limit_denominator(10 ** 6)
    return U.vscale(v, fr)


def interpolate(fun, deg, nvar, rng=None):
    """exact interpolation of a polynomial value-function of nvar coordinates of degree <= deg (<= 2) on the principal
    lattice {0, e_i, 2 e_i, e_i + e_j}: returns {monomial exponent tuple: value}"""
    assert deg <= 2
    z = [0] * nvar
    def at(*idx):
        x = [0] * nvar
        for i in idx:
            x[i] += 1
        return fun([complex(v) for v in x])
    out = {}
    f0 = at()
    if not U.viszero(f0):
        out[tuple(z)] = f0
    if deg == 0:
        return out
    half = (Fraction(1, 2), Fraction(0))
    f1 = [at(i) for i in range(nvar)]
    cii = [U.vzero()] * nvar
    if deg == 2:
        f2 = [at(i, i) for i in range(nvar)]
        cii = [U.vscale(U.vadd(U.vsub(f2[i], U.vscale(f1[i], (Fraction(2), Fraction(0)))), f0), half) for i in range(nvar)]
    bi = [U.vsub(U.vsub(f1[i], f0), cii[i]) for i in range(nvar)]
    for i in range(nvar):
        e = list(z); e[i] = 1
        if not U.viszero(bi[i]):
            out[tuple(e)] = bi[i]
        if deg == 2 and not U.viszero(cii[i]):
            e2 = list(z); e2[i] = 2
            out[tuple(e2)] = cii[i]
    if deg == 2:
        for i in range(nvar):
            for j in range(i + 1, nvar):
                v = at(i, j)
                for t in (f0, bi[i], bi[j], cii[i], cii[j]):
                    v = U.vsub(v, t)
                if not U.viszero(v):
                    e = list(z); e[i] = 1; e[j] = 1
                    out[tuple(e)] = v
    return out


# ------------------------------------------------------------------------------------------------ entries
def entry(T, gx, La, Lb, out_types, out_vecs, D, rng, log=print, label=''):
    """kernel on twisted inputs (La, Lb: lists of (coef, zeta, deriv, type)) and twisted outputs: for every choice of
    output vectors (out_vecs: list per output of [(name, vector)]) and both orderings, the exact polynomial in the null
    coordinates of the two output momenta (variables c1_w, c2_w, w in hol + exact).  Returns rows."""
    wl = T.hol + T.exa
    rows = []
    for choice in itertools.product(*out_vecs):
        for order in ((0, 1), (1, 0)):
            def fun(x, choice=choice, order=order):
                q = []
                for j in range(2):
                    qq = np.zeros(4, complex)
                    for a, w in enumerate(wl):
                        qq = qq + x[4 * j + a] * T.nv[w]
                    q.append(qq)
                legs = [(out_types[j], choice[j][1], q[j]) for j in range(2)]
                legs = [legs[i] for i in order]
                return multi_value(La, Lb, legs, gx)['fin']
            poly = interpolate(fun, D, 8, rng)
            # the coordinate c_w(q) multiplies the derivative along n_{-w}: holomorphic iff -w is a u'-type weight
            def lab(w):
                mw = (-w[0], -w[1])
                return ('h%d' % T.hol.index(mw)) if mw in T.hol else ('e%d' % T.exa.index(mw))
            names = ['c%d%s' % (j + 1, lab(w)) for j in range(2) for w in wl]
            terms = {'*'.join('%s^%d' % (names[k], e[k]) if e[k] > 1 else names[k] for k in range(8) if e[k]) or '1': U.vstr(v)
                     for e, v in poly.items()}
            hol_idx = [k for k in range(8) if names[k][2] == 'h']
            exact_free = all(all(e[k] == 0 for k in range(8) if k not in hol_idx) for e in poly)
            rows.append({'entry': label, 'outputs': [c[0] for c in choice], 'order': [out_types[i] for i in order],
                         'poly': terms, 'free_of_exact_directions': exact_free})
            log('    %s %s order %s: %s' % (label, [c[0] for c in choice], order, terms))
    return rows


def letters(T, which, idx=None):
    if which == 'beta':
        return [(1, T.beta_z[idx], None, 'lamP')]
    if which == 'gamma':
        return [(1, T.gam_w[idx], None, 'AW')]
    if which == 'b':
        return [(c, z, f, 'AV') for (c, z, f) in b_letters(T)]
    if which == 'c':
        return [(1, np.array([1 + 0j]), None, 'c')]
    if which == 'bstar':
        return [(1, T.xi, None, 'lamsP')]
    raise ValueError(which)


def outs(T, which):
    if which == 'dc':
        return [('dc%d' % a, v) for a, v in enumerate(T.dc_out)]
    if which == 'gamma':
        return [('g%d' % a, v) for a, v in enumerate(T.gam_out)]
    if which == 'beta':
        return [('b%d' % a, v) for a, v in enumerate(T.beta_out)]
    raise ValueError(which)


OUT_TYPE = {'dc': 'lamM', 'gamma': 'AW', 'beta': 'lamP'}


def run(log=print, seed=3, which=('betagamma', 'betabeta', 'bstar_c', 'bgamma', 'betab_gg', 'betab_bdc', 'gammagamma')):
    T = Twist()
    gx = QC.odd_gen(T.xi)
    rng = np.random.default_rng(seed)
    res = {'xi_weight': [T.u, T.z], 'rows': {}}
    spec = {
        'betagamma': (('beta', 0), ('gamma', 0), ('dc', 'dc'), 1),
        'betabeta': (('beta', 0), ('beta', 1), ('gamma', 'dc'), 1),
        'bstar_c': (('c', None), ('bstar', None), ('dc', 'dc'), 0),
        'bgamma': (('b', None), ('gamma', 0), ('gamma', 'dc'), 1),
        'betab_gg': (('beta', 0), ('b', None), ('gamma', 'gamma'), 2),
        'betab_bdc': (('beta', 0), ('b', None), ('beta', 'dc'), 1),
        'gammagamma': (('gamma', 0), ('gamma', 1), ('gamma', 'beta'), 1),
    }
    for nm in which:
        (wa, ia), (wb, ib), (o1, o2), D = spec[nm]
        La, Lb = letters(T, wa, ia), letters(T, wb, ib)
        log('  TQ3 entry %s: letters (%s%s, %s%s), outputs (%s, %s), degree %d' % (nm, wa, '' if ia is None else ia, wb, '' if ib is None else ib, o1, o2, D))
        res['rows'][nm] = entry(T, gx, La, Lb, (OUT_TYPE[o1], OUT_TYPE[o2]), [outs(T, o1), outs(T, o2)], D, rng, log=log, label=nm)
    return res, T
