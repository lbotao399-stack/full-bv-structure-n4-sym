"""Finite transfer endpoint T: one-loop graph engine (independent numerics for the K-route).

Imports the Step-1 core of engine A (engine_a/r2_core.py: field layout c, A, A^*, lambda, lambda^*, c^*; su(2) in an
orthonormal basis of (-,-); exact cubic/quartic vertex functions scaled by 2 / 8; propagator numerators Pnum(k, line),
special-edge kernel numerators Khat(line); parity diagonal SDIAG).  Nothing from step2 engine outputs is read.

Kernels at finite T (momentum space, scalar factors):
   P(eps,T)(k)  = int_eps^T dt e^{-t k^2} Pnum(k)                 (edges 'AV','AP','gh','l')
   (K_eps-K_0)(k) = (e^{-eps k^2}-1) Khat = -int_0^eps ds k^2 e^{-s k^2} Khat   (edges 'sAV','sAP','sl','sc')
   T d/dT P(eps,T)(k) = T e^{-T k^2} Pnum(k)                         (edge pinned at t = T)
One-loop formulas (verified in 0 dimensions, check0d_special_triangle.py / check0d_quartic_bubble.py):
   flow, 1PI:          sum_n 1/(2n) sum_{orderings} str[ P S h_1 P S h_2 ... ]
   special edge, 1PI:  1/2 sum_{orderings} str[ (K_eps-K_0) S h_1 P S h_2 ... ]
(S = parity, h = Hessian with left derivatives, odd parameters factored to the left; sign of the odd parameters in
cycle order relative to the canonical order (b1, b2, gamma) as in engine A).

Insertion of S^(1)(v) for a plane-wave vector field v = u e^{iP.y}: Hessian built from ch02 eq:c2-S1even,
   S^(1)(v) = int ( L_v A + 1/4 div v A|_W , A^* ) - ( v.lambda, lambda^* ) - ( v.c, c^* ),
   v.lambda = v(lambda) + 1/4 rho(d(v,-)) lambda + 3/8 div v lambda,  rho(q(a^b)) = rho(a)rho(b) - (a,b),
on plane waves: d -> i k, d(v^n) -> i P u^n; the functional is S^(1)(v) gamma with gamma the odd ghost of the
generator (gamma placed last).  K_v (v = |y|^2 b - 2 (b.y) y) is reached through y^a y^b = -d_{P_a} d_{P_b} e^{iPy}|_{P=0}:
the P-jet to second order is computed exactly with hyper-dual numbers (P = e1 d1 + e2 d2, e1^2 = e2^2 = 0).
Units: hbar/(16 pi^2) (the factor of the one-loop term and (16 pi^2)^{-1} of the Gaussian integral stripped);
divide by C_ADJ = tr_g(ad_X ad_Y)/(X,Y) = -2 for the unit hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))."""
import sys, os, itertools, math, warnings
import numpy as np
np.seterr(all="ignore"); warnings.filterwarnings("ignore")
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'engine_a'))
import r2_core as C                         # noqa: E402

D = C.D
NT = C.NTOT
SD = C.SDIAG.real.copy()
SARR = C.SARR; SBARR = C.SBARR
E4 = np.eye(4)
PLINES = ['AV', 'AP', 'gh', 'l']
SLINES = ['sAV', 'sAP', 'sl', 'sc']
PARITY = dict(C.PARITY); PARITY['g'] = 1      # gamma: odd ghost of the even generator


# ------------------------------------------------------------------ hyper-dual numbers / matrices
class HD:
    """a + b e1 + c e2 + d e1 e2 ; components are numpy arrays or None (= 0)"""
    __slots__ = ('a', 'b', 'c', 'd')

    def __init__(self, a, b=None, c=None, d=None):
        self.a, self.b, self.c, self.d = a, b, c, d

    @staticmethod
    def _mul(x, y, op):
        if x is None or y is None:
            return None
        return op(x, y)

    @staticmethod
    def _add(*xs):
        xs = [x for x in xs if x is not None]
        if not xs:
            return None
        out = xs[0]
        for x in xs[1:]:
            out = out + x
        return out

    def mat(self, o):
        m = lambda x, y: HD._mul(x, y, np.matmul)
        return HD(self.a @ o.a,
                  HD._add(m(self.a, o.b), m(self.b, o.a)),
                  HD._add(m(self.a, o.c), m(self.c, o.a)),
                  HD._add(m(self.a, o.d), m(self.b, o.c), m(self.c, o.b), m(self.d, o.a)))

    def scale_rows(self, v):
        f = lambda x: None if x is None else v[:, None] * x
        return HD(f(self.a), f(self.b), f(self.c), f(self.d))

    def comps(self):
        z = lambda x: 0.0 if x is None else x
        return [z(self.a), z(self.b), z(self.c), z(self.d)]


def hd_vec(re, d1=None, d2=None):
    return (np.asarray(re, dtype=float), None if d1 is None else np.asarray(d1, float), None if d2 is None else np.asarray(d2, float))


def hd_vec_add(x, y, sy=1.0):
    f = lambda a, b: a if b is None else (sy * b if a is None else a + sy * b)
    return (x[0] + sy * y[0], f(x[1], y[1]), f(x[2], y[2]))


# ------------------------------------------------------------------ line blocks
def _line_support(L):
    kk = np.array([0.31, -0.77, 1.13, 0.29])
    N = C.Pnum(kk, L) if L in PLINES else C.Khat(L)
    rows = np.nonzero(np.abs(N).sum(axis=1))[0]
    cols = np.nonzero(np.abs(N).sum(axis=0))[0]
    return rows, cols


SUPPORT = {L: _line_support(L) for L in PLINES + SLINES}


class Line:
    """numerator of an edge restricted to its support: affine in the edge momentum"""
    def __init__(self, L):
        self.L = L
        self.rows, self.cols = SUPPORT[L]
        ix = np.ix_(self.rows, self.cols)
        if L in PLINES:
            z = np.zeros(4)
            self.N0 = C.Pnum(z, L)[ix]
            self.Nm = np.array([C.Pnum(E4[m], L)[ix] - self.N0 for m in range(4)])
        else:
            self.N0 = C.Khat(L)[ix].real.astype(complex)
            self.Nm = np.zeros((4,) + self.N0.shape, complex)
        self.special = L in SLINES

    def at(self, k):   # k hyper-dual vector
        a = self.N0 + np.tensordot(k[0], self.Nm, 1)
        b = None if k[1] is None or not np.any(self.Nm) else np.tensordot(k[1], self.Nm, 1)
        c = None if k[2] is None or not np.any(self.Nm) else np.tensordot(k[2], self.Nm, 1)
        return HD(a, b, c, None)


LINES = {L: Line(L) for L in PLINES + SLINES}


# ------------------------------------------------------------------ insertions
class Insertion:
    """Hessian of an insertion, affine in (kj, kl, p): h = h0 + kj.Hj + kl.Hl + p.Hp (full NT x NT, unscaled).
    kind: 'bg' (one background leg of I: cubic vertices), 'bg2' (two background A legs: quartic vertex),
    'S1' (S^(1)(u e^{iPy}), p = P), 'zero'."""
    def __init__(self, kind, **kw):
        self.kind = kind
        z = np.zeros(4)
        if kind == 'bg':
            t, pol = kw['t'], kw['pol']
            f = lambda kj, kl, p: C.hessian([(t, pol, p)], kj, kl) / 2.0
        elif kind == 'bg2':
            (t1, pol1), (t2, pol2) = kw['legs']
            f = lambda kj, kl, p: C.hessian([(t1, pol1, z), (t2, pol2, z)], kj, kl) / 8.0
        elif kind == 'S1':
            u = np.asarray(kw['u'], float)
            f = lambda kj, kl, p: s1_hessian(u, kj, kl, p)
        else:
            raise ValueError(kind)
        self.h0 = f(z, z, z)
        self.Hj = np.array([f(E4[m], z, z) - self.h0 for m in range(4)])
        self.Hl = np.array([f(z, E4[m], z) - self.h0 for m in range(4)])
        self.Hp = np.array([f(z, z, E4[m]) - self.h0 for m in range(4)])
        # affinity check at a random point
        kj, kl, p = np.array([.3, -1.2, .7, 2.]), np.array([-.4, .9, 1.1, -.6]), np.array([1.3, .2, -.8, .5])
        lin = self.h0 + np.tensordot(kj, self.Hj, 1) + np.tensordot(kl, self.Hl, 1) + np.tensordot(p, self.Hp, 1)
        assert np.allclose(lin, f(kj, kl, p)), 'insertion not affine'
        self.cache = {}

    def block(self, rows, cols):
        key = (tuple(rows), tuple(cols))
        if key not in self.cache:
            ix = np.ix_(rows, cols)
            self.cache[key] = (self.h0[ix], self.Hj[(slice(None),) + ix], self.Hl[(slice(None),) + ix], self.Hp[(slice(None),) + ix])
        return self.cache[key]

    def at(self, rows, cols, kj, kl, p):
        h0, Hj, Hl, Hp = self.block(rows, cols)
        a = h0 + np.tensordot(kj[0], Hj, 1) + np.tensordot(kl[0], Hl, 1) + np.tensordot(p[0], Hp, 1)
        def eps_part(i):
            terms = []
            if kj[i] is not None: terms.append(np.tensordot(kj[i], Hj, 1))
            if kl[i] is not None: terms.append(np.tensordot(kl[i], Hl, 1))
            if p[i] is not None: terms.append(np.tensordot(p[i], Hp, 1))
            return HD._add(*terms)
        return HD(a, eps_part(1), eps_part(2), None)


def s1_hessian(u, kj, kl, P):
    """H_jl = (-1)^{|j||l|} * coefficient of the odd parameters (j, l, gamma) in S^(1)(u e^{iPy}) gamma, legs j at kj,
    l at kl (incoming), P incoming at the vertex.  Only field-antifield blocks."""
    H = np.zeros((NT, NT), complex)
    uP = float(np.dot(u, P[:4]))
    sP = np.einsum('m,mab->ab', np.asarray(P[:4], complex), SARR[:4])
    sbu = np.einsum('m,mab->ab', np.asarray(u, complex), SBARR[:4])
    Ml = lambda kla: (1j * np.dot(u, kla) + 1j / 8 * uP) * np.eye(16) + 1j / 4 * (sP @ sbu)   # v.lambda = Ml lambda
    for a in range(D):
        for M in range(10):
            for N in range(10):
                if M < 4 and N < 4:
                    vj = 1j * np.dot(u, kj) * (M == N) + u[M] * 1j * P[N]        # A=e_M at kj, A^*=e_N
                    vl = 1j * np.dot(u, kl) * (M == N) + u[M] * 1j * P[N]        # A=e_M at kl
                elif M >= 4 and N >= 4 and M == N:
                    vj = 1j * np.dot(u, kj) + 0.25j * uP
                    vl = 1j * np.dot(u, kl) + 0.25j * uP
                else:
                    continue
                iA = C.OFF['A'] + M * D + a; iS = C.OFF['As'] + N * D + a
                H[iA, iS] += vj          # j = A (even), l = A^*: sign (+)
                H[iS, iA] += vl          # j = A^*, l = A (at kl)
        Mj, Mlk = Ml(kj), Ml(kl)
        for al in range(16):
            for be in range(16):
                il = C.OFF['l'] + al * D + a; isl = C.OFF['ls'] + be * D + a
                H[il, isl] += -Mj[be, al]        # j = lambda at kj, l = lambda^*: -(v.lambda)_be lambda^*_be
                H[isl, il] += -Mlk[be, al]       # j = lambda^*, l = lambda at kl
        ic = C.OFF['c'] + a; ics = C.OFF['cs'] + a
        H[ic, ics] += -1j * np.dot(u, kj)
        H[ics, ic] += -1j * np.dot(u, kl)
    return H


# ------------------------------------------------------------------ cycles
def perm_sign(seq):
    s_ = 1; a = list(seq)
    for i in range(len(a)):
        for j in range(len(a) - 1 - i):
            if a[j] > a[j + 1]:
                a[j], a[j + 1] = a[j + 1], a[j]; s_ = -s_
    return s_


class Cycle:
    """one cyclic sequence of insertions with line types; line i enters insertion i (j-leg), leaves insertion i-1.
    ins: list of (Insertion, momentum_shift_HD, odd_label or None); momentum shift = incoming momentum at the insertion.
    Line momenta k_0 = k, k_{i+1} = k_i - Pin_i  (all incoming)."""
    def __init__(self, ins, lines, pref):
        self.ins = ins; self.lines = [LINES[L] for L in lines]; self.Lnames = list(lines)
        self.n = len(ins); self.pref = pref
        odd = [lab for (_, _, lab) in ins if lab is not None]
        self.sign = perm_sign(odd)
        # shifts q_i = a_i + beta_i P  (P hyper-dual direction coefficients stored as scalars beta_i via d1/d2 parts)
        self.q = [hd_vec(np.zeros(4))]
        for i in range(self.n - 1):
            self.q.append(hd_vec_add(self.q[-1], ins[i][1], -1.0))
        self.special = self.lines[0].special
        assert all(not L.special for L in self.lines[1:])

    def value(self, k):
        n = self.n
        ks = [hd_vec_add(hd_vec(k), q) for q in self.q]
        X = None
        for i in range(n):
            Li, Lnext = self.lines[i], self.lines[(i + 1) % n]
            phi = Li.at(ks[i])
            neg = (-ks[i][0], None if ks[i][1] is None else -ks[i][1], None if ks[i][2] is None else -ks[i][2])
            h = self.ins[i][0].at(Li.cols, Lnext.rows, neg, ks[(i + 1) % n], self.ins[i][1])
            Mi = phi.mat(h.scale_rows(SD[Li.cols]))
            X = Mi if X is None else X.mat(Mi)
        r0 = self.lines[0].rows
        out = []
        for comp in X.comps():
            if isinstance(comp, float):
                out.append(0.0)
            else:
                out.append(np.sum(SD[r0] * np.diag(comp)))
        return np.array(out, complex) * (self.sign * self.pref)

    def nonzero(self, trials=2):
        rng = np.random.default_rng(5)
        for _ in range(trials):
            if np.abs(self.value(rng.normal(size=4))).max() > 1e-12:
                return True
        return False


# ------------------------------------------------------------------ polynomials in k (sparse monomial basis)
def stirling1_matrix(Dm):
    """B[a, r] = coefficient of x^r in binom(x, a)"""
    B = np.zeros((Dm + 1, Dm + 1))
    for a in range(Dm + 1):
        poly = np.array([1.0])
        for j in range(a):
            poly = np.convolve(poly, np.array([-j, 1.0]))
        poly = poly / math.factorial(a)
        B[a, :len(poly)] = poly
    return B


def simplex(Dm):
    return [al for al in itertools.product(range(Dm + 1), repeat=4) if sum(al) <= Dm]


class Poly:
    """sum_c coef[c] k^{exps[c]} ; coef may carry a trailing axis (several polynomials at once)"""
    def __init__(self, exps, coef):
        self.exps = np.asarray(exps, int).reshape(-1, 4); self.coef = np.asarray(coef, complex)

    @staticmethod
    def fit(f, Dm, ncomp):
        pts = simplex(Dm)
        vals = {al: f(np.array(al, float)) for al in pts}
        # forward differences
        dif = {}
        for al in pts:
            s_ = np.zeros(ncomp, complex)
            for be in itertools.product(*[range(a + 1) for a in al]):
                c = (-1) ** (sum(al) - sum(be)) * np.prod([math.comb(al[j], be[j]) for j in range(4)])
                s_ = s_ + c * vals[be]
            dif[al] = s_
        B = stirling1_matrix(Dm)
        coef = {}
        for al, dv in dif.items():
            if np.abs(dv).max() == 0:
                continue
            for r in itertools.product(*[range(a + 1) for a in al]):
                c = np.prod([B[al[j], r[j]] for j in range(4)])
                if c != 0:
                    coef[r] = coef.get(r, 0) + c * dv
        exps = list(coef.keys()) or [(0, 0, 0, 0)]
        P_ = Poly(exps, np.array([coef.get(e, np.zeros(ncomp)) for e in exps]))
        return P_, vals

    def eval(self, X):
        X = np.atleast_2d(X)
        mon = np.ones((X.shape[0], len(self.exps)))
        for j in range(4):
            mon = mon * X[:, j:j + 1] ** self.exps[None, :, j]
        return mon @ self.coef

    def deriv(self, j):
        e = self.exps.copy(); m = e[:, j]
        keep = m > 0
        e = e[keep]; c = self.coef[keep] * m[keep].reshape((-1,) + (1,) * (self.coef.ndim - 1))
        e[:, j] -= 1
        return Poly(e if len(e) else np.zeros((1, 4), int), c if len(e) else np.zeros((1,) + self.coef.shape[1:]))

    def dirderiv(self, d):
        out = None
        for j in range(4):
            if d[j] != 0:
                t = self.deriv(j).scale(d[j])
                out = t if out is None else out.add(t)
        return out if out is not None else Poly(np.zeros((1, 4), int), np.zeros((1,) + self.coef.shape[1:]))

    def scale(self, s_):
        return Poly(self.exps, self.coef * s_)

    def add(self, o):
        d = {}
        for e, c in zip(map(tuple, self.exps), self.coef): d[e] = d.get(e, 0) + c
        for e, c in zip(map(tuple, o.exps), o.coef): d[e] = d.get(e, 0) + c
        ex = list(d.keys())
        return Poly(ex, np.array([d[e] for e in ex]))

    def lap(self):
        out = None
        for j in range(4):
            t = self.deriv(j).deriv(j)
            out = t if out is None else out.add(t)
        return out

    def mul_ksq(self, s_=1.0):
        d = {}
        for e, c in zip(map(tuple, self.exps), self.coef):
            for j in range(4):
                e2 = list(e); e2[j] += 2; e2 = tuple(e2)
                d[e2] = d.get(e2, 0) + s_ * c
        ex = list(d.keys())
        return Poly(ex, np.array([d[e] for e in ex]))

    def degree(self):
        nz = np.abs(self.coef).reshape(len(self.exps), -1).max(axis=1) > 1e-13 * max(1.0, np.abs(self.coef).max())
        return int(self.exps[nz].sum(axis=1).max()) if nz.any() else 0

    def prune(self, tol=1e-13):
        m = np.abs(self.coef).reshape(len(self.exps), -1).max(axis=1)
        keep = m > tol * max(1.0, m.max())
        if not keep.any():
            keep[0] = True
        return Poly(self.exps[keep], self.coef[keep])


# ------------------------------------------------------------------ Gaussian + Schwinger assembly
class GraphPoly:
    """numerator jets of one cycle along (d1, d2): N(k, P = e1 d1 + e2 d2) -> 4 polynomials in k; plus the shift data"""
    def __init__(self, cyc, d1, d2, verify=True):
        self.cyc = cyc; self.d1 = np.asarray(d1, float); self.d2 = np.asarray(d2, float)
        Dm = 2 * cyc.n
        poly, vals = Poly.fit(cyc.value, Dm, 4)
        if verify:
            rng = np.random.default_rng(11)
            for _ in range(2):
                kk = rng.normal(size=4) * 1.7
                dv = cyc.value(kk); pv = poly.eval(kk)[0]
                assert np.allclose(dv, pv, rtol=1e-9, atol=1e-9 * max(1, np.abs(dv).max())), ('poly fit', dv, pv)
        if cyc.special:
            poly = poly.mul_ksq(-1.0)            # (e^{-eps k^2} - 1) = -int_0^eps ds k^2 e^{-s k^2}, k = k_0
        self.poly = poly.prune()
        # shifts: q_i = a_i + beta_i P ; extract a_i and beta_i (beta from the d1 part)
        self.a = np.array([q[0] for q in cyc.q])
        bet = []
        for q in cyc.q:
            if q[1] is None:
                bet.append(0.0)
            else:
                b_ = q[1] @ self.d1 / (self.d1 @ self.d1) if self.d1 @ self.d1 else 0.0
                assert np.allclose(q[1], b_ * self.d1) and (q[2] is None or np.allclose(q[2], b_ * self.d2))
                bet.append(b_)
        self.beta = np.array(bet)
        # polynomial list for the hyper-dual assembly (all components share exps)
        N = [Poly(self.poly.exps, self.poly.coef[:, c]) for c in range(4)]
        dd1, dd2 = self.d1, self.d2
        flist = {'N0': N[0], 'N1': N[1], 'N2': N[2], 'N12': N[3],
                 'd1N0': N[0].dirderiv(dd1), 'd2N0': N[0].dirderiv(dd2),
                 'd1N2': N[2].dirderiv(dd1), 'd2N1': N[1].dirderiv(dd2),
                 'd12N0': N[0].dirderiv(dd1).dirderiv(dd2)}
        deg = self.poly.degree()
        self.M = deg // 2
        # stack Laplacian powers
        self.names = list(flist)
        allp = []
        for nm in self.names:
            p_ = flist[nm]
            for m in range(self.M + 1):
                allp.append(p_)
                p_ = p_.lap()
        # merge into one coefficient matrix over a common exponent set
        exset = {}
        for p_ in allp:
            for e in map(tuple, p_.exps):
                if e not in exset: exset[e] = len(exset)
        coefm = np.zeros((len(exset), len(allp)), complex)
        for j, p_ in enumerate(allp):
            for e, c in zip(map(tuple, p_.exps), p_.coef):
                coefm[exset[e], j] += c
        self.stack = Poly(list(exset.keys()), coefm)

    def integrate(self, t, w, chunk=4096):
        """t: (npts, n) Schwinger parameters per edge (edge 0 = special s if special), w: (npts,) weights.
        returns hyper-dual integral [G, d1.dG, d2.dG, d1 d2 : d^2 G] (units hbar/(16 pi^2))"""
        out = np.zeros(4, complex)
        nm = self.names; M = self.M; nf = len(nm)
        for s0 in range(0, len(w), chunk):
            tt = t[s0:s0 + chunk]; ww = w[s0:s0 + chunk]
            tau = tt.sum(axis=1); u = tt / tau[:, None]
            K0 = u @ self.a
            kap = u @ self.beta
            q0 = u @ (self.a ** 2).sum(axis=1) - (K0 ** 2).sum(axis=1)
            wv = u @ (self.beta[:, None] * self.a) - kap[:, None] * K0
            lam = u @ self.beta ** 2 - kap ** 2
            vals = self.stack.eval(-K0).reshape(len(ww), nf, M + 1)
            fac = np.array([1.0 / (4.0 ** m * math.factorial(m)) for m in range(M + 1)])
            taupow = tau[:, None] ** (-np.arange(M + 1))[None, :]
            E = np.einsum('pfm,pm,m->pf', vals, taupow, fac)
            Ed = {name: E[:, i] for i, name in enumerate(nm)}
            E0 = Ed['N0']
            E1 = Ed['N1'] - kap * Ed['d1N0']
            E2 = Ed['N2'] - kap * Ed['d2N0']
            E12 = Ed['N12'] - kap * Ed['d1N2'] - kap * Ed['d2N1'] + kap ** 2 * Ed['d12N0']
            X0 = np.exp(-tau * q0)
            wd1 = wv @ self.d1; wd2 = wv @ self.d2
            X1 = X0 * (-2 * tau * wd1); X2 = X0 * (-2 * tau * wd2)
            X12 = X0 * (4 * tau ** 2 * wd1 * wd2 - 2 * tau * lam * (self.d1 @ self.d2))
            g = ww * tau ** -2.0
            out += np.array([np.sum(g * E0 * X0), np.sum(g * (E1 * X0 + E0 * X1)), np.sum(g * (E2 * X0 + E0 * X2)),
                             np.sum(g * (E12 * X0 + E1 * X2 + E2 * X1 + E0 * X12))])
        return out


# ------------------------------------------------------------------ quadrature grids
def gl(n):
    x, w = np.polynomial.legendre.leggauss(n)
    return x, w


def log_grid(lo, hi, per_decade=10):
    """Gauss-Legendre in log t on [lo, hi], panels of at most one decade"""
    edges = [lo]
    e = 10.0 ** np.floor(np.log10(lo) + 1)
    while e < hi * (1 - 1e-12):
        if e > lo * (1 + 1e-12):
            edges.append(e)
        e *= 10.0
    edges.append(hi)
    x, w = gl(per_decade)
    ts, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        la, lb = np.log(a), np.log(b)
        lt = 0.5 * (lb - la) * x + 0.5 * (lb + la)
        ts.append(np.exp(lt)); ws.append(0.5 * (lb - la) * w * np.exp(lt))
    return np.concatenate(ts), np.concatenate(ws)


def lin_grid(lo, hi, n):
    x, w = gl(n)
    return 0.5 * (hi - lo) * x + 0.5 * (hi + lo), 0.5 * (hi - lo) * w


# ------------------------------------------------------------------ two independent injected momenta (arity two)
class GraphPoly2(GraphPoly):
    """as GraphPoly, but the e1 and e2 parts of the line shifts are independent: q_i = a_i + beta_i P1 + gam_i P2,
    P1 = e1 d1, P2 = e2 d2 (P1^2 = P2^2 = 0).  For beta = gam this reduces to GraphPoly."""
    def __init__(self, cyc, d1, d2, verify=True):
        self.cyc = cyc; self.d1 = np.asarray(d1, float); self.d2 = np.asarray(d2, float)
        Dm = 2 * cyc.n
        poly, vals = Poly.fit(cyc.value, Dm, 4)
        if verify:
            rng = np.random.default_rng(11)
            kk = rng.normal(size=4) * 1.3
            dv = cyc.value(kk); pv = poly.eval(kk)[0]
            assert np.allclose(dv, pv, rtol=1e-8, atol=1e-8 * max(1, np.abs(dv).max())), ('poly fit', dv, pv)
        if cyc.special:
            poly = poly.mul_ksq(-1.0)
        self.poly = poly.prune()
        self.a = np.array([q[0] for q in cyc.q])
        b1, b2 = [], []
        for q in cyc.q:
            b1.append(0.0 if q[1] is None else (q[1] @ self.d1) / (self.d1 @ self.d1))
            b2.append(0.0 if q[2] is None else (q[2] @ self.d2) / (self.d2 @ self.d2))
            assert q[1] is None or np.allclose(q[1], b1[-1] * self.d1)
            assert q[2] is None or np.allclose(q[2], b2[-1] * self.d2)
        self.beta = np.array(b1); self.gam = np.array(b2)
        N = [Poly(self.poly.exps, self.poly.coef[:, c]) for c in range(4)]
        flist = {'N0': N[0], 'N1': N[1], 'N2': N[2], 'N12': N[3],
                 'd1N0': N[0].dirderiv(self.d1), 'd2N0': N[0].dirderiv(self.d2),
                 'd1N2': N[2].dirderiv(self.d1), 'd2N1': N[1].dirderiv(self.d2),
                 'd12N0': N[0].dirderiv(self.d1).dirderiv(self.d2)}
        self.M = self.poly.degree() // 2
        self.names = list(flist)
        allp = []
        for nm in self.names:
            p_ = flist[nm]
            for m in range(self.M + 1):
                allp.append(p_); p_ = p_.lap()
        exset = {}
        for p_ in allp:
            for e in map(tuple, p_.exps):
                if e not in exset: exset[e] = len(exset)
        coefm = np.zeros((len(exset), len(allp)), complex)
        for j, p_ in enumerate(allp):
            for e, c in zip(map(tuple, p_.exps), p_.coef):
                coefm[exset[e], j] += c
        self.stack = Poly(list(exset.keys()), coefm)

    def integrate(self, t, w, chunk=4096):
        out = np.zeros(4, complex)
        nm = self.names; M = self.M; nf = len(nm)
        for s0 in range(0, len(w), chunk):
            tt = t[s0:s0 + chunk]; ww = w[s0:s0 + chunk]
            tau = tt.sum(axis=1); u = tt / tau[:, None]
            K0 = u @ self.a
            k1 = u @ self.beta; k2 = u @ self.gam
            q0 = u @ (self.a ** 2).sum(axis=1) - (K0 ** 2).sum(axis=1)
            w1 = u @ (self.beta[:, None] * self.a) - k1[:, None] * K0
            w2 = u @ (self.gam[:, None] * self.a) - k2[:, None] * K0
            l12 = u @ (self.beta * self.gam) - k1 * k2
            vals = self.stack.eval(-K0).reshape(len(ww), nf, M + 1)
            fac = np.array([1.0 / (4.0 ** m * math.factorial(m)) for m in range(M + 1)])
            taupow = tau[:, None] ** (-np.arange(M + 1))[None, :]
            E = np.einsum('pfm,pm,m->pf', vals, taupow, fac)
            Ed = {name: E[:, i] for i, name in enumerate(nm)}
            E0 = Ed['N0']; E1 = Ed['N1'] - k1 * Ed['d1N0']; E2 = Ed['N2'] - k2 * Ed['d2N0']
            E12 = Ed['N12'] - k1 * Ed['d1N2'] - k2 * Ed['d2N1'] + k1 * k2 * Ed['d12N0']
            X0 = np.exp(-tau * q0)
            wd1 = w1 @ self.d1; wd2 = w2 @ self.d2
            X1 = X0 * (-2 * tau * wd1); X2 = X0 * (-2 * tau * wd2)
            X12 = X0 * (4 * tau ** 2 * wd1 * wd2 - 2 * tau * l12 * (self.d1 @ self.d2))
            g = ww * tau ** -2.0
            out += np.array([np.sum(g * E0 * X0), np.sum(g * (E1 * X0 + E0 * X1)), np.sum(g * (E2 * X0 + E0 * X2)),
                             np.sum(g * (E12 * X0 + E1 * X2 + E2 * X1 + E0 * X12))])
        return out


def gamma_pair(z1, z2):
    """Gamma(z1,z2)^M = z1^T sb[M] z2 for z1, z2 in Sigma (ten components)"""
    return np.array([z1 @ SBARR[M] @ z2 for M in range(10)])


def s1eta_hessian(zeta, kj, kl, P=None):
    """quadratic part of S^(1) along the Killing spinor zeta e^{iPy}, zeta in Sigma ([ch02 eq:c2-S1Q] with Q_1 -> zeta):
       -(Gamma(eps,lambda), A^*) + 1/2 (rho(dA) eps, lambda^*),  rho(q(dA)) on plane waves = i (s[k] sb[a] - (k,a)).
    Even functional (the ghost of an odd generator is even); H_jl = (-1)^{|j||l|} coefficient(j, l)."""
    H = np.zeros((NT, NT), complex)
    zeta = np.asarray(zeta, complex)
    sbz = np.array([SBARR[M] @ zeta for M in range(10)])          # (M, beta): (sb[M] zeta)_beta
    def T2vec(k, M):    # (i/2)(s[k] sb[e_M] - k_M) zeta  in Sigma
        sk = np.einsum('m,mab->ab', np.asarray(k[:4], complex), SARR[:4])
        km = k[M] if M < 4 else 0.0
        return 0.5j * (sk @ (SBARR[M] @ zeta) - km * zeta)
    for a in range(D):
        for M in range(10):
            for be in range(16):
                il = C.OFF['l'] + be * D + a; iS = C.OFF['As'] + M * D + a
                H[il, iS] += sbz[M, be]           # j = lambda, l = A^*:  (-1)^{1}( - (sb[M] zeta)_be )
                H[iS, il] += -sbz[M, be]          # j = A^*, l = lambda
            vj = T2vec(kj, M); vl = T2vec(kl, M)
            for be in range(16):
                iA = C.OFF['A'] + M * D + a; isl = C.OFF['ls'] + be * D + a
                H[iA, isl] += vj[be]              # j = A (at kj), l = lambda^*
                H[isl, iA] += vl[be]              # j = lambda^*, l = A (at kl)
    return H


def s1eta_const_hessian(eta):
    """the last term of [ch02 eq:c2-S1S]: -(rho(A|_W) eta, lambda^*), eta in Sigma^*; legs A|_W, lambda^*"""
    H = np.zeros((NT, NT), complex)
    for a in range(D):
        for M in range(4, 10):
            v = -(SARR[M] @ np.asarray(eta, complex))
            for be in range(16):
                iA = C.OFF['A'] + M * D + a; isl = C.OFF['ls'] + be * D + a
                H[iA, isl] += v[be]; H[isl, iA] += v[be]
    return H


def s2_hessian(z1, z2, kj, kl, P=None):
    """S^(2)(eps1, eps2) of [ch02 eq:c2-S2] for spinors z1, z2 in Sigma (pointwise values), quadratic:
       1/4 (Gamma(z1,z2), Gamma(lambda^*,lambda^*)) - 1/2 (z1,lambda^*)(z2,lambda^*) - (iota_{Gamma(z1,z2)} A, c^*).
    Gamma(lambda^*,lambda^*)^M = lambda^*^T s[M] lambda^* (pairing on Sigma^*), (z,lambda^*) = z^T lambda^*.
    Even functional; H_jl = (-1)^{|j||l|} coefficient."""
    H = np.zeros((NT, NT), complex)
    G = gamma_pair(np.asarray(z1, complex), np.asarray(z2, complex))
    Ms = np.einsum('M,Mab->ab', G, SARR)                             # sum_M G^M s[M]  (on Sigma^*)
    z1 = np.asarray(z1, complex); z2 = np.asarray(z2, complex)
    B = 0.25 * (Ms + Ms.T) - 0.5 * (np.outer(z1, z2) + np.outer(z2, z1))   # lambda^*: lambda^*^T B lambda^* ... second derivative
    for a in range(D):
        for al in range(16):
            for be in range(16):
                i1 = C.OFF['ls'] + al * D + a; i2 = C.OFF['ls'] + be * D + a
                H[i1, i2] += B[al, be]
        for M in range(10):
            iA = C.OFF['A'] + M * D + a; ic = C.OFF['cs'] + a
            H[iA, ic] += -G[M]; H[ic, iA] += -G[M]
    return H


class InsertionGen(Insertion):
    """Insertion from an arbitrary Hessian function f(kj, kl, P) (affine)"""
    def __init__(self, f):
        self.kind = 'gen'
        z = np.zeros(4)
        self.h0 = f(z, z, z)
        self.Hj = np.array([f(E4[m], z, z) - self.h0 for m in range(4)])
        self.Hl = np.array([f(z, E4[m], z) - self.h0 for m in range(4)])
        self.Hp = np.array([f(z, z, E4[m]) - self.h0 for m in range(4)])
        kj, kl, p = np.array([.3, -1.2, .7, 2.]), np.array([-.4, .9, 1.1, -.6]), np.array([1.3, .2, -.8, .5])
        lin = self.h0 + np.tensordot(kj, self.Hj, 1) + np.tensordot(kl, self.Hl, 1) + np.tensordot(p, self.Hp, 1)
        assert np.allclose(lin, f(kj, kl, p)), 'insertion not affine'
        self.cache = {}
