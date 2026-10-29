"""kn_engine: independent one-loop engine for the K-route at a finite transfer endpoint T (Step 2, K_1 numerics).

Imports ONLY the Step-1 core of engine A (engine_a/r2_core.py): field layout (c, A, A^*, lambda, lambda^*, c^*) with
su(2) in an orthonormal basis of (-,-), vertex functions of I (cubic scaled by 2, quartic by 8), propagator
numerators Pnum(k, line) of ch03 prop:c3-P, the kernel Khat(line) of K_t, the parity diagonal SDIAG.
Nothing from step2 engine outputs and nothing from the partial files of the earlier run is used.

One-loop formulas (checked by brute force in 0 dimensions, kn_0d_cycles.py, every vertex quadratic in the
fluctuation so that only cycles occur; H = Hessian with left derivatives, S = parity, str = supertrace):
  flow cycle over n distinct insertions:  (1/(2n)) sum_{sigma in S_n} str[ P S H_s1 P S H_s2 ... P S H_sn ]
  special-edge cycle (defect (1/2){J,J}_{K'} closed by P):  (1/2) sum_{sigma in S_n} str[ K' S H_s1 P S H_s2 ... ]
  (K' = K_eps - K_0 on line 0, entering insertion s1).
Line i enters insertion i through its j-leg (incoming momentum -k_i) and leaves insertion i-1 through its l-leg
(incoming k_i); k_{i+1} = k_i - Pin_i.  Odd background parameters are factored to the left; the value is the
coefficient of the canonical product of odd parameters (perm sign).

Kernels at a finite endpoint T (momentum space; ch08 def:c8-PT):
  P(eps,T)(k) = Pnum(k) int_eps^T dt e^{-t k^2};   K'(k) = Khat (e^{-eps k^2} - 1)  (t_0 = eps minus t_0 = 0).
Gaussian loop integral (exact): int d^4k/(2pi)^4 N(k) e^{-sum t_i (k+q_i)^2}
   = (16 pi^2)^{-1} tau^{-2} e^{-tau V} sum_m (Delta^m N)(k_0)/(4^m m! tau^m),  tau = sum t_i, k_0 = -sum u_i q_i,
   V = sum u_i q_i^2 - |sum u_i q_i|^2.  Units: hbar/(16 pi^2) (the factor hbar and (16 pi^2)^{-1} stripped).
Injected momenta (plane-wave insertions u e^{iP.y}, spinors zeta e^{iP.x}) enter as jets: P = sum_v z_v w_v,
the jet algebra truncated by an allowed exponent set; position-dependent generators are reached by
y^a -> -i d/dP_a at P = 0."""
import sys, os, itertools, math
from fractions import Fraction
import numpy as np

np.seterr(all='ignore')
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
sys.path.insert(0, os.path.join(ROOT, 'engine_a'))
import r2_core as C  # noqa: E402

D = C.D
NT = C.NTOT
SD = C.SDIAG.real.copy()
PAR = dict(C.PARITY)
SARR = C.SARR
SBARR = C.SBARR
PLINES = ['AV', 'AP', 'gh', 'l']
SLINES = ['sAV', 'sAP', 'sl', 'sc']
E4 = np.eye(4)


# ============================================================== jet algebra
class Jet:
    """truncated polynomial algebra in abstract variables z_0..z_{nvar-1}; basis = allowed exponents (closed under
    divisors); coefficients are monomial coefficients (value = sum_e c_e z^e)."""

    def __init__(self, nvar, allowed):
        self.nvar = nvar
        self.basis = [tuple(int(x) for x in e) for e in allowed]
        assert self.basis[0] == (0,) * nvar
        self.index = {e: i for i, e in enumerate(self.basis)}
        self.nb = len(self.basis)
        tab = []
        for i, ei in enumerate(self.basis):
            for j, ej in enumerate(self.basis):
                s = tuple(a + b for a, b in zip(ei, ej))
                if s in self.index:
                    tab.append((i, j, self.index[s]))
        self.tab = tab
        self.lin = []
        for v in range(nvar):
            e = tuple(1 if k == v else 0 for k in range(nvar))
            self.lin.append(self.index.get(e))
        self.order = [sum(e) for e in self.basis]
        self.maxorder = max(self.order)
        # shift tables: z^{e_b} * X  -> list of (src j, dst k)
        self.shift_tab = []
        for b, eb in enumerate(self.basis):
            lst = []
            for j, ej in enumerate(self.basis):
                s = tuple(x + y for x, y in zip(eb, ej))
                if s in self.index:
                    lst.append((j, self.index[s]))
            self.shift_tab.append(lst)

    def shift(self, b, X):
        out = np.zeros_like(X)
        for j, k in self.shift_tab[b]:
            out[k] = X[j]
        return out

    # scalar jets: arrays (nb, ...)
    def mul(self, A, B):
        sh = np.broadcast(A[0], B[0]).shape
        out = np.zeros((self.nb,) + sh, complex)
        for i, j, k in self.tab:
            out[k] = out[k] + A[i] * B[j]
        return out

    def exp0(self, X):
        """exp of a jet with zero constant term (truncated)"""
        out = np.zeros_like(X)
        out[0] = 1.0
        term = out.copy()
        for n in range(1, self.maxorder + 1):
            term = self.mul(term, X) / n
            out = out + term
        return out


def jet_order2(nvar=4):
    allowed = [e for e in itertools.product(range(3), repeat=nvar) if sum(e) <= 2]
    allowed.sort(key=lambda e: (sum(e), [-x for x in e]))
    return Jet(nvar, allowed)


def jet_hyperdual():
    return Jet(2, [(0, 0), (1, 0), (0, 1), (1, 1)])


def jet_trivial():
    return Jet(0, [()])


def jet_order1(nvar=4):
    allowed = [tuple([0] * nvar)] + [tuple(1 if k == v else 0 for k in range(nvar)) for v in range(nvar)]
    return Jet(nvar, allowed)


class JMat:
    """matrix-valued jet: dict basis-index -> ndarray (batch, r, c) or (r, c)"""

    def __init__(self, jet, comps):
        self.jet = jet
        self.c = {i: m for i, m in comps.items() if m is not None}

    def matmul(self, o):
        out = {}
        for i, j, k in self.jet.tab:
            a = self.c.get(i)
            if a is None:
                continue
            b = o.c.get(j)
            if b is None:
                continue
            p = np.matmul(a, b)
            out[k] = p if k not in out else out[k] + p
        return JMat(self.jet, out)


# ============================================================== polynomials in k (4 variables)
def simplex(Dg):
    pts = [a for a in itertools.product(range(Dg + 1), repeat=4) if sum(a) <= Dg]
    pts.sort(key=lambda a: (sum(a), a))
    return pts


_TM_CACHE = {}


def fit_matrix(Dg):
    """matrix Tm with monomial coefficients = Tm @ (values on simplex(Dg)) for polynomials of degree <= Dg"""
    if Dg in _TM_CACHE:
        return _TM_CACHE[Dg]
    pts = simplex(Dg)
    idx = {p: i for i, p in enumerate(pts)}
    n = len(pts)
    # Newton forward differences ND[alpha, beta]
    ND = np.zeros((n, n))
    for a in pts:
        ia = idx[a]
        for b in itertools.product(*[range(x + 1) for x in a]):
            c = (-1) ** (sum(a) - sum(b)) * np.prod([math.comb(a[j], b[j]) for j in range(4)])
            ND[ia, idx[b]] += c
    # Stirling-type conversion: C(x, a) = sum_r s(a, r) x^r
    smax = Dg
    S1 = [[Fraction(0)] * (smax + 1) for _ in range(smax + 1)]
    for a in range(smax + 1):
        poly = [Fraction(1)]
        for j in range(a):
            new = [Fraction(0)] * (len(poly) + 1)
            for r, cf in enumerate(poly):
                new[r + 1] += cf
                new[r] += -j * cf
            poly = new
        for r, cf in enumerate(poly):
            S1[a][r] = cf / math.factorial(a)
    ST = np.zeros((n, n))
    for r in pts:
        for a in pts:
            if all(a[j] >= r[j] for j in range(4)):
                v = Fraction(1)
                for j in range(4):
                    v *= S1[a[j]][r[j]]
                ST[idx[r], idx[a]] = float(v)
    Tm = ST @ ND
    _TM_CACHE[Dg] = (pts, Tm)
    return pts, Tm


class Poly:
    """sum_i coef[i] k^{exps[i]} ; coef has trailing component axes"""

    def __init__(self, exps, coef):
        self.exps = np.asarray(exps, int).reshape(-1, 4)
        self.coef = np.asarray(coef, complex)

    def eval(self, X):
        X = np.atleast_2d(np.asarray(X, float))
        mx = int(self.exps.max()) if self.exps.size else 0
        pw = np.ones((X.shape[0], 4, mx + 1))
        for p in range(1, mx + 1):
            pw[:, :, p] = pw[:, :, p - 1] * X
        mon = np.ones((X.shape[0], len(self.exps)))
        for j in range(4):
            mon = mon * pw[:, j, :][:, self.exps[:, j]]
        return np.tensordot(mon, self.coef, axes=(1, 0))

    def _from_dict(self, d):
        if not d:
            return Poly(np.zeros((1, 4), int), np.zeros((1,) + self.coef.shape[1:], complex))
        ex = list(d.keys())
        return Poly(ex, np.array([d[e] for e in ex]))

    def deriv(self, j):
        d = {}
        for e, c in zip(map(tuple, self.exps), self.coef):
            if e[j] > 0:
                e2 = list(e); e2[j] -= 1; e2 = tuple(e2)
                d[e2] = d.get(e2, 0) + e[j] * c
        return self._from_dict(d)

    def lap(self):
        d = {}
        for e, c in zip(map(tuple, self.exps), self.coef):
            for j in range(4):
                if e[j] >= 2:
                    e2 = list(e); e2[j] -= 2; e2 = tuple(e2)
                    d[e2] = d.get(e2, 0) + e[j] * (e[j] - 1) * c
        return self._from_dict(d)

    def degree(self, tol=1e-9):
        m = np.abs(self.coef).reshape(len(self.exps), -1).max(axis=1)
        if not np.any(m > tol * max(1.0, m.max() if m.size else 1.0)):
            return 0
        return int(self.exps[m > tol * max(1.0, m.max())].sum(axis=1).max())


def fit_poly(values, Dg):
    """values: array (npts, ...) on simplex(Dg) -> Poly"""
    pts, Tm = fit_matrix(Dg)
    sh = values.shape[1:]
    coef = Tm @ values.reshape(len(pts), -1)
    keep = np.abs(coef).max(axis=1) > 1e-10 * max(1.0, np.abs(coef).max())
    if not keep.any():
        keep[0] = True
    ex = np.array(pts)[keep]
    return Poly(ex, coef[keep].reshape((-1,) + sh))


# ============================================================== lines
def _support(L):
    kk = np.array([0.31, -0.77, 1.13, 0.29])
    N = C.Pnum(kk, L) if L in PLINES else C.Khat(L)
    rows = np.nonzero(np.abs(N).sum(axis=1))[0]
    cols = np.nonzero(np.abs(N).sum(axis=0))[0]
    return rows, cols


SUPPORT = {L: _support(L) for L in PLINES + SLINES}


class Line:
    def __init__(self, L):
        self.L = L
        self.rows, self.cols = SUPPORT[L]
        ix = np.ix_(self.rows, self.cols)
        if L in PLINES:
            z = np.zeros(4)
            self.N0 = C.Pnum(z, L)[ix].astype(complex)
            self.Nm = np.array([C.Pnum(E4[m], L)[ix] - C.Pnum(z, L)[ix] for m in range(4)]).astype(complex)
        else:
            self.N0 = C.Khat(L)[ix].astype(complex)
            self.Nm = np.zeros((4,) + self.N0.shape, complex)
        self.special = L in SLINES
        self.Scols = SD[self.cols]


LINES = {L: Line(L) for L in PLINES + SLINES}


# ============================================================== insertions
class Insertion:
    """affine Hessian h(kj, kl, pin) = h0 + kj.Hj + kl.Hl + pin.Hp (full NT x NT, scaled by `scale`).
    tags: set of labels (background legs / generator arguments) carried by this insertion;
    odd: list of odd-parameter labels carried, in the order they stand in the vertex function."""

    def __init__(self, fn, scale, tags, odd, name, pin_const=None, pin_lin=None):
        z = np.zeros(4)
        self.name = name
        self.scale = scale
        self.tags = tuple(tags)
        self.odd = list(odd)
        self.h0 = fn(z, z, z)
        self.Hj = np.array([fn(E4[m], z, z) - self.h0 for m in range(4)])
        self.Hl = np.array([fn(z, E4[m], z) - self.h0 for m in range(4)])
        self.Hp = np.array([fn(z, z, E4[m]) - self.h0 for m in range(4)])
        kj = np.array([.3, -1.2, .7, 2.]); kl = np.array([-.4, .9, 1.1, -.6]); p = np.array([1.3, .2, -.8, .5])
        lin = self.h0 + np.tensordot(kj, self.Hj, 1) + np.tensordot(kl, self.Hl, 1) + np.tensordot(p, self.Hp, 1)
        assert np.allclose(lin, fn(kj, kl, p)), 'insertion %s not affine' % name
        # injected momentum pin(z) = pin_const + sum_v z_v pin_lin[v]
        self.pin_const = np.zeros(4) if pin_const is None else np.asarray(pin_const, float)
        self.pin_lin = pin_lin          # array (nvar, 4) or None (set per run)
        self.cache = {}

    def block(self, cols, rows):
        key = (tuple(cols), tuple(rows))
        if key not in self.cache:
            ix = np.ix_(cols, rows)
            self.cache[key] = (self.h0[ix], self.Hj[(slice(None),) + ix], self.Hl[(slice(None),) + ix],
                               self.Hp[(slice(None),) + ix])
        return self.cache[key]

    def nonzero_between(self, Lin, Lout):
        h0, Hj, Hl, Hp = self.block(LINES[Lin].cols, LINES[Lout].rows)
        return bool(np.abs(h0).sum() + np.abs(Hj).sum() + np.abs(Hl).sum() + np.abs(Hp).sum() > 0)


def bg_insertion(t, pol, mom_const, tag, mom_lin=None, quartic_with=None):
    """Hessian of I at one background leg (cubic part, scale 2) or two A legs (quartic, scale 8)"""
    pol = np.asarray(pol, complex)
    if quartic_with is None:
        def fn(kj, kl, p):
            return C.hessian([(t, pol, p)], kj, kl)
        ins = Insertion(fn, 2.0, [tag], [tag] if PAR[t] else [], 'bg:%s' % tag, mom_const, mom_lin)
    else:
        (t2, pol2, tag2) = quartic_with
        pol2 = np.asarray(pol2, complex)
        z = np.zeros(4)

        def fn(kj, kl, p):
            return C.hessian([(t, pol, z), (t2, pol2, z)], kj, kl, fluct=['A'])
        ins = Insertion(fn, 8.0, [tag, tag2], [], 'bg4:%s%s' % (tag, tag2), mom_const, mom_lin)
    return ins


# ---------------------------------------------------------------- S^(1)(v) for v = u e^{iP.y} (ch02 eq:c2-S1even)
def s1v_hessian(u, kj, kl, P):
    """8 x Hessian of S^(1)(u e^{iP.y}) gamma (gamma = odd ghost of the even generator, placed last in the
    canonical order).  S^(1)(v) = int (L_v A + 1/4 div v A|_W, A^*) - (v.lambda, lambda^*) - (v.c, c^*),
    v.lambda = v(lambda) + 1/4 rho(d(v,-)) lambda + 3/8 div v lambda;  on plane waves d -> i k,
    v(f) -> i (u.k) f, div v -> i (u.P), d(v,-) -> i P ^ u, rho(q(P^u)) = rho(P) rho(u) - (P.u).
    H_jl = (-1)^{|j||l|} (coefficient of the odd parameters in leg order (j, l, gamma)); leg j has incoming kj,
    leg l incoming kl."""
    u = np.asarray(u, float)
    H = np.zeros((NT, NT), complex)
    uP = float(np.dot(u, P[:4]))
    sP = np.einsum('m,mab->ab', np.asarray(P[:4], complex), SARR[:4])     # rho(P) : Sigma^* -> Sigma
    sbu = np.einsum('m,mab->ab', np.asarray(u, complex), SBARR[:4])       # rho(u) : Sigma -> Sigma^*

    def Mlam(k):        # 8 * (v.lambda) operator on Sigma at lambda momentum k
        return (8j * np.dot(u, k) + 1j * uP) * np.eye(16) + 2j * (sP @ sbu)
    for a in range(D):
        for M in range(10):
            for N in range(10):
                if M < 4 and N < 4:
                    vj = 8j * np.dot(u, kj) * (M == N) + 8j * u[M] * P[N]
                    vl = 8j * np.dot(u, kl) * (M == N) + 8j * u[M] * P[N]
                elif M >= 4 and N >= 4 and M == N:
                    vj = 8j * np.dot(u, kj) + 2j * uP
                    vl = 8j * np.dot(u, kl) + 2j * uP
                else:
                    continue
                iA = C.OFF['A'] + M * D + a
                iS = C.OFF['As'] + N * D + a
                H[iA, iS] += vj           # j = A (even) at kj, l = A^* (odd): (-1)^0
                H[iS, iA] += vl           # j = A^*, l = A at kl: (-1)^0
        Mj, Ml = Mlam(kj), Mlam(kl)
        for al in range(16):
            for be in range(16):
                il = C.OFF['l'] + al * D + a
                isl = C.OFF['ls'] + be * D + a
                H[il, isl] += -Mj[be, al]     # j = lambda (odd) at kj, l = lambda^* (even): -(v.lambda, lambda^*)
                H[isl, il] += -Ml[be, al]     # j = lambda^*, l = lambda at kl
        ic = C.OFF['c'] + a
        ics = C.OFF['cs'] + a
        H[ic, ics] += -8j * np.dot(u, kj)     # -(v.c, c^*), c at kj
        H[ics, ic] += -8j * np.dot(u, kl)
    return H


def s1v_insertion(u, tag='g', pin_lin=None):
    u = np.asarray(u, float)

    def fn(kj, kl, p):
        return s1v_hessian(u, kj, kl, p)
    return Insertion(fn, 8.0, [tag], [tag], 'S1v', None, pin_lin)


# ============================================================== cycles
def perm_sign(seq, canon):
    pos = [canon.index(x) for x in seq]
    s_ = 1
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            if pos[i] > pos[j]:
                s_ = -s_
    return s_


class Cycle:
    """insertions in cyclic order with line types; line i enters insertion i (j-leg) and leaves insertion i-1.
    Line momenta k_i(z) = k + a_i + sum_v z_v B_i[v]; k_0 = k; k_{i+1} = k_i - pin_i."""

    def __init__(self, ins, lines, jet, canon_odd, special):
        self.ins = ins
        self.Ls = list(lines)
        self.lines = [LINES[L] for L in lines]
        self.n = len(ins)
        self.jet = jet
        self.special = special
        assert (self.lines[0].special == special) and all(not L.special for L in self.lines[1:])
        odd = [lab for I in ins for lab in I.odd]
        self.sign = perm_sign(odd, canon_odd) if odd else 1
        self.scale = float(np.prod([I.scale for I in ins]))
        nv = jet.nvar
        a = [np.zeros(4)]
        B = [np.zeros((nv, 4))]
        for i in range(self.n - 1):
            I = ins[i]
            a.append(a[-1] - I.pin_const)
            B.append(B[-1] - (I.pin_lin if I.pin_lin is not None else np.zeros((nv, 4))))
        self.a = np.array(a)            # (n, 4)
        self.B = np.array(B)            # (n, nv, 4)
        # momentum conservation check: k_n = k_0
        tot = sum(I.pin_const for I in ins)
        totl = sum((I.pin_lin if I.pin_lin is not None else np.zeros((nv, 4))) for I in ins)
        assert np.allclose(tot, 0) and np.allclose(totl, 0), 'momentum not conserved in cycle'
        self.deg = sum(1 for L in self.Ls if L in ('gh', 'l')) + self.n     # numerator degree bound in k

    def factors(self, K):
        """jet-matrix factors N_{L_i}(k_i) S h_i for a batch of loop momenta K (npts, 4)"""
        jet = self.jet
        out = []
        n = self.n
        npts = K.shape[0]
        for i in range(n):
            Li = self.lines[i]
            Ln = self.lines[(i + 1) % n]
            I = self.ins[i]
            kq = K + self.a[i][None, :]                       # (npts, 4)  k_i const part
            kn = K + self.a[(i + 1) % n][None, :]
            # line numerator jet
            Nc = Li.N0[None] + np.tensordot(kq, Li.Nm, axes=(1, 0))          # (npts, r, c)
            Nj = {0: Nc}
            for v in range(jet.nvar):
                bv = self.B[i][v]
                if np.any(bv) and np.any(Li.Nm):
                    Nj[jet.lin[v]] = np.tensordot(bv, Li.Nm, axes=(0, 0))[None]
            # insertion jet: kj = -k_i, kl = k_{i+1}
            h0, Hj, Hl, Hp = I.block(Li.cols, Ln.rows)
            hc = h0[None] + np.tensordot(-kq, Hj, axes=(1, 0)) + np.tensordot(kn, Hl, axes=(1, 0)) \
                + np.tensordot(I.pin_const, Hp, axes=(0, 0))[None]
            hj = {0: hc}
            for v in range(jet.nvar):
                bv_in = self.B[i][v]
                bv_out = self.B[(i + 1) % n][v]
                pv = I.pin_lin[v] if I.pin_lin is not None else np.zeros(4)
                m = np.tensordot(-bv_in, Hj, axes=(0, 0)) + np.tensordot(bv_out, Hl, axes=(0, 0)) \
                    + np.tensordot(pv, Hp, axes=(0, 0))
                if np.any(m):
                    hj[jet.lin[v]] = m[None]
            Sc = Li.Scols
            hjS = {k: Sc[None, :, None] * m for k, m in hj.items()}
            F = JMat(jet, Nj).matmul(JMat(jet, hjS))
            out.append(F)
        return out

    def trace_values(self, K):
        """str of the cycle product at loop momenta K: array (npts, nb) (raw: sign/scale/prefactor not applied)"""
        Fs = self.factors(K)
        X = Fs[0]
        for F in Fs[1:]:
            X = X.matmul(F)
        r0 = self.lines[0].rows
        Sr = SD[r0]
        out = np.zeros((K.shape[0], self.jet.nb), complex)
        for kidx, M in X.c.items():
            d = np.einsum('pii->pi', np.broadcast_to(M, (K.shape[0],) + M.shape[-2:]))
            out[:, kidx] = d @ Sr
        return out

    def numerator(self, verify=True):
        Dg = self.deg
        pts, _ = fit_matrix(Dg)
        K = np.array(pts, float)
        vals = self.trace_values(K)
        poly = fit_poly(vals, Dg)
        if verify:
            rng = np.random.default_rng(17)
            Kt = rng.normal(size=(3, 4)) * 1.3
            dv = self.trace_values(Kt)
            pv = poly.eval(Kt)
            err = np.abs(dv - pv).max() / max(1.0, np.abs(dv).max())
            assert err < 1e-8, ('polynomial fit failed', err)
        return poly

    def is_zero(self, trials=2):
        rng = np.random.default_rng(5)
        K = rng.normal(size=(trials, 4))
        v = self.trace_values(K)
        return np.abs(v).max() < 1e-11


# ============================================================== Gaussian + Schwinger assembly
class Assembled:
    """precomputed derivative polynomials of one cycle numerator for the Gaussian integral with jets"""

    def __init__(self, cyc, poly):
        self.cyc = cyc
        jet = cyc.jet
        self.jet = jet
        self.n = cyc.n
        deg = poly.degree()
        self.M = deg // 2
        # list of (comp b, multi-index alpha (tuple of k-directions), m) -> stacked poly
        entries = []
        polys = []
        for b in range(jet.nb):
            ordb = jet.order[b]
            pb = Poly(poly.exps, poly.coef[:, b])
            if np.abs(pb.coef).max() < 1e-14:
                continue
            nder = jet.maxorder - ordb
            lp = pb
            for m in range(self.M + 1):
                entries.append((b, (), m)); polys.append(lp)
                if nder >= 1:
                    d1 = [lp.deriv(a) for a in range(4)]
                    for a in range(4):
                        entries.append((b, (a,), m)); polys.append(d1[a])
                    if nder >= 2:
                        for a in range(4):
                            for c in range(a, 4):
                                entries.append((b, (a, c), m)); polys.append(d1[a].deriv(c))
                lp = lp.lap()
        exset = {}
        for p_ in polys:
            for e in map(tuple, p_.exps):
                if e not in exset:
                    exset[e] = len(exset)
        coefm = np.zeros((max(1, len(exset)), len(polys)), complex)
        for j, p_ in enumerate(polys):
            for e, c in zip(map(tuple, p_.exps), p_.coef):
                coefm[exset[e], j] += c
        self.stack = Poly(list(exset.keys()) or [(0, 0, 0, 0)], coefm)
        self.entries = entries

    def integrand(self, t):
        """t: (npts, n) Schwinger times of the n lines (t_0 may be 0). Returns jets (nb, npts) of the Gaussian
        integral (units hbar/(16 pi^2) before prefactors): tau^{-2} e^{-tau V} sum_m (Delta^m N)/(4^m m! tau^m)."""
        jet = self.jet
        cyc = self.cyc
        npts = t.shape[0]
        tau = t.sum(axis=1)
        u = t / tau[:, None]
        K0 = u @ cyc.a                                    # (npts, 4)
        k0 = -K0
        nv = jet.nvar
        # Ub[v] = sum_i u_i B_i[v]   (npts, nv, 4)
        Ub = np.einsum('pi,ivm->pvm', u, cyc.B)
        dk = -Ub                                          # dk/dz_v
        V0 = u @ (cyc.a ** 2).sum(axis=1) - (K0 ** 2).sum(axis=1)
        # V - V0 = 2 sum_v z_v w_v + sum_{v,w} z_v z_w lam_vw
        wv = np.einsum('pi,im,ivm->pv', u, cyc.a, cyc.B) - np.einsum('pm,pvm->pv', K0, Ub)
        lam = np.einsum('pi,ivm,iwm->pvw', u, cyc.B, cyc.B) - np.einsum('pvm,pwm->pvw', Ub, Ub)
        Xj = np.zeros((jet.nb, npts), complex)
        for v in range(nv):
            if jet.lin[v] is not None:
                Xj[jet.lin[v]] += -tau * 2 * wv[:, v]
        for v in range(nv):
            for w in range(nv):
                e = tuple((1 if k == v else 0) + (1 if k == w else 0) for k in range(nv))
                if e in jet.index:
                    Xj[jet.index[e]] += -tau * lam[:, v, w]
        Xexp = jet.exp0(Xj)
        # polynomial values at k0
        vals = self.stack.eval(k0)                        # (npts, nentries)
        E = np.zeros((jet.nb, npts), complex)
        fac = [1.0 / (4.0 ** m * math.factorial(m)) for m in range(self.M + 1)]
        # dk jets: dkj[a] = jet with linear comps
        dkj = np.zeros((4, jet.nb, npts), complex)
        for v in range(nv):
            if jet.lin[v] is not None:
                for a in range(4):
                    dkj[a, jet.lin[v]] = dk[:, v, a]
        dkdk = {}
        if jet.maxorder >= 2:
            for a in range(4):
                for c in range(a, 4):
                    dkdk[(a, c)] = jet.mul(dkj[a], dkj[c])
        one = np.zeros((jet.nb, npts), complex); one[0] = 1.0
        inner = {}
        for col, (b, alpha, m) in enumerate(self.entries):
            coef = vals[:, col] * fac[m] * tau ** (-m)
            if len(alpha) == 0:
                Dl = one
            elif len(alpha) == 1:
                Dl = dkj[alpha[0]]
            else:
                a, c = alpha
                Dl = dkdk[(a, c)] * (0.5 if a == c else 1.0)
            inner[b] = inner.get(b, 0) + Dl * coef[None]
        for b, X in inner.items():
            E = E + jet.shift(b, X)
        F = jet.mul(E, Xexp) * (tau ** -2.0 * np.exp(-tau * V0))[None]
        return F


# ============================================================== quadrature
def gl(n):
    return np.polynomial.legendre.leggauss(n)


def log_panels(lo, hi, per_decade):
    """Gauss-Legendre in log t on [lo, hi] with panels of at most one decade (aligned at powers of 10)"""
    edges = [lo]
    e = 10.0 ** (np.floor(np.log10(lo) + 1e-12) + 1)
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


def tensor_grid(axes):
    """axes: list of (t, w) per dimension -> (pts (N, d), weights (N,))"""
    ts = np.meshgrid(*[a[0] for a in axes], indexing='ij')
    ws = np.meshgrid(*[a[1] for a in axes], indexing='ij')
    pts = np.stack([x.ravel() for x in ts], axis=1)
    w = np.prod(np.stack([x.ravel() for x in ws], axis=1), axis=1)
    return pts, w


def integrate_cycle(asm, eps, T, per_decade=10, mode='value', chunk=20000):
    """Schwinger integral of one assembled cycle.  P lines: t in [eps, T]; special line 0: K_eps (t_0 = eps)
    minus K_0 (t_0 = 0).  mode 'value': the cycle value; mode 'TdT': T d/dT of it (sum over P lines pinned at T).
    Returns jets (nb,) in units hbar/(16 pi^2), raw (without sign, prefactor, scale)."""
    cyc = asm.cyc
    n = cyc.n
    Pidx = list(range(1, n)) if cyc.special else list(range(n))
    ax = log_panels(eps, T, per_decade)
    jetnb = asm.jet.nb
    total = np.zeros(jetnb, complex)

    def run(fixed, free):
        # fixed: dict line -> value; free: list of lines integrated over [eps, T]
        acc = np.zeros(jetnb, complex)
        if free:
            pts, w = tensor_grid([ax] * len(free))
        else:
            pts, w = np.zeros((1, 0)), np.ones(1)
        for s0 in range(0, len(w), chunk):
            pp = pts[s0:s0 + chunk]; ww = w[s0:s0 + chunk]
            t = np.zeros((len(ww), n))
            for li, val in fixed.items():
                t[:, li] = val
            for j, li in enumerate(free):
                t[:, li] = pp[:, j]
            F = asm.integrand(t)
            acc = acc + F @ ww
        return acc

    if mode == 'value':
        if cyc.special:
            total = run({0: eps}, Pidx) - run({0: 0.0}, Pidx)
        else:
            total = run({}, Pidx)
    elif mode == 'TdT':
        for pin in Pidx:
            free = [x for x in Pidx if x != pin]
            if cyc.special:
                total = total + T * (run({0: eps, pin: T}, free) - run({0: 0.0, pin: T}, free))
            else:
                total = total + T * run({pin: T}, free)
    else:
        raise ValueError(mode)
    return total


def prefactor(cyc):
    if cyc.special:
        return 0.5 * cyc.sign / cyc.scale
    return cyc.sign / (2.0 * cyc.n * cyc.scale)


# ============================================================== enumeration
def enumerate_cycles(insertions, jet, canon_odd, special, lines_P=PLINES, lines_S=SLINES):
    """all orderings of the given distinct insertions and all line assignments (line 0 special if `special`)
    with nonvanishing adjacent blocks"""
    out = []
    n = len(insertions)
    for perm in itertools.permutations(range(n)):
        ins = [insertions[i] for i in perm]

        def rec(i, acc):
            if i == n:
                if ins[n - 1].nonzero_between(acc[n - 1], acc[0]):
                    yield list(acc)
                return
            choices = (lines_S if (special and i == 0) else lines_P)
            for L in choices:
                if i > 0 and not ins[i - 1].nonzero_between(acc[i - 1], L):
                    continue
                yield from rec(i + 1, acc + [L])
        for Ls in rec(0, []):
            cyc = Cycle(ins, Ls, jet, canon_odd, special)
            out.append((perm, Ls, cyc))
    return out


# ============================================================== tau-analytic integration of flow cycles
from scipy import special as _sps


def _geom(asm, u):
    cyc = asm.cyc
    K0 = u @ cyc.a
    Ub = np.einsum('pi,ivm->pvm', u, cyc.B)
    V0 = u @ (cyc.a ** 2).sum(axis=1) - (K0 ** 2).sum(axis=1)
    V0 = np.maximum(V0, 0.0)
    wv = np.einsum('pi,im,ivm->pv', u, cyc.a, cyc.B) - np.einsum('pm,pvm->pv', K0, Ub)
    lam = np.einsum('pi,ivm,iwm->pvw', u, cyc.B, cyc.B) - np.einsum('pvm,pwm->pvw', Ub, Ub)
    return K0, Ub, V0, wv, lam


def tau_coeffs(asm, u):
    """F(tau u) = tau^{-2} e^{-tau V0(u)} sum_j G[:, j] tau^{j + jmin}; returns (G (nb, nj, npts), jmin, V0)"""
    jet = asm.jet
    npts = u.shape[0]
    nv = jet.nvar
    K0, Ub, V0, wv, lam = _geom(asm, u)
    k0 = -K0
    dk = -Ub
    X1 = np.zeros((jet.nb, npts), complex)
    for v in range(nv):
        if jet.lin[v] is not None:
            X1[jet.lin[v]] += -2 * wv[:, v]
    for v in range(nv):
        for w in range(nv):
            e = tuple((1 if k == v else 0) + (1 if k == w else 0) for k in range(nv))
            if e in jet.index:
                X1[jet.index[e]] += -lam[:, v, w]
    Xp = [np.zeros((jet.nb, npts), complex)]
    Xp[0][0] = 1.0
    for pp in range(1, jet.maxorder + 1):
        Xp.append(jet.mul(Xp[-1], X1) / pp)
    vals = asm.stack.eval(k0)
    fac = [1.0 / (4.0 ** m * math.factorial(m)) for m in range(asm.M + 1)]
    dkj = np.zeros((4, jet.nb, npts), complex)
    for v in range(nv):
        if jet.lin[v] is not None:
            for a in range(4):
                dkj[a, jet.lin[v]] = dk[:, v, a]
    dkdk = {}
    if jet.maxorder >= 2:
        for a in range(4):
            for c in range(a, 4):
                dkdk[(a, c)] = jet.mul(dkj[a], dkj[c])
    Em = [np.zeros((jet.nb, npts), complex) for _ in range(asm.M + 1)]
    one = np.zeros((jet.nb, npts), complex); one[0] = 1.0
    inner = {}
    for col, (b, alpha, m) in enumerate(asm.entries):
        coef = vals[:, col] * fac[m]
        if len(alpha) == 0:
            Dl = one
        elif len(alpha) == 1:
            Dl = dkj[alpha[0]]
        else:
            a, c = alpha
            Dl = dkdk[(a, c)] * (0.5 if a == c else 1.0)
        key = (b, m)
        inner[key] = inner.get(key, 0) + Dl * coef[None]
    for (b, m), X in inner.items():
        Em[m] = Em[m] + jet.shift(b, X)
    jmin = -asm.M
    jmax = jet.maxorder
    G = np.zeros((jet.nb, jmax - jmin + 1, npts), complex)
    for m in range(asm.M + 1):
        for pp in range(jet.maxorder + 1):
            G[:, pp - m - jmin, :] += jet.mul(Em[m], Xp[pp])
    return G, jmin, V0


def incgam(alpha, V, a, b):
    """int_a^b tau^alpha e^{-V tau} d tau  (alpha integer, V >= 0, 0 < a <= b), vectorized"""
    V = np.asarray(V, float); a = np.asarray(a, float); b = np.asarray(b, float)
    out = np.zeros(np.broadcast(V, a, b).shape)
    V, a, b = np.broadcast_arrays(V, a, b)
    small = V * b < 0.05
    if np.any(small):
        Vs, as_, bs = V[small], a[small], b[small]
        acc = np.zeros(Vs.shape)
        ta = np.ones(Vs.shape)            # (-V a)^k / k!
        tb = np.ones(Vs.shape)            # (-V b)^k / k!
        tv = np.ones(Vs.shape)            # (-V)^k / k!  (only used for the log term)
        pa = as_ ** (alpha + 1.0)
        pb = bs ** (alpha + 1.0)
        for k in range(0, 30):
            e = alpha + k + 1
            if e == 0:
                acc = acc + tv * np.log(bs / as_)
            else:
                acc = acc + (pb * tb - pa * ta) / e
            ta = ta * (-Vs * as_) / (k + 1)
            tb = tb * (-Vs * bs) / (k + 1)
            tv = tv * (-Vs) / (k + 1)
        out[small] = acc
    big = ~small
    if np.any(big):
        Vb, ab, bb = V[big], a[big], b[big]
        if alpha >= 0:
            def G(x):
                s_ = np.zeros(x.shape); t_ = np.ones(x.shape)
                for k in range(alpha + 1):
                    s_ = s_ + t_
                    t_ = t_ * x / (k + 1)
                return math.factorial(alpha) * np.exp(-x) * s_
            out[big] = (G(Vb * ab) - G(Vb * bb)) / Vb ** (alpha + 1)
        elif alpha == -1:
            out[big] = _sps.exp1(Vb * ab) - _sps.exp1(Vb * bb)
        else:
            n = -alpha
            out[big] = ab ** (1 - n) * _sps.expn(n, Vb * ab) - bb ** (1 - n) * _sps.expn(n, Vb * bb)
    return out


def _ratio_grid(n, L, per_decade):
    """points rho in R^{n-1}, rho_k <= 0, sum rho_k >= L, with GL weights (panels of one decade in each rho_k,
    nested so that the constraint is respected exactly).  Returns (rho (N, n-1), w (N,))"""
    dec = np.log(10.0)
    x, wq = gl(per_decade)

    def panels(lo, hi):
        if hi - lo <= 1e-14:
            return np.zeros(0), np.zeros(0)
        nd = max(1, int(np.ceil((hi - lo) / dec - 1e-9)))
        edges = np.linspace(lo, hi, nd + 1)
        ts, ws = [], []
        for e0, e1 in zip(edges[:-1], edges[1:]):
            ts.append(0.5 * (e1 - e0) * x + 0.5 * (e1 + e0)); ws.append(0.5 * (e1 - e0) * wq)
        return np.concatenate(ts), np.concatenate(ws)
    if n - 1 == 1:
        r, w = panels(L, 0.0)
        return r[:, None], w
    pts, wts = [], []

    def rec(k, prefix, wpre, used):
        lo = L - used
        r, w = panels(lo, 0.0) if k < n - 2 else panels(L - used, 0.0)
        for ri, wi in zip(r, w):
            if k == n - 2:
                pts.append(prefix + [ri]); wts.append(wpre * wi)
            else:
                rec(k + 1, prefix + [ri], wpre * wi, used + ri)
    rec(0, [], 1.0, 0.0)
    return np.array(pts), np.array(wts)


def integrate_flow_simplex(asm, eps, T, per_decade=10, chunk=20000):
    """flow cycle (all lines P(eps,T)): sum over orderings u_s1 <= ... <= u_sn of
    int drho Jac W^{-n} sum_j G_j I(n-3+j, V0, eps/u_min, T/u_max).  Returns raw jets (nb,)"""
    cyc = asm.cyc
    n = cyc.n
    assert not cyc.special
    L = np.log(eps / T)
    rho, w = _ratio_grid(n, L, per_decade)
    r = np.exp(rho)                                   # (N, n-1): r_k = u_{s_k}/u_{s_{k+1}}
    # w_k = prod_{j>=k} r_j  (k = 1..n-1), w_n = 1
    wk = np.ones((len(w), n))
    for k in range(n - 2, -1, -1):
        wk[:, k] = wk[:, k + 1] * r[:, k]
    W = wk.sum(axis=1)
    jac = np.ones(len(w))
    for j in range(1, n - 1):
        jac = jac * r[:, j] ** j
    jac = jac * np.prod(r, axis=1)                    # d r = r d rho
    meas = w * jac / W ** n
    total = np.zeros(asm.jet.nb, complex)
    for perm in itertools.permutations(range(n)):
        u = np.zeros((len(w), n))
        for k in range(n):
            u[:, perm[k]] = wk[:, k] / W
        umin = u[:, perm[0]]; umax = u[:, perm[n - 1]]
        for s0 in range(0, len(w), chunk):
            sl = slice(s0, s0 + chunk)
            G, jmin, V0 = tau_coeffs(asm, u[sl])
            a_ = eps / umin[sl]; b_ = T / umax[sl]
            ok = b_ > a_
            acc = np.zeros((asm.jet.nb, int(ok.sum())), complex)
            for jj in range(G.shape[1]):
                Iv = incgam(n - 3 + jj + jmin, V0[sl][ok], a_[ok], b_[ok])
                acc = acc + G[:, jj, ok] * Iv[None]
            total = total + acc @ meas[sl][ok]
    return total
