"""Realisation-parametrised tree-level operations on plane-wave densities with position dependence (both routes):
a functional with position degree <= 1 is PFn(w0, [w_mu]) = int(w0 + x^mu w_mu); its canonical data on a configuration
legs_at(r) (last leg momentum -sum p - r) are alpha = w0 + sum_mu (-i) d/dt w_mu(r = t e_mu)|_0 and
beta_mu = w_mu(r = 0) (Def c11-canonical; the same rule as the injected-momentum evaluation of the cycles, eb_pos).
Q (Prop c1-Q symbol of the route's realisation, Koszul rule of eb_qb.apply_Q) and the bracket {G, F} with a
translation-invariant quadratic G (derivation M_G(p) = -K0 H_G(p,-p)^T, Pi for odd G; cb_tree, gate GB-tree) act
leg by leg on the plane waves, so they commute with the momentum injection."""
import numpy as np
from fractions import Fraction as Fr
import fitlib as FL

OFF = {'A': 0, 'lam': 10, 'c': 26, 'As': 27, 'lams': 37, 'cs': 53}
DIM = {'A': 10, 'lam': 16, 'c': 1, 'As': 10, 'lams': 16, 'cs': 1}
PAR = {'A': 0, 'lam': 1, 'c': 1, 'As': 1, 'lams': 0, 'cs': 0}
SPECIES = ['A', 'lam', 'c', 'As', 'lams', 'cs']
SPEC_OF = [s for s in SPECIES for _ in range(DIM[s])]
PARV = np.array([PAR[s] for s in SPEC_OF])
PI = np.diag([(-1) ** int(x) for x in PARV]).astype(complex)
QT = {'c': 'A', 'A': 'As', 'As': 'cs', 'lam': 'lams'}
TV = (0, 1, -1, 2, -2)


def lagrange_deriv0(ts):
    ts = [Fr(t) for t in ts]
    ws = []
    for i, ti in enumerate(ts):
        others = [tj for j, tj in enumerate(ts) if j != i]
        den = Fr(1)
        for tj in others:
            den *= (ti - tj)
        tot = Fr(0)
        for k in range(len(others)):
            pr = Fr(1)
            for jj, tj in enumerate(others):
                if jj != k:
                    pr *= (-tj)
            tot += pr
        ws.append(tot / den)
    return ws


WD = lagrange_deriv0(TV)


class PFn:
    def __init__(self, w0=None, wmu=None):
        self.w0, self.wmu = w0, wmu


def canon(fn, legs_at):
    """canonical data [alpha, beta_0..beta_3] as Gaussian rational pairs."""
    zero = [Fr(0)] * 4
    ar, ai = Fr(0), Fr(0)
    if fn.w0 is not None:
        a = FL.to_fr(fn.w0(legs_at(zero)))
        ar += a[0]; ai += a[1]
    betas = [(Fr(0), Fr(0))] * 4
    if fn.wmu is not None:
        betas = []
        for mu in range(4):
            bz = (Fr(0), Fr(0))
            for t, w in zip(TV, WD):
                r = [Fr(t) if m == mu else Fr(0) for m in range(4)]
                v = FL.to_fr(fn.wmu[mu](legs_at(r)))
                ar += w * v[1]; ai += -w * v[0]
                if t == 0:
                    bz = v
            betas.append(bz)
    return [(ar, ai)] + betas


class Tree:
    def __init__(self, qmat, K0):
        self.qmat, self.K0 = qmat, K0
        self._hp = {}

    def qleg(self, leg):
        t = leg['sp']
        if t not in QT:
            return None
        s = QT[t]
        v = np.zeros(54, complex)
        v[OFF[t]:OFF[t] + DIM[t]] = leg['u']
        w = self.qmat(leg['p']) @ v
        u = w[OFF[s]:OFF[s] + DIM[s]].copy()
        return None if not np.any(u) else dict(sp=s, u=u, p=leg['p'])

    def Qd(self, dens):
        def f(legs):
            tot = 0j
            for j, leg in enumerate(legs):
                ql = self.qleg(leg)
                if ql is None:
                    continue
                nb = sum(PAR[l['sp']] for l in legs[:j])
                new = list(legs); new[j] = ql
                tot += (-1) ** nb * dens(new)
            return tot
        return f

    def Q(self, fn):
        return PFn(None if fn.w0 is None else self.Qd(fn.w0), None if fn.wmu is None else [self.Qd(w) for w in fn.wmu])

    def hessian(self, g, kr, kc):
        H = np.zeros((54, 54), complex)
        for J in range(54):
            sJ = SPEC_OF[J]; uJ = np.zeros(DIM[sJ], complex); uJ[J - OFF[sJ]] = 1
            lJ = dict(sp=sJ, u=uJ, p=list(kr))
            for K in range(54):
                sK = SPEC_OF[K]; uK = np.zeros(DIM[sK], complex); uK[K - OFF[sK]] = 1
                v = g([lJ, dict(sp=sK, u=uK, p=list(kc))])
                if v != 0:
                    H[J, K] = (-1) ** (PARV[J] * PARV[K]) * v
        return H

    def hpoly(self, g):
        """H_g(k,-k) as an exact quadratic polynomial in k (15 evaluations; checked at a random k)."""
        key = id(g)
        if key in self._hp:
            return self._hp[key][1]
        H = lambda k: self.hessian(g, k, [-x for x in k])
        E4 = np.eye(4, dtype=int)
        H0 = H([0, 0, 0, 0])
        Hp = [H(list(E4[m])) for m in range(4)]; Hm = [H(list(-E4[m])) for m in range(4)]
        a = [(Hp[m] - Hm[m]) / 2 for m in range(4)]
        b = {(m, m): (Hp[m] + Hm[m]) / 2 - H0 for m in range(4)}
        for m in range(4):
            for n in range(m + 1, 4):
                b[(m, n)] = H(list(E4[m] + E4[n])) - H0 - a[m] - a[n] - b[(m, m)] - b[(n, n)]

        def ev(k):
            k = [complex(x) for x in k]
            out = H0 + sum(k[m] * a[m] for m in range(4))
            for (m, n), B_ in b.items():
                out = out + k[m] * k[n] * B_
            return out
        kt = [2, -1, 3, 1]
        if not np.allclose(ev(kt), H(kt)):
            memo = {}

            def ev(k):
                kk = tuple(complex(x) for x in k)
                if kk not in memo:
                    memo[kk] = H(list(k))
                return memo[kk]
        self._hp[key] = (g, ev)
        return ev

    def brd(self, g, gpar, f):
        cache = {}
        hp = self.hpoly(g)

        def h(legs):
            tot = 0j
            for j, leg in enumerate(legs):
                key = tuple(complex(x) for x in leg['p'])
                if key not in cache:
                    M = -self.K0 @ hp(leg['p']).T
                    cache[key] = M if gpar == 0 else PI @ M
                v = np.zeros(54, complex); v[OFF[leg['sp']]:OFF[leg['sp']] + DIM[leg['sp']]] = leg['u']
                img = cache[key] @ v
                nb = sum(PAR[l['sp']] for l in legs[:j])
                sg = (-1) ** (nb * (1 - gpar))
                for s in SPECIES:
                    u = img[OFF[s]:OFF[s] + DIM[s]]
                    if np.any(u):
                        new = list(legs); new[j] = dict(sp=s, u=u.copy(), p=leg['p'])
                        tot += sg * f(new)
            return tot
        return h

    def br(self, g, gpar, fn):
        """{G, F} for translation-invariant quadratic G (density g, parity gpar) and position-dependent F (PFn)."""
        return PFn(None if fn.w0 is None else self.brd(g, gpar, fn.w0),
                   None if fn.wmu is None else [self.brd(g, gpar, w) for w in fn.wmu])


def k0_from_pairing(MP):
    """K_0 = -(fibre pairing)^{-1} on the 54-slot space (engine B ce_core.k0_matrix with its pairing matrix)."""
    Mpinv = np.linalg.inv(MP)
    M = np.zeros((54, 54), complex)
    for m in range(10):
        M[OFF['A'] + m, OFF['As'] + m] = -1; M[OFF['As'] + m, OFF['A'] + m] = -1
    for i in range(16):
        for j in range(16):
            M[OFF['lam'] + i, OFF['lams'] + j] = -Mpinv[j, i]
            M[OFF['lams'] + i, OFF['lam'] + j] = -Mpinv[i, j]
    M[OFF['c'], OFF['cs']] = -1; M[OFF['cs'], OFF['c']] = -1
    return M


def psum(*terms):
    """linear combination of PFn's: terms = (coef, PFn)."""
    def comb(getter):
        fs = [(c, getter(f)) for c, f in terms if getter(f) is not None]
        if not fs:
            return None
        return lambda legs: sum((c * g(legs) for c, g in fs), 0j)
    w0 = comb(lambda f: f.w0)
    has = [f for c, f in terms if f.wmu is not None]
    wmu = None
    if has:
        wmu = []
        for mu in range(4):
            fs = [(c, f.wmu[mu]) for c, f in terms if f.wmu is not None]
            wmu.append(lambda legs, fs=fs: sum((c * g(legs) for c, g in fs), 0j))
    return PFn(w0, wmu)


# ------------------------------------------------------------------ position degree two
PAIRS2 = [(m, n) for m in range(4) for n in range(m, 4)]
TV5 = (0, 1, -1, 2, -2)
W1 = lagrange_deriv0(TV5)


def lagrange_deriv2_0(ts):
    """weights for f''(0) exactly for polynomials of degree < len(ts)."""
    import sympy as sp
    t = sp.symbols('t')
    ws = []
    for i, ti in enumerate(ts):
        Lp = 1
        for j, tj in enumerate(ts):
            if j != i:
                Lp *= (t - tj) / sp.Integer(ti - tj)
        ws.append(Fr(str(sp.diff(Lp, t, 2).subs(t, 0))))
    return ws


W2 = lagrange_deriv2_0(TV5)


class PFn2:
    """int(w0 + x^mu w_mu + sum_{mu <= nu} x^mu x^nu w_{mu nu}); wmn: dict (mu,nu) mu <= nu -> density."""

    def __init__(self, w0=None, wmu=None, wmn=None):
        self.w0, self.wmu, self.wmn = w0, wmu, wmn


def _ex(v):
    return (Fr(v[0]), Fr(v[1])) if isinstance(v, tuple) else FL.to_fr(v)


def _dir_vals(f, legs_at, d):
    """values of f(legs_at(t d)) for t in TV5 (Gaussian pairs; f may return exact pairs)."""
    return [_ex(f(legs_at([Fr(t) * x for x in d]))) for t in TV5]


def canon2(fn, legs_at):
    """canonical data [alpha, beta_0..3, gamma_{mu<=nu} (10)] of Phi = alpha delta + (-i) beta^mu d_mu delta
    + (-i)^2 gamma^{mu nu} d_mu d_nu delta (sum over ordered mu, nu with gamma symmetric; stored for mu <= nu as
    the coefficient of x^mu x^nu), derived from the distribution identity g(P) d d delta = ...:
       alpha = f0 - i sum d_mu f_mu - sum_{mu,nu} d_mu d_nu f_{mu nu},  beta_r = f_r - 2i sum_mu d_mu f_{mu r},
       gamma_{mu nu} = w_{mu nu}(0),   f_{mu mu} = w_{mu mu}, f_{mu nu} = f_{nu mu} = w_{mu nu}/2 (mu < nu)."""
    E4 = [[Fr(int(i == m)) for i in range(4)] for m in range(4)]
    z = (Fr(0), Fr(0))
    def add(a, b, c=(Fr(1), Fr(0))):
        return (a[0] + c[0] * b[0] - c[1] * b[1], a[1] + c[0] * b[1] + c[1] * b[0])
    alpha = z; beta = [z] * 4; gamma = [z] * 10
    zero = [Fr(0)] * 4
    if fn.w0 is not None:
        alpha = add(alpha, _ex(fn.w0(legs_at(zero))))
    if fn.wmu is not None:
        for mu in range(4):
            v = _dir_vals(fn.wmu[mu], legs_at, E4[mu])
            d1 = (sum(w * x[0] for w, x in zip(W1, v)), sum(w * x[1] for w, x in zip(W1, v)))
            alpha = add(alpha, d1, (Fr(0), Fr(-1)))
            beta[mu] = add(beta[mu], v[0])
    if fn.wmn is not None:
        for k, (m, n) in enumerate(PAIRS2):
            if (m, n) not in fn.wmn:
                continue
            f = fn.wmn[(m, n)]
            half = Fr(1) if m == n else Fr(1, 2)
            vm = _dir_vals(f, legs_at, E4[m])
            gamma[k] = add(gamma[k], vm[0])
            d1m = (sum(w * x[0] for w, x in zip(W1, vm)), sum(w * x[1] for w, x in zip(W1, vm)))
            if m == n:
                # beta_m += -2i d_m f_mm ; alpha += - d_m^2 f_mm
                beta[m] = add(beta[m], d1m, (Fr(0), Fr(-2)))
                d2 = (sum(w * x[0] for w, x in zip(W2, vm)), sum(w * x[1] for w, x in zip(W2, vm)))
                alpha = add(alpha, d2, (Fr(-1), Fr(0)))
            else:
                vn = _dir_vals(f, legs_at, E4[n])
                d1n = (sum(w * x[0] for w, x in zip(W1, vn)), sum(w * x[1] for w, x in zip(W1, vn)))
                # f_mn = f_nm = w/2: beta_n += -2i d_m (w/2) ; beta_m += -2i d_n (w/2)
                beta[n] = add(beta[n], d1m, (Fr(0), Fr(-1)))
                beta[m] = add(beta[m], d1n, (Fr(0), Fr(-1)))
                # alpha += -(d_m d_n + d_n d_m)(w/2) = - d_m d_n w ; d_m d_n = (D2_{m+n} - D2_m - D2_n)/2
                dmn = [x + y for x, y in zip(E4[m], E4[n])]
                vs = _dir_vals(f, legs_at, dmn)
                D2 = lambda vv: (sum(w * x[0] for w, x in zip(W2, vv)), sum(w * x[1] for w, x in zip(W2, vv)))
                a, b, c = D2(vs), D2(vm), D2(vn)
                mixed = ((a[0] - b[0] - c[0]) / 2, (a[1] - b[1] - c[1]) / 2)
                alpha = add(alpha, mixed, (Fr(-1), Fr(0)))
    return [alpha] + beta + gamma


def tree_apply2(op, fn):
    mp = lambda d: None if d is None else op(d)
    return PFn2(mp(fn.w0), None if fn.wmu is None else [op(w) for w in fn.wmu],
                None if fn.wmn is None else {k: op(w) for k, w in fn.wmn.items()})
