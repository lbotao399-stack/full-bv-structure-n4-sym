"""Engine B -- eps -> 0 asymptotics of Schwinger integrals by Hepp sectors.

Input: an integral over t in [eps, inf)^n of  sum_nu h_nu(t) exp(-phi(t))  (|x|^2 = 1 by scaling), where h_nu is
rational and homogeneous of degree -(n + nu) and phi is rational, positive and homogeneous of degree -1.
Scale integration (t = S tau):  int_{S >= eps/min tau} S^(n-1) S^(-n-nu) exp(-phi(tau)/S) dS
                               = phi(tau)^(-nu) * lowergamma(nu, phi(tau) * min(tau) / eps).
Hepp sector sigma (ordering t_sigma(1) <= ... <= t_sigma(n)):  t_sigma(k) = S * prod_{j>=k} y_j, y in [0,1]^(n-1),
Jacobian S^(n-1) prod_j y_j^(j-1), min(tau) = prod_j y_j.  In each sector every polynomial with positive
coefficients is (monomial) x (polynomial with positive constant term); the integrand becomes
    prod_j y_j^(a_j) * R(y) * lowergamma(nu, Z(y)/eps),   Z = phi * prod_j y_j = y^p * Ztilde(y).
For a_j >= 0 for all j the sector has a limit (lowergamma -> Gamma(nu)).  A sector with a single a_{j*} = -1
contributes  -(Gamma(nu)/p_{j*}) * int R(0,y') dy'  to the coefficient of log eps and
    Gamma(nu) int [R(y) - R(0,y')]/y_{j*} dy + (1/p_{j*}) int R(0,y') [Gamma(nu) log Ztilde(0,y') + C_nu] dy'
to the constant, with C_nu = -Gamma(nu) psi(nu)  (int_0^Y lowergamma(nu,z) dz/z = Gamma(nu) log Y + C_nu + o(1)).
Exact statements: see notes.txt, section H.
"""
import itertools
import math
import numpy as np
import sympy as sp


class Poly:
    """Polynomial in t_1..t_n with exponent matrix E (m x n) and float coefficients c (m)."""

    def __init__(self, expr, tsyms):
        P = sp.Poly(sp.expand(expr), *tsyms)
        terms = P.terms()
        self.E = np.array([list(e) for e, _ in terms], dtype=np.int64)
        self.c = np.array([float(cf) for _, cf in terms], dtype=np.float64)
        self.cexact = [sp.nsimplify(cf) for _, cf in terms]
        self.positive = bool(np.all(self.c > 0))

    def sector(self, sigma, dominance=True):
        """y-exponent matrix in Hepp sector sigma (tuple of edge indices, smallest first)."""
        Es = self.E[:, list(sigma)]
        Y = np.cumsum(Es, axis=1)[:, :-1]  # exponent of y_j = sum_{k<=j} alpha_sigma(k), j=1..n-1
        m = Y.min(axis=0)
        if dominance:
            # dominance: the componentwise minimum must be a monomial of the polynomial
            if not np.any(np.all(Y == m[None, :], axis=1)):
                raise RuntimeError('no dominant monomial in sector %s' % (sigma,))
        return Y - m[None, :], m


def eval_reduced(Yred, c, y):
    """sum_i c_i prod_j y_j^Yred_ij at points y (P x d)."""
    out = np.zeros(y.shape[0])
    for i in range(len(c)):
        term = np.full(y.shape[0], c[i])
        for j in range(Yred.shape[1]):
            e = Yred[i, j]
            if e:
                term = term * y[:, j] ** e
        out += term
    return out


class Term:
    """One nu-term:  const * N(t) / prod_f f(t)^e_f * lowergamma(nu, phi*min/eps), phi = cphi * Pphi/Qphi."""

    def __init__(self, gexpr, nu, phiexpr, tsyms):
        self.nu = nu
        self.n = len(tsyms)
        g = sp.factor(sp.together(gexpr))
        num, den = sp.fraction(g)
        cn, nfac = sp.factor_list(num)
        cd, dfac = sp.factor_list(den)
        self.const = float(sp.nsimplify(cn / cd))
        self.const_exact = sp.nsimplify(cn / cd)
        Nexpr = sp.Integer(1)
        for f, e in nfac:
            Nexpr *= f ** e
        self.N = Poly(Nexpr, tsyms)
        self.dens = [(Poly(f, tsyms), int(e)) for f, e in dfac]
        for P, e in self.dens:
            if not P.positive:
                raise RuntimeError('denominator factor without positive coefficients')
        ph = sp.factor(sp.together(phiexpr))
        pn, pd = sp.fraction(ph)
        cpn, pnf = sp.factor_list(pn)
        cpd, pdf = sp.factor_list(pd)
        self.cphi = float(sp.nsimplify(cpn / cpd))
        self.phinum = [(Poly(f, tsyms), int(e)) for f, e in pnf]
        self.phiden = [(Poly(f, tsyms), int(e)) for f, e in pdf]

    def sector(self, sigma):
        n = self.n
        a = np.arange(n - 1, dtype=np.int64)  # Jacobian prod y_j^(j-1), 0-based j -> j
        NY, mN = self.N.sector(sigma, dominance=False)
        a = a + mN
        dens = []
        for P, e in self.dens:
            Y, m = P.sector(sigma)
            a = a - e * m
            dens.append((Y, P.c, e))
        p = np.ones(n - 1, dtype=np.int64)
        phn, phd = [], []
        for P, e in self.phinum:
            Y, m = P.sector(sigma)
            p = p + e * m
            phn.append((Y, P.c, e))
        for P, e in self.phiden:
            Y, m = P.sector(sigma)
            p = p - e * m
            phd.append((Y, P.c, e))
        return dict(a=a, p=p, NY=NY, Nc=self.N.c, dens=dens, phn=phn, phd=phd)

    def R(self, sd, y, shift=None):
        """const * Ntilde / prod dens~^e * prod_j y_j^(a_j) (with a_{shift} raised by one)."""
        val = self.const * eval_reduced(sd['NY'], sd['Nc'], y)
        for Y, c, e in sd['dens']:
            val = val / eval_reduced(Y, c, y) ** e
        a = sd['a'].copy()
        if shift is not None:
            a[shift] += 1
        for j in range(len(a)):
            if a[j]:
                val = val * y[:, j] ** a[j]
        return val

    def logZt(self, sd, y, jstar):
        """log of Z / y_{j*}^{p_{j*}} (requires y_{j*} = 0 in y for the limit)."""
        val = np.full(y.shape[0], math.log(self.cphi))
        for Y, c, e in sd['phn']:
            val += e * np.log(eval_reduced(Y, c, y))
        for Y, c, e in sd['phd']:
            val -= e * np.log(eval_reduced(Y, c, y))
        for j in range(len(sd['p'])):
            if j != jstar and sd['p'][j]:
                val += sd['p'][j] * np.log(y[:, j])
        return val

    def Z(self, sd, y):
        val = self.cphi * np.ones(y.shape[0])
        for Y, c, e in sd['phn']:
            val = val * eval_reduced(Y, c, y) ** e
        for Y, c, e in sd['phd']:
            val = val / eval_reduced(Y, c, y) ** e
        for j in range(len(sd['p'])):
            if sd['p'][j]:
                val = val * y[:, j] ** sd['p'][j]
        return val


def gl_nodes(N, power=1):
    """Gauss-Legendre on [0,1] in u with y = u^power: returns y-nodes and weights dy."""
    x, w = np.polynomial.legendre.leggauss(N)
    u = 0.5 * (x + 1.0)
    wu = 0.5 * w
    y = u ** power
    wy = wu * power * u ** (power - 1)
    return y, wy


def tensor_grid(d, N, power, chunk=None):
    y1, w1 = gl_nodes(N, power)
    grids = np.meshgrid(*([y1] * d), indexing='ij')
    wg = np.meshgrid(*([w1] * d), indexing='ij')
    Y = np.stack([g.ravel() for g in grids], axis=1)
    W = np.prod(np.stack([g.ravel() for g in wg], axis=1), axis=1)
    return Y, W


def C_nu(nu):
    return -math.gamma(nu) * float(sp.polygamma(0, nu))


def C_nu_parts(nu):
    """C_nu = -Gamma(nu) psi(nu) = Gamma(nu) (EulerGamma - H_{nu-1}):  (rational part, coefficient of EulerGamma)."""
    H = sum(1.0 / k for k in range(1, nu))
    return -math.gamma(nu) * H, math.gamma(nu)


def asymptotics(terms, n, N=20, power=2, sectors=None, verbose=False):
    """Return (a, b_rat, b_gamma, info):  I(eps) = a log eps + b_rat + b_gamma * EulerGamma + o(1).
    (b_gamma = -a identically, since the EulerGamma part of C_nu is Gamma(nu) EulerGamma.)"""
    d = n - 1
    Yg, Wg = tensor_grid(d, N, power)
    Yh, Wh = tensor_grid(d - 1, N, power) if d > 1 else (np.zeros((1, 0)), np.ones(1))
    a_tot = 0.0
    b_rat = 0.0
    b_gam = 0.0
    info = {'divergent_sector_terms': 0, 'sectors': 0, 'min_a': 99, 'p_jstar': set(), 'divergent_sets': set()}
    perms = list(itertools.permutations(range(n))) if sectors is None else sectors
    for sigma in perms:
        info['sectors'] += 1
        for T in terms:
            sd = T.sector(sigma)
            a = sd['a']
            info['min_a'] = min(info['min_a'], int(a.min()))
            neg = np.where(a < 0)[0]
            G = math.gamma(T.nu)
            if len(neg) == 0:
                b_rat += G * np.dot(Wg, T.R(sd, Yg))
                continue
            if len(neg) > 1 or a[neg[0]] < -1:
                raise RuntimeError('sector %s: exponents %s not handled' % (sigma, a))
            js = int(neg[0])
            info['divergent_sector_terms'] += 1
            info['divergent_sets'].add(tuple(sorted(sigma[:js + 1])))
            ps = sd['p'][js]
            info['p_jstar'].add(int(ps))
            if ps <= 0:
                raise RuntimeError('sector %s: p_j* = %d' % (sigma, ps))
            Rfull = T.R(sd, Yg, shift=js)
            Y0 = Yg.copy()
            Y0[:, js] = 0.0
            R0 = T.R(sd, Y0, shift=js)
            b_rat += G * np.dot(Wg, (Rfull - R0) / Yg[:, js])
            Yf = np.insert(Yh, js, 0.0, axis=1)
            Rf = T.R(sd, Yf, shift=js)
            Lf = T.logZt(sd, Yf, js)
            cr, cg = C_nu_parts(T.nu)
            IR0 = np.dot(Wh, Rf)
            a_tot += -(G / ps) * IR0
            b_rat += (1.0 / ps) * (G * np.dot(Wh, Rf * Lf) + cr * IR0)
            b_gam += (1.0 / ps) * cg * IR0
    info['p_jstar'] = sorted(info['p_jstar'])
    info['divergent_sets'] = sorted(info['divergent_sets'])
    return a_tot, b_rat, b_gam, info


def lower_gamma(nu, z):
    """Lower incomplete gamma for nu = 2, 3 (vectorised, series for small z)."""
    z = np.asarray(z, dtype=np.float64)
    out = np.empty_like(z)
    small = z < 0.5
    zs = z[small]
    s = np.zeros_like(zs)
    term = np.ones_like(zs)  # (-z)^k / k!
    for k in range(16):
        s += term * zs ** nu / (k + nu)
        term = term * (-zs) / (k + 1)
    out[small] = s
    zl = z[~small]
    e = np.exp(-np.minimum(zl, 700.0))
    if nu == 2:
        out[~small] = 1.0 - (1.0 + zl) * e
    elif nu == 3:
        out[~small] = 2.0 - (2.0 + 2.0 * zl + zl * zl) * e
    else:
        raise NotImplementedError
    return out


def graded_nodes(nodes_per_panel=4, levels=8, power=5):
    """Composite Gauss-Legendre in u on panels [0,2^-levels], [2^-levels, 2^-levels+1], ..., [1/2, 1], y = u^power."""
    xg, wg = np.polynomial.legendre.leggauss(nodes_per_panel)
    edges = [0.0] + [2.0 ** (-k) for k in range(levels, -1, -1)]
    us, ws = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        us.append(a + (b - a) * 0.5 * (xg + 1.0))
        ws.append((b - a) * 0.5 * wg)
    u = np.concatenate(us)
    w = np.concatenate(ws)
    return u ** power, w * power * u ** (power - 1)


def graded_grid(d, **kw):
    y1, w1 = graded_nodes(**kw)
    grids = np.meshgrid(*([y1] * d), indexing='ij')
    wg = np.meshgrid(*([w1] * d), indexing='ij')
    Y = np.stack([g.ravel() for g in grids], axis=1)
    W = np.prod(np.stack([g.ravel() for g in wg], axis=1), axis=1)
    return Y, W


def _accumulate(res, eps_list, nu, W, R, Z, zcut=60.0):
    """res[i] += sum W R lowergamma(nu, Z/eps_i), using lowergamma = Gamma(nu) for Z/eps > zcut."""
    G = math.gamma(nu)
    base = G * np.dot(W, R)
    order = np.argsort(Z)
    Zs = Z[order]
    WRs = (W * R)[order]
    for i, e in enumerate(eps_list):
        k = np.searchsorted(Zs, zcut * e)
        corr = 0.0
        if k > 0:
            corr = np.dot(WRs[:k], G - lower_gamma(nu, Zs[:k] / e))
        res[i] += base - corr


def finite_eps(terms, n, eps_list, nodes_per_panel=4, levels=8, power=5, panel_nodes=8, extra=40.0,
               chunk=400000, plain_N=None):
    """I(eps) for each eps in eps_list (same Hepp representation, no asymptotic expansion).
    Divergent direction y_{j*} = exp(-u), u in [0, log(1/eps_min) + extra], composite Gauss-Legendre panels of
    width 1; other directions graded composite Gauss-Legendre (graded_nodes)."""
    eps_list = np.asarray(eps_list, dtype=np.float64)
    d = n - 1
    if plain_N is None:
        Yg, Wg = graded_grid(d, nodes_per_panel=nodes_per_panel, levels=levels, power=power)
        Yh, Wh = graded_grid(d - 1, nodes_per_panel=nodes_per_panel, levels=levels, power=power)
    else:
        Yg, Wg = tensor_grid(d, plain_N, power)
        Yh, Wh = tensor_grid(d - 1, plain_N, power)
    Umax = math.log(1.0 / eps_list.min()) + extra
    npan = int(math.ceil(Umax))
    xg, wg = np.polynomial.legendre.leggauss(panel_nodes)
    uu = np.concatenate([k + 0.5 * (xg + 1.0) for k in range(npan)])
    wu = np.concatenate([0.5 * wg for k in range(npan)])
    res = np.zeros(len(eps_list))
    for sigma in itertools.permutations(range(n)):
        for T in terms:
            sd = T.sector(sigma)
            a = sd['a']
            neg = np.where(a < 0)[0]
            if len(neg) == 0:
                for c0 in range(0, Yg.shape[0], chunk):
                    Y = Yg[c0:c0 + chunk]
                    _accumulate(res, eps_list, T.nu, Wg[c0:c0 + chunk], T.R(sd, Y), T.Z(sd, Y))
                continue
            js = int(neg[0])
            P = Yh.shape[0]
            step = max(1, chunk // P)
            for k0 in range(0, len(uu), step):
                ub = uu[k0:k0 + step]
                wb = wu[k0:k0 + step]
                Y = np.repeat(Yh, len(ub), axis=0)
                Y = np.insert(Y, js, np.tile(np.exp(-ub), P), axis=1)
                W = np.repeat(Wh, len(ub)) * np.tile(wb, P)
                _accumulate(res, eps_list, T.nu, W, T.R(sd, Y, shift=js), T.Z(sd, Y))
    return res


def fit_log(eps_list, vals, nterms=6, logsq=False):
    """Least-squares fit vals = a log eps + b + eps (c log eps + d) + eps^2 (e log eps + f) + ..."""
    L = np.log(eps_list)
    cols = [L, np.ones_like(L)]
    k = 1
    while len(cols) < nterms:
        cols.append(eps_list ** k * L)
        if len(cols) < nterms:
            cols.append(eps_list ** k)
        if logsq and len(cols) < nterms:
            cols.append(eps_list ** k * L * L)
        k += 1
    A = np.stack(cols, axis=1)
    coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
    return coef
