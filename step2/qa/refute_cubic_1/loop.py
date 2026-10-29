"""One-loop words, Schwinger integrals at finite eps, and assembly of the order-hbar part of
W(P(eps,infty), I - I^CT[eps]) on plane-wave legs (refute_cubic_1; independent code).

Formulas (checked on the toy model, toy_conventions.py):
  W_1^{1PI}(e) = sum_n 1/(2n) STr(M^n),  M = P Pi H(e),  H_{bc} = (d_b d_c I)(e)
  order-hbar W(e) = [W_1^{1PI} - I^CT](e + chi(e)),  chi^a = sum_b P^{ab} (d_b I_3)(e)
Edge factor: P(eps,infty)(q) = Phat(q) * int_eps^infty dt e^{-t |q|^2}.
"""
import itertools
import warnings
import numpy as np
from scipy.special import exp1

import ft

warnings.filterwarnings('ignore')

INV16PI2 = 1.0 / (16 * np.pi ** 2)


class Insertion:
    """a background insertion into the loop: cubic (one leg) or quartic (two legs)."""

    def __init__(self, legs):
        self.legs = legs
        self.p = sum(L.p for L in legs)
        terms = ft.CUBIC if len(legs) == 1 else ft.QUARTIC
        z = np.zeros(4)
        H0, m0 = ft.hessian(terms, legs, z, z)
        self.mono = m0
        self.H0 = H0
        self.Hb = []
        self.Hc = []
        for mu in range(4):
            e = np.zeros(4)
            e[mu] = 1.0
            Hb, mb = ft.hessian(terms, legs, e, z)
            Hc, mc = ft.hessian(terms, legs, z, e)
            self.Hb.append(Hb - H0)
            self.Hc.append(Hc - H0)
        self.zero = (m0 is None) or (not np.any(H0) and not any(np.any(x) for x in self.Hb + self.Hc))
        # linearity check
        kb = np.array([0.3, -0.7, 0.2, 0.5])
        kc = np.array([-0.4, 0.1, 0.9, -0.6])
        Hd, _ = ft.hessian(terms, legs, kb, kc)
        if not self.zero:
            assert np.max(abs(Hd - self.H(kb, kc))) < 1e-10, "Hessian not linear in momenta"

    def H(self, kb, kc):
        out = self.H0.copy()
        for mu in range(4):
            if kb[mu]:
                out = out + kb[mu] * self.Hb[mu]
            if kc[mu]:
                out = out + kc[mu] * self.Hc[mu]
        return out

    def m(self, kin):
        q = kin + self.p
        return ft.Phat(q) @ ft.PI @ self.H(-q, kin)


def words_for(legs):
    """all ordered words of insertions whose leg sets partition the legs, with weight 1/(2n)."""
    L = len(legs)
    idx = list(range(L))
    blocks = []
    # set partitions of idx into blocks of size 1 or 2
    def parts(rem):
        if not rem:
            yield []
            return
        a = rem[0]
        for rest in parts(rem[1:]):
            yield [(a,)] + rest
        for j in range(1, len(rem)):
            b = rem[j]
            r2 = rem[1:j] + rem[j + 1:]
            for rest in parts(r2):
                yield [(a, b)] + rest
    out = []
    for P in parts(idx):
        for order in itertools.permutations(P):
            out.append(list(order))
    return out


class OneLoop:
    """precomputed polynomial numerators for the 1PI one-loop multilinear part on legs"""

    def __init__(self, legs, rng):
        self.legs = legs
        self.items = []   # (monomial, weight*sign, n, offsets s_e, coeffs, fiterr)
        cache = {}
        for word in words_for(legs):
            ins = []
            for blk in word:
                if blk not in cache:
                    cache[blk] = Insertion([legs[i] for i in blk])
                ins.append(cache[blk])
            if any(x.zero for x in ins):
                continue
            n = len(ins)
            # monomial: matrix order = last applied leftmost
            seq = []
            for x in reversed(ins):
                seq += list(x.mono)
            sg, mono = ft.canon(seq)
            if sg == 0:
                continue
            ops = [(x.m, x.p) for x in ins]
            c, err = ft.fit_poly(ops, rng, maxdeg=n)   # degree <= number of edges (see header)
            offs = []
            acc = np.zeros(4)
            for x in ins:
                acc = acc + x.p
                offs.append(acc.copy())
            self.items.append((mono, sg / (2.0 * n), n, np.array(offs), c, ft.lap_coeffs(c), err))

    def max_fit_err(self):
        # relative to the largest numerator of this set of words (identically vanishing words have noise only)
        scale = max([np.abs(it[4]).max() for it in self.items] + [1e-300])
        return max([it[6] * np.abs(it[4]).max() / scale for it in self.items] + [0.0])

    def value(self, eps, quad):
        res = {}
        for mono, w, n, offs, c, cl, err in self.items:
            v = w * quad.integrate(n, offs, c, cl, eps)
            res[mono] = res.get(mono, 0) + v
        return res


def polyval(c, K):
    return ft.poly_basis(K) @ c


class Quad:
    """Schwinger integrals at finite eps after the Gaussian loop integral:
    int d^4k/(2pi)^4 N(k) prod_e int_eps^infty dt_e exp(-t_e |k + s_e|^2)."""

    def __init__(self, nt=80, ny=(60, 80, 60), h=None):
        self.nt = nt
        self.ny = ny
        # tanh-sinh nodes on (0,1)
        hstep = 6.0 / nt if h is None else h
        tt = np.arange(-nt, nt + 1) * hstep
        x = 0.5 * (1 + np.tanh(0.5 * np.pi * np.sinh(tt)))
        w = 0.5 * 0.5 * np.pi * np.cosh(tt) / np.cosh(0.5 * np.pi * np.sinh(tt)) ** 2 * hstep
        keep = (x > 0) & (x < 1) & (w > 0)
        self.ts_x, self.ts_w = x[keep], w[keep]
        # 1-x computed accurately
        self.ts_1mx = 0.5 * (1 - np.tanh(0.5 * np.pi * np.sinh(tt)))[keep]

    def ynodes(self, eps, umax):
        """nodes for u = umax*exp(-y): returns u, weight (including du = u dy)"""
        ye = np.log(umax / eps)
        segs = [(0.0, max(ye - 6, 0.5)), (max(ye - 6, 0.5), ye + 6), (ye + 6, ye + 45)]
        us, ws = [], []
        for (a, b), n in zip(segs, self.ny):
            g, gw = np.polynomial.legendre.leggauss(n)
            y = 0.5 * (b - a) * g + 0.5 * (b + a)
            wy = 0.5 * (b - a) * gw
            u = umax * np.exp(-y)
            us.append(u)
            ws.append(wy * u)
        return np.concatenate(us), np.concatenate(ws)

    def integrate(self, n, offs, c, cl, eps):
        if n == 1:
            # single edge with zero offset (seagull): int dt t^-2 (Nc + N2/t), Gaussian 1/(16 pi^2 t^2)
            assert np.allclose(offs[-1], 0)
            Nc = polyval(c, np.zeros((1, 4)))[0]
            N2 = 0.25 * polyval(cl, np.zeros((1, 4)))[0]
            return INV16PI2 * (Nc / eps + N2 / (2 * eps ** 2))
        if n == 2:
            s1 = offs[0]
            q2 = s1 @ s1
            u, wu = self.ynodes(eps, 0.5)
            tot = 0.0
            for x, w in ((u, wu), (1 - u, wu)):
                kcen = -(x[:, None] * s1[None, :])  # k_c = -sum x_e s_e, s_2 = 0
                D = u * (1 - u) * q2               # = x(1-x)|s_1|^2, computed from the small variable u
                T0 = eps / u
                Nc = polyval(c, kcen)
                N2 = 0.25 * polyval(cl, kcen)
                z = T0 * D
                f = Nc * exp1(z) + N2 * (np.exp(-z) / T0 - D * exp1(z))
                tot = tot + np.sum(w * f)
            return INV16PI2 * tot
        if n == 3:
            s = offs  # s_1, s_2, s_3 (= 0)
            u, wu = self.ynodes(eps, 1.0 / 3.0)
            t, wt = self.ts_x, self.ts_w
            omt = self.ts_1mx
            tot = 0.0
            for imin in range(3):
                others = [j for j in range(3) if j != imin]
                U = u[:, None]
                X = np.zeros((len(u), len(t), 3))
                X[:, :, imin] = U
                X[:, :, others[0]] = U + (1 - 3 * U) * t[None, :]
                X[:, :, others[1]] = U + (1 - 3 * U) * omt[None, :]
                W = (wu[:, None] * wt[None, :]) * (1 - 3 * U)
                Xf = X.reshape(-1, 3)
                kcen = -(Xf @ s)
                # D = sum_{i<j} x_i x_j |s_i - s_j|^2  (no cancellation near the corners)
                D = (Xf[:, 0] * Xf[:, 1] * ((s[0] - s[1]) @ (s[0] - s[1])) + Xf[:, 0] * Xf[:, 2] * ((s[0] - s[2]) @ (s[0] - s[2]))
                     + Xf[:, 1] * Xf[:, 2] * ((s[1] - s[2]) @ (s[1] - s[2])))
                T0 = eps / Xf[:, imin]
                Nc = polyval(c, kcen)
                N2 = 0.25 * polyval(cl, kcen)
                z = T0 * D
                f = Nc * np.exp(-z) / D + N2 * exp1(z)
                tot = tot + np.sum(W.reshape(-1) * f)
            return INV16PI2 * tot
        raise ValueError


# ---------------------------------------------------------------- local functionals as plane-wave values
def local(T2f, slots, legs, coef=1.0):
    return {m: coef * v for m, v in ft.multilinear([(slots, T2f)], legs).items()}


def addto(d, e, s=1.0):
    for k, v in e.items():
        d[k] = d.get(k, 0) + s * v
    return d


def ICT2(legs, eps, mass=9.0 / 4.0, lA=-1.0, ll=4.0, lc=-1.0):
    out = {}
    addto(out, local(ft.T2_AA, ('A', 'A'), legs, ft.UNIT * mass / eps))
    addto(out, local(ft.T2_dAdA, ('A', 'A'), legs, ft.UNIT * lA * np.log(eps)))
    addto(out, local(ft.T2_ldl, ('l', 'l'), legs, ft.UNIT * ll * np.log(eps)))
    addto(out, local(ft.T2_dcAs, ('c', 'As'), legs, ft.UNIT * lc * np.log(eps)))
    return out


def F2long(legs, coef=27.0 / 16.0):
    return local(ft.T2_dsAdsA, ('A', 'A'), legs, ft.UNIT * coef)


MULT3 = {0: 6.0, 1: 10.0, 2: 1.0, 3: 1.0, 4: 1.0}   # claimed log multiples of the cubic vertices (thm c4-log)


def ICT3(legs, eps, mult=MULT3):
    out = {}
    for i, (slots, Tf) in enumerate(ft.CUBIC):
        addto(out, {m: ft.UNIT * mult[i] * np.log(eps) * v for m, v in ft.multilinear([(slots, Tf)], legs).items()})
    return out


class ThreeLeg:
    """order-hbar three-leg part of W(P(eps,infty), I - I^CT - F) on one leg set."""

    def __init__(self, legs, rng):
        self.legs = legs
        self.ol3 = OneLoop(legs, rng)
        self.pairs = []
        for (j, m) in [(0, 1), (0, 2), (1, 2)]:
            l = [i for i in range(3) if i not in (j, m)][0]
            chis = ft.tree_field([legs[j], legs[m]], 0.0)   # fibre without e^{-eps q^2}; times exp below
            q = legs[j].p + legs[m].p
            ols = [(OneLoop([legs[l], ch], rng), ch) for ch in chis]
            self.pairs.append((q, legs[l], ols))

    def components(self, eps, quad, with_F=False):
        comp = {}
        comp['1PI'] = self.ol3.value(eps, quad)
        se, ct2, f2 = {}, {}, {}
        for q, Ll, ols in self.pairs:
            sc = np.exp(-eps * q @ q)
            for ol, ch in ols:
                addto(se, ol.value(eps, quad), sc)
                addto(ct2, ICT2([Ll, ch], eps), -sc)
                if with_F:
                    addto(f2, F2long([Ll, ch]), -sc)
        comp['1PR_selfenergy'] = se
        comp['tree_ICT2'] = ct2
        comp['tree_F2'] = f2
        comp['ICT3'] = {m: -v for m, v in ICT3(self.legs, eps).items()}
        return comp

    def max_fit_err(self):
        e = self.ol3.max_fit_err()
        for q, Ll, ols in self.pairs:
            for ol, ch in ols:
                e = max(e, ol.max_fit_err())
        return e
