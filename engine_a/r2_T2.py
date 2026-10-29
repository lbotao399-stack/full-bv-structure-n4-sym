"""T2 (round 2, engine A): one-loop QME defect at the UV boundary, quadratic sector (external c, A|_V).

Route (a) [graphs with one D_eps vertex]: D^(0) = (1/2)({I,I}_eps - {I,I}_0) = (1/2){I,I}_{K'}, K' = K_eps - K_0;
its one-loop self-contraction is (verified in 0 dimensions, r2_verify_special0d.py)
    d_P (1/2){I,I}_{K'} |_{loop} = (1/2) str(K' S H P S H),
K'(k) = (e^{-eps k^2} - 1) Khat = -int_0^eps dtau k^2 e^{-tau k^2} Khat (special edge, tau in [0,eps]),
P(k) = int_eps^infty dt e^{-t k^2} Phat(k).  D^(1) = Delta_eps I - Q I^CT[eps] - {I, I^CT[eps]}_eps; in this sector only
-Q I^CT[eps] (Delta_eps I = 0: tr ad = 0; {I_3, I^CT_2}_eps is cubic).
Bubble with special edge (momentum k) and P edge (momentum k - P1): after the exact Gaussian integral and tau = eps a,
t = eps b, (a+b) = T', u = b/T' (region u in [1/2,1], T' in [1/u, 1/(1-u)]):
    eps^{-1}: int_{1/2}^1 (2u-1) (Delta M)(u P1)/4 du
    eps^0  : int_{1/2}^1 log(u/(1-u)) [ M(u P1) - u(1-u) P1^2 (Delta M)(u P1)/4 ] du
    no log eps, M(k) = -k^2 str(Khat S h_1(-k,k-P1) Phat(k-P1) S h_2(-(k-P1),k)).
Route (b) [cross-check]: Q W_1 at finite eps equals the same total (QP(eps,infty) = K_eps); in this sector
Q W_1 = Q(A|_V self-energy) (the (c,A^*) self-energy is proportional to p and is annihilated).  The finite part of every
self-energy bubble is computed exactly from int_0^1 du int_{eps/min(u,1-u)}^infty dT [T^{-1} G(u) + T^{-2} c] e^{-T u(1-u)s}.
Units: hbar/(16 pi^2) before division by c_adj; reported in u = hbar c_adj/(16 pi^2)."""
import json, sys, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import r2_core as C
from r2_core import Q2
from r2_graphs import lap, cached

E4 = [np.eye(4, dtype=int)[m] for m in range(4)]
uu = sp.symbols('u', positive=True)
LOG2 = sp.log(2)
rng = np.random.default_rng(31)

def interp_poly(vals):
    """vals: list of (u (Fraction), Q2) -> sympy polynomial in u with Q(i) coefficients (Lagrange, exact)"""
    poly = 0
    for i, (ui, vi) in enumerate(vals):
        term = sp.Rational(vi.re) + sp.I * sp.Rational(vi.im)
        for j, (uj, _) in enumerate(vals):
            if j != i:
                term = term * (uu - sp.Rational(uj)) / (sp.Rational(ui) - sp.Rational(uj))
        poly += term
    return sp.expand(poly)

def Rat(x):
    return sp.Rational(Fraction(x).numerator, Fraction(x).denominator)

class Ins:
    def __init__(self, bgs):
        z = np.zeros(4, dtype=complex)
        self.h0 = C.hessian(bgs, z, z)
        self.Hj = [C.hessian(bgs, E4[m].astype(complex), z) - self.h0 for m in range(4)]
        self.Hl = [C.hessian(bgs, z, E4[m].astype(complex)) - self.h0 for m in range(4)]
    def __call__(self, kj, kl):
        h = self.h0.copy()
        for m in range(4):
            if kj[m]: h = h + kj[m] * self.Hj[m]
            if kl[m]: h = h + kl[m] * self.Hl[m]
        return h

def special_bubbles(gam, a, p):
    """route (a): returns {(Ls, Lp): {'epsinv': sympy, 'fin': sympy}} in units hbar/(16pi^2) (vertex scale removed)"""
    ext = {'c': ('c', gam, p), 'A': ('A', a, -p)}
    ins = {nm: Ins([e]) for nm, e in ext.items()}
    out = {}
    for (s1, s2) in (('c', 'A'), ('A', 'c')):
        P1 = ext[s1][2]
        assert np.all(P1 % 6 == 0)
        for Ls in C.SLINES:
            Kh = C.Khat(Ls)
            for Lp in C.PLINES:
                def strpart(k, Kh=Kh, Lp=Lp, s1=s1, s2=s2, P1=P1):
                    k2 = k - P1
                    X = C.mm(Kh, C.SDIAG[:, None] * ins[s1](-k, k2))
                    X = C.mm(X, C.Pnum(k2, Lp))
                    X = C.mm(X, C.SDIAG[:, None] * ins[s2](-k2, k))
                    return np.sum(C.SDIAG * np.diag(X))
                def M(k, strpart=strpart):
                    return -int(k @ k) * strpart(k)
                Mc = cached(M)
                # quick zero test
                test = [Mc(np.array([1, 2, -1, 3])), Mc(np.array([-2, 1, 1, 0])), Mc(P1 // 2)]
                if all(abs(t) == 0 for t in test):
                    continue
                us = [Fraction(j, 6) for j in range(7)]
                Mv = [(ui, Q2.of(Mc((P1 * ui.numerator) // ui.denominator))) for ui in us]
                Mpoly = interp_poly(Mv[::2])                 # u = 0,1/3,2/3,1 (cubic)
                for ui, vi in Mv:
                    assert sp.expand(Mpoly.subs(uu, sp.Rational(ui)) - (sp.Rational(vi.re) + sp.I * sp.Rational(vi.im))) == 0
                Lv = [(ui, Q2.of(lap(Mc, (P1 * ui.numerator) // ui.denominator))) for ui in us[::3]]
                Lpoly = interp_poly(Lv[:2])
                assert sp.expand(Lpoly.subs(uu, sp.Rational(Lv[2][0])) - (sp.Rational(Lv[2][1].re) + sp.I * sp.Rational(Lv[2][1].im))) == 0
                s2v = int(P1 @ P1)
                epsinv = sp.integrate((2 * uu - 1) * Lpoly / 4, (uu, sp.Rational(1, 2), 1))
                fin = sp.integrate(sp.expand_log(sp.log(uu) - sp.log(1 - uu)) * sp.expand(Mpoly - uu * (1 - uu) * s2v * Lpoly / 4),
                                   (uu, sp.Rational(1, 2), 1))
                pref = sp.Rational(1, 2) * sp.Rational(1, 4)    # (1/2) str(K'SHPSH); two cubic insertions scaled by 2
                key = (Ls, Lp)
                e = out.setdefault(key, {'epsinv': 0, 'fin': 0, 'orders': []})
                e['epsinv'] = sp.nsimplify(sp.expand(e['epsinv'] + pref * epsinv))
                e['fin'] = sp.expand(sp.simplify(e['fin'] + pref * fin))
                e['orders'].append(s1 + s2)
    return out

def selfenergy_finite(a1, a2, p):
    """route (b) ingredients: per A|_V-A|_V self-energy class, exact eps^{-1}, log eps and finite parts (with the nonlocal
    -(gamma_E + log p^2) * [log coefficient] separated) of the cross coefficient W_1(a1 e^{ipx} + a2 e^{-ipx})."""
    from r2_graphs import GraphSum
    G = GraphSum([('A', 'V', a1, p), ('A', 'V', a2, -p)])
    out = {}
    for (seq, Ls, sgn) in G.enumerate():
        n = len(seq)
        N, qs = G.numerator(seq, Ls)
        Nc = cached(N)
        key = str([(Ls[i], seq[i]) for i in range(n)])
        pref = Fraction(sgn, 2 * n) * Fraction(1, 2 ** sum(1 for b in seq if len(b) == 1) * 8 ** sum(1 for b in seq if len(b) == 2))
        if n == 1:
            ent = {'epsinv': Q2.of(Nc(np.zeros(4, dtype=int))), 'log': 0, 'fin': 0, 'nonlocal': 0}
            ent = {'epsinv': sp.Rational(ent['epsinv'].re) + sp.I * sp.Rational(ent['epsinv'].im), 'log': 0, 'fin': 0, 'nonlocal': 0}
        else:
            P = -qs[1]
            vals = [(Fraction(j, 2), Q2.of(Nc((P * j) // 2))) for j in range(3)]
            Gp = interp_poly(vals)
            cc = Q2.of(lap(Nc, 0 * P)) * Fraction(1, 4)
            c_ = sp.Rational(cc.re) + sp.I * sp.Rational(cc.im)
            s = int(P @ P)
            Mx = sp.Max(uu, 1 - uu)
            intG = sp.integrate(Gp, (uu, 0, 1))
            logc = -intG + c_ * s / 6
            # int G log M  and int u(1-u) log M  (split at 1/2)
            def intlogM(f):
                return sp.integrate(f.subs(uu, uu) * sp.log(1 - uu), (uu, 0, sp.Rational(1, 2))) + sp.integrate(f * sp.log(uu), (uu, sp.Rational(1, 2), 1))
            fin = -intlogM(Gp) - c_ * s / 6 + c_ * s * intlogM(uu * (1 - uu))
            ent = {'epsinv': c_ / 4, 'log': logc, 'fin': sp.expand(sp.simplify(fin)), 'nonlocal': logc}
            # J = epsinv/eps + log*log(eps) + fin + (-(gamma_E + log s)) * nonlocal_coefficient ... with nonlocal = -int G + c s/6 = log
        e = out.setdefault(key, {'epsinv': 0, 'log': 0, 'fin': 0, 'nonlocal': 0})
        for k_ in e:
            e[k_] = sp.expand(e[k_] + sp.Rational(pref.numerator, pref.denominator) * ent[k_])
    return out

def struct_cA(gam, a, p):
    return int(round(np.einsum('a,m,ma->', gam.real, p, a[:4].real)))

def run():
    res = {}
    samples = []
    for it in range(3):
        p = rng.integers(-1, 2, size=4) * 6
        if not p.any(): p[1] = 6
        gam = rng.integers(-2, 3, size=3).astype(complex)
        a = np.zeros((10, 3), complex); a[:4] = rng.integers(-2, 3, size=(4, 3))
        if struct_cA(gam, a, p) == 0:
            a[0, 0] += 1
        S = struct_cA(gam, a, p); s = int(p @ p)
        sb_ = special_bubbles(gam, a, p)
        samples.append((S, s, sb_))
        print('sample', it, 'S=', S, 's=', s, flush=True)
    # per graph: epsinv = x_e * i * S ; fin = x_f * i * p^2 * S   (units hbar/16pi^2 -> divide by c_adj)
    graphs = {}
    for (S, s, sb_) in samples:
        for key, v in sb_.items():
            xe = sp.nsimplify(sp.expand(v['epsinv'] / (sp.I * S * C.C_ADJ)))
            xf = sp.expand(sp.simplify(v['fin'] / (sp.I * s * S * C.C_ADJ)))
            g = graphs.setdefault(key, {'eps_inv_samples': [], 'fin_samples': [], 'orders': v['orders']})
            g['eps_inv_samples'].append(xe); g['fin_samples'].append(xf)
    out = []
    tot_e = 0; tot_f = 0
    for key, g in graphs.items():
        assert all(sp.simplify(x - g['eps_inv_samples'][0]) == 0 for x in g['eps_inv_samples']), g
        assert all(sp.simplify(x - g['fin_samples'][0]) == 0 for x in g['fin_samples']), g
        tot_e += g['eps_inv_samples'][0]; tot_f += g['fin_samples'][0]
        out.append({'special_edge': key[0], 'P_line': key[1], 'orders': g['orders'],
                    'eps_inv_coeff_of_i(gamma,iota_p a)': str(g['eps_inv_samples'][0]),
                    'finite_coeff_of_i p^2 (gamma,iota_p a)': str(sp.nsimplify(g['fin_samples'][0], [LOG2]))})
        print(out[-1], flush=True)
    print('TOTAL route (a): eps^-1', sp.nsimplify(tot_e), ' finite', sp.simplify(tot_f), flush=True)
    return out, sp.nsimplify(tot_e), sp.simplify(tot_f), samples

if __name__ == '__main__':
    out, te, tf, samples = run()
    json.dump({'graphs': out, 'eps_inv_total': str(te), 'finite_total': str(tf)}, open('logs/r2_T2_routeA.json', 'w'), indent=1)
