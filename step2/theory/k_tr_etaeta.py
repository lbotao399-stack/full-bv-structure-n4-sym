# k_tr_etaeta.py -- arity two along (eta1, eta2) in Sigma^* x Sigma^*: two floating S^(1)(eta_a) joined by a
# lambda-lambda edge, at a finite transfer endpoint T (k_transfer.tex, Thm thm:ktr-flowlog (arity two)).
#
# Momentum-space symbol of the chain c --[A^*|_V c]-- S^(1)(eta_1) --[lambda lambda]-- S^(1)(eta_2) --[A^*|_V c]-- c
# as in step2/qa/design_audit_2/eta_pair_ir.py (same Euclidean sigma matrices, same Wick-sign assembly
# Total = C1[z,z'] - C1[z',z] + (1<->2)), with every scalar factor 1/|q|^2 of P(eps,infty) replaced by the scalar
# factor of P(eps,T) in the infrared region, h_T(|q|^2) = (1 - e^{-T|q|^2})/|q|^2 (eps -> 0 there; the eps-end only
# matters at |q| ~ eps^{-1/2}).  f, g are the attachments of the rest of the cycle (Taylor data f0, f1, g0, g1).
# Checks:
#   (i)  at finite T the integrand is bounded at q -> 0, so J(T) = int_{|p|<1} d^4p Total_T(p) converges;
#   (ii) J(T2) - J(T1) -> (1/2) * (coefficient of log R) * log(T2/T1), coefficient of log R = 2 pi^2 * 6i v1.(f0 g1 - f1 g0),
#        i.e. coefficient of log T = 6 pi^2 i v1.(f0 g1 - f1 g0)   [half of the large-distance coefficient];
#   (iii) for v1 = Gamma(eta1,eta2)|_V = 0 the log T coefficient vanishes.
import json
import numpy as np

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory/ktr_out/k_tr_etaeta.json'
rng = np.random.default_rng(20261029)
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]], complex); s3 = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2, dtype=complex)
sig = [s1, s2, s3, 1j*I2]
sigb = [-s1, -s2, -s3, 1j*I2]


def sphere_nodes(ns=20, nphi=20):
    xs, ws = np.polynomial.legendre.leggauss(ns); s = 0.5*(xs + 1); w = 0.5*ws
    ph = 2*np.pi*np.arange(nphi)/nphi
    S, P1, P2 = np.meshgrid(s, ph, ph, indexing='ij'); W = np.meshgrid(w, ph, ph, indexing='ij')[0]/nphi**2
    ca, sa = np.sqrt(1 - S), np.sqrt(S)
    pts = np.stack([ca*np.cos(P1), ca*np.sin(P1), sa*np.cos(P2), sa*np.sin(P2)], -1).reshape(-1, 4)
    return pts, W.reshape(-1)


NODES, WTS = sphere_nodes()


def hT(s, T):
    """scalar factor (1-e^{-T s})/s and its derivative in s; T = inf gives 1/s."""
    if np.isinf(T):
        return 1/s, -1/s**2
    em = np.exp(-T*s)
    h = -np.expm1(-T*s)/s
    dh = (T*s*em + np.expm1(-T*s))/s**2
    return h, dh


def mats(p, F, G, T):
    """F = (f0, f1), G = (g0, g1): linear Taylor data.  Left factor d_{q_a}(q_m F(q) h(q^2)) at q = -p,
    middle p_k h(p^2), right d_{p_b}(p_n G(p) h(p^2))."""
    f0, f1 = F; g0, g1 = G
    N = p.shape[0]; pp = (p**2).sum(1)
    q = -p
    h, dh = hT(pp, T)
    Fq = f0 + q@f1; Gp = g0 + p@g1
    L = np.zeros((N, 2, 2), complex); R = np.zeros((N, 2, 2), complex)
    for a in range(4):
        for m in range(4):
            la = (a == m)*h*Fq + q[:, m]*f1[a]*h + q[:, m]*Fq*2*q[:, a]*dh
            L += la[:, None, None]*(sig[a]@sigb[m])[None]
            rb = (a == m)*h*Gp + p[:, m]*g1[a]*h + p[:, m]*Gp*2*p[:, a]*dh
            R += rb[:, None, None]*(sigb[m]@sig[a])[None]
    Mid = sum((p[:, k]*h)[:, None, None]*sig[k][None] for k in range(4))
    return L, Mid, R


def I_or(p, e1, eb2, F, G, T):
    L, Mid, R = mats(p, F, G, T)
    return np.einsum('ai,nab,bi->n', e1, L@Mid@R, eb2)


def total(p, eta, F, G, T):
    e1, eb1, e2, eb2 = eta
    out = 0
    for (ea, ebb) in [(e1, eb2), (e2, eb1)]:
        out = out + I_or(p, ea, ebb, F, G, T) - I_or(p, ea, ebb, G, F, T)
    return out


def bil(e1, eb2, m):
    return np.einsum('ai,ab,bi->', e1, sig[m], eb2)


def rnd():
    return rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4))


def J(T, eta, F, G, rmin):
    """int_{rmin<|p|<1} d^4p Total = 2 pi^2 int r^3 <Total>_{S^3} dr, log-radius Gauss-Legendre per decade."""
    xr, wr = np.polynomial.legendre.leggauss(24)
    edges = np.arange(np.floor(np.log10(rmin)), 0.0 + 1e-12, 0.5)
    val = 0
    for lo, hi in zip(edges[:-1], edges[1:]):
        tt = (lo + (xr + 1)/2*(hi - lo))*np.log(10); wt = wr/2*(hi - lo)*np.log(10)
        for ti, wi in zip(tt, wt):
            r = np.exp(ti)
            val += wi*r**4*(total(r*NODES, eta, F, G, T)*WTS).sum()
    return 2*np.pi**2*val


res = {}
e1, eb1, e2, eb2 = rnd(), rnd(), rnd(), rnd()
eta = (e1, eb1, e2, eb2)
v1 = np.array([-1j*(bil(e2, eb1, c) + bil(e1, eb2, c)) for c in range(4)])
f0, g0 = 0.7 - 0.2j, -0.4 + 0.9j
f1 = rng.normal(size=4) + 1j*rng.normal(size=4); g1 = rng.normal(size=4) + 1j*rng.normal(size=4)
F, G = (f0, f1), (g0, g1)
pred_logT = 6*np.pi**2*1j*(v1@(f0*g1 - f1*g0))
# (i) boundedness at p -> 0 for finite T
T = 1e4
small = total(1e-6*NODES, eta, F, G, T)
res['i_bounded_at_p_to_0_T=1e4'] = {'max_abs_|Total|_at_|p|=1e-6': float(np.abs(small).max()),
                                     'compare_T=inf_at_|p|=1e-6': float(np.abs(total(1e-6*NODES, eta, F, G, np.inf)).max())}
# (ii) J(T) for increasing T
Js = {}
for Tv in [1e2, 1e4, 1e6]:
    Js[Tv] = J(Tv, eta, F, G, rmin=1e-3/np.sqrt(Tv))
res['ii_J(T)'] = {'T=%g' % k: repr(complex(v)) for k, v in Js.items()}
d1 = (Js[1e4] - Js[1e2])/np.log(1e2); d2 = (Js[1e6] - Js[1e4])/np.log(1e2)
res['ii_logT_slope'] = {'(J(1e4)-J(1e2))/log100': repr(complex(d1)), '(J(1e6)-J(1e4))/log100': repr(complex(d2)),
                        'prediction_6 pi^2 i v1.(f0 g1 - f1 g0)': repr(complex(pred_logT)),
                        'rel_dev_last': float(abs(d2 - pred_logT)/abs(pred_logT))}
# (ii') independence of the lower radius (convergence at p -> 0): J with rmin 10x smaller
Jb = J(1e4, eta, F, G, rmin=1e-4/np.sqrt(1e4))
res['ii_rmin_independence_T=1e4'] = {'J(rmin=1e-3/sqrtT)': repr(complex(Js[1e4])), 'J(rmin=1e-4/sqrtT)': repr(complex(Jb)),
                                     'abs_diff': float(abs(Jb - Js[1e4]))}
# (iii) pair with v1 = 0
basis = []
for idx in range(16):
    E = np.zeros(16, complex); E[idx] = 1
    e2t, eb2t = E[:8].reshape(2, 4), E[8:].reshape(2, 4)
    basis.append([-1j*(bil(e2t, eb1, c) + bil(e1, eb2t, c)) for c in range(4)])
Amat = np.array(basis).T
_, sv, Vh = np.linalg.svd(Amat)
ker = Vh[4:].conj().T
vec = ker@(rng.normal(size=12) + 1j*rng.normal(size=12))
e2k, eb2k = vec[:8].reshape(2, 4), vec[8:].reshape(2, 4)
etak = (e1, eb1, e2k, eb2k)
v1k = np.array([-1j*(bil(e2k, eb1, c) + bil(e1, eb2k, c)) for c in range(4)])
Jk = {Tv: J(Tv, etak, F, G, rmin=1e-3/np.sqrt(Tv)) for Tv in [1e2, 1e4]}
res['iii_v1=0'] = {'|v1|': float(np.abs(v1k).max()),
                   '(J(1e4)-J(1e2))/log100': repr(complex((Jk[1e4] - Jk[1e2])/np.log(1e2))),
                   'generic_pair_slope_for_scale': repr(complex(d2))}
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
