"""Calibration: arity-zero special bubbles (c, A|_V) at T -> large, compared with engine A Step-1 (logs/r2_T2_routeA.json):
per-graph eps^-1 and finite coefficients; total eps^-1 = -9/2 i(gamma, iota_p a), finite = -27/8 i p^2 (gamma, iota_p a)
(units hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))).  Also T d/dT of the same graphs at finite T (should be O(eps))."""
import json, time, itertools
import numpy as np
from kt_engine import *
t0 = time.time()
p = np.array([6., 0., -6., 6.]) / 6.0 * 1.0          # |p|^2 = 3
gam = np.array([1., -2., 0.]).astype(complex)
a = np.zeros((10, D), complex); a[:4] = np.array([[1, 0, -1], [2, 1, 0], [0, -1, 2], [-1, 2, 1]])
S = np.einsum('a,m,ma->', gam, p, a[:4]).real; s2 = p @ p
z = hd_vec(np.zeros(4))
ins = {'c': (Insertion('bg', t='c', pol=gam), hd_vec(p), 0), 'A': (Insertion('bg', t='A', pol=a), hd_vec(-p), None)}
print('insertions built', round(time.time() - t0, 1), 's', flush=True)
cycles = []
for order in (('c', 'A'), ('A', 'c')):
    for Ls in SLINES:
        for Lp in PLINES:
            cyc = Cycle([ins[order[0]], ins[order[1]]], [Ls, Lp], 0.5)
            if cyc.nonzero():
                cycles.append(((order, Ls, Lp), GraphPoly(cyc, np.zeros(4), np.zeros(4))))
print('nonzero special bubbles:', [c[0] for c in cycles], round(time.time() - t0, 1), 's', flush=True)
def G_full(gp, eps, T, nsig=12, per=12):
    sg, sw = lin_grid(0, eps, nsig); tg, tw = log_grid(eps, T, per)
    tt = np.array([[x, y] for x in sg for y in tg]); ww = np.array([x * y for x in sw for y in tw])
    return gp.integrate(tt, ww)[0]
def G_TdT(gp, eps, T, nsig=12):
    sg, sw = lin_grid(0, eps, nsig)
    tt = np.array([[x, T] for x in sg]); ww = np.array([x * T for x in sw])
    return gp.integrate(tt, ww)[0]
CADJ = C.C_ADJ
res = {'p': p.tolist(), 'S=(gamma,iota_p a)': S, 'per_graph': {}}
Tbig = 1e9
eps_list = [1e-3, 1e-4, 1e-5]
tot = {e: 0 for e in eps_list}
for key, gp in cycles:
    vals = {e: G_full(gp, e, Tbig) for e in eps_list}
    for e in eps_list: tot[e] += vals[e]
    # Richardson: G = c/eps + d + O(eps)  -> from two smallest eps
    e1, e2 = eps_list[-2], eps_list[-1]
    cm1 = (vals[e1] - vals[e2]) / (1 / e1 - 1 / e2); c0 = vals[e2] - cm1 / e2
    res['per_graph'][str(key)] = {'eps_inv/(i S C_ADJ)': complex(cm1 / (1j * S * CADJ)).real,
                                  'finite/(i p^2 S C_ADJ)': complex(c0 / (1j * s2 * S * CADJ)).real}
    print(key, res['per_graph'][str(key)], flush=True)
e1, e2 = eps_list[-2], eps_list[-1]
cm1 = (tot[e1] - tot[e2]) / (1 / e1 - 1 / e2); c0 = tot[e2] - cm1 / e2
res['total'] = {'eps_inv/(i S C_ADJ)': complex(cm1 / (1j * S * CADJ)).real, 'finite/(i p^2 S C_ADJ)': complex(c0 / (1j * s2 * S * CADJ)).real,
                'engineA_step1': {'eps_inv': -4.5, 'finite': -27 / 8}}
print('TOTAL', res['total'], flush=True)
# T d/dT at finite T: should be O(eps)
res['TdT'] = {}
for T in [1.0, 10.0, 100.0]:
    for eps in [1e-2, 1e-3, 1e-4]:
        v = sum(G_TdT(gp, eps, T) for _, gp in cycles)
        res['TdT']['T=%g eps=%g' % (T, eps)] = complex(v / (1j * S * CADJ)).real
print('TdT', res['TdT'], flush=True)
res['seconds'] = round(time.time() - t0, 1)
json.dump(res, open('calib_arity0.json', 'w'), indent=1)
