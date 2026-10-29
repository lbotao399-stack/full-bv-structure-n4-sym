# k_tr_crosscheck.py -- cross-checks of k_transfer.tex against the independent finite-T graph engine of
# step2/qa/k_transfer_numeric (kt_engine.py on engine_a/r2_core.py: su(2), explicit Schwinger integrals, T d/dT by
# pinning one edge at t = T).  That engine was written before this file and shares no code with cubic_core.py.
#
# (1) Reproducibility: rerun one flow sector and one defect sector with the engine's finer grid and compare with the
#     stored values at common (eps, T).
# (2) Flow (arity one along K_b, 1PI two-point, no special edge): the T -> infinity limit of T d/dT of the
#     translation-invariant coefficient kap0 must equal (3/(8 pi^2)) hbar Y_W, i.e. 6 Y_W in the engine's unit
#     hbar/(16 pi^2), with the closed forms of Y_W (k_tr_zeromode.py, cubic_core conventions):
#        (A_Wperp, A*_Wperp): Y_W = -2 i C (b.p)/|p|^2 (u, u*)
#        (A_V, A*_V):         Y_W = i C/|p|^2 [ -2 (b.p)(u.u*) + (b.u)(p.u*) ]
#        (c, c*):             Y_W = i C (b.p)/|p|^2 (u_c u_c*)
#     C = tr_g(ad_X ad_Y)/(X,Y) = -2 (su(2), orthonormal), colour pairing of the legs = 1 (both legs along T_0).
#     The first-derivative and second-derivative coefficients kap1, kap2 must have T d/dT -> 0 (the log T term is
#     translation invariant).
# (3) Defect (one special edge K_eps - K_0 on the cycle): T d/dT kap = O(eps) at fixed T, and T * (T d/dT kap) tends to
#     a constant as T -> infinity, i.e. no log T at fixed eps: Q Y_W(eps) = 0 in every two-point sector
#     (k_tr_defectlog.py).
import sys
import json
import time
import numpy as np

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
NUM = ROOT + '/step2/qa/k_transfer_numeric'
OUT = ROOT + '/step2/theory/ktr_out/k_tr_crosscheck.json'
sys.path.insert(0, NUM)
import kv_driver as kv  # noqa: E402

res = {'engine': 'step2/qa/k_transfer_numeric/kt_engine.py + kv_driver.py (independent of cubic_core.py)'}


def load(mode, name):
    return json.load(open('%s/kv_%s_%s_main.json' % (NUM, mode, name)))


def c(x):
    return complex(x[0], x[1])


# ------------------------------------------------------------------ (1) reproducibility
t0 = time.time()
rr = {}
for mode, name, eps_l, T_l in [('flow', 'AWAsW', [1e-6], [1e2, 1e6]), ('defect', 'AWAW', [1e-4, 1e-6], [1e2, 1e5])]:
    out = kv.compute(mode, name, eps_l, T_l, per=15, nsig=12, log=lambda s: None)
    st = load(mode, name)
    cmp = {}
    for key, v in out['values'].items():
        a = c(v['TdT_kap0']); b_ = c(st['values'][key]['TdT_kap0'])
        cmp[key] = {'rerun_finer_grid': repr(a), 'stored': repr(b_), 'rel_dev': abs(a - b_)/abs(b_)}
    rr['%s_%s' % (mode, name)] = cmp
res['1_reproducibility'] = rr
res['1_seconds'] = round(time.time() - t0, 1)

# ------------------------------------------------------------------ (2) flow: log T coefficient vs zero-mode trees
C = -2.0
flow = {}
for name in ['AWAsW', 'AVAsV', 'ccs', 'lls']:
    d = load('flow', name)
    p = np.array(d['p1']); bb = np.array(d['b'], float)
    pol1 = np.array(d['pol1']); pol2 = np.array(d['pol2'])
    V = d['values']
    Tmax = max(v['T'] for v in V.values()); emin = min(v['eps'] for v in V.values())
    v = V['eps=%g,T=%g' % (emin, Tmax)]
    num = c(v['TdT_kap0'])
    k1 = max(abs(c(x)) for x in v['TdT_kap1']); k2 = max(abs(c(x)) for r in v['TdT_kap2'] for x in r)
    if name == 'AWAsW':
        u, us = pol1[4:10, 0], pol2[4:10, 0]
        Y = -2j*C*(bb @ p)/(p @ p)*(u @ us)
    elif name == 'AVAsV':
        u, us = pol1[:4, 0], pol2[:4, 0]
        Y = 1j*C/(p @ p)*(-2*(bb @ p)*(u @ us) + (bb @ u)*(p @ us))
    elif name == 'ccs':
        Y = 1j*C*(bb @ p)/(p @ p)*1.0
    else:
        Y = None
    ent = {'eps': emin, 'T': Tmax, 'TdT_kap0_engine': repr(num), 'max|TdT kap1|': k1, 'max|TdT kap2|': k2}
    if Y is not None:
        ent['prediction_6_Y_W'] = repr(6*Y)
        ent['rel_dev'] = abs(num - 6*Y)/abs(6*Y)
    else:
        ent['prediction_6_Y_W'] = ('not compared: the engine uses the spinor realisation engine_a/cliff.py, the closed form '
                                   '(rho(ip) rho(b) u, u*) C/|p|^2 is in the realisation of classical/n4_model.py')
    # T-approach: T d/dT kap0 at T = 1e2, 1e4, 1e6 (eps min): converges to a constant (log T coefficient)
    ent['TdT_kap0_vs_T'] = {('T=%g' % T): repr(c(V['eps=%g,T=%g' % (emin, T)]['TdT_kap0'])) for T in [1.0, 1e2, 1e4, 1e6]
                            if 'eps=%g,T=%g' % (emin, T) in V}
    flow[name] = ent
res['2_flow_logT_coefficient'] = flow

# ------------------------------------------------------------------ (3) defect: O(eps) and no log T
defect = {}
for name in ['AWAW', 'AVAV', 'll', 'cAs']:
    d = load('defect', name)
    V = d['values']
    rows = {}
    for T in [1.0, 1e2, 1e4, 1e6]:
        r = {}
        for e in [1e-2, 1e-4, 1e-6]:
            v = V['eps=%g,T=%g' % (e, T)]
            r['eps=%g' % e] = {'TdT_kap0/eps': repr(c(v['TdT_kap0'])/e), 'T*TdT_kap0/eps': repr(T*c(v['TdT_kap0'])/e)}
        rows['T=%g' % T] = r
    a = T = None
    big = [c(V['eps=1e-06,T=%g' % T]['TdT_kap0'])*T/1e-6 for T in [1e4, 1e5, 1e6]]
    defect[name] = {'table': rows,
                    'T*TdT_kap0/eps at eps=1e-6, T=1e4,1e5,1e6': [repr(x) for x in big],
                    'no_logT': bool(abs(big[2] - big[1]) < 1e-3*abs(big[2])),
                    'proportional_to_eps': bool(abs(c(V['eps=1e-06,T=100']['TdT_kap0'])/1e-6
                                                    - c(V['eps=1e-05,T=100']['TdT_kap0'])/1e-5)
                                                < 1e-3*abs(c(V['eps=1e-06,T=100']['TdT_kap0'])/1e-6))}
res['3_defect_TdT'] = defect
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps({k: res[k] for k in res if k != '3_defect_TdT'}, indent=1))
print(json.dumps({n: {k: v for k, v in defect[n].items() if k != 'table'} for n in defect}, indent=1))
