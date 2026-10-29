"""Delta(T, eps) = G(eps, T) - G(eps, T_ref = 1) = int_1^T (T' d/dT' G) dT'/T' for the K_v kernel kap0 (and the lower
jets), from T d/dT on a grid of 4 points per decade (Simpson in log T), for the 1PI defect sums (a) and flow sums (b).
The eps-singular parts (eps^-1, log eps) are T-independent and cancel in Delta; lim_{eps->0} Delta(T, eps) is the
T-dependence of the eps -> 0 limit."""
import sys, json, time
import numpy as np
import kv_driver as KD


def simpson_cum(f, x):
    """cumulative Simpson on an odd number of equally spaced points; returns values at even indices"""
    out = [0.0]; h = x[1] - x[0]
    for j in range(2, len(x), 2):
        out.append(out[-1] + h / 3 * (f[j - 2] + 4 * f[j - 1] + f[j]))
    return np.array(out)


if __name__ == '__main__':
    mode = sys.argv[1]; secs = sys.argv[2].split(',')
    eps_list = [1e-3, 1e-4, 1e-5] if mode == 'defect' else [1e-4, 1e-6]
    T_list = [10 ** (j / 4) for j in range(0, 25)]
    res = {}
    for s in secs:
        t0 = time.time()
        out = KD.compute(mode, s, eps_list, T_list, per=10, nsig=8, log=lambda x: None)
        V = out['values']
        lx = np.log(np.array(T_list))
        sec = {}
        for e in eps_list:
            k0 = np.array([complex(*V['eps=%g,T=%g' % (e, T)]['TdT_kap0']) for T in T_list])
            k1 = np.array([[complex(*x) for x in V['eps=%g,T=%g' % (e, T)]['TdT_kap1']] for T in T_list])
            D0 = simpson_cum(k0, lx)
            D1 = np.array([simpson_cum(k1[:, a], lx) for a in range(4)]).T
            sec['eps=%g' % e] = {'T': [T_list[j] for j in range(0, 25, 2)],
                                 'Delta_kap0': [[z.real, z.imag] for z in D0],
                                 'max|Delta_kap1|': [float(np.abs(r).max()) for r in D1]}
        res[s] = sec
        print(mode, s, 'done', round(time.time() - t0, 1), 's', flush=True)
        for e in eps_list:
            d = sec['eps=%g' % e]['Delta_kap0']
            print('   eps=%g  Delta_kap0(T=1e2, 1e4, 1e6) =' % e, ['%.4e%+.4ei' % tuple(d[j]) for j in (4, 8, 12)])
    json.dump(res, open('delta_T_%s.json' % mode, 'w'), indent=1)
