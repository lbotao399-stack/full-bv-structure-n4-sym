"""Independent evaluation of the log T coefficients (zero-mode tree), without Schwinger integrals.

For a cycle with K_v floating between two lines (A|_V and c-A^*|_V), the coefficient of log T in the K_v kernel kap0
(units hbar/(16 pi^2)) is the S^3 average of |k|^4 times the degree -4 part of
   -C^n_ab d_{P_a} d_{P_b} str[ ... phi_i(k) S h_K(-k, k-P; u = e_n) phi_{i+1}(k-P) ... ] / (k^2 (k-P)^2) |_{P=0},
all other factors of the cycle (the 'rest', containing the external legs and, for the defect, the special edge
(e^{-eps q^2} - 1) Khat at its momentum q) frozen at loop momentum 0 and P = 0: the tree with the two freed legs at zero
momentum.  Reason: int d^4k/(2pi)^4 f(khat)/k^4 over T^{-1/2} < |k| < O(1) = <f>_{S^3} log T / (16 pi^2) + O(1).
Compared with lim_{T->inf} T dW/dT of kv_driver (per structure)."""
import sys, json, itertools
import numpy as np
from kt_engine import *
import kv_driver as KD


def sphere_nodes(ns=16, nphi=16):
    xs, ws = np.polynomial.legendre.leggauss(ns); s = 0.5 * (xs + 1); w = 0.5 * ws
    ph = 2 * np.pi * np.arange(nphi) / nphi
    S, P1, P2 = np.meshgrid(s, ph, ph, indexing='ij'); W = np.meshgrid(w, ph, ph, indexing='ij')[0] / nphi ** 2
    ca, sa = np.sqrt(1 - S), np.sqrt(S)
    pts = np.stack([ca * np.cos(P1), ca * np.sin(P1), sa * np.cos(P2), sa * np.sin(P2)], -1).reshape(-1, 4)
    return pts, W.reshape(-1)


NODES, WTS = sphere_nodes(5, 8)   # exact for polynomials of degree <= 9 on S^3 (integrand: degree <= 4)


def hd_scalar_mul(x, y):
    return [x[0] * y[0], x[0] * y[1] + x[1] * y[0], x[0] * y[2] + x[2] * y[0], x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0]]


def ir_value(sec, order, L, pref, u, d1, d2, kvec, eps):
    """HD value of the cycle with the K segment at (k, P hyper-dual) and the rest frozen at k=0, P=0"""
    ins = sec.ins(u, d1, d2)
    n = len(order)
    iK = order.index('K')
    cyc = Cycle([ins[o] for o in order], list(L), pref)
    zero = hd_vec(np.zeros(4))
    # line momenta: line iK enters K with k, line iK+1 leaves with k - P; rest lines at their k=0, P=0 values
    ks = []
    for i in range(n):
        q = cyc.q[i]
        if i == iK:
            ks.append(hd_vec(kvec))
        elif i == (iK + 1) % n:
            ks.append((np.asarray(kvec, float), -np.asarray(d1, float), -np.asarray(d2, float)))
        else:
            ks.append(hd_vec(q[0] - cyc.q[iK][0]))        # loop momentum such that line iK carries 0 -> shift
    X = None
    for i in range(n):
        Li, Ln = cyc.lines[i], cyc.lines[(i + 1) % n]
        insert, pin, _ = cyc.ins[i]
        if i == iK:
            kin, kout, pv = ks[i], ks[(i + 1) % n], pin
        else:
            kin = (ks[i][0], None, None) if i != (iK + 1) % n else (np.zeros(4), None, None)
            kout = (ks[(i + 1) % n][0], None, None) if (i + 1) % n != iK else (np.zeros(4), None, None)
            pv = (pin[0], None, None)
        phi = Li.at(ks[i] if i in (iK, (iK + 1) % n) else (ks[i][0], None, None))
        neg = (-kin[0], None if kin[1] is None else -kin[1], None if kin[2] is None else -kin[2])
        h = insert.at(Li.cols, Ln.rows, neg, kout, pv)
        Mi = phi.mat(h.scale_rows(SD[Li.cols]))
        if Li.special:   # special edge scalar factor at its frozen momentum
            qs = ks[i][0]
            Mi = HD(Mi.a * (np.exp(-eps * qs @ qs) - 1), None if Mi.b is None else Mi.b * (np.exp(-eps * qs @ qs) - 1),
                    None if Mi.c is None else Mi.c * (np.exp(-eps * qs @ qs) - 1), None if Mi.d is None else Mi.d * (np.exp(-eps * qs @ qs) - 1))
        elif i not in (iK, (iK + 1) % n):
            qq = ks[i][0] @ ks[i][0]
            fac = 1.0 / qq          # rest P edge at T -> infinity, eps -> 0: 1/q^2
            Mi = HD(Mi.a * fac, None if Mi.b is None else Mi.b * fac, None if Mi.c is None else Mi.c * fac, None if Mi.d is None else Mi.d * fac)
        X = Mi if X is None else X.mat(Mi)
    r0 = cyc.lines[0].rows
    val = [0.0 if isinstance(c, float) else np.sum(SD[r0] * np.diag(c)) for c in X.comps()]
    val = [v * cyc.sign * pref for v in val]
    # scalar factors 1/k^2 and 1/(k-P)^2 of the two adjacent lines
    k = np.asarray(kvec, float); k2 = k @ k
    g1 = [1 / k2, 0, 0, 0]
    g2 = [1 / k2, 2 * (k @ d1) / k2 ** 2, 2 * (k @ d2) / k2 ** 2, 8 * (k @ d1) * (k @ d2) / k2 ** 3 - 2 * (d1 @ d2) / k2 ** 2]
    return np.array(hd_scalar_mul(hd_scalar_mul(val, g1), g2))


def ir_kap0(sec, order, L, pref, eps, nodes=NODES, wts=WTS):
    tot = 0
    for pt, w in zip(nodes, wts):
        J = {}
        for run in KD.RUNS_B0:
            u, d1, d2 = KD.run_dirs(run[1], run[0])
            J[run] = ir_value(sec, order, L, pref, u, d1, d2, pt, eps)
        k0, _, _ = KD.kernels(J)
        tot += w * k0
    return tot


if __name__ == '__main__':
    out = {}
    for mode, name, eps in [('flow', 'AWAsW', 0.0), ('flow', 'lls', 0.0), ('flow', 'ccs', 0.0), ('flow', 'AVAsV', 0.0),
                            ('defect', 'AWAW', 1e-6), ('defect', 'cAs', 1e-6), ('defect', 'AVAV', 1e-6), ('defect', 'll', 1e-6)]:
        sec = KD.Sector(mode, name)
        num = json.load(open('kv_%s_%s_main.json' % (mode, name)))
        V = num['values']
        epsn = min(v['eps'] for v in V.values()); Tn = max(v['T'] for v in V.values())
        res = []
        for si, (order, L, pref) in enumerate(sec.structs):
            iK = order.index('K'); n = len(order)
            adj = {L[iK], L[(iK + 1) % n]}
            if adj != {'AV', 'gh'}:
                continue
            pred = ir_kap0(sec, order, L, pref, eps if mode == 'defect' else 0.0)
            numv = complex(*V['eps=%g,T=%g' % (epsn, Tn)]['per_structure_TdT'][str(si)]['kap0'])
            res.append({'structure': ['-'.join(order), '-'.join(L)], 'IR_prediction': [pred.real, pred.imag],
                        'numeric_TdT_at_T=%g_eps=%g' % (Tn, epsn): [numv.real, numv.imag]})
            print(mode, name, res[-1], flush=True)
        tp = sum(complex(*r['IR_prediction']) for r in res)
        out['%s_%s' % (mode, name)] = {'floating_structures': res, 'sum_IR_prediction': [tp.real, tp.imag],
                                       'numeric_total_TdT': V['eps=%g,T=%g' % (epsn, Tn)]['TdT_kap0']}
        print('  SUM prediction', tp, ' numeric total', V['eps=%g,T=%g' % (epsn, Tn)]['TdT_kap0'], flush=True)
    json.dump(out, open('ir_check.json', 'w'), indent=1)
