"""log T coefficients from the exact integrand at T = infinity (no Schwinger integrals), general cycles.

For a one-loop cycle with loop momentum k, the integrand at T = infinity is F(k) = str[ prod_i phi_i(k_i) S h_i ] *
prod_{P edges} 1/k_i^2 * (e^{-eps k_s^2} - 1) (special edge), with the hyper-dual P-jets of the inserted momenta.
Near k = 0 (the region of the floating insertions), F(r khat) = sum_j r^j f_j(khat).  With the transfer endpoint T the
propagators carry (1 - e^{-T k^2}) and
   coefficient of log T in int d^4k/(2pi)^4 F = <f_{-4}>_{S^3} / (16 pi^2),
provided the more singular terms average to zero (f_{-5} odd).  In the units hbar/(16 pi^2) the coefficient is <f_{-4}>.
f_{-4} is extracted from g(r) = r^4 F(r khat) = f_{-5}/r + f_{-4} + f_{-3} r + f_{-2} r^2 at r = h, 2h, 4h, 8h.
Validated against the Schwinger numerics of kv_driver (arity one, K_v) and applied to the arity-two eta-eta chains."""
import sys, json, itertools, time
import numpy as np
from kt_engine import *
import kv_driver as KD
import eta_driver as ED
from ir_check import sphere_nodes


def hd_inv_sq(k):
    """1/|k|^2 for a hyper-dual momentum k = (re, v1, v2) -> [f, d1 f, d2 f, d1 d2 f]"""
    re = k[0]; v1 = np.zeros(4) if k[1] is None else k[1]; v2 = np.zeros(4) if k[2] is None else k[2]
    r2 = re @ re
    return np.array([1 / r2, -2 * (re @ v1) / r2 ** 2, -2 * (re @ v2) / r2 ** 2,
                     8 * (re @ v1) * (re @ v2) / r2 ** 3 - 2 * (v1 @ v2) / r2 ** 2])


def hd_special(k, eps):
    """e^{-eps |k|^2} - 1 as a hyper-dual scalar"""
    re = k[0]; v1 = np.zeros(4) if k[1] is None else k[1]; v2 = np.zeros(4) if k[2] is None else k[2]
    e = np.exp(-eps * (re @ re))
    a1 = -2 * eps * (re @ v1); a2 = -2 * eps * (re @ v2); a12 = -2 * eps * (v1 @ v2)
    return np.array([e - 1, e * a1, e * a2, e * (a12 + a1 * a2)])


def hd_mul(x, y):
    return np.array([x[0] * y[0], x[0] * y[1] + x[1] * y[0], x[0] * y[2] + x[2] * y[0],
                     x[0] * y[3] + x[1] * y[2] + x[2] * y[1] + x[3] * y[0]])


def full_integrand(cyc, k, eps):
    val = cyc.value(k)                      # numerator, HD
    for i, L in enumerate(cyc.lines):
        ki = hd_vec_add(hd_vec(k), cyc.q[i])
        fac = hd_special(ki, eps) if L.special else hd_inv_sq(ki)
        val = hd_mul(val, fac)
    return val


def soft_points(cyc):
    """distinct loop momenta at which some line is soft (real parts of -q_i at P = 0)"""
    pts = []
    for q in cyc.q:
        c = -np.asarray(q[0], float)
        if not any(np.allclose(c, x) for x in pts):
            pts.append(c)
    return pts


def f_minus4(cyc, khat, eps, h=2e-3):
    """sum over the soft points k* of the r^-4 (and r^-5) Laurent coefficients of F(k* + r khat)"""
    rs = h * np.array([1., 2., 4., 8.])
    Mx = np.array([[1 / r, 1, r, r * r] for r in rs])
    f4 = 0; f5 = 0
    for kc in soft_points(cyc):
        G = np.array([r ** 4 * full_integrand(cyc, kc + r * khat, eps) for r in rs])   # g(r) = A/r + B + C r + D r^2
        coef = np.linalg.solve(Mx, G)
        f4 = f4 + coef[1]; f5 = f5 + coef[0]
    return f4, f5


def logT_K(sec, eps, nodes, wts):
    """arity one (K_v): per structure log T coefficient of kap0"""
    out = []
    for si, (order, L, pref) in enumerate(sec.structs):
        acc = 0; acc5 = 0
        for pt, w in zip(nodes, wts):
            J = {}; J5 = {}
            for run in KD.RUNS_B0:
                u, d1, d2 = KD.run_dirs(run[1], run[0])
                ins = sec.ins(u, d1, d2)
                cyc = Cycle([ins[o] for o in order], list(L), pref)
                J[run], J5[run] = f_minus4(cyc, pt, eps)
            acc += w * KD.kernels(J)[0]; acc5 += w * KD.kernels(J5)[0]
        out.append((si, acc, acc5))
    return out


def logT_eta(sec, eps, nodes, wts, only=None, arange=range(4)):
    """arity two (eta1, eta2): per structure log T coefficient of kap0 = -sum_ab (e1 e2 component of G_ab);
    arange restricts the first index a (partial sums, added afterwards)"""
    out = []
    cls = [ED.classify(o, L) for (o, L, _) in sec.structs]
    for si, (order, L, pref) in enumerate(sec.structs):
        if only is not None and cls[si] not in only:
            continue
        acc = 0; acc5 = 0
        for pt, w in zip(nodes, wts):
            for a in arange:
                for b in range(4):
                    ins = sec.eins(a, b)
                    cyc = Cycle([ins[o] for o in order], list(L), pref)
                    f4, f5 = f_minus4(cyc, pt, eps)
                    acc += -w * f4[3]; acc5 += -w * f5[3]
        out.append((si, cls[si], acc, acc5))
    return out


if __name__ == '__main__':
    what = sys.argv[1]
    nodes, wts = sphere_nodes(int(sys.argv[2]) if len(sys.argv) > 2 else 5, int(sys.argv[3]) if len(sys.argv) > 3 else 8)
    t0 = time.time()
    res = {}
    if what == 'K':
        for mode, name, eps in [('flow', 'AWAsW', 0.0), ('defect', 'AWAW', 1e-3)]:
            sec = KD.Sector(mode, name)
            r = logT_K(sec, eps, nodes, wts)
            tot = sum(x[1] for x in r)
            res['%s_%s' % (mode, name)] = {'eps': eps, 'per_structure': [['-'.join(sec.structs[si][0]), '-'.join(sec.structs[si][1]),
                                           [v.real, v.imag], [v5.real, v5.imag]] for si, v, v5 in r], 'total': [tot.real, tot.imag]}
            print(mode, name, 'total log T coefficient', tot, round(time.time() - t0, 1), 's', flush=True)
    else:
        mode, name, eps = sys.argv[4], sys.argv[5], float(sys.argv[6])
        only = tuple(sys.argv[7].split(',')) if len(sys.argv) > 7 and sys.argv[7] != 'all' else None
        arange = [int(x) for x in sys.argv[8].split(',')] if len(sys.argv) > 8 else range(4)
        sec = ED.EtaSector(mode, name)
        r = logT_eta(sec, eps, nodes, wts, only=only, arange=arange)
        tot = sum(x[2] for x in r); totc = sum(x[2] for x in r if x[1] == 'chain')
        res['%s_%s_eps=%g' % (mode, name, eps)] = {
            'eps': eps, 'Gamma(eta1,eta2)|_V': [[x.real, x.imag] for x in ED.V1],
            'per_structure': [['-'.join(sec.structs[si][0]), '-'.join(sec.structs[si][1]), c, [v.real, v.imag], [v5.real, v5.imag]] for si, c, v, v5 in r],
            'total': [tot.real, tot.imag], 'total_chain': [totc.real, totc.imag]}
        print(mode, name, eps, 'total log T coefficient', tot, ' chains', totc, round(time.time() - t0, 1), 's', flush=True)
        for si, c, v, v5 in r:
            if abs(v) > 1e-9 * max(1, abs(tot)) or abs(v5) > 1e-6:
                print('   ', '-'.join(sec.structs[si][0]), '-'.join(sec.structs[si][1]), c, v, ' f-5 avg', v5)
    tag = what if what == 'K' else '%s_%s_%s_%g_%s_a%s' % (what, sys.argv[4], sys.argv[5], float(sys.argv[6]),
                                                              sys.argv[7] if len(sys.argv) > 7 else 'all', (sys.argv[8] if len(sys.argv) > 8 else '0123').replace(',', ''))
    json.dump(res, open('ir_eta_%s.json' % tag, 'w'), indent=1)
    print('seconds', round(time.time() - t0, 1))
