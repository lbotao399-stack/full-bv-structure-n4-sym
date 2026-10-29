"""Zero-mode kernel of the composite A|_V--c line along K_b (normalisation of the log T coefficient).
Composite segment of a cycle: N_AV(l) S h_S(-l, l-P; P) N_gh(l-P) / (|l|^2 |l-P|^2)   (T = infinity, eps -> 0),
h_S = Hessian of S^(1)(u e^{iP.y}) gamma (kn_engine.s1v_hessian / 8).  The K_b combination of the plane-wave jets
(X0 = sum_{m,a,c} C^m_{ac} d_a d_c at P = 0) of this segment is homogeneous of degree -4 in l:
  G(l^) = |l|^4 X0[segment],  and int d^4l/(2pi)^4 G(l^)/|l|^4 (IR cutoff l^2 ~ 1/T) = <G>_{S^3} log T /(16 pi^2).
In engine units (hbar/(16 pi^2) stripped) the segment's log T coefficient is therefore <G>_{S^3}.
Claim to test (ch08 prop:c8-Tsafe, design thm:o-ir (ii)): <G> = g0 Z_b with Z_b[(A_V mu, colour a), (c, colour a')] =
b_mu delta_{a a'} (the constant modes b (x) (-,-)^{-1} and 1), and g0 = 6 = 16 pi^2 x 3/(8 pi^2)."""
import json, os
import numpy as np
import kn_engine as E

OUT = os.path.dirname(os.path.abspath(__file__))
C = E.C
b = np.array([1.0, 0.0, -1.0, 2.0])
I4 = np.eye(4)
jet = E.jet_order2(4)
LAV = E.LINES['AV']; LGH = E.LINES['gh']
Cm = np.zeros((4, 4, 4))
for m in range(4):
    for a in range(4):
        for c in range(4):
            Cm[m, a, c] = (a == c) * b[m] - b[a] * (c == m) - b[c] * (a == m)


BLOCKS = [E.s1v_insertion(I4[m], 'g', pin_lin=I4).block(LAV.cols, LGH.rows) for m in range(4)]


def segment_X0(l):
    """X0 combination of the composite segment [N_AV/|l|^2] S h_S [N_gh(l-P)/|l-P|^2] (h_S = Hessian/8), exact P-jets"""
    l = np.asarray(l, float)
    l2 = float(l @ l)
    # scalar jet of 1/|l-P|^2 = (1/l^2)[1 + 2 l.P/l^2 - P^2/l^2 + 4 (l.P)^2/l^4 + O(P^3)]
    inv = np.zeros(jet.nb, complex)
    inv[0] = 1.0 / l2
    for v in range(4):
        inv[jet.lin[v]] += 2 * l[v] / l2 ** 2
    for v in range(4):
        for w in range(4):
            e = tuple((1 if k == v else 0) + (1 if k == w else 0) for k in range(4))
            inv[jet.index[e]] += (-(v == w) / l2 + 4 * l[v] * l[w] / l2 ** 2) / l2
    tot = 0
    for m in range(4):
        h0, Hj, Hl, Hp = BLOCKS[m]
        hc = h0 + np.tensordot(-l, Hj, 1) + np.tensordot(l, Hl, 1)
        comps = {0: LAV.Scols[:, None] * hc}
        for v in range(4):
            comps[jet.lin[v]] = LAV.Scols[:, None] * (-Hl[v] + Hp[v])
        SH = E.JMat(jet, comps)
        NAV = E.JMat(jet, {0: LAV.N0 / l2})
        Ngc = LGH.N0 + np.tensordot(l, LGH.Nm, 1)
        dg = {0: Ngc}
        dg.update({jet.lin[v]: -LGH.Nm[v] for v in range(4)})
        NGH = E.JMat(jet, dg)
        M = NAV.matmul(SH).matmul(NGH)
        out = {}
        for i, j, k in jet.tab:
            if i in M.c and abs(inv[j]) > 0:
                out[k] = out.get(k, 0) + M.c[i] * inv[j]
        for a in range(4):
            for c in range(4):
                e = tuple((1 if k == a else 0) + (1 if k == c else 0) for k in range(4))
                val = out.get(jet.index[e], 0)
                tot = tot + Cm[m, a, c] * (2 * val if a == c else val)
    return tot / 8.0


def sphere_average(f, n=16):
    xc, wc = np.polynomial.legendre.leggauss(n)
    acc = 0; wt = 0
    for x1, w1 in zip(xc, wc):              # chi in (0, pi): cos chi = x1, measure sin^2 chi d chi = sqrt(1-x1^2) dx1
        chi = np.arccos(x1)
        for x2, w2 in zip(xc, wc):          # theta: cos theta = x2, measure sin theta d theta = dx2
            th = np.arccos(x2)
            for k in range(n):
                ph = 2 * np.pi * (k + 0.5) / n
                lhat = np.array([np.cos(chi), np.sin(chi) * np.cos(th), np.sin(chi) * np.sin(th) * np.cos(ph),
                                 np.sin(chi) * np.sin(th) * np.sin(ph)])
                w = w1 * np.sqrt(1 - x1 ** 2) * w2 * (2 * np.pi / n)
                acc = acc + w * f(lhat)
                wt += w
    return acc / wt


if __name__ == '__main__':
    # segment_X0 returns X0 of [N_AV S h N_gh /(|l-P|^2)] divided by 8 and by |l|^2 (first propagator): degree -4.
    # check homogeneity
    lh = np.array([0.3, -0.5, 0.7, 0.4]); lh /= np.linalg.norm(lh)
    g1 = segment_X0(1e-2 * lh) * (1e-2) ** 4
    g2 = segment_X0(1e-3 * lh) * (1e-3) ** 4
    avg = sphere_average(lambda u: segment_X0(u), n=int(os.environ.get("KN_NS", "24")))
    Zb = np.zeros_like(avg)
    rows = LAV.rows; cols = LGH.cols
    # Z_b on (A_V rows) x (gh cols): rows index A-slot (M*D + a), cols c-slot (a') and A^*_V slots
    Dc = C.D
    for ia, r_ in enumerate(rows):
        off = r_ - C.OFF['A']; Mdir, a = off // Dc, off % Dc
        for jc, c_ in enumerate(cols):
            if C.OFF['c'] <= c_ < C.OFF['c'] + Dc and (c_ - C.OFF['c']) == a and Mdir < 4:
                Zb[ia, jc] = b[Mdir]
    nz = np.abs(Zb) > 0
    ratio = avg[nz] / Zb[nz]
    rest = np.abs(avg[~nz]).max()
    res = {'homogeneity_|G(0.01)-G(0.001)|/|G|': float(np.max(np.abs(g1 - g2)) / max(1e-300, np.max(np.abs(g1)))),
           'sphere_average_over_Zb_on_support': [[z.real, z.imag] for z in np.unique(np.round(ratio, 10))],
           'max_|<G>|_off_support': float(rest),
           'g0_expected_from_design': 6.0,
           'note': 'g0 = <G>/Z_b in engine units (hbar/(16pi^2)); design: coefficient 3/(8 pi^2) = 6/(16 pi^2) times the cut tree'}
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(OUT, 'kn_zeromode_kernel.json'), 'w'), indent=1)
