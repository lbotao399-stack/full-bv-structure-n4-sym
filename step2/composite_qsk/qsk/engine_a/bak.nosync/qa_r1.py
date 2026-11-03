"""Engine A, part Q, gate G-R1 on the scalar bilinears: Cor c8-homotopy (i) at order hbar,
   [s, Lambda_1(xi)] + [L_1^(0), Lambda_0(xi)] = 0  on O = int omega (iota_w1 A, iota_w2 A)  (sO = 0),
evaluated in the two-field sector (c, lam) on plane-wave legs (exact):
   A = s Q_1(O)|_(c,lam)     (from the gauge-variant trace term a3 F3: s A_V = -dc; F1, F2 are s-closed, a4 = 0),
   B = L_1^(0)({S^(1)(xi), O}) = L_1^(0)(-F1)   (arity-zero two-body kernels of the pairs (A_perp, lam_+) of F1),
   C = {S^(1)(xi), L_1^(0)(O)}  (arity-zero kernel of (A_perp, A_perp) on (c, A_V) with A_V <- xi.lam = -Gamma(xi,lam)).
The graded commutator of the odd L_1^(0) with the odd (parity-changing) {S^(1)(xi), -} is an anticommutator; the test
reports A + e1 B + e2 C for e1, e2 = +-1 and both Koszul readings of the derivation on the two output orderings."""
from fractions import Fraction
import itertools
import numpy as np
import qa_s1 as S1, qa_core as QC, qa_util as U

SB = S1.SB


def run(xi, a3, npts=3, seed=77, log=print):
    R = U.Rng(seed)
    g0 = QC.zero_gen()
    res = []
    for ip in range(npts):
        w1, w2 = R.gv(6), R.gv(6)
        ww = complex(np.dot(w1, w2))
        for lt in ('lamM', 'lamP'):
            qc, ql = R.mom(), R.mom()
            l = R.vec(lt)
            legs = [('c', np.array([1 + 0j]), qc), (lt, l, ql)]
            # A
            A = U.vzero()
            if lt == 'lamM':
                z = ww * sum(-1j * qc[m] * np.dot(SB[m] @ xi, l) for m in range(4))
                A = U.vscale(a3, U.gq(z))
            # B: L_1^(0)(-F1)
            B = U.vzero()
            for wa, wb in ((w1, w2), (w2, w1)):
                zl = sum(wb[i] * (SB[4 + i] @ xi) for i in range(6))
                Xa = ('AW', wa, None); Xb = ('lamP', zl, None)
                B = U.vadd(B, U.both_orders(Xa, Xb, legs, g0)['fin'])
            B = U.vneg(B)
            # C: {S, L_1^(0) O}: A_V leg <- -Gamma(xi, l)|_V at momentum ql (only for lam_- legs: Gamma|_V pairs S_+ W with S_- W*)
            Xa = ('AW', w1, None); Xb = ('AW', w2, None)
            pv = np.array([-(xi @ SB[m] @ l) for m in range(4)])
            Cs = {}
            for kos in ('c_first_minus', 'A_first_minus'):
                Ct = U.vzero()
                for order in ((0, 1), (1, 0)):
                    lg = [('c', np.array([1 + 0j]), qc), ('AV', pv, ql)]
                    lg = [lg[i] for i in order]
                    v = QC.pair_value(Xa, Xb, lg, g0)['fin']
                    c_first = order == (0, 1)
                    sg = -1 if ((kos == 'c_first_minus') == c_first) else 1
                    Ct = U.vadd(Ct, U.vscale(v, (Fraction(sg), Fraction(0))))
                Cs[kos] = Ct
            row = {'point': ip, 'lam': lt, 'A': U.vstr(A), 'B': U.vstr(B), 'C': {k: U.vstr(v) for k, v in Cs.items()}}
            combos = {}
            for e1 in (1, -1):
                for e2 in (1, -1):
                    for kos, C in Cs.items():
                        tot = U.vadd(A, U.vadd(U.vscale(B, (Fraction(e1), Fraction(0))), U.vscale(C, (Fraction(e2), Fraction(0)))))
                        combos['e1=%d,e2=%d,%s' % (e1, e2, kos)] = U.viszero(tot)
            row['vanishing_combos'] = [k for k, v in combos.items() if v]
            res.append(row)
            log('  G-R1 point %d %s: A=%s | B=%s | C=%s | vanishing: %s' % (ip, lt, row['A'], row['B'], row['C'], row['vanishing_combos']))
    return res
