"""Engine A, part Q, gate G-R1 on the scalar bilinears (rewritten in the revision of 2026-10-30, findings Q5, Q8).

Cor c8-homotopy (i) at order hbar, with Lambda^(1)(x) = {S^(1)(x),-} + Q_1 (Def c8-components, |O| = 0), Lambda^(0) =
s + hbar L_1^(0), both odd as operators once the constant spinor xi is stripped, so the graded commutator is an
anticommutator; on O = int omega (iota_w1 A, iota_w2 A) with sO = 0:
      s Q_1(O) + L_1^(0)({S^(1)(xi),O}) + {S^(1)(xi), L_1^(0)(O)} = 0 .
Evaluated on the two-field sector (c, lam_-) with one consistent Grassmann bookkeeping:
  * an engine value with ordered legs (Y_1, Y_2) is the coefficient of theta_1 theta_2 for the chains visiting the
    outputs in that order (Def c1-vertex); the functional's coefficient in the fixed order (c, lam) is
    v(c,lam) - v(lam,c) for two odd outputs and v(c,A_V) + v(A_V,c) for one even output;
  * {S^(1)(xi), O} = +F1 (Cor c7-Q; Lemma qs-arrows with the corrected sign: A -> +Gamma(xi,lam));
  * {S^(1)(xi), -} is an odd derivation; on a functional of (c, A_V) it passes c and replaces A_V by
    +Gamma(xi,lam)|_V: {S, g(c, A_V)} = -g(c, Gamma(xi,lam)|_V).
Terms (coefficient of theta_c theta_lam):
  A = s(a3 F3),  F3 = (w1,w2) int omega (rho(A|_V) xi, lam),  s: A_V -> -dc;
  B = L_1^(0)(F1)  (arity-zero two-body kernels of the pairs (A_perp, lam_+) of F1, odd-odd antisymmetrised);
  C = {S^(1)(xi), L_1^(0)(O)}  (arity-zero kernel of (A_perp, A_perp) on (c, A_V), symmetrised, then the derivation).
Finite part psi with L_1^(0) = 0 on the family {O, F1} (Def c10-hyp (ii), Rem c5-cocycle), two-field sector:
  psi(O) = kappa (w1,w2) int omega (A_V, A_V)  with  L_1^(0)(O) + s psi(O) = 0;
  psi(F1) = kappa' F3                           with  L_1^(0)(F1) + s psi(F1) = 0;
the component changes by [{S^(1)(xi),-}, psi] (Prop c8-laws (b)): Q_1(O) -> Q_1(O) + 2 kappa F3 - kappa' F3, so
a3 -> a3' = a3 + 2 kappa - kappa', which vanishes exactly when G-R1 holds (the shifted Q_1(O) is s-closed)."""
from fractions import Fraction
import numpy as np
import qa_s1 as S1, qa_core as QC, qa_util as U

SB = S1.SB


def run(xi, a3, npts=3, seed=77, log=print):
    R = U.Rng(seed)
    g0 = QC.zero_gen()
    rows = []
    kappas, kappas_p = [], []
    for ip in range(npts):
        while True:
            w1, w2 = R.gv(6), R.gv(6)
            if np.dot(w1, w2) != 0:
                break
        tw = complex(np.dot(w1, w2))
        qc, ql = R.mom(), R.mom()
        l = R.vec('lamM')
        cvec = np.array([1 + 0j])
        Gam = np.array([xi @ SB[m] @ l for m in range(4)])            # Gamma(xi, l)|_V, corrected sign +
        # A = s(a3 F3) on (c, lam): F3 with A_V -> -dc = -i qc c
        z = U.gq(tw * sum(-1j * qc[m] * np.dot(SB[m] @ xi, l) for m in range(4)))
        A = U.vscale(a3, z)
        # B = L_1^(0)(F1), F1 = sum over (wa, wb) of (iota_wa A, (rho(wb) xi, lam))
        b1, b2 = U.vzero(), U.vzero()
        legs = [('c', cvec, qc), ('lamM', l, ql)]
        for wa, wb in ((w1, w2), (w2, w1)):
            zl = sum(wb[i] * (SB[4 + i] @ xi) for i in range(6))
            Xa, Xb = ('AW', wa, None), ('lamP', zl, None)
            b1 = U.vadd(b1, QC.pair_value(Xa, Xb, legs, g0)['fin'])
            b2 = U.vadd(b2, QC.pair_value(Xa, Xb, legs[::-1], g0)['fin'])
        B = U.vsub(b1, b2)
        # L_1^(0)(O) on (c, A_V) with polarisation pol: U(pol) = u(c,A_V) + u(A_V,c)
        Xa, Xb = ('AW', w1, None), ('AW', w2, None)

        def Ufun(pol):
            u1 = QC.pair_value(Xa, Xb, [('c', cvec, qc), ('AV', pol, ql)], g0)['fin']
            u2 = QC.pair_value(Xa, Xb, [('AV', pol, ql), ('c', cvec, qc)], g0)['fin']
            return U.vadd(u1, u2)
        C = U.vneg(Ufun(Gam))
        tot = U.vadd(A, U.vadd(B, C))
        # psi: kappa from L1(O) on a generic A_V polarisation, against s psi_O = -2 (w1,w2) (i qc . pol) on (c, A_V)
        pol = np.array([1, 2, -1, 1], dtype=complex)
        Ug = Ufun(pol)
        spsiO = U.gq(-2 * tw * 1j * np.dot(qc, pol))
        kap = U.vneg(U.gdiv(Ug, spsiO))
        # check that L1(O) is proportional to (dc, A_V) and not to c d.A_V: second polarisation
        pol2 = np.array([0, 1, 3, -2], dtype=complex)
        Ug2 = Ufun(pol2)
        kap2 = U.vneg(U.gdiv(Ug2, U.gq(-2 * tw * 1j * np.dot(qc, pol2)))) if np.dot(qc, pol2) != 0 else kap
        # kappa': L1(F1) + kappa' sF3 = 0 on (c, lam)
        kapp = U.vneg(U.gdiv(B, z))
        kappas.append(U.vstr(kap)); kappas_p.append(U.vstr(kapp))
        a3p = U.vadd(a3, U.vsub(U.vscale(kap, (Fraction(2), Fraction(0))), kapp))
        row = {'point': ip, 'w1.w2': str(tw), 'q_c': [int(x) for x in qc], 'q_lam': [int(x) for x in ql],
               'A = s(a3 F3)': U.vstr(A), 'B = L1(F1) = b(c,lam) - b(lam,c)': U.vstr(B),
               'C = {S, L1(O)} = -U(c,A_V -> +Gamma)': U.vstr(C), 'A+B+C': U.vstr(tot), 'pass': U.viszero(tot),
               'kappa (psi(O) = kappa (w1,w2) int omega (A_V,A_V))': U.vstr(kap),
               'kappa from a second A_V polarisation': U.vstr(kap2), 'L1(O) proportional to (dc,A_V)': U.vstr(kap) == U.vstr(kap2),
               "kappa' (psi(F1) = kappa' F3)": U.vstr(kapp), "a3' = a3 + 2 kappa - kappa'": U.vstr(a3p)}
        rows.append(row)
        log('  G-R1 point %d: A=%s | B=%s | C=%s | A+B+C=%s | kappa=%s kappa\'=%s a3\'=%s' % (
            ip, row['A = s(a3 F3)'], row['B = L1(F1) = b(c,lam) - b(lam,c)'], row['C = {S, L1(O)} = -U(c,A_V -> +Gamma)'],
            row['A+B+C'], row['kappa (psi(O) = kappa (w1,w2) int omega (A_V,A_V))'], row["kappa' (psi(F1) = kappa' F3)"],
            row["a3' = a3 + 2 kappa - kappa'"]))
    summary = {'all_pass': all(r['pass'] for r in rows),
               'kappa_constant': len(set(kappas)) == 1, 'kappa': kappas[0],
               'kappa_prime_constant': len(set(kappas_p)) == 1, 'kappa_prime': kappas_p[0],
               'L1O_is_dc_AV': all(r['L1(O) proportional to (dc,A_V)'] for r in rows),
               'a3_shifted': rows[0]["a3' = a3 + 2 kappa - kappa'"],
               'a3_shifted_zero_everywhere': all(r["a3' = a3 + 2 kappa - kappa'"] == '0' for r in rows)}
    return rows, summary
