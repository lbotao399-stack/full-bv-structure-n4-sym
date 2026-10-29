"""Audit #3: intermediate equalities in the proof of delta_theory prop:d-psi (ii), exact rational arithmetic in the
finite model A of step2/theory/delta_toy.py (imported read-only; nothing written there).

The proof displays
   Q P*MP* - P*MP* Q = k M P* - P* M k = Q^GF (k M f - f M k),
   k M f - f M k = -4 (k H f' + k f - f H k') = -4 (k (H f)' - f H k') = 0 .
With [M,Q^GF] = -4 Q^GF and [M,g(H)] = -4 H g'(H) one finds instead
   k M P* - P* M k = Q^GF (k M f - f M k - 4 k f),     k M f - f M k = -4 (k H f' - f H k'),
so that the sum is -4 Q^GF (k f + k H f' - f H k') = -4 Q^GF (k (Hf)' - f H k') = 0 (Hf = k).
The conclusion [Q, P*MP*] = 0 is right; the two displayed intermediate equalities are each false, the errors cancel.
"""
import json, os, sys
from fractions import Fraction as Fr
from math import factorial
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory')
import delta_toy as T


def poly_Z(coeffs, Z, n):
    """sum_j coeffs[j] Z^j as a matrix"""
    out = T.mzero(n); P = T.mid(n)
    for cj in coeffs:
        out = T.madd(out, P, 1, cj); P = T.mmul(P, Z)
    return out


def main():
    m = T.model_A(T.FEYN)
    N, eps, Tt = m['N'], m['eps'], m['T']
    nz = N + 1
    Z = T.mzero(nz)
    for a in range(N):
        Z[a + 1][a] = Fr(1)
    # functions of H = 1 x Z as power series truncated at Z^N
    cf = lambda n: (Fr(-1) ** n) * (eps ** n - Tt ** n) / factorial(n)          # coefficient of z^n in e^{-eps z} - e^{-T z}
    k_c = [cf(n) for n in range(nz)]                                             # k(z)
    f_c = [cf(n + 1) for n in range(nz)]                                         # f(z) = k(z)/z
    fp_c = [cf(n + 2) * (n + 1) for n in range(nz)]                              # f'(z)
    kp_c = [cf(n + 1) * (n + 1) for n in range(nz)]                              # k'(z)
    I6 = T.mid(6)
    lift = lambda c: T.kron(I6, poly_Z(c, Z, nz))
    k, f, fp, kp = lift(k_c), lift(f_c), lift(fp_c), lift(kp_c)
    H, M, G, Q, Ps = m['H'], m['M'], m['G'], m['Q'], m['Pstar']
    mm, ad = T.mmul, T.madd
    res = {}
    res['|P* - Q^GF f|'] = T.mmax(ad(Ps, mm(G, f), 1, -1))
    res['|k - H f|'] = T.mmax(ad(k, mm(H, f), 1, -1))
    res['|[M,f] + 4 H f\'|'] = T.mmax(ad(ad(mm(M, f), mm(f, M), 1, -1), mm(H, fp), 1, 4))
    res['|[M,Q^GF] + 4 Q^GF|'] = T.mmax(ad(ad(mm(M, G), mm(G, M), 1, -1), G, 1, 4))
    PMP = mm(mm(Ps, M), Ps)
    lhs = ad(mm(Q, PMP), mm(PMP, Q), 1, -1)
    res['|Q P*MP* - P*MP* Q|'] = T.mmax(lhs)
    X1 = ad(mm(mm(k, M), Ps), mm(mm(Ps, M), k), 1, -1)                         # k M P* - P* M k
    kMf_fMk = ad(mm(mm(k, M), f), mm(mm(f, M), k), 1, -1)
    X2 = mm(G, kMf_fMk)                                                          # as displayed
    X3 = mm(G, ad(kMf_fMk, mm(k, f), 1, -4))                                     # corrected
    res['|Q P*MP* - P*MP* Q - (k M P* - P* M k)|'] = T.mmax(ad(lhs, X1, 1, -1))
    res['|k M P* - P* M k - Q^GF(kMf - fMk)|  (displayed, should be 0 if right)'] = T.mmax(ad(X1, X2, 1, -1))
    res['|k M P* - P* M k - Q^GF(kMf - fMk - 4kf)|  (corrected)'] = T.mmax(ad(X1, X3, 1, -1))
    Y2 = ad(ad(mm(k, mm(H, fp)), mm(k, f)), mm(f, mm(H, kp)), 1, -1)             # kHf' + kf - fHk'
    Y3 = ad(mm(k, mm(H, fp)), mm(f, mm(H, kp)), 1, -1)                           # kHf' - fHk'
    res['|kMf - fMk + 4(kHf\' + kf - fHk\')|  (displayed)'] = T.mmax(ad(kMf_fMk, Y2, 1, 4))
    res['|kMf - fMk + 4(kHf\' - fHk\')|  (corrected)'] = T.mmax(ad(kMf_fMk, Y3, 1, 4))
    res['|Q^GF (k f + kHf\' - fHk\')|  (the quantity that vanishes)'] = T.mmax(mm(G, Y2))
    res['|Q^GF k f| (size of the compensating term)'] = T.mmax(mm(G, mm(k, f)))
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'check_psi_proof_steps.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
