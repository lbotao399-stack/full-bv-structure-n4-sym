"""Diagnostic of the verification claimed by step2/theory/delta_toy.py (the theory note's toy model).

This file does not verify anything about the theory with delta_toy.py; it only measures how strongly the
toy's own tests probe the composite-line identities.  It imports delta_toy read-only (no file is
written by the import) and re-runs three of its identities, in its own conventions and on its own model A,
(1) on its own four sparse random test functionals and (2) on every monomial of degree <= 3:
  [d_P, C_P(S_M,-)] = d_Psi, Psi = (P*MP* x 1)K_0          (its check '[d_P,C_P(S_M,-)] - d_Psi')
  d_Psi Q - Q d_Psi = 0                                     (its check 'd_Psi Q - Q d_Psi')
  conjugation formula with (h^2/2) d_Psi                    (its check 'conjugation formula')
"""
import json, os, random, sys, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
THEORY = os.path.abspath(os.path.join(HERE, '..', '..', 'theory'))
sys.path.insert(0, THEORY)
import delta_toy as dt   # read-only import

def main():
    out = {}
    m = dt.model_A(dt.FEYN); alg = m['alg']; K0 = m['K0']
    Q, M, P, Ps = m['Q'], m['M'], m['P'], m['Pstar']
    hb = Fr(2, 7)
    rng = random.Random(1); cf = lambda r: Fr(r.randint(-3, 3) or 1, r.randint(1, 3))
    tests = [alg.rand([1, 2, 3], 5, rng, p, cf) for p in (0, 1, 0, 1)]   # exactly the toy's tests (seed 1)
    SM = dt.quad_for_derivation(alg, K0, M, 0)
    Qf = lambda f: alg.D(Q, f)
    PMP = dt.mmul(dt.mmul(Ps, M), Ps); Psi = dt.kernel_of(PMP, K0)
    Cp = lambda f: alg.C(P, SM, f)
    def comp(f):
        return alg.add(alg.dphi(P, Cp(f)), Cp(alg.dphi(P, f)), coeffs=[1, -1])
    # all monomials of degree <= 3
    n = alg.n; dense = []
    for d in (1, 2, 3):
        for mono in itertools.combinations_with_replacement(range(n), d):
            mm, s = alg.sortmono(mono)
            if mm is not None:
                dense.append({mm: Fr(1)})
    out['toy: number of variables'] = n
    out['toy: its test functionals, number of monomials each'] = [len(t) for t in tests]
    out['toy tests: max |d_Psi F| over its 4 test functionals'] = str(max(alg.dist(alg.dphi(Psi, F), {}) for F in tests))
    out['toy tests: max |[d_P,C_P(S_M,-)] F| over its 4 test functionals'] = str(max(alg.dist(comp(F), {}) for F in tests))
    out['dense (all monomials deg<=3): number'] = len(dense)
    out['dense: max |[d_P,C_P(S_M,-)]F - d_Psi F|  (toy claims 0)'] = str(max(alg.dist(comp(F), alg.dphi(Psi, F)) for F in dense))
    out['dense: max |[d_P,C_P(S_M,-)]F - 2 d_Psi F|'] = str(max(alg.dist(comp(F), alg.scale(2, alg.dphi(Psi, F))) for F in dense))
    out['dense: max |d_Psi Q F - Q d_Psi F|  (toy claims 0)'] = str(max(alg.dist(alg.dphi(Psi, Qf(F)), Qf(alg.dphi(Psi, F))) for F in dense))
    out['dense: max |d_Psi Q F + Q d_Psi F|'] = str(max(alg.dist(alg.add(alg.dphi(Psi, Qf(F)), Qf(alg.dphi(Psi, F))), {}) for F in dense))
    def conj(F, c):
        eF = alg.exp_dphi(P, hb, F)
        rhs = alg.add(alg.mul(SM, eF), alg.scale(hb, alg.mul(alg.dphi(P, SM), eF)), alg.scale(hb, alg.C(P, SM, eF)),
                      alg.scale(c * hb * hb / 2, alg.dphi(Psi, eF)))
        return alg.dist(alg.exp_dphi(P, hb, alg.mul(SM, F)), rhs)
    out['dense: conjugation formula residual with (h^2/2) d_Psi  (toy claims 0)'] = str(max(conj(F, 1) for F in dense))
    out['dense: conjugation formula residual with h^2 d_Psi'] = str(max(conj(F, 2) for F in dense))
    return out

if __name__ == '__main__':
    res = main()
    with open(os.path.join(HERE, 'diag_theory_toy_results.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
