"""Audit 2b: large-|y| order of a single insertion of S^(1)(x) on a cycle (ch08 Prop c8-Tsafe; theory issue CH08 #1).

The insertion point y carries the term of S^(1)(x) with both legs on propagator components: (L_v A + (1/4) div v A|_perp, A*)
for even x.  Its A leg meets an A x A edge G(y - z1), its A*|_V leg an A*|_V x c edge d G(y - z2) (Prop c3-P, eps -> 0,
T = infinity, G(r) = 1/(4 pi^2 |r|^2)).  Insertion kernel (A leg along W_perp, so L_v A = v(A) on components):
    K_nu(y) = [v(y).d_y G(y - z1) + w(y) G(y - z1)] * d_nu G(y - z2),
  Delta: v = 2y, w = 2 (weight term of Example c2-dil);  translation: v = b, w = 0;  special conformal: v = |y|^2 b - 2(b.y)y,
  w = (1/4) div v = -2 (b.y).
The decay exponent is measured exactly by rational arithmetic along rays y = R n at R = 10^3 .. 10^6 (fixed generic z1, z2,
b, n): log|K(2R)/K(R)|/log 2 -> -order.  Also the sphere average of the leading term is reported (Monte Carlo free:
exact average of the leading homogeneous term computed from its closed form).  Pure Python fractions.
"""
import json, os, math
from fractions import Fraction as Fr

z1 = [Fr(1, 3), Fr(-2, 5), Fr(1, 7), Fr(3, 4)]
z2 = [Fr(-1, 2), Fr(2, 3), Fr(-3, 5), Fr(1, 9)]
b = [Fr(2, 7), Fr(1, 3), Fr(-1, 5), Fr(3, 11)]
dirs = [[Fr(1), Fr(2), Fr(-2), Fr(1)], [Fr(3), Fr(-1), Fr(1), Fr(2)], [Fr(0), Fr(1), Fr(4), Fr(-1)]]

def sq(v):
    return sum(x * x for x in v)

def G(r):                      # without the 1/(4 pi^2)
    return 1 / sq(r)

def dG(r):                     # d_nu of 1/|r|^2 = -2 r_nu/|r|^4
    s = sq(r)
    return [-2 * x / (s * s) for x in r]

def kernel(y, case):
    r1 = [y[i] - z1[i] for i in range(4)]; r2 = [y[i] - z2[i] for i in range(4)]
    g1 = G(r1); dg1 = dG(r1); dg2 = dG(r2)
    if case == 'Delta':
        v = [2 * x for x in y]; w = Fr(2)
    elif case == 'translation':
        v = b[:]; w = Fr(0)
    else:
        yy = sq(y); by = sum(b[i] * y[i] for i in range(4))
        v = [yy * b[i] - 2 * by * y[i] for i in range(4)]; w = -2 * by
    A = sum(v[i] * dg1[i] for i in range(4)) + w * g1
    return [A * dg2[nu] for nu in range(4)]

out = {}
for case in ('translation', 'Delta', 'special conformal'):
    rec = []
    for n in dirs:
        exps = []
        for k in (3, 4, 5, 6):
            R = Fr(10) ** k
            K1 = kernel([R * x for x in n], case); K2 = kernel([2 * R * x for x in n], case)
            nrm1 = math.sqrt(float(sq(K1))); nrm2 = math.sqrt(float(sq(K2)))
            exps.append(round(math.log(nrm2 / nrm1) / math.log(2), 6))
        rec.append(exps)
    out[case] = {'measured exponent of |K| along three rays, R = 1e3..1e6': rec}
# leading terms in closed form (y = R n, |n| = 1):
#   Delta: [2 y.d(1/|y|^2) + 2/|y|^2] (-2 y_nu/|y|^4) = [(-4 + 2)/R^2](-2 n_nu/R^3) = 4 n_nu / R^5   (odd, |y|^-5)
#   translation: (-2 b.n/R^3)(-2 n_nu/R^3) = 4 (b.n) n_nu / R^6                                     (|y|^-6)
#   special conformal: [v.d(1/|y|^2) + w/|y|^2] = [-2(b.n)R^3 + 4 (b.n) R^3 - 2 (b.n) R]/R^4 ... leading (b.n)(2-2)/R = 0 at
#   order 1/R; with z1 the next order gives |y|^-2 for the A-leg factor, times |y|^-3: |y|^-5?  -> measured above.
out['closed-form leading terms (|n| = 1)'] = {
    'Delta': '4 n_nu / |y|^5 (odd in y; absolutely integrable at large |y| in R^4)',
    'translation': '4 (b.n) n_nu / |y|^6',
}
here = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(here, 'insertion_decay.json'), 'w') as f:
    json.dump(out, f, indent=1)
print(json.dumps(out, indent=1))
