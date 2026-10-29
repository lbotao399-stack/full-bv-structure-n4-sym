"""Engine A, gate G8f: a fermion family.  The su(1|1) sector {Z, psi}: Z = A_{I=5} + i A_{I=6} (engine directions 8, 9),
psi the component of lambda_+ with weights (1/2, 1/2, 1/2) under the three rotation planes (4,5), (6,7), (8,9) of W
and weight +1/2 or -1/2 under the planes (0,1), (2,3) of V (exact null vector of (1/2) s[M] sb[N] - t, sympy; weights +i/2, -i/2 on (0,1), (2,3); chirality checked).
Planar two-site kernels (units u, Lemma lem:tl-colour, adjacent pair -> (-1)^r (value/2) Y_1 ... Y_r):
    K(Z,psi) -> alpha (Z psi) + beta (psi Z),   K(psi,Z) -> alpha' (psi Z) + beta' (Z psi),
    K(Z,Z) -> (Z Z) coefficient,  K(psi,psi) -> (psi psi) coefficient.
Beisert's su(1|1) density H_12 = 1 - Pi (graded permutation) predicts, with lc^natural = sum_pairs K - u sum n_F
(n_F(Z) = 2, n_F(psi) = 4, half of each letter's n_F to each of its two pairs), lc^natural|_pair = -2u H_12:
    (Z,psi): alpha - 3 = -2, beta = 2;   (Z,Z): coefficient - 2 = 0;   (psi,psi): coefficient - 4 = -2 * 2 = -4 -> 0."""
import json, os, sys
from fractions import Fraction
import numpy as np
import sympy as sp
import ea_core as E, ea_ops as O

SRC = os.path.dirname(os.path.abspath(__file__))


def psi_vector(sv=+1):
    S = lambda M, N: sp.Matrix(E.SG[M]) * sp.Matrix(E.SGB[N]) / 2
    to = lambda X: sp.Matrix(X.shape[0], X.shape[1], lambda i, j: sp.nsimplify(complex(X[i, j]).real) + sp.I * sp.nsimplify(complex(X[i, j]).imag))
    planes = [(0, 1), (2, 3), (4, 5), (6, 7), (8, 9)]
    target = [sv * sp.I / 2, -sv * sp.I / 2, -sp.I / 2, -sp.I / 2, -sp.I / 2]
    A = sp.Matrix.vstack(*[to(np.asarray(S(*pl)).astype(complex)).T - t * sp.eye(16) for pl, t in zip(planes, target)])
    ns = A.nullspace()
    return ns


def as_gauss(v):
    v = [sp.nsimplify(x) for x in v]
    den = sp.ilcm(*[sp.fraction(sp.re(x))[1] for x in v] + [sp.fraction(sp.im(x))[1] for x in v])
    w = [complex(sp.N(x * den)) for x in v]
    return np.array([complex(round(z.real), round(z.imag)) for z in w])


def coeff_matrix(ta, tb, za, zb, types, moms):
    """planar kernel (units u) on two ordered outputs, as a matrix over unit output polarisations"""
    n0, n1 = (O.slotvec_to_stripped(types[0], np.ones(1)).size if False else None), None
    rows0 = {'AW': 6, 'lamP': 16}[types[0]]; rows1 = {'AW': 6, 'lamP': 16}[types[1]]
    km = np.zeros((rows0, rows1), dtype=object)
    for i in range(rows0):
        for j in range(rows1):
            p0 = np.zeros(rows0, complex); p0[i] = 1
            p1 = np.zeros(rows1, complex); p1[j] = 1
            legs = [(types[0], p0, moms[0]), (types[1], p1, moms[1])]
            v = O.pair_kernel_on_legs((ta, za, O.ONE), (tb, zb, O.ONE), legs)
            km[i, j] = (v[0] / 2, v[1] / 2)            # (-1)^2 / 2
    return km


def project(km, a, b):
    """km == c * a (x) b exactly?  returns c (Gaussian rational) or None"""
    c = None
    for i in range(km.shape[0]):
        for j in range(km.shape[1]):
            w = complex(a[i]) * complex(b[j])
            val = complex(float(km[i, j][0]), float(km[i, j][1]))
            if abs(w) > 0 and c is None:
                c = val / w
    ok = True
    for i in range(km.shape[0]):
        for j in range(km.shape[1]):
            val = complex(float(km[i, j][0]), float(km[i, j][1]))
            if abs(val - c * complex(a[i]) * complex(b[j])) > 1e-12:
                ok = False
    return c, ok


def run():
    ns = psi_vector(+1)
    assert len(ns) >= 1
    zeta = as_gauss(ns[0])
    chir = np.asarray(E.G5S) @ zeta
    assert np.allclose(chir, zeta), 'psi must lie in S_+ (x) W (lambda_+)'
    z = np.zeros(6, complex); z[4] = 1; z[5] = 1j              # Z = A_8 + i A_9 (engine directions 8, 9)
    r1, r2 = np.array([1, 2, -1, 3]), np.array([2, -1, 3, 1])
    out = {'psi_nullspace_dim': len(ns), 'psi_vector': [str(x) for x in zeta]}
    k_A = coeff_matrix('AW', 'lamP', z, zeta, ('AW', 'lamP'), (r1, r2))
    k_S = coeff_matrix('AW', 'lamP', z, zeta, ('lamP', 'AW'), (r1, r2))
    a, oka = project(k_A, z, zeta); b, okb = project(k_S, zeta, z)
    kZZ = coeff_matrix('AW', 'AW', z, z, ('AW', 'AW'), (r1, r2)); cZZ, okZZ = project(kZZ, z, z)
    kPP = coeff_matrix('lamP', 'lamP', zeta, zeta, ('lamP', 'lamP'), (r1, r2)); cPP, okPP = project(kPP, zeta, zeta)
    out.update({'K(Z,psi)->(Z psi)': str(a), 'K(Z,psi)->(psi Z)': str(b), 'rank-one': [oka, okb],
                'K(Z,Z)->(Z Z)': str(cZZ), 'K(psi,psi)->(psi psi)': str(cPP), 'rank-one ZZ, psi psi': [okZZ, okPP]})
    lcZp = (a - 3, b); lcZZ = cZZ - 2; lcPP = cPP - 4
    out['lc_natural_pair'] = {'(Z,psi)': '%s (Z psi) + %s (psi Z)' % (lcZp[0], lcZp[1]), '(Z,Z)': str(lcZZ),
                              '(psi,psi)': str(lcPP)}
    out['target_-2u(1-Pi)'] = {'(Z,psi)': '-2 (Z psi) + 2 (psi Z)', '(Z,Z)': '0', '(psi,psi)': '-4'}
    out['pass'] = bool(oka and okb and abs(lcZp[0] + 2) < 1e-12 and abs(lcZp[1] - 2) < 1e-12 and abs(lcZZ) < 1e-12
                       and abs(lcPP + 4) < 1e-12)
    return out


if __name__ == '__main__':
    r = run()
    json.dump(r, open(os.path.join(SRC, '..', 'out.nosync', 'su11_gate.json'), 'w'), indent=1, default=str)
    print(json.dumps(r, indent=1, default=str))
