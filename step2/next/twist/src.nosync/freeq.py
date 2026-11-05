"""Free twisted differential Q^xi_free = Q + xi. + (S^(2) term) on single letters, engine-A conventions.
Fields e = (c, A_0..9, A*_0..9, lam_1..16, lam*_1..16, c*) (54, single colour).  Letters = row vectors w (functional
w.e on a plane wave of momentum p); Q acts on letters by w -> w M(p).  M(p) = op_Q(p) + Xi(p) + kappa S2."""
import env
import numpy as np, sympy as sp, itertools
from fractions import Fraction
import qa_e0 as E0, qa_s1 as S1, qa_twist as TW, ea_core as E
NS = E.NS; O = E.OFS

def Xi_op(xi, p):
    gen = S1.quad_block(S1.odd_terms(xi, ('T1', 'T2', 'T3', 'T4')))
    Xb = E0.block_operator(gen, p)
    Xf = E0.op_xi_fields(xi, p)
    fld = [O['A'] + i for i in range(10)] + [O['l'] + a for a in range(16)] + [O['c']]
    sub_b = Xb[np.ix_(fld, fld)]; sub_f = Xf[np.ix_(fld, fld)]
    nz = np.abs(sub_f) > 1e-9
    rr = np.unique(np.round(sub_b[nz] / sub_f[nz], 9))
    assert len(rr) == 1, rr
    return Xb / complex(rr[0])

def S2_op(xi):
    M = np.zeros((NS, NS), complex)
    M[O['l']:O['l'] + 16, O['ls']:O['ls'] + 16] = np.outer(xi, xi)
    return M

def kappa_S2(xi):
    """unique kappa with (Q + Xi + kappa S2)^2 = 0 on plane waves"""
    rng = np.random.default_rng(1)
    ks = []
    for _ in range(3):
        p = rng.integers(-3, 4, 4).astype(complex)
        Qm = E0.op_Q(p); X = Xi_op(xi, p); S = S2_op(xi)
        A = X @ X; B = Qm @ S + S @ Qm
        assert np.allclose(Qm @ X + X @ Qm, 0)
        nz = np.abs(B) > 1e-9
        r = np.unique(np.round(-A[nz] / B[nz], 9))
        assert len(r) == 1, r
        assert np.allclose(A + r[0] * B, 0)
        ks.append(complex(r[0]))
    assert np.allclose(ks, ks[0])
    return ks[0]

def M_free(xi, p, kap):
    return E0.op_Q(p) + Xi_op(xi, p) + kap * S2_op(xi)

if __name__ == '__main__':
    T = TW.Twist()
    kap = kappa_S2(T.xi)
    print('kappa_S2 =', kap)
    p = np.array([1, -2, 3, 1], complex)
    Mm = M_free(T.xi, p, kap)
    print('nilpotent:', np.allclose(Mm @ Mm, 0))
