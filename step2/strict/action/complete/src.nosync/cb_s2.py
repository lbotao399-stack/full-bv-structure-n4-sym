"""Route B: the arity-two vertex S^(2)(Q1,Q2) of Definition c2-S2 (constant spinors Q1, Q2 in Sigma, or the values of
conformal Killing spinors at the injected position), in engine B's realisation:
   S^(2)(Q1,Q2) = int ( 1/4 (Gamma(Q1,Q2), Gamma(lam*,lam*)) - 1/2 (Q1,lam*)(Q2,lam*) - (iota_{Gamma(Q1,Q2)} A, c*) ),
Gamma(Q1,Q2)^M = Q1^T MP R_M Q2 (Sigma), Gamma(w1,w2)^M = w1^T RB_M^T MP w2 (Sigma^*), (Q,w) = Q^T MP w, colour pairing
(c0,c1); iota_v A = sum_M v^M A_M.  Two forms: numpy density (tree level) and vertex TERMS (loop letter)."""
import numpy as np
import sympy as sp
import cb_common as B
Q, V, C = B.Q, B.V, B.C
R, RB, MP = B.R, B.RB, B.MP


def gam(Q1, Q2):
    return np.array([Q1 @ MP @ R[M] @ Q2 for M in range(10)])


def gam_star(w1, w2):
    return np.array([w1 @ RB[M].T @ MP @ w2 for M in range(10)])


def s2_dens(Q1, Q2, parts=('ll', 'Ac')):
    G12 = gam(Q1, Q2)
    terms = []
    if 'll' in parts:
        L = lambda s: 0.25 * (G12 @ gam_star(s[0]['u'], s[1]['u'])) - 0.5 * (Q1 @ MP @ s[0]['u']) * (Q2 @ MP @ s[1]['u'])
        terms.append((1, Q.dens(1, ('lams', 'lams'), L)))
    if 'Ac' in parts:
        terms.append((1, Q.dens(-1, ('A', 'cs'), lambda s: (G12 @ s[0]['u']) * s[1]['u'][0])))
    return Q.dsum(*terms)


def translation_dens(b, sectors=('V', 'perp', 'lam')):
    """S(P_b) two-leg sector parts: (b.dA, A*) [V: A_V and c; perp: A_perp], -(b.d lam, lam*), -(b.d c, c*)."""
    bd = lambda leg: 1j * np.dot(Q.pv(leg), b)
    out = []
    if 'V' in sectors:
        out.append((1, Q.dens(1, ('A', 'As'), lambda s: bd(s[0]) * np.dot(s[0]['u'][:4], s[1]['u'][:4]))))
        out.append((1, Q.dens(-1, ('c', 'cs'), lambda s: bd(s[0]) * s[0]['u'][0] * s[1]['u'][0])))
    if 'perp' in sectors:
        out.append((1, Q.dens(1, ('A', 'As'), lambda s: bd(s[0]) * np.dot(s[0]['u'][4:], s[1]['u'][4:]))))
    if 'lam' in sectors:
        out.append((1, Q.dens(-1, ('lam', 'lams'), lambda s: bd(s[0]) * (s[0]['u'] @ MP @ s[1]['u']))))
    return Q.dsum(*out)


def S2xi_dens(xi):
    T = Q.xi_terms(xi)
    return Q.dsum(*[(1, T[k]) for k in ('T1', 'T2', 'T3', 'T4')])


def _sym(v):
    return [sp.nsimplify(complex(x).real) + sp.I * sp.nsimplify(complex(x).imag) for x in v]


def s2_terms(Q1, Q2, tag, parts=('ll', 'Ac')):
    """vertex TERMS (vertices.py convention) of S^(2)(Q1,Q2) at numeric spinors; returns list of names."""
    Q1s, Q2s = sp.Matrix(_sym(Q1)), sp.Matrix(_sym(Q2))
    Mp = C.Mp
    G12 = [(Q1s.T * Mp * C.Rm[M] * Q2s)[0] for M in range(10)]
    names = []
    if 'll' in parts:
        def L_ll(s):
            w1, w2 = sp.Matrix(s[0]['u']), sp.Matrix(s[1]['u'])
            gs = sum((G12[M] * (w1.T * C.Rbm[M].T * Mp * w2)[0] for M in range(10)), 0)
            return sp.Rational(1, 4) * gs - sp.Rational(1, 2) * (Q1s.T * Mp * w1)[0] * (Q2s.T * Mp * w2)[0]
        V.TERMS['S2ll_' + tag] = (sp.Integer(1), ('lams', 'lams'), L_ll, ('k', (0, 1))); names.append('S2ll_' + tag)
    if 'Ac' in parts:
        def L_Ac(s):
            return sum((G12[M] * s[0]['u'][M] for M in range(10)), 0) * s[1]['u'][0]
        V.TERMS['S2Ac_' + tag] = (sp.Integer(-1), ('A', 'cs'), L_Ac, ('k', (0, 1))); names.append('S2Ac_' + tag)
    return names
