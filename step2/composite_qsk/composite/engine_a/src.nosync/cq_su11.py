"""Engine A, part C, family F2 (su(1|1)): two-site densities of the composite kernels on {Z, psi}.
Z = A_4 + i A_5 (letter zZ, leg eZ), psi the component of lambda_+ with SO(6) weights (1/2,1/2,1/2) and SO(4) weights
(+i/2, -i/2) on the planes (0,1), (2,3) (ea_su11.psi_vector(+1)).  On-shell legs: complex null Gaussian momenta q,
the psi leg polarisation u(q) in span{psi(+1), psi(-1)} solving the Weyl equation rho(q) u = 0 (sb[mu] q_mu u = 0).
Densities (units u, adjacent pair, planar factor 1/2): for an ordered pair (X, Y) and an ordered output word (Y1, Y2),
   d0[(X,Y)][word]  = translation-invariant part (constant),
   dP[(X,Y)][word]  = c1 (i beta.q_{Y1}) + c2 (i beta.q_{Y2})   (translation kernel along beta, on shell),
normalised by the letter-on-leg values of the output word (zZ.eZ = 2, zeta.u(q))."""
import itertools, time
from fractions import Fraction
import numpy as np
import sympy as sp
import cq_gval as GV, cq_eval as V
import ea_core as E
import ea_su11 as S11

ONE = {(0, 0, 0, 0): 1}
LOGSYM = [sp.Integer(1), sp.Symbol('log2'), sp.Symbol('log3'), sp.Symbol('log5'), sp.Symbol('log7')]
zZ = np.zeros(6, complex); zZ[4] = 1; zZ[5] = 1j          # Z = A_8 + i A_9 (as ea_su11)
eZ = np.zeros(6, complex); eZ[4] = 1; eZ[5] = -1j


def psis():
    p = S11.as_gauss(S11.psi_vector(+1)[0]); m = S11.as_gauss(S11.psi_vector(-1)[0])
    assert np.allclose(np.asarray(E.G5S) @ p, p) and np.allclose(np.asarray(E.G5S) @ m, m)
    return p, m


PSI, PSIM = psis()


_UCACHE = {}


def onshell_u(q, zeta=None):
    """an on-shell lambda_+ polarisation: u in S_+ (x) W (ea_core BASIS lamP) with (sum_mu q_mu sb[mu]) u = 0 and
    zeta.u != 0 (zeta the psi letter vector); exact Gaussian-integer vector"""
    zeta = PSI if zeta is None else zeta
    key = (tuple(complex(x) for x in q), tuple(complex(x) for x in zeta))
    if key in _UCACHE:
        return _UCACHE[key]
    W = sum(complex(q[mu]) * np.asarray(E.SGB[mu]) for mu in range(4))
    B = np.column_stack(E.BASIS['lamP']).astype(complex)
    M = W @ B
    Ms = sp.Matrix(M.shape[0], M.shape[1], lambda i, j: sp.Integer(int(round(M[i, j].real))) + sp.I * int(round(M[i, j].imag)))
    ns = Ms.nullspace()
    best = None
    for v in ns:
        den = sp.ilcm(*[sp.fraction(sp.nsimplify(sp.re(x)))[1] for x in v] + [sp.fraction(sp.nsimplify(sp.im(x)))[1] for x in v])
        x = np.array([complex(sp.N(c * den)) for c in v])
        u = B @ x
        u = np.array([complex(round(z.real), round(z.imag)) for z in u])
        if abs(complex(np.dot(zeta, u))) > 1e-9:
            best = u; break
    if best is None:
        _UCACHE[key] = None
        return None
    assert np.allclose(W @ best, 0) and np.allclose(np.asarray(E.G5S) @ best, best)
    _UCACHE[key] = best
    return best


LET = {'Z': ('AW', zZ, ONE), 'psi': ('lamP', PSI, ONE)}


def leg(name, q):
    if name == 'Z':
        return ('AW', eZ, q)
    return ('lamP', onshell_u(q), q)


def lol(name, q):
    """letter-on-leg value of an output letter"""
    if name == 'Z':
        return complex(np.dot(zZ, eZ))
    u = onshell_u(q)
    return 0j if u is None else complex(np.dot(PSI, u))


NULLS = [np.array(v) for v in [(1, 1j, 0, 0), (0, 0, 1, 1j), (1, 0, 1j, 0), (1 + 1j, 1 - 1j, 0, 0), (1, 1, 1j, 1j),
                               (1, 1j, 1, -1j), (2, 2j, 1, 1j), (1, -1j, 1j, 1), (0, 0, 1, -1j), (1, 0, -1j, 0),
                               (0, 1, 1j, 0), (1, -1j, 1, 1j), (2, 1, 2j, 1j), (1, 2, 2j, 1j)]]


def check_nulls():
    return all(abs(complex(np.dot(q, q))) < 1e-12 for q in NULLS)


WORDS = {('Z', 'psi'): [('Z', 'psi'), ('psi', 'Z')], ('psi', 'Z'): [('psi', 'Z'), ('Z', 'psi')],
         ('Z', 'Z'): [('Z', 'Z')], ('psi', 'psi'): [('psi', 'psi')]}


def densities(beta=(1, 2, -1, 3), kappa=-1, npts=4, log=print):
    out = {}
    beta = np.array(beta)
    good = [q for q in NULLS if abs(lol('psi', q)) > 1e-9]
    kin = [(good[i], good[j]) for i, j in itertools.combinations(range(len(good)), 2)][:npts + 2]
    for pair, words in WORDS.items():
        Xa, Xb = LET[pair[0]], LET[pair[1]]
        out[str(pair)] = {}
        for w in words:
            d0s = []; rows = []; rhs = []; vvnz = False
            for (q1, q2) in kin:
                legs = [leg(w[0], q1), leg(w[1], q2)]
                norm = lol(w[0], q1) * lol(w[1], q2)
                a, b = GV.kernel(Xa, Xb, legs, 'D0')
                vvnz |= not V.viszero(b)
                tot = GV.combine((a, b), kappa)
                d0s.append(sp.nsimplify(0) + sum(_sym(tot[i]) * LOGSYM[i] for i in range(V.NB)) / 2 / _symc(norm))
                a, b = GV.kernel(Xa, Xb, legs, ('Pv', tuple(beta)))
                vvnz |= not V.viszero(b)
                tot = GV.combine((a, b), kappa)
                val = sum(_sym(tot[i]) * LOGSYM[i] for i in range(V.NB)) / 2 / _symc(norm)
                rows.append([_symc(1j * np.dot(beta, q1)), _symc(1j * np.dot(beta, q2))]); rhs.append(sp.expand(val))
            d0s = [sp.simplify(x) for x in d0s]
            const = all(sp.simplify(x - d0s[0]) == 0 for x in d0s)
            A = sp.Matrix(rows); Y = sp.Matrix(rhs)
            sol = (A.H * A).inv() * (A.H * Y)
            sol = sp.simplify(sol)
            resid = sp.simplify(A * sol - Y)
            ok = all(x == 0 for x in resid)
            out[str(pair)][str(w)] = {'d0': str(d0s[0]), 'd0_constant': const, 'dP_c1': str(sol[0]), 'dP_c2': str(sol[1]),
                                      'dP_fit_exact': ok, 'vertex_vertex_nonzero': vvnz}
            log('   su11 pair %s word %s: d0 = %s (const %s), dP = %s q_1 + %s q_2 (exact %s)' % (
                pair, w, d0s[0], const, sol[0], sol[1], ok))
    return out


def _sym(z):
    return sp.Rational(z[0]) + sp.I * sp.Rational(z[1])


def _symc(z):
    z = complex(z)
    return sp.nsimplify(round(z.real, 9)) + sp.I * sp.nsimplify(round(z.imag, 9))
