"""Engine A, part C, family F1: twist two tr(D_n^k Z D_n^{N-k} Z) (n = (1, i, 0, 0), Z = A_4 + i A_5 of W_R^perp in engine
directions), two-field part, gauge completion to O(A_n) (ea_spectra.covariant_twist_op).  Matrices (units u) in the
basis B_k = tr(d_n^k Z d_n^{N-k} Z), k = 0..N/2:
  M0   : sum C^(2),0_Delta (translation-invariant part)            P7(c): -2/3 (4 + N) * identity
  Mn   : sum C^(2)_{P_n} into spin N+1 (beta = n)                   P7(b): -1/3 * i n.(q1+q2) B_k
  Mbeta: sum C^(2)_{P_beta}, beta generic, fitted on the Lorentz ansatz (beta.q_i)(n.q)^N, (beta.n)(q.q)(n.q)^{N-1}
Fits are exact (sympy, Gaussian rationals per log-basis component) with a residual check on extra points."""
import itertools, math, time
from fractions import Fraction
import numpy as np
import sympy as sp
import cq_ops as CO, cq_eval as V
import ea_spectra as ES

zZ = ES.zZ; eZ = ES.eZ; NV = ES.NV
ONE = {(0, 0, 0, 0): 1}


def basis_val(k, N, legs):
    return CO.trace_on_legs([('AW', zZ, ES.nullpoly(k)), ('AW', zZ, ES.nullpoly(N - k))], legs)


def to_sym(z):
    return sp.Rational(z[0]) + sp.I * sp.Rational(z[1])


def fit(A, Ys):
    """exact least squares A M = Y for every log-basis component; returns (M components, max residual)"""
    As = sp.Matrix(A)
    AH = As.H
    inv = (AH * As).inv()
    out = []; res = 0
    for Y in Ys:
        Yc = sp.Matrix(Y)
        Mq = sp.simplify(inv * (AH * Yc))
        r = sp.simplify(As * Mq - Yc)
        res = max([res] + [abs(complex(x)) for x in r])
        out.append(Mq)
    return out, res


def run(N, npts_extra=2, seed=11, kinds=('D0', 'Pn'), kappa=-1, log=print):
    rng = np.random.default_rng(seed + 7 * N)
    ks = list(range(N // 2 + 1))
    ops = [ES.covariant_twist_op(k, N) for k in ks]
    cache = {}
    out = {'N': N}
    for kind in kinds:
        t0 = time.time()
        if kind == 'D0':
            nunk = len(ks)
            def bfun(q1, q2):
                return [complex(basis_val(j, N, [('AW', eZ, q1), ('AW', eZ, q2)])) for j in ks]
            kd = 'D0'
        elif kind == 'Pg':
            # generic beta: ansatz (i beta.q_i)(i n.q1)^a (i n.q2)^{N-a}; (beta.n) s(q) (i n.q1)^a (i n.q2)^{N-1-a}
            beta = np.array([1, 2, -1, 3])
            nunk = 2 * (N + 1) + 3 * N
            def bfun(q1, q2):
                a1 = 1j * complex(np.dot(NV, q1)); a2 = 1j * complex(np.dot(NV, q2))
                b1 = 1j * complex(np.dot(beta, q1)); b2 = 1j * complex(np.dot(beta, q2)); bn = complex(np.dot(beta, NV))
                f = [b * a1 ** a * a2 ** (N - a) for b in (b1, b2) for a in range(N + 1)]
                for sv in (complex(np.dot(q1, q2)), complex(np.dot(q1, q1)), complex(np.dot(q2, q2))):
                    f += [bn * sv * a1 ** a * a2 ** (N - 1 - a) for a in range(N)]
                return f
            kd = ('Pv', tuple(beta))
        elif kind == 'Pn':
            # image in spin N+1: basis (i n.q1)^a (i n.q2)^{N+1-a}, a = 0..N+1
            nunk = N + 2
            def bfun(q1, q2):
                a1 = 1j * complex(np.dot(NV, q1)); a2 = 1j * complex(np.dot(NV, q2))
                return [a1 ** a * a2 ** (N + 1 - a) for a in range(N + 2)]
            kd = ('Pv', tuple(NV))
        npts = nunk + npts_extra
        pts = []
        while len(pts) < npts:
            pts.append((rng.integers(-3, 4, 4), rng.integers(-3, 4, 4)))
        A = []; Y = [[[] for _ in range(V.NB)] for _ in ks]
        vv_nonzero = [False]
        for (q1, q2) in pts:
            A.append([sp.nsimplify(complex(x).real) + sp.I * sp.nsimplify(complex(x).imag) for x in bfun(q1, q2)])
            legs = [('AW', eZ, q1), ('AW', eZ, q2)]
            for ik, op in enumerate(ops):
                tl, tv = CO.planar_on_legs(op, legs, kd, cache=cache)
                if not V.viszero(tv):
                    vv_nonzero[0] = True
                tot = V.vadd(tl, V.vscale(tv, (Fraction(kappa), Fraction(0))))
                for b in range(V.NB):
                    Y[ik][b].append(to_sym(tot[b]))
        # fit per operator, per basis component
        Ms = []; res = 0
        for ik in range(len(ks)):
            comps, r = fit(A, [[[y] for y in Y[ik][b]] for b in range(V.NB)])
            res = max(res, r)
            Ms.append(comps)
        # matrix: column ik = image of O_k ; rows = basis functions
        mats = []
        for b in range(V.NB):
            mats.append([[Ms[ik][b][j] for ik in range(len(ks))] for j in range(nunk)])
        out[kind] = {'matrix_by_logbasis': {V.Q.BASIS_NAMES[b]: [[str(x) for x in row] for row in mats[b]] for b in range(V.NB)},
                     'fit_residual': float(res), 'points': npts, 'seconds': round(time.time() - t0, 1),
                     'vertex_vertex_part_nonzero': vv_nonzero[0]}
        log('  twist two N=%d %s done %.1fs residual %g' % (N, kind, time.time() - t0, res))
    return out
