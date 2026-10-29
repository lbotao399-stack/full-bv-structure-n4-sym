"""Engine A, part C, family F2 (fermion sl(2)): tr(D_n^k psi D_n^{N-k} psi), one R component psi (the su(1|1) letter),
n a null direction aligned with psi, gauge completion to O(A_n) (cq_gops.covariant_op).  Two-field part, on shell
(null Gaussian momenta, Weyl-equation polarisations u(q) with psi.u(q) != 0).  Basis B_j = tr(d^j psi d^{N-j} psi),
j < N/2 (B_j = -B_{N-j}).  Matrices (units u):
  L    : lc^natural = sum_pairs K^log - u sum_letters n_F  (n_F(lambda) = 4)      target eig -8 h(N+1), N odd
  M0   : sum C^(2),0_Delta                                                          P7(c): -2/3 (4 + N)
  Pn   : sum C^(2)_{P_n} - (-1/3) i n.(q1+q2) (two-field part), residual            P7(b): 0
  D1   : 4 L + M0 + 2/3 (4 + N)  (theta = -1/3, Lemma lem:cc-onebody)                target eig -32 h(N+1)."""
import itertools, math, time
from fractions import Fraction
import numpy as np
import sympy as sp
import cq_gops as GO, cq_gval as GV, cq_eval as V, cq_su11 as S11, cq_ops as CO

LOGSYM = S11.LOGSYM


def nullpoly(n, k):
    """(i n.p)^k as {exp: complex}"""
    out = {(0, 0, 0, 0): 1}
    for _ in range(k):
        new = {}
        for e, c in out.items():
            for mu in range(4):
                if n[mu] == 0:
                    continue
                e2 = list(e); e2[mu] += 1; e2 = tuple(e2)
                new[e2] = new.get(e2, 0) + c * 1j * n[mu]
        out = new
    return {e: complex(round(c.real), round(c.imag)) for e, c in out.items() if abs(c) > 1e-12}


def letters(n):
    psi = lambda m: ('lamP', S11.PSI, nullpoly(n, m))
    An = lambda m: ('AV', np.asarray(n, dtype=complex), nullpoly(n, m))
    return psi, An


def kin_points(npts, seed=1):
    good = [q for q in S11.NULLS if abs(S11.lol('psi', q)) > 1e-9]
    pairs = list(itertools.combinations(range(len(good)), 2))
    rng = np.random.default_rng(seed)
    rng.shuffle(pairs)
    return [(good[i], good[j]) for i, j in pairs[:npts]]


def symv(v):
    return sum((sp.Rational(v[i][0]) + sp.I * sp.Rational(v[i][1])) * LOGSYM[i] for i in range(V.NB))


def symc(z):
    z = complex(z)
    return sp.nsimplify(round(z.real, 9)) + sp.I * sp.nsimplify(round(z.imag, 9))


def run(N, n=(1, 1j, 0, 0), kinds=('L', 'M0', 'Pn'), kappa=-1, log=print, npts_extra=2):
    psi, An = letters(n)
    ks = list(range((N + 1) // 2))            # j < N/2
    ops = [GO.covariant_op(psi, k, N, An) for k in ks]
    nb = len(ks)
    pts = kin_points(nb + npts_extra + 2)
    legsets = [[('lamP', S11.onshell_u(q1), q1), ('lamP', S11.onshell_u(q2), q2)] for q1, q2 in pts]
    A = sp.Matrix([[symc(GO.trace_graded([psi(j), psi(N - j)], lg)) for j in ks] for lg in legsets])
    out = {'N': N, 'n': [str(x) for x in n]}
    caches = {}
    for kind in kinds:
        t0 = time.time()
        Y = []
        for lg, (q1, q2) in zip(legsets, pts):
            row = []
            for ik, op in enumerate(ops):
                if kind == 'L':
                    val = GO.planar_graded(op, lg, lambda a, b, o: GV.logkernel(a, b, o), cache=caches.setdefault('L', {}))
                    y = symv(val) - 8 * symc(GO.op_graded([op[0]], lg))
                elif kind == 'M0':
                    a, b = GO.planar_graded(op, lg, lambda a, b, o: GV.kernel(a, b, o, 'D0'), cache=caches.setdefault('M0', {}))
                    y = symv(GV.combine((a, b), kappa))
                elif kind == 'Pn':
                    a, b = GO.planar_graded(op, lg, lambda a, b, o: GV.kernel(a, b, o, ('Pv', tuple(n))), cache=caches.setdefault('Pn', {}))
                    y = symv(GV.combine((a, b), kappa)) - sp.Rational(-1, 3) * symc(1j * np.dot(n, q1 + q2)) * symc(GO.op_graded([op[0]], lg))
                row.append(sp.expand(y))
            Y.append(row)
        Y = sp.Matrix(Y)
        if kind == 'Pn':
            res = max([abs(complex(sp.N(x.subs({s: 1 for s in LOGSYM[1:]})))) for x in Y] + [0])
            out[kind] = {'max_abs_residual_vs_-1/3_i_n.dq': res, 'pass': res == 0, 'seconds': round(time.time() - t0, 1)}
            log('  fermion sl2 N=%d Pn residual %s (%.1fs)' % (N, res, time.time() - t0))
            continue
        M = sp.simplify((A.H * A).inv() * (A.H * Y))
        resid = sp.simplify(A * M - Y)
        rmax = max([abs(complex(sp.N(x.subs({s: 1 for s in LOGSYM[1:]})))) for x in resid] + [0])
        ev = M.eigenvals()
        out[kind] = {'matrix': [[str(M[i, j]) for j in range(nb)] for i in range(nb)], 'fit_residual': rmax,
                     'eigenvalues': {str(k): v for k, v in ev.items()}, 'seconds': round(time.time() - t0, 1)}
        out.setdefault('_mats', {})[kind] = M
        log('  fermion sl2 N=%d %s: eig %s residual %s (%.1fs)' % (N, kind, ev, rmax, time.time() - t0))
    if 'L' in out.get('_mats', {}) and 'M0' in out['_mats']:
        D1 = 4 * out['_mats']['L'] + out['_mats']['M0'] + sp.Rational(2, 3) * (4 + N) * sp.eye(nb)
        ev = D1.eigenvals()
        target = sorted([-32 * sum(sp.Rational(1, j) for j in range(1, s + 2)) for s in range(1, N + 1, 2)])
        got = sorted([sp.nsimplify(e) for e, m in ev.items() for _ in range(m)], key=lambda z: float(sp.re(z)))
        out['D1'] = {'eigenvalues': [str(x) for x in got], 'target_-32h(n+1)_n_odd<=N': [str(x) for x in target],
                     'pass': [sp.simplify(a - b) == 0 for a, b in zip(got, target)] == [True] * len(target) and len(got) == len(target)}
        log('  fermion sl2 N=%d D1 eig %s target %s' % (N, got, target))
    out.pop('_mats', None)
    return out
