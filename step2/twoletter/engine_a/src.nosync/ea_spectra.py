"""Engine A, task T2 / gate G8: planar spectra from the two-body log kernel.
(i) su(2) sector Tr(Z^{L-M} X^M): two-site matrix of the planar kernel computed by the engine, cyclic chain, exact
    characteristic polynomials of lc^natural = sum_l K_{l,l+1} - u sum_letters n_F, in units of -2u (E_H normalisation).
(ii) sl(2) sector tr(d_n^k Z d_n^{N-k} Z), n null (n = (1, i, 0, 0)): lc^natural on the symmetric polynomials,
    exact matrix, eigenvalues vs -8 h(S) units (B03c twist two)."""
import itertools, math, json, time
from fractions import Fraction
import numpy as np
import sympy as sp
import ea_core as E, ea_ops as O

zZ = np.array([1, 1j, 0, 0, 0, 0]); zX = np.array([0, 0, 1, 1j, 0, 0])
eZ = np.array([1, -1j, 0, 0, 0, 0]); eX = np.array([0, 0, 1, -1j, 0, 0])   # extractors: (zZ, eZ) = 2


def two_site_su2():
    """M[(o1,o2),(a,b)] planar kernel value in units u (normalised letters), states Z=0, X=1"""
    q = [np.array([1, 0, 2, -1]), np.array([0, 1, -1, 1])]
    Z = {0: zZ, 1: zX}; X = {0: eZ, 1: eX}
    M = [[Fraction(0)] * 4 for _ in range(4)]
    for a, b in itertools.product(range(2), repeat=2):
        for o1, o2 in itertools.product(range(2), repeat=2):
            legs = [('AW', X[o1], q[0]), ('AW', X[o2], q[1])]
            v = O.pair_kernel_on_legs(('AW', Z[a], O.ONE), ('AW', Z[b], O.ONE), legs)
            assert v[1] == 0
            M[2 * o1 + o2][2 * a + b] = v[0] / 2 / 4          # planar factor (-1)^2/2, extractor normalisation 2 x 2
    return M


def cyclic_spectrum(M2, L, Mmag):
    """lc^natural on cyclic states of length L with Mmag letters X, in units u: sum over links of the two-site matrix
    minus 2 per letter (n_F(A) = 2).  Returns (matrix on cyclic orbits, charpoly of E = -lc^natural/(2u))."""
    states = [s for s in itertools.product(range(2), repeat=L) if sum(s) == Mmag]
    def canon(s):
        return min(tuple(s[r:] + s[:r]) for r in range(L))
    orbits = sorted(set(canon(s) for s in states))
    # cyclic-invariant vectors: sum over the orbit; the matrix of H on orbit sums
    idx = {o: i for i, o in enumerate(orbits)}
    n = len(orbits)
    H = [[Fraction(0)] * n for _ in range(n)]
    for o in orbits:
        orb = set(tuple(o[r:] + o[:r]) for r in range(L))
        vec = {s: Fraction(1) for s in orb}
        out = {}
        for s, c in vec.items():
            for l in range(L):
                a, b = s[l], s[(l + 1) % L]
                for o1, o2 in itertools.product(range(2), repeat=2):
                    m = M2[2 * o1 + o2][2 * a + b]
                    if m == 0:
                        continue
                    t = list(s); t[l] = o1; t[(l + 1) % L] = o2; t = tuple(t)
                    out[t] = out.get(t, 0) + c * m
            out[s] = out.get(s, 0) - 2 * L * c
        # express out (cyclic invariant) in orbit sums: coefficient of the canonical rep
        for t, c in out.items():
            if t == canon(t):
                H[idx[t]][idx[o]] += c
    # E = lc^natural / (-2u)
    Em = sp.Matrix(n, n, lambda i, j: sp.Rational(H[i][j]) / (-2))
    x = sp.Symbol('E')
    cp = sp.factor(Em.charpoly(x).as_expr())
    return Em, cp


def harmonic(n):
    return sum((Fraction(1, k) for k in range(1, n + 1)), Fraction(0))


def nullpoly(k):
    """(i n.p)^k, n = (1, i, 0, 0): dict exp -> Gaussian integer"""
    out = {}
    for j in range(k + 1):
        c = (1j ** k) * math.comb(k, j) * (1j ** j)
        e = (k - j, j, 0, 0)
        out[e] = out.get(e, 0) + complex(round(c.real), round(c.imag))
    return {e: c for e, c in out.items() if c != 0}


NV = np.array([1, 1j, 0, 0])      # null direction n; A_n = (n, A|_V)


def covariant_twist_op(k, N):
    """tr(D_n^k Z D_n^{N-k} Z) to first order in A_n, D = d + [A, -]:
    D^k Z = d^k Z + sum_i C(k, i+1) [d^i A_n, d^{k-1-i} Z] + O(A^2)"""
    Zk = lambda m: ('AW', zZ, nullpoly(m))
    An = lambda m: ('AV', NV, nullpoly(m))
    op = [(1, [Zk(k), Zk(N - k)])]
    for i in range(k):
        c = math.comb(k, i + 1)
        op.append((c, [An(i), Zk(k - 1 - i), Zk(N - k)]))
        op.append((-c, [Zk(k - 1 - i), An(i), Zk(N - k)]))
    for i in range(N - k):
        c = math.comb(N - k, i + 1)
        op.append((c, [Zk(k), An(i), Zk(N - k - 1 - i)]))
        op.append((-c, [Zk(k), Zk(N - k - 1 - i), An(i)]))
    return op


def twist_two(N, seed=3, covariant=True):
    """matrix of lc^natural (units u) on the two-field parts B_k = tr(d^k Z d^{N-k} Z) of the gauge-invariant
    operators tr(D^k Z D^{N-k} Z); the two-field part of lc gets the pairs of the two-letter term (r = 2 outputs) and
    the pairs (A_n, Z) of the three-letter terms (r = 1 output).  Eigenvalues vs -8 h(S), S even <= N."""
    rng = np.random.default_rng(seed + N)
    ks = list(range(N // 2 + 1))
    configs = []
    while len(configs) < len(ks) + 2:
        q1 = rng.integers(-2, 3, 4); q2 = rng.integers(-2, 3, 4)
        configs.append((q1, q2))
    def zq(z):
        return sp.Rational(z[0]) + sp.I * sp.Rational(z[1])
    def Bval(k, q1, q2):
        legs = [('AW', eZ, q1), ('AW', eZ, q2)]
        v = O.trace_on_legs([('AW', zZ, nullpoly(k)), ('AW', zZ, nullpoly(N - k))], legs)
        return sp.Integer(int(round(v.real))) + sp.I * int(round(v.imag))
    def lcnat_val(k, q1, q2):
        legs = [('AW', eZ, q1), ('AW', eZ, q2)]
        op = covariant_twist_op(k, N) if covariant else [(1, [('AW', zZ, nullpoly(k)), ('AW', zZ, nullpoly(N - k))])]
        v = O.lc_planar_on_legs(op, legs)
        return zq(v) - 4 * Bval(k, q1, q2)
    A = sp.Matrix([[Bval(k, *c) for k in ks] for c in configs])
    Y = sp.Matrix([[lcnat_val(k, *c) for k in ks] for c in configs])
    # exact least squares (A has full column rank): normal equations, then exact residual
    AH = A.H
    Mq = sp.simplify((AH * A).inv() * (AH * Y))
    resid = sp.simplify(A * Mq - Y)
    maxres = max([abs(complex(x)) for x in resid] + [0.0])
    Mq = Mq.applyfunc(sp.nsimplify)
    ev = Mq.eigenvals()
    expected = sorted([-8 * harmonic(S) for S in range(0, N + 1, 2)])
    got = [sp.nsimplify(sp.simplify(e)) for e, mlt in ev.items() for _ in range(mlt)]
    try:
        got = sorted(got, key=lambda z: float(z))
        ok = [sp.Rational(e) for e in expected] == got
    except TypeError:
        ok = False
    return {'N': N, 'matrix_units': [[str(Mq[i, j]) for j in range(len(ks))] for i in range(len(ks))],
            'eigenvalues_units': [str(g) for g in got], 'expected_-8h(S)_S_even<=N': [str(e) for e in expected],
            'fit_residual': maxres, 'pass': bool(ok)}


if __name__ == '__main__':
    M2 = two_site_su2()
    print('two-site', [[str(x) for x in r] for r in M2])
    for L in range(2, 7):
        for Mm in range(0, L // 2 + 1):
            Em, cp = cyclic_spectrum(M2, L, Mm)
            print(L, Mm, cp)
    for N in range(0, 7):
        t0 = time.time()
        r = twist_two(N)
        print(N, r['eigenvalues_units'], r['expected_-8h(S)_S_even<=N'], r['pass'], r['fit_residual'], round(time.time() - t0, 1))
