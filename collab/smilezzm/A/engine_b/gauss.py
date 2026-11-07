"""Engine B -- symbolic Gaussian integration over internal vertex positions.

A graph has two external vertices 'x' (position x, |x|^2 = X2) and 'o' (the origin) and internal vertices
0..k-1 (positions z_i in R^4, integrated).  Every edge e = (u, v, t_e) carries the heat kernel
    Kt(u - v) = (4 pi t_e)^(-2) exp(-|u - v|^2 / (4 t_e)),
the Schwinger-time integrand of D_eps(r) = int_eps^inf Kt(r) dt.  A derivative d/du_m of the kernel of edge e
gives the vector factor  -(u - v)_m / (2 t_e).  The numerator is a sum of terms coeff * prod of dot products of
such vector factors (here at most one dot product per term).

gaussian(...) returns (pref, num_x2, num_1, phi) with
    int prod_i d^4 z_i  prod_e Kt_e * numerator  =  pref * (X2 * num_x2 + num_1) * exp(-X2 * phi),
all sympy expressions in the Schwinger times.
"""
import sympy as sp

PI = sp.pi


def vec_deriv(u, v, t):
    """Vector factor of d/du of the kernel of an edge (u, v, t):  -(u - v)/(2 t).
    Returned as (coefficients over positions, tag); the tag (edge, end) lets gaussian() add the contact term
    of two derivatives acting on the same edge:  d_u d_u K = [w w/(4t^2) - 1/(2t)] K,
    d_u d_v K = -[w w/(4t^2) - 1/(2t)] K  (w = u - v), i.e. -dim/(2t) resp. +dim/(2t) in a dot product."""
    d = {}
    for w, c in ((u, -sp.Rational(1, 2) / t), (v, sp.Rational(1, 2) / t)):
        if w == 'o':
            continue
        d[w] = d.get(w, 0) + c
    return (d, (frozenset([u, v]), t, u))


def gaussian(k, edges, terms, dim=4):
    """k internal vertices; edges = [(u, v, t)], u, v in {'x','o',0..k-1};
    terms = [(coeff, [(vecA, vecB)])] (list of dot-product pairs; [] = constant 1)."""
    M = sp.zeros(k, k)
    B = sp.zeros(k, 1)
    C = sp.Integer(0)
    for (u, v, t) in edges:
        c = 1 / t
        ends = [u, v]
        ints = [w for w in ends if isinstance(w, int)]
        exts = [w for w in ends if not isinstance(w, int)]
        if len(ints) == 2:
            i, j = ints
            M[i, i] += c
            M[j, j] += c
            M[i, j] -= c
            M[j, i] -= c
        elif len(ints) == 1:
            i = ints[0]
            M[i, i] += c
            if exts[0] == 'x':
                B[i] += c
                C += c
        else:
            if set(exts) == {'x', 'o'}:
                C += c
    detM = sp.factor(sp.together(M.det()))
    Minv = M.inv() if k > 0 else M
    Minv = Minv.applyfunc(lambda e: sp.factor(sp.together(e)))
    mbar = (Minv * B).applyfunc(lambda e: sp.factor(sp.together(e)))
    phi = sp.factor(sp.together((C - (B.T * Minv * B)[0, 0] if k > 0 else C) / 4))
    pref = sp.Integer(1)
    for (_, _, t) in edges:
        pref *= (4 * PI * t) ** (-2)
    pref *= ((4 * PI) ** k / detM) ** (dim // 2) if k > 0 else 1

    def mean_coeff(vec):
        s = vec.get('x', 0)
        for i in range(k):
            s += vec.get(i, 0) * mbar[i]
        return s

    num_x2 = sp.Integer(0)
    num_1 = sp.Integer(0)
    for coeff, pairs in terms:
        if len(pairs) == 0:
            num_1 += coeff
            continue
        if len(pairs) != 1:
            raise NotImplementedError('only one dot product per term')
        a, b = pairs[0]
        ta = tb = None
        if isinstance(a, tuple):
            a, ta = a
        if isinstance(b, tuple):
            b, tb = b
        if ta is not None and tb is not None and ta[0] == tb[0] and ta[1] == tb[1]:
            num_1 += coeff * (-1 if ta[2] == tb[2] else 1) * sp.Integer(dim) / (2 * ta[1])
        num_x2 += coeff * mean_coeff(a) * mean_coeff(b)
        cov = 0
        for i in range(k):
            for j in range(k):
                cov += a.get(i, 0) * b.get(j, 0) * 2 * Minv[i, j]
        num_1 += coeff * dim * cov
    return (sp.factor(pref), sp.factor(sp.together(num_x2)), sp.factor(sp.together(num_1)), phi)
