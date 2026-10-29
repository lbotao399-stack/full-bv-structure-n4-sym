"""Twist two (sl(2) sector) from the referee's OWN kernels: O_k = tr(D^k Z D^{N-k} Z), D = n.d + [n.A, -], n null.
lc^natural restricted to the ZZ part (two-letter ZZ kernel + r = 1 outputs of the pairs (A_n, Z), (Z, A_n) in the
O(A_n) terms of the covariant derivatives) -> matrix on {O_k}; Beisert: eig = -8 h(S) u, S even <= N."""
import sys, json, itertools, time
import sympy as sp
from fractions import Fraction
from ref_core import *
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 6
n = [1, sp.I, 0, 0]
nb = [1, -sp.I, 0, 0]                                # n.nb = 2
rho1, rho2 = sp.symbols('rho1 rho2')
r1 = [rho1 * x / 2 for x in nb]; r2 = [rho2 * x / 2 for x in nb]   # n.r_j = rho_j
z = [0] * 8 + [1, sp.I]
zbar = [0] * 8 + [sp.Rational(1, 2), -sp.I / 2]     # selector: value of the output letter Z
nsel = [1, sp.I] + [0] * 8                          # letter n.A
def dn(k):
    """letter factor (i n.p)^k"""
    return lambda p: (sp.I * sum(n[m] * p[m] for m in range(4))) ** k
def coeffs2(expr):
    P = sp.Poly(sp.expand(expr), rho1, rho2)
    return {(mon[0], mon[1]): c for mon, c in P.terms()}
def coeffs1(expr):
    P = sp.Poly(sp.expand(expr), rho1)
    return {mon[0]: c for mon, c in P.terms()}
# output (n.r)^c <-> (-i n.d)^c on the output: Z_c coefficient picks (-i)^c
t0 = time.time()
KZZ = {}
for a in range(NMAX + 1):
    for b in range(NMAX + 1 - a):
        tot = 0
        for verts, rs in [([('AAAA', zbar, zbar)], [[r1[m] + r2[m] for m in range(4)]]),
                          ([('AAA', zbar), ('AAA', zbar)], [r1, r2])]:
            res = kernel_chain('A', z, 'A', z, verts, rs, fa=dn(a), fb=dn(b))
            tot += res.get('log', 0)
        KZZ[(a, b)] = {cd: sp.nsimplify(c * (-sp.I) ** (cd[0] + cd[1])) for cd, c in coeffs2(tot).items()}
print('KZZ done', time.time() - t0)
KAZ, KZA = {}, {}
for i in range(NMAX):
    for m in range(NMAX - i):
        res = kernel_chain('A', nsel, 'A', z, [('AAA', zbar)], [r1], fa=dn(i), fb=dn(m)).get('log', 0)
        KAZ[(i, m)] = {c: sp.nsimplify(v * (-sp.I) ** c) for c, v in coeffs1(res).items()}
        res = kernel_chain('A', z, 'A', nsel, [('AAA', zbar)], [r1], fa=dn(m), fb=dn(i)).get('log', 0)
        KZA[(m, i)] = {c: sp.nsimplify(v * (-sp.I) ** c) for c, v in coeffs1(res).items()}
print('KAZ/KZA done', time.time() - t0)
# planar engine on cyclic words of letters ('Z', k) / ('A', i); kernel value c (units hbar/16pi^2) -> (-1)^r c/2 u
def kern(X, Y):
    if X[0] == 'Z' and Y[0] == 'Z':
        return [([('Z', c), ('Z', d)], v) for (c, d), v in KZZ.get((X[1], Y[1]), {}).items()]
    if X[0] == 'A' and Y[0] == 'Z':
        return [([('Z', c)], v) for c, v in KAZ.get((X[1], Y[1]), {}).items()]
    if X[0] == 'Z' and Y[0] == 'A':
        return [([('Z', c)], v) for c, v in KZA.get((X[1], Y[1]), {}).items()]
    return []
def canon(word):
    """cyclic word of bosonic letters -> canonical rotation"""
    rots = [tuple(word[i:] + word[:i]) for i in range(len(word))]
    return min(rots)
def lc_planar(op):
    """op: {word(tuple): coeff}; returns planar lc on ZZ-only outputs, units u"""
    out = {}
    for word, cw in op.items():
        Lw = len(word)
        for p in range(Lw):
            q = (p + 1) % Lw
            if Lw == 2 and p == 1:
                # L = 2: second adjacency = wrap term (S = full) of the pair (word[0], word[1])
                for outs, v in kern(word[0], word[1]):
                    nw = canon(list(reversed(outs)))
                    if all(x[0] == 'Z' for x in nw):
                        out[nw] = out.get(nw, 0) + cw * v / 2
                continue
            if q == 0:
                # wrap pair (X_1, X_L): S = full term, outputs reversed, coefficient +c/2
                a, b = 0, Lw - 1
                for outs, v in kern(word[a], word[b]):
                    nw = canon(list(word[1:b]) + list(reversed(outs)))
                    if all(x[0] == 'Z' for x in nw):
                        out[nw] = out.get(nw, 0) + cw * v / 2
                continue
            for outs, v in kern(word[p], word[q]):
                r = len(outs)
                nw = canon(list(word[:p]) + outs + list(word[q + 1:]))
                if all(x[0] == 'Z' for x in nw):
                    out[nw] = out.get(nw, 0) + cw * (-1) ** r * v / 2
    return out
def binom(a, b):
    return sp.binomial(a, b)
def Dm_terms(mp):
    """D^{m'} Z = Z_{m'} + sum_{0<=i<=j<m'} C(j,i) [A_i, Z_{m'-1-i}] + O(A^2): list of (coef, [letters]) for the
    commutator expanded: [A_i, Z_m] = A_i Z_m - Z_m A_i"""
    terms = [(1, [('Z', mp)])]
    for j in range(mp):
        for i in range(j + 1):
            c = binom(j, i); m = mp - 1 - i
            terms.append((c, [('A', i), ('Z', m)]))
            terms.append((-c, [('Z', m), ('A', i)]))
    return terms
results = {}
for N in range(0, NMAX + 1):
    basis = sorted(set(canon([('Z', k), ('Z', N - k)]) for k in range(N + 1)))
    ks = list(range(N // 2 + 1))
    Mrows = []
    for k in ks:
        op = {}
        for c1, w1 in Dm_terms(k):
            for c2, w2 in Dm_terms(N - k):
                w = w1 + w2
                if sum(1 for x in w if x[0] == 'A') > 1:
                    continue
                cw = canon(w)
                op[cw] = op.get(cw, 0) + c1 * c2
        lc = lc_planar(op)
        # lc^natural = lc - u sum n_F on the ZZ part
        zz = canon([('Z', k), ('Z', N - k)])
        lc[zz] = lc.get(zz, 0) - 4
        Mrows.append([sp.nsimplify(lc.get(bv, 0)) for bv in basis])
    # O_k has ZZ part tr(Z_k Z_{N-k}) (+ the same for N-k): basis vector index of O_k
    Bmat = sp.zeros(len(ks), len(basis))
    for r_, k in enumerate(ks):
        Bmat[r_, basis.index(canon([('Z', k), ('Z', N - k)]))] = 1
    Mt = sp.Matrix(Mrows) * Bmat.T.inv() if Bmat.shape[0] == Bmat.shape[1] else None
    ev = Mt.eigenvals()
    harm = lambda S: sum(sp.Rational(1, j) for j in range(1, S + 1))
    pred = sorted([-8 * harm(S) for S in range(0, N + 1, 2)])
    got = sorted(sum([[k] * v for k, v in ev.items()], []))
    results[N] = dict(eig=[str(x) for x in got], beisert=[str(x) for x in pred], match=(got == pred))
    print('N =', N, 'eig lc^natural (units u):', got, ' -8h(S):', pred, ' match', got == pred)
json.dump(dict(results=results, KZZ={str(k): {str(a): str(b) for a, b in v.items()} for k, v in KZZ.items()},
               KAZ={str(k): {str(a): str(b) for a, b in v.items()} for k, v in KAZ.items()},
               KZA={str(k): {str(a): str(b) for a, b in v.items()} for k, v in KZA.items()}),
          open('out_twist2.json', 'w'), indent=1)
print('time', time.time() - t0)
