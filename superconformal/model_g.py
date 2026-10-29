"""
model_g.py -- the functionals S, S^(1), S^(2) of model.py for an arbitrary gauge algebra g given by
  * a basis T_1..T_d with rational structure constants [T_j, T_k] = sum_l f_jk^l T_l,
  * a diagonal invariant pairing (T_i, T_j) = k_i delta_ij (k_i rational, nonzero).
Fields A, lam, c carry components along T_i; antifields A^*, lam^*, c^* carry DUAL components
(phi^*_i = (phi^*, T_i)), so that the BV bracket is the plain contraction of engine.bracket and
  (X, phi^*) = X^i phi^*_i,   (X, Y) = k_i X^i Y^i,   (phi^*, psi^*) = k_i^{-1} phi^*_i psi^*_i.
For g = su(2), k = (1,1,1) this reproduces model.py term by term (checked in exhaustive_tower.py).
Clifford data are those of model.py.
"""
import os, sys, itertools
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
import engine as E
import model as Mo

V4, WP, NS = Mo.V4, Mo.WP, Mo.NS
SIG, SIGB, SIG_NZ, SIGB_NZ = Mo.SIG, Mo.SIGB, Mo.SIG_NZ, Mo.SIGB_NZ
matvec = Mo.matvec
e_mu = Mo.e_mu


class Gauge:
    def __init__(self, name, dim, FL, k):
        self.name, self.dim, self.FL, self.k = name, dim, FL, k
        self.kinv = [1 / x for x in k]

    def lie_tests(self):
        """Jacobi, invariance of the pairing, nondegeneracy; exact."""
        f = {}
        for (j, k, l), v in self.FL:
            f[(j, k, l)] = f.get((j, k, l), Fr(0)) + v
        d = self.dim
        F = lambda j, k, l: f.get((j, k, l), Fr(0))
        anti = all(F(i, j, l) == -F(j, i, l) for i in range(d) for j in range(d) for l in range(d))
        jac = all(sum(F(j, k, m) * F(i, m, l) + F(k, i, m) * F(j, m, l) + F(i, j, m) * F(k, m, l) for m in range(d)) == 0
                  for i, j, k, l in itertools.product(range(d), repeat=4))
        # invariance: ([T_i,T_j],T_l) + (T_j,[T_i,T_l]) = 0  <=>  f_ij^l k_l + f_il^j k_j = 0
        inv = all(F(i, j, l) * self.k[l] + F(i, l, j) * self.k[j] == 0 for i in range(d) for j in range(d) for l in range(d))
        return {'dim': d, 'antisymmetric': anti, 'jacobi': jac, 'pairing_invariant': inv,
                'pairing_nondegenerate': all(x != 0 for x in self.k), 'pairing_diagonal': [str(x) for x in self.k]}


def su2():
    return Gauge('su2', Mo.DIMG, list(Mo.FL), [Fr(1)] * Mo.DIMG)


def su3():
    """basis T_1..T_7 orthonormal (compact Gell-Mann, [T_a,T_b] = f_abc T_c), T_8' = sqrt(3) T_8:
    f'_ab^c = s_a s_b f_abc / s_c with s = (1,..,1,sqrt 3) is rational; k = (1,..,1,3)."""
    import sympy.combinatorics  # noqa: F401  (used by n4_model.structure_constants)
    CL = os.path.join(os.path.dirname(HERE), 'classical')
    sys.path.insert(0, CL)
    import n4_model as nm
    sys.path.pop(0)
    d, f = nm.structure_constants('su3')
    s = [sp.Integer(1)] * 7 + [sp.sqrt(3)]
    FL = []
    for (a, b, c), v in f.items():
        w = sp.nsimplify(s[a] * s[b] * v / s[c])
        assert w.is_rational, (a, b, c, w)
        if w != 0:
            FL.append(((a, b, c), Fr(int(w.p), int(w.q))))
    return Gauge('su3', d, FL, [Fr(1)] * 7 + [Fr(3)])


class Fields(Mo.Fields):
    def __init__(self, V, G):
        self.V = V
        self.G = G
        self.dimg = G.dim

    def bracket_g(self, X, Y, i):
        out = {}
        for (j, k, l), val in self.G.FL:
            if l == i:
                E.padd_into(out, E.pmul(X(j), Y(k)), (val, Fr(0)))
        return out


def build_S(V, G):
    Fd = Fields(V, G)
    S = {}
    for M in range(10):
        for N in range(10):
            if M == N:
                continue
            for i in range(G.dim):
                f = Fd.F(M, N, i)
                E.padd_into(S, E.pmul(f, f), (Fr(-1, 8) * G.k[i], Fr(0)))
    for M in range(10):
        for (a, b, s) in SIG_NZ[M]:
            for i in range(G.dim):
                E.padd_into(S, E.pmul(Fd.lam(a, i), Fd.Dlam(M, b, i)), E.cmul(s, (Fr(1, 2) * G.k[i], Fr(0))))
    for M in range(10):
        for i in range(G.dim):
            Dc = {}
            if M in V4:
                E.padd_into(Dc, Fd.c(i, e_mu(M)))
            E.padd_into(Dc, Fd.bracket_g(lambda j: Fd.A(M, j), lambda k: Fd.c(k), i))
            E.padd_into(S, E.pmul(Dc, Fd.As(M, i)))
    for a in range(NS):
        for i in range(G.dim):
            cl = Fd.bracket_g(lambda j: Fd.c(j), lambda k: Fd.lam(a, k), i)
            E.padd_into(S, E.pmul(cl, Fd.lams(a, i)))
    for i in range(G.dim):
        cc = Fd.bracket_g(lambda j: Fd.c(j), lambda k: Fd.c(k), i)
        E.padd_into(S, E.pmul(cc, Fd.cs(i)), (Fr(1, 2), Fr(0)))
    return S


def killing_spinor(V, Q, eta):
    return Mo.killing_spinor(V, Q, eta)


def S1_odd(V, G, Q, eta, aS=Fr(-1)):
    """identical to model.S1_odd (antifield-linear: no pairing factor), with g = G."""
    Fd = Fields(V, G)
    eps = killing_spinor(V, Q, eta)
    out = {}
    for M in range(10):
        for (a, b, s) in SIG_NZ[M]:
            if not eps[b]:
                continue
            for i in range(G.dim):
                E.padd_into(out, E.pmul(E.pmul(eps[b], Fd.lam(a, i)), Fd.As(M, i)), E.cneg(s))
    for M in range(10):
        for N in range(M + 1, 10):
            vec = [dict() for _ in range(NS)]
            for (b, cc, s2) in SIG_NZ[N]:
                if not eps[cc]:
                    continue
                for a in range(NS):
                    s1 = SIGB[M][a][b]
                    if E.ciszero(s1):
                        continue
                    E.padd_into(vec[a], eps[cc], E.cmul(s1, s2))
            for i in range(G.dim):
                f = Fd.F(M, N, i)
                if not f:
                    continue
                for a in range(NS):
                    if vec[a]:
                        E.padd_into(out, E.pmul(E.pmul(vec[a], f), Fd.lams(a, i)), (Fr(1, 2), Fr(0)))
    if any(not E.ciszero(e) for e in eta):
        for w in WP:
            t = matvec(SIGB[w], eta)
            for a in range(NS):
                if E.ciszero(t[a]):
                    continue
                for i in range(G.dim):
                    E.padd_into(out, E.pmul(Fd.A(w, i), Fd.lams(a, i)), E.cmul(t[a], (aS, Fr(0))))
    return out


def S1_even(V, G, vf, R, kl=Fr(1, 4), wl=Fr(3, 8), wA=Fr(1, 4)):
    Fd = Fields(V, G)
    D = Mo.even_data(V, vf, R)
    out = {}
    for i in range(G.dim):
        for M in range(10):
            dA = {}
            for nu in V4:
                if vf[nu]:
                    E.padd_into(dA, E.pmul(vf[nu], Fd.A(M, i, e_mu(nu))))
            if M in V4:
                for nu in V4:
                    if D['dv'][M][nu]:
                        E.padd_into(dA, E.pmul(Fd.A(nu, i), D['dv'][M][nu]))
            else:
                for (a_, b_), cval in R.items():
                    if b_ == M:
                        E.padd_into(dA, Fd.A(a_, i), cval)
                if D['div']:
                    E.padd_into(dA, E.pmul(D['div'], Fd.A(M, i)), (wA, Fr(0)))
            if dA:
                E.padd_into(out, E.pmul(dA, Fd.As(M, i)))
        for a in range(NS):
            dl = {}
            for nu in V4:
                if vf[nu]:
                    E.padd_into(dl, E.pmul(vf[nu], Fd.lam(a, i, e_mu(nu))))
            for (M, N), om in D['Om'].items():
                for (b, cc, s2) in SIG_NZ[N]:
                    s1 = SIGB[M][a][b]
                    if E.ciszero(s1):
                        continue
                    E.padd_into(dl, E.pmul(om, Fd.lam(cc, i)), E.cmul((kl, Fr(0)), E.cmul(s1, s2)))
            if D['div']:
                E.padd_into(dl, E.pmul(D['div'], Fd.lam(a, i)), (wl, Fr(0)))
            if dl:
                E.padd_into(out, E.pmul(dl, Fd.lams(a, i)), (Fr(-1), Fr(0)))
        dc = {}
        for nu in V4:
            if vf[nu]:
                E.padd_into(dc, E.pmul(vf[nu], Fd.c(i, e_mu(nu))))
        if dc:
            E.padd_into(out, E.pmul(dc, Fd.cs(i)), (Fr(-1), Fr(0)))
    return out


def S2_odd(V, G, e1, e2):
    """[ESW (34)] pointwise: 1/4 (Gamma(e1,e2), Gamma(lam^*,lam^*)) - 1/2 (e1,lam^*)(e2,lam^*) - iota_Gamma A c^*;
    the two lam^* lam^* terms use the inverse pairing k_i^{-1} on dual components."""
    Fd = Fields(V, G)
    Gm = {M: Mo.Gamma_poly(e1, e2, M) for M in range(10)}
    out = {}
    for M in range(10):
        if not Gm[M]:
            continue
        for i in range(G.dim):
            gl = {}
            for (a, b, s) in SIGB_NZ[M]:
                E.padd_into(gl, E.pmul(Fd.lams(a, i), Fd.lams(b, i)), s)
            if gl:
                E.padd_into(out, E.pmul(Gm[M], gl), (Fr(1, 4) * G.kinv[i], Fr(0)))
            E.padd_into(out, E.pmul(E.pmul(Gm[M], Fd.A(M, i)), Fd.cs(i)), (Fr(-1), Fr(0)))
    for i in range(G.dim):
        p1, p2 = {}, {}
        for a in range(NS):
            if e1[a]:
                E.padd_into(p1, E.pmul(e1[a], Fd.lams(a, i)))
            if e2[a]:
                E.padd_into(p2, E.pmul(e2[a], Fd.lams(a, i)))
        E.padd_into(out, E.pmul(p1, p2), (Fr(-1, 2) * G.kinv[i], Fr(0)))
    return out
