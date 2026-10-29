"""Finite-dimensional graded model for Part D (dilatation, arity one).

Checks, in an odd-symplectic graded vector space E with Q^2 = 0, a gauge-fixing operator Q^GF
(Feynman type with (Q^GF)^2 != 0 and Landau type with (Q^GF)^2 = 0), H = [Q,Q^GF], heat kernels
K_t = (e^{-tH} x 1)K_0, propagator P(eps,T) and an even operator M preserving the pairing with
[M,Q] = 0, every algebraic identity used in step2/theory/delta_theory.tex.

Model A (exact, Fractions): fibre (c, A_L, A_T, A*_L, A*_T, c*) tensor R[z]/(z^{N+1}),
  Q = Qf x 1, Q^GF = Qgf x z, H = 1 x z (nilpotent), M = 1 x (-4 z d/dz + 2N):
  [M,Q] = 0, [M,Q^GF] = -4 Q^GF, [M,H] = -4H exactly, so the scaling identities hold at finite T.
Model B (floats, 1e-12): fibre with momentum parameter s (H = s^2 on every slot), two copies,
  M = a generic even pairing-preserving operator commuting with Q (found as a null space);
  H invertible, so the homotopy Psi = [Q, Phi] with Phi = k^{-1} P* Psi is available.
Conventions (fixed here and used in the text):
  functionals = polynomials in coordinates x^i, left derivatives d_i;
  D_X = sum X^k_j x^j d_k (derivation induced by X), Q on functionals := D_Q;
  d_phi = (1/2) sum phi^{ij} d_i d_j; kernel of an operator A: (A x 1)K_0, components (A K_0)^{ij};
  Delta_t = -d_{K_t}; {f,g}_t = Delta_t(fg) - Delta_t(f) g - (-1)^{|f|} f Delta_t(g);
  S_X = the quadratic functional with {S_X, -} = D_X.
Run: python3 delta_toy.py  (writes delta_toy_results.json next to this file).
"""
import itertools, json, os, random
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------- super polynomials
class Alg:
    def __init__(self, par, zero=Fr(0)):
        self.par = list(par); self.n = len(par); self.zero = zero

    def clean(self, f, tol=0.0):
        return {m: c for m, c in f.items() if (abs(c) > tol if tol else c != 0)}

    def add(self, *fs, coeffs=None):
        out = {}
        for k, f in enumerate(fs):
            a = 1 if coeffs is None else coeffs[k]
            for m, c in f.items():
                out[m] = out.get(m, self.zero) + a * c
        return self.clean(out)

    def scale(self, a, f):
        return self.clean({m: a * c for m, c in f.items()})

    def sortmono(self, m):
        m = list(m); sign = 1
        for i in range(1, len(m)):
            j = i
            while j > 0 and m[j - 1] > m[j]:
                if self.par[m[j - 1]] and self.par[m[j]]:
                    sign = -sign
                m[j - 1], m[j] = m[j], m[j - 1]; j -= 1
        for i in range(1, len(m)):
            if m[i] == m[i - 1] and self.par[m[i]]:
                return None, 0
        return tuple(m), sign

    def mul(self, f, g):
        out = {}
        for m1, c1 in f.items():
            for m2, c2 in g.items():
                m, s = self.sortmono(m1 + m2)
                if m is None: continue
                out[m] = out.get(m, self.zero) + s * c1 * c2
        return self.clean(out)

    def parity(self, f):
        ps = {sum(self.par[i] for i in m) % 2 for m in f}
        assert len(ps) <= 1, "inhomogeneous parity"
        return ps.pop() if ps else 0

    def x(self, i):
        return {(i,): Fr(1) if isinstance(self.zero, Fr) else 1.0}

    def d(self, i, f):
        out = {}
        for m, c in f.items():
            if i not in m: continue
            k = m.index(i); mult = m.count(i)
            sign = (-1) ** (self.par[i] * sum(self.par[j] for j in m[:k]))
            mm = m[:k] + m[k + 1:]
            out[mm] = out.get(mm, self.zero) + sign * mult * c
        return self.clean(out)

    def D(self, X, f):  # derivation induced by the matrix X (X[k][j]: e_j -> sum_k X[k][j] e_k)
        out = {}
        for k in range(self.n):
            dk = self.d(k, f)
            if not dk: continue
            for j in range(self.n):
                if X[k][j] != 0:
                    out = self.add(out, self.mul(self.x(j), dk), coeffs=[1, X[k][j]])
        return out

    def dphi(self, phi, f):
        out = {}
        for i in range(self.n):
            for j in range(self.n):
                if phi[i][j] != 0:
                    out = self.add(out, self.d(i, self.d(j, f)), coeffs=[1, phi[i][j] / 2])
        return out

    def lap(self, K, f):
        return self.scale(-1, self.dphi(K, f))

    def br(self, K, f, g):
        pf = self.parity(f)
        return self.add(self.lap(K, self.mul(f, g)), self.mul(self.lap(K, f), g),
                        self.mul(f, self.lap(K, g)), coeffs=[1, -1, -((-1) ** pf)])

    def C(self, phi, f, g):  # d_phi(fg) - d_phi(f) g - f d_phi(g)
        return self.add(self.dphi(phi, self.mul(f, g)), self.mul(self.dphi(phi, f), g),
                        self.mul(f, self.dphi(phi, g)), coeffs=[1, -1, -1])

    def exp_dphi(self, phi, hbar, f, nmax=12):
        out, term = dict(f), dict(f)
        for k in range(1, nmax):
            term = self.scale(hbar / k, self.dphi(phi, term))
            if not term: break
            out = self.add(out, term)
        return out

    def dist(self, f, g):
        h = self.add(f, g, coeffs=[1, -1])
        return max([abs(float(c)) for c in h.values()] + [0.0])

    def rand(self, deg_list, nterms, rng, par_target, coeff):
        f = {}
        while len(f) < nterms:
            dgr = rng.choice(deg_list)
            m = tuple(rng.randrange(self.n) for _ in range(dgr))
            if sum(self.par[i] for i in m) % 2 != par_target: continue
            m2, s = self.sortmono(m)
            if m2 is None: continue
            f[m2] = coeff(rng)
        return self.clean(f)

# ---------------------------------------------------------------- matrices (lists, exact or float)
def mzero(n, z=Fr(0)): return [[z] * n for _ in range(n)]
def mid(n, one=Fr(1)):
    I = mzero(n, one * 0)
    for i in range(n): I[i][i] = one
    return I
def mmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]
def madd(A, B, a=1, b=1): return [[a * A[i][j] + b * B[i][j] for j in range(len(A[0]))] for i in range(len(A))]
def mscale(a, A): return [[a * x for x in r] for r in A]
def mT(A): return [list(r) for r in zip(*A)]
def kron(A, B):
    n, m = len(A), len(B)
    return [[A[i // m][j // m] * B[i % m][j % m] for j in range(n * m)] for i in range(n * m)]
def mmax(A): return max(abs(float(x)) for r in A for x in r)
def comm(A, B, graded_sign=1):  # AB - graded_sign*BA
    return madd(mmul(A, B), mmul(B, A), 1, -graded_sign)
def expm_nil(Z, t, N):  # e^{-tZ} for nilpotent Z with Z^{N+1} = 0, exact
    out, term = mid(len(Z), Fr(1)), mid(len(Z), Fr(1))
    for k in range(1, N + 1):
        term = mscale(Fr(-t) / k, mmul(term, Z)); out = madd(out, term)
    return out
def kernel_of(A, K0): return mmul(A, K0)              # (A x 1)K_0
def op_of(phi, K0inv): return mmul(phi, K0inv)
def one_tensor(A, phi, par, parA):                     # (1 x A)phi, component (i,j) = sum_l phi[i][l] (-1)^{|A| p_i} A[j][l]
    n = len(phi)
    return [[sum(phi[i][l] * A[j][l] for l in range(n)) * ((-1) ** (parA * par[i])) for j in range(n)] for i in range(n)]
def minv_exact(A):
    n = len(A); M = [list(r) + [Fr(int(i == j)) for j in range(n)] for i, r in enumerate(A)]
    for c in range(n):
        p = next(r for r in range(c, n) if M[r][c] != 0); M[c], M[p] = M[p], M[c]
        pv = M[c][c]; M[c] = [x / pv for x in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c]; M[r] = [a - f * b for a, b in zip(M[r], M[c])]
    return [r[n:] for r in M]

# ---------------------------------------------------------------- fibre model
# basis: 0 c (odd), 1 A_L (even), 2 A_T (even), 3 A*_L (odd), 4 A*_T (odd), 5 c* (even)
FPAR = [1, 0, 0, 1, 1, 0]
def fibre(s, a1, a3, kc, kL, kT, b1, one=Fr(1)):
    z = one * 0
    Q = mzero(6, z); G = mzero(6, z); K0 = mzero(6, z)
    a2, b2 = one * Fr(-1, 2) if isinstance(one, Fr) else -0.5, -2 * one
    Q[1][0] = a1 * s; Q[4][2] = a2 * s * s; Q[5][3] = a3 * s
    G[1][3] = b1 * one; G[2][4] = b2; G[0][1] = s / a1; G[3][5] = s / a3
    for (i, j, k) in [(0, 5, kc), (1, 3, kL), (2, 4, kT)]:
        K0[i][j] = K0[j][i] = k * one
    return Q, G, K0

def adjoint_defects(Q, G, K0, par):
    dQ = madd(kernel_of(Q, K0), one_tensor(Q, K0, par, 1))          # (Q x 1 + 1 x Q)K_0
    dG = madd(kernel_of(G, K0), one_tensor(G, K0, par, 1), 1, -1)   # (G x 1 - 1 x G)K_0
    return mmax(dQ), mmax(dG)

def find_fibre_signs():
    sols = []
    for a1, a3, kc, kL, kT in itertools.product([1, -1], repeat=5):
        for b1 in [-2, 0]:
            Q, G, K0 = fibre(Fr(1), Fr(a1), Fr(a3), kc, kL, kT, Fr(b1))
            if adjoint_defects(Q, G, K0, FPAR) == (0.0, 0.0):
                sols.append((a1, a3, kc, kL, kT, b1))
    return sols

FEYN = (1, 1, 1, -1, 1, -2)    # (a1, a3, kc, kL, kT, b1): b1 = -2 is the Feynman arrow A*_L -> A_L
LANDAU = (1, 1, 1, -1, 1, 0)   # same data without that arrow: (Q^GF)^2 = 0

def solve_exact(rows, rhs):
    """least-squares-free exact solve of an overdetermined consistent system; returns solution or None"""
    n = len(rows[0]); M = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv, r = [], 0
    for c in range(n):
        p = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if p is None: continue
        M[r], M[p] = M[p], M[r]; pv = M[r][c]; M[r] = [x / pv for x in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f * b for a, b in zip(M[i], M[r])]
        piv.append(c); r += 1
    if any(all(x == 0 for x in row[:-1]) and row[-1] != 0 for row in M): return None
    sol = [Fr(0)] * n
    for i, c in enumerate(piv): sol[c] = M[i][-1]
    return sol

def quad_for_derivation(alg, K0, X, parX):
    """S_X quadratic with {S_X, x^k} = D_X x^k for all k (exact). Returns polynomial or None."""
    want = (parX + 1) % 2
    monos = [(i, j) for i in range(alg.n) for j in range(i, alg.n)
             if (alg.par[i] + alg.par[j]) % 2 == want and not (i == j and alg.par[i])]
    cols = [alg.br(K0, {m: Fr(1)}, alg.x(k)) for m in monos for k in [None]] if False else None
    rows, rhs = [], []
    brs = {m: [alg.br(K0, {m: Fr(1)}, alg.x(k)) for k in range(alg.n)] for m in monos}
    for k in range(alg.n):
        target = alg.D(X, alg.x(k))
        keys = set(target) | set().union(*[set(brs[m][k]) for m in monos])
        for key in keys:
            rows.append([brs[m][k].get(key, Fr(0)) for m in monos]); rhs.append(target.get(key, Fr(0)))
    sol = solve_exact(rows, rhs)
    if sol is None: return None
    return alg.clean({m: c for m, c in zip(monos, sol)})

def deriv_matrix(alg, Cmap):
    """matrix Y with Cmap(x^k) = sum_j Y[k][j] x^j (Cmap linear on linear functionals)"""
    Y = mzero(alg.n)
    for k in range(alg.n):
        for m, c in Cmap(alg.x(k)).items():
            assert len(m) == 1; Y[k][m[0]] = c
    return Y

def kernel_from_second_order(alg, Op):
    """phi with Op = d_phi on quadratics: Op(x^a x^b) = (1/2)(phi^{ba} + (-1)^{p_a p_b} phi^{ab})"""
    phi = mzero(alg.n)
    for a in range(alg.n):
        for b in range(alg.n):
            v = Op(alg.mul(alg.x(a), alg.x(b)))
            c = v.get((), Fr(0))
            if a == b and alg.par[a]: continue
            phi[a][b] = c * (-1) ** (alg.par[a] * alg.par[b]) if a != b else c
    # symmetrise consistently with d_phi(x^a x^b) for a != b: value = (1/2)(phi^{ba} + (-1)^{pp} phi^{ab}) = (-1)^{pp} phi^{ab}
    return phi

# ---------------------------------------------------------------- model A (exact scaling model)
def model_A(signs, N=3, eps=Fr(1, 3), T=Fr(5, 2)):
    a1, a3, kc, kL, kT, b1 = signs
    Qf, Gf, Kf = fibre(Fr(1), Fr(a1), Fr(a3), kc, kL, kT, Fr(b1))
    n = N + 1
    Z = mzero(n); Kz = mzero(n); Mz = mzero(n); I = mid(n)
    for a in range(N): Z[a + 1][a] = Fr(1)                  # z^a -> z^{a+1}
    for a in range(n): Kz[a][N - a] = Fr(1); Mz[a][a] = Fr(-4 * a + 2 * N)
    par = [FPAR[i // n] for i in range(6 * n)]
    m = dict(par=par, N=N, eps=eps, T=T, n=6 * n,
             Q=kron(Qf, I), G=kron(Gf, Z), H=kron(mid(6), Z), M=kron(mid(6), Mz), K0=kron(Kf, Kz))
    m['K0inv'] = minv_exact(m['K0'])
    Ek = lambda t: kron(mid(6), expm_nil(Z, t, N))
    m['Kt'] = lambda t: kernel_of(Ek(t), m['K0'])
    # P(eps,T) = int_eps^T G e^{-tH} dt = Gf x (z int e^{-tz}) = Gf x (e^{-eps z} - e^{-T z}) (since Z G-factor)
    gz = madd(expm_nil(Z, eps, N), expm_nil(Z, T, N), 1, -1)    # e^{-eps z} - e^{-T z}
    m['k'] = kron(mid(6), gz)                                    # Q P* + P* Q
    # P* = int G e^{-tH} dt = Gf x (Z int_eps^T e^{-tZ} dt) = Gf x (e^{-eps z} - e^{-T z})
    m['Pstar'] = kron(Gf, gz)
    # dictionary: the kernel of an odd operator X is (Pi X x 1)K_0 (Pi = parity operator); this is the
    # convention in which [d_P, Q] = d_{K_eps} - d_{K_T} for Q = D_Q (checked below)
    m['Pi'] = [[(Fr(1) if par[i] == 0 else Fr(-1)) if i == j else Fr(0) for j in range(6 * n)] for i in range(6 * n)]
    m['P'] = kernel_of(mmul(m['Pi'], m['Pstar']), m['K0'])
    # exact derivative kernel: -4(eps d/deps + T d/dT)P = Gf x (4 eps z e^{-eps z} - 4 T z e^{-T z})
    dg = madd(mmul(Z, expm_nil(Z, eps, N)), mmul(Z, expm_nil(Z, T, N)), 4 * eps, -4 * T)
    m['scaleP'] = kernel_of(mmul(m['Pi'], kron(Gf, dg)), m['K0'])
    m['Ek'] = Ek
    m['alg'] = Alg(par)
    return m

# ---------------------------------------------------------------- checks in model A
def run_A(signs, seed=1, label=''):
    m = model_A(signs); alg = m['alg']; par = m['par']; K0 = m['K0']
    Q, G, H, M, P, Ps = m['Q'], m['G'], m['H'], m['M'], m['P'], m['Pstar']
    eps, T = m['eps'], m['T']; hb = Fr(2, 7)
    rng = random.Random(seed); cf = lambda r: Fr(r.randint(-3, 3) or 1, r.randint(1, 3))
    tests = [alg.rand([1, 2, 3], 5, rng, p, cf) for p in (0, 1, 0, 1)]
    res = {}
    Keps, KT = m['Kt'](eps), m['Kt'](T)
    Qf = lambda f: alg.D(Q, f)
    # S_M and S_Q (free action) as quadratics generating D_M and D_Q through {,} with K_0
    SM = quad_for_derivation(alg, K0, M, 0); SQ = quad_for_derivation(alg, K0, Q, 1)
    res['S_M exists'] = SM is not None; res['S_free exists'] = SQ is not None
    res['{S_M,F} - D_M F'] = max(alg.dist(alg.br(K0, SM, F), alg.D(M, F)) for F in tests)
    res['{S_free,F} - QF'] = max(alg.dist(alg.br(K0, SQ, F), Qf(F)) for F in tests)
    res['Q S_M'] = alg.dist(Qf(SM), {})
    # [d_phi, D_X] = -d_{(Xx1+1xX)phi}  (X = M even, phi = P);  for odd X = Q: [d_phi,Q] = d_{(Qx1+1xQ)phi} (graded)
    MP = madd(kernel_of(M, P), one_tensor(M, P, par, 0))
    res['[D_M,d_P] + d_{(Mx1+1xM)P}'] = max(alg.dist(alg.add(alg.D(M, alg.dphi(P, F)), alg.dphi(P, alg.D(M, F)),
                                                   alg.dphi(MP, F), coeffs=[1, -1, 1]), {}) for F in tests)
    Kd = madd(Keps, KT, 1, -1)
    res['d_P Q - Q d_P - (d_{K_eps} - d_{K_T})'] = max(alg.dist(alg.add(alg.dphi(P, Qf(F)), Qf(alg.dphi(P, F)),
                                              alg.dphi(Kd, F), coeffs=[1, -1, -1]), {}) for F in tests)
    # (Mx1+1xM)P(eps,T) = -4(eps d/deps + T d/dT)P(eps,T)
    res['(Mx1+1xM)P + 4(eps d_eps + T d_T)P'] = mmax(madd(MP, m['scaleP'], 1, -1))
    # C_P(S_M,-) is the derivation D_Y; identify Y
    Y = deriv_matrix(alg, lambda f: alg.C(P, SM, f))
    PM = mmul(Ps, M)
    res['C_P(S_M,-) = D_Y: Y - P*M'] = mmax(madd(Y, PM, 1, -1))
    res['C_P(S_M,F) - D_{P*M}F'] = max(alg.dist(alg.C(P, SM, F), alg.D(PM, F)) for F in tests)
    # second-order operator [d_P, C_P(S_M,-)] = d_Psi, Psi = kernel of P* M P*
    PMP = mmul(mmul(Ps, M), Ps); Psi = kernel_of(PMP, K0)   # even operator: no parity twist
    res['[d_P,C_P(S_M,-)] - d_Psi'] = max(alg.dist(alg.add(alg.dphi(P, alg.C(P, SM, F)), alg.C(P, SM, alg.dphi(P, F)),
                                                coeffs=[1, -1]), alg.dphi(Psi, F)) for F in tests)
    res['[d_P,d_Psi]'] = max(alg.dist(alg.dphi(P, alg.dphi(Psi, F)), alg.dphi(Psi, alg.dphi(P, F))) for F in tests)
    # conjugation formula e^{h d_P}(S_M F) = (S_M + h d_P S_M + h C_P(S_M,-) + (h^2/2) d_Psi) e^{h d_P} F
    def rhs(F):
        eF = alg.exp_dphi(P, hb, F)
        return alg.add(alg.mul(SM, eF), alg.scale(hb, alg.mul(alg.dphi(P, SM), eF)),
                       alg.scale(hb, alg.C(P, SM, eF)), alg.scale(hb * hb / 2, alg.dphi(Psi, eF)))
    res['conjugation formula'] = max(alg.dist(alg.exp_dphi(P, hb, alg.mul(SM, F)), rhs(F)) for F in tests)
    # Psi = (1/2){M,(P*)^2}, (P*)^2 = (Q^GF)^2 f(H)^2, [Q, P*MP*] = 0
    P2 = mmul(Ps, Ps)
    res['P*MP* - (1/2){M,P*^2}'] = mmax(madd(PMP, madd(mmul(M, P2), mmul(P2, M)), 1, Fr(-1, 2)))
    res['[Q,P*MP*]'] = mmax(comm(Q, PMP)); res['|(Q^GF)^2|'] = mmax(mmul(G, G)); res['|Psi|'] = mmax(Psi)
    res['d_Psi Q - Q d_Psi'] = max(alg.dist(alg.dphi(Psi, Qf(F)), Qf(alg.dphi(Psi, F))) for F in tests)
    # bracket at scale t: {S_M,F}_t = D_{e^{-tH}M}F
    for nm, t in (('eps', eps), ('T', T)):
        EM = mmul(m['Ek'](t), M)
        res['{S_M,F}_%s - D_{e^{-tH}M}F' % nm] = max(alg.dist(alg.br(m['Kt'](t), SM, F), alg.D(EM, F)) for F in tests)
    # homotopy identity: Q C_P(S_M,F) + C_P(S_M,QF) = {S_M,F}_eps - {S_M,F}_T   (Q = D_Q, S_M odd, d_P S_M = 0 here)
    res['d_P S_M'] = float(alg.dphi(P, SM).get((), 0)); res['d_{K_eps-K_T} S_M'] = float(alg.dphi(Kd, SM).get((), 0))
    res['Q C_P(S_M,F) + C_P(S_M,QF) - ({S_M,F}_eps - {S_M,F}_T)'] = max(alg.dist(
        alg.add(Qf(alg.C(P, SM, F)), alg.C(P, SM, Qf(F))),
        alg.add(alg.br(Keps, SM, F), alg.br(KT, SM, F), coeffs=[1, -1])) for F in tests)
    # transfer identity (Q + h Delta_T) e^{h d_P} = e^{h d_P} (Q + h Delta_eps)
    res['transfer identity'] = max(alg.dist(alg.add(Qf(alg.exp_dphi(P, hb, F)), alg.scale(hb, alg.lap(KT, alg.exp_dphi(P, hb, F)))),
                                            alg.exp_dphi(P, hb, alg.add(Qf(F), alg.scale(hb, alg.lap(Keps, F))))) for F in tests)
    # scaling Ward identity at operator level: D_M e^{h d_P} F = e^{h d_P} D_M F + 4 h (eps d_eps + T d_T) applied as d_{-scaleP}
    res['Ward: D_M e^{hd_P} - e^{hd_P} D_M - h d_{-(Mx1+1xM)P} e^{hd_P}'] = max(alg.dist(
        alg.D(M, alg.exp_dphi(P, hb, F)),
        alg.add(alg.exp_dphi(P, hb, alg.D(M, F)), alg.scale(-hb, alg.dphi(MP, alg.exp_dphi(P, hb, F))))) for F in tests)
    res['label'] = label
    return {k: (float(v) if isinstance(v, Fr) else v) for k, v in res.items()}


# ---------------------------------------------------------------- model B (floats): homotopy and large distance
def model_B(s, signs, eps=0.37, T=np.inf, copies=2, seed=0):
    a1, a3, kc, kL, kT, b1 = signs
    Qf, Gf, Kf = fibre(float(s), float(a1), float(a3), kc, kL, kT, float(b1), one=1.0)
    Qf, Gf, Kf = (np.array(X, dtype=float) for X in (Qf, Gf, Kf))
    I = np.eye(copies); Q = np.kron(I, Qf); G = np.kron(I, Gf); K0 = np.kron(I, Kf)
    par = FPAR * copies; n = len(par); Pi = np.diag([(-1.0) ** p for p in par])
    H = Q @ G + G @ Q
    h = s * s
    fT = (np.exp(-eps * h) - (0.0 if np.isinf(T) else np.exp(-T * h))) / h
    kT_ = h * fT
    Ps = G * fT; k = np.eye(n) * kT_
    # even pairing-preserving operators commuting with Q: null space of the linear constraints
    basis = []
    for i in range(n):
        for j in range(n):
            if par[i] == par[j]:
                E = np.zeros((n, n)); E[i, j] = 1; basis.append(E)
    rows = []
    for E in basis:
        c1 = (E @ Q - Q @ E).ravel(); c3 = (E @ G - G @ E).ravel()   # covariant generator: commutes with Q, Q^GF
        c2 = (E @ K0 + np.array(one_tensor(E.tolist(), K0.tolist(), par, 0))).ravel()
        rows.append(np.concatenate([c1, c2, c3]))
    A = np.array(rows).T
    u, sv, vt = np.linalg.svd(A); null = vt[np.sum(sv > 1e-10):]
    rng = np.random.default_rng(seed)
    M = sum(c * sum(v * E for v, E in zip(vec, basis)) for c, vec in zip(rng.normal(size=len(null)), null))
    return dict(Q=Q, G=G, K0=K0, par=par, Pi=Pi, H=H, Ps=Ps, k=k, M=M, h=h, eps=eps, T=T, n=n, null_dim=len(null))

def run_B(signs, label):
    out = {'label': label}
    m = model_B(1.3, signs)
    Q, G, M, Ps, k, K0, Pi, par = m['Q'], m['G'], m['M'], m['Ps'], m['k'], m['K0'], m['Pi'], m['par']
    out['dim of even pairing-preserving operators commuting with Q and Q^GF'] = m['null_dim']
    out['|[M,Q^GF]|'] = float(np.abs(M @ G - G @ M).max())
    out['H - s^2'] = float(np.abs(m['H'] - m['h'] * np.eye(m['n'])).max())
    Psi = Ps @ M @ Ps
    out['|Psi| = |P*MP*|'] = float(np.abs(Psi).max())
    out['[Q,Psi]'] = float(np.abs(Q @ Psi - Psi @ Q).max())
    Phi = np.linalg.inv(k) @ Ps @ Psi
    out['Psi - (Q Phi + Phi Q)'] = float(np.abs(Psi - (Q @ Phi + Phi @ Q)).max())
    # functional level: d_Psi = d_Phi Q + Q d_Phi with kernels (Psi x 1)K_0 (even), (Pi Phi x 1)K_0 (odd)
    alg = Alg(par, zero=0.0); rng = random.Random(11)
    Psik = (Psi @ K0).tolist(); Phik = (Pi @ Phi @ K0).tolist(); Ql = Q.tolist()
    Qd = lambda f: alg.D(Ql, f)
    err = 0.0
    for p in (0, 1):
        F = alg.rand([2, 3], 6, rng, p, lambda r: r.uniform(-1, 1))
        err = max(err, alg.dist(alg.dphi(Psik, F), alg.add(alg.dphi(Phik, Qd(F)), Qd(alg.dphi(Phik, F)))))
    out['d_Psi - (d_Phi Q + Q d_Phi) on functionals'] = err
    # large-distance (small H = s^2) behaviour of the composite line and of its primitive
    table = []
    for s in (1.0, 0.1, 0.01, 0.001):
        row = {'s': s}
        for T in (np.inf, 100.0):
            mb = model_B(s, signs, T=T)
            Psi_b = mb['Ps'] @ mb['M'] @ mb['Ps']
            nk = mb['k'][0, 0]
            row['|(P*)^2| T=%s' % T] = float(np.abs(mb['Ps'] @ mb['Ps']).max())
            row['|Phi| T=%s' % T] = float(np.abs(mb['Ps'] @ Psi_b).max() / nk) if nk != 0 else None
        row['|(Q^GF)^2|'] = float(np.abs(mb['G'] @ mb['G']).max()); row['|(Q^GF)^3|'] = float(np.abs(mb['G'] @ mb['G'] @ mb['G']).max())
        table.append(row)
    out['small-momentum table'] = table
    return out

def main():
    results = {'model A (exact), Feynman': run_A(FEYN, label='Feynman (Q^GF)^2 != 0'),
               'model A (exact), Landau': run_A(LANDAU, label='Landau (Q^GF)^2 = 0'),
               'model B (float), Feynman': run_B(FEYN, 'Feynman'),
               'model B (float), Landau': run_B(LANDAU, 'Landau')}
    with open(os.path.join(HERE, 'delta_toy_results.json'), 'w') as f:
        json.dump(results, f, indent=1)
    bad = []
    for k, r in results.items():
        for kk, v in r.items():
            if isinstance(v, float) and kk.startswith(('|', 'dim')) is False and 'table' not in kk and abs(v) > 1e-12:
                bad.append((k, kk, v))
    print(json.dumps(results, indent=1)); print('NONZERO RESIDUALS:', bad)

if __name__ == '__main__':
    main()
