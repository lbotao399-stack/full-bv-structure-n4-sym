"""Engine B explicit realization (code-only bases).
4d: two-component blocks sigma: S_+ -> S_-, sigmabar: S_- -> S_+.
6d: W = C^4, W_R^perp realised on Wedge^2: antisymmetric 4x4 blocks Sig^i: W -> W^*, Sigb^i: W^* -> W ('t Hooft symbols).
10d: Gamma^mu = gamma4^mu (x) 1_8, Gamma^{4+i} = gamma5 (x) gamma6^i ;  Sigma = S_+(x)W + S_-(x)W^*, Sigma^* = S_-(x)W + S_+(x)W^*.
Lie algebra: su(2), su(3) with deliberately non-orthonormal bases and a rescaled invariant pairing."""
import itertools
import sympy as sp

I = sp.I
Z = sp.zeros
def eye(n): return sp.eye(n)
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]])
sig = [eye(2), I * s1, I * s2, I * s3]        # S_+ -> S_-
sigb = [eye(2), -I * s1, -I * s2, -I * s3]    # S_- -> S_+

def block(a, b, c, d):
    return sp.Matrix(sp.BlockMatrix([[a, b], [c, d]]))

gamma4 = [block(Z(2, 2), sigb[m], sig[m], Z(2, 2)) for m in range(4)]   # on (S_+, S_-)
gamma5 = sp.diag(1, 1, -1, -1)

def thooft(selfdual):
    out = []
    for a in range(3):
        M = Z(4, 4)
        for m in range(3):
            for n in range(3):
                M[m, n] = sp.LeviCivita(a, m, n)
        s = 1 if selfdual else -1
        M[a, 3] = s; M[3, a] = -s
        out.append(M)
    return out
eta = thooft(True); etab = thooft(False)
Sig = [eta[0], eta[1], eta[2], I * etab[0], I * etab[1], I * etab[2]]      # W -> W^*  (antisymmetric)
Sigb = [-eta[0], -eta[1], -eta[2], I * etab[0], I * etab[1], I * etab[2]]  # W^* -> W
gamma6 = [block(Z(4, 4), Sigb[i], Sig[i], Z(4, 4)) for i in range(6)]      # on (W, W^*)
gamma7 = sp.diag(*([1] * 4 + [-1] * 4))

def kron(A, B):
    return sp.Matrix(sp.BlockMatrix([[A[i, j] * B for j in range(A.cols)] for i in range(A.rows)]))

Gam = [kron(gamma4[m], eye(8)) for m in range(4)] + [kron(gamma5, gamma6[i]) for i in range(6)]
Gam11 = kron(gamma5, gamma7)
# index of the 32-dim space: a*8+A, a in 0..3 (0,1 = S_+, 2,3 = S_-), A in 0..7 (0..3 = W, 4..7 = W^*)
SIG = [a * 8 + A for a in (0, 1) for A in range(4)] + [a * 8 + A for a in (2, 3) for A in range(4, 8)]
SIGS = [a * 8 + A for a in (2, 3) for A in range(4)] + [a * 8 + A for a in (0, 1) for A in range(4, 8)]
R = [Gam[M].extract(SIGS, SIG) for M in range(10)]    # rho(e_M): Sigma -> Sigma^*
Rb = [Gam[M].extract(SIG, SIGS) for M in range(10)]   # rho(e_M): Sigma^* -> Sigma

def pairing():
    """(s, s*) = s^T Mp s*;  require Mp R_M symmetric and Rb_M^T Mp symmetric for all M."""
    xs = sp.symbols('m0:256'); Mp = sp.Matrix(16, 16, xs)
    eqs = []
    for M in range(10):
        A = Mp * R[M]; eqs += list(A - A.T)
        B = Rb[M].T * Mp; eqs += list(B - B.T)
    sol = sp.solve(eqs, xs, dict=True)[0]
    Ms = Mp.subs(sol)
    free = sorted(Ms.free_symbols, key=lambda s: s.name)
    assert len(free) == 1, free
    Ms = Ms.subs(free[0], 1)
    # normalise so that the first nonzero entry is +1
    nz = [v for v in Ms if v != 0]
    return Ms / nz[0]

def checks(Mp):
    out = {}
    ok = all((Gam[M] * Gam[N] + Gam[N] * Gam[M]) == 2 * (1 if M == N else 0) * eye(32) for M in range(10) for N in range(10))
    out['clifford_10d'] = ok
    out['Gamma11_anticommutes'] = all(Gam11 * Gam[M] + Gam[M] * Gam11 == Z(32, 32) for M in range(10))
    out['Sigma_is_+1'] = all(Gam11[i, i] == 1 for i in SIG) and all(Gam11[i, i] == -1 for i in SIGS)
    out['Sig_antisymmetric'] = all(S.T == -S for S in Sig) and all(S.T == -S for S in Sigb)
    out['pairing_symmetry'] = all((Mp * R[M]).T == Mp * R[M] and (Rb[M].T * Mp).T == Rb[M].T * Mp for M in range(10))
    out['pairing_nondegenerate'] = Mp.det() != 0
    # Assumption 3.1 (10d Fierz) with random rational spinors
    import random
    random.seed(3)
    Qs = [sp.Matrix([random.randint(-3, 3) for _ in range(16)]) for _ in range(3)]
    def Gm(Q1, Q2):
        return [(Q1.T * Mp * R[M] * Q2)[0] for M in range(10)]
    def rho_of(v):
        return sum((v[M] * R[M] for M in range(10)), Z(16, 16))
    tot = Z(16, 1)
    for (a, b, c) in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]:
        tot += rho_of(Gm(Qs[a], Qs[b])) * Qs[c]
    out['assumption_3_1'] = tot == Z(16, 1)
    # chirality: rho(e_mu) mu<4 exchanges S_+xW <-> S_-xW ; rho(e_i) keeps S_+ (W -> W^*)
    return out

# ---------------- Lie algebras ----------------
def su2_data(scale=sp.Rational(3, 2)):
    t = [-I / 2 * s for s in (s1, s2, s3)]
    basis = [t[0], t[0] + t[1], t[2] + 2 * t[0]]        # non-orthonormal on purpose
    return lie_data(basis, scale)

def su3_data(scale=sp.Rational(5, 7)):
    l = [None] * 8
    l[0] = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 0]]); l[1] = sp.Matrix([[0, -I, 0], [I, 0, 0], [0, 0, 0]])
    l[2] = sp.Matrix([[1, 0, 0], [0, -1, 0], [0, 0, 0]]); l[3] = sp.Matrix([[0, 0, 1], [0, 0, 0], [1, 0, 0]])
    l[4] = sp.Matrix([[0, 0, -I], [0, 0, 0], [I, 0, 0]]); l[5] = sp.Matrix([[0, 0, 0], [0, 0, 1], [0, 1, 0]])
    l[6] = sp.Matrix([[0, 0, 0], [0, 0, -I], [0, I, 0]]); l[7] = sp.diag(1, 1, -2) / sp.sqrt(3)
    t = [-I / 2 * x for x in l]
    basis = [t[0] + t[3], t[1], t[2] - t[7], t[3], t[4] + 2 * t[5], t[5], t[6], t[7]]
    return lie_data(basis, scale)

def lie_data(basis, scale):
    """pairing (X,Y) = -2*scale*tr(XY) (positive definite on the compact form); structure constants; c_adj."""
    n = len(basis)
    kap = sp.Matrix(n, n, lambda a, b: sp.nsimplify(sp.simplify(-2 * scale * (basis[a] * basis[b]).trace())))
    kinv = kap.inv()
    def coords(X):
        v = sp.Matrix([sp.simplify(-2 * scale * (X * basis[b]).trace()) for b in range(n)])
        return kinv * v
    f = [[coords(basis[a] * basis[b] - basis[b] * basis[a]) for b in range(n)] for a in range(n)]  # [t_a,t_b] = f[a][b][c] t_c
    # ftil(a,b,c) = ([t_a,t_b],t_c)
    ft = sp.MutableDenseNDimArray.zeros(n, n, n)
    for a in range(n):
        for b in range(n):
            v = f[a][b]
            for c in range(n):
                ft[a, b, c] = sp.nsimplify(sum(v[d] * kap[d, c] for d in range(n)))
    # c_adj: sum_a [t_a,[t^a, X]] = c_adj X
    dual = [sum((kinv[a, b] * basis[b] for b in range(n)), Z(basis[0].rows, basis[0].rows)) for a in range(n)]
    X = basis[0]
    Y = sum(((ta * (td * X - X * td) - (td * X - X * td) * ta) for ta, td in zip(basis, dual)), Z(X.rows, X.rows))
    ratio = None
    for i in range(X.rows):
        for j in range(X.cols):
            if X[i, j] != 0:
                ratio = sp.nsimplify(sp.simplify(Y[i, j] / X[i, j])); break
        if ratio is not None:
            break
    assert sp.simplify(Y - ratio * X) == Z(X.rows, X.rows)
    return dict(n=n, basis=basis, kap=kap, kinv=kinv, ft=ft, cadj=ratio, f=f)

if __name__ == '__main__':
    Mp = pairing()
    print('pairing matrix nonzero pattern:', [(i, j, Mp[i, j]) for i in range(16) for j in range(16) if Mp[i, j] != 0][:20])
    print(checks(Mp))
    for name, L in (('su2', su2_data()), ('su3', su3_data())):
        print(name, 'dim', L['n'], 'c_adj =', L['cadj'])
