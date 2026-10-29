"""Checks of the derived items [D] of step2/literature/ip_targets.tex (inner product / Q^dagger = S targets).

Run:  python3 check_ip_targets.py > check_ip_targets_output.txt
Small exact (sympy) and double-precision (numpy) computations; no external input files.

Blocks of SW (17): rows/columns ordered (1) = 2-dim, (2) = 2-dim | (3) = 4-dim; Delta = diag(1_2, -1_2 | 0_4).
Project table (ch02): (2,1) translations, (1,2) special conformal, (3,1) = S_+ (x) W (Q^i_alpha),
(2,3) = S_- (x) W^* (Qbar_{alphadot i}), (1,3) = S_+ (x) W^* (S^alpha_i), (3,2) = S_- (x) W (Sbar^{alphadot i}).

C1  B03 Table 1: hermiticity S_2 = Q_2^dagger (adjoint of B03 (3.9), structure i -> s_i structure i) as conditions on gamma_k.
C2  B03 (4.13) (corrected entry, n4_untwisted D5) is self-adjoint for the B03 Sec.2.1 tree scalar product.
C3  one-loop hermiticity identity  X#_1 = X_1^{dag0} + [X#_0, M],  M = G0^{-1} G1, and its frame covariance.
C4  conformal sl(2) toy (radial quantisation, P uncorrected): K_1 = [K_0, M], 2 D_1 = [[K_0, M], P_0].
C5  supermatrix conjugations X -> H X^dagger H^{-1} on gl(4|4): anti-automorphism, block map, positivity pattern,
    and transport of the Beisert oscillator conjugation (hep-th/0407277 (D.16)-(D.19)) by the supertranspose.
C6  inversion R(x) = x/(x,x): pushforward of the conformal Killing fields; v -> -R_*(conj v) is an antilinear
    anti-automorphism of the vector-field bracket fixing Delta = 2E.
C7  DFS (2.6): sum_p [[T^a,T^p],T^p] = N T^a for SU(N), tr T^a T^b = delta/2; unit conversions.
"""
import itertools
import numpy as np
import sympy as sp

rng = np.random.default_rng(20261028)
out = []


def say(*a):
    s = ' '.join(str(x) for x in a)
    print(s)


# ----------------------------------------------------------------------------------------------------------
say('=== C1  B03 Table 1: S_2 = (Q_2)^dagger')
a1, b1 = sp.symbols('alpha1 beta1', real=True)
gr = sp.symbols('g1r g2r g3r g4r', real=True)
gi = sp.symbols('g1i g2i g3i g4i', real=True)
gam = [gr[k] + sp.I * gi[k] for k in range(4)]
q = sp.Rational(1, 4)
Q2 = [q*a1**2 - sp.I*gam[2]/2 + sp.I*gam[3]/2, sp.I*gam[2]/2 + sp.I*gam[3]/2, -q*a1**2 - sp.I*gam[0], q*a1**2 + sp.I*gam[1]]
S2 = [q*a1**2 + sp.I*gam[2]/2 - sp.I*gam[3]/2, -sp.I*gam[2]/2 - sp.I*gam[3]/2, -q*a1**2 + sp.I*gam[0], q*a1**2 - sp.I*gam[1]]
# phases: Q_2 carries e^{+i beta1}, S_2 carries e^{-i beta1}; conj(e^{i beta1}) = e^{-i beta1} for real beta1.
for signs in itertools.product([1, -1], repeat=4):
    eqs = []
    for k in range(4):
        d = sp.expand((sp.conjugate(Q2[k]) - signs[k]*S2[k]).subs(a1, 1))   # alpha1 = 1 (alpha1 != 0 rescales)
        eqs += [sp.re(d), sp.im(d)]
    const = [e for e in eqs if e.free_symbols == set()]
    if any(c != 0 for c in const):
        say('  signs', signs, '-> no solution (gamma-independent equation', [c for c in const if c != 0], ')')
        continue
    sol = sp.solve([e for e in eqs if e.free_symbols], list(gr) + list(gi), dict=True)
    say('  signs', signs, '->', sol if sol else 'no solution')
say('  gamma-independent sums Q2[i]+S2[i] =', [sp.simplify(x + y) for x, y in zip(Q2, S2)])
gQ0 = {gr[0]: 0, gr[1]: 0, gr[2]: 0, gr[3]: 0,
       gi[0]: a1**2/4, gi[1]: a1**2/4, gi[2]: -a1**2/4, gi[3]: a1**2/4}
say('  frame Q_2 = 0 (n4_untwisted D2): Q2 =', [sp.simplify(c.subs(gQ0)) for c in Q2],
    ' gamma imaginary -> not a hermitian frame')

# ----------------------------------------------------------------------------------------------------------
say('\n=== C2  B03 (4.13) self-adjoint for the B03 Sec.2.1 tree scalar product')
g, b2 = sp.symbols('g beta2', real=True)
a1p, a3, eps = sp.symbols('alpha1 alpha3 epsilon', positive=True)
X = sp.symbols('X')
# convention of check_beisert_short_states.py: H O_i = sum_j M_ij O_j
M = sp.Matrix([[3 + eps*g**4, -2*sp.sqrt(2)*sp.exp(sp.I*b2)*a1p**3*g**3],
               [-9*sp.sqrt(2)*X*a1p**3*g**3, 3 + 6*a1p**2*g**2 - 18*a1p**4*g**4 + 12*a1p*a3*g**4]])
# tree norms: O_1 = eps^{ab}|psi_a psi_b> = 2 eps^{12}|psi_1 psi_2>, |psi_1 psi_2> norm 1 (n = 1);
# O_2 = eps^{abc}|phi_a phi_b phi_c> = 3(|123> - |132>), two cyclic classes of norm 1.
G11, G22 = 4, 9*2
say('  tree norms <O1|O1> =', G11, ' <O2|O2> =', G22)
Gm = sp.diag(G11, G22)
# self-adjointness: <O_i, H O_j> = <H O_i, O_j>  <=>  sum_k M_jk G_ik = sum_k conj(M_ik) G_kj
for Xv, name in [(sp.exp(-sp.I*b2), 'corrected (n4_untwisted D5)'), (1, 'printed e^{i b2} e^{-i b2}')]:
    Mv = M.subs(X, Xv)
    lhs = (Gm * Mv.T)
    rhs = (Mv.conjugate() * Gm)
    diff = sp.simplify(lhs - rhs)
    say('  M21 ', name, ': G M^T - conj(M) G =', diff)
c_prot = sp.sqrt(2)/3*sp.exp(sp.I*b2)*a1p   # protected: O1 + c g O2
d_kon = -3/sp.sqrt(2)*sp.exp(-sp.I*b2)*a1p  # Konishi descendant: O2 + d g O1
ip = sp.simplify(sp.expand(sp.conjugate(c_prot*g)*G22*1 + sp.conjugate(1)*G11*(d_kon*g)))
say('  <O1 + c g O2 , O2 + d g O1>_tree at order g =', ip)

# ----------------------------------------------------------------------------------------------------------
say('\n=== C3  one-loop hermiticity identity and frame covariance (random matrices, n = 7)')
n = 7


def rand_c(n):
    return rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))


A = rand_c(n)
G0 = A.conj().T @ A + n*np.eye(n)          # hermitian positive
G1h = rand_c(n)
G1 = G1h + G1h.conj().T                      # hermitian
X0, X1 = rand_c(n), rand_c(n)
G0i = np.linalg.inv(G0)


def dag0(Y):
    return G0i @ Y.conj().T @ G0


Mm = G0i @ G1
# exact X#(h) = G(h)^{-1} X(h)^dag G(h), first-order coefficient by complex-free central difference + Richardson
def Xsharp(h):
    G = G0 + h*G1
    Xh = X0 + h*X1
    return np.linalg.solve(G, Xh.conj().T @ G)


hs = 1e-4
der = (8*(Xsharp(hs) - Xsharp(-hs)) - (Xsharp(2*hs) - Xsharp(-2*hs)))/(12*hs)
pred = dag0(X1) + (dag0(X0) @ Mm - Mm @ dag0(X0))
say('  |d/dh X#(0) - (X1^{dag0} + [X#_0, M])| =', float(np.max(np.abs(der - pred))))
say('  |X#_0 - X0^{dag0}| =', float(np.max(np.abs(Xsharp(0.0) - dag0(X0)))))
# frame change T = 1 + h Z: X' = T X T^{-1}, G' = T^{-dag} G T^{-1}
Z = rand_c(n)


def first_order(f, h=1e-4):
    return (8*(f(h) - f(-h)) - (f(2*h) - f(-2*h)))/(12*h)


Xp1 = first_order(lambda h: (np.eye(n) + h*Z) @ (X0 + h*X1) @ np.linalg.inv(np.eye(n) + h*Z))
Gp1 = first_order(lambda h: np.linalg.inv(np.eye(n) + h*Z).conj().T @ (G0 + h*G1) @ np.linalg.inv(np.eye(n) + h*Z))
Mp = G0i @ Gp1
say('  X1\' - (X1 + [Z, X0])                 :', float(np.max(np.abs(Xp1 - (X1 + Z @ X0 - X0 @ Z)))))
say('  M\'  - (M - Z - Z^{dag0})             :', float(np.max(np.abs(Mp - (Mm - Z - dag0(Z))))))
lhs_new = first_order(lambda h: np.linalg.solve(G0 + h*Gp1, (X0 + h*Xp1).conj().T @ (G0 + h*Gp1)))
pred_new = dag0(Xp1) + (dag0(X0) @ Mp - Mp @ dag0(X0))
say('  identity in the new frame             :', float(np.max(np.abs(lhs_new - pred_new))))
say('  X#_1\' - (X#_1 + [Z, X#_0])           :', float(np.max(np.abs(lhs_new - (der + Z @ dag0(X0) - dag0(X0) @ Z)))))

# ----------------------------------------------------------------------------------------------------------
say('\n=== C4  conformal sl(2) toy: [K,P] = 2D, [D,P] = P, P^dag = K, lowest weight Delta0 + h*gam, P uncorrected')
h, gm, D0 = sp.symbols('h gamma Delta0', positive=True)
L = 6
Dl = D0 + h*gm
P = sp.zeros(L, L)
K = sp.zeros(L, L)
D = sp.zeros(L, L)
for m in range(L):
    D[m, m] = Dl + m
    if m + 1 < L:
        P[m + 1, m] = 1                      # P|m> = |m+1>
    if m >= 1:
        K[m - 1, m] = m*(2*Dl + m - 1)       # K|m> = m(2 Delta + m - 1)|m-1>
Gt = sp.diag(*[sp.factorial(m)*sp.rf(2*Dl, m) for m in range(L)])   # <m|m> = m! (2 Delta)_m
# hermiticity of the toy at finite h: K = G^{-1} P^dag G
say('  K - G^{-1} P^T G =', sp.simplify(K - Gt.inv()*P.T*Gt) == sp.zeros(L, L))
G0t = Gt.subs(h, 0)
G1t = sp.diff(Gt, h).subs(h, 0)
Mt = sp.simplify(G0t.inv()*G1t)
K0, K1 = K.subs(h, 0), sp.diff(K, h).subs(h, 0)
P0, D1 = P, sp.diff(D, h).subs(h, 0)
say('  M (diagonal)       =', [sp.simplify(Mt[m, m]) for m in range(L)])
say('  K_1 == [K_0, M]    :', sp.simplify(K1 - (K0*Mt - Mt*K0)) == sp.zeros(L, L))
comm = sp.simplify((K0*Mt - Mt*K0)*P0 - P0*(K0*Mt - Mt*K0))
# truncation: the last row/column of P*K products are cut at level L-1; compare on levels 0..L-2
say('  2 D_1 == [[K_0,M],P_0] on levels 0..L-2 :',
    sp.simplify((comm - 2*D1)[:L-1, :L-1]) == sp.zeros(L-1, L-1))
say('  level-one ratio M_1 - M_0 = gamma/Delta0 :', sp.simplify(Mt[1, 1] - Mt[0, 0] - gm/D0) == 0)

# ----------------------------------------------------------------------------------------------------------
say('\n=== C5  conjugations X -> H X^dagger H^{-1} on gl(4|4)')
par = np.array([0, 0, 0, 0, 1, 1, 1, 1])
blk = np.array([1, 1, 2, 2, 3, 3, 3, 3])


def sbr(Xa, pa, Xb, pb):
    return Xa @ Xb - (-1)**(pa*pb) * Xb @ Xa


def homog(p):
    Y = rand_c(8)
    mask = np.array([[(par[i] + par[j]) % 2 == p for j in range(8)] for i in range(8)])
    return Y*mask


def conj(Hm):
    Hi = np.linalg.inv(Hm)
    return lambda Y: Hm @ Y.conj().T @ Hi


Delta = np.diag([1, 1, -1, -1, 0, 0, 0, 0]).astype(complex)
sgn = np.diag([1, 1, 1, 1, -1, -1, -1, -1]).astype(complex)


def str_(Y):
    return np.trace(sgn @ Y)


def delta_comp(Y):
    return (str_(Delta @ Y)/str_(Delta @ Delta)).real


def E(i, j):
    Y = np.zeros((8, 8), complex)
    Y[i, j] = 1
    return Y


cands = {'eta=diag(-1_2, 1_2 | 1_4)': np.diag([-1, -1, 1, 1, 1, 1, 1, 1]).astype(complex),
         'eta=diag( 1_2,-1_2 | 1_4)': np.diag([1, 1, -1, -1, 1, 1, 1, 1]).astype(complex),
         'eta=1 (plain dagger)     ': np.eye(8, dtype=complex)}
reps = {'(2,1) P': E(2, 0), '(1,2) K': E(0, 2), '(3,1) Q': E(4, 0), '(2,3) Qbar': E(2, 4),
        '(1,3) S': E(0, 4), '(3,2) Sbar': E(4, 2)}
for name, Hm in cands.items():
    c = conj(Hm)
    err = 0.0
    for pa, pb in itertools.product([0, 1], repeat=2):
        Ya, Yb = homog(pa), homog(pb)
        err = max(err, np.max(np.abs(c(sbr(Ya, pa, Yb, pb)) - sbr(c(Yb), pb, c(Ya), pa))))
    lin = np.max(np.abs(c(1j*E(4, 0)) + 1j*c(E(4, 0))))
    pattern = {k: round(delta_comp(sbr(v, int(k.split()[1] in ('Q', 'Qbar', 'S', 'Sbar')), c(v),
                                       int(k.split()[1] in ('Q', 'Qbar', 'S', 'Sbar')))), 4) for k, v in reps.items()}
    say('  ', name, ': max|[X,Y]^c - [Y^c,X^c]| =', float(err), ' antilinear:', float(lin) < 1e-14,
        ' Delta^c = Delta:', np.allclose(c(Delta), Delta))
    say('      Delta-component of [X, X^c} :', pattern)
# random block-diagonal hermitian H with the signature of eta=diag(-1_2,1_2|1_4)
Bb = [rand_c(2), rand_c(2), rand_c(4)]
Hr = np.zeros((8, 8), complex)
Hr[0:2, 0:2] = -(Bb[0].conj().T @ Bb[0] + np.eye(2))
Hr[2:4, 2:4] = Bb[1].conj().T @ Bb[1] + np.eye(2)
Hr[4:8, 4:8] = Bb[2].conj().T @ Bb[2] + np.eye(4)
c = conj(Hr)
odd_ok = all(delta_comp(sbr(v, 1, c(v), 1)) < 0 for v in [homog(1) * (np.outer(blk == 3, blk == 1)) for _ in range(5)])
say('   general H (block-diag hermitian, signature -,+|+): {Q,Q^c} Delta-component < 0 for random Q in (3,1):', odd_ok)
# transport of Beisert's oscillator conjugation: (J_A^B)^dag = eta_A eta_B J_B^A, eta_B = (+1_2 a, -1_2 b | +1_4 c),
# Beisert J-matrices are the supertranspose of the SW/project matrices.


def st(Y):
    Aa, Bb_, Cc, Dd = Y[:4, :4], Y[:4, 4:], Y[4:, :4], Y[4:, 4:]
    Z_ = np.zeros_like(Y)
    Z_[:4, :4], Z_[:4, 4:], Z_[4:, :4], Z_[4:, 4:] = Aa.T, Cc.T, -Bb_.T, Dd.T
    return Z_


def st_inv(Y):
    Aa, Bb_, Cc, Dd = Y[:4, :4], Y[:4, 4:], Y[4:, :4], Y[4:, 4:]
    Z_ = np.zeros_like(Y)
    Z_[:4, :4], Z_[:4, 4:], Z_[4:, :4], Z_[4:, 4:] = Aa.T, -Cc.T, Bb_.T, Dd.T
    return Z_


etaB = np.diag([1, 1, -1, -1, 1, 1, 1, 1]).astype(complex)
cB = conj(etaB)
cSW = conj(cands['eta=diag(-1_2, 1_2 | 1_4)'])
err_st = 0.0
for p in (0, 1):
    Y = homog(p)
    err_st = max(err_st, np.max(np.abs(st(cB(st_inv(Y))) - cSW(Y))))
errst = 0.0
for pa, pb in itertools.product([0, 1], repeat=2):
    Ya, Yb = homog(pa), homog(pb)
    errst = max(errst, np.max(np.abs(st(sbr(Ya, pa, Yb, pb)) - (-1)**(pa*pb)*sbr(st(Yb), pb, st(Ya), pa))))
say('   st([X,Y}) = (-1)^{|X||Y|} [st Y, st X}  (st reverses brackets up to the Koszul sign):', float(errst))
say('   st(D_B) with D_B = diag(1/2,-1/2|0) in the (a,b|c) blocks equals Delta/2:',
    np.allclose(st(np.diag([.5, .5, -.5, -.5, 0, 0, 0, 0]).astype(complex)), Delta/2))
say('   st o (Beisert conjugation) o st^{-1} == conjugation with eta=diag(-1_2,1_2|1_4):', float(err_st) < 1e-12)

# ----------------------------------------------------------------------------------------------------------
say('\n=== C6  inversion pushforward on conformal Killing fields of R^4')
x = sp.symbols('x1:5', real=True)
y = sp.symbols('y1:5', real=True)
xs = sp.Matrix(x)
r2 = sum(xi**2 for xi in x)
Rmap = xs / r2
Jac = Rmap.jacobian(xs)


def push(vfun):
    """pushforward of the vector field v(x) under R: (R_* v)(y) = dR_{x}(v(x)) at x = R^{-1}(y) = y/(y,y)."""
    vx = sp.Matrix(vfun(x))
    w = Jac * vx
    ys2 = sum(yi**2 for yi in y)
    sub = {x[i]: y[i]/ys2 for i in range(4)}
    return sp.simplify(w.subs(sub))


def as_y(vfun):
    return sp.Matrix(vfun(y))


e1 = [1, 0, 0, 0]
P_e = lambda z: e1
K_e = lambda z: [sum(zi**2 for zi in z)*e1[i] - 2*sum(e1[k]*z[k] for k in range(4))*z[i] for i in range(4)]
Ef = lambda z: list(z)
om = sp.Matrix([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]])
Rot = lambda z: list(om*sp.Matrix(z))
say('  R_*(P_e) == K_e           :', sp.simplify(push(P_e) - as_y(K_e)) == sp.zeros(4, 1))
say('  R_*(K_e) == P_e           :', sp.simplify(push(K_e) - as_y(P_e)) == sp.zeros(4, 1))
say('  R_*(E)   == -E            :', sp.simplify(push(Ef) + as_y(Ef)) == sp.zeros(4, 1))
say('  R_*(rotation) == rotation :', sp.simplify(push(Rot) - as_y(Rot)) == sp.zeros(4, 1))


def vbr(v, w, z):
    """vector-field bracket [v,w] = v(w) - w(v) (manuscript def:c2-even)."""
    V, W = sp.Matrix(v(z)), sp.Matrix(w(z))
    return list(W.jacobian(z)*V - V.jacobian(z)*W)


zz = sp.Matrix(y)
pairs = [(P_e, K_e), (Ef, P_e), (Ef, K_e), (Rot, K_e), (P_e, Rot)]
ok = True
for v, w in pairs:
    # dagger(v) = -R_*(v) for real v; anti-automorphism: dagger([v,w]) = [dagger(w), dagger(v)]
    lhs = -push(lambda z: vbr(v, w, z))
    dv, dw = -push(v), -push(w)
    rhs = sp.Matrix(vbr(lambda z: list(dw.subs({y[i]: z[i] for i in range(4)})),
                        lambda z: list(dv.subs({y[i]: z[i] for i in range(4)})), list(y)))
    ok = ok and sp.simplify(lhs - rhs) == sp.zeros(4, 1)
say('  v -> -R_*(v) reverses brackets on (P,K),(E,P),(E,K),(rot,K),(P,rot):', ok)
say('  hence Delta = 2E -> 2E, P_e -> -K_e, K_e -> -P_e, rotation -> -rotation (real fields)')

# ----------------------------------------------------------------------------------------------------------
say('\n=== C7  colour identity of DFS (2.6) and unit conversions')
for N in (2, 3, 4):
    # generalized Gell-Mann basis, tr T^a T^b = delta/2
    gens = []
    for i in range(N):
        for j in range(i + 1, N):
            S_ = np.zeros((N, N), complex); S_[i, j] = S_[j, i] = 0.5; gens.append(S_)
            A_ = np.zeros((N, N), complex); A_[i, j] = -0.5j; A_[j, i] = 0.5j; gens.append(A_)
    for k in range(1, N):
        Dg = np.zeros((N, N), complex)
        Dg[:k, :k] = np.eye(k); Dg[k, k] = -k
        gens.append(Dg/np.sqrt(2*k*(k + 1)))
    err = 0.0
    for Ta in gens:
        cas = sum((Ta @ Tp - Tp @ Ta) @ Tp - Tp @ (Ta @ Tp - Tp @ Ta) for Tp in gens)
        err = max(err, np.max(np.abs(cas - N*Ta)))
    say('  SU(%d): max|sum_p [[T^a,T^p],T^p] - N T^a| =' % N, float(err))
gN = sp.symbols('g N', positive=True)
unit = -gN[0]**2*gN[1]/(32*sp.pi**2)
say('  unit (su(N), (X,Y) = -(4 hbar/g^2) tr XY) =', unit, '; -24 unit =', sp.simplify(-24*unit))
say('  B04 (2.20): 6 g_B^2 with g_B^2 = g^2 N/(8 pi^2) =', sp.simplify(6*gN[0]**2*gN[1]/(8*sp.pi**2)))
gstd = sp.symbols('g_std', positive=True)
lhs24 = -(3*gN[1]/(4*sp.pi**2))*3*(gN[1]**2 - 1)/(4*sp.pi**2)**2       # -gamma_1 x (24) prefactor
rhs27 = -9*gN[1]*(gN[1]**2 - 1)/(4*(4*sp.pi**2)**3)
say('  BKRS: order-g^2 log coefficient from (24) / printed (27) =', sp.simplify(lhs24/rhs27),
    '(consistent with g^2 in (27) = g_std^2/4)')

# ----------------------------------------------------------------------------------------------------------
say('\n=== C8  conjugation on the odd parameters (dictionary entries SPN-sigmaE, Gamma, pairing, brackets)')
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
sig = [s1, s2, s3, 1j*np.eye(2)]            # sigma_E^m[alpha, alphadot]
sigb = [-s1, -s2, -s3, 1j*np.eye(2)]        # sigmabar_E^m[alphadot, alpha]
say('  sigmabar_E = -sigma_E^dagger :', all(np.allclose(sb, -sg.conj().T) for sg, sb in zip(sig, sigb)))


def GamV_xi(xa, xab, ya, yab):
    """Gamma^m(xi_1, xi_2) = i( xibar_2 sigmabar^m xi_1 + xibar_1 sigmabar^m xi_2 ); xi = (xi_alpha^i, xibar_{alphadot i})."""
    return np.array([1j*(np.einsum('ui,uv,vi->', yab, sigb[m], xa) + np.einsum('ui,uv,vi->', xab, sigb[m], ya)) for m in range(4)])


def GamV_eta(ea, eab, fa, fab):
    """Gamma^m(eta_1, eta_2) = -i( eta_2 sigma^m etabar_1 + eta_1 sigma^m etabar_2 ); eta = (eta^alpha_i, etabar^{alphadot i})."""
    return np.array([-1j*(np.einsum('vi,vu,ui->', fa, sig[m], eab) + np.einsum('vi,vu,ui->', ea, sig[m], fab)) for m in range(4)])


def dag(xa, xab, f_und, f_dot):
    """xi -> eta with eta^alpha_i = -f_und conj(xi_alpha^i), etabar^{alphadot i} = f_dot conj(xibar_{alphadot i})."""
    return -f_und*xa.conj(), f_dot*xab.conj()


def rnd(shape):
    return rng.normal(size=shape) + 1j*rng.normal(size=shape)


x1a, x1b, x2a, x2b = rnd((2, 4)), rnd((2, 4)), rnd((2, 4)), rnd((2, 4))
for fu, fd in [(1.0, 1.0), (2.0, 0.5), (2.0, 1.0), (1.0, -1.0)]:
    e1 = dag(x1a, x1b, fu, fd)
    e2 = dag(x2a, x2b, fu, fd)
    lhs = GamV_eta(e2[0], e2[1], e1[0], e1[1])            # Gamma(xi_2^c, xi_1^c)|_V
    rhs = np.conj(GamV_xi(x1a, x1b, x2a, x2b))            # conj Gamma(xi_1, xi_2)|_V
    pair = np.einsum('ai,ai->', x1a, e1[0]) - np.einsum('ai,ai->', x1b, e1[1])   # (xi, xi^c)
    say('  factors (undotted, dotted) =', (fu, fd), ': Gamma(xi2^c,xi1^c)|_V == conj Gamma(xi1,xi2)|_V :',
        np.allclose(lhs, rhs), ' (xi1, xi1^c) =', np.round(pair, 6))
say('  => [xi_1,xi_2]^c = [xi_2^c, xi_1^c] with P_v^c = -K_{conj v} (C6) and prop:c2-special forces')
say('     (undotted factor) x (dotted factor) = 1; (xi, xi^c) < 0 for positive factors.')
