"""Audit (QFT structure): is the engines' translation rescaling P_1 an inner automorphism / universal?

Checks (exact or float with explicit tolerance; one process, < 0.2 GB):
 A. Colour-dipole identity on gauge singlets: sum_{j != i} T_i.T_j = -C_A on every gauge-invariant tensor of
    adjoint letters (su(2), su(3); lengths 2, 3, 4).  Hence a two-letter kernel kappa T_i.T_j (rho_i + rho_j)
    (rho_i = translation of letter i) acts on every singlet as -kappa C_A rho_0(P) and as 0 on a single letter:
    this reproduces the engine pattern (-1/3 unit on the scalar bilinears, 0 on letters) without being an inner
    automorphism.
 B. Flavour structure: a pair kernel (T_i.T_j)(kap_id + kap_perm P^flav_ij)(rho_i + rho_j) gives
    lambda = -C_A (kap_id + kap_perm) on every scalar bilinear (20' and Konishi: symmetric flavour) but
    lambda = -C_A (kap_id - kap_perm) on B_w-type trilinears (f_abc colour, antisymmetric flavour).
    The bilinear data (compare_delta.json: equal -1/3 on 20' and K) fix only kap_id + kap_perm.
 C. The inner-direction frame change U = 1 + c unit w acts on ALL local functionals (it is the field
    redefinition by phi ~ S^(1)(Delta), Thm ip-laws (a)): change of the P component on weight w is
    c (w - (w+1)) = -c for every operator, letters included.  With the engine values (P_1 = -(1/3+f) on
    bilinears, -f on letters), no c removes both.
 D. Arity two at (P_b, xi), [P_b, xi] = 0: with P-rescalings r_P^F on families F and a classical map
    rho_0(xi): F -> F' (nonzero), the order-unit part of [rho(P), rho(xi)] is (r_P^{F'} - r_P^F) rho_0(P) rho_0(xi),
    so r_P^{F'} = r_P^F: K -> O_xi (tau != 0) and O_zeta -> B_w (kappa_B != 0) force one P-rescaling on K, O_zeta,
    B_w.  With B: kap_perm = 0 is a prediction of the master equation, not a consequence of the bilinear data.
"""
import json, itertools
import numpy as np
import sympy as sp

OUT = {}

# ---------------------------------------------------------------- colour algebras (adjoint rep, real basis)
def su2():
    f = np.zeros((3, 3, 3))
    for p in itertools.permutations(range(3)):
        f[p] = np.linalg.det(np.eye(3)[list(p)])
    return f

def su3():
    # Gell-Mann structure constants f_abc (totally antisymmetric), basis T_a = lambda_a/2
    f = np.zeros((8, 8, 8))
    vals = {(0, 1, 2): 1.0, (0, 3, 6): 0.5, (0, 4, 5): -0.5, (1, 3, 5): 0.5, (1, 4, 6): 0.5,
            (2, 3, 4): 0.5, (2, 5, 6): -0.5, (3, 4, 7): np.sqrt(3) / 2, (5, 6, 7): np.sqrt(3) / 2}
    for (a, b, c), v in vals.items():
        for p in itertools.permutations(range(3)):
            idx = tuple([a, b, c][i] for i in p)
            f[idx] = v * np.linalg.det(np.eye(3)[list(p)])
    return f

def adj_gens(f):
    # (T_a)_{bc} = -i f_abc  (Hermitian); T_a.T_a = C_A
    return np.array([-1j * f[a] for a in range(f.shape[0])])

def apply_on_slot(T, t, slot):
    # (T t)_{..a..} = sum_b T_{a b} t_{..b..} on the given slot
    t2 = np.tensordot(T, t, axes=([1], [slot]))
    return np.moveaxis(t2, 0, slot)

def dipole(Tg, t, i, j):
    out = np.zeros_like(t, dtype=complex)
    for a in range(Tg.shape[0]):
        out += apply_on_slot(Tg[a], apply_on_slot(Tg[a], t, j), i)
    return out

def singlets(f, L):
    n = f.shape[0]
    d = None
    res = []
    if L == 2:
        res.append(('delta', np.eye(n)))
    if L == 3:
        res.append(('f', f.astype(complex)))
        if n == 8:
            # d_abc of su(3): d_abc = 2 tr({T_a,T_b} T_c) with T = lambda/2 (built numerically)
            lam = gell_mann()
            dd = np.zeros((8, 8, 8))
            for a in range(8):
                for b in range(8):
                    for c in range(8):
                        dd[a, b, c] = np.real(np.trace((lam[a] @ lam[b] + lam[b] @ lam[a]) @ lam[c])) / 4
            res.append(('d', dd.astype(complex)))
    if L == 4:
        e = np.eye(n)
        res.append(('delta_delta', np.einsum('ab,cd->abcd', e, e).astype(complex)))
        res.append(('ff', np.einsum('abe,cde->abcd', f, f).astype(complex)))
    return res

def gell_mann():
    l = np.zeros((8, 3, 3), dtype=complex)
    l[0][0, 1] = l[0][1, 0] = 1
    l[1][0, 1] = -1j; l[1][1, 0] = 1j
    l[2][0, 0] = 1; l[2][1, 1] = -1
    l[3][0, 2] = l[3][2, 0] = 1
    l[4][0, 2] = -1j; l[4][2, 0] = 1j
    l[5][1, 2] = l[5][2, 1] = 1
    l[6][1, 2] = -1j; l[6][2, 1] = 1j
    l[7] = np.diag([1, 1, -2]) / np.sqrt(3)
    return l

A = {}
for name, f in [('su2', su2()), ('su3', su3())]:
    Tg = adj_gens(f)
    CA = np.real(np.trace(sum(Tg[a] @ Tg[a] for a in range(Tg.shape[0])))) / Tg.shape[1]
    rows = []
    for L in (2, 3, 4):
        for tname, t in singlets(f, L):
            if np.linalg.norm(t) < 1e-12:
                continue
            # gauge invariance check
            gi = max(np.linalg.norm(sum(apply_on_slot(Tg[a], t, s) for s in range(L))) for a in range(Tg.shape[0]))
            per_letter = []
            for i in range(L):
                m = sum(dipole(Tg, t, i, j) for j in range(L) if j != i)
                lam = np.vdot(t, m) / np.vdot(t, t)
                resid = np.linalg.norm(m - lam * t) / np.linalg.norm(t)
                per_letter.append((float(np.real(lam)), float(resid)))
            rows.append({'L': L, 'tensor': tname, 'gauge_invariance_residual': float(gi),
                         'sum_j!=i Ti.Tj eigenvalue per letter': [p[0] for p in per_letter],
                         'eigen_residuals': [p[1] for p in per_letter],
                         'equals_minus_C_A': bool(all(abs(p[0] + CA) < 1e-10 and p[1] < 1e-10 for p in per_letter))})
    A[name] = {'C_A': float(CA), 'singlets': rows}
OUT['A_colour_dipole_on_singlets'] = A
OUT['A_conclusion'] = ('a two-letter kernel kappa T_i.T_j (rho_i+rho_j) acts on every gauge singlet as -kappa C_A rho_0(P), '
                       'on a single letter as 0 (no pair): the engine pattern (bilinears -1/3, letters 0) is a two-letter '
                       'effect, not an inner automorphism (which acts on letters with the same coefficient)')

# ---------------------------------------------------------------- B: flavour exchange component of the pair kernel
kid, kperm, CA_s = sp.symbols('kappa_id kappa_perm C_A')
# symmetric-flavour bilinear (20' or Konishi): P^flav_12 = +1, colour delta: T1.T2 = -C_A
lam_bilinear = -CA_s * (kid + kperm)
# B_w-type trilinear: colour f_abc (T1.(T2+T3) = -C_A), flavour antisymmetric c_[IJK]: P^flav_1j = -1
lam_B = -CA_s * (kid - kperm)
# numerical confirmation with su(2) colour f = eps, 3 flavours eps_IJK (antisymmetric) for B
f = su2(); Tg = adj_gens(f)
epsF = su2()  # eps_IJK on 3 flavours
tB = np.einsum('abc,IJK->aIbJcK', f, epsF).astype(complex)   # slots (a1 I1 a2 I2 a3 I3)
def flav_swap(t, i, j):
    # swap flavour indices of letters i, j (slots 2i+1, 2j+1)
    axes = list(range(t.ndim)); axes[2 * i + 1], axes[2 * j + 1] = axes[2 * j + 1], axes[2 * i + 1]
    return np.transpose(t, axes)
def dip(t, i, j):
    out = np.zeros_like(t)
    for a in range(3):
        out += apply_on_slot(Tg[a], apply_on_slot(Tg[a], t, 2 * j), 2 * i)
    return out
num = {}
for (ki, kp) in [(1.0, 0.0), (0.0, 1.0), (0.7, 0.3)]:
    m = sum(dip(ki * tB + kp * flav_swap(tB, 0, j), 0, j) for j in (1, 2))
    lam = np.vdot(tB, m) / np.vdot(tB, tB)
    num[str((ki, kp))] = {'lambda_B': float(np.real(lam)), 'predicted': float(-2.0 * (ki - kp)),
                          'residual': float(np.linalg.norm(m - lam * tB))}
# bilinear, su(2) colour delta, symmetric flavour (identity on 2 flavours, traceless part and trace)
for tag, F in [('trace', np.eye(3)), ('traceless', np.diag([1.0, -1.0, 0.0]))]:
    tK = np.einsum('ab,IJ->aIbJ', np.eye(3), F).astype(complex)
    for (ki, kp) in [(1.0, 0.0), (0.0, 1.0)]:
        m = dip(ki * tK + kp * flav_swap(tK, 0, 1), 0, 1)
        lam = np.vdot(tK, m) / np.vdot(tK, tK)
        num['bilinear_' + tag + str((ki, kp))] = {'lambda': float(np.real(lam)), 'predicted': float(-2.0 * (ki + kp))}
OUT['B_flavour_exchange'] = {'lambda_bilinears': str(lam_bilinear), 'lambda_B_w': str(lam_B),
                             'numerical_su2_C_A=2': num,
                             'conclusion': 'the scalar-bilinear value fixes kappa_id + kappa_perm only; P_1(B_w) needs '
                                           'kappa_id - kappa_perm, not determined by the bilinear data'}

# ---------------------------------------------------------------- C: inner-direction frame change on letters
c, fP = sp.symbols('c f_P')
weights = {'A (scalar or vector letter)': sp.Integer(1), 'lambda': sp.Rational(3, 2), 'c ghost': sp.Integer(0),
           'scalar bilinear': sp.Integer(2), 'O_zeta': sp.Rational(5, 2), 'B_w, Fbar_w': sp.Integer(3)}
# change of the coefficient r_P of {S1(P),O} under U = 1 + c unit w: c (w_O - w_PO) = -c (w_PO = w_O + 1)
chg = {k: str(sp.simplify(c * (w - (w + 1)))) for k, w in weights.items()}
P1_letters_after = sp.simplify(-fP - c).subs(c, -(sp.Rational(1, 3) + fP))
P1_bilin_after = sp.simplify(-(sp.Rational(1, 3) + fP) - c).subs(c, -(sp.Rational(1, 3) + fP))
OUT['C_inner_frame'] = {'change_of_r_P_on_weight_w': chg,
                        'engine_values_before': {'scalar bilinears': '-(1/3 + f_P)', 'letters': '-f_P'},
                        'after_c=r_P(bilinears)': {'scalar bilinears': str(P1_bilin_after),
                                                   'letters': str(sp.simplify(P1_letters_after))},
                        'conclusion': 'no inner-direction frame has P_1 = 0 on letters and on bilinears; the frame '
                                      'P_1 = 0 can hold at most on the gauge-invariant cohomology, and only if the '
                                      'two-letter coefficient is universal there'}

# ---------------------------------------------------------------- D: [P, xi] = 0 forces equal P-rescalings
u, rF, rFp = sp.symbols('u r_F r_Fp')
P0F, P0Fp, X0 = sp.symbols('P0F P0Fp X0', commutative=False)
# rho(P) on F, F': (1 + u r) P0 ; rho(xi): F -> F' : X0 + u X1 with [P0, X0] = 0 classically (P0Fp X0 = X0 P0F)
X1 = sp.Symbol('X1', commutative=False)
comm = sp.expand((1 + u * rFp) * P0Fp * (X0 + u * X1) - (X0 + u * X1) * (1 + u * rF) * P0F)
comm1 = sp.expand(comm.coeff(u, 1)).subs({P0Fp * X0: X0 * P0F})
# translation covariance of the one-loop xi component: P0Fp X1 = X1 P0F
comm1 = sp.expand(comm1.subs({P0Fp * X1: X1 * P0F}))
OUT['D_arity_two_P_xi'] = {'order_unit_part_of_[rho(P),rho(xi)]': str(sp.factor(comm1)),
                           'vanishing_requires': 'r_P(F\') = r_P(F) whenever the classical xi-map F -> F\' is nonzero',
                           'instances': ['K -> O_xi (tau = 2, Prop ip-tree)', 'O_zeta -> B_w (kappa_B = 1/4, Prop ip-weight3)'],
                           'consequence': 'one P-rescaling on K, O_zeta, B_w; with B this predicts kappa_perm = 0 and the '
                                          'same coefficient for the lambda-A and lambda-lambda pair kernels (O_zeta, Fbar_w)'}

ok = (all(r['equals_minus_C_A'] for g in A.values() for r in g['singlets'])
      and all(abs(v['lambda_B'] - v['predicted']) < 1e-10 for k, v in num.items() if k.startswith('('))
      and all(abs(v['lambda'] - v['predicted']) < 1e-10 for k, v in num.items() if k.startswith('bilinear'))
      and OUT['C_inner_frame']['after_c=r_P(bilinears)']['letters'] == '1/3')
OUT['all_checks_pass'] = bool(ok)
json.dump(OUT, open(__file__.replace('.py', '.json'), 'w'), indent=1)
print(json.dumps({'all_checks_pass': ok, 'B': OUT['B_flavour_exchange']['lambda_B_w'],
                  'C_letters_after': OUT['C_inner_frame']['after_c=r_P(bilinears)'],
                  'D': OUT['D_arity_two_P_xi']['order_unit_part_of_[rho(P),rho(xi)]']}, indent=1))
