"""Route B (engine B), PART Q: slot space, the terms of S^(1)(x) written anew, and the operators x. on plane waves.

Reused verified code (read only): step2/composite_qsk/composite/engine_b/src.nosync/ce_core.py (54-component slot space,
edge matrices N, k0, qgf, q) and its vendor/ folder (Step-1 engine_b vertices.py, loops.py, realization.py: Clifford
realisation, action terms on plane waves, Hessian blocks).

New here: the terms of S^(1)(xi) (manuscript eq c2-S1Q, Table qs-terms), of S^(1)(eta) (eq c2-S1S, the intrinsic term
s_0), of S^(1)(v) for even generators (Def c7-vdot) as plane-wave vertex functions in the conventions of vertices.py
(Def c1-vertex: coefficient of the odd parameters in leg order), and the operators xi., eta., v. as 54 x 54 symbols.

Realisation facts used (realization.py): Sigma components 0..7 = S_+ (x) W (spinor a in {0,1}, W index A in 0..3,
component a*4+A), 8..15 = S_- (x) W^*; Sigma^* components 0..7 = S_- (x) W (lamsM), 8..15 = S_+ (x) W^* (lamsP).
Gamma(Q1,Q2)_M = (rho(e_M) Q1, Q2) = Q1^T MR_M Q2, MR_M = Mp R_M symmetric; (s, s*) = s^T Mp s*.
The labels lamP = components 0..7 are the -1 eigenspace of rho(e0 e1 e2 e3) (composite engine B revision note); with
the dictionary map of engine indices (0,1,2,3) -> (4,1,2,3) this is the +1 eigenspace of rho(e1 e2 e3 e4), i.e.
S_+ = undotted for the orientation dx^1 dx^2 dx^3 dx^4 of step2/notation/dictionary.json (checked in E0)."""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
QSK = os.path.abspath(os.path.join(HERE, '..', '..'))          # step2/composite_qsk/qsk
ROOT = os.path.abspath(os.path.join(QSK, '..'))                 # step2/composite_qsk
CEB = os.path.join(ROOT, 'composite', 'engine_b', 'src.nosync')
VEND = os.path.join(CEB, 'vendor')
for p in (VEND, CEB):
    if p not in sys.path:
        sys.path.insert(0, p)
import sympy as sp
import ce_core as C

V = C.V
L = C.L
I = sp.I
NT = C.NT
OFF = C.OFF
PARITY = C.PARITY
SPEC_OF = C.SPEC_OF
SLOT = C.SLOT
Mp = C.Mp
Mpinv = C.Mpinv
Rm = C.Rm            # rho(e_M): Sigma -> Sigma^*
Rbm = C.Rbm          # rho(e_M): Sigma^* -> Sigma
MR = [Mp * Rm[M] for M in range(10)]
PI = C.PI

XI = list(sp.symbols('xi0:16'))        # generic element of Sigma (components of the realisation)
ETA = list(sp.symbols('eta0:16'))      # generic element of Sigma^*


def slot_of(i):
    """slot name of a component index of the 54-dim slot space."""
    s = SPEC_OF[i]
    r = i - OFF[s]
    if s == 'A':
        return 'AV' if r < 4 else 'AW'
    if s == 'As':
        return 'AsV' if r < 4 else 'AsW'
    if s == 'lam':
        return 'lamP' if r < 8 else 'lamM'
    if s == 'lams':
        return 'lamsM' if r < 8 else 'lamsP'
    return s


# ------------------------------------------------------------------ Clifford helpers
def rho_two(Fmn, xi):
    """rho(F) xi for a two-form with components Fmn[M][N] (antisymmetric): (1/2) sum_{MN} F_MN rho(e_M) rho(e_N) xi."""
    v = sp.zeros(16, 1)
    xv = sp.Matrix(xi)
    for M in range(10):
        for N in range(10):
            if Fmn[M][N] != 0:
                v += sp.Rational(1, 2) * Fmn[M][N] * (Rbm[M] * (Rm[N] * xv))
    return v


def Gam(Q1, Q2):
    q1 = sp.Matrix(Q1); q2 = sp.Matrix(Q2)
    return [sp.expand((q1.T * MR[M] * q2)[0]) for M in range(10)]


# ------------------------------------------------------------------ terms of S^(1)(xi), xi in Sigma (eq c2-S1Q)
def _dA(a):
    """(dA)_MN of a plane-wave leg a (d along R^4 only)."""
    F = [[0] * 10 for _ in range(10)]
    for M in range(10):
        for N in range(10):
            v = 0
            if M < 4:
                v += I * a['p'][M] * a['u'][N]
            if N < 4:
                v -= I * a['p'][N] * a['u'][M]
            F[M][N] = v
    return F


def register_xi_terms(xi, tag='x'):
    """register T1+T2 (-(Gamma(xi,lam),A*)), T3+T4 (1/2 (rho(dA) xi, lam*)), T5 (1/2 (rho(1/2[A^A]) xi, lam*)) in
    vertices.TERMS under the names S1lam_<tag>, S1F2_<tag>, S1F3_<tag>; xi a list of 16 entries (sympy)."""
    xv = sp.Matrix(xi)
    MRxi = [(xv.T * MR[M]) for M in range(10)]          # row vectors: Gamma(xi, .)_M

    def L_lam(s):      # slots (lam, As): sum_M Gamma(xi,lam)_M A*_M
        l, As = s
        u = sp.Matrix(l['u'])
        return sum(((MRxi[M] * u)[0] * As['u'][M] for M in range(10)), 0)

    def L_F2(s):       # slots (A, lams): (rho(dA) xi, lam*)
        a, ls = s
        w = rho_two(_dA(a), xi)
        return (w.T * Mp * sp.Matrix(ls['u']))[0]

    RRxi = [[(Rbm[M] * (Rm[N] * xv)) for N in range(10)] for M in range(10)]

    def L_F3(s):       # slots (A, A, lams), colour ([c1,c2], c3): sum_{MN} u_M u'_N (rho_M rho_N xi, lam*)
        a, b, ls = s
        lv = Mp * sp.Matrix(ls['u'])
        tot = 0
        for M in range(10):
            if a['u'][M] == 0:
                continue
            for N in range(10):
                if b['u'][N] == 0 or M == N:
                    continue
                tot += a['u'][M] * b['u'][N] * (RRxi[M][N].T * lv)[0]
        return tot
    V.TERMS['S1lam_' + tag] = (sp.Integer(-1), ('lam', 'As'), L_lam, ('k', (0, 1)))
    V.TERMS['S1F2_' + tag] = (sp.Rational(1, 2), ('A', 'lams'), L_F2, ('k', (0, 1)))
    V.TERMS['S1F3_' + tag] = (sp.Rational(1, 4), ('A', 'A', 'lams'), L_F3, ('f', (0, 1, 2)))
    return ['S1lam_' + tag, 'S1F2_' + tag], ['S1F3_' + tag]


def register_s0_term(eta, tag='e'):
    """s_0 = -(rho(A|_perp) eta, lam*) for eta in Sigma^* (eq c2-S1S): slots (A, lams), A|_perp only."""
    ev = sp.Matrix(eta)

    def L_s0(s):
        a, ls = s
        w = sp.zeros(16, 1)
        for M in range(4, 10):
            if a['u'][M] != 0:
                w += a['u'][M] * (Rbm[M] * ev)
        return (w.T * Mp * sp.Matrix(ls['u']))[0]
    V.TERMS['S1s0_' + tag] = (sp.Integer(-1), ('A', 'lams'), L_s0, ('k', (0, 1)))
    return ['S1s0_' + tag]


def register_translation_term(b, tag='P'):
    """S^(1)(P_b) = int ((v.A, A*) - (v.lam, lam*) - (v.c, c*)), v = b constant: v.f = b.d f (Def c7-vdot)."""
    def bd(x):
        return sum((I * b[m] * x['p'][m] for m in range(4)), 0)

    def L_A(s):
        a, As = s
        return bd(a) * sum((a['u'][M] * As['u'][M] for M in range(10)), 0)

    def L_l(s):
        l, ls = s
        return bd(l) * (sp.Matrix(l['u']).T * Mp * sp.Matrix(ls['u']))[0]

    def L_c(s):
        c, cs = s
        return bd(c) * c['u'][0] * cs['u'][0]
    V.TERMS['SPA_' + tag] = (sp.Integer(1), ('A', 'As'), L_A, ('k', (0, 1)))
    V.TERMS['SPl_' + tag] = (sp.Integer(-1), ('lam', 'lams'), L_l, ('k', (0, 1)))
    V.TERMS['SPc_' + tag] = (sp.Integer(-1), ('c', 'cs'), L_c, ('k', (0, 1)))
    return ['SPA_' + tag, 'SPl_' + tag, 'SPc_' + tag]


# ------------------------------------------------------------------ the operators on plane waves (Lemma qs-arrows)
def xi_op_fields(k, xi):
    """field arrows of xi. on a plane wave e^{ikx}: lam -> -Gamma(xi,lam) (slot A), A -> -(1/2) rho(dA) xi (slot lam);
    returned as a 54 x 54 matrix (row = target component, column = source component), antifield arrows empty."""
    X = sp.zeros(NT, NT)
    xv = sp.Matrix(xi)
    for M in range(10):
        row = (xv.T * MR[M])
        for j in range(16):
            if row[j] != 0:
                X[OFF['A'] + M, OFF['lam'] + j] = -row[j]
    for N in range(10):
        leg = dict(u=[1 if M == N else 0 for M in range(10)], p=list(k))
        w = rho_two(_dA(leg), xi)
        for j in range(16):
            if w[j] != 0:
                X[OFF['lam'] + j, OFF['A'] + N] = sp.expand(-sp.Rational(1, 2) * w[j])
    return X


def complete_antifield_arrows(Xf, k, pairs, name='y'):
    """solve (X (x) 1 + 1 (x) X) K_0 = 0 for the antifield blocks; pairs = list of (target species, source species)
    of the unknown antifield arrows; entries linear polynomials in k with unknown coefficients."""
    mk = [-x for x in k]
    unknowns = []
    X = Xf.copy()
    for (ts, ss) in pairs:
        for i in range(V.DIM[ts]):
            for j in range(V.DIM[ss]):
                cs = sp.symbols(f'{name}_{ts}_{ss}_{i}_{j}_0:5')
                unknowns += list(cs)
                X[OFF[ts] + i, OFF[ss] + j] = cs[0] + sum((cs[1 + m] * k[m] for m in range(4)), 0)
    K0 = C.k0_matrix()
    Xm = X.subs({k[m]: mk[m] for m in range(4)}, simultaneous=True)
    E = (X * K0 + PI * K0 * Xm.T).applyfunc(sp.expand)
    eqs = []
    for v in E:
        if v != 0:
            P = sp.Poly(v, *k)
            eqs += [c for c in P.coeffs()]
    sol = sp.solve(eqs, unknowns, dict=True)
    assert len(sol) == 1
    X = X.subs(sol[0]).subs({u: 0 for u in unknowns})
    return X.applyfunc(sp.expand), sol[0]


def graded_comm(A, B, pa, pb):
    return (A * B - (-1) ** (pa * pb) * B * A).applyfunc(sp.expand)


def symbol_derivative(D, k, mu):
    """sigma^mu(D) = [D, x^mu] on plane waves: x^mu e^{ikx} = -i d/dk_mu e^{ikx}, hence [D, x^mu] = -i dD/dk_mu."""
    return D.applyfunc(lambda e: sp.expand(-I * sp.diff(e, k[mu])))
