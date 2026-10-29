"""QA compare_delta (Part D4): reassemble the arity-one two-point obstruction Ob(Delta) = Q(R)_sym from engine A's
per-graph two-point data (step2/engine_a/delta_results.json, task1) and compare with engine B's exact Ob(Delta).

Canonical form (both engines): a two-leg kernel K = alpha(p) delta(Sigma) + beta^mu(p) d_mu delta(Sigma), field leg at p
(held fixed), derivative along the second (antifield / Q-image) leg, p2 = -p on the support.
Q on the second leg: polarisation b -> h(p2) a2  =>  alpha_Q = alpha h(-p) - beta^mu (d_mu h)(-p),  beta_Q = beta h(-p).
Bosonic swap (both legs A): alpha'(p) = alpha(-p)^T - (div beta)^T(-p), beta'(p) = beta(-p)^T; kernel normalisation
K_sym = (K + K')/2 (the same normalisation as the target local functional's kernel, e.g. int(dA,dA) <-> |p|^2(a,b)).
Conventions used (engine B's Q_basis, ch01): Q int(A,A^*) = -1/2 int(dA,dA) - int(dc,A^*), Q int(c^*,c) = -int(dc,A^*),
Q int(lambda^*,lambda) = -int(lambda, dslash lambda).  L2 = log|p|^2 + gamma + log 2, L3 = log|p|^2 + gamma + log 3.
"""
import json
import sympy as sp

P = sp.symbols('p0:4', real=True)
p = sp.Matrix(P)
pp = sum(x * x for x in P)
c2, c3 = sp.symbols('c2 c3')          # gamma + log 2, gamma + log 3 (kept symbolic: must drop out)
L2 = sp.log(pp) + c2
L3 = sp.log(pp) + c3
I4 = sp.eye(4)
E = [sp.Matrix([1 if i == m else 0 for i in range(4)]) for m in range(4)]
R = sp.Rational
out = {}


def subs_neg(expr):
    return expr.subs({x: -x for x in P}, simultaneous=True)


def div(beta):
    """beta: list of 4 matrices (functions of p); returns sum_mu d/dp_mu beta^mu"""
    return sum((beta[m].diff(P[m]) for m in range(4)), sp.zeros(*beta[0].shape))


def Q_second_leg(alpha, beta, h):
    """h: function q -> matrix (Q-image of the second leg at momentum q)"""
    q = sp.symbols('q0:4', real=True)
    hq = h(sp.Matrix(q))
    at = {q[i]: -P[i] for i in range(4)}
    h0 = hq.subs(at, simultaneous=True)
    dh = [hq.diff(q[m]).subs(at, simultaneous=True) for m in range(4)]
    aQ = alpha * h0 - sum((beta[m] * dh[m] for m in range(4)), sp.zeros(*alpha.shape))
    bQ = [beta[m] * h0 for m in range(4)]
    return aQ, bQ


def bos_sym(alpha, beta):
    db = div(beta)
    a_sw = subs_neg(alpha).T - subs_neg(db).T
    b_sw = [subs_neg(beta[m]).T for m in range(4)]
    return (alpha + a_sw) / 2, [(beta[m] + b_sw[m]) / 2 for m in range(4)]


def simp(M):
    return M.applyfunc(lambda e: sp.simplify(sp.expand(e)))


# ---------------- A|_Wperp sector (flavour pairing (a,b) -> scalar 1x1 structure) ----------------
def run_AP(alpha_c, beta_c):
    alpha = sp.Matrix([[alpha_c]]); beta = [sp.Matrix([[beta_c * P[m]]]) for m in range(4)]
    aQ, bQ = Q_second_leg(alpha, beta, lambda q: sp.Matrix([[-R(1, 2) * sum(x * x for x in q)]]))
    aS, bS = bos_sym(aQ, bQ)
    return sp.simplify(aS[0] / pp), [sp.simplify(b[0]) for b in bS]

ext_AP = run_AP(R(43, 2) - 8 * L2, R(43, 2) - 8 * L2)
pi_AP = run_AP(R(8, 3), R(10, 3))
tot_AP = run_AP(R(145, 6) - 8 * L2, R(149, 6) - 8 * L2)
out['A_Wperp'] = {
    'ext_leg Q_sym / int(dA,dA)': str(ext_AP[0]), 'Pi_int Q_sym / int(dA,dA)': str(pi_AP[0]),
    'total Q_sym / int(dA,dA)': str(tot_AP[0]), 'beta_sym': [str(x) for x in tot_AP[1]],
    'engine_B_exact': '-11/3', 'agree': bool(sp.simplify(tot_AP[0] + R(11, 3)) == 0 and all(x == 0 for x in tot_AP[1]))}

# ---------------- A|_V sector ----------------
Lo = p * p.T / pp
def RV_matrices(with_F0):
    # engine A pieces (task1 'A^*_V,A_V'), summed here (not the hand-typed total)
    gV_ab = R(43, 4) - 4 * L2
    gV_Lo = 4 * L2 - 4 + (-R(27, 4) if with_F0 else 0)        # -F_0 = -(27/16) int(d^*A_V,d^*A_V) in W_1^ren
    al = 2 * (gV_ab * I4 + gV_Lo * Lo)                        # ext leg on W_1^AV: alpha = 2 g_V
    be = [2 * P[m] * (gV_ab * I4 + gV_Lo * Lo) for m in range(4)]   # beta = 2 g_V p
    al += (2 + 2 * L2) * Lo                                   # ext leg on the ghost self-energy
    be = [be[m] + (2 * L2 - 2) * P[m] * Lo for m in range(4)]
    # Pi_int (log eps part removed: L3e -> L3)
    al += (R(13, 6) + L3) * I4 + 2 * Lo
    be = [be[m] + R(10, 3) * P[m] * I4 + (L3 - R(1, 2)) * E[m] * p.T - 2 * P[m] * Lo for m in range(4)]
    return al, be

def total_alpha_ab_coefficient():
    # (a.b) coefficient of alpha of the sum of the pieces, without logs: 2*43/4 + 13/6
    return R(43, 2) + R(13, 6)

hV = lambda q: -R(1, 2) * (sum(x * x for x in q) * sp.eye(4) - q * q.T)
resV = {}
for wf0 in (False, True):
    al, be = RV_matrices(wf0)
    aQ, bQ = Q_second_leg(al, be, hV)
    aS, bS = bos_sym(aQ, bQ)
    aS = simp(aS)
    # decompose on x1 (|p|^2 I - p p^T) + x2 p p^T
    x1 = sp.simplify(aS[0, 1] / (-P[0] * P[1])) if True else None
    x1 = sp.simplify(-aS[0, 1] / (P[0] * P[1]))
    x2 = sp.simplify(aS[0, 1] / (P[0] * P[1]) + x1)
    resid = simp(aS - (x1 * (pp * sp.eye(4) - p * p.T) + x2 * p * p.T))
    bres = [simp(b) for b in bS]
    resV['with_F0' if wf0 else 'without_F0'] = {
        'coefficient int(dA_V,dA_V)': str(x1), "coefficient int(d^*A_V,d^*A_V)": str(x2),
        'alpha_residual_zero': all(sp.simplify(e) == 0 for e in resid),
        'beta_sym_zero': all(all(sp.simplify(e) == 0 for e in b) for b in bres)}
out['A_V'] = resV
out['A_V']['engine_B_exact'] = '-11/3 int(dA_V,dA_V) + 0 int(d^*A_V,d^*A_V)'
out['A_V']['engine_A_hand_total_alpha_(a.b)_constant'] = {'reported': '70/3', 'sum_of_reported_pieces': str(total_alpha_ab_coefficient())}

# ---------------- (c, A^*_V) sector ----------------
# R_V with field leg A_V(p) -> Q-image -dc: a -> -i p gamma (scalar gamma); R_c with c^*(p2) -> Q-image i (p2.b)
g = sp.Symbol('gamma_pol'); bvec = sp.Matrix(sp.symbols('b0:4'))
def c_sector(with_F0):
    al, be = RV_matrices(with_F0)
    a_vec = -sp.I * p * g
    alV = (a_vec.T * al * bvec)[0]
    beV = [(a_vec.T * be[m] * bvec)[0] for m in range(4)]
    # R_c total (task1 'c^*,c'): alpha = -8/3, beta = p (7/6 - 2 L2 - L3), structure (gamma, w)
    al_c = -R(8, 3); be_c = [P[m] * (R(7, 6) - 2 * L2 - L3) for m in range(4)]
    # w -> h(p2) = i (p2.b): alpha_Q = al_c h(-p) - be_c . dh, beta_Q = be_c h(-p)
    h0 = -sp.I * (p.T * bvec)[0]; dh = [sp.I * bvec[m] for m in range(4)]
    alc = g * (al_c * h0 - sum(be_c[m] * dh[m] for m in range(4)))
    bec = [g * be_c[m] * h0 for m in range(4)]
    a_tot = sp.simplify(sp.expand(alV + alc)); b_tot = [sp.simplify(sp.expand(beV[m] + bec[m])) for m in range(4)]
    coeff = sp.simplify(a_tot / (sp.I * g * (p.T * bvec)[0]))
    return coeff, b_tot
cs = {}
for wf0 in (False, True):
    coeff, bt = c_sector(wf0)
    cs['with_F0' if wf0 else 'without_F0'] = {'coefficient int(dc,A^*)': str(coeff), 'beta': [str(x) for x in bt]}
out['c_AsV'] = cs
out['c_AsV']['engine_B_exact'] = '-14/3'
out['c_AsV']['note'] = ('int(dc,A^*) <-> i(gamma,p.b); log cancellation between the A_V and c contributions fixes their '
                        'relative sign; F_0 = (27/16) unit int(d^*A_V,d^*A_V) (Step 1) must be in W_1^ren')

# ---------------- lambda sector ----------------
# kernel (u, X v), u = field (Sigma), v = antifield (Sigma^*); structures alpha*(u,v), beta = f p^mu (u,v) + g (u, rho(e^mu) rho(p) v)
# Q: v -> kappa i rho(p2) u2 (kappa = -1: Q int(lambda^*,lambda) = -int(lambda,dslash lambda));
# rules: rho(a)rho(b)+rho(b)rho(a) = 2(a,b) (engine A's log-eps identity), (u1,rho(x)u2) = (u2,rho(x)u1);
# fermionic swap K_sym = (K - K')/2 gives beta_sym = 0 and
#   alpha_sym = -(i kappa/2)[2 alpha - 3 f - 6 g - 2|p|^2 (f' + g')] (u1, rho(p) u2),   int(lambda,dslash lambda) <-> -i(u1,rho(p)u2)
# so the coefficient of int(lambda,dslash lambda) is  (kappa/2)[2 alpha - 3 f - 6 g - 2|p|^2 (f'+g')].
s = sp.Symbol('s', positive=True)     # s = |p|^2
L2s = sp.log(s) + c2; L3s = sp.log(s) + c3
def lam_coeff(al, f, gg, kappa=-1):
    expr = 2 * al - 3 * f - 6 * gg - 2 * s * (sp.diff(f, s) + sp.diff(gg, s))
    return sp.simplify(sp.expand(sp.Rational(kappa, 2) * expr))
lam = {
    'ext_leg (D_M Sigma_lambda)': str(lam_coeff(6 * (4 * L2s - 4), 4 * (4 * L2s - 4), 0)),
    'Pi_int': str(lam_coeff(-R(8, 3), -R(1, 3) - 2 * L3s, L3s - R(1, 2))),
    'total': str(lam_coeff(24 * L2s - R(80, 3), 16 * L2s - 2 * L3s - R(49, 3), L3s - R(1, 2))),
    'engine_B_exact': '47/3',
    'expected ext-leg value 4 x (log eps coefficient of I^CT_2 in the sector) = 4 x 4': '16'}
out['lambda'] = lam

# ---------------- F(Delta): check Q F = Ob and the field-counting expectation ----------------
# Q images (engine B Q_basis): Q int(A,A^*) = -1/2 dAdA - dcA*, Q int(lambda^*,lambda) = -lam dslash lam, Q int(c^*,c) = -dcA*
def QF(xA, xl, xc):
    return {'dAdA': -R(1, 2) * xA, 'lamdlam': -xl, 'dcAs': -xA - xc}
FB = QF(R(22, 3), -R(47, 3), -R(8, 3))
FC = QF(2 * R(11, 3), -4 * R(11, 3), -R(11, 3))
out['F_Delta'] = {
    'Q F_B (engine B translation-invariant F)': {k: str(v) for k, v in FB.items()},
    'Ob(Delta) engine B': {'dAdA': '-11/3', 'lamdlam': '47/3', 'dcAs': '-14/3'},
    'Q((11/3) field-counting)': {k: str(v) for k, v in FC.items()},
    'translation-invariant Q-closed directions x1 int(A,A^*) + x2 int(lambda^*,lambda) + x3 int(c^*,c)':
        'Q = -x1/2 dAdA - x2 lamdlam - (x1+x3) dcAs = 0 only for x1 = x2 = x3 = 0: the translation-invariant F is unique',
    'F proportional to the field-counting functional': False}

json.dump(out, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/compare_delta/ob_reassembly.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
