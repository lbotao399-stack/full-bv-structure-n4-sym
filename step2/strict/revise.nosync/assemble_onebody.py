"""Reviser: one-body parts -{hbar F_1^(1)(x), O} on F0 from the action-side solutions, and the assembled components.
Inputs (read, not recomputed): action/complete/results.json, operator engine results (two-body F0 entries as recorded
in operator/engine_a/results.json and engine_b/results.json; checked here by string search), action engines.
Units u; L = log(3/2). Output: revise.nosync/onebody_F0.json"""
import json, os, sympy as sp
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
L = sp.log(sp.Rational(3, 2))
tD, tP, tM, tR, txp, txm, tep, tem, tK = sp.symbols('t_Delta t_P t_M t_R t_xi+ t_xi- t_eta+ t_eta- t_K')
a, c = sp.symbols('a c')   # outer, inner directions
comp = json.load(open(os.path.join(ROOT, 'action/complete/results.json')))
# --- action-side particular solutions on the A_perp arrow (two legs) ---
p_xi = L - sp.Rational(79, 12)          # coefficient of T1+T2 in F1(xi)|_2, both chiralities, both routes
fp = comp['F1(eta_+)|2']['F1_particular_B']
assert sp.nsimplify(sp.sympify(fp['x.T1(rho(e)eta)'].replace('log', 'sp.log'), locals={'sp': sp}) - p_xi) == 0
# eta: coefficient of x.T2(rho eta) in the primed parametrisation = p_xi + t_eta (uniform)
# --- two-body Q_1 on F0 (both engines, both chiralities): a1..a4 ---
a1, a2, a3, a4 = -2*L, 4*L, 4*L, 0
for f in ['operator/engine_a/notes.txt', 'operator/engine_b/notes.txt']:
    s = open(os.path.join(ROOT, f)).read()
    assert '-2L F1 + 4L F2 + 4L F3' in s, f
# --- one-body on O = int omega (iota_w1 A, iota_w2 A): only the A*_perp terms of F1 act: T2 in F1(xi) ---
# {T2, O} = {S^(1)(xi), O} = F1 (the A_perp arrow is T2 alone); hence -{F1(xi),O} = -(p_xi + t_xi) F1
cQ = {'+': sp.simplify(a1 - (p_xi + txp)), '-': sp.simplify(a1 - (p_xi + txm))}
cS = {'+': sp.simplify(a1 - (p_xi + tep)), '-': sp.simplify(a1 - (p_xi + tem))}
# dilatation on F0 (engines): sum Upsilon^0_Delta = 40/3 (traceless), -104/3 (trace); one-body -{F1(Delta),O} = -(44/3 + 4 t_Delta) O
D1_tl = sp.Rational(40, 3) - sp.Rational(44, 3) - 4*tD
D1_tr = -sp.Rational(104, 3) - sp.Rational(44, 3) - 4*tD
p1 = -(sp.Rational(1, 3) + tP)          # P_1 = p1 u d_b O
# action-side computed relations (two legs, two routes)
rel = {tD: -sp.Rational(1, 3)}
sumxi = tP + sp.Rational(27, 4) + 2*L    # t_xi+ + t_xi-
sumxe = tD + sp.Rational(27, 4) + 2*L    # t_xi+ + t_eta+
inv1 = sp.simplify(cQ['+'] + cQ['-'] - p1)          # depends on t_xi+ + t_xi- - t_P only
inv1v = sp.simplify(inv1.subs(txm, sumxi - txp))
inv2 = sp.simplify(cQ['+'] + cS['+'] - D1_tl/4)     # depends on t_xi+ + t_eta+ - t_Delta only
inv2v = sp.simplify(inv2.subs(tep, sumxe - txp))
out = {
 'units': 'u = -g^2 N/(32 pi^2); L = log(3/2)',
 'Q1_F0_traceless': {'form': 'Q_1(O) = c_Q {S^(1)(xi),O}  (F2,F3,F4 vanish on traceless: (w1,w2) = 0)', 'c_Q(xi+)': str(cQ['+']), 'c_Q(xi-)': str(cQ['-'])},
 'Q1_F0_general': {'two_body': 'a1 F1 + a2 F2 + a3 F3 + a4 F4', 'a': [str(a1), str(a2), str(a3), str(a4)],
                    'one_body': '-(L - 79/12 + t_xi) F1', 'total_F1_coefficient': str(cQ['+']).replace('t_xi+', 't_xi'),
                    'arity_zero': 'L_1^(0)(O) = -(2 + 4L) u (w1,w2) int omega (dc, A_V) at psi = 0'},
 'S1_F0_traceless': {'form': 'S_1(O) = c_S {S^(1)(eta),O}; position part x_O^mu Q_1 along rho(e_mu) eta; intrinsic part 0 on F0',
                     'c_S(eta+)': str(cS['+']), 'c_S(eta-)': str(cS['-'])},
 'D1_F0': {'traceless': str(D1_tl), 'trace': str(D1_tr), 'at_t_Delta=-1/3': [str(D1_tl.subs(rel)), str(D1_tr.subs(rel))]},
 'P1_F0': str(p1) + ' u d_b O',
 'K1_F0_onebody': ('-{F1(K_b)|_2, O} = (44/3) u int omega (b,x)(iota_w1 A, iota_w2 A) - u theta^K_perp {S(K_b)_perp, O} '
                   '- sum_j kappa^K_j {G^K_j, O}; x-linear part = -(b,x_O) x (one-body of D^0_1 from F1part(Delta)) (translation covariance)'),
 'free_direction_invariants_at_psi0': {
     'c_Q(xi+) + c_Q(xi-) - p_1': str(inv1), 'value': str(inv1v),
     'c_Q(xi+) + c_S(eta+) - D_1(traceless)/4': str(inv2), 'value_2': str(inv2v)},
 'checks': {'inv1 independent of t_xi, t_P': inv1v.free_symbols == set(), 'inv2 independent of t_xi, t_eta, t_Delta': inv2v.free_symbols == set()},
}
# outer/inner invariance of the relations used
shift = {txp: txp + a - c, txm: txm - a - c, tep: tep - a + c, tem: tem + a + c, tP: tP - 2*c, tD: tD}
out['checks']['t_xi+ + t_xi- - t_P invariant'] = sp.simplify((txp + txm - tP).subs(shift, simultaneous=True) - (txp + txm - tP)) == 0
out['checks']['t_xi+ + t_eta+ - t_Delta invariant'] = sp.simplify((txp + tep - tD).subs(shift, simultaneous=True) - (txp + tep - tD)) == 0
json.dump(out, open(os.path.join(ROOT, 'revise.nosync/onebody_F0.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
