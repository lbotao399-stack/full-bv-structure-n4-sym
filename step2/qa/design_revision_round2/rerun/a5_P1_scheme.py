"""Audit script a5: what P_1 = 0 costs, given the composite-line finding of engine B
(step2/engine_b/delta_results.json, part2.coefficient.normalisation_of_the_classical_part and
part2.translation_part_check: in the translation-covariant scheme the translation component of the scalar
bilinears is (1 - unit/3) D_P O, i.e. P_1(O) = -(1/3) unit {S^(1)(P_b),O}; letters: 0 per
step2/theory/delta_predictions.json translation_component).

Scheme psi' = sum_n c_n (E + w) P_n = sum_n (c_n/2) {S^(1)(Delta),-} P_n, P_n the projector on n-letter
local functionals (well defined on functionals int omega f, unlike 'c_n x Weyl weight of the density', since
int (d omega) f = - int omega (d f) relates densities of different weight).  Def c8-renorm allows it
(degree and weight preserving, C^bullet-linear).  By Prop c8-laws (b) the component along x changes by
[s,psi']^(1)(x):
  even x (preserves the letter number):  (c_n/2) e_x delta_x O,      e_x = ad_Delta eigenvalue
  odd x:  (c_n/2) e_x delta_x O + ((c_n - c_{n+1})/2) {S^(1)(Delta), (delta_x O)_{n+1}}
using [delta_x, delta_Delta] = e_x delta_x (from {S^(1)(x),S^(1)(Delta)} = S^(1)([x,Delta]), a2 toy signs).
"""
import json, os
import sympy as sp

c1, c2, c3, u = sp.symbols('c1 c2 c3 unit')
e = {'P_b': -2, 'K_b': 2, 'Delta': 0, 'so(V)': 0, 'so(W)': 0, 'Sigma (Q)': -1, 'Sigma^* (S)': 1}
res = {}
# P_1 = 0 on letters (n=1) and scalar bilinears (n=2): -c_n = -(P_1 coefficient)
P1 = {1: 0, 2: sp.Rational(-1, 3) * u}
sol = {c1: -P1[1] * 0, c2: None}
# shift along P on n-letter operators: (c_n/2)(-2) = -c_n; need P_1 + shift = 0 -> c_n = P1_n
cs = {1: P1[1], 2: P1[2]}
res['required_c_n'] = {str(k): str(v) for k, v in cs.items()}
tab = {}
for n in (1, 2):
    row = {}
    for x, ex in e.items():
        row[x] = str(sp.simplify(cs[n] * ex / 2)) + ' x delta_x O'
    tab[f'n={n}'] = row
res['component_shifts_in_the_P1=0_scheme_relative_to_covariant_scheme'] = tab
res['odd_x_extra_letter_raising_term'] = '((c_n - c_{n+1})/2) {S^(1)(Delta), (delta_x O)_{n+1}}: zero for the scalar bilinears (delta_xi of A|_Wperp is linear in lambda), nonzero for O_zeta (delta_xi O_zeta = kappa_B B_w, three letters): shifts c_B\', not c_K'
res['consequence_for_thm_o_KD_targets'] = {
    'P_1 = 0 scheme': 'K_1 = (1/4) div v gamma_O O (thm:o-KD, if its other hypotheses hold)',
    'translation-covariant scheme (P_1 = -(1/3) unit delta_P on bilinears)': 'K_1 = (1/4) div v gamma_O O + (1/3) unit {S^(1)(K_b),O} modulo s-exact (shift -c_2 delta_K with c_2 = -(1/3) unit)',
    'Q_1, S_1 on bilinears': 'differ between the two schemes by -(c_2/2) delta_xi = +(1/6) unit delta_xi and +(c_2/2) delta_eta = -(1/6) unit delta_eta',
    'D_1, M_1, R_1': 'unchanged (e_x = 0)'}
res['reflection_compatibility'] = 'psi\' commutes with Theta_n for real c_n (E and w are reflection invariant; {S^(1)(Delta),a}^Theta = {S^(1)(Delta),a^Theta})'
res['pair_kernel_conjecture'] = ('if the composite remainder acts on pairs of letters as -(1/3) unit (d_i + d_j), then on n-letter '
                                 'operators R = -((n-1)/3) unit delta_P, c_n = -(n-1)/3 unit: reproduces c_1 = 0, c_2 = -1/3; '
                                 'conjecture, to be tested by E11 on a three-letter operator')
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a5_P1_scheme.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
