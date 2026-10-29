#!/usr/bin/env python3
"""cubic_assemble.py -- builds cubic_maps.json from cubic_maps_numeric.json (exact plane-wave results) and
cubic_census.json; solves Ob_3 = Q F_3 + {I_3, F_2} symbolically (sympy)."""
import json
import sympy as sp

D = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory/'
num = json.load(open(D + 'cubic_maps_numeric.json'))
cen = json.load(open(D + 'cubic_census.json'))
R = num['results']
SEC_OF = {'G1': 'c_AV_AV', 'G2': 'c_AV_AV', 'G3': 'c_AW_AW', 'G4': 'c_c_AsV', 'G5': 'c_lam_lam'}
GB = ['G1', 'G2', 'G3', 'G4', 'G5']


def vec(entry):
    out = []
    for g in GB:
        out.append(sp.Rational(entry[SEC_OF[g]][g]))
        assert entry[SEC_OF[g]]['_residual_zero']
    return out


QH = sp.Matrix([vec(R[f'QH{j}']) for j in range(1, 7)]).T          # 5 x 6, columns = Q H_j
BR = {X: sp.Matrix(vec(R['{I3,' + X + '}'])) for X in ('D1', 'B1', 'B2', 'B3', 'B4', 'mass', 'massH')}
TREE = sp.Matrix(vec(num['trees']['tree_Ob2int_I3']))

a = sp.symbols('alpha1:6')
b = sp.symbols('beta1:5')
h = sp.symbols('h1:7')
F2br = sp.Rational(27, 16) * BR['D1'] + b[0] * BR['B1'] + b[1] * BR['B2'] + b[2] * BR['B3'] + b[3] * BR['B4']
eqs = list(QH * sp.Matrix(h) + F2br - sp.Matrix(a))
sol = sp.solve(eqs[:1] + eqs[2:], [h[0], h[2], h[3], h[4]], dict=True)[0]   # G2 row has no h
g2_condition = sp.simplify(eqs[1])
kernel = QH.nullspace()
rankQ = QH.rank()
imageQ_plus_beta = (QH.row_join(BR['B1']).row_join(BR['B2']).row_join(BR['B3']).row_join(BR['B4'])).rank()

# s-closure directions: Q X3 + {I3, X2} = 0
dirs = {
    'beta1: int(dA|_V,dA|_V) completed by H1 [cubic part of int(F_{A|_V},F_{A|_V})]': (sp.Matrix([1, 0, 0, 0, 0, 0]), 'B1'),
    'beta2: int(dA|_Wperp,dA|_Wperp) completed by 2 H2 [cubic part of int(d_A A|_Wperp, d_A A|_Wperp)]':
        (sp.Matrix([0, 2, 0, 0, 0, 0]), 'B2'),
    'beta3: int(lambda,dslash lambda) completed by -2 H5 (in span H; the cubic part int(lambda, rho([A,lambda])) of '
    'int(lambda,dslash_A lambda) has the same Q-image -2 G5)': (sp.Matrix([0, 0, 0, 0, -2, 0]), 'B3'),
    'beta4: int(dc,A^*) completed by H3+H4+H5+(1/2)H6 [cubic part of the BRST (antifield) term of S]':
        (sp.Matrix([0, 0, 1, 1, 1, sp.Rational(1, 2)]), 'B4'),
}
dir_check = {k: str(list(QH * v + BR[X])) for k, (v, X) in dirs.items()}

particular = {str(k): str(v) for k, v in sol.items()}
family = ('F_3/unit = (alpha4/4 - alpha1/2) H1 - alpha4 H3 - 2 alpha3 H4 + alpha5 H5 + beta1 H1 + 2 beta2 H2 - 2 beta3 H5 '
          '+ beta4 (H3 + H4 + H5 + H6/2) + mu1 (H1 - 4 H3 - 2 H6) + mu2 (H2 - 2 H4),  mu1, mu2 free')
# verify the family
mu1, mu2 = sp.symbols('mu1 mu2')
Fvec = (sp.Matrix([a[3] / 4 - a[0] / 2, 0, -a[3], -2 * a[2], a[4], 0]) + b[0] * dirs[list(dirs)[0]][0]
        + b[1] * dirs[list(dirs)[1]][0] + b[2] * dirs[list(dirs)[2]][0] + b[3] * dirs[list(dirs)[3]][0]
        + mu1 * sp.Matrix([1, 0, -4, 0, 0, -2]) + mu2 * sp.Matrix([0, 1, 0, -2, 0, 0]))
resid = sp.simplify(QH * Fvec + F2br - sp.Matrix(a))
family_check = [str(x) for x in resid]

kappa = sp.Rational(-27, 8)
out = {
    'unit': 'hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y))',
    'conventions': {
        'Q_on_functionals': 'Q F = left d/dr F(e + r Qe); Qe: c -> -dc (A|_V), A -> -1/2 d^*dA (A^*), A^*|_V -> d^*A^*|_V '
                            '(c^*), lambda -> dslash lambda (lambda^*)  [Prop c1-Q]',
        'bracket': '{f,g} = Delta_0(fg) - Delta_0(f)g - (-1)^{|f|} f Delta_0(g), Delta_0 = sum d_x d_xi, x in {A, lambda^*, '
                   'c^*}  [Def c1-bracket]; checks: {S_2,-} = Q, field-count lemma, quartic CME',
        'G_basis': {'G1': 'int([c,A|_V], d^*d A|_V)', 'G2': 'int([c,A|_V], d d^* A|_V)',
                    'G3': 'int([c,A|_Wperp], d^*d A|_Wperp)', 'G4': 'int([c,dc], A^*)',
                    'G5': 'int([c,lambda], dslash lambda)'},
        'H_basis': {'H1': 'int(dA|_V,[A|_V ^ A|_V])', 'H2': 'int(dA|_Wperp,[A|_V,A|_Wperp])',
                    'H3': 'int([A|_V,c],A^*|_V)', 'H4': 'int([A|_Wperp,c],A^*|_Wperp)',
                    'H5': 'int([c,lambda],lambda^*)', 'H6': 'int([c,c],c^*)'},
        'quadratic': {'D1': 'int(d^*A|_V,d^*A|_V)', 'B1': 'int(dA|_V,dA|_V)', 'B2': 'int(dA|_Wperp,dA|_Wperp)',
                      'B3': 'int(lambda,dslash lambda)', 'B4': 'int(dc,A^*)', 'mass': 'int(A,A)',
                      'massH': 'int(A, H A), H = d^*d + dd^* componentwise',
                      'Ob2int': 'int(d^*dc, d^*A|_V)  (Ob_2 = -27/8 unit Ob2int)'},
        'plane_wave_values': {'G1': 'f r0 [(p2^2-p1^2)(a1.a2) - (a1.p2)(a2.p2) + (a1.p1)(a2.p1)]',
                              'G2': 'f r0 [(a1.p2)(a2.p2) - (a1.p1)(a2.p1)]', 'G3': 'f r0 (w1.w2)(p2^2 - p1^2)',
                              'G4': '-i f (p2.b) r0 r1 r2', 'G5': '-i f (lambda1, rho(p0) lambda2) r0 r1 r2',
                              'f': '([X0,X1],X2)'},
    },
    'sector_dimensions': {
        '(c,A|_V,A|_V)': {'weight4_even': 2, 'basis': ['G1', 'G2'], 'parity_odd': 0,
                          'reason': 'K antisymmetric in legs 1,2: (a1.a2)(p1^2-p2^2), (a1.p1)(a2.p1)-(a1.p2)(a2.p2); '
                                    'eps(a1,a2,p1,p2) is symmetric'},
        '(c,A|_Wperp,A|_Wperp)': {'weight4_even': 1, 'basis': ['G3'], 'parity_odd': 0},
        '(c,lambda,lambda)': {'weight4_even': 1, 'basis': ['G5'], 'parity_odd_candidates': 1,
                              'odd_candidate': 'int([c,lambda], dslash rho(v1 v2 v3 v4) lambda) = '
                                               'i f (lambda1, rho(p2-p1) rho(vol) lambda2): excluded by the reflection '
                                               '(one direction of V) x (one direction of W_R^perp) in Spin(10)'},
        '(c,c,A^*|_V)': {'weight4_even': 1, 'basis': ['G4'], 'parity_odd': 0},
        'weight2_structures': 'int([c,A|_V],A|_V), int([c,A|_Wperp],A|_Wperp): identically zero (verified numerically)',
        'other_degree1_three_leg_sets': '(c,A|_V,A|_Wperp), (c,A,lambda), (c,c,lambda^*), (c,c,A^*|_Wperp): no invariant',
    },
    'matrices': {
        'Q_H_in_G_basis(rows G1..G5, cols H1..H6)': [[str(x) for x in QH.row(i)] for i in range(5)],
        'bracket_I3_in_G_basis': {X: [str(x) for x in v] for X, v in BR.items()},
        'Q_G_degree2_sector(c,c,A|_V)': {'QG1': '0 (identically: [c,dc] = d[c,c]/2, d^*d^* = 0)', 'QG2':
                                         'int([c,dc], dd^*A|_V) + int([c,A|_V], dd^*dc) != 0', 'QG3': '0', 'QG4': '0',
                                         'QG5': '0'},
        '{I3,Ob2int}': 'Q G2',
        '{I3,Ob_2}': '-(27/8) unit Q G2',
        'tree(Ob2int, I3; P(0,infinity)) in G basis': [str(x) for x in TREE],
        'tree(Ob_2, I3; P(0,infinity)) in units': [str(kappa * x) for x in TREE],
    },
    'rank': {'rank Q on span H': int(rankQ), 'image Q span H + beta brackets': int(imageQ_plus_beta),
             'missing direction': 'G2 (reached only by {I3, int(d^*A|_V,d^*A|_V)} = 2 G2)',
             'kernel of Q on span H': ['H1 - 4 H3 - 2 H6', 'H2 - 2 H4'],
             'kernel_check': [str(list(v)) for v in kernel]},
    'solution': {
        'equation': 'Ob_3 = Q F_3 + {I_3, F_2},  Ob_3 = unit sum alpha_i G_i,  F_3 = unit sum h_j H_j',
        'G2_component_condition': f'alpha2 = 27/8  (from {g2_condition} = 0)',
        'particular_solution(h2=h6=0)': particular,
        'family': family,
        'family_residual': family_check,
        's_closure_directions(Q X3 + {I3,X2})': dir_check,
    },
    'wess_zumino': {
        'statement': 'Q Ob_3 + {I_3, Ob_2} = (alpha2 - 27/8) unit Q G2, Q G2 != 0: WZ <=> alpha2 = 27/8',
        'alpha1_alpha4': 'unconstrained (Q G1 = Q G4 = 0 identically; Q G3 = Q G5 = 0)',
        'automatic': 'yes, as a consequence of Q^2 = 0 on Q W(P(eps,infty),J), Q I_3 = 0, Q Ob_2 = 0 and the two limit '
                     'statements (I-b + III_{eps-0} -> 0, I-d + II-b -> tree(Ob_2,I_3;P(0,inf))); not from linear '
                     'algebra: Q + beta-brackets is not surjective onto span G (G2 missing). Solvability <=> WZ <=> '
                     'alpha2 = 27/8.',
    },
    'log_and_inverse_eps_bookkeeping(units of unit)': {
        'III_eps^{-1}': '0: {I_3, int(A,A)} = 0 identically',
        'II-a_log_per_row': {'(c,A|_V,A|_V)': {'(A|_V)^3 x6': '-3 G1', 'A|_V c A^*|_V x1': '+1/2 G1', 'total': '-5/2 G1'},
                             '(c,A|_Wperp,A|_Wperp)': {'A|_V (A|_Wperp)^2 x6': '-3 G3', 'A|_Wperp c A^*|_Wperp x1': '+1/2 G3',
                                                       'total': '-5/2 G3'},
                             '(c,lambda,lambda)': {'lambda lambda A|_V x10': '+10 G5', 'c lambda lambda^* x1': '-G5',
                                                   'total': '+9 G5'},
                             '(c,c,A^*|_V)': {'A|_V c A^*|_V x1': '+G4', 'c c c^* x1': '-G4', 'total': '0'}},
        'III_K0_log_per_row': {'(c,A|_V,A|_V)': {'(A|_V)^3 with -int(dc,A^*)': '+1/2 G1',
                                                 'A|_V c A^*|_V with -int(dA|_V,dA|_V)': '+2 G1', 'total': '+5/2 G1'},
                               '(c,A|_Wperp,A|_Wperp)': {'A|_V (A|_Wperp)^2 with -int(dc,A^*)': '+1/2 G3',
                                                         'A|_Wperp c A^*|_Wperp with -int(dA|_Wperp,dA|_Wperp)': '+2 G3',
                                                         'total': '+5/2 G3'},
                               '(c,lambda,lambda)': {'lambda lambda A|_V with -int(dc,A^*)': '-G5',
                                                     'c lambda lambda^* with 4 int(lambda,dslash lambda)': '-8 G5',
                                                     'total': '-9 G5'},
                               '(c,c,A^*|_V)': {'A|_V c A^*|_V with -int(dc,A^*)': '-G4', 'c c c^* with -int(dc,A^*)': '+G4',
                                                'total': '0'}},
        'II-a + III_K0 (log)': '0 in every sector (= -(s s Y)_3 with Y = 2int(A,A^*) - 4int(lambda^*,lambda) - int(c^*,c))',
        'I-a': 'eps^{-1}: weight-2 forms, zero after the sum over identical legs; log eps: 0 graph by graph (power counting)',
        'lim III_{K_eps-K_0}': {'(c,A|_V,A|_V)': '+9/2 (G1+G2)', '(c,A|_Wperp,A|_Wperp)': '+9/2 G3',
                                '(c,lambda,lambda)': '0', '(c,c,A^*|_V)': '0'},
        'lim I-b': {'(c,A|_V,A|_V)': '-9/2 (G1+G2)', '(c,A|_Wperp,A|_Wperp)': '-9/2 G3', '(c,lambda,lambda)': '0',
                    '(c,c,A^*|_V)': '0'},
        'lim (I-d + II-b)': {'(c,A|_V,A|_V)': '27/8 G1 - 27/8 G2', '(c,A|_Wperp,A|_Wperp)': '27/8 G3',
                             '(c,lambda,lambda)': '-27/4 G5', '(c,c,A^*|_V)': '-27/4 G4'},
        'total three-leg limit of Q W^(1)': 'Ob_3 + lim(I-d + II-b)',
    },
    'checks': num['checks'],
    'wz_numeric': num['wz_degree2'],
    'tree_identity_check': num['trees']['check_Q_tree_identity'],
    'census': cen,
}
json.dump(out, open(D + 'cubic_maps.json', 'w'), indent=1)
print('QH =', QH)
print('brackets', {k: list(v) for k, v in BR.items()})
print('G2 condition', g2_condition, 'particular', particular)
print('family residual', family_check)
print('dir checks', dir_check)
print('kernel', [list(v) for v in kernel], 'rank', rankQ, imageQ_plus_beta)
print('tree', list(TREE))
