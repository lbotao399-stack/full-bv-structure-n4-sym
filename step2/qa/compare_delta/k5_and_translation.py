"""QA compare_delta (Part D4, bilinear component), two checks.

(A) Engine B's K5 fit (alpha = 16 T1 - 16 T2, beta = (2/3) k T1) re-evaluated on its raw exact samples
    (delta_results.json part2.K5.raw), and engine A's pieces mapped onto K5:
      K5 = 4 x (log-eps coefficient of d_delta I^CT) + {Pi_int|_P, O} + Q(Y1a),
    with Pi_int|_P = (5/3) S_M|_P - (2/3) int(A_P,A^*_P) (both engines), Q(Y1a) = -2 D_M O (y1a_symmetrize.py),
    in plane-wave (vertex-value) normalisation: O <-> (T1, 0), D_M O <-> (-4 T1, -2 k T1).

(B) Independent analytic evaluation of the translation closed form (engine B):
      [(d_CE+Q)O[infty]](translation) - (D_P O)[infty] = lim (hbar/2) d_{K''}(Q O_W)_tree,
    the one-loop graph with the composite A_V--c line (translation letter), the A_V(A_P)^2 vertex X and the vertex
    {I,O}_eps - {I,O}.  Ward identity at X (longitudinal composite line l_nu): l.(2 p1 + l) = (l+p1)^2 - |p1|^2.
    Calibration by Pi_int (the same loop with Y's A^* leg external): beta(Pi_int) = C pi^2 p [1 + r], r = regulator
    remnant of the cancelled A_P line.  The two K_eps assignments of {Y,O}_eps carry opposite colour sign
    (([X,c],Y) = -([Y,c],X)); assignment 1 is {Pi_int^T, O}, assignment 2 is the Q-image of Y1a_T.
    All eps -> 0 limits are Gaussian integrals after l = u/sqrt(eps):
      int d^4u u_mu u_nu u^-4 e^{-a u^2} = (delta_mu nu / 4) pi^2 / a.
"""
import json
import sympy as sp

R = sp.Rational
out = {}

# ---------------- (A) engine B raw K5 samples ----------------
B = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_b/delta_results.json'))
raw = B['part2']['K5']['raw']['ApAp']
chk = []
for smp in raw:
    w1 = [sp.Rational(x) for x in smp['w1']]; w2 = [sp.Rational(x) for x in smp['w2']]
    u1 = [sp.Rational(x) for x in smp['u1']]; u2 = [sp.Rational(x) for x in smp['u2']]
    k = [sp.Rational(x) for x in smp['k']]
    d = lambda a, b: sum(x * y for x, y in zip(a, b))
    T1 = d(w1, u1) * d(w2, u2) + d(w1, u2) * d(w2, u1)
    T2 = d(w1, w2) * d(u1, u2)
    a0 = sp.Rational(smp['res']['0']['alpha']['1'])
    bvals = [sp.Rational(v['1']) for v in smp['res']['0']['beta']]
    logs_zero = all(v == '0' for e in [smp['res']['0']['alpha']] + smp['res']['0']['beta'] for kk_, v in e.items() if kk_ != '1')
    sing_zero = all(v == '0' for e in [smp['res']['-1']['alpha']] + smp['res']['-1']['beta'] for v in e.values())
    chk.append({'T1': str(T1), 'T2': str(T2), 'alpha': str(a0), '16T1-16T2': str(16 * T1 - 16 * T2),
                'beta': [str(x) for x in bvals], '(2/3) k T1': [str(R(2, 3) * x * T1) for x in k],
                'alpha_ok': a0 == 16 * T1 - 16 * T2, 'beta_ok': bvals == [R(2, 3) * x * T1 for x in k],
                'eps^-1 and log parts zero': bool(logs_zero and sing_zero)})
out['engineB_K5_raw_vs_fit'] = chk

# engine A pieces -> K5 (traceless: log coefficient L = 4; per unit of T1)
kk = sp.Symbol('k')
O = (1, 0); DMO = (-4, -2)                        # (alpha/T1, beta/(k T1))
def comb(*terms):
    return tuple(sum(c * v[i] for c, v in terms) for i in range(2))
L = 4
Pi_bracket = comb((R(5, 3), DMO), (-R(2, 3) * 2, O))        # {S_M,O} = D_M O, {int(A,A^*),O} = 2 O
QY1a = comb((-2, DMO))
K5_A = comb((4 * L, O), (1, Pi_bracket), (1, QY1a))
out['K5_from_engineA_pieces_traceless'] = {'alpha/T1': str(K5_A[0]), 'beta/(k T1)': str(K5_A[1]),
                                          'engineB': {'alpha/T1': '16', 'beta/(k T1)': '2/3'},
                                          'agree': K5_A == (16, R(2, 3))}
QY1a_unsym = comb((-1, DMO))
K5_A_unsym = comb((4 * L, O), (1, Pi_bracket), (1, QY1a_unsym))
out['K5_with_engineA_half_normalised_Q(Y1a)=-D_M O'] = {'alpha/T1': str(K5_A_unsym[0]), 'beta/(k T1)': str(K5_A_unsym[1])}

# ---------------- (B) analytic translation closed form ----------------
u, a = sp.symbols('u a', positive=True)
J = lambda aa: sp.integrate(2 * sp.pi ** 2 * u ** 3 * u ** -2 * sp.exp(-aa * u ** 2), (u, 0, sp.oo))   # int d^4u u^-2 e^{-a u^2}
assert sp.simplify(J(3) - sp.pi ** 2 / 3) == 0
# vectors as coefficients on (p1, p2) [all results are linear combinations of the leg momenta]
# Pi_int (legs: field at p, A^* external): first Ward term int l D^2 e^{-eps(l+p)^2}, D(l) = e^{-eps l^2}/l^2:
#   eps^{-1/2} int u u^-4 e^{-2u^2} e^{-(u + sqrt(eps) p)^2}  ->  -2 (1/4) J(3) p
rem_Pi = -2 * R(1, 4) * J(3) / sp.pi ** 2               # coefficient of p (units pi^2)
main_Pi = 1                                               # -|p|^2 int l/(l^4 (l+p)^2) = +pi^2 p
out['Pi_int_remnant'] = {'remnant/p': str(rem_Pi), 'total/p': str(main_Pi + rem_Pi),
                         'meaning': 'beta(Pi_int) = (5/6) C p: with beta(Pi_int) = (10/3) p (both engines) C = 4'}
# assignment 2 (K_eps between O's second leg and Y's A^*, Y's A_P external at p2):
#   first Ward term: eps^{-1/2} int u u^-4 e^{-2u^2} e^{-(u+sqrt(eps)p1)^2} e^{-(u - sqrt(eps)p2)^2} -> (1/4) J(4) (-2 p1 + 2 p2)
#   second: -|p1|^2 int l/(l^4 (l+p1)^2) = + pi^2 p1
A2 = (R(1, 4) * J(4) / sp.pi ** 2 * (-2) + 1, R(1, 4) * J(4) / sp.pi ** 2 * 2)          # (coef p1, coef p2)
# assignment 1 = {Pi_int^T, O}: Pi_int value at external field momentum p1, opposite colour sign
A1 = (-(main_Pi + rem_Pi), 0)
G = (A2[0] + A1[0], A2[1] + A1[1])
# plane-wave value of a kernel c1 p1 + c2 p2 (alpha-type, bosonic legs): (c1 + c2)(p1 + p2) = -(c1 + c2) k
pw = lambda v: -(v[0] + v[1])                               # coefficient of k
ratio_G = sp.nsimplify(pw(G) / pw(A1))                     # G / {Pi_int^T,O}
ratio_A2 = sp.nsimplify(pw(A2) / pw(A1))
pi_T = R(5, 3)                                              # {Pi_int^T, O} = (5/3) D_P O  (Pi_int^T = (5/3) S^(1)(P)|_P)
out['translation_closed_form'] = {
    'assignment2_kernel (p1,p2) coefficients': [str(x) for x in A2],
    'assignment1_kernel (p1,p2) coefficients': [str(x) for x in A1],
    'G = assignment1 + assignment2 (Khat part; K_0 parts cancel)': [str(x) for x in G],
    'G / {Pi_int^T,O}': str(ratio_G), 'Q(Y1a_T) / {Pi_int^T,O}': str(ratio_A2),
    'translation component - (D_P O)[infty] (F(P) = 0)': str(ratio_G * pi_T) + ' D_P O',
    'Q(Y1a_T)': str(ratio_A2 * pi_T) + ' D_P O',
    'engineB': '-1/3 D_P O', 'engineA_corrected (5/3 - 2)': '-1/3 D_P O'}
json.dump(out, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/compare_delta/k5_and_translation.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
