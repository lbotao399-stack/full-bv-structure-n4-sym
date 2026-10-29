"""Assemble step2/engine_b/delta_results.json (Part D4, engine B) from the raw outputs of the delta_* scripts."""
import json, os
import sympy as sp

def L(f):
    return json.load(open(f)) if os.path.exists(f) else None

tests = L('delta_tests.json'); lt = L('delta_local_tests.json'); xe = L('delta_Xexact.json'); qb = L('delta_Qbasis.json')
cl = L('delta_classical.json'); k5raw = L('delta_K5_ApAp_s31.json'); k5cls = L('delta_K5_class_s0_mu0.json')
fb1 = L('delta_fit_B1.json'); ra = L('delta_routeA_Ap_n40.json'); rc = L('delta_routeA_c_n40.json')
ra2 = L('delta_routeA_Ap_n56.json'); rc2 = L('delta_routeA_c_n56.json')
idT = L('delta_identity_T.json'); flam = L('delta_fit_lam.json')
s1 = {k: L(f'delta_part1_sing_{k}.json') for k in ('Aperp_c', 'AV', 'lam')}

R = {}
R['statement'] = ('Engine B, Part D4 of DESIGN_STEP2.md: the dilatation along psl(4|4) at one loop, quadratic sector, N = 4, Feynman-gauge '
                  'Q^GF, UV boundary eps only (upper Schwinger limit infinity), no dimensional regularisation, exact rational arithmetic '
                  'for every number marked exact.')
R['units'] = 'unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)); every two-leg value multiplies the colour pairing (Z_1, Z_2) of the legs'
R['conventions'] = {
 'S_M': 'S_M = S^(1)(Delta), Delta = 2E: int (2E(A)+2A, A*) - int (2E(lambda)+3lambda, lambda*) - int (2E(c), c*)  (ex:c2-dil); M = v. (def:c7-vdot)',
 'E_in_momentum_space': 'on Fourier components E -> -(4 + p.d/dp); v.A -> -(6 + 2p.d/dp), v.lambda -> -(5 + 2p.d/dp), v.c -> -(8 + 2p.d/dp), '
                        'v.A* -> -(2 + 2p.d/dp), v.lambda* -> -(3 + 2p.d/dp), v.c* -> -2p.d/dp',
 'position_dependence': 'S_M = S_w + S_E, S_w = 2 int(A,A*) - 3 int(lambda,lambda*), S_E = sum_mu (-i d/dkappa_mu)|_0 T_mu[kappa], '
                        'T_mu[kappa] = 2 int e^{i kappa.x}[(d_mu A,A*) - (d_mu lambda,lambda*) - (d_mu c,c*)] (= 2 S^(1)(translation) injecting kappa)',
 'canonical_form': 'a two-leg functional on (leg 1 at p, leg 2 at p_2) is alpha(p) delta(p+p_2) + beta^mu(p) d/dp_2^mu delta(p+p_2) '
                   '<-> int dvol [(leg2, alpha(-i d) leg1) + x^mu (leg2, i beta_mu(-i d) leg1)]  (x placed on leg 2 = antifield leg in Part 1)',
 'ghost_convention': 'right ghost: J = I - I^CT[eps] - F_0 + (S_M - I^CT_Delta[eps] - F(Delta)) xi, component = coefficient left of xi '
                     '(the cochain evaluation of Def c2-dce); then [(d_CE+Q)O[infty]](Delta) = Q(arity-one part of O[infty]) with no sign',
 'sign_of_D_M': '{S_M, -} = + (derivation induced by e -> v.e).  Test: {S_M, int(A,A)} = %s int(A,A) (expected -4 = 2(2-4)); '
                '{S_M, -1/4 int(dA,dA)} = %s; Q S_M = 0 on all four two-leg sectors (machine, exact).' % (
                    lt['{S_M, int(A,A)} / int(A,A)'][0] if lt else '?', lt['{S_M, -1/4 int(dA,dA)}'][0] if lt else '?'),
 'D_M_on_canonical_two_point': 'for w_f + w_a = 4: alpha -> 2 p.d_p alpha, beta -> 2 (p.d_p - 1) beta;  on translation-invariant F(p) with '
                               'leg weights w_1, w_2: F -> 2(w_1 + w_2 - 4 + p.d_p) F',
 'word_formula': 'P only: W_1 coefficient = sum_words (1/2n) rsign colour STr[P Pi h ...]; one special edge Khat = K_eps - K_0: '
                 'hbar Delta_Khat W + 1/2{W,W}_Khat (1PI part) = sum_words sum_pos (1/2n) rsign (-1)^{odd legs before pos} colour STr; '
                 'S_M is a letter with one odd leg xi placed last (toy-verified formula of engine_b/toy_r2.py)',
}
# ---------------------------------------------------------------- Part 1
P1 = {}
P1['sectors'] = {
 'occurring': ['(A|_Wperp, A*|_Wperp)', '(A|_V, A*|_V)', '(lambda, lambda*)', '(c, c*)'],
 'excluded': '(A*|_V, c), (A, A), (lambda, lambda), (c, A*): wrong degree (the arity-one functional has degree -1: coordinate degrees '
             'A, lambda: 0; A*, lambda*: -1; c: +1; c*: -2); (A|_V, A*|_Wperp) and spinor-vector mixings: Spin(4) x Spin(6) covariance '
             '(Delta commutes with both).'}
P1['graph_classes'] = {
 '(i) S_M on the cycle (1PI)': 'only the block (A|_V, A*|_V) of S_M can sit on an internal line (lambda*, c*, A*|_Wperp have no propagator '
     'component): the internal S_M turns a line into an A|_V - c line.  One graph per sector: (A*_perp,A_perp): AcA*[A*_perp] - AAA[A_perp], lines '
     'A_perp and A_V-(S_M)-c;  (A*_V,A_V): two graphs (AcA*-AAA with A_V and A_V-(S_M)-c lines; AcA*-AcA* with c-A* and A_V-(S_M)-c lines);  '
     '(lambda*,lambda): c lambda lambda* - lambda lambda A, lines lambda and A_V-(S_M)-c;  (c*,c): ccc* - AcA*, lines c-A* and A_V-(S_M)-c.  '
     'No graph with a quartic vertex (the A_V - c line cannot close on AAAA).',
 '(ii) S_M on an external leg': 'B_P(S_M, W^(2)_1loop(eps) - I^CT_2[eps]) = S_M joined by P(eps,infty) to the renormalised one-loop two-point '
     'function: (A*_perp,A_perp) <- Sigma_perp; (A*_V,A_V) <- Sigma_VV and (via the c-A* edge) Sigma_(c,A*); (lambda*,lambda) <- Sigma_lambda; '
     '(c*,c) <- Sigma_(A*,c).  Finite as eps -> 0 (I^CT_2 tree included).',
 'closed form (D1(a))': 'W\'_J(S_M) = S_M + hbar d_P S_M + B_P(S_M, W) + (hbar/2) d_{K\'\'} W + (1/4) B_{K\'\'}(W,W), K\'\' the kernel of '
     '[d_P,[d_P,S_M]]; P*M P* = (M P*^2 + P*^2 M)/2 since [M,P*] = 4 eps L_eps* commutes with P*; at two legs, order hbar: class (ii) = '
     'B_P(S_M, W^(2)), class (i) = (hbar/2) d_{K\'\'} W_0^(4) (tadpoles vanish by colour).',
}
P1['eps_expansion'] = {
 'eps^-1': '0 in every sector (superficial degree 0 for alpha, 1 for beta: no power divergence)',
 'log eps (exact; I^CT[eps](Delta) on the quadratic sector = unit log eps X)': {
   'A|_Wperp': 'X = 0',
   'A|_V': 'alpha = (u_f,u_a)_V, beta^mu = u_f^mu (p,u_a)  <->  X_V = int (d(iota_x A|_V), A*|_V)',
   'c': 'alpha = 0, beta^mu = -p^mu u_c u_c*  <->  X_c = -int (iota_x dc, c*) = -int (E(c), c*)',
   'lambda': 'alpha = 0, beta^mu = -(rho(e_mu) rho(p) u_f, u_a) (both chiralities, fit %s)  <->  X_lambda = -int (rho(x) dslash lambda, lambda*)' % (
       flam['solution'] if flam else '?'),
   'total': 'X = int(d(iota_x A|_V), A*|_V) - int(iota_x dc, c*) - int(rho(x) dslash lambda, lambda*)',
   'Q-closed': 'Q X = 0 exactly on (A_V,A_V), (A_perp,A_perp), (c,A*_V), (lambda,lambda) (machine); Q X_V and Q X_c separately nonzero on '
               '(c,A*_V) and cancel',
   'Q-exact': 'X = Q( int(iota_x A|_V, c*) - (1/2) int(rho(x) lambda*, lambda*) ) exactly (machine, all sectors)',
   'raw_samples': s1,
 },
 'finite parts at generic momenta': {
   'class (i), A|_Wperp': 'no log; numerical limits alpha -> a0 (u_f,u_a), beta -> b0 p (u_f,u_a), i.e. class (i)_perp -> '
       'a0 int(A_perp,A*_perp) + b0 int(E(A_perp),A*_perp) = (5/3) S_M|_perp - (2/3) int(A_perp,A*_perp) for a0 = 8/3, b0 = 10/3: numerically a0 = 2.66644 (eps = 1e-7) -> 2.6666 +- 0.0002, b0 = 3.33326 -> 3.3333 +- 0.0001 (remainder O(eps log eps), Gauss-Legendre 56^3 in log u, eps = 1e-5, 1e-6, 1e-7, Schwinger representation without Taylor expansion); raw: see routeA',
   'class (i), c': 'beta = -p (u_c u_c*)[log(eps|p|^2) + c_c] + o(1): slope per log eps -0.9999 (eps 1e-6 -> 1e-7, exact value -1), c_c = 2.509 +- 0.001 (numerical; not identified); alpha -> a_c u_c u_c*, a_c = -2.6666 +- 0.0001 (-8/3)',
   'class (ii)': 'B_{P(0,infty)}(S_M, Sigma_ren), Sigma_ren = lim (W^(2)_1loop - I^CT_2) (Step-1 two-point functions, constant terms not '
       'recomputed here); its Q-image is exact: Q B_P(S_M, Sigma) = -B_{K_eps}(S_M, Sigma) -> {S_M, Sigma_ren} = D_M Sigma_ren = 4 x (the log eps '
       'coefficient of I^CT_2 in the sector), e.g. -4 int(dA_perp,dA_perp)',
   'routeA_raw': {'Aperp_n40': ra, 'c_n40': rc, 'Aperp_n56': ra2, 'c_n56': rc2},
 },
}
P1['scaling_Ward_identity'] = {
 'identity': 'D_M W(P(eps,infty),J) = W\'_J(D_M J) + 4 eps d/deps|_J W(P(eps,infty),J), D_M = {S_M,-} = + derivation of M; proof: '
             '[D_M, d_phi] = -d_{(M x 1 + 1 x M) phi} and (M x 1 + 1 x M) P(eps,infty) = 4 eps L_eps = -4 eps dP/deps (prop:c7-dil). '
             'With the opposite sign of D_M the right side would carry -4 eps d/deps.',
 'two-point check (arity zero)': 'W^(2) = Sigma - I^CT_2: D_M Sigma = 4 eps d/deps Sigma (homogeneity Sigma(sp, eps/s^2) = s^{4 - w_1 - w_2} Sigma), '
     'W\'(D_M J) = -D_M I^CT_2 = +9 unit eps^-1 int(A,A) (the eps^-1 coefficient; the log eps part of I^CT_2 has weight 4 and is D_M-closed)',
 'two-point check (arity one, class (i))': 'D_M acts on (alpha, beta) as (2p.d_p alpha, 2(p.d_p - 1) beta); numerically alpha(2p, eps/4) = alpha(p, eps) '
     'and beta(2p, eps/4) = 2 beta(p, eps) to 2e-5 relative (routeA scaling_check), i.e. D_M(class i) = 4 eps d/deps(class i); the extra '
     'term W\'\'(D_M J, S_M) = B_P(S_M, 9 unit eps^-1 int(A,A)) is the class-(ii) counterterm tree.',
}
R['part1'] = P1
# ---------------------------------------------------------------- arity-one obstruction (D2 input)
D2 = {}
D2['formula'] = ('Q W\'_J(S_M) = W\'_J(D_xi) + W\'\'_J(D_0, S_M) (UV-boundary identity, arity one); at two legs, order hbar: D_xi flow = '
                 '1PI special-edge graphs with S_M on the cycle (B1, exact masters) + {I^CT_2, S_M} term = +9 unit eps^-1 int(A,A) (D_M I^CT_2 = -9 eps^-1 '
                 'int(A,A)) - Q I^CT_Delta (= 0, Q X = 0) + Khat bridges (O(eps)) + S_M - P - (arity-zero two-leg defect, -> 0 after F_0).  '
                 'Ob(Delta)^(2) = lim [B1 + 9 unit eps^-1 int(A,A)].')
D2['B1_fits'] = fb1
D2['A_perp_result'] = ('Ob(Delta)|_(A_perp,A_perp) = -(11/3) unit int(dA_perp, dA_perp), no x-part, no log 2 / log(3/2) parts '
                       '(exact, 3 samples incl. two independent implementations of the eps-expansion); eps^-1 part of B1 = -9 unit int(A,A) '
                       'cancels the K_0 term exactly.')
D2['route_A_crosscheck'] = ('Independent: Q(class (i)) + Q(class (ii)) at eps -> 0 with class (i)_perp -> a0 int(A,A*) + b0 int(E A,A*) (numerical) '
                            'and Q int(A_perp,A*_perp) = -1/2 int(dA,dA), Q int(E A_perp,A*_perp) = +1/2 int(dA,dA) (exact): '
                            'Q(i) = -(1/2)(a0 - b0) int(dA,dA) = +0.33 int(dA,dA) (numerics: a0 - b0 = -0.667), Q(ii) -> D_M Sigma_ren = -4 int(dA,dA); '
                            'sum -3.67 = -11/3 (agreement to the numerical precision).  The translation part (b0) drops out of Q(i) on (A_perp, A_perp) '
                            'because int(a.dA, d*dA) is a total derivative; equivalently B1 with the translation letter vanishes identically '
                            '(identity test, delta_identity_T.json).')
D2['all_sectors'] = ('Ob(Delta) on the quadratic sector (exact; per sector 1-3 samples, residual-free fits, no x-dependent part, no log 2 / log(3/2) '
                     'constants, eps^-1 parts: -9 unit int(A,A) on (A,A) cancelling the K_0 term, 0 elsewhere): '
                     'unit [ -(11/3) int(dA,dA) (all ten components; coefficient of int(d*A_V,d*A_V) = 0) + (47/3) int(lambda, dslash lambda) '
                     '- (14/3) int(dc, A*) ]  (lambda: 2 samples, residual-free fit to the single structure int(lambda,dslash lambda); beta = 0)')
D2['F(Delta)_translation_invariant'] = ('unique translation-invariant quadratic solution of Q F = Ob(Delta): F(Delta) = unit [ (22/3) int(A,A*) '
                     '- (47/3) int(lambda*,lambda) - (8/3) int(c*,c) ] = unit [ (11/3) FC - int(lambda*,lambda) + int(c*,c) ], FC = 2 int(A,A*) - 4 int(lambda*,lambda) '
                     '- int(c*,c) (using Q int(A_perp,A*_perp) = -1/2 int(dA_perp,dA_perp), Q int(A_V,A*_V) = -1/2 int(dA_V,dA_V) - int(dc,A*), '
                     'Q int(lambda*,lambda) = -int(lambda,dslash lambda), Q int(c*,c) = -int(dc,A*)).  Q-closed directions (not fixed by the two-leg QME): '
                     'S_M|_perp, S_M|_lambda, S_M|_V + S_M|_c, X; they change the translation-invariant coefficients only through their weight terms.')
D2['F(Delta)_perp'] = ('Q F(Delta) = Ob(Delta) on (A_perp,A_perp): F(Delta)|_perp = unit [a int(A_perp,A*_perp) + b int(E(A_perp), A*_perp)] with '
                       'a - b = 22/3; b is the Q-closed direction int((E+1)A_perp, A*_perp) = (1/2) S_M|_perp, not fixed by the two-leg QME.  '
                       'Translation-invariant choice: b = 0, F = (22/3) unit int(A_perp,A*_perp) = (11/3) unit x (Wperp part of the field-counting functional).  '
                       'Choice cancelling the x-dependent part of K5 (below): b = -2/3, F|_perp = unit[(11/3) FC|_perp - (1/3) S_M|_perp], FC = 2 int(A,A*) '
                       '- 4 int(lambda*,lambda) - int(c*,c).')
R['arity_one_obstruction_quadratic'] = D2
# ---------------------------------------------------------------- Part 2
P2 = {}
P2['classical'] = {'statement': '{S_M, O} = 2(E + w_O) O verified in canonical form: alpha = -4 v, beta^mu = -2 k^mu v (v = vertex value of O, '
                                'k = momentum of omega), i.e. 2(E+w_O)O_k = -2(2 + k.d/dk)O_k; at order hbar^0 [(d_CE+Q)O[infty]](Delta) = '
                                '-B_{K_eps}(S_M,O) flowed -> {S_M,O}[infty]', 'samples': cl}
P2['decomposition'] = ('[(d_CE+Q)O[infty]](Delta) = d_xi^R d_delta lim [flow of D_{K_0} + hbar Delta_Khat W + 1/2{W,W}_Khat] (UV-boundary identity, '
    'D = (d_CE+Q)J + 1/2{J,J}_eps + hbar Delta_eps J, {,}_eps = {,}_{K_0} + {,}_Khat).  On (A_perp,A_perp), order hbar: '
    'K_0 part = (D_M O)[infty] - {I^CT_Delta + F(Delta), O} + O-P-[Ob(Delta) - Q F(Delta)] (+ S_M-P-(cocycle defect of O), zero on A_perp by Spin(6)); '
    'Khat part = K5 (1PI graphs with O, S_M and Khat on the cycle) + O-P-(B1) [inside Ob] + bridges (-> 0).  Hence '
    'LHS - 2((E+w_O)O)[infty] = lim [K5 - {F(Delta),O} - {I^CT_Delta,O}] with {I^CT_Delta, O} = 0 (X_perp = 0), given Q F(Delta) = Ob(Delta).')
P2['K5'] = {'fit': 'alpha = 16 T1 - 16 T2, beta^mu = (2/3) k^mu T1; T1 = (w1,u1)(w2,u2)+(w1,u2)(w2,u1) (= vertex value of O), T2 = (w1,w2)(u1,u2); '
                   'eps^-1 parts 0, log(3/2), log(4/3) parts 0; 3 samples, 15 equations, 6 unknowns, exact',
            'alpha_part': '16 T1 - 16 T2 = 4 x (log eps coefficient of partial_delta I^CT[eps]) = 4 [4 O - 2 (w1,w2) int omega (A_perp,A_perp)]',
            'beta_part': '(2/3) k T1 = -(1/3) x (beta part of 2(E+w_O)O): the translation component of K5 is -(1/3) unit {S^(1)(P), O}',
            'raw': k5raw}
P2['translation_part_check'] = {'statement': 'closed form for a translation letter T (T commutes with Q, Q^GF, H, P*; the kernel of P* T P* is Q-closed): '
    'LHS_P - (D_P O)[infty] = lim (hbar/2) d_{K\'\'}(Q O_W)_tree, i.e. only graphs with T on an A_V - c line between two P edges and Khat '
    'between O and an AcA* vertex; all other K5 graphs must cancel', 'result': k5cls['totals'] if k5cls else None,
    'meaning': 'subset = full translation part (844 on sample 0, = (2/3) k_0 T1), complement = 0 exactly: internal consistency of K5.  '
               'The nonzero translation part comes from the scale-eps BRST variation of O ({I,O}_eps - {I,O}) joined to the A_V - c line '
               '(which exists because (Q^GF)^2 != 0).'}
P2['coefficient'] = {
 'formula': 'with F(Delta)|_perp = a int(A,A*) + b int(E A, A*), a - b = 22/3: LHS - 2((E+w_O)O)[infty] = -(1/3 + b/2) 2(E+w_O)O + 2 c_O O with '
            'alpha: 16 T1 - 16 T2 - 2a T1 + 4b T1, beta: (2/3 + b) k T1;  2 c_O(T1 part) = 44/3 - 2(a - b) = 0 for every admissible b',
 'traceless': '0',
 'trace': '-24',
 'trace_detail': 'O_tr = int omega (A_perp, A_perp): T1 -> 2(u1,u2), T2 -> 6(u1,u2), v_tr = 2(u1,u2): alpha = 0 x T1 - 16 x 6 (u1,u2) = -96 (u1,u2) = -48 v_tr -> '
                 '2 c = -48, c = -24 (units); = 3 g^2 N/(4 pi^2) for su(N) with unit = -g^2 N/(32 pi^2)',
 'normalisation_of_the_classical_part': 'the factor of 2(E+w_O)O is 1 - (1/3 + b/2) unit: = 1 for b = -2/3 (F(Delta) contains -(1/3) unit S_M|_perp); '
                 '= 1 - unit/3 for the translation-invariant choice b = 0 (then the translation component of O is (1 - unit/3) D_P O as well)',
 'Step-1 scaling route prediction': '0 (traceless), -24 (trace): reproduced',
}
P2['not_computed'] = ['the (A_V, A_V) and (c, A*_V) sectors of the trace component (K5 with A_V legs, the cocycle-defect tree S_M-P-(cocycle '
                      'obstruction of O_tr), and the finite counterterm g (w1,w2) int omega (A_V,A_V) of Remark c5-cocycle)',
                      'the cubic arity-one obstruction that would fix b (the S_M|_perp component of F(Delta))']
R['part2'] = P2
R['comparison_with_theory'] = ('step2/theory/delta_predictions.json did not exist when this file was written; no comparison made. '
  'Step-1 scaling-route values (STORY.md, ch06 Theorem thm:c6-dil): coefficient 0 (traceless), -24 (trace): reproduced (part2.coefficient).')
R['checks'] = {'regressions': tests, 'Q_basis': qb, 'X_exactness': xe, 'local_tests': lt, 'identity_T': idT}
R['machine'] = 'one python process at a time, RSS < 0.2 GB, no Wolfram kernel, no files outside the project folder'
json.dump(R, open('delta_results.json', 'w'), indent=1, default=str)
print('written', len(json.dumps(R, default=str)))
