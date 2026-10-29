"""Assemble step2/engine_a/delta_results.json (Part D4, engine A) from the logs of d4_*.py."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
L = lambda f: json.load(open(os.path.join(HERE, 'logs', f)))
signs = L('d4_signs_su2.json')
pi = {s: L('d4_pi_%s_su2_p1201.json' % s) for s in ('cs_c', 'AsV_AV', 'ls_l')}
se = L('d4_se_su2.json')
ward = L('d4_y1a_ward_su2.json'); ward_gl3 = L('d4_y1a_ward_gl3.json')
L2 = '(log|p|^2 + gamma + log 2)'; L3 = '(log|p|^2 + gamma + log 3)'
L2e = '(log(eps|p|^2) + gamma + log 2)'; L3e = '(log(eps|p|^2) + gamma + log 3)'
res = {
 'role': 'engine A, Part D4 of DESIGN_STEP2.md (raw results written before reading step2/theory)',
 'units': 'unit = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)); colour universality checked (su(2) and gl(3) -> sl(3) give identical numbers)',
 'conventions': {
  'S_M': 'S^(1)(Delta) = int((2E(A)+2A,A^*) - (2E(lambda)+3lambda,lambda^*) - (2E(c),c^*)), M = v. for v = 2E',
  'D_M': 'derivation F(e) -> d/ds F(e + s M e); D_M F = 2 s d/ds|_{s=1} sigma_s F; E -> -(4 + p.d/dp) on Fourier components',
  'bracket_sign': '{S_M, -} = + D_M (measured exactly: {S_X,F}/(D_X F) = 1 for F = int(A,A), int(dA,dA), int(lambda,dslash lambda), int(dc,A^*))',
  'scaling_ward_identity': 'D_M W(P,J) = W_J\'(D_M J) + 4 eps d/deps|_J W(P,J)  (sign +; checked on every two-point function below)',
  'canonical_form': 'two-point functional = alpha(p) delta(Sigma) + beta^mu(p) (d_mu delta)(Sigma), field leg first, its momentum p fixed, Sigma = total momentum (the x-explicit S_M produces the d delta part)',
  'odd_parameters': 'coefficient of (xi, r_odd leg); Hessian of xi S_M = -xi Hess(S_M)',
  'Q_on_A*_leg': 'A^*|_P -> -(1/2) d^*d A|_P (s = +1, measured)',
  'tree_constants_measured': signs,
  'warning': 'Q C_P(S_X,F) = -D_X F for constant pairing-preserving X but Q C_P(S_M,F) = +D_M F (|p|^-2 M_{A^*} |p|^2 = M_A); component normalised so that its classical part is +2(E+w_O)O',
  'notation': 'L2 := log|p|^2 + gamma_E + log 2, L3 := log|p|^2 + gamma_E + log 3; eps-dependence only through eps|p|^2 (log eps coefficient = log|p|^2 coefficient in every entry)'},
 'task1_two_point_arity_one': {
  'sectors_that_occur': ['(A^*|_Wperp, A|_Wperp)', '(A^*|_V, A|_V)', '(lambda^*, lambda)', '(c^*, c)'],
  'why_only_these': 'd/dxi W has degree -1: two-point parts pair a field with an antifield; (A^*|_V, c) has degree 0 and belongs to W itself',
  'graph_classes': {
   'internal insertion Pi_int': 'S_M on the cycle; P has components A--A, c--A^*|_V, lambda--lambda only, so S_M sits on a cycle only through (MA,A^*) as a composite A|_V--c line, kernel P* M P* = (1/2){M,(P*)^2}, nonzero only because (Q^GF)^2 != 0',
   'external leg': 'S_M with one leg external, the other joined by P(eps,infty) to W_1^(2) = Sigma(eps) - I^CT_2 (per graph: Sigma and the counterterm vertex separately)'},
  'self_energies_W1_plane_wave_ren': {
   'A_Wperp': '(a,b) |p|^2 (43/8 - 2 L2)  [Sigma(eps) = (a,b)(9/2 eps^-1 + |p|^2(43/8 - 2 L2e))]',
   'A_V': '(a.b)|p|^2 (43/8 - 2 L2) + (p.a)(p.b)(2 L2 - 2)  [+ 9/2 eps^-1 (a.b)]',
   'lambda': '8 i (1 - L2) (u1, rho(p) u2)  (= (4 L2 - 4) x plane wave of int(lambda, dslash lambda))',
   'ghost (c,A^*_V)': 'i (gamma, p.b) (1 - L2)',
   'raw_fits': se},
  'A^*_Wperp,A_Wperp': {
   'structure': '(a,b) = sum_w (a_w, b_w); S_M itself = (alpha, beta) = (2, 2p)(a,b)',
   'external_leg_graph S_M-P-Sigma': 'alpha = 2g, beta = 2g p, g = 9/(eps|p|^2) - 4 log eps - 4 log|p|^2 - 4 gamma - 4 log 2 + 7/4 + O(eps log eps)',
   'external_leg_counterterm S_M-P-(-I^CT_2)': 'g = -9/(eps|p|^2) + 4 log eps + 9',
   'external_leg_total': 'g = 43/4 - 4 L2 (finite; the eps^-1/|p|^2 terms are nonlocal and cancel)',
   'Pi_int': 'alpha = 8/3, beta = 10/3 p  (finite, local; = (5/3) S_M|_P - (2/3) int(A_P,A^*_P)); recomputed at p = (2,-1,1,3) and with the field leg dependent (-32/3, -10/3 p): consistent',
   'total_eps^-1': '0', 'total_log_eps': '0', 'arity_one_counterterm I^CT_Delta|_P': '0',
   'total_finite R': 'alpha = 145/6 - 8 L2, beta = (149/6 - 8 L2) p',
   'Q(R)_sym (two-point obstruction along Delta, F_Delta = 0)': '-(11/3) unit int(dA|_P, dA|_P)  [= D_M W_1^ren (-4) + Q(Pi_int)_sym (+1/3)]',
   'F_Delta|_P solving the arity-one master equation': '(22/3) unit int(A|_P, A^*|_P) + f_S S_M|_P (f_S free: Q(S_M)_sym = 0)'},
  'A^*_V,A_V': {
   'structures': 'Lo := (p.a)(p.b)/|p|^2; beta: p^mu (a.b), a^mu (p.b), b^mu (p.a), p^mu Lo, eps^{mu m n s}(a_m,b_n)p_s (coefficient 0)',
   'external_leg (S_M A^*_V external, A_V leg on W_1^AV)': 'alpha = 2g_V, beta = 2 g_V p, g_V = (a.b)(43/4 - 4 L2) + Lo(4 L2 - 4)  [per graph: (a.b)(9/(eps|p|^2) - 9 + 43/4 - 4 L2e) + Lo(4 L2e - 4); counterterm graph -(a.b)(9/(eps|p|^2) - 9) + 4 log eps ((a.b) - Lo)]',
   'external_leg (S_M A_V external, A^*_V leg on the ghost self-energy)': 'alpha = (2 + 2 L2) Lo, beta = (2 L2 - 2) p Lo',
   'Pi_int': 'alpha = (a.b)(13/6 + L3e) + 2 Lo,  beta = 10/3 p (a.b) + (L3e - 1/2) a^mu (p.b) - 2 p Lo',
   'Pi_int_raw': pi['AsV_AV'],
   'total_eps^-1': '0',
   'total_log_eps': 'alpha = (a.b), beta = a^mu (p.b)  <=>  log eps * int(iota_x A|_V, d^*A^*|_V)',
   'total_finite R': 'alpha = (a.b)(70/3 - 8 L2 + L3) + Lo(10 L2 - 4), beta = p(a.b)(149/6 - 8 L2) + p Lo (10 L2 - 12) + a^mu(p.b)(L3 - 1/2)'},
  'lambda^*,lambda': {
   'structures': '(u,v), (u, chi v), p^mu(u,v), (u, rho(e^mu) rho(p) v) (chi = 4d chirality; chi coefficients all 0)',
   'external_leg': 'alpha = 6 phi, beta = 4 phi p, phi = 4 L2 - 4 [per graph phi = 4 L2e - 4, counterterm graph -4 log eps]',
   'Pi_int': 'alpha = -8/3, beta = p (-1/3 - 2 L3e) + (u, rho(e^mu) rho(p) v)(L3e - 1/2)',
   'Pi_int_raw': pi['ls_l'],
   'total_eps^-1': '0',
   'total_log_eps': 'beta = -2 p^mu (u,v) + (u, rho(e^mu) rho(p) v) = -(u, rho(p) rho(e^mu) v)',
   'total_finite R': 'alpha = 24 L2 - 80/3, beta = p(16 L2 - 2 L3 - 49/3) + (u, rho(e^mu) rho(p) v)(L3 - 1/2)'},
  'c^*,c': {
   'structure': '(gamma, w)',
   'external_leg': 'alpha = 0, beta = (2 - 2 L2) p [per graph (2 - 2 L2e) p, counterterm 2 log eps p]',
   'Pi_int': 'alpha = -8/3, beta = p (-5/6 - L3e)',
   'Pi_int_raw': pi['cs_c'],
   'total_eps^-1': '0', 'total_log_eps': 'beta = -p  <=>  -log eps int(E c, c^*)',
   'total_finite R': 'alpha = -8/3, beta = p (7/6 - 2 L2 - L3)'},
  'arity_one_counterterm_two_point': 'I^CT_Delta[eps] two-point = log eps [ int(iota_x A|_V, d^*A^*|_V) - int(E c, c^*) + (lambda term: beta = -(u, rho(p) rho(e^mu) v)) ] unit; no eps^-1; no A|_Wperp term. The (A_V, c) part equals log eps * Q int(iota_x A|_V, c^*): Q-exact, hence Q-closed as required by the log-eps part of the arity-one defect ({S_M, Y_2} = D_M Y_2 = 0).',
  'checks': {
   'Ward identity': 'every entry depends on eps only through eps|p|^2 (log eps coefficient = log|p|^2 coefficient; eps^-1 terms ~ 1/(eps|p|^2)); this is D_M = 4 eps d/deps|_J per graph; on W_1^(2)(A_Wperp): D_M Sigma = -18 eps^-1 - 8|p|^2 = 4 eps d/deps Sigma (the opposite sign fails); with the counterterm: D_M(Sigma - I^CT) = -D_M I^CT + 4 eps d/deps|_J (Sigma - I^CT), -D_M I^CT_2 = 9 eps^-1 unit int(A,A)',
   'momentum independence': 'A_Wperp sector at p = (1,2,0,1), (2,-1,1,3), (1,1,1,-1)',
   'dependent-leg swap (d delta machinery)': 'Pi_int(A_P): (8/3, 10/3 p) vs (-32/3, -10/3 p): exact',
   'regression': 'W_1^(2) eps^-1 and log eps coefficients = Step-1 I^CT in all four sectors'}},
 'task2_component': {
  'operator': 'O = int omega (iota_w1 A, iota_w2 A), omega = e^{i(k,x)} dvol; bilinear (A|_Wperp, A|_Wperp) part',
  'route': 'Lemma c7-image: [(d_CE+Q)O[infty]](Delta) = Q(arity-one part of O[infty] along Delta) (normalised: classical part +2(E+w_O)O); equivalent to the UV-boundary graphs (a)-(y) + counterterm + defect insertions by the UV-boundary identity (thm:c3-uv). Graph classes of the arity-one part (complete at one loop): B0 = O-P-S_M; B-S = S_M-P-(one-loop O-flow incl. dressed legs, O loops, d_delta I^CT trees); Y1a = 1PI triangle (O, composite S_M line, X = A_V(A_P)^2, Y = ([A,c],A^*)); Y1b = O-P-Pi_int; counterterm trees O-P-(I^CT_Delta + F_Delta).',
  'classical part': 'Q(B0) = +D_M O = 2(E+w_O)O: verified (canonical form alpha = -4, beta = -2k on the O structure)',
  'pieces': {
   'Q(B0 + B-S)': 'D_M O[infty]_0 = 2 s d/ds sigma_s O[infty] = 2(E+w_O)O[infty] + {-16 (traceless), -64 (trace)} unit O  [thm:c6-CS with Step-1 log coefficients 4/-8 and tree insertion 4 lambda per leg (rechecked)]',
   'Q(Y1b + counterterms) = Q C_P(Z,O)': 'Z = Pi_int - F_Delta fixed by the arity-one master equation: z_1 = -8 unit => +16 unit O (= -Q(4 unit tree insertion of 2int(A,A^*) - 4int(lambda^*,lambda) - int(c^*,c) into O)), plus z_2 D_M O (z_2 = S_M-component of Z, scheme freedom: z_2 = 5/3 if F_Delta has no S_M component)',
   'Q(Y1a)': 'alpha_Q = p2.beta, beta_Q = -(1/2)|p2|^2 beta with beta = unit T_O [8 p2/|p2|^2 + 8 |p1|^2 That(p1,p2)], alpha(Y1a) = 0 (to 1e-14)',
   'Y1a fits': {'su2': ward, 'gl3': ward_gl3}},
  'result': 'C = 2(E+w_O)O[infty] + 2 c unit O + z_2 unit 2(E+w_O)O + Q(Y1a),  c = 0 (traceless), c = -24 (trace)',
  'coefficient_multiplying_O': {'traceless': '0', 'trace': '-24',
     'meaning': 'coefficient c in C = 2(E + w_O + c unit)O[infty] + (rescaling z_2) + Q(Y1a); Q(Y1a) is an additional Q-exact one-particle-irreducible term'},
  'comparison_with_eq:c6-step2': 'C - [2 s d/ds sigma_s O[infty] - Q(4 unit tree(field-counting -> O[infty]))] = z_2 unit D_M O + Q(Y1a). (c6-step2) holds up to the scheme rescaling z_2 and the term Q(Y1a) != 0. Y1a exists only because (Q^GF)^2 != 0; with a gauge fixing satisfying (Q^GF)^2 = 0 the internal insertions vanish and (c6-step2) is exact.',
  'Q(Y1a) structure': 'Ward identity at X: bubble part (local after Q) + |p1|^2 x triangle. Symmetrised local part of the bubble: -4 O - D_M O. The |p1|^2 part is not a pure equation-of-motion term (That ~ p1/|p1|^2 as |p1|^2 -> 0), so this split does not define an O-coefficient; using it anyway would give -2 (traceless) and -26 (trace), contradicting the protection of the traceless bilinear.',
  'IR / well-definedness': 'Y1a UV and IR convergent at generic momenta (composite line ~ l^-3 at l -> 0, derivative part ~ l^-3); single Feynman-parameter monomials of the r-derivative terms are spuriously IR divergent after elimination of u_1, resolved by homogenisation; the x-explicit S_M makes all graph values distributions (delta, d delta), resolved by the canonical form; no IR scale used.',
  'locality_of_component': 'C is not of the form (local functional)[infty] at one loop: Q(Y1a) is nonlocal (logs/dilogs of k^2, |p1|^2, |p2|^2 ratios) and no local functional of weight 2 and degree -1 exists to absorb it; same mechanism for translations (beta(Y1a) is -i times the translation internal-insertion graph).'},
 'files': ['d4_signs.py', 'd4_treecheck.py', 'd4_cycles.py', 'd4_integrals.py', 'd4_two_point.py', 'd4_pi_sectors.py',
           'd4_se_sectors.py', 'd4_y1a.py', 'd4_y1a_ward.py', 'd4_notes.txt', 'logs/']}
json.dump(res, open(os.path.join(HERE, 'delta_results.json'), 'w'), indent=1)
print('written')
