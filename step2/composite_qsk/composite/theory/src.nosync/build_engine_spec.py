"""Build composite/theory/engine_spec.json (composite kernels C2_Delta, C2_Pb; Conjecture P7; P_1).
Reads census_sp.json, ha_predictions.json, ha_su11.json from this folder."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
THEORY = os.path.dirname(HERE)
census = json.load(open(os.path.join(HERE, 'census_sp.json')))
ha = json.load(open(os.path.join(HERE, 'ha_predictions.json')))
su11 = json.load(open(os.path.join(HERE, 'ha_su11.json')))

field_pairs = {}
for key in ['AV|AV', 'AV|AW', 'AV|lamP', 'AV|lamM', 'AW|AW', 'AW|lamP', 'AW|lamM', 'lamP|lamP', 'lamP|lamM', 'lamM|lamM']:
    v = census['pairs'][key]
    field_pairs[key] = {
        'W0': v['W0'],
        'leg_sets_W0_Pb_and_position_part': {s: [d['k_min'], d['chains'], d['w_out']] for s, d in v['by_letter_derivatives']['0']['leg_sets_Pb_and_position'].items()},
        'leg_sets_W0_translation_invariant_part': {s: [d['k_min'], d['chains'], d['w_out']] for s, d in v['by_letter_derivatives']['0']['leg_sets_translation_invariant'].items()},
        'counts_n0..3_Pb': [len(v['by_letter_derivatives'][n]['leg_sets_Pb_and_position']) for n in '0123'],
        'counts_n0..3_translation_invariant': [len(v['by_letter_derivatives'][n]['leg_sets_translation_invariant']) for n in '0123'],
    }

spec = {
 'title': 'Engine specification: composite kernels C^(2)_Delta, C^(2)_{P_b} on all local operators; Conjecture P7 in testable form; P_1',
 'theory': 'step2/composite_qsk/composite/theory/composite_theory.tex (+ _build/composite_theory.pdf); builds on '
           'step2/twoletter/theory/twoletter_theory.tex and manuscript Ch7-Ch10',
 'conventions': {
   'manuscript': 'Ch1 (fields, action, Q), Ch3 (Q^GF, K_t, P(eps,T)), Ch4 (special edge K_eps - K_0, Prop c4-defect1; defect flow Prop c4-defectflow), '
                 'Ch6 (Lemma c6-bvexact), Ch7 (v., Prop c7-cov, Prop c7-dil), Ch8 (Def c8-components, Thm c8-L, Cor c8-homotopy), '
                 'Ch9 (composite kernel Def/Prop c9-composite, Thm c9-quadratic, Rem c9-Pi, Thm c9-D1, Prop c9-triangle, Thm c9-P1, Prop c9-scheme), '
                 'Ch10 (pair kernels, Thm c10-D1, Prop c10-P1, Prop c10-scalarsector, Cor c10-P7)',
   'regulator': 'Schwinger cutoff t_e >= eps on every edge P(eps,T); special edge t in [0, eps]; composite line = two segments, each t >= eps; '
                'T = infinity at generic momenta (Prop c8-Tsafe), results T-independent (Thm c8-L (ii)); NO dimensional regularisation in any role; no IR regulator other than T',
   'unit': 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2) for su(N); planar adjacent rule of Lemma c10-colour (iv): a term a*hbar/(16 pi^2) (1 (x) ad_{Y_r}..ad_{Y_1})(-,-)^{-1} gives (-1)^r a/2 units on an adjacent pair',
   'dilatation': 'v = 2E; v. = 2E + 2w on a slot of Weyl weight w (w = 0,1,3/2,3,5/2,4 for c,A,lambda,A*,lambda*,c*); {S^(1)(Delta),M} = v.M; v.O = 2(E + w_O)O on int omega(D)',
   'translations': 'v. = partial_b on every slot; {S^(1)(P_b), int omega D} = int omega partial_b D',
   'two_component': 'Def c8-twocomp / Prop c8-twocomp: sigma^m = (sigma^1,sigma^2,sigma^3,i), <A_m A_n> = 2(-,-)^-1 delta_mn, <A_ij A_kl> = 4(-,-)^-1 eps_ijkl, '
                    '<c A*_m> = -(-,-)^-1 i p_m, <lambda_alpha^i lambdabar_betadot j> = -(-,-)^-1 p_m sigma^m delta^i_j, each times e^{-eps p^2}/p^2 at T = infinity; '
                    'step2/notation/dictionary.json for engine index maps',
 },
 'definitions': {
   'square_kernel': 'X := ((P*)^2 (x) 1) K_0; components c(x)A|_V and A|_V(x)c only; at T = inf (Q^GF)^2 = 2 d^* on A*|_V times e^{-2 eps p^2} p^-4 (a doubled edge, momentum degree -3)',
   'composite_kernels': 'K_x := (P* v_x. (x) 1 + 1 (x) P* v_x.) P = 2 (P* v_x. P* (x) 1) K_0 = (v_x.(x)1 - 1(x)v_x.) X;  K_Delta = 2(v.(x)1)X + 4(eps d_eps + T d_T)X;  K_{P_b} = 2(d_b (x) 1)X  (Lemma lem:cc-kernels)',
   'internal_vertex': 'in a graph K_x is S^(1)(x) with both legs on edges P: only int(v_x.A|_V, A*|_V) contributes; leg A|_V on an A-A edge, leg A*|_V on an A*|_V-c edge. '
                      'Along Delta: int((2E+2)A|_V, A*|_V) = int x^mu s_mu + int s_0 with s_mu = 2(d_mu A|_V, A*|_V) (twice the internal vertex density of S^(1)(P_mu)) and s_0 = 2(A|_V, A*|_V). Along P_b: int(d_b A|_V, A*|_V)',
   'special_edge': 'K_eps - K_0 = -int_0^eps (H (x) 1) K_t dt = (e^{-eps p^2} - 1) x components of K_0 (field slot with its antifield slot: A(x)A*, lambda_+(x)lambda*_+ [S+(x)W with S+(x)W*], lambda_-(x)lambda*_- [S-(x)W* with S-(x)W], c(x)c*, and transposes); numerator -|l+q|^2 e^{-t|l+q|^2}, t in [0,eps]',
   'Z_x': 'Z_x(X_a,X_b;eps) := minus (hbar/2) d_{K_x} of the amputated chains X_a - ... - X_b with k >= 1 vertices of I, exactly one edge = K_x, exactly one other edge = special edge, all remaining edges P(eps,infty); signs of Lemma c10-sign along the chain from X_a (Def def:cc-sp). By Prop prop:cc-master, for sO = 0, Q(composite graphs of O) = sum of Z_x with trees',
   'C2_x': 'C^(2)_x(X_a,X_b) := coefficient of eps^0 of Z_x(X_a,X_b;eps): local, two-body (Thm thm:cc-local); this is the two-body composite kernel. Ch10/twoletter object tilde C^(2)_Delta := Q-image of renormalised composite cycles through the pair (nonlocal); C^(1)_Delta(X) = {Pi_Delta, X} (nonlocal on A|_V, lambda, c). On cocycles sum tilde C^(2) + sum C^(1) = sum C^(2)',
   'canonical_form': 'int omega D with all derivatives on fields; omega-linear M(int omega D) = int omega (x^b M^b(D) + M^0(D)); plane waves: field momenta independent, omega (and the injected r) absorbs the total momentum; value alpha delta + beta^b d_b delta <-> M^0 = alpha, M^b = i beta^b (Def def:cc-canonical)',
   'position_identity': '(C^(2)_Delta)^b = 2 C^(2)_{P_b} graph by graph; same for ct^(1), I^CT(1), C_Delta (Prop prop:cc-position)',
   'ct1_x': 'ct^(1)_x(M) = hbar Sing_eps of the plain composite cycles through pairs of M (one K_x line, other edges P); lc^(1)_x its log eps coefficient; degree d_ab - 1',
   'L1_x': 'L^(1)(x) := log eps coefficient of I^CT_psl[eps]^(1)(x) (composite loops of I); quadratic: L^(1)(Delta) = u Q G_Delta, L^(1)(P_b) = (1/2) u Q G_{P_b}; no eps^-m terms at all (Thm thm:cc-local (ii))',
   'G': 'G_Delta = int(iota_x A|_V, c*) - (1/2) int(rho(x) lambda*, lambda*);  G_{P_b} = int(iota_b A|_V, c*) - (1/2) int(rho(b) lambda*, lambda*); two-component G = int x_m(A_m,c*) + i int x_m(lambda*^alpha_i, sigma^m_{alpha alphadot} lambdabar*^{alphadot i})',
   'free_directions': "hbar F_1^(1)(Delta) = 4uF + F'_Delta, F'_Delta = u[-2/3 int(A,A*) + 1/3 int(lambda*,lambda) + 4/3 int(c*,c)] + u sum_d theta_d S_d + u Q(kappa_V int(iota_x A|_V,c*) + kappa_lambda int(rho(x)lambda*,lambda*)) + cubic...; "
                       "hbar F_1^(1)(P_b) = u sum_d theta^P_d S_d(P_b) + u Q(kappa^P ...) + ...; covariant scheme: theta^P = theta, kappa^P = kappa",
   'N': "N(int omega X_1..X_n) = sum_a (2 + |alpha_a|) int omega X_1..X_n for field letters X_a = d^{alpha_a} phi_a (canonical form); at theta = -1/3 the one-body term -{F'_Delta, X}^0 = (2/3) u (2 + |alpha|) X for every field slot (Lemma lem:cc-onebody)",
 },
 'theorems_to_use': {
   'locality': 'Thm thm:cc-local: for sO = 0, C_x(O) = sum_{a<b} +- C^(2)_x(X_a,X_b) (no one-body term); no log eps in Z_x; poles: sum z_m = s(eps^-m part of ct^(1)_x(O)); log covariance s lc^(1)_x(O) = -{L^(1)(x), O}',
   'D1_local': "Cor cor:cc-D1local: D_1(O) = sum_pairs [4 K^log + C^(2)_Delta] + sum_letters [-4u n_F X - {F'_Delta, X}];  P_1(O) = sum_pairs C^(2)_{P_b} - sum_letters {hbar F_1^(1)(P_b), X}",
   'arity_two': 'Thm thm:cc-arity2: position part of D_1 is s-exact on cocycles; Cor cor:cc-theta: 2 C_{P_b} - {(F_Delta)^b, -} exact; equal thetas => C_{P_b} ~ u theta {S^(1)(P_b),-}; covariant scheme: D_1^b = 2 P_1, arity two <=> P_1 ~ 0',
   'P7': 'Conj conj:cc-P7: (a) cubic arity-one forces theta_perp = theta_lambda = theta_V; (b) sum C^(2)_{P_b} ~ -(1/3) u {S^(1)(P_b),-}; (c) sum C^(2),0_Delta ~ -(2/3) u N, on the cohomology of each family. Thm thm:cc-P7: (a)-(c) => Xi_Delta ~ -(1/3+theta) u v., then theta = -1/3, D_1 ~ 4 lc^natural, Spec D_1 = -8 u E_H',
 },
 'census': {
   'file': 'step2/composite_qsk/composite/theory/src.nosync/census_sp.json (src.nosync/census_sp.py, base step2/twoletter/theory/checks_census.py)',
   'summary': '45 of 55 unordered slot pairs carry special-edge composite cycles; the 10 without: pairs of two letters among A*|_perp, lambda*_+, lambda*_-, c*. k = 1 exactly for pairs with a letter A|_V and for (c, antifield). 0 output-bound violations',
   'field_pairs': field_pairs,
   'leg_set_rules': 'translation-invariant part along Delta: chains with w_out <= w_ab; position part along Delta and P_b part: w_out <= w_ab + 1; n_out <= w_ab + 1 + d_ab; outputs carry the pair degree d_ab',
 },
 'closed_forms': {
   'simplex': 'Thm thm:cc-simplex: k+2 Schwinger times (composite = two segments), R_sigma = {alpha_sigma smallest}; 16pi^2 [eps^-m] = (1/m) int_{R_sigma} [(min_{f!=sigma} alpha_f)^m - alpha_sigma^m] F_m; 16pi^2 [eps^0] = int_{R_sigma} F_0 log(min_{f!=sigma} alpha_f / alpha_sigma); no log eps; '
              'F_m = sum_{j >= m+k} (Lap^j num)(centre)/(j! 4^j) (-Phi)^{j-k-m}/(j-k-m)!',
   'genfun': 'Prop prop:cc-genfun: Bessel generating function in y = y_b - y_a with index m+k: (-1)^{m+k} e^{i(y,centre)} (|y|/(2 sqrt Phi))^{m+k} I_{m+k}(|y| sqrt Phi), chain numerator with l -> -i d_y',
   'k1': 'Prop prop:cc-k1: beta = alpha_sigma in [0,1/3]; W_0(beta) = (1-beta) log((1-beta)/(2 beta)) - (1 - 3 beta); W_m(beta) = (2/m)[(((1-beta)/2)^{m+1} - beta^{m+1})/(m+1) - beta^m((1-beta)/2 - beta)]; r-linear part: factor (1-beta)/2; values in Q + Q log2 + Q log3',
   'exact_basis': 'k >= 2: exact sympy integration over ordered subsimplices {alpha_sigma <= alpha_f <= others} (homogenise, Laplace representation); else >= 40-digit numerics + integer relation on {1, log2, log3, log5, pi^2, Li_2(1/2), Li_2(1/3), Li_2(1/4)} confirmed by a second evaluation',
 },
 'tasks': [
   {'id': 'TC0', 'what': 'completeness of the Q-closed quadratic local functionals of degree -1, Weyl weight 4 (position <= 1) in each sector (expected: S_perp; S_lambda, Q int(rho(x)lambda*,lambda*); S_V, Q int(iota_x A|_V,c*)); same along P_b with weight 5, no position. Pure linear algebra on plane waves.', 'output': 'tc0_directions.json'},
   {'id': 'TC1', 'what': 'regression: GC1-GC4 (quadratic L^(1), composite part, Pi_Delta|perp, composite triangle 8/8, bilinear C^(2) = -(1/3) u v. - (4/3) u and C^(2)_{P_b} = -(1/3) u d_b, position identity)', 'output': 'tc1_regression.json'},
   {'id': 'TC2', 'what': 'C^(2)_{P_b} and C^(2)_Delta (canonical: position and translation-invariant parts) for the 10 field pairs, letter derivatives n <= 3 (n <= 2 for route B), all leg sets of census_sp.json within the bound, general colour position and planar adjacent form; plus the pairs (d^a A|_V, d^b X) with ONE output (gauge completion: leg sets A|_perp of (A|_V,A|_perp), lambda of (A|_V,lambda), A|_V of (A|_V,A|_V)); the plain composite cycles and ct^(1)_x for the same pairs', 'output': 'tc2_kernels.json'},
   {'id': 'TC3', 'what': 'families F0-F3 of Table tab:cc-families: matrices of sum C^(2)_{P_b} and sum C^(2),0_Delta on a basis of the two-field part, projected on shell (|p|^2 = 0 scalars, Weyl equation gauginos, transverse on-shell gauge legs), descendants kept; test Conjecture P7 (b), (c); D_1 at theta = -1/3 and its spectrum vs -8 u E_H (ha_predictions.json, ha_su11.json, twist two -32 u h(S))', 'output': 'tc3_families.json'},
   {'id': 'TC4', 'what': 'pole and log gates GC5, GC6 on the families: Z_x poles = s(ct^(1)_x poles); s lc^(1)_x(O) = -{L^(1)(x),O}; lc^(1)_Delta(O) - u{G_Delta,O} s-closed', 'output': 'tc4_covariance.json'},
   {'id': 'TC5', 'what': 'cubic arity-one test (GC13): rank of {I_3,S_perp}, {I_3,S_lambda} modulo Q(cubic local functionals of degree -1, weight 4, position degree <= 1); if rank 2, compute the cubic composite part (three-leg composite loops of I with one special edge) and solve for (delta_perp, delta_lambda); same along P_b', 'output': 'tc5_cubic.json'},
   {'id': 'TC6', 'what': 'two routes (GC8): route of the special edges (all cycles through the pair with one special edge and one S^(1)(Delta), Lemma c8-transfer) on the bilinears and twist two S\' <= 2, compared with sum [4 K^log + C^(2)_Delta] (K^log from step2/twoletter/results/kernels_final.json)', 'output': 'tc6_routes.json'},
   {'id': 'TC7', 'what': 'locality demonstration (GC7, route B): on twist two S\' = 2 and on the su(1|1) pair (Z,psi), evaluate sum tilde C^(2)_Delta and sum C^(1)_Delta separately at >= 3 kinematic points (nonlocal: log|p|^2 terms), and their sum against sum C^(2)_Delta', 'output': 'tc7_locality.json'},
 ],
 'families': {
   'F0': {'what': 'scalar bilinears; planar scalar sector; su(2) L <= 8', 'expect': 'C^(2)_Delta(bilinear) = -(1/3) u v.O - (4/3) u O (traceless and trace); per adjacent pair -(1/6) u v.(X_a X_b) - (2/3) u X_a X_b; C_{P_b} = -(1/3) u d_b; D_1 = 0, -48 u (bilinears) at theta = -1/3; D_1 = -8 u * spectrum of sum(1 - P + K/2)'},
   'F1': {'what': "twist two tr((iota_n d_A)^k Z (iota_n d_A)^{S'-k} Z), n = (0,0,1,i) or (1,i,0,0), Z = iota_w A, (w,w) = 0, S' <= 8 (route B S' <= 4)", 'expect': "Spec D_1 = {-32 u h(S): S even <= S'}; matrix D_1 = -8u(H'_12 + H'_21) [B03c (3.14)] on shell in basis D^kZ D^{S'-k}Z/(k!(S'-k)!); (b): sum C^(2)_{P_b} ~ -(1/3) u d_b; (c): sum C^(2),0 ~ -(2/3) u (4 + S')"},
   'F2': {'what': 'su(1|1) tr(Z^{L-M} psi^M) L <= 6 (Z = A_5 + iA_6 / Phi_34, psi = lambda component of so(6) weights (1/2,1/2,1/2)); fermion sl(2) D^k psi (one R component); two components psi^(3), psi^(4) with tr(Z F_11)', 'expect': {'su11_E_H': {k: v for k, v in su11.items() if k.startswith('L')}, 'fermion_sl2_E_H': ha['fermion_Psi1_c4'], 'two_component_E_H': ha['fermion_Psi1_c3_c4'], 'D_1': '-8 u E_H'}},
   'F3': {'what': 'gluon sl(2) D^k F_11 (F_11 = self-dual F, alpha=beta=1; two-field part d A|_V, three-field [A,A]); mixed (F_11, psi)', 'expect': {'gluon_E_H': ha['gluon_F11'], 'mixed_E_H': ha['mixed_F11_Psi1'], 'D_1': '-8 u E_H; tr(F_11 F_11): -48 u'}},
   'harmonic_action_source': '[B03c] (4.9), (4.27)-(4.29); implementation src.nosync/ha_predict.py (checks: twist two 4h(S), scalars 0^20 and 6)',
 },
 'routes': {
   'A': {'name': 'Schwinger/simplex exact (extends step2/twoletter/engine_a)',
         'method': ['chain numerator from the vertex functions of step2/twoletter/engine_a (ea_core, ea_kern; vendored Step-1 r2_core, cliff); new edge types: special edge (K_0 components, t in [0,eps], numerator -|l+q|^2) and composite line (two segments A-A and A*|_V-c, internal vertex of Lemma lem:cc-kernels (iv))',
                    'position along Delta by the injected momentum r at the internal vertex (absorbed by omega): alpha = G^{s_0} - i d_r G^{s_mu}|_0 (r-derivative on the integrand), position part from G^{s_mu}(0) = 2 x (P_mu graph)',
                    'Gaussian loop integral, Thm thm:cc-simplex with exact Dirichlet-log integrals over ordered subsimplices of R_sigma; letter derivatives via Prop prop:cc-genfun; k = 1 cross-check Prop prop:cc-k1',
                    'colour: (1 (x) ad...)(-,-)^-1 stripping as engine_a; planar adjacent form'],
         'output_dir': 'step2/composite_qsk/composite/engine_a/ (code and caches in src.nosync/, out.nosync/)'},
   'B': {'name': 'momentum-space numerator algebra + Feynman parameters + numerical eps-extrapolation (extends step2/twoletter/engine_b and step2/engine_b delta_*)',
         'method': ['composite line in kernel form K_Delta = 2(v.(x)1)X + 4 eps d_eps X with v. = -2(4 + p.d_p) + 2w on Fourier components (Euler operator on the ghost end), K_{P_b} = 2(d_b (x) 1)X',
                    'special edge as the explicit factor (e^{-eps p^2} - 1); Feynman parameters for the other edges; one-dimensional quadrature in the special Schwinger time',
                    'eps in {1e-2, 5e-3, 2e-3, 1e-3, 5e-4} (>= 5 values), >= 3 kinematic points per leg set; fit sum_m eps^-m z_m + z_0 + eps(a + b log eps); polynomial ansatz of Prop prop:cc-bound; integer relation for rational/log2/log3 coefficients',
                    'canonical form: alpha delta + beta d delta with field momenta fixed (engine_b delta_* convention, swap formula alpha(p2) = alpha(-p1-k) - d.beta)'],
         'output_dir': 'step2/composite_qsk/composite/engine_b/ (code and caches in *.nosync/)'},
   'independence': 'A: vertex form + injected momentum; B: kernel form + Euler operator; neither reads the other before the comparison phase',
 },
 'gates': [
   'GC1 quadratic L^(1)(Delta) = u Q G_Delta, no eps^-1 (Thm c9-quadratic (i)); L^(1)(P_b) = (1/2) u Q G_{P_b}',
   'GC2 composite part u[1/3 int(dA,dA) - 1/3 int(lambda,dslash lambda) - 2/3 int(dc,A*)]; Pi_Delta|perp = u[5/3 S_perp - 2/3 int(A_perp,A*_perp)]; composite triangle coefficients 8, 8; quadratic composite part along P_b = 0',
   'GC3 bilinears: C^(2)_Delta(A_perp,A_perp)|_{A_perp A_perp} = -(1/3) u v.O - (4/3) u O, traceless and trace; C^(2)_{P_b} = -(1/3) u d_b O; one chain A_perp -K- ([A_perp,c],A*_perp) -K_Delta- A_V(A_perp)^2 -P- A_perp',
   'GC4 position identity (C^(2)_Delta)^b = 2 C^(2)_{P_b} graph by graph',
   'GC5 no log eps in Z_x; pole identity z_m = s(ct^(1)_{x,m}) summed over pairs on cocycles',
   'GC6 log covariance s lc^(1)_x(O) = -{L^(1)(x), O}; cubic part of L^(1)(Delta) = u {I_3, G_Delta}',
   'GC7 locality: sum tilde C^(2) + sum C^(1) (each nonlocal) = sum C^(2) on cocycles',
   'GC8 two routes: special-edge route (Thm c9-D1 first row on bilinears: 4 lc - 4/3 u O - 1/3 u {S,O}) = sum [4 K^log + C^(2)_Delta]',
   'GC9 Conjecture P7 (b) per family',
   'GC10 Conjecture P7 (c) per family',
   'GC11 spectra at theta = -1/3 vs -8 u E_H (harmonic action)',
   'GC12 colour realisations su(2), su(3), gl(3) of the planar adjacent rule',
   'GC13 cubic arity-one rank test; completeness of quadratic Q-closed directions (TC0)',
 ],
 'outputs_expected': {
   'engine_a': ['kernels_C2.json (per pair/leg set/chain: z_m, C2 position + translation-invariant, colour general + planar)', 'families.json', 'gates.json', 'PROGRESS.txt', 'notes.txt'],
   'engine_b': ['kernels_C2.json (same schema)', 'families.json', 'gates.json', 'PROGRESS.txt', 'notes.txt'],
   'schema_note': 'coefficients as strings of exact rationals in the basis {1, log2, log3}; momenta and polarisations as in step2/twoletter/engine_a kernels.json; units u',
 },
 'machine': {'memory': 'each process < 1.5 GB RSS, one heavy process at a time per agent', 'tmp': 'never /tmp; write under step2/composite_qsk/composite/; code and caches in *.nosync subfolders (iCloud eviction); write PROGRESS.txt incrementally', 'language': 'pure Python (fractions, sympy, numpy, mpmath) or C; cluster-portable; no Wolfram, no cluster, no browser', 'files': 'subagents write .json/.txt/.tex/.py, no .md'},
 'checks_done': {
   'census_sp': 'src.nosync/census_sp.py -> census_sp.json: 45 pair types, degree/parity/Z4 asserted on every chain, 0 bound violations',
   'harmonic_action': 'src.nosync/ha_predict.py -> ha_predictions.json (twist two 0,6,25/3,49/5 = 4h(S); scalars 0^20, 6), ha_su11.py -> ha_su11.json',
   'tex': 'composite_theory.tex compiled with pdflatex (_build/composite_theory.pdf), no undefined references',
 },
}
out = os.path.join(THEORY, 'engine_spec.json')
json.dump(spec, open(out, 'w'), indent=1)
json.dump(spec, open(os.path.join(HERE, 'engine_spec.json'), 'w'), indent=1)
print('written', out, len(json.dumps(spec)))
