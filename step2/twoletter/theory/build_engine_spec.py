"""Builds engine_spec.json (twoletter_theory.tex, Section sec:tl-engine) from census.json and the theory.
Run after checks_census.py."""
import json, os
from fractions import Fraction as Fr
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
C = json.load(open(os.path.join(HERE, 'census.json')))
W = {'c': Fr(0), 'AV': Fr(1), 'AW': Fr(1), 'lamP': Fr(3, 2), 'lamM': Fr(3, 2), 'AsV': Fr(3)}
DEG = {'c': 1, 'AV': 0, 'AW': 0, 'lamP': 0, 'lamM': 0, 'AsV': -1}

kernels = {}
for key, legsets in C['output_leg_sets_by_pair_weight'].items():
    a, b = key.split('|')
    W0 = W[a] + W[b]
    explicit = [s for s in C['chains_listed'][key] if int(s.split('|')[0].split('=')[1]) <= 3]
    kernels[key] = dict(
        letters=[a, b], pair_weight_without_derivatives=str(W0), degree=DEG[a] + DEG[b],
        output_bound='n_out <= W + %d' % (DEG[a] + DEG[b]),
        chains_k_le_3_wout_le_5=explicit,
        chain_counts_by_k_and_wout=C['chain_counts'][key],
        leg_sets_by_pair_weight={w: dict(n=v['letter_derivatives'], log=v['log'], eps_max_m=v['eps_max_m'])
                                 for w, v in legsets.items() if Fr(w) <= 6},
    )

spec = dict(
    title='Engine specification: universal two-letter kernels and D_1 on all local operators (R8)',
    theory='step2/twoletter/theory/twoletter_theory.tex (+ _build/twoletter_theory.pdf)',
    conventions=dict(
        manuscript='Ch1 (fields, action, Q), Ch3 (Q^GF, K_t, P(eps,T)), Ch4 (I^CT, unit, 27/16), Ch5 (O(eps)), '
                   'Ch6 (dilatation), Ch8 (Def c8-quant, def:c8-renorm, components), Ch9 (composite kernel, D_1, P_1)',
        regulator='per-edge Schwinger cutoff t_e >= eps on every edge P(eps,T); T = infinity (generic momenta) or '
                  'finite T (any momenta): counterterms identical (Prop prop:tl-map (iii)); NO dimensional regularisation',
        unit='hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2) for su(N), (X,Y) = -(4 hbar/g^2) tr(XY)',
        two_component='Def c8-twocomp; step2/notation/dictionary.json; every reported kernel also in components'),
    definitions=dict(
        counterterm_map='d_delta I^CT[eps](M) = sum_m eps^-m ct_m(M) + log eps lc(M) = hbar Sing_eps d_{P_e} M (Thm thm:tl-decomp)',
        two_letter_kernel='K_eps(X_a,X_b) = hbar Sing_eps [(d^alpha (x) d^beta) P_e^{ab}(x,y)]_{y=x}, P_e = sum_k P (I''_e P)^k',
        signs='sigma_ab = (-1)^{|X_b|(|X_{a+1}|+...+|X_{b-1}|)}; component of P in the order (a,b) (Def c3-P)',
        colour='(1 (x) w)(-,-)^{-1}, w = ad_{Y_r}...ad_{Y_1} (Y_1 nearest X_a); su(N) single trace: adjacent -> (-1)^r N tr(M Y_1..Y_r) '
               '+ multi-trace; general position -> multi-trace only; across traces -> merged trace tr(M1 Y_{S^c}^-> M2 Y_S^<-), no N; '
               'wrap-around pair (X_1,X_n) -> +N tr(X_2..X_{n-1} Y_r..Y_1) (Lemma lem:tl-colour)',
        lc_natural='lc(M) - u {2 int(A,A*) - 4 int(lambda*,lambda) - int(c*,c), M}; commutes with s = Q + {I,-}',
        D1='D_1 = 4 lc^natural - {F_Delta prime, -} + C_Delta  (Thm thm:tl-D1); F_Delta prime with theta_perp, theta_lambda, theta_V, kappa_V, kappa_lambda (Def def:tl-H)'),
    kernels=kernels,
    weight_cutoffs=dict(
        explicit_tables='all 21 pair types, pair weight W <= 6 (leg sets above); counts to W <= 7 (census.json)',
        closed_forms='k = 0: Taylor coefficients of P(eps,infty) (Lemma lem:tl-k0), all W; '
                     'k = 1: Bessel I_0/I_m generating functions and the closed form (eq:tl-bubbleclosed) in the derivative '
                     'numbers (a,b,j), all W; k = 2: I_1/I_{m+1} with Dirichlet integrals on the 2-simplex; '
                     'general k: eq:tl-genfun with k+1 edges; k <= w_ab + d_ab',
        derivative_generating_function='letters point-split phi_a(x+y_a) phi_b(x+y_b); factor e^{i(y_ab,l)} e^{i(y_b,R)}, y_ab = y_b-y_a; '
                                       'kernel on d^alpha phi_a d^beta phi_b = alpha! beta! [y_a^alpha y_b^beta]'),
    tasks=[
        dict(id='T1', what='K^log and ct_m for every pair type and every chain, W <= 6, general colour position '
                           '(colour element (1 (x) w)(-,-)^{-1}) and planar adjacent form', output='kernels_W6.json'),
        dict(id='T2', what='lc^natural on gauge-invariant families: scalar bilinears (gate), length-2 sl(2) twist-two '
                           'tr(D^k Z D^{n-k} Z) for n <= 8, su(2) Tr(Z^{L-M}X^M) for L <= 6, Konishi multiplet members '
                           'of weight <= 4 with fermions and field strengths', output='spectra.json'),
        dict(id='T3', what='composite two-letter kernels C2_Delta, C_{P_b} and the one-body C1_Delta beyond the planar scalar '
                           'sector (where Prop prop:tl-scalarsector gives C2_Delta = -u v.(X_a X_b), C2_Pb = -u d_b(X_a X_b)); '
                           'method: special-edge evaluation of Thm c9-D1 (cycles through the pair with one special edge and the '
                           'composite line or S^(1)(Delta) on a letter edge), checked on the bilinears against the table of '
                           'Thm c9-D1; composite leg sets in census_composite.json; first targets: (d^a Z, d^b Z), (A_V, Z), '
                           '(lambda_+, lambda_-)', output='composite_kernels.json', status='not computed'),
        dict(id='T4', what='Ward gates W1 (Prop prop:tl-W1, every letter type) and W2 (Cor cor:tl-ward) on all pairs '
                           'with W <= 4', output='ward_gates.json', status='W1 done (both engines); W2 done in engine A on field-letter pairs (G3-W2)'),
        dict(id='T5', what='W2 with s on antifield outputs and antifield letters: (c,A*_V), (A_V,lambda_+-), derivative (lambda,c), (A_V,A_V) at W_0+1 with three outputs'),
        dict(id='T6', what='Conjecture P7 (conj:tl-P7) on the families of T2/T3: twist two (d^a Z, d^b Z), (A_V, Z), a fermion pair, a field-strength family'),
        dict(id='T7', what='two-body kernels K^(x) of Q_1, S_1, K_1 (Prop prop:tl-arity1): census with edges P, one special edge K_eps-K_0 or composite line, vertices of I and at most one S^(1)(x); bound Prop prop:tl-bound1')],
    routes=dict(
        A=dict(name='Schwinger-parameter Gaussian generating functions with exact simplex integrals',
               method=['chain numerator N(l) from Hessians at plane-wave backgrounds (engine A conventions)',
                       'Gaussian loop integral: (16 pi^2 tau^2)^-1 sum_j num_j tau^-j, num_j = (Lap^j N)(c)/(j! 4^j)',
                       'eps^-m: (1/m) sum_e int_{R_e} alpha_e^m F_m; log: - int_Delta F_0 (Thm thm:tl-simplex), exact '
                       'Dirichlet integrals over regions R_e (cone from barycentre, volume 1/E)',
                       'letter derivatives through the Bessel generating functions (Prop prop:tl-genfun) or directly as '
                       'polynomial numerators',
                       'special edges: Lemma lem:tl-special (no log; finite part -int_{R_s} F_0 log(alpha_s/min_{f!=s} alpha_f))'],
               library='step2/twoletter/theory/tl_simplex.py (cycle_singular, cycle_singular_totaltime, '
                       'special_edge_singular, special_edge_finite, F_m, genfun_F, bubble_log_closed)',
               reuse=['step2/engine_a/s2_core.py: Pnum (propagator numerators), Khat (K_0), canon/V (vertex functions, '
                      'cubic x2, quartic x8 scale), hessian, strace (verified: checks_engines_results.json A1-A4)',
                      'engine_a/cliff.py (10d Clifford blocks), engine_a/r2_core.py']),
        B=dict(name='explicit momentum-space numerator algebra with Feynman parameters and exact one-dimensional integrals, '
                    'plus numerical eps-extrapolation at sample points',
               method=['full numerator in momentum space (sympy), edges e^{-eps q^2}/q^2',
                       'singular part from the UV boundary: eps d/deps = - sum_e (edge e pinned at t = eps), pinned graphs '
                       'Taylor-expanded in the external momenta with t_f = eps u_f, u_f in [1,oo), exact one-dimensional '
                       'u-integrals after Feynman parametrisation of the remaining edges',
                       'numerical check at 2-3 sample momenta: residual r(eps) = eps dw/deps + sum k a_k eps^-k - b at '
                       'eps = 1e-2..1e-4 (40 digits for E=2, 25 for E=3), required to decrease like eps log eps'],
               reuse=['engine_b/vertices.py (Lagrangian terms, vertex), engine_b/loops.py (derive_propagators: inversion '
                      'of the free Hessian, prop_numerator, gaussian_moments, boundary_singular)',
                      'step2/engine_b/fa.py (functionals on plane waves, contraction B_K, K_0, special edge, Q)',
                      'verified: checks_engines_results.json B1 (A = 2/q^2, lambda, ghost -i q_m)']),
        independence='A and B do not read each other outputs before the comparison phase; A uses the simplex formula, B the '
                     'pinned-edge Taylor expansion and numerical eps residuals'),
    gates=[
        'G1 scalar bilinears: lc = 4 O (traceless), -8 O (trace); graph table of Prop c5-table; eps^-1 vacuum term',
        'G2 Callan-Symanzik values: lc^natural = 0, -12 O; letters -2,-2,-4,-1 units (Thm c6-dil)',
        'G3 Ward W1 (Prop prop:tl-W1) for A_perp, lambda, c, A_V, A*_perp, A*_V, lambda*, c* letters; G3-W2: W2 (Rem rem:tl-W2gate)',
        'G4 colour: su(2), su(3), gl(3) realisations give the same rational multiples; adjacency rule (checks_colour.py)',
        'G5 shape: log eps identical for per-edge and total-time cutoffs on every computed cycle; eps^-m differ',
        'G6 k = 0: pure powers eps^{-(w_a+w_b)/2}, no log, no finite part (Lemma lem:tl-k0)',
        'G7 output bound n_out <= W + degree and weights W - 2m (Prop prop:tl-bound) on every kernel',
        'G8 spectrum: eig(lc^natural) = -2 E_H units on planar single traces (predictions.json P5, P6)',
        'G9 D_1 on bilinears and planar scalar sector = 4 lc^natural - (1/3+theta_perp){S1(Delta),O} (Prop prop:tl-scalarsector); P_1 = -(1/3+theta_P){S1(P_b),O}',
        'G8f su(1|1) two-site density: lc^natural = -2u(1 - Pi) on (Z,psi), (Z,Z), (psi,psi)',
        'G10 two-route agreement: exact rationals equal between routes A and B for every kernel coefficient'],
    machine=dict(memory='each process < 1.5 GB RSS, one heavy process at a time', tmp='never /tmp; write under step2/twoletter/',
                 language='pure Python (fractions, sympy, numpy, mpmath) or C; cluster-portable; no Wolfram'),
    checks_done=dict(census='checks_census.py -> census.json (0 degree mismatches, 0 bound violations)',
                     census_composite='checks_census_composite.py -> census_composite.json (one composite c-A_V line, W_0)',
                     simplex='checks_simplex.py -> checks_simplex_results.json (ALL PASS)',
                     colour='checks_colour.py -> checks_colour_results.json (ALL PASS)',
                     engines='checks_engines.py -> checks_engines_results.json (all_pass)',
                     predictions='build_predictions.py -> predictions.json'),
)
with open(os.path.join(HERE, 'engine_spec.json'), 'w') as f:
    json.dump(spec, f, indent=1)
print('engine_spec.json written; pair types:', len(kernels))
