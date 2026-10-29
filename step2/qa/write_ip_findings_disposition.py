"""Writes step2/qa/ip_findings_disposition.json: disposition of the second-round audit of the inner-product design."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
T = 'step2/theory/'
D = []
def add(i, lens, sev, loc, claim, verdict, action, evidence):
    D.append({'id': i, 'lens': lens, 'severity': sev, 'location': loc, 'claim': claim, 'verdict': verdict,
              'action': action, 'evidence': evidence})

add('AS-1', 'algebra-signs', 'major', 'ch08 def:c8-2pt with thm:c8-adjoint; design rem:ip-spin, prop:ip-tree, konishi, twenty, weight3, thm:ip-logic (i); engine_plan_ip two_point, signs_recorded, G0.1, open_items[0]; tp def:t-2pt',
    'The printed coefficient of t1 t2 is (-1)^{|O1||O2|} times the expectation; the Ward identity of thm:c8-adjoint holds for the expectation (coefficient of t2 t1); the fermionic norms are negative in the printed form; fixing the sign by positivity hides a convention error.',
    'accepted',
    'New lem:ip-ward (Grassmann order and Ward signs) and rem:ip-inputs item; prop:ip-spin computes the sign (+24 in the t2 t1 form, -24 in the t1 t2 form or with spinor factor +rho(x)/|x|); all relations use the coefficient of t2 t1. tp eq:t-def changed to d_t1 d_t2 (coefficient of t2 t1), proof of thm:t-ward reads the t2 t1 coefficient, prop:t-census-zeta sign fixed to +; two_point_census.json regenerated. engine_plan_ip.json: conventions.two_point = coefficient of t2 t1, signs_recorded without "fixed by positivity", G0.1 a genuine gate, open item removed. CH08 corrections def:c8-2pt (two entries) and thm:c8-adjoint proposed.',
    T + 'inner_product_checks.json C12: coefficient of t1 t2 = -1, of t2 t1 = +1 for odd O; Ward table (t2t1 form = hypothesis in all four cases; t1t2 form opposite for odd Y); independent realisation of Def c8-twocomp: radial O_zeta ratio +1 (t2t1 form), -1 with s_lambda = +1; separated value 24 (rho(x)xi,chi)/(16pi^4|x|^6). Derivation: E(delta(FG)) = 0 for a derivation of total parity |Y|. Agrees with step2/qa/ip_audit2_algebra-signs/a3_results.json G1, G2 and conj prop:RP.')
add('AS-2', 'algebra-signs', 'minor', 'thm:ip-P (ii)',
    'The odd-generator consequences fail when D_1 != 0: the order-hbar part of [rho(xi),rho(eta)] contains the D_1 term; on the Konishi pair r_xi + r_xidag = r_Delta - 12.',
    'accepted',
    'thm:ip-P (ii) restated: additive relations for even generators and for odd generators on D_1 = 0 families; with D_1 != 0 the dilatation part of [x,y] contributes its one-loop component; on the Konishi pair the rescalings of xi and xi^dagger add up to the rescaling of Delta minus 12 (prop:ip-konishi).',
    T + 'inner_product_checks.json C6 (MC on the Konishi pair: q + s = k/2 = -12); qa/ip_audit2_algebra-signs/a4_results.json R1; hand derivation in the proof of thm:ip-P (ii).')
add('AS-3', 'algebra-signs', 'minor', 'proof of thm:ip-P (v)',
    '[Delta,P_b] = -2 P_b, not +2.', 'accepted',
    'Proof of (v) now reads [rho(Delta),rho(P_b)] = rho([Delta,P_b]) = -2 rho(P_b); conclusion unchanged.',
    'qa/ip_audit2_algebra-signs/a1_results.json 3_psl44 (ad_Delta on P = -2); ch02; consistency with the proof of (ii).')
add('AS-4', 'algebra-signs', 'minor', 'rem:ip-radius',
    'The radial Gram for the sphere of radius R is R^{-2(w+kappa)} times the unit-sphere Gram; signs of the identifications follow.', 'accepted',
    'rem:ip-radius rewritten (radius written as a second unit of length l\', no symbol R): factor (l\'/l)^{-2(w+kappa)}; basis change 1 + c unit w <-> log(l\'/l) = -c unit; inner direction of thm:ip-laws (a) <-> -2c unit; unit-of-length change (e) <-> radius l^2/l\'.',
    T + 'inner_product_checks.json C13 (lim (R^2/x^2)^Delta c |R^2 x/x^2|^{-2 Delta} = c R^{-2 Delta}; log R = -c unit); qa/ip_audit2_algebra-signs/a3_results.json G5.')
add('AS-5', 'algebra-signs', 'minor', 'thm:ip-laws (a)',
    '"+-4a with the sign of ch01 def:c1-bracket" leaves open a sign already fixed by {S^(1)(Delta), f(0)} = 2 w f(0): the shift is +4a.', 'accepted',
    'thm:ip-laws (a): the field rescaling is normalised as phi = a tr_g(ad ad)/(16 pi^2 (X,Y)) int(A_perp,A*_perp), so hbar{phi,O} = 2 a unit O on the scalar bilinears, U = 1 + 2a unit, shift +4a (in units), with the convention stated; proof gives {S^(1)(Delta),(iota_a A,iota_b A)(0)} = 4 (iota_a A,iota_b A)(0) from the term int(2A_perp,A*_perp) of S^(1)(Delta).',
    'S^(1)(Delta) two-component display (thm:ip-laws) and {S^(1)(Delta),f(0)} = 2 w f(0) (ch02, compare_delta.json convention_alignment: {S_M,-} = +D_M); tp remark after cor:t-checks (+4a).')
add('AS-6', 'algebra-signs', 'minor', 'inner_product_checks.py C8',
    'C8 coded U = +conj(kappa_B) n_B/n_zeta and never asserted it; prop:ip-weight3 (i) states the minus sign.', 'accepted',
    'C8 now derives U from the matrices (entry (O_zeta,B) of conj(X)^T Gram = Gram Y with X_{B zeta} = -kappa_B - u c_B\', complex kappa_B) and asserts U = -conj(kappa_B) n_B/n_zeta; docstring corrected.',
    T + 'inner_product_checks.json C8: order_sqrt_unit_relation_U = -n_Bp*u*conjugate(kappa_Bc)/n_zeta, U_equals_minus_conj_kappaB_nB_over_nzeta = true.')
add('AS-7', 'algebra-signs', 'minor', 'ch08 thm:c8-L (iv), def:c8-components, thm:c8-adjoint',
    '|O| in Lambda^(1)(x)O = (-1)^{|O|(1+|x|)}{S^(1)(x),O} must be the degree (ghost number); in the Ward identity |F|, |Y| are total parities.', 'accepted',
    'rem:ip-inputs item "Gradings"; lem:ip-ward and thm:ip-logic (i) state total parities. CH08 corrections thm:c8-L (iv) (|O| = degree of Def c2-dce, written without the new symbol gh(O) proposed by the auditor) and thm:c8-adjoint (total parities).',
    'Def c2-dce: "|g| the degree of g in C^bullet"; the (1,0) cochain bracket gives exactly (-1)^{|O|(1+|x|)} with the degree; hand computation with the identity of cor:c8-homotopy (ii) on an odd O of degree 0 (total-parity reading gives +Lambda([x,y])).')
add('AS-8', 'algebra-signs', 'minor', 'prop:ip-weight3 two-component paragraph',
    'Sigmabar^{Iij} is undefined in ch08 (new notation); the formulas are right only for Sigmabar = Sigma^dagger.', 'accepted',
    'Replaced by A^{ij} = (1/2) eps^{ijkl} A_kl = A_I conj(Sigma^I_ij) and c_IJK = conj(Sigma^I_ik) Sigma^J_kj z^j Sigma^K_il z^l (also in the names paragraph); C11 recoded with conj(Sigma) and the identity (1/2) eps Sigma = conj(Sigma) asserted; flavour sum 48 unchanged.',
    T + 'inner_product_checks.json C11 (thooft_blocks true with conj(Sigma); flavour sum 48); qa/ip_audit2_algebra-signs/a1_results.json thooft_conj_equals_half_eps.')
add('QS-1', 'qft-structure', 'major', 'thm:ip-P (iii), (iv); rem:ip-inputs; rem:ip-radius; prop:ip-Kline; lem:ip-W; engine_plan_ip IP6',
    'Nothing shows that P_1 is an inner automorphism or a universal rescaling on all test families; the composite effect is a two-letter kernel (bilinears -1/3, letters 0); flavour-exchange and letter-type parts can differ; the master equation at (P_b, xi) forces one rescaling on families joined by classical xi-maps; universality is a prediction.', 'accepted',
    'thm:ip-P (iii) restated as proposed (scalar bilinears computed; MC at (P_b,xi) forces the same rescaling on families reached by classical xi-maps; P_1 = 0 on the gauge-invariant test families exactly when the composite components of O_zeta, the 20, Fbar_z, B_z equal -(1/3) unit {S^(1)(P_b),-}, a prediction; letters keep +(1/3) unit). New prop:ip-universal (colour-dipole lemma with proof; flavour/letter-type caveat; MC at (P_b,xi)). rem:ip-Pneeds lists what the test needs from P_1. engine_plan_ip.json: E11 extended to the composite A_V--c translation components on O_zeta, 20, Fbar_z, B_z; IP6 gate "equal to the bilinear value".',
    'qa/compare_delta.json comparisons_part2 [9] and manuscript_and_design_predictions (B reports no such graphs for letters); delta_theory prop:d-P (ii); qa/ip_audit2_qft-structure/p1_universality.json A-D; proof of prop:ip-universal (i): on a singlet sum_j ad^{(j)} = 0, so sum_{j!=i} ad^{(i)}ad^{(j)} = -(ad^{(i)})^2 = -tr_g(ad ad)/(X,Y).')
add('QS-2', 'qft-structure', 'major', 'thm:ip-P (v); Status; engine_plan_ip IP6 and prerequisites',
    'The -1/3 multiplies the scalar-sector part S^(1)(Delta)|_perp, coefficient f_S = z + offset (cubic arity-one equation); r_Delta = 0 must hold separately on the other families.', 'accepted',
    'thm:ip-P (v) restated in terms of the coefficient of S^(1)(Delta)|_{W_R^perp} in hbar F_1^(1)(Delta) (= coefficient of the whole S^(1)(Delta) plus the offset of the cubic arity-one equation), prediction -1/3 units, and vanishing of the total coefficient of {S^(1)(Delta),-} on the other test families; z fixed by (K_b,P_b\') (delta_theory prop:d-z), (Delta,P_b) constrains the same data. Status and IP6 updated; prerequisites: cubic arity-one equation along Delta, R_O on O_zeta, Fbar_z, B_z.',
    'qa/compare_delta.json comparisons_part1 [18] (free directions S_M|_perp, S_M|_lambda, S_M|_V + S_M|_c) and part2 [7], [8] (exact for f_S = -1/3); delta_theory ex:d-quadratic (cubic equation fixes S^(1)(Delta)|_perp), eq:d-D1 (R_O), prop:d-z.')
add('QS-3', 'qft-structure', 'major', 'thm:ip-laws (f); covariance sentence; summary',
    'Covariance under a change of gauge fixing is asserted, not proved; in the minimal complex no smooth family joins the Feynman operator ((Q^GF)^2 != 0) to a nilpotent gauge fixing with a Laplacian.', 'accepted',
    'thm:ip-laws (f) restated as a hypothesis (wording of the finding); covariance "(a)-(e), and (f) under its hypothesis"; rem:ip-inputs item "Gauge fixing"; tp prop:t-laws (gauge) and cor:t-inv (i) marked as hypothesis; engine_plan_ip.json new stage IP8 (composite components in a nilpotent gauge fixing of the enlarged complex, or a proof that the Feynman-gauge composite effects are of types (a)-(c)).',
    'odd_components_design prop:o-laws (c) is a proof plan; delta_theory rem:d-gf and thm:d-inv (hypothesis); engine_plan.json open item "gauge-fixing law (c): proof plan only".')
add('QS-4', 'qft-structure', 'minor', 'ch08 cor:c8-gram',
    'On a family mixing lengths the hbar-expansion puts the one-loop length-decreasing component into Y_0; report as a Chapter 8 issue.', 'accepted',
    'CH08 corrections cor:c8-gram (statement: "of one length"; proof: graded expansion appended). thm:ip-graded states that the hbar-expansion of cor:c8-gram applies on a family of one length.',
    'qa/ip_audit2_qft-structure/ch08_gram_mixed_length.json; ' + T + 'inner_product_checks.json C1.')
add('QS-5', 'qft-structure', 'minor', 'prop:ip-orders; engine_plan_ip IP3; tp prop:t-census-w3',
    'The order-one-half B_w-Fbar_w/F_n graphs have one cubic vertex (8 half-edges, even) and vanish by colour, not by an odd half-edge count.', 'accepted',
    'prop:ip-orders statement and proof: one cubic vertex, 3+2+3 = 8 half-edges; vanishing by colour (delta_ac f_abc = 0 with the B self-contraction) or by field content (no cubic vertex (A_perp,A_perp,lambda)); odd counts 3+2+3+3, 3+2+4 at order one. tp prop:t-census-w3, two_point_census.json (builder edited and rerun) and engine_plan_ip IP0/IP3 corrected; new gate G0.7.',
    T + 'inner_product_checks.json C14 (enumeration over all cubic vertex types: 6 connected graphs for the Yukawa vertex of the right chirality, all with a self-contraction; surviving 0); qa/ip_audit2_qft-structure/bf_order_census.json.')
add('QS-6', 'qft-structure', 'minor', 'thm:ip-logic (i); tp prop:t-K',
    'The K Ward identity is invoked under a hypothesis (existence of the arity-two graphs) that fails at T = infinity; the proof needs only the T-independent local representative L_1^(2)(eta1,eta2) at finite T.', 'accepted',
    'thm:ip-logic (i) and thm:ip-origin (ii) cite the finite-endpoint local representative (k_transfer thm:ktr-defect (v)) and state that the arity-two flows grow as log T when Gamma(eta1,eta2)|_V != 0; tp prop:t-K hypothesis replaced accordingly.',
    T + 'k_transfer.json T_to_infinity: safe includes L_1^(2)(eta1,eta2); not_safe includes W^(2), O[T]^(2) along (eta1,eta2) with v1 != 0; k_transfer thm:ktr-defect (v).')
add('QS-7', 'qft-structure', 'minor', 'cor:ip-D; test table T5; engine_plan_ip IP5',
    'On the chosen families T5 (D_1 = D_1^adj) is automatic at weights 2, 5/2 and equivalent to T3 (iii) at weight 3.', 'accepted',
    'cor:ip-D adds the sentence of the finding; T5 moved to a family with one-loop mixing within one length (new prop:ip-mixing: dimension-four 20\' of [BERS], weight-four singlets) and to IP7; table marks the old T5 as a consequence.',
    T + 'inner_product_checks.json C17 (Gram_tree D_1 Hermitian at leading order exactly at the orthogonality solution); qa/ip_audit2_qft-structure/t5_redundancy.json.')
add('QS-8', 'qft-structure', 'minor', 'prop:ip-tree; rem:ip-spin; rem:ip-bps',
    'T0 and T3 (iv) check moduli only and T3 (iv) uses the radial adjointness at order one half as input.', 'accepted',
    'prop:ip-tree (b) "passes in modulus", sign gates listed in (c) (positivity computed, not chosen; sign of the classical MC identity); rem:ip-bps wording replaced as proposed ("a check of the moduli ..., given the radial adjointness at order one half and the master equation at (xi,xi^dagger); signs are tested by the gates of prop:ip-tree (c) and IP0"); Status says "in modulus".',
    'qa/ip_audit2_qft-structure/t3_tree_inputs.json (kappa_B = -1/4, flavour 48, factor 6; only moduli enter k); C11.')
add('QS-9', 'qft-structure', 'minor', 'rem:ip-twobody; engine_plan_ip IP7',
    'Where Gram_{1/2} != 0 the order-one relation contains products of two order-one-half structures, i.e. three- and four-letter components.', 'accepted',
    'Last sentence of rem:ip-twobody replaced as proposed; prop:ip-orders notes the nonzero order-one-half entries at weights 5/2 and 7/2; IP7 includes the three- and four-letter components.',
    'eq:ip-graded terms [X_{1/2}^adj, M_{1/2}] and M_{1/2}[Y_0, M_{1/2}]; tp prop:t-orderg (O_zeta with the trilinear of {S^(1)(xi),F_n}).')
add('TP-1', 'targets-plan', 'major', 'prop:ip-konishi, prop:ip-twenty, thm:ip-logic (ii),(iv), rem:ip-tested, thm:ip-laws (d); engine_plan_ip IP5 and prerequisites; engine_plan.json E3',
    'The real parts of the T1, T2 radial lines are not decisive: the free multiple of S^(1)(xi) (G_B restriction, inner R_xi) shifts them by 2 Re; only reflection compatibility (S_+ (x) W <-> S_- (x) W^*) fixes it, and nothing computes along S_- (x) W^*.', 'accepted',
    'New prop:ip-realpart (free multiple, shifts of the lines, reflection compatibility + MC at (xi,xi^{dagger_n}) fix the real part: universal real rescaling of the supercharges = half the translation rescaling; consequence for the test). New prop:ip-reflline (T1n, T2n: reflection lines along xi^{dagger_n} = rho(n) conj(xi) in S_- (x) W^* with Q_1 along S_- (x) W^* (E4n, E5n with the E12 machinery) and the separated-point two-point functions). Alternative recorded (E3: counterterm along S_- (x) W^* as the Theta_n-image, then the arity-two equation at (xi,xi^{dagger_n})). Without either: T2 radial line fixes Re and only T1 - T2 and T3 test (stated in prop:ip-realpart (iv), prop:ip-konishi, table, IP5). CH08 correction of the last sentence of cor:c8-gram.',
    T + 'inner_product_checks.json C16 (real reflection-compatible rescalings: s3 = (s1+s2)/2, r_Q = r_P/2; real inner Delta allowed, real B excluded; imaginary: r_P = 0, imaginary B); qa/ip_audit2_targets-plan/ip_audit2_checks.json A1, A2; engine_plan.json E3 and open item on Re a; odd rem:o-alpha, rem:o-GBdagger.')
add('TP-2', 'targets-plan', 'major', 'engine_plan_ip IP0, IP5, IP6; lem:ip-W, thm:ip-origin, thm:ip-P (iii), prop:ip-konishi',
    'Frame bookkeeping inconsistent: the radial Gram equals the two-point coefficient only in the frame P_1 = 0; the inner redefinition leaves the two-point function unchanged; mixing frames gives a residual of G_B pattern; IP6 double-counts (multiplying G1 by 1 + 2c unit w); routes (a), (b) are not interchangeable.', 'accepted',
    'New prop:ip-frames: (a) inner direction: matrices U^{-1}XU, two-point unchanged, radial Gram of the engines\' frame = (1 - 2c unit w) x two-point; (b) finite part c x weight: same matrices on jets at the origin, two-point x (1 + c unit (w1+w2)), hypothesis of lem:ip-W not established, not used for thm:ip-origin; (c) conversion rule and the mixed-frame residual 4 c_1 - 5 c_2 (-c universal). thm:ip-P (iii) keeps route (a) for the identification and states that (b) changes the two-point function; thm:ip-laws (a) explains why the radial Gram changes although the two-point function does not. engine_plan_ip.json: IP6 runs before IP5, computes per-family remainders, converts matrices, keeps G0, G1 = two-point coefficients; the "G1 multiplied by 1 + 2c unit w" step deleted (conventions.frame warns against it); test set states that all data are in the frame P_1 = 0.',
    T + 'inner_product_checks.json C15 (sl(2) model: frame-E relation exact with Gram_E = (1 - 2 c unit w) x two-point; mixed frames fail; previous IP6 rule leaves residual 4 c w = 2 c K_0 at entry (0,1); Konishi pair residual 4 c_K - 5 c_zeta, -c universal); qa/ip_audit2_targets-plan/ip_audit2_checks.json A1, A3.')
add('TP-3', 'targets-plan', 'major', 'engine_plan_ip IP1 method, IP2, IP3, independence clauses; design Engine plan section',
    'The engines are independent in code but not in method: both prescribed engine A\'s Beta/Dirichlet route and a shared census; the finite-eps quadrature confirms only to O(eps log eps).', 'accepted',
    'engine_plan_ip.json: constraints.independence requires code AND method independence and own census; IP1 method split (engine A position-space Beta/Dirichlet; engine B momentum-space whole-graph Schwinger-Gaussian with its masters.py and the Fourier lemma, or several-eps extrapolation; engine B re-derives lemmas t-reduce, t-local, t-ft and the five-line contact claim); numeric confirmation described as O(eps log eps); IP2, IP3 "as IP1 with the engine-specific routes"; comparison protocol compares censuses graph by graph. Design section "Engine plan" updated.',
    'engine_b/method_notes.json (engine B used closed-form masters J3/J2 and finite-eps whole-graph numerics in Parts C/D); two_point_example_run.txt (J ratio 0.999912 at eps = 1e-4); grep: no cross-imports.')
add('TP-4', 'targets-plan', 'minor', 'thm:ip-P (v), Status, engine_plan_ip IP6',
    'The engine input is f_S (coefficient of S_M|_perp), not z; cite (K_b,P_b\') for z (delta_theory prop:d-z); add the Delta composite remainders of O_zeta, 20, Fbar_w to IP6.', 'accepted',
    'Merged with QS-2: (v) in terms of the S^(1)(Delta)|_{W_R^perp} coefficient (and the S^(1)(Delta)|_lambda, S^(1)(Delta)|_V + S^(1)(Delta)|_c coefficients for the families with lambda letters); z cited from (K_b,P_b\'); IP6 tasks include the Delta composite remainders; the use of (Delta,P_b) kept as a second constraint (bracket term +-2 S^(1)(P_b) != 0, as noted by QS-2).',
    'qa/compare_delta.json part1 [18], part2 [7]; delta_theory prop:d-z, eq:d-Fsol.')
add('TP-5', 'targets-plan', 'minor', 'test-set table, cor:ip-D, prop:ip-Kline; engine_plan_ip IP3 outputs, IP5 covariance runs',
    'The test set is not minimal: T5 redundant, MC lines of T1/T2 repeat E9, T4 is the kinematic identity and the MC at (K,P), m_F enters no relation, covariance runs cannot fail.', 'accepted',
    'Table column "status" marks gates (T0 (a), the T1/T2 gates = E9, T4) and decisive tests; T5 moved to a mixing family (prop:ip-mixing); m_F recorded as a by-product in IP3 (enters no relation at relative order unit); IP5 covariance runs split into bookkeeping runs (cannot fail, stated) and one genuine run with the finite part inserted in the IP1/IP2 and E5/E8 graphs; rem:ip-tested lists the decisive content.',
    'qa/ip_audit2_targets-plan/ip_audit2_checks.json A5; ' + T + 'inner_product_checks.json C17; ex:ip-scalars (D_1 diagonal and real); engine_plan.json E9.')
add('TP-6', 'targets-plan', 'minor', 'prop:ip-orders; tp prop:t-census-w3; engine_plan_ip IP3',
    'The order-one-half class is the one-cubic-vertex tree (8 half-edges) and vanishes by field content; IP3 records it as zero without checking.', 'accepted',
    'Merged with QS-5: the proof gives both mechanisms (the Yukawa topology vanishes by colour, the topology with two scalar legs at the vertex by field content); IP0/IP3 gate G0.7 checks the one-cubic-vertex trees graph by graph.',
    T + 'inner_product_checks.json C14; qa/ip_audit2_targets-plan/ip_audit2_checks.json A4.')
add('TP-7', 'targets-plan', 'minor', 'literature table; ip_targets Target 11 and check_ip_targets.py C7; engine_plan_ip IP7',
    'Five target records wrong or incomplete: [Ryz] ratio no cross-check; T1 protection from [Mi97],[DO02]; BKRS anomaly/orthogonality are statements; g^2 in BKRS (27) = 4 g_std^2; hep-th/0307015v5 regulator missing.', 'accepted',
    '(a) [Ryz] row: statement only, fails in the UV-boundary scheme, no cross-check (design table; ip_targets.tex Target 10 project object; build_ip_targets_json.py IP10 and open item resolved). (b) T1 cites [Mi97], [DO02]. (c) BKRS row: coefficient quoted from their ref. [8], orthogonality from the U_B(1) bonus symmetry (statements) in the design table and ip_targets Target 12 regulator field (tex and json). (d) C7 annotation and Target 11 text: g^2 in (27) = 4 g_std^2 (footnote 3), check_ip_targets.py rerun. (e) [B03c] hep-th/0307015v5 added with inputs and regulator (point splitting, Sec. 3, App. B as corrected in v5; n4_untwisted B03c), statement only, in the design table, IP7 and the ip_targets.json open items.',
    'qa/ip_audit2_targets-plan/lit_spotcheck.json; ip_sources/hep-th_0104016/paper.txt lines 196 (footnote 3: "replace g there, with 2g"), 223, 262-263 ("as expected on the basis of the so-called U_B(1) bonus symmetry [21]"), 403-425 ((24), (27)); A6 (m_20\' = 27/2 + 4 log(3/2), log coefficient 0); n4_untwisted.json B03c (regulator point splitting).')
add('TP-8', 'targets-plan', 'minor', 'engine_plan_ip cluster_portability, IP1 method, IP7',
    'The prescribed code is not cluster-portable: 11 step2/engine_a files hard-code absolute paths; IP7 has no truncation or RSS estimate.', 'accepted',
    'engine_plan_ip.json: file-relative paths, an explicit module bundle and a copied-directory gate required; the 11 files are listed as a known defect to be fixed before cluster runs (engine code itself not edited in this design task); IP7 gets a per-job truncation (one job per pair of letter types and derivative order, first pass <= 4) and a pilot job measuring RSS (< 1.5 GB gate).',
    'grep: /Users/libotao/Desktop/N4_Costello_BV_2026-10-27 occurs in step2/engine_a/{d4_cycles,d4_pi_sectors,d4_signs,d4_se_sectors,d4_two_point,d4_y1a,s2_core,s2_ga,s2_trees,s2_make_results,s2_verify0d}.py; none in step2/engine_b (qa/ip_audit2_targets-plan/ip_audit2_checks.json A9).')
add('TP-9', 'targets-plan', 'minor', 'props ip-tree, konishi, twenty, weight3, Kline; thm:ip-P; rem:ip-radius',
    'The test-set formulas display symbols not defined in Chapter 8 (q, s, m_K, n_zeta, tau, sigma, kappa_B, c_K, b_F, b_B, U, V, k, r_x, s_a, z, f, R).', 'accepted',
    'The test set is restated with Chapter-8 objects only: the equations are identities between pairings (<,>_rad,tree, <,>_rad,1, <,>_n) of the operators and components Q_1, S_1 (at eta = conj xi), D_1, K_1 and the classical brackets (eq:ip-radline, eq:ip-mcline, eq:ip-konishi, eq:ip-reflline, prop:ip-weight3 (i)-(iv), eq:ip-konishiclass); the coefficients are described in words; rescalings are "the rescaling of x"; the radius is a second unit of length; kappa is that of ch08 prop:c8-primary; the W-polarisation is written z as in ch08 lem:c8-ev (B_z, Fbar_z = B_w, Fbar_w of [odd] at w = z). The code names are kept in engine_plan_ip.json "code_names" with their defining equations. Operator labels O_zeta, B_z, Fbar_z (defined by formulas, from [tp]/[odd]) are kept; the Konishi primary is written as the trace bilinear in full.',
    'DESIGN.md rule 2; DEFINITIONS_v2 item 11; grep of the revised tex for the listed symbols in displays (none outside code-name references).')
add('TP-10', 'targets-plan', 'minor', 'thm:ip-P (i), rem:ip-inputs; ch08 def:c8-radial, thm:c8-adjoint, prop:c8-primary',
    'The design amends the P_1 = 0 hypothesis of Chapter 8 but reports no CH08 issue; def:c8-radial cites prop:c7-cov for P_1 = 0, contradicted on the scalar bilinears.', 'accepted',
    'CH08 corrections def:c8-radial (text of the finding), thm:c8-adjoint and prop:c8-primary (same qualification) in step2/qa/ch08_corrections_from_ip.json and in the structured result.',
    'qa/compare_delta.json comparisons_part2 [9]; delta_theory prop:d-P; prop:ip-frames.')
doc = {'file': 'step2/qa/ip_findings_disposition.json', 'date': '2026-10-29',
       'design_files': [T + 'inner_product_design.tex', T + 'engine_plan_ip.json', T + 'inner_product_checks.py', T + 'inner_product_checks.json',
                        T + 'two_point_design.tex', T + 'two_point_census.json (via two_point_census_build.py)',
                        'step2/literature/ip_targets.tex', 'step2/literature/ip_targets.json (via build_ip_targets_json.py)',
                        'step2/literature/check_ip_targets.py (+ output)'],
       'chapter8_corrections': 'step2/qa/ch08_corrections_from_ip.json',
       'counts': {'accepted': sum(1 for d in D if d['verdict'] == 'accepted'), 'rejected': sum(1 for d in D if d['verdict'] == 'rejected'), 'total': len(D)},
       'findings': D}
json.dump(doc, open(os.path.join(HERE, 'ip_findings_disposition.json'), 'w'), indent=1, ensure_ascii=False)
print(doc['counts'])
