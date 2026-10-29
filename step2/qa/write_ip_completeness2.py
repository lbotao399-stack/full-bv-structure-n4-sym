#!/usr/bin/env python3
"""Writes step2/qa/ip_completeness2.json (completeness critic, second pass, of the inner-product design
against DESIGN_STEP2 R0 and R5).  Every CH08 'old_text' is asserted to occur exactly once in
manuscript/ch08_definitions.tex.  Checks referenced as [K1]..[K4] are in ip_completeness2_checks.py/json."""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CH08 = open(os.path.join(ROOT, "manuscript", "ch08_definitions.tex")).read()
UNIT = r"\frac{\hbar\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{16\pi^2(X,Y)}"

items = []
def item(**kw):
    items.append(kw)

# ---------------------------------------------------------------- R5, high
item(id="IPC2-01", area="R5 (b), (c)", severity="high",
     what="The one-loop two-point coefficients that feed the radial Gram (Example ip-scalars, prop:ip-konishi, IP1-IP3) are computed with the arity-zero finite counterterm (27/16) unit int(d*A|_V, d*A|_V) only ('minimal scheme'). Through relative order one unit the two-point function of two operators of equal length depends on the two-leg part of F^(0) only (a cubic order-hbar vertex cannot join two equal-length operators: 2L+3 half-edges; a vertex with four or more legs adds a cycle, relative order unit^2), but that two-leg part is not a free choice: the cubic obstruction (ch09 thm:c9-ob3, alpha_1 = alpha_3 = 11/4 + 19/4 log(3/2), alpha_5 = -9 log(3/2)) ties the quadratic beta_k to the cubic h_j (thm:c9-F3), and the classes (beta_2, h_4) = (1/4, 1) and (beta_3, h_5) = (-1/2, 1) of C/V (odd prop:o-arity0 (iii)) are fixed only by the arity-one equation along xi (E3, two legs) and E3b. The minimal choice is the representative beta = 0 with cubic part -(11/8 + 19/8 log 3/2) H1 - (11/2 + 19/2 log 3/2) H4 - 9 log(3/2) H5 (admissible at arity zero); whether it is the arity-zero part of an equivariant solution is undecided.",
     why_it_matters="beta_2 shifts the one-loop radial constants of the trace bilinear and of the 20' by -8 beta_2 (tp remark 'Scheme dependence of the result'), those of O_zeta and of the 20 by -4 beta_2 (one A|_Wperp line) plus a beta_3 term (lambda line). The shifts move the T1 and T2 radial lines by equal amounts, the pattern of the real part of prop:ip-realpart, and leave T1 - T2 and T3 unchanged: a wrong F^(0) violates hypothesis (H1) of the Ward identities and would be absorbed into the 'fixed' real part instead of being detected. The same F^(0) must be used in E1, E3, E4, E5, E8 and in IP1-IP3.",
     next_action="(1) In engine_plan_ip.json conventions_to_fix_in_code.two_point replace 'F^(0) = (27/16) unit int(d*A|_V,d*A|_V)' by 'the two-leg part of the arity-zero part of the E3/E3b solution F_psl(4|4) (beta_1..beta_4 recorded with their C/V class)'. (2) Until E3 is done, carry beta_2, beta_3 symbolically in IP1, IP2, IP3 (they enter linearly through two-leg insertions on the A|_Wperp and lambda lines) and add to IP5 the beta-dependence of every line (T1, T2, T1n, T2n, T3). (3) In inner_product_design.tex Example ip-scalars and prop:ip-konishi state the beta-dependence of the constants.",
     two_component="beta_2 int(dA|_Wperp, dA|_Wperp) = (beta_2/4) int d^4x (d_m A^{ij}, d_m A_{ij}); beta_3 int(lambda, dslash lambda) = 2i beta_3 int d^4x (lambdabar_{adot i}, sigmabar^{m adot a} d_m lambda_a^i) (ch09 rem:c9-twocomp3); the forced term int(d*A|_V, d*A|_V) = int d^4x (d_m A_m, d_n A_n).",
     evidence=["manuscript/ch09_cubic_and_dilatation.tex thm:c9-ob3, thm:c9-F3, eq:c9-F2rep", "step2/theory/odd_components_design.tex prop:o-arity0", "step2/theory/two_point_design.tex prop:t-exist (classes), prop:t-laws (F0), remark after cor:t-checks", "step2/theory/engine_plan_ip.json conventions_to_fix_in_code.two_point", "step2/theory/inner_product_design.tex ex:ip-scalars"])

item(id="IPC2-02", area="R5 (c); R0 (sign rules on fermionic operators)", severity="high",
     what="Two QA files propose incompatible corrections of the same formula of ch08 thm:c8-L (iv), Lambda^(1)(x)O = (-1)^{|O|(1+|x|)}{S^(1)(x),O}: step2/qa/ch08_corrections_from_ip.json (AS-7) reads |O| as the degree (ghost number), step2/qa/ch08_corrections_from_design.json (thm:c8-L entry and the outside-chapter entry on ch02 def:c2-dce) reads |O| (and |Lambda'|) as the total parity. With the total-parity reading, Lambda^(1)(x) = (-1)^F {S^(1)(x),-} for even x, and x -> -Lambda^(1)(x) is not a representation on odd operators (the graded commutator gives +Lambda^(1)([x,y]) on odd O), against cor:c8-homotopy (ii) [K4]; the design audit's evidence (a1b, a1c) concerns products and brackets of functional-valued cochains, a different formula.",
     why_it_matters="O_zeta, the 20 and Fbar_z are odd; the T2, T3 lines and the Ward sign table (lem:ip-ward) use Lambda^(1)(x) = {S^(1)(x),-} on degree-0 functionals of either parity (engine_plan_ip conventions_to_fix_in_code.representation). Applying the design correction to ch08 would flip the sign of every classical component on fermionic operators.",
     next_action="One derivation, checked by a script: the (1,0) component of the cochain bracket of def:c2-dce with Koszul signs taken in the total parity, applied to an odd arity-zero O, and the resulting exponent; then (a) the representation property on an odd O of degree 0 (as in [K4]) and (b) the design audit's tests a1b, a1c on the same sign rule. Apply to ch08 only the reconciled text (proposed in ch08_issues CH08-7).",
     two_component="on O_zeta = int omega (lambda_a^j, A_{ji}) zeta^{a i} (odd, degree 0): Lambda^(1)(xi) O_zeta = {S^(1)(xi), O_zeta} with S^(1)(xi) of ch08 remark 'Comparison with the holomorphic twist'.",
     evidence=["step2/qa/ch08_corrections_from_ip.json (thm:c8-L (iv))", "step2/qa/ch08_corrections_from_design.json (thm:c8-L; outside_chapter_8 ch02 def:c2-dce)", "step2/qa/odd_rev_audit_algebra/a1b_product_rule.json, a1c_bracket_rule.json", "step2/qa/ip_completeness2_checks.json K4"])

item(id="IPC2-03", area="R5 (a), (c): P_1 and the frame", severity="high",
     what="The translation rescaling on families with lambda letters is presented as a prediction with ONE number f_P (coefficient of S^(1)(P_b) in hbar F_1^(1)(P_b)): thm:ip-P (iii), prop:ip-universal, IP6 gate 'c_F equal ... to the scalar-bilinear value -1/3 - f_P'. But, exactly as along Delta (compare_delta.json part1 [18]: free Q-closed directions S_M|_perp, S_M|_lambda, S_M|_V + S_M|_c), the sector parts of S^(1)(P_b) in (A|_Wperp, A*|_Wperp), (lambda, lambda*), (A|_V, A*|_V) + (c, c*) are separately Q-closed at two legs (Q int(d_b lambda, lambda*) and Q int(d_b A, A*) are total derivatives), so hbar F_1^(1)(P_b) has three free sector coefficients at two legs. A sector-dependent choice gives P_1 on O_zeta of the form (number) d_b on the A|_Wperp letter + (another number) d_b on the lambda letter, not a multiple of {S^(1)(P_b),O_zeta}. Moreover the one-letter composite kernels (composite line on a leg) are letter-type dependent along Delta (delta_theory predictions: Pi_Psi entries -2/3 int(A|_V,A*|_V), +1/3 int(lambda*,lambda), +4/3 int(c*,c); on A|_Wperp (5/3) S_M|_perp - (2/3) int(A|_Wperp, A*|_Wperp)); prop:ip-universal (i) treats only the two-letter colour-dipole kernel.",
     why_it_matters="The master equation at (P_b, xi) on the cohomology forces one translation rescaling on families joined by classical xi-maps; with three sector coefficients this is not a prediction but the equation that fixes f_P(lambda) - f_P(A|_Wperp) (with the cubic arity-one equation along P_b). The frame P_1 = 0 on O_zeta, the 20, Fbar_z, B_z, the conversion U = 1 + c_F unit w of prop:ip-frames and the T2/T3 lines depend on it; if a single f_P is imposed the IP6 gate can fail without any anomaly.",
     next_action="(1) In E11/IP6 record the three sector coefficients of S^(1)(P_b) in hbar F_1^(1)(P_b) (code names f_P_perp, f_P_lam, f_P_Vc) and the one- and two-letter composite translation kernels per letter type (A|_Wperp, lambda, A|_V). (2) Add the cubic arity-one equation along P_b (sectors (lambda, lambda, A|_Wperp), (A|_Wperp, A|_Wperp, A|_V), (lambda, lambda, A|_V) and the ghost sectors, complete composite sums) to the plan; it fixes the differences of the sector coefficients up to Q-exact terms. (3) Restate the IP6 gate: the MC at (P_b, xi) on (trace bilinear, O_xi) determines f_P_lam - f_P_perp; then equality of the rescalings on B_z (three letters) and Fbar_z is a genuine test. (4) Correct prop:ip-universal (i) to include one-letter kernels.",
     two_component="sector parts of S^(1)(P_b): int d^4x (b^m d_m A_n, A*_n), (1/4) int d^4x (b^m d_m A^{ij}, A*_{ij}), and the lambda sector built from (b^m d_m lambda_a^i, lambda^{*a}_i) - (b^m d_m lambdabar_{adot i}, lambdabar^{*adot i}) with the overall sign of the lambda sector of S^(1)(Delta) in the display of ip thm:ip-laws; P_1 on O_zeta = int omega (lambda_a^j, A_{ji}) zeta^{a i} of the form int omega [c_A (lambda_a^j, b^m d_m A_{ji}) + c_lam (b^m d_m lambda_a^j, A_{ji})] zeta^{a i}, a multiple of the classical translation iff c_A = c_lam.",
     evidence=["step2/qa/compare_delta.json comparisons_part1_two_point[18], comparisons_part2_bilinear_component[9]", "step2/theory/delta_theory.tex prop:d-P, section 'Predictions'", "step2/theory/delta_predictions.json", "step2/theory/inner_product_design.tex thm:ip-P (iii), prop:ip-universal, prop:ip-frames", "step2/theory/engine_plan_ip.json IP6 gates"])

item(id="IPC2-04", area="R5 (a), (c): CH08 frame bookkeeping", severity="high",
     what="The reality condition of the rescaling that removes P_1 is stated oppositely in the QA record: k_judgement.tex cor:kj-gram and k_judgement.json open issue 2 and its proposed CH08 text for def:c8-radial say the finite part (a/2){S^(1)(Delta),-} is reflection compatible exactly when conj(a) = -a; ch08_corrections_from_design.json (outside-chapter entry), ch09 prop:c9-scheme, ip thm:ip-laws (a), prop:ip-realpart (iii) and checks C16 say: for real a. [K2] confirms real: {S^(1)(Delta),-} acts on smeared densities as the real operator 2(x.d + w), which commutes with Theta_n; the sign algebra of thm:c8-fields (i), (ii) with the ch01 symmetry {g,f} = (-1)^{|f||g|}{f,g} gives {S^(1)(Delta),O}^{Theta_n} = {S^(1)(Delta),O^{Theta_n}} for both parities of O. Three versions of the def:c8-radial correction exist (ip: keep P_1 = 0 as a frame condition; design: drop it, finite part with projection pi_F, real p; k_judgement: conj(a) = -a).",
     why_it_matters="With the k_judgement version the frame P_1 = 0 would be unreachable in a reflection-compatible scheme and cor:c8-gram would carry an extra term; the T1, T2 lines and the K line T4 would be assembled in an inconsistent frame.",
     next_action="Withdraw the k_judgement/k_transfer statements (k_transfer.tex prop:ktr-rescale (ii), k_judgement cor:kj-gram, open issue 2, its def:c8-radial proposal); adopt for ch08 the design text for def:c8-radial (real p) together with the frame qualifiers of CH08-3 and CH08-4 below.",
     evidence=["step2/qa/k_judgement.json open_issues[1], ch08_issues[3]", "step2/qa/ch08_corrections_from_design.json outside_chapter_8[3], entries def:c8-radial", "manuscript/ch09_cubic_and_dilatation.tex prop:c9-scheme", "step2/theory/inner_product_checks.json C16", "step2/qa/ip_completeness2_checks.json K2"])

item(id="IPC2-05", area="R5 (c): K line", severity="high",
     what="The proposed CH08 text for the third relation of cor:c8-gram with P_1 = p{S^(1)(P_b),-} (ch08_corrections_from_design.json, k_judgement.json), 'K_1 along K_bbar = pbar{S^(1)(K_bbar),-} + hbar[{S^(1)(K_bbar),-}, Gram_tree^{-1} Gram_1]', is correct only if Gram is the radial Gram of def:c8-inner in that scheme, which is (1 - 2pw) times the two-point coefficient on Weyl weight w (prop:ip-frames (a)). With the two-point coefficients inserted (the engines' natural data) the rescaling term is -p{S^(1)(K_bbar),-}, since [{S^(1)(K_b),-}, weight] = {S^(1)(K_b),-} [K3, checks C15].",
     why_it_matters="IP5/T4 compare the K-route K_1 (E10c) with the K line; mixing the frame formula with two-point Gram data gives a residual 2p{S^(1)(K_bbar),-}, of the size of the effect being tested (-1/3 units on the bilinears).",
     next_action="Apply CH08-4 (below); in engine_plan_ip IP5/T4 state which Gram enters the K line.",
     two_component="{S^(1)(K_b),-} on a scalar density of weight w at x: v^m d_m + (w/4) d_m v^m with v^m = x^2 b^m - 2(b.x) x^m, d_m v^m = -8 b.x.",
     evidence=["step2/theory/inner_product_design.tex thm:ip-P (iv), prop:ip-frames", "step2/theory/inner_product_checks.json C15", "step2/qa/ip_completeness2_checks.json K3"])

item(id="IPC2-06", area="R5 (b), (c): gauge invariance of the renormalised test operators", severity="high",
     what="L_1(O)^(0) = 0 (one-loop BRST closure) on the test operators is assumed everywhere (tp thm:t-ward (H2), ch08 prop:c8-2pt (i), ip thm:ip-logic) and never computed. The cutoff breaks gauge invariance and the operator counterterms have gauge-variant output sectors that no census lists: for the trace bilinear the outputs (A|_V, A|_V) and (c, A*|_V) (compare_delta.json open[2], ch09 rem:c9-next); for O_zeta the weight-5/2 Lorentz-, SU(4)- and colour-covariant gauge-variant structure int omega (lambda, rho(A|_V) zeta) (superficial degree 0: log eps and finite parts; triangles O_zeta - (lambda lambda A|_V) - Yukawa and O_zeta - (A|_V A|_Wperp A|_Wperp) - Yukawa); analogous structures for the 20, Fbar_z, B_z.",
     why_it_matters="(H2) is what makes the two-point function descend to the cohomology of Lambda^(0) and the Ward identities hold; the finite parts psi that restore L_1^(0) = 0 change Q_1, S_1 (E1, E5, E8) and must be reflection compatible. These sectors do not enter IP1/IP2 at relative order one unit (P joins A|_V only to A|_V), so the omission is invisible in the two-point engines and appears only in the components.",
     next_action="New stage (IP0b, run with E5): for each test family compute all output sectors of partial_delta I^CT[eps](O) allowed by weight, Lorentz, SU(4), colour (log eps and the eps^0 part of the arity-zero defect), hence hbar L_1(O)^(0); choose psi with L_1^(0) = 0 and check psi(O^{Theta_n}) = psi(O)^{Theta_n}; record psi for E1, E5, E8.",
     two_component="int omega (lambda, rho(A|_V) zeta) = i int omega (lambdabar_{adot i}, A_m) sigmabar^{m adot a} zeta_a^i for zeta = (zeta_a^i, 0), from (lambda, eta) = lambda_a^i eta^a_i - lambdabar_{adot i} etabar^{adot i} and (rho(v)xi)^{adot i} = -i v_m sigmabar^{m adot a} xi_a^i (ch08 eq:c8-pairing).",
     evidence=["step2/theory/two_point_design.tex prop:t-census-zeta (OC lists legs (lambda, A|_Wperp) only)", "step2/qa/compare_delta.json open[2]", "manuscript/ch09_cubic_and_dilatation.tex rem:c9-next"])

item(id="IPC2-07", area="R0, R5 (c): arity-two equivariant QME", severity="high",
     what="No arity-two obstruction Ob^(2)(x1,x2) or arity-two finite counterterm hbar F_1^(2) has been computed along any pair. The test uses its consequences: the K Ward identity through (eta1, eta2) (tp prop:t-K, hence D_1 = D_1^adj in ip thm:ip-origin (ii)), theta = -1/3 through (Delta, P_b) (ch09 thm:c9-kappa, ip thm:ip-P (v)), z through (K_b, P_b') (delta_theory prop:d-z), the translation-supercharge link through (xi, xi') (prop:ip-realpart (iii)), the gates through (xi, eta) (E9); k_judgement open: L_1^(2)(eta1, eta2), L_1^(2)(K_b, P_c) only shown T-independent, not computed.",
     why_it_matters="R0 asks for the structure in arities zero, one and two; every 'gate' of the test set and the self-adjointness of D_1 are conditional on it.",
     next_action="Define the arity-two defect graphs (one special edge, two S^(1) insertions or one S^(2)) and compute Ob^(2) in the two-leg sectors for (Delta, P_b), (K_b, P_b') and (eta1, eta2) at finite T, (xi, xi'), (xi, eta) at T = infinity with the orientation sum in one integrand (engine_plan.json T_infinity_safe); solve for F_1^(2); then E9.",
     evidence=["step2/qa/k_judgement.json open_issues[4]", "step2/theory/engine_plan.json E9, T_infinity_safe", "manuscript/ch09_cubic_and_dilatation.tex thm:c9-kappa hypotheses (i), (ii)"])

item(id="IPC2-08", area="R0 (scheme freedom); R5 (c)", severity="high",
     what="The classification of the degree-0, weight-4 s-closed local cochains (hypothesis (H4) of odd thm:o-twist, used in ch08 thm:c8-invariants (iii), delta_theory prop:d-z, odd prop:o-arity0 (iv)) is unproved (E3b not run). It decides whether the three classes of C/V, the Yukawa direction H8+ + H8- and the SO(6)-singlet non-commutator (A|_Wperp)^4 directions are fixed by the equivariant QME or are genuine scheme freedom.",
     why_it_matters="It is the statement 'the scheme freedom is the coupling direction, G_B and field redefinitions' on which the invariance of the spectrum, of [Q^xi_1] and of the test lines rests; it also fixes the beta-classes of IPC2-01.",
     next_action="Run E3b: for each class of C/V solve {S^(1)(x1), G^(0)} = (-1)^{|x1|}(Q + {I,-}) G^(1)(x1) for x1 = xi, eta, P_b, R in so(Wperp) with G^(1) in the arity-one primitives of degree -1, first checking the four-leg gauge completion; report the rank.",
     evidence=["step2/theory/odd_components_design.tex prop:o-arity0, thm:o-twist (H4)", "step2/theory/engine_plan.json E3 (E3b), open_items"])

item(id="IPC2-09", area="R0, R5 (c): status of the operator-side inputs", severity="high",
     what="None of the one-loop components entering the test has been computed: Q_1 (E1, E4, E5, E6), S_1 (E8), K_1 on operators (E10c), Q_1 along S_- x W^* (E4n, E5n), the arity-two (xi, eta) component (E9), P_1 and M_1 on the lambda families (E11 extended, IP6), D_1 on the lambda families from the counterterms (Z(O_zeta) of IP2). The only evaluated test content is T0 (b) and T3 (iv) in modulus.",
     why_it_matters="R5 (c) cannot be executed; R0's 'action on local operators' is open except D_1 on the scalar bilinears and P_1 on the bilinears.",
     next_action="Order of runs (dependencies found here): E0 (signs tau, sigma, kappa_B) -> E3 + E3b (fix the two-leg F^(0), F^(1)(xi)) -> IP0b (L_1^(0), psi) -> E11/IP6 with sector coefficients (IPC2-03) and the cubic arity-one equations along Delta and P_b -> IP1-IP3 with that F^(0) and E1, E4, E5, E8 with the same scheme -> arity-two (IPC2-07) and E9, E10c -> IP5; two engines independent in code and method at each stage.",
     evidence=["step2/engine_a/, step2/engine_b/ (no Q_1, S_1, K_1 result files)", "step2/theory/engine_plan.json", "step2/theory/engine_plan_ip.json prerequisites_from_other_plans"])

# ---------------------------------------------------------------- R5, medium
item(id="IPC2-10", area="R5 (b)", severity="medium",
     what="Hermiticity of the one-loop radial (and reflection) form is neither stated nor proved; only the tree Gram is gated as Hermitian positive (G0.1).",
     why_it_matters="thm:ip-logic (iv) ('the relation for X^dagger is the adjoint of that for X'), cor:ip-D and the meaning of def:c8-inner as an inner product use it; a non-Hermitian Gram_1 would signal a non-reflection-compatible scheme or a failed Ward identity along P_n - K_n.",
     next_action="Add the proposition of CH08-6 (proof: W(P,J)^{Theta_n} = W(P,J^{Theta_n}), value at the zero field conjugated, anti-multiplicativity and Theta_n^2 = 1, then the Ward identity along P_n - K_n moves exp(-pi/4 delta^hbar) to the second slot); add an engine gate 'Gram_1 Hermitian' on every family with more than one basis element (T5 family, the weight-3 block) and 'one-loop constants real' on the Schur families.",
     evidence=["manuscript/ch08_definitions.tex def:c8-inner, thm:c8-fields", "step2/theory/inner_product_design.tex thm:ip-logic (iv)"])

item(id="IPC2-11", area="R5 (c), (d): coverage of the test set", severity="medium",
     what="(a) No test family contains an A|_V letter or a field strength: the Feynman-gauge connection (not reflection positive, conj prop:RP), the forced (27/16) counterterm, the composite A|_V--c effects on A|_V letters and the ghost loops never enter a Gram_1 of the test. (b) The graded relation eq:ip-graded is never tested: Gram_1/2 = 0 on all test families (prop:ip-orders), while families of equal weight and lengths L, L+1 with nonzero order-one-half entries exist (weight-four colour-singlet scalars: length-two bilinears in d A|_Wperp, F_A, lambda d lambda against the length-three (lambda, rho([A|_Wperp, -]) lambda)-type operators through one Yukawa vertex). (c) R5 (d) 'normalisation fixed by symmetry' is not implemented: the R-current (Noether current of S^(1)(R), R in so(Wperp)) and the stress tensor have normalisations fixed by their charges, which through the radial lines along the 20' multiplet fix the 20' constant; engine_plan_ip lists the R-current only as an 'extension'.",
     why_it_matters="R5 (c) asks for every generator and a test of the independently computed G1; the gauge sector, the mixed-length structure and the protected normalisation are where the Feynman-gauge and length effects live.",
     next_action="Add IP9: the weight-3 R-current family (gauge invariant, vector, letters A|_Wperp d A|_Wperp and lambdabar lambda, with the cubic [A|_V, -] term of D) with its charge Ward identity (contact term with S^(1)(R)) and the radial line from the 20' to it; extend T5/IP7 at weight four to evaluate eq:ip-graded with Gram_1/2 != 0; optionally int(F_A, F_A) in the 20' multiplet.",
     two_component="R-current: J_m^{[IJ]} built from A_I d_m A_J - A_J d_m A_I and lambdabar_{adot i} sigmabar_m^{adot a} lambda_a^j with the generators of su(4) (Srednicki-type bilinear lambda^dagger sigmabar lambda).",
     evidence=["step2/theory/inner_product_design.tex sec:ip-tests, prop:ip-orders", "step2/theory/engine_plan_ip.json open_items (R-current extension)", "step2/theory/conjugation.tex prop:RP (A|_V not reflection positive)"])

item(id="IPC2-12", area="R5 (c): engine-plan consistency", severity="medium",
     what="engine_plan.json E11 asks to 'test the pair conjecture p = -(n-1)/3 unit on operators with n letters (two-body kernel -(1/3)(d_i + d_j))', which contradicts ip prop:ip-universal (i): a one-loop two-letter kernel carries the colour factor sum_a ad^(i)_{e_a} ad^(j)_{e_a}, whose sum over pairs on a gauge-invariant operator is -tr_g(ad_X ad_Y)/(X,Y) per letter, so the translation coefficient is independent of n.",
     why_it_matters="B_z (three letters) discriminates: -1/3 (colour dipole) against -2/3 (pair conjecture); the two plans would read the same engine output as confirmation of different predictions.",
     next_action="Replace the E11 conjecture by the colour-dipole statement with the one-letter caveat of IPC2-03; make B_z an explicit IP6 gate.",
     evidence=["step2/theory/engine_plan.json E11 tasks[0]", "step2/theory/inner_product_design.tex prop:ip-universal (i)", "step2/qa/ip_audit2_qft-structure/p1_universality.json"])

item(id="IPC2-13", area="R0, R5 (c): arity-one equations still missing", severity="medium",
     what="Along Delta the arity-one equation is solved at two legs only (the offset between f_S and z needs the cubic arity-one equation, thm:ip-P (v)); along P_b at two legs only (sector coefficients, IPC2-03); along K_b only the 1PI special-edge graphs were summed (k_judgement open[2]: one-particle-reducible graphs, -Q I^CT(K_b), -{S^(1)(K_b), I^CT}_eps not added, eps^{-1} first-derivative parts 18i b, -108i b units not shown to cancel); along xi (E3) only two legs are planned; along so(V) the composite kernel 2(v.(P*)^2 x 1)K_0 is nonzero (rotation v) and M_1 is assumed zero (prop:c8-primary, thm:c8-adjoint, eq:ip-mcline) without computation.",
     why_it_matters="Hypothesis (H1) of every Ward identity of the test and R0's arity-one structure.",
     next_action="Add to the plans: cubic arity-one equations along Delta and P_b; complete the K_b two-leg obstruction; M_1 on the scalar and spinor families (composite sums; on a spinor family check separately the orbital and spin parts, whose relative rescaling is excluded only by [M,M] = M); extend E3 to three legs.",
     evidence=["step2/qa/k_judgement.json open_issues[2]", "step2/theory/inner_product_design.tex thm:ip-P (v)", "manuscript/ch09_cubic_and_dilatation.tex rem:c9-next"])

item(id="IPC2-14", area="R5 (b): identification of the two D_1", severity="medium",
     what="thm:ip-origin uses D_1 twice: in lem:ip-W as the component of -Lambda^(1)(Delta) on the family (def:ip-family) and in eq:ip-H through ch08 prop:c8-2pt (iii), where D_1 is the scaling-route coefficient. They agree under eq:c6-step2, which the engines find holds only up to -(1/3 + f_S) unit {S^(1)(Delta),-} (compare_delta part2 [8]); a residual rescaling r of Delta leaves a term r w log(x,x) in eq:ip-origin and the origin limit does not exist.",
     why_it_matters="The existence of the radial Gram at the origin and the identification Gram = two-point coefficient silently assume f_S = -1/3 (ch09 thm:c9-kappa, conditional on the arity-two equation along (Delta, P_b)).",
     next_action="State the hypothesis in lem:ip-W/thm:ip-origin ('the coefficient of {S^(1)(Delta),-} in D_1 vanishes on the family') and add the check as an IP6 gate on every family (already listed as a prediction in thm:ip-P (v)).",
     evidence=["step2/qa/compare_delta.json comparisons_part2_bilinear_component[8]", "step2/theory/inner_product_design.tex lem:ip-W, thm:ip-origin", "manuscript/ch08_definitions.tex prop:c8-2pt (iii)"])

item(id="IPC2-15", area="R5 (b): gauge fixing", severity="medium",
     what="Covariance under a change of gauge-fixing operator is a hypothesis (thm:ip-laws (f), tp prop:t-laws (gauge)); IP8 is planned but has no construction: the nilpotent gauge fixing with Laplacian H needs an enlarged complex (delta_theory rem:d-gf), whose propagator, counterterms and composite effects (expected zero) are not set up.",
     why_it_matters="The Feynman-gauge composite effects (P_1 = -1/3 units, the Delta offset) are exactly the frame data of the test; comparing with a gauge fixing where they vanish is the cleanest check that the frame analysis and the invariant content are right.",
     next_action="Specify the enlarged complex (a trivial pair with its antifields), its Q^GF with (Q^GF)^2 = 0 and H the Laplacian on every slot, check [C07, Def 4.3.3] in the finite-dimensional momentum-space model of ch09 rem:c9-methodD, then recompute the scalar-bilinear two-point constants, P_1 and D_1 there; compare T3 and T1 - T2.",
     evidence=["step2/theory/delta_theory.tex rem:d-gf, thm:d-inv", "step2/theory/engine_plan_ip.json IP8"])

item(id="IPC2-16", area="R5: endpoint and large distance in IP6", severity="medium",
     what="IP6 evaluates composite translation and dilatation remainders at T = infinity on new families but does not require complete composite sets; single composite graphs at T = infinity have Q-exact but nonlocal images and the composite line is O(|p|^{-4}) on c* at small momentum (delta_theory prop:d-T (ii), (iii)); engine A's earlier 'nonlocal remainder' came from this.",
     why_it_matters="A partial set at T = infinity produces spurious nonlocal remainders or log-divergent cycles; the allowed procedures are the complete set of rem:d-Y1a at T = infinity or the Q-primitive decomposition at finite T.",
     next_action="In engine_plan_ip IP6 add: 'complete composite sets (delta_theory rem:d-Y1a) or finite T with eps -> 0 first; kernel 2(v.(P*)^2 x 1)K_0 (factor 2)'.",
     evidence=["step2/theory/delta_theory.tex prop:d-T, rem:d-Y1a", "step2/qa/compare_delta.json comparisons_part2_bilinear_component[4]", "step2/qa/ch08_corrections_from_delta.json (audD1 kernel factor 2)"])

item(id="IPC2-17", area="R5 (a): choice of the reflection hyperplane", severity="low",
     what="The radial conjugation of def:c8-radial is built with a unit vector n; its independence of n at order hbar (modulo Lambda^(0)-exact terms) is not stated; it follows from the rotation covariance of the regulator only if the one-loop rotation components vanish on the family (M_1 = 0).",
     why_it_matters="The radial Gram must not depend on an auxiliary choice; a dependence would reveal a nonzero M_1.",
     next_action="State the lemma (finite rotations act exactly on the regularised theory; with M_1 = 0 the transports for n and Rn are conjugate) and add a gate: radial Gram of one family computed with two n.",
     evidence=["manuscript/ch08_definitions.tex def:c8-radial, prop:c8-C"])

item(id="IPC2-18", area="R0: arity three and higher", severity="low",
     what="The vanishing of one-loop components of arity >= 3 is proved only along constant xi (odd lem:o-actiontrunc); along position-dependent generators (Sigma^*, K_b, Delta, so(V)) no bound is given, although ch08 thm:c8-L (iv) and def:c8-twisted use all arities.",
     why_it_matters="Completeness of the arity-two truncation of the MC equation used in the test and in R0.",
     next_action="Weight/degree/x-degree census of local cochains of degree |O|+1-k with the position of weight -1 (ch08 correction 'def:c8-quant' of the design file), for k = 3, 4 on the test families.",
     evidence=["step2/theory/odd_components_design.tex lem:o-actiontrunc and the paragraph after it"])

item(id="IPC2-19", area="R0 (R3, R8)", severity="medium",
     what="(a) R3: the transport of the untwisted Q^xi_1 to the holomorphic letters (ESW Thm 10.9) and the normalisation map to [BK] (3.51), (3.64)-(3.70) are not constructed (engine_plan.json open item 1). (b) R8: the two-letter kernels on all local operators and the spectrum of D_1 against hep-th/0307015v5 (IP7) are not started; IP7 also needs the two-letter kernels of Q_1, S_1 for the kernel form of cor:c8-gram.",
     why_it_matters="R0 includes the action on all local operators and the twisted consequence; the inner product on all operators (IP7) depends on the same kernels.",
     next_action="Schedule after IPC2-09: pilot job of IP7 (pair (A|_Wperp, A|_Wperp), derivative order 2) on the cluster with the file-relative paths fixed; construct the ESW Thm 10.9 transport on letters in code (dictionary.json).",
     evidence=["step2/theory/engine_plan.json open_items[0]", "step2/theory/engine_plan_ip.json IP7"])

item(id="IPC2-20", area="R1, R9, notation", severity="low",
     what="(a) R9 physics-notation callouts are missing in two_point_design.tex, conjugation.tex and odd_components_design.tex (present in inner_product_design.tex). (b) Symbols not defined in Chapter 8 are displayed in the design: O_zeta, O_xi, O_chi, B_z, Fbar_z, F_n, Bbar_n (from odd eq:o-BF), Gram_1/2, X_1/2, Y_1/2, <,>_{rad,tree}, <,>_{rad,1}, <,>_{n,tree}, <,>_{n,1}; two_point_design.tex writes the finite part of operators as varphi (ch08: psi; ch08 uses phi for the degree -1 generator of a field redefinition). (c) Chapter 8 gives no two-component form of Theta_I (CH08-5).",
     why_it_matters="Hard rule 'no new notation beyond Chapter 8' and directive R1/R9 for every report.",
     next_action="Decide (user) whether a definitions block of the inner-product chapter may name the test operators and the order-one-half parts; otherwise write them out; rename varphi -> psi in two_point_design.tex; add the physics-notation remarks.",
     evidence=["grep 'Physics notation' step2/theory/*.tex", "step2/theory/two_point_design.tex eq:t-J"])

item(id="IPC2-21", area="R5 (b): sign of the coupling direction in tp", severity="low",
     what="two_point_design.tex prop:t-laws (F0) still says that the coupling direction multiplies every edge by 1 - a and every vertex of I by 1 + a; with J = I - I^CT - F the signs are opposite (ip thm:ip-laws (c), engine_plan_ip open item).",
     why_it_matters="The coupling direction shifts Gram_1 by a L Gram_tree on length L; a sign error propagates to the length-changing entries.",
     next_action="Correct tp prop:t-laws (F0) and its proof.",
     evidence=["step2/theory/two_point_design.tex prop:t-laws", "step2/theory/inner_product_design.tex thm:ip-laws (c)"])

item(id="IPC2-22", area="R5 (b), (d): single-source constants", severity="low",
     what="The 20' and Konishi one-loop constants 27/2 + 4 log(3/2) and 27/2 + 4 log(3/2) + 24(gamma_E - 3 log 2) rest on one derivation (two_point_scalar_example.py, numerical confirmation only to O(eps log eps)); the lemmas are re-derived (C10), the graph sums are not yet reproduced by a second engine; the 'Gaiotto-Lee-Wu' twisted inner-product target was not found and HSW was read in abstract and introduction only.",
     why_it_matters="Every radial line uses these constants (with the beta-dependence of IPC2-01).",
     next_action="IP1 by both engines with different methods, per class; record the literature items as not used.",
     evidence=["step2/theory/two_point_example_results.json", "step2/literature/ip_targets.json not_found, sources (HSW)"])

# ---------------------------------------------------------------- CH08 issues (format of the task)
ch08 = []
def ch(id_, label, problem, evidence, old, new):
    assert CH08.count(old) == 1, (label, CH08.count(old))
    ch08.append(dict(id=id_, label=label, problem=problem, evidence=evidence, old_text=old, replacement_text=new,
                     line="CH08: %s: %s; evidence: %s; proposed fix: replace '%s' by '%s'" % (label, problem, evidence, old, new)))

ch("CH08-1", "def:c8-2pt",
   "the parenthetical describes the arity-zero finite counterterm as the 27/16 term 'plus s-closed choices'; the cubic part of the obstruction (alpha_1, alpha_3, alpha_5 != 0) forces cubic or quadratic terms that are not s-closed, and the quadratic coefficients that enter the two-point function at relative order one unit are fixed by the arity-one equations",
   "ch09 thm:c9-ob3, thm:c9-F3; odd prop:o-arity0; tp prop:t-exist (only two-leg insertions of order hbar enter at relative order one unit for equal lengths)",
   r"(at one loop $\tfrac{27}{16}" + UNIT + r"\int\mathrm{dvol}\,(d^*A|_V,d^*A|_V)$ plus $s$-closed choices)",
   r"(at one loop the arity-zero part of a solution $\hbar F_1$ of Theorem \ref{thm:c8-obstruction} (ii): in the quadratic sector $\tfrac{27}{16}" + UNIT + r"\int\mathrm{dvol}\,(d^*A|_V,d^*A|_V)$ plus $Q$-closed quadratic terms, tied to the cubic terms by the cubic part of the obstruction and fixed, modulo the arity-zero parts of $s$-closed cochains, by the arity-one components of $sF_1=\mathrm{Ob}$; through relative order one unit the two-point function of two operators of equal length depends on its quadratic part only)")

ch("CH08-2", "def:c8-2pt",
   "the two-point cochain is defined at T = infinity, but its components along K_b and along (eta1, eta2) with Gamma(eta1,eta2)|_V != 0 contain log T (the flows W^(1)(K_b), O[T]^(1)(K_b), W^(2), O[T]^(2) have no limit)",
   "k_judgement.json flows_have_logT, T_to_infinity_table; ch08 prop:c8-Tsafe",
   r"with the global ghosts the same formula with $I_{\mathfrak{psl}(4|4)}$, $F_{\mathfrak{psl}(4|4)}$, $I^{CT}_{\mathfrak{psl}(4|4)}[\varepsilon]$ defines the two-point cochain.",
   r"with the global ghosts the same formula with $I_{\mathfrak{psl}(4|4)}$, $F_{\mathfrak{psl}(4|4)}$, $I^{CT}_{\mathfrak{psl}(4|4)}[\varepsilon]$ and $P(\varepsilon,T)$ in place of $P(\varepsilon,\infty)$ defines the two-point cochain at scale $T$; its components of arity at most one along $\Sigma$, $\Sigma^*$, $\Delta$, the translations and the rotations, and of arity two along $(\xi,\xi')$ and $(\xi,\eta)$, have limits as $T\to\infty$ (Proposition \ref{prop:c8-Tsafe}), while those along $K_b$ and along $(\eta_1,\eta_2)$ with $\Gamma(\eta_1,\eta_2)|_V\neq0$ contain $\log T$ and are used at finite $T$ only.")

ch("CH08-3", "def:c8-inner",
   "the limit at the origin is asserted without its condition, and the identification of the radial Gram with the conjugated two-point coefficient holds only in a scheme with P_1 = 0; with P_1 = p{S^(1)(P_b),-} on a family of weight w the radial Gram is 1 - 2pw times that coefficient",
   "inner_product_design.tex thm:ip-origin (i), prop:ip-frames (a); inner_product_checks.json C3, C15; ip_completeness2_checks.json K3",
   r"On operators at the origin, the radial inner product is the limit as the supports shrink to the origin;",
   r"On operators at the origin, the radial inner product is the limit as the supports shrink to the origin, which exists on a family of one classical Weyl weight exactly when the one-loop dilatation matrix of the family is self-adjoint for the tree part; in a scheme with $P_1=0$ on the family it equals the two-point coefficient of Proposition \ref{prop:c8-2pt} (iii) with the first family conjugated, at the separation $\ell$, and in a scheme with $P_1=p\,\{S^{(1)}(P_b),-\}$ on a family of Weyl weight $w$, $p$ real, it is $1-2pw$ times that coefficient;")

ch("CH08-4", "cor:c8-gram (third relation)",
   "the relation for K_1 is given for P_1 = 0 only; the proposed extensions (design file, k_judgement) do not say that Gram is the radial Gram of that scheme, and with the two-point coefficients the rescaling term has the opposite sign",
   "inner_product_design.tex thm:ip-P (iv), prop:ip-frames; inner_product_checks.json C15; ip_completeness2_checks.json K3",
   r"With $P_1=0$ the third relation is the statement $K=P^\dagger$: $K_1$ along $K_{\bar b}$ equals $\hbar\big[\{S^{(1)}(K_{\bar b}),-\},\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}\big]$.",
   r"With $P_1=0$ the third relation is the statement $K=P^\dagger$: $K_1$ along $K_{\bar b}$ equals $\hbar\big[\{S^{(1)}(K_{\bar b}),-\},\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}\big]$. For $P_1=p\,\{S^{(1)}(P_b),-\}$ on the family, $p$ real, it reads $K_1$ along $K_{\bar b}=p\,\{S^{(1)}(K_{\bar b}),-\}+\hbar\big[\{S^{(1)}(K_{\bar b}),-\},\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}\big]$, with the radial Gram matrix of Definition \ref{def:c8-inner} in that scheme, which is $1-2pw$ times the two-point coefficient on Weyl weight $w$; with the ratio of the two-point coefficients in place of $\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}$ the first term is $-p\,\{S^{(1)}(K_{\bar b}),-\}$, since $\{S^{(1)}(K_b),-\}$ lowers the Weyl weight by one.")

ch("CH08-5", "prop:c8-C",
   "two omissions: Theta_I is not a symmetry of the regularised theory (K_t, P(eps,T) and the eps^{-1} counterterm int(A,A) are not inversion covariant), which is why def:c8-radial goes through Theta_n; and there is no two-component form of Theta_I (needed by engine B and by every two-component report)",
   "conjugation.tex cor:scheme (last sentence); ip_completeness2_checks.json K1 (two-component Theta_I on lambda, lambdabar, lambda*, lambdabar* equal to the index-free -(x,x)^{-2} rho(x) conj(lambda(I x)), (x,x)^{-3} rho(x) conj(lambda*(I x)); convention gate: the same code reproduces the printed Theta_n of def:c8-Theta; Theta_I^2 = 1 on lambda, lambda*, A_m)",
   r"$\Theta_I$ satisfies (i) of Theorem \ref{thm:c8-fields} with $\dagger$ in place of $\dagger_n$, the scalar kinetic density changing by a total derivative.",
   r"$\Theta_I$ satisfies (i) of Theorem \ref{thm:c8-fields} with $\dagger$ in place of $\dagger_n$, the scalar kinetic density changing by a total derivative; it does not satisfy (ii) and (iii) there, since $K_t$, $P(\varepsilon,T)$ and the counterterm $\tfrac94\varepsilon^{-1}" + UNIT + r"\int(A,A)$ are not inversion covariant, and the one-loop radial conjugation is therefore built from $\Theta_n$ (Definition \ref{def:c8-radial}). Two-component: $(\Theta_IA)_m(x)=(x,x)^{-1}\big(\delta_{mk}-2x_mx_k/(x,x)\big)\overline{A_k(I(x))}$, $(\Theta_IA)_{ij}(x)=(x,x)^{-1}\overline{A^{ij}(I(x))}$, $(\Theta_I\lambda)_\alpha{}^i(x)=-i(x,x)^{-2}x_m\sigma^m_{\alpha\dot\alpha}\overline{\bar\lambda_{\dot\alpha i}(I(x))}$, $(\Theta_I\bar\lambda)_{\dot\alpha i}(x)=-i(x,x)^{-2}x_m\overline{\lambda_\alpha{}^i(I(x))}\sigma^m_{\alpha\dot\alpha}$, $(\Theta_I\lambda^*)^\alpha{}_i(x)=-i(x,x)^{-3}x_m\overline{\bar\lambda^{*\dot\alpha i}(I(x))}\bar\sigma^{m\dot\alpha\alpha}$, $(\Theta_I\bar\lambda^*)^{\dot\alpha i}(x)=-i(x,x)^{-3}x_m\bar\sigma^{m\dot\alpha\alpha}\overline{\lambda^{*\alpha}{}_i(I(x))}$, $(\Theta_Ic)(x)=\overline{c(I(x))}$, and $\Theta_IA^*$, $\Theta_Ic^*$ as $\Theta_IA$, $\Theta_Ic$ with the weights $3$, $4$ in place of $1$, $0$ and the factor $-1$.")

ch("CH08-6", "def:c8-inner",
   "the Hermiticity of the reflection and radial forms at one loop is not stated, although thm:c8-adjoint and cor:c8-gram are read as adjointness for an inner product and the generation argument (ip thm:ip-logic (iv)) uses it",
   "thm:c8-fields (ii), (iii); ip lem:ip-ward (order of factors of the two-point convention); the Ward identity along P_n - K_n of thm:c8-adjoint",
   r"and $Y^{\rm adj}:=\mathrm{Gram}_{\rm tree}^{-1}\,\overline{Y}^{\,T}\mathrm{Gram}_{\rm tree}$ is the adjoint for the tree inner product.",
   r"and $Y^{\rm adj}:=\mathrm{Gram}_{\rm tree}^{-1}\,\overline{Y}^{\,T}\mathrm{Gram}_{\rm tree}$ is the adjoint for the tree inner product. In a reflection-compatible scheme $\langle b,a\rangle_n=\overline{\langle a,b\rangle_n}$ modulo $O(\hbar^2)$, since $W(P(\varepsilon,T),J)^{\Theta_n}=W(P(\varepsilon,T),J^{\Theta_n})$, the value of $F^{\Theta_n}$ at the zero field is the conjugate of that of $F$, $(FG)^{\Theta_n}=G^{\Theta_n}F^{\Theta_n}$ and $\Theta_n^2=1$ (Theorem \ref{thm:c8-fields}); with the Ward identity along $P_n-K_n$, which moves $\exp(-\tfrac\pi4\delta^\hbar)$ to the second argument, the same holds for $\langle-,-\rangle_{\rm rad}$, so $\mathrm{Gram}_{\rm tree}$ and $\mathrm{Gram}_{\rm 1}$ are Hermitian.")

ch("CH08-7", "thm:c8-L (iv)",
   "|O| in Lambda^(1)(x)O = (-1)^{|O|(1+|x|)}{S^(1)(x),O} must be the degree for x -> -Lambda^(1)(x) to be a representation on odd O (with the total parity it gives +Lambda^(1)([x,y]) on odd O); the design audit's total-parity rule concerns products of functional-valued cochains; the two QA proposals (ip: degree, design: total parity for the same exponent) are incompatible and must be reconciled by one derivation before either is applied",
   "ch08_corrections_from_ip.json (thm:c8-L (iv)); ch08_corrections_from_design.json (thm:c8-L); odd_rev_audit_algebra a1b, a1c; ip_completeness2_checks.json K4",
   r"$\Lambda^{(1)}(x)O=(-1)^{|O|(1+|x|)}\{S^{(1)}(x),O\}$, $\Lambda^{(2)}(x_1,x_2)O=(-1)^{|O|(|x_1|+|x_2|)}\{S^{(2)}(x_1,x_2),O\}$, $\Lambda^{(k\ge3)}=0$.",
   r"$\Lambda^{(1)}(x)O=(-1)^{|O|(1+|x|)}\{S^{(1)}(x),O\}$, $\Lambda^{(2)}(x_1,x_2)O=(-1)^{|O|(|x_1|+|x_2|)}\{S^{(2)}(x_1,x_2),O\}$, $\Lambda^{(k\ge3)}=0$, where in these two formulas $|O|$ is the degree of $O$ (its ghost number), so that $\Lambda^{(1)}(x)=\{S^{(1)}(x),-\}$ on functionals of degree $0$ of either parity and $x\mapsto-\Lambda^{(1)}(x)$ is a representation on the cohomology also on odd operators; products of functional-valued cochains carry the Koszul signs of the total parity.")

ch("CH08-8", "def:c8-radial",
   "P_1 = 0 is attributed to prop:c7-cov and the reality condition of the finite part that removes a translation rescaling is recorded with opposite signs in the QA files (k_judgement: conj(a) = -a; design, ch09 prop:c9-scheme, ip: real); the definition of delta^hbar does not need P_1 = 0",
   "compare_delta.json comparisons_part2_bilinear_component[9]; ch09 thm:c9-P1, prop:c9-scheme; ch08_corrections_from_design.json def:c8-radial and outside_chapter_8[3]; ip_completeness2_checks.json K2 ({S^(1)(Delta),-} = 2(x.d + w) commutes with Theta_n; compatibility iff the coefficient is real)",
   r"Assume the local representatives along $P_n$ and $K_n$ exist, with $P_1=0$ (Proposition \ref{prop:c7-cov}), and put",
   r"Assume the local representatives along $P_n$ and $K_n$ exist (Proposition \ref{prop:c7-cov} gives the covariance of the regulator, not $P_1=0$; on a family on which $P_1=p\,\{S^{(1)}(P_b),-\}$, the finite part $\tfrac p2\{S^{(1)}(\Delta),-\}$ of the operators of the family gives $P_1=0$ there and is compatible with the reflection of Theorem \ref{thm:c8-fields} exactly when $p$ is real, since $\{S^{(1)}(\Delta),-\}$ acts on a smeared density of weight $w$ as $2(x\cdot\partial+w)$, which commutes with $\Theta_n$), and put")

doc = {
    "file": "step2/qa/ip_completeness2.json",
    "date": "2026-10-29",
    "role": "completeness critic (second pass) of the inner-product design against DESIGN_STEP2 R0 (complete one-loop BV realisation of psl(4|4)) and R5 (one-loop inner product and the test Q^dagger = S); no design or manuscript file edited",
    "read": ["DESIGN.md (rule 4 amended)", "DEFINITIONS_v2.md", "STORY.md", "DESIGN_STEP2.md R0-R9", "HANDOFF_STEP2_2026-10-28.md",
             "manuscript/ch08_definitions.tex (in full)", "manuscript/ch09_cubic_and_dilatation.tex (thm:c9-ob3, thm:c9-F3, thm:c9-kappa, thm:c9-P1, prop:c9-scheme, rem:c9-next)",
             "manuscript/ch01_classical.tex def:c1-bracket, rem:c1-real",
             "step2/theory/inner_product_design.tex (in full)", "step2/theory/engine_plan_ip.json (in full)", "step2/theory/two_point_design.tex (sections 1-4, census)",
             "step2/theory/conjugation.tex (prop:RP, prop:ops, thm:adjoint, thm:inv, cor:scheme, thm:oneloop)", "step2/theory/delta_theory.tex (def:d-PiPsi..prop:d-T, predictions; read-only)",
             "step2/theory/odd_components_design.tex (prop:o-arity0, prop:o-Kclaim, lem:o-actiontrunc; read-only)", "step2/theory/engine_plan.json (E3, E6, E8-E12, T_infinity_safe, open_items)",
             "step2/theory/inner_product_checks.json", "step2/literature/ip_targets.json (sources, targets, not_found, open_items)", "step2/literature/n4_untwisted.json (not_read, B03c record)",
             "step2/qa/ip_findings_disposition.json", "step2/qa/ch08_corrections_from_ip.json", "step2/qa/ch08_corrections_from_design.json", "step2/qa/ch08_corrections_from_delta.json",
             "step2/qa/compare_delta.json", "step2/qa/k_judgement.json, k_judgement.tex (cor:kj-gram)", "step2/qa/odd_rev_audit_algebra/a1*.json"],
    "unit": "hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) ('unit' is a code name only)",
    "names": "xi in Sigma, eta in Sigma^* spinor parameters (ch08); Q_1, S_1, D_1, K_1, P_1, M_1, R_1 one-loop components (def:c8-components); Q the linear BV differential; code names f_P_perp, f_P_lam, f_P_Vc, beta_k, h_j, kappa_k, c_F only in this file",
    "checks": {"script": "step2/qa/ip_completeness2_checks.py", "output": "step2/qa/ip_completeness2_checks.json",
               "K1": "two-component Theta_I (lambda, lambdabar, lambda*, lambdabar*) equals the index-free form; the same code reproduces ch08's printed Theta_n; Theta_I^2 = 1 on lambda, lambda*, A_m",
               "K2": "the classical dilatation 2(x.d + w) commutes with Theta_n; c{S^(1)(Delta),-} is reflection compatible iff c is real; {S^(1)(Delta),O}^Theta = {S^(1)(Delta),O^Theta} for both parities",
               "K3": "-pK + [K, M_tp] = pK + [K, M_tp - 2pW] in a lowest-weight sl(2) module, not = pK + [K, M_tp]",
               "K4": "Lambda^(1)(x) = ad: representation property holds; (-1)^F ad: fails on odd vectors",
               "resources": "0.9 s, 83 MB RSS, one process"},
    "items": items,
    "ch08_issues": ch08,
    "not_found_missing": {
        "dimensional_regularisation": "no step of the design needs DR; the Fourier lemma of tp lem:t-ft continues in the power of |q| at fixed dimension 4 (allowed); DR literature is used as statements only",
        "IR_regulators": "none other than the transfer endpoint T; the only large-distance hazard is the partial composite sets at T = infinity (IPC2-16)",
        "two_component_versions": "present in inner_product_design.tex for every statement; missing for Theta_I in Chapter 8 (CH08-5)"
    },
}
out = os.path.join(HERE, "ip_completeness2.json")
with open(out, "w") as fh:
    json.dump(doc, fh, indent=1, ensure_ascii=False)
print(len(items), "items,", len(ch08), "CH08 issues ->", out)
