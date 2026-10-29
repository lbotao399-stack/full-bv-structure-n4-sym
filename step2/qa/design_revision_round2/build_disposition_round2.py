"""Append the round-2 dispositions to step2/qa/design_revision_disposition.json (key round_2)."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
P = os.path.join(HERE, '..', 'design_revision_disposition.json')
d = json.load(open(P))
D = 'step2/theory/odd_components_design.tex'
EP = 'step2/theory/engine_plan.json'
CH = 'step2/qa/ch08_corrections_from_design.json'
R = 'step2/qa/design_revision_round2/'
items = []
def it(i, lens, sev, loc, verdict, action, evidence):
    items.append({'id': i, 'lens': lens, 'severity': sev, 'location': loc, 'verdict': verdict, 'action': action, 'evidence': evidence})

it('R2-A1', 'algebra-scheme', 'major', 'thm:o-KD, thm:o-invariants (iv), prop:o-S (ii), prop:o-Kconsistent (iii), freedom tables; E8, E10c, E10e, E11',
   'accepted (refined)',
   D + ': new prop:o-P1 (translation component, P_1 = 0 as a renormalisation condition, values -(1/3+f_P) on bilinears, -f_P on letters, R_1 = 0 in every so(W)-covariant scheme, M_1 on the composite line); prop:o-laws (a) exact direction kappa S^(1)(Delta), (b) family finite parts (c/2) ad_{S1(Delta)} pi_F; thm:o-KD with general P_1 = p ad (K_1 == 1/4 div v gamma O - p ad_{S1(K_b)} O) and values in the P_1 = 0 and translation-covariant schemes; prop:o-S (ii) with -p ad_{S1(eta)}; thm:o-invariants (iv) invariant combination; both tables gain the rows exact direction, f_P, z, composite remainders, family finite parts. ' + EP + ': E11 moved before E8/E10 with the three-letter test, E8/E10c targets in both schemes, E10e with pbar. Refinement: the projection pi_n on n-letter functionals of the fix is replaced by a family projection pi_F commuting with ad_{S1(Delta)}, ad_{S1(P_b)}; with c_n != c_{n+1} the polynomial-degree projection does not commute with Q+{I,-} on gauge-invariant operators with covariant derivatives and would give them L_1(O)^(0) != 0; c_n = -(n-1)/3 is kept as the E11 test.',
   'engine_b/delta_results.json part2 (coefficient.normalisation_of_the_classical_part, translation_part_check subset 844 = full); qa/compare_delta.json [9] (three routes); ' + R + 'scheme_shift_round2.json (kappa ad_{S1(Delta)}: P -2 kappa, K +2 kappa, so(V) 0, Delta 0; exact); odd_rev_audit_algebra a3 part 4, a5 (rerun identical)')
it('R2-A2', 'algebra-scheme', 'major', 'prop:o-arity0; E3b; targets F0',
   'accepted with correction',
   D + ': prop:o-arity0 restated (three-leg part in span{H1..H6,H7,H8+,H8-}; closed space dim 9, 8 even; ker Q dim 5; span of the six dim 6 (5 even, 1 odd); quotient dim 3, all even, representatives (b1,H1)=(1,1), (b2,H4)=(1/4,1), (b3,H5)=(-1/2,1); Yukawa and H5+H7 expressed in them); two-component H7, H8+-; E3b rewritten; targets F0 corrected. Correction to the finding: the parity-even quotient is 3-dimensional, not 2 (the audit subtracted the full rank 6 of the six from the parity-even dimension 8, but one of the six, s(N_lp) - s(N_lm) = -2(H8+ - H8-), is parity odd); the Yukawa direction and H5+H7 are not two extra directions beyond the Part-C family: modulo the span of the six every closed direction is a Part-C direction (C = V + Part-C).',
   R + 'cubic_quotient.py/.json (exact rationals, inputs transcribed from cubic_structure.tex prop:s2c-QH, cor:s2c-QI3, prop:s2c-br as in audit a4); odd_rev_audit_algebra a4 rerun identical')
it('R2-A3', 'algebra-scheme', 'major', 'CH08 def:c8-radial', 'accepted (merged with R2-G6 and the K-judgement proposal; reflection condition corrected)',
   CH + ' entry def:c8-radial: drop "with P_1 = 0 (Prop c7-cov)", append the scheme-condition paragraph with the family finite part (p/2) ad_{S1(Delta)} pi_F, its effect (K_1 +p, Q_1 -p/2, S_1 +p/2), reflection compatibility for real p, and the exact direction kappa S^(1)(Delta) (-2kappa, 2kappa, -kappa, kappa)',
   'engine B, compare_delta, scheme_shift_round2.json; reflection: ad_{S1(Delta)} = 2(x.d + w) on densities, real and r_n-invariant')
it('R2-A4', 'algebra-scheme', 'minor', 'CH08 thm:c8-adjoint', 'accepted in modified form',
   CH + ' entry thm:c8-adjoint: the hypothesis P_1 = M_1 = R_1 = 0 is dropped (the proof uses only reflection compatibility and Cor c8-homotopy for the full Lambda^(1)); the R_1 statement of the finding is in design prop:o-P1 (iv); prop:c8-primary gets P_1 = M_1 = R_1 = 0 instead',
   'proof of thm:c8-adjoint; k_transfer cor:ktr-gram; k_judgement cor:kj-gram')
it('R2-A5', 'algebra-scheme', 'minor', 'CH08 cor:c8-gram', 'accepted (merged with the K-judgement formula)',
   CH + ' entry cor:c8-gram: K_1 along K_bbar = pbar ad_{S1(K_bbar)} + hbar[ad_{S1(K_bbar)}, Gram_tree^{-1} Gram_1] for P_1 = p ad; design prop:o-Kconsistent (iii)',
   'k_transfer cor:ktr-gram; ch08 thm:c8-dag (i)')
it('R2-A6', 'algebra-scheme', 'minor', 'prop:o-GB (ii) proof; CH08 prop:c8-B proof', 'accepted',
   D + ': prop:o-GB (ii) proof rewritten via G_B^(1)(P_b) = 0 (induction on the degree in the position, translation-invariant phi^(0), then the four-term contradiction along S_+ (x) W); ' + CH + ' entry prop:c8-B',
   'odd_rev_audit_algebra a3 part 2 (S+W alone: c = 1, all a_f equal; both blocks: none; translation: c = 0), a7 (x-dependent degree -2 weight 9/2 contents); rerun identical')
it('R2-A7', 'algebra-scheme', 'minor', 'CH08 thm:c8-L (iv), lem:c8-ev', 'accepted',
   D + ': lem:o-ev states the total-parity rule and the product intertwining; thm:o-twist (H1) removed; ' + CH + ' entries thm:c8-L, lem:c8-ev; ch02 def:c2-dce noted (outside chapter 8)',
   'odd_rev_audit_algebra a1 (T1-T4, T6), a1b (V1 total parity 10/10), a1c (V1 8/8); rerun identical')
it('R2-A8', 'algebra-scheme', 'minor', 'CH08 def:c8-quant', 'accepted (extended)',
   CH + ' entry def:c8-quant: Weyl weight 4 with the weight convention for the position (-1) and the arguments (half the ad_Delta eigenvalue); the same convention in def:c8-renorm and Weyl weight 4 for phi in prop:c8-laws (a); design def:o-names (weights), prop:o-laws (a)',
   'design thm:o-invariants (ii) proof (S^(1)(P_b) shifts D_1 by 2 hbar ad_{S1(P_b)})')
it('R2-A9', 'algebra-scheme', 'minor', 'CH08 remark after thm:c8-twist; design thm:o-twist (iv)', 'accepted',
   D + ': thm:o-twist (iv) with C S^xi + Im Q^xi_0 (N S^xi = Q^xi_0 S^(1)(R_xi), r_xi(G_B) = -Q^xi_0 S^(1)(R_xi)), comparison remark updated; ' + CH,
   'odd_rev_audit_algebra a2 ({psi,f} = {f,psi}: 20/20)')
it('R2-A10', 'algebra-scheme', 'minor', 'CH08 remark (existence along K); design def:o-K, prop:o-Kclaim, E10', 'accepted (merged)',
   D + ': rem:o-Kcomposite (composite line K.(P*)^2 + [P*,K.]P*, complete composite sums, R_{O,K}); E10c requires complete composite sums; ' + CH + ' entry for the remark (merged with the K-judgement text)',
   'delta_theory rem:d-Y1a, prop:d-T (iii); k_judgement lem:kj-composite')
it('R2-A11', 'algebra-scheme', 'minor', 'def:o-B; targets F5', 'accepted',
   D + ' def:o-B: bracket with diag(0_4,1_4), equivalently minus the bracket with diag(1_4,0_4) (as ch08 def:c8-B); step2/literature/q1_targets.tex F5 and build_q1_targets.py F5 (json regenerated)',
   'odd_rev_audit_algebra a3 part 1 (block table: B = ad diag(0_4,1_4) = -ad diag(1_4,0_4); ad Delta eigenvalues reproduced)')
it('R2-A12', 'algebra-scheme', 'minor', 'prop:o-GB (iii), eq:o-twistlaws, thm:o-invariants (iii), tables, rem:o-alpha; prop:o-laws (a)', 'accepted',
   D + ': hbar inserted in every law (L_1 changes multiplied by hbar in the components); prop:o-laws (a) with (-1)^{1+|x_1|}',
   'ch08 prop:c8-laws ("-a hbar"); shuffle sign of ch02 def:c2-dce for the (1,0) shuffle; odd_rev_audit_algebra a3 part 4')
it('R2-A13', 'algebra-scheme', 'minor', 'caption of the D_1..R_1 table; arity-one counterterm row', 'accepted',
   D + ': caption: int(c*, iota_v A_V) exists for every even generator with a vector field on V (translations, rotations, Delta, K); none for so(W_R^perp); row entries "s-exact on gauge-invariant O"',
   'odd_rev_audit_algebra a7 (content (A,c*) degree -2 weight 5 along P_b)')
it('R2-A14', 'algebra-scheme', 'minor', 'thm:o-invariants (ii); D_1 column', 'accepted',
   D + ': weight 4 imposed on the finite counterterms (def:o-names, prop:o-laws (a)); thm:o-invariants (ii) proved for every weight-4 phi^(0), position-dependent or not ({S^(1)(Delta),-} = 2(x.d + 4) on a weight-4 density, a total derivative), with the counterexample S^(1)(P_b); D_1 column gains the z row (fixed at arity two, delta_theory prop:d-z) and the composite remainder R_O',
   'delta_theory thm:d-D3, prop:d-z; engine A z_2, engine B b (part2.coefficient)')
it('R2-A15', 'algebra-scheme', 'minor', 'thm:o-twist (H1); engine_plan open_items, E0', 'accepted',
   D + ': (H1) replaced by the statement of the sign rule with citation of ch08 lem:c8-ev and thm:c8-L (iv) (CH08 correction); ' + EP + ': E0 task and conventions.evaluation_along_xi updated, open item removed',
   'odd_rev_audit_algebra a1, a1b')
it('R2-A16', 'algebra-scheme', 'minor', 'thm:o-twist (ii) proof', 'accepted', D + ': (Q^xi_0)^2 = -1/2 {{S^xi,S^xi},-}',
   'odd_rev_audit_algebra a2 (ad_f ad_f = -1/2 ad_{f,f}: 8/8)')
it('R2-A17', 'algebra-scheme', 'minor', 'remark (comparison with the twisted computation)', 'accepted',
   D + ': "the twisted BV action [BK, (3.53)] has no coupling constant; loops are counted by hbar [BK, (1.3)]"',
   'odd_rev_audit_algebra/lit_check/bk_pdf/bk_v3.txt: no occurrence of "coupling"; (1.3) is the hbar expansion of Q, (3.53) the Chern-Simons-type action on C^{2|3} (checked in this round)')
it('R2-A18', 'algebra-scheme', 'minor', 'lem:o-Tflow (iv)', 'accepted', D + ': T in (0,infinity), and T = infinity in the components of prop:o-Tsafe (i)', 'ch08 prop:c8-Tindep (i)')
it('R2-A19', 'algebra-scheme', 'minor', 'lem:o-actiontrunc, prop:o-trunc', 'accepted',
   D + ': hypothesis that the Sigma-components carry no explicit position dependence, with its justification (I^CT by prop:c7-cov; the schemes used); x-dependent contents at k = 3, 4 recorded in the proof',
   'odd_rev_audit_algebra a7 (x-dependent contents (lam*,c*)x at degree -3, (c*,c*)x^2 at degree -4)')
it('R2-A20', 'algebra-scheme', 'minor', 'def:o-scale, targets F0 notation', 'accepted',
   D + ': F^(k)_psl for arities, F_psl = hbar F_1 + ... (ch08); targets F0/F1 now use F^(0)_psl, two-leg / three-leg part in words', 'ch08 def:c8-quant, def:c8-2pt')
it('R2-A21', 'algebra-scheme', 'minor', 'rem:o-alpha; targets F1', 'accepted',
   D + ' rem:o-alpha and targets F1: the a of F1 (added to the action) is minus the a of prop:o-laws (e)', 'J = I_psl - F_psl - I^CT: F -> F + hbar a G_B changes the twisted action by -a hbar (S^(1)(xi) - S^(2)(xi,xi))')

it('R2-G1', 'graphs-T', 'major', 'thm:o-KD, prop:o-S (ii), freedom table; E8, E10c, E10e, E11', 'accepted (merged with R2-A1; one formula corrected)',
   'as R2-A1; the relations K_1 == 1/4 div v gamma O - c ad_{S1(K_b)}, S_1(eta) == Q_1(rho(x)eta) - c ad_{S1(eta)} + U_eta and the E10c targets adopted, with two-component forms. Not adopted: the finite part "psi\' = c x (Weyl weight)": multiplication by c times the weight of the density does not commute with integration by parts in omega and is not a map of local functionals; the design uses (c/2) ad_{S1(Delta)} pi_F, which has the same component shifts',
   R + 'scheme_shift_round2.json (check (4): c x density weight does not commute with derivatives; ad_{S1(Delta)} does); odd_rev_audit_graphs scheme_shift_check.json (rerun identical)')
it('R2-G2', 'graphs-T', 'major', 'lem:o-legs, prop:o-ir, thm:o-ir, prop:o-Tsafe; engine_plan placement, large_distance, T_infinity_safe', 'accepted',
   D + ': lem:o-legs (even-even adjacency through A_V (x) A_V, two ghost ends; odd pairs of opposite chirality; two-component chain); prop:o-ir (degree d1+d2-2); thm:o-ir (v) (orientation sum, which pairs converge, (K_b,x) at finite T); prop:o-Tsafe (i), (iii); ' + EP,
   'odd_rev_audit_graphs census_ir.json, koszul_orientation.json (relative sign -1, leading term 0, first moment = prediction); rerun identical')
it('R2-G3', 'graphs-T', 'major', 'thm:o-classes, eq:o-L1graphs; engine_plan formula, E6, E12', 'accepted',
   D + ': class (T0) added; eq:o-L1graphs now subtracts the trees of hbar L_1(O)^(0) with one S^(1)(x1) and adds the arity-one tree flow of the counterterm of the cochain sO; reduces to the former formula for gauge-invariant O with L_1^(0) = 0; ' + EP + ' graph_class_reference.T0/formula, E6',
   'Lemma c8-rep, thm:c8-L (iv): L(O)[T] at arity one = (sO)[T]^(1) + hbar L_1(O)^(1)[T]_tree + trees of hbar L_1(O)^(0) with S^(1)(x1); letter int omega (X, iota_v A_V): L_1^(0) a multiple of int omega (X, iota_v dc), tree with one lambda leg (same leg count as L_1^(1)(xi)); delta_theory thm:d-D3 subtracts the same tree')
it('R2-G4', 'graphs-T', 'major', 'prop:o-laws (a), tables, caption, prop:o-arity0, prop:o-GB (ii), thm:o-invariants, actiontrunc; E2, E3b, E4, E5', 'accepted except the sign of the exact direction',
   D + ': exact direction kappa S^(1)(Delta) in prop:o-laws (a), both tables, thm:o-invariants (iii), (iv), with the shifts (P, K, xi, eta) = (-2kappa, +2kappa, -kappa, +kappa) hbar x classical; caption corrected; prop:o-GB (ii) repaired through the translation component (the two-block repair of the finding needs an unproved lemma on position-dependent invariants and is not adopted); lem:o-actiontrunc and prop:o-trunc restricted to x-independent Sigma-components; E3b stated for translation-covariant cochains. Rejected: the shift pattern (2a, -2a, a, -a) of the finding has the wrong overall sign (it is the sign of -(s phi)^(1)(x) with the shuffle sign dropped)',
   R + 'scheme_shift_round2.json (exact: -2 kappa along P, +2 kappa along K); odd_rev_audit_algebra a3 part 4 (same); k_transfer prop:ktr-rescale (i) (psi\' = (a/2) ad_{S1(Delta)}: P_1 -a, K_1 +a, Q_1 -a/2, S_1 +a/2) and delta_theory prop:d-P (ii) agree')
it('R2-G5', 'graphs-T', 'major', 'engine_plan E10c, E10d, E0, E9', 'accepted',
   EP + ': E10c exact procedure (eps -> 0 at fixed T = eps^0 coefficient of the integrand with the P-edge upper limit extended, per-graph proof of the needed Taylor order, exact or >= 40 digits + PSLQ, cross-check at T|p|^2 = 1, 4, 16); E0 canonical form with gamma^{mn} d_m d_n delta (test on {S^(1)(K_b), int(A,A)}); E9 both orientations in one integrand; finite-T propagator in every K job',
   'k_judgement thm:kj-Tindep (bound O(eps) at fixed T); census_ir.py (S^3 average 2 eta.xi per orientation); engine conventions alpha delta + beta d delta only (engine_a/engine_b delta_results conventions)')
it('R2-G6', 'graphs-T', 'major', 'CH08 def:c8-radial', 'accepted (merged; see R2-A3)', CH + ' entry def:c8-radial (family finite part (p/2) ad_{S1(Delta)} pi_F instead of c x Weyl weight)', 'as R2-A3, R2-G1')
it('R2-G7', 'graphs-T', 'minor', 'CH08 prop:c8-B', 'accepted in substance (repair by the translation route)', CH + ' entry prop:c8-B', 'as R2-A6')
it('R2-G8', 'graphs-T', 'minor', 'CH08 prop:c8-Tsafe', 'accepted (merged with the K-judgement coefficient and the composite-sum sentence)', CH + ' entry prop:c8-Tsafe', 'census_ir, koszul_orientation; k_judgement prop:kj-composite')
it('R2-G9', 'graphs-T', 'minor', 'lem:o-Tflow (iv), def:o-scale, prop:o-unique', 'accepted', D + ': T ranges (0,infinity), T = infinity in the convergent components (prop:o-Tsafe (i))', 'ch08 prop:c8-Tindep; thm:o-ir (ii)-(v)')
it('R2-G10', 'graphs-T', 'minor', 'prop:o-uv, eq:o-defect, thm:o-classes (E1)', 'accepted', D + ': defect vertex and classes stated at psi = 0, the extra vertex delta hbar((d_CE+Q)psi(O) + {I_psl,psi(O)}_eps) displayed, its effect is the law (b)', 'ch08 def:c8-renorm')
it('R2-G11', 'graphs-T', 'minor', 'def:o-B', 'accepted (= R2-A11)', D + ' def:o-B', 'odd_rev_audit_algebra a3 part 1')

it('KJ-1', 'K judgement', 'decision', 'prop:o-Kclaim (placeholder), def:o-K, thm:o-ir (ii)-(iv), prop:o-Tsafe, E10', 'integrated',
   D + ': the placeholder replaced by Proposition "Existence and independence of the endpoint along K" (route (A); proof from k_transfer thm:ktr-defect and k_judgement thm:kj-Tindep, numbers from the validated engine kn), def:o-K with the order of limits, thm:o-ir (ii)-(iv) with the established coefficients and trees, rem:o-Kcomposite, prop:o-Tsafe (ii)-(iii); ' + EP + ' E10a-E10g',
   'step2/qa/k_judgement.json/.tex; step2/theory/k_transfer.tex')
it('KJ-2', 'K judgement', 'decision', 'thm:o-KD, prop:o-Kconsistent', 'integrated',
   D + ': thm:o-KD general P_1 = p (d = 0, m = 0, alpha = -p, kappa = gamma - p w), predicted f_S = -1/3 noted for the dilatation workflow; prop:o-Kconsistent (i) (K_b,P_b\') and (P_b\',Delta), (ii) (eta1,eta2), (iii) Cor c8-gram with pbar',
   'k_transfer thm:ktr-relations, thm:ktr-KD, cor:ktr-gram; k_tr_relations.json rerun identical in k_judgement')
it('KJ-3', 'K judgement', 'dispute', 'k_transfer prop:ktr-rescale (ii); k_judgement cor:kj-gram, open issue 2, CH08 proposal for def:c8-radial',
   'rejected (reflection condition)',
   'the design and ' + CH + ' state reflection compatibility of the finite part (a/2) ad_{S1(Delta)} for real a; consequently a real rescaling p is removable in a reflection-compatible scheme',
   'ad_{S1(Delta)} acts on densities of weight w as 2(x.d + w), real coefficients, commuting with the linear isometry r_n and preserving weights; Theta_n antilinear => (kappa ad_{S1(Delta)} O)^Theta = conj(kappa) ad_{S1(Delta)} O^Theta; in the bracket computation of ktr, S^(1)(Delta)^Theta = -S^(1)(Delta) is compensated by the exchange sign {g,f} = (-1)^{|f||g|}{f,g} of the ch01 bracket (odd_rev_audit_algebra a2: 60/60), which ktr omitted; a6: ad_Delta commutes with dagger_n; a5 agrees')
it('KJ-4', 'K judgement', 'evidence', 'citations of k_transfer prop:ktr-values (ii), prop:ktr-Tsafe rows 5-6',
   'noted', D + ' cites the K-judgement values (validated engine kn, third realisation) for every number; the partial files k_tr_etaeta.py, k_tr_kernel.py, k_tr_zeromode.py, ktr_out/, qa/k_transfer_numeric/ were reused only through the K judgement, which reran them (theory scripts bit for bit, numerics on subsets)',
   'step2/qa/k_transfer_numeric/PARTIAL_FILES_NOTE.txt; k_judgement.json argument item 1')
it('KJ-5', 'K judgement', 'decision', 'CH08 proposals of the K judgement', 'adopted (merged)',
   CH + ' entries def:c8-K, thm:c8-K, prop:c8-Tsafe, remark after thm:c8-K, def:c8-radial (reflection condition corrected, KJ-3), thm:c8-adjoint, prop:c8-primary, cor:c8-gram', 'k_judgement.json open_issues')

d['round_2'] = {
    'date': '2026-10-29',
    'input': 'auditor findings (lenses algebra-scheme: 21, graphs-T: 11) on the revised design, the K judgement (step2/qa/k_judgement.json), alignment with manuscript/ch08_definitions.tex',
    'files_changed': [D, D.replace('odd_components_design.tex', '_build_odd/odd_design.pdf'), EP, 'step2/literature/q1_targets.tex',
                      'step2/literature/build_q1_targets.py', 'step2/literature/q1_targets.json', 'step2/literature/build_q1targets/q1_targets.pdf',
                      CH, 'step2/qa/design_revision_disposition.json (this key)'],
    'new_checks': [R + 'cubic_quotient.py/.json', R + 'scheme_shift_round2.py/.json', R + 'rerun/ (audit scripts a1-a7, census_ir, koszul_orientation, scheme_shift_check rerun: outputs byte-identical)',
                   R + 'build_ch08_corrections.py (asserts every quoted ch08 text occurs exactly once)', R + 'update_engine_plan.py'],
    'backups': R + 'backup/ (pre-round-2 copies of the design, plan, targets, disposition)',
    'counts': {'findings': 32, 'accepted': 26, 'accepted_with_correction_or_partly': 6, 'rejected': 0, 'K_judgement_items': 5, 'ch08_corrections': 14},
    'items': items,
    'open_after_round_2': [
        'E11: composite remainders along P_b and so(V) on a three-letter operator; M_1 rescaling; scheme record (translation covariant / P_1 = 0) before E8, E10',
        'E10c, E10f: K_1 on letters and bilinears from the K-route graphs at finite T; completion of the arity-one obstruction along K_b',
        'E3b: the three directions of C/V admit no s-closed psl(4|4)-completion (hypothesis (H4))',
        'f_S = -1/3 predicted by d = 0 (dilatation workflow)',
        'transport Phi from ESW Thm 10.9 (rem:o-transport); gauge-fixing law (c): proof plan']}
json.dump(d, open(P, 'w'), indent=1, ensure_ascii=False)
print(len(items), 'items')
