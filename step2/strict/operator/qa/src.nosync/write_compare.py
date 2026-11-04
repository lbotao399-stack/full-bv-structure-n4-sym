"""Assemble qa/compare.json from the QA comparator outputs (cmp_tables, cmp_families, recheck_rowsA, toy_special_cycle)
and the hand derivations recorded below."""
import json, os
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/operator'
Q = f'{ROOT}/qa/out.nosync'
ct = json.load(open(f'{Q}/cmp_tables.json'))
cf = json.load(open(f'{Q}/cmp_families.json'))
ra = json.load(open(f'{Q}/recheck_rowsA.json'))
ty = json.load(open(f'{Q}/toy_special_cycle.json'))

def summ_fam(d):
    return {S: {'sum_Upsilon0_Delta_entrywise_equal': v.get('sum_Upsilon0_Delta_equal_entrywise'),
                'D1_eigenvalues_A': v.get('D1_eigs_A'), 'D1_eigenvalues_B': v.get('D1_eigs_B'),
                'agree': bool(v.get('sum_Upsilon0_Delta_equal_entrywise') and v.get('eigs_agree'))} for S, v in d.items()}

compare = {
 'file': 'step2/strict/operator/qa/compare.json',
 'date': '2026-11-04',
 'role': 'operator-side referee/comparator: engine A (route A) vs engine B (route B), entry by entry, with independent QA checks',
 'units': 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); -32 h(J) u = (g^2 N/pi^2) h(J)',
 'inputs_read': ['engine_a/results.json, gates.json, notes.txt, out.nosync/tables_compact.json, out.nosync/partial_*.json (via results)',
                 'engine_b/results.json, gates.json, notes.txt, out.nosync/tables/*.json',
                 'theory/operator_theory.tex (all 846 lines)', 'manuscript ch10 (Thm c10-brst, eq c10-finitesp)',
                 'ch11 (Def c11-defect, Thm c11-classes, eq c11-pole) from the size-identical scratchpad copy (manuscript file is iCloud-dataless)',
                 'step2/composite_qsk/results/qsk_summary.txt; qsk engine_b notes (GitHub copy, size-identical)',
                 'step2/literature/q1_targets.json'],
 'code': ['qa/src.nosync/cmp_tables.py', 'qa/src.nosync/cmp_families.py', 'qa/src.nosync/recheck_rowsA.py',
          'qa/src.nosync/toy_special_cycle.py', 'qa/src.nosync/write_compare.py'],
 'A_vs_B': {
  'F0_scalar_bilinears': {k: ({'A': v.get('A'), 'B': v.get('B'), 'agree': v['agree']} if 'agree' in v else
                              {kk: vv for kk, vv in v.items()}) for k, v in cf['F0'].items()},
  'F0_further': {
     'line_class_C0_Delta': {'A': '-8/3 (same-order word)', 'B': '-8/3 O', 'agree': True},
     'special_class': {'A': 'exch 16, trace -8 per word (sum 16 O - 16 T)', 'B': '16 O - 16 (w1,w2)(Y0,Y1)', 'agree': True},
     'Upsilon0_eta_singlet': {'A': 'zero on the colour singlet (single-lambda leg set only)', 'B': '0 (both chiralities)', 'agree': True},
     'Upsilon0_K_singlet': {'A': '0', 'B': '0', 'agree': True},
     'Q1_traceless_coefficient': {'A': '79/12 - 3 log(3/2) - theta^xi_perp', 'B': 'not assembled (two-body only)',
                                  'QA': '79/12 - 3 log(3/2) - theta^xi_perp: two-body -2L F1 plus one-body -(L - 79/12 + theta^xi_perp) F1 (T2 and S^perp_2 act on A_perp letters as the classical arrow); equals Part Q c_Q = -2L - theta - p, p = L - 79/12', 'agree': True},
     'one_body_per_letter_Delta': {'A': 'A (w=1+k): 2k/3 - 20/3; lambda (w=3/2+k): 2k/3 - 44/3', 'B': '-4 n_X + c_X - 2 theta (w_X + k), c_A = 2/3, c_lam = 1/3, c_c = 4/3',
                                   'QA': 'F_comp = -(1/3) N_A - (1/12) N_lam - (4/3) N_c from s(N_A, N_lam, N_c) = ([-1,0,-2],[0,4,0],[0,0,1]) and the composite part u[1/3,-1/3,-2/3] = Ob1(Delta)|_2 - 4 s(uN); with {N_A,A} = 2A, {N_lam,lam} = 4 lam, {N_c,c} = c this gives c_A = 2/3, c_lam = 1/3, c_c = 4/3: both engines reproduced', 'agree': True}},
  'F1_su2': {'density': {'A': '-8u(1 - P + K/2) on SO(6)', 'B': '-8u(1 - P) on {Z,X}', 'agree': True,
                         'note': 'K annihilates the su(2) sector; same Minahan-Zarembo density'},
             'spectra_L_le_8': {'agree_with_independent_Heisenberg_A': cf['F1_su2_independent']['all_agree'],
                                'agree_with_independent_Heisenberg_B': cf['F1_su2_independent']['all_agree']}},
  'F2_su11': {'two_site_densities': {'A': '(Z,Z) 20/3; (Z,psi) same 8/3, exchanged 8; (psi,psi) -4/3',
                                     'B': 'L2_M1 sector Upsilon0 = 64/3 = 2(8/3 + 8); L2_M0 40/3 = 2(20/3)', 'agree': True},
              'spectra_L_le_6': {'A_and_B_agree_with_independent_free_fermion_diagonalisation': cf['F2_su11_independent']['all_agree']}},
  'F3_twist_two_S_le_6': summ_fam(cf['F3_twist_two']),
  'F4_fermion_sl2_S_le_5': summ_fam(cf['F4_fermion_sl2']),
  'F3_F4_line_class_and_Pn': {'A': 'G-D1sum: sum Upsilon0 = 4 lc^nat + 16 - (2/3)(4+S\'); Upsilon_Pn = -1/3 d_n',
                              'B': 'line class = -(2/3)(4+S\') 1; special class = 4 lc; Upsilon_Pn = -1/3 d_n',
                              'QA': 'B Upsilon_Pn matrices equal -1/3 times the d_n matrix built here on every step (F3 0..6, F4 1..5)', 'agree': True},
  'F5_F6_F7': {'status': 'engine A only (route B did not compute F5-F7): no G-AB possible',
               'QA_recomputation': cf['R7_matrix_identities'].get('F5', {}).get('eigenvalues_recomputed'),
               'QA_recomputation_F7': cf['R7_matrix_identities'].get('F7', {}).get('eigenvalues_recomputed')},
  'tables_F8': {
     'overlap_rows_by_generator_n_le_1': ct['overlap_rows'],
     'superficial_degree_D_mismatches_on_P_rows': len(ct['D_mismatch']),
     'G_D1pair_verdict_disagreements': ct['conj_disagree'],
     'G_D1pair_disagreement_reading': 'all 11 are route-B rows with a missing class (EP timed out or ESW above the memory cap); engine A has the full row and the pair form holds there',
     'F_x_nonzero_sets_identical': all(x['A_F_nonzero'] == x['B_F_nonzero'] for x in ct['F_nonzero_A_vs_B']),
     'F_x_nonzero_rows': [x['key'] for x in ct['F_nonzero_A_vs_B']],
     'B_rows_absent_in_A_D0': '60 rows with D_EP = 0 and no ESW block: the Delta intrinsic kernel is identically zero (degree -1); A omits them',
     'numeric_values': 'not comparable row by row: A evaluates at seeded Gaussian-integer letters in the Chapter 2 realisation, B symbolically or at its own points in the 54-slot basis; the families F0-F4 are the value-level comparison'},
  'gates': {'A': '86/86 pass', 'B': '70/82 pass; 12 fail = 10 with the theory F_x added (finding QA-F1) + 2 reduced G-pole gates whose failing rows are line-placement rows (QA-F8) and ESW-incomplete rows (QA-F6)'}},
 'QA_recheck_of_engine_A_rows_from_raw_numbers': {
     'pair_form': 'eps0 = 4 lc + line (D0), = line (P, M0), = 0 (R) recomputed exactly on all 5280 rows: ' + str(ra['counts']),
     'Ward_pole_form': 'total poles = -4 m ct_m (D0) / 0 (P, M0, R) holds on every pole entry without line placements; on the 204 pole entries with line placements the stored total pole differs, i.e. the line class carries poles (consistent with B-F3; full eq c11-pole untested)',
     'shifted_eq_c11_pole_counterexamples': ra['shifted_pole_identity_counterexamples']},
 'independent_QA_checks': {
   'spectra_vs_psu224_harmonic_action': cf['spectra_vs_independent_harmonic_targets'],
   'selection_rule_used': 'two-site psu(2,2|4) modules V_J with E_12 = 2h(J), L = 2 cyclic => D_1 = -32 h(J) u, J even; F3: J = S, F4: J = S+1, F6: J = S+2, F7: J in {S+1, S+2}; spectrum at S\' = union over S <= S\' (descendants)',
   'R7_matrix_identities': cf['R7_matrix_identities'],
   'prop_op_tfix_by_hand': 't_Delta = 0; t(S_+W*) = -alpha, t(S_-W) = -beta, t_P = alpha + beta, t_K = -(alpha + beta); t = aB + (c/2)lambda with a = (alpha-beta)/2, c = -(alpha+beta) reproduces all blocks; Rem op-tc-main factors 1-(a-c/2), 1-(-a-c/2), 1-(-a+c/2), 1-(a+c/2), 1+c, 1-c verified',
   'weight_count_Thm_op_master_iv': 'total weight 4 - lambda_x/2 - 2m at m = 1: 2,3,2,2,5/2,3/2,1 < 7/2 - p_x = 5/2,7/2,5/2,7/2,7/2,5/2,3/2 along Delta,P,M,R,Sigma,Sigma*,K: verified',
   'composite_line_Delta_Q_closed': '[Q, P* D. P*] = -4 Q^GF (E (g + H g\') - g H E\') = 0 since g + H g\' = E\' and gH = E: verified by hand at every eps, T',
   'toy_special_edge_cycle': ty}}

inv = {
 'D_1': {'representative': 'sum_{a<b} +- [eps^0] D_Delta(X_a,X_b) - {hbar F_1^(1)(Delta), O}; hbar F_1^(1)(Delta) = 4uN + F_comp + u(theta S_perp + theta_lam S_lam + theta_V S_V)(Delta) + kappa_V, kappa_lam (s-exact) [+ G_B, coupling direction]; theta = theta_lam = theta_V = -1/3 fixed by the arity-two equation along (Delta, P_b)',
         'scheme_independent': 'spectra on Q_0-cohomology: 0 (half-BPS bilinear, tr Z^L, j = 0 twist two), -48 u (Konishi), -32 h(J) u on F3-F7, -8u E_Heisenberg on F1, -8u E_su(1|1) on F2; D_1 commutes with the classical action (R7)',
         'artefacts': 'individual pair kernels and their log 2, log 3 parts; the split special (4 lc) / line (C_Delta) / one-body; 40/3 and -104/3 alone; matrix entries in a given basis'},
 'P_1': {'representative': '-(1/3 + theta_P) u d_b O on every family (two-body -1/3 from the line class, special class 0)',
         'scheme_independent': 'none alone: theta_P is the c-direction of t = aB + (c/2)lambda; P_1 = 0 is reachable (Def c8-radial). Invariant: the (R2) relation c_Q + c_Qbar = c_P on cohomology',
         'artefacts': '-1/3 itself'},
 'M_1, R_1': {'representative': 'M_1: two-body 0 on A_perp letters, one-body -(1/3) u spin_Om (A_V; lambda coefficient flagged 1/3 vs -1/12) - theta_M M.; R_1 = -theta_R R. (Upsilon_R = 0 on every row)',
              'scheme_independent': 'zero: the restriction of delta to so(V) + so(W_R^perp) is a coboundary (R8, Whitehead)', 'artefacts': 'all numbers'},
 'Q_1': {'representative': 'F0, psi = 0: two-body -2L F1 + 4L F2 + 4L F3 + 0 F4 (L = log 3/2), one-body -{(L - 79/12)(T1+T2) + theta^xi_V S^V_2 + theta^xi_perp S^perp_2, O} [+ a B]; traceless: c_Q = 79/12 - 3L - theta^xi_perp; cocycle scheme removes 4L F3',
         'scheme_independent': 'c_Q + c_S (same chirality block) tied to (1/2)(xi,eta) D_1 by (R3); c_Q + c_Qbar - c_P by (R2); the class [Q^xi_1] mod Im Q^xi_0 and C S^xi for the twist (not yet computed)',
         'artefacts': 'c_Q alone (theta^xi_perp free, outer a B), every log(3/2) term, 79/12, the F2/F3 entries individually'},
 'S_1': {'representative': 'two-body: Q_1 along rho(x_O) eta (position part, G-pos graph by graph) + Upsilon^0_eta (zero on the scalar singlet); one-body: NOT AVAILABLE (no particular solution of the action-side QME along eta)',
         'scheme_independent': 'only through (R3) with Q_1 and D_1, and (R6); untested', 'artefacts': 'everything individually'},
 'K_1': {'representative': 'two-body: v(x_O) P_1 - (b,x_O) D^0_1 + M^0_1 + K^0_1 (G-pos); Upsilon^0_K = 0 on the singlet; one-body NOT AVAILABLE',
         'scheme_independent': 'through (R5) with D_1, M_1, P_1 and through K = P^dagger (Thm op-c), both untested', 'artefacts': 'individual entries'},
 'L_1^(0)': {'representative': '-(2 + 4L) u (w1,w2) int omega (dc, A_V) at psi = 0',
             'scheme_independent': 'zero: removed by psi(O) = -(1 + 2L) u (w1,w2) int omega (A_V, A_V)', 'artefacts': 'the whole coefficient'}}
compare['invariant_content'] = inv
compare['not_comparable_or_not_done'] = [
  'F5, F6, F7: route A only', 'F8 n >= 2: route A only (B n <= 1); n = 3 not run in A; antifield pairs n >= 2 in neither',
  'xi, eta tables: A n <= 2 (xi), n <= 1 (eta); B n = 0 only', 'K_b tables: A n <= 1 intrinsic, B OP4/OP8 samples',
  'G-twist (F9), G-KP, G-KPbdry, G-Cid, G-typeI, G-logcov, full G-pole, two-sided (R2)-(R6), M_1 and generic-b P_1 on F3-F7, Q_1/S_1/K_1 on F1-F7: neither route']
os.makedirs(f'{ROOT}/qa', exist_ok=True)
json.dump(compare, open(f'{ROOT}/qa/compare.json', 'w'), indent=1, default=str)
print('written', len(json.dumps(compare)))
