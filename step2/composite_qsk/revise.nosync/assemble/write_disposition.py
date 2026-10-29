"""Writes step2/composite_qsk/disposition.json (reviser's decision on every referee finding)."""
import json, os
W = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/composite_qsk'
RV = os.path.join(W, 'revise.nosync')
wc = json.load(open(os.path.join(RV, 'word_map', 'word_compare.json')))
bc = {'agree': wc['bosonic_summary']['agree'], 'rows': wc['bosonic_summary']['rows'], 'agree_p1': wc['bosonic_summary']['agree_with_vv_sign_+1']}
ge = json.load(open(os.path.join(W, 'composite', 'engine_a', 'gates_ext.json'))) if os.path.exists(os.path.join(W, 'composite', 'engine_a', 'gates_ext.json')) else {}
ext_done = {k: v['pass'] for k, v in ge.items()}
nf, naf, nfp = wc['fermionic_summary']['rows'], wc['fermionic_summary']['agree'], wc['fermionic_summary']['agree_with_vv_sign_+1']

D = {
 'meta': {'role': 'reviser', 'date': '2026-10-30', 'referee_file': 'composite/qa/findings.json',
          'qsk': 'theory Q, both QSK engines and the QSK referee returned nothing: no QSK finding exists to dispose of (qsk_final.json)'},
 'findings': [
  {'id': 0, 'severity': 'major', 'claim': 'P7(a) and completeness of the quadratic Q-closed directions untested; D_1 = 4 lc^natural conditional',
   'decision': 'accepted (partly fixed, partly reformulated)',
   'reason': 'Checked: GC13 and TC0 were not run by either engine; Def cc-free had "possibly further Q-closed quadratic directions".',
   'actions': [
     'TC0 computed exactly (revise.nosync/tc0/tc0.py, tc0.json): parity-even Spin(4)xSpin(6)-invariant quadratic degree -1 weight-4 functionals, '
     'position degree <= 1: dims (perp, lambda, V+c) = (2, 3, 6); Q-closed (1, 2, 2); Q-exact (0, 1, 1); cohomology = span{S_perp, S_lambda, S_V}; '
     'kappa_V, kappa_lambda directions exact. Without parity one more class int(lambda*, rho(v1v2v3v4) lambda) (free-gaugino chiral rotation), '
     'excluded by the reflection used in Thm c9-quadratic (iii). New Prop prop:cc-tc0; Def cc-free amended.',
     'New Prop prop:cc-theta-sector: under the hypotheses of Thm cc-arity2, (b) on a single-sector family forces theta_diamond = -1/3; with the '
     'computed (b) on twist two (perp), fermion sl(2) (lambda), gluon sl(2) (V) this gives (a); (a) is a consequence and the cubic arity-one '
     'equation must admit delta = 0 (consistency condition for the existence of the quantisation). New Cor cor:cc-noA: D_1 = 4 lc^natural without assumption (a).',
     'GC13 (cubic rank test with the cubic composite part) not computed: needs the composite loops of I with three legs, a new graph class; '
     'recorded as open in the status table.'],
   'status_after': 'TC0 closed; (a) reduced to (b) + existence of the arity-two quantisation; GC13 open'},
  {'id': 1, 'severity': 'major', 'claim': 'Thm cc-local marked proved; proof is a sketch; direct consequences not computed',
   'decision': 'accepted',
   'reason': 'Checked: steps (C), (D) of the proof are asserted; GC1, GC2, GC5c, GC6, GC7, GC8 are "not run" in both gates.json. The referee is right that a '
             'uniform one-body term t u v. would be absorbed into theta and stay invisible in the spectra.',
   'actions': ['Proof labelled "Outline"; new Rem rem:cc-localstatus lists what is checked (no log eps numerically, GC3 via Prop cc-routes, scalar sector, '
               'spectra) and what is not (GC5c, GC6, GC7, GC8, Pi_Delta on A_V, lambda, c); states that theta = -1/3 rests on (iv).',
               'Status table: Thm cc-local "outline; consequences checked indirectly"; Cor cc-D1local, Prop cc-routes "conditional on Thm cc-local (iv)".'],
   'status_after': 'downgraded; GC5c/GC6/GC7/GC8 remain open (not run in this revision)'},
  {'id': 2, 'severity': 'minor', 'claim': 'Def cc-sp read literally gives Z_x the opposite sign to Prop cc-master',
   'decision': 'accepted, fixed',
   'reason': 'Checked against manuscript Def c1-bracket (Delta_0 = -d_{K_0}) and Prop c4-defect1: {f,g}_eps - {f,g} = -(f, g joined by K_eps - K_0); '
             'so Q(composite graphs) = +(hbar/2) d_{K_x} of trees with one edge K_eps - K_0, while Def cc-sp said "minus".',
   'actions': ['Def cc-sp now defines Z_x := (hbar/2) d_{K_x} d_{K_eps - K_0}(chains) with the Koszul sign s(c,sigma) = (-1)^{n(<sigma)+n(<c)+[c<sigma]}; '
               'Prop cc-master states the bracket-defect = -(edge K_eps - K_0); proof of Thm cc-local adjusted.'],
   'status_after': 'closed'},
  {'id': 3, 'severity': 'minor', 'claim': "engine A's vertex-vertex sign is -1; tables correct only read with -1",
   'decision': 'accepted, fixed and derived',
   'reason': 'Derivation (engine_a/notes.txt sec. 2b): A\'s chain T = L_1 S h_1 ... L_n S h_n S L_0 S puts the parity matrix on both sides of the last line; '
             'the special edge joins slots of opposite parity, so at X_b it carries the extra -1 = (-1)^[c<sigma] of s(c,sigma); on interior lines '
             'the factor must be explicit: vertex-vertex sign = -(VV_a) = -1.',
   'actions': ['engine A defaults set to -1 (cq_eval, cq_gval.VVSIGN, cq_twist, cq_fsl2, cq_su11); kernels_C2.json rows carry value = letter - vertex_vertex; '
               'run_all.py rerun: 26/26 gates (617 s, 241 MB).',
               'Affected comparison rerun word by word (revise.nosync/word_map/word_compare.py; A word = B word with reversed indices): %d/%d bosonic off-shell rows agree with B (%d/%d with sign +1); A*_V c and A*_W c leg sets vanish.' % (bc['agree'], bc['rows'], bc['agree_p1'], bc['rows']),
               'Fermionic rows (new, see finding 4): %d/%d agree with sign -1; %d/%d would agree with +1.' % (naf, nf, nfp, nf)],
   'status_after': 'closed'},
  {'id': 4, 'severity': 'minor', 'claim': "engine A's fermionic sign conventions chosen to give the P7 answer; not an independent test",
   'decision': 'accepted; resolved by derivation and by a direct fermionic comparison',
   'reason': 'Checked engine_a/notes.txt sec. 2: R2 and the intrinsic sign were selected by requiring su(1|1) universality.',
   'actions': ['Derivation: R2 = (-1)^{n(<sigma)+n(<c)} is s(c,sigma) without the mutual-order factor, which A\'s chain convention supplies (finding 3); the '
               'intrinsic pair sign is the two-letter engine\'s (ea_kern.kernel). The selection therefore coincides with the first-principles rule.',
               'Exact intertwiner between A\'s and B\'s realisations of Sigma (revise.nosync/ferm_map/spinor_map.py): unique up to sign, S^T (Mp R_B[M]) S = sb_A[M] '
               'for all ten M with identical vector labels; it maps A\'s lamP (+1 eigenspace of rho(e0e1e2e3)) onto B\'s lamM: the engines orient V oppositely.',
               'Fermionic off-shell comparison (ferm_compare.py, chirality-pure letters, letters as covectors, outputs as vectors, all leg sets with <= 2 outputs '
               'of the nine pairs with a gaugino letter or output): %d/%d rows agree exactly (rational and log parts).' % (naf, nf)],
   'status_after': 'closed'},
  {'id': 5, 'severity': 'minor', 'claim': 'several family results rest on one engine only',
   'decision': 'accepted, partly fixed',
   'reason': 'Checked against both gates.json: gluon S\'=2, fermion S\'=4,5, (F11,psi), psi3psi4/ZF11 were B only.',
   'actions': ['engine_a/ext_families.py: fermion sl(2) S\'=4, 5, gluon S\'=2, twist two S\'=5, 6 in route A. Outcome: %s.' % (json.dumps(ext_done) if ext_done else 'see results_ext.json'),
               '(F11,psi) and psi3psi4/ZF11 remain route B only (route A has no code for these two families); generic b on non-collinear legs beyond twist two not run.'],
   'status_after': 'reduced; two mixed families single-route'},
  {'id': 6, 'severity': 'minor', 'claim': 'P7(b) along n not independent of P7(c)',
   'decision': 'accepted, fixed',
   'reason': 'Verified: Lemma cc-canonical with Prop cc-position gives [C^(2),0_Delta, d_b] = 2 C^(2)_{P_b} pair by pair, and [N, d_b] = d_b; so (c) on F and d_bF implies (b) along b.',
   'actions': ['New Rem rem:cc-P7dep (1) after Conjecture cc-P7.'], 'status_after': 'closed'},
  {'id': 7, 'severity': 'minor', 'claim': '"(b) and (c) contain no scheme parameter" hides regulator dependence',
   'decision': 'accepted, fixed',
   'reason': 'Verified: manuscript Rem c10-scheme1 (d) and Prop c10-shape: finite parts of special-edge cycles depend on the cutoff shape.',
   'actions': ['Sentence replaced by Rem rem:cc-P7dep (3): identities at the regulator of Def c8-PT; scheme-independent statement = spectrum at the arity-two-fixed coefficients.'],
   'status_after': 'closed'},
  {'id': 8, 'severity': 'minor', 'claim': 'spec deliverables incomplete (TC2, B kernels_C2.json, families.json, plain cycles and ct^(1), gl(3))',
   'decision': 'accepted in part',
   'reason': 'Checked: as stated.',
   'actions': ['gl(3) realisation of the planar adjacent rule: exact gl identity (Casimir = P/kappa) r = 1,2,3, orthonormal and rescaled bases: pass (revise.nosync/colour_gl3).',
               'Route-B values at route-A table points (22 bosonic from the QA runs, 26 fermionic from this revision) written to composite/engine_b/kernels_C2_at_A_points.json; agreement summary and the family matrices of both routes in results/composite_final.json (closed_forms.off_shell_tables, families).',
               'Not done (rejected for this revision, cost): TC2 tables with letter derivatives n <= 3 in general colour position for all 45 pair types; tabulation of the plain composite cycles and ct^(1)_x (implemented in route A, kinds Plog/Peps, untabulated).'],
   'status_after': 'partly open'},
  {'id': 9, 'severity': 'minor', 'claim': "engine B's vendored realization.py, loops.py not verified against the originals",
   'decision': 'accepted; functional regression done, byte diff still impossible',
   'reason': 'The originals engine_b/realization.py, loops.py are still dataless; a 90 s copy attempt timed out (Operation timed out).',
   'actions': ['Fermion sl(2) S\'=3 rerun with the vendored modules (revise.nosync/b_vendor_check, 83 s, 109 MB): all 40 family entries and all 109 density entries '
               'identical to the 20:27 output made with the original modules.',
               'The exact intertwiner of finding 4 reproduces B\'s Clifford data from the vendored realization.py: Clifford, pairing and rho on Sigma* checks pass.'],
   'status_after': 'functionally verified; byte diff pending the iCloud download'},
 ],
}
json.dump(D, open(os.path.join(W, 'disposition.json'), 'w'), indent=1)
print('disposition written')
