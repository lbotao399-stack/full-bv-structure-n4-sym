"""QA: assemble qa/compare.json from the QA work files."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
W = os.path.join(HERE, '..')
QA = os.path.join(W, '..')
ld = lambda f: json.load(open(os.path.join(W, f)))
ce = ld('compare_entries.json')
off = ld('offshell_B_at_A.json')
mom = ld('qa_mom_checks.json')
tw = ld('qa_mom_twist.json')
ha = ld('qa_harmonic_two_site.json')
su = ld('qa_su11.json')
tg = ld('qa_targets_B.json')
twsum = {S: {kind: {k: v['agree'] for k, v in d.items()} for kind, d in r.items() if kind != 'seconds'} for S, r in tw['by_S'].items()}
out = {
 'meta': {
  'role': 'PART C adversarial referee and comparator (QA)',
  'date': '2026-10-30',
  'inputs': ['composite/engine_a/results.json, gates.json, kernels_C2.json', 'composite/engine_b/results.json, gates.json',
             'composite/theory/composite_theory.tex, engine_spec.json, src.nosync/census_sp.json'],
  'qa_code': 'composite/qa/work.nosync/scripts/ (qa_compare.py, b_at_a_points.py, qa_offshell_analyse.py, qa_mom.py, qa_mom_chains.py, '
             'qa_mom_checks.py, qa_mom_twist.py, qa_harmonic.py, qa_su11.py, qa_targets.py, qa_assemble.py); raw outputs composite/qa/work.nosync/*.json',
  'engine_code_use': 'engine B source copied read-only to qa/work.nosync/b_src and run at engine A table points (python -B, no writes to engine folders); '
                     'no engine code used in the QA third method or in the harmonic action',
  'units': 'u = hbar tr_g(ad ad)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); raw chain values in hbar/(16 pi^2) times C_w = (1 (x) ad_{Y_w1}...ad_{Y_wr})(-,-)^{-1}; '
           'planar: (-1)^r c/2 u per adjacent pair (L >= 3), sum_w c_w u on a two-site trace',
 },
 'summary': {
  'engine_disagreements_found': 0,
  'overlapping_entries_compared': len(ce['entries']),
  'overlapping_entries_agree': sum(1 for e in ce['entries'] if e['agree'] is True),
  'single_engine_entries': sum(1 for e in ce['entries'] if e['agree'] is None),
  'offshell_bosonic_tables': 'engine B evaluated at every bosonic row of engine A kernels_C2.json (AV|AV: AV, AV AV, AW AW, AsV c; AV|AW: AW, AV AW, AsW c; AW|AW: AW AW; '
                             '24 rows, C_{P_beta} and C^0_Delta): exact agreement, log2/log3 coefficients included, iff A\'s vertex-vertex sign kappa = -1; '
                             'kappa = +1 fails on every row with a nonzero vertex-vertex part',
  'kappa_resolution': 'kappa = -1 (engine A convention VV_b); in the leg sets AsV c, AsW c A has letter part = vertex-vertex part (pure log2, log3) and B gives 0, consistent only with kappa = -1',
  'third_method': 'QA momentum-space Taylor expansion with exact tensor reduction and radial integrals, own Feynman rules from Ch1/Ch3 and the kernel form of Lemma lem:cc-kernels (i); '
                  'one overall factor F0 = -hbar/2 (sign of the odd skeleton) fixed once; predictions: (AV,AW)->AW at both A points (C_Pb, C^0) and its eps^-1 pole (= B symbolic pole), GC3 (-1/3, -8/3 per colour word), '
                  'AW|AW off-shell table (2 points), twist two S\' <= 4 at generic complex null momenta: P7 (c) and P7 (b) along n and along generic b: all exact',
  'harmonic_action': 'independent implementation of [B03c] (4.27)-(4.29): every spectral target of Prop prop:cc-targets and Table tab:cc-su11 reproduced; '
                     'engine D_1 matrices at theta = -1/3 equal -8 x (H12 + H21) at matrix level in the engines\' bases (twist two S\'<=4, fermion S\'<=5, gluon S\'<=2, (F11,psi) S\'<=4)',
 },
 'entries_static': ce['entries'],
 'offshell_B_at_A_points': off,
 'third_method_scalar_sector': {'k1_and_bilinear_and_AWAW_table': mom, 'twist_two_family': {'kinematics': tw['kinematics'], 'agree_by_S': twsum,
                                'detail': tw['by_S']}},
 'k1_closed_forms_and_pole': ld('qa_k1_closedforms.json'),
 'harmonic_action_two_site': ha,
 'su11_harmonic_action': su,
 'B_family_targets_recomputed': tg,
 'not_compared': [
  'fermionic off-shell tables of engine A (AV|lamP, AW|lamP, lamP|lamM, ...): the Step-1 spinor bases of the two engines differ; no map constructed',
  'poles z_m (eps^-1 coefficients) between engines beyond the k=1 (AV,AW)->AW pole (QA third method = B)',
  'gluon S\'=2, fermion S\'=4,5, (F11,psi) S\'<=4, psi psi\'/Z F11 S\'<=2: engine B only',
  'generic-beta P7(b) and density universality incl. lambda_+ lambda_-: engine A only',
  'GC1, GC2, GC5c, GC6, GC7, GC8, GC12 (gl(3)), GC13, TC0, TC2: run by neither engine',
 ],
}
json.dump(out, open(os.path.join(QA, 'compare.json'), 'w'), indent=1)
print('compare.json written', os.path.getsize(os.path.join(QA, 'compare.json')))
