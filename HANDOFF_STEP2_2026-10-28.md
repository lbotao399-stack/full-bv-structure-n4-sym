# Step 2 handoff — 2026-10-28 (paused for the 5h usage limit)

Binding spec: DESIGN_STEP2.md (Parts C, D; user directives R0–R7). DESIGN.md rule 4 amended (transfer endpoint T allowed;
core quantities T-independent; T -> infinity where safe). Memory: n4-costello-uv-boundary-program.md.

## Workflows (all stopped by me at the pause; resume = same script + resumeFromRunId, completed agents replay from cache)
Scripts in ~/.claude/projects/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/workflows/scripts/

| run | script | status at pause |
|---|---|---|
| wf_20406333-59a | n4-step2-cubic-and-dilatation-wf_20406333-59a.js | 11 agents done (cubic engines A/B, compare_cubic, Δ theory/engines in progress) |
| wf_c6617485-557 | n4-step2-q1-targets-notation-design-wf_c6617485-557.js | COMPLETE (Q1 targets, dictionary 118/118, odd design + audits) |
| wf_f56a6c01-b4b | n4-step2-inner-product-design-wf_f56a6c01-b4b.js | 3 foundation agents done (ip_targets, conjugation, two_point_design); synthesis/audit/revise not run |
| wf_c35ebbf2-79d | n4-step2-odd-design-revision-and-k-transfer-wf_c35ebbf2-79d.js | 1 agent done; revision + K-at-finite-T theory/numerics in progress |

Resume (same session only): Workflow({scriptPath: "<script>", resumeFromRunId: "<run>"}).
In a new session: the journals are in .../subagents/workflows/<run>/journal.jsonl; rerun the script fresh (outputs on disk
are kept; agents overwrite their own files) or hand-author a continuation from the journal.
Note: wf_f56a6c01 PRE still says "no IR scale"; R7 in DESIGN_STEP2.md supersedes it (edit PRE before rerunning).
Killed at pause: step2/qa/k_transfer_numeric/ir_eta.py (partial).

## Outputs on disk (step2/)
- engine_a/, engine_b/: cubic_results.json, Δ (D4) partial runs; qa/compare_cubic.json, masters_indep, plugback_*.
- literature/: sweep, bgkwwy, gkw, n4_untwisted, q1_targets (11 targets), ip_targets (inner product).
- notation/: esw_to_weinberg_srednicki.tex/.pdf, dictionary.json, check_dictionary.py (118/118).
- theory/: cubic_structure.tex, odd_components_design.tex (AUDIT NOTICE at top), engine_plan.json,
  conjugation.tex, two_point_design.tex, two_point_census.json, k_tr_*.py (partial).
- qa/design_audit_corrections.json (36 flagged items: outer derivation B freedom; K-route log terms incl. eta-eta
  P-edge chains; F1 = alpha-conjugation; c_K scheme dependent; Zwiebel 1111.0083 prior planar 2->1 S, K).

## Next steps
1. Finish wf_20406333 (cubic obstruction α1..α5, F_3 family, D1–D4) and review.
2. Finish wf_c35ebbf2 (design revision; K_1 at finite T, T-independence; judge).
3. Finish wf_f56a6c01 (inner product synthesis, Q† = S test, audits).
4. Engines: Q_1, S_1, K_1 with scheme freedom; inner product; Q† = S test.
5. Chapters (FORMAT.md, theorem form) with Weinberg–Srednicki versions; PDF; reader.

## Resumed 2026-10-28 evening
- User: reports go into the reader (append chapters, never touch learner highlights); cluster allowed for pure compute;
  core = clarify the definitions first; 2-loop+ ahead. Ultracode off (no Workflow runs without opt-in).
- Wrote manuscript/ch08_definitions.tex (definitions chapter, all orders in hbar, Weinberg–Srednicki two-component
  forms); independent referee review (step2/qa/ch08_review/findings.json: 0 fatal, 6 major, 27 minor, 1 markup) —
  all applied; v1 kept as step2/qa/ch08_review/ch08_definitions_v1_reviewed.tex.
- Added to the live reader: revision 50 -> 51, 413 new blocks with IDs c8-NNNNN (+ sentence fragments), all 2350
  earlier blocks byte-identical, 17 learner-highlight blocks intact. Tools: reader_build/prep_new_chapter.py,
  reader_build/merge_into_live.py (use --replace-prefix c8 for later revisions of chapter 8; it refuses if learner
  formatting on a changed block would be lost). Backups: reader/qa/backup_port8787_rev000050_*.
- Paused workflows (see table above) were not resumed.

## Continuations launched 2026-10-28 (late evening), aligned with Chapter 8
Finished results of the paused runs were extracted to step2/qa/resume_2026-10-28/wf{1,2,3}_done.json and are read by
the new agents (the old runs were not resumed; their unfinished agents are re-run in new scripts).
| run | script | content |
|---|---|---|
| wf_86256db2-042 | n4-step2-dilatation-theory-continued-wf_6eae9d6e-1fe.js | Part D theory (D1-D3, composite-line remainder Q(Y1a), P_1), 3 audits, engine-vs-theory comparison, cubic refuter #2, draft manuscript/ch09_cubic_and_dilatation.tex |
| wf_ca3e017f-d02 | n4-step2-k-transfer-continued-wf_ca3e017f-d02.js | K_1 at finite T (theory + numerics), 2 refuters, judge, final fix of odd_components_design.tex; CH08 corrections -> step2/qa/ch08_corrections_from_design.json |
| wf_14c0ce14-cc9 | n4-step2-inner-product-continued-wf_14c0ce14-cc9.js | inner-product synthesis on Chapter 8, 3 refuters, revision, critic; CH08 corrections -> step2/qa/ch08_corrections_from_ip.json |
(The first launch wf_6eae9d6e-1fe was stopped at once: placeholder args.)
Already established before the pause: cubic obstruction Ob_3 = unit[(11/4+19/4 log(3/2))(G1+G3) + 27/8 G2 - 9 log(3/2) G5]
(engines A = B, refuter #1 confirms alpha_3, alpha_4); Konishi 0 / -24 confirmed by both Delta engines; open: Feynman-gauge
composite A_V--c lines give a nonlocal Q-exact remainder (Delta) and a translation component -(1/3) unit (bilinears).

## 2026-10-29 status
- wf_86256db2 (dilatation) DONE: cubic sector fully confirmed (refuter #2); Part D theory partly corrected by audits;
  draft manuscript/ch09_cubic_and_dilatation.tex (referee review running, step2/qa/ch09_review/); 34 CH08 items in
  step2/qa/ch08_corrections_from_delta.json (main: D_1 = 2c O - (1/3+theta) units {S^(1)(Delta),O} with theta a scheme
  parameter fixed by the (Delta,P_b) component; P_1 = 0 needs a scheme choice; Prop c8-Tsafe proof wording; Lemma c8-transfer
  conventions). Apply after the other two runs report.
- Physics-notation supplements (user 2026-10-29): manuscript/supplements/key.tex + chNN.tex; tools reader_build/prep_supplements.py,
  insert_supplements.py (inserts after anchors, never edits existing blocks). ch03-05 drafted; ch00-02, ch06-08 running.

## 2026-10-29 (night) reader update
- ch09 refereed (step2/qa/ch09_review/findings.json: 0 fatal/major, 19 minor) — all applied (pre-fix copy ch09_before_fixes.tex).
- ch08 consolidated from the 4 correction streams (step2/qa/ch08_consolidated/: disposition.json 31 applied / 21 merged /
  13 rejected; ch08_before_consolidation.tex). Labels unchanged. ch07 rem:c7-internal overstatement ("no one-loop
  correction") fixed in manuscript + live reader (reader_build/fix_c7_internal.py; backup ch07_before_internal_fix.tex).
- Live reader 8787: rev 58 -> 59 (c8 replaced) -> 60 (c9 appended, IDs c9-NNNNN) -> 61 (ch07 remark). 19 learner-formatting
  blocks byte-identical. Scratch copy: scratchpad/reader_scratch on port 8899.
- reader_build/chapter_ids.py: live IDs = b-NNNNN for chapters 0-7, cN-NNNNN for chapters >= 8 (used by prep_new_chapter.py
  [new usage: CHAP LIVE_IDS OUT], prep_supplements.py, insert_supplements.py). Re-run insert_supplements.py after any
  --replace-prefix merge (it re-places all ph* blocks at their anchors).
- Supplements (key, ch00-ch09; 112 blocks, 675 reader blocks, IDs ph*) reviewed (step2/qa/physics_supp/review/) and inserted: live rev 62. Re-run insert_supplements.py 8787 after any supplement edit or chapter replacement.

## 2026-10-29 R8 workflow launched (user go-ahead)
- run wf_c061cef0-9fc, script ~/.claude/projects/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/workflows/scripts/n4-two-letter-counterterm-kernels-wf_c061cef0-9fc.js
- Theory -> engines A (Schwinger/Gaussian) || B (momentum/Feynman + numeric) -> compare -> 2 refuters -> revise -> ch10 draft
  (manuscript/ch10_two_letter_kernels.tex, supplements/ch10.tex). Outputs under step2/twoletter/. Reader not touched by it.
- Scope: counterterms of all local operators (two-letter kernels, output bound, covariance reduction) + D_1/P_1 two-body.
  Q_1, S_1, K_1 two-body kernels from the defect graphs = next workflow.
- 09:15 first run wf_c061cef0-9fc died (network ECONNRESET, all 8 agents; theory left census.json, tl_simplex.py, checks).
  Relaunched as wf_f959070c-6dd (same script + resume note for theory, PROGRESS.md logs, stop guards).
- 2026-10-29 16:xx wf_f959070c-6dd COMPLETE: theory (step2/twoletter/theory/twoletter_theory.tex, 22 pp), engines A/B agree
  389/389 rows; numbers referee: no wrong number (7920+560+3833 entries rechecked); theory referee: 4 major (reduction prop,
  D_1 on all operators not established, P7 unproven, Q1/S1/K1 two-body absent) — all accepted and addressed by the reviser
  (results/kernels_final.json, qa/disposition.json). Draft manuscript/ch10_two_letter_kernels.tex + supplements/ch10.tex.
  Open: composite kernels C2_Delta, C_{P_b} beyond the planar scalar sector (D_1 there conditional on P7); Q1/S1/K1 two-body
  kernels; W2 with antifields; field-strength families. ch10 referee running -> step2/qa/ch10_review/.
- ch10 refereed (step2/qa/ch10_review/findings.json: 2 major, 15 minor, 7 markup — all applied; pre-fix copies in pre/).
  Live reader rev 63 (c10 appended) -> 64 (126 supplements re-placed). GitHub pushed via API (scratchpad api_push.py):
  dab670b (ch10 + step2/twoletter, big engine_b/kernels*.json gitignored), fd67eca (README). Local .git is behind origin.

## 2026-10-29 composite + Q/S/K workflow (user go-ahead)
- run wf_5e5a15e1-2e4, script .../workflows/scripts/n4-composite-and-odd-two-body-kernels-wf_5e5a15e1-2e4.js; root step2/composite_qsk/
- theory C (C2_Delta, C1_Delta, C_{P_b}, P7 test) || theory Q (Q_1, S_1, K_1 two-body, twist check) -> composite engines A||B
  -> QSK engines A||B || composite referee -> QSK referee -> revise -> manuscript/ch11_composite_and_odd_kernels.tex + supplements/ch11.tex
- 2026-10-30 wf_5e5a15e1-2e4 done for PART C: P7 (b),(c) exact on all families (twist two S'<=6, su(1|1) L<=6, fermion sl(2)
  S'<=5, gluon S'<=2, psi psi'/Z F11, (F11,psi)); arity two forces theta_perp = theta_lambda = theta_V = -1/3 (conditional on
  arity-two quantisation + Thm cc-local, outline proof); D_1 = 4 lc^nat, P_1 = 0 in the translation-covariant scheme;
  spectra -8u E_H. Results step2/composite_qsk/results/composite_final.json; draft manuscript/ch11_composite_and_odd_kernels.tex
  (11.1-11.5 final draft, 11.6 structure only). PART Q failed (network) -> rerun wf_a7e26f82-275 (script copy
  step2/composite_qsk/partq_workflow.js, retries on dropped agents), ends by replacing Section 11.6.
- ch00_overview.tex had been evicted by iCloud (dataless); restored from the byte-identical GitHub copy (sha 7df6c9c).
- ch11 referee (11.1-11.5) running -> step2/qa/ch11_review/.
- ch11 refereed (27 findings: 3 major overclaims, all applied; pre-fix copies step2/qa/ch11_review/pre/). Live reader rev 65
  (c11 appended) -> 66 (136 supplements). GitHub a1f594b (ch11 + step2/composite_qsk), d3e5156 (README/CONTRIBUTING).
  reader_build prep scripts now put ~/.npm/_npx/4b4c857f6efdfb61/node_modules first in NODE_PATH (the Documents copy of
  markdown-it is partly evicted); regenerated supplement blocks byte-identical to the live ones.
- When Part Q (wf_a7e26f82-275) rewrites Section 11.6: referee it, then merge_into_live.py 8787 ch11_blocks.json
  --replace-prefix c11 and re-run insert_supplements.py 8787; push via scratchpad api_push.py.
- 2026-10-30 user directive R10 (strict Costello route; DESIGN_STEP2.md). User approved the strict-route workflow, to be
  launched WHEN Part Q (wf_a7e26f82-275) finishes: script step2/strict_workflow.js (copy of scratchpad strict.js);
  root step2/strict/ (action/ and operator/ sides); ends by rewriting ch09 §9.2-9.5, ch10 §10.8, all of ch11, and writing
  manuscript/ch12_strict_construction.tex + supplements/ch12.tex.
- 2026-10-30 Part Q DONE (wf_a7e26f82-275): referee verdict fail (fatal: engine A rounding bug zeroed twist beta entries;
  sign of Lemma qs-arrows wrong; route B missed a T5 placement) -> all fixed by the reviser, 79/79 gates
  (results/qsk_final.json, qsk/disposition.json). Section 11.6 rewritten (24 statements; NOT yet refereed, NOT in reader).
  Open: TQ1 (action side along xi), F_x on D>=2 rows, Upsilon_xixi, R2-R6, class-level twist comparison.
- Strict-route workflow launched: wf_5bfa6f93-8f4 (script scratchpad strict.js = step2/strict_workflow.js).
  Reader/GitHub update for ch09-ch12 after it finishes + referee.

## 2026-10-30 14:xx PAUSED at user request
- Strict-route run wf_5bfa6f93-8f4 STOPPED during the operator engines. Done and cached: theory-action, theory-operator,
  engine-action-A/B, refute-action (pass_with_fixes). Action-side results so far (two engines agree exactly):
  Ob^(1)(Delta)|_2 regression OK; Ob^(1)(xi)|_2 = u[(79/24 - 1/2 log(3/2))(g1+g3) + 27/8 g2]; F_1^(1)(xi)|_2 =
  (log(3/2) - 79/12)u(T1+T2) + theta^xi_V S^V_2 + theta^xi_perp S^perp_2; Ob^(2)(Delta,P_b)|_2 class -2/3 u per sector
  => theta = -1/3 COMPUTED at two legs; (Delta,xi) arity two residual class -1/3 u[S(xi)]. Referee: printed WZ signs in
  the action theory wrong; open: K_b, eta solution, S^(2)-pairs ({Q,Q},{Q,S}), (K_b,P_b'), (xi,P_b), 3-4 leg classes.
- Script amended (user 2026-10-30): invariant content (class mod Q_0-exact/homotopy, scheme artefacts named), one-body
  parts from F_1^(1)(x), class-level twist comparison via Phi = high priority; new agent 'action-complete' for the open
  action items. RESUME: Workflow({scriptPath: "<scratchpad>/strict.js" (copy: step2/strict_workflow.js),
  resumeFromRunId: "wf_5bfa6f93-8f4"}). Operator engines restart (their files in step2/strict/operator/engine_*/ stay).

## 2026-11-03 resumed (user)
- Disk 98% full (11 GB free); iCloud had evicted 10 manuscript files: 8 restored from byte-identical GitHub copies,
  ch11 + supplements/ch11 downloaded by reading them. Strict-route run resumed: wf_5bfa6f93-8f4 (task wov0e7nzm).
- GitHub push script now persistent: tools.nosync/github_api_push.py LIST MSG (scratchpad copy was lost).
  Backup push of qsk results + 11.6 draft + step2/strict in progress.
- User 2026-10-31: cluster allowed for this project when helpful (budget shared with the fortuity session, which owns
  the cluster jobs; coordinate). Candidates: operator engines to high weight, eta intrinsic kernel (memory), 2-loop.
- 2026-11-04 strict-route run wf_5bfa6f93-8f4 COMPLETE (12 agents). theta = -1/3 computed twice from action-level arity two
  ((Delta,P_b) and (Delta,xi)); F_1^(1) along Delta, xi, eta, K_b on two legs; invariant relations
  t_xi+ + t_xi- - t_P = 27/4 + 2log(3/2) = t_xi+ + t_eta+ - t_Delta, t_M = t_Delta, t_R = t_Delta + 1/3; D_1 from defect
  cycles alone on F2-F7 = Beisert spectra; finite-special-vertex double count found and removed; WZ signs corrected.
  Chapters 9-11 rewritten (rewrite_log.json), ch12_strict_construction.tex new. Open: (K_b,P_b'), (eta,eta'), 3-4-leg
  classes, Q/S/M/K_1 beyond F0, G-Cid/G-typeI (locality), Phi. Referee of ch09-12 running -> step2/qa/ch09_12_review/.
- 2026-11-04 ch09-12 refereed (31 findings: 8 major overclaims/leftovers, all applied; pre-fix copies
  step2/qa/ch09_12_review/pre/). Live reader rev 71 (c9/c10/c11 replaced, c12 appended, 156 supplements);
  reader_build/update_ch9_12.sh PORT does the whole sequence. README/CONTRIBUTING updated; GitHub push of ch09-12 +
  step2/strict running (tools.nosync/github_api_push.py).

## 2026-11-05 (user: "continue, don't stop")
- Three background agents (local compute): step2/next/action2 ([K,P], (eta,eta'), (xi-,eta+-), three/four-leg classes),
  step2/next/operator2 (Q/S/M/K_1 on F1-F7, generic-b P_1, route B on F5-F7, locality gates, (R2)-(R6)),
  step2/next/twist (transport map Phi, class-level BK comparison, Upsilon_xixi).
- Supplements ch01/ch03: [X] slot notation replaced by (Qe)_X; functional side sA_m = -d_m c (substitution rule of
  Prop c9-maps proof; common +d_m c is c -> -c). Reader rev 115 (targeted reader_build/update_one_supplement.py keeps
  learner highlights; the user now highlights supplement blocks, so insert_supplements.py refuses -> use the targeted
  script). The reader server at 8787 is started on demand by the library (8840); start_reader.py --port 8787 if down.
- 2026-11-05 later: user directive R11 (DESIGN_STEP2.md). New agents (outputs in ~/N4_BV_DR): freedom/ (classification of
  renormalization freedom, meaning of Q_1 etc.), literature/ (critical high-loop DR study), theory/ (BV formulation of
  DR for N=4 with constant ghosts, one-loop Q_1/K_1/D_1 plan). Survey of the user's earlier projects -> step2/next/survey/.
  Twist agent told to use the user's map (All-loop article Appendix B, reader :8778 from b-01003).
  Diagnosis given to user: Step-1 Q^GF violates Costello axioms (3),(5) (Ch.3) -> composite-line layer; survey checks
  whether earlier projects used a Q^GF with (Q^GF)^2 = 0.
- 2026-11-05 evening: user R12 -> rebuild on (Q^GF)^2 = 0. Stopped step2/next/action2, operator2 (old-GF data kept).
  Gauge-fixing design agent -> ~/N4_BV_HK2/gf/ (formulations A non-minimal (cbar,b) / B first-order; machine-checked axioms;
  equivariant extension; [x.,Q^GF] per generator; migration plan; draft revised Chapter 3). Twist agent continues.
- 2026-11-05 twist agent DONE (step2/next/twist/: twist.tex/pdf, results.json, src.nosync/qa2.py for Upsilon_xixi):
  Phi = projection onto single-letter H(Q + xi. + 1/2 xi xi^T), BK normalisation via the user's All-loop map;
  class of Q_1^xi mod Im Q_0^xi agrees with BK (3.64)-(3.69) in all 8 pair channels with ONE common constant (-i/4 engine
  units); log(3/2) cancels; Konishi entry nonzero (arity-zero kernel); needs arity-0, arity-2 (Upsilon_xixi) and the
  lambda* term of b. Manuscript Table c11-twistentries incomplete as a class statement (fix in rebuild).
  PROVISIONAL: old gauge fixing; engine B not run on arity-0/2; absolute constant open. -> gate for the rebuild.
- 2026-11-05 DR literature study DONE: ~/N4_BV_DR/literature/ (dr_survey.txt/json 174 works, synthesis.txt + PDF in build/).
  Key: no ">=6-loop" spectral number is a regulated-QFT computation (QSC, bootstrap, 4d integrands); DRED SUSY identities
  verified only partially (1-4 loops, beta-function equalities, few identities); NO composite-operator SUSY Ward identity
  ever evaluated in DRED/SSDR/BMHV; SM 4-loop g_s term reading-point dependent, R=3 via Weyl consistency. Q1-Q9 testable
  BV-DR questions in synthesis.txt (forwarded to the BV-DR theory agent).
- 2026-11-05 BV-DR theory note (first pass) DONE: ~/N4_BV_DR/theory/bv_dr.tex/pdf (25 pp), engine_spec.json (26 gates,
  routes A residues/QAP, B direct Green functions). DRED (Stoeckinger Q4S, V10 -> Q4S + W_perp) and BMHV-type (J = delta,
  keeps SO(4)_R); non-minimal (cbar,b), Psi = int(cbar, d.A + xi_GF b) Feynman at xi_GF = 1; QME = CME in DR, defect =
  evanescent breaking (QAP); D_1 = 2 gamma via RG; K_1 = -(v1,x_O) D_1 + K_1^0; relations fail exactly by Ob - sF_1.
  Agent resumed to add the literature Q1-Q9 targets. Downloaded papers in theory/src_pdf (never push).
- 2026-11-05 freedom note DONE: ~/N4_BV_DR/freedom/freedom.pdf (10 pp, 17/17 checks). Ch12 "invariant relations" with
  log(3/2) are NOT invariant; Thm c8-invariants (iii) understated (whole D_1 matrix canonical); Thm c12-opalg (iv) fixed by
  Prop 4.3. Corrections agent -> ~/N4_BV_DR/manuscript_fixes/fixes.json (apply, then reader c8/c11/c12 + GitHub).
- 2026-11-05 BV-DR theory final (27 pp; invariants I1-I6 from the freedom note; Q1-Q9 targets; 29 gates; in DR the
  one-loop BRST correction on operators is zero, so the twisted Konishi entry sits in Q_1 + arity-two). DR engines
  launched: ~/N4_BV_DR/engines.nosync/route_a (residues/QAP) and route_b (direct Green functions, no QAP).
- 2026-11-05/06 NEW GAUGE FIXING DONE: ~/N4_BV_HK2/gf/ (gf.pdf, ch03_revised.tex, checks.py 130/130, migration_plan.json
  84 items). Formulation A: antighost cbar, cbar* (weight 2), S' = S - 1/4 int(cbar*,cbar*), Psi = int(cbar, d*A|_V);
  Q^GF: A*->-2A*, lambda*->-dslash lambda*, cbar*->cbar* (into c), c*->-c* (into cbar); Feynman gauge forced by (Q^GF)^2=0;
  all 5 axioms; old Q^GF = compression via chain map j: (cbar, A*) -> A* - d cbar => arity-zero results (I^CT, 27/16, Ob3,
  F3) are pullbacks. Composite lines vanish along Delta, P, M, R; remain along K (4 iota_b A_V-c) and Q, S (lambda-cbar,
  identical to old). Recompute: arity one along Delta, P, M, K; theta, t_x; D1/P1/M1/K1 representatives; arity two with
  V-part brackets; log T along K; Theta_n on cbar. Referee -> ~/N4_BV_HK2/gf_review/; engine A -> ~/N4_BV_HK2/engine_a/.
