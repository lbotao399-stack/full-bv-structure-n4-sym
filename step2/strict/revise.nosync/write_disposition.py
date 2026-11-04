"""Reviser: disposition of every QA finding (action QA F1-F12, operator QA QA-F1..QA-F15)."""
import json, os
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
A = 'action/qa/findings.json'; O = 'operator/qa/findings.json'
D = []
def add(src, fid, sev, decision, check, action, status):
    D.append(dict(source=src, id=fid, severity=sev, decision=decision, own_check=check, action_taken=action, status_after_revision=status))

add(A,'F1','major','accept',
 'Re-derived by hand: Def c1-bracket gives {f,g} = Delta0(fg) - Delta0(f)g - (-1)^|f| f Delta0(g) = (-1)^|f| {f,g}_std, with {,}_std the standard odd bracket. '
 'A degree-1 derivation D of {,}_std (Q, {I,-}, ad_F for |F| = 0) satisfies D{a,b}_std = {Da,b}_std + (-1)^{|a|+1}{a,Db}_std, hence D{a,b} = -{Da,b} - (-1)^{|a|}{a,Db} for the bracket of Def c1-bracket. '
 'Forcing s^2 = 0 on G = sF (F = F^(0) only, then F = F^(1) only) gives (sG)^(1)(x) = Lambda G1(x) + {G0,S(x)} and (sG)^(2)(x,y) = Lambda G2 - G1([x,y]) + (-1)^{|x|}{G1(x),S(y)} - {S(x),G1(y)} + {G0,S2(x,y)}, '
 'i.e. the WZ identities Lambda G1 = -{G0,S(x)}, Lambda G2 = G1([x,y]) - (-1)^{|x|}{G1(x),S(y)} + {S(x),G1(y)} - {G0,S2}. The degree-0 formulas (eq at-eqs) are unchanged. '
 'Closed log combination: Lambda(lc1 - {uN,S(x)}) = -{u Lambda N,S} + {u Lambda N,S} = 0. Rerun here: qa/src.nosync/q_bv_toy.py (C rule True, standard rule False; s(sF) = 0 with the corrected formula, fails with the printed one).',
 'action_theory parts 2 and 3/4 rewritten: Lemma at-sformulas (ii),(iii) and proof, Prop at-log, Rem at-logexact, Thm at-eq1 (iv),(v), Prop at-xi2 (v),(vi), Thm at-eq2 (ii),(iii); derivation rule stated in the proof; recompiled.',
 'fixed')
add(A,'F2','major','accept',
 'Confirmed by the independent action-completion run (action/complete, both routes, results Ob2(Delta,xi+-)): solvable exactly for theta^Delta_d = -1/3 in every sector; F1^(2)(Delta,xi)|_2 has the log(4/3) position structures in both routes; G-pole2, G-F1Delta pass. The engine-A gate demanded closedness of Ob2 and dropped the position part.',
 'Flag lifted; engine-A gate G-WZ2 (Delta,xi) replaced in gates_revised.json by the corrected full-functional gate (pass, from action/complete and qa/q_dxi_fit.json). Reported in action_final.json.',
 'fixed')
add(A,'F3','major','accept in part',
 'action/complete computed: S^(2) vertex, (xi+,xi-), (xi,xi\'), (xi,P_b), (Delta,xi), eta+- arity one with cohomology, (xi+,eta+-), K_b arity one with cohomology; 73 gates, two routes. Still not computed: (K_b,P_b\'), (eta,eta\'), (xi-,eta+-), all three/four-leg classes (S1-cubic, Hypothesis at-H, GC13), cubic couplings, G-log numerically, G-T numerically, the 4 extra closed directions of F1(K_b)|_2.',
 'Completion results taken into action_final.json and into action_theory (new Section: computed solutions and relations). Remaining items listed as open.',
 'partly fixed; remainder open')
add(A,'F4','minor','accept',
 'Checked engine_a/out.nosync/task_T6.json: (AW,lam) coefficient "-3/2 + (-3)*log2 + (3)*log3" = -3/2 + 3 log(3/2); notes.txt line 76 prints -3/2 - 3 log(3/2). T2: notes give 116/3 + 6 log2 - 3 log3 = 116/3 + log(64/27) and 116/3 + 8 log2 - 4 log3 = 116/3 + log(256/81); the earlier summary printed log(8/27), log(16/81).',
 'engine_a/notes.txt corrected (T6 line); compact forms added.', 'fixed')
add(A,'F5','minor','accept',
 'qa/q_mspin_replay.json: theory/candidate = 1/3 exactly on 3/3 lam-lam configurations of engine B data. Engine A value -1/12 int(rho(Om)lam,dslash lam) is 1/3 x (1/4) in magnitude, consistent with S(M)_spin lambda part (1/4) omega rho rho lam; the realisation sign map between the two Clifford realisations is not needed because each route fits within its own realisation.',
 'engine_b/notes.txt: lambda-sector fit recorded; action_final.json states F1part(M_omega) = 1/3 u S(M_omega)_spin in all sectors.', 'fixed')
add(A,'F6','minor','accept',
 'Rerun here: qa/src.nosync/q_lcN.py with eb_qb.apply_Q gives Q int(lam,lam*) = -int(lam,dslash lam) and, with (lam*,lam) = +(lam,lam*) (lam* even, Def c1-pairing), Q N = [-1, 4, -1] = engine-A plain-cycle lc0/u.',
 'engine_b/notes.txt G-sign convention text corrected.', 'fixed')
add(A,'F7','minor','accept',
 'The inner direction shifts F1(P_b) by 2c u S(P_b) and F1(Delta) by 0, so translation covariance of F1 (Prop at-pos applied to the solution) identifies t_P = t_Delta. action/complete computed (xi+,eta+) and (Delta,xi): t_M = t_Delta, t_R = t_Delta + 1/3, t_xi+ + t_eta+ - t_Delta = 27/4 + 2L, t_xi+ + t_xi- - t_P = 27/4 + 2L, sector uniformity of t_P, t_xi, t_eta.',
 'Thm at-fix table completed (M_omega, R_r, (Delta,xi), (xi,P_b), (xi,eta)); covariant solution set stated (outer direction only); computed relations inserted.', 'fixed')
add(A,'F8','minor','accept',
 'The omega-linear (lam,lam)_1 space at weight 4: rho^(3) = rho_[mu rho_nu rho_rho] with Gamma symmetric makes int(lam, omega rho^(3) d lam) a total derivative, so the space is one-dimensional, spanned by Q S(M)_spin|_lam; engine B fit (finding F5) gives the uniform 1/3.',
 'Prop at-two (iii) extended to the lambda sector.', 'fixed')
add(A,'F9','minor','accept',
 'A* is odd and the fibre pairing on V_10 (x) g symmetric: (A*_perp,A*_perp) = sum A*^a_I A*^a_I = 0 identically.',
 'Vacuous step removed from the proof of Thm at-theta; the zero leg sets marked in the text (tables are machine generated, entries annotated).', 'fixed')
add(A,'F10','minor','accept',
 'gates.json of both engines: B "G-H" is the two-leg S1 along Sigma; engine A G-logex was redefined after the sign finding; G-pos structural; G-log not run.',
 'gates_revised.json written for both action engines (B G-H renamed G-S1xi2; G-logex documented with the corrected sign; G-pos marked structural; G-log listed as not run); engine_spec gate texts G-WZ1, G-WZ2, G-logex corrected in theory note Sec. engines.', 'fixed (G-log open)')
add(A,'F11','minor','accept',
 'Intrinsic parts with the last-leg-eliminated density are not invariant under integration by parts; (Delta,xi) is an example (intrinsic log(4/3) on T1 is part of Q of a position functional).',
 'action/complete reports eta and K_b with position parts (canonical data alpha, beta, gamma). Three-leg Delta fits remain intrinsic-only and are marked as not usable for cohomology.', 'fixed for eta, K_b; open for three legs')
add(A,'F12','minor','accept',
 'Ob1(Delta)|_2 = u[-11/3, 47/3, -14/3] and Q int(A,A*) = -1/2 int(dA,dA) - int(dc,A*)|_V, Q int(lam*,lam) = -int(lam,dslash lam), Q int(c*,c) = -int(dc,A*) give F1part(Delta)|_2 = u[22/3 int(A,A*) - 47/3 int(lam*,lam) - 8/3 int(c*,c)] (A: -11/3 = 22/3 (-1/2); lam: 47/3; c: -14/3 = -22/3 + 8/3); recomputed by action/complete as gate G-F1Delta in both routes.',
 'Explicit two-leg solutions along Delta, P_b, xi, eta, K_b, M_omega, R written in action_theory (new section) and action_final.json.', 'fixed (three legs open)')

add(O,'QA-F1','major','accept',
 'Re-derived: a bridge = divergent sub-cycle (sum_m eps^-m a_m + a_0 + o(1)) times a special tree edge (e^{-eps|p|^2}-1) = sum_j (-eps|p|^2)^j/j!; its eps^0 part sum_m ((-1)^m/m!)|p|^{2m} a_m is nonzero, and the counterterm partner -{f, eps^-m ct_m}^sigma_eps supplies exactly its negative. Hence bridge + partner = (e^{-eps H}-1)(cycle - Sing cycle) -> 0 and the finite special vertex must not be added again. Data: engine A 646/646 D>=2 rows pass with [eps^0]D, all 68 F_x != 0 single-output rows fail with F_x added; engine B same 13 rows with F_x != 0; arity zero: the vacuum-pole term of eq op-master0 is the partner of the bridges through the action self-energy and cancels the same way.',
 'operator_theory: Def op-defect (iii)-(v) (Upsilon_x := [eps^0] D_x, Upsilon_0 := [eps^0] special-edge cycles), proof of Thm op-master (Br) step regrouped bridge+partner, eq op-master0 last term deleted, Prop op-master2 and the pole identity without j-shift terms; recompiled. Manuscript items (Def c11-defect (iii)-(iv), Thm c11-classes, eq c11-pole, eq c10-finitesp, qsk_summary A1, twoletter build_final) listed for the author: manuscript/ not edited (not permitted to this role). Numbers unchanged: both engines already use [eps^0] D; F0 arity-zero entries have no D>=2 pair and no antifield letter.',
 'fixed in step2/strict; manuscript propagation pending (author)')
add(O,'QA-F2','major','accept in part',
 'action/complete now provides F1(eta)|_2 (both chiralities, both routes) and F1(K_b)|_2 (position degree <= 2, both routes).',
 'One-body parts of S_1 and K_1 assembled here on F0 (revise.nosync/assemble_onebody.py): S_1(O)|_traceless = c_S {S^(1)(eta),O} with c_S = 79/12 - 3L - t_eta; K_1 one-body = -{F1(K_b)|_2,O} with the computed x-linear part and the unfixed closed directions explicit. T-independence of F1(K_b) is structural (eps^0 simplex coefficients, Lemma at-special), no numerical finite-T run. S_1, K_1 on F1-F7 and (R3), (R5), (R6) remain open.',
 'partly fixed; open beyond F0')
add(O,'QA-F3','major','accept',
 'Confirmed from both engines\' gates and open lists.',
 'No new operator engine runs in this revision (the operator engine A run is resumable, n=3 tables and F5 at S\'=4 pending). Listed in operator_final.json under open.', 'open')
add(O,'QA-F4','major','accept',
 'The transport map Phi needs the BK normalisation of the letters, the field map (3.51) and the representative of b, which no source gives (q1_targets open items); the (xi,xi) arity-two operator kernels (needed for the bc and bb entries, k = 2) are not enumerated by any engine.',
 'Class-level check that needs no normalisation of Phi done: Q_1 along a rank-one xi on the chiral-primary bilinears Tr gamma^(I gamma^J) = int omega (iota_w1 A, iota_w2 A) with w_a in the 3-plane of w vanishes identically (Q_0^xi O = 0 and (w1,w2) = 0 kill c_Q F1 and the F2, F3 terms), in every scheme: this is the n = 2 case of [BK, Sec. 3.4.1, (3.77)]. The full class comparison with [BK (3.7), (3.59)-(3.70)] remains open.',
 'open (partial class-level check passed)')
add(O,'QA-F5','minor','accept',
 'One global sign per route is calibrated (G-reg / G-pole); after calibration every other entry is a prediction (F1-F7 spectra, Heisenberg, su(1|1)).',
 'Recorded; derivation of the global factor from Def c11-defect (ii) and Lemma c10-sign not done.', 'open')
add(O,'QA-F6','minor','accept',
 'qa/out.nosync/cmp_tables.json conj_disagree: 11 rows, all with B has.EP = False or has.ESW = False; engine A complete there and the pair form holds.',
 'operator/engine_b/gates_revised.json: incomplete rows marked "no verdict"; B-F3 restricted to the (A_V,X)->X line-class rows; the 8 gates "with the theory finite vertex F_x" reclassified as tests of the withdrawn definition (evidence for QA-F1).', 'fixed')
add(O,'QA-F7','minor','accept', 'Both gates labelled structural in the engines.', 'Recorded; direct finite-eps Schwinger integration not run.', 'open')
add(O,'QA-F8','minor','accept', 'qa/out.nosync/recheck_rowsA.json: poles satisfy the Ward form on every entry without line placements, fail on all 204 entries with line placements.', 'Recorded; the full pole identity with ct^(1)(sO) and spectator terms remains to be computed.', 'open')
add(O,'QA-F9','minor','accept',
 'Along Delta a cycle with one Y = -P(eps,2eps) edge has no log eps, so the edge identity alone cannot produce 4 lc; the factor comes from the end actions.',
 'Mechanism sentence of Rem op-pair (2) removed; the remark now states only the edge identity and that the pair form rests on G-D1pair.', 'fixed')
add(O,'QA-F10','minor','accept (resolved by new action data)',
 'action/complete: the arity-two equations (xi+,xi+\') and (xi-,xi-\') force theta^xi_V = theta^xi_perp (both routes), and with it the eta sector difference is fixed (3L - 3/2 in the parametrisation of action/complete 2.4), so along Sigma and Sigma* every free direction of the solution acts as -t_x rho_0(x).',
 'Prop op-tfix: hypothesis now a consequence of the computed action-side equations; statement amended.', 'fixed')
add(O,'QA-F11','minor','accept', 'Engine B OP7 uses one numeric point on F0; exact matrix identities exist on F3, F4 (B) and F5, F7 (A).', 'Recorded.', 'open')
add(O,'QA-F12','minor','accept', 'Spectra are basis-normalisation independent; off-diagonal entries 8, 64 scale with the normalisation of Z\'F_11.', 'Recorded; a priori normalisation of Z\' (Def c11-families) not re-read here (ch11 is iCloud-evicted; the GitHub copy has a different size and was not used).', 'open')
add(O,'QA-F13','minor','accept', 'Only route-B OP8 configurations; action side K_b T-independence structural.', 'Recorded.', 'open')
add(O,'QA-F14','minor','accept (resolved)',
 'phi_c: action/complete gate G-F1Delta (both routes) uses F1part(Delta)|_2 = 22/3 int(A,A*) - 47/3 int(lam*,lam) - 8/3 int(c*,c) and reproduces Ob1(Delta)|_2 exactly, so phi_c = 8/3 (14/3 is the coefficient of int(dc,A*) in Ob1, not of int(c*,c) in F1). M_Om spin: 1/3 in all sectors (engine B fit + QA replay, finding F5).',
 'operator_final.json uses phi_c = 8/3 and the uniform 1/3.', 'fixed')
add(O,'QA-F15','minor','accept', 'Engine B notes 7.', 'Recorded; caps 1250/1400 MB kept for any rerun; no heavy process run in this revision exceeded 1.5 GB.', 'recorded')

out = dict(file='step2/strict/disposition.json', date='2026-11-04', role='reviser (strict route R10)',
           counts=dict(total=len(D), accepted=sum(1 for d in D if d['decision'].startswith('accept')), rejected=0),
           dispositions=D)
json.dump(out, open(os.path.join(ROOT, 'disposition.json'), 'w'), indent=1)
print(out['counts'])
