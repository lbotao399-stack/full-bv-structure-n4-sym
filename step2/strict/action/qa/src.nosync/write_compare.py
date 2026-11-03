import json, os
QA = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
L = 'log(3/2)'
E = []
def e(item, A, B, QAres, status, note=''):
    E.append(dict(item=item, engine_A=A, engine_B=B, qa_third=QAres, status=status, note=note))

e('G-reg0, two legs: Ob^(0)|_2', '-27/8 u int(d*dc,d*A_V)', '-27/8 u int(d*dc,d*A_V)', '', 'agree',
  'regression to Thm c4-obstruction / Ch9 table; exact in both')
e('G-reg0, two-leg special-edge pole', '-9/2 eps^-1 int(dc,A_V) (= -Q I^CT_2)', '(9/4) u eps^-1 Q int(A,A)',
  'consistent: Q int(A,A) = -2 int(dc,A_V) in both normalisations', 'agree')
e('G-reg0, three legs: Ob^(0)|_3', '(11/4 + 19/4 L)(G1+G3) + 27/8 G2 + 0 G4 - 9 L G5, L = log(3/2)',
  'same; parity-odd coefficient 0', '', 'agree', 'regression to Thm c9-ob3 (eq c9-ob3)')
e('T1: hbar Ob^(1)(Delta)|_2', 'u[-11/3 int(dA,dA) (A_V and A_perp) + 0 int(d*A_V,d*A_V) + 47/3 int(lam,dslash lam) - 14/3 int(dc,A*)]; position part 0; 5 further sectors zero',
  'same; parity-odd (lam,dslash rho(vol) lam) coefficient 0; position part 0', '', 'agree',
  'regression to Thm c9-quadratic (ii) reproduced by two independent lineages (A: two-letter engine A + Part-C mode model; B: Step-1 engine B vertices + ce_core/ce_int)')
e('G-pole along Delta (two legs)', '-9 int(A,A) = (9/4){int(A,A),S(Delta)}', '-9 u int(A,A); 0 in (lam,lam), (c,A*)', '', 'agree')
e('G-P2: Ob^(1)(P_b)|_2', '0 (direct P_b cycles and position part of Delta)', '0 (position part of Delta cycles)', '', 'agree')
e('X-scal (cross-check only): Ob^(1)(Delta)|_2 = 4 lc^(0) + composite',
  'lc^(0)|_2 = u[-1, 4, -1] from plain cycles; remainder u[1/3,-1/3,-2/3]',
  'special-edge placements u[-4,-4,16,-4] = 4 lc^(0); line placements u[1/3,1/3,-1/3,-2/3]',
  'Q N with engine-B Q and N of Thm c4-log, (lam*,lam) = (lam,lam*) (lam* even): Q N = [-1, 4, -1] = u^-1 lc^(0) of engine A (qa/out.nosync/q_lcN.txt)',
  'agree', 'B states (lam*,lam) = -(lam,lam*); with that reading Q N|_lam = -4, contradicting A; see findings (minor, no effect on B numbers)')
e('T4: Ob^(1)(xi_+-)|_2', '(79/24 - 1/2 %s)(g1+g3) + 27/8 g2 + 0 g4, both chiralities' % L, 'same, both chiralities; c-lam* sector 0', '', 'agree')
e('G-pole along xi', '(9/2)[int(Gamma_V,A_V) + int(Gamma_perp,A_perp)] = (9/4){int(A,A),T1+T2}', '(9/2) u int(A,Gamma(xi,lam)) = (9/4)u{int(A,A),S(xi)}', '', 'agree')
e('S1 along xi: Q images', 'QT1 = -QT3 = -1/2 g1, QT2 = -QT4 = -1/2 g3, QT0 = -g2-g4; {int(d*A_V,d*A_V),T1} = 2 g2', 'same; Q g2 = -Q g4; {F1^(0),T1} = 27/8 g2',
  '', 'agree', 'B found two Koszul rules passing its G-sign; A (Q = {S_free,-} in the mode model) fixes q1 = -1/2, i.e. B\'s chosen left rule')
e('S1 along xi: two-leg cohomology', 'H^0 = 0; ker Q at degree -1 = span{T1+T3, T2+T4}; no exact', 'same (dims 5/4, rank 3, H^0 0, H^-1 2)', '', 'agree')
e('S2 along xi: residual and particular solution', 'R^(1)(xi)|_2 = (79/24 - 1/2 %s)(g1+g3); F1^(1)(xi)|_2 = (%s - 79/12)(T1+T2) + u(theta^xi_V S^V_2 + theta^xi_perp S^perp_2), r0 = 0' % (L, L),
  'same', '', 'agree', 'theta^xi_V, theta^xi_perp remain free (not fixed by any computed equation)')
e('G-WZ1 along xi (two legs)', 'holds as Lambda Ob^(1) = -{Ob^(0),S(xi)} (378i vs -378i); fails with printed sign', 'holds with the same sign (residual 0 vs 6156)',
  'derived from Def c2-dce: (sG)^(1)(x) = Lambda G^(1)(x) + {G^(0),S(x)} for |G| = 1; toy BV algebra (q_bv_toy.py) confirms the derivation rule Lambda{a,b} = -{Lambda a,b} - (-1)^{|a|}{a,Lambda b}',
  'agree', 'theory note Lemma at-sformulas (ii) wrong (finding F1)')
e('B1: Ob^(2)(Delta,P_b)|_2, classes', 'o_perp = o_lam = o_V = -2/3 u; no other class; position part 0; no poles', 'o = -2/3 u in every block; other class (V block) 0; logs cancel in the class', '',
  'agree', 'theta^Delta_d = -1/3 (d = perp, lam, V); X-theta agrees (cross-check only). Normalisation pinned independently of X-theta by the corrected WZ2 along (Delta,xi) (below)')
e('B1: explicit representatives', 'perp -2/3 S(P_b)_perp; lam 2/3 int(d_b lam,lam*) + (-3/2+log(4/3)) int(rho(b) dslash lam,lam*); V -2/3 int(d_bA_V,A*_V) + (3/2-log(4/3)) int(d_mu(b.A_V),A*_mu) + (-5/6+log(4/3)) int(d_b c,c*)',
  'not reported (classes only, p<=1 bases: V 20/9/7/2, perp 3/1/0/1, lam 7/3/2/1)', 'QA arithmetic: A representatives = -2/3 S(P_b)_d + (3/2 - log(4/3)) Q int(iota_b A_V,c*) [V], + (-3/4+log 2-1/2 log 3) Q int(rho(b)lam*,lam*) [lam]: consistent with A\'s exact directions',
  'A only (representatives); classes agree')
e('T3: Ob^(1)(R_r)', '0 on 7 two-leg sectors and 3 three-leg sectors', '0 on 4 two-leg and the same 3 three-leg sectors', '', 'agree')
e('T3: Ob^(1)(M_omega)|_2, V sector', '1/3 Q int Om_{mu nu}(A_mu,A*_nu)', '+1/3 u Q int omega_{mu nu}(A_mu,A*_nu) (fitted, 3 configs)', '', 'agree')
e('T3: Ob^(1)(M_omega)|_2, lambda sector', '-1/12 int(rho(Om) lam, dslash lam), rho(Om) = sum Om_{mu nu} s[mu] sb[nu]',
  'asserted +1/3 u Q int((1/4) omega rho rho lam, lam*) but never fitted (raw lam-lam values stored, candidate list empty)',
  'replay of B\'s seed-23 stream (reproduces B\'s stored V candidates exactly) and fit of B\'s raw lam-lam values: theory = (1/3) Q(spinor spin part of S(M_omega)) exactly in 3/3 configurations (qa/q_mspin_replay.json)',
  'agree in magnitude (1/3 x 1/4 = 1/12); B claim now substantiated', 'sign map between A\'s s[mu]sb[nu] basis and B\'s RB[mu]R[nu] basis not established; both engines give the uniform statement F1^(1)(M_omega)|_2 = 1/3 u S(M_omega)_spin + closed in the V sector, QA confirms it in the lam sector for B')
e('T6: Ob^(1)(eta_+)|_2 intrinsic', '(A_V,lam): -5/2, -5/2; (A_perp,lam): -3/2 - 3 log2 + 3 log3 = -3/2 + 3 log(3/2) [results.json; notes/summary print -3/2 - 3 log(3/2)]; (c,lam*): 5/2; position part nonzero/nonzero/zero',
  'not computed',
  'second evaluation with engine B code (q_eta.py, q_eta_fit.py): S(eta) = sum_mu x^mu S(rho(e_mu)eta) + int s0, sign of s0 fixed by Q S_2(eta) = 0 (residual 0 for the printed sign of c11-Seta, 356 for the opposite); intrinsic coefficients on A\'s bases realised in B\'s conventions: -5/2, -5/2; -3/2 - 3 log2 + 3 log3; 5/2; no poles; position part nonzero in the A_V and A_perp sectors, zero in (c,lam*)',
  'agree (A and QA)', 'intrinsic parts only (canonical density, last leg eliminated); no solution along eta, no S1-eta cohomology; eta_- not recomputed')
e('T2: Ob^(1)(Delta)|_3 on classical cubic structures', '8/3 (four ghost structures); 116/3 + 6 log2 - 3 log3 [summary prints log(8/27)]; 116/3 + 8 log2 - 4 log3 [summary prints log(16/81)]; (A_V)^3, A_V A_perp^2 unfitted; position parts nonzero',
  'not computed', '', 'A only, partial', 'intrinsic parts with nonzero position parts: convention dependent, not usable as cohomology data (see F2 lesson)')
e('T5: Ob^(1)(xi)|_3', 'seeded tables; no poles', 'not computed', '', 'A only, tables')
e('T8: lc^(1)(xi)|_2, G-logex', 'lc^(1)(xi)|_2 = {uN,S_2(xi)} = 2u(T1+T2-T3-T4); lc^(1) - {uN,S} = 0; l_xi = 0, c = a = 0',
  'not computed', 'consistent with the corrected closed combination lc^(1)(x) - {uN,S(x)} (Def c2-dce + C-bracket rule)', 'A only')
e('B2: Ob^(2)(Delta,xi_+-)|_2', 'intrinsic (p = 0): (8 + log(4/3)) T1 + 8 T2 - 26/3 (T3+T4) + 0 T0, both chiralities; flagged: G-WZ2 fails',
  'not computed',
  'third evaluation with engine B\'s cycle code (q_dxi.py, both chiralities, 65 equations each): intrinsic part identical to A (T1 8 + log(4/3), T2 8, T0 0, T3 -26/3, T4 -26/3); poles 0; position part nonzero, purely log(4/3); rational position part 0; the whole log(4/3) part (intrinsic + position) equals log(4/3) Q[int x^mu (Gamma(xi,lam)^mu,c*) + int x^mu (A*_mu,(xi,lam*))] exactly (rank 4, unique). Hence hbar Ob^(2)(Delta,xi)|_2 = u[8(T1+T2) - 26/3(T3+T4) + log(4/3) Q(psi1+psi2)]',
  'agree on intrinsic part; A\'s G-WZ2 failure resolved (artefact)',
  'corrected WZ2 (Q Ob^(2)(Delta,xi) = Ob^(1)([Delta,xi]) - {Ob^(1)(Delta),S(xi)} + {S(Delta),Ob^(1)(xi)}) holds: Q Ob^(2)|_2 = -25/3 u (g1+g3) = -{Ob^(1)(Delta)|_2,S_2(xi)}; printed sign gives +25/3 (fails)')
e('G-tDelta (two legs): t_Delta from (Delta,xi)', 'not computed (value flagged)', 'not computed',
  'R^(2)(Delta,xi)|_2 at theta = 0 = Ob^(2) - {F1_part(Delta),S_2(xi)} = -1/3 u S_2(xi) + log(4/3) u Q(psi1+psi2), with F1_part(Delta)|_2 = u[22/3 int(A,A*) - 47/3 int(lam*,lam) - 8/3 int(c*,c)] (Def c11-free) and {F1_part(Delta),S_2(xi)} = 25/3 u (T1+T2-T3-T4); hence rho(Delta,xi) = -xi/3, t_Delta = -1/3 in the V and perp sectors and both chiralities; with theta^Delta_d = -1/3 the (Delta,xi) equation is solved on two legs by F1^(2)(Delta,xi)|_2 = log(4/3) u (psi1+psi2)',
  'QA result: pass', 'agrees with t_Delta = -1/3 from (Delta,P_b)')
e('(xi,P_b) two-leg class', 'not computed', 'not computed', 'position part of the (Delta,xi) cycles = 2 Ob^(2)(P_b,xi)|_2 (Prop at-pos) is log(4/3) times a Q-exact functional: class 0, but Ob^(2)(xi,P_b)|_2 != 0 (needs F1^(2)(xi,P_b) prop. log(4/3))', 'QA inference (via Prop at-pos), not a direct computation')
for pair in ['(xi_+,xi_-)', '(xi_+,xi\'_+)', '(xi,eta)', '(eta,eta\')', '(K_b,P_b\')', '(M,P), (M,xi), (R,xi), (Delta,eta), (Delta,K), ... (pairs absent from the theory table)']:
    e('arity two %s' % pair, 'not computed', 'not computed', '', 'neither')
e('T7: Ob^(1)(K_b) (finite T), S1-K, G-T', 'not computed', 'not computed', '', 'neither')
e('S1 along eta, solution along eta', 'not computed', 'not computed', '', 'neither')
e('G-H (Hypothesis at-H at three legs), S1-cubic, GC13, G-zero', 'not computed', 'not computed (B\'s gate named G-H is the two-leg S1 along Sigma)', '', 'neither')
e('G-log (numerical finite-eps extrapolation)', 'not run', 'not run (structural)', 'Prop at-small (ii) closed forms verified by Monte Carlo, n_t = 2..5 (n_t = 4 at 4e-4 relative, 40M samples)', 'neither (closed forms QA-checked)')
e('G-pos', 'structural (by construction)', 'structural (beta = cycle with S(d_mu v))', '', 'untested as an independent identity')
out = {'meta': {'role': 'QA comparator, action side', 'date': '2026-10-30', 'unit': 'u = -g^2 N/(32 pi^2)',
                'qa_code': 'qa/src.nosync: q_bv_toy.py, q_mspin_replay.py, q_dxi.py, q_dxi_fit.py, q_eta.py, q_eta_fit.py, q_lcN.py, q_small.py, write_compare.py, write_findings.py',
                'qa_data': 'qa/q_dxi_fit.json, qa/q_dxi_fit_minus.json, qa/q_eta_fit.json, qa/q_mspin_replay.json, qa/out.nosync/'},
       'entries': E}
json.dump(out, open(os.path.join(QA, 'compare.json'), 'w'), indent=1)
print(len(E), 'entries')
