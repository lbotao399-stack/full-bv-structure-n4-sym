# Step 1 design spec — N=4 SYM, Costello formalism, UV boundary only

Project root: /Users/libotao/Desktop/N4_Costello_BV_2026-10-27
Notation ledger (READ FIRST): sources/ledger_SYNTHESIS.md (+ per-source ledgers in sources/).
Primary-source texts already downloaded (read-only):
- [C07] Costello arXiv:0706.1533v3: /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/costello0706/c.txt (and .pdf)
- [CG2] Costello–Gwilliam Vol. 2 draft: .../scratchpad/cg2/cg2.txt
- [ESW] Elliott–Safronov–Williams arXiv:2002.10517: .../scratchpad/esw_pdf/esw.pdf
- [SW] Saberi–Williams arXiv:1910.04120: .../scratchpad/sw_src/page.txt
(prefix: /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407)

## Hard rules (user)
1. Primary sources only. Do NOT read, cite, reuse or compare with any other local project/file on this machine (all earlier N=4 / QED / SQED / supergraph work is excluded). Standard textbook facts may be used only as sanity checks, never as inputs.
2. No new notation. Every displayed symbol must come from the ledger (with its source) or be written out in full. Spinors index-free: S_+, S_-, Sigma, W, W^*, Gamma, rho, (-,-). Explicit component bases (gamma matrices, spinor indices, Lie-algebra bases) live ONLY inside code; reported formulas are index-free.
3. Euclidean R^4, complexified fields, Weyl spinors.
4. UV boundary only. Propagator P(eps,infty) = \int_eps^infty (Q^{GF} (x) 1) K_t dt. On R^4, K_infty = 0. Never introduce an IR scale (no T, no L, no masses as regulators). External momenta generic and nonzero when evaluating graphs.
   AMENDED by the user 2026-10-28: the standard Costello transfer endpoint T IS allowed (P(eps,T) = \int_eps^T (Q^{GF} (x) 1) K_t dt, effective action and observables at scale T, [C07]); it is the endpoint of the transfer, not a physical IR regulator. Requirements: every core UV quantity (counterterms, obstructions, finite counterterms, local representatives L_1(O) and the components Q_1, S_1, D_1, K_1, P_1, M_1, R_1) must be independent of T and T must not enter the core theory; wherever T -> infinity is safe (proved), take it (as in Step 1: P(eps,infty), K_infty = 0). Still forbidden: masses, position cutoffs R or any other IR regulator in a definition (R only as a diagnostic).
5. No dimensional regularization in any role.
6. Machine: 24 GB Mac shared with the user. Each process < 2 GB; at most 2 heavy processes at a time per agent; never write to /tmp (use the project folder); kill every process you start; Wolfram: at most one kernel, check `pgrep -fl WolframKernel` first, never kill a kernel you did not start.
7. Deliver derivations and exact results; no fitting to expected answers.

## Formalism conventions (fixed)
- Weight e^{+S/hbar}; S = (1/2)<e,Qe> + I, real slice with <x,Qx> negative definite [C07 CS1].
- QME: QI + (1/2){I,I} + hbar Delta I = 0; Delta_eps = -partial_{K_eps}; partial_phi = (1/2) sum ... [C07 §6.4]; [partial_phi, Q] = partial_{Q phi}; Q P(eps,infty) = K_eps - K_infty = K_eps [C07 Lemma 10.2.1 at T=infty].
- Flow: W(P, I) = hbar log( e^{hbar partial_P} e^{I/hbar} ) [CG2 §7.2.7] (= [C07] Gamma(P,S)).
- Counterterms: I - I^{CT}[eps], I^{CT}[eps] in O_loc (x) A_{<0}, one loop A_{<0} = span{eps^{-1}, log eps} [C07 Thm B, CG2 §13.7.2 Sing_eps]. Constructed as the singular part of the eps -> 0 expansion, graph by graph (locality is checked, not assumed).
- Bracket {f,g} [C07 Lemma 4.3.2 = ESW Def 1.15]; scale-eps bracket {,}_eps defined from Delta_eps by the [CG2 §7.2.8] product rule.
- Observables [CG2 Def 8.3.1.1, 8.3.1.3]: square-zero parameter delta of degree -|O|. For a local O, apply [C07 Thm B] to the family I + delta O; write O(eps) := O - partial_delta I^{CT}[eps] (O(eps) is the user's own symbol). Then
  O[infty] := lim_{eps->0} partial_delta W(P(eps,infty), I - I^{CT}[eps] + delta O(eps)).
  At scale infty, Delta_infty = 0 and {,}_infty = 0, so the CG2 differential on observables is (d_CE + Q) O[infty] (d_CE = Chevalley–Eilenberg differential of the symmetry algebra when global ghosts are present [CG2 Prop 13.4.0.1, ESW Def 1.26–1.27]).
- UV-boundary identity: (d_CE + Q) e^{hbar partial_{P(eps,infty)}} = e^{hbar partial_{P(eps,infty)}} (d_CE + Q + hbar Delta_eps). Every anomaly / every singular coefficient therefore sits at t = eps.

## Fields (ESW, dimensionally reduced 10d SYM; index-free)
- 4d N=4 = dimensional reduction of 10d N=(1,0) SYM [ESW §10.3, Prop 2.20, §2.6]. Keep ESW's own 10d symbols: the gauge field A is a section over R^4 of V_10^vee (x) g with V_10 = V (x) ... = V (+) W_R^perp (W_R^perp = 6 transverse directions); its V-part is the 4d connection, its W^perp-part is the six scalars (write them as iota_w A for w in W_R^perp, or "the W^perp component of A"). Gaugino lambda in Pi Sigma_10 (x) g, Sigma_10 = S_+^{10} -> S_+ (x) W (+) S_- (x) W^* (ESW §10, W of dim 4 for N=4).
- Action: ESW (31) reduced: S_gauge = \int dvol( -1/4 (F_A,F_A) + 1/2 (lambda, dslash_A lambda) + (d_A c, A^*) + ([c,lambda],lambda^*) + 1/2([c,c],c^*) ), F_A = dA + 1/2[A ^ A] with d only along R^4. (-,-) is the nondegenerate invariant pairing on g; coupling sits in (-,-) (gap G12 option a). Color factors in results: write them with [-,-] and (-,-) only, e.g. tr_g(ad_X ad_Y) written in full; no new symbol for a Casimir.
- Antifields: A^*, lambda^*, c^* [ESW Rem 1.36]. Degrees: c in -1, A, lambda in 0, A^*, lambda^* in +1, c^* in +2 [ESW Rem 1.33].

## Gauge fixing (user decision 2026-10-27): second-order, Feynman-gauge Q^GF
Minimal complex (c, A, A^*, c^*, lambda, lambda^*). Q^GF: A^* -> A (identity via the metric on V_10, all ten components), A -> d^* A (4d part only, into the c slot), c^* -> d c^* (into the A^* slot), lambda^* -> (Dirac adjoint) so that H = [Q,Q^GF] = (positive Laplacian) (x) id on EVERY component. Signs to be fixed so that the A–A covariance is positive. Facts to prove (by machine and by hand): H Laplacian; [Q^GF,H]=0; P(eps,infty) graded symmetric; Q P(eps,infty) = K_eps; (Q^GF)^2 != 0 (state explicitly: this is the one axiom of [C07 Def 4.3.3] that fails; Im Q^GF is not Lagrangian). Ghost line = c – A^* propagator (no antighost).

## Superconformal symmetry
g = psl(4|4) (complexified Euclidean, [SW (7),(8),(17)], [ESW Def 1.26–1.27]). Global ghosts = coordinates on g[1]; components written as cochain evaluations S^{(k)}_g(X_1,...,X_k) [ESW Def 1.27]. Dilatation generator: fix ONE normalization from [SW (17)] (state weights) and keep it. Component of a differential along a generator X: evaluate the cochain on X, e.g. [(d_CE + Q)O[infty]](X).

## Step-1 computations
(a) All one-loop 1PI two-point graphs with superficial degree >= 0 (A along V, A along W^perp, lambda, c–A^*, and antifield two-point sectors such as (A^*, c)); exact eps -> 0 singular parts (eps^{-1}, log eps) as local polynomials in external momentum; I^{CT}[eps] on the quadratic sector, index-free.
(b) Scalar bilinears O = ( iota_{w1} A, iota_{w2} A ) (gauge invariant, w1,w2 in W_R^perp; decomposes into traceless symmetric and trace) and single letters; all one-loop graphs with one O insertion and degree >= 0; partial_delta I^{CT}[eps] -> O(eps).
(c) Dilatation component [(d_CE+Q)O[infty]](Delta) at one loop, two routes: (i) scaling route via eps∂_eps of the counterterms; (ii) direct UV-boundary graphs. Gauge-invariant results must be gauge-parameter independent in the sense available (Feynman gauge only here) and the traceless symmetric part must come out 0.
(d) Selection rules for all other superconformal components (translations, rotations, R: regulator-covariant -> 0; Q, S, K: allowed output structures by weight/R/Lorentz/ghost degree).
