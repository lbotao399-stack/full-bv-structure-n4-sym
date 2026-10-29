# Step 2 (first installment) — binding spec

User request (2026-10-28): "Start Step 2: the equivariant one-loop QME along the superconformal ghosts via K_eps
boundary graphs, beginning with the dilatation Delta and the cubic gauge sector."

Read first, in this order: DESIGN.md (hard rules, unchanged), DEFINITIONS_v2.md (binding corrections), STORY.md,
sources/ledger_SYNTHESIS.md, then manuscript/ch03_uv_boundary.tex, ch04_counterterms.tex (section "The quantum master
equation at the UV boundary"), ch06_dilatation.tex, ch07_selection_roadmap.tex (section "The road ahead").
Hard rules recap: primary sources only ([C07],[CG2],[ESW],[SW]); NO earlier project on this machine except this
folder; no new notation (write objects out; reuse the manuscript's names); index-free in reported formulas; Euclidean;
UV boundary only (no T, no IR scale, no masses); NO dimensional regularisation in any role; exact results, no fitting.
Machine: 24 GB shared Mac. Each agent: at most ONE heavy process at a time, each < 1.5 GB RSS; never /tmp (use the
project folder); kill what you start; Wolfram: only if unavoidable, at most one kernel, check `pgrep -fl WolframKernel`
first, never kill a kernel you did not start. Subagents cannot write .md files: write .tex / .json / .py / .txt.

Unit everywhere: hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) (code-level name "unit"; never a new symbol in text).
Output folder: step2/ with subfolders theory/, engine_a/, engine_b/, numeric/, qa/. Engines A and B are independent:
neither reads the other's step2 outputs before the comparison phase. Both may import their own Step-1 code
(engine_a/, engine_b/).

## Part C — cubic sector of the arity-zero one-loop obstruction

Interaction J = I - I^CT[eps] (Step-1 I^CT, Theorem thm:c4-quadratic + thm:c4-log). Defect vertex at one loop
(prop:c4-defect1): D(eps) = 1/2({I,I}_eps - {I,I}) + hbar Delta_eps I - Q I^CT[eps] - {I, I^CT[eps]}_eps.
Q W(P(eps,infty),J) = connected graphs with exactly one D(eps) vertex, other vertices I, edges P(eps,infty)
(prop:c4-defectflow). Special edge = K_eps - K_0 (field-antifield components only), Schwinger time in [0,eps].

Sectors (degree 1, three legs, the complete list by weight/Lorentz/colour): (c, A|_V, A|_V), (c, A|_Wperp, A|_Wperp),
(c, lambda, lambda), (c, c, A*|_V). Weight 4 (two derivatives for (c,A,A); one for the others). Colour: only the
antisymmetric invariant (X_0,[X_1,X_2]) (tr of three ad's in the adjoint is proportional to it).

Basis of the target (each a local functional; (-,-) = pairing of g tensored with the metric; d^*d, dd^* on one-forms
of V, d^*d on the scalars iota_w A componentwise):
  G1 = int dvol ([c, A|_V], d^*d A|_V)
  G2 = int dvol ([c, A|_V], d d^* A|_V)
  G3 = int dvol ([c, A|_Wperp], d^*d A|_Wperp)
  G4 = int dvol ([c, dc], A^*)            (same pairing as the term (dc, A^*) of the action)
  G5 = int dvol ([c, lambda], dslash lambda)   (same pairing as the kinetic term 1/2(lambda, dslash lambda))
Report Ob_3 = unit * sum_i alpha_i G_i, alpha_i exact; prove/verify that the residual outside span{G_i} is zero by
evaluation on random plane waves (symmetric parts that vanish identically must vanish numerically too).

Graph classes (one loop, three legs, one D vertex):
  I-a  special edge ON the loop: triangles (three cubic vertices) and bubbles (one quartic + one cubic vertex), the
       other loop edges P(eps,infty). 1PI.
  I-b  special edge as a bridge between a cubic vertex with two external legs and a regulated one-loop two-point
       subgraph (bubbles + seagull, all edges P(eps,infty)). Pair it with the K_eps - K_0 part of -{I_3, I^CT_2}_eps:
       (e^{-eps|p|^2} - 1) * (renormalised self-energy) -> 0. Verify the cancellation explicitly.
  II-a -Q I^CT_3 (log part only; Step-1 cubic log coefficients).
  III  -{I_3, I^CT_2} (K_0 part) and its K_eps - K_0 part (paired with I-b).
  I-d, II-b  one-particle-reducible: two-leg defect subgraph (Step-1 special bubbles + -Q I^CT_2) joined by
       P(eps,infty) to a cubic vertex: these converge to the tree of the quadratic obstruction Ob_2 with I_3 (not part
       of the local Ob_3). Verify that their eps -> 0 limit is exactly that tree (no extra local remainder).
Expected checks: eps^{-1} total = 0 in each sector; log eps of I-a = Q Y_3 + {I_3, Y_2} = 0 where Y = log-coefficient
of I^CT (Y = s(unit*(2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c))) by thm:c4-log, s = Q + {I,-}).
Integrals: scaled Schwinger variables tau_1 in [0,1] (special), tau_2, tau_3 in [1,infty); the eps^0 coefficient is the
order-|p|^2 (resp. |p|^1) Taylor coefficient of the Gaussian-integrated integrand; prove convergence of each integral
used. Exact evaluation (rational + possibly log 2, log 3, pi^2, Li_2 values): exact symbolic integration preferred;
otherwise >= 40-digit numerics + integer relation on an explicit basis, confirmed by a second evaluation.

Exactness: finite counterterm F = F_2 + F_3 with F_2 = (27/16) unit int(d^*A|_V, d^*A|_V) + (Q-closed quadratic terms
beta_1 int(dA|_V,dA|_V) + beta_2 int(dA|_Wperp,dA|_Wperp) + beta_3 int(lambda,dslash lambda) + beta_4 int(dc,A^*)),
F_3 in span of
  H1 = int dvol (dA|_V, [A|_V ^ A|_V]),  H2 = int dvol (d A|_Wperp, [A|_V, A|_Wperp]) (one-form part of d on scalars),
  H3 = int dvol ([A|_V, c], A^*|_V),  H4 = int dvol ([A|_Wperp, c], A^*|_Wperp),  H5 = int dvol ([c,lambda],lambda^*),
  H6 = int dvol ([c,c], c^*).
Equation: Ob_3 = Q F_3 + {I_3, F_2}. Solve exactly; state the solution family (beta-directions = cubic parts of the
s-closed completions int(F_A,F_A), int(d_A phi, d_A phi), int(lambda, dslash_A lambda), BRST term). Check the
Wess-Zumino closure Q Ob_3 + {I_3, Ob_2} = 0 and state precisely whether it is automatic.

## Part D — the arity-one equivariant master equation along the dilatation

S_M := S^(1)(Delta) (ex:c2-dil), quadratic; M := v. for v = 2E (def:c7-vdot, lem:c7-weights): even, [M,Q] = 0,
M preserves <,>, [M,Q^GF] = -4 Q^GF, [M,H] = -4H, (M(x)1 + 1(x)M) K_t = -4 t dK_t/dt,
(M(x)1 + 1(x)M) P(eps,infty) = 4 eps L_eps = -4 eps dP/deps (prop:c7-dil). {S_M, -} = derivation induced by M (one
sign fixed by the bracket conventions: fix it explicitly in code and text).
Targets:
 D1 (theorem) Exact identities, valid at finite eps for every even J at least cubic mod hbar:
    (a) graphs with one S_M vertex: d/dxi W(P(eps,infty), J + xi S_M) in closed form (S_M quadratic, so the conjugation
        e^{hbar d_P} S_M e^{-hbar d_P} terminates): terms S_M, a constant, the derivation along e -> P*Me applied to
        W(P,J), and a second-order operator with kernel built from P* M P*; state each with signs;
    (b) Q of (a) via Q P* + P* Q = e^{-eps H} (K_infty = 0) and [Q,M] = 0;
    (c) [D_M, d_phi] = -d_{(M(x)1+1(x)M)phi}; D_M W(P,J) = W'_J(D_M J) + 4 eps d/deps|_J W(P,J) (check sign).
 D2 (theorem) The arity-one counterterm I^CT_{psl(4|4)}[eps](Delta) at one loop, the arity-one defect along Delta
    (ch07 "The road ahead", displayed formula), its eps -> 0 limit Ob(Delta) (modulo constants), and a finite local
    arity-one counterterm F(Delta) of degree -1 with Ob(Delta) = (Q + {I,-}) F(Delta) + {S_M, F_0} (F_0 = arity-zero
    finite counterterm). Expected: no trace-anomaly class (coefficient along int(F_A,F_A) proportional to the coupling
    log coefficient, which is 0), F(Delta) a multiple of unit*(2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c)).
    Derive, do not assume.
 D3 (theorem) The identification (eq:c6-step2) for every homogeneous local operator O at one loop, i.e.
    [(d_CE+Q) O[infty]](Delta) = 2 s d/ds|_{s=1} sigma_s(O[infty]) - Q(4 unit * tree insertion of
    (2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c)) into O[infty]) + (terms built from (Q+{I,-})O),
    or the corrected statement if it is false; hence the dilatation coefficients of thm:c6-dil as components of the
    equivariant cochain.
 D4 (numeric verification in N=4) For the quadratic sector: evaluate the one-loop graphs with one S_M vertex
    (S_M on each internal edge, on an external leg) for the two-point functions of A|_Wperp (and A|_V, lambda, c-A^*),
    with M in momentum space (E acts on Fourier components as -(4 + p.d/dp)), at several eps and momenta; compare with
    the closed forms of D1/D2 (singular parts, finite parts, eps -> 0 limits). Then the direct UV-boundary evaluation of
    the bilinear part of [(d_CE+Q)O[infty]](Delta) for O = int omega (iota_{w1}A, iota_{w2}A) (traceless and trace),
    i.e. graphs with one O vertex, one S_M vertex and one K_eps edge (ch06 figure (a)-(y)) plus the counterterm and
    defect graphs; compare with 2(E + w_O + coefficient) O: coefficients 0 (traceless), -24 (trace), in units.

## User directives added 2026-10-28 (binding for all later work)
R1. Working notation is free, but EVERY report to the user carries a second version in Weinberg–Srednicki
    two-component spinor notation (undotted/dotted indices, epsilon, sigma/sigmabar as in Srednicki ch. 34-36;
    Lorentzian base entries S-SPN-*, S-SYM-* of ~/.claude/skills/notation-unify/STANDARD.md, Weinberg via
    notation-db/weinberg-qtf-1995.md). The Euclidean adaptation is DERIVED from this project's own realisation of
    ESW (explicit, verified); the Euclidean/conformal blocks SIG/REP/CONF of STANDARD.md come from the user's blocked
    SUSY note and are NOT used. Dictionary: step2/notation/ (built per the notation-unify workflow P1-P4).
R2. No twist anywhere in the computation: the theory stays the full untwisted N=4 BV theory with all psl(4|4) ghosts.
    A twisted statement may only appear as a consequence (components of the untwisted equivariant cochain restricted
    to a nilpotent supercharge), never as an input.
R3. Reproduce the conclusions about Q_1 (the one-loop correction to the action of the twisting supercharge on local
    operators) of the holomorphic-twist literature, from the untwisted one-loop components along that supercharge
    (arity one, and arity two along it where S^(2) contributes), and display the renormalisation freedom explicitly:
    how each finite counterterm / scheme choice (arity zero and arity one) changes Q_1, and what is scheme independent.
    Name clash: in the manuscript Q_1 in Sigma and Q_2 in Sigma^* are spinor parameters (Prop prop:c7-susy); the
    user's Q_1, S_1, D_1, K_1 are one-loop (order hbar) components. Reports must never mix the two.
R4. Compute the one-loop components S_1 (special supercharges), D_1 (dilatation; Part D), K_1 (special conformal),
    and the remaining ones (P, M, R: covariance realised), on letters and on gauge-invariant local operators.

## User directives added 2026-10-28, second message (binding)
R0. Program goal, in the user's words: understand completely the one-loop BV structure of N=4 SYM, i.e. the
    realisation of the symmetries. Concretely: the one-loop psl(4|4)-equivariant BV structure in all arities
    (arity zero, one, two), its obstruction theory and scheme freedom, and its action on local operators (all
    components Q_1, S_1, D_1, K_1, P_1, M_1, R_1 and the arity-two parts). Every Step-2 deliverable is judged
    against this.
R5. The CFT inner product on gauge-invariant local operators up to one loop, as a test of Q^dagger = S.
    (a) Conjugation. The antilinear anti-automorphism X -> X^dagger of psl(4|4) of radial quantisation
        (Q^dagger = S, P^dagger = K, D^dagger = D, rotations and R-symmetry in their compact real forms),
        derived in the realisation of Chapter 2 and stated in two-component form; its realisation on fields
        (inversion composed with the antilinear real-slice conjugation); its relation to the Osterwalder-Schrader
        reflection in a hyperplane (a linear isometry, so a candidate exact symmetry of K_t, P(eps,infty), I and
        the rotation-covariant counterterms) through a conformal map generated by P and K. The one-loop
        realisation of the conjugation is then fixed by the reflection and by P_1, K_1, D_1.
    (b) Inner product. From the two-point functions of the observables O_1[infty], O_2[infty] (ch05) at one
        loop: tree part and one-loop part (code names G0, G1; never new symbols in text), including the
        order-g tree-level entries between operators of different length; separation dependence fixed by D_1;
        what is canonical and what depends on the scheme (arity-zero finite counterterms, finite parts of
        operators, arity-one/two finite counterterms, gauge fixing).
    (c) Test. G(X a, b) = G(a, X^dagger b) through one loop for X = Q (and every generator), with Q_1, S_1,
        K_1, D_1 from R3/R4 and G1 computed independently from two-point graphs; covariance of the test
        under every scheme change (finite parts of operators must be compatible with the conjugation).
    (d) Scheme-independent targets: D_1 symmetric for the tree inner product; order-g^2 behaviour of protected
        two-point functions with the normalisation fixed by symmetry; every literature statement, with its
        regulator recorded (statements from DR-based papers are statements only).
    Evaluation: position-space UV-boundary graphs, exact; no dimensional regularisation in any role; no IR scale.
R6. (User decision 2026-10-28) K_1 at large distance: first try to prove that the total log R coefficient of the
    K-route graphs (arity one along K_v, and arity two along (eta1,eta2) in Sigma* x Sigma*) cancels. If it does not
    cancel, define K_1 through K = P^dagger with the one-loop inner product of R5 (P_1 = 0), with no IR cutoff.
    An eps -> 0 before R -> infinity prescription is NOT adopted.
R7. (User correction 2026-10-28, supersedes the "no IR scale" wording in R5 and in R6) The standard Costello
    transfer endpoint T is allowed (DESIGN.md rule 4, amended): P(eps,T), effective action and observables at
    scale T. Core UV quantities (counterterms, obstructions, finite counterterms, local representatives L_1(O),
    components Q_1..K_1) must be T-independent; T must not enter the core theory; take T -> infinity wherever it
    is safe. Consequence for K_1 (replaces R6): define K_1 by the K-route graphs at finite T and prove that the
    local representative is T-independent (the large-distance log terms of the defect graphs carry O(eps)
    coefficients, so they should vanish as eps -> 0 at fixed T); then T -> infinity where safe. K = P^dagger with
    the one-loop inner product is a cross-check, and the definition only if T-independence fails.
R8. (User insight 2026-10-28) Finiteness of the one-loop data on local operators. The full Q_1, S_1, D_1, ... involve
    the infinite tower of covariant derivatives, even in the 1/16-BPS sector. But at one loop a graph meets a local
    operator in at most two letters (a third would need a second independent loop momentum): O on the cycle uses
    exactly two of its letters; O on a tree attached to the cycle uses one. Hence the one-loop counterterm of operators
    and every one-loop component L_1 act as (one-body derivation by finitely many action-level local functionals:
    field counting, finite counterterms, arity-one counterterms) + (a two-body, second-order operator on the jet algebra:
    O -> sum over pairs of letters of a universal kernel(pair) x spectators). The outputs are bounded as well: for a pair
    of total Weyl weight w the singular parts have weight w - 2m (eps^{-m}) and w (log eps), so a graph with enough extra
    external legs is UV finite and contributes neither a counterterm nor a defect; the number of output fields is
    bounded by the weight, and gauge covariance (the BRST/QME relations of the counterterms) organises the extra legs
    into covariant derivatives and field strengths, so the independent data are the minimal cycles. At two loops at most
    three letters enter; each loop order needs finitely many kernel types.
    Task: compute the one-loop counterterm d_delta I^CT[eps](O) for ALL local operators O (all letters with arbitrary
    derivatives, gauge-variant ones included), as universal two-letter kernels in closed form (generating functions in
    the momenta of the two letters; Schwinger/Feynman-parameter integrals; eps^{-m} and log eps parts), exact, no DR,
    with the output bound and the covariance reduction proved; then D_1 on all local operators (check: Beisert's
    complete one-loop dilatation operator, hep-th/0307015, harmonic action, scheme-independent spectrum, e.g. twist-two
    harmonic sums), and the two-body kernels of Q_1, S_1, K_1 from the defect graphs (cf. [BK, (3.7)]: the twisted Q_1 is
    a second-order operator on pairs of letters). Engines cluster-portable (pure Python/C).
R9. (User, 2026-10-29) Every report also carries physics notation as a supplement to each part (core stays abstract):
    S_+ = undotted / S_- = dotted spinors; W, W^* = 4, 4bar of SU(4)_R, W_R^perp = vector of SO(6)_R = antisymmetric 6
    of SU(4) (A_ij = A_I Sigma^I_ij); rho = gamma matrices (10d chiral Gamma^M; 4d sigma, sigmabar; Weinberg/Srednicki
    gamma); Pi = parity shift (Grassmann-odd fields); 10d spinor indices with upper/lower kept distinct; ghost numbers
    and statistics; su(N) normalisation. In the reader: "Physics notation" callouts inserted after the definitions and
    sections they translate; the existing blocks are never edited.
