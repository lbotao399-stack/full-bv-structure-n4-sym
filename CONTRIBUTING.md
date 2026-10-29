# Contributing

This file is the working brief for collaborators — people or AI agents (e.g. a Claude Code session) — who continue the
project. Contributions come by **fork and pull request**; no direct write access to this repository is needed.
An agent can be started by giving it this file as its instructions.

You are joining an ongoing research project as an independent collaborator. Repository (public):
https://github.com/lbotao399-stack/full-bv-structure-n4-sym — "Full BV structure of N=4 SYM", work in progress.

## 0. Setup

1. Fork the repository to your account (`gh repo fork lbotao399-stack/full-bv-structure-n4-sym --clone`). Work only in your
   fork, on a branch `collab/<your-github-handle>/<track>` (tracks in §4). Do not push to the upstream repository.
2. Read, in this order, before doing anything else:
   - `README.md`, `DESIGN.md` and `DESIGN_STEP2.md`, the design rules and research directives R0–R9;
   - `HANDOFF_STEP2_2026-10-28.md`, the running log of status and open items. Paths in it are absolute paths on the
     owner's machine; read them relative to the repository root;
   - `manuscript/FORMAT.md` and the chapters `manuscript/ch00_overview.tex` … `ch10_two_letter_kernels.tex`.
     Chapter 8 contains the binding definitions;
   - the material of your track, listed in §4.
3. Python 3 with `fractions`, `sympy`, `numpy` and `mpmath` is enough. Node.js is needed only for the KaTeX check in §5.
   Do not use Mathematica.

## 1. What the project is

The project uses Costello's BV formalism to understand how the symmetry algebra acts on local operators in N=4 SYM, and
in ordinary QFTs. The setting has three ingredients:
- a UV boundary ε: a heat-kernel proper-time cutoff on every edge, with propagator P(ε,T) = ∫_ε^T (Q^GF⊗1)K_t dt;
- the transfer to scale T;
- the quantum master equation, with the superconformal algebra psl(4|4) built in through constant ghosts.

The targets are:
- the one-loop components Q₁, S₁, K₁, D₁, P₁, M₁, R₁ of the symmetry action on local operators, with the
  obstruction theory and the renormalisation freedom explicit;
- the CFT (radial) inner product to one loop, and the test Q† = S.

The holomorphic twist is never an input. It is only a check.

## 2. Non-negotiable rules

1. **No dimensional regularisation in any role**: not as a scheme, not as an intermediate tool, not as a cross-check.
   Literature results obtained with DR may be compared only at the level of scheme-independent quantities, such as
   spectra.
2. **T is the endpoint of Costello's transfer, not an IR cutoff.** There is no other IR regulator. Core UV quantities
   (counterterms, obstructions, local representatives, K₁) must be T-independent. Use T = ∞ where safe and finite T
   where Def c8-K requires it.
3. **Notation.** Use the manuscript's index-free notation:
   - V, W_R^⊥, Σ = S₊⊗W ⊕ S₋⊗W*, ρ, Γ;
   - the letters A|_V, A|_⊥, λ, c and the antifields A*, λ*, c*;
   - the unit u = ħ tr_g(ad ad)/(16π²(X,Y)) = −g²N/(32π²).

   Every reported formula also gets its Weinberg–Srednicki two-component form (Def c8-twocomp,
   `step2/notation/dictionary.json`). Introduce no notation without an explicit definition, and never rename existing
   objects.
4. **Evidence standard.** A numerical claim enters a chapter only after all of the following:
   - two methodologically independent engines agree entry by entry, for example Schwinger/simplex generating functions
     versus momentum-space Feynman-parameter algebra with numerical ε-extrapolation;
   - an adversarial referee has recomputed a sample by a third method;
   - every gate (known values: Konishi γ = 3g²N/(4π²), 0 for the 20′, the Chapter 9/10 tables, scheme-independent
     spectra) is recorded with its value.

   Exact arithmetic wherever possible. Code must be pure Python or C and runnable end to end from one `run_all.py`.
5. **Style** (FORMAT.md). Theorem form throughout: numbered Definition/Lemma/Proposition/Theorem/Corollary/Remark,
   with framed proofs. Formula chains carry the argument, with short connecting prose. State results affirmatively.
   Each chapter opens with 2–4 sentences on the problem and why it is not trivial. No audit notes, script names or file
   paths inside chapters. Only the markup listed in FORMAT.md is allowed (the reader converter and KaTeX depend on it).
6. **Third-party material.** Cite papers as [Author/arXiv id, eq.]. Never commit papers, their PDFs, full texts or
   sources.
7. **Do not edit existing chapters (ch00–ch19) or the supplements of other chapters.** Propose corrections as a
   corrections file (§5).

## 3. Current state (late October 2026)

- **Chapters 0–7**: classical BV data, superconformal action, UV boundary, one-loop counterterms (finite counterterm
  27/16), local operators, the dilatation component, and the components along the other generators.
- **Chapter 8**: definitions of the one-loop equivariant structure (obstruction, local representative L(M) and its
  components, invariants, scheme laws, conjugations, inner product).
- **Chapter 9**: the cubic obstruction
  Ob₃ = u[(11/4 + 19/4 log 3/2)(G₁+G₃) + 27/8 G₂ − 9 log(3/2) G₅], the composite A_V–c line in Feynman gauge, and
  D₁ = 2cO − (1/3+θ)u{S⁽¹⁾(Δ),O}.
- **Chapter 10**: the counterterms of all local operators as universal two-letter kernels, with closed forms. They are
  checked against the su(2) Heisenberg spectra and the twist-two harmonic sums −8h(S)u.
- **Chapter 11, §11.1–11.5**: the composite kernels, and D₁, P₁ on all operators (results in
  `step2/composite_qsk/results/composite_final.json`).
- **In progress by the owner's Claude** (do not duplicate): the two-body kernels of Q₁, S₁, K₁ and the comparison with
  the holomorphic twist, which will replace §11.6 (`step2/composite_qsk/qsk/`).

## 4. Tracks

Claim one track before starting: open an issue titled `Claim: <track letter> — <your handle>` in the upstream
repository. Tracks are listed in priority order.

**Track A — the one-loop inner product and the test Q† = S (directive R5).**
- **Read:** Chapter 8 (Def c8-2pt, Def c8-radial, Def c8-inner, Thm c8-adjoint, Cor c8-gram);
  `step2/theory/inner_product_design.tex`, `two_point_design.tex`, `conjugation.tex`;
  `step2/theory/engine_plan_ip.json` (stages IP0–IP8); and the gap list `step2/qa/ip_completeness2.json` (22 items).
- **Do:**
  - Build the engines for IP0–IP4 and IP6: tree and one-loop two-point coefficients, radial Grams, descendants and the
    Shapovalov form, and frame changes.
  - Re-derive independently the minimal-scheme Gram constants for the 20′, unit·(27/2 + 4 log 3/2), and for Konishi,
    unit·(27/2 + 4 log 3/2 + 24(γ_E − 3 log 2)). Each currently rests on one derivation.
  - Prove or refute Hermiticity of the one-loop radial form.
- The assembly IP5 (tests T0–T4) needs Q₁ and S₁ from Chapter 11. Prepare it so that it runs as soon as those are
  pushed.

**Track B — the arity-zero and arity-two master equations.**
- **The minimal scheme:** the one-loop radial constants use the "minimal" two-leg finite counterterm (the β = 0
  representative). The cubic obstruction ties the quadratic coefficients β_k to the cubic h_j: the classes
  (β₂,h₄) and (β₃,h₅) are fixed only by the equations E3/E3b. Determine whether the minimal representative extends
  equivariantly, and classify the degree-0, weight-4 s-closed local cochains (hypothesis H4, used in Thm c8-invariants
  (iii)).
- **The arity-two obstruction:** compute the arity-two obstruction Ob⁽²⁾(x₁,x₂) and the finite counterterm
  ħF₁⁽²⁾ along (Δ, P_b) and (K_b, P_b′). Decide whether they fix θ and the other free directions of ħF₁⁽¹⁾(Δ) (Def
  c10-hyp).
- **Read:** Chapter 4, Chapter 9 (Thm c9-F3, Thm c9-wz, Thm c9-quadratic) and Chapter 10 §10.8.

**Track C — rotations, R-symmetry and higher arity.**
- **so(V):** the composite line along so(V) is nonzero (Rem c7-internal) and has never been evaluated. Compute M₁ and
  M̄₁ on all operators with the two-letter machinery of Chapter 10.
- **R-symmetry:** R₁ is not computed either.
- **Higher arity:** prove or refute the vanishing of the one-loop components of arity ≥ 3 along position-dependent
  generators (Σ*, K_b).

**Track D — preparation for two loops.**
- At two loops, at most three letters of an operator enter a graph. Develop the structure theory: the decomposition
  into three-, two- and one-letter kernels, the census and the output bound, the BRST relations, and the closed-form
  method for two-loop simplex integrals with the per-edge cutoff.
- Then build a census engine that is cluster-portable (C or pure Python, no network, no numpy dependence for the
  heavy part).

## 5. Delivery

Deliver by pull request from your fork to `lbotao399-stack/full-bv-structure-n4-sym:main`. Open it as a draft PR early
and push incrementally. Mark it ready for review only when the checklist below holds.

**Where things go**
- All work goes in `collab/<handle>/<track>/`:
  - `theory/` (a theorem-form .tex that compiles with pdflatex, with no errors and no undefined references);
  - `engine_a/`, `engine_b/` (code + `run_all.py` + `results.json` + `gates.json`);
  - `qa/` (comparison table, referee findings, disposition of each finding);
  - `PROGRESS.txt`.
- The draft chapter goes in `manuscript/ch2N_<slug>.tex`. Chapter numbers 20–29 are reserved for collaborators: take
  the next free one, and prefix all labels with `c2N-`. Its physics-notation supplement goes in
  `manuscript/supplements/ch2N.tex`. Supplement format:
  - lines `% ANCHOR: end_of_section <n>` or `% ANCHOR: after_label <label>`, each followed by
    `\begin{remark}[Physics notation: <topic>] … \end{remark}`;
  - no `\label`, no equation environments, no `\item[`;
  - copy the pattern of `manuscript/supplements/ch10.tex`.
- Corrections to existing chapters go in `collab/<handle>/corrections_chNN.json`, one file per chapter. Each entry is
  {label, old_text, replacement_text, problem, evidence}. `old_text` must occur exactly once in the chapter (check it
  with a count), and `evidence` points to a script or result in your folder.
- `collab/<handle>/<track>/DELIVERY.json` is the manifest. It lists:
  - every claim with its status: proved / computed (two engines + referee) / conjectured / open;
  - every gate with value and target;
  - the engine agreement table (entries compared, max difference);
  - the referee verdict and findings with their disposition;
  - the command that reproduces everything, with runtime and peak memory.
- `collab/<handle>/<track>/report.pdf` is a formula-only PDF: formulas, logic symbols, tables, flowcharts and graph
  pictures, with no prose. All explanation goes in the PR description.

**Checklist before "ready for review"**
1. `python3 reader_build/convert.py` (from the repository root) reports `"unresolved": []` with your chapter present.
2. Every math string of your chapter and its supplement renders with KaTeX: use the vendored
   `reader/notion-reader/app/vendor/katex/katex.min.js`, displayMode for equations, `throwOnError: true`,
   `strict: 'ignore'`. Adapt the pattern in `step2/qa/ch10_review/katex_check.py`; it contains absolute paths from the
   owner's machine, so replace them.
3. Both engines run from `run_all.py`, agree, and all gates pass. The referee is independent: a separate session or
   agent with no access to the engines' code.
4. No dimensional regularisation, no IR regulator other than T, no new undefined notation, no third-party papers
   committed.
5. The PR description, in plain prose, gives: what was established, what is conjectured, what is open, and anything
   that contradicts an existing chapter (with the corrections file).

**Rhythm**
- One PR per milestone, for example "Track A: IP0–IP2". Keep PRs reviewable: a chapter draft of tens of pages is fine,
  but do not bundle unrelated tracks.
- Before opening a PR, rebase on upstream `main`, which may have gained Chapter 11 or other chapters, and re-run the
  checklist.
- If you find an error in an existing chapter, open an issue at once with the corrections file. Do not wait for the
  PR.
