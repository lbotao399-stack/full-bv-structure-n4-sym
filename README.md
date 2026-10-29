# Full BV structure of N=4 SYM

> **Status: work in progress, still exploratory and far from complete.** The definitions, statements and numbers here
> are being worked out. Some will change, and several central questions are open. Treat everything as a research notebook,
> not a finished result.

## Aim

The project uses Costello's BV formalism (renormalisation with a UV boundary ε, the transfer to scale T, and the
quantum master equation) to understand how the symmetry algebra acts on **local operators** in N=4 super Yang–Mills,
and in ordinary quantum field theories more generally. The focus is on supersymmetry and conformal symmetry:

- the one-loop corrections to the action of the superconformal generators psl(4|4) on local operators: the supercharges
  and conformal supercharges **Q_loop, S_loop**, the special conformal transformations **K_loop**, the dilatation
  **D_loop**, and the translations, rotations and R-symmetry, computed directly, **without twisting**;
- the obstruction theory and the renormalisation freedom that make this action well defined, up to homotopy;
- the **CFT inner product** (radial quantisation) to one loop, and the test **Q† = S** that it provides.

The holomorphic twist is not an input. Its one-loop Q₁ should come out as a special case, with the renormalisation
freedom shown explicitly.

## Approach

1. **Classical BV data** (Ch. 1–2): N=4 SYM in Euclidean signature in index-free 10d notation
   (Σ = S₊⊗W ⊕ S₋⊗W*, ρ, Γ), with the superconformal algebra built in as an equivariant (L∞) action through constant
   ghosts of psl(4|4).
2. **UV boundary** (Ch. 3–4): heat-kernel propagator P(ε,T) = ∫_ε^T (Q^GF⊗1)K_t dt. The one-loop counterterms come
   from the ε-singular parts, and the finite counterterm is fixed by the quantum master equation. The transfer
   endpoint T is the Wilsonian scale of Costello's transfer, not an IR cutoff: the core UV data are T-independent.
3. **Local operators** (Ch. 5–7): renormalised operators O(ε), their flow to scale T, and the components of the
   one-loop symmetry action. Examples: the dilatation component reproduces the Konishi anomalous dimension and the
   vanishing for the 20′.
4. **Definitions of the one-loop equivariant structure** (Ch. 8): the obstruction Ob and its cancellation sF₁ = Ob,
   the local representative L(M) of d_T M[T], its components (Q₁, S₁, D₁, K₁, P₁, M₁, R₁), invariants, scheme laws,
   the conjugations (radial †, Osterwalder–Schrader reflection), and the inner product.
5. **The cubic obstruction and the dilatation cochain** (Ch. 9): the three-leg obstruction,
   Ob₃ = unit·[(11/4 + 19/4 log 3/2)(G₁+G₃) + 27/8 G₂ − 9 log(3/2) G₅], its Wess–Zumino consistency,
   the composite A–c line, which exists because (Q^GF)² ≠ 0 in Feynman gauge, and
   D₁ = 2cO − (1/3+θ)·unit·{S⁽¹⁾(Δ),O} with a scheme parameter θ.
6. **In progress**: one-loop counterterms of *all* local operators. At one loop a graph meets a local operator in at
   most two letters, so the full infinite tower of covariant derivatives reduces to finitely many universal two-letter
   kernels with bounded outputs. D₁ follows on all operators, and the result is compared with Beisert's complete one-loop
   dilatation operator. Next: the two-body kernels of Q₁, S₁, K₁, and the inner-product test Q† = S.

Each chapter carries "Physics notation" supplements that translate the index-free notation into the
Weinberg–Srednicki two-component notation (dotted and undotted spinors, SU(4)/SO(6) R-symmetry indices, γ-matrices).

## Layout

| Path | Content |
|---|---|
| `manuscript/` | LaTeX chapters `ch00`–`ch09` (theorem form), `supplements/` (physics-notation remarks), `FORMAT.md` |
| `DESIGN.md`, `DESIGN_STEP2.md` | binding design rules and research directives |
| `HANDOFF_STEP2_*.md`, `STORY.md` | running log of what has been done and what is open |
| `step2/` | Step-2 work: theory notes (`theory/`), independent computation engines (`engine_a/`, `engine_b/`), audits (`qa/`), two-letter kernels (`twoletter/`) |
| `classical/`, `superconformal/`, `engine_a/`, `engine_b/`, `compare/`, `qa/` | Step-1 checks: classical algebra, superconformal closure, one-loop counterterm engines and audits |
| `diagrams/` | graph census and diagrams |
| `reader_build/` | scripts that turn the manuscript into an interactive reader |

## Conventions and rules of the project

- Regularisation is the heat-kernel proper-time cutoff ε on every edge. **Dimensional regularisation is not used in any
  role**, and there is no IR regulator besides the transfer endpoint T.
- Each numerical claim is computed by at least two independent engines and checked by adversarial refuters before it
  enters the manuscript.
- Conventions are fixed in the manuscript itself. New notation is introduced only with a definition.

## Contributing

Questions, corrections and contributions are welcome: open an issue, or a pull request against `manuscript/`
or the computation folders. Please keep to the conventions above.
