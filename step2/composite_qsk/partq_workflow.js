export const meta = {
  name: 'n4-odd-and-k-two-body-kernels',
  description: 'Part Q rerun: one-loop two-body kernels of Q_1, S_1, K_1 on all local operators, twist check, two engines, referee, revision, Chapter 11 section 11.6',
  phases: [
    { title: 'Theory', detail: 'defect graphs, decomposition, census, closed forms, renormalisation freedom, twist check' },
    { title: 'Engines', detail: 'two independent engines' },
    { title: 'Refute', detail: 'adversarial referee and comparator' },
    { title: 'Revise', detail: 'apply findings, final results' },
    { title: 'Write', detail: 'replace Section 11.6 of Chapter 11 and its supplement' },
  ],
}

const PRE = `Project: /Users/libotao/Desktop/N4_Costello_BV_2026-10-27 — N=4 super Yang–Mills in Costello's BV formalism on Euclidean R^4, heat-kernel proper-time UV boundary epsilon on every edge (P(eps,T) = int_eps^T (Q^GF (x) 1) K_t dt), transfer endpoint T (allowed; it is the endpoint of Costello's transfer, NOT an IR cutoff; core UV data are T-independent; T = infinity where safe, finite T where Def c8-K requires it), constant psl(4|4) ghosts.
Binding spec: DESIGN.md, DESIGN_STEP2.md (R0–R9; R8 is the two-letter programme). Source of truth: manuscript/ch01..ch10 (*.tex) and manuscript/FORMAT.md. Key places: Ch7 (components, Prop c7-cov, Rem c7-internal), Ch8 (Def c8-components, Thm c8-L, Def c8-K, Thm c8-K, Prop c8-Tsafe, conjugations, Def c8-radial), Ch9 (composite A|_V–c line, Prop c9-composite, Thm c9-D1, Rem c9-Pi, Thm c9-P1, Prop c9-scheme), Ch10 (two-letter kernels: Thm c10-D1, Def c10-hyp, Prop c10-spec, Cor c10-P7, Prop c10-bound1, Rem c10-next). Previous workflow outputs to build on: step2/twoletter/ (theory/twoletter_theory.tex, engine_a, engine_b — reuse their verified code and conventions —, results/kernels_final.json, qa/). Earlier design for the odd components and K: step2/theory/odd_components_design.tex, engine_plan.json, k_transfer.tex/.json, step2/qa/k_judgement.json; twist targets step2/literature/q1_targets.(tex|json), bgkwwy, gkw; notation dictionary step2/notation/dictionary.json.
Notation: the manuscript's index-free notation (V, W_R^perp, Sigma = S_+ (x) W + S_- (x) W^*, rho, Gamma; letters A|_V, A|_perp, lambda, c, antifields), with the Weinberg–Srednicki two-component form of every reported formula. No new notation beyond what you define explicitly; name objects by their role. Unit u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2).
Hard rules: never use dimensional regularisation in any role; no IR regulator other than the transfer endpoint T; do not read or reference earlier N=4/QED/SQED/SUSY projects elsewhere on this machine; at most one heavy process of yours at a time, each < 1.5 GB RSS; never write to /tmp; pure Python (fractions, sympy, numpy, mpmath) or C, cluster-portable; no Wolfram, no cluster, no browser; do not edit manuscript/ or reader/ unless told; subagents may not write .md report files (use .txt or .json).
Environment: the disk is nearly full and iCloud evicts files in this Desktop folder minutes after they are written (reading an evicted file can hang): keep code and caches in subfolders named *.nosync, and if a read stalls, move on and report it. The network drops occasionally: write results and a progress log (PROGRESS.txt) incrementally so that an interruption loses little. A single response is capped at 128k output tokens: write long files in pieces.
Workflow root: step2/composite_qsk/. Part C (composite kernels) is finished: results in step2/composite_qsk/results/composite_final.json, theory composite/theory/composite_theory.tex, engines composite/engine_a, composite/engine_b, and manuscript/ch11_composite_and_odd_kernels.tex Sections 11.1-11.5. If a repository file is evicted by iCloud (ls -lO shows dataless) and reading hangs, read the identical copy from GitHub with gh api repos/lbotao399-stack/full-bv-structure-n4-sym/contents/<path> -H 'Accept: application/vnd.github.raw' when its size equals the local size; never overwrite local files.`


const SUM = {
  type: 'object',
  properties: {
    summary: { type: 'string' },
    files: { type: 'array', items: { type: 'string' } },
    key_results: { type: 'array', items: { type: 'string' } },
    open_issues: { type: 'array', items: { type: 'string' } },
  },
  required: ['summary', 'files', 'key_results', 'open_issues'],
}
const AUDIT = {
  type: 'object',
  properties: {
    verdict: { type: 'string', enum: ['pass', 'pass_with_fixes', 'fail'] },
    summary: { type: 'string' },
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          severity: { type: 'string', enum: ['fatal', 'major', 'minor'] },
          claim: { type: 'string' }, evidence: { type: 'string' }, fix: { type: 'string' },
        },
        required: ['severity', 'claim', 'evidence', 'fix'],
      },
    },
  },
  required: ['verdict', 'summary', 'findings'],
}


const RESUME = '\nRESUME NOTE: an earlier attempt of this exact task was cut off by a network failure. Reuse what is already in your output folder (read PROGRESS.txt first) instead of starting over.'
async function run(prompt, opts) {
  let r = await agent(prompt, opts)
  for (let k = 0; !r && k < 2; k++) { log(`retry ${opts.label}`); r = await agent(prompt + RESUME, { ...opts, label: opts.label + '-retry' + (k + 1) }) }
  return r
}

phase('Theory')
const thQ = await run(`${PRE}

ROLE: theorist, PART Q (two-body kernels of Q_1, S_1, K_1). Output folder: step2/composite_qsk/qsk/theory/.
Goal: the one-loop components of the superconformal action on ALL local operators along the supercharges Q (S_+ (x) W and S_- (x) W^*), the conformal supercharges S, and the special conformal K_b, as universal one-body + two-body kernels on letters (R8: at one loop at most two letters of O enter), with the renormalisation freedom explicit, and the comparison with the holomorphic twist (R3: the twisted Q_1 of [BK, (3.7)] is a second-order operator on pairs of letters; targets in step2/literature/q1_targets.(tex|json)).
Deliver in theorem form, step2/composite_qsk/qsk/theory/qsk_theory.tex (pdflatex into _build; no errors, no undefined refs), written section by section:
1. The defect graphs: for x = Q, S, K_b, which one-loop graphs with one vertex S^(1)(x) and one vertex O give the component L_1(O)(x) (Def c8-components, Thm c8-L): S^(1)(x) on the cycle (does not commute with the regulator: Q^GF and H are not covariant under supersymmetry or K; derive the commutators and the resulting special/composite edges), S^(1)(x) on a leg, and the arity-one counterterms; reuse step2/theory/odd_components_design.tex and k_transfer.tex after checking them against Ch8 (they predate Chapters 8–10).
2. Decomposition into one-body + two-body kernels on letters, with Koszul signs and colour structure as in Ch10; census of pair types and topologies; output bound (weights and ghost number; for S and K the explicit x-dependence of the Killing data changes the weight count — state it).
3. Closed-form method for the singular and finite parts (the components are finite limits of graphs with special edges: use the special-edge lemma and simplex formulas of Ch10), at T = infinity for Q and S, at finite T for K_b with the explicit log T term of Def c8-K (Prop c8-Tsafe: coefficient (3/8 pi^2) hbar times the cut tree), and the K-route vs K = P^dagger relation as a cross-check.
4. Renormalisation freedom: which finite parts of the operators and of the arity-one counterterms change Q_1, S_1, K_1 and how (Prop c8-laws, Prop c9-scheme), and the invariant content (the class up to homotopy, Cor c8-homotopy).
5. The twist check: restrict Q_1 along a twisting supercharge Q^4_- (u (x) z with z = e_4, Ch8 twist remark) and compare with the twisted Q_1 of [BK, (3.7)] / the q1_targets entries (with the renormalisation freedom of the twist made explicit); and the expected relation S_1 = -(Q_1)^adj + Gram terms (Cor c8-gram) as a later inner-product test.
6. Engine spec qsk/theory/engine_spec.json for two methodologically independent routes extending step2/twoletter/engine_a and engine_b (and the composite engines of part C if needed: those run before the QSK engines; their outputs will be in step2/composite_qsk/composite/engine_a, engine_b), with gates: Q^2-type closure relations ([Q_1, Q_0] + [Q_0, Q_1] terms from Thm c8-L / Cor c8-homotopy), the twist comparison, {Q, S} relations at one loop through D_1, and K_1 T-independence.
Return the summary object.`, { label: 'theory-qsk', phase: 'Theory', schema: SUM })
if (!thQ) { log('theory failed; stopping'); return { thQ } }

const ENG = (part, name, route, other, extra) => `${PRE}

ROLE: engine ${name} for PART ${part}. Output folder: step2/composite_qsk/${part === 'C' ? 'composite' : 'qsk'}/engine_${name.toLowerCase()}/ (code and caches in src.nosync/, out.nosync/). Do NOT read step2/composite_qsk/${part === 'C' ? 'composite' : 'qsk'}/engine_${other}/ (independence).
Read the theory of this part (step2/composite_qsk/${part === 'C' ? 'composite/theory/composite_theory.tex and engine_spec.json' : 'qsk/theory/qsk_theory.tex and engine_spec.json'}) and implement ${route}. ${extra}
Compute every item of the engine spec with exact coefficients (rationals, pi powers, logs as they arise), explicit tables to the spec's weight cutoffs and closed forms where the method gives them (verified against the tables). Record every gate (value, target, pass/fail) in gates.json; results in results.json (machine-readable, exact strings); method and conventions in notes.txt. One script run_all.py runs everything and prints resources. Return the summary object.`


phase('Engines')
const engQ = await parallel([
  () => run(ENG('Q', 'A', 'route A (simplex/Schwinger generating functions, extending step2/twoletter/engine_a and step2/composite_qsk/composite/engine_a)', 'b', 'The composite engines of part C are finished (step2/composite_qsk/composite/engine_a, engine_b); use them where the theory needs composite edges.'), { label: 'engine-Q-A', phase: 'Engines', schema: SUM }),
  () => run(ENG('Q', 'B', 'route B (explicit momentum-space numerator algebra with Feynman parameters plus numerical eps-extrapolation, extending step2/twoletter/engine_b and step2/composite_qsk/composite/engine_b)', 'a', 'The composite engines of part C are finished (step2/composite_qsk/composite/engine_a, engine_b); use them where the theory needs composite edges.'), { label: 'engine-Q-B', phase: 'Engines', schema: SUM }),
])
if (!engQ[0] && !engQ[1]) { log('both engines failed; stopping'); return { thQ, engQ } }

phase('Refute')
const refQ = await run(`${PRE}

ROLE: adversarial referee and comparator of PART Q (two-body kernels of Q_1, S_1, K_1). Output folder: step2/composite_qsk/qsk/qa/. Default to finding problems.
Engine summaries: A ${JSON.stringify(engQ[0])}; B ${JSON.stringify(engQ[1])}.
1. Compare engines A and B entry by entry; resolve disagreements by your own independent computation.
2. Recompute independently a sample: each generator type (Q along S_+ (x) W and S_- (x) W^*, S, K_b), each topology, fermion and ghost/antifield pairs, derivative letters; the finite-T log T term of K_1 and its T-independence after epsilon -> 0.
3. Check the closure gates yourself (the one-loop relations from Thm c8-L / Cor c8-homotopy, {Q,S} through D_1), the renormalisation-freedom statements, and the twist comparison against [BK, (3.7)] / q1_targets (the twist must not be an input: it is a check).
4. Attack the theory (qsk/theory/qsk_theory.tex) against Ch7–Ch10.
Write compare.json and findings.json. Return the audit object.`, { label: 'refute-qsk', phase: 'Refute', schema: AUDIT })

phase('Revise')
const revised = await run(`${PRE}

ROLE: reviser of PART Q. You may edit step2/composite_qsk/qsk/ and step2/composite_qsk/results/qsk_final.json; not manuscript/, not reader/.
Inputs: theory ${JSON.stringify(thQ)}; engines ${JSON.stringify(engQ)}; referee ${JSON.stringify(refQ)} (full findings in step2/composite_qsk/qsk/qa/findings.json).
Check every finding yourself; accept or reject each with a reason (step2/composite_qsk/qsk/disposition.json). Fix the theory (recompile, 0 errors, 0 undefined refs) and engines for accepted findings; rerun what changed and the affected comparisons. Produce step2/composite_qsk/results/qsk_final.json (exact coefficients, closed forms, renormalisation freedom explicit, gates with values, twist comparison) and step2/composite_qsk/results/qsk_summary.txt (formula-centred, both notations, proved vs computed vs open). Return the summary object.`, { label: 'revise-qsk', phase: 'Revise', schema: SUM })
if (!revised) { log('revision failed; stopping'); return { thQ, engQ, refQ } }

phase('Write')
const chapter = await run(`${PRE}

ROLE: author of the Q_1, S_1, K_1 part of Chapter 11. You may edit manuscript/ch11_composite_and_odd_kernels.tex ONLY in its last section (Section 11.6, the components along Q, S and K, currently structure only) and append to manuscript/supplements/ch11.tex (do not change its existing blocks, which may already have been refereed). Sections 11.1-11.5 are final: do not change them.
Inputs: step2/composite_qsk/qsk/theory/qsk_theory.tex (final), step2/composite_qsk/results/qsk_final.json and qsk_summary.txt, revision summary ${JSON.stringify(revised)}.
Replace Section 11.6 by the computed results in the style of Sections 11.1-11.5 and FORMAT.md: numbered statements labelled c11- (keep every existing c11- label of Section 11.6 that other text cites, grep for them; new labels must not clash), proofs, tables, closed forms, the renormalisation freedom, the twist comparison stated as what it is, open items stated affirmatively, a Remark with the Weinberg-Srednicki two-component forms; no audit narratives; cross-references only to existing labels. Append physics-notation remarks for the new statements to supplements/ch11.tex in the same anchor format (supercharges Q^i_alpha, Qbar_alphadot i, S, K_m; the twisted Q_1 in the usual language).
Verify: python3 reader_build/convert.py from the project root ("unresolved": []), KaTeX on every math string of Chapter 11 (pattern step2/qa/ch10_review/katex_check.py adapted), python3 reader_build/prep_supplements.py (0 KaTeX errors, no problems). Do not touch reader/. Return the summary object.`, { label: 'write-ch11-sec6', phase: 'Write', schema: SUM })

return { thQ, engQ, refQ, revised, chapter }
