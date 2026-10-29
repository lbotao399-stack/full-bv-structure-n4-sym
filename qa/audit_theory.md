# Theory audit (adversarial, lens: theory)

Inputs read: `DESIGN.md`, `sources/ledger_SYNTHESIS.md`, `classical/classical_spec.md` + `checks.json`, `superconformal/spec.md` + `verify.json` (+ `verify.py` generator code), `engine_a/results.json` (+ `faces.py` header), `engine_b/results.json` (+ `route2.py` header), `compare/compare.json`, `diagrams/manifest.json` (route labels only).

Independent checks (own code, exact sympy, no engine code imported):
- `qa/audit_theory/superpoly.py`: super-polynomial algebra with left derivatives.
- `qa/audit_theory/toy_gf.py` → `toy_gf.json`: a one-mode finite-dimensional copy of the minimal complex. The fields are $c$; $A$ (two $V$ components and one $W_{\mathbb R}^\perp$ component); $A^*$; $c^*$; $\lambda$ (2); $\lambda^*$ (2). The operators are replaced as follows: $d$ by a column with $|d|^2=h$; $d^*d$ on $A$ by $h-dd^{T}$ on $V$ and by $h$ on $W_{\mathbb R}^\perp$; $\not d$ by an antisymmetric matrix squaring to $-h$. All conventions are those of DESIGN, the ledger and classical_spec §3.
- `qa/audit_theory/toy_acyclic.py` → `toy_acyclic.json`: Q-cohomology of the polynomial functionals in the same toy.

Verdict scale:
- **confirmed**: the claim holds as stated, including by independent check.
- **conditional**: the algebra is right, but it rests on an unproved input, which is named.
- **refuted**: the argument has a gap or the statement is false; the exact gap is given.
- **open**: nothing was computed.

---

## 0. Verdict table

| # | Claim | Verdict | Exact gap (details in §1–§6) |
|---|---|---|---|
| 1 | $H=[Q,Q^{GF}]$ = positive Laplacian ⊗ id on every component | confirmed | — |
| 2 | $[Q^{GF},H]=0$ | confirmed | trivial, since $H$ is scalar |
| 3 | $P(\varepsilon,\infty)$ graded symmetric | confirmed | rests on $(Q^{GF}\otimes1-1\otimes Q^{GF})K_0=0$ (verified) |
| 4 | $QP(\varepsilon,\infty)=K_\varepsilon$ | conditional | exact algebraically; on $\mathbb R^4$ it uses $K_\infty:=0$ (G1, no source) |
| 5 | $(Q^{GF})^2\neq0$ | confirmed | $A^*\mapsto2d^*(A^*\vert_V)[c]$, $c^*\mapsto-2dc^*[A\vert_V]$ |
| 6 | DESIGN: "$(Q^{GF})^2\ne0$ is **the one** axiom of [C07 Def 4.3.3] that fails" | **refuted** | 4.3.3(5) also fails. $\operatorname{Im}Q^{GF}$ is not isotropic. [C07 Lemma 6.5.1], [C07 Lemma 6.6.2] and CS1 (negative definite on $(\operatorname{Im}Q^{GF})_{\mathbb R}$) fail. [CG2 §7.2.1] also requires square-zero |
| 7 | UV-boundary identity $(d_{CE}+Q)e^{\hbar\partial_P}=e^{\hbar\partial_P}(d_{CE}+Q+\hbar\Delta_\varepsilon)$ | confirmed (identity); **refuted** (stated corollary) | The identity holds with $+\hbar\Delta_\varepsilon$, $\Delta_\varepsilon=-\partial_{K_\varepsilon}$, and without $(Q^{GF})^2=0$. The corollary "every singular coefficient therefore sits at $t=\varepsilon$" does not follow from it |
| 8 | $O(\varepsilon):=O-\partial_\delta I^{CT}[\varepsilon]$ | **refuted** (as sourced) | The sign is consistent with CS4. [C07 Thm B] needs "at least cubic mod ħ", which fails for $\delta O$ with $O$ bilinear or a single letter. Quantum closedness of $O$ is not checked |
| 9 | $O[\infty]:=\lim\partial_\delta W(P(\varepsilon,\infty),I-I^{CT}[\varepsilon]+\delta O(\varepsilon))$ | **refuted** | Built from $I$, not $I_{\mathfrak g}$. Then $O[\infty]$ evaluated on any $X$ is $0$, so $[(d_{CE}+Q)O[\infty]](X)\equiv0$. Equivariant counterterms are absent. [CG2 Def 8.3.1.3, 8.4.0.1] (family over scales, support) are not satisfied |
| 10 | "At scale ∞, $\Delta_\infty=0$, $\{,\}_\infty=0$, so the differential is $(d_{CE}+Q)$" | **refuted** (as a contentful statement) | Formally true given $K_\infty=0$. But then $P(0,\infty)$ is a contracting homotopy and the scale-∞ complex is acyclic mod constants (verified in degrees 1–3). [CG2 Prop 13.4.0.1] does not cover scale ∞ or global $\mathfrak g$ quantum |
| 11 | Dilatation route (i) (scaling) | conditional → **refuted** as a derivation | The scaling algebra is correct. Assumption (A) is checked only at quadratic order (independently re-verified) and at the ghost vertex. Assumption (B) is asserted, not derived. Identifying it with $[(d_{CE}+Q)O[\infty]](\Delta)$ is not derived |
| 12 | Dilatation route (ii) and "two routes agree" | **refuted** | The engines' route (ii) is $\varepsilon\,d/d\varepsilon$ of the same graphs ($L_\varepsilon$-pinned), so it is not independent. DESIGN's route (ii) ($K_\varepsilon$ edges with an $S^{(1)}_{\mathfrak g}(\Delta)$ vertex) was evaluated by nobody |
| 13 | Sign / normalization of the dilatation | **refuted** (compliance, prose) | The engines use [SW (3)] $\Delta=-E$; DESIGN mandates [SW (17)] ($=2E$). Engine A's prose has the wrong sign for the special-edge sum. The absolute sign of one-generator components is a convention, not fixed in DESIGN |
| 14 | Selection rules $P,M,R$: "exactly $(X.O)[\infty]$, no one-loop correction" | **refuted** (as "exactly") | Argued from regulator covariance only. The one-loop Ward identity $Q([O[\infty]](X))=(X.O)[\infty]$ is not derived |
| 15 | Selection rules $Q,\bar Q,S,\bar S,K$ | open | Lists of allowed structures only. The "$K$ fixed by $[K,P]$" step is not a consequence of the trivial-module CE complex. IR finiteness of $x^2$-weighted internal insertions is not addressed |
| 16 | Superconformal tower: algebra, $S^{(1)}$, $S^{(2)}$, $S^{(3)}=0$, arities 1–4 | confirmed with caveats | Random-dense evidence, not a proof ("probability 1" is false). The psl(4\|4) identification uses necessary conditions only. Only $\mathfrak{su}(2)$ was run. $S^{(2)}$ uniqueness is open |
| 17 | Tower downstream claims: "one-generator component involves $S^{(1)}(X)$ only; its sign is independent of the $d_{CE}$ convention" | **refuted** | At one loop the ħ-corrections and arity-1 counterterms enter. The sign of $S^{(1)}$ is tied to the $d_{CE}$ convention by the arity-2 equation |

---

## 1. Feynman-type $Q^{GF}$ (rows 1–6)

**1a. Independent reproduction** (`toy_gf.json`). With $Q:=\{S_{free},-\}$ and $S_{free}=-\tfrac14(dA,dA)+\tfrac12(\lambda,\not d\lambda)+(dc,A^*)$:
$$
Qc=-dc\,[A\vert_V],\quad QA=-\tfrac12 d^*dA\,[A^*],\quad QA^*=d^*(A^*\vert_V)\,[c^*],\quad Q\lambda=\not d\lambda\,[\lambda^*]\qquad(\text{all four signs of classical\_spec §5}).
$$
Take the DESIGN ansatz $Q^{GF}A^*=\alpha A^*[A]$, $Q^{GF}A=\beta d^*(A\vert_V)[c]$, $Q^{GF}c^*=\gamma dc^*[A^*\vert_V]$, $Q^{GF}\lambda^*=\delta\not d\lambda^*[\lambda]$. Then
$$
H=[Q,Q^{GF}]=h\cdot\mathrm{id}\iff(\alpha,\beta,\gamma,\delta)=(-2,-1,1,-1)\ \ (\text{unique}),\qquad [Q^{GF},H]=0,
$$
$$
(Q^{GF})^2A^*=2d^*(A^*\vert_V)[c],\qquad (Q^{GF})^2c^*=-2dc^*[A\vert_V],
$$
$$
(Q\otimes1+1\otimes Q)K_0=0,\qquad (Q^{GF}\otimes1-1\otimes Q^{GF})K_0=0,\qquad (Q^{GF}\otimes1+1\otimes Q^{GF})K_0\neq0,
$$
$$
P(\varepsilon,\infty)\ \text{even and graded symmetric},\qquad (Q\otimes1+1\otimes Q)P(\varepsilon,\infty)=K_\varepsilon .
$$
On $\mathbb R^4$ the last identity is $QP(\varepsilon,T)=K_\varepsilon-K_T$ [C07 Lemma 10.2.1] with $K_\infty:=\lim_{T\to\infty}K_T=0$. That limit is an input (ledger G1); no source supplies $K_\infty$ on $\mathbb R^4$.

**1b. Row 6 is false.** All of the following fail together with $(Q^{GF})^2=0$:
- **[C07 Def 4.3.3(5)].** $\mathscr E=\operatorname{Im}Q\oplus\operatorname{Im}Q^{GF}\oplus\operatorname{Ker}H$ fails. In the field theory $Qc=-dc=Q^{GF}(\tfrac12dc[A^*])$ lies in both images (classical_spec §6). In the toy, $\dim(\operatorname{Im}Q\cap\operatorname{Im}Q^{GF})=1$.
- **Isotropy.** $\operatorname{Im}Q^{GF}$ is not isotropic. In the toy, $\langle Q^{GF}A^*,Q^{GF}c^*\rangle\neq0$, i.e. $\langle Q^{GF}x,Q^{GF}y\rangle=\pm\langle x,(Q^{GF})^2y\rangle\ne0$.
- **[C07 Lemma 6.5.1].** "$P(0,\infty)$ inverts $Q:\operatorname{Im}Q^{GF}\to\operatorname{Im}Q$" fails. In the toy, $P(0,\infty)Q\neq1$ on $\operatorname{Im}Q^{GF}$. What survives is only the homotopy identity $P(0,\infty)Q+QP(0,\infty)=1$, which needs $H=[Q,Q^{GF}]$ alone.
- **CS1, [C07 Lemma 6.6.2].** The Gaussian $\int_{(\operatorname{Im}Q^{GF})_{\mathbb R}}e^{\frac12\langle x,Qx\rangle/\hbar}$ is not defined: $\operatorname{Im}Q^{GF}$ contains all of $A$, and $\langle x,Qx\rangle$ vanishes on pure-gauge $A$. The positivity statement in classical_spec §7 holds only for the "Feynman completion" $-(N_{AA}+N')^{-1}$. That is an observation, not a source formula, and it is not $e^{\frac12\langle x,Qx\rangle}$ restricted to $L$.
- **[CG2 §7.2.1].** "Gauge fixing operator" requires square-zero, so the [CG2] definitions of theory and observable (7.2.9.1, 8.3.1.3) do not apply literally.

What does **not** need $(Q^{GF})^2=0$, all re-verified in the toy:
- $QP=K_\varepsilon-K_\infty$;
- $[\partial_P,Q]=\partial_{QP}$;
- the RG/QME compatibility [C07 Lemma 10.2.2];
- the UV-boundary identity (§2).

**1c. A consequence nobody tracked.** For a quadratic $S^{(1)}_{\mathfrak g}(X)$ (field × antifield), an insertion with both legs internal produces the chain $P^{AA}\to S^{(1)}_{\mathfrak g}(X)\to P^{A^*c}$. As an operator this is $Q^{GF}H^{-1}\,X\,Q^{GF}H^{-1}$, which equals $X(Q^{GF})^2H^{-2}\neq0$ when $X$ commutes with $Q^{GF}$. This vanishes for square-zero $Q^{GF}$ and is nonzero here. Any census of graphs with an $S^{(1)}_{\mathfrak g}(X)$ vertex (DESIGN route (ii), §5) must contain these $A\to c$ converting insertions.

---

## 2. UV-boundary identity (row 7)

Derivation chain [C07 §6.4, Lemma 10.2.1]:
$$
[\partial_P,Q]=\partial_{QP}=\partial_{K_\varepsilon}\ \Rightarrow\ Qe^{\hbar\partial_P}=e^{\hbar\partial_P}\big(Q-\hbar\partial_{K_\varepsilon}\big)=e^{\hbar\partial_P}(Q+\hbar\Delta_\varepsilon).
$$
- $[\partial_P,d_{CE}]=0$, because $P$ carries no ghost and $d_{CE}$ acts only on cochain arguments.
- **Machine.** `toy_gf.json` verifies $[\partial_{P(\varepsilon,\infty)},Q]=\partial_{K_\varepsilon}$ and $Qe^{\hbar\partial_P}F=e^{\hbar\partial_P}(Q+\hbar\Delta_\varepsilon)F$ on random mixed-parity polynomials. Here $Q=\{S_{free},-\}$, the [C07] contraction $\partial_\phi=\frac12\sum\partial_{\phi''}\partial_{\phi'}$ is used, and $(Q^{GF})^2\neq0$. The variant with $-\hbar\Delta_\varepsilon$ fails.
- **Verdict on the identity:** confirmed.

The DESIGN corollary "every anomaly / every singular coefficient therefore sits at $t=\varepsilon$" is **refuted as a consequence**:
1. The identity localises the action of **$Q$** on $K_\varepsilon$ edges. The localisation of **$\varepsilon$-derivatives** is a different and elementary fact, which produces $L_\varepsilon$, not $K_\varepsilon$:
$$G(\varepsilon)=\int_{[\varepsilon,\infty)^{E}}f\ \Rightarrow\ \varepsilon\tfrac{dG}{d\varepsilon}=-\textstyle\sum_e\varepsilon f\vert_{t_e=\varepsilon}.$$
   The dilatation is not $Q$.
2. Applying the identity to $e^{(I_{\mathfrak g}-I_{\mathfrak g}^{CT})/\hbar}$ produces $\{,\}_\varepsilon$, while the classical master equation holds with $\{,\}_0$. The tree remainders $\{V_1,V_2\}_\varepsilon-\{V_1,V_2\}_0$ multiply singular subgraphs and leave finite terms. For $x$-dependent $S^{(1)}_{\mathfrak g}(\Delta),S^{(1)}_{\mathfrak g}(K),S^{(1)}_{\mathfrak g}(S)$ these remainders are where the regulator breaks the symmetry, since the smearing $e^{-\varepsilon H}$ does not commute with $x$-weights. They must be part of any route-(ii) evaluation.

---

## 3. $O(\varepsilon)$ and $O[\infty]$ (rows 8–9)

**3a. Sign.** For the family $I+\delta O$, CS4 gives
$$(I+\delta O)^R=I+\delta O-I^{CT}_{\rm family}=(I-I^{CT}\vert_{\delta=0})+\delta\big(O-\partial_\delta I^{CT}\big),$$
which is consistent with DESIGN.

**3b. Source coverage (refuted as sourced).**
- [C07 Thm B] assumes $S$ "even, at least cubic mod ħ". $\delta O$ is quadratic for the bilinears and linear for the letters. The family version over $\mathbb C[\delta]/\delta^2$ is [CG2] text after Def 8.3.1.3 ("satisfies the axioms ... except locality"); it is not a theorem.
- The extension is straightforward ($\delta^2=0$ keeps the graph sums finite), but it is not stated or proved in the project.
- Univalent $O$ (letters) gives 1PR graphs $O$–$P$–(loop) and tadpoles. Their finiteness relies on the quadratic $I^{CT}$ and on colour vanishing (engines), and is not covered by Thm B.

**3c. Closedness not checked.** $O(\varepsilon)$ removes singular parts only. Whether $O$ lifts to a closed quantum observable needs the one-loop observable obstruction, i.e. the $K_\varepsilon$ graphs with one $O$ vertex. That obstruction may require a finite correction to $O$. It was not computed. The one-loop gauge QME (the obstruction of [C07 Prop 12.0.3]) was not computed either. Without it, $Q+\{I[\infty],-\}_\infty+\hbar\Delta_\infty$ is not shown to be a differential on the observables built from $I-I^{CT}$.

**3d. Hidden scale.** The scheme $\mathscr A_{<0}=\mathrm{span}\{\varepsilon^{-1},\log\varepsilon\}$ fixes the unit of length inside $\log\varepsilon$. This is not an IR regulator, but it is the only source of dilatation non-covariance that route (i) uses. It should be stated.

**3e. $O[\infty]$: refuted (definition gap).**
- With trivial coefficients, $(d_{CE}F)(X)=0$ for a cochain $F$ with no ghost argument ([ESW Def 1.26] formula; superconformal/spec §2 item 5). Hence
$$[(d_{CE}+Q)O[\infty]](X)=\pm\,Q\big([O[\infty]](X)\big).$$
- DESIGN defines $O[\infty]$ from $I-I^{CT}[\varepsilon]$, which is ghost-independent, so $[O[\infty]](X)=0$ and **$[(d_{CE}+Q)O[\infty]](X)\equiv0$ for every $X$, the dilatation included.**
- Step (c) as defined therefore computes zero. The engines' numbers come from a different construction: geometric scaling plus an ħ-corrected arity-1 action (assumption (B)).
- Required definition:
$$O[\infty]:=\lim_{\varepsilon\to0}\partial_\delta W\big(P(\varepsilon,\infty),\,I_{\mathfrak g}-I_{\mathfrak g}^{CT}[\varepsilon]+\delta O(\varepsilon)\big),\qquad I_{\mathfrak g}=I+\textstyle\sum_{k\ge1}S^{(k)}_{\mathfrak g}+(\hbar\text{-corrections})\ \ [\text{ESW Rem 1.29}].$$
  Here $I_{\mathfrak g}^{CT}[\varepsilon]$ includes the arity-1 counterterms, i.e. the singular parts of graphs with one $S^{(1)}_{\mathfrak g}(X)$ vertex. None of these were computed.
- [CG2 Def 8.3.1.3] defines an observable as a family $\{O[\Phi]\}$ over parametrices with proper support. Only the scale-∞ member is defined, since finite $L$ is forbidden, and $P(0,\infty)$ is not properly supported.
- The support condition [CG2 Def 8.4.0.1] is omitted. §4 shows why it matters.

---

## 4. "At scale ∞ the differential is $(d_{CE}+Q)$" (row 10)

1. **IR.** $\Delta_\infty F=\lim_{t\to\infty}(-\partial_{K_t}F)=0$ needs dominated convergence on the nonlocal $F$, which contains trees and loops with $P(0,\infty)\sim|p|^{-2}$. For trees with generic momenta this is fine. A general argument is absent (G1/G4). For $x^2$-weighted vertices see §6b.
2. **Acyclicity (main point).** With $K_\infty=0$:
$$QP(0,\infty)+P(0,\infty)Q=K_0=\mathrm{id}\quad\text{on }\mathscr E\quad(\text{toy: verified}).$$
   Hence $\big(\mathscr O(\mathscr E)/\mathbb C,\ Q\big)$ is acyclic. Machine: $\dim H=0$ in polynomial degrees 1, 2, 3 (`toy_acyclic.json`: ranks 6/12, 36/72, 146/292). Two consequences follow:
   - every $Q$-closed $O[\infty]$ is $Q$-exact up to a constant;
   - $[(d_{CE}+Q)O[\infty]](\Delta)=\pm Q([O[\infty]](\Delta))$ is $Q$-exact by construction.

   So a statement "$[(d_{CE}+Q)O[\infty]](\Delta)=-(E+w_O+a\,\cdot)O[\infty]$" holds **modulo $Q$-exact terms for every $a$**. The anomalous coefficient is not a cohomological invariant of the scale-∞ complex. It is fixed only as an identity between specific representatives, i.e. by the **locality at scale 0** of the equivariant data ($S^{(1)}_{\mathfrak g}(\Delta)$ and its ħ-correction). That locality is precisely the unverified assumption (B). The support condition of [CG2 Def 8.4.0.1] is what would make the complex non-trivial, and DESIGN omits it.
3. **Source coverage.** [CG2 Prop 13.4.0.1] is stated for the local symmetry algebra $\mathcal L$ and finite parametrices. [CG2] has no $T=\infty$, and its quantum global→local statement (§13.1.3) was not found (ledger). Using it at scale ∞ with $\mathfrak g=\mathfrak{psl}(4|4)$ is a definition, not a sourced result.

---

## 5. Dilatation: the two routes and their signs (rows 11–13)

**5a. Route (i): correct steps** (re-derived). The scaling map is $\phi(x)\mapsto\lambda^{w}\phi(\lambda x)$, and $(\varepsilon^{-1},\log\varepsilon)$ coefficients are those of $I^{CT}[\varepsilon]$:
$$W(P(\varepsilon,\infty),F)\circ\text{scaling}=W(P(\lambda^2\varepsilon,\infty),F\circ\text{scaling}),$$
$$(I-I^{CT}[\varepsilon])\circ\text{scaling}=I-I^{CT}[\lambda^2\varepsilon]+2\hbar\log\lambda\cdot(\log\varepsilon\text{ coefficient}).$$
- The $\varepsilon^{-1}$ coefficients have dimension 2 and drop out. Only $\log\varepsilon$ coefficients survive.
- The tree identity used by engine A, (insertion of $Q\Psi$ into $O$) = $Q$(insertion of $\Psi$) $-\{\Psi,O\}_\varepsilon$ for $QO=0$, follows from $[\partial_P,Q]=\partial_{K_\varepsilon}$ (§2).

**5b. Assumption (A), quadratic level: independently confirmed** (`toy_gf.json`):
$$Q\big(a\textstyle\int(A,A^*)-l\int(\lambda^*,\lambda)-g\int(c^*,c)\big)=-\tfrac a2(dA,dA)+l(\lambda,\not d\lambda)+(g-a)(dc,A^*),$$
so the engines' log coefficients $(-1,\,4,\,-1)$ give the unique weights $(a,l,g)=(2,4,1)$. The quadratic part of the coupling direction is itself $Q$-exact:
$$Q\big(\tfrac12\textstyle\int(A,A^*)-\tfrac12\int(\lambda^*,\lambda)-\tfrac12\int(c^*,c)\big),$$
which reproduces compare's one-parameter ambiguity.

**5c. Gaps in route (i).**
- **(A) beyond quadratic.**
  - Neither $(Q+\{I,-\})$-closedness nor exactness of the $\log\varepsilon$ coefficient is shown for the $AAA$, $\lambda\lambda A$, $A\vert_V(A\vert_{W^\perp})^2$, $c\lambda\lambda^*$, $cc\,c^*$ and quartic vertices.
  - The ghost vertex fixes the coupling direction to $0$ only inside $\mathrm{span}\{\text{exact},\text{coupling}\}$. Other closed, non-exact classes allowed by a non-SUSY regulator (separate Yukawa or quartic directions, double-trace quartics) are not excluded by any computation.
  - The weights in 5b predict the cubic and quartic log coefficients exactly: AAA $=6\times$ tree, $\lambda\lambda A=10\times$, quartic $A^4=8\times$, $c\lambda\lambda^*=cc\,c^*=1\times$. These are cheap falsification tests and none was run.
  - Note on robustness: the coefficient along $O$ uses only the quadratic weight of $A$ (verified) and the ghost vertex. The numbers $0$ and $-24$ survive an extra non-exact cubic or quartic class, but the equivariant QME along $\Delta$, i.e. the absence of a one-loop dilatation anomaly, would not.
- **(B).** "$S^{(1)}_{\mathfrak g}(\Delta)$ is corrected by $2\hbar\times\Psi$ so that the equivariant QME holds along $\Delta$" is asserted. The arity-1 equivariant obstruction was never computed. (B) follows from (A) only via the unproven identification below.
- **(C) Identification.** The claim is
$$[(d_{CE}+Q)O[\infty]](\Delta)=\pm Q([O[\infty]](\Delta))=(\text{geometric dilatation of }O[\infty])+(\text{tree insertion of the ħ-correction}).$$
  This is the one-loop dilatation Ward identity. Its proof is exactly DESIGN route (ii), which was not evaluated.

**5d. Route (ii) is not a second route.**
- Engine A (`faces.py`) and engine B (`route2.py`) both replace edges by $\frac{d}{d\lambda}P(\varepsilon/\lambda^2,\infty)=2\varepsilon L_\varepsilon$, i.e. they compute $\varepsilon\frac{d}{d\varepsilon}$ of the same graphs. That carries exactly the information already in the $\log\varepsilon$ coefficients used by route (i). They then assemble with the same (A) and (B).
- The agreement checks the $\varepsilon$-asymptotics, not the dilatation formula.
- DESIGN's route (ii) is the set of $K_\varepsilon$-edge graphs with an $S^{(1)}_{\mathfrak g}(\Delta)$ vertex plus the $\partial_\delta I^{CT}$–$K_\varepsilon$ trees (21+4 in the FeynArts census). It is **unevaluated** (compare.json "unresolved" #1). Beyond those 21+4, its census must also contain:
  1. $K_\varepsilon$ edges not adjacent to $S^{(1)}_{\mathfrak g}(\Delta)$;
  2. the $\{,\}_\varepsilon-\{,\}_0$ tree remainders of §2;
  3. the internal $S^{(1)}_{\mathfrak g}(\Delta)$ insertions $\propto(Q^{GF})^2$ of §1c.

  The manifest's route-(ii) labels show $S^{(1)}$ attached through $K_{AA^*}$ only; completeness is not demonstrated.

**5e. Signs and normalisation.**
- **Normalisation.** The engines use [SW (3)] $\Delta=-E$, acting as $-(E+w)$ on a weight-$w$ field. DESIGN requires [SW (17)], and superconformal/spec uses $\Delta=2E$. The conversion $\Delta_{(17)}=-2\Delta_{(3)}$ (ledger X13) appears only in compare. `results.json` of both engines is in the non-mandated normalisation.
- **Engine A prose sign.** Engine A's prose says "sum over special edges of a graph = $2\varepsilon\,d/d\varepsilon$ of the graph". From $\frac{d}{d\lambda}P(\varepsilon/\lambda^2,\infty)\vert_{\lambda=1}=2\varepsilon L_\varepsilon$:
$$\textstyle\sum_e G[e\mapsto2\varepsilon L_\varepsilon]=-2\varepsilon\frac{dG}{d\varepsilon}.$$
  The code (face value $=-\int g$) uses the correct sign; the prose is wrong. Engine B is consistent: $-2\varepsilon\,dP/d\varepsilon=2\varepsilon L_\varepsilon$ and $\varepsilon\,dG/d\varepsilon=-\sum(\text{pinned})$.
- **Relative signs.** The equivariant compensation enters $I$ and $O$ by the same rule (engine B: "$\Delta.F-Q(\text{correction})$" annihilates $I_{\mathfrak g}$ and acts on $O$). The relative sign between classical scaling and anomalous term is therefore convention-free. The value "traceless $=0$" is then a non-trivial check: the opposite compensation sign gives $8-16-8=-16$.
- **Absolute sign.** The absolute sign of a one-generator component depends on:
  1. the resolution of ESW's internal sign conflict (Def 1.26 printed vs (32)–(34)/Lemma 3.6, §6c);
  2. the placement of the ghost when "evaluating the cochain on $X$". DESIGN does not fix it.
- **Units (sanity only, not an input).** I re-derived the coupling normalisation from the Hodge-normalised $-\frac14(F_A,F_A)$: $(-,-)=\frac{2\hbar}{g^2}\times$ (unit form in a basis with Casimir $C_A$), hence $\hbar\times$(scalar of $X\mapsto\sum_a[t_a,[t^a,X]]$) $=-g^2C_A/2$. With this, trace $-24$ maps to $+3g^2C_A/(4\pi^2)$, consistent with compare.
- **Gauge parameter.** "Gauge-parameter independent in the sense available" cannot be tested: there is no gauge parameter in this $Q^{GF}$.

---

## 6. Selection rules and superconformal tower (rows 14–17)

**6a. $P,M,R$ (row 14).** $X$ commutes with $H$ and $Q^{GF}$, so $X.P(\varepsilon,\infty)=0$, $X.K_\varepsilon=0$, and $I^{CT}$ is $X$-invariant by uniqueness in Thm B. These facts are necessary. The claimed equality
$$Q([O[\infty]](X))=(X.O)[\infty]\quad\text{exactly}$$
also needs the arity-1 counterterms to vanish (or to be compatible) and the one-loop Ward identity. That identity includes the internal insertions of §1c, which exist only because $(Q^{GF})^2\ne0$. Neither is derived. The claim is plausible, but not established as "exactly".

**6b. $Q,\bar Q,S,\bar S,K$ (row 15): open.**
- Only structure-allowed lists were produced; nothing was computed.
- **$K$.** Engine B says the anomalous part is "fixed by $[K_\mu,P_\nu]=2M_{\mu\nu}+2\delta_{\mu\nu}E$" (and uses [SW (4)], not the mandated (17)). With trivial coefficients, the arity-2 component of $(d_{CE}+Q)^2O[\infty]=0$ only gives
$$[(d_{CE}+Q)O[\infty]]([K,P])=\mp Q\big([(d_{CE}+Q)O[\infty]](K,P)\big),$$
  and the $K$-component itself never enters. The step is a representation-theory heuristic, not a consequence of the complex. Engines A and B also disagree in wording.
- **IR (risk).** Consider an internal $S^{(1)}_{\mathfrak g}(K)$ insertion between $P^{AA}\sim|x|^{-2}$ and $P^{A^*c}\sim|x|^{-3}$, with the vertex weight $(x,x)u-2(u,x)x$ and one derivative. Large-$|x|$ power counting gives $\int d^4x\,|x|^{-4}$, i.e. a logarithmic IR divergence on $\mathbb R^4$ with $K_\infty=0$. For $\Delta$ and $S$ (linear in $x$) the count gives $|x|^{-5}$, which converges. No cancellation has been shown, and the "generic external momenta" rule does not control internal lines.

**6c. Tower (rows 16–17).**
- **Supported.**
  - The algebra: $30|32$; super-Jacobi on 84 type triples; conformal-Killing closure; [SW (17)] block grading; vanishing $\mathrm{str}(\mathrm{ad}\,\mathrm{ad})$.
  - Unique $S^{(1)}$ coefficients $(\frac14,\frac38,\frac14,-1)$.
  - Arity-2 closure for 31 ordered pairs; arity 3 for 24 patterns.
  - $S^{(3)}=0$; arity $\ge4$ structurally zero.
  - Negative controls fail as they should.
- **Refuted rigor claim.** `verify.py` states that the identity "holds identically with probability 1 (Schwartz–Zippel)".
  - Generator components are drawn from finite sets. For even generators the set is $\{\pm1,\pm2\}$ (`ri(nz=True)`), so $|S|=4$. For odd generators it is $\{-2..2\}+i\{-1..1\}$, so $|S|=15$.
  - Each arity-2 and arity-3 pattern is a single draw (one key per pattern in `verify.json`).
  - Schwartz–Zippel bounds the false-pass probability of a nonzero defect of degree $k$ by $k/|S|$. For example, it is $\le1/2$ for an even–even arity-2 pair.
  - The identities are multilinear, so exact verification on bases ($\dim=30|32$) is possible and was not done.
- **psl(4\|4).** The identification uses necessary conditions. There is no explicit isomorphism to the [SW (17)] supermatrices; the block rule tests only the grading.
- **Sign independence refuted.** spec §2 item 5 / §5 say the one-generator sign does not depend on the $d_{CE}$ convention. But spec §2 item 4 itself says that the printed [ESW Def 1.26] is the CE differential of the opposite algebra, via $x\mapsto-x$, $S^{(k)}\mapsto(-1)^kS^{(k)}$. Under the printed convention the MC-consistent arity-1 term is therefore $-S^{(1)}$ of [ESW (32)–(33)], and the classical part $\{S^{(1)}_{\mathfrak g}(X),O\}$ flips sign. The sign is a convention that must be stated together with $d_{CE}$.
- **"Involves $S^{(1)}(X)$ only" refuted at one loop.** The one-generator component also contains:
  - the ħ-corrections of the arity-1 action (assumption (B) for $\Delta$; unknown for $K,S,\bar S,Q,\bar Q$, all broken by the non-SUSY, non-conformal regulator);
  - the arity-1 counterterms.
- **Open.**
  - Uniqueness of $S^{(2)}$, i.e. the homotopy (acknowledged).
  - Gauge algebras other than $\mathfrak{su}(2)$ (acknowledged).
  - The quantum equivariant extension (ledger G19).

---

## 7. Minimal list of what would close the refutations

1. Redefine $O[\infty]$ with $I_{\mathfrak g}-I_{\mathfrak g}^{CT}[\varepsilon]$ (row 9). Compute the arity-1 counterterms along $\Delta$.
2. Evaluate DESIGN route (ii) completely, including the three classes in §5d. Compare it with route (i). This is the only independent check of (B) and (C).
3. Compute the cubic and quartic $\log\varepsilon$ coefficients and test them against the predictions in §5c.
4. Compute the one-loop gauge obstruction ([C07 Prop 12.0.3] at $T=\infty$) and the observable obstruction for $O$.
5. Restate all dilatation numbers in the [SW (17)] normalisation, with the $d_{CE}$ convention and the ghost-placement convention stated.
6. Replace the random-dense tower checks by basis-wise checks (multilinearity), or report the Schwartz–Zippel bounds.
7. Correct DESIGN's "the one axiom that fails" (§1b), and state the support condition of [CG2 Def 8.4.0.1] for scale-∞ observables (§4.2).
