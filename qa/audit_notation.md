# Notation and sources audit

Auditor lens: notation and sources. Reference: `DESIGN.md` hard rules 1–5, `sources/ledger_SYNTHESIS.md` (§5b collisions, §6 gaps).
Scope: every `.md`, `.tex`, `.json` under `classical/`, `engine_a/`, `engine_b/`, `superconformal/`, `compare/` (no `.tex` exists). For the "other local projects" check, `.py` and `.log` files were grepped as well.
Severity: **H** = breaks a DESIGN hard rule, or a ledger symbol is reused for a different object inside one reported entry. **M** = a new or colliding symbol in a reported file. **L** = cosmetic, or a code-level artefact.

## 0. Verdict

| check | result |
|---|---|
| references to or reuse of other local projects | **PASS**: every absolute path points inside the project (`compare/compare_tables.py:6`, `compare/make_compare.py:3`, `*/logs`, `superconformal/out/log*.txt`); there are no imports from outside the project. FeynArts appears only through the in-project `diagrams/manifest.json` and `diagrams/data/census_comparison.json` (`compare/compare_tables.py:9,11`). A grep of `diagrams/` finds no external paths or project names. |
| IR scale (T, L, masses, upper t-cutoff) | **PASS**: `T` appears only in the ledger/DESIGN citation phrase "[C07 Lemma 10.2.1 at T=∞]" (`classical/classical_spec.md:161`, `engine_a/results.json:2`). $K_\infty=0$ is stated (`classical_spec.md:142,241`). No `L`, mass or finite upper Schwinger limit appears in md/json/py. |
| dimensional regularization | **PASS**: none in md/json/py. |
| Q-hat / $Q_1$-style shorthand | **PASS for Q**: $Q_1,Q_2,\dots$ are always [ESW (28)] spinors (`superconformal/spec.md:36`, `classical_spec.md:17-18`). **FAIL elsewhere**: `O_1` (engine_b), `W_1` (engines), `Delta_(17)` (compare). See V27, V34, V49. |
| symbols not in the ledger or not written out | **FAIL**: 14 H, 24 M, 12 L (§2). |
| explicit spinor indices / gamma components / Lie bases in reported text | **FAIL**: V1, V2, V23, V30, V42, V43. Code-level JSON is listed separately (§3). |
| DESIGN–ledger conflicts needing a decision | 4 (§1) |
| citations not in the ledger | 5 locations, all verified in the ESW text (§4) |

## 1. DESIGN versus ledger conflicts (orchestrator decision; engines follow DESIGN)

| # | symbol | conflict | occurrences | compliant resolution |
|---|---|---|---|---|
| D1 | $W$ | DESIGN uses $W(P,I)$ for the RG flow [CG2 §7.2.7]. Ledger §5b reserves $W$ for the [ESW] multiplicity space and names the flow [C07] $\Gamma$. | Both meanings in `classical/classical_spec.md:19` vs `:212-217`; `classical/census.md:1`; `engine_a/results.json:2`, `engine_b/results.json:2,582` (`W_1`) | Keep $W(\cdot,\cdot)$ only with its two arguments and amend ledger §5b, **or** write $\Gamma(P(\varepsilon,\infty),I)$ [C07 Def 6.6.1]. Never write bare `W_1`; write "the $\hbar^1$ part of $W(P(\varepsilon,\infty),I)$". |
| D2 | $\mathfrak g$ | DESIGN uses $\mathfrak g$ for the gauge algebra (Fields) and also sets $\mathfrak g=\mathfrak{psl}(4\vert4)$ (Superconformal). Ledger §5b has no $\mathfrak g$ entry. | `superconformal/spec.md:7-10,34,40,106-133` ($\mathfrak{psl}(4\vert4)$) vs `:206,223` (gauge $\mathfrak{su}(2)$) in one document | $\mathfrak g$ = gauge algebra only. Write the symmetry algebra out: $S_{\mathfrak{psl}(4\vert4)}$, $S^{(k)}_{\mathfrak{psl}(4\vert4)}$, $C^\bullet(\mathfrak{psl}(4\vert4),\mathcal O_{\rm loc})$. Add a §5b row. |
| D3 | $P,M,\bar M,K,R,Q,\bar Q,S,\bar S$ unindexed | Ledger §5b: unindexed $S$ = action/interaction, $Q$ = [C07] differential, $K_{t}$ = kernel, two-slot $P$ = propagator, superscript $R$ = renormalised; generators "always carry indices". DESIGN rule 2 forbids spinor indices. | superconformal, engine_a, engine_b, verify JSON | Index-free generator names built from ledger symbols: translation by $u\in V$; $\mathfrak{so}(V)$ rotation; $\Delta$; special conformal field $x\mapsto(x,x)u-2(u,x)x$; element of $\mathfrak{so}(W_{\mathbb R}^\perp)$; $Q_1\in S_+\otimes W$ or $S_-\otimes W^*$ ($\subset\Sigma$); $Q_2\in S_-\otimes W$ or $S_+\otimes W^*$ ($\subset\Sigma^*$, realised as $x\mapsto\rho(x)Q_2$). This is the form `superconformal/spec.md:54-61` already uses. |
| D4 | $S$ | DESIGN: $S=\tfrac12\langle e,Qe\rangle+I$ (full action, the ESW sense). Ledger §5b/X1: unindexed $S$ = [C07] interaction. | all files (consistently DESIGN) | No file change. Amend ledger §5b/X1 to the DESIGN reading. |

## 2. Violations in reported files

### 2.1 `superconformal/spec.md`

| id | sev | line | found | compliant replacement |
|---|---|---|---|---|
| V1 | H | 33, 44-50, 56-61, 97, 210-218 | [SW] spinor-index generator names $Q^i_\alpha,\bar Q_{\dot\alpha i},\bar S^{\dot\alpha i},S^\alpha_i$ and block labels $(\dot\alpha,\alpha)$, $(i,\alpha)$ appear in reported tables. Unindexed generator $S$ appears in the same row as the action $S$: `:210` "$\{S,S^{(1)}(x)\}=0$ \| $P,M,\Delta,K,R,Q,S$"; `:211` "$-1$ from $S$"; `:217` "$S$-compensator", "$S^{(2)}=0$ on $QS$". | Use the D3 names, e.g. `:210` "$\{S,S^{(1)}(x)\}=0$ for $x$ a translation, rotation, $\Delta$, special conformal field, $\mathfrak{so}(W_{\mathbb R}^\perp)$ element, $Q_1\in\Sigma$, $Q_2\in\Sigma^*$"; `:211` "$-1$ from $Q_2\in\Sigma^*$". Keep the [SW (8)] names only as one cross-reference column, labelled as a quotation. |
| V2 | H | 206 | Lie-algebra basis in reported text: "$(T_i,T_j)=\delta_{ij}$" | "$\mathfrak g=\mathfrak{su}(2)$ with $(X,Y)=-\tfrac12\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)$" (basis in `superconformal/model.py:10`; $\operatorname{tr}\,\mathrm{ad}\,\mathrm{ad}=-2\,\mathrm{id}$ in `classical/checks.json:33`) |
| V3 | M | 42, 80-84, 197 | $r\in\mathfrak{so}(W_{\mathbb R}^\perp)$ collides with [SW (17)] $r$ (ledger §4, SC4) | "the $\mathfrak{so}(W_{\mathbb R}^\perp)$ component of $[Q_1+\rho(x)Q_2,\,Q_3+\rho(x)Q_4]$", defined by the `:82` pairing |
| V4 | M | 69 vs 82 | $w$ is a vector field in $[v,w]=v(w)-w(v)$ but $w\in W_{\mathbb R}^\perp$ in `:82` (same §1) | $[v,v']=v(v')-v'(v)$ |
| V5 | L | 31, 49 | vector-index forms $E=x^\mu\partial_\mu$, $u^\mu K_\mu$ | "$E$ = Euler vector field, $E(x)=x$ [SW Prop 2.2 (3)]"; "the [SW (3)] special conformal field along $u$" |
| V6 | L | 42 | $y\mapsto ry$: $y$ collides with the second position $y$ (`classical_spec.md:138`) | "the linear vector field on $W_{\mathbb R}^\perp$ given by an element of $\mathfrak{so}(W_{\mathbb R}^\perp)$" |
| V7 | L | 42, 54 | $\mathfrak g_{\bar0},\mathfrak g_{\bar1}$ (not in ledger) | "the even (odd) part of $\mathfrak{psl}(4\vert4)$" |
| V8 | L | 99 | $\text{block}(a,b)$, Iverson bracket $[b=c]$ | write "$\subset\text{block}(a,d)$ if $b=c$, plus $\text{block}(c,b)$ if $d=a$" |

### 2.2 `classical/classical_spec.md`

| id | sev | line | found | compliant replacement |
|---|---|---|---|---|
| V9 | M | 234 | $\Delta_O$ for the dimension of $O$. Ledger §5b: subscripted $\Delta$ = BV Laplacian, unsubscripted = dilatation. | "the dimension of $O$" in words |
| V10 | M | 232-234; also `census.md:7,20,50,59,61` | $\omega$ = superficial degree, colliding with the [ESW] BV pairing $\omega$ (ledger §5b) | "superficial degree" (column header "degree") |
| V11 | M | 194-210 | Polarisation symbols collide: $v\in\Sigma^*\otimes\mathfrak g$ vs $v\in V_{10}$ (`:16`); $w\in\mathfrak g$ vs $w\in W_{\mathbb R}^\perp$ (`:12`); $\gamma\in\mathfrak g$ vs graph $\gamma$ (`:217`, same §9); $\alpha\in V_{10}\otimes\mathfrak g$ vs [SW] index / [CG2] anomaly $\alpha$; $a$ vs [ESW §10.2] field $a$ | Evaluate on the field letters with momentum labels, e.g. "$A\,c\,A^*$: $r_cr_*\,([A(p_1),c(p_2)],A^*(p_3))$", "$c\,c\,c^*$: $r_1r_2\,([c(p_1),c(p_2)],c^*(p_3))$" |
| V12 | M | 164-167 | $N$, $N_{\lambda\lambda}$, $N_{AA}$, $N'$ (new) | "the graded Hessian of $\tfrac12\langle e,Qe\rangle$ on the $\lambda$ (resp. $A$) sector"; "the Hessian of $-\tfrac14(d^*A\vert_V,d^*A\vert_V)$" |
| V13 | M | 140, 163-167 | $\tilde K_t(p)$, $\tilde P^{AA},\tilde P^{\lambda\lambda},\tilde P^{cA^*},\tilde P^{A^*c}$ (not in ledger) | "the $A\otimes A$ (…) component of $P(\varepsilon,\infty)$ in momentum space", as in the `:148-153` table |
| V14 | M | 171-177; `census.md:11` | $I_3,I_4$: single-index shorthand colliding with [CG2 Def 7.2.9.1] $I_{i,k}$ ($\hbar$ order, degree) | $I_{0,3}$, $I_{0,4}$ (or "cubic / quartic part of $I$") |
| V15 | M | 12, 34, 80-85, 105-106 and every later use; `census.md`, `spec.md:27,141-165` | $A\vert_{W^\perp}$: DESIGN prescribes $\iota_wA$ or "the $W_{\mathbb R}^\perp$ component of $A$", and $W^\perp$ drops the $\mathbb R$, colliding with the ESW $W$ | $\iota_wA$, $w\in W_{\mathbb R}^\perp$; at minimum $A\vert_{W_{\mathbb R}^\perp}$ |
| V16 | L | 29 | $F$ for [ESW Rem 1.36] $F_{\rm gauge}$, colliding with $F_A$ | $F_{\rm gauge}$ |
| V17 | L | 70 | $\mathrm{dvol}_{10}$ (not in ledger) | "the volume form of $V_{10}$" |
| V18 | L | 103, 115, 122, 132 | slot notation "$[A\vert_V]$", colliding with $[-,-]$ and the shift $[1]$ | "(in the $A\vert_V$ slot)" |
| V19 | L | 224-228 | $\mathrm{str}$, while [CG2 §13.7.2] has $\mathrm{STr}$ (engine_b uses STr) | $\mathrm{STr}$ |
| V20 | L | 106 | momentum $p$ used before its definition at `:138` | move the Fourier convention before §5 |

### 2.3 `classical/census.md`

| id | sev | line | found | compliant replacement |
|---|---|---|---|---|
| V21 | M | 10-13, 22-35, 52-59 | Vertex labels $V_3^{VVV},V_3^{VWW},V_4^{\dots},Y^V,Y^W,G^V,G^W$ are new. $V$ (vertex) collides with $V$ ($\mathbb R^4$ directions). Superscript $W$ means $W_{\mathbb R}^\perp$ but collides with ESW $W$ ($\dim 4$). $G$ collides with ESW $G_R,G_W$. | Name each vertex by its term: "$(A\vert_V)^3$ part of $-\tfrac12(dA,\tfrac12[A\wedge A])$", "$\lambda\lambda\,\iota_wA$ part of $\tfrac12(\lambda,(\not d_A-\not d)\lambda)$", "$A\vert_V\,c\,A^*\vert_V$ part of $([A,c],A^*)$", … |

### 2.4 `engine_a/results.json`

| id | sev | line | found | compliant replacement |
|---|---|---|---|---|
| V23 | H | 2, 754, 803-804, 836; 805 | `c_adj` (DESIGN: "no new symbol for a Casimir"); basis `t_a`, `t^a` in reported text (`:2`); `f_abb` (`:805`) | unit $\hbar\,\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)/\big((X,Y)\,16\pi^2\big)$ (as at `compare/compare.json:572`); `f_abb` → "$\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X)=0$" |
| V24 | H | 826-837, 841-844 | Dilatation in $\Delta=-E$ ([SW (3)]). DESIGN fixes [SW (17)], i.e. $\Delta=2E$ (`superconformal/spec.md:48-52`). Weights "Q +1/2, S −1/2, K −1" (`:843-844`) contradict `spec.md:52` (Q −1, S +1, K +2). | Restate in $\Delta=2E$: component $=2(E+\text{dim}\,O+a\,u)O$, where $u$ is the V23 unit written out. The coefficient of $u\,O$ is $2a$: traceless 0, trace −48, $\iota_wA$ −8, $A\vert_V$ (transverse) −8, $\lambda$ −16, $c$ −4. Weights: Q, Q̄ −1; S, S̄ +1; K +2. |
| V25 | H | 826, 836 | `w` = Weyl weight inside the entry that uses `iota_w A` (`:830`); `D_lambda`, `lambda^w phi(lambda x)` (λ = gaugino, φ = ESW §10.2 field); `gamma_O` (γ = [C07] graph); `Psi` ([C07] Thm A $\Psi_{q,r}$); `O_{omega-dot}` undefined; `F` vs $F_A$ | "Weyl weight" in words; scale factor written as a word or a non-ledger letter, e.g. "$x\mapsto\ell x$"; `Psi` → $2\int\mathrm{dvol}(A,A^*)-4\int\mathrm{dvol}(\lambda^*,\lambda)-\int\mathrm{dvol}(c^*,c)$; `gamma_O` → "coefficient of $u\,O$ in $[(d_{CE}+Q)O[\infty]](\Delta)$" |
| V26 | H | 16-520 (all `eps_inv`/`log_eps` keys), 757-790 | Labels `B1`–`B7` are undefined in the file; they are defined only in code (`engine_a/sector_a.py:59-97`) | written-out functionals, §5 R1 |
| V27 | M | 2 | `W_1`; `S_free`; `Phi(x)` ([C07] $\Phi$); `y_a y_b` (basis) | "the $\hbar^1$ part of $W(P(\varepsilon,\infty),I)$"; $\tfrac12\langle e,Qe\rangle$; "a field"; "two coordinate functionals" |
| V28 | M | 610, 636, 662, 686, 712, 738, 751, 793-821, 843 | `int omega` uses $\omega$ as the [CG2 Lemma 10.1.1.3] background density, uncited; ledger §5b reserves $\omega$ for the ESW pairing | $\int\mathrm{dvol}$, or cite [CG2 Lemma 10.1.1.3] at each use and add the meaning to §5b |
| V29 | M | 804, 842-844 | `O_ww`; `w_O` (collides with $w\in W_{\mathbb R}^\perp$ and $w_\gamma$) | $\sum$ over $W_{\mathbb R}^\perp$ written as the trace with the inverse metric of $W_{\mathbb R}^\perp$; "the dimension of $O$" |
| V30 | M | 843-844 | `Q^i_alpha` (spinor index); rep labels `6 x 4 (or 4bar)`; unindexed `S, Sbar`, `K` | $Q_1\in S_+\otimes W$; "$W_{\mathbb R}^\perp\otimes W$ (or $\otimes W^*$)"; D3 names |

### 2.5 `engine_b/results.json`

| id | sev | line | found | compliant replacement |
|---|---|---|---|---|
| V31 | H | 2; all 22 `units` lines (24 … 432); 515, 526, 536, 584 | `c_adj`, basis `t_a`, `t^a`, `t^` | as V23 |
| V32 | H | 506, 574-578 | $\Delta=-E$; "Q dimension +1/2, K −1, S −1/2" | as V24 |
| V33 | H | 515 | `g^2 C_A`, "standard coupling", "Killing-normalised": new coupling and Casimir symbols (textbook sanity check) | delete, or state it with $\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)$ and $(-,-)$ only |
| V34 | H | 510, 535 | `R` = the local functional $2\int(A,A^*)+\dots$ (collides with SW $R$, ESW $R$ and superscript $R$; R/C as operator names are banned by the user); `O_1` ($Q_1$-style subscript shorthand); `I^(1)(Delta)`; `xi-linear` (ξ collides with [C07] $\xi_i$, `classical_spec.md:46-55`; ghost coordinates have no source name); `D_lambda^*` | `R` → the functional written out (V25); `O_1` → "the tree-level local functional $-2\hbar\,(\cdots)$" written out; `I^(1)(Delta)` → $S^{(1)}_{\mathfrak{psl}(4\vert4)}(\Delta)$ plus its $\hbar$ correction [ESW Def 1.27, DESIGN]; `xi-linear …` → "the component along $\Delta$ (arity 1) of the equivariant QME" |
| V35 | H | 569 | `Sigma` = self-energy (collides with ESW $\Sigma$) | "the one-loop two-point part of $W(P(\varepsilon,\infty),I)$" |
| V36 | M | 521, 530, 535-537, 569, 587 | `b_k`, `B_k`, `b_4`, `b_6`, `b'_O`, `k_X`, `b_X`, `(d_delta I^CT)_log`, `d^cross_{P(0,infty)}(B_k,O)` | coefficients written out; `d^cross_P(B,O)` → $\partial_{P(0,\infty)}(BO)-(\partial_{P(0,\infty)}B)O-B\,\partial_{P(0,\infty)}O$ ([C07 §6.4] symbols only) |
| V37 | M | 508-537, 575 | `gamma`, `gamma_O`, `gamma_letter`, `kappa`, `kappa_classical_Delta_on_functional` | words, as V25 |
| V38 | M | 526, 536, 569 | `u = hbar c_adj/16pi^2`; `u` is also the Schwinger ratio `u_f` (`:2`) | unit written out (V23) |
| V39 | M | 2, 582, 596-687 | `Pi` = parity operator (ESW $\Pi$ = parity shift) | $(-1)^{\vert\cdot\vert}$ (`classical_spec.md:224`) |
| V40 | M | 588, 598-653 | `C0`, `C1`, Feynman parameter `x`: $x$ collides with position $x$ and [C07] $x_i$; Feynman parameters are in no source | Schwinger times $t_e$ [C07 F8] with the ratio written out, $t_1/(t_1+t_2)$ |
| V41 | M | 2, 582 | `W_1`, `N(q)`, `N_AA`, `N_{lam lam}`, `Phi(x)` | as V27, V12 |
| V42 | M | 2 | Explicit gamma components in a reported string: `sigma = (1, i sigma_k)`, `Gamma^mu = gamma^mu (x) 1, Gamma^{4+i} = gamma_5 (x) gamma_6^i`, `'t Hooft blocks` | move to the `engine_b/realization.py` docstring; keep only "explicit realization: engine_b/realization.py" |
| V43 | M | 574-578 | `P_mu`, `M_munu`, `R^i_j`, `K_mu`, `x_mu`, `delta_munu`, `Q^i_alpha`, `Qbar_{alphadot i}`, `S^alpha_i`, `Sbar^{alphadot i}`, `{K,Q} ~ S` | D3 index-free names; the [SW (4)] relation stated for vector fields |
| V44 | L | 358-432, 471-502; 494; 502; 574 | `int omega` (V28); `O_tr` (defined in place); `zeta` (defined in place, acceptable); letter test element `X` vs generator `X` (`:574`) | V28; trace written out; letter element written out |
| V45 | L | 556 | fit constants `a,b,c,d` (`c` vs ghost) | numeric text only; rename in prose |

### 2.6 `compare/compare.json` (and `compare/compare_tables.json`, same labels)

| id | sev | line | found | compliant replacement |
|---|---|---|---|---|
| V46 | H | 3, 827 | "compare/report.md writes them out": `compare/report.md` does not exist. The written-out forms exist only at `:572-581`. | Write `compare/report.md`, or cite `compare.json:572-581` |
| V47 | H | 565-567, 817-819 | `konishi`: `g^2 C_A`, `N`, "standard" | as V33 |
| V48 | H | 804-815 | $\Delta=-E$ values; conversion documented at `:815` but values not restated | restate per V24 |
| V49 | M | 572, 799-800, 805-815, 392-418, 478-503, 561 | `u`, `a` (collides with ESW §10.2 $a$), `w_O`, `O_ww`, `O_tr`, `Delta_(17)`, `Delta_(3)`, `gamma` keys, `B_kin`, `S_free`, `d_{A^*}` | written out; `Delta_(17)` → "$\Delta$ of [SW (17)]" |
| V50 | M | 262-330, 582-590 | `B1`–`B7` as table labels. Compliant only next to `:573-581`; never inside formulas (cf. `engine_a/bv_exact.py:1` "L_2 = -(B1+B4)+4B6-B7") | §5 R1 |
| V51 | L | 740-831 | `int omega` | V28 |

## 3. Code-level JSON (allowed as code output; must not be copied into reports)

- `classical/checks.json`: 13-27 (gamma labels `IIIX…`, `sig^mu`, `Om4`, `e_mu`); 71-160 (`kappa` = scale of $(-,-)$, which DESIGN keeps symbol-free; `N_AA`, `N_GF`, `N_L`); 168-171 (`h_t`, `G_eps`, `Delta_+`, colliding with $\Delta$); 175-180 (`theta1 theta2`, θ = [SW] superspace coordinate, vs $r_j$ in `classical_spec.md:195`); 187 (`W_1`).
- `classical/census.json`: labels `g3[VVV]`, `y[W]`, `gh[V]`, `AW`; 117-221 `f(X,i,i)` (Lie index).
- `engine_a/logs/*.json`: `Psi`, `zA/zc/zl`, `L2`, `gll`, component tuples `('A',4,1)`, `p0..p3`.
- `engine_b/report_data.json`, `report_data_O.json`, `two_point_all.json`, `o_graphs.json`, `route2_numeric.json`: `u0 v0 w1_0 p0`, `x`, `C0/C1`, `Pi`, `I_3`. The name `report_data` suggests report input; if it is typeset, apply V39/V40.
- `superconformal/verify.json`, `out/*.json`: `S1/S2/S3`, `kl/wl/wA/aS/cR`, `D#0`; key `/arity1_invariance/{S, S^(1)(x)} = 0/S#0` uses `S` both as the action and as a generator (D3).

## 4. Sources and citations

| id | item | status |
|---|---|---|
| S1 | [ESW Def 1.55, Prop 1.56] (`classical_spec.md:21,29`; `spec.md:146`), [ESW Lemma 3.4] (`classical_spec.md:72`), [ESW Lemma 3.6 item 3] (`spec.md:109`), [ESW Lemma 3.7] (`spec.md:129`) | present in the ESW text (`esw.txt` l.1842, 1867, 2647, 2684, 2768); absent from all ledgers → add to `sources/ledger_esw-sym-bv.md` |
| S2 | `classical_spec.md:21` cites [ESW Def 1.55] for $d^*$ and dvol | Def 1.55 defines dimensional reduction only; mark $d^*$ and dvol as "written out" |
| S3 | `spec.md:65` "This fixes ledger G20" | fixed by the orientation of the code frame (`:63`); label it as a program convention, not a source statement |
| S4 | `spec.md:107,111-121`: sign deviations from the printed [ESW Def 1.26] (antisymmetry, extra overall sign of $d_{CE}$) | documented with machine evidence; add as a ledger §5a X-entry |
| S5 | Ledger §5b omissions | add rows for $\mathfrak g$ (D2), $\omega$ (ESW §2.1 two-forms used at `spec.md:24,47`; CG2 background density used by the engines), $r$ (V3), $x$ ([C07] $x_i$ / position / Feynman parameter), $W$ (D1) |
| S6 | out of lens, recorded | `compare.json:827` unresolved[4]: engine A two-point and O coefficients are float64 fits with rational reconstruction, not exact arithmetic |

## 5. Compliant restatements

**R1. Basis functionals** (from `engine_a/sector_a.py:59-97`, which agrees with `compare.json:573-581`). Pairings are Hodge $\otimes(-,-)$:
- B1 $=\int\mathrm{dvol}\,(dA\vert_V,dA\vert_V)$
- B2 $=\int\mathrm{dvol}\,(d^*A\vert_V,d^*A\vert_V)$
- B3 $=\int\mathrm{dvol}\,(A\vert_V,A\vert_V)$
- B4 $=\int\mathrm{dvol}\,(dA\vert_{W_{\mathbb R}^\perp},dA\vert_{W_{\mathbb R}^\perp})$
- B5 $=\int\mathrm{dvol}\,(A\vert_{W_{\mathbb R}^\perp},A\vert_{W_{\mathbb R}^\perp})$
- B6 $=\int\mathrm{dvol}\,(\lambda,\not d\lambda)$
- B7 $=\int\mathrm{dvol}\,(dc,A^*\vert_V)$

The quadratic counterterm, written with no new symbols:
$$
I^{CT}[\varepsilon]\big\vert_{\rm quadratic}=\frac{\hbar\,\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{(X,Y)\,16\pi^2}\Big[\tfrac94\,\varepsilon^{-1}\big(\text{B3}+\text{B5}\big)+\log\varepsilon\,\big(-\text{B1}-\text{B4}+4\,\text{B6}-\text{B7}\big)\Big]
$$
Each B$k$ is to be substituted by its integral above.

**R2. Dilatation in the DESIGN normalisation ([SW (17)], $\Delta=2E$; ledger X13: $\Delta_{(17)}=-2\Delta_{(3)}$).** $[(d_{CE}+Q)O[\infty]](\Delta)$ equals $2\big(E+\text{dim}\,O\big)O$ plus $\dfrac{\hbar\,\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{(X,Y)\,16\pi^2}\,k\,O$, with $k$:
- traceless symmetric $(\iota_{w_1}A,\iota_{w_2}A)$: 0
- trace: −48
- $\iota_wA$: −8
- $A\vert_V$ (transverse): −8
- $\lambda$: −16
- $c$: −4

The letters are Feynman-gauge only.
