# Classical psl(4|4)-equivariant BV action of reduced 4d N=4 SYM (Euclidean, complexified)

Sources: [C07], [CG2], [ESW], [SW] as keyed in `sources/ledger_SYNTHESIS.md`; the theory is `classical/classical_spec.md` (cited as [CS §n]).
Machine checks: `superconformal/verify.json` (section names in brackets), produced by `superconformal/verify.py`. The explicit realization is in `model.py`, `algebra.py` and `structure.py`, and the sign conventions are in `conventions.py`. Explicit bases (10d gamma blocks, su(2), spinor components) live only in the code.

**Result.** The action is
$$
S_{\mathfrak g}=S+S^{(1)}_{\mathfrak g}+S^{(2)}_{\mathfrak g},\qquad S^{(k\ge3)}_{\mathfrak g}=0,
$$
with $S^{(2)}_{\mathfrak g}$ nonzero only on pairs of odd generators. The equation $d_{CE}S_{\mathfrak g}+\frac12\{S_{\mathfrak g},S_{\mathfrak g}\}=0$ is machine-verified at every arity and for every pair and triple of generator types (§6).

---

## 0. Symbols

All symbols are ledger symbols. The table only recalls them.

| symbol | meaning | source |
|---|---|---|
| $V$, $W_{\mathbb R}^\perp$, $V_{10}=V\oplus W_{\mathbb R}^\perp$ | $\mathbb R^4$ directions, six transverse directions, 10d space | [ESW §2.6], DESIGN, [CS §1] |
| $x\in V$ | position on $\mathbb R^4$ | [CS §7] |
| $\Sigma,\Sigma^*$, $\rho$, $(-,-)$, $\Gamma$ | 10d semi-spin reps, Clifford action, pairing $\Sigma\otimes\Sigma^*\to\mathbb C$, Γ-pairing on $\mathrm{Sym}^2\Sigma$ and on $\mathrm{Sym}^2\Sigma^*$ | [ESW §2.1, (28)] |
| $q$ | antisymmetrisation $\wedge^\bullet V_{10}\to\mathrm{Cl}(V_{10})$, $q(v_1\wedge v_2)=v_1v_2-(v_1,v_2)$ | [ESW (29)] |
| $\rho(\omega)$, for a 2-form $\omega$ | $\rho(q(\omega))$, as in $\rho(F_A)$ | [ESW (33)] |
| $S_\pm$, $W$, $W^*$ | $\Sigma\cong S_+\otimes W\oplus S_-\otimes W^*$ | [ESW §10, §10.3] |
| $\not d$ | Dirac operator $\Gamma(\Sigma)\to\Gamma(\Sigma^*)$ | [ESW §2.1] |
| $A,\lambda,c,A^*,\lambda^*,c^*$, $A\vert_V$, $A\vert_{W^\perp}$, $\iota_wA$ | fields, antifields, restrictions | [ESW Def 3.2, Rem 1.36], [CS §1–2] |
| $S$, $F_A$, $d_A$, $\not d_A$, $\{-,-\}$ | action, curvature, covariant derivative, Dirac operator, bracket | [ESW (31)], [CS §3–4] |
| $(v,-)$ | metric dual one-form of a vector field $v$ | [ESW §2.1] |
| $L_v$, $v(\cdot)$, $\iota_v$ | Lie derivative, derivative along $v$, contraction | [ESW (32), (34)] |
| $E=x^\mu\partial_\mu$ | Euler vector field on $V$ | [SW Prop 2.2 (3)] |
| $\operatorname{div}v$ | divergence on $V$: $L_v\,\mathrm{dvol}=(\operatorname{div}v)\,\mathrm{dvol}$ | written out |
| $P,M,\bar M,\Delta,K,R,Q,\bar Q,S,\bar S$ | generator names | [SW (8), (17)] |
| $S^{(k)}_{\mathfrak g}(X_1,\dots,X_k)$, $d_{CE}$, $C^\bullet(\mathfrak g,\mathcal O_{\rm loc})$ | cochain components, CE differential, CE complex | [ESW Def 1.26, 1.27] |

Supercharge spinors are always written with a subscript, $Q_1,Q_2,\dots$ (ledger §5b: unsubscripted $Q$ is the [C07] differential). Following [ESW (28)], where $Q_1,Q_2\in\Sigma^*$ is allowed, $Q_i$ may be an element of $\Sigma$ or of $\Sigma^*$. Each formula states which one is meant. Generator-type labels $Q,\bar Q,S,\bar S$ in tables are the [SW (8)] names.

---

## 1. The algebra $\mathfrak g=\mathfrak{psl}(4|4)$: complexified Euclidean realization

**Even part.** $\mathfrak g_{\bar 0}$ is the conformal Killing vector fields on $V$ plus $\mathfrak{so}(W_{\mathbb R}^\perp)$. Each element is one vector field on $V_{10}$. The $\mathfrak{so}(W_{\mathbb R}^\perp)$ elements are the linear vector fields $y\mapsto ry$ on $W_{\mathbb R}^\perp$.

| [SW] name | vector field | [SW (17)] block | $\mathrm{ad}_\Delta$ weight |
|---|---|---|---|
| $P$ | $x\mapsto u$, $u\in V$ (translations) | $(\dot\alpha,\alpha)$ | $-2$ |
| $M$, $\bar M$ | $x\mapsto\omega x$, $\omega\in\mathfrak{so}(V)$; $M$ acts only on $S_+$, $\bar M$ only on $S_-$ | $(\alpha,\alpha)$, $(\dot\alpha,\dot\alpha)$ | $0$ |
| $\Delta$ | $2E$ | $\mathrm{diag}(1,-1\vert 0)$ | $0$ |
| $K$ | $x\mapsto (x,x)u-2(u,x)x$, $u\in V$ ($=u^\mu K_\mu$ of [SW (3)]) | $(\alpha,\dot\alpha)$ | $+2$ |
| $R$ | $\mathfrak{so}(W_{\mathbb R}^\perp)$ | $(i,i)$ | $0$ |

**Dilatation normalisation (one choice, ledger X13).** $\Delta$ is the $\Delta$ of [SW (17)], i.e. the vector field $2E$. Then $\mathrm{ad}_\Delta$ has eigenvalue $-2$ on $P$, $+2$ on $K$, $-1$ on $Q,\bar Q$, $+1$ on $S,\bar S$, and $0$ on $M,\bar M,R$ [alg: ad_Delta_weights]. In terms of the other [SW] normalisations, $\Delta=-2\times$ the $\Delta=-E$ of [SW (3)], and $\Delta=2\times$ the bosonic part of [SW (12)].

**Odd part.** $\mathfrak g_{\bar1}=\Pi(\Sigma\oplus\Sigma^*)$. An element $Q_1\in\Sigma$ is realised by the constant spinor $Q_1$. An element $Q_2\in\Sigma^*$ is realised by the $\Sigma$-valued function $x\mapsto\rho(x)Q_2$. Together these span the conformal Killing spinors $x\mapsto Q_1+\rho(x)Q_2$.

| [SW] name | subspace | [SW (17)] block |
|---|---|---|
| $Q^i_\alpha$ | $S_+\otimes W\subset\Sigma$ | $(i,\alpha)$ |
| $\bar Q_{\dot\alpha i}$ | $S_-\otimes W^*\subset\Sigma$ | $(\dot\alpha,i)$ |
| $\bar S^{\dot\alpha i}$ | $S_-\otimes W\subset\Sigma^*$ | $(i,\dot\alpha)$ |
| $S^\alpha_i$ | $S_+\otimes W^*\subset\Sigma^*$ | $(\alpha,i)$ |

- $S_+\otimes W$ is the $+1$ eigenspace of $\rho(q(\mathrm{dvol}))$ on $\Sigma$. Here $\mathrm{dvol}$ is the volume form of $V$ for the orientation of the orthonormal frame of `classical/n4_model.py`; the opposite orientation exchanges the names $Q\leftrightarrow\bar Q$, $S\leftrightarrow\bar S$, $M\leftrightarrow\bar M$.
- $S_-\otimes W$ is the $-1$ eigenspace of $\rho(q(\mathrm{dvol}))$ on $\Sigma^*$. It is $\rho(V)\cdot(S_+\otimes W)$.
- This fixes ledger G20: the [SW] R-index of $Q^i_\alpha$ is $W$.

**Brackets** (derived; all machine-checked):
$$
[v,w]=v(w)-w(v)\qquad(\text{vector fields on }V_{10}),
$$
$$
[v,\,Q_1+\rho(x)Q_2]=v(Q_1+\rho(x)Q_2)+\tfrac14\rho(d(v,-))(Q_1+\rho(x)Q_2)-\tfrac18(\operatorname{div}v)(Q_1+\rho(x)Q_2).
$$
The right-hand side is again of the form $Q_1'+\rho(x)Q_2'$.

For odd elements $Q_1+\rho(x)Q_2$ and $Q_3+\rho(x)Q_4$:
$$
[Q_1+\rho(x)Q_2,\;Q_3+\rho(x)Q_4]=\Big(x\mapsto\Gamma\big(Q_1+\rho(x)Q_2,\,Q_3+\rho(x)Q_4\big)\big|_V\Big)\;+\;r,
$$
where $r\in\mathfrak{so}(W_{\mathbb R}^\perp)$ is given by
$$
(w,r\,w')=(\rho(q(w\wedge w'))Q_1,Q_4)+(\rho(q(w\wedge w'))Q_3,Q_2),\qquad w,w'\in W_{\mathbb R}^\perp .
$$
- The coefficient $1$ of $r$ is the unique solution of super-Jacobi on $(Q,Q,S)$ [alg: cR_from_super_Jacobi].
- The $W_{\mathbb R}^\perp$-components of $\Gamma$ drop out of the bracket [ESW §2.6, Prop 2.20]. They reappear in $S^{(2)}$ (§4).

Special cases (machine-confirmed):
- $[Q_1,Q_2]=\Gamma(Q_1,Q_2)|_V$ for $Q_1,Q_2\in\Sigma$. This is [ESW Def 2.5] reduced.
- $[Q_1,Q_2]=-\big(x\mapsto(x,x)u-2(u,x)x\big)$ with $u=\Gamma(Q_1,Q_2)|_V$, for $Q_1,Q_2\in\Sigma^*$.
- For $Q_1\in\Sigma$ and $Q_2\in\Sigma^*$, $[Q_1,Q_2]=(x\mapsto\Gamma(Q_1,\rho(x)Q_2)|_V)+r$. Its dilatation part is $\tfrac12(Q_1,Q_2)\,\Delta$, and the rest is a rotation plus $r$.
- $[u,Q_1]=\rho(u)Q_1$ for $Q_1\in\Sigma^*$ and $u$ a translation; $[u,Q_1]=0$ for $Q_1\in\Sigma$.
- $[(x,x)u-2(u,x)x,\ Q_1]=\rho(u)Q_1\in\Sigma^*$ for $Q_1\in\Sigma$, and $=0$ for $Q_1\in\Sigma^*$.
- $[\Delta,Q_1]=-Q_1$ for $Q_1\in\Sigma$, and $+Q_1$ for $Q_1\in\Sigma^*$.

**Identification with $\mathfrak{psl}(4|4)$** [alg]:
- dimensions $30|32$;
- super-Jacobi for all 84 unordered type triples $\{P,M,\Delta,K,R,Q,S\}^3$;
- the $V$-part of every odd–odd bracket is a conformal Killing vector;
- the supermatrix block rule of [SW (17), (19)] holds for every pair of the ten named subspaces, i.e. $[\text{block}(a,b),\text{block}(c,d)]\subset\text{block}(a,d)\,[b=c]+\text{block}(c,b)\,[d=a]$;
- the Killing form $\mathrm{str}(\mathrm{ad}_X\mathrm{ad}_Y)$ vanishes, as it must for $\mathfrak{psl}(n|n)$.

---

## 2. Conventions of the cochain complex (locked by machine)

- $C^\bullet(\mathfrak g,\mathcal O_{\rm loc})=\bigoplus_k\mathrm{Hom}(\wedge^k\mathfrak g,\mathcal O_{\rm loc})[-k]$ [ESW Def 1.26].
- Antisymmetry: $f(\dots,x_i,x_{i+1},\dots)=(-1)^{|x_i||x_{i+1}|+1}f(\dots,x_{i+1},x_i,\dots)$. The source prints $|x_1||x_2|$ [sic].
- The value $S^{(k)}_{\mathfrak g}(x_1,\dots,x_k)$ has parity $k+\sum|x_i|$.
- $\{-,-\}$ is the bracket of [CS §3] ([C07 Lemma 4.3.2], $\langle x_i,\xi_i\rangle=+\int\mathrm{dvol}(x_i,\xi_i)$). It equals ESW's bracket: it reproduces [ESW Lemma 3.6 item 3] term by term, and the opposite sign does not [conventions].
- Cochain bracket: the [ESW Def 1.26] formula, with $|g|$ read as the degree of $g$ in $C^\bullet(\mathfrak g,\mathcal O_{\rm loc})$ ($=0$ for every component of $S_{\mathfrak g}$).
- $d_{CE}$ is the [ESW Def 1.26] formula with $|f|$ read as the degree in $C^\bullet(\mathfrak g,\mathcal O_{\rm loc})$, and with one extra overall sign:
$$
(d_{CE}f)(x_1,\dots,x_n)=\sum_{i<j}(-1)^{|x_i|\sum_{p<i}|x_p|+|x_j|\sum_{p<j,p\ne i}|x_p|+i+j+|f|+1}f([x_i,x_j],x_1,\dots,\widehat{x_i},\dots,\widehat{x_j},\dots,x_n).
$$

Why the extra sign [conventions]:
1. Read $|f|,|g|$ as either the degree of the multilinear map or the degree in $C^\bullet$. Only the reading $|f|=|g|=$ degree in $C^\bullet$ makes the printed formulas equivalent to a generating-function model on $\mathcal O(\mathfrak g[1])\otimes\mathcal O_{\rm loc}$ ($\mathfrak{gl}(2|1)$, all parity patterns, arity $\le3$).
2. In that reading the printed $d_{CE}$ gives $(d_{CE}S^{(1)})(Q_1,Q_2)=-S^{(1)}(\Gamma(Q_1,Q_2))$. [ESW Lemma 3.6 item 2] states $+S^{(1)}(\Gamma(Q_1,Q_2))$.
3. With ESW's bracket and ESW (33), (34), the printed Lemma 3.6 holds by machine, and the variant with $-S^{(1)}(\Gamma)$ fails.
4. The formula with the extra sign is the CE differential of $\mathfrak g$. The printed one is that of the opposite algebra, isomorphic via $x\mapsto-x$, $S^{(k)}\mapsto(-1)^kS^{(k)}$. Both pass the generating-function test.
5. Arity-1 components of $d_{CE}$ vanish, so the one-generator components used downstream do not depend on this sign.

The resulting component equations (machine-derived from the formulas above, `conventions.component_forms()`) are as follows. Write $\{S^{(k)},S^{(l)}\}$ for the bracket of values, $S=S^{(0)}_{\mathfrak g}$, and $S^{(k)}=S^{(k)}_{\mathfrak g}$.

- arity 1: $\{S,S^{(1)}(x)\}=0$;
- arity 2: $\{S,S^{(2)}(x_1,x_2)\}+S^{(1)}([x_1,x_2])+(-1)^{|x_1|}\{S^{(1)}(x_1),S^{(1)}(x_2)\}=0$;
- arity 3, $(x_1,x_2,x_3)$ = (even, odd, odd):
$$\{S,S^{(3)}\}+S^{(2)}([x_1,x_2],x_3)+S^{(2)}([x_1,x_3],x_2)+S^{(2)}([x_2,x_3],x_1)+\{S^{(1)}(x_1),S^{(2)}(x_2,x_3)\}-\{S^{(1)}(x_2),S^{(2)}(x_1,x_3)\}-\{S^{(1)}(x_3),S^{(2)}(x_1,x_2)\}=0;$$
- arity 3, (odd, odd, odd): the same with all signs $+$. This is [ESW Lemma 3.7]'s cyclic form;
- arity 4, (odd)$^4$:
$$\{S,S^{(4)}\}+\textstyle\sum_{i<j}S^{(3)}([x_i,x_j],\dots)-\sum_i\{S^{(1)}(x_i),S^{(3)}(\dots)\}+\{S^{(2)}(x_1,x_2),S^{(2)}(x_3,x_4)\}+\{S^{(2)}(x_1,x_3),S^{(2)}(x_2,x_4)\}+\{S^{(2)}(x_1,x_4),S^{(2)}(x_2,x_3)\}=0.$$

ESW Thm 3.3's four equations are the case $\mathfrak g=\mathfrak A$ of these.

---

## 3. $S^{(1)}_{\mathfrak g}$ (Noether functionals)

**Even generators** ($v$ a conformal Killing vector field on $V$, or $v\in\mathfrak{so}(W_{\mathbb R}^\perp)$ as a vector field on $V_{10}$):
$$
S^{(1)}(v)=\int\mathrm{dvol}\Big(\big(L_vA+\tfrac14(\operatorname{div}v)\,A|_{W^\perp},\,A^*\big)-\big(v.\lambda,\lambda^*\big)-\big(v(c)\big)c^*\Big),
$$
$$
v.\lambda=v(\lambda)+\tfrac14\rho(d(v,-))\lambda+\tfrac38(\operatorname{div}v)\lambda .
$$
- $A$ is regarded as a one-form on $V_{10}$ constant along $W_{\mathbb R}^\perp$ [ESW Def 1.55, Prop 1.56], and $L_v$ is its Lie derivative on $V_{10}$.
- For $v$ on $V$: $L_v(A|_V)$ is the Lie derivative of the one-form, and $L_v(\iota_wA)=v(\iota_wA)$.
- For $v=(y\mapsto ry)\in\mathfrak{so}(W_{\mathbb R}^\perp)$: $\iota_w(L_vA)=\iota_{rw}A$.
- $\operatorname{div}v$ is $0$ for $P,M,R$; it is $8$ for $\Delta=2E$, and $-8(u,x)$ for $x\mapsto(x,x)u-2(u,x)x$.
- For $v\in\mathfrak{iso}(V)\oplus\mathfrak{so}(W_{\mathbb R}^\perp)$ this is [ESW (32)]. ESW's "so(V) action on Σ" is $\tfrac14\rho(d(v,-))$.
- The three coefficients $\tfrac14$ (spin), $\tfrac38$ (Weyl weight $\tfrac32$ of $\lambda$) and $\tfrac14$ (Weyl weight $1$ of $A|_{W^\perp}$) are the unique solutions of $\{S,S^{(1)}(v)\}=0$ [arity1_invariance: derived_coefficients]. The Weyl weight of $A|_V$ (1) and of $c$ (0) is automatic.
- Explicitly, $S^{(1)}(\Delta)=\int\mathrm{dvol}\big((2E(A)+2A,A^*)-(2E(\lambda)+3\lambda,\lambda^*)-(2E(c))c^*\big)$, with $E(A)$ componentwise.

**Odd generators, $Q_1\in\Sigma$** [ESW (33), reduced, Prop 2.20]:
$$
S^{(1)}(Q_1)=\int\mathrm{dvol}\Big(-(\Gamma(Q_1,\lambda),A^*)+\tfrac12(\rho(F_A)Q_1,\lambda^*)\Big).
$$

**Odd generators, $Q_2\in\Sigma^*$** (conformal Killing spinor $x\mapsto\rho(x)Q_2$):
$$
S^{(1)}(Q_2)=\int\mathrm{dvol}\Big(-(\Gamma(\rho(x)Q_2,\lambda),A^*)+\tfrac12(\rho(F_A)\rho(x)Q_2,\lambda^*)-(\rho(A|_{W^\perp})Q_2,\lambda^*)\Big).
$$
- Since $\not d(\rho(x)Q_2)=4Q_2$, the last term is $-\tfrac14(\rho(A|_{W^\perp})\not d(\rho(x)Q_2),\lambda^*)$. In this form one formula covers both cases on $x\mapsto Q_1+\rho(x)Q_2$.
- The coefficient $-1$ is the unique solution of $\{S,S^{(1)}(Q_2)\}=0$ [arity1_invariance].
- $\rho(A|_{W^\perp})$ is Clifford multiplication by the $W_{\mathbb R}^\perp$-component of $A$, with the $\mathfrak g$ factor kept.

**$c^*$ terms in $S^{(1)}$:** only $-(v(c))c^*$ for even $v$. The functionals $S^{(1)}(Q_1)$ ($Q_1\in\Sigma$ or $\Sigma^*$) contain no $c^*$.

---

## 4. $S^{(2)}_{\mathfrak g}$ (quadratic-antifield terms) and $S^{(k\ge3)}_{\mathfrak g}$

$S^{(2)}(x_1,x_2)=0$ if $x_1$ or $x_2$ is even. On odd generators, $S^{(2)}$ is [ESW (34)] evaluated pointwise on the conformal Killing spinors. For $Q_1,Q_2$, each in $\Sigma$ or in $\Sigma^*$, replace every $Q_i\in\Sigma^*$ by the $\Sigma$-valued function $x\mapsto\rho(x)Q_i$ in
$$
S^{(2)}(Q_1,Q_2)=\int\mathrm{dvol}\Big(\tfrac14\big(\Gamma(Q_1,Q_2),\Gamma(\lambda^*,\lambda^*)\big)-\tfrac12(Q_1,\lambda^*)(Q_2,\lambda^*)-\iota_{\Gamma(Q_1,Q_2)}A\,c^*\Big).
$$
- $\Gamma(Q_1,Q_2)$ is the full $V_{10}$-valued Γ-pairing (ten components).
- Its $W_{\mathbb R}^\perp$-components give the scalar-dependent part of the gauge parameter.

**Field-dependent gauge terms.** The $c^*$ terms of $S^{(2)}$ are exactly $-\iota_{\Gamma(Q_1,Q_2)}A\,c^*$, with gauge parameter $\iota_{\Gamma}A=\iota_{\Gamma|_V}(A|_V)+\iota_{\Gamma|_{W^\perp}}A$:
- for $Q_1,Q_2\in\Sigma$: as in ESW (constant $\Gamma$);
- for $Q_1\in\Sigma$, $Q_2\in\Sigma^*$: $\Gamma(Q_1,\rho(x)Q_2)$, linear in $x$;
- for $Q_1,Q_2\in\Sigma^*$: $\Gamma(\rho(x)Q_1,\rho(x)Q_2)$, with $V$-part $-((x,x)u-2(u,x)x)$, $u=\Gamma(Q_1,Q_2)|_V$, and $W_{\mathbb R}^\perp$-part $-(x,x)\,\Gamma(Q_1,Q_2)|_{W^\perp}$.

**$S^{(3)}=0$.** Every arity-3 equation holds with $S^{(3)}=0$ (§6). All higher arities then vanish:
- $S^{(2)}$ contains only $A$, $\lambda^*$ and $c^*$, whose bracket partners $A^*,\lambda,c$ are absent, so $\{S^{(2)},S^{(2)}\}=0$ identically;
- the arity-4 equation reduces to this;
- arity $\ge5$ has no terms.

**Not addressed:** uniqueness of $S^{(2)}$ modulo $\{S,-\}$-closed degree-2 local functionals, i.e. the choice of homotopy.

---

## 5. Physical reading (for downstream use)

- **Even–odd pairs.** $S^{(2)}=0$, so the conformal and R transformations, with the Weyl weights above, commute with the supercharges off shell. The arity-2 equation is $S^{(1)}([v,Q_1])=-\{S^{(1)}(v),S^{(1)}(Q_1)\}$.
- **Odd–odd pairs.** $\{S^{(1)}(Q_1),S^{(1)}(Q_2)\}$ equals $S^{(1)}([Q_1,Q_2])+\{S,S^{(2)}(Q_1,Q_2)\}$. Here $[Q_1,Q_2]$ is the conformal Killing vector plus the $\mathfrak{so}(W_{\mathbb R}^\perp)$ rotation $r$ of §1. The term $\{S,S^{(2)}\}$ contains:
  - the fermion equation of motion: the $\lambda^*\lambda^*$ terms;
  - the field-dependent gauge transformation with parameter $\iota_{\Gamma(Q_1,Q_2)}A$: the $c^*$ term.
- **One-generator components downstream.** $[(d_{CE}+Q)O[\infty]](X)$ involves $S^{(1)}(X)$ only. Its sign does not depend on the $d_{CE}$ convention of §2.

---

## 6. Machine verification (`verify.json`)

Realisation: $\mathfrak g=\mathfrak{su}(2)$ gauge algebra with $(T_i,T_j)=\delta_{ij}$, all jets, exact Gaussian rationals. Each check is an identity of local functionals modulo total derivatives: all Euler–Lagrange derivatives of the density vanish identically. Generators are random and dense, with small Gaussian-integer components and fixed seeds; all identities are multilinear in the generators.

| equation | generator types verified | status |
|---|---|---|
| $\{S,S^{(1)}(x)\}=0$ | $P,M,\Delta,K,R,Q,S$ (2 random each) | holds [arity1_invariance] |
| coefficient derivations | spin $\tfrac14$ (from $M$ and from $R$); $(\tfrac38,\tfrac14)$ from $\Delta$; $(\tfrac14,\tfrac38,\tfrac14)$ from $K$; $-1$ from $S$ | unique solutions [arity1_invariance] |
| arity 2 (closure up to gauge and EOM, with $S^{(2)}$) | all 28 unordered pairs of $\{P,M,\Delta,K,R,Q,S\}$, plus the reorderings $SQ$, $KP$, $SK$ (31 checks) | all hold [arity2_closure] |
| arity 3 (Jacobi level, $S^{(3)}=0$) | all 15 (even, odd, odd) type patterns; $QQQ$, $QQS$, $QSS$, $SSS$; reorderings $QKS$, $SQ\Delta$; patterns with $\le1$ odd argument ($PKQ$, $KKS$, $MR\Delta$: structurally zero) (24 checks) | all hold [arity3_jacobi] |
| arity 3 non-vacuity | for $QQQ$, $QQS$, $QSS$, $SSS$ the cyclic sum cancels already at density level (pointwise Fierz identity, [ESW Assumption 3.1]); each single term $\{S^{(1)}(x_i),S^{(2)}(x_j,x_k)\}$ is nonzero | [arity3_nonvacuity] |
| arity 4 | $QSQS$, $SSQQ$, $KQSS$: every term vanishes identically ($\{S^{(2)},S^{(2)}\}=0$ structurally, $S^{(3)}=0$) | all hold [arity4] |
| explicit component equations | arity 2, 3, 4 equations of §2 generated from the adopted sign formulas | [component_equations] |
| negative controls | the S-compensator coefficient set to $0$; wrong λ Weyl weight; spin coefficient $-\tfrac14$; no Weyl weight on $A\vert_{W^\perp}$; $r=0$ in $[Q,S]$; $S^{(2)}$ without the $W_{\mathbb R}^\perp$-components of Γ (on $QQ$ and $QS$); $S^{(2)}=0$ on $QS$ | each must fail [arity1, arity2 negative_controls] |
| algebra | super-Jacobi (84 type triples), CKV closure, [SW (17)] block rule, $\mathrm{ad}_\Delta$ weights, vanishing Killing form, dimension $30|32$ | [algebra] |
| conventions | readings of [ESW Def 1.26]; ESW bracket sign (Lemma 3.6 item 3); ESW Lemma 3.6 as printed; graded symmetry of $\{-,-\}$ | [conventions] |

Open items (not verified):
- uniqueness of $S^{(2)}$, i.e. the space of homotopies;
- random-dense checks were done for $\mathfrak g=\mathfrak{su}(2)$ only. The identities use only the Jacobi identity and invariance of $(-,-)$ (any $\mathfrak g$ by construction), but other gauge algebras were not run.
