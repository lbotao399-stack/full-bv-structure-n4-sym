# Classical BV data — 4d N=4 SYM on Euclidean R^4 (dimensional reduction of [ESW (31)])

Sources: [C07], [CG2], [ESW], [SW] as keyed in `sources/ledger_SYNTHESIS.md`. Machine checks: `classical/checks.json` (section numbers §1–§9 below refer to it), produced by `classical/checks.py` from the explicit realization `classical/n4_model.py`, `classical/action.py` and `classical/jets.py`. Explicit bases (gamma blocks, Lie-algebra basis, spinor components) occur only in the code.

---

## 1. Spaces

| symbol | meaning | source |
|---|---|---|
| $V_{10}=V\oplus W_{\mathbb R}^\perp$ | 10d Euclidean space; $V$ = the $\mathbb R^4$ directions, $W_{\mathbb R}^\perp$ = 6 transverse directions | [ESW §2.6]; DESIGN |
| $A\vert_V$, $A\vert_{W^\perp}$, $\iota_wA$ | restriction of $A\in V_{10}^\vee\otimes\mathfrak g$ to $V$, to $W_{\mathbb R}^\perp$; $\iota_wA=A(w)$, $w\in W_{\mathbb R}^\perp$ | DESIGN |
| $\mathfrak g$, $(-,-)$, $[-,-]$ | complexified gauge Lie algebra, nondegenerate invariant symmetric pairing (positive definite on the compact real form; carries the coupling), bracket | [ESW §3.1]; G12(a) |
| $(-,-)^{-1}$ | inverse pairing, element of $\mathfrak g\otimes\mathfrak g$ | written out |
| $\Sigma=S_+^{10}$, $\Sigma^*=S_-^{10}$ | 16-dim semi-spin reps of $\mathfrak{so}(V_{10})$, $M=\Sigma\oplus\Sigma^*$ | [ESW §2.1, §3.1] |
| $\rho(v)$ | Clifford action, $\rho(v)\rho(v')+\rho(v')\rho(v)=2(v,v')$, $\rho(v):\Sigma\to\Sigma^*$ and $\Sigma^*\to\Sigma$ | [ESW §2.1] |
| $(-,-):\Sigma\otimes\Sigma^*\to\mathbb C$ | $(\rho(v)Q_1,Q_2)=(Q_1,\rho(v)Q_2)$ | [ESW §2.1] |
| $\Gamma$ | $(\rho(v)Q_1,Q_2)=(v,\Gamma(Q_1,Q_2))$ | [ESW (28)] |
| $\Sigma\to S_+\otimes W\oplus S_-\otimes W^*$ | 4d decomposition, $\dim W=4$, $W\cong S_+^6$ of $\mathrm{Spin}(6;\mathbb C)$ | [ESW §10, §10.3] |
| $\not d=\rho\circ d$ | Dirac operator $\Gamma(\Sigma)\to\Gamma(\Sigma^*)$; the same formula on $\Gamma(\Sigma^*)\to\Gamma(\Sigma)$ is written $\not d\vert_{\Sigma^*}$ | [ESW §2.1] |
| $d$, $d^*$, $\mathrm{dvol}$ | exterior derivative along $\mathbb R^4$ only ($\iota_w d=0$), its formal adjoint, Euclidean volume form | [ESW Def 1.55] |

Machine facts (§1): $\rho(v)\rho(v')+\rho(v')\rho(v)=2(v,v')$ on $\Sigma$ and $\Sigma^*$; $\rho(v)$ symmetric for $(-,-)$; [ESW Assumption 3.1] $\rho(\Gamma(Q_1,Q_2))Q_3+\text{cyclic}=0$ holds; 4d chirality on $\Sigma$ has eigenspaces $8+8$; $\rho(v)$, $v\in V$, exchanges the 4d chiralities ($S_+\otimes W\leftrightarrow S_-\otimes W$ type); $\rho(w)$, $w\in W_{\mathbb R}^\perp$, preserves them ($S_+\otimes W\leftrightarrow S_+\otimes W^*$ type). Euclidean signature forces one non-real Clifford generator; it is placed in $W_{\mathbb R}^\perp$, so $\not d$ on $\mathbb R^4$ is real.

---

## 2. Graded space of fields $\mathscr E$

$\mathscr E=F\oplus F^![-1]$, $F=\Omega^0(\mathbb R^4;\mathfrak g[1])\oplus\Omega^0(\mathbb R^4;V_{10}^\vee\otimes\mathfrak g)\oplus\Omega^0(\mathbb R^4;\Pi\Sigma\otimes\mathfrak g)$, restricted to fields constant along $W_{\mathbb R}^\perp$ [ESW Def 3.2, Rem 1.36, Def 1.55, Prop 1.56, Prop 2.20]. Duals are identified through $(-,-)$ on $\mathfrak g$, the metric of $V_{10}$, $(-,-):\Sigma\otimes\Sigma^*\to\mathbb C$ and $\mathrm{dvol}$.

| field | fibre | $\mathbb Z$-degree [ESW Rem 1.33] | $\mathbb Z/2$ | total parity | $\mathrm{Spin}(4)\times\mathrm{Spin}(6)$ content |
|---|---|---|---|---|---|
| $c$ | $\mathfrak g$ | $-1$ | even | odd | scalar |
| $A$ | $V_{10}^\vee\otimes\mathfrak g$ | $0$ | even | even | $A\vert_V$: vector of $\mathrm{Spin}(4)$; $A\vert_{W^\perp}$: vector of $\mathrm{Spin}(6)$ (six scalars) |
| $\lambda$ | $\Pi\Sigma\otimes\mathfrak g$ | $0$ | odd | odd | $S_+\otimes W\oplus S_-\otimes W^*$ |
| $A^*$ | $V_{10}\otimes\mathfrak g$ | $+1$ | even | odd | as $A$ |
| $\lambda^*$ | $\Pi\Sigma^*\otimes\mathfrak g$ | $+1$ | odd | even | $S_-\otimes W\oplus S_+\otimes W^*$ |
| $c^*$ | $\mathfrak g$ | $+2$ | even | even | scalar |

Koszul signs use the total parity [ESW §1.1, ledger X10].

---

## 3. Pairing, bracket, $K_0$

Normalisation [C07 §10.1, §2.1]: $\langle x_i,\xi_i\rangle=+\int\mathrm{dvol}\,(x_i,\xi_i)$ with $x_i$ the even member and $\xi_i$ the odd member of each canonical pair,
$$
(x_i,\xi_i)\in\{(A,A^*),\ (\lambda^*,\lambda),\ (c^*,c)\},\qquad \langle\xi_i,x_i\rangle=-\langle x_i,\xi_i\rangle ,
$$
all other pairings zero; $\langle,\rangle$ is odd, graded antisymmetric [C07 Def 4.3.1; ESW Def 1.13, 1.16].
$$
\Delta_0=\sum\partial_{x_i}\partial_{\xi_i}\ (\text{contracted with }(-,-)^{-1}),\qquad \{f,g\}=\Delta_0(fg)-\Delta_0(f)g-(-1)^{|f|}f\Delta_0(g)\quad[\text{C07 §10.1, Lemma 4.3.2}],
$$
$$
K_0=-\sum\,(x_i\otimes\xi_i+\xi_i\otimes x_i)\otimes(-,-)^{-1}\ \ (\text{supported on the diagonal of }\mathbb R^4\times\mathbb R^4)\quad[\text{C07 §10.1}].
$$
Machine facts (§2): for this normalisation $\{S,S\}=0$ holds exactly ($\mathfrak g=\mathfrak{su}(2)$, all jets). The three other relative signs of the $(\lambda,\lambda^*)$ and $(c,c^*)$ pairs against $(A,A^*)$ violate the CME. The overall sign of $\langle,\rangle$ is not fixed by the CME. It is fixed here by $\langle x_i,\xi_i\rangle=+\int\mathrm{dvol}\,(x_i,\xi_i)$.

---

## 4. Action

$$
S=\tfrac12\langle e,Qe\rangle+I=\int_{\mathbb R^4}\mathrm{dvol}\Big(-\tfrac14(F_A,F_A)+\tfrac12(\lambda,\not d_A\lambda)+(d_Ac,A^*)+([c,\lambda],\lambda^*)+\tfrac12([c,c],c^*)\Big)\quad[\text{ESW (31), Def 1.55, Prop 2.20}],
$$
$$
F_A=dA+\tfrac12[A\wedge A],\qquad d_A=d+[A,-],\qquad \not d_A=\rho\circ d_A,
$$
- $d$ acts along $\mathbb R^4$ only. $d_A$ and $\rho\circ d_A$ use all ten components of $A$.
- $(F_A,F_A)$ is the Hodge pairing on $\wedge^2V_{10}^\vee$ tensored with $(-,-)$, i.e. $(F_A,F_A)\,\mathrm{dvol}_{10}=(F_A\wedge *F_A)$.
- Every pairing is read with its left entry first, as printed in (31).
- The normalisation of $(F_A,F_A)$ is fixed by [ESW (33), Lemma 3.4] (§3): $\{S,S^{(1)}(Q)\}=0$ holds for the reduced $\mathfrak{su}(2)$ theory with
$$
S^{(1)}(Q)=\int\mathrm{dvol}\big(-(\Gamma(Q,\lambda),A^*)+\tfrac12(\rho(F_A)Q,\lambda^*)\big).
$$
  It holds only with coefficient $-\frac14$ of the Hodge pairing. The coefficients $+\frac14$ and $-\frac14\times2$ (i.e. $-\frac14\sum_{\text{all }M,N}$) fail, and it also fails with the uniform bracket.

Reduced form, $A=A\vert_V+A\vert_{W^\perp}$:
$$
(F_A,F_A)=(F_{A|_V},F_{A|_V})+(d_{A|_V}A\vert_{W^\perp},d_{A|_V}A\vert_{W^\perp})+(\tfrac12[A\vert_{W^\perp}\wedge A\vert_{W^\perp}],\tfrac12[A\vert_{W^\perp}\wedge A\vert_{W^\perp}]),
$$
$$
\not d_A\lambda=\not d\lambda+\rho([A\vert_V,\lambda])+\rho([A\vert_{W^\perp},\lambda]),\qquad (d_Ac,A^*)=(dc,A^*\vert_V)+([A,c],A^*).
$$
Here $\rho([A,\lambda])$ means Clifford multiplication by the $V_{10}^\vee$ factor of $A$ together with $[-,-]$ on the $\mathfrak g$ factors, i.e. $(\not d_A-\not d)\lambda$.

Free part:
$$
\tfrac12\langle e,Qe\rangle=\int\mathrm{dvol}\Big(-\tfrac14(dA,dA)+\tfrac12(\lambda,\not d\lambda)+(dc,A^*)\Big).
$$
Real slice [C07 §6.6, CS1] (§6): real $A$ with respect to the compact real form, on which $(-,-)$ is positive definite. Then
$$
\langle A,QA\rangle=-\tfrac12\int\mathrm{dvol}\,(dA,dA)\le0,
$$
negative definite modulo $\operatorname{Im}Q$ (eigenvalues, in units of $(-,-)$: $-\frac12|p|^2$ with multiplicity 9, and $0$ with multiplicity 1, per colour).

---

## 5. $Q$ on every field

$Q$ is the linear part of $\{S,-\}$, i.e. $\{\tfrac12\langle e,Qe\rangle,-\}$ (machine-derived, §4). It acts on $\mathscr E$, with targets in square brackets:
$$
Q c=-dc\ [A\vert_V],\qquad QA=-\tfrac12\,d^*dA\ [A^*],\qquad QA^*=d^*(A^*\vert_V)\ [c^*],\qquad Q\lambda=\not d\lambda\ [\lambda^*],\qquad Qc^*=0,\quad Q\lambda^*=0 .
$$
- $d^*dA\vert_{W^\perp}=d^*d(A\vert_{W^\perp})$, the Laplacian on the six scalars.
- In momentum space $\partial\to ip$: $QA=-\frac12\big(|p|^2-p\,\iota_p\big)A$ on $V$ and $-\frac12|p|^2A$ on $W_{\mathbb R}^\perp$.
- $Q^2=0$ (§2, §4), and $(Q\otimes1+1\otimes Q)K_0=0$ (§4).

---

## 6. Gauge fixing $Q^{GF}$, $H$

$Q^{GF}$ has the form fixed in DESIGN. Its coefficients are the unique solution of $H=[Q,Q^{GF}]=d^*d\otimes\mathrm{id}$ (§4):
$$
Q^{GF}A^*=-2A^*\ [A],\qquad Q^{GF}A=-d^*(A\vert_V)\ [c],\qquad Q^{GF}c^*=dc^*\ [A^*\vert_V],\qquad Q^{GF}\lambda^*=-\not d\vert_{\Sigma^*}\lambda^*\ [\lambda],\qquad Q^{GF}c=Q^{GF}\lambda=0 .
$$
- $Q^{GF}A^*$ is placed in the $A$ slot via the metric of $V_{10}$, all ten components. The factor $-2$ is forced by $-\frac14$ with the Hodge $(F_A,F_A)$.
- $H=[Q,Q^{GF}]=d^*d\otimes\mathrm{id}$ on all 54 components per colour (§4). $d^*d$ on $\Omega^0(\mathbb R^4)$ is the positive Laplacian.
- $[Q^{GF},H]=0$ (§4).
- $(Q^{GF})^2\neq0$ (§4):
$$
(Q^{GF})^2A^*=2\,d^*(A^*\vert_V)\ [c],\qquad (Q^{GF})^2c^*=-2\,dc^*\ [A\vert_V],
$$
  and it is zero on $c,A,\lambda,\lambda^*$.
- $(Q^{GF}\otimes1-1\otimes Q^{GF})K_0=0$, i.e. $Q^{GF}$ has the adjointness opposite to $Q$ [C07 Lemma 6.3.1, Def 4.3.3(1) with §4.3; ledger X8] (§4).

Status of [C07 Def 4.3.3]:
- (1) holds;
- (2) holds;
- (3) fails;
- (4) holds;
- (5) fails as well: $Qc=-dc=Q^{GF}(\tfrac12dc\ [A^*])\in\operatorname{Im}Q\cap\operatorname{Im}Q^{GF}$ (§4). $\operatorname{Im}Q^{GF}$ is not Lagrangian. On $\mathbb R^4$, $\operatorname{Ker}H$ also falls under G1.

---

## 7. $K_t$ and $P(\varepsilon,\infty)$

Fourier convention: a field is $\int\frac{d^4p}{(2\pi)^4}e^{ip\cdot x}(\cdots)(p)$ and a kernel on $\mathbb R^4\times\mathbb R^4$ is $\int\frac{d^4p}{(2\pi)^4}e^{ip\cdot(x-y)}(\cdots)(p)$; $x,y\in\mathbb R^4$ are positions. The heat kernel is
$$
K_t=-(4\pi t)^{-2}e^{-|x-y|^2/4t}\sum(x_i\otimes\xi_i+\xi_i\otimes x_i)\otimes(-,-)^{-1},\qquad \tilde K_t(p)=-e^{-t|p|^2}(\cdots),\qquad K_t\star e=e^{-tH}e\quad[\text{C07 §6.5}],
$$
with the identity of the fibre in each slot: the metric of $V_{10}$ for $A\otimes A^*$, the canonical element of $\Sigma\otimes\Sigma^*$ for $\lambda\otimes\lambda^*$. $K_\infty=0$ on $\mathbb R^4$.
$$
P(\varepsilon,\infty)=\int_\varepsilon^\infty(Q^{GF}\otimes1)K_t\,dt\qquad[\text{C07 §6.5}],\qquad \int_\varepsilon^\infty(4\pi t)^{-2}e^{-|x|^2/4t}\,dt=\frac{1-e^{-|x|^2/4\varepsilon}}{4\pi^2|x|^2}\quad(§7).
$$
Every nonzero component (§5). The component in the row "$X\otimes Y$" is the value of $\partial_P$ on the product of a coordinate functional of $X$ at $x$ and one of $Y$ at $y$, in that order (§9):

| component | position space | momentum space |
|---|---|---|
| $A\otimes A$ | $2\,(-,-)^{-1}\otimes(\text{metric of }V_{10})^{-1}\ \dfrac{1-e^{-\lvert x-y\rvert^2/4\varepsilon}}{4\pi^2\lvert x-y\rvert^2}$ | $2\,(-,-)^{-1}\otimes(\text{metric})^{-1}\ \dfrac{e^{-\varepsilon\lvert p\rvert^2}}{\lvert p\rvert^2}$ |
| $c\otimes A^*$ (only $A^*\vert_V$) | $-(-,-)^{-1}\ d_x\dfrac{1-e^{-\lvert x-y\rvert^2/4\varepsilon}}{4\pi^2\lvert x-y\rvert^2}$ | $-i\,(-,-)^{-1}\ p\ \dfrac{e^{-\varepsilon\lvert p\rvert^2}}{\lvert p\rvert^2}$ |
| $A^*\otimes c$ (only $A^*\vert_V$) | $-(-,-)^{-1}\ d_x\dfrac{1-e^{-\lvert x-y\rvert^2/4\varepsilon}}{4\pi^2\lvert x-y\rvert^2}$ | $-i\,(-,-)^{-1}\ p\ \dfrac{e^{-\varepsilon\lvert p\rvert^2}}{\lvert p\rvert^2}$ |
| $\lambda\otimes\lambda$ | $(-,-)^{-1}\otimes\not d_x\vert_{\Sigma^*}\dfrac{1-e^{-\lvert x-y\rvert^2/4\varepsilon}}{4\pi^2\lvert x-y\rvert^2}$ | $i\,(-,-)^{-1}\otimes\rho(p)\vert_{\Sigma^*\to\Sigma}\ \dfrac{e^{-\varepsilon\lvert p\rvert^2}}{\lvert p\rvert^2}$ |

The kernel $\rho(p)\vert_{\Sigma^*\to\Sigma}$ is read as an element of $\Sigma\otimes\Sigma$ through $(-,-):\Sigma\otimes\Sigma^*\to\mathbb C$. All other components ($A\otimes A^*$, $\lambda\otimes\lambda^*$, $c\otimes c^*$, $c\otimes A^*\vert_{W^\perp}$, …) vanish. Explicitly
$$
d_x\frac{1-e^{-|x-y|^2/4\varepsilon}}{4\pi^2|x-y|^2}=(x-y)\,\frac{\frac{|x-y|^2}{2\varepsilon}e^{-|x-y|^2/4\varepsilon}-2\big(1-e^{-|x-y|^2/4\varepsilon}\big)}{4\pi^2|x-y|^4},\qquad \frac{1-e^{-|x|^2/4\varepsilon}}{4\pi^2|x|^2}\Big|_{x=0}=\frac1{16\pi^2\varepsilon}\quad(§7).
$$
Facts (§5–§7):
- Graded symmetry: exchanging the two tensor factors (and $x\leftrightarrow y$, i.e. $p\to-p$) multiplies each component by the Koszul sign of the two factors.
- $QP(\varepsilon,\infty)=K_\varepsilon$, with $Q=Q\otimes1+1\otimes Q$, exactly [C07 Lemma 10.2.1 at $T=\infty$]. The position-space core is $d^*d\int_\varepsilon^\infty=(4\pi\varepsilon)^{-2}e^{-|x|^2/4\varepsilon}$.
- $(Q\otimes1+1\otimes Q)K_t=0$ and $\frac{d}{dt}K_t+(Q\otimes1+1\otimes Q)L_t=0$ [C07 §6.5].
- Positivity: $|p|^2\tilde P^{AA}$ has eigenvalue $2(-,-)^{-1}$ (all ten directions), positive.
- Hessian inversion at $\varepsilon\to0$, with $N$ the graded Hessian of $\frac12\langle e,Qe\rangle$ and covariance $=-N^{-1}$:
  - $\lambda$: $\tilde P^{\lambda\lambda}=-N_{\lambda\lambda}^{-1}$ exactly.
  - $A$: $\tilde P^{AA}(-N_{AA})=$ the projector killing $\operatorname{Im}Q$ (the longitudinal part of $A\vert_V$). Equivalently $\tilde P^{AA}=-(N_{AA}+N')^{-1}$, with $N'$ the Hessian of $-\frac14(d^*A\vert_V,d^*A\vert_V)$ (Feynman gauge, derived, not added).
  - $c$–$A^*$: $\tilde P^{cA^*},\tilde P^{A^*c}=-$(Hessian on $c\oplus Q^{GF}(c^*)$)$^{-1}$, pushed forward.

---

## 8. Interaction vertices ($I=I_3+I_4$, index-free)

$$
I_3=\int\mathrm{dvol}\Big(-\tfrac12(dA,\tfrac12[A\wedge A])+\tfrac12(\lambda,(\not d_A-\not d)\lambda)+([A,c],A^*)+([c,\lambda],\lambda^*)+\tfrac12([c,c],c^*)\Big),
$$
$$
I_4=\int\mathrm{dvol}\Big(-\tfrac14(\tfrac12[A\wedge A],\tfrac12[A\wedge A])\Big).
$$
Component types, using the restrictions of §1:

| term | component types |
|---|---|
| $-\frac12(dA,\frac12[A\wedge A])$ | $(A\vert_V)^3$; $A\vert_V\,(A\vert_{W^\perp})^2$. There is no $(A\vert_{W^\perp})^3$ and no $(A\vert_V)^2A\vert_{W^\perp}$ |
| $I_4$ | $(A\vert_V)^4$; $(A\vert_V)^2(A\vert_{W^\perp})^2$; $(A\vert_{W^\perp})^4$ |
| $\frac12(\lambda,(\not d_A-\not d)\lambda)$ | $\lambda\lambda A\vert_V$ (gauge); $\lambda\lambda A\vert_{W^\perp}$ (Yukawa) |
| $([A,c],A^*)$ | $A\vert_V\,c\,A^*\vert_V$; $A\vert_{W^\perp}\,c\,A^*\vert_{W^\perp}$ |
| $([c,\lambda],\lambda^*)$ | $c\lambda\lambda^*$ |
| $\frac12([c,c],c^*)$ | $cc\,c^*$ |

---

## 9. Momentum-space Feynman rules

**Vertex functions.** $I$ is evaluated on plane-wave test fields $\sum_j r_ju_je^{ip_j\cdot x}$, all momenta incoming:
- $r_j$ is an odd element of an auxiliary nilpotent ring for odd legs [C07 §6.4.1] and $1$ for even legs;
- the polarisations $u_j$ commute.

The vertex function is the coefficient of $r_1r_2\cdots$ (odd legs in leg order) in the multilinear part, with $(2\pi)^4\delta(\sum p_j)$ stripped (§8). Polarisations are $a\in V_{10}^\vee\otimes\mathfrak g$, $u\in\Sigma\otimes\mathfrak g$, $v\in\Sigma^*\otimes\mathfrak g$, $\gamma,w\in\mathfrak g$, $\alpha\in V_{10}\otimes\mathfrak g$. In the table:
- $\iota_pa=a(p)\in\mathfrak g$, which sees only $a\vert_V$;
- $[a_1,a_2]$ contracts the $V_{10}^\vee$ factors with the metric and brackets the $\mathfrak g$ factors;
- $([x\wedge y],[z\wedge w])$ is the Hodge pairing, with $[x\wedge y](e,e')=[x(e),y(e')]-[x(e'),y(e)]$.

| vertex | vertex function |
|---|---|
| $AAA$ | $-\frac i2\big[(\iota_{p_2-p_1}a_3,[a_1,a_2])+(\iota_{p_3-p_2}a_1,[a_2,a_3])+(\iota_{p_1-p_3}a_2,[a_3,a_1])\big]$ |
| $AAAA$ | $-\frac12\big[([a_1\wedge a_2],[a_3\wedge a_4])+([a_1\wedge a_3],[a_2\wedge a_4])+([a_1\wedge a_4],[a_2\wedge a_3])\big]$ |
| $\lambda\lambda A$ | $r_1r_2\,(u_1,\rho([a,u_2]))$ |
| $A\,c\,A^*$ | $r_cr_*\,([a,\gamma],\alpha)$ |
| $c\,\lambda\,\lambda^*$ | $r_cr_\lambda\,([\gamma,u],v)$ |
| $c\,c\,c^*$ | $r_1r_2\,([\gamma_1,\gamma_2],w)$ |

**Expansion of $W$.** [C07 Def 6.6.1, §15.4; CG2 §7.2.7; ledger F8]
$$
W(P(\varepsilon,\infty),I)=\hbar\log\big(e^{\hbar\partial_{P(\varepsilon,\infty)}}e^{I/\hbar}\big),\qquad \partial_P=\tfrac12\sum\partial_{P''}\partial_{P'}\quad[\text{C07 §6.4}],
$$
$$
W=\sum_{\gamma\ \text{connected}}\hbar^{\,b_1(\gamma)+\sum_v(\text{vertex }\hbar\text{-degree})}\ \frac{1}{|\operatorname{Aut}\gamma|}\,w_\gamma .
$$
- Vertices carry $I$ (weight $\hbar^{-1}$ each) and edges carry $\hbar P(\varepsilon,\infty)$. A counterterm vertex $\hbar\,I^{CT}$ has $\hbar$-degree 1.
- Edge rule: two contracted coordinate functionals, brought adjacent by Koszul reordering, become $\hbar$ times the component of $P(\varepsilon,\infty)$ in the same order (§9).
- $|\operatorname{Aut}\gamma|$ counts automorphisms fixing the external legs. Vertex functions are the full multilinear forms of the table, and the sum runs over all component labels.
- One loop (§9, verified in a 0-dim model with mixed parities and odd backgrounds): the background-field part of $W$ at order $\hbar^1$ is
$$
-\frac\hbar2\,\mathrm{str}\log\big(1-P\circ(-1)^{|\cdot|}\circ\partial^2I\big)=\hbar\sum_{n\ge1}\frac1{2n}\,\mathrm{str}\big[(P\circ(-1)^{|\cdot|}\circ\partial^2I)^n\big].
$$
  - $\partial^2I$ is the second derivative of $I$ in the fluctuation, with the right-hand slot differentiated first.
  - $(-1)^{|\cdot|}$ is the parity operator.
  - $\mathrm{str}$ is the supertrace $(-1)^{|\cdot|}$; it gives the fermion-loop and ghost-loop signs.
  - The forms without the parity operator, or with the plain trace, fail.
  - $n=1$ reproduces $\partial_PI$.

**Superficial degree.** With dimensions $A:1$, $\lambda:\frac32$, $c:0$, $A^*:3$, $\lambda^*:\frac52$, $c^*:4$ (all vertices have dimension 4, and the $c$/$A^*$ split is irrelevant for degree-0 functionals):
- a 1PI graph has $\omega=4-\sum_{\text{ext}}\dim$;
- with one insertion of an operator of dimension $\Delta_O$, $\omega=\Delta_O-\sum_{\text{ext}}\dim$.

---

## 10. Open items for downstream use

- $(Q^{GF})^2\neq0$, and [C07 Def 4.3.3(5)] fails too (§6). The [C07 Lemma 6.5.1] identity $P(0,\infty)=Q^{-1}$ on $\operatorname{Im}Q^{GF}$ therefore holds only sector-wise (§7, Hessian inversion). $QP(\varepsilon,\infty)=K_\varepsilon$ holds exactly.
- G1: non-compact $\mathbb R^4$. $K_\infty=0$ is used; no IR regulator is introduced.
- The coupling sits in $(-,-)$, and every propagator carries $(-,-)^{-1}$.
