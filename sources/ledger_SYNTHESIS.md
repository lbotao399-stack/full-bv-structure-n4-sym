# Consolidated notation ledger: one-loop BV quantization of 4d N=4 SYM in Costello's formalism

**Program scope:**
- Euclidean $\mathbb R^4$, with complexified fields.
- Weyl spinors.
- All superconformal ghosts.
- One loop.
- Local $O(\varepsilon)$.
- UV boundary $\varepsilon$ only. The upper Schwinger limit is fixed at $\infty$, and the symbol $T$ is not used by the program.

**Source keys:**
- **[C07]** Costello, arXiv:0706.1533v3. The book was not read.
- **[CG2]** Costello–Gwilliam, *Factorization Algebras in QFT* Vol. 2, authors' draft of 27 Oct 2020. Book page = PDF page − 10. The source prints the BV Laplacian as $\triangle$.
- **[ESW]** Elliott–Safronov–Williams, arXiv:2002.10517v1. Printed page = PDF page.
- **[SW]** Saberi–Williams, arXiv:1910.04120v3. The ESW ledger read v1 §3.1 p.13 for (7)–(8); the labels are the same.

**Markers:**
- **[derived]**: a one-line check made from the ledgers. It is not in the source.
- **[sic]**: the source's own typo.
- **Program column**:
  - "T=∞": the source formula with its upper scale set to $\infty$. [C07] allows $0<\varepsilon<T\le\infty$.
  - "finite-T": the item needs a finite IR scale, so the program does not use it.

## 1. FORMALISM [C07]

| Symbol | Meaning | Source, location | Verbatim formula | Program |
|---|---|---|---|---|
| $M,E,\mathscr E$ | compact manifold; complex super vector bundle; fields | §1, §4.1 | $\mathscr{E}=\Gamma(M,E)$; $\mathscr E\otimes\mathscr E=\Gamma(M^2,E\boxtimes E)$ (completed projective $\otimes$, fn. 1) | compactness is assumed throughout (G1) |
| $\langle\,,\rangle_M,\ \langle\,,\rangle$ | odd symplectic pairing | Def 4.3.1 | $\langle\,,\rangle_M:E\otimes E\to\mathrm{Densities}_M$, odd, antisymmetric, $E\to E^\vee\otimes\mathrm{Densities}_M$ iso; induces $\langle\,,\rangle:\mathscr E\otimes_{\mathbb C}\mathscr E\to\mathbb C$ | used |
| $Q$ | linear BV differential | §1, §2.1, §4.3 | odd, square zero, "skew self adjoint" (§1, §2.1); "odd linear elliptic differential operator … self adjoint", $[Q,Q]=0$ (§4.3) | used; X8 |
| action | free part + interaction | §1, §4.3 | $\frac{1}{2}\langle e,Qe\rangle+S(e)$, $S\in\mathscr O_l(\mathscr E)$ at least cubic; CME $QS+\tfrac12\{S,S\}=0$ | $S$ = interaction (X1) |
| $Q^{GF}$, $L$ | gauge fixing; gauge-fixing subspace | §1.1, Def 4.3.3 | $L=\operatorname{Im}Q^{GF}$. Def 4.3.3: (1) skew self adjoint (§1.1: "self adjoint"); (2) order $\le1$; (3) $[Q^{GF},Q^{GF}]=0$; (4) $H:=[Q,Q^{GF}]$ generalised Laplacian; (5) $\mathscr E=\operatorname{Im}Q\oplus\operatorname{Im}Q^{GF}\oplus\operatorname{Ker}H$ | needed for N=4 (G15) |
| $H$ | generalised Laplacian | §1.1, 4.3.3(4), §1.9 | $H\overset{\text{def}}{=}[Q,Q^{GF}]$; §1.9: $H=[Q^{GF},Q]$; $\sigma(H)\in\Gamma(M,\mathrm{Sym}^2TM\otimes\mathrm{End}\,E)$ positive definite metric $\otimes\,\mathrm{id}_E$; $H^*(\mathscr E,Q)=\operatorname{Ker}H$ | used |
| $\mathscr O(\mathscr E)$ | functionals; families | Def 4.1.1, 4.1.2 | $\mathscr{O}(\mathscr{E})=\prod_{n\geq0}\operatorname{Hom}(\mathscr{E}^{\otimes n},\mathbb{C})_{S_n}$; $\mathscr O(\mathscr E,C^\infty(X))=\prod_{n\ge0}\operatorname{Hom}(\mathscr E^{\otimes n},C^\infty(X))_{S_n}$; $\mathscr O(\mathscr E,\mathbb C[[\hbar]])=\varprojlim\mathscr O(\mathscr E)\otimes\mathbb C[\hbar]/\hbar^n$ | used |
| $\mathscr O_l(\mathscr E)$ | local functionals | Def 4.2.2 | each $f_n$ factors $\mathscr E^{\otimes n}\to\mathrm{Densities}(M)\xrightarrow{\int_M}\mathbb C$, the first map polydifferential; "not an algebra" | used |
| $\{f,g\}$ | odd bracket | Lemma 4.3.2, §10.1 | defined for $f\in\mathscr O_l$, $g\in\mathscr O$; for $\dim M=0$: $\{f,g\}=\Delta_0(fg)-\Delta_0(f)g-(-1)^{\lvert f\rvert}f\Delta_0(g)$ | used; [C07] has no scale-dependent bracket |
| $\partial_e$, $Qf$ | derivations | §6.2 | $f(x+\varepsilon e)=f(x)+\varepsilon\partial_e f$; $f(x+\varepsilon Qx)=f(x)+\varepsilon(Qf)(x)$; $[\partial_e,Q]=\partial_{Qe}$ | here $\varepsilon$ is a formal parameter (X16) |
| $K\star$ | convolution | §6.3, Lemma 6.3.1 | $K\star e=(-1)^{\lvert e\rvert}\sum K'\otimes\langle K'',e\rangle$; $(QK)\star e=[Q,K\star]e$; $K$ (anti)symmetric $\iff K\star$ (skew) self adjoint | used |
| $\partial_\phi$ | order-two contraction | §6.4 | $\partial_{\phi}=\frac{1}{2}\sum\partial_{\phi''}\partial_{\phi'}=\frac{1}{2}\sum(-1)^{\lvert\phi'\rvert\lvert\phi''\rvert}\partial_{\phi'}\partial_{\phi''}$, $\phi\in\mathrm{Sym}^2\mathscr E$; $[\partial_\phi,Q]=\partial_{Q\phi}$; $[\partial_\phi,e^\vee]=-\partial_{\phi\star e}$, $e^\vee(e')=\langle e',e\rangle$ | used; X2, X3 |
| $K_t,\ L_t$ | heat kernel | §6.5 | $K_t\star e=e^{-tH}e$; $L_t=(Q^{GF}\otimes1)K_t$; $(Q\otimes1+1\otimes Q)K_t=0$; $\frac{d}{dt}K_t+(Q\otimes1+1\otimes Q)L_t=0$; $K_t+dtL_t$ closed in $\mathscr E\otimes\mathscr E\otimes\Omega^*(\mathbb R_{>0})$ | $t\in(\varepsilon,\infty)$ |
| $P(\varepsilon,T)$, $P_\varepsilon^T$ | propagator | §1.2, §6.5, §10.2 | $P(\varepsilon,T)=\int_\varepsilon^T L_t\,dt=\int_\varepsilon^T(Q^{GF}\otimes1)K_t\,dt\in\mathscr E\otimes\mathscr E$, $0<\varepsilon<T\le\infty$ | T=∞: only $P(\varepsilon,\infty)$. It lies in $\mathscr E\otimes\mathscr E$ "because $Q^{GF}$ annihilates the 0-eigenspace" (§6.5, compact $M$) |
| $P(0,\infty)$ | full propagator | §1.2, Lemma 6.5.1 | $\int_0^\infty(Q^{GF}\otimes1)K_t\,dt$, distributional; represents $Q^{-1}:\operatorname{Im}Q\to\operatorname{Im}Q^{GF}$ | $\varepsilon\to0$ limit |
| $K_\infty$ | — | proof of 10.3.2 | element of $\operatorname{Ker}H\otimes\operatorname{Ker}H$ representing the identity on $\operatorname{Ker}H$ | T=∞; G1, G4 |
| $\Delta_T$ | BV Laplacian | §1.8, §10.1, §2.1 | $\Delta_T=-\partial_{K_T}$, odd, $\Delta_T^2=0$; $\Delta_0$ ill-defined for $\dim M>0$; $\dim M=0$: $K_0=-\sum_i(x_i\otimes\xi_i+\xi_i\otimes x_i)$, $\Delta_0=-\partial_{K_0}=\sum\partial_{x_i}\partial_{\xi_i}$; §2.1: $\Delta=\sum\partial_{x_i}\partial_{\eta_i}$, $\langle x_i,\eta_i\rangle=1$ | program uses $\Delta_\varepsilon=-\partial_{K_\varepsilon}$ and $\Delta_\infty=-\partial_{K_\infty}$ |
| $\Gamma(\Phi,S)$ | RG flow (the book's $W$) | Def 6.6.1, §1.2 | $\Gamma(\Phi,S)=\hbar\log(\exp(\hbar\partial_\Phi)\exp(S/\hbar))\in\mathscr O(\mathscr E,\mathbb C[[\hbar]])$, $\Phi$ even, $S$ even, at least cubic mod $\hbar$; §1.2 [sic]: $\Gamma(P(\varepsilon,T),S)=\hbar\log(e^{\hbar\partial_{P(\varepsilon,\infty)}}e^{S/\hbar})$ | $\Gamma(P(\varepsilon,\infty),S)$ |
| $\Gamma_{i,k}$; order | $\hbar^i$ coefficient, degree $k$ | §1.4, §7.1, proof of Thm B | $\Gamma(P(\varepsilon,T),S)=\sum_{i,k\geq0}\hbar^i\Gamma_{i,k}(P(\varepsilon,T),S)$; lexicographic on $(i,k)$; $A_{\le(I,K)}=\sum_{(i,k)\le(I,K)}\hbar^iA_{i,k}$ | one loop: $i=1$ |
| $\mathscr A$ | $\varepsilon$-functions | Def 7.1.1 | finite spans of $\int_{U\subset(\varepsilon,1)^n}\frac{F(t_1,\dots,t_n)^{1/2}}{G(t_1,\dots,t_n)^{1/2}}dt_1\dots dt_n$ and $\int_{U\subset(\varepsilon,1)^{n-1}}\frac{F(t_1,\dots,t_n=\varepsilon)^{1/2}}{G(t_1,\dots,t_n=\varepsilon)^{1/2}}dt_1\cdots dt_{n-1}$, $F,G\in\mathbb Z_{\ge0}[t_1,\dots,t_n]\setminus\{0\}$, $U$ cut out by $t_i^l>t_j$, $l\in\mathbb Z$ | used |
| $\mathscr A_{\ge0},\ \mathscr A_{<0}$ | renormalisation scheme | Def 1.4.1, 8.0.2 | $\mathscr A_{\ge0}$: the $\varepsilon\to0$ limit exists; $\mathscr A=\mathscr A_{<0}\oplus\mathscr A_{\ge0}$; "singular part" = projection onto $\mathscr A_{<0}$ | used |
| $\mathscr A_{\le0},\ \mathscr A_{>0}$ | — | §12 | $\mathscr A_{\le0}=\mathscr A_{<0}\oplus\mathbb C$; $\mathscr A_{>0}$: the limit exists and $=0$; $\mathscr A=\mathscr A_{\le0}\oplus\mathscr A_{>0}$ | used (obstruction) |
| $S^{CT},\ S^R$ | counterterms | §8 | $S^R(\hbar,\varepsilon)=S(\hbar)-S^{CT}(\hbar,\varepsilon)=S-\sum_{i>0,k\geq0}\hbar^iS^{CT}_{i,k}(\varepsilon)$ | one loop: $S^{CT}_{1,k}$ |
| $\Gamma^R$ | renormalised flow | §10.3 | $\Gamma^R(P(0,T),S)\overset{def}{=}\lim_{\varepsilon\to0}\Gamma(P(\varepsilon,T),S^R)$ | T=∞: $\Gamma^R(P(0,\infty),S)$ |
| $Z,\ Z^R,\ d\mu$ | partition function; Gaussian measure | §1, §1.2, §1.5, §6.6 | $Z(S,\hbar,a)=\int_{x\in L}\exp\left(\frac{1}{2\hbar}\langle x,Qx\rangle+\frac1\hbar S(x+a)\right)dx$ $=\lim_{\varepsilon\to0}(e^{\hbar\partial_{P(\varepsilon,\infty)}}e^{S/\hbar})(a)$ up to $\det Q$; $Z^R=\exp(\Gamma^R(P(0,\infty),S)/\hbar)\in\mathscr O(\mathscr E)((\hbar))$; $\int_{x\in(\operatorname{Im}Q^{GF})_{\mathbb R}}\exp(\tfrac12\langle x,Qx\rangle/\hbar)d\mu=1$, with $\langle x,Qx\rangle$ negative definite on the real slice | T=∞ is native here |
| $\delta$ | odd parameter | (12.0.1) | "odd parameter" | X17, G6 |
| $O(S),\ O_{(I,K)}$ | obstruction | Prop 12.0.3 | see F7 | T=∞ (G3, G4) |
| $w_\gamma,\ \operatorname{Aut}\gamma$ | graph weights | §15.4 | $\Gamma_{i,k}(P(\varepsilon,T),S)(\alpha)=\sum_\gamma\frac{1}{\operatorname{Aut}\gamma}w_\gamma(\alpha)$; the source does not give $w_\gamma$ | G16 |
| $S^{eff}(T)$; $\mathrm{GF},\mathrm{BV}_{\mathscr E}$ | effective actions; simplicial sets | Def 9.0.7; §11 | — | finite-T / homotopy; not used |

### Key statements [C07] in UV-boundary form

**F1. Lemma 6.5.1.** $\alpha\mapsto\int_0^\infty Q^{GF}e^{-tH}\alpha\,dt$ inverts $Q:\operatorname{Im}Q^{GF}\to\operatorname{Im}Q$. In the proof, $\int_0^\infty Q^{GF}e^{-tH}Q\alpha\,dt=\alpha-\pi\alpha$, with $\pi:\mathscr E\to\operatorname{Ker}H$.

**F2. Lemma 6.6.2 ($M$ a point).** $\Gamma(P(0,\infty),S)(a)=\hbar\log\int_{x\in(\operatorname{Im}Q^{GF})_{\mathbb R}}\exp\left(\tfrac12\langle x,Qx\rangle/\hbar+S(x+a)/\hbar\right)d\mu$. For $\dim M>0$ the source takes this as the definition (§7.1).

**F3. Theorem A(1).** $\Gamma_{i,k}(P(\varepsilon,T),S(\varepsilon,\hbar))(e)\simeq\sum_{r\ge0}f_r(\varepsilon)\Phi_r(T,e)$, with $f_r\in\mathscr A$ and $\Phi_r\in\mathscr O(\mathscr E,C^\infty((0,\infty)))$. The meaning of $\simeq$ is given in §15.3. Parts (2) and (3), the $T\to0$ expansion into $\Psi_{q,r}\in\mathscr O_l$, are finite-T.

**F4. Theorem B.**
- Statement: let $S\in\mathscr O_l(\mathscr E,\mathbb C[[\hbar]])$ be even and at least cubic mod $\hbar$. Then there is a unique $S^{CT}(\varepsilon,\hbar)=\sum_{i\ge1,k\ge0}\hbar^iS^{CT}_{i,k}(s)$ [sic] with:
  - $S^{CT}_{i,k}(\varepsilon)\in\mathscr O_l(\mathscr E,\mathscr A_{<0})$, homogeneous of degree $k$;
  - $\lim_{\varepsilon\to0}\Gamma(P(\varepsilon,T),S-S^{CT})$ exists.
- Construction: $S^{CT}_{I,K}=$ singular part of $\Gamma_{I,K}(P(\varepsilon,T),S-S^{CT}_{<(I,K)})$.
- §1.5 order: $S^{CT}_{1,1}$, then $S^{CT}_{1,2}$ from $S-\hbar S^{CT}_{1,1}$, and so on. Lemma 8.0.5: $S^{CT}_{0,k}=0$. §1.5: the counterterms are independent of $T$.
- After the statement: $\Gamma^R(P(0,\infty),S)$ "is a renormalised functional integral".
- The proofs of $T$-independence and locality (Lemmas 8.0.3–8.0.4) use finite $T,T'$ and $T\to0$. See G2.

**F5. Lemma 10.2.1 and Lemma 10.2.2.**
- Lemma 10.2.1: $QP_\varepsilon^T=-K_T+K_\varepsilon$. In the proof, $QP(\varepsilon,T)$ is the kernel of $\alpha\mapsto Q\int_\varepsilon^TQ^{GF}e^{-tH}\alpha+\int_\varepsilon^TQ^{GF}e^{-tH}Q\alpha=-e^{-TH}\alpha+e^{-\varepsilon H}\alpha$.
  - T=∞ substitution: $QP(\varepsilon,\infty)=-K_\infty+K_\varepsilon$.
- Lemma 10.2.2: $(Q+\hbar\Delta_{T_1})\exp(S/\hbar)=0\iff(Q+\hbar\Delta_{T_2})\exp(\Gamma(P(T_1,T_2),S)/\hbar)=0$.
- Identities used in the proof:
  - $[\partial_{P(T_1,T_2)},Q]=\partial_{QP(T_1,T_2)}=\partial_{K_{T_1}}-\partial_{K_{T_2}}=\Delta_{T_2}-\Delta_{T_1}$;
  - $(Q+\hbar\Delta_{T_2})\exp(\hbar\partial_{P(T_1,T_2)})\exp(S/\hbar)=\exp(\hbar\partial_{P(T_1,T_2)})(Q+\hbar\Delta_{T_1})\exp(S/\hbar)$.

**F6. Def 10.3.1 and Lemma 10.3.2.**
- Def 10.3.1: $(Q+\hbar\Delta_T)\exp\left(\hbar^{-1}\Gamma^R(P(0,T),S)\right)=0$ "for some (or equivalently, all) $T>0$".
- Lemma 10.3.2: $\Gamma^R(P(0,\infty),S)|_{\operatorname{Ker}H}$ satisfies the QME on $\operatorname{Ker}H=H^*(\mathscr E,Q)$.

**F7. Obstruction.**
- (12.0.1): $\Gamma\big(P(\varepsilon,T)+\delta\Delta_\varepsilon,S^R+\delta QS^R+\delta d_{DR}S^R\big)=\hbar\log\big((1+\delta Q+\delta d_{DR}+\delta\Delta_T)\exp(\hbar\partial_{P(\varepsilon,T)})\exp(S^R/\hbar)\big)$. It follows from $[\partial_{P(\varepsilon,T)},Q+d_{DR}]=\Delta_T-\Delta_\varepsilon$.
- Prop 12.0.3:
  - There is a unique odd $O(S)\in\mathscr O_l(\mathscr E,\mathscr A_{\le0}\otimes\Omega^*(\Delta^n)\otimes\mathbb C[[\hbar]])$ with $\lim_{\varepsilon\to0}\frac{d}{d\delta}\langle\Gamma(P(\varepsilon,T)+\delta\Delta_\varepsilon,S^R+\delta QS^R+\delta d_{DR}S^R+\delta O(S))\rangle=0$.
  - $O_{(I,K)}(S_{<(I,K)}+\hbar^IS_{(I,K)})=\hbar^IQS_{(I,K)}+\hbar^Id_{DR}S_{(I,K)}+O_{(I,K)}(S_{<(I,K)})$.
  - RQME $\iff O(S)=0$.
- Lemma 12.0.4: if $O_{<(I,K)}(S)=0$, then $O_{(I,K)}(S)$ is $\varepsilon$-independent, lies in $\mathscr O_l(\mathscr E,\Omega^*(\Delta^n))$, and $(Q+d_{DR})O_{(I,K)}=0$. The lift is $S+\hbar^IS'$ with $(Q+d_{DR})S'=-O_{(I,K)}(S)$.
- Lemma 14.4.2(3): $O_{(I,K)}(S+\hbar^IS')=O_{(I,K)}(S)+(Q+d_{DR})S'$ [sic: the $\hbar^I$ is missing].

**F8. Graphs (§15.4–15.5).**
- Each vertex degree equals its $\hbar$ power. Degree-0 vertices are at least trivalent. $\sum$(vertex degrees) $+\,b_1=i$.
  - [derived] For $i=1$: either $b_1=1$ with all vertex degrees $0$, or $b_1=0$ with exactly one degree-1 vertex.
- General class: $f_\gamma(t_e,\alpha)=\int_{M^{V(\gamma)}}\otimes S_v\big(\alpha\otimes\otimes_{e\in E(\gamma)}K_{t_e}\big)$ and $F_\gamma(\varepsilon,T,\alpha)=\int_{t_e\in(\varepsilon,T)^{E(\gamma)}}f_\gamma\prod dt_e$.
- Remark: the variant propagator $P(\varepsilon,T)+\delta K_\varepsilon$ has one special edge carrying $K_\varepsilon$.

**F9. Yang–Mills (§5.2).**
- First-order formulation [MZ96, Wit04, Cos06] on a compact oriented 4-manifold, with $d_AF(A)_+=0$.
- The dga $\mathcal B$ by degree:
  - degree 0: $\Omega^0(\mathrm{End}V)$;
  - degree 1: $\Omega^1\oplus\Omega^2_+$;
  - degree 2: $\Omega^2_+\oplus\Omega^3$;
  - degree 3: $\Omega^4$.
- Arrows of $Q$:
  - $d_A$ and $-[F(A)_+,\,]$ out of $\Omega^0$;
  - $d_{A+}$ and $-[F(A)_+,\,]$ out of $\Omega^1$;
  - $\mathrm{Id}$ and $-d_A$ out of $\Omega^2_+$;
  - $-[F(A)_+,\,]$ and $-d_A$ into $\Omega^4$.
- $\mathscr E=\mathcal B[-1]$ [sic], $S(e)=\mathrm{Tr}\,a^3$ [sic]. $\mathscr E^{-1}$ = ghosts, $\mathscr E^0$ = fields, $\mathscr E^1$ = antifields, $\mathscr E^2$ = antighosts.
- $Q^{GF}$ is the Hermitian adjoint of $Q$.
- The author does not see how to fit the standard (second-order) formulation.

**F10. Flat $\mathbb R^n$ (§15.2).**
- $K_t=Ct^{-n/2}e^{-\lVert x-y\rVert^2/4t}d(x_1-y_1)\wedge\cdots\wedge d(x_n-y_n)\otimes I_{\mathfrak g}$.
- Lemma 15.2.2: $P(\varepsilon,T)=\left(\int_{\lVert x-y\rVert^2/4T}^{\lVert x-y\rVert^2/4\varepsilon}u^{n/2-1}e^{-u}du\right)\pi^*\mathrm{Vol}_{S^{n-1}}\otimes I_{\mathfrak g}$, up to an overall constant.
- Prop 15.2.5 fixes $0<T<\infty$.

**Finite-T only, not used:**
- the semigroup $\Gamma(P(T_2,T_3),\Gamma(P(T_1,T_2),f))=\Gamma(P(T_1,T_3),f)$;
- Lemma 9.0.6, Def 9.0.7, Theorem C;
- Theorem A(2)–(3).

### Sign conventions [C07]

- **CS1. Weight.**
  - $e^{+S/\hbar}$, with integrand $\exp(\frac1{2\hbar}\langle x,Qx\rangle+\frac1\hbar S(x+a))$.
  - The real slice is chosen with $\langle x,Qx\rangle$ negative definite (§6.6).
  - §5.1 example: $\exp\big(\hbar^{-1}\int_M\tfrac12\phi\Delta\phi+\tfrac12m^2\phi^2+\phi^4\big)$, with $H(\phi e_i)=(\Delta\phi+m^2\phi)e_i$ and $Q^{GF}(\phi e_1)=\phi e_0$.
- **CS2. QME.**
  - $(Q+\hbar\Delta)e^{S/\hbar}=0\iff QS+\tfrac12\{S,S\}+\hbar\Delta S=0$, with $\Delta_T=-\partial_{K_T}$.
  - §1.11: $\{\Phi_{CS},\Phi_{CS}\}+\hbar\Delta\Phi_{CS}=0$ [sic].
- **CS3. Bracket.** $\{f,g\}=\Delta_0(fg)-\Delta_0(f)g-(-1)^{\lvert f\rvert}f\Delta_0(g)$.
  - [derived] This gives $\{g,f\}=(-1)^{\lvert f\rvert\lvert g\rvert}\{f,g\}$.
- **CS4. Counterterm.** $S^R=S-S^{CT}$, with $S^{CT}$ starting at $\hbar^1$.
- **CS5. Kernel identities.**
  - $QP(\varepsilon,T)=-K_T+K_\varepsilon$;
  - $[\partial_P,Q]=\partial_{QP}$;
  - $[\partial_{P(\varepsilon,T)},Q+d_{DR}]=\Delta_T-\Delta_\varepsilon$.
- **CS6. Normalisation and signs.**
  - The $\tfrac12$ sits inside $\partial_\phi$.
  - Convolution carries the sign $(-1)^{\lvert e\rvert}$.
  - $\int\exp(\frac12\langle x,Qx\rangle/\hbar)d\mu=1$.
- **CS7. Adjointness.** See X8.
- **CS8. Grading.**
  - Mainly $\mathbb Z/2$, with $[1]$ the parity shift (§1.1).
  - §2: $\mathfrak g[1]$ has $\mathfrak g$ in degree $-1$, and $E=\mathfrak g[1]\oplus V\oplus V^\vee[-1]\oplus\mathfrak g^\vee[-2]$.
  - $S$ is even; $\Delta_T$ and $O(S)$ are odd; families enter as $+d_{DR}$.

## 2. OBSERVABLES / SYMMETRIES [CG2]

| Symbol | Meaning | Location | Verbatim formula | Program |
|---|---|---|---|---|
| $\mathscr E,\langle-,-\rangle,Q$ | free BV fields, pairing, linear differential | Def 7.2.1.1 | "a graded antisymmetric map of vector bundles $\langle-,-\rangle_{loc}:E\otimes E\to \mathrm{Dens}(M)$ of cohomological degree $-1$"; "$Q:\mathscr E\to\mathscr E$ of cohomological degree 1 that is skew self adjoint" | same objects as [C07] |
| $Q^{GF}$, $D$ | gauge fixing; generalised Laplacian | §7.2.1 | "square-zero differential operator of cohomological degree $-1$ … self adjoint … $D=[Q,Q^{GF}]$ is a generalized Laplacian" | [C07] calls this $H$ (X11) |
| $K_t$ | heat kernel | §7.2.3 | $\frac{d}{dt}K_t+(D\otimes 1)K_t=0$; $K_0$ is the kernel of the identity; on non-compact $M$ "there is more than one heat kernel" | G1 |
| $\Phi$; $\Phi<\Psi$ | parametrix (the scale); its order | Def 7.2.4.2, §7.2.4 | symmetric, degree 1, proper support, $([Q,Q^{GF}]\otimes 1)\Phi-K_{\mathrm{Id}}$ smooth; $\Phi<\Psi$ "if the support of $\Phi$ is contained in the support of $\Psi$" | finite scale; not used |
| $P(\Phi)$, $K_\Phi$ | propagator; smooth kernel | §7.2.5 | $P(\Phi)=\tfrac12 Q^{GF}\Phi$; $K_\Phi=K_{\mathrm{Id}}-QP(\Phi)$; $K_\Phi-K_\Psi=QP(\Psi)-QP(\Phi)$; for $\Phi=\int_0^L K_t\,dt$: $P(\Phi)=P(0,L)=\int_0^L(Q^{GF}\otimes1)K_t\,dt$ and $K_\Phi=K_L$ | X3 |
| $P^L_\epsilon$ | cut-off propagator (βγ only) | §14.2.2 | $P^L_\epsilon=\int_{t=\epsilon}^L\partial_z\frac{1}{4\pi t}e^{-\lvert z-z'\rvert^2/4t}$ | not used; X4 |
| $\partial_\alpha$ | order-2 contraction | §7.2.7 | "the unique order 2 differential operator that is given by pairing with $\alpha$ on $\mathrm{Sym}^2\mathscr E^\vee$ and that is zero on $\mathrm{Sym}^{\le1}\mathscr E^\vee$" | X2, X3 |
| $W(\alpha,-)$; $W_0$ | RG flow; tree level | §7.2.7 | $F\mapsto\hbar\log\left(e^{\hbar\partial_\alpha}e^{F/\hbar}\right)$ on $\mathscr O^+(\mathscr E)[[\hbar]]$; $W_0(P(\Phi),I)=W((P(\Phi),I)\bmod\hbar$ (sic) | same form as [C07] $\Gamma$ (X1) |
| $\triangle_\Phi$, $\{-,-\}_\Phi$ | scale-$\Phi$ Laplacian and bracket | §7.2.8 | $\triangle_\Phi=\partial_{K_\Phi}$, $(\triangle_\Phi)^2=0$; $\{I,J\}_\Phi=\triangle_\Phi(IJ)-(\triangle_\Phi I)J-(-1)^{\lvert I\rvert}I\triangle_\Phi J$ | X2, X5 |
| $\mathscr O_{P,sm},\ \mathscr O^+_{P,sm},\ \mathscr O_{loc}$ | functional classes | §7.2.6, Def 7.1.0.2 | proper support; smooth first derivative; "+" = at least cubic mod $\hbar$; local = $\int_M(D_1\phi)\cdots(D_n\phi)\,d\mu$; all modulo constants | G9 |
| $I[\Phi]$, $I_{i,k}[\Phi]$ | effective interaction | Def 7.2.9.1 | (i) $I[\Phi]\in\mathscr O^+_{P,sm}(\mathscr E_c)[[\hbar]]$; (ii) $I[\Phi]=W(P(\Phi)-P(\Psi),I[\Psi])$; (iii) $(Q+\hbar\triangle_\Phi)e^{I[\Phi]/\hbar}=0$; (iv) locality via $\mathrm{Supp}(\Phi)^{3i+k}$ | finite scale |
| $\mathrm{Obs}^q_\Phi(M)$, $\widehat Q_\Phi$ | BD algebra at scale $\Phi$ | Def 8.2.1.1, §8.3.1 | $(\mathscr O(\mathscr E)[[\hbar]],\,Q+\hbar\triangle_\Phi+\{I[\Phi],-\}_\Phi)$; $\widehat Q_\Phi=Q+\{I[\Phi],-\}_\Phi+\hbar\triangle_\Phi$ | G5 |
| $W^\Phi_\Psi$; $\delta$ | linearised flow; square-zero parameter | Def 8.3.1.1 | $I[\Phi]+\delta W^\Phi_\Psi(f)=W(P(\Phi)-P(\Psi),I[\Psi]+\delta f)$; $\delta$ "square-zero parameter of cohomological degree $-i$" for $f$ of degree $i$ | X17, G6 |
| $O[\Phi]$, $\mathrm{Obs}^q(M)$ | global quantum observables | Def 8.3.1.3 | $O[\Phi]\in\mathrm{Obs}^q_\Phi(M)=\mathscr O(\mathscr E)[[\hbar]]$, $W^\Phi_\Psi O[\Psi]=O[\Phi]$; $\widehat Q(O)[\Phi]=QO[\Phi]+\{I[\Phi],O[\Phi]\}_\Phi+\hbar\triangle_\Phi O[\Phi]$ | G5 |
| support; $\mathrm{Obs}^q(U)$ | — | Def 8.4.0.1 | "for each $(i,k)$, there exists a compact subset $K\subset U^k$ and a parametrix $\Phi$, such that $\mathrm{Supp}\,O_{i,k}[\Psi]\subset K$ for all parametrices $\Psi\le\Phi$" | — |
| $\mathrm{Obs}^q(p)$ | local operators | Def 10.1.0.1 | "the limit of $\mathrm{Obs}^q(V)$ as $V$ runs over the open sets containing $p$" | G5 |
| $\mathcal M$, $\mathcal M[1]$ | elliptic $L_\infty$ algebra; fields | Def 4.2.0.4, §4.4 | invariant pairing of degree $-3$; "the graded vector bundle $L[1]$ on $M$ has a $-1$ symplectic form" | — |
| $C^*(\mathfrak g),\ C^*_{red},\ \mathrm d_{\mathfrak g}$ | Chevalley–Eilenberg cochains | §11.2.2, §12.2.4 | "$C^*(\mathfrak g)=\mathscr O(B\mathfrak g)$"; $\mathscr O(\mathcal M[1])\otimes\mathscr O(\mathfrak g[1])$ | superconformal ghosts (G7) |
| $\mathcal L$, $\mathrm d_{\mathcal L}$, $\mathcal L\ltimes\mathcal M$ | local symmetry algebra | §12.2, Def 12.2.2.1 | "$0\to\mathcal M\to\mathcal L\ltimes\mathcal M\to\mathcal L\to0$ consists of maps of $L_\infty$ algebras" | — |
| $\mathrm{Act}(\mathfrak g,\mathcal M)$ | global action complex | Def 12.2.4.1(iii) | $C^*_{red}(\mathfrak g)\otimes\mathscr O_{loc}(B\mathcal M)[-1]$; $\mathrm d_{\mathfrak g}S^{\mathfrak g}+\mathrm d_{\mathcal M}S^{\mathfrak g}+\tfrac12\{S^{\mathfrak g},S^{\mathfrak g}\}=0$ | G7 |
| $\mathrm{Act},\ \mathrm{InnerAct}(\mathcal L,\mathcal M)$ | local action complexes | Def 12.2.2.4 | InnerAct = kernel of $C^*_{red,loc}(\mathcal L\oplus\mathcal M)[-1]\to C^*_{red,loc}(\mathcal M)[-1]$; Act = InnerAct modulo $C^*_{red,loc}(\mathcal L)$ | — |
| $\mathrm{Act}^q_{P,sm}$, $\mathrm{InnerAct}^q_{P,sm}$ | quantum versions | §13.2.2–3 | $\mathscr O^+_{P,sm}(\mathcal L[1]\oplus\mathcal M[1])/\mathscr O_{P,sm}(\mathcal L[1])[[\hbar]]$; $\mathscr O^+_{P,sm}(\mathcal L[1]\oplus\mathcal M[1])[[\hbar]]$ | G9 |
| $I^{\mathcal L}[\Phi]$ | equivariant effective interaction | Def 13.2.2.1 | (ii) $\mathrm d_{\mathcal L}I^{\mathcal L}[\Phi]+QI^{\mathcal L}[\Phi]+\tfrac12\{I^{\mathcal L}[\Phi],I^{\mathcal L}[\Phi]\}_\Phi+\hbar\triangle_\Phi I^{\mathcal L}[\Phi]=0$ | G19 |
| $S^{eq},S^{inner},S^{\mathfrak g},S^{\mathcal L},S^{tot}$ | equivariant actions | §11.2.2, §11.3.1, 12.2.4.1, 12.2.2.5 | $S^{tot}=S^{\mathcal L}+S$ | — |
| $\alpha$ | anomaly | Prop 13.3.2.1 | $\alpha=\mathrm d_{\mathcal L}\tilde I^{\mathcal L}[\Phi]+Q\tilde I^{\mathcal L}[\Phi]+\tfrac12\{\tilde I^{\mathcal L}[\Phi],\tilde I^{\mathcal L}[\Phi]\}_\Phi+\hbar\triangle_\Phi\tilde I^{\mathcal L}[\Phi]$, $\Phi$-independent, $\mathrm d_{\mathcal L}\alpha=0$, $\alpha\in C^*_{red,loc}(\mathcal L)$ of degree 1 | G19 |
| $\mathrm{Obs}^q_{\mathcal L}(U)$ | equivariant observables | Prop 13.4.0.1 | $W(P(\Phi)-P(\Psi),I^{\mathcal L}[\Psi]+\epsilon O[\Psi])=I^{\mathcal L}[\Phi]+\epsilon O[\Phi]$, $\epsilon$ of degree $-k$; $(\mathrm dO)[\Phi]=\mathrm d_{\mathcal L}O[\Phi]+QO[\Phi]+\{I^{\mathcal L}[\Phi],O[\Phi]\}_\Phi+\hbar\triangle_\Phi O[\Phi]$ | X16, G6 |
| $U^{BD}(\mathfrak g)$, $\mathbb U^{BD}_\alpha\mathcal L$ | enveloping BD algebras | Def 11.3.1.2, §13.1.2 | $U^{BD}(\mathfrak g)=C_*(\mathfrak g)[[\hbar]]=\mathrm{Sym}^*(\mathfrak g[-1])[[\hbar]]$, differential $d_{\mathfrak g}+\hbar d_{CE}$ | not used |
| $I^{CT}[\epsilon]$, $\mathrm{Sing}_\epsilon$ | one-loop counterterm | §13.7.2 | see below | closest UV-only formula (G4) |
| $\mathsf c_{ij}$ | $so(n)$ ghosts | §12.3 Example | $\mathrm d\Psi_i=\mathsf c_{ij}\Psi_j$, $\mathrm d\mathsf c_{ij}=\mathsf c_{ij}\mathsf c_{jk}$ (sic); degree 1 | only ghost coordinates named in any source |
| $T^{(n)}(\mathscr E,I)$ | simplicial set of theories mod $\hbar^{n+1}$ | §7.4, §7.5 | not a scale | X11 |

### Key statements [CG2]

- **Background-field route (Lemma 10.1.1.3 and its proof).**
  - "$\mathrm{Obs}^{cl}(0)$ are canonically identified with translation-invariant local functionals of both $\mathscr E$ and $\Omega^n(\mathbb R^n)$ that depend linearly on the background field."
  - Proof: $\tilde S=S+\int\omega O(x)$, $I^{b.g.}[\Psi]=I[\Psi]+I'[\Psi]$, $O[\Psi]=\frac{\partial}{\partial\delta_0}I'[\Psi]$; "Here we do not ask that the quantum master equation hold."
- **Family over $\mathbb C[\delta]/\delta^2$ (after Def 8.3.1.3).** $\{I[\Phi]+\delta O[\Phi]\}$ satisfies the axioms of a family of theories over $\mathbb C[\delta]/\delta^2$ except locality.
- **Differential in components (Lemma 8.4.0.2).** $\widehat QO_{i,k}[\Phi]=QO_{i,k}[\Phi]+\sum_{a+b=i,\ r+s=k+2}\{I_{a,r}[\Phi],O_{b,s}[\Phi]\}_\Phi+\Delta_\Phi O_{i-1,k+2}[\Phi]$.
- **Mod $\hbar$ (Lemma 8.6.1.2).** $O[\Phi]=\lim_{\Gamma\to0}W^\Phi_\Gamma(O)\bmod\hbar$, a sum over trees with one special vertex labelled by $O$ and every other vertex labelled by $I_0$.
- **Products (§8.5, 10.3.2.1).** $m(O,O')[\Phi]=\lim W^\Phi_\Psi(O[\Psi]*O'[\Psi])$.
- **Global to local symmetry.**
  - Lemma 12.2.4.2: $\mathfrak g$-actions are canonically homotopy equivalent to $\Omega^*_X\otimes\mathfrak g$-actions.
  - Lemma 12.2.4.3: if there is no gauge symmetry, the action factors through $\Omega^{\le1}_X\otimes\mathfrak g$.
  - The quantum analogue is announced in §13.1.3 but was not found.
- **Equivariant QMEs.**
  - §11.3.1: $\mathrm d_{\mathfrak g}S^{eq}+\hbar\triangle S^{eq}+\tfrac12\{S^{eq},S^{eq}\}=0$.
  - (12.2.2.2): $(\mathrm d_{\mathcal L}+\mathrm d_{\mathcal M})S^{\mathcal L}+\tfrac12\{S^{\mathcal L},S^{\mathcal L}\}=0$, i.e. $\mathrm d_{\mathcal L}S^{tot}+\tfrac12\{S^{tot},S^{tot}\}=0$.
  - Def 13.2.1.1: $\mathrm d_{\mathcal L}I+QI+\{I,I\}_\Phi+\hbar\triangle_\Phi I=0$ (no ½).
- **Obstruction complexes (§13.3.1).**
  - $\mathrm{Act}(\mathcal L,\mathcal M)$ when the non-equivariant quantization is fixed;
  - $C^*_{red,loc}(\mathcal L\ltimes\mathcal M\mid\mathcal L)$ for the full equivariant problem;
  - $C^*_{red,loc}(\mathcal L\ltimes\mathcal M)$ or $\mathrm{InnerAct}(\mathcal L,\mathcal M)$ for inner actions.
- **§13.7.2 (the only explicit ε-counterterm in the book).**
  - $I^{(1)}_1[\epsilon,L](\rho_1)=\int_{t=\epsilon}^L\int_{x\in X}\mathrm{STr}_E\,\rho_1(x)\,(\not\partial_xK_t(x,y))_{y=x}\,\mathrm dt$.
  - $I^{CT}[\epsilon](\rho_1)=\mathrm{Sing}_\epsilon I^{(1)}_1[\epsilon,L](\rho_1)$.
  - $I^{(1)}_1[L]=\lim_{\epsilon\to0}(I^{(1)}_1[\epsilon,L]-I^{CT}[\epsilon])$.
  - $\mathrm d_{dR}I^{(1)}_1[L]+\triangle_LI^{(1)}_0=0$, and $\mathrm d_{dR}I^{(1)}[\epsilon,L]=\triangle_\epsilon I^{(1)}_0-\triangle_LI^{(1)}_0$.
  - $\alpha^{(1)}(\rho_0)=\lim_{\epsilon\to0}\left(\triangle_\epsilon I^{(1)}_0(\rho_0)-I^{CT}[\epsilon](\mathrm d\rho_0)\right)$.
  - The counterterm term vanishes as "a purely singular function of $\epsilon$ (… a linear combination of $\log\epsilon$ and $\epsilon^{-k}$ for $k>0$)".
  - Typos on the same page: "$I^{(1)}[L]=I^{(1)}_0[L]+\hbar I^{(1)}_0[L]$" (sic) and "$\triangle_L=K_L$" (sic).

### Sign conventions [CG2]

- **CG1. Weight.**
  - The formal definitions use $e^{+I[\Phi]/\hbar}$.
  - The §8.3.1 Remark uses $\langle O\rangle=\frac1{Z_S}\int O\,e^{-S/\hbar}$ and $e^{(-S+\hbar\delta O)/\hbar}$, with $O$ at weight $\hbar^0$. Def 8.3.1.1 instead has $\delta O$ inside $e^{(I+\delta O)/\hbar}$. The source does not reconcile these.
  - Finite dimensions: volume form $e^Sdx$, $\mathrm{div}_S=\triangle+\{S,-\}$.
- **CG2. QME.** $QI+\hbar\triangle_\Phi I+\tfrac12\{I,I\}_\Phi=0$.
- **CG3. ½ in the equivariant QME.** Present in §11.3.1, 13.2.2.1(ii) and 13.3.2.1; absent in 13.2.1.1.
- **CG4. Counterterm.** $\lim_{\epsilon\to0}(I[\epsilon,L]-I^{CT}[\epsilon])$, with $I^{CT}=\mathrm{Sing}_\epsilon$ (pure $\log\epsilon$ and $\epsilon^{-k}$). Shown only for a one-loop interaction term.
- **CG5. Kernels.**
  - $K_\Phi=K_{\mathrm{Id}}-QP(\Phi)$; $Q(P(\Phi)-P(\Psi))=K_\Psi-K_\Phi$; $Q(P(0,L)-P(0,\epsilon))=K_\epsilon-K_L$.
  - §14.2.2: $\bar\partial P(\epsilon,L)=K_L-K_\epsilon$ (opposite sign).
- **CG6. Bracket.**
  - §7.2.8 form as in the table above.
  - §6.3: $\mathrm{div}(ab)=(\mathrm{div}\,a)b+(-1)^{\lvert a\rvert}a(\mathrm{div}\,b)+(-1)^{\lvert a\rvert}\{a,b\}$.
- **CG7. Degrees.**
  - Field pairing has degree $-1$; the $\mathcal M$ pairing has degree $-3$; fields are $\mathcal M[1]$; $C^*(\mathcal L)\cong\mathscr O(\mathcal L[1])$.
  - $\Phi$ has degree 1, $P(\Phi)$ degree 0; $\delta$ (or $\epsilon$) has degree $-i$ for $O$ of degree $i$.
  - $\hbar$ has weight $-1$ (§6.4); $\mathsf c_{ij}$ has degree 1.
  - "$\mathbb Z/2\times\mathbb Z$ graded … Both gradings contribute to signs".
- **CG8. Enveloping BD algebra.**
  - Def 11.3.1.2 uses $\mathrm{Sym}(\mathfrak g[1])$, while §13.1.2 uses $\mathrm{Sym}^*(\mathfrak g[-1])$.
  - The differentials are written as $d(xy)=d_{\mathfrak g}(x)y\pm x\,d_{\mathfrak g}(y)+\hbar[x,y]$ and as $d_{\mathfrak g}+\hbar d_{CE}$.

## 3. N=4 FIELDS / SUSY [ESW], with [SW §4.1]

| Symbol | Meaning | Location | Verbatim formula | Program |
|---|---|---|---|---|
| $V_{\mathbb R},V$ | spacetime and its complexification | ESW §2.2 p.36 | real, $\dim n$, oriented, nondegenerate symmetric form; $V=V_{\mathbb R}\otimes_{\mathbb R}\mathbb C$ | Euclidean $\mathbb R^4$ |
| $\mathrm{Cl}(V)$ | Clifford algebra | §2.1 p.35 | $v_1v_2+v_2v_1=2(v_1,v_2)$ | G11 |
| $M=M^+\oplus M^-,\ \rho(v)$, $(-,-)$ | graded Clifford module; action; pairing | §2.1 p.35 | $(-,-)\colon M^+\otimes M^-\to\mathbb C$, $(\rho(v)Q_1,Q_2)=(Q_1,\rho(v)Q_2)$; "$M^+=\Sigma$ and $M^-=\Sigma^*$" | — |
| $\Gamma$ | Γ-pairing | (28) p.35 | $\Gamma\colon \mathrm{Sym}^2(\Sigma)\to V,\ \mathrm{Sym}^2(\Sigma^*)\to V$ by $(\rho(v)Q_1,Q_2)=(v,\Gamma(Q_1,Q_2))$; "typically unique up to a scale" (p.36) | G11; X12 |
| $\wedge^2V\to\mathfrak{so}(V)$; $q$ | two-forms; antisymmetrisation | §2.1, (29) p.35 | $\omega\mapsto(w\mapsto-2\iota_{(w,-)}\omega)$; $q(v_1\wedge v_2)=v_1v_2-(v_1,v_2)$ | — |
| $\not d,\ \not d_A$ | Dirac operators | §2.1 p.35, §3.1 p.44 | $\not d\colon\Gamma(N,\Sigma)\to\Gamma(N,\Sigma^*)$; $\not d_A$ "obtained from $\Gamma$" | — |
| $S;\ S_\pm$ | spin ($n$ odd) and semi-spin ($n$ even) reps | §2.2 p.36 | — | X12 |
| $\Sigma$; $\mathcal N$ | spinorial rep; number of supersymmetries | Def 2.4 p.36 | $\Sigma=S_+\otimes W_+\oplus S_-\otimes W_-$; "$\mathcal N=\dim(W)$" ($n\equiv0,1,3,4\bmod8$) | — |
| 4d $\Sigma$ | 4d odd part | §10 p.68 | $\Sigma\cong S_+\otimes W\oplus S_-\otimes W^*$; $S_\pm$ the 2-dim semi-spin reps of $\mathrm{Spin}(4,\mathbb C)\cong SL(2,\mathbb C)\times SL(2,\mathbb C)$, with symplectic pairings $S_\pm\otimes S_\pm\to\mathbb C$ | G10 |
| N=4 $W$, $G_R$ | 4d N=4 R-symmetry | §10.3 p.75 | $G_R=\mathrm{Spin}(6;\mathbb C)\cong SL(4,\mathbb C)$; $W\cong S_+^6$, $W^*\cong S_-^6$ | G20 |
| $\mathrm{Spin}(4,\mathbb R)$ | rotations | §10.3 p.75 | $\mathrm{Spin}(4,\mathbb R)\cong SU(2)_+\times SU(2)_-$ | — |
| $\mathfrak A$ | supertranslation algebra | Def 2.5 p.36 | $\mathfrak A=\Pi\Sigma\oplus V$, "only nontrivial bracket is given by $\Gamma$" | — |
| $\mathfrak t_{\mathcal N}$ | complexified supertranslations | SW §4.1 | $\mathfrak{t}_{\mathcal{N}}=\mathbb{C}^{4}\oplus\Pi(S_{+}\otimes\mathbb{C}^{\mathcal{N}}\oplus S_{-}\otimes\mathbb{C}^{\mathcal{N}})$, $\Gamma:S_{+}\otimes S_{-}\to\mathbb{C}^{4}$ (normalisation not given) | X12 |
| SUSY algebra; $ISO(V_{\mathbb R})$ | — | Def 2.6, p.37 | $(\mathfrak{so}(V)\oplus\mathfrak g_R)\ltimes\mathfrak A$; $ISO(V_{\mathbb R})=\mathrm{Spin}(V_{\mathbb R})\ltimes V_{\mathbb R}$ | — |
| $\mathfrak g$ | gauge Lie algebra | §3.1 p.44; Part II p.51 | "a Lie algebra equipped with a nondegenerate symmetric bilinear pairing"; complexified | G12 |
| $\Pi,\ \lvert x\rvert$ | parity shift; total degree | §1.1 p.17 | $\Pi E$ = $E$ in odd degree; $\lvert x\rvert\in\mathbb Z/2\mathbb Z$ total degree | — |
| $E,\mathcal E,\mathcal E_c,E^!$ | bundle, sections, compact support, dual | §1.1, Intro p.3 | $E^!=E^*\otimes\mathrm{Dens}_M$; "$E\cong E^![-1]$" | — |
| $\mathcal O_{\rm loc},\ \mathcal O^+_{\rm loc}$ | local functionals; at least cubic | §1.1 | $\mathcal O_{\rm loc}(\mathcal E)\subset\mathcal O(\mathcal E_c)$ | — |
| $\omega$ | BV pairing | Def 1.13, 1.16 | $\omega\colon E\otimes E\to\mathrm{Dens}_M[-1]$, graded skew symmetric; $\int_M\omega(e_0,Q_{BV}e_1)=(-1)^{\lvert e_0\rvert}\int_M\omega(Q_{BV}e_0,e_1)$ | matches [C07] $\langle\,,\rangle$ (X1, X8) |
| $X_F,\ \{-,-\}$ | Hamiltonian vector field; bracket | Def 1.15 p.21 | $X_F:=\omega\circ d_{\rm dR}(F)$; $\{F,G\}=X_F(G)$; $\{F,G\}=(-1)^{\lvert F\rvert\lvert G\rvert}\{G,F\}$ | X5 |
| $Q_{BV},S,I$ | BV differential; action; interaction | p.21, Def 1.16 | $S=\frac12\int_M\omega(e,Q_{BV}e)+I$; $Q_{BV}I+\frac12\{I,I\}=0$; CME $\{S,S\}=0$ | X1 |
| $A,\lambda,c$ | BRST fields | Def 3.2 p.44 | $A\in\Omega^1(V_{\mathbb R};\mathfrak g)$, $\lambda\in\Gamma(V_{\mathbb R};\Pi\Sigma\otimes\mathfrak g)$, $c\in\Omega^0(V_{\mathbb R};\mathfrak g[1])$ | — |
| $F_A$ | curvature | §3.1 p.44 | $F_A=\mathrm dA+\frac12[A\wedge A]$ | — |
| $A^*,\lambda^*,c^*$ | antifields | §3.1, Rem 1.36 | "$T^*[-1]F_{\rm gauge}=F_{\rm gauge}\oplus F^!_{\rm gauge}[-1]$"; antifield of $\phi$ is $\phi^*$ | X9 |
| $C^\bullet(\mathfrak g,\mathcal O_{\rm loc}),\ d_{\rm CE}$ | equivariant complex | Def 1.26 p.25 | $\bigoplus_n\mathrm{Hom}(\wedge^n\mathfrak g,\mathcal O_{\rm loc}(\mathcal E))[-n]$; $d_{\rm CE}$ with [Saf17] signs (formula in Def 1.26); $(d_{\rm CE}+Q_{BV})^2=0$ | G7 |
| $S^{(k)}_{\mathfrak g},\ S^{(k)}_{\mathfrak A}$ | $L_\infty$-action components | Def 1.27, §2.5 | $S_{\mathfrak g}=\sum_{k\ge0}S^{(k)}_{\mathfrak g}$, (a) $S^{(0)}=S$, (b) $S^{(k\ge1)}$ "at least quadratic in the fields", (c) $d_{\rm CE}S_{\mathfrak g}+\frac12\{S_{\mathfrak g},S_{\mathfrak g}\}=0$ | G7, G8, G9 |
| $I_{\mathfrak g}$ | — | Rem 1.29 | $I_{\mathfrak g}=I+\sum_{k\ge1}S^{(k)}_{\mathfrak g}$, $(d_{\rm CE}+Q_{BV})I_{\mathfrak g}+\frac12\{I_{\mathfrak g},I_{\mathfrak g}\}=0$ | — |
| $Q$ (ESW) | square-zero supercharge | Def 2.16 | "a nonzero element $Q\in\Sigma$ such that $\Gamma(Q,Q)=0$"; $S^Q=S+\sum_{k\ge1}S^{(k)}_{\mathfrak A}(Q,\dots,Q)$ | twist only; X11 |
| $\alpha,\phi$; $\eta$ | compatible $U(1)$, twisting map; supergravity ghost | Defs 2.18–2.19; Intro p.6 | $\eta$ "is bosonic: it lives in bidegree $(-1,1)$" | X11 |
| $W_{\mathbb R},W,\Gamma_W,G_W$ | dimensional-reduction data | §2.6 p.42 | $V_{\mathbb R}=W_{\mathbb R}\oplus W_{\mathbb R}^\perp$; $\Gamma_W\colon\mathrm{Sym}^2(\Sigma)\xrightarrow{\Gamma_V}V\to W$; $G_W=G_V\times\mathrm{Spin}(W_{\mathbb R}^\perp)$ | X11 |
| $R$ (ESW) | chiral-matter $\mathfrak g$-rep | §3, §10.1 | "Let $R$ be a complex $\mathfrak g$-representation." | X11 |

### Key formulas [ESW]

- **(31) p.44.** $S_{\rm gauge}=\int_{V_{\mathbb R}}\mathrm{dvol}\left(-\frac14(F_A,F_A)+\frac12(\lambda,\not d_A\lambda)+(\mathrm d_Ac,A^*)+([c,\lambda],\lambda^*)+\frac12([c,c],c^*)\right)$. The pairing on $\mathfrak g$ is implicit thereafter.
- **(32).** $S^{(1)}_{\rm gauge}(v)=\int\mathrm{dvol}\,((L_vA,A^*)-(v.\lambda,\lambda^*)-(v.c)c^*)$ for $v\in\mathfrak{iso}(V)$.
- **(33).** $S^{(1)}_{\rm gauge}(Q)=\int\mathrm{dvol}\left(-(\Gamma(Q,\lambda),A^*)+\frac12(\rho(F_A)Q,\lambda^*)\right)$.
- **(34).** $S^{(2)}_{\rm gauge}(Q_1,Q_2)=\int\mathrm{dvol}\left(\frac14(\Gamma(Q_1,Q_2),\Gamma(\lambda^*,\lambda^*))-\frac12(Q_1,\lambda^*)(Q_2,\lambda^*)-\iota_{\Gamma(Q_1,Q_2)}Ac^*\right)$.
- **Thm 3.3.** $d_{\rm CE}(S_{\rm gauge,\mathfrak A})+\frac12\{S_{\rm gauge,\mathfrak A},S_{\rm gauge,\mathfrak A}\}=0$. This is equivalent to four equations:
  - $\{S,S^{(1)}\}=0$;
  - $\{S,S^{(2)}\}+d_{\rm CE}S^{(1)}+\frac12\{S^{(1)},S^{(1)}\}=0$;
  - $d_{\rm CE}S^{(2)}+\{S^{(1)},S^{(2)}\}=0$;
  - $\{S^{(2)},S^{(2)}\}=0$.
- **Lemma 3.6, item 2.** $(d_{\rm CE}S^{(1)})(Q_1,Q_2)=S^{(1)}(\Gamma(Q_1,Q_2))$.
- **Props 2.1–2.3.**
  - Prop 2.1: $\Gamma(Q_1,\rho(X)Q_2)+\Gamma(Q_2,\rho(X)Q_1)=-2\iota_{\Gamma(Q_1,Q_2)}X$.
  - Prop 2.2: $\not d\rho(\Gamma(Q_1,\lambda))Q_2=\rho(d\Gamma(Q_1,\lambda))Q_2+(Q_1,\not d\lambda)Q_2$.
  - Prop 2.3: $\not d(\rho(X)Q)=\rho(\mathrm dX)Q+(-1)^{n(1+p)}\rho(*\mathrm d*X)Q$.
- **Assumption 3.1.** $\rho(\Gamma(Q_1,Q_2))Q_3+\rho(\Gamma(Q_2,Q_3))Q_1+\rho(\Gamma(Q_3,Q_1))Q_2=0$. It holds in the cases 4d $\mathcal N=1$ with $\Sigma=S_+\oplus S_-$ and 10d $\mathcal N=(1,0)$ with $\Sigma=S_+$.
- **Prop 2.20 and §10.3 (route to N=4).** The dimensional reduction has a unique $\mathfrak A'$-supersymmetric structure with $p^*S^{(i)}_{V_{\mathbb R}}=S^{(i)}_{W_{\mathbb R}}$. "The 4d $\mathcal N=4$ super Yang–Mills theory is obtained by dimensional reduction from the 10d $\mathcal N=(1,0)$ super Yang–Mills."
- **Thm 10.9 proof (N=2 view).** N=4 = N=2 YM with matter $U=T^*\mathfrak g$. Fields (§10.2):
  - $A$;
  - $a,\widetilde a\in\Omega^0(\mathfrak g)$;
  - $(\lambda_+,\lambda_-)\in\Omega^0(\Pi S_+\otimes W\otimes\mathfrak g\oplus\Pi S_-\otimes W^*\otimes\mathfrak g)$;
  - $\phi\in\Omega^0(U\otimes W)$;
  - $(\psi_-,\psi_+)\in\Omega^0(\Pi S_+\otimes U\oplus\Pi S_-\otimes U)$;
  - $c\in\Omega^0(\mathfrak g)[1]$.

  Here $\dim W=2$ and $G_R=SL(2;\mathbb C)\times GL(1,\mathbb C)$, with weights "2 on $a$, −2 on $\widetilde a$ and ±1 on $\lambda_\pm,\psi_\pm$".
- **Signature.**
  - "we always consider complexified bundles of fields".
  - "$\mathcal N=1$ super Yang–Mills theory on $M=\mathbb R^4$ with the Euclidean metric".
  - Rem 4.1: square-zero supercharges need the complexified setting.

### Sign conventions [ESW]

- **ES1. Weight.** No path integral and no $e^{\pm S/\hbar}$.
- **ES2. Grading.** Bigrading $\mathbb Z\times\mathbb Z/2$: "$d$ has bidegree (1, 0) while $Q$ has bidegree (0, 1)". $Q_{BV}$ has degree +1.
- **ES3. Ghost degree.** As a section the ghost has degree $-1$ ($c\in\Omega^0(\mathfrak g[1])$), which forces $\deg A^*=+1$ and $\deg c^*=+2$. Coordinate functions carry the opposite sign.
- **ES4. SUSY ghosts.** No field symbol is introduced. Supersymmetry enters as cochains $\mathrm{Hom}(\wedge^n\mathfrak A,-)[-n]$, i.e. functions on $\mathfrak A[1]$. $S^{(k)}(Q,\dots,Q)$ has degree $-k$ and α-weight $k$.
- **ES5. Bracket.** $\{F,G\}=X_F(G)$, with symmetry $(-1)^{\lvert F\rvert\lvert G\rvert}$. The bracket has degree $-n$, which is $+1$ for $n=-1$.
- **ES6. CME.** $\{S,S\}=0$; $Q_{BV}I+\frac12\{I,I\}=0$; $d_{\rm CE}S_{\mathfrak g}+\frac12\{S_{\mathfrak g},S_{\mathfrak g}\}=0$.
- **ES7. Clifford conventions.** $+2(v_1,v_2)$; $\mathfrak{so}$ map $\omega\mapsto(w\mapsto-2\iota_{(w,-)}\omega)$.
- **ES8. Internal inconsistencies.**
  - Thm 2.10 prints a repeated term; Assumption 3.1 is cyclic.
  - p.21 says $S$ has "cohomological degree $1-n$", against degree 0 in Defs 1.16 and 1.27.

## 4. SUPERCONFORMAL [SW]

| Symbol | Meaning | Location | Verbatim formula |
|---|---|---|---|
| $\mathfrak{so}(p+1,q+1)$ | conformal algebra, signature $(p,q)$ | Def 2.1 | Euclidean $\mathbb R^4$ not written out |
| $P_\mu,M_{\mu\nu},\Delta,K_\mu,E$ | conformal vector fields | Prop 2.2, (3) | $P_{\mu}=\frac{\partial}{\partial x^\mu},\ M_{\mu\nu}=x_{\mu}\frac{\partial}{\partial x^\nu}-x_{\nu}\frac{\partial}{\partial x^\mu},\ \Delta=-E,\ K_{\mu}=\lvert x\rvert^{2}\frac{\partial}{\partial x^\mu}-2x_{\mu}E$, $E=x^{\mu}\partial_{\mu}$ |
| $D$ | used in (4), not defined | (4) | [derived] (4) matches (3) only if $D\equiv\Delta=-E$ |
| $g_{\mu\nu}$; $\epsilon_{\alpha\beta}$ | metric; spinor form | §2.1 opening | $\epsilon_{\alpha\beta}$ "$\mathfrak{su}(2)$-invariant alternating form"; raising convention not given |
| $\mathfrak{sconf}(4,\mathcal N)=\mathfrak c(4,\mathcal N)$ | superconformal algebra | (7) | $\mathfrak{su}(2,2\lvert\mathcal{N})$ for $\mathcal N\neq4$; $\mathfrak{psu}(2,2\lvert 4)$ for $\mathcal{N}=4$ |
| $\mathfrak{sl}(4\lvert\mathcal N)$, $\mathfrak{psl}(4\lvert4)$ | complexified versions | after (7) | for $\mathcal N=4$, quotient by the one-dimensional center of $\mathfrak{sl}(k\lvert k)$ |
| $P_{\alpha\dot\alpha},Q^i_\alpha,\bar Q_{\dot\alpha i},M_{\alpha\beta},\Delta,R^i_j,\bar M_{\dot\alpha\dot\beta},S^\alpha_i,\bar S^{\dot\alpha i},K^{\alpha\dot\alpha}$ | generators | (8) | vertical position = conformal weight; horizontal position = (#chiral − #antichiral) spinor indices; parity = $\mathbb Z$-grading mod 2 |
| matrix presentation | blocks ordered $(\alpha,\dot\alpha\,\lvert\,i)$ | Rem 2.9, (17) | $\left[\begin{array}{cc\lvert c}M^{\alpha}_{\beta}&K^{\alpha\dot{\beta}}&S^{\alpha}_{j}\\P_{\alpha\dot{\beta}}&\bar{M}_{\dot{\alpha}}^{\dot{\beta}}&\bar{Q}_{\dot{\alpha}j}\\\hline Q_{\alpha}^{i}&\bar{S}^{\dot{\alpha}i}&R^{i}_{j}\end{array}\right]$; $\Delta=\mathrm{diag}(1,-1\,\lvert\,0)$ (block form in source); $r=\frac{1}{\mathcal{N}-4}\mathrm{diag}(\mathcal N,\mathcal N\,\lvert\,4)$ |
| $e_{\mu\nu}$ | elementary supermatrices | Thm 2.10 proof, (19) | $[e_{0+},e_{\mu\nu}]=\delta_{+\mu}e_{0\nu}-(-)^{\lvert\mu\nu\rvert}\delta_{\nu 0}e_{\mu+}$ |
| $\pm$; $0,i$ | index splittings under $SO(4)\to U(2)$ and $U(\mathcal N)\to U(1)\times U(\mathcal N-1)$ | (18) | holomorphic supercharge $Q^0_+=e_{0+}$ |
| $x_{\alpha\dot\alpha}$ | coordinate on the open cell of $\mathrm{Fl}(2;4)$ | Rem 2.3, (5)–(6) | ${\rm Fl}(2;4)=GL(4,\mathbb{C})/P$, representative $\begin{bmatrix}1&0\\x_{\alpha\dot{\alpha}}&1\end{bmatrix}$ |
| $V,S_+,R,\Pi$ | chiral superspace | Prop 2.4, (9) | $\mathbb{C}^{4\lvert 2\mathcal{N}}\cong V\oplus\Pi(S_{+}\otimes R)$, $R$ = defining rep of $U(\mathcal N)$ |
| superspace cells | — | Prop 2.8, (16) | left-chiral superspace open dense in $\mathrm{Fl}(2\lvert0;4\lvert\mathcal N)$, representative with $\theta^\alpha_i$; right-chiral in $\mathrm{Fl}(2\lvert\mathcal N;4\lvert\mathcal N)$, representative with $\bar\theta^{\dot\alpha i}$ |
| N=1 vector fields | — | (11)–(13) | $P_{\alpha\dot{\alpha}}=\partial_{\alpha\dot{\alpha}}$, $M_{\alpha\beta}=x_{\alpha\dot{\alpha}}\partial_{\beta\dot{\alpha}}+\theta_{\alpha}\partial_{\beta}+(\alpha\leftrightarrow\beta)$, $\bar{M}_{\dot{\alpha}\dot{\beta}}=x_{\alpha\dot{\alpha}}\partial_{\alpha\dot{\beta}}+(\dot{\alpha}\leftrightarrow\dot{\beta})$, $\Delta=E+\frac{1}{2}\theta^{\alpha}\partial_{\alpha}$, $R=\theta^{\alpha}\partial_{\alpha}$, $K^{\alpha\dot{\alpha}}=\lvert x\rvert^{2}\partial^{\alpha\dot{\alpha}}-2x^{\alpha\dot{\alpha}}E+\theta^{\alpha}x^{\beta\dot{\alpha}}\partial_{\beta}$, $Q_{\alpha}=\partial_{\alpha}$, $\bar{Q}_{\dot{\alpha}}=\theta^{\alpha}\partial_{\alpha\dot{\alpha}}$, $S^{\alpha}=-x^{\alpha\dot{\alpha}}\theta^{\beta}\partial_{\beta\dot{\alpha}}+\theta^{2}\partial^{\alpha}$, $\bar{S}^{\dot{\alpha}}=x^{\alpha\dot{\alpha}}\partial_{\alpha}$ |
| $D_{\alpha i},\bar D^j_{\dot\alpha}$ | used, not defined | (10) | $[X,D_{\alpha i}]\sim D_{\alpha i},\quad[X,\bar{D}_{\dot{\alpha}}^{j}]\sim D_{\dot{\alpha}}^{j}$ (sic) |
| $\varepsilon_i$ | odd coordinates of degree $-1$ | §3.1.1 | $\mathcal{O}(X^{d\lvert m})=\mathcal{O}(X)[\varepsilon_{1},\ldots,\varepsilon_{m}]$ |
| $\mathcal G_{\mathcal N}$, $\mathfrak X_{\mathcal N}$ | twisted currents, vector fields | Defs 3.1, 3.2 | twisted theory only |
| N=4 twisted fields, hCS action | twisted theory | Prop 4.5, (29), §4.1 | $\underline{\underline{A}}\in\Omega^{0,\bullet}(\mathbb{C}^{2\lvert 3},\mathfrak{g})$; $\int_{\mathbb{C}^{2\lvert 3}}{\rm d}^{2}z\,{\rm d}^{3}\varepsilon\,\left(\langle\underline{A}\wedge\bar\partial\underline{A}\rangle+\frac{1}{6}\langle\underline{A}\wedge[\underline{A},\underline{A}]\rangle\right)$ |
| $\underline\partial$, $\mathfrak X^{div}_4$, $I(\underline\xi)$ | twisted symmetry and current | Def 4.8, Prop 4.9 | $\frac{1}{2}\int_{\mathbb{C}^{2\lvert 3}}{\rm d}^{2}z{\rm d}^{3}\varepsilon\,\langle\underline{A}L_{\underline{\xi}}\underline{A}\rangle$ |
| $Q^{GF},P_{\epsilon<L},K_\epsilon,W,I[L],\Delta_L,\{,\}_L,\Theta$ | holomorphic renormalisation | §7.2–7.3 (N=2) | $Q^{GF}=\bar\partial^{*}\otimes 1$; $I[L]+I_{S}[L]\overset{\rm def}{=}\lim_{\epsilon\to 0}W(P_{\epsilon<L},I+I_{S})$ |

### Key statements [SW]

- **After (7).** $Q^i_\alpha$ and $\bar S^{\dot\alpha i}$ form a chiral spinor of $\mathcal O(4,2)$, equivalent to the antifundamental of $\mathfrak{su}(2,2)$, in the fundamental of the R-symmetry. $\bar Q_{\dot\alpha i}$ and $S^\alpha_i$ are in the fundamental of $\mathfrak{su}(2,2)$ and the antifundamental of the R-symmetry. The source follows [DO09] (not read).
- **Prop 2.4.** $\mathfrak{sl}(4\lvert\mathcal N)$ acts on $\mathbb C^{4\lvert2\mathcal N}$, through $\mathfrak{psl}(4\lvert4)$ when $\mathcal N=4$. Explicit vector fields are given only for $\mathcal N=1$.
- **Thm 2.10.** $H(\mathfrak c(4),Q)=\mathfrak{psl}(3\lvert3)$ (twist only).

### Sign conventions [SW]

- **SC1. Conformal brackets (4).**
  - $[D,P_\mu]=P_\mu$, $[D,K_\mu]=-K_\mu$.
  - $[M_{\mu\nu},K_\rho]=g_{\rho\nu}K_\mu-g_{\mu\rho}K_\nu$, $[M_{\mu\nu},P_\rho]=g_{\rho\nu}P_\mu-g_{\mu\rho}P_\nu$.
  - $[K_\mu,P_\nu]=2M_{\mu\nu}+2g_{\mu\nu}E$.
  - $[M_{\mu\nu},M_{\rho\sigma}]=g_{\mu\sigma}M_{\nu\rho}+g_{\nu\rho}M_{\mu\sigma}-g_{\mu\rho}M_{\nu\sigma}-g_{\nu\sigma}M_{\mu\rho}$.
- **SC2. Dilatation weights** (see X13).
  - (12) [derived]: $[\Delta,P]=-P$ and $[\Delta,Q_\alpha]=-\tfrac12Q_\alpha$.
  - (17) [derived]: $P:-2$, $K:+2$, $Q,\bar Q:-1$, $S,\bar S:+1$.
- **SC3. N=1 bracket [derived].** $\{Q_\alpha,\bar Q_{\dot\alpha}\}=\partial_{\alpha\dot\alpha}=P_{\alpha\dot\alpha}$, with unit coefficient and no $i$.
- **SC4. R-charge [derived].**
  - $\mathrm{str}\,r=0$.
  - $[r,Q]=-Q$, $[r,\bar S]=-\bar S$, $[r,\bar Q]=+\bar Q$, $[r,S]=+S$.
  - $r$ is singular at $\mathcal N=4$.
- **SC5. QME (§7.2).** $\bar\partial I[L]+\frac{1}{2}\{I[L],I[L]\}_{L}+\hbar\Delta_{L}I[L]=0$.
- **SC6. Anomaly.**
  - $\Theta[L]=\hbar^{-1}\big(\bar\partial(I[L]+I_{S}[L])+\frac{1}{2}\{I[L]+I_{S}[L],I[L]+I_{S}[L]\}_{L}+\hbar\Delta_{L}(I[L]+I_{S}[L])\big)$.
  - $\hbar\Theta=\lim_{L\to 0}\lim_{\epsilon\to 0}\sum_{\Gamma\in{\rm Wheel}_{3},e}W_{\Gamma,e}(P_{\epsilon<L},K_{\epsilon},I+I_{S})$.
- **SC7. Renormalisation.** No counterterms (Lemma 7.5). The weight $e^{\pm S/\hbar}$ and the definition of $W$ are not stated; they are deferred to [Wil18] and [Cos11].
- **SC8. Degrees.** $\varepsilon_i$ has degree $-1$; $\underline A\in\Omega^{0,\bullet}[1]$, $\underline B\in\Omega^{2,\bullet}$; $\hbar=2\pi i$ (Rem 7.12, VOA only).

## 5. CONFLICTS

### 5a. Sign and normalisation conflicts

| # | Item | Sources | Resolution |
|---|---|---|---|
| X1 | Names of action, flow, scale and counterterm | [C07]: $S$ = interaction, $\Gamma(P,S)$, $T$, $S^{CT}$, $S^R=S-S^{CT}$. [CG2]: $I[\Phi]$, $W$, $\Phi$/$L$, $I^{CT}[\epsilon]$. [SW]: $I[L]$, $W(P_{\epsilon<L},\cdot)$. [ESW]: $S$ = full action $\frac12\int\omega(e,Q_{BV}e)+I$, $I$ = interaction | Program formulas use [C07] symbols with $P(\varepsilon,\infty)$. When ESW content is imported: ESW $I$ ↔ [C07] $S$, ESW $Q_{BV}$ ↔ [C07] $Q$, ESW $\int\omega$ ↔ [C07] $\langle\,,\rangle$ (both are $\frac12\langle e,Qe\rangle$ + interaction). CG/SW formulas are quoted in their own symbols |
| X2 | BV Laplacian vs kernel | [C07] $\Delta_T=-\partial_{K_T}$; [CG2] $\triangle_\Phi=\partial_{K_\Phi}$. Both have $QP(\varepsilon,L)=K_\varepsilon-K_L$ and both give $[Q,\partial_{P(T_1,T_2)}]=\Delta_{T_1}-\Delta_{T_2}$ ([CG2] with $\Phi\leftrightarrow T_1$, $\Psi\leftrightarrow T_2$). [derived] So the sign relating $[Q,\partial_\alpha]$ to $\partial_{Q\alpha}$ in [CG2] is opposite to [C07]'s $[\partial_\phi,Q]=\partial_{Q\phi}$; the [CG2] ledger does not record it | Use the [C07] chain: $\Delta_\varepsilon=-\partial_{K_\varepsilon}$, $[\partial_\phi,Q]=\partial_{Q\phi}$, Lemma 10.2.1. Do not take explicit signs from [CG2]'s $\partial_\alpha$; [CG2] $\triangle$ formulas are cited only at the level of the operator in $Q+\hbar\Delta$ |
| X3 | Where the ½ sits | [C07]: $\partial_\phi=\frac12\sum\partial\partial$, $P=\int(Q^{GF}\otimes1)K_t$. [CG2]: $\partial_\alpha$ is a pairing on $\mathrm{Sym}^2$ (normalisation not recorded), $P(\Phi)=\tfrac12Q^{GF}\Phi$, equated with $\int_0^L(Q^{GF}\otimes1)K_tdt$; whether $Q^{GF}$ acts on one or both factors is not recorded | [C07] |
| X4 | Kernel identity sign | [C07] $QP(\varepsilon,T)=-K_T+K_\varepsilon$ = [CG2] §7.2.5/§13.2.1. [CG2] §14.2.2: $\bar\partial P(\epsilon,L)=K_L-K_\epsilon$ | [C07]; [CG2] §14.2.2 is the βγ example and is not used |
| X5 | Bracket | [C07] and [CG2] §7.2.8: $\{f,g\}=\Delta(fg)-\Delta(f)g-(-1)^{\lvert f\rvert}f\Delta(g)$, [derived] symmetry $(-1)^{\lvert f\rvert\lvert g\rvert}$, which equals [ESW] Def 1.15. [CG2] §6.3 differs by $(-1)^{\lvert a\rvert}$ | [C07] Lemma 4.3.2 bracket. [C07] has no scale-ε bracket; if one is needed, use the [CG2] §7.2.8 form with [C07]'s $\Delta_\varepsilon$ |
| X6 | ½ in the (equivariant) QME | ½ in [C07] CS2, [CG2] (except Def 13.2.1.1), [SW] §7.2, [ESW] Def 1.16/1.27. Missing in [C07] §1.11 and [CG2] Def 13.2.1.1 | ½ everywhere; the two omissions are source typos |
| X7 | Exponent weight | [C07] $e^{+S/\hbar}$, with real slice $\langle x,Qx\rangle$ negative definite and the §5.1 integrand with $+\tfrac12\phi\Delta\phi$. [CG2] formal $e^{+I/\hbar}$ vs Remark $e^{-S/\hbar}$. [ESW] and [SW] state none | [C07] $e^{+S/\hbar}$. The relative sign against ESW's $-\frac14(F,F)$ is G14 |
| X8 | Adjointness | [C07]: $Q$ "skew self adjoint" (§1, §2.1) vs "self adjoint" (§4.3); $Q^{GF}$ "self adjoint" (§1.1) vs "skew self adjoint" (4.3.3). [CG2]: $Q$ skew self adjoint, $Q^{GF}$ self adjoint. [ESW] Def 1.16: $\int\omega(e_0,Q_{BV}e_1)=(-1)^{\lvert e_0\rvert}\int\omega(Q_{BV}e_0,e_1)$ | For $Q$, the explicit [ESW] Def 1.16 formula. For $Q^{GF}$, fix its meaning through [C07] Lemma 6.3.1 once a kernel symmetry is chosen (G15) |
| X9 | Ghost degree: section vs coordinate | Sections in degree $-1$: [C07] §2 ($\mathfrak g[1]$), [C07] §5.2 $\mathscr E^{-1}$, [ESW] Rem 1.33 ($c$; $A^*$ at +1, $c^*$ at +2). Coordinates in degree $+1$: [CG2] $\mathsf c_{ij}$, $C^*(\mathfrak g)=\mathscr O(\mathfrak g[1])$. [C07] §5.2 $\mathscr E=\mathcal B[-1]$ [sic] | Field degrees follow [ESW] Rem 1.33 / [C07] §2. Superconformal ghosts are coordinates on $\mathfrak g[1]$ ([CG2] §11.2.2, [ESW] Rem 1.28): ghosts of even generators are odd, ghosts of odd generators are even ([ESW] "$\eta$ is bosonic") |
| X10 | Grading | [C07] mainly $\mathbb Z/2$, $\mathbb Z$ in §2. [ESW] $\mathbb Z\times\mathbb Z/2$. [CG2] "$\mathbb Z/2\times\mathbb Z$ … Both gradings contribute to signs" | [ESW] bigrading, Koszul sign by total degree $\lvert x\rvert$ |
| X12 | Normalisation of $\Gamma$ and of $\{Q,\bar Q\}$ | [ESW] (28): fixed by $(\rho(v)Q_1,Q_2)=(v,\Gamma(Q_1,Q_2))$, "unique up to a scale". [SW] §4.1: $\Gamma:S_+\otimes S_-\to\mathbb C^4$, unnormalised. [SW] (11)/(13) [derived]: $\{Q_\alpha,\bar Q_{\dot\alpha}\}=P_{\alpha\dot\alpha}$ | Open (G11) |
| X13 | Dilatation normalisation ([SW] internal) | (3)/(4): $[D,P]=+P$ with $D\equiv\Delta=-E$. (12): $[\Delta,P]=-P$. (17): $P$ has weight $-2$. [derived from the SC2 entries] On the bosonic generators, ad of (17)'s $\Delta$ = $2\times$ ad of (12)'s $\Delta$ = $-2\times$ ad of (4)'s $D$ | Use the (17) supermatrix presentation, the only one that covers $Q^i,S,R$ at $\mathcal N=4$, and state weights relative to it. (3)/(4) are cited only with the factor above |
| X14 | Signature and real form | [SW] (7) Lorentzian $\mathfrak{psu}(2,2\lvert4)$, Def 2.1 $\mathfrak{so}(p+1,q+1)$. [ESW] Euclidean, complexified; Rem 4.1 | Complexified $\mathfrak{psl}(4\lvert4)$, $\mathrm{Spin}(4,\mathbb C)\cong SL(2,\mathbb C)\times SL(2,\mathbb C)$, Euclidean metric; no real form used |
| X15 | YM formulation | [C07] §5.2: first-order on $\mathcal B$; the standard formulation does not fit its gauge-fixing type. [ESW] (31): second-order $-\frac14(F_A,F_A)$ | Unresolved (G15) |
| X17 | Square-zero parameters | [C07] $\delta$ odd (12.0.1). [CG2] $\delta$ of degree $-i$ (Def 8.3.1.1), $\epsilon$ of degree $-k$ (Prop 13.4.0.1). [SW] $\delta$ of degree one, $=\varepsilon_3$ | $\delta$ = [C07] obstruction parameter. An observable parameter in the same formula has no clean source symbol (G6) |

### 5b. Symbol collisions

| Symbol | Meanings by source | Program use |
|---|---|---|
| $Q$ | [C07]/[CG2] linear BV differential; [ESW] square-zero $Q\in\Sigma$, with $Q_{BV}$ the differential; [SW] supercharges $Q^i_\alpha$ | $Q$ = [C07]; supercharges always carry indices, $Q^i_\alpha,\bar Q_{\dot\alpha i}$ ([SW]) |
| $W$ | [CG2] RG flow; [ESW] multiplicity space, and §2.6 a spacetime subspace; [SW] graph weight $W_{\Gamma,e}$ | [ESW] multiplicity space only; the RG flow is [C07] $\Gamma$ |
| $\Delta,\triangle,D,H$ | [C07] $\Delta_T$ BV Laplacian, $H$ generalised Laplacian, §15.5 $D$ odd operator commuting with $H$; [CG2] $\triangle_\Phi$, $D$ generalised Laplacian; [SW] $\Delta$ dilatation, $D$ in (4), $\Delta_L$ BV Laplacian, $D_{\alpha i}$ | BV Laplacian always subscripted ($\Delta_\varepsilon,\Delta_\infty$); unsubscripted $\Delta$ = dilatation, only in [SW] formulas; $H$ for the generalised Laplacian |
| $K$ | [C07]/[CG2] heat kernel $K_t,K_\varepsilon,K_\infty,K_\Phi$, convolution $K\star$; [SW] $K_\mu,K^{\alpha\dot\alpha}$ | subscript $t,\varepsilon,\infty$ = kernel; vector or spinor index = special conformal |
| $P$ | [C07] propagator $P(\varepsilon,T)$; [SW] $P_\mu,P_{\alpha\dot\alpha}$, parabolic subgroup $P$ | two-slot argument = propagator; index = translation |
| $\Gamma$ | [C07] $\Gamma(\Phi,S)$; [ESW] $\Gamma(Q_1,Q_2)$; [SW] $\Gamma:S_+\otimes S_-\to\mathbb C^4$; [CG2] parametrix in Lemma 8.6.1.2 | distinguished by argument type (propagator and functional vs two spinors) |
| $M$ | [C07] manifold; [ESW] Clifford module $M^\pm$, manifold; [SW] $M_{\mu\nu},M_{\alpha\beta}$ | indexed = rotation |
| $L$, $\mathcal L$ | [C07] $L=\operatorname{Im}Q^{GF}$, $L_t$; [CG2]/[SW] scale $L$; [CG2] $\mathcal L$ symmetry algebra; [ESW] $L$ $L_\infty$ algebra | [C07] $L$, $L_t$; $\mathcal L$ in [CG2] statements |
| $E$ | [C07]/[CG2]/[ESW] bundle; [SW] Euler field $x^\mu\partial_\mu$ | — |
| $S$ | [C07] interaction; [ESW] full action, spin rep $S$, $S_\pm$; [SW] $S^\alpha_i$, $S_\pm$; [C07] $S^{CT},S^R$ | $S$ unindexed = [C07] interaction; $S_\pm$ = semi-spin; $S^\alpha_i,\bar S^{\dot\alpha i}$ = generators |
| $R$ | [ESW] matter rep; [SW] defining rep of $U(\mathcal N)$, $R^i_j$, N=1 $R=\theta\partial$, $r$; [C07] $S^R,\Gamma^R$ | superscript $R$ = renormalised; $R^i_j$ = R-symmetry |
| $\rho$ | [ESW] Clifford action; [SW] twisting homomorphism; [CG2] $\rho_0,\rho_1$ background fields | [ESW] |
| $\alpha$ | [SW] spinor index; [CG2] anomaly, and argument of $\partial_\alpha$; [ESW] compatible $U(1)$; [C07] test element | spinor index; anomaly only in [CG2] statements |
| $\omega$ | [ESW] BV pairing; [CG2] background field in $\int\omega O$, and $\Omega^n$ forms | [ESW] pairing ↔ [C07] $\langle\,,\rangle$ |
| $O$ | [C07] obstruction $O(S),O_{(I,K)}$; [CG2] observable $O[\Phi]$; program $O(\varepsilon)$ | obstruction always written with a functional argument or $(I,K)$; observable with $\varepsilon$ |
| $I,K$ | [C07] multi-index $(I,K)$; [CG2]/[ESW]/[SW] interaction $I$ | — |
| $\varepsilon,\epsilon$ | [C07] UV cutoff, plus the §6.2 formal parameter; [CG2] cutoff (§13.7) and square-zero parameter (§13.4); [SW] $\varepsilon_i$ odd coordinates, $\epsilon$ cutoff; [ESW] Koszul sign in Def 1.26 | $\varepsilon$ = UV cutoff only |
| $T$ | [C07] scale; [CG2] $T^{(n)}$ simplicial set | not used |
| $c$ | [ESW] gauge ghost; [CG2] $\mathsf c_{ij}$ $so(n)$ ghosts | [ESW] $c$ for the gauge ghost |
| $\Sigma$ | [ESW] spinorial rep ($=M^+$ in general; in 4d $S_+\otimes W\oplus S_-\otimes W^*$) | [ESW] |
| $\Phi$ | [C07] even kernel in $\Gamma(\Phi,S)$, $\Phi_r$; [CG2] parametrix | [C07] |

## 6. GAPS

| # | Gap | Options | Closest source |
|---|---|---|---|
| G1 | $P(\varepsilon,\infty)$ on non-compact $\mathbb R^4$ (massless). [C07] needs compact $M$ (heat kernel, $\operatorname{Ker}H$, 4.3.3(5)); the $m=0$ case violates the axioms (§5.1); $\mathbb R^n$ appears only in Chern–Simons with $\mathscr E_c$ and $0<T<\infty$ (Prop 15.2.5). [CG2] has no $T=\infty$ and notes non-unique heat kernels | (a) compact $M$ as in [C07] §4.3, so that $P(\varepsilon,\infty)$, $K_\infty$ and Lemma 10.3.2 hold as stated (superconformal data on $M$ not in any source); (b) $\mathbb R^4$, $\mathscr E_c$, flat $K_t$ of [C07] §15.2 with lower limit $\lVert x-y\rVert^2/4T\to0$ in a Lemma 15.2.2-type formula (the source fixes $T<\infty$) | [C07] §15.2, Lemma 15.2.2 |
| G2 | Locality and $T$-independence of $S^{CT}$ without finite $T$ ([C07] Lemmas 8.0.3–8.0.4 use $T,T'$ and $T\to0$; [CG2] §13.7.2 uses $\mathrm{Sing}_\epsilon$ of $\int_\epsilon^L$) | (a) compact $M$: [C07] Thm B, counterterm $T$-independent (§1.5), evaluated with $P(\varepsilon,\infty)$; (b) direct singular-part extraction from the $\varepsilon$-expansion of the one-loop graphs on $\mathbb R^4$; no source proves locality in this setting | [C07] Thm B, Lemma 8.0.4 |
| G3 | Sign and normalisation of the UV-boundary special edge. [C07] (12.0.1) puts the operator $\Delta_\varepsilon$ in the kernel slot, while §15.3/15.5 speak of a propagator $\delta K_\varepsilon$, relation unstated; ħ is dropped in 10.2.2 and §12; the proof of 12.0.4 writes "$P(\varepsilon,T)+\Delta_\varepsilon$" | Fix it from [C07] Lemma 10.2.1 at T=∞, $QP(\varepsilon,\infty)=-K_\infty+K_\varepsilon$, with $\Delta_\varepsilon=-\partial_{K_\varepsilon}$ and $[\partial_\phi,Q]=\partial_{Q\phi}$ | [C07] (12.0.1), Prop 12.0.3; [SW] $W_{\Gamma,e}(P_{\epsilon<L},K_\epsilon,\cdot)$ |
| G4 | QME/obstruction at T=∞. [C07] Def 10.3.1 at $T=\infty$ needs $K_\infty\in\operatorname{Ker}H\otimes\operatorname{Ker}H$ (compact); no source gives a $K_\infty$ on $\mathbb R^4$ | (a) $O(S)=0$ via [C07] Prop 12.0.3 with $T=\infty$; (b) the [CG2] §13.7.2 form $\alpha^{(1)}(\rho_0)=\lim_{\epsilon\to0}(\triangle_\epsilon I^{(1)}_0(\rho_0)-I^{CT}[\epsilon](\mathrm d\rho_0))$, the only formula in any source whose final expression has no IR scale | [CG2] §13.7.2; [C07] Prop 12.0.3 |
| G5 | Composite operators $O(\varepsilon)$ with counterterms. No source gives an explicit one. [CG2]: abstract $O[\Phi]$ (Def 8.3.1.3); background field $\int\omega O$ with existence "without the QME" (Lemma 10.1.1.3); tree level only mod $\hbar$ (Lemma 8.6.1.2). [CG2] Remark puts $O$ at weight $\hbar^0$, Def 8.3.1.1 puts $\delta O$ inside $e^{(I+\delta O)/\hbar}$ | (a) [C07] Thm B applied to $S+\delta O$ over $\mathbb C[\delta]/\delta^2$ ([CG2] Def 8.3.1.3 text: a family over $\mathbb C[\delta]/\delta^2$, locality holds for local $O$). Under [C07] CS4 the counterterm enters with a minus sign; the form "$O+O^{CT}(\varepsilon)$" is ledger phrasing, not a source formula. (b) The [CG2] background-field route of Lemma 10.1.1.3 | [CG2] Def 8.3.1.3, Lemma 10.1.1.3; [C07] Thm B |
| G6 | Two square-zero parameters in one formula (obstruction and observable) | [C07] $\delta$ with [CG2] Def 8.3.1.1 $\delta$ (clash), or with [CG2] Prop 13.4.0.1 $\epsilon$ (clashes with the cutoff); no source has both | none |
| G7 | Superconformal ghosts and the $L_\infty$ action of $\mathfrak{psl}(4\lvert4)$ on untwisted N=4. [ESW] covers only $\mathfrak A$; [SW] gives generators but no action on fields; no source names ghost coordinates (except [CG2] $\mathsf c_{ij}$) | (a) [ESW] Def 1.27 with $\mathfrak g=\mathfrak{psl}(4\lvert4)$ (super Lie algebras allowed): $S_{\mathfrak g}\in C^\bullet(\mathfrak g,\mathcal O_{\rm loc})$, $d_{\rm CE}S_{\mathfrak g}+\frac12\{S_{\mathfrak g},S_{\mathfrak g}\}=0$, cochain arguments labelled by the [SW] (8) generators; (b) [CG2] Def 12.2.4.1 $\mathrm{Act}(\mathfrak g,\mathcal M)$ and Lemma 12.2.4.2 (local $\Omega^*\otimes\mathfrak g$) | [ESW] (32)–(34), Thm 3.3; [CG2] Def 12.2.4.1 |
| G8 | Quadratic-antifield terms ($S^{(2)}$ and possibly $S^{(k\ge3)}$) for the superconformal action. Only the supertranslation $S^{(2)}$ (34) exists ($\lambda^*\lambda^*$ and $\iota_{\Gamma(Q_1,Q_2)}Ac^*$); the tower stops at $k=2$ for $\mathfrak A$ (Thm 3.3) | Solve the Thm 3.3 tower for $C^\bullet(\mathfrak{psl}(4\lvert4),\mathcal O_{\rm loc})$, starting from the reduced (33)–(34). No source covers pairs involving $S,\bar S,K,\Delta,R$ | [ESW] (34), Thm 3.3 |
| G9 | Ghost-dependent terms quadratic in fields vs [C07] Thm B. [C07] needs $S$ at least cubic mod $\hbar$, with families only over $C^\infty(X)$, $\Omega^*(\Delta^n)$, $\mathbb C[[\hbar]]$, $\mathscr A$. [ESW] Def 1.27(b): $S^{(k)}$ is at least quadratic | (a) the [CG2] counting $\mathscr O^+_{P,sm}(\mathcal L[1]\oplus\mathcal M[1])$, with ghosts counted as fields ($\mathrm{Act}^q$); (b) [C07] families over the dga $C^*(\mathfrak g)$, not in [C07] | [CG2] Def 13.2.2.1, §13.2.2 |
| G10 | Explicit component bases for $S_+,S_-$. [ESW]: abstract 2-dim semi-spin reps with symplectic pairings. [SW]: indices $\alpha,\dot\alpha$, $\epsilon_{\alpha\beta}$ with no raising convention, $x_{\alpha\dot\alpha}$; no $\epsilon_{\dot\alpha\dot\beta}$. [SW] (9)+(16) put $\theta^\alpha_i$ on $\Pi(S_+\otimes R)$, so $S_+$ carries undotted indices there | (a) [SW] index labels, with the raising convention and the dotted alternating form still to be fixed; (b) [ESW] index-free | [SW] §2.1, (9), (16); [DO09] cited, not read |
| G11 | Explicit Γ-matrices (4d and 10d) and the map $x^\mu\leftrightarrow x_{\alpha\dot\alpha}$ | (a) [ESW] (28) with a chosen $(-,-)$ and $v_1v_2+v_2v_1=2(v_1,v_2)$; (b) the [SW] unit-coefficient $\{Q_\alpha,\bar Q_{\dot\alpha}\}=P_{\alpha\dot\alpha}$ | [ESW] (28), Props 2.1–2.3; [SW] (11)–(13) |
| G12 | Normalisation of the coupling: absent everywhere ([ESW] leaves the pairing on $\mathfrak g$ implicit, with $F_A=\mathrm dA+\frac12[A\wedge A]$; [C07] §5.2 has $\mathrm{Tr}$ only) | (a) in the invariant pairing on $\mathfrak g$; (b) in the connection, through $F_A$ | [ESW] (31) |
| G13 | Explicit N=4 action in 4d Weyl components (six scalars, four Weyl fermions, antifields, $S^{(1)},S^{(2)}$). [ESW] gives no reduced formulas and no symbol for the six scalars | (a) [ESW] Prop 2.20 reduction of (31)–(34) via $\mathrm{Spin}(10)\supset\mathrm{Spin}(4)\times\mathrm{Spin}(6)$, $\Sigma=S_+^{10}\to S_+\otimes W\oplus S_-\otimes W^*$, $W\cong S_+^6$; (b) [ESW] Thm 10.9, N=2 + $T^*\mathfrak g$ matter, with the §10.2 fields $A,a,\widetilde a,\lambda_\pm,\phi,\psi_\pm,c$ | [ESW] §10.3, (31)–(34), Prop 2.20 |
| G14 | Euclidean real structure: reality conditions, conjugation on $\Sigma$, sign of $\int\mathrm{dvol}$ relative to $e^{+S/\hbar}$ | (a) fully complexified fields ([ESW]), with the [C07] real slice ($\langle x,Qx\rangle$ negative definite) used only for Gaussian convergence; (b) none other in the sources | [C07] §6.6 |
| G15 | $Q^{GF}$ for the N=4 BV complex satisfying [C07] Def 4.3.3 (fermion sector included) | (a) [C07] §5.2 first-order YM with $Q^{GF}$ the Hermitian adjoint of $Q$, plus N=4 matter; (b) the [ESW] (31) second-order complex with a Hermitian-adjoint $Q^{GF}$, which [C07] states it cannot fit | [C07] §5.2; [SW] $Q^{GF}=\bar\partial^*\otimes1$ (holomorphic only) |
| G16 | Explicit $w_\gamma$, heat-kernel constants, a basis of $\mathscr A$, explicit counterterm ε-functions | One-loop graph class of [C07] §15.4–15.5 ($f_\gamma,F_\gamma$) with the flat $K_t$ of §15.2 | [CG2] §13.7.2 ($\mathrm{STr}_E\rho_1(\not\partial_xK_t)_{y=x}$); [C07] F8, F10 |
| G17 | Action of $K,\Delta$ (and $S,\bar S$) on fields, Weyl weights, conformal Killing data | Extend [ESW] (32) $S^{(1)}(v)$ from $\mathfrak{iso}(V)$ to conformal Killing $v$, and read $S$ from the [SW] (13) vector fields (N=1 superspace only) | [SW] Prop 2.2 (3), (13); [ESW] (32) |
| G18 | Full $\mathfrak{psl}(4\lvert4)$ brackets with one fixed normalisation; the Euclidean real form $\mathfrak{so}(5,1)$ is not written | Supercommutators of the [SW] (17) supermatrices, (19) | [SW] Rem 2.9 |
| G19 | Superconformal anomaly in UV-only form. [CG2] Prop 13.3.2.1 $\alpha$ is $\Phi$-indexed, though $\Phi$-independent; [SW] $\Theta$ uses $\lim_{L\to0}\lim_{\epsilon\to0}$ | $\alpha$ with $\Phi$ replaced by $P(\varepsilon,\infty)$ and the $\varepsilon$-boundary edge of G3 (not in any source) | [CG2] Prop 13.3.2.1, §13.7.2, §13.3.1 |
| G20 | Identification of [ESW] $W$ ($\cong S_+^6$ of $\mathrm{Spin}(6;\mathbb C)\cong SL(4,\mathbb C)$) with the [SW] R-index $i$ / $R$ ($Q^i_\alpha$ "fundamental"; $R$ is the coordinate rep in (9)): $W$ or $W^*$ is not stated | — | [SW] after (7); [ESW] §10.3 |
| G21 | Explicit one-loop QME or anomaly formula beyond the abstract $O_{(I,K)}$, and two-loop specifics | — | [C07] Prop 12.0.3; [CG2] §13.7.2 |

No source uses dimensional regularisation as a regulator: [C07] mentions it only as a hypothetical (§1.6, §9), [CG2] and [SW] use a heat-kernel cutoff, and [ESW] does no quantization.

Working files from the input ledgers:
- /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/costello0706/0706.1533v3.pdf
- /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/cg2/cg2.txt
- /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/esw_pdf/esw.pdf
- /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/sw_src/page.txt