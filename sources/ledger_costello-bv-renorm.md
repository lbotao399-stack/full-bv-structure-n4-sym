# Costello, arXiv:0706.1533: notation ledger

## SOURCE
- K. J. Costello, "Renormalisation and the Batalin-Vilkovisky formalism", arXiv:0706.1533v3 [math.QA], 29 Jun 2007 (v1 11 Jun, v2 17 Jun 2007), 88 pp. The abs page lists no journal reference.
- URLs read:
  - https://arxiv.org/abs/0706.1533
  - https://arxiv.org/pdf/0706.1533v3 (full text, fetched by WebFetch, text layer extracted with pdftotext). Copies: /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/costello0706/0706.1533v3.pdf and /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/costello0706/c.txt
  - https://ar5iv.labs.arxiv.org/html/0706.1533 (LaTeX cross-check; the ar5iv rendering stops at heading 11.2).
- The book (AMS Math. Surveys Monogr. 170, 2011) could not be reached:
  - http://www.math.northwestern.edu/~costello/renormalization.pdf: TLS certificate rejected.
  - https://sites.math.northwestern.edu/~costello/renormalization.pdf: 404.
  - The AMS postpub page: 403.
  - Nothing below comes from the book.
- Formulas are transcribed exactly, source typos included, marked [sic]. Prose of statements is closely paraphrased.
- Dictionary to the book's names:
  - This paper uses $S$ (not $I$), $\Gamma(P,S)$ (not $W(P,I)$), scale $T$ (not $L$), $\Delta_T$, $S^{CT}$, and "singular part".
  - The letter $L$ here is the gauge-fixing subspace $L=\operatorname{Im}Q^{GF}$.

## SYMBOLS

| Symbol | Meaning | Location | Verbatim formula |
|---|---|---|---|
| $M,E,\mathscr E$ | compact manifold, complex super vector bundle, fields | §1, §4.1 | $\mathscr{E}=\Gamma(M,E)$; $\mathscr E\otimes\mathscr E=\Gamma(M^2,E\boxtimes E)$ (completed projective $\otimes$, fn. 1) |
| $\langle\,,\rangle_M$, $\langle\,,\rangle$ | odd symplectic structure | Def 4.3.1 | $\langle\,,\rangle_M:E\otimes E\to\mathrm{Densities}_M$; odd, antisymmetric, with $E\to E^\vee\otimes\mathrm{Densities}_M$ an iso; induces $\langle\,,\rangle:\mathscr E\otimes_{\mathbb C}\mathscr E\to\mathbb C$ |
| $Q$ | linear differential | §1, §2.1, §4.3 | odd, square zero, "skew self adjoint" (§1, §2.1). §4.3 fixed data: "odd linear elliptic differential operator … self adjoint", $[Q,Q]=0$ [conflict] |
| action | — | §1, §4.3 | $\frac{1}{2}\langle e,Qe\rangle+S(e)$, $S\in\mathscr O_l(\mathscr E)$ at least cubic; CME $QS+\tfrac12\{S,S\}=0$ |
| $Q^{GF}$, $L$ | gauge fixing | §1.1, Def 4.3.3 | $L=\operatorname{Im}Q^{GF}$ |
| $H$ | generalised Laplacian | §1.1, 4.3.3(4), §1.9 | $H\overset{\text{def}}{=}[Q,Q^{GF}]$; §1.9 writes $H=[Q^{GF},Q]$ |
| $\mathscr O(\mathscr E)$ | functionals | Def 4.1.1 | $\mathscr{O}(\mathscr{E})=\prod_{n\geq0}\operatorname{Hom}(\mathscr{E}^{\otimes n},\mathbb{C})_{S_n}$ |
| base/families | parameter dependence | Def 4.1.2 | $\mathscr O(\mathscr E,C^\infty(X))=\prod_{n\ge0}\operatorname{Hom}(\mathscr E^{\otimes n},C^\infty(X))_{S_n}$; $\mathscr O(\mathscr E,\mathbb C[[\hbar]])=\varprojlim\mathscr O(\mathscr E)\otimes\mathbb C[\hbar]/\hbar^n$ |
| $\mathscr O_l(\mathscr E)$ | local functionals | Def 4.2.2 | each $f_n$ factors $\mathscr E^{\otimes n}\to\mathrm{Densities}(M)\xrightarrow{\int_M}\mathbb C$, the first map polydifferential. Families: $\mathrm{PolyDiff}(\mathscr E^{\otimes n},\mathrm{Densities}(M))\otimes\Gamma(X,F)$. $\mathscr O_l(\mathscr E,\mathbb C[[\hbar]])=\varprojlim\mathscr O_l(\mathscr E)\otimes\mathbb C[\hbar]/\hbar^n$. "not an algebra" |
| $\{f,g\}$ | odd Poisson bracket | Lemma 4.3.2, §10.1 | well defined for $f\in\mathscr O_l$, $g\in\mathscr O$. For $\dim M=0$: $\{f,g\}=\Delta_0(fg)-\Delta_0(f)g-(-1)^{\lvert f\rvert}f\Delta_0(g)$ |
| $\partial_e$, $Qf$ | derivations | §6.2 | $f(x+\varepsilon e)=f(x)+\varepsilon\partial_e f$; $f(x+\varepsilon Qx)=f(x)+\varepsilon(Qf)(x)$; $[\partial_e,Q]=\partial_{Qe}$ |
| $K\star$ | convolution | §6.3 | $K\star e=(-1)^{\lvert e\rvert}\sum K'\otimes\langle K'',e\rangle$; $(QK)\star e=[Q,K\star]e$ |
| $\partial_\phi$ | order-two operator | §6.4 | $\partial_{\phi}=\frac{1}{2}\sum\partial_{\phi''}\partial_{\phi'}=\frac{1}{2}\sum(-1)^{\lvert\phi'\rvert\lvert\phi''\rvert}\partial_{\phi'}\partial_{\phi''}$, $\phi\in\mathrm{Sym}^2\mathscr E$; $[\partial_\phi,Q]=\partial_{Q\phi}$; Lemma 6.4.1: $[\partial_\phi,e^\vee]=-\partial_{\phi\star e}$ with $e^\vee(e')=\langle e',e\rangle$ |
| $K_t$, $L_t$ | heat kernel | §6.5 | $K_t\star e=e^{-tH}e$; $L_t=(Q^{GF}\otimes1)K_t$; $(Q\otimes1+1\otimes Q)K_t=0$; $\frac{d}{dt}K_t+(Q\otimes1+1\otimes Q)L_t=0$; $K_t+dtL_t$ closed in $\mathscr E\otimes\mathscr E\otimes\Omega^*(\mathbb R_{>0})$ |
| $P(\varepsilon,T)$, $P_\varepsilon^T$ | propagator | §1.2, §6.5, §10.2 | $P(\varepsilon,T)=\int_\varepsilon^T L_t\,dt=\int_\varepsilon^T(Q^{GF}\otimes1)K_t\,dt\in\mathscr E\otimes\mathscr E$ for $0<\varepsilon<T\le\infty$ |
| $P(0,\infty)$ | full propagator | §1.2, Lemma 6.5.1 | $\int_0^\infty(Q^{GF}\otimes1)K_t\,dt$, distributional; represents $Q^{-1}:\operatorname{Im}Q\to\operatorname{Im}Q^{GF}$ |
| $K_\infty$ | — | proof of 10.3.2 | element of $\operatorname{Ker}H\otimes\operatorname{Ker}H$ representing the identity on $\operatorname{Ker}H$ |
| $\Delta_T$ | BV Laplacian | §1.8, §10.1, §2.1 | $\Delta_T=-\partial_{K_T}$, odd, $\Delta_T^2=0$; $\Delta_0$ ill-defined for $\dim M>0$. $\dim M=0$: $K_0=-\sum_i(x_i\otimes\xi_i+\xi_i\otimes x_i)$, $\Delta_0=-\partial_{K_0}=\sum\partial_{x_i}\partial_{\xi_i}$. §2.1: $\Delta=\sum\partial_{x_i}\partial_{\eta_i}$, $\langle x_i,\eta_i\rangle=1$ |
| $\Gamma(\Phi,S)$ | RG flow (book's $W$) | Def 6.6.1, §1.2 | $\Gamma(\Phi,S)=\hbar\log(\exp(\hbar\partial_\Phi)\exp(S/\hbar))\in\mathscr O(\mathscr E,\mathbb C[[\hbar]])$, with $\Phi$ even and $S$ even, at least cubic mod $\hbar$. §1.2 [sic]: $\Gamma(P(\varepsilon,T),S)=\hbar\log(e^{\hbar\partial_{P(\varepsilon,\infty)}}e^{S/\hbar})$ |
| $\Gamma_{i,k}$ | $\hbar^i$, degree $k$ | §1.4, §7.1 | $\Gamma(P(\varepsilon,T),S)=\sum_{i,k\geq0}\hbar^i\Gamma_{i,k}(P(\varepsilon,T),S)$ |
| order | — | proof of Thm B | lexicographic on $(i,k)$; $A_{\le(I,K)}=\sum_{(i,k)\le(I,K)}\hbar^iA_{i,k}$ |
| $\mathscr A$ | $\varepsilon$-functions | Def 7.1.1 | finite spans of $\int_{U\subset(\varepsilon,1)^n}\frac{F(t_1,\dots,t_n)^{1/2}}{G(t_1,\dots,t_n)^{1/2}}dt_1\dots dt_n$ and $\int_{U\subset(\varepsilon,1)^{n-1}}\frac{F(t_1,\dots,t_n=\varepsilon)^{1/2}}{G(t_1,\dots,t_n=\varepsilon)^{1/2}}dt_1\cdots dt_{n-1}$, with $F,G\in\mathbb Z_{\ge0}[t_1,\dots,t_n]\setminus\{0\}$ and $U$ cut out by $t_i^l>t_j$, $l\in\mathbb Z$ |
| $\mathscr A_{\ge0}$, $\mathscr A_{<0}$ | scheme | Def 1.4.1, 8.0.2 | $\mathscr A_{\ge0}$: the $\varepsilon\to0$ limit exists; $\mathscr A=\mathscr A_{<0}\oplus\mathscr A_{\ge0}$; "singular part" = projection onto $\mathscr A_{<0}$ |
| $\mathscr A_{\le0}$, $\mathscr A_{>0}$ | — | §12 | $\mathscr A_{\le0}=\mathscr A_{<0}\oplus\mathbb C$; $\mathscr A_{>0}$: limit exists and $=0$; $\mathscr A=\mathscr A_{\le0}\oplus\mathscr A_{>0}$ |
| $S^{CT}$, $S^R$ | counterterms | §8 | $S^R(\hbar,\varepsilon)=S(\hbar)-S^{CT}(\hbar,\varepsilon)=S-\sum_{i>0,k\geq0}\hbar^iS^{CT}_{i,k}(\varepsilon)$ |
| $\Gamma^R$ | renormalised flow | §10.3 | $\Gamma^R(P(0,T),S)\overset{def}{=}\lim_{\varepsilon\to0}\Gamma(P(\varepsilon,T),S^R)$ |
| $Z$, $Z^R$ | — | §1, §1.2, §1.5 | $Z(S,\hbar,a)=\int_{x\in L}\exp\left(\frac{1}{2\hbar}\langle x,Qx\rangle+\frac1\hbar S(x+a)\right)dx$; $=\lim_{\varepsilon\to0}(e^{\hbar\partial_{P(\varepsilon,\infty)}}e^{S/\hbar})(a)$ up to $\det Q$; $Z^R=\exp(\Gamma^R(P(0,\infty),S)/\hbar)\in\mathscr O(\mathscr E)((\hbar))$ |
| $d\mu$ | — | §6.6 | $\int_{x\in(\operatorname{Im}Q^{GF})_{\mathbb R}}\exp(\tfrac12\langle x,Qx\rangle/\hbar)d\mu=1$, with $\langle x,Qx\rangle$ negative definite on the real slice |
| $S^{eff}(T)$ | effective actions | Def 9.0.7 | see K7 |
| $O(S)$, $O_{(I,K)}$ | obstruction | Prop 12.0.3 | see K11 |
| $w_\gamma$, $\operatorname{Aut}\gamma$ | graph weights | §15.4 | see K12 |
| $\mathrm{GF}$, $\mathrm{BV}_{\mathscr E}$, $\mathrm{BV}_{H^*(\mathscr E,Q)}$ | simplicial sets | §11 | see K10 |

## KEY STATEMENTS

**K1. Def 4.3.3 (gauge fixing condition).** $Q^{GF}:\mathscr E\to\mathscr E$ is an odd linear differential operator with:
1. skew self adjoint (§1.1 says "self adjoint");
2. order $\le1$;
3. $[Q^{GF},Q^{GF}]=0$;
4. $H:=[Q,Q^{GF}]$ is a generalised Laplacian [BGV92]: $\sigma(H)\in\Gamma(M,\mathrm{Sym}^2TM\otimes\mathrm{End}\,E)$ is a positive definite metric on $T^*M$ tensored with $\mathrm{id}_E$;
5. $\mathscr E=\operatorname{Im}Q\oplus\operatorname{Im}Q^{GF}\oplus\operatorname{Ker}H$.

Consequence (§1.9): $H^*(\mathscr E,Q)=\operatorname{Ker}H$, finite dimensional.

**K2. Lemma 6.5.1.** $\alpha\mapsto\int_0^\infty Q^{GF}e^{-tH}\alpha\,dt$ inverts $Q:\operatorname{Im}Q^{GF}\to\operatorname{Im}Q$. In the proof, $\int_0^\infty Q^{GF}e^{-tH}Q\alpha\,dt=\alpha-\pi\alpha$, where $\pi:\mathscr E\to\operatorname{Ker}H$. §6.5: $P(\varepsilon,\infty)\in\mathscr E\otimes\mathscr E$ because $Q^{GF}$ annihilates the 0-eigenspace.

**K3. Lemma 6.6.2 ($M$ a point).** $\Gamma(P(0,\infty),S)(a)=\hbar\log\int_{x\in(\operatorname{Im}Q^{GF})_{\mathbb R}}\exp\left(\tfrac12\langle x,Qx\rangle/\hbar+S(x+a)/\hbar\right)d\mu$. For $\dim M>0$ the source takes this as the definition (§7.1).

**K4. RG semigroup (§1.3).** $\Gamma(P(T_2,T_3),\Gamma(P(T_1,T_2),f))=\Gamma(P(T_1,T_3),f)$. It holds for $f$ at least cubic mod $\hbar$ and $0<\varepsilon<T\le\infty$, and the flow is invertible.

**K5. Theorem A (§7.1).**
1. $\Gamma_{i,k}(P(\varepsilon,T),S(\varepsilon,\hbar))(e)\simeq\sum_{r\ge0}f_r(\varepsilon)\Phi_r(T,e)$, with $f_r\in\mathscr A$ and $\Phi_r\in\mathscr O(\mathscr E,C^\infty((0,\infty)))$.
2. $\Phi_r(T,e)=\sum g_q(T)\Psi_{q,r}(e)$, with $\Psi_{q,r}\in\mathscr O_l(\mathscr E)$.
3. For $k>0$ the germs of $\Psi_{q,r}$ at $x$ depend only on the germ of $(E,\langle\,,\rangle,Q,Q^{GF},S)$.

Meaning of $\simeq$ (§15.3): $\sup_{T\in K}\|\Gamma_{i,k}(P(\varepsilon,T),S)(\alpha)-\sum_{r=0}^Rf_r(\varepsilon)\Phi_r(T,\alpha)\|_l^{\Delta^d}<\varepsilon^{R+1}\|\alpha\|_m$ ("for all T sufficiently small" [sic]).

**K6. Theorem B (§8).** Let $S\in\mathscr O_l(\mathscr E,\mathbb C[[\hbar]])$ be even and at least cubic mod $\hbar$. Then there is a unique $S^{CT}(\varepsilon,\hbar)=\sum_{i\ge1,k\ge0}\hbar^iS^{CT}_{i,k}(s)$ [sic: $(s)$ for $(\varepsilon)$] with:
1. $S^{CT}_{i,k}(\varepsilon)\in\mathscr O_l(\mathscr E,\mathscr A_{<0})$, homogeneous of degree $k$;
2. $\lim_{\varepsilon\to0}\Gamma(P(\varepsilon,T),S-S^{CT})$ exists.

Further properties:
- For $k>0$ the germ of $S^{CT}_{i,k}$ at $x$ is determined by the germ of the data.
- §1.5: the counterterms are independent of $T$.
- Construction: $S^{CT}_{I,K}=$ singular part of $\Gamma_{I,K}(P(\varepsilon,T),S-S^{CT}_{<(I,K)})$.
- Lemma 8.0.3 ($T$-independence) uses $P(\varepsilon,T')-P(\varepsilon,T)=P(T,T')\in\mathscr E\otimes\mathscr E$.
- Lemma 8.0.4: locality follows from $T$-independence plus the $T\to0$ expansion of Theorem A.
- Lemma 8.0.5: $S^{CT}_{0,k}=0$, i.e. no tree-level counterterms.
- §1.5 lists the first steps: $S^{CT}_{1,1}$, then $S^{CT}_{1,2}$ from $S-\hbar S^{CT}_{1,1}$, then $S^{CT}_{2,0}$ from $S-\hbar\sum_kS^{CT}_{1,k}$.
- After the statement: $\Gamma^R(P(0,\infty),S)$ "is a renormalised functional integral".

**K7. Lemma 9.0.6, Def 9.0.7, Theorem C.**
- Lemma 9.0.6: $\Gamma^R(P(0,T'),S)=\Gamma(P(T,T'),\Gamma^R(P(0,T),S))$.
- Def 9.0.7: $S^{eff}(T)\in\mathscr O(\mathscr E,\mathbb C[[\hbar]])$ for each $T>0$, at least cubic mod $\hbar$, with:
  1. $S^{eff}(T_2)=\Gamma(P(T_1,T_2),S^{eff}(T_1))$;
  2. there is $\Phi\in\mathscr O_l(\mathscr E,C^\infty(0,\infty)\otimes\mathbb C[[\hbar]])$ with $\lim_{T\to0}(S^{eff}(T)-\Phi(T))=0$.
- Theorem C: for each $\mathscr A_{<0}$, the map $S\mapsto\{\Gamma^R(P(0,T),S)\mid T\in\mathbb R_{>0}\}$ is a bijection.

**K8. Lemmas 1.8.1, 10.2.1, 10.2.2.**
- Lemma 10.2.1: $QP_\varepsilon^T=-K_T+K_\varepsilon$ ($Q$ = tensor-product differential). In the proof, $QP(\varepsilon,T)$ is the kernel of $\alpha\mapsto Q\int_\varepsilon^TQ^{GF}e^{-tH}\alpha+\int_\varepsilon^TQ^{GF}e^{-tH}Q\alpha=-e^{-TH}\alpha+e^{-\varepsilon H}\alpha$.
- Lemma 10.2.2: $(Q+\hbar\Delta_{T_1})\exp(S/\hbar)=0\iff(Q+\hbar\Delta_{T_2})\exp(\Gamma(P(T_1,T_2),S)/\hbar)=0$.
- Proof identities:
  - $[\partial_{P(T_1,T_2)},Q]=\partial_{QP(T_1,T_2)}=\partial_{K_{T_1}}-\partial_{K_{T_2}}=\Delta_{T_2}-\Delta_{T_1}$;
  - $[Q,\exp(\partial_{P(T_1,T_2)})]=\exp(\partial_{P(T_1,T_2)})(\Delta_{T_1}-\Delta_{T_2})$ [sic, $\hbar$ omitted];
  - $(Q+\hbar\Delta_{T_2})\exp(\hbar\partial_{P(T_1,T_2)})\exp(S/\hbar)=\exp(\hbar\partial_{P(T_1,T_2)})(Q+\hbar\Delta_{T_1})\exp(S/\hbar)$.

**K9. Def 10.3.1 (renormalised QME), Lemma 10.3.2, quantisation.**
- Def 10.3.1: $(Q+\hbar\Delta_T)\exp\left(\hbar^{-1}\Gamma^R(P(0,T),S)\right)=0$. This is independent of $T$ by 10.2.2; Def 1.8.2 says "for some (or equivalently, all) $T>0$".
- Lemma 10.3.2: $\Gamma^R(P(0,\infty),S)|_{\operatorname{Ker}H}$ satisfies the QME on $\operatorname{Ker}H=H^*(\mathscr E,Q)$.
- Quantisation (§1.11, §10.3): $S=S_0+\sum_{i>0}\hbar^iS_i$ solving the renormalised QME. Equivalently, a system $\{S^{eff}(T)\}$ satisfying the QME that converges mod $\hbar$ to $S_0$ as $T\to0$. Only Chern-Simons is quantised, and only modulo constants.

**K10. Families and homotopies.**
- Def 11.1.1: $Q^{GF}:\mathscr E\otimes C^\infty(\Delta^n)\to\mathscr E\otimes C^\infty(\Delta^n)$, $C^\infty(\Delta^n)$-linear and first order, with:
  - $H_0=[Q\otimes1,Q^{GF}]$ a smooth family of generalised Laplacians;
  - $\mathscr E\otimes C^\infty(\Delta^n)=\operatorname{Im}(Q\otimes1)\oplus\operatorname{Im}Q^{GF}\oplus\operatorname{Ker}H_0$;
  - $H=[Q+d_{DR},Q^{GF}]=H_0+H_1$, with $H_1=[d_{DR},Q^{GF}]$ and $K_t\in\mathscr E\otimes\mathscr E\otimes\Omega^*(\Delta^n)$.
- §11.2: $S\in\mathscr O_l(\mathscr E,\Omega^*(\Delta^n)[[\hbar]])$ satisfies $(d_{DR}+Q+\hbar\Delta_T)\exp(\hbar^{-1}\Gamma^R(P(0,T),S))=0$.
- Def 11.2.2: $\mathrm{BV}_{\mathscr E}[n]$ consists of pairs ($Q^{GF}\in\mathrm{GF}[n]$, $S$ solving that equation). $\mathrm{BV}_{H^*}[n]$: $(d_{DR}+\hbar\Delta_{H^*(\mathscr E,Q)})e^{S/\hbar}=0$.
- §1.9 homotopy: $\left(Q+dt\frac{d}{dt}+\hbar\Delta\right)e^{(S_t+dtS'_t)/\hbar}=0$.
- §2.3: $d_{DR}F+\tfrac12\{F,F\}+\hbar\Delta F=0$; with $F=A(t)+dtB(t)$ this gives $\tfrac12\{A,A\}+\hbar\Delta A=0$ and $\frac{d}{dt}A+\{A,B\}+\hbar\Delta B=0$.
- Theorem 11.3.1: $\mathrm{BV}_{\mathscr E}\to\mathrm{GF}$ is a fibration.

**K11. Obstruction theory (§12–14).**
- §12, first equation: RQME $\iff\lim_{\varepsilon\to0}(Q+d_{DR}+\Delta_T)\exp(\hbar\partial_{P(\varepsilon,T)})\exp(S^R/\hbar)=0$ [sic].
- (12.0.1): $\Gamma\big(P(\varepsilon,T)+\delta\Delta_\varepsilon,S^R+\delta QS^R+\delta d_{DR}S^R\big)=\hbar\log\big((1+\delta Q+\delta d_{DR}+\delta\Delta_T)\exp(\hbar\partial_{P(\varepsilon,T)})\exp(S^R/\hbar)\big)$. Here $\delta$ is an odd parameter, and the identity follows from $[\partial_{P(\varepsilon,T)},Q+d_{DR}]=\Delta_T-\Delta_\varepsilon$.
- Prop 12.0.3:
  1. There is a unique odd $O(S)\in\mathscr O_l(\mathscr E,\mathscr A_{\le0}\otimes\Omega^*(\Delta^n)\otimes\mathbb C[[\hbar]])$ with $\lim_{\varepsilon\to0}\frac{d}{d\delta}\langle\Gamma(P(\varepsilon,T)+\delta\Delta_\varepsilon,S^R+\delta QS^R+\delta d_{DR}S^R+\delta O(S))\rangle=0$.
  2. $O_{(I,K)}(S_{<(I,K)}+\hbar^IS_{(I,K)})=\hbar^IQS_{(I,K)}+\hbar^Id_{DR}S_{(I,K)}+O_{(I,K)}(S_{<(I,K)})$.
  3. The germ of $O_{i,k}$ is local in the data.

  RQME $\iff O(S)=0$. Definition: $O_{(I,K)}(S,T)=-$projection onto $\mathscr A_{\le0}$ of $\frac{d}{d\delta}\langle\Gamma_{(I,K)}(\cdots+\delta O_{<(I,K)}(S))\rangle$. It is $T$-independent, hence local.
- Lemma 12.0.4: if $O_{<(I,K)}(S)=0$, then $O_{(I,K)}(S)$ is $\varepsilon$-independent, lies in $\mathscr O_l(\mathscr E,\Omega^*(\Delta^n))$, and $(Q+d_{DR})O_{(I,K)}(S)=0$. The lift is $S+\hbar^IS'$ with $(Q+d_{DR})S'=-O_{(I,K)}(S)$.
- Lemma 14.4.2(3): $O_{(I,K)}(S+\hbar^IS')=O_{(I,K)}(S)+(Q+d_{DR})S'$.
- Theorem D: for triangular $Q$, $\Gamma(M,\mathrm{BV})\to R\Gamma(M,\mathrm{BV})$ is a weak equivalence.

**K12. Graphs (§15.4–15.5).**
- $\Gamma_{i,k}(P(\varepsilon,T),S)(\alpha)=\sum_\gamma\frac{1}{\operatorname{Aut}\gamma}w_\gamma(\alpha)$, summed over connected graphs with $k$ legs.
- Each vertex has a degree in $\mathbb Z_{\ge0}$ equal to its $\hbar$ power. Degree-0 vertices are at least trivalent. $\sum$(vertex degrees) $+\,b_1=i$.
- The formula for $w_\gamma$ is not given.
- General class: $f_\gamma(t_e,\alpha)=\int_{M^{V(\gamma)}}\otimes S_v\big(\alpha\otimes\otimes_{e\in E(\gamma)}K_{t_e}\big)$ and $F_\gamma(\varepsilon,T,\alpha)=\int_{t_e\in(\varepsilon,T)^{E(\gamma)}}f_\gamma(t_e,\alpha)\prod dt_e$. The propagator is $\int_\varepsilon^T(D\otimes1)K_t\,dt$ with $D$ odd and commuting with $H$.
- Remark: the variant propagator $P(\varepsilon,T)+\delta K_\varepsilon$ has one special edge carrying $K_\varepsilon$.

**K13. Yang-Mills (§5.2).**
- First-order formulation [MZ96, Wit04], details in [Cos06]. The author states he does not see how to fit the standard formulation (no gauge fixing of the required type).
- Setting: $M$ compact oriented 4-manifold, $d_AF(A)_+=0$.
- The dga $\mathcal B$ by degree:
  - degree 0: $\Omega^0(\mathrm{End}V)$;
  - degree 1: $\Omega^1\oplus\Omega^2_+$;
  - degree 2: $\Omega^2_+\oplus\Omega^3$;
  - degree 3: $\Omega^4$.
- Arrows ($Q$), read from the diagram:
  - $d_A:\Omega^0\to\Omega^1$ and $-[F(A)_+,\,]:\Omega^0\to\Omega^2_+$;
  - $d_{A+}:\Omega^1\to\Omega^2_+$ and $-[F(A)_+,\,]:\Omega^1\to\Omega^3$;
  - $\mathrm{Id}:\Omega^2_+\to\Omega^2_+$ and $-d_A:\Omega^2_+\to\Omega^3$;
  - $-[F(A)_+,\,]:\Omega^2_+\to\Omega^4$ and $-d_A:\Omega^3\to\Omega^4$.
- A trace of degree $-3$ gives an odd symmetric pairing $\mathrm{Tr}\,ab$. $\mathscr E=\mathcal B[-1]$ [sic], $S(e)=\mathrm{Tr}\,a^3$ [sic].
- $\mathscr E^{-1}=\mathcal B^0$ ghosts, $\mathscr E^0=\mathcal B^1$ fields, $\mathscr E^1=\mathcal B^2$ antifields, $\mathscr E^2=\mathcal B^3$ antighosts.
- $Q^{GF}$ is the Hermitian adjoint of $Q$, using a metric on $M$ and a Hermitian metric on $V$.

**K14. Flat $\mathbb R^n$ (§15.2, Chern-Simons).**
- $\mathscr E_c=\Omega^*_c(\mathbb R^n)\otimes\mathfrak g[1]$ and $P(\varepsilon,T)=\int_\varepsilon^Td^*K_t\,dt$.
- $K_t=Ct^{-n/2}e^{-\|x-y\|^2/4t}d(x_1-y_1)\wedge\cdots\wedge d(x_n-y_n)\otimes I_{\mathfrak g}$.
- Lemma 15.2.2 (up to an overall constant): $P(\varepsilon,T)=\left(\int_{\|x-y\|^2/4T}^{\|x-y\|^2/4\varepsilon}u^{n/2-1}e^{-u}du\right)\pi^*\mathrm{Vol}_{S^{n-1}}\otimes I_{\mathfrak g}$.
- Prop 15.2.5 fixes $0<T<\infty$.

## SIGN CONVENTIONS
1. **Weight.** Weight $e^{+S/\hbar}$, integrand $\exp(\frac1{2\hbar}\langle x,Qx\rangle+\frac1\hbar S(x+a))$. A real slice with $\langle x,Qx\rangle$ negative definite is chosen (§6.6). In the $\phi^4$ example (§5.1) the integrand is $\exp\big(\hbar^{-1}\int_M\tfrac12\phi\Delta\phi+\tfrac12m^2\phi^2+\phi^4\big)$, with $H(\phi e_i)=(\Delta\phi+m^2\phi)e_i$ and $Q^{GF}(\phi e_1)=\phi e_0$.
2. **QME and Laplacian.**
   - QME: $(Q+\hbar\Delta)e^{S/\hbar}=0\iff QS+\tfrac12\{S,S\}+\hbar\Delta S=0$.
   - $\Delta_T=-\partial_{K_T}$.
   - On cohomology, §2.3 has $\Delta e^{f/\hbar}=0$, and §1.11 has $\{\Phi_{CS},\Phi_{CS}\}+\hbar\Delta\Phi_{CS}=0$ [sic, no $\frac12$].
3. **Bracket.** $\{f,g\}=\Delta_0(fg)-\Delta_0(f)g-(-1)^{\lvert f\rvert}f\Delta_0(g)$.
4. **Counterterms.** $S^R=S-S^{CT}$, with $S^{CT}$ starting at $\hbar^1$ and $S^{CT}_{0,k}=0$.
5. **Propagator identities.** $QP(\varepsilon,T)=-K_T+K_\varepsilon$; $[\partial_P,Q]=\partial_{QP}$; $[\partial_{P(\varepsilon,T)},Q+d_{DR}]=\Delta_T-\Delta_\varepsilon$.
6. **Gaussian normalisations.** The $\frac12$ sits inside $\partial_\phi$, and $\Gamma=\hbar\log(e^{\hbar\partial_\Phi}e^{S/\hbar})$. Convolution carries the sign $(-1)^{\lvert e\rvert}$. Measure: $\int\exp(\frac12\langle x,Qx\rangle/\hbar)d\mu=1$; $Z$ holds "up to … the determinant of $Q$".
7. **Adjointness.** $Q$ is skew self adjoint (§1, §2.1) vs self adjoint (§4.3). $Q^{GF}$ is self adjoint (§1.1) vs skew self adjoint (4.3.3). Lemma 6.3.1: $K$ is (anti)symmetric $\iff K\star$ is (skew) self adjoint.
8. **Grading.**
   - Mainly $\mathbb Z/2$ (super); $[1]$ is parity shift (§1.1).
   - §2: $\mathfrak g[1]$ has $\mathfrak g$ in degree $-1$, and $E=\mathfrak g[1]\oplus V\oplus V^\vee[-1]\oplus\mathfrak g^\vee[-2]$.
   - $S$ is even; $\Delta_T$ and $O(S)$ are odd.
   - Families enter as $+d_{DR}$.

## GAPS
1. The book's symbols ($W(P,I)$, $I$, $I^{CT}$, $L$ as a scale, $P(\varepsilon,L)$, $\Delta_L$, $\{\,,\}_L$, "purely singular", $\mathscr O_{loc}$) are not in this paper. The book itself could not be fetched.
2. There is no scale-dependent bracket $\{\,,\}_T$. The only brackets are the $\dim M=0$ formula and the local/non-local bracket of Lemma 4.3.2.
3. $M$ must be compact (heat kernel, $\operatorname{Ker}H$, the decomposition in 4.3.3(5)).
   - The massless case violates the axioms (§5.1, $m=0$).
   - $\mathbb R^n$ appears only in the local Chern-Simons model with $\mathscr E_c$ and $0<T<\infty$.
   - $P(\varepsilon,\infty)$ is treated only for compact $M$ (§1.2, §6.5, §7.1, Lemma 10.3.2).
4. Where the proofs use finite scales: $T$-independence and locality (Lemmas 8.0.3–8.0.4, Prop 12.0.3) use finite $T,T'$ and $T\to0$ asymptotics.
5. Not given in the source:
   - the explicit $w_\gamma$;
   - Feynman-rule normalisations;
   - the heat-kernel constants $C,C',C''$;
   - an explicit basis of $\mathscr A$, or the explicit $\varepsilon$-functions that appear in counterterms.
6. Theory content:
   - None of the program's field content appears: no fermions, Weyl spinors, Euclidean spinor conventions, SUSY, N=4, or superconformal ghosts.
   - Yang-Mills appears only in the first-order, background-field form on a compact 4-manifold.
   - Yang-Mills is renormalised but not quantised; no QME solution is constructed for it.
7. Absent topics:
   - composite operators or observables $O(\varepsilon)$ and their counterterms;
   - explicit one-loop QME or anomaly formulas beyond the abstract $O_{(I,K)}$;
   - two-loop specifics.

   Parameter dependence is available only through $C^\infty(X)$, $\Omega^*(\Delta^n)$, $\mathbb C[[\hbar]]$ and $\mathscr A$-valued coefficients.
8. Source inconsistencies:
   - §1.2 defines $\Gamma(P(\varepsilon,T),S)$ using $P(\varepsilon,\infty)$.
   - $\hbar$ is dropped in the proof of 10.2.2 and in the first equation of §12.
   - (12.0.1) puts the operator $\Delta_\varepsilon$ in the kernel slot, while §15.3/15.5 speak of a "propagator $\delta K_\varepsilon$"; the relation, including its sign, is not stated.
   - The proof of 12.0.4 writes "$P(\varepsilon,T)+\Delta_\varepsilon$" and "$\delta T\Delta_T$".
   - Lemma 14.4.2(3) lacks the $\hbar^I$ that appears in 12.0.3(2).
   - Theorem B writes $S^{CT}_{i,k}(s)$.
   - §5.2 writes $\mathscr E=\mathcal B[-1]$ (against §2's shift convention) and $S(e)=\mathrm{Tr}\,a^3$.
9. No dimensional regularisation is used. DR appears only as a hypothetical reparametrisation in §1.6 and §9. Morava's $\int_0^\infty t^z(Q^{GF}\otimes1)K_t\,dt$ is mentioned in §1.4, unproven and not used.