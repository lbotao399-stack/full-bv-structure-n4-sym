# NOTATION LEDGER: Saberi–Williams, arXiv:1910.04120

## SOURCE

- I. Saberi, B. R. Williams, *Superconformal algebras and holomorphic field theories (or: Towards derived enhancements of superconformal algebras, I)*, arXiv:1910.04120v3 [math-ph], 4 Nov 2021. Versions: v1 9 Oct 2019, v2 28 Feb 2020, v3 4 Nov 2021. The abstract page lists no journal reference.
- Read at: https://arxiv.org/abs/1910.04120 and https://arxiv.org/html/1910.04120v3. The math was taken from the LaTeXML `alttext`; the extracted text is at `/private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/sw_src/page.txt`.
- Why this source was kept: it treats 4d $\mathcal N=4$ explicitly, in eq. (7), §2.1, Prop. 2.4, Thm. 2.10, Prop. 4.5, Prop. 4.9 and §6.5. So the fallback clause does not apply and no other source was used.
- Transcription notes:
  - Formulas are exact LaTeX. Only rendering macros were expanded: `\smash{\underline{X}}`→`\underline{X}`, `\uprightintop\displaylimits`→`\int`, `\partialderivative{}{y}`→`\frac{\partial}{\partial y}`.
  - The prose of theorems is paraphrased, but the formulas inside them are exact. This follows the copyright limit on long verbatim prose.
  - Anything marked **[derived]** is my own one-line check from the source's formulas. It is not in the source.

## SYMBOLS

| symbol | meaning | location | defining formula (verbatim) |
|---|---|---|---|
| $\mathfrak{so}(p+1,q+1)$ | conformal algebra, $d>2$, signature $(p,q)$ | Def. 2.1 | (Euclidean $\mathbb R^4$ is not written out separately) |
| $P_\mu,M_{\mu\nu},\Delta,K_\mu,E$ | conformal vector fields on $\mathbb R^{p,q}$ | Prop. 2.2, eq. (3) | $P_{\mu}=\frac{\partial}{\partial x^\mu},\ M_{\mu\nu}=x_{\mu}\frac{\partial}{\partial x^\nu}-x_{\nu}\frac{\partial}{\partial x^\mu},\ \Delta=-E,\ K_{\mu}=\lvert x\rvert^{2}\frac{\partial}{\partial x^\mu}-2x_{\mu}E$, with $E=x^{\mu}\partial_{\mu}$ |
| $D$ | used in (4), never defined | eq. (4) | **[derived]** (4) matches the vector fields of (3) only if $D\equiv\Delta=-E$ |
| $g_{\mu\nu}$ | metric (vector indices are raised and lowered with it) | §2.1 opening, (4) | — |
| $\epsilon_{\alpha\beta}$ | raises and lowers spinor indices | §2.1 opening | "$\mathfrak{su}(2)$-invariant alternating form"; the raising convention is not given |
| $\mathfrak{sconf}(d,\mathcal N)$, also $\mathfrak c(d,\mathcal N)$ | superconformal algebra | eq. (7) | $\mathfrak{sconf}(d,\mathcal{N})=\mathfrak{su}(2,2\lvert\mathcal{N}),\ d=4,\mathcal N\neq4;\ \mathfrak{psu}(2,2\lvert 4),\ d=4,\ \mathcal{N}=4$ (the 3, 5, 6d cases are omitted here) |
| $\mathfrak{sl}(4\lvert\mathcal N)$, $\mathfrak{psl}(4\lvert4)$ | complexified 4d superconformal algebra | §2.1, after (7) | for $\mathcal N=4$ one passes to the quotient because $\mathfrak{sl}(k\lvert k)$ has a one-dimensional center |
| $P_{\alpha\dot\alpha},Q^i_\alpha,\bar Q_{\dot\alpha i},M_{\alpha\beta},\Delta,R^i_j,\bar M_{\dot\alpha\dot\beta},S^\alpha_i,\bar S^{\dot\alpha i},K^{\alpha\dot\alpha}$ | generators of $\mathfrak{sl}(4\lvert\mathcal N)$ | diagram (8) | vertical position = conformal weight; horizontal position = (number of chiral) − (number of antichiral) spinor indices; parity = $\mathbb Z$-grading mod 2 |
| $x_{\alpha\dot\alpha}$ | coordinate on the open cell of $\mathrm{Fl}(2;4)$ | Rem. 2.3, (5)–(6) | ${\rm Fl}(2;4)=GL(4,\mathbb{C})/P,\ P=\begin{bmatrix}*&*\\0&*\end{bmatrix}$; coset representative $\begin{bmatrix}1&0\\x_{\alpha\dot{\alpha}}&1\end{bmatrix}$ |
| $V,S_+,R,\Pi$ | vector rep of $\mathrm{Spin}(4)$, chiral spinor, defining rep of $U(\mathcal N)$, parity shift | Prop. 2.4, (9) | $\mathbb{C}^{4\lvert 2\mathcal{N}}\cong V\oplus\Pi(S_{+}\otimes R)$ (chiral superspace) |
| $D_{\alpha i},\bar D^j_{\dot\alpha}$ | used, not defined | (10) | $[X,D_{\alpha i}]\sim D_{\alpha i},\quad[X,\bar{D}_{\dot{\alpha}}^{j}]\sim D_{\dot{\alpha}}^{j}$ (the missing bar on the right is in the source) |
| $\mathcal N=1$ chiral-superspace vector fields | — | (11)–(13) | $P_{\alpha\dot{\alpha}}=\partial_{\alpha\dot{\alpha}}$, $M_{\alpha\beta}=x_{\alpha\dot{\alpha}}\partial_{\beta\dot{\alpha}}+\theta_{\alpha}\partial_{\beta}+(\alpha\leftrightarrow\beta)$, $\bar{M}_{\dot{\alpha}\dot{\beta}}=x_{\alpha\dot{\alpha}}\partial_{\alpha\dot{\beta}}+(\dot{\alpha}\leftrightarrow\dot{\beta})$; $\Delta=E+\frac{1}{2}\theta^{\alpha}\partial_{\alpha}$, $R=\theta^{\alpha}\partial_{\alpha}$, $K^{\alpha\dot{\alpha}}=\lvert x\rvert^{2}\partial^{\alpha\dot{\alpha}}-2x^{\alpha\dot{\alpha}}E+\theta^{\alpha}x^{\beta\dot{\alpha}}\partial_{\beta}$; $Q_{\alpha}=\partial_{\alpha}$, $\bar{Q}_{\dot{\alpha}}=\theta^{\alpha}\partial_{\alpha\dot{\alpha}}$, $S^{\alpha}=-x^{\alpha\dot{\alpha}}\theta^{\beta}\partial_{\beta\dot{\alpha}}+\theta^{2}\partial^{\alpha}$, $\bar{S}^{\dot{\alpha}}=x^{\alpha\dot{\alpha}}\partial_{\alpha}$ |
| $\mathrm{Fl}(m_i\lvert n_i;m\lvert n)$ | super flag variety | Def. 2.6, (14) | $0\subset V_{1}\subset\cdots\subset V_{k}\subset\mathbb{C}^{m\lvert n},\ \dim V_{i}=m_{i}\lvert n_{i}$ |
| chiral-superspace cells | — | Prop. 2.8, (16) | left-chiral superspace is open dense in $\mathrm{Fl}(2\lvert0;4\lvert\mathcal N)$; right-chiral in $\mathrm{Fl}(2\lvert\mathcal N;4\lvert\mathcal N)$; full superspace $\mathbb C^{4\lvert4\mathcal N}$ compactifies to $\mathrm{Fl}(2\lvert0,2\lvert\mathcal N;4\lvert\mathcal N)$; representatives $\left[\begin{smallmatrix}1&0&0\\x^{\alpha\dot\alpha}&1&0\\\theta^\alpha_i&0&0\end{smallmatrix}\right]$, $\left[\begin{smallmatrix}1&0&0\\x^{\alpha\dot\alpha}&1&\bar\theta^{\dot\alpha i}\\0&0&0\end{smallmatrix}\right]$ |
| matrix presentation | $\mathfrak{sl}(4\lvert\mathcal N)$ blocks, ordered $(\alpha,\dot\alpha\,\lvert\,i)$ | Rem. 2.9, (17) | $\left[\begin{array}{cc\lvert c}M^{\alpha}_{\beta}&K^{\alpha\dot{\beta}}&S^{\alpha}_{j}\\P_{\alpha\dot{\beta}}&\bar{M}_{\dot{\alpha}}^{\dot{\beta}}&\bar{Q}_{\dot{\alpha}j}\\\hline Q_{\alpha}^{i}&\bar{S}^{\dot{\alpha}i}&R^{i}_{j}\end{array}\right];\ \Delta=\left[\begin{array}{cc\lvert c}1&0&0\\0&-1&0\\\hline 0&0&0\end{array}\right],\ r=\frac{1}{\mathcal{N}-4}\left[\begin{array}{cc\lvert c}\mathcal{N}&0&0\\0&\mathcal{N}&0\\\hline 0&0&4\end{array}\right]$ |
| $e_{\mu\nu}$; $\mu$ | elementary matrices; $\mu$ is any basis index of $\mathbb C^{4\lvert\mathcal N}$ | Thm. 2.10 proof, (19) | $[e_{0+},e_{\mu\nu}]=\delta_{+\mu}e_{0\nu}-(-)^{\lvert\mu\nu\rvert}\delta_{\nu 0}e_{\mu+}$ |
| $\pm$; $0$, $i$ | $SO(4)\to U(2)$: the left-chiral index becomes $\pm$. $U(\mathcal N)\to U(1)\times U(\mathcal N-1)$: the R-index becomes $0$, $i$ | Thm. 2.10, (18) | holomorphic supercharge $Q^0_+$, identified with $e_{0+}$ |
| $\mathfrak t_{\mathcal N}$, $S_\pm$, $\Gamma$ | complexified supertranslation algebra | §4.1 | $\mathfrak{t}_{\mathcal{N}}=\mathbb{C}^{4}\oplus\Pi(S_{+}\otimes\mathbb{C}^{\mathcal{N}}\oplus S_{-}\otimes\mathbb{C}^{\mathcal{N}})$, $\Gamma:S_{+}\otimes S_{-}\to\mathbb{C}^{4}$ (normalization not given) |
| holomorphic supercharge | odd, square-zero $Q\in\mathfrak t_{\mathcal N}$ with $\mathrm{im}[Q,-]\subset\mathbb C^4$ of rank two | §4.1 | — |
| $\rho$ | twisting homomorphism; $Q$ has weight $+1$ | Rem. 4.4 ($\mathcal N=2$ only) | $\rho:\mathbb{C}^{\times}\to G_{R}$ |
| $\varepsilon_i$, $X^{d\lvert m}$ | odd coordinates of cohomological degree $-1$ (treated algebraically) | §3.1.1 | $\mathcal{O}(X^{d\lvert m})=\mathcal{O}(X)[\varepsilon_{1},\ldots,\varepsilon_{m}]$ |
| $\mathcal G_{\mathcal N}$ | holomorphic currents | Def. 3.1, (25)–(26) | $\mathcal{G}_{\mathcal{N}}^{\bullet}(X)=\Omega^{0,\bullet}(X,\mathfrak{g}\otimes_{\mathbb{C}}A)$, $A=\mathbb{C}[\varepsilon_{1},\ldots,\varepsilon_{\mathcal{N}-1}]$ |
| $\mathfrak X_{\mathcal N}$ | extended holomorphic vector fields | Def. 3.2 | $\mathfrak{X}_{\mathcal{N}}^{\bullet}(X)=\left(\Omega^{0,\bullet}(X,TX)\otimes_{\mathbb{C}}A\right)\bowtie\left(\Omega^{0,\bullet}(X)\otimes_{\mathbb{C}}{\rm Der}(A)\right)$, differential $\bar\partial$ |
| $\underline A,\underline B$ | holomorphic BF fields | Def. 3.4 | $\underline{A}\in\Omega^{0,\bullet}(X,\mathfrak{h})[1]$, $\underline{B}\in\Omega^{2,\bullet}(X,\mathfrak{h}^{\vee})$; $S=\int_X\langle\underline{B}\wedge\bar\partial\underline{A}\rangle_{\mathfrak{h}}+\frac{1}{2}\int_X\langle\underline{B}\wedge[\underline{A},\underline{A}]\rangle_{\mathfrak{h}}$ |
| hCS | holomorphic Chern–Simons | Def. 3.4 | $S(\underline{A})=\frac{1}{2}\int_X\langle\underline{A}\wedge\bar\partial\underline{A}\rangle+\frac{1}{6}\int_X\langle\underline{A}\wedge[\underline{A},\underline{A}]\rangle$ |
| $\mathcal N=4$ twisted fields | — | Prop. 4.5, (29) | $\underline{\underline{A}}\in\Omega^{0,\bullet}(\mathbb{C}^{2},\mathfrak{g})[\varepsilon_{1},\varepsilon_{2},\varepsilon_{3}]\cong\Omega^{0,\bullet}(\mathbb{C}^{2\lvert 3},\mathfrak{g})$; pairing $\int_{\mathbb{C}^{0\lvert 3}}{\rm d}^{3}\varepsilon\,\langle\underline{X}_{1},\underline{X}_{2}\rangle$; ${\rm d}^{2}z={\rm d}z_{1}{\rm d}z_{2}$ |
| $\mathcal N=4$ hCS action | — | §4.1 | $\int_{\mathbb{C}^{2\lvert 3}}{\rm d}^{2}z\,{\rm d}^{3}\varepsilon\,\left(\langle\underline{A}\wedge\bar\partial\underline{A}\rangle+\frac{1}{6}\langle(\underline{A}\wedge[\underline{A},\underline{A}])\rangle\right)$ |
| $\delta$ | degree-one parameter; identified with $\varepsilon_3$ | §4.1 | $\underline{A}+\delta\underline{B}\cong\Omega^{0,\bullet}(\mathbb{C}^{2},\mathfrak{g}[\varepsilon_{1},\varepsilon_{2}])[1][\delta]$ |
| $\underline\partial$, $\mathfrak X^{div}_4$ | super-divergence; divergence-free vector fields on $\mathbb C^{2\lvert3}$ | Def. 4.8 | $\underline{\partial}=\sum_{i=1,2}\frac{\partial}{\partial z_{i}}\frac{\partial}{\partial(\partial_{z_{i}})}-\sum_{a=1,2,3}\frac{\partial}{\partial\varepsilon_{a}}\frac{\partial}{\partial(\partial_{\varepsilon_{a}})}$, complex $\mathfrak X_4\xrightarrow{\underline\partial}\Omega^{0,\bullet}(\mathbb C^{2\lvert3})$ in degrees $0,1$ |
| embedding | $\mathfrak X_3\hookrightarrow\mathfrak X_4^{div}$ | after Prop. 4.9 | $\underline{\zeta}\mapsto\underline{\zeta}+(\underline{\partial}\zeta)\varepsilon_{3}\frac{\partial}{\partial\varepsilon_{3}}$ |
| $I(\underline\xi)$ | local Noether current | (27); Prop. 4.9 | $I(\underline{\xi})=\int_X\langle\underline{B}\wedge L_{\underline{\xi}}\underline{A}\rangle$; for $\mathcal N=4$: $\frac{1}{2}\int_{\mathbb{C}^{2\lvert 3}}{\rm d}^{2}z{\rm d}^{3}\varepsilon\,\langle\underline{A}L_{\underline{\xi}}\underline{A}\rangle$ |
| $\mathcal O_{\rm loc}$, $\{-,-\}$, ${\rm d}_{CE}$ | local functionals, BV bracket, Chevalley–Eilenberg differential | §3.3 | ${\rm d}_{CE}I+\{S,I\}+\frac{1}{2}\{I,I\}=0$ |
| $Q^{GF}$, $P_{\epsilon<L}$, $K_\epsilon$, $W$, $I[L]$, $\Delta_L$, $\{,\}_L$, $\Theta$ | gauge fixing, propagator, heat kernel, graph weight, effective action | §7.2–7.3 (for $\mathcal N=2$) | $Q^{GF}=\bar\partial^{*}\otimes 1$; $I[L]+I_{S}[L]\overset{\rm def}{=}\lim_{\epsilon\to 0}W(P_{\epsilon<L},I+I_{S})$ |

## KEY STATEMENTS

All of these are paraphrases; the formulas are exact.

- **(7) and following text.** The complexified algebra is $\mathfrak{sl}(4\lvert\mathcal N)$, or $\mathfrak{psl}(4\lvert4)$ when $\mathcal N=4$.
  - $Q^i_\alpha$ and $\bar S^{\dot\alpha i}$ together form a chiral spinor of $\mathcal O(4,2)$, equivalent to the antifundamental of $\mathfrak{su}(2,2)$. They are in the fundamental of the R-symmetry.
  - $\bar Q_{\dot\alpha i}$ and $S^\alpha_i$ are in the fundamental of $\mathfrak{su}(2,2)$ and the antifundamental of the R-symmetry.
  - The source says it generally follows [DO09] (Dolan–Osborn).
- **Prop. 2.4.** $\mathfrak{sl}(4\lvert\mathcal N)$ acts by supervector fields on $\mathbb C^{4\lvert2\mathcal N}$. For $\mathcal N=4$ the action factors through $\mathfrak{psl}(4\lvert4)$. Explicit vector fields are given only for $\mathcal N=1$, in (11)–(13).
- **Thm. 2.10.** The cohomology of $\mathfrak c(\mathcal N)$ with respect to a holomorphic supercharge is $\mathfrak{sl}(3\lvert\mathcal N-1)$ for $\mathcal N=1,2$, and $\mathfrak{psl}(3\lvert3)$ for $\mathcal N=4$.
  - Proof: with $Q=e_{0+}$, $\ker(\mathrm{ad}_Q)$ is spanned by the $e_{\mu\nu}$ with $\mu\neq+$, $\nu\neq0$, together with $e_{00}+e_{++}$. $\mathrm{im}(\mathrm{ad}_Q)$ is spanned by the $e_{0\nu}$ and $e_{\mu+}$.
  - Survivors, as listed in the text: $P^-_{\dot\alpha}$, $K^{\dot\alpha}_+$, $\bar M_{\dot\alpha\dot\beta}$, traceless $R_i^j$, $\bar Q^i_{\dot\alpha}$, $Q^i_-$, $\bar S^{\dot\alpha i}$, $S^-_i$, plus two central generators.
- **Cor. 2.11, (21)–(24).** $\mathfrak{sl}(3\lvert\mathcal N-1)$ acts on $\mathbb C^{2\lvert\mathcal N-1}$ by holomorphic vector fields:
  - $p_{i}=\frac{\partial}{\partial z_i},\ m_{ij}=z_{i}\frac{\partial}{\partial z_j},\ k_{i}=z_{i}e$, with $e=z_i\partial/\partial z_i=\mathrm{tr}(m)$.
  - Weight generator: $z_{i}\frac{\partial}{\partial z_i}+\frac{1}{2}\varepsilon\frac{\partial}{\partial\varepsilon}$.
  - Odd generators: $\frac{\partial}{\partial\varepsilon},\ \varepsilon\frac{\partial}{\partial z_i}$ (positive weight) and $z_i\frac{\partial}{\partial\varepsilon},\ \varepsilon e$ (negative weight). For $\mathcal N=2$, $\hat e=z_i\partial/\partial z_i+\varepsilon\partial/\partial\varepsilon$ replaces $e$ in $k_i$.
  - The R-symmetry acts by $\varepsilon_i\partial/\partial\varepsilon_j$.
- **Prop. 4.5.** The holomorphic twist of $\mathcal N=4$ SYM with ordinary $\mathfrak g$ is equivalent to holomorphic BF theory with values in $\mathfrak g[\varepsilon_1,\varepsilon_2]$. If $\mathfrak g$ has a non-degenerate invariant pairing, it is also equivalent to hCS with values in $\mathfrak g[\varepsilon_1,\varepsilon_2,\varepsilon_3]$. The source cites [Bau11, Cos13, ES19, ESW20] for this.
- **Prop. 4.9.** For any holomorphic $Q$, twisted $\mathcal N=4$ SYM has a classical symmetry by $\mathfrak X_4^{div}$. Its Noether current is given in the table above.
- **§3.3.** The Noether map is one quote: it "defines a map of Lie algebras from $\mathfrak{X}_{\mathcal{N}}$ to the space of local functionals". This is equivalent to the MC equation listed in the table.
- **§6.5.** $\mathfrak X_3$ is a symmetry of twisted $\mathcal N=4$ SYM. $z_2\partial_{\varepsilon_1}$ is divergence-free. Localization to the small $\mathcal N=4$ vertex algebra is only expected; no proof is given.
- **Off-shell closure.**
  - For the untwisted $\mathcal N=4$ SYM, the source says nothing about whether the SUSY or superconformal action closes off shell.
  - Rem. 7.10 is a generic remark: when SUSY is off-shell via auxiliary fields, twisting can be done at the BRST level. After auxiliaries are eliminated, the twisting supercharge may depend on interaction terms.
  - On the twisted theory, the symmetry is a strict Lie map into $\mathcal O_{\rm loc}[-1]$, at the classical level only.

## SIGN CONVENTIONS

- **Conformal brackets, (4), verbatim.**
  - $[D,P_\mu]=P_\mu$, $[D,K_\mu]=-K_\mu$.
  - $[M_{\mu\nu},K_\rho]=g_{\rho\nu}K_\mu-g_{\mu\rho}K_\nu$, $[M_{\mu\nu},P_\rho]=g_{\rho\nu}P_\mu-g_{\mu\rho}P_\nu$.
  - $[K_\mu,P_\nu]=2M_{\mu\nu}+2g_{\mu\nu}E$.
  - $[M_{\mu\nu},M_{\rho\sigma}]=g_{\mu\sigma}M_{\nu\rho}+g_{\nu\rho}M_{\mu\sigma}-g_{\mu\rho}M_{\nu\sigma}-g_{\nu\sigma}M_{\mu\rho}$.
  - **[derived]** These are the Lie brackets of the vector fields (3) with $D=\Delta=-E$. The $E$ in $[K,P]$ equals $-\Delta$.
- **Dilatation normalization is not consistent across the paper.**
  - (3)/(4): $\mathrm{ad}_\Delta P=+P$.
  - (12): $\Delta=E+\tfrac12\theta\partial$. **[derived]** This gives $[\Delta,P]=-P$ and $[\Delta,Q_\alpha]=-\tfrac12Q_\alpha$.
  - (17): **[derived]** $\mathrm{ad}_\Delta e_{rc}=(d_r-d_c)e_{rc}$ gives weights $P:-2$, $K:+2$, $Q,\bar Q:-1$, $S,\bar S:+1$.
- **$\mathcal N=1$ supersymmetry bracket.** **[derived]** From (11) and (13), $\{Q_\alpha,\bar Q_{\dot\alpha}\}=\partial_{\alpha\dot\alpha}=P_{\alpha\dot\alpha}$, with unit coefficient and no $i$. No $\mathcal N$-general $\{Q,\bar Q\}$, $\{Q,S\}$ or $\{S,\bar S\}$ is written; only the $e_{\mu\nu}$ relation (19) is given.
- **R-symmetry.** **[derived]** $\mathrm{str}\,r=0$, and $[r,Q]=-Q$, $[r,\bar S]=-\bar S$, $[r,\bar Q]=+\bar Q$, $[r,S]=+S$. At $\mathcal N=4$, $r$ is singular, $\mathbb 1\in\mathfrak{sl}(4\lvert4)$ is central, and that centre is quotiented out.
- **Degrees and BV.**
  - $\varepsilon_i$ has degree $-1$. With the other twisting map $\rho_2$ ($\mathcal N=2$), $\varepsilon$ would have degree $+1$ instead.
  - The field space carries a $(-1)$-shifted symplectic pairing.
  - $\underline A\in\Omega^{0,\bullet}[1]$ and $\underline B\in\Omega^{2,\bullet}$.
  - Pairings are graded symmetric (BF) or graded skew-symmetric (hCS).
- **QME (§7.2).** $\bar\partial I[L]+\frac{1}{2}\{I[L],I[L]\}_{L}+\hbar\Delta_{L}I[L]=0$, so $+\hbar\Delta$.
- **Anomaly.**
  - $\Theta[L]=\hbar^{-1}\big(\bar\partial(I[L]+I_{S}[L])+\frac{1}{2}\{I[L]+I_{S}[L],I[L]+I_{S}[L]\}_{L}+\hbar\Delta_{L}(I[L]+I_{S}[L])\big)$.
  - $\hbar\Theta=\lim_{L\to 0}\lim_{\epsilon\to 0}\sum_{\Gamma\in{\rm Wheel}_{3},e}W_{\Gamma,e}(P_{\epsilon<L},K_{\epsilon},I+I_{S})$.
- **Counterterms.** None. The holomorphic theory is one-loop finite (Lemma 7.5): $\lim_{\epsilon\to0}W(P_{\epsilon<L},\cdot)$ exists.
- **Not stated in the source:** the exponent weight $e^{\pm S/\hbar}$, the definition of $W$ (deferred to [Wil18], [Cos11]), and the sign of $I^{CT}$.
- **$\hbar$.** Set to $\hbar=2\pi i$ for vertex algebras (Rem. 7.12).
- **Notation clash.** `\varepsilon` is an odd coordinate; `\epsilon` is the UV cutoff.
- **Regularization.** No dimensional regularization is used; the regulator is a propagator cutoff with a heat kernel.

## GAPS

1. There is no BV action of untwisted 4d $\mathcal N=4$ SYM in this source: no fields, no Weyl-spinor component action, no Euclidean real form, and no explicit $\mathfrak{so}(5,1)$.
2. There is no action of $\mathfrak{psl}(4\lvert4)$ on the $\mathcal N=4$ fields or BV complex. Explicit supervector fields exist only for $\mathcal N=1$ chiral superspace, (11)–(13).
3. There are no ghosts for the superconformal algebra and no ghost notation. Only ${\rm d}_{CE}$ and the Noether map $I$ appear, and only for the twisted $\mathfrak X_{\mathcal N}$.
4. Conformal Killing data: the vector fields (3) are given, but the conformal Killing equation and the Weyl-weight action on fields are not written.
5. No full set of $\mathfrak{psl}(4\lvert4)$ (anti)commutators with normalizations is given. The dilatation normalization conflicts between (3)/(4), (12) and (17).
6. There is no statement on off-shell closure for $\mathcal N=4$, and no $L_\infty$ corrections for the untwisted theory.
7. The renormalization uses $P_{\epsilon<L}$ with finite $L$, then $L\to0$ for the anomaly. Propagators $P(\epsilon,\infty)$, composite operators $O(\epsilon)$, counterterms $I^{CT}(\epsilon)$ and two-loop terms are absent. The theory is one-loop exact and holomorphic only.
8. The exponent convention and the definition of $W(P,I)$ are deferred to [Wil18] (Williams, *Renormalization for holomorphic field theories*, 2018) and [Cos11] (Costello, *Renormalization and effective field theory*). The $\mathcal N=4$ twist is deferred to [Cos13] (Costello, *Notes on supersymmetric and holomorphic field theories in dimensions 2 and 4*) and [ESW20] (Elliott–Safronov–Williams, *A Taxonomy of Twists of Supersymmetric Yang–Mills Theory*). These are citations listed in the source; I did not read them.
9. Possible source typos: in Def. 3.2, $[X\otimes a,X'\otimes a']=[X,Y]\otimes aa'$ has $Y$ where $X'$ should be. Prop. 4.5's hCS action has no $\tfrac12$ on the kinetic term, while Def. 3.4 has one. In (10) the right-hand side is unbarred.