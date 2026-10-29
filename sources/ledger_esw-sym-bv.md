# NOTATION LEDGER: Elliott–Safronov–Williams, "A Taxonomy of Twists of Supersymmetric Yang–Mills Theory"

## SOURCE

- C. Elliott, P. Safronov, B. R. Williams, *A Taxonomy of Twists of Supersymmetric Yang–Mills Theory*, arXiv:2002.10517. Version v1 (24 Feb 2020) is the only version on the abstract page. arXiv DOI 10.48550/arXiv.2002.10517. No journal reference is listed.
- What was read: https://arxiv.org/pdf/2002.10517v1 (98 pp.), fetched with WebFetch. I checked it by text extraction and by rendered page images. `arxiv.org/html/2002.10517v1` comes back cut off before §2.1, and ar5iv failed. Page numbers below are the printed page numbers, which equal the PDF page numbers.
- Supplementary reading, limited to the item noted in GAPS: Saberi–Williams [SW19a], arXiv:1910.04120v1, §3.1 p.13. ESW cite [EWY18] (arXiv:1702.05973) only for non-SUSY YM asymptotic freedom (p.6), not for N=4 components, so I did not read it.

## SYMBOLS

| symbol | meaning | location | verbatim defining formula / statement |
|---|---|---|---|
| $V_{\mathbb R},\ V$ | real spacetime vector space and its complexification | §2.2 p.36 | "Let $V_{\mathbb R}$ be a finite-dimensional real vector space of dimension $n=\dim_{\mathbb R}(V_{\mathbb R})$ equipped with an orientation and a nondegenerate symmetric bilinear form. Denote by $V=V_{\mathbb R}\otimes_{\mathbb R}\mathbb C$ its complexification." |
| $\mathrm{Cl}(V)$ | Clifford algebra | §2.1 p.35 | $v_1v_2+v_2v_1=2(v_1,v_2)$ |
| $M=M^+\oplus M^-,\ \rho(v)$ | $\mathbb Z/2$-graded Clifford module and Clifford action | §2.1 p.35 | $(-,-)\colon M^+\otimes M^-\to\mathbb C$, $\ (\rho(v)Q_1,Q_2)=(Q_1,\rho(v)Q_2)$; "From now on we denote $M^+=\Sigma$ and $M^-=\Sigma^*$." |
| $\Gamma$ | $\Gamma$-pairings | (28) p.35 | $\Gamma\colon \mathrm{Sym}^2(\Sigma)\to V,\ \Gamma\colon\mathrm{Sym}^2(\Sigma^*)\to V$ by $(\rho(v)Q_1,Q_2)=(v,\Gamma(Q_1,Q_2))$ |
| $\wedge^2V\to\mathfrak{so}(V)$ | two-forms acting on spinors | §2.1 p.35 | $\omega\mapsto(w\mapsto-2\iota_{(w,-)}\omega)$ |
| $q$ | antisymmetrization map | (29) p.35 | $q(v_1\wedge v_2)=v_1v_2-(v_1,v_2)$ |
| $\not d,\ \not d_A$ | Dirac operator, and twisted Dirac operator "obtained from $\Gamma$" | §2.1 p.35, §3.1 p.44 | $\not d\colon\Gamma(N,\Sigma)\to\Gamma(N,\Sigma^*)$ |
| $S;\ S_\pm$ | spin rep ($n$ odd); semi-spin reps ($n$ even) | §2.2 p.36 | — |
| $\Sigma$ | spinorial representation | Def 2.4 p.36 | "a finite sum of spin or semi-spin representations of $\mathfrak{so}(V)$"; "$\Sigma=S\otimes W$ … $\Sigma=S_+\otimes W_+\oplus S_-\otimes W_-$, where $W$ is a finite dimensional multiplicity space" |
| $\mathcal N$ | number of supersymmetries | §2.2 pp.36–37 | "If $n\equiv0,1,3,4\pmod 8$, we let $\mathcal N=\dim(W)$." |
| $\mathfrak A$ | supertranslation Lie algebra | Def 2.5 p.36 | "the $\mathfrak{so}(V)$-equivariant super Lie algebra $\mathfrak A=\Pi\Sigma\oplus V$ whose only nontrivial bracket is given by $\Gamma$" |
| $G_R,\ \mathfrak g_R$ | R-symmetry group and its Lie algebra | §2.2 p.37 | "A Lie group $G_R$, the group of R-symmetries, which acts on $\Sigma$ by $\mathfrak{so}(V)$-equivariant automorphisms preserving $\Gamma$." |
| supersymmetry algebra | — | Def 2.6 p.37 | $(\mathfrak{so}(V)\oplus\mathfrak g_R)\ltimes\mathfrak A$ |
| $ISO(V_{\mathbb R})$ | Poincaré group | p.37 | $ISO(V_{\mathbb R})=\mathrm{Spin}(V_{\mathbb R})\ltimes V_{\mathbb R}$ |
| 4d $\Sigma$ | odd part in 4d | §10 p.68 | $\Sigma\cong S_+\otimes W\oplus S_-\otimes W^*$; "$S_+,S_-$ are the 2-dimensional semi-spin representations of $\mathrm{Spin}(4,\mathbb C)\cong SL(2,\mathbb C)\times SL(2,\mathbb C)$ … The semi-spin representations carry symplectic pairings $S_\pm\otimes S_\pm\to\mathbb C$." "$\mathcal N=\dim(W)=1,2,4$" |
| N=4 $W$, $G_R$ | 4d N=4 R-symmetry | §10.3 p.75 | "$G_R=\mathrm{Spin}(6;\mathbb C)\cong SL(4,\mathbb C)$ under which $W\cong S_+^6$ is the positive six-dimensional semi-spin representation and $W^*\cong S_-^6$." |
| $\mathrm{Spin}(4,\mathbb R)$ | 4d rotation group | §10.3 p.75 | $\mathrm{Spin}(4,\mathbb R)\cong SU(2)_+\times SU(2)_-$ |
| $W_{\mathbb R},W,\Gamma_W,G_W$ | data for dimensional reduction. **Name clash:** here $W$ is a spacetime subspace, not the multiplicity space | §2.6 p.42 | $V_{\mathbb R}=W_{\mathbb R}\oplus W_{\mathbb R}^\perp$, $W=W_{\mathbb R}\otimes_{\mathbb R}\mathbb C$; $\Gamma_W\colon\mathrm{Sym}^2(\Sigma)\xrightarrow{\Gamma_V}V\to W$; $G_W=G_V\times\mathrm{Spin}(W_{\mathbb R}^\perp)$ |
| $R$ | **a $\mathfrak g$-representation (chiral matter)**, not an R-symmetry rep | §3 p.43, §10.1 p.68 | "Let $R$ be a complex $\mathfrak g$-representation." |
| $\mathfrak g$ | gauge Lie algebra | §3.1 p.44; Part II p.51 | "a Lie algebra equipped with a nondegenerate symmetric bilinear pairing"; Part II: "a complex Lie algebra $\mathfrak g$ … the complexified Lie algebra of the gauge group" |
| $\Pi,\ \lvert x\rvert$ | parity shift; total degree | §1.1 p.17 | "We will write $\Pi E$ to denote $E$ placed in odd $\mathbb Z/2\mathbb Z$ degree. For an element $x$ we denote by $\lvert x\rvert\in\mathbb Z/2\mathbb Z$ the total degree." |
| $E,\ \mathcal E,\ \mathcal E_c,\ E^!$ | bundle; smooth sections; compactly supported sections; $E^*\otimes\mathrm{Dens}_M$ | §1.1 p.17; Intro p.3 | $E^!=E^*\otimes\mathrm{Dens}_M$; "$E\cong E^![-1]$" |
| $\mathcal O_{\rm loc}(\mathcal E),\ \mathcal O^+_{\rm loc}(\mathcal E)$ | local functionals; those that are at least cubic | §1.1 pp.17–18 | $\mathcal O_{\rm loc}(\mathcal E)\subset\mathcal O(\mathcal E_c)$ |
| $\omega$ | BV pairing | Def 1.13 p.20; Def 1.16 p.22 | $\omega\colon L\otimes L\to\mathrm{Dens}_M[n-2]$ (strict $n$-shifted); free BV theory: $\omega\colon E\otimes E\to\mathrm{Dens}_M[-1]$ |
| $X_F,\ \{-,-\}$ | Hamiltonian vector field; BV bracket | Def 1.15 p.21 | $X_F:=\omega\circ d_{\rm dR}(F)\in\mathrm{Fun}_{\rm loc}(\mathcal E,\mathcal E[-n])$; $\{F,G\}=X_F(G)$ |
| $Q_{\rm BV},\ S,\ I$ | BV differential; action; interaction | p.21 | $S=\frac12\int_M\omega(e,Q_{BV}e)+I$ |
| $A,\lambda,c$ | BRST fields | §3.1 p.44 | $A\in\Omega^1(V_{\mathbb R};\mathfrak g)$, $\lambda\in\Gamma(V_{\mathbb R};\Pi\Sigma\otimes\mathfrak g)$, $c\in\Omega^0(V_{\mathbb R};\mathfrak g[1])$ |
| $F_A$ | curvature | §3.1 p.44 | $F_A=\mathrm dA+\frac12[A\wedge A]$ |
| $A^*,\lambda^*,c^*$ | antifields | §3.1 p.44; Rem 1.36 p.28 | "$T^*[-1]F_{\rm gauge}=F_{\rm gauge}\oplus F^!_{\rm gauge}[-1]$"; "Given a field $\phi$ the corresponding antifield will be denoted by $\phi^*$." |
| $C^\bullet(\mathfrak g,\mathcal O_{\rm loc}(\mathcal E)),\ d_{\rm CE}$ | equivariant (Chevalley–Eilenberg) complex | Def 1.26 p.25 | see KEY STATEMENTS |
| $S_{\mathfrak g}^{(k)},\ S_{\mathfrak A}^{(k)}$ | $L_\infty$-action components | Def 1.27 p.26; §2.5 p.41 | $S_{\mathfrak A}=S+\sum_{k\ge1}S^{(k)}_{\mathfrak A}\in C^\bullet(\mathfrak A,\mathcal O_{\rm loc}(\mathcal E))$ |
| $Q$ (square-zero) | square-zero supercharge | Def 2.16 p.41 | "a nonzero element $Q\in\Sigma$ such that $\Gamma(Q,Q)=0$" |
| $\alpha,\ \phi$ | compatible $U(1)$; twisting homomorphism | Defs 2.18–2.19 p.42 | — |
| $\eta$ | supergravity "ghost for supertranslations" (Intro only) | p.6 | $\eta\in\Gamma(M,(P_{\rm Spin}\times P_R)\times^{\mathrm{Spin}(n)\times G_R}\Sigma)$ |

There are no 2-component index notations ($\alpha,\dot\alpha$), no $\sigma^\mu$ or $\gamma^\mu$ matrices, and no symbol for an R-symmetry index anywhere in ESW.

## KEY STATEMENTS (verbatim)

**Def 1.14 (p.20).** "A classical BV field theory (or, simply, classical field theory) is a local formal moduli problem on the spacetime manifold $M$ equipped with a strict $(-1)$-shifted symplectic structure."

**Def 1.15 and the text after it (p.21).** $\{F,G\}=X_F(G)$. "The Poisson bracket enjoys the graded skew symmetry property $\{F,G\}=(-1)^{|F||G|}\{G,F\}$ as well as the graded Jacobi identity." Also: $\{-,-\}\colon\mathcal O_{\rm loc}(\mathcal M)\times\mathcal O_{\rm loc}(\mathcal M)\to\mathcal O_{\rm loc}(\mathcal M)[-n]$.

**Action and CME (p.21).** "…the action functional $S\in\mathcal O_{\rm loc}(BL)$ of cohomological degree $1-n$ such that $S=\frac12\int_M\omega(e,Q_{BV}e)+I$, where $e\in L$ and $I$ … is at least cubic. … classical master equation $\{S,S\}=0$."

**Def 1.16 (p.22).** Free BV theory:
- $Q_{BV}\colon E\to E[1]$, "even … of cohomological degree +1", with $Q_{BV}^2=0$ and elliptic.
- $\omega\colon E\otimes E\to\mathrm{Dens}_M[-1]$, fiberwise nondegenerate, graded skew symmetric, with $\int_M\omega(e_0,Q_{BV}e_1)=(-1)^{|e_0|}\int_M\omega(Q_{BV}e_0,e_1)$.
- Interaction: "$I\in\mathcal O^+_{\rm loc}(\mathcal E)$ of cohomological degree zero satisfying the Maurer-Cartan equation $Q_{BV}I+\frac12\{I,I\}=0$".

**Def 1.26 (p.25).**
- $C^\bullet(\mathfrak g,\mathcal O_{\rm loc}(\mathcal E))=\bigoplus_n\mathrm{Hom}(\wedge^n\mathfrak g,\mathcal O_{\rm loc}(\mathcal E))[-n]$.
- Antisymmetry: $f(x_1,\dots,x_i,x_{i+1},\dots,x_n)=(-1)^{|x_1||x_2|+1}f(x_1,\dots,x_{i+1},x_i,\dots,x_n)$.
- Differential, "following the sign conventions of [Saf17]": $(d_{\rm CE}f)(x_1,\dots,x_n)=\sum_{i<j}(-1)^{|x_i|\sum_{p=1}^{i-1}|x_p|+|x_j|\sum_{p=1,p\ne i}^{j-1}|x_p|+i+j+|f|}f([x_i,x_j],x_1,\dots,\widehat{x_i},\dots,\widehat{x_j},\dots,x_n)$.
- "degree +1 BV bracket": $\{f,g\}(x_1,\dots,x_{k+l})=\sum_{\sigma\in S_{k,l}}\mathrm{sgn}(\sigma)(-1)^{\epsilon+\epsilon_1}\{f(x_{\sigma(1)},\dots,x_{\sigma(k)}),g(x_{\sigma(k+1)},\dots,x_{\sigma(k+l)})\}$, where $\epsilon_1=|g|k+\sum_{i=1}^k|x_{\sigma(i)}|(l+|g|)$ and "$\epsilon$ is the usual Koszul sign".
- p.26: $(Q_{BV}f)(x_1,\dots,x_n)=Q_{BV}f(x_1,\dots,x_n)$, $(d_{\rm CE}+Q_{BV})^2=0$.

**Def 1.27 (p.26).** "An action of a super Lie algebra $\mathfrak g$ on $(E,\omega,S)$ is an element $S_{\mathfrak g}=\sum_{k\ge0}S^{(k)}_{\mathfrak g}\in C^\bullet(\mathfrak g,\mathcal O_{\rm loc}(\mathcal E))$ of cohomological degree zero", subject to three conditions:
- (a) $S^{(0)}_{\mathfrak g}=S$.
- (b) For $k\ge1$, $S^{(k)}_{\mathfrak g}(x_1,\dots,x_k)$ "is at least quadratic in the fields".
- (c) $d_{\rm CE}S_{\mathfrak g}+\frac12\{S_{\mathfrak g},S_{\mathfrak g}\}=0$.

**Rem 1.28 (p.26).** "the operator $d_{\rm CE}+\{S_{\mathfrak g},-\}$ defines a differential on the graded vector space $\mathcal O(\mathfrak g[1]\oplus\mathcal E)$ … maps $\mathfrak g^{\otimes k}\otimes\mathcal E^{\otimes\ell}\to\mathcal E$ combining to give $\mathcal E$ the structure of an $L_\infty$-module for $\mathfrak g$."

**Rem 1.29 (p.26).** $I_{\mathfrak g}=I+\sum_{k\ge1}S^{(k)}_{\mathfrak g}$, with $(d_{\rm CE}+Q_{BV})I_{\mathfrak g}+\frac12\{I_{\mathfrak g},I_{\mathfrak g}\}=0$.

**Rem 1.33 (p.27).** "the bundle $L[1]$ is concentrated in $\mathbb Z$-degrees −1 and 0. In this case, sections in degree 0 are thought of as physical fields, and ghosts – sections in degree −1 – are thought of as generators of the infinitesimal gauge symmetry."

**Def 1.34 (p.27).** The associated BV theory is $T^*[-1]\mathcal M$ with action $S_{\rm BRST}+S_{\rm anti}$.

**Def 2.7 (p.37).** "A classical field theory $(\mathcal E,S,\omega)$ is supersymmetric if $E\to M$ is an $ISO(V_{\mathbb R})\times G_R$-equivariant vector bundle and the infinitesimal strict action of the translation Lie algebra $V$ on the classical theory is extended to a $\mathrm{Spin}(V_{\mathbb R})\times G_R$-equivariant $L_\infty$ action of the supertranslation Lie algebra $\mathfrak A$ on the classical theory."

**Prop 2.1 (p.35).** $\Gamma(Q_1,\rho(X)Q_2)+\Gamma(Q_2,\rho(X)Q_1)=-2\iota_{\Gamma(Q_1,Q_2)}X$ for $X\in\wedge^2(V)$.

**Prop 2.2 (p.36).** $\not d\rho(\Gamma(Q_1,\lambda))Q_2=\rho(d\Gamma(Q_1,\lambda))Q_2+(Q_1,\not d\lambda)Q_2$.

**Prop 2.3 (p.36).** $\not d(\rho(X)Q)=\rho(\mathrm dX)Q+(-1)^{n(1+p)}\rho(*\mathrm d*X)Q$ for $X\in\Omega^p(N)$.

**Assumption 3.1 (p.43).** $\rho(\Gamma(Q_1,Q_2))Q_3+\rho(\Gamma(Q_2,Q_3))Q_1+\rho(\Gamma(Q_3,Q_1))Q_2=0$. The relevant cases are listed:
- "(4d $\mathcal N=1$) $\dim(V)=4$ and $\Sigma=S_+\oplus S_-$"
- "(10d $\mathcal N=(1,0)$) $\dim(V)=10$ and $\Sigma=S_+$"

**Def 3.2 (p.44).**
- Bundle: $F_{\rm gauge}=\Omega^1(V_{\mathbb R};\mathfrak g)\oplus\Gamma(V_{\mathbb R};\Pi\Sigma\otimes\mathfrak g)\oplus\Omega^0(V_{\mathbb R};\mathfrak g[1])$.
- Bracket: $[c,A+\lambda+c']=[c,A]+[c,\lambda]+[c,c']$.
- BRST action: $S_{\rm BRST}(A,\lambda)=\int_{V_{\mathbb R}}\mathrm{dvol}\left(-\frac14(F_A,F_A)+\frac12(\lambda,\not d_A\lambda)\right)$.

**(31) p.44, full BV action.**
$$S_{\rm gauge}=\int_{V_{\mathbb R}}\mathrm{dvol}\left(-\frac14(F_A,F_A)+\frac12(\lambda,\not d_A\lambda)+(\mathrm d_Ac,A^*)+([c,\lambda],\lambda^*)+\frac12([c,c],c^*)\right).$$
"To simplify the notation, the pairing on $\mathfrak g$ from now on will be implicit."

**(32) p.44, Poincaré Hamiltonian.** $S^{(1)}_{\rm gauge}(v)=\int_{V_{\mathbb R}}\mathrm{dvol}\,((L_vA,A^*)-(v.\lambda,\lambda^*)-(v.c)c^*)$ for $v\in\mathfrak{iso}(V)$, "where $v.\lambda$ contains both a derivative and the $\mathfrak{so}(V)$ action on $\Sigma$."

**(33)–(34) p.45, supersymmetry action.**
$$S^{(1)}_{\rm gauge}(Q)=\int_{V_{\mathbb R}}\mathrm{dvol}\left(-(\Gamma(Q,\lambda),A^*)+\frac12(\rho(F_A)Q,\lambda^*)\right)$$
$$S^{(2)}_{\rm gauge}(Q_1,Q_2)=\int_{V_{\mathbb R}}\mathrm{dvol}\left(\frac14(\Gamma(Q_1,Q_2),\Gamma(\lambda^*,\lambda^*))-\frac12(Q_1,\lambda^*)(Q_2,\lambda^*)-\iota_{\Gamma(Q_1,Q_2)}Ac^*\right).$$
$S^{(2)}$ is the non-closure term: it is quadratic in $\lambda^*$, and it contains a field-dependent gauge transformation through the $c^*$ term.

**Thm 3.3 (p.45).** "$S_{\rm gauge,\mathfrak A}=S_{\rm gauge}+S^{(1)}_{\rm gauge}+S^{(2)}_{\rm gauge}\in C^\bullet(\mathfrak A,\mathcal O_{\rm loc}(\mathcal E_{\rm gauge}))$ satisfies the classical master equation $d_{\rm CE}(S_{\rm gauge,\mathfrak A})+\frac12\{S_{\rm gauge,\mathfrak A},S_{\rm gauge,\mathfrak A}\}=0$." It decomposes into four equations:
- $\{S_{\rm gauge},S^{(1)}_{\rm gauge}\}=0$
- $\{S_{\rm gauge},S^{(2)}_{\rm gauge}\}+d_{\rm CE}S^{(1)}_{\rm gauge}+\frac12\{S^{(1)}_{\rm gauge},S^{(1)}_{\rm gauge}\}=0$
- $d_{\rm CE}S^{(2)}_{\rm gauge}+\{S^{(1)}_{\rm gauge},S^{(2)}_{\rm gauge}\}=0$
- $\{S^{(2)}_{\rm gauge},S^{(2)}_{\rm gauge}\}=0$

**Lemma 3.6 proof, item 2 (p.46).** $(d_{\rm CE}S^{(1)}_{\rm gauge})(Q_1,Q_2)=S^{(1)}_{\rm gauge}(\Gamma(Q_1,Q_2))$.

**Defs 2.16–2.17 (p.41).**
- Twisted action: $S^Q=S+\sum_{k\ge1}S^{(k)}_{\mathfrak A}(Q,\dots,Q)$. The Q-twisted theory is "$\mathbb Z/2\mathbb Z$-graded".
- Adjectives: "topological if $d=n$", "holomorphic if $n$ is even and $d=n/2$".

**After Def 2.18 (p.42).** "the element $S^{(k)}_{\mathfrak A}(Q,\dots,Q)$ has α-weight $k$. But it also has cohomological degree $-k$."

**Prop 2.20 (p.42).** The dimensional reduction "has a unique $\mathfrak A'$-supersymmetric structure, compatible with the supersymmetry on $V_{\mathbb R}$ in the sense that $p^*S^{(i)}_{V_{\mathbb R}}=S^{(i)}_{W_{\mathbb R}}$."

**4d N=4 (§10.3 p.75).** "The 4d $\mathcal N=4$ super Yang–Mills theory is obtained by dimensional reduction from the 10d $\mathcal N=(1,0)$ super Yang–Mills."
- The only route the source gives to N=4 components is the proof of Thm 10.9 (p.76): "The 4d $\mathcal N=4$ super Yang–Mills theory viewed as a $\mathcal N=2$ theory is 4d $\mathcal N=2$ Yang–Mills theory with matter valued in $U=T^*\mathfrak g$."
- The N=2 BRST fields (§10.2 p.71) are:
  - $A\in\Omega^1(M;\mathfrak g)$
  - $a,\widetilde a\in\Omega^0(M;\mathfrak g)$
  - $(\lambda_+,\lambda_-)\in\Omega^0(M;\Pi S_+\otimes W\otimes\mathfrak g\oplus\Pi S_-\otimes W^*\otimes\mathfrak g)$
  - $\phi\in\Omega^0(M;U\otimes W)$
  - $(\psi_-,\psi_+)\in\Omega^0(M;\Pi S_+\otimes U\oplus\Pi S_-\otimes U)$
  - $c\in\Omega^0(M;\mathfrak g)[1]$

  Here $W$ is 2-dimensional, $G_R=SL(2;\mathbb C)\times GL(1,\mathbb C)$, and the $GL(1,\mathbb C)$ weights are "2 on $a$, −2 on $\widetilde a$ and ±1 on $\lambda_\pm,\psi_\pm$".
- The 4d N=1 analogue (§10.1 p.69): $(\lambda_+,\lambda_-)\in\Omega^0(M;\Pi S_+\otimes\mathfrak g\oplus\Pi S_-\otimes\mathfrak g)$, $(\overline\phi,\phi)\in\Omega^0(M;R\oplus R^*)$, $(\psi_-,\psi_+)\in\Omega^0(M;\Pi S_+\otimes R^*\oplus\Pi S_-\otimes R)$, ghost "$A_0\in\Omega^0(M;\mathfrak g)[1]$".
- Kapustin–Witten (§10.3.5 p.79): "$W\cong S_+\oplus S_-$", $Q=uQ_++ivQ_-$, $s=-u^2/v^2$.

**Signature (Intro p.4; §4.1 p.52; §10.1 p.68).**
- "we always consider complexified bundles of fields."
- "$\mathcal N=(1,0)$ super Yang–Mills theory on $M=\mathbb R^{10}$ with the Euclidean metric"
- "$\mathcal N=1$ super Yang–Mills theory on $M=\mathbb R^4$ with the Euclidean metric"
- Rem 4.1: "It is important here that we are working in the complexified setting. … in this real supersymmetry algebra there are no square-zero supercharges … which only exists in the split signature 5 + 5 … twists of maximal super Yang–Mills theory cannot be compatible with any unitary structure."

## SIGN CONVENTIONS

1. **Exponent weight.** None. No path integral and no $e^{\pm S/\hbar}$ appear.
2. **Grading.** Bigrading $\mathbb Z\times\mathbb Z/2$. Intro p.4: "$d$ has bidegree (1, 0) while $Q$ has bidegree (0, 1)". $Q_{BV}$ has degree +1.
3. **Ghost degree is −1 as a section (Rem 1.33).** $c\in\Omega^0(\mathfrak g[1])$. For $S$ to have degree 0 this forces $\deg A^*=+1$ and $\deg c^*=+2$, since $E=F\oplus F^![-1]$. This is opposite in sign to the physics ghost number of the field symbol $c$; the coordinate functions carry the opposite sign.
4. **Supersymmetry "ghost".** No field symbol is introduced. Supersymmetry enters only through cochains $\mathrm{Hom}(\wedge^n\mathfrak A,-)[-n]$, which are functions on $\mathfrak A[1]$ (Rem 1.28). The twist sets the argument to $Q$, and $S^{(k)}(Q,\dots,Q)$ has degree $-k$. The supergravity ghost $\eta$ (Intro p.6) "is bosonic: it lives in bidegree $(-1,1)$".
5. **Bracket.** $\{F,G\}=X_F(G)$, with $X_F=\omega\circ d_{\rm dR}F$. Stated symmetry: $\{F,G\}=(-1)^{|F||G|}\{G,F\}$. The bracket has degree $-n$, which is +1 for $n=-1$.
6. **Classical master equation.** $\{S,S\}=0$; $Q_{BV}I+\frac12\{I,I\}=0$; equivariant form $d_{\rm CE}S_{\mathfrak g}+\frac12\{S_{\mathfrak g},S_{\mathfrak g}\}=0$.
7. **BV action.** $S=\frac12\int\omega(e,Q_{BV}e)+I$. YM part $-\frac14(F_A,F_A)+\frac12(\lambda,\not d_A\lambda)$, with antifield couplings $+(\mathrm d_Ac,A^*)+([c,\lambda],\lambda^*)+\frac12([c,c],c^*)$.
8. **Clifford relation and $\Gamma$.** $v_1v_2+v_2v_1=+2(v_1,v_2)$; $\Gamma$ is defined by $(\rho(v)Q_1,Q_2)=(v,\Gamma(Q_1,Q_2))$; $\mathfrak{so}$ identification $\omega\mapsto(w\mapsto-2\iota_{(w,-)}\omega)$.
9. **Curvature.** $F_A=\mathrm dA+\frac12[A\wedge A]$.
10. **Source-internal inconsistencies, recorded and not resolved.**
    - Thm 2.10 (p.38) prints $\rho(\Gamma(Q_1,Q_2))Q_3+\rho(\Gamma(Q_2,Q_3))Q_1+\rho(\Gamma(Q_2,Q_3))Q_1=0$, with a repeated term. Assumption 3.1 has the cyclic form.
    - p.21 says $S$ has "cohomological degree $1-n$". That equals 2 for $n=-1$, while Def 1.16 and Def 1.27 require degree 0.
    - Def 2.7 uses the notation $(\mathcal E,S,\omega)$, while §2.5 uses $(E,S,\omega)$.

## GAPS (needed by the program, absent from ESW)

1. **No quantization.** Rem 1.7: "only relevant for quantization, which we do not consider in this paper". The source has no $\hbar$-weighted QME and no BV Laplacian $\Delta$ or $\Delta_P$. It also has no propagator $P(\epsilon,\infty)$ or heat kernel, no $W(P,I)$, no counterterm sign convention ($I-I^{CT}$ or $I+\dots$), no UV regulator, and no composite operators $O(\epsilon)$. Costello's framework [Cos11] is cited only in prose (pp.6, 21).
2. **No superconformal algebra.** Only the supertranslation algebra $\mathfrak A$ and $(\mathfrak{so}(V)\oplus\mathfrak g_R)\ltimes\mathfrak A$ appear. There are no $S$, $\bar S$, $K$ or dilatation generators, so no ghosts for them, and no $L_\infty$ action of $\mathfrak{psl}(4|4)$.
   - Supplementary, [SW19a] §3.1 p.13, eq. (7): $\mathfrak c(4,4)=\mathfrak{psu}(2,2|4)$. "In our considerations, we will always complexify, and thus deal with the complex Lie algebras $\mathfrak{sl}(4|\mathcal N)$ or $\mathfrak{psl}(4|4)$."
   - Generator labels in (8): $P_{\alpha\dot\alpha},Q^i_\alpha,\bar Q_{\dot\alpha i},M_{\alpha\beta},\Delta,R^i_j,\bar M_{\dot\alpha\dot\beta},S_i^\alpha,\bar S^{\dot\alpha i},K^{\alpha\dot\alpha}$.
   - [SW19a] gives no BV/$L_\infty$ action on untwisted N=4 fields.
3. **No untwisted 4d N=4 component action.** The source never writes the six scalars, the four Weyl fermions, or an explicit reduced action with antifields and $S^{(1)},S^{(2)}$. They are defined only implicitly, through Prop 2.20 applied to (31)–(34) with $\mathrm{Spin}(10)\supset\mathrm{Spin}(4)\times\mathrm{Spin}(6)$.
   - The reduction formulas for $\Sigma=S_+^{10}$ and $V^{10}$ are not written.
   - $S_\pm$ appear only as abstract 2-dimensional representations, with no index notation and no explicit $\Gamma$ in components.
4. **No explicit Euclidean real structure.** Spacetime is $\mathbb R^n$ with the Euclidean metric, but all fields are complexified. No reality conditions or conjugation on $\Sigma$ are given, there is no choice of $\int\mathrm{dvol}$ sign relative to a Euclidean weight, and no positivity statement.
5. **Normalizations not fixed.** The pairing on $\mathfrak g$ (no Killing-form normalization), the coupling constant (none appears), and the scale of $\Gamma$ ("typically unique up to a scale", p.36) are all left open.
6. **No scalar coupling.** 4d N=4 in ESW has no superpotential and no $\theta$-term (Intro p.4: "We do not consider superpotential, mass, or Fayet–Iliopoulos terms … we ignore all topological terms").
7. **Shift convention $V[n]$ not stated explicitly.** It is fixed only implicitly, by Rem 1.33 together with $c\in\Omega^0(\mathfrak g[1])$ sitting in degree −1.

Working files: the fetched PDF and extracted text are at /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/esw_pdf/esw.pdf and /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/esw_pdf/esw.txt. The [SW19a] PDF is at /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/sw_pdf/sw.pdf.