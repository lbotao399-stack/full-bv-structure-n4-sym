# Notation ledger: Costello and Gwilliam, *Factorization Algebras in Quantum Field Theory, Vol. 2*

## SOURCE

- Kevin Costello and Owen Gwilliam, *Factorization algebras in quantum field theory, Volume 2*. Authors' draft, title page "Volume 2 (27 October 2020)", 474 PDF pages. The published edition is Cambridge University Press 2021 (New Mathematical Monographs 41); this ledger does not use it.
- URL read: https://bpb-us-e2.wpmucdn.com/websites.umass.edu/dist/8/48667/files/2025/04/factorization2.pdf, linked from https://websites.umass.edu/ogwilliam.
- The old link https://people.math.umass.edu/~gwilliam/factorization2.pdf could not be reached.
- How it was read: WebFetch saved the PDF, the text was extracted with pdftotext, and every formula below was checked against rendered page images.
- Book page = PDF page − 10.
- Typography: the source writes the BV Laplacian as $\triangle$. In the transcriptions below, "(sic)" means the source's own typo.

## SYMBOLS

| symbol | meaning | location | verbatim defining formula |
|---|---|---|---|
| $\mathscr E$, $\langle-,-\rangle$, $Q$ | free BV fields, pairing, linear differential | Def 7.2.1.1 | "a graded antisymmetric map of vector bundles $\langle-,-\rangle_{loc}:E\otimes E\to \mathrm{Dens}(M)$ of cohomological degree $-1$"; "$Q:\mathscr E\to\mathscr E$ of cohomological degree 1 that is skew self adjoint" |
| $Q^{GF}$ | gauge-fixing operator | §7.2.1 | "square-zero differential operator of cohomological degree $-1$ … self adjoint … $D=[Q,Q^{GF}]$ is a generalized Laplacian" |
| $K_t$ | heat kernel | §7.2.3 | $\frac{d}{dt}K_t+(D\otimes 1)K_t=0$; $K_0$ is the kernel of the identity |
| $\Phi$ | parametrix (the "scale") | Def 7.2.4.2 | symmetric, degree 1, proper support, $([Q,Q^{GF}]\otimes 1)\Phi-K_{\mathrm{Id}}$ smooth |
| $\Phi<\Psi$ | partial order on parametrices | §7.2.4 | "if the support of $\Phi$ is contained in the support of $\Psi$" |
| $P(\Phi)$, $K_\Phi$ | propagator and smooth kernel | §7.2.5 | $P(\Phi)=\tfrac12 Q^{GF}\Phi$; $K_\Phi=K_{\mathrm{Id}}-QP(\Phi)$; $K_\Phi-K_\Psi=QP(\Psi)-QP(\Phi)$ |
| $P(0,L)$ | heat-kernel propagator | §7.2.5 | if $\Phi=\int_0^L K_t\,dt$ then $P(\Phi)=P(0,L)=\int_0^L(Q^{GF}\otimes1)K_t\,dt$ and $K_\Phi=K_L$ |
| $P(\epsilon,L)$, $P^L_\epsilon$ | cut-off propagator | §14.2.2 only (βγ example) | $P^L_\epsilon=\int_{t=\epsilon}^L\partial_z\frac{1}{4\pi t}e^{-\lvert z-z'\rvert^2/4t}$ |
| $\partial_\alpha$ | order-2 contraction operator | §7.2.7 | "the unique order 2 differential operator that is given by pairing with $\alpha$ on $\mathrm{Sym}^2\mathscr E^\vee$ and that is zero on $\mathrm{Sym}^{\le1}\mathscr E^\vee$" |
| $W(\alpha,-)$ | renormalization group (RG) flow | §7.2.7 | $F\mapsto\hbar\log\left(e^{\hbar\partial_\alpha}e^{F/\hbar}\right)$, on $\mathscr O^+(\mathscr E)[[\hbar]]$ (at least cubic mod $\hbar$) |
| $W_0$ | tree-level part | §7.2.7 | $W_0(P(\Phi),I)=W((P(\Phi),I)\bmod\hbar$ (sic) |
| $\triangle_\Phi$ | BV Laplacian at scale $\Phi$ | §7.2.8 | $\triangle_\Phi=\partial_{K_\Phi}$; $(\triangle_\Phi)^2=0$ |
| $\{-,-\}_\Phi$ | BV bracket at scale $\Phi$ | §7.2.8 | $\{I,J\}_\Phi=\triangle_\Phi(IJ)-(\triangle_\Phi I)J-(-1)^{\lvert I\rvert}I\triangle_\Phi J$ |
| $\mathscr O_{P,sm}$, $\mathscr O^+_{P,sm}$, $\mathscr O_{loc}$ | functional classes | §7.2.6, Def 7.1.0.2 | proper support; smooth first derivative; "+" means at least cubic mod $\hbar$; local means $\int_M(D_1\phi)\cdots(D_n\phi)\,d\mu$; all modulo constants |
| $I[\Phi]$ | scale-$\Phi$ effective interaction | Def 7.2.9.1 | see KEY STATEMENTS |
| $I_{i,k}[\Phi]$ | $k$-th Taylor component of the $\hbar^i$ coefficient | Def 7.2.9.1(iv) | — |
| $\mathrm{Obs}^q_\Phi(M)$, $\widehat Q_\Phi$ | BD algebra at scale $\Phi$ | Def 8.2.1.1, §8.3.1 | $(\mathscr O(\mathscr E)[[\hbar]],\,Q+\hbar\triangle_\Phi+\{I[\Phi],-\}_\Phi)$; $\widehat Q_\Phi=Q+\{I[\Phi],-\}_\Phi+\hbar\triangle_\Phi$ |
| $W^\Phi_\Psi$ | linearized RG flow on observables | Def 8.3.1.1 | $I[\Phi]+\delta W^\Phi_\Psi(f)=W(P(\Phi)-P(\Psi),I[\Psi]+\delta f)$ |
| $\delta$ | square-zero parameter | Def 8.3.1.1 | "square-zero parameter of cohomological degree $-i$" (for $f$ of degree $i$) |
| $O[\Phi]$, $\mathrm{Obs}^q(M)$, $\mathrm{Obs}^q(U)$ | quantum observables | Def 8.3.1.3, Def 8.4.0.1 | see KEY STATEMENTS |
| $\mathrm{Obs}^q(p)$ | quantum point observables ("local operators") | Def 10.1.0.1 | "the limit of $\mathrm{Obs}^q(V)$ as $V$ runs over the open sets containing $p$" |
| $\mathcal M$, $\mathcal M[1]$ | classical theory and its fields | Def 4.2.0.4, §4.4 | $\mathcal M$ is an elliptic $L_\infty$ algebra with invariant pairing of degree $-3$; "the graded vector bundle $L[1]$ on $M$ has a $-1$ symplectic form" |
| $C^*(\mathfrak g)$, $C^*_{red}(\mathfrak g)$, $d_{\mathfrak g}$ | Chevalley–Eilenberg cochains and differential | §11.2.2, §12.2.4 | "$C^*(\mathfrak g)=\mathscr O(B\mathfrak g)$"; "$\mathrm d_{\mathfrak g}$ is the Chevalley-Eilenberg differential"; $\mathscr O(\mathcal M[1])\otimes\mathscr O(\mathfrak g[1])$ (§11.2.2) |
| $\mathcal L$, $\mathrm d_{\mathcal L}$ | local $L_\infty$ algebra of symmetries | §12.2, §13.2.1 | "$\mathrm d_{\mathcal L}$ indicates the Chevalley differential on $\mathscr O(\mathcal L[1])$" |
| $\mathcal L\ltimes\mathcal M$ | semidirect product encoding an action | Def 12.2.2.1 | "$0\to\mathcal M\to\mathcal L\ltimes\mathcal M\to\mathcal L\to0$ consists of maps of $L_\infty$ algebras" |
| $\mathrm{Act}(\mathfrak g,\mathcal M)$ | dg Lie algebra controlling actions of a global $\mathfrak g$ | Def 12.2.4.1(iii) | $C^*_{red}(\mathfrak g)\otimes\mathscr O_{loc}(B\mathcal M)[-1]$ |
| $\mathrm{Act}(\mathcal L,\mathcal M)$, $\mathrm{InnerAct}(\mathcal L,\mathcal M)$ | same, for a local $\mathcal L$ | Def 12.2.2.4 | InnerAct is the kernel of $C^*_{red,loc}(\mathcal L\oplus\mathcal M)[-1]\to C^*_{red,loc}(\mathcal M)[-1]$; Act is InnerAct modulo the central subcomplex $C^*_{red,loc}(\mathcal L)$ |
| $\mathrm{Act}^q_{P,sm}$, $\mathrm{InnerAct}^q_{P,sm}$ | quantum versions | §13.2.2, §13.2.3 | $\mathscr O^+_{P,sm}(\mathcal L[1]\oplus\mathcal M[1])/\mathscr O_{P,sm}(\mathcal L[1])[[\hbar]]$; $\mathscr O^+_{P,sm}(\mathcal L[1]\oplus\mathcal M[1])[[\hbar]]$ |
| $I^{\mathcal L}[\Phi]$ | equivariant effective interaction | Def 13.2.2.1 | see KEY STATEMENTS |
| $S^{eq}$, $S^{inner}$, $S^{\mathfrak g}$, $S^{\mathcal L}$, $S^{tot}$ | equivariant action functionals | §11.2.2, §11.3.1, 12.2.4.1, 12.2.2.5 | $S^{tot}=S^{\mathcal L}+S$ |
| $\alpha$ | obstruction to an inner action (anomaly) | Prop 13.3.2.1 | "$\alpha\in\mathscr O_{loc}(\mathcal L[1])[[\hbar]]=C^*_{red,loc}(\mathcal L)$ of cohomological degree 1" |
| $\mathrm{Obs}^q_{\mathcal L}(U)$ | equivariant quantum observables | Prop 13.4.0.1 | $C^*(\mathcal L)$-module |
| $U^{BD}(\mathfrak g)$, $U^{BD}_\alpha$, $\mathbb U^{BD}_\alpha\mathcal L$ | enveloping BD algebras | Def 11.3.1.2, §13.1.2 | $U^{BD}(\mathfrak g)=C_*(\mathfrak g)[[\hbar]]=\mathrm{Sym}^*(\mathfrak g[-1])[[\hbar]]$, differential $d_{\mathfrak g}+\hbar d_{CE}$ |
| $I^{CT}[\epsilon]$, $\mathrm{Sing}_\epsilon$ | counterterm | §13.7.2 proof | $I^{CT}[\epsilon](\rho_1)=\mathrm{Sing}_\epsilon I^{(1)}_1[\epsilon,L](\rho_1)$ |
| $\mathsf c_{ij}$ | ghost operators of $so(n)$ | §12.3 Example | $\mathrm d\Psi_i=\mathsf c_{ij}\Psi_j$, $\mathrm d\mathsf c_{ij}=\mathsf c_{ij}\mathsf c_{jk}$ (sic); "in cohomological degree 1" |
| $T^{(n)}(\mathscr E,I)$ | simplicial set of theories mod $\hbar^{n+1}$ | §7.4, §7.5, §13.3.1 | not a scale |

## KEY STATEMENTS

**Def 7.2.9.1** (a quantum field theory consists of):
- (i) "$I[\Phi]\in\mathscr O^+_{P,sm}(\mathscr E_c)[[\hbar]]$" for every parametrix $\Phi$.
- (ii) RG flow: $I[\Phi]=W(P(\Phi)-P(\Psi),I[\Psi])$.
- (iii) The $\Phi$-QME: $(Q+\hbar\triangle_\Phi)e^{I[\Phi]/\hbar}=0$. "Equivalently, $QI[\Phi]+\hbar\triangle_\Phi I[\Phi]+\tfrac12\{I[\Phi],I[\Phi]\}_\Phi$." The source omits "$=0$" here.
- (iv) Locality: "$I_{i,k}(e_1,\dots,e_k)=0$ unless, for all $1\le r<s\le k$, $\mathrm{Supp}(e_r)\times\mathrm{Supp}(e_s)\subset\mathrm{Supp}(\Phi)^{3i+k}$".

**Def 7.2.9.2** (quantization of a classical theory): $I[\Phi]=I_0[\Phi]\bmod\hbar$, where $I_0[\Phi]=W_0(P(\Phi),I_0)$. Equivalently, "$\lim_{\Phi\to0}I[\Phi]-I_0\bmod\hbar=0$".

**§7.2.8** (RG flow preserves the QME): if $F$ satisfies the $\Phi$-QME, then $W(P(\Psi)-P(\Phi),F)$ satisfies the $\Psi$-QME, because $[Q,\partial_{P(\Phi)}-\partial_{P(\Psi)}]=\triangle_\Psi-\triangle_\Phi$.

**Thm 7.5.0.1** (obstruction theory): $T^{(n+1)}\to T^{(n)}\to DK(\mathscr O_{loc}(\mathscr E)[1],Q+\{I,-\})$ is a homotopy fibre sequence. When the obstruction class vanishes, the lifts are "a torsor for the simplicial Abelian group $DK(\mathscr O_{loc}(\mathscr E))[1]$".

**Def 8.3.1.3** (global quantum observables): "A global observable $O$ of cohomological degree $i$ is an assignment to every parametrix $\Phi$ of an element $O[\Phi]\in\mathrm{Obs}^q_\Phi(M)=\mathscr O(\mathscr E)[[\hbar]]$ of cohomological degree $i$ such that $W^\Phi_\Psi O[\Psi]=O[\Phi]$." The differential is
$$\widehat Q(O)[\Phi]=\widehat Q_\Phi(O[\Phi])=QO[\Phi]+\{I[\Phi],O[\Phi]\}_\Phi+\hbar\triangle_\Phi O[\Phi].$$
The text that follows says that $\{I[\Phi]+\delta O[\Phi]\}$ satisfies the axioms "needed to define a family of quantum field theories over the base ring $\mathbb C[\delta]/\delta^2$. The only axiom which is not satisfied is the locality axiom."

**Def 8.4.0.1** (support): "$O$ is supported on $U$ if, for each $(i,k)$, there exists a compact subset $K\subset U^k$ and a parametrix $\Phi$, such that $\mathrm{Supp}\,O_{i,k}[\Psi]\subset K$ for all parametrices $\Psi\le\Phi$."

**Lemma 8.4.0.2 proof** (the differential in components):
$$\widehat QO_{i,k}[\Phi]=QO_{i,k}[\Phi]+\sum_{a+b=i,\ r+s=k+2}\{I_{a,r}[\Phi],O_{b,s}[\Phi]\}_\Phi+\Delta_\Phi O_{i-1,k+2}[\Phi].$$

**Lemma 8.6.1.2 proof** (classical to quantum, mod $\hbar$): $O[\Phi]=\lim_{\Gamma\to0}W^\Phi_\Gamma(O)\bmod\hbar$. This is a sum over trees "with one special vertex … labelled by $O$, and every other vertex … labelled by the classical interaction $I_0$".

**Product of observables** (§8.5): $m(O,O')[\Phi]=\lim W^\Phi_\Psi(O[\Psi]*O'[\Psi])$. For point observables (10.3.2.1): $(O_1(0)\cdot O_2(x))[\Phi]=\lim_{\Psi\to0}W^\Phi_\Psi(I[\Psi],\tilde O_1(0)[\Psi]\tilde O_2(x)[\Psi])$.

**Prop 10.1.1.1 and Lemma 10.1.1.3** (local operators):
- Prop 10.1.1.1: "$\mathrm{Obs}^q(p)$ … is flat over $\mathbb C[[\hbar]]$, and it reduces modulo $\hbar$ to … $\mathrm{Obs}^{cl}(p)$."
- Lemma 10.1.1.3: "$\mathrm{Obs}^{cl}(0)$ are canonically identified with translation-invariant local functionals of both $\mathscr E$ and $\Omega^n(\mathbb R^n)$ that depend linearly on the background field."
- Proof: $\tilde S=S+\int\omega O(x)$, then $I^{b.g.}[\Psi]=I[\Psi]+I'[\Psi]$ and $O[\Psi]=\frac{\partial}{\partial\delta_0}I'[\Psi]$. Existence is cited from Costello (2011b), Thm 13.4.3 of Ch. 2 and Thm 3.1.2 of Ch. 3: "Here we do not ask that the quantum master equation hold."

**§11.3.1** (equivariant QME, finite-dimensional): "$\mathrm d_{\mathfrak g}S^{eq}+\hbar\triangle S^{eq}+\tfrac12\{S^{eq},S^{eq}\}=0$." In addition, $S^{eq}$ reduces to $S$ modulo the maximal ideal of $C^*(\mathfrak g)$. The inner version has the same equation for $S^{inner}\in(C^*(\mathfrak g)\otimes_{\mathbb C}\mathscr O(X)_q)/\mathbb C[[\hbar]]$.

**Def 11.3.1.1**: $\mathrm{Act}^q(\mathfrak g,\mathcal M)=(\mathscr O^q(\mathcal M[1])[-1]\otimes C^*_{red}(\mathfrak g))/C^*_{red}(\mathfrak g)[[\hbar]][-1]$, with differential $\mathrm d_{\mathfrak g}+\{S,-\}+\hbar\triangle$. The obstruction to an inner action lies in $H^1(C^*_{red}(\mathfrak g))[[\hbar]]$. The obstruction to lifting from mod $\hbar^n$ to mod $\hbar^{n+1}$ lies in "$H^2(\mathrm{Act}(\mathfrak g,\mathcal M))$".

**Def 12.2.4.1** (global $\mathfrak g$, classical): the source's (iii) reads "An degree one element $S^{\mathfrak g}$ in … $\mathrm{Act}(\mathfrak g,\mathcal M)=C^*_{red}(\mathfrak g)\otimes\mathscr O_{loc}(B\mathcal M)[-1]$ that satisfies $\mathrm d_{\mathfrak g}S^{\mathfrak g}+\mathrm d_{\mathcal M}S^{\mathfrak g}+\tfrac12\{S^{\mathfrak g},S^{\mathfrak g}\}=0$". This is read as "a family of classical field theories over the dg ring $C^*(\mathfrak g)$ that reduces to the original … modulo the maximal ideal $C^{>0}(\mathfrak g)$."

**Lemma 12.2.4.2** (global to local symmetry): actions of $\mathfrak g$ are canonically homotopy equivalent to actions of $\Omega^*_X\otimes\mathfrak g$. **Lemma 12.2.4.3**: if $\mathcal M$ is concentrated in degrees 1 and 2 (no gauge symmetry), the action factors through $\Omega^{\le1}_X\otimes\mathfrak g$.

**Lemma 12.2.2.5 and (12.2.2.2)** (local $\mathcal L$, classical): $(\mathrm d_{\mathcal L}+\mathrm d_{\mathcal M})S^{\mathcal L}+\tfrac12\{S^{\mathcal L},S^{\mathcal L}\}=0$, equivalently $\mathrm d_{\mathcal L}S^{tot}+\tfrac12\{S^{tot},S^{tot}\}=0$.
- Def 12.2.3.1 (inner action): $\mathrm d_{\mathcal L}S^{\mathcal L}+\tfrac12\{S^{\mathcal L},S^{\mathcal L}\}=0$ in InnerAct.
- Lemma 12.2.3.4: if no component of the structure maps sends $\mathcal L^{\otimes n}$ to $\mathcal M$, the action extends canonically to an inner action.

**Def 13.2.1.1** (equivariant QME, effective): "$\mathrm d_{\mathcal L}I+QI+\{I,I\}_\Phi+\hbar\triangle_\Phi I=0$." No ½ here.

**Def 13.2.2.1** (equivariant quantum field theory): $I^{\mathcal L}[\Phi]\in\mathrm{Act}^q_{P,sm}(\mathcal L,\mathcal M)$ satisfying four conditions.
- (i) $W(P(\Phi)-P(\Psi),I^{\mathcal L}[\Psi])=I^{\mathcal L}[\Phi]$.
- (ii) $\mathrm d_{\mathcal L}I^{\mathcal L}[\Phi]+QI^{\mathcal L}[\Phi]+\tfrac12\{I^{\mathcal L}[\Phi],I^{\mathcal L}[\Phi]\}_\Phi+\hbar\triangle_\Phi I^{\mathcal L}[\Phi]=0$.
- (iii) Locality, with the bounds of Def 7.2.9.1.
- (iv) Restricted to $\mathcal M[1]$, it gives $I[\Phi]$.
- Remark: this is "a family of theories over the dg base ring $C^*(\mathcal L(X))$".
- Def 13.2.3.1 (inner action): the same four conditions, in $\mathrm{InnerAct}^q_{P,sm}$.

**Prop 13.3.2.1** (anomaly): with $\tilde I^{\mathcal L}[\Phi]$ any lift,
$$\alpha=\mathrm d_{\mathcal L}\tilde I^{\mathcal L}[\Phi]+Q\tilde I^{\mathcal L}[\Phi]+\tfrac12\{\tilde I^{\mathcal L}[\Phi],\tilde I^{\mathcal L}[\Phi]\}_\Phi+\hbar\triangle_\Phi\tilde I^{\mathcal L}[\Phi].$$
$\alpha$ is independent of $\Phi$ and $\mathrm d_{\mathcal L}\alpha=0$. Lifting to an inner action "is clearly equivalent to making $\alpha$ exact in $C^*_{red,loc}(\mathcal L)[[\hbar]]$".

**§13.3.1** (equivariant obstruction theory): with the quantization fixed non-equivariantly, the obstruction-deformation complex is $\mathrm{Act}(\mathcal L,\mathcal M)$. For the full equivariant problem it is $C^*_{red,loc}(\mathcal L\ltimes\mathcal M\mid\mathcal L)$. For inner actions it is $C^*_{red,loc}(\mathcal L\ltimes\mathcal M)$ or $\mathrm{InnerAct}(\mathcal L,\mathcal M)$ respectively.

**Prop 13.4.0.1 proof** (equivariant observables):
- An element of degree $k$ is a family $O[\Phi]$ on $\mathcal L(X)[1]\oplus\mathcal M(X)[1]$ such that, for "$\epsilon$ … a square-zero parameter of cohomological degree $-k$", $W(P(\Phi)-P(\Psi),I^{\mathcal L}[\Psi]+\epsilon O[\Psi])=I^{\mathcal L}[\Phi]+\epsilon O[\Phi]$.
- Differential: $(\mathrm dO)[\Phi]=\mathrm d_{\mathcal L}O[\Phi]+QO[\Phi]+\{I^{\mathcal L}[\Phi],O[\Phi]\}_\Phi+\hbar\triangle_\Phi O[\Phi]$.
- Module structure: $(\alpha\cdot O)[\Phi]=\alpha\cdot(O[\Phi])$, with $\mathrm d(\alpha\cdot O)=(\mathrm d\alpha)\cdot O+\alpha\cdot\mathrm dO$.

**Thm 13.1.2.2 and Thm 13.5.1.1** (quantum Noether): there is a map $\mathbb U^{BD}_\alpha(\mathcal L_c)\to\mathrm{Obs}^q$, which refines to $\mathrm{Rees}\,\mathbb U^{\hbar}_\alpha\mathcal L\to\mathrm{Obs}^q$ over $\mathbb C[[\hbar]]$.

**§13.7.2 proof** (the only explicit ε-counterterm in the book):
- $I^{(1)}_1[\epsilon,L](\rho_1)=\int_{t=\epsilon}^L\int_{x\in X}\mathrm{STr}_E\,\rho_1(x)\,(\not\partial_xK_t(x,y))_{y=x}\,\mathrm dt$.
- $I^{CT}[\epsilon](\rho_1)=\mathrm{Sing}_\epsilon I^{(1)}_1[\epsilon,L](\rho_1)$.
- $I^{(1)}_1[L](\rho_1)=\lim_{\epsilon\to0}\left(I^{(1)}_1[\epsilon,L](\rho_1)-I^{CT}[\epsilon](\rho_1)\right)$.
- QME at this order: $\mathrm d_{dR}I^{(1)}_1[L]+\triangle_LI^{(1)}_0=0$, with $\mathrm d_{dR}I^{(1)}[\epsilon,L]=\triangle_\epsilon I^{(1)}_0-\triangle_LI^{(1)}_0$.
- Obstruction: $\alpha^{(1)}(\rho_0)=\lim_{\epsilon\to0}\left(\triangle_\epsilon I^{(1)}_0(\rho_0)-I^{CT}[\epsilon](\mathrm d\rho_0)\right)$.
- The counterterm term "must vanish because it is a purely singular function of $\epsilon$ (… a linear combination of $\log\epsilon$ and $\epsilon^{-k}$ for $k>0$)."
- Source typos on the same page: "$I^{(1)}[L]=I^{(1)}_0[L]+\hbar I^{(1)}_0[L]$" (sic, the second term should be $I^{(1)}_1$) and "$\triangle_L=K_L$" (sic).

## SIGN CONVENTIONS

1. **Exponent weight.**
   - The formal definitions use $e^{+I[\Phi]/\hbar}$, both in the QME and in $W=\hbar\log(e^{\hbar\partial_\alpha}e^{F/\hbar})$.
   - The motivating Remark in §8.3.1 uses $\langle O\rangle=\frac1{Z_S}\int O\,e^{-S/\hbar}$ and $e^{(-S+\hbar\delta O)/\hbar}$. There $O$ comes with weight $\hbar^0$, while in Def 8.3.1.1 it is $\delta O$ inside $e^{(I+\delta O)/\hbar}$.
   - The source does not reconcile these.
   - In finite dimensions (§6.2) the volume form is $e^S dx$ and $\mathrm{div}_S=\triangle+\{S,-\}$.
2. **QME.** $QI+\hbar\triangle_\Phi I+\tfrac12\{I,I\}_\Phi=0$: plus sign on $\hbar\triangle$, factor ½ on the bracket. There is no separate free action; $Q$ carries the quadratic part.
3. **Equivariant QME: the ½ is inconsistent.** It appears in §11.3.1, 13.2.2.1(ii) and 13.3.2.1, but is missing in 13.2.1.1.
4. **Counterterm sign.** The finite part is $\lim_{\epsilon\to0}(I[\epsilon,L]-I^{CT}[\epsilon])$, with $I^{CT}=\mathrm{Sing}_\epsilon$. The counterterm is pure $\log\epsilon$ and $\epsilon^{-k}$, with no finite part. This is shown only for a one-loop interaction term, not for observables.
5. **Kernel identities.**
   - $K_\Phi=K_{\mathrm{Id}}-QP(\Phi)$ and $Q(P(\Phi)-P(\Psi))=K_\Psi-K_\Phi$ (§13.2.1(ii)).
   - Substituting the §7.2.5 heat-kernel parametrices gives $Q(P(0,L)-P(0,\epsilon))=K_\epsilon-K_L$.
   - Conflict: §14.2.2 states $\bar\partial P(\epsilon,L)=K_L-K_\epsilon$, the opposite sign. It is not resolved in the source.
6. **Bracket vs Laplacian: the sign differs between chapters.**
   - §7.2.8: $\{I,J\}_\Phi=\triangle_\Phi(IJ)-(\triangle_\Phi I)J-(-1)^{\lvert I\rvert}I\triangle_\Phi J$.
   - §6.3: $\mathrm{div}(ab)=(\mathrm{div}\,a)b+(-1)^{\lvert a\rvert}a(\mathrm{div}\,b)+(-1)^{\lvert a\rvert}\{a,b\}$.
7. **Degrees.**
   - Field pairing: degree $-1$. $\mathcal M$ pairing: degree $-3$. Fields are $\mathcal M[1]$, and $C^*(\mathcal L)\cong\mathscr O(\mathcal L[1])$ as graded algebras.
   - $\Phi$ has degree 1; $P(\Phi)$ has degree 0 (§7.2.7 takes $\alpha$ of degree 0).
   - $\delta$ (or $\epsilon$) has degree $-i$ for $O$ of degree $i$.
   - In the filtered BD algebra, $\hbar$ has weight $-1$ (§6.4).
   - Ghosts $\mathsf c_{ij}$ have degree 1.
   - "$\mathbb Z/2\times\mathbb Z$ graded … Both gradings contribute to signs" (§12.2.2 Example).
   - $S^{\mathfrak g}$ is "degree one" in the shifted Lie algebra (12.2.4.1); $S^{\mathcal L}$ is "of cohomological degree 0" as a functional (12.2.3.1).
8. **Two inconsistencies in the enveloping BD algebra.**
   - Def 11.3.1.2 uses $\mathrm{Sym}(\mathfrak g[1])$, while §13.1.2 uses $\mathrm{Sym}^*(\mathfrak g[-1])$.
   - Def 11.3.1.2 writes the differential as $d(xy)=d_{\mathfrak g}(x)y\pm x\,d_{\mathfrak g}(y)+\hbar[x,y]$, while §13.1.2 writes it as $d_{\mathfrak g}+\hbar d_{CE}$.

## GAPS

1. **No $P(\epsilon,\infty)$, no $T=\infty$, no IR-free formulation.**
   - The scale is a parametrix $\Phi$ with proper support, and all limits are $\Phi\to0$.
   - The only heat-kernel propagators are $P(0,L)$ (§7.2.5) and $P^L_\epsilon$ (§14.2.2).
   - "On a non-compact manifold $M$, there is more than one heat kernel" (§7.2.3).
   - The source has no statement allowing the upper limit to go to infinity.
2. **No explicit $O(\epsilon)=O+O^{CT}(\epsilon)$ for composite operators.**
   - Quantum point observables exist abstractly: background-field coupling $\int\omega O$ (Lemma 10.1.1.3), with existence cited from Costello (2011b), without the QME.
   - Explicit tree-level lifting is given only mod $\hbar$ (Lemma 8.6.1.2).
   - The only explicit ε-counterterm is §13.7.2, for an interaction term linear in a background field.
   - No two-loop formulas.
3. **No supersymmetry, no superconformal algebra, no Weyl spinors, no Euclidean spinor conventions.**
   - Only holomorphically twisted N=1, N=2 and Kapustin–Witten-twisted N=4, at the classical cotangent level (§4.6.4–4.6.6).
   - No explicit treatment of a Lie superalgebra $\mathfrak g$ with odd generators.
4. **Ghosts of a global symmetry.**
   - The source works only through $C^*(\mathfrak g)=\mathscr O(B\mathfrak g)$, $d_{\mathfrak g}$ and $C^*_{red}$.
   - It names no coordinates on $\mathfrak g[1]$, except $\mathsf c_{ij}$ in the $so(n)$ chiral fermion example.
   - There is no explicit formula for $d_{\mathfrak g}$ on generators beyond $\mathrm d\mathsf c_{ij}=\mathsf c_{ij}\mathsf c_{jk}$ (sic, as printed).
5. **Quantum global-to-local symmetry is announced but not found.** §13.1.3 says the quantum analogue of Lemma 12.2.4.2 will be proved "below". A text search found no separate statement.
6. **Notation collisions the program must avoid.**
   - The source's $\epsilon$ means two things: a square-zero parameter (§13.4) and the UV cut-off (§13.7).
   - $T^{(n)}$ is a simplicial set of theories, not an IR scale.
7. **No dimensional regularization.** The text search found no use of dimensional regularization. The regulator is the ε cut-off with $\mathrm{Sing}_\epsilon$ subtraction.

Extracted text, used only as a reading aid: /private/tmp/claude-501/-Users-libotao-Desktop/c75e7fb6-0585-46ee-8195-29049c94d407/scratchpad/cg2/cg2.txt