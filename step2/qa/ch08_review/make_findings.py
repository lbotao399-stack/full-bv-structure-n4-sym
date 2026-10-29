#!/usr/bin/env python3
"""Writes findings.json for the referee report on manuscript/ch08_definitions.tex."""
import json, os

HERE = os.path.dirname(os.path.abspath(__file__))

F = []


def add(id_, location, severity, claim, evidence, replacement):
    F.append(dict(id=id_, location=location, severity=severity, claim=claim, evidence=evidence, replacement=replacement))


# ------------------------------------------------------------------ major
add('M1', 'prop:c8-Tindep (i), line 108; proof (i), line 114', 'major',
    r"The singular part of every graph with edges $P(\varepsilon,T)$ is independent of $T\in(0,\infty]$.",
    r"True for one-loop graphs (all singular coefficients come from the region where the Schwinger times of the single cycle are small). "
    r"False for unsubtracted graphs with a divergent subgraph at higher order: the region where the subgraph times are small and the outer times are large "
    r"gives (singular coefficient of the subgraph) x (function of T); only the subtracted sums are T-independent. The chapter claims all orders "
    r"(Remark 'The definitions at higher order'). The design's version (odd_components_design lem:o-Tflow (iv)) is one-loop only. "
    r"At T=infinity the graphs with S^(1)(K_v) floating on a cycle diverge at large distance for every eps (Prop c8-Tsafe), so their singular part is undefined there.",
    r"\item (i) At order $\hbar$ the singular part of every graph with edges $P(\varepsilon,T)$ is independent of $T\in(0,\infty)$, and of $T=\infty$ in the components whose graphs converge at large distance. At every order in $\hbar$ the counterterms $I^{CT}_{\mathfrak{psl}(4|4)}[\varepsilon]$ constructed with $P(\varepsilon,T)$ are the same for all $T\in(0,\infty)$ [C07, §1.5]."
    r"  || Proof (i): At order $\hbar$, in the Schwinger representation of Lemma \ref{lem:c4-local} the coefficients of $\varepsilon^{-k}$ and $\log\varepsilon$ come from the region in which the Schwinger times of the cycle sum to at most a fixed $\tau_0<T$; there the integration domain does not involve $T$, and the complementary region gives a function with a limit as $\varepsilon\to0$. At higher order, if $\lim_{\varepsilon\to0}W(P(\varepsilon,T),J)$ exists, so does $\lim_{\varepsilon\to0}W(P(\varepsilon,T'),J)=W\big(P(T,T'),\lim_{\varepsilon\to0}W(P(\varepsilon,T),J)\big)$, the edges $P(T,T')$ being smooth, and uniqueness of the counterterms in the scheme [C07, Thm B] identifies the counterterms at $T$ and at $T'$.")

add('M2', 'thm:c8-invariants (iii), line 268; proof (iii), line 273', 'major',
    r"The spectrum of the operator induced by $\Lambda^{(1)}(\Delta)$ on classes of fixed classical weight is invariant under (a), (b) and (c).",
    r"(c) of prop:c8-laws is an arbitrary s-closed G. Under it Lambda^(0) changes by -a hbar ad_{G^(0)} and Lambda^(1)(Delta) by -a hbar ad_{G^(1)(Delta)}; a closed equivariant deformation other than the coupling and G_B can change the spectrum (it changes the theory). "
    r"The proof treats only conjugations, the coupling direction (Lambda^(1) unchanged) and G_B (B Delta=0, G_B^(0)=0). Whether H^0 is spanned by S-S^(2), G_B and Im s is not established. "
    r"This repeats the overclaim flagged in design_audit_corrections.json (design_audits[0] item 2, correction 3: restrict to (a),(b), coupling, G_B, or add hypothesis E3b/H4).",
    r"\item (iii) The spectrum of the operator induced by $\Lambda^{(1)}(\Delta)$ on the cohomology classes of a fixed classical Weyl weight (the one-loop anomalous dimensions) is invariant under (a), under (b), under the coupling direction and under the outer direction $G_B$; under (c) with an arbitrary $s$-closed $G$ it is invariant when every $s$-closed local cochain of degree $0$ and Weyl weight $4$ lies in $\mathbb C(S-S^{(2)})+\mathbb C\,G_B+\operatorname{Im}s$."
    r"  || Proof (iii): Conjugation preserves the spectrum; the coupling direction rescales $\Lambda^{(0)}$ and leaves $\Lambda^{(1)}$ unchanged; for $G_B$, $G_B^{(0)}=0$ and $B\Delta=0$; under the stated hypothesis a closed $G$ is a combination of these and of $s\phi$, which is (a).")

add('M3', 'prop:c8-2pt (i), line 463; proof (i), line 474', 'major',
    r"The two-point function of a $Q+\{I,-\}$-exact operator with a closed one vanishes.",
    r"At one loop the BRST Ward identity gives 0 = <L(A)^(0), O_2> +- <A, L(O_2)^(0)> with L(A)^(0) = (Q+{I,-})A + hbar L_1(A)^(0) (two_point_design cor:t-exact). "
    r"Hence <(Q+{I,-})A, O_2> = -hbar <L_1(A)^(0), O_2>_tree, of relative order one unit and nonzero in general; it also needs L(O_2)^(0)=0 through order hbar (quantum closedness), "
    r"whereas Def c8-2pt only assumes (Q+{I,-})O_a=0, plus the QME and disjoint supports. The proof given covers only the first sentence (constants, Q-exact terms).",
    r"\item (i) Constant terms of $O_a(\varepsilon)$ and $Q$-exact terms of the flow do not contribute. Assume \eqref{eq:c8-qme} modulo $\hbar^2$, and local representatives with $L(O_2)^{(0)}=0$ modulo $\hbar^2$. Then $\langle L(A)^{(0)},O_2\rangle=0$ through relative order one unit for every local $A$ of degree $-1$ with a local representative, where $L(A)^{(0)}=(Q+\{I,-\})A+\hbar L_1(A)^{(0)}$; the two-point function descends to the cohomology of $\Lambda^{(0)}$."
    r"  || add to the proof: The last statement is the arity-zero, $t_1t_2$, $e=0$ part of the flow of the defect vertex of $J+t_1A(\varepsilon)+t_2O_2(\varepsilon)$: the left side vanishes by the first statement, and the right side is $\langle L(A)^{(0)},O_2\rangle\pm\langle A,L(O_2)^{(0)}\rangle$.")

add('M4', 'cor:c8-gram, line 506 (covariance under finite parts); proof line 510', 'major',
    r"Under $a\mapsto(1+\hbar\psi')a$: $X_1\mapsto X_1+[\psi',X_0]$ and $G_{\rm tree}^{-1}G_1\mapsto G_{\rm tree}^{-1}G_1+\psi'+\psi'^{\rm adj}$, and both sides change covariantly.",
    r"With the conventions of the proof ($\langle a,b\rangle=\bar a^TGb$, $\bar X^TG=GY$, $Y^{\rm adj}=G^{-1}\bar Y^TG$, columns), the change is $G\mapsto\bar g^TGg$, $X\mapsto g^{-1}Xg$, $g=1+\hbar\psi'$ "
    r"(consistent with prop:c8-laws (b): L conjugated as g^{-1}Lg). First order: $X_1\mapsto X_1+[X_0,\psi']$, $M\mapsto M+\psi'+\psi'^{\rm adj}$. "
    r"Numerically (checks_algebra.py C2-C7): exact expansion gives X_1+[X_0,psi]; the printed pair (X_1+[psi,X_0], M+psi+psi^adj) violates the relation (defect 1.8e2), the corrected pair preserves it (defect 3e-14). "
    r"IP-P1 of ip_targets.json states the other consistent pair (X_1+[Z,X_0], M-Z-Z^adj, i.e. Z=-psi); ch08 combined the X-rule of one convention with the M-rule of the other. "
    r"Related: prop:c8-2pt (iii) uses the row convention D_1(O_a)=sum_b (D_1)_{ab}O_b, whereas def:c8-inner and this proof use columns.",
    r"Under a finite part $a\mapsto(1+\hbar\psi')a$ of the operators, with $\psi'$ reflection compatible, $G\mapsto\overline{(1+\hbar\psi')}^{\,T}G\,(1+\hbar\psi')$ and $X\mapsto(1+\hbar\psi')^{-1}X(1+\hbar\psi')$ for the matrix of every generator, so $X_1\mapsto X_1+[X_0,\psi']$ and $G_{\rm tree}^{-1}G_1\mapsto G_{\rm tree}^{-1}G_1+\psi'+\psi'^{\rm adj}$, and the relation is preserved."
    r"  || in def:c8-inner add: Matrices act on coefficient columns, $YO_b=\sum_aY_{ab}O_a$, and $\langle a,b\rangle_{\rm rad}=\overline a^{\,T}(G_{\rm tree}+\hbar G_1)b$."
    r"  || in prop:c8-2pt (iii) replace '$D_1(O_a)=\sum_b(D_1)_{ab}O_b$ ... $-\tfrac14(D_1C+CD_1^{\,T})$' by '$D_1(O_b)=\sum_a(D_1)_{ab}O_a$ ... $-\tfrac14(D_1^{\,T}C+CD_1)$'.")

add('M5', 'def:c8-twocomp (lines 14-35) and prop:c8-twocomp (line 39)', 'major',
    r"The formulas of Def c8-twocomp realise $\rho$, $(-,-)$, $\Gamma$ of Def c1-spaces exactly: Clifford relation, symmetry, [ESW, Assumption 3.1], chirality.",
    r"Def c8-twocomp gives rho(v) only on Sigma^* for v in V, Gamma only along V, and the 't Hooft blocks only through sum_I Sig^I_ij Sig^I_kl = 2 eps_ijkl (which fixes them only up to SO(6)). "
    r"rho on Sigma follows from rho(e_m)^2=1 (checks_twocomp T4), but rho(e_I), Gamma|_W (G2),(G4) and the explicit blocks are missing; the V-data are compatible with the sign-flipped W_perp action (T5). "
    r"These data enter the ten-dimensional Clifford relation, ESW Assumption 3.1 (T1, T9 pass only with (R3),(R4) of the dictionary), the term (1/2)eps_jikl w^j((u lambda^i),A^{*kl}) and the commutator term of S^(1)(xi) in Remark c8-twist, S^(2), and A^{ij}, w^{ij}. "
    r"With the dictionary data every printed formula checks (T1-T3, T6-T9).",
    r"Insert after \eqref{eq:c8-pairing}: With $\Sigma^I_{ij}=(\eta^1,\eta^2,\eta^3,i\bar\eta^1,i\bar\eta^2,i\bar\eta^3)_{ij}$, $\eta^a_{bc}=\bar\eta^a_{bc}=\epsilon_{abc}$, $\eta^a_{b4}=-\eta^a_{4b}=\delta_{ab}$, $\bar\eta^a_{b4}=-\bar\eta^a_{4b}=-\delta_{ab}$, $w_{ij}=w_I\Sigma^I_{ij}$, $w^{ij}=\tfrac12\varepsilon^{ijkl}w_{kl}$ for $w\in W_{\mathbb R}^\perp$:"
    r" \[ (\rho(v)\xi)^\alpha{}_i=i\,v_m\bar\xi_{\dot\alpha i}\bar\sigma^{m\dot\alpha\alpha},\quad (\rho(v)\xi)^{\dot\alpha i}=-i\,v_m\bar\sigma^{m\dot\alpha\alpha}\xi_\alpha{}^i,\quad (\rho(w)\xi)^\alpha{}_j=\varepsilon^{\alpha\beta}w_{ji}\xi_\beta{}^i,\quad (\rho(w)\xi)^{\dot\alpha j}=\varepsilon^{\dot\alpha\dot\beta}w^{ji}\bar\xi_{\dot\beta i},\quad (\rho(w)\eta)_\alpha{}^k=-w^{kj}\varepsilon_{\alpha\beta}\eta^\beta{}_j,\quad (\rho(w)\eta)_{\dot\beta k}=-w_{kj}\varepsilon_{\dot\beta\dot\alpha}\bar\eta^{\dot\alpha j} \]"
    r" \[ \Gamma(\xi_1,\xi_2)_{kl}=2\varepsilon_{jikl}\xi_{1\alpha}{}^j\varepsilon^{\alpha\beta}\xi_{2\beta}{}^i-2\big(\bar\xi_{1\dot\alpha k}\varepsilon^{\dot\alpha\dot\beta}\bar\xi_{2\dot\beta l}-\bar\xi_{1\dot\alpha l}\varepsilon^{\dot\alpha\dot\beta}\bar\xi_{2\dot\beta k}\big),\qquad \Gamma(\eta_1,\eta_2)^{kl}=-2\varepsilon^{klpj}\varepsilon_{\alpha\beta}\eta_1{}^\beta{}_j\eta_2{}^\alpha{}_p+2\big(\varepsilon_{\dot\beta\dot\alpha}\bar\eta_1^{\dot\alpha l}\bar\eta_2^{\dot\beta k}-\varepsilon_{\dot\beta\dot\alpha}\bar\eta_1^{\dot\alpha k}\bar\eta_2^{\dot\beta l}\big), \]"
    r" with $\Gamma_{kl}=\Gamma^I\Sigma^I_{kl}$, $\Gamma^{kl}=\tfrac12\varepsilon^{klij}\Gamma_{ij}$ and $(\Gamma,u)=\Gamma^mu_m+\tfrac14\Gamma_{ij}u^{ij}$.")

add('M6', 'thm:c8-L (iv) line 204; cor:c8-homotopy; def:c8-twisted; thm:c8-twist hypotheses line 307 and proof line 311', 'major',
    r"$L\circ L=0$ is the Maurer--Cartan equation $d_{CE}\Lambda+\Lambda d_{CE}+\Lambda^2=0$; by Lemma c8-ev, $r_\xi$ maps it to $r_\xi(\Lambda)^2=0$.",
    r"The product of operator-valued cochains (Lambda^2) and the signs of the C^bullet-linear extension (L(cM)=(d_CE c)M +- cL(M)) are never defined; Lemma c8-ev proves compatibility of r_xi with d_CE and with the cochain bracket only, not with the product, "
    r"so r_xi(Lambda^2)=r_xi(Lambda)^2 and the coefficients -1/2, -1/6 in Def c8-twisted rest on an unstated convention. This is hypothesis (H1) of the design (thm:o-twist) and design_audit_corrections.json design_audits[0] item 1 / correction 2; ch08 drops it. "
    r"With the shuffle sign of def:c2-dce applied to composition, all-equal odd arguments give (-1)^{kl}, and a_k a_l = a_{k+l} binom(k+l,k) (-1)^{kl} (checks_algebra I1) gives multiplicativity, so the gap closes by a definition. "
    r"Thm c8-twist also uses thm:c8-L (iv), which assumes eq:c8-qme, but does not assume it.",
    r"In thm:c8-L (iv), after 'with values in the operators on $\mathscr O_l(\mathscr E)[[\hbar]]$' add: where $(\Lambda\Lambda')(x_1,\dots,x_n)$ is the sum over $(k,l)$-shuffles $\sigma$ of $\operatorname{sgn}(\sigma)\cdot(\text{Koszul sign of }\sigma)\cdot(-1)^{|\Lambda'|k+\sum_{i=1}^k|x_{\sigma(i)}|(l+|\Lambda'|)}\,\Lambda^{(k)}(x_{\sigma(1)},\dots,x_{\sigma(k)})\circ\Lambda'^{(l)}(x_{\sigma(k+1)},\dots,x_{\sigma(n)})$, the sign rule of Definition \ref{def:c2-dce}, and the same rule defines the products of field-independent cochains with local cochains."
    r"  || In lem:c8-ev add: 'and the product of operator-valued cochains with the composition of their values' after 'the bracket of values'; proof: 'the same shuffle count applies to the product'."
    r"  || thm:c8-twist: 'Assume \eqref{eq:c8-qme}, that the local representatives exist on a space of local functionals stable under $\Lambda$, and that the sums defining $r_\xi(\Lambda)$ are finite there.'")

# ------------------------------------------------------------------ minor
add('m1', 'def:c8-B, line 143', 'minor',
    r"$B$ is induced by the bracket with $\operatorname{diag}(1_4,0_4)$ in $\mathfrak{gl}(4|4)$.",
    r"With the block table of ch02 (S_+(x)W = block (3,1), S_-(x)W = (3,2), S_-(x)W^* = (2,3), S_+(x)W^* = (1,3); the same convention reproduces ad_Delta = -2,+2,-1,-1,+1,+1), "
    r"[diag(1_4,0_4),X] = -X on (3,1),(3,2) and +X on (1,3),(2,3), the opposite of B (checks_algebra A1-A3); B = ad diag(0_4,1_4). The same slip is in the design and in design_audit_corrections.json. The eigenvalues and all uses of B are correct.",
    r"it is induced by the bracket with $\operatorname{diag}(0_4,1_4)$ in $\mathfrak{gl}(4|4)$, equivalently minus the bracket with $\operatorname{diag}(1_4,0_4)$.")

add('m2', 'def:c8-ob line 123; thm:c8-obstruction (ii) line 131, proof line 138', 'minor',
    r"Ob is the local cochain whose tree-level flow equals the order-$\hbar$ part of the limit; $\hbar F_1$ is a quantisation exactly when $sF_1=\mathrm{Ob}$.",
    r"The order-hbar part carries the factor hbar (in ch04 eq:c4-obstruction the limit is -27/8 hbar tr/(16 pi^2 (X,Y)) int(...)), while replacing F=0 by hbar F_1 changes the defect by -hbar sF_1. "
    r"As written the condition should read hbar sF_1 = Ob. Signs agree with ch04 (F enters as I - F, QF = Ob with F = 27/16 unit int(d^*A|_V,d^*A|_V)).",
    r"def:c8-ob: 'whose tree-level flow with edges $P(0,T)$ equals the coefficient of $\hbar$ in this limit'; proof of (ii): 'the coefficient of $\hbar$ of the obstruction becomes $\mathrm{Ob}-sF_1$.'")

add('m3', 'thm:c8-obstruction (iv) line 133; proof line 138', 'minor',
    r"Two classes of $H^0$ are always present: the coupling direction $S-S^{(2)}$ and $G_B$.",
    r"The proof shows S - S^(2) is s-closed, but the claim that its class is nonzero is never proved (prop:c8-B shows only that G_B is not in C(S-S^(2)) + Im s). "
    r"The gap closes with Step-1 results: (s phi)^(0) = {S, phi^(0)}; phi^(0) of degree -1 and weight 4 is a combination of field-counting functionals (lem:c4-enumeration); {S, field count} = S would need an additive weight with 2a=3a=4a=1 on the (A|_V)^2,3,4 terms (cor:c4-coupling proof). "
    r"This needs [g,g] != 0.",
    r"add to the proof of (iv): The class of $S-S^{(2)}$ is nonzero: in arity zero $(s\phi)^{(0)}=\{S,\phi^{(0)}\}$, with $\phi^{(0)}$ of degree $-1$ and Weyl weight $4$ a combination of the field-counting functionals (Lemma \ref{lem:c4-enumeration}), and $\{S,\phi^{(0)}\}=S$ is excluded by the computation of Corollary \ref{cor:c4-coupling}.")

add('m4', 'def:c8-2pt, line 455', 'minor',
    r"$W(P(\varepsilon,\infty),I-F_0-I^{CT}[\varepsilon]+t_1O_1(\varepsilon)+t_2O_2(\varepsilon))$",
    r"F_0 is defined nowhere in ch08; with the chapter's own indexing F_psl = hbar F_1 + hbar^2 F_2 + ... (def:c8-quant), F_0 reads as the order-hbar^0 term, which is zero; that reading drops the forced 27/16 counterterm. "
    r"In two_point_design def:t-2pt, F_0 is the arity-zero finite counterterm. The square-zero parameters t_a should also be specified as of degree 0 and of the parity of O_a, as in the design.",
    r"W\big(P(\varepsilon,\infty),\,I-F^{(0)}_{\mathfrak{psl}(4|4)}-I^{CT}[\varepsilon]+t_1O_1(\varepsilon)+t_2O_2(\varepsilon)\big) ... with $F^{(0)}_{\mathfrak{psl}(4|4)}$ the arity-zero part of the finite counterterm (at one loop $\tfrac{27}{16}\frac{\hbar\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{16\pi^2(X,Y)}\int\mathrm{dvol}\,(d^*A|_V,d^*A|_V)$ plus $s$-closed choices) and square-zero parameters $t_1,t_2$ of degree $0$ and of the parities of $O_1,O_2$.")

add('m5', 'cor:c8-gram lines 504-506', 'minor',
    r"$S_1=Q_1^{\rm adj}+[S_0,G_{\rm tree}^{-1}G_1]$, $K_1=P_1^{\rm adj}+[K_0,\dots]$; under $aG_B$ 'the first relation holds exactly when $a$ is imaginary'.",
    r"The relation (X^dag)_1 = X_1^adj + [(X^dag)_0, M] is correct (checks_algebra C1). In the 'in particular' form, S, K must be the conjugate generators X^dag, which carry signs: (Q^i_alpha)^dag = -S^alpha_i, (Qbar)^dag = +Sbar, P_b^dag = -K_{bbar} (T16). "
    r"Read as the components of def:c8-components along the unit generators, the relation has the wrong sign for the Q^i_alpha/S^alpha_i pair. "
    r"Q_1 here is the one-loop part of the matrix of -Lambda^(1)(Q), i.e. minus the def:c8-components functional; the relation is homogeneous, so this is consistent once stated. "
    r"Outer direction: C9 confirms the shift preserves the relation iff conj(a) = -a; and G_B is anti-invariant under the conjugation (T22: B anticommutes with dagger and dagger_n), so aG_B is reflection compatible iff a is imaginary. "
    r"Under the corollary's hypotheses (reflection-compatible scheme) a is already imaginary, so 'holds exactly when' should read 'is preserved exactly when'. Separately, K_0 also names the kernel of the identity (line 348).",
    r"in particular, with $X_0+\hbar X_1$ the matrix of $-\Lambda^{(1)}(X)$ and $Y_0+\hbar Y_1$ that of $-\Lambda^{(1)}(X^\dagger)$ for $X=Q^i_\alpha$ ($X^\dagger=-S^\alpha{}_i$), $X=\bar Q_{\dot\alpha i}$ ($X^\dagger=\bar S^{\dot\alpha i}$), $X=P_b$ ($X^\dagger=-K_{\bar b}$) and $X=\Delta$: $Y_1=X_1^{\rm adj}+[Y_0,G_{\rm tree}^{-1}G_1]$. ... The shift $F_{\mathfrak{psl}(4|4)}\mapsto F_{\mathfrak{psl}(4|4)}+\hbar aG_B$ changes $X_1$ by $-aX_0$ and $Y_1$ by $+aY_0$ for $X=Q^i_\alpha$ (Proposition \ref{prop:c8-laws}); it preserves the first relation exactly when $\bar a=-a$, which is also the condition for $aG_B$ to be reflection compatible, $B$ anticommuting with $\dagger_n$: the conjugation fixes the real part of the outer direction.")

add('m6', 'thm:c8-dag (iii) line 383; proof line 388', 'minor',
    r"Positivity, in a representation in which $-\tfrac12\Delta$ is self-adjoint with a positive eigenvalue on a vector annihilated by $\Sigma^*$, selects the sign $+$.",
    r"The computation ||Qu||^2 = <u,[Q^natural,Q]u> uses (1) a representation pi with pi(X^natural) = pi(X)^adj for the conjugation being tested, and (2) the rotation parts of [Q^natural,Q] (so(V), so(W_perp)) acting trivially on u. "
    r"Neither is stated. The design (conjugation.tex prop:sign (iii)) assumes u invariant under so(V)+so(W_perp); without it the rotation terms contribute (alternatively, sum over an orthonormal basis of S_+(x)W).",
    r"Positivity of the norm, in a representation $\pi$ with $\pi(X^\natural)=\pi(X)^{\rm adj}$ for the conjugation $X^\natural$ in question, in which $-\tfrac12\Delta$ acts as a self-adjoint operator with a positive eigenvalue on a vector annihilated by $\Sigma^*$ and invariant under $\mathfrak{so}(V)\oplus\mathfrak{so}(W_{\mathbb R}^\perp)$, selects the sign $+$; ... (proof) For $Q\in S_+\otimes W$ and a unit vector $u$ annihilated by $\Sigma^*$, invariant under the rotations, with $\Delta$-eigenvalue $-2d$, $d>0$,")

add('m7', 'lem:c8-cc proof, line 363', 'minor',
    r"Two such maps differ by an equivariant automorphism, a positive scalar on each of the irreducible modules $\Sigma$, $\Sigma^*$, fixed to $1$ by involutivity and positivity.",
    r"With independent scalars a on Sigma and b on Sigma^*, involutivity gives b = 1/conj(a) and positivity a > 0, which leaves a > 0 free. "
    r"checks_twocomp T14: c_a = (2 conj on Sigma, (1/2) conj on Sigma^*) is involutive and positive but violates conj(rho(v)Q) = rho(vbar) conj(Q). "
    r"Compatibility with rho, which exchanges Sigma and Sigma^*, forces a = b, and then a = 1. The design's proof has the same gap.",
    r"two such maps differ by a linear automorphism of $\Sigma\oplus\Sigma^*$ commuting with every $\rho(v)$, hence by one scalar $a$ on the irreducible Clifford module $\Sigma\oplus\Sigma^*$ ($\rho(v)$ exchanges $\Sigma$ and $\Sigma^*$, so the scalars on the two summands agree); involutivity gives $|a|=1$ and positivity $a>0$, so $a=1$.")

add('m8', 'thm:c8-L (iv) line 204', 'minor',
    r"Modulo $\hbar$, $\Lambda^{(1)}(x)=\{S^{(1)}(x),-\}$ (on all local functionals).",
    r"By def:c2-dce the (1,0)-shuffle of {S_psl, O} carries (-1)^{|g|k + |x_1|(l+|g|)} with |g| = |O|, k=1, l=0, i.e. (-1)^{|O|(1+|x|)}: for odd O and even x, Lambda^(1)(x)O = -{S^(1)(x),O}. "
    r"The design states the formula for |O| = 0 only. The ghost letter int omega (X,c) (degree 1) is one of the operators of ch06.",
    r"Modulo $\hbar$, $L=s$, $\Lambda^{(0)}=Q+\{I,-\}$, $\Lambda^{(1)}(x)O=(-1)^{|O|(1+|x|)}\{S^{(1)}(x),O\}$, $\Lambda^{(2)}(x_1,x_2)O=(-1)^{|O|(|x_1|+|x_2|)}\{S^{(2)}(x_1,x_2),O\}$, $\Lambda^{(k\ge3)}=0$.")

add('m9', 'def:c8-PT, line 71', 'minor',
    r"For every $T$ the components of $P(\varepsilon,T)$ and of its derivatives are bounded by the corresponding components of $P(0,\infty)$.",
    r"checks_algebra E1/E2: at T = infinity the ratios |d^k P(eps,infty)|/|d^k P(0,infty)| (radial, k <= 4) stay <= 1, but at finite T the second radial derivative exceeds the P(0,infty) one by up to 1.1% "
    r"and the third by 6% (eps = 1, T = 50). Algebraically f'' = e^{-u}(6/r^4)(1 + u + (2/3)u^2), u = r^2/4T, for the T-term at large r. Dominated convergence needs only a constant multiple.",
    r"and for every $T$ they are bounded by a constant multiple, independent of $\varepsilon$ and $T$, of the corresponding components of $P(0,\infty)$.")

add('m10', 'def:c8-quant line 92; remark (higher order) line 520', 'minor',
    r"Counterterms: 'order by order in $\hbar$ and graph by graph, ... coefficients in the span of the negative powers of $\varepsilon$ and of the powers of $\log\varepsilon$ (Definition c4-ict)'.",
    r"def:c4-ict is the one-loop, non-equivariant definition. Beyond one loop, (i) products eps^{-m} log^k eps occur (the chapter's own higher-order remark says so), and they are outside the stated span; "
    r"(ii) 'graph by graph' singular parts of unsubtracted graphs are nonlocal; the counterterms are defined recursively with lower-order counterterm vertices [C07, Thm B].",
    r"... the singular parts of the flow, order by order in $\hbar$ with the counterterms of lower order inserted, local cochains with coefficients in the span of the functions $\varepsilon^{-m}\log^k\varepsilon$, $(m,k)\neq(0,0)$ (Definition \ref{def:c4-ict} at order $\hbar$, [C07, Thm B] at higher order).")

add('m11', 'thm:c8-fields (iii) proof, line 420', 'minor',
    r"At every order the counterterms are reflection invariant: 'the singular part is unique ... (Proposition c4-unique); induction on k'.",
    r"prop:c4-unique is the order-hbar uniqueness statement. Beyond one loop the uniqueness of the counterterms is [C07, Thm B], adapted to R^4 with upper scale T.",
    r"the singular part is unique (Proposition \ref{prop:c4-unique} at order $\hbar$, [C07, Thm B] beyond) and the functions $\varepsilon^{-m}\log^k\varepsilon$ are real; induction on $k$.")

add('m12', 'thm:c8-fields (ii), line 414', 'minor',
    r"$\Theta_n$ is compatible with ... $\{-,-\}_t$ and the flow.",
    r"The proof of thm:c8-adjoint needs the explicit bracket rule {f,g}_t^Theta = (-1)^{|f|+|g|}{g^Theta,f^Theta}_t (conjugation.tex prop:ops), and the anti-multiplicativity (FG)^Theta = G^Theta F^Theta. "
    r"I re-derived the rule from Delta^Theta = (-1)^{|F|}Delta(F^Theta), the anti-multiplicativity and the definition of the bracket. With it {S^(1)(X),a}^Theta = (-1)^{1+|a||X|}{S^(1)(X^dag),a^Theta} follows, as used.",
    r"append to (ii): and $\{f,g\}_t^{\Theta_n}=(-1)^{|f|+|g|}\{g^{\Theta_n},f^{\Theta_n}\}_t$, $(FG)^{\Theta_n}=G^{\Theta_n}F^{\Theta_n}$.")

add('m13', 'prop:c8-C line 433', 'minor',
    r"$\Theta_I$ is the inversion of the fields with the Weyl weights ... and the spinor factor $-\rho(x)/|x|$ on $\lambda$.",
    r"The description omits the factors on the antifields: in conjugation.tex def:Theta, Theta_I has -1 on A^* and on c^*, and +rho(x)/|x| on lambda^*. "
    r"These are needed for Theta_I = (C^{-1})^* Theta_n C^* and for S^(1)(X)^{Theta_I} = (-1)^{|X|}S^(1)(X^dag).",
    r"where $\Theta_I$ is the inversion of the fields with the Weyl weights $1,0,\tfrac32,3,\tfrac52,4$ of $A,c,\lambda,A^*,\lambda^*,c^*$, the spinor factors $-\rho(x)/|x|$ on $\lambda$ and $\rho(x)/|x|$ on $\lambda^*$, and the factors $-1$ on $A^*$ and $c^*$,")

add('m14', 'def:c8-radial, line 446', 'minor',
    r"Inside the definition: for a scalar primary with $D_1(O)=2\gamma_OO$, $K_1(O)=-2(v,x)\gamma_OO$, $O^{\Theta_I,\hbar}=(x,x)^{-(w_O+\gamma_O)}\overline f(I(x))$.",
    r"This is a claim, not a definition, and no proof is given. It needs the rotation components to vanish on O (M_1=R_1=0, as in conjugation.tex cor:weight), holds modulo Lambda^(0)-exact terms, and uses the cocycle v -> (1/4)div v. "
    r"The values are correct: div K_b = -8(b,x); the weight w_O+gamma_O matches D_1 = 2 gamma_O. gamma_O is a new symbol, which DEFINITIONS_v2 item 11 forbids.",
    r"Move the last sentence to a proposition: 'Proposition. Assume moreover $M_1(O)=R_1(O)=0$. For the density $O_x$ at $x$ of a scalar primary $\int\omega f$ of Weyl weight $w_O$ with $D_1(O)=2cO$ and $K_1(O)=\tfrac14(\operatorname{div}K_v)\,cO=-2(v,x)\,cO$, $c$ a number of order one unit, $(O_x)^{\Theta_I,\hbar}=(x,x)^{-(w_O+c)}\overline f(I(x))$ modulo $\Lambda^{(0)}$-exact terms and $O(\hbar^2)$. Proof: $\Lambda^{(1)}(v)$ acts on $O_x$ as the classical action on a density of weight $w_O+c$, since $v\mapsto\tfrac14\operatorname{div}v$ is a $1$-cocycle; the conformal factors of $C^{-1}$, $r_n$, $C$ multiply to that of $I$.'")

add('m15', 'prop:c8-2pt (iii), lines 465-469', 'minor',
    r"Scaling identity and the log form 'for scaling operators of common weight $w$'.",
    r"The log form is correct (checks_algebra J1), and for D_1 = 2 gamma it gives |x|^{-2w}C(1 - gamma log|x|^2). "
    r"It needs Lorentz-scalar operators (translation and rotation covariance give a function of |x-y|) and a finite counterterm of Weyl weight 4 with weight-preserving psi: two_point_design thm:t-scaling assumes F_0 of weight 4. def:c8-quant imposes no weight condition on F.",
    r"(iii) Scaling identity: assume the finite counterterm has Weyl weight $4$ and $\psi$ preserves the Weyl weight. ... so for Lorentz-scalar operators of common weight $w$ with ...")

add('m16', 'prop:c8-Tsafe, line 326', 'minor',
    r"Consequently ... $Q_1$, $S_1$, $D_1$, $P_1$, $M_1$, $R_1$ and the twisted $Q^\xi_1$ may be computed at $T=\infty$.",
    r"(1) Q^xi_1 contains the arity-two and higher components along xi, i.e. graphs with two or more floating S^(1)(xi); the single-insertion decay does not cover them. "
    r"The pair analysis (odd_components_design thm:o-ir (iv): P-edge pair (xi,xi') of degree -1) is needed. "
    r"(2) thm:c8-L (ii) passes to T = infinity only when the graph sums of M[T] and L(M)[T] converge, and these contain the divergent K-components; the conclusion requires the component-wise version (the component along x involves only components of arity <= 1 along x). "
    r"(3) Minor: so(W_perp) has no floating insertion (its A^* leg is transverse and has no propagator), which is why it converges.",
    r"... converge as $T\to\infty$, and so do the arity-two graphs along $(\xi,\xi')$ and $(\xi,\eta)$: for two floating odd insertions joined by a $\lambda$ edge the integrand has degree $-1$ at large distance in $\mathbb R^8$ for $(\xi,\xi')$ and $-1$ after the sum over the two orientations for $(\xi,\eta)$. Since the component of \eqref{eq:c8-rep} along $x$ involves only the components of arity at most one along $x$, Theorem \ref{thm:c8-L} (ii) applies component by component, and $Q_1$, $S_1$, $D_1$, $P_1$, $M_1$, $R_1$ and the twisted $Q^\xi_1$ may be computed at $T=\infty$.")

add('m17', 'remark (comparison with the holomorphic twist), line 315', 'minor',
    r"[BK] is a statement about $[Q^\xi_1]$ transported by [ESW, Thm 10.9]; 'by Theorem c8-twist this class is scheme independent'.",
    r"(1) The transport of [Q^xi_1] to BK's letters by ESW Thm 10.9 is an identification still to be derived: design_audit_corrections design_audits[0] item 5 (UNVERIFIED); ESW organise the fields in their own U(2)-equivariant symbols. "
    r"Agreement also requires the two quantisations (BK: holomorphic gauge, no counterterms) to differ by an element of C S^xi + C N_{R_xi}S^xi + Im Q^xi_0 (design thm:o-twist (iv)). "
    r"(2) 'scheme independent' overstates thm:c8-twist, which covers (a), (b), the coupling direction and G_B.",
    r"The one-loop differential of the holomorphic twist [BK, (3.64)--(3.70)] is to be compared with the class $[Q^\xi_1]$ on the cohomology of $Q^\xi_0$, after its transport to the holomorphic letters through the equivalence of [ESW, Thm 10.9] and provided the two quantisations differ by an element of $\mathbb C\,S^\xi+\mathbb C\,\{S^{(1)}(R_\xi),S^\xi\}+\operatorname{Im}Q^\xi_0$; by Theorem \ref{thm:c8-twist} this class is invariant under the scheme changes listed there, while the individual entries of $Q_1$ are not.")

add('m18', 'lem:c8-ev statement line 284 and proof line 293', 'minor',
    r"'For $\xi=u\otimes w\in S_+\otimes W$ ... every $\xi$ satisfies $[\xi,\xi]=0$'; '$(2c+1)\{S^{(1)}(\xi),S^{(1)}(\xi)\}$, a nonzero functional'.",
    r"Garbled quantifier: every xi in S_+(x)W has [xi,xi] = Gamma|_V = 0, but the ten-component Gamma(xi,xi) vanishes only at rank one (checks_twocomp T19, T20). Rank one is needed for the two-component S^(2)(xi,xi); {S^xi,S^xi} = 0 holds for every xi with [xi,xi] = 0. "
    r"The nonvanishing of {S^(1)(xi),S^(1)(xi)}, needed for 'only c = -1/2', is asserted without argument.",
    r"Every $\xi\in S_+\otimes W$ satisfies $[\xi,\xi]=0$. For $\xi=u\otimes w$ (two-component $\xi_\alpha{}^i=u_\alpha w^i$, $\bar\xi=0$), $\Gamma(\xi,\xi)=0$ in all ten components, and ... ; proof: '... is $(2c+1)\{S^{(1)}(\xi),S^{(1)}(\xi)\}$; this equals $(2c+1)\{S,S^{(2)}(\xi,\xi)\}$, whose part linear in $\lambda$ is $-(2c+1)\int\mathrm{dvol}\,\big((\xi,\lambda^*),(\xi,{\not d}\lambda)\big)$ up to the sign of the pair, with $(\xi,{\not d}\lambda)=i\,u_\alpha w^i\partial_m\bar\lambda_{\dot\alpha i}\bar\sigma^{m\dot\alpha\alpha}\neq0$.'")

add('m19', 'def:c8-renorm, line 167', 'minor',
    r"At $T=\infty$ and $\psi=0$ this is Definition c7-component.",
    r"def:c7-component has no finite counterterm F, while M[T] here contains -F_psl.",
    r"At $T=\infty$, $\psi=0$ and $F_{\mathfrak{psl}(4|4)}=0$ this is Definition \ref{def:c7-component};")

add('m20', 'remark (what is computed), line 243', 'minor',
    r"'... graphs with one defect vertex of Lemma c8-transfer applied to $e^{(J+\delta O(\varepsilon))/\hbar}$, as in Proposition c6-Leps along $\Delta$.'",
    r"prop:c6-Leps is the L_eps boundary form of the scaling route (eps d/deps of the counterterms), not the K_eps defect graphs. The defect formula along Delta is eq:c6-uvobs, in the last section of ch06.",
    r"as in \eqref{eq:c6-uvobs} along $\Delta$.")

add('m21', 'prop:c8-B proof, line 155', 'minor',
    r"Averaging over the compact form of $\mathrm{Spin}(4)\times SU(4)$, $\phi^{(0)}$ is a combination $\sum_fa_f\int(f,f^*)$.",
    r"phi^(0) must also be averaged over the compact form of the gauge group, so that its colour structure is the invariant pairing; the equation G_B = a(S-S^(2)) + s phi is gauge covariant, so this averaging is allowed. "
    r"The remainder (four linear terms with differences a_{lambda,-}-a_{A_V}, a_{lambda,+}-a_{A_W}, a_{A_V}-a_{lambda,+}, a_{A_W}-a_{lambda,-}, summing to 0) is correct.",
    r"averaging over the compact forms of $\mathrm{Spin}(4)\times SU(4)$ and of the gauge group,")

add('m22', 'thm:c8-K, line 340', 'minor',
    r"If a local representative along $K_v$ exists at one finite $T_0$ ...",
    r"Local representatives are defined for the whole cochain (def:c8-rep), and thm:c8-L (ii), used in the proof, is stated for the whole L(M).",
    r"If the local representative $L(O)$ exists at one finite $T_0$, then $K_1(O)$ is the same local functional at every finite $T$, ...")

add('m23', 'def:c8-inner line 479; def:c8-2pt', 'minor',
    r"Radial inner product for $a,b$ supported in the open unit ball; Gram matrix on operators at the origin.",
    r"a^{Theta_I} is compactly supported only if supp a avoids the origin. def:c8-2pt needs compactly supported omega with disjoint supports. Operators at the origin require a limit (the factor (x,x)^{w_O+gamma} at infinity), which is not stated.",
    r"for $a,b$ supported in the open unit ball minus the origin ... On operators at the origin, the radial inner product is the limit, as the support of $a$ shrinks to the origin, written in a basis, a Gram matrix ...")

add('m24', 'whole chapter: notation', 'minor',
    r"Names: '$Q$ alone is the linear BV differential'; spinor parameters are $\xi$, $\eta$.",
    r"(1) lem:c8-cc and thm:c8-dag (iii) use Q, Q_3, Q_4 for spinors; ch02 writes supercharge spinors with subscripts. "
    r"(2) K_0 is the kernel of the identity (line 348, special edge K_eps - K_0) and also the classical matrix of K (lines 504-506). "
    r"(3) w is an element of W_perp in iota_w A (line 34) and of W in xi = u (x) w (line 284); the design reserves w for W. "
    r"(4) v is the vector field of K_b (line 34) and the constant vector labelling K_v, P_v (lines 218, 446). "
    r"(5) New symbols contrary to DEFINITIONS_v2 item 11 and DESIGN_STEP2 R5(b) ('G0, G1 code names; never new symbols in text'): gamma_O, O_1/O_2/O_a, G_tree, G_1, X_0, X_1, R_n.",
    r"Write $\xi,\xi_3,\xi_4$ for the spinors in Lemma \ref{lem:c8-cc} and Theorem \ref{thm:c8-dag} (iii); write the classical matrices as 'the matrix of $-\{S^{(1)}(X^\dagger),-\}$'; write $w'$, or 'the element $w\in W$', in $\xi=u\otimes w$; write $K_b$ with $b\in V$ throughout ($K_1(O):=\hbar L_1(O)^{(1)}(K_b)$, $-2(b,x)$); write the tree and one-loop Gram matrices in words.")

add('m25', 'citations [BK] (line 315) and [Srednicki] (lines 12, 47)', 'minor',
    r"[BK, (3.64)--(3.70)], [Srednicki, §§34--36], [Srednicki, (69.16), (95.77)].",
    r"Neither key occurs in ch00-ch07 or in the source ledgers, and the manuscript has no reference list; the other keys are introduced in ch00/ch01 by author. "
    r"BK = K. Budzik, J. Kulp, arXiv:2512.07771 (literature/gkw.tex). The equation range (3.64)-(3.70), with (3.70) = (1.10), agrees with q1_targets.json.",
    r"at first use: '[BK] (K. Budzik and J. Kulp, arXiv:2512.07771)', '[Srednicki] (M. Srednicki, Quantum Field Theory, Cambridge University Press, 2007)'.")

add('m26', 'opening paragraph (line 3) and defensive wording (lines 3, 114, 151, 155, 289, 293, 315, 340, 463)', 'minor',
    r"FORMAT.md: chapters open with 2-4 sentences; no defensive wording ('only', 'does not', 'not a ...').",
    r"The opening paragraph has five sentences. Defensive phrases: 'entering only through' (3), 'does not involve T' (114), 'for no y' (151), 'since no local functional ... exists' (155), 'only c = -1/2 does' (289), 'contains only [xi,xi]=0' (293), "
    r"'is never used as an input' (315), 'enters no local representative' (340), 'do not contribute' (463).",
    r"e.g. (3) 'with the global ghosts of $\mathfrak{psl}(4|4)$ entering through the classical cochain $S^{(1)}+S^{(2)}$'; (151) '$B$ is an outer derivation of $\mathfrak{psl}(4|4)$'; (289) 'and $c=-\tfrac12$ is the unique value for which $S+S^{(1)}(\xi)+c\,S^{(2)}(\xi,\xi)$ solves it'; (315) 'The twisted theory appears as a consequence:'; (340) 'such a term is part of the flow of Lemma \ref{lem:c8-flow}, and the local representative is $T$-independent'; merge sentences 4 and 5 of the opening.")

add('m27', 'thm:c8-adjoint hypotheses, line 484', 'minor',
    r"'the two-point cochain satisfies the Ward identity of the equivariant quantum master equation' (for every $X$).",
    r"The proof uses the identity <Lambda^(1)(Y)F,G> = -(-1)^{|Y||F|}<F,Lambda^(1)(Y)G> for Y = X^{dag_n}, at T = infinity, for all generators. "
    r"two_point_design thm:t-ward derives it only for x in Sigma+Sigma^*, Delta, P, so(V), so(W_perp) (large-distance hypothesis H3), and the K-direction via the (eta_1,eta_2) arity-two identity, whose graphs diverge logarithmically for v_1 != 0 (design_audits[1]). "
    r"As a hypothesis this is acceptable, but it should be stated as the explicit identity.",
    r"... and the two-point function satisfies, for every $Y$ in $\mathfrak{psl}(4|4)$, the Ward identity $\langle\Lambda^{(1)}(Y)F,G\rangle=-(-1)^{|Y||F|}\langle F,\Lambda^{(1)}(Y)G\rangle$ on the observables considered.")

# ------------------------------------------------------------------ markup
add('k1', 'lines 108, 109, 130-133, 201-204, 228-229, 252-254, 266-268, 377, 382, 383, 413-415, 463-465 (27 items)', 'markup',
    r"\item[(i)] ... (optional argument of \item).",
    r"FORMAT.md allows '\begin{itemize} \item ... \end{itemize}' only; ch00-ch07 never use \item[...]. "
    r"reader_build/convert.py splits on r'\\item\b' and leaves the optional argument, so the reader shows '* [(i)] ...'.",
    r"\item (i) ...   (label as text at the start of the item; same for (ii), (iii), (iv), (a), (b), (c)).")

checked_ok = [
    'All 70 \\ref/\\eqref labels resolve (ch01-ch08); no duplicate labels; all LaTeX commands are KaTeX-supported (only \\item[...] is outside FORMAT).',
    'Def c8-twocomp: Wick map v^0=-iv^4, x^4=it, d_0=i d_4; sigma^m=(sigma,i), sigmabar^m=(-sigma,i); sigma^mu v_mu = sigma^m v_m; {sigma,sigmabar}=-2 delta; epsilon conventions; (R2); Gamma|_V on Sigma x Sigma and Sigma^* x Sigma^* (T6, T7); Killing spinor eq:c8-killing (T8); with the dictionary W_perp data: Clifford, symmetry, ESW Assumption 3.1, rho(e1e2e3e4)=+1 on S_+(x)W (T1-T3, T9).',
    'Prop c8-twocomp: kinetic action (-1/8(F,F), -1/16(DA^ij,DA_ij), i(lambdabar,sigmabar D lambda), (Dc,A^*)), (X,Y)=2hbar X^aY^a, F=-igF^aT^a, S/hbar = Wick image of Srednicki, e^{iS_M}=e^{S/hbar}; Dirac column rho(v)=-i gamma_S v = gamma_W v; propagator components 2,4 eps_ijkl,-i p_m,-p.sigma (dictionary prop:prop).',
    'Def c8-PT scalar factors; Lemma c8-transfer both identities and proof (Delta=-d_K, [d_P,Q]=d_{QP}); Prop c8-Tindep (ii).',
    'Def c8-quant and Thm c4-finite: F enters as I - F; the 27/16 coefficient and sign agree with Step 1 (s F = Ob at arity zero).',
    'Thm c8-obstruction (iv): the family mu S + S^(1) + mu^{-1} S^(2) solves eq:c2-mc (arity-n component scales as mu^{2-n} by lem:c2-components, S^(k>=3)=0); derivative S - S^(2).',
    'B: eigenvalues +1 on S_+(x)W, S_-(x)W and -1 on S_-(x)W^*, S_+(x)W^*; two-component (upper minus lower W indices); brackets inside each eigenspace vanish (T18); outerness; G_B components; G_B s-closed, degree 0, weight 4; independence from the coupling class (weight count: no degree -2, weight 9/2 functional; four differences of a_f summing to 0).',
    'Prop c8-laws: (a) [s,ad_phi] = -ad_{s phi} (checks_bracket B4) and conjugation g^{-1}Lg; (b) L_1 -> L_1 + [s,psi]; (c) L_1 -> L_1 - a ad_G, coupling components (0: -a(Q+{I,-}), 1: none, 2: +a{S^(2),-}), G_B arity one -a b{S^(1)(x),-}; the Q_1/S_1 shifts on the four blocks.',
    'Thm c8-L (i)-(iii) arguments; Lemma c8-flow ((d+Q+hbar Delta_T)(delta M e^{I/hbar}) = +-delta (d_T M) e^{I/hbar}).',
    'Cor c8-homotopy (ii): the displayed identity holds for all parities (checks_bracket B2, B3: ad_{f,g}=(-1)^{|f|+1}[ad_f,ad_g] for the ch01 bracket; B5: the exact term is [ad_S, ad_{S^(2)}]); x -> -Lambda^(1)(x) is a representation on cohomology.',
    'Thm c8-invariants (i), (ii): conjugations, coupling direction ((1-a hbar) Lambda^(0)), composition with e^{-a hbar B}, factors 1 -+ a hbar on the blocks.',
    'Lemma c8-ev: shuffle sign (-1)^{kl}; a_k a_l = a_{k+l} binom(k+l,k)(-1)^{kl}, a_0..a_4 = 1,1,-1/2,-1/6,1/24 (I1); c = -1/2 via {S,S^(2)(xi,xi)} = {S^(1)(xi),S^(1)(xi)} and {S^(1)(xi),S^(2)(xi,xi)}=0; S^(2)(xi,xi) = -1/2 int(u w lambda^*, u w lambda^*) = (F3) of the dictionary; rank-one condition (T19, T20).',
    'Def c8-twisted coefficients; Thm c8-twist: r_xi(S - S^(2)) = S^xi - S^(1)(xi) + S^(2)(xi,xi), r_xi(G_B) = S^(1)(xi) - S^(2)(xi,xi), N S^xi = -S^(1)(xi) + S^(2)(xi,xi) from lem:c2-components (arity 2 and 3 with (R_xi,xi,xi)); the resulting shifts -aQ_0^xi + [Q_0^xi, a N] and [Q_0^xi, -aN].',
    'Remark c8-twist two-component S^(1)(xi): all five terms agree with (F1) of the dictionary specialised to xi = u (x) w (by hand); index-free form = eq:c2-S1Q.',
    'Prop c8-Tsafe: decay |y|^{-6}, |y|^{-5}, |y|^{-4}; int_{|y|^2<T} |y|^{-4} d^4y = pi^2 log T; [K_v,P_v\'] = (v,v\')Delta + rotation (ch07 rem:c7-KDelta).',
    'Def c8-K, Thm c8-K and Remark c8-K are consistent with directive R7 (K_1 at finite T, T-independent; K = P^dag as cross-check).',
    'Lemma c8-cc: two-component conjugation, conj(rho(e_M)Q) = rho(e_M) conj(Q) for all ten real e_M, involutive, (xi, conj xi) = sum |xi|^2, Gamma(conj,conj) = conj Gamma (T10-T13); chirality mapping.',
    'Def c8-dag / Thm c8-dag (i): (xi + rho(x) eta)^dag = -conj eta - rho(x) conj xi; P_b^dag = -K_{bbar}, K^dag = -P, Delta^dag = Delta, omega^dag = -conj omega, xi^dag = -conj xi, eta^dag = -conj eta; two-component and unit elements (Q^i_alpha)^dag = -S^alpha_i, (Qbar)^dag = +Sbar, (S)^dag = -Q, (Sbar)^dag = +Qbar (T15, T16); matches conjugation.tex prop:dag2 (ip_targets IP7 uses the label convention Q^dag = S, a renaming S -> -S).',
    'Thm c8-dag (ii): X -> -X^dag (even), -i X^dag (odd) is an involutive antilinear automorphism given [X,Y]^dag = [Y^dag,X^dag]; fixed even part and maximal compact so(4)+R i Delta+so(6). (iii): the family e^{tB}((+-1)^{|X|}X^dag), t real, and the sign + from positivity (up to m6).',
    'Def c8-Theta: (xi + rho(x)eta)^{dag_n} = rho(n) conj xi - rho(x) rho(n) conj eta; two-component xi^{dag_n} (T17), Theta_n lambda, Theta_n lambdabar, Theta_n lambda^*, Theta_n lambdabar^*, Theta_n A_ij = conj A^ij (by hand from (R1),(R2) and the conjugation).',
    'Thm c8-fields (i)-(iii): statements agree with conjugation.tex thm:fields, prop:ops, thm:inv; 27/16 term reflection invariant.',
    'Prop c8-C: C is the time -pi/4 flow of P_n - K_n (D1), maps the half space into the unit ball (D2), C r_n C^{-1} = I (D3), Ad(C) = C_* = exp(+pi/4 ad_{P_n-K_n}) (D4); proof of X^dag = Ad(C)((Ad(C)^{-1}X)^{dag_n}) on the even part.',
    'Def c8-radial: K_1 = (1/4) div K_v gamma O = -2(v,x) gamma O and the weight (x,x)^{-(w+gamma)}; order hbar^0 = O^{Theta_I}.',
    'Prop c8-2pt (ii), (iii): existence and large-distance convergence (sigma > 0 for every vertex set), scaling identity and log form (J1); D_1 = -48 units for the trace bilinear = 2 x (-24) of thm:c6-dil; unit = -g^2N/(32 pi^2), gamma = 3g^2N/(4 pi^2) (K1, su(2), su(3), su(4)).',
    'Thm c8-adjoint: sign chain of the proof (both identities) given the bracket rule of m12 and the Ward identity; the second identity via e^{t delta} Lambda(X) e^{-t delta} = Lambda(e^{-t ad_Y}X) and Ad(C)^{-1} = exp(-pi/4 ad).',
    'Cor c8-gram: the order-hbar relation (X^dag)_1 = (X_1)^adj + [(X^dag)_0, G^{-1}G_1] (C1); the outer direction forces conj(a) = -a (C9) and G_B is anti-invariant under the conjugations (T22); K_1 = [K_0, G^{-1}G_1] when P_1 = 0; D_1 tree self-adjoint since D_0 is a multiple of the identity at fixed weight.',
]

out = dict(
    file='manuscript/ch08_definitions.tex',
    date='2026-10-28',
    summary=dict(fatal=0, major=sum(f['severity'] == 'major' for f in F), minor=sum(f['severity'] == 'minor' for f in F),
                 markup=sum(f['severity'] == 'markup' for f in F)),
    checks=['step2/qa/ch08_review/checks_algebra.py -> checks_algebra.json',
            'step2/qa/ch08_review/checks_bracket.py -> checks_bracket.json',
            'step2/qa/ch08_review/checks_twocomp.py -> checks_twocomp.json'],
    findings=F,
    checked_correct=checked_ok,
)
json.dump(out, open(os.path.join(HERE, 'findings.json'), 'w'), indent=1, ensure_ascii=False)
print(json.dumps(out['summary']))
