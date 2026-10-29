#!/usr/bin/env python3
"""Consolidated Chapter 8 revision: list of exact (old, new) replacements.
Each entry: id, target label, old, new, sources (list of 'file#index').
apply.py applies them in order to manuscript/ch08_definitions.tex, asserting each old occurs exactly once."""

E = []


def add(eid, label, old, new, sources):
    E.append(dict(id=eid, label=label, old=old, new=new, sources=sources))


D = 'step2/qa/ch08_corrections_from_delta.json'
G = 'step2/qa/ch08_corrections_from_design.json'
P = 'step2/qa/ch08_corrections_from_ip.json'
C = 'step2/qa/ip_completeness2.json'

# ---------------------------------------------------------------- Def c8-PT, Lemma c8-transfer
add('R01', 'def:c8-PT',
    r"""every component of $P(\varepsilon,T)$ and of its derivatives decays as $e^{-|x-y|^2/4T}$, and for every $T$ they are bounded by a constant multiple, independent of $\varepsilon$ and $T$, of the corresponding components of $P(0,\infty)$.
\end{definition}""",
    r"""every component of $P(\varepsilon,T)$ and of its derivatives decays as $e^{-|x-y|^2/4T}$, and for every $T$ they are bounded by a constant multiple, independent of $\varepsilon$ and $T$, of the corresponding components of $P(0,\infty)$. Here $\partial_\phi$ is read as in Definition \ref{def:c3-P}, its value on the product of a coordinate functional of the first factor and one of the second, in this order, being the corresponding component of $\phi$, and a linear operator $X$ on $\mathscr E$ acts on functionals as the derivation induced by $\ell\mapsto\ell\circ X$ [C07, §6.2]; then $[\partial_\phi,X]=\partial_{(X\otimes1+1\otimes X)\phi}$ as a graded commutator for kernels and operators of either parity (Lemma \ref{lem:c9-ops}).
\end{definition}""",
    [D + '#17', D + '#18'])

add('R02', 'lem:c8-transfer (proof)',
    r"""By [C07, §6.4], $[\partial_\phi,Q]=\partial_{Q\phi}$; with \eqref{eq:c8-PT},""",
    r"""By [C07, §6.4] and Definition \ref{def:c8-PT}, $[\partial_\phi,Q]=\partial_{(Q\otimes1+1\otimes Q)\phi}$; with \eqref{eq:c8-PT},""",
    [D + '#17', D + '#18'])

# ---------------------------------------------------------------- Def c8-quant (weight)
add('R03', 'def:c8-quant',
    r"""A finite counterterm is a local cochain $F_{\mathfrak{psl}(4|4)}=\hbar F_1+\hbar^2F_2+\cdots$ of degree $0$, independent of $\varepsilon$.""",
    r"""A finite counterterm is a local cochain $F_{\mathfrak{psl}(4|4)}=\hbar F_1+\hbar^2F_2+\cdots$ of degree $0$ and Weyl weight $4$, independent of $\varepsilon$; here the explicit position has weight $-1$ and an argument $x$ has weight one half of its $\mathrm{ad}_\Delta$-eigenvalue (Remark \ref{rem:c2-dilatation}), so that a value on $x_1,\dots,x_k$ has weight $4$ minus the weights of the arguments (Proposition \ref{prop:c7-select}).""",
    [G + '#8'])

# ---------------------------------------------------------------- Thm c8-obstruction (iv) proof
add('R04', 'thm:c8-obstruction (proof of (iv))',
    r"""in arity zero $(s\phi)^{(0)}=\{S,\phi^{(0)}\}$, with $\phi^{(0)}$ of degree $-1$ and Weyl weight $4$ a combination of the field-counting functionals (Lemma \ref{lem:c4-enumeration}), and $\{S,\phi^{(0)}\}=S$ is excluded by the computation of Corollary \ref{cor:c4-coupling}.""",
    r"""in arity zero $(s\phi)^{(0)}=\{S,\phi^{(0)}\}$; $\{S,-\}$ preserves the degree in the explicit position, so $\{S,\phi^{(0)}\}=S$ requires the same identity for the part of $\phi^{(0)}$ of degree zero in the position, which is translation invariant of degree $-1$ and Weyl weight $4$, a combination of the field-counting functionals (Lemma \ref{lem:c4-enumeration}), and that identity is excluded by the computation of Corollary \ref{cor:c4-coupling}.""",
    [G + '#5'])

# ---------------------------------------------------------------- Prop c8-B proof
add('R05', 'prop:c8-B (proof)',
    r"""Independence: if $G_B=a(S-S^{(2)})+s\phi$, the arity-one component along $\xi\in S_+\otimes W$ reads $S^{(1)}(\xi)=\pm\{S^{(1)}(\xi),\phi^{(0)}\}$, since the space of local functionals of degree $-2$ and Weyl weight $\tfrac92$ is zero;""",
    r"""Independence: suppose $G_B=a(S-S^{(2)})+s\phi$ with $\phi$ local of degree $-1$. Along a translation, $0=G_B^{(1)}(P_b)=(s\phi)^{(1)}(P_b)$, so $\{S^{(1)}(P_b),\phi^{(0)}\}$ is $(Q+\{I,-\})$-exact; since $\{S^{(1)}(P_b),-\}$ acts on the explicit polynomial dependence on the position as $-\partial_b$ and $Q+\{I,-\}$ preserves the degree in the position, induction on that degree replaces $\phi$ by $\phi-s\chi$, with $\chi$ local of degree $-2$ and $\phi^{(0)}$ translation invariant, and $s\phi$ is unchanged; this excludes $\phi^{(0)}=S^{(1)}(\Delta)$, which is invariant under the compact groups and satisfies $\{S^{(1)}(\xi),S^{(1)}(\Delta)\}=S^{(1)}(\xi)$ for all $\xi\in\Sigma$, while $\{S^{(1)}(P_b),S^{(1)}(\Delta)\}=-2S^{(1)}(P_b)$. The part of degree zero in the position of the arity-one component along $\xi\in S_+\otimes W$ then reads $S^{(1)}(\xi)=\pm\{S^{(1)}(\xi),\phi^{(0)}\}$, since the space of translation-invariant local functionals of degree $-2$ and Weyl weight $\tfrac92$ is zero;""",
    [G + '#4'])

# ---------------------------------------------------------------- Def c8-renorm (weight)
add('R06', 'def:c8-renorm',
    r"""A scheme for local operators is a linear map $\psi$ of local cochains, preserving the degree and the Weyl weight,""",
    r"""A scheme for local operators is a linear map $\psi$ of local cochains, preserving the degree and the Weyl weight (counted as in Definition \ref{def:c8-quant}),""",
    [G + '#8'])

# ---------------------------------------------------------------- Thm c8-L (iv) and (ii) proof
add('R07', 'thm:c8-L (iv), product rule',
    r"""the sign rule of Definition \ref{def:c2-dce}; the same rule defines the products of field-independent cochains with local cochains.""",
    r"""the sign rule of Definition \ref{def:c2-dce}; the same rule defines the products of field-independent cochains with local cochains, and, with the composition replaced by the product of values and $|\Lambda'|$ read as the total parity (degree plus fermion number modulo $2$), the product of any two functional-valued cochains, in particular $e^{J/\hbar}$ in \eqref{eq:c8-IT}; this product is associative and graded commutative.""",
    [G + '#6', C + '#7'])

add('R08', 'thm:c8-L (iv), sign of the components',
    r"""$\Lambda^{(1)}(x)O=(-1)^{|O|(1+|x|)}\{S^{(1)}(x),O\}$, $\Lambda^{(2)}(x_1,x_2)O=(-1)^{|O|(|x_1|+|x_2|)}\{S^{(2)}(x_1,x_2),O\}$, $\Lambda^{(k\ge3)}=0$.""",
    r"""$\Lambda^{(1)}(x)O=(-1)^{|O|(1+|x|)}\{S^{(1)}(x),O\}$, $\Lambda^{(2)}(x_1,x_2)O=(-1)^{|O|(|x_1|+|x_2|)}\{S^{(2)}(x_1,x_2),O\}$, $\Lambda^{(k\ge3)}=0$, where in these two formulas $|O|$ is the degree of $O$ in the sense of Definition \ref{def:c2-dce} (its ghost number), so that $\Lambda^{(1)}(x)=\{S^{(1)}(x),-\}$ on functionals of degree $0$ of either parity and $x\mapsto-\Lambda^{(1)}(x)$ is a representation on the cohomology also on odd operators.""",
    [P + '#2', C + '#7', G + '#6'])

add('R09', 'thm:c8-L (proof of (ii))',
    r"""and the flow intertwines $d_{T_0}$ with $d_T$; for $T<T_0$ the flow is inverted.""",
    r"""and the flow intertwines $d_{T_0}$ with $d_T$; for $T<T_0$ the flow is inverted. At $T=\infty$ the flow with $P(T_0,\infty)$ is the limit of the convergent graph sums, and the intertwining is the first identity of \eqref{eq:c8-transfer} with $(T_0,\infty)$ in place of $(\varepsilon,T)$, in which $\Delta_\infty=0$ (Definition \ref{def:c3-delta}).""",
    [D + '#27', D + '#33'])

# ---------------------------------------------------------------- Remark What is computed
add('R10', 'Remark (What is computed)',
    r"""Chapter 6 computed the dilatation coefficients through the scaling covariance of the regulator, which give $D_1$ under the identification \eqref{eq:c6-step2}: for the traceless scalar bilinears they vanish, and for the trace bilinear $D_1(O)$ is $-48$ units times $O$, the value $2\times(-24)$ of Theorem \ref{thm:c6-dil}.""",
    r"""Chapter 6 computed the dilatation coefficients through the scaling covariance of the regulator. On the part of the scalar bilinears with two legs $A|_{W_{\mathbb R}^\perp}$, Theorem \ref{thm:c9-D1} gives $D_1(O)$ equal to twice the coefficient of Theorem \ref{thm:c6-dil} times $O$ plus $-(\tfrac13+\theta)$ units times $\{S^{(1)}(\Delta),O\}$, where $\theta$ units is the coefficient of $\int(2E(A|_{W_{\mathbb R}^\perp})+2A|_{W_{\mathbb R}^\perp},A^*|_{W_{\mathbb R}^\perp})$ in $\hbar F_1^{(1)}(\Delta)$, a scheme parameter that the two-point part of $sF_1=\mathrm{Ob}$ leaves free. When the one-loop quantisation and the local representatives exist through arity two along $(\Delta,P_b)$, the component of Corollary \ref{cor:c8-homotopy} (ii) along $(\Delta,P_b)$ gives $\theta=-\tfrac13$ (Theorem \ref{thm:c9-kappa}); then $D_1(O)$ vanishes for the traceless scalar bilinears and is $-48$ units times $O$ for the trace bilinear, the value $2\times(-24)$ of Theorem \ref{thm:c6-dil}, and \eqref{eq:c6-step2} holds on that part (Corollary \ref{cor:c9-ident}).""",
    [D + '#12', D + '#6', D + '#23', D + '#28'])

# ---------------------------------------------------------------- Prop c8-laws (a)
add('R11', 'prop:c8-laws (a)',
    r"""\item (a) Field redefinition: $F_{\mathfrak{psl}(4|4)}\mapsto F_{\mathfrak{psl}(4|4)}+\hbar\,s\phi$, $\phi$ local of degree $-1$:""",
    r"""\item (a) Field redefinition: $F_{\mathfrak{psl}(4|4)}\mapsto F_{\mathfrak{psl}(4|4)}+\hbar\,s\phi$, $\phi$ local of degree $-1$ and Weyl weight $4$:""",
    [G + '#8'])

# ---------------------------------------------------------------- Lemma c8-ev
add('R12', 'lem:c8-ev',
    r"""and the product of operator-valued cochains of Theorem \ref{thm:c8-L} (iv) with the composition of their values.""",
    r"""and the products of operator-valued and of functional-valued cochains of Theorem \ref{thm:c8-L} (iv) with the composition, resp. the product, of their values; hence $r_\xi W(P(\varepsilon,T),J)=W(P(\varepsilon,T),r_\xi J)$.""",
    [G + '#7'])

# ---------------------------------------------------------------- Remark twist
add('R13', 'Remark (Comparison with the holomorphic twist)',
    r"""provided the two quantisations differ by an element of $\mathbb C\,S^\xi+\mathbb C\,\{S^{(1)}(R_\xi),S^\xi\}+\operatorname{Im}Q^\xi_0$""",
    r"""provided the two quantisations differ by an element of $\mathbb C\,S^\xi+\operatorname{Im}Q^\xi_0$, which contains $\{S^{(1)}(R_\xi),S^\xi\}=Q^\xi_0S^{(1)}(R_\xi)$ and $r_\xi(G_B)|_{q=1}=-Q^\xi_0S^{(1)}(R_\xi)$""",
    [G + '#9'])

# ---------------------------------------------------------------- Prop c8-Tsafe statement
add('R14', 'prop:c8-Tsafe (statement)',
    r"""whose bracket is a special conformal field.
\end{proposition}""",
    r"""whose bracket is a special conformal field. Two floating insertions of even generators are adjacent through an $A|_V\otimes A|_V$ edge; the chain has both ends on ghost legs, its degree at large distance in $\mathbb R^8$ is $d_1+d_2-2$ for one orientation ($d=0,1,1,2$ for translations, $\mathfrak{so}(V)$, $\Delta$ and the special conformal fields) and one less after the sum over the two orientations. The arity-two graphs along pairs of translations, rotations of $V$ and $\Delta$ converge after this sum, and so do those along $(\eta_1,\eta_2)$ with $\Gamma(\eta_1,\eta_2)|_V=0$; those containing a special conformal field are defined at finite $T$. The coefficient of $\log T$ of a graph with $S^{(1)}(K_b)$ floating on its cycle is $\tfrac{3}{8\pi^2}\hbar$ times the tree obtained by cutting the insertion out of the cycle and evaluating the freed legs on the constant modes $A|_V=b\otimes X_a$ and $c=X^a$, summed over a basis $\{X_a\}$ of $\mathfrak g$ with dual basis $\{X^a\}$ for $(-,-)$. Along $\Delta$, the translations and $\mathfrak{so}(V)$ the insertions on a cycle are composite $A|_V$--$c$ lines as well; their graphs converge at $T=\infty$, single composite graphs have nonlocal $Q$-images, and the components are read off from the complete sets of composite graphs.
\end{proposition}""",
    [G + '#10'])

# ---------------------------------------------------------------- Prop c8-Tsafe proof
add('R15', 'prop:c8-Tsafe (proof, edge decay)',
    r"""The coefficient of the insertion grows as $1$, $|y|$, $|y|^2$ for constant, linear and quadratic vector fields or Killing spinors, and the two adjacent edges decay as $|y|^{-3}$ each (one derivative, Proposition \ref{prop:c3-P});""",
    r"""On a cycle the terms of $S^{(1)}(x)$ with both legs on components of the propagator are $(v.A,A^*)$ for even $x$, with $v.A=L_vA+\tfrac14(\operatorname{div}v)A|_{W_{\mathbb R}^\perp}$ (Definition \ref{def:c7-vdot}), and $-(\Gamma(\epsilon(y),\lambda),A^*)$ for odd $x$, $\epsilon(y)=\xi+\rho(y)\eta$. The $A^*|_V$ leg meets an $A^*|_V\otimes c$ component and the $\lambda$ leg a $\lambda\otimes\lambda$ component, both of order $|y|^{-3}$; the $A$ leg meets an $A\otimes A$ component, of order $|y|^{-2}$ (Proposition \ref{prop:c3-P}). A vector field of degree $k$ enters through the derivative along $v$, of order $|y|^k$ acting on the $A\otimes A$ component, and through $dv$ and $\operatorname{div}v$, of order $|y|^{k-1}$, so the insertion is of order $|y|^{k-6}$; a Killing spinor of degree $g$ gives $|y|^{g-6}$;""",
    [D + '#21', D + '#16', D + '#10'])

add('R16', 'prop:c8-Tsafe (proof, arity two and K_b)',
    r"""for $(\xi,\eta)$ the two orientations traverse the $\lambda$ edge in opposite directions and their leading terms cancel. For $K_b$ the product is of order $|y|^{-4}$ with a nonzero average over the sphere, and $\int_{|y|^2<T}|y|^{-4}d^4y$ grows as $\pi^2\log T$.""",
    r"""for $(\xi,\eta)$ the two orientations traverse the $\lambda$ edge in opposite directions and their leading terms cancel. For two floating even insertions the $A\otimes A$ edge between them carries both derivatives and is of order $|y|^{-4}$, the two $A^*|_V\otimes c$ edges are of order $|y|^{-3}$, and the two orientations of the chain enter with opposite signs, so their leading terms cancel. For $K_b$ the product is of order $|y|^{-4}$ with a nonzero average over the sphere; this term comes from the antisymmetric part of $dv$ acting on $A|_V$, while on $A|_{W_{\mathbb R}^\perp}$ the leading terms of the derivative along $v$ and of $\tfrac14\operatorname{div}v$ cancel. $\int_{|y|^2<T}|y|^{-4}d^4y$ grows as $\pi^2\log T$; $T\partial_T$ of the scalar integral of the composite line tends to $-\tfrac{3}{16\pi^2}b$, and with the fibre factors $2(-,-)^{-1}$ of $A\otimes A$ and $-(-,-)^{-1}$ of $A^*|_V\otimes c$ (Proposition \ref{prop:c8-twocomp}) it gives the coefficient $\tfrac{3}{8\pi^2}$ of the cut tree. After the integration over the insertion point the two edges form the kernel of $P(\varepsilon,T)\star\,v.\,P(\varepsilon,T)\star$, with $P(\varepsilon,T)\star=\int_\varepsilon^TQ^{GF}e^{-tH}dt$ the operator represented by $P(\varepsilon,T)$; along $\Delta$ it is of order $|p|^{-3}$ at small momentum ($|x-y|^{-1}$ in position space at $T=\infty$), and along translations of order $|p|^{-2}$ (Proposition \ref{prop:c9-composite} (iv)), so the one-loop integrals with external legs at generic momenta converge.""",
    [D + '#21', G + '#10'])

# ---------------------------------------------------------------- Def c8-K, Thm c8-K, Remark existence along K
add('R17', 'def:c8-K',
    r"""the arity-two components along $(\eta_1,\eta_2)$ and $(K_b,P_{b'})$ are defined in the same way at finite $T$.
\end{definition}""",
    r"""the arity-two components along $(\eta_1,\eta_2)$ and $(K_b,P_{b'})$ are defined in the same way at finite $T$. The limit $\varepsilon\to0$ is taken at fixed $T$; at a fixed $\varepsilon$ the graphs of the defect carry a term of order $\varepsilon\log T$ for the scalar bilinears and along $(\eta_1,\eta_2)$ with $\Gamma(\eta_1,\eta_2)|_V\neq0$, so the defect graphs are evaluated at finite $T$.
\end{definition}""",
    [G + '#11'])

add('R18', 'thm:c8-K',
    r"""The observables $O[T]$ themselves may keep a term in $\log T$ in their component along $K_b$; such a term is part of the flow of Lemma \ref{lem:c8-flow}, and the local representative is independent of $T$.""",
    r"""The observables keep, in their component along $K_b$, the term $\tfrac{3}{8\pi^2}\hbar\log T$ times $\sum_a\partial_{A|_V=b\otimes X_a}\partial_{c=X^a}$ of the order-$\hbar^0$ part of $O[\infty]$, with $\{X_a\}$ a basis of $\mathfrak g$ and $\{X^a\}$ the dual basis for $(-,-)$; this tree is nonzero for the scalar bilinears, the term is generated by the flow of Lemma \ref{lem:c8-flow}, its image under $Q$ vanishes for gauge-invariant $O$, and the local representative is independent of $T$. The effective interaction $I_{\mathfrak{psl}(4|4)}[T]$ keeps the same term with the tree of $W(P(0,\infty),I)$, which is nonzero in every two-point sector.""",
    [G + '#12'])

add('R19', 'Remark (The existence along K)',
    r"""In the graphs of the defect along $K_b$ the logarithm in $T$ multiplies a tree that contains the special edge $K_\varepsilon-K_0$ at a generic momentum $p$, and therefore carries the factor $e^{-\varepsilon|p|^2}-1$; it vanishes as $\varepsilon\to0$ at fixed $T$. The existence of the local representative at fixed $T$ is thus the locality of the $\varepsilon\to0$ limit of these graphs, the same statement as along the other generators. Two consistency relations constrain $K_1$: the component of $L\circ L=0$ at $(K_b,P_{b'})$, which ties $K_1$ to $D_1$ through $[K_b,P_{b'}]=(b,b')\Delta+(\text{rotation})$ (Remark \ref{rem:c7-KDelta}), and the relation $K=P^\dagger$ of Corollary \ref{cor:c8-gram}.""",
    r"""In a graph of the defect along $K_b$ with $S^{(1)}(K_b)$ floating on its cycle, the coefficient of $\log T$ at a fixed $\varepsilon$ is $\tfrac{3}{8\pi^2}\hbar$ times the cut tree of Proposition \ref{prop:c8-Tsafe}; this tree contains the special edge $K_\varepsilon-K_0$ at a generic momentum $p$ and carries the factor $e^{-\varepsilon|p|^2}-1$. Summed over the graphs, the coefficient is $\tfrac{3}{8\pi^2}\hbar$ times $Q$ of $\sum_a\partial_{A|_V=b\otimes X_a}\partial_{c=X^a}$ of the order-$\hbar^0$ part of $\partial_\delta W(P(\varepsilon,\infty),I+\delta O)$: it vanishes identically in the two-point sectors of the action and is of order $\varepsilon$, and nonzero, for the scalar bilinears. Along $(\eta_1,\eta_2)\in\Sigma^*\times\Sigma^*$ with $\Gamma(\eta_1,\eta_2)|_V\neq0$ the coefficient multiplies the first moment of the cut tree and carries the gradient $-2\varepsilon p\,e^{-\varepsilon|p|^2}$ of the special-edge factor. In every case it vanishes as $\varepsilon\to0$ at fixed $T$. At a fixed $T$, the derivative $T\partial_T$ of every graph with the special edge on its cycle is of order $\varepsilon$, because the edge pinned at $T$ confines the loop momentum and the special edge carries $e^{-\varepsilon|q|^2}-1$; the limit $\varepsilon\to0$ is therefore the same at every $T$, and the existence of the local representative at fixed $T$ is the locality of that limit, the same statement as along the other generators. Every insertion of $S^{(1)}(K_b)$ on a cycle is a composite $A|_V$--$c$ line with kernel $P(\varepsilon,T)\star\,v.\,P(\varepsilon,T)\star=v.\,(P(\varepsilon,T)\star)^2+[P(\varepsilon,T)\star,v.]\,P(\varepsilon,T)\star$, $v$ the vector field of $K_b$; the first term exists because $(Q^{GF})^2\neq0$, and the existence of $L(O)$ along $K_b$ is a statement about the complete set of composite graphs (Proposition \ref{prop:c8-Tsafe}). Three relations constrain $K_1$ on the cohomology of $\Lambda^{(0)}$: the component of $L\circ L=0$ at $(K_b,P_{b'})$, which ties $K_1$ to $D_1$, $P_1$ and $M_1$ through $[K_b,P_{b'}]=(b,b')\Delta+(\text{rotation})$ (Remark \ref{rem:c7-KDelta}); the component at $(\eta_1,\eta_2)$, which expresses $K_1$ along $K_{\Gamma(\eta_1,\eta_2)|_V}$ through $S_1$ and the classical transformations along $\eta_1$, $\eta_2$; and the relation $K=P^\dagger$ of Corollary \ref{cor:c8-gram}.""",
    [G + '#13'])

# ---------------------------------------------------------------- Prop c8-C (Theta_I)
add('R20', 'prop:c8-C',
    r"""$\Theta_I$ satisfies (i) of Theorem \ref{thm:c8-fields} with $\dagger$ in place of $\dagger_n$, the scalar kinetic density changing by a total derivative.""",
    r"""$\Theta_I$ satisfies (i) of Theorem \ref{thm:c8-fields} with $\dagger$ in place of $\dagger_n$, the scalar kinetic density changing by a total derivative; (ii) and (iii) there hold for $\Theta_n$ and fail for $\Theta_I$, since $K_t$, $P(\varepsilon,T)$ and the counterterm $\tfrac94\varepsilon^{-1}\frac{\hbar\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{16\pi^2(X,Y)}\int(A,A)$ change under the inversion, and the one-loop radial conjugation is therefore built from $\Theta_n$ (Definition \ref{def:c8-radial}). Two-component: $(\Theta_IA)_m(x)=(x,x)^{-1}\big(\delta_{mk}-2x_mx_k/(x,x)\big)\overline{A_k(I(x))}$, $(\Theta_IA)_{ij}(x)=(x,x)^{-1}\overline{A^{ij}(I(x))}$, $(\Theta_I\lambda)_\alpha{}^i(x)=-i(x,x)^{-2}x_m\sigma^m_{\alpha\dot\alpha}\overline{\bar\lambda_{\dot\alpha i}(I(x))}$, $(\Theta_I\bar\lambda)_{\dot\alpha i}(x)=-i(x,x)^{-2}x_m\overline{\lambda_\alpha{}^i(I(x))}\sigma^m_{\alpha\dot\alpha}$, $(\Theta_I\lambda^*)^\alpha{}_i(x)=-i(x,x)^{-3}x_m\overline{\bar\lambda^{*\dot\alpha i}(I(x))}\bar\sigma^{m\dot\alpha\alpha}$, $(\Theta_I\bar\lambda^*)^{\dot\alpha i}(x)=-i(x,x)^{-3}x_m\bar\sigma^{m\dot\alpha\alpha}\overline{\lambda^{*\alpha}{}_i(I(x))}$, $(\Theta_Ic)(x)=\overline{c(I(x))}$, and $\Theta_IA^*$, $\Theta_Ic^*$ as $\Theta_IA$, $\Theta_Ic$ with the weights $3$, $4$ in place of $1$, $0$ and the factor $-1$.""",
    [C + '#4'])

# ---------------------------------------------------------------- Def c8-radial
add('R21', 'def:c8-radial (assumption)',
    r"""Assume the local representatives along $P_n$ and $K_n$ exist, with $P_1=0$ (Proposition \ref{prop:c7-cov}), and put""",
    r"""Assume the local representatives along $P_n$ and $K_n$ exist, and put""",
    [G + '#0', C + '#7', P + '#4', D + '#14', D + '#31', D + '#3', D + '#8', D + '#24'])

add('R22', 'def:c8-radial (P_1 as a renormalisation condition)',
    r"""modulo $\Lambda^{(0)}$-exact terms and $O(\hbar^2)$. At order $\hbar^0$ this is $O^{\Theta_I}$ (Proposition \ref{prop:c8-C}); at order $\hbar$ it depends on $K_1$ through $\Lambda^{(1)}(K_n)$.""",
    r"""modulo $\Lambda^{(0)}$-exact terms and $O(\hbar^2)$. At order $\hbar^0$ this is $O^{\Theta_I}$ (Proposition \ref{prop:c8-C}); at order $\hbar$ it depends on $P_1$ and $K_1$ through $\Lambda^{(1)}(P_n)$ and $\Lambda^{(1)}(K_n)$. The condition $P_1=0$ is a renormalisation condition. Proposition \ref{prop:c7-cov} gives the covariance of the regulator; the translations also enter cycles through the composite kernel $2\big(\partial_b(P(\varepsilon,\infty)\star)^2\otimes1\big)K_0$, present because $(Q^{GF})^2\neq0$, and on the part of the scalar bilinears with two legs $A|_{W_{\mathbb R}^\perp}$ these graphs give $P_1(O)=-(\tfrac13+\theta_P)$ units times $\{S^{(1)}(P_b),O\}$, with $\theta_P$ the coefficient in the finite counterterm along $P_b$ of Theorem \ref{thm:c9-P1}. On a family stable under $\Lambda$ on which $P_1=a_P\{S^{(1)}(P_b),-\}$ for a number $a_P$, the finite part $\hbar\psi'=\tfrac{a_P}2\{S^{(1)}(\Delta),-\}$ of the operators of the family (Definition \ref{def:c8-renorm}) gives $P_1=0$ there, leaves $D_1$, $M_1$, $R_1$ unchanged and changes the components along $K_b$, $\Sigma$, $\Sigma^*$ by $a_P$, $-\tfrac{a_P}2$, $\tfrac{a_P}2$ times the classical transformation (Proposition \ref{prop:c9-scheme}); it is compatible with the reflection of Theorem \ref{thm:c8-fields} exactly when $a_P$ is real, since $\{S^{(1)}(\Delta),-\}$ acts on a local operator of Weyl weight $w$ as $2(E+w)$ (Lemma \ref{lem:c9-SDelta} (i)), which commutes with $\Theta_n$. A field redefinition generated by a multiple of $S^{(1)}(\Delta)$ (Proposition \ref{prop:c8-laws} (a)) changes these components by the same multiples of the classical transformation on every operator, so the condition $P_1=0$ is imposed family by family through the finite parts of the operators.""",
    [G + '#0', C + '#7', D + '#14', D + '#31', D + '#19'])

# ---------------------------------------------------------------- Prop c8-primary
add('R23', 'prop:c8-primary',
    r"""Assume moreover $M_1(O)=R_1(O)=0$.""",
    r"""Assume moreover $P_1(O)=M_1(O)=R_1(O)=0$, the first in a scheme with $P_1=0$ on the operators considered (Definition \ref{def:c8-radial}).""",
    [G + '#2', P + '#5'])

# ---------------------------------------------------------------- Def c8-2pt
add('R24', 'def:c8-2pt (derivative order)',
    r"""\big\langle O_1,O_2\big\rangle:=\lim_{\varepsilon\to0}\partial_{t_2}\partial_{t_1}W""",
    r"""\big\langle O_1,O_2\big\rangle:=\lim_{\varepsilon\to0}\partial_{t_1}\partial_{t_2}W""",
    [P + '#0'])

add('R25', 'def:c8-2pt (convention and finite counterterm)',
    r"""the coefficient of $t_1t_2$ at the zero field, with $F^{(0)}_{\mathfrak{psl}(4|4)}$ the arity-zero part of the finite counterterm (at one loop $\tfrac{27}{16}\frac{\hbar\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{16\pi^2(X,Y)}\int\mathrm{dvol}\,(d^*A|_V,d^*A|_V)$ plus $s$-closed choices)""",
    r"""the coefficient of $t_2t_1$ at the zero field, i.e. $(-1)^{|O_1||O_2|}$ times the coefficient of $t_1t_2$ ($|O_a|$ the total parity), which at tree level is the free expectation of $O_1O_2$ with $O_1$ on the left, the form in which the Ward identity of Theorem \ref{thm:c8-adjoint} holds, with $F^{(0)}_{\mathfrak{psl}(4|4)}$ the arity-zero part of the finite counterterm (at one loop the arity-zero part of a solution $\hbar F_1$ of Theorem \ref{thm:c8-obstruction} (ii): in the quadratic sector $\tfrac{27}{16}\frac{\hbar\operatorname{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{16\pi^2(X,Y)}\int\mathrm{dvol}\,(d^*A|_V,d^*A|_V)$ plus $Q$-closed quadratic terms, tied to the cubic terms by the cubic part of the obstruction and fixed, modulo the arity-zero parts of $s$-closed cochains, by the arity-one components of $sF_1=\mathrm{Ob}$; through relative order one unit the two-point function of two operators of equal length depends on its quadratic part only)""",
    [P + '#1', C + '#0'])

add('R26', 'def:c8-2pt (two-point cochain at scale T)',
    r"""with the global ghosts the same formula with $I_{\mathfrak{psl}(4|4)}$, $F_{\mathfrak{psl}(4|4)}$, $I^{CT}_{\mathfrak{psl}(4|4)}[\varepsilon]$ defines the two-point cochain.""",
    r"""with the global ghosts the same formula with $I_{\mathfrak{psl}(4|4)}$, $F_{\mathfrak{psl}(4|4)}$, $I^{CT}_{\mathfrak{psl}(4|4)}[\varepsilon]$ and $P(\varepsilon,T)$ in place of $P(\varepsilon,\infty)$ defines the two-point cochain at scale $T$; its components of arity at most one along $\Sigma$, $\Sigma^*$, $\Delta$, the translations and the rotations, and of arity two along $(\xi,\xi')$ and $(\xi,\eta)$, have limits as $T\to\infty$ (Proposition \ref{prop:c8-Tsafe}), while those along $K_b$ and along $(\eta_1,\eta_2)$ with $\Gamma(\eta_1,\eta_2)|_V\neq0$ contain $\log T$ and are used at finite $T$ only.""",
    [C + '#1'])

# ---------------------------------------------------------------- Prop c8-2pt (iii)
add('R27', 'prop:c8-2pt (iii)',
    r"""(the component of Definition \ref{def:c8-components} along $\Delta$ under the identification \eqref{eq:c6-step2})""",
    r"""(on the part of the scalar bilinears with two legs $A|_{W_{\mathbb R}^\perp}$, the component of Definition \ref{def:c8-components} along $\Delta$ when the one-loop quantisation and the local representatives exist through arity two along $(\Delta,P_b)$, Corollary \ref{cor:c9-ident}; the scaling identity is an arity-zero statement and holds with the coefficients of Theorem \ref{thm:c6-dil})""",
    [D + '#13', D + '#7', D + '#23', D + '#29'])

# ---------------------------------------------------------------- Def c8-inner
add('R28', 'def:c8-inner (limit at the origin)',
    r"""On operators at the origin, the radial inner product is the limit as the supports shrink to the origin;""",
    r"""On operators at the origin, the radial inner product is the limit as the supports shrink to the origin, which exists on a family of one classical Weyl weight exactly when the one-loop dilatation matrix of the family is self-adjoint for the tree part; in a scheme with $P_1=0$ on the family it equals the two-point coefficient of Proposition \ref{prop:c8-2pt} (iii) with the first family conjugated, at the separation $\ell$, and in a scheme with $P_1=a_P\{S^{(1)}(P_b),-\}$ on a family of Weyl weight $w$, $a_P$ real, it is $1-2a_Pw$ times that coefficient;""",
    [C + '#2'])

add('R29', 'def:c8-inner (Hermiticity)',
    r"""and $Y^{\rm adj}:=\mathrm{Gram}_{\rm tree}^{-1}\,\overline{Y}^{\,T}\mathrm{Gram}_{\rm tree}$ is the adjoint for the tree inner product.""",
    r"""and $Y^{\rm adj}:=\mathrm{Gram}_{\rm tree}^{-1}\,\overline{Y}^{\,T}\mathrm{Gram}_{\rm tree}$ is the adjoint for the tree inner product. In a reflection-compatible scheme $\langle b,a\rangle_n=\overline{\langle a,b\rangle_n}$ modulo $O(\hbar^2)$, since $W(P(\varepsilon,T),J)^{\Theta_n}=W(P(\varepsilon,T),J^{\Theta_n})$, the value of $F^{\Theta_n}$ at the zero field is the conjugate of that of $F$, $(FG)^{\Theta_n}=G^{\Theta_n}F^{\Theta_n}$ and $\Theta_n^2=1$ (Theorem \ref{thm:c8-fields}); with the Ward identity along $P_n-K_n$, which moves $\exp(-\tfrac\pi4\delta^\hbar)$ to the second argument, the same holds for $\langle-,-\rangle_{\rm rad}$, so $\mathrm{Gram}_{\rm tree}$ and $\mathrm{Gram}_{\rm 1}$ are Hermitian.""",
    [C + '#5'])

# ---------------------------------------------------------------- Thm c8-adjoint
add('R30', 'thm:c8-adjoint',
    r"""the local representatives exist on the operators considered, with $P_1=M_1=R_1=0$; and the two-point function satisfies, for every $Y\in\mathfrak{psl}(4|4)$, the Ward identity $\langle\Lambda^{(1)}(Y)F,G\rangle=-(-1)^{|Y||F|}\langle F,\Lambda^{(1)}(Y)G\rangle$ on the observables considered.""",
    r"""the local representatives exist on the operators considered; and the two-point function of Definition \ref{def:c8-2pt} satisfies, for every $Y\in\mathfrak{psl}(4|4)$, the Ward identity $\langle\Lambda^{(1)}(Y)F,G\rangle=-(-1)^{|Y||F|}\langle F,\Lambda^{(1)}(Y)G\rangle$ on the observables considered, where $|F|$, $|Y|$ are total parities.""",
    [G + '#1', P + '#3', D + '#8', D + '#14'])

# ---------------------------------------------------------------- Cor c8-gram
add('R31', 'cor:c8-gram (one length)',
    r"""Under the hypotheses of Theorem \ref{thm:c8-adjoint}, on a family of operators at the origin closed under the generators considered,""",
    r"""Under the hypotheses of Theorem \ref{thm:c8-adjoint}, on a family of operators at the origin of one length, closed under the generators considered,""",
    [P + '#6'])

add('R32', 'cor:c8-gram (third relation)',
    r"""With $P_1=0$ the third relation is the statement $K=P^\dagger$: $K_1$ along $K_{\bar b}$ equals $\hbar\big[\{S^{(1)}(K_{\bar b}),-\},\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}\big]$.""",
    r"""In a scheme with $P_1=0$ on the family (Definition \ref{def:c8-radial}) the third relation is the statement $K=P^\dagger$: $K_1$ along $K_{\bar b}$ equals $\hbar\big[\{S^{(1)}(K_{\bar b}),-\},\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}\big]$. For $P_1=a_P\{S^{(1)}(P_b),-\}$ on the family, $a_P$ real, it reads $K_1$ along $K_{\bar b}=a_P\{S^{(1)}(K_{\bar b}),-\}+\hbar\big[\{S^{(1)}(K_{\bar b}),-\},\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}\big]$, since $(\{S^{(1)}(P_b),-\})^{\rm adj}=-\{S^{(1)}(K_{\bar b}),-\}$, with the radial Gram matrix of Definition \ref{def:c8-inner} in that scheme, $1-2a_Pw$ times the two-point coefficient on Weyl weight $w$; with the ratio of the two-point coefficients in place of $\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\rm 1}$ the first term is $-a_P\{S^{(1)}(K_{\bar b}),-\}$, since $\{S^{(1)}(K_b),-\}$ lowers the Weyl weight by one.""",
    [C + '#3', G + '#3', D + '#8', D + '#14'])

add('R33', 'cor:c8-gram (outer direction)',
    r"""which is also the condition for $aG_B$ to be reflection compatible, $B$ anticommuting with $\dagger_n$: the conjugation fixes the real part of the outer direction.""",
    r"""which is also the condition for $aG_B$ to be reflection compatible, $B$ anticommuting with $\dagger_n$: reflection compatibility (Theorem \ref{thm:c8-fields}) fixes the real part of the outer direction, since the arity-one finite counterterm along $S_-\otimes W^*$ must be the $\Theta_n$-image of the one along $S_+\otimes W$.""",
    [P + '#8'])

add('R34', 'cor:c8-gram (proof, two lengths)',
    r"""For the outer direction, the change $-aX_0$ requires the change $-\bar aY_0$ on the right side, while Proposition \ref{prop:c8-laws} gives $+aY_0$; hence $\bar a=-a$.""",
    r"""For the outer direction, the change $-aX_0$ requires the change $-\bar aY_0$ on the right side, while Proposition \ref{prop:c8-laws} gives $+aY_0$; hence $\bar a=-a$. On a family containing lengths $L$ and $L+1$ the tree norms of the two lengths differ by one power of the unit; in a basis normalised by the tree norms the relation is expanded in the square root of the unit: $Y_{\frac12}=X_{\frac12}^{\rm adj}+\big[Y_0,\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\frac12}\big]$ and $Y_1=X_1^{\rm adj}+\big[Y_0,\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_1\big]+\big[X_{\frac12}^{\rm adj},\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\frac12}\big]-\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\frac12}\big[Y_0,\mathrm{Gram}_{\rm tree}^{-1}\mathrm{Gram}_{\frac12}\big]$, with $X_{\frac12}$ the classical length-increasing part and $\mathrm{Gram}_{\frac12}$ the entries between the two lengths of order one half.""",
    [P + '#7'])
