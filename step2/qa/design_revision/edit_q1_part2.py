p='/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/literature/q1_targets.tex'
s=open(p).read()
R=[]
def r(a,b): R.append((a,b))
r(r"""(c) Project counterpart: $O(\varepsilon)$ has a unique singular part and a free finite part (DEFINITIONS\_v2 item 5); finite counterterms in $C^\bullet(\mathfrak{psl}(4|4),\mathscr O_l(\mathscr E))$ (ch07, Step 2). Section~\ref{sec:freedom}.
\end{target}""",
  r"""(c) Project counterpart: $O(\varepsilon)$ has a unique singular part and a free finite part (DEFINITIONS\_v2 item 5); finite counterterms in $C^\bullet(\mathfrak{psl}(4|4),\mathscr O_l(\mathscr E))$ (ch07, Step 2). Section~\ref{sec:freedom}.
\end{target}

\begin{target}[the rank-$(1,1)$ supercharge: no loop corrections]\label{t:ht}
(a) \Sl\ [CLW] (3.1): the holomorphic-topological supercharge $Q_{HT}=Q^4_-+Q_{\dot-3}$. [CLW] \S3.1: ``all the loop diagrams naively contributing to the BRST differential vanish identically'' (the $(H+T)$-Laman graphs with $H=1$, $T=2$, via the vanishing of the top forms $\omega^{H,T}_\Gamma$ on $3$-Laman graphs), so $Q_{\rm int}$ truncates at tree level to $Q_0$ and obeys the Leibniz rule, to all orders. [CLW] footnote 5: the leading loop correction $Q_1$ computed in [32]$=$[Z05] ``takes the form $[Q_0,-]$''. [CLW] (3.17)--(3.19), [BK] footnote 6, (1.11): in the 2d chiral algebra the one-loop correction $Q^{2d}_1$ is nonzero, $Q_1^{2d}(C^A(\theta)C^B(\theta'))=-2f^{ABC}(\theta^1-\theta'^1)(\theta^2-\theta'^2)\partial C^C$ ((3.19)).

(c) \Dl\ In the project $Q_{HT}$ is $x=\xi+\bar\xi'$ with $\xi=u\otimes w$ and $\bar\xi'\in S_-\otimes W^*$ of rank one with $\bar\xi'_{\dot\alpha i}w^i=0$ (the commutant of Target~\ref{t:psl33} (c)); $\Gamma(x,x)=0$ in all ten components, and the evaluation along $x$ and the twisted action of the design (lem:o-ev) apply. Target: $[Q_1^x]=0$ on $H(Q_0^x)$, from the components along $\xi$, $\bar\xi'$ and $(\xi,\bar\xi')$.

(d) Scheme independent: the class statement $[Q_1^x]=0$ (Lemma~\ref{lem:inv}); the pair values vanish only in the HT scheme of [CLW]. This addresses open item 6 (the [Z05] tension) and is an R3 reproduction check.
\end{target}""")
r(r"their $s$-exact part is a field-type rescaling $N_\beta$, acting as $Q_1\mapsto Q_1+[Q_0,N_\beta]$; the coupling direction (total action) acts as $Q_1\mapsto Q_1+aQ_0$, since $Q_0$ does not depend on the coupling when the coupling sits in $(-,-)$.",
  r"their $s$-exact part is a field-type rescaling $N_\beta$, acting as $Q_1\mapsto Q_1+[Q_0,N_\beta]$; the coupling direction (total action, $G=S-S^{(2)}$, design prop:o-laws (d), $F\mapsto F+\hbar aG$) acts as $Q_1\mapsto Q_1-aQ_0-a[N_{R_\xi},Q_0]$: the $aQ_0$ part is the overall rescaling of the twisted action, the second part the R-rotation of F1; untwisted, it changes the one-loop BRST transformation by $-a(Q+\{I,-\})$, which vanishes on gauge-invariant operators.")
r(r"by Schur's lemma for the covariant subalgebra the $s$-closed non-exact choice is $a\,S^{(1)}(\xi)$, a renormalisation of the supercharge, acting as $Q_1\mapsto Q_1+aQ_0$; $s$-exact choices act as $[Q_0,-]$.",
  r"by Schur's lemma for the covariant subalgebra the closed choice along $\xi$ is $a\,S^{(1)}(\xi)$, a renormalisation of the supercharge. With the correlated F2 rescaling, $S+(1+a\hbar)S^{(1)}(\xi)-\frac12(1+2a\hbar)S^{(2)}(\xi,\xi)$ is $S^\xi$ composed with the R-rotation $\alpha(1+a\hbar)^{-1}$ of [ESW, Def 2.18] at first order ($S$ and $\omega$ are $\alpha$-invariant, $S^{(k)}(\xi,\dots,\xi)$ has $\alpha$-weight $k$), so F1 acts by conjugation: $Q_1\mapsto Q_1+a[N_\alpha,Q_0]=Q_1+aQ_0-a[N_{\rm gh},Q_0]$, a $[Q_0,\varphi]$-type shift, zero on $c$ and $\gamma$ and the term $\frac12\varepsilon_{IJK}[\gamma^J,\gamma^K]$ on $\beta_I$ up to the sign of the letter map (design rem:o-alpha). In the untwisted equivariant setting $a\,S^{(1)}(\xi)$ is the restriction to $\xi$ both of the outer cocycle $G_B$ (F5) and of the exact cocycle of the inner R-rotation $R_\xi$; $s$-exact choices act as $[Q_0,-]$.")
r(r"\item[F4.] Choice of representatives and homotopy data in the transport ([BGKWWY] (D.14) $h_0$; [ESW] Thm.~10.9 eliminations; [BGKWY22] (B.6) propagator changes): $Q_1\mapsto Q_1+[Q_0,\varphi_h]$.",
  r"""\item[F4.] Choice of representatives and homotopy data in the transport ([BGKWWY] (D.14) $h_0$; [ESW] Thm.~10.9 eliminations; [BGKWY22] (B.6) propagator changes): $Q_1\mapsto Q_1+[Q_0,\varphi_h]$.
\item[F5.] \Dl\ The outer derivation $B$ of $\mathfrak{psl}(4|4)$ (induced by $\mathrm{ad}\,\mathrm{diag}(1_4,0_4)$ in $\mathfrak{gl}(4|4)$; $+1$ on $S_+\otimes W$ and $S_-\otimes W$, $-1$ on $S_-\otimes W^*$ and $S_+\otimes W^*$, $0$ on the even part): $G_B=\frac d{dt}(e^{tB})^*S_{\mathfrak{psl}(4|4)}$ is $s$-closed and not $s$-exact (design prop:o-GB). Untwisted: $Q_1\mapsto Q_1-a\,\delta_\xi$ on $S_+\otimes W$, $+a\,\delta_{\bar\xi}$ on $S_-\otimes W^*$; $S_1\mapsto S_1+a\,\delta_\eta$ on $S_+\otimes W^*$, $-a\,\delta_{\bar\eta}$ on $S_-\otimes W$; $D_1,K_1,P_1,M_1,R_1$ unchanged. Twisted: $Q_1\mapsto Q_1+a[N_{R_\xi},Q_0]$, of the form of Lemma~\ref{lem:inv}.""")
r(r"Rescalings of the letters at tree level conjugate $[Q_1]$ by an invertible map of $H(Q_0)$ and leave kernel, image and rank unchanged.",
  r"Rescalings of the letters at tree level conjugate $[Q_1]$ by an invertible map of $H(Q_0)$ and leave kernel, image and rank unchanged. The shifts of F0, F1 and F5, $-aQ_0-a[N_{R_\xi},Q_0]$, $a[N_\alpha,Q_0]$ and $a[N_{R_\xi},Q_0]$, are of this form, since $[N,Q_0]=-[Q_0,N]$ for an even derivation $N$.")
r(r"the coefficient of the Konishi entry before the orthogonality normalisation (Target~\ref{t:konishi});",
  r"the coefficient $c_K$ of the Konishi entry, in every normalisation of the weight-$3$ operators, orthogonality to the protected operator included (Target~\ref{t:konishi});")
r(r"(Target~\ref{t:cgamma}, a selection-rule consequence); the Konishi coefficient in the orthogonality normalisation.",
  r"(Target~\ref{t:cgamma}, a selection-rule consequence); $[Q_1^x]=0$ for the rank-$(1,1)$ supercharge (Target~\ref{t:ht}).")
for a,b in R:
    n=s.count(a); assert n==1,(n,a[:80]); s=s.replace(a,b)
open(p,'w').write(s); print('ok',len(R))
