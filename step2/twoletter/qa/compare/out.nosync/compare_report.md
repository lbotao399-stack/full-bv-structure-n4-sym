# R8 two-letter kernels: engine A vs engine B (comparator)

Data: `compare.json` (mirror `out.nosync/compare.json`). Code: `src.nosync/` (`cmp_lib.py`, `cmp_W0.py`, `cmp_deriv.py`, `check_census.py`, `check_predictions.py`, `check_twist2.py`, `check_p6b.py`, `check_closed.py`, `check_shape_numeric.py`, `build_compare.py`). Unit: $u=\hbar\,\mathrm{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)/(16\pi^2(X,Y))=-g^2N/(32\pi^2)$. No dimensional regularisation, $T=\infty$, per-edge cutoff $t_e\ge\varepsilon$.

## 0. Verdict

| gate | content | result |
|---|---|---|
| C1 | $W_0$, 21 pair types, every leg set / ordered word / block structure / kind / component / monomial | 389/389 rows; 76 902 numbers (word totals) + 77 286 (blocks); $\max\lvert A-B\rvert\le10^{-12}$ |
| C1n | negative controls (wrong conventions) | 188, 144, 30 rows fail |
| C2 | engine B $W_0{+}1$ (all pairs), $W_0{+}2$ ($A_V A_V$, $A_V A_\perp$, $A_\perp A_\perp$), $\le2$ outputs, vs engine A on demand | 164/164 |
| C3 | census, output bound, degrees, ghost number | 316 nonzero kernels, 15 641 monomials, 0 violations |
| C4–C6 | P1, P2, P6 from the agreed tables (own assembly) | pass |
| C7–C10, C14 | twist two: two_site.json vs tables; own assembly; $=$ engine A matrices; P6b matrix level | pass |
| C11–C13 | $k=0$ lemma; eq. tl-bubbleclosed; P8 by quadrature at finite $\varepsilon$ | pass |

Engine wrong on any compared entry: none.

## 1. Basis and ordering dictionary (engine B $\to$ engine A)

$$A^{(A)}_M=A^{(B)}_{D2B(M)},\quad D2B=(1,2,3,0)\ (M\le3),\qquad A^{(A)}_{M}=A^{(B)}_{I=M-3}\ (4\le M\le8),\qquad A^{(A)}_9=-A^{(B)}_{I=6}$$
$$\lambda_A=U\,B_S\,\lambda_B,\qquad \lambda^*_A=U^*B_{S^*}\lambda^*_B,\qquad \rho_A(e_M)|_\Sigma=\bar s_M,\ \ \bar s_M U=U^*\rho_D(\pm e_M)$$
$$\zeta^{(B)}_{\rm letter}=(U B_S)^{T}\,\mathrm{BASIS}_{i_a},\qquad \lambda^{(B)}_{\rm output}=(UB_S)^{-1}\mathrm{BASIS}_{i_o}$$
$$\text{word}_A=(Y_1,\dots,Y_k)\ \ (\mathrm{ad}_{Y_k}\cdots\mathrm{ad}_{Y_1})\ \longleftrightarrow\ \text{word}_B=(\ell(Y_k),\dots,\ell(Y_1)),\qquad \mathrm K_A=\operatorname{sgn}(\pi_{\rm odd})\,\mathrm K_B$$

$\pi_{\rm odd}$: permutation from increasing label order (B) to chain order from letter $a$ (A) of the odd outputs.

Two-component (Weinberg–Srednicki, Euclidean realisation D of `step2/notation`): engine B strings are already two-component, $A_0^{(B)}=A_{m=4}$, $A^{(B)}_{k}=A_{m=k}$, $\bar\lambda_{\dot\alpha i}=\epsilon_{\dot\alpha\dot\beta}\bar\lambda^{\dot\beta}{}_i$ ($B_S$), $A_{ij}=A_I\Sigma^I_{ij}$.

## 2. Resolved disagreements

| id | symptom | resolution | engine wrong |
|---|---|---|---|
| R-1 | identity map on $V$: 188/389 rows | $D2B=(1,2,3,0)$ | none |
| R-2 | intertwiner of check_dictionary A1–A3 with $M\to M$: 144/389 rows | engine A: $\rho\vert_\Sigma=\bar s_M$ (r2_core Yukawa $\lambda^T\bar s_M\lambda A_M$), $\bar s_9=-s_9$ $\Rightarrow$ $e_9\to-e_9$ | none (notation layer, F-1) |
| R-3 | no Koszul sign: 30/389 rows ($c\vert c\to cc$: $+2$ vs $-2$) | $\mathrm K_A=\operatorname{sgn}(\pi_{\rm odd})\mathrm K_B$ | none |
| R-4 | odd-pair Hessian sign (B: $+1$ empirical; A: $(-1)^{\lvert X_a\rvert+\lvert X_b\rvert}$ derived) | all odd-pair entries agree; W1 fixes it in both | none |
| R-5 | comparator: 5/164 derivative rows | zero direction $\Rightarrow$ empty polynomial read as $1$ by `ea_kern.kernel`; resampled | none |
| R-6 | comparator: sign of $R$ | $iR=\partial$ on the outputs, $R=p_a+p_b=\sum_j r_j$ | none |

## 3. Census and output bound (C3)

$$\deg_{r}[\log\varepsilon]\mathrm K=W-w_{\rm out},\qquad \deg_r[\varepsilon^{-m}]\mathrm K=W-w_{\rm out}-2m,\qquad n_{\rm out}\le W+d_{ab},\qquad \textstyle\sum_{\rm out}\mathrm{gh}=d_{ab}$$

Nonzero outside the predicted sets: none. Predicted, zero in both engines:
$$\mathrm K(A_\perp,A_\perp)\big|_{A_V}=0\ \ (\text{graph (b)}),\qquad \mathrm K(A_\perp,A_\perp)\big|_{A_VA_V}=\underbrace{(+2)}_{\rm triangle}+\underbrace{(-2)}_{\rm bubble}=0,$$
$$\mathrm K^{(0)}\big|_{(\rm none)}\propto\varepsilon^{-W/2}\ \text{only}\ \Rightarrow\ [\varepsilon^{-1}]=0\ \text{at}\ W=3\ (\lambda_+\lambda_-,\ cA^*_V,\ A A\ \text{with one derivative}),\ W=4.$$

## 4. Predictions against the agreed kernels

P1 (own contraction, $\mathcal O=\int\omega(\iota_{w_1}A,\iota_{w_2}A)$, units $u$):
$$(e)=2\mathcal O-2(w_1,w_2){\textstyle\int}(A_\perp,A_\perp),\quad (f)=2\mathcal O,\quad (c)=-2(w_1,w_2){\textstyle\int}(A_V,A_V),\quad (d)=+2(w_1,w_2){\textstyle\int}(A_V,A_V),\quad (b)=0$$
$$\mathrm{lc}=4\mathcal O\ (\text{traceless}),\qquad \mathrm{lc}=4\mathcal O-12\mathcal O=-8\mathcal O\ (\text{trace})$$
Two-component: $\iota_wA=\tfrac14 w^{ij}A_{ij}$, $(w_1,w_2)=\tfrac14 w_1^{ij}w_{2\,ij}$, $\int(A_\perp,A_\perp)=\tfrac14\int(A^{ij},A_{ij})$, $A^{ij}=\tfrac12\epsilon^{ijkl}A_{kl}$.

P2: $\mathrm{lc}^{\rm natural}=\mathrm{lc}-u\{F,\mathcal O\}$, $\{F,\mathcal O\}=4\mathcal O$:
$$\mathrm{lc}^{\rm natural}(20')=0,\qquad \mathrm{lc}^{\rm natural}(\mathcal K)=-12u\,\mathcal K,\qquad \gamma_{\mathcal K}=-24u=\frac{3g^2N}{4\pi^2}.$$

P6 (own assembly from $\mathrm K(A_\perp,A_\perp)|_{A_\perp A_\perp}$, colour reading $(-1)^r\tfrac{\rm value}2u$):
$$\mathrm{tr}(\cdots A_IA_J\cdots)\mapsto u\big[2\,\mathrm{tr}(\cdots A_JA_I\cdots)-\delta_{IJ}\mathrm{tr}(\cdots A_KA_K\cdots)\big],\qquad \mathrm{lc}^{\rm natural}=-2u\sum_l\big(1-P_{l,l+1}+\tfrac12K_{l,l+1}\big)$$
$$(A_{ij},A_{kl})\mapsto u\big[2A_{kl}A_{ij}-\tfrac12\epsilon_{ijkl}A^{mn}A_{mn}\big]$$
$$\mathrm{Tr}(Z^{L-M}X^M),\ L\le8:\ \ \det\!\big(E-\mathrm{lc}^{\rm natural}/(-2u)\big)=\text{P6 charpolys (23/23)},\qquad Z=-\tfrac14(\bar\Sigma^1+i\bar\Sigma^2)^{ij}A_{ij},\ X=-\tfrac14(\bar\Sigma^3+i\bar\Sigma^4)^{ij}A_{ij}$$

P5 twist two, $\mathcal O_k=\mathrm{tr}(D^kZ\,D^{N-k}Z)$, $D=n\!\cdot\!\partial+[n\!\cdot\!A,-]$, $n=(1,i,0,0)$; own assembly of the $ZZ$ part:
$Z_a:=(n\!\cdot\!\partial)^aZ$, $A_i:=(n\!\cdot\!\partial)^i(n\!\cdot\!A)$; $h$, $\mathrm{AZ}$, $\mathrm{ZA}$ from two_site.json (units $u$; entries with $a+b\le2$ re-derived from the agreed tables, C7):
$$D^{m'}Z\ni\sum_{0\le i\le j<m'}\tbinom{j}{i}\,[A_i,Z_{m'-1-i}],\qquad G(m')=\{(i,m'-1-i,\tbinom ji)\}$$
$$\begin{aligned}\mathrm{lc}^{\rm natural}(\mathcal O_k)\big|_{ZZ}=&\sum_{(x,y)\in\{(k,N-k),(N-k,k)\}}\sum_{c,d}h_{(c,d)\leftarrow(x,y)}\,\mathrm{tr}(Z_cZ_d)-4\,\mathrm{tr}(Z_kZ_{N-k})\\&+\sum_{(i,m,c)\in G(k)}c\Big[(\mathrm{AZ}_{i,m}-\mathrm{ZA}_{i,m})\,\mathrm{tr}(Z_{i+m+1}Z_{N-k})+(\mathrm{ZA}_{i,N-k}-\mathrm{AZ}_{i,N-k})\,\mathrm{tr}(Z_mZ_{i+N-k+1})\Big]\\&+\sum_{(i,m,c)\in G(N-k)}c\Big[(\mathrm{AZ}_{i,m}-\mathrm{ZA}_{i,m})\,\mathrm{tr}(Z_kZ_{i+m+1})+(\mathrm{ZA}_{i,k}-\mathrm{AZ}_{i,k})\,\mathrm{tr}(Z_{i+k+1}Z_m)\Big]\end{aligned}$$
$$\operatorname{eig}\mathrm{lc}^{\rm natural}=-8h(S)\,u,\ S\ \text{even}\le N\le8:\ 0,-12,-\tfrac{50}3,-\tfrac{98}5,-\tfrac{761}{35};\qquad \gamma=2\operatorname{eig}=\frac{g^2N}{2\pi^2}h(S)\ \ [\text{B03c (2.31)}]$$
$$ZZ\ \text{part without the}\ O(A_n)\ \text{terms}:\ N=2:\{-\tfrac83,0\},\ N=4:\{-\tfrac{16}5,-\tfrac83,0\}\ \ (=2\cdot\tfrac{2}{j+1}-4)$$
Own matrices $=$ engine A matrices (kernels.json), $N=0,\dots,8$, entrywise. Two-component: $n_{\alpha\dot\alpha}=n_m\sigma^m_{E\,\alpha\dot\alpha}=\kappa_\alpha\tilde\kappa_{\dot\alpha}$, $n\!\cdot\!D=-\tfrac12\kappa_\alpha\tilde\kappa_{\dot\alpha}D^{\dot\alpha\alpha}$.

P6b (matrix level, new): in the letter basis $D^kZ/k!\leftrightarrow(a^\dagger)^k|0\rangle$,
$$\mathrm{lc}^{\rm natural}\big|_{\rm twist\,two}=-2\big(H'_{12}+H'_{21}\big)u,\qquad H'_{12}(a_1^\dagger)^k(a_2^\dagger)^{n-k}|00\rangle:\ \big(h(k)+h(n-k)\big)\delta_{kk'}-\frac{1-\delta_{kk'}}{\lvert k-k'\rvert}\quad(N\le8).$$

P3/P7: $D_1=4\mathrm{lc}^{\rm natural}-\{F'_\Delta,-\}+C_\Delta$; $4\mathrm{lc}^{\rm natural}=0,\,-48u$ on the bilinears (both); composite parts are Ch9 inputs in both engines; P7 untested.

P4: W1 passes in both; W2 untested.

## 5. Closed forms

$$16\pi^2S_\varepsilon(v)=\sum_{m\ge1}\varepsilon^{-m}\frac{(-1)^{m+1}}{m!}\Big(\frac{\lvert v\rvert^2}4\Big)^{m-1},\qquad \sum\frac{s^\alpha u^\beta}{\alpha!\beta!}\mathrm K^{(0)}=2(X_a\!\cdot\!X_b)\,16\pi^2S_\varepsilon(u-s),\qquad \lvert v\rvert^2=-\tfrac12v_{\alpha\dot\alpha}v^{\dot\alpha\alpha}$$
(C11: engine B (none) sets, $W_0,W_0{+}1,W_0{+}2$.)

Quartic bubble (own derivation, total-time cutoff; $\log\varepsilon$ coefficient shape independent):
$$16\pi^2[\log\varepsilon]=-N_V\sum_j\frac{(\lvert v\rvert^2\lvert R\rvert^2/4)^j}{(j!)^2}\int_0^1d\alpha\,\alpha^j(1-\alpha)^je^{i(1-\alpha)(u,R)+i\alpha(s,R)}$$
$$=-N_V\sum_{a,b,j}\big(i(s,R)\big)^a\big(i(u,R)\big)^b\lvert u-s\rvert^{2j}\Big(\frac{\lvert R\rvert^2}4\Big)^j\frac{(a+j)!\,(b+j)!}{a!\,b!\,(j!)^2\,(a+b+2j+1)!},\qquad iR=\partial\big|_{\rm outputs}$$
$$(s,R)=-\tfrac12s_{\alpha\dot\alpha}R^{\dot\alpha\alpha},\qquad \lvert R\rvert^2=-\tfrac12R_{\alpha\dot\alpha}R^{\dot\alpha\alpha}$$
(C12: 1350 checks against engine A quartic blocks, $n\le3$.)

P8 (C13), $E=2$, $N=\lvert l\rvert^4$, $\lvert p\rvert^2=1$:
$$16\pi^2I(\varepsilon)=\iint_{t_i\ge\varepsilon}\frac{dt_1dt_2}{T^2}\Big[\frac6{T^2}+\frac{6c^2\lvert p\rvert^2}{T}+c^4\lvert p\rvert^4\Big]e^{-t_1t_2\lvert p\rvert^2/T},\quad c=\tfrac{t_2}T$$
$$\varepsilon\partial_\varepsilon I+\frac{2}{4\varepsilon^2}+\frac{\lvert p\rvert^2}{8\varepsilon}=8.3\cdot10^{-5},\,8.3\cdot10^{-6},\,8.3\cdot10^{-7},\,8.3\cdot10^{-8}\quad(\varepsilon=10^{-3..-6})\ \Rightarrow\ [\log\varepsilon]=0.$$

## 6. Findings

- F-1: translating engine-A fermionic components to two-component form with $U$ of `check_dictionary.py` needs $A^{(A)}_9=-A_{I=6}$ (A1 "vector map $M\to M$" holds for $\rho|_\Sigma=s_M$, engine A uses $\bar s_M$).
- F-2: engine B, identical outputs: functional $=\big(\prod_t n_t!\big)^{-1}\times$ labelled multilinear form (B states $1/2$; used only for $n=2$).
- F-3: two_site.json `D1_eigenvalues_units` strings omit $u$.
- F-4: P6b holds at matrix level for twist two, $N\le8$ (Sec. 4).
- F-5: expected zeros (Sec. 3).

## 7. Open

P7 / T3 composite kernels; W2; leg sets with $\ge3$ outputs at $W>W_0$ and $W_0{+}3\dots6$ (no engine-B tables); fermion and field-strength spectra; P6b beyond two sites.
