# Audit — numbers lens: three decisive coefficients

Script (fresh, exact sympy; no project code imported or copied): `qa/audit_numbers_src/audit_numbers.py`
Runs: `ALG=su2` → `qa/audit_numbers_src/audit_numbers_results_su2.json` (2 s, 81 MB); `ALG=su3` → `qa/audit_numbers_src/audit_numbers_results_su3.json` (5 s, 79 MB). Both give identical coefficients.
Compared against: `compare/compare.json` (`final_tables`, `coefficients`).

Unit for every coefficient below:
$$
\frac{\hbar\,\mathrm{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{(X,Y)\,16\pi^2}.
$$

## Inputs used (DESIGN, classical_spec §4, §7)

- weight $e^{+S/\hbar}$, $S$ = [ESW (31)], $(F_A,F_A)$ = Hodge pairing $\otimes(-,-)$;
- every line $P(\varepsilon,\infty)=\int_\varepsilon^\infty(Q^{GF}\otimes1)K_t\,dt$, $H$ = positive Laplacian on every component (Feynman gauge), $A$–$A$ covariance positive;
- one loop: $-\tfrac\hbar2\log\det$ of the $A$ block, $+\tfrac\hbar2\log\det(\not d_A)$ ($\lambda$ Pfaffian), $+\hbar\log\det$ of the $c$–$A^*$ block;
- $I^{CT}[\varepsilon]$ = $\{\varepsilon^{-1},\log\varepsilon\}$ part of $W(P(\varepsilon,\infty),I)$ at $\hbar^1$; $\partial_\delta I^{CT}[\varepsilon]$ = same for $\partial_\delta W(P(\varepsilon,\infty),I+\delta O)$, 1PI (1PR graphs cancel against the $\delta^0$ counterterm);
- singular parts: $\int_\varepsilon^\infty\!\int_\varepsilon^\infty dt\,dt'\,(t+t')^{-3}=\tfrac1{4\varepsilon}$, tadpole $P(\varepsilon,\infty)$ at coincident points $\propto\tfrac1{16\pi^2\varepsilon}$; $\log\varepsilon$ from the scale-invariant part over the full simplex.

## Results

$$
I^{CT}[\varepsilon]\supset\frac{\hbar\,\mathrm{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{(X,Y)\,16\pi^2}\Big[\tfrac94\,\varepsilon^{-1}\!\int\!\mathrm{dvol}\big((A|_V,A|_V)+(A|_{W^\perp},A|_{W^\perp})\big)-\log\varepsilon\!\int\!\mathrm{dvol}\big((dA|_V,dA|_V)+(dA|_{W^\perp},dA|_{W^\perp})\big)+0\cdot\log\varepsilon\!\int\!\mathrm{dvol}\,(d^*A|_V,d^*A|_V)\Big]
$$

$$
O=\int\mathrm{dvol}\,(\iota_{w_1}A,\iota_{w_2}A):\quad
\partial_\delta I^{CT}[\varepsilon]=\frac{\hbar\,\mathrm{tr}_{\mathfrak g}(\mathrm{ad}_X\mathrm{ad}_Y)}{(X,Y)\,16\pi^2}\log\varepsilon\Big[4\,O-2\,(w_1,w_2)\!\int\!\mathrm{dvol}\sum_w(\iota_wA,\iota_wA)+0\cdot(w_1,w_2)\!\int\!\mathrm{dvol}\,(A|_V,A|_V)\Big]+\frac{2\hbar\dim\mathfrak g\,(w_1,w_2)}{16\pi^2\varepsilon}\int\mathrm{dvol}
$$
($w$ over an orthonormal basis of $W_{\mathbb R}^\perp$, $\dim W_{\mathbb R}^\perp=6$.)

| # | coefficient | audit | compare.json | verdict |
|---|---|---|---|---|
| 1 | $\log\varepsilon$ of $\int\mathrm{dvol}(dA\vert_{W^\perp},dA\vert_{W^\perp})$ in $I^{CT}[\varepsilon]$ | $-1$ | $-1$ (`ICT_quadratic.B4.log_eps`) | agree |
| 1′ | $\varepsilon^{-1}$ of $\int\mathrm{dvol}(A\vert_{W^\perp},A\vert_{W^\perp})$ | $9/4$ | $9/4$ (`B5`) | agree |
| 2a | traceless ($(w_1,w_2)=0$): $\log\varepsilon$ of $\partial_\delta I^{CT}[\varepsilon]$, multiplying $O$ | $4$ | $4$ (`partial_delta_ICT."O log"`) | agree |
| 2b | trace $\sum_w(\iota_wA,\iota_wA)$: $\log\varepsilon$ coefficient multiplying the trace | $4+6\cdot(-2)=-8$ | $4,\,-2$ → $-8$ (`O_eps.trace`) | agree |
| 2c | $(w_1,w_2)\int\mathrm{dvol}(A\vert_V,A\vert_V)$ mixing | $0$ | $0$ | agree |
| 2d | field-independent $\varepsilon^{-1}$, per $\hbar\dim\mathfrak g\,(w_1,w_2)/16\pi^2$ | $2$ (trace: $12$) | $2$ (trace: $12$) | agree |
| 3 | $\varepsilon^{-1}$ of $\int\mathrm{dvol}(A\vert_V,A\vert_V)$ in $I^{CT}[\varepsilon]$ | $9/4$ | $9/4$ (`B3.eps_inv`) | agree |
| 3′ | $\log\varepsilon$ of $\int\mathrm{dvol}(dA\vert_V,dA\vert_V)$; of $\int\mathrm{dvol}(d^*A\vert_V,d^*A\vert_V)$ | $-1$; $0$ | $-1$; $0$ | agree |

Per graph (same unit; `final_tables.ICT_two_point_per_graph`, `partial_delta_ICT_scalar_bilinear.per graph`):

| graph (census) | audit $\varepsilon^{-1}$ | audit $\log\varepsilon$ | compare.json |
|---|---|---|---|
| $A\vert_V$ ext, $A\vert_V$ lines, cubic² (1) | $-9/16$ | $19/24$; $-1/8$ | same |
| $A\vert_V$ ext, $A\vert_{W^\perp}$ lines, cubic² (2) | $-3/4$ | $-1/2$; $0$ | same |
| $A\vert_V$ ext, quartic, $A\vert_V$ line (3) | $3/2$ | $0$ | same |
| $A\vert_V$ ext, quartic, $A\vert_{W^\perp}$ line (4) | $3$ | $0$ | same |
| $A\vert_V$ ext, $c$–$A^*$ loop (5) | $1/16$ | $1/24$; $1/8$ | same |
| $A\vert_V$ ext, $\lambda$ loop (6) | $-1$ | $-4/3$; $0$ | same |
| $A\vert_V$–$A\vert_{W^\perp}$ ext, all loops (7) | $0$ | $0$ | same |
| $A\vert_{W^\perp}$ ext, mixed lines, cubic² (8) | $-1/4$ | $1$ | same |
| $A\vert_{W^\perp}$ ext, quartic, $A\vert_V$ line (9) | $2$ | $0$ | same |
| $A\vert_{W^\perp}$ ext, quartic, $A\vert_{W^\perp}$ line (10) | $5/2$ | $0$ | same |
| $A\vert_{W^\perp}$ ext, $\lambda$ loop (11) | $-2$ | $-2$ | same |
| $O$ vacuum (O1) | $2$ | — | same |
| $O$ + quartic, $A\vert_V$ ext (O3) | — | $A\vert_V$-mix $-2$ | same |
| $O$ + cubic², $A\vert_V$ ext (O4) | — | $A\vert_V$-mix $+2$ | same |
| $O$ + quartic, $A\vert_{W^\perp}$ ext (O5) | — | $O$: $2$; trace-mix: $-2$ | same |
| $O$ + cubic², $A\vert_{W^\perp}$ ext (O6) | — | $O$: $2$; trace-mix: $0$ | same |

## Internal checks of the audit script

| check | result |
|---|---|
| cubic/quartic vertex functions vs direct expansion of $-\tfrac14(F_A,F_A)$ on plane waves | 31/31 |
| $\rho(v)\rho(v')+\rho(v')\rho(v)=2(v,v')$, all ten directions, exact Gaussian integers; trace over $\Sigma$ of four Clifford factors | exact |
| singular-part rule vs high-precision quadrature (remainder after subtracting $\varepsilon^{-1}$, $\log\varepsilon$ terms) | $-0.3072350,\,-0.3072507,\,-0.3072523$ at $\varepsilon=10^{-4},10^{-5},10^{-6}$ |
| $\log\varepsilon$ part transverse: $d^*$-term $-1/8$ (graph 1) $+1/8$ (graph 5) | $0$ |
| $A\vert_V$–$A\vert_{W^\perp}$ mixing; off-colour entries | $0$; $0$ |
| $\mathfrak{su}(2)$ vs $\mathfrak{su}(3)$ | identical |
| two-point coefficients vs normalisation of $(F_A,F_A)$ relative to $(-,-)$ | independent |
| insertion coefficients vs that normalisation | inversely proportional |

Normalisation dependence of row 2: with $(F_A,F_A)$ read as twice the Hodge pairing, rows 2a–2d become $2,\,-4,\,0,\,1$. The Hodge reading is the one compatible with [ESW (33)] at free level: $\delta A=-\Gamma(Q,\lambda)$ in $-\tfrac14 c\,(F_A,F_A)$ gives $c$ times the functional that $\delta\lambda=\tfrac12\rho(F_A)Q$ in $\tfrac12(\lambda,\not d\lambda)$ gives with coefficient $\tfrac12$ (Clifford reduction of $\rho(v)\rho(F_A)$), so $|c|=1$ for the Hodge pairing (hand derivation, not machine-checked here; classical_spec §4 reports the machine check).

Sanity (textbook, not input): doubled $\log\varepsilon$ coefficients per loop reproduce the Feynman-gauge pattern $19/12,\,-11/6$ (gauge), $1/12,\,1/6$ (ghost), $-8/3$ ($\lambda$), $-1$ ($A|_{W^\perp}$ lines), total $-2$; traceless: $4$ = $\log\varepsilon$ ratio of the $A|_{W^\perp}$ kinetic counterterm to its tree term ($-1/(-\tfrac14)$) (no anomalous scaling); trace: $-8\neq4$.

## Verdict

(1) agree, (2) agree, (3) agree — exact equality with compare.json; no disagreement found.
Scope note: the $\varepsilon^{-1}$ values ($9/4$, $2$) are specific to the DESIGN prescription (lower cutoff $\varepsilon$ on every line of every graph); they are confirmed for that prescription only.
