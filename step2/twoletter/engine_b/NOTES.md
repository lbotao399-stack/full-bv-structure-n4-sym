# Engine B (route B): method and conventions

Folder: `step2/twoletter/engine_b/`. Independent of `step2/twoletter/engine_a/` (never read). Reused, from outside
the workflow root: `engine_b/vertices.py` (action terms on plane waves, Hessian blocks with colour keys),
`engine_b/loops.py` (propagator numerators by inversion of the free Hessian), `engine_b/realization.py` (component
realisation). No dimensional regularisation; T = infinity at generic momenta; no Wolfram, no cluster.

## Definition computed

For a two-letter vertex M = int omega (zeta_a . phi_a)(zeta_b . phi_b) (slot order a, b):

    ct(M) = (hbar/2) sum_{ordered vertex sequences} STr[ P Pi h_1 P Pi h_2 ... P Pi H_M ]   (singular part)

h_v the Hessian of I at the outputs (cubic: one output, quartic: two), P(q) = N(q) e^{-eps q^2}/q^2. This is the
delta-linear part of the one-loop formula W_1 = (hbar/2) sum (1/n) STr[(P Pi H)^n] (the convention with which
engine_b reproduced Ch5), i.e. hbar Sing_eps d_{P_e} M (Thm thm:tl-decomp). K(X_a,X_b) := ct(int omega X_a X_b),
sigma_ab = 1. For odd M the sign of H_M relative to `hessian_blocks` is +1: the naive left-derivative argument gave -1,
which reverses all three W1 statements for (A_perp,c), (A_V,c) (gates.json G3); with +1 they hold, together with the
even-M W1 gates (lambda,c), (c,c), which do not involve the choice.

## Route B steps

1. Numerator algebra in momentum space, sparse exact polynomials over Q(i) (sympy `ring`, QQ_I): vector propagation
   through P, Pi h_1, ..., P and contraction with Pi H_M; edge j carries q - cur_j.
2. Colour: each vertex block carries a key; cubic -> ad_E, quartic K12 -> ad_E1 ad_E2 - ad_E2 ad_E1, K13 -> -ad_E2 ad_E1,
   K14 -> -ad_E1 ad_E2; the chain gives (1 (x) W_k...W_1)(-,-)^{-1}; reported as sum_w L_w C_w,
   C_w = (1 (x) ad_{Y_w1}...ad_{Y_wr})(-,-)^{-1}, X_a on the first factor.
3. Loop integral: for each loop monomial q^mu: Laplacian powers, q -> c(alpha) = sum alpha_j cur_j, Phi(alpha);
   num_j/(j! 4^j), F_m. Pinned edges: eps d/deps w = -sum_e (graph with t_e = eps); with t_f = eps u_f,
   [eps^-m] = (1/m) sum_e I_e(F_m), [log eps] = -sum_e I_e(F_0),
   I_e(alpha^n) = int_{[1,oo)^{E-1}} prod u_f^{n_f} tau_u^{-(m+|n|+E)} du.
4. Feynman parametrisation of the remaining edges (u_f = U beta_f): exact U-integral (polynomial in 1/(1+min beta)),
   then exact one-dimensional beta-integrals region by region; log(3/2) (E=3), log(4/3) (E=4) appear and cancel in every
   pinned integral (asserted, gate G5). E <= 4 implemented (k <= 3 vertices).
5. Numerical check (mpmath, 30 digits): eps dw/deps by quadrature over the pinned faces, without expansion, at random
   rational sample data; residual r(eps) = eps w' + sum m a_m eps^-m - b at eps = 1e-2, 1e-3, 1e-4.

## Units and notation of the output strings

Coefficients in units hbar/(16 pi^2). Xa[comp], Xb[comp]: letter components; Y<j>[comp]: output polarisations;
p<j>_m: output momenta (i p_m = d_m on that output). Odd outputs: Grassmann parameters in front in increasing label
order. Components (realization.py, Euclidean two-component): A_m (m=0..3) = A|_V; A_I1..A_I6 = A|_{W_R^perp};
lam_a<alpha>^<i> = lambda_+ (S_+ (x) W); lambar_ad<alphadot>_<i> = lambda_- (S_- (x) W^*); lams^a<alpha>_<i> =
lambda^*_+; lambars^ad<alphadot><i> = lambda^*_-; As = A^*; cs = c^*; sigma = (1, i sigma_k), sigmabar = (1, -i sigma_k).
For su(N): a planar adjacent pair with coefficient c and word w gives (-1)^r c/2 units u = -g^2 N/(32 pi^2), with the
pair replaced by Y_wr ... Y_w1 (Lemma lem:tl-colour (iii)); the colour contractions used are verified in
colour_rules.py (G4).

## D_1 and P_1

D_1 = 4 lc^natural - {F'_Delta,-} + C_Delta (Thm thm:tl-D1). Engine B computes the two-body log kernels, hence
4 lc^natural = 4 (sum over pairs K^log - u sum_letters n_F); the composite parts (C2_Delta, Pi_Delta, C_{P_b}) are not
computed by route B: on the bilinears they are taken from Ch9 (composite triangle -2u{S1(Delta),O}, one-body
u[-8X + (5/3-theta) v.X] on A_perp), giving D_1 = 4 lc^natural - (1/3+theta) u {S1(Delta),O}; on the spectrum families
the reported D_1 eigenvalues are 4 eig(lc^natural) - (1/3+theta) u (classical dilatation), the second term being the
theory prediction P7, not an engine output. P_1 two-body part: -(1/3+theta_P) u {S1(P_b),-} on the bilinears (Ch9 input).

## Files

tl_integrals.py, tl_engine.py, kernels_W0.py (kernels_W0.json), gates.py (gates.json, two_site.json),
colour_rules.py, export.py (kernels.json), run_all.py, PROGRESS.md, logs log_*.txt.


## Revision 2026-10-29 (qa/disposition.json)

Sign of the odd-pair kernels (derivation, no number changes). The kernel is fixed by twoletter_theory Lemma
lem:tl-sign: K(X_a,X_b) = hbar Sing <X_a(x) X_b(y)>_e at y = x, with <phi_a phi_b>_e = P_e^{phi_a phi_b} read with the
first factor in the slot of X_a; in a supertrace with the Hessian of M instead of that of delta M the factor
d_j d_l (delta M) = (-1)^{|delta|(|j|+|l|)} delta d_j d_l M, |delta| = |M| = |X_a| + |X_b|, enters. Engine A realises
this as (-1)^{|X_a|+|X_b|} times its operator-Hessian supertrace (derived, ea_kern.kernel). Route B's Hessian
convention (vertices.hessian_blocks with the parity factor on the row index, M_template) with ODD_M_SIGN = +1 gives the
same K on every odd-pair entry (qa/compare R-4, rerun results/compare_rerun: 389/389 W_0 rows); so ODD_M_SIGN = +1 is
the value of Lemma lem:tl-sign in route B's convention, and the W1 gates (A_perp,c), (A_V,c) are consistency checks of
it rather than the source of the sign.

Functionals from labelled multilinear forms: a labelled output form with n_t identical outputs of type t becomes a
functional after division by prod_t n_t! (not 1/2 in general); the spectra used n = 2 only.

D_1 entries of two_site.json: rewritten as numbers times the unit u, conditional on prediction P7 at theta = -1/3; the
classical part -(1/3+theta) u {S^(1)(Delta),-} is kept separately; gamma = 2 eig(lc^natural) recorded.

W2: not run independently in route B. Engine B equals engine A entrywise on all compared tables, so engine A's W2 gate
(engine_a/out.nosync/w2_gate.json) holds on those entries; recorded in gates.json as a derived gate.
