"""Audit 2b (Part D theory, D2): cross-engine test of thm:d-D2 (ii) in every quadratic sector.

Claim tested: Ob(Delta)|_quad = 4 x (log-eps coefficient of I^CT)|_quad + Q Pi_Psi, with Pi_Psi the renormalised
class (i) (composite cycle) of the arity-one two-point functional.  Engine A computed class (i) ("Pi_int") in all four
sectors (A_perp, A_V, c, lambda) including its nonlocal parts; engine B computed Ob(Delta) directly from special-edge
graphs.  Here Q is applied (exactly, symbolically) to engine A's class (i) and compared with Ob_B - 4 C_log.

Conventions (checked against engine B's stated identities, see checks at the end):
  two-leg functional F = int dvol (leg2, [alpha(-i d) + x^rho i beta_rho(-i d)] leg1), leg1 = field at momentum p fixed
  (this is engine A's "alpha delta + beta^mu d_mu delta, field leg first, its momentum p fixed").
  Q on functionals = D_Q (directional derivative along Qe, odd parameter moved to the left), with the arrows of
  Prop c1-Q: c -> -dc (A_V slot), A -> -(1/2) d*dA (A* slot), A* -> d*(A*|_V) (c* slot), lambda -> dslash lambda.
  L := log(eps |p|^2) + gamma + log 3 (engine A's L3e), a function of s = |p|^2 with dL/ds = 1/s.
Writes q_of_engineA_classi.json next to this file.
"""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
p = sp.Matrix(sp.symbols('p0:4', real=True))
s = (p.T * p)[0]
Lc = sp.Symbol('Lconst')
L = sp.log(s) + Lc          # L3e with log eps + gamma + log 3 absorbed in Lconst (only dL/ds = 1/s matters)
I4 = sp.eye(4)
e = [I4[:, r] for r in range(4)]
res = {}

def d_dp(expr, r):
    return sp.diff(expr, p[r])

def neg(expr):
    return expr.subs({p[i]: -p[i] for i in range(4)}, simultaneous=True)

# ------------------------------------------------------------------ (A_V, A*_V) and (c, c*) class (i), engine A
Lo = (p * p.T) / s
alpha_V = (sp.Rational(13, 6) + L) * I4 + 2 * Lo
beta_V = [sp.Rational(10, 3) * p[r] * I4 + (L - sp.Rational(1, 2)) * (p * e[r].T) - 2 * p[r] * Lo for r in range(4)]
alpha_c = sp.Rational(-8, 3)
beta_c = [p[r] * (sp.Rational(-5, 6) - L) for r in range(4)]

# (c, A*_V) output.  From Pi_V: A_V leg <- dc (symbol i p): alpha'' = alpha_V (i p), beta'' = beta_V (i p).
# From Pi_c: c* leg <- d*A*: int(d*A*, O c) = int(A*, d(O c)): alpha' = i p alpha_c + i beta_c, beta' = i p beta_c.
alpha_cA = alpha_V * (sp.I * p) + sp.I * p * alpha_c + sp.I * sp.Matrix(beta_c)
beta_cA = [beta_V[r] * (sp.I * p) + sp.I * p * beta_c[r] for r in range(4)]
alpha_cA = sp.simplify(alpha_cA)
beta_cA = [sp.simplify(b) for b in beta_cA]
target_cA = sp.Rational(2, 3) * sp.I * p          # -(2/3) int(dc,A*) = +(2/3) int(A*, dc)
res['(c,A*_V) output: alpha - target'] = str(sp.simplify(alpha_cA - target_cA).T)
res['(c,A*_V) output: beta'] = [str(b.T) for b in beta_cA]

# (A_V, A_V) output.  A*_V leg <- -(1/2) d*dA: K = ell o O_V, ell(k) = -(1/2)(|k|^2 - k k^T) (self-adjoint).
ell = -sp.Rational(1, 2) * (s * I4 - p * p.T)
dell = [sp.Matrix(4, 4, lambda i, j: d_dp(ell[i, j], r)) for r in range(4)]
alpha_K = ell * alpha_V + sum((dell[r] * beta_V[r] for r in range(4)), sp.zeros(4, 4))
beta_K = [ell * beta_V[r] for r in range(4)]

def symmetrise_even(alpha_K, beta_K):
    """canonical form (p1 = p fixed) of the quadratic functional int(A, K A), both orderings of two plane waves:
       alpha_tot = alpha_K(p) + alpha_K(-p)^T - (d_rho beta_K^rho)(-p)^T ; beta_tot^rho = beta_K^rho(p) + beta_K^rho(-p)^T
       (matrices act a1 -> paired with a2)."""
    dbeta = sum((sp.Matrix(4, 4, lambda i, j: d_dp(beta_K[r][i, j], r)) for r in range(4)), sp.zeros(4, 4))
    a_tot = alpha_K + neg(alpha_K).T - neg(dbeta).T
    b_tot = [beta_K[r] + neg(beta_K[r]).T for r in range(4)]
    return sp.simplify(a_tot), [sp.simplify(b) for b in b_tot]

aVV, bVV = symmetrise_even(alpha_K, beta_K)
target_VV = 2 * sp.Rational(1, 3) * (s * I4 - p * p.T)      # (1/3) int(dA_V, dA_V)
res['(A_V,A_V) output: alpha - target'] = str(sp.simplify(aVV - target_VV))
res['(A_V,A_V) output: beta'] = [str(sp.simplify(b)) for b in bVV]

# ------------------------------------------------------------------ (A_perp, A*_perp): alpha = 8/3, beta = 10/3 p
alpha_P = sp.Rational(8, 3)
beta_P = [sp.Rational(10, 3) * p[r] for r in range(4)]
ellP = -sp.Rational(1, 2) * s
aK = ellP * alpha_P + sum(d_dp(ellP, r) * beta_P[r] for r in range(4))
bK = [ellP * beta_P[r] for r in range(4)]
dbK = sum(d_dp(bK[r], r) for r in range(4))
aP_tot = sp.simplify(aK + neg(aK) - neg(dbK))
bP_tot = [sp.simplify(bK[r] + neg(bK[r])) for r in range(4)]
res['(A_perp,A_perp) output: alpha / (2 |p|^2)'] = str(sp.simplify(aP_tot / (2 * s)))
res['(A_perp,A_perp) output: beta'] = [str(b) for b in bP_tot]

# ------------------------------------------------------------------ (lambda, lambda*): engine A structures
# Euclidean Clifford algebra rho(e_mu) = gamma_mu (4x4); fibre transpose = reversal, realised by M -> C^-1 M^T C
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -sp.I], [sp.I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); s0 = sp.eye(2)
sig = [s1, s2, s3, sp.I * s0]
sigb = [-s1, -s2, -s3, sp.I * s0]
def blk(a, b):
    M = sp.zeros(4, 4); M[0:2, 2:4] = a; M[2:4, 0:2] = b; return M
gam = [blk(sig[m], sigb[m]) for m in range(4)]          # {g_m, g_n} = -2 delta (Euclidean sigma of ch08)
gam = [sp.I * g for g in gam]                              # now {g_m, g_n} = +2 delta
assert all(sp.simplify(gam[m] * gam[n] + gam[n] * gam[m] - 2 * (1 if m == n else 0) * sp.eye(4)) == sp.zeros(4, 4)
           for m in range(4) for n in range(4))
Cs = sp.Matrix(4, 4, sp.symbols('C0:16'))
eqs = []
for m in range(4):
    eqs += list(gam[m].T * Cs - Cs * gam[m])
sol = sp.solve(eqs, list(Cs), dict=True)[0]
C = Cs.subs(sol)
free = [x for x in C.free_symbols]
C = C.subs({free[0]: 1}).subs({x: 0 for x in free[1:]})
assert C.det() != 0
Ci = C.inv()
def rev(M):          # fibre transpose (reversal anti-automorphism, rho symmetric)
    return sp.simplify(Ci * M.T * C)
assert all(sp.simplify(rev(gam[m]) - gam[m]) == sp.zeros(4, 4) for m in range(4))
rp = sum((p[m] * gam[m] for m in range(4)), sp.zeros(4, 4))
# engine A: beta^mu = p^mu(-1/3 - 2L)(u,v) + (u, rho(e^mu) rho(p) v)(L - 1/2), u = field lambda, v = antifield;
# (u, rho(e)rho(p) v) = (rho(p) rho(e) u, v): operator on the field
alpha_l = sp.Rational(-8, 3) * sp.eye(4)
beta_l = [p[r] * (sp.Rational(-1, 3) - 2 * L) * sp.eye(4) + (L - sp.Rational(1, 2)) * rp * gam[r] for r in range(4)]
# Q on the lambda* leg: D_Q F = -int(O lambda, dslash lambda) = -int(lambda, K lambda), K = dslash o O (symbol i rho(k))
iK = sp.I * rp
alpha_Kl = iK * alpha_l + sum((sp.I * gam[r] * beta_l[r] for r in range(4)), sp.zeros(4, 4))
beta_Kl = [iK * beta_l[r] for r in range(4)]
# contributing part for odd fields with rho symmetric: (1/2)(K - K^T); K^T of alpha(-i d) + x^r i beta_r(-i d) has
# canonical form alpha(-p)^T + (d_r beta_r)(-p)^T,  beta_r(-p)^T  (^T = rev); for odd legs the second ordering enters with -1
dbl = sum((sp.Matrix(4, 4, lambda i, j: d_dp(beta_Kl[r][i, j], r)) for r in range(4)), sp.zeros(4, 4))
a_anti = sp.simplify((alpha_Kl - (rev(neg(alpha_Kl)) - rev(neg(dbl)))) / 2)
b_anti = [sp.simplify((beta_Kl[r] - rev(neg(beta_Kl[r]))) / 2) for r in range(4)]
# functional = -int(lambda, K_anti lambda); target -(1/3) int(lambda, dslash lambda) <-> K_anti = (1/3) i rho(p)
res['(lambda,lambda) output: K_anti - (1/3) i rho(p)'] = str(sp.simplify(a_anti - sp.Rational(1, 3) * iK))
res['(lambda,lambda) output: x-part of K_anti'] = [str(sp.simplify(b)) for b in b_anti]

# ------------------------------------------------------------------ engine-B F(Delta) and Q images used above
# Q int(A_perp,A*_perp) = -(1/2) int(dA,dA); Q int(A_V,A*_V) = -(1/2) int(dA_V,dA_V) - int(dc,A*);
# Q int(lambda*,lambda) = -int(lambda,dslash lambda); Q int(c*,c) = -int(dc,A*)  (engine B, re-derived here with D_Q)
a, b, c_ = sp.Rational(22, 3), sp.Rational(-47, 3), sp.Rational(-8, 3)      # F = a int(A,A*) + b int(l*,l) + c int(c*,c)
QF = {'int(dA,dA)': -a / 2, 'int(lambda,dslash lambda)': -b, 'int(dc,A*)': -a - c_}
res['Q F_B(Delta) (units)'] = {k: str(v) for k, v in QF.items()}
res['Ob_B(Delta) (units)'] = {'int(dA,dA)': '-11/3', 'int(lambda,dslash lambda)': '47/3', 'int(dc,A*)': '-14/3'}
res['4 C_log|quad (units)'] = {'int(dA,dA)': '-4', 'int(lambda,dslash lambda)': '16', 'int(dc,A*)': '-4'}
res['Q Pi (engine A class i), from this script'] = {'int(dA_perp,dA_perp)': '1/3 (if alpha/(2|p|^2) = 1/3)',
                                                     'int(dA_V,dA_V)': '1/3 (if alpha - target = 0)',
                                                     'int(dc,A*)': '-2/3 (if alpha - target = 0)',
                                                     'int(lambda,dslash lambda)': '-1/3 (if K_anti - target = 0)'}
with open(os.path.join(HERE, 'q_of_engineA_classi.json'), 'w') as f:
    json.dump(res, f, indent=1)
print(json.dumps(res, indent=1))
