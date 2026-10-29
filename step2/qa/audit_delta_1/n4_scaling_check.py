# Auditor #1, Part D: exact (sympy) checks of the N=4 kernel statements used by D1(c),
# in the conventions of Prop c3-P (position space kernels on R^4) and Lemma c7-weights.
# M = v. for v = 2E: on a slot of Weyl weight w (in units of E), M = 2E + 2w, E = componentwise x.d/dx.
# On a kernel k(x-y) in the slot pair (s1,s2): (M(x)1 + 1(x)M) k = (2 u.grad_u + 2w1 + 2w2) k,  u = x - y.
import sympy as sp, json, os

u = sp.symbols('u1:5', real=True)
eps, t, r2s = sp.symbols('epsilon t r2', positive=True)
r2 = sum(ui**2 for ui in u)

def udotgrad(F):
    return sum(ui * sp.diff(F, ui) for ui in u)

f = (1 - sp.exp(-r2 / (4 * eps))) / (4 * sp.pi**2 * r2)          # scalar core of P(eps,infty), Prop c3-P
k = (4 * sp.pi * t)**-2 * sp.exp(-r2 / (4 * t))                   # scalar core of K_t
res = {}

def z(expr):
    return sp.simplify(sp.expand(expr)) == 0

# P(eps,infty) components (weights in units of E; Lemma c7-weights: A 1, lambda 3/2, c 0, A* 3, lambda* 5/2, c* 4)
# A(x)A : scalar f ; c(x)A*|_V and A*|_V(x)c : grad f ; lambda(x)lambda : rho(grad) f (componentwise = grad f)
res["P_AA"] = z((2 * udotgrad(f) + 2 * (1 + 1) * f) - (-4 * eps * sp.diff(f, eps)))
res["P_cAstar"] = all(z((2 * udotgrad(sp.diff(f, ui)) + 2 * (0 + 3) * sp.diff(f, ui)) - (-4 * eps * sp.diff(sp.diff(f, ui), eps))) for ui in u)
res["P_lamlam"] = all(z((2 * udotgrad(sp.diff(f, ui)) + 2 * (sp.Rational(3, 2) + sp.Rational(3, 2)) * sp.diff(f, ui)) - (-4 * eps * sp.diff(sp.diff(f, ui), eps))) for ui in u)
# and -4 eps dP/deps = 4 eps L_eps : dP/deps = -L_eps (P = int_eps^infty L_t dt)
res["dP_deps_equals_minus_L"] = z(sp.diff(f, eps) + (4 * sp.pi * eps)**-2 * sp.exp(-r2 / (4 * eps)))
# wrong sign control
res["control_P_AA_plus_sign_fails"] = not z((2 * udotgrad(f) + 4 * f) - (+4 * eps * sp.diff(f, eps)))
# K_t components: field-antifield pairs, weights add to 4 (in units of E) => 2*4 = 8
res["K_t"] = z((2 * udotgrad(k) + 8 * k) - (-4 * t * sp.diff(k, t)))
# L_t on A(x)A: (M(x)1+1(x)M)L_t = -4 L_t - 4 t dL_t/dt
res["L_t_AA"] = z((2 * udotgrad(k) + 4 * k) - (-4 * k - 4 * t * sp.diff(k, t)))
# pairing preservation <-> weights add to div(2E)=8 (in units of 2E-weights: 2w1+2w2 = 8)
w = {"A": 1, "lam": sp.Rational(3, 2), "c": 0, "Astar": 3, "lamstar": sp.Rational(5, 2), "cstar": 4}
res["pair_weights_sum_8"] = all(2 * w[a] + 2 * w[b] == 8 for a, b in [("A", "Astar"), ("lam", "lamstar"), ("c", "cstar")])
# momentum space: E acts on Fourier components as -(4 + p.d/dp)  (f(x) = int d^4p/(2pi)^4 e^{i(p,x)} f~(p))
p = sp.symbols('p1:5', real=True)
aa = sp.symbols('a', positive=True)
gx = sp.exp(-aa * r2)
# 4d Fourier transform factorises; compute FT[g] and FT[u.grad g] by 1d integrals
def ft1(expr_u, ui, pi_):
    return sp.integrate(expr_u * sp.exp(-sp.I * pi_ * ui), (ui, -sp.oo, sp.oo))
g1 = [sp.exp(-aa * ui**2) for ui in u]
G1 = [sp.simplify(ft1(g1[i], u[i], p[i])) for i in range(4)]
xg1 = [sp.simplify(ft1(u[i] * sp.diff(g1[i], u[i]), u[i], p[i])) for i in range(4)]
Ghat = sp.Mul(*G1)
FT_Eg = sum(xg1[i] * sp.Mul(*[G1[j] for j in range(4) if j != i]) for i in range(4))
res["E_momentum_minus_4_minus_pdp"] = z(FT_Eg - (-(4 * Ghat + sum(p[i] * sp.diff(Ghat, p[i]) for i in range(4)))))
print(json.dumps({k_: bool(v) for k_, v in res.items()}, indent=1))
json.dump({k_: bool(v) for k_, v in res.items()}, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "n4_scaling_results.json"), "w"), indent=1)
