# Audit #2: remaining finite checks for the K_1 part of the design.
# M1 rank of span{eta1 eta2 + eta2 eta1 : Gamma(eta1,eta2)|_V = 0} in Sym^2 Sigma^* (design: 132 = 136 - 4), using the
#    two-component v1^m = -i(eta2 sig^m etabar1 + eta1 sig^m etabar2) (design Def def:o-K).
# M2 1/4 div is a 1-cocycle on conformal Killing fields of R^4: 1/4 div[v,w] = v(1/4 div w) - w(1/4 div v).
# M3 [K_{v1}, P_{v2}] = (v1,v2) 2E + (x -> 2(v1,x)v2 - 2(v2,x)v1) with [v,w] = v(w) - w(v) (ch07 rem:c7-KDelta).
# M4 the values of Thm thm:o-KD from gamma = 0,-24,-4,-4,-8,-2 and K_1 = -2(b.x) gamma.
import numpy as np, sympy as sp, json
rng = np.random.default_rng(5)
s1 = np.array([[0,1],[1,0]], complex); s2 = np.array([[0,-1j],[1j,0]], complex); s3 = np.array([[1,0],[0,-1]], complex)
sig = [s1, s2, s3, 1j*np.eye(2)]
def bil(e, eb, m): return np.einsum('ai,ab,bi->', e, sig[m], eb)
def split(v): return v[:8].reshape(2, 4), v[8:].reshape(2, 4)
def v1(x, y):
    ex, ebx = split(x); ey, eby = split(y)
    return np.array([-1j*(bil(ey, ebx, m)+bil(ex, eby, m)) for m in range(4)])
iu = np.triu_indices(16)
rows = []
for _ in range(400):
    x = rng.normal(size=16)+1j*rng.normal(size=16)
    A = np.array([v1(x, np.eye(16)[j]) for j in range(16)]).T   # 4 x 16, linear in eta2
    _, _, Vh = np.linalg.svd(A); ker = Vh[4:].conj().T
    y = ker@(rng.normal(size=12)+1j*rng.normal(size=12))
    assert np.abs(v1(x, y)).max() < 1e-10
    Sm = np.outer(x, y)+np.outer(y, x)
    rows.append(Sm[iu])
rk = np.linalg.matrix_rank(np.array(rows), tol=1e-8)
out = {'M1_rank': int(rk), 'M1_expected': 132}
# M2, M3
x = sp.symbols('x1:5'); b = sp.symbols('b1:5'); c = sp.symbols('c1:5')
X = sp.Matrix(x)
def div(v): return sum(sp.diff(v[i], x[i]) for i in range(4))
def brk(v, w): return sp.Matrix([sum(v[n]*sp.diff(w[k], x[n]) - w[n]*sp.diff(v[k], x[n]) for n in range(4)) for k in range(4)])
def K(bb): return sp.Matrix([sum(xi**2 for xi in x)*bb[k] - 2*sum(bb[i]*x[i] for i in range(4))*x[k] for k in range(4)])
def P(bb): return sp.Matrix(list(bb))
M = sp.Matrix(4, 4, lambda i, j: sp.Symbol('m%d%d' % (min(i, j), max(i, j)))*(1 if i < j else (-1 if i > j else 0)))
def Rot(): return M*X
def Dil(): return 2*X
gens = {'K_b': K(b), 'K_c': K(c), 'P_b': P(b), 'P_c': P(c), 'R': Rot(), 'D': Dil()}
ok = True
for n1, v in gens.items():
    for n2, w in gens.items():
        lhs = div(brk(v, w))/4
        rhs = sum(v[n]*sp.diff(div(w)/4, x[n]) - w[n]*sp.diff(div(v)/4, x[n]) for n in range(4))
        ok &= sp.simplify(lhs - rhs) == 0
out['M2_quarter_div_cocycle'] = bool(ok)
KP = brk(K(b), P(c))
pred = sum(b[i]*c[i] for i in range(4))*2*X + sp.Matrix([2*sum(b[i]*x[i] for i in range(4))*c[k] - 2*sum(c[i]*x[i] for i in range(4))*b[k] for k in range(4)])
out['M3_[K,P]'] = bool(sp.simplify(KP - pred) == sp.zeros(4, 1))
out['M3_quarter_div_K'] = str(sp.simplify(div(K(b))/4))
gam = {'traceless': 0, 'trace': -24, 'iota_w A': -4, 'A_V': -4, 'lambda': -8, 'c': -2}
out['M4_K1_coeff_of_(b.x)unit'] = {k: -2*v for k, v in gam.items()}
print(json.dumps(out, indent=1))
json.dump(out, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/design_audit_2/misc_checks.json', 'w'), indent=1)
