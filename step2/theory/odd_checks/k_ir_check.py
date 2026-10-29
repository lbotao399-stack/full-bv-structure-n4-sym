# Large-distance check of the floating S^(1)(K) insertion on an A|_V -> c line (odd_components_design.tex, Prop. prop:o-ir).
# Tails (ch03 Prop c3-P, |x-y| >> sqrt(eps)): A(x)A(y): 2 D(x-y) delta_{mn}; A*_V(x) c(y): -d_x D(x-y); D(r)=1/(4 pi^2 r^2).
# Vertex (L_v A|_V, A*_V), v(y)=(y,y)b-2(b,y)y. Integrand (propagator factors 2 and -1 omitted):
#   F_k(y) = sum_m [L_v(D(.-za) e_k)]_m(y) d_{y_m} D(y-zb),  [L_v w]_m = v^n d_n w_m + w_n d_m v^n.
import sympy as sp, numpy as np
y = sp.symbols('y1:5', real=True); bb = sp.symbols('b1:5', real=True)
c = 1/(4*sp.pi**2)
r2 = sum(t**2 for t in y)
yv = sp.Matrix(y); b = sp.Matrix(bb)
v = r2*b - 2*(b.dot(yv))*yv
D0 = c/r2
lead = []
for k in range(4):
    tot = 0
    for m in range(4):
        Lv = (1 if m == k else 0)*sum(v[n]*sp.diff(D0, y[n]) for n in range(4)) + D0*sp.diff(v[k], y[m])
        tot += Lv*sp.diff(D0, y[m])
    lead.append(sp.together(sp.expand(tot)))
print('F_0 at za=zb=0 :', sp.factor(lead[0]))
# angular average on S^3 of a polynomial numerator: <prod y_i^{2a_i}> = Gamma(2) prod Gamma(a_i+1/2) / (Gamma(1/2)^4 Gamma(2+sum a_i))
def sphere_avg_poly(num):
    P = sp.Poly(sp.expand(num), *y); out = 0
    for mon, co in P.terms():
        if any(e % 2 for e in mon): continue
        a = [e//2 for e in mon]
        out += co*sp.gamma(2)*sp.prod([sp.gamma(sp.Rational(2*ai+1,2)) for ai in a])/(sp.gamma(sp.Rational(1,2))**4*sp.gamma(2+sum(a)))
    return sp.nsimplify(sp.simplify(out))
num0, den0 = sp.fraction(lead[0])
# den0 = const * r2^p ; evaluate on |y|=1
den_on_sphere = sp.simplify(den0.subs({y[0]:1, y[1]:0, y[2]:0, y[3]:0}))
avg0 = sp.simplify(sphere_avg_poly(num0)/den_on_sphere)
print('S^3-average of |y|^4 F_0 :', avg0, '   ratio to c^2 b_1 :', sp.simplify(avg0/(c**2*bb[0])))
# numerical check with anchored points za, zb != 0: average of |y|^4 F_k over S^3 at radius R -> -3 c^2 b_k
za = np.array([0.3,-0.2,0.5,0.1]); zb = np.array([-0.4,0.25,0.15,-0.3]); bnum = np.array([0.7,-1.1,0.4,0.9])
cn = 1/(4*np.pi**2)
def Fnum(Y):
    da = Y-za; db = Y-zb; ra2 = (da**2).sum(-1); rb2 = (db**2).sum(-1); yy = (Y**2).sum(-1); by = Y@bnum
    vY = yy[:,None]*bnum[None,:] - 2*by[:,None]*Y
    gradDa = -2*cn*da/ra2[:,None]**2; Da = cn/ra2; gradDb = -2*cn*db/rb2[:,None]**2
    # d_m v^n = 2 y_m b^n - 2 b_m y^n - 2 (b.y) delta_mn
    out = np.zeros((Y.shape[0],4))
    vdotgradDa = (vY*gradDa).sum(-1)
    for k in range(4):
        term1 = vdotgradDa*gradDb[:,k]
        dmvk = 2*Y*bnum[k] - 2*bnum[None,:]*Y[:,k][:,None] - 2*by[:,None]*np.eye(4)[k][None,:]
        term2 = Da*(dmvk*gradDb).sum(-1)
        out[:,k] = term1+term2
    return out
rng = np.random.default_rng(1)
G = rng.normal(size=(400000,4)); U = G/np.linalg.norm(G,axis=1)[:,None]
for Rr in [10.,100.,1000.]:
    vals = Fnum(Rr*U)*Rr**4
    print('R=%7.0f  <|y|^4 F_k> / (c^2 b_k) =' % Rr, np.round(vals.mean(0)/(cn**2*bnum),4))
