# Audit #2: independent check of design eq:o-Kir (Prop prop:o-ir) by a second method, and calibration of the
# momentum-space IR method used in eta_pair_ir.py.
# (1) position space: y = Y/u, u -> 0, exact rational anchors za, zb (two different sets), leading u^4 coefficient,
#     exact S^3 average from <prod y_i^(2a_i)> = prod (2a_i-1)!! / (4*6*...*(2n+2)).
# (2) momentum space: B_K = int dp/(2pi)^4 [FT(L_v alpha)](-p) . Ghat(p), alpha_m = delta_mk (D*f), G_m = (d_m D)*g,
#     FT conventions hhat(p)=int e^{-ipy} h, (y_a h)^ = i d_{p_a} hhat, (d_a h)^ = i p_a hhat, D^ = 1/p^2.
#     fhat(q) = f0 + f1.q, ghat(p) = g0 + g1.p (Taylor data of the rest at the IR point); the degree -4 part is the
#     f0*g0 part. The coefficient of log R (|y|<R <-> |p|>1/R) is (2pi^2/(2pi)^4) <p^4 I_-4>; position space gives
#     2pi^2 (-3 b_k)(4pi^2)^-2 f0 g0 = -3 b_k f0 g0/(8 pi^2).
import sympy as sp, json
from fractions import Fraction

out = {}
def dfact(n): return 1 if n <= 0 else n*dfact(n-2)
def sphere_avg(poly, vars_):
    P = sp.Poly(sp.expand(poly), *vars_); res = 0
    for mon, co in P.terms():
        if any(e % 2 for e in mon): continue
        n = sum(mon)//2
        den = 1
        for j in range(n): den *= (4+2*j)
        res += co*sp.Integer(sp.prod([dfact(e-1) for e in mon]))/den
    return sp.nsimplify(res)

# ---------- (1) position space ----------
ys = sp.symbols('y1:5', real=True); Y = sp.symbols('Y1:5', real=True); u = sp.symbols('u', positive=True)
b = [sp.Rational(7,10), sp.Rational(-11,10), sp.Rational(2,5), sp.Rational(9,10)]
cst = 1/(4*sp.pi**2)
def D(x): return cst/sum(xi**2 for xi in x)
vv = [sum(s**2 for s in ys)*b[n] - 2*sum(b[i]*ys[i] for i in range(4))*ys[n] for n in range(4)]
res_pos = {}
for label, za, zb in [('anchors0', [0,0,0,0], [0,0,0,0]),
                      ('anchorsA', [sp.Rational(3,10),sp.Rational(-1,5),sp.Rational(1,2),sp.Rational(1,10)], [sp.Rational(-2,5),sp.Rational(1,4),sp.Rational(3,20),sp.Rational(-3,10)]),
                      ('anchorsB', [1,2,-1,0], [0,-3,1,2])]:
    Da = D([ys[i]-za[i] for i in range(4)]); Db = D([ys[i]-zb[i] for i in range(4)])
    avgs = []
    for k in range(4):
        tot = 0
        for m in range(4):
            Lv = (1 if m == k else 0)*sum(vv[n]*sp.diff(Da, ys[n]) for n in range(4)) + Da*sp.diff(vv[k], ys[m])
            tot += Lv*sp.diff(Db, ys[m])
        Fu = sp.together(tot.subs({ys[i]: Y[i]/u for i in range(4)}))
        # leading behaviour: u^4 * F at u=0
        lead = sp.cancel(sp.together(Fu/u**4)).subs(u, 0)
        lead = sp.cancel(lead)
        R2 = sum(s**2 for s in Y)
        poly = sp.cancel(lead*R2**4)
        avgs.append(sphere_avg(poly, Y))
    res_pos[label] = {'S3_avg_of_|y|^4F_k_over_(4pi^2)^-2': [str(sp.nsimplify(a/cst**2)) for a in avgs],
                      'b': [str(x) for x in b]}
out['position_space'] = res_pos
out['position_space']['expected'] = '-3 b_k'

# ---------- (2) momentum space calibration ----------
p = sp.symbols('p1:5', real=True); q = sp.symbols('q1:5', real=True)
f0, g0 = sp.symbols('f0 g0')
P2 = lambda x: sum(s**2 for s in x)
def ymul(a, expr): return sp.I*sp.diff(expr, q[a])
cal = []
for k in range(4):
    alpha = [(f0/P2(q) if m == k else 0) for m in range(4)]
    def FT_Lv(m):
        tot = 0
        for n in range(4):
            h = sp.I*q[n]*alpha[m]
            tot += sum(b[n]*ymul(j, ymul(j, h)) for j in range(4)) - 2*sum(b[j]*ymul(j, ymul(n, h)) for j in range(4))
        for n in range(4):
            h = alpha[n]
            tot += 2*b[n]*ymul(m, h) - 2*b[m]*ymul(n, h) - (2*sum(b[j]*ymul(j, h) for j in range(4)) if m == n else 0)
        return tot
    Ghat = [sp.I*p[m]*g0/P2(p) for m in range(4)]
    I = sum(FT_Lv(m).subs({q[i]: -p[i] for i in range(4)})*Ghat[m] for m in range(4))
    polyK = sp.cancel(sp.together(I)*P2(p)**4)
    avgK = sphere_avg(polyK, p)
    logcoef = sp.simplify(2*sp.pi**2/(2*sp.pi)**4*avgK)
    target = sp.simplify(-3*b[k]/(8*sp.pi**2)*f0*g0)
    cal.append({'k': k+1, 'log_R_coefficient': str(logcoef), 'position_space_prediction': str(target),
                'agree': bool(sp.simplify(logcoef - target) == 0)})
out['momentum_calibration'] = cal
out['design_factor_6*2pi^2*(4pi^2)^-2'] = str(sp.nsimplify(6*2*sp.pi**2/(4*sp.pi**2)**2))
print(json.dumps(out, indent=1))
json.dump(out, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/design_audit_2/k_tail_and_calibration.json', 'w'), indent=1)
