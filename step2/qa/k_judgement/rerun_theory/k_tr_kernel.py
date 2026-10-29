# k_tr_kernel.py -- the floating special-conformal insertion at a finite transfer endpoint T (k_transfer.tex).
#
# Objects (manuscript names; code names only inside this file):
#   D_T(r) = int_eps^T (4 pi t)^{-2} e^{-r^2/4t} dt = (e^{-r^2/4T} - e^{-r^2/4eps})/(4 pi^2 r^2)  (scalar factor of P(eps,T)),
#   v(y)   = (y,y) b - 2 (b,y) y  (special conformal field along b = v_1),
#   F_k(z_a,z_b;eps,T) = int d^4y sum_m [L_v(D_T(.-z_a) e_k)]_m(y) d_{y_m} D_T(y-z_b),
#   the composite of the floating term (L_v A|_V, A^*|_V) of S^(1)(K_v) between an A(x)A line ending at z_a and an
#   A^*|_V (x) c line ending at z_b, with the fibre factors 2(-,-)^{-1} and -(-,-)^{-1} removed [design eq:o-Kir].
# Parts:
#   A  exact: T dF_k/dT at z_a = z_b = 0 as a closed function of T/eps and its limit T -> infinity;
#      exact F_k(T) - F_k(T_0).
#   B  z_a, z_b != 0: T dF_k/dT by exact Gaussian moments in y and a one-dimensional Schwinger integral; the limit is
#      independent of z_a, z_b (the zero-mode statement).
#   C  non-invariance of K_t and P(eps,T) under the special conformal field: exact symbolic identities.
#   D  the special edge: the factor e^{-eps|p|^2} - 1 at a generic momentum, |.| <= eps |p|^2.
import json
import sympy as sp
import mpmath as mp
import time
T0clock = time.time()

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_judgement/rerun_theory/k_tr_kernel.json'
res = {}

# ------------------------------------------------------------------ A: exact, z_a = z_b = 0
u, T, e, tau = sp.symbols('u T epsilon tau', positive=True)
D = (sp.exp(-u/(4*T)) - sp.exp(-u/(4*e)))/(4*sp.pi**2*u)          # D_T as a function of u = r^2
TD = sp.exp(-u/(4*T))/(16*sp.pi**2*T)                              # T dD_T/dT = T (4 pi T)^{-2} e^{-u/4T}
# radial reduction (sphere average <y_k y_l> = r^2 delta_kl/4):
#   sum_m [L_v(f e_k)]_m d_m g = (v.df)(d_k g) + f (d_m v^k)(d_m g),  v.y = -r^2 (b.y),  (d_m v^k) y_m = 2 r^2 b_k - 4(b.y) y_k
#   sphere average: b_k [ -(r^2/4) f' g' + r f g' ];  d^4y = 2 pi^2 r^3 dr;  u = r^2:  F_k = 2 pi^2 b_k int du [-(1/2)u^3 D_u^2 + u^2 D D_u]
Du, TDu = sp.diff(D, u), sp.diff(TD, u)
dens_T = -sp.Rational(1, 2)*u**3*2*Du*TDu + u**2*(TD*Du + D*TDu)
TdF = sp.simplify(2*sp.pi**2*sp.integrate(sp.expand(dens_T), (u, 0, sp.oo), conds='none'))
TdF_tau = sp.simplify(sp.factor(TdF.subs(T, tau*e)))
lim = sp.limit(TdF, T, sp.oo)
res['A_TdFdT_over_b_k_at_z0'] = {'closed_form_in_tau=T/eps': str(TdF_tau), 'limit_T_to_infinity': str(sp.simplify(lim)),
                                  'expected_half_of_logR_coefficient_-3/(8pi^2)': str(sp.Rational(-3, 16)/sp.pi**2),
                                  'agree': bool(sp.simplify(lim + sp.Rational(3, 16)/sp.pi**2) == 0)}
# exact F_k(T) - F_k(T0) = int_{T0}^{T} (T' dF/dT') dT'/T'
asym = sp.simplify(sp.series(TdF_tau, tau, sp.oo, 3).removeO())
res['A_large_tau_expansion_of_TdF/dT'] = str(asym)
# the eps -> 0 limit at fixed T of the same quantity (UV end): the T-dependence of the composite survives eps -> 0
res['A_eps_to_0_at_fixed_T'] = str(sp.simplify(sp.limit(TdF, e, 0)))

print('A done', time.time()-T0clock, flush=True)
# ------------------------------------------------------------------ B: z_a, z_b != 0
y = sp.symbols('y1:5', real=True)
Y = sp.symbols('Y1:5', real=True)
tf, tg, sig2 = sp.symbols('t_f t_g sigma2', positive=True)
bvec = [sp.Rational(7, 10), sp.Rational(-11, 10), sp.Rational(2, 5), sp.Rational(9, 10)]
za = [sp.Rational(3, 10), sp.Rational(-1, 5), sp.Rational(1, 2), sp.Rational(1, 10)]
zb = [sp.Rational(-2, 5), sp.Rational(1, 4), sp.Rational(3, 20), sp.Rational(-3, 10)]
yy = sum(t**2 for t in y)
by = sum(bvec[i]*y[i] for i in range(4))
v = [yy*bvec[n] - 2*by*y[n] for n in range(4)]
dv = [[sp.diff(v[k], y[m]) for k in range(4)] for m in range(4)]   # dv[m][k] = d_m v^k


import numpy as np
from scipy import integrate as si
bn = np.array([float(q) for q in bvec]); zan = np.array([float(q) for q in za]); zbn = np.array([float(q) for q in zb])
# 4-d Gauss-Hermite with 3 nodes per axis: exact for polynomials of degree <= 5 in each variable
gh_x, gh_w = np.polynomial.hermite_e.hermegauss(3)          # weight e^{-x^2/2}
GX = np.array(np.meshgrid(gh_x, gh_x, gh_x, gh_x, indexing='ij')).reshape(4, -1).T
GW = np.prod(np.array(np.meshgrid(gh_w, gh_w, gh_w, gh_w, indexing='ij')).reshape(4, -1), 0)/(2*np.pi)**2


def pair_integral(k, tfv, tgv):
    """int d^4y sum_m [L_v(G_tf(.-z_a) e_k)]_m d_m G_tg(.-z_b), exactly (Gaussian moments, polynomial degree 4)."""
    alpha = 1/(4*tfv) + 1/(4*tgv)
    y0 = (zan/tfv + zbn/tgv)/(1/tfv + 1/tgv)
    beta = ((zan - zbn)**2).sum()/(4*(tfv + tgv))
    s = np.sqrt(1/(2*alpha))
    Yp = y0[None, :] + s*GX
    yy_ = (Yp**2).sum(1); by_ = Yp@bn
    vv = yy_[:, None]*bn[None, :] - 2*by_[:, None]*Yp                 # v(y)
    df = -(Yp - zan[None, :])/(2*tfv); dg = -(Yp - zbn[None, :])/(2*tgv)
    vdf = (vv*df).sum(1)
    # d_m v^k = 2 y_m b_k - 2 b_m y_k - 2 (b.y) delta_mk
    dvk = 2*Yp*bn[k] - 2*bn[None, :]*Yp[:, k:k+1] - 2*by_[:, None]*(np.arange(4) == k)[None, :]
    poly = vdf*dg[:, k] + (dvk*dg).sum(1)
    pref = (4*np.pi*tfv)**-2*(4*np.pi*tgv)**-2*(np.pi/alpha)**2*np.exp(-beta)
    return pref*(poly*GW).sum()


resB = {}
for k in range(4):
    vals = {}
    for Tval in [1e2, 1e4, 1e6, 1e8]:
        epsv = 1.0
        g = lambda w: np.exp(w)*Tval*(pair_integral(k, Tval, np.exp(w)) + pair_integral(k, np.exp(w), Tval))
        val = si.quad(g, np.log(epsv), np.log(Tval), limit=400, epsabs=0, epsrel=1e-11)[0]
        vals['T=%g' % Tval] = val
    resB['k=%d' % (k + 1)] = {'TdF/dT': vals, 'prediction_-3b_k/(16pi^2)': float(-3*bvec[k]/(16*sp.pi**2))}
# control at z_a = z_b = 0 against the exact closed form of part A
zan0, zbn0 = zan.copy(), zbn.copy()
zan[:] = 0; zbn[:] = 0
ctrl = {}
for Tval in [1e2, 1e4]:
    g = lambda w: np.exp(w)*Tval*(pair_integral(0, Tval, np.exp(w)) + pair_integral(0, np.exp(w), Tval))
    num = si.quad(g, 0.0, np.log(Tval), limit=400, epsabs=0, epsrel=1e-11)[0]
    exact = float(TdF.subs({T: Tval, e: 1})*bvec[0])
    ctrl['T=%g' % Tval] = {'gauss_moment_route': num, 'closed_form_A': exact}
zan[:] = zan0; zbn[:] = zbn0
res['B_z_nonzero_eps=1'] = {'b': [str(x) for x in bvec], 'z_a': [str(x) for x in za], 'z_b': [str(x) for x in zb],
                            'values': resB, 'control_z=0_vs_part_A': ctrl}

print('B done', time.time()-T0clock, flush=True)
# ------------------------------------------------------------------ C: non-invariance under K_v
# Exact identities, checked at random rational points with 40-digit arithmetic (the remainders are derived by hand
# in k_transfer.tex, Lemma lem:ktr-noninv; nothing is fitted).
import random
random.seed(20261028)
mp.mp.dps = 40
x = sp.symbols('x1:5', real=True)
w_ = sp.symbols('w', real=True)
t = sp.symbols('t', positive=True)
r = [x[i] - y[i] for i in range(4)]
rr = sum(ri**2 for ri in r)
bm = sum(bvec[i]*(x[i] + y[i])/2 for i in range(4))      # (b, midpoint); (1/8 div v)(midpoint) = -(b, midpoint)
br_ = sum(bvec[i]*r[i] for i in range(4))
xx = sum(t_**2 for t_ in x)
bx = sum(bvec[i]*x[i] for i in range(4))
vx = [xx*bvec[n] - 2*bx*x[n] for n in range(4)]
vy = v
divx = sum(sp.diff(vx[n], x[n]) for n in range(4))
divy = sum(sp.diff(vy[n], y[n]) for n in range(4))
Gt = lambda tt: (4*sp.pi*tt)**-2*sp.exp(-rr/(4*tt))
DT = (sp.exp(-rr/(4*T)) - sp.exp(-rr/(4*e)))/(4*sp.pi**2*rr)


def maxres(expr, nsample=6):
    f = sp.lambdify(list(x) + list(y) + [T, e, t, w_], expr, 'mpmath')
    m = mp.mpf(0)
    for _ in range(nsample):
        args = [mp.mpf(random.randint(-9, 9))/7 for _ in range(8)] + [mp.mpf(random.randint(20, 60))/10, mp.mpf(random.randint(1, 9))/10,
                                                                         mp.mpf(random.randint(1, 30))/10, mp.mpf(random.randint(0, 8))/2]
        m = max(m, abs(f(*args)))
    return m


resC = {}
G = Gt(t)
lhs = sum(vx[n]*sp.diff(G, x[n]) for n in range(4)) + w_/4*divx*G + sum(vy[n]*sp.diff(G, y[n]) for n in range(4)) + (4 - w_)/4*divy*G
rhs = (-bm)*(-4*t*sp.diff(G, t)) + (4 - 2*w_)*br_*G
resC['C1_heat_kernel_weight_w'] = {'identity': '(v.(x)1+1(x)v.)K_t|_{field weight w at x, antifield weight 4-w at y} = (1/8 div v)(midpoint) (-4 t dK_t/dt) + (4-2w)(b,x-y) K_t',
                                   'max_residual_40digits': mp.nstr(maxres(lhs - rhs), 5)}
Pphi = 2*DT
lhs2 = sum(vx[n]*sp.diff(Pphi, x[n]) for n in range(4)) + sp.Rational(1, 4)*divx*Pphi \
    + sum(vy[n]*sp.diff(Pphi, y[n]) for n in range(4)) + sp.Rational(1, 4)*divy*Pphi
rhs2 = (-bm)*(4*e*2*Gt(e) - 4*T*2*Gt(T))
resC['C2_P_scalar'] = {'identity': '(v.(x)1+1(x)v.)P(eps,T)|_{A_Wperp(x)A_Wperp} = (1/8 div v)(midpoint)(4 eps L_eps - 4 T L_T)',
                       'max_residual_40digits': mp.nstr(maxres(lhs2 - rhs2), 5)}
m3 = mp.mpf(0)
for m_ in range(4):
    for n_ in range(4):
        Pmn = lambda a, c: (2*DT if a == c else 0)
        val = sum(vx[q]*sp.diff(Pmn(m_, n_), x[q]) for q in range(4)) + sum(sp.diff(vx[k], x[m_])*Pmn(k, n_) for k in range(4)) \
            + sum(vy[q]*sp.diff(Pmn(m_, n_), y[q]) for q in range(4)) + sum(sp.diff(vy[k], y[n_])*Pmn(m_, k) for k in range(4))
        pred = (-bm)*(4*e*(2*Gt(e) if m_ == n_ else 0) - 4*T*(2*Gt(T) if m_ == n_ else 0)) + 4*(r[m_]*bvec[n_] - r[n_]*bvec[m_])*DT
        m3 = max(m3, maxres(val - pred, 3))
resC['C3_P_vector'] = {'identity': '(v.(x)1+1(x)v.)P(eps,T)|_{A_m(x)A_n} = (1/8 div v)(midpoint)(4 eps L_eps - 4 T L_T)_{mn} + 4((x-y)_m b_n - (x-y)_n b_m) D_T',
                       'max_residual_40digits': mp.nstr(m3, 5)}
m4 = mp.mpf(0)
for m_ in range(4):
    Pc = lambda a: -sp.diff(DT, x[a])
    val = sum(vx[q]*sp.diff(Pc(m_), x[q]) for q in range(4)) \
        + sum(vy[q]*sp.diff(Pc(m_), y[q]) for q in range(4)) + divy*Pc(m_) - sum(sp.diff(vy[m_], y[n_])*Pc(n_) for n_ in range(4))
    Lc = lambda tt, a: -sp.diff(Gt(tt), x[a])
    pred = (-bm)*(4*e*Lc(e, m_) - 4*T*Lc(T, m_)) - sum((2*r[m_]*bvec[q] + bvec[m_]*r[q])*sp.diff(DT, x[q]) for q in range(4))
    m4 = max(m4, maxres(val - pred, 3))
resC['C4_P_ghost_line'] = {'identity': '(v.(x)1+1(x)v.)P(eps,T)|_{c(x)A^*_m(y)} = (1/8 div v)(midpoint)(4 eps L_eps - 4 T L_T)_{c A^*_m} - (2(x-y)_m b_q + b_m (x-y)_q) d_{x_q} D_T  [fibre factor -(,)^{-1} removed]',
                           'max_residual_40digits': mp.nstr(m4, 5)}
res['C_noninvariance'] = resC

# ------------------------------------------------------------------ D: special edge factor
pp, ee = sp.symbols('p2 eps', positive=True)
res['D_special_edge'] = {'factor': 'e^{-eps|p|^2}-1', 'bound': '|e^{-eps|p|^2}-1| <= eps|p|^2',
                         'series': str(sp.series(sp.exp(-ee*pp) - 1, ee, 0, 3))}

json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
