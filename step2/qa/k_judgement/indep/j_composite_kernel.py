# j_composite_kernel.py -- judge's independent check of the composite A|_V--c line along K_b at a finite endpoint T.
#
# Written from the definitions only (ch03 Prop c3-P table, ch08 Def c8-PT, ch07 Def c7-vdot); no code of the theory
# agent (step2/theory/k_tr_*.py) or of the numerics agent (step2/qa/k_transfer_numeric/) is imported.
#
# Objects.  D(r) = (e^{-r^2/4T} - e^{-r^2/4eps})/(4 pi^2 r^2) is the scalar factor of P(eps,T); v(y) = |y|^2 b - 2(b,y)y.
#   F_k(z_a,z_b) = int d^4y sum_m [L_v(D(.-z_a) e_k)]_m(y) d_{y_m} D(y-z_b),   (L_v w)_m = v^n d_n w_m + w_n d_m v^n.
# Composite line (judge's derivation, see k_judgement.tex): with the A(x)A component 2(-,-)^{-1} D and the A^*|_V(x)c
# component -(-,-)^{-1} d_x D of ch03 Prop c3-P (value of d_P on the product of the coordinate functionals in the order
# of the slots), the kernel CL in the slots (A|_V at z_a, c at z_b) defined by [d_P,[d_P,S]] = d_{2 CL},
# S = int (L_v A|_V, A^*|_V), is  CL_k(z_a,z_b) = 2 * (-1) * (-,-)^{-1} F_k(z_a,z_b) = -2 (-,-)^{-1} F_k.
#
# Checks (all independent implementations):
#  (1) T dF_k/dT at z_a = z_b = 0 by direct 4-d quadrature of the Cartesian integrand (no sphere-average identity):
#      radial Gauss-Legendre in log r, exact angular quadrature on S^3 (Hopf coordinates), at tau = T/eps = 2, 10, 100.
#  (2) the same quantity in closed form by the judge's own radial reduction (sympy), compared with (1) and with the
#      theory's closed form -(tau-1)(3tau^2+8tau+3)/(16pi^2(tau+1)^3); limit tau -> infinity.
#  (3) z_a, z_b != 0: T dF_k/dT = T int_eps^T dt [Phi_k(T,t) + Phi_k(t,T)], with Phi_k(t_f,t_g) the Gaussian integral of the
#      polynomial integrand computed EXACTLY by Isserlis moments in sympy (not by quadrature); 1-d integral in log t by
#      mpmath; T up to 1e8; limit compared with -3 b_k/(16 pi^2).
#  (4) the orientation: CL limit = -2 x (limit of F_k) = +3/(8 pi^2) b_k (-,-)^{-1}.
import json
import time
import itertools
import numpy as np
import sympy as sp
import mpmath as mp

t_start = time.time()
OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_judgement/indep/j_composite_kernel.json'
res = {}
bvec = np.array([0.7, -1.1, 0.4, 0.9])


# ------------------------------------------------------------------ (1) direct 4-d quadrature at z = 0
def Dfun(r2, eps, T):
    return (np.exp(-r2/(4*T)) - np.exp(-r2/(4*eps)))/(4*np.pi**2*r2)


def dD_dr2(r2, eps, T):
    # d/d(r^2) of D
    a = (-np.exp(-r2/(4*T))/(4*T) + np.exp(-r2/(4*eps))/(4*eps))/(4*np.pi**2*r2)
    return a - Dfun(r2, eps, T)/r2


def TdD(r2, T):
    # T d/dT D = e^{-r^2/4T}/(16 pi^2 T)
    return np.exp(-r2/(4*T))/(16*np.pi**2*T)


def dTdD_dr2(r2, T):
    return -np.exp(-r2/(4*T))/(64*np.pi**2*T**2)


def s3_nodes(n1=8, n2=8):
    """Hopf coordinates y = (cos a cos p1, cos a sin p1, sin a cos p2, sin a sin p2), measure sin a cos a da dp1 dp2;
    with s = sin^2 a in [0,1]: measure (1/2) ds dp1 dp2, total 2 pi^2.  Exact for polynomials of low degree."""
    xs, ws = np.polynomial.legendre.leggauss(n1)
    s = 0.5*(xs + 1); w = 0.5*ws
    ph = 2*np.pi*np.arange(n2)/n2
    pts, wts = [], []
    for si, wi in zip(s, w):
        ca, sa = np.sqrt(1 - si), np.sqrt(si)
        for p1 in ph:
            for p2 in ph:
                pts.append([ca*np.cos(p1), ca*np.sin(p1), sa*np.cos(p2), sa*np.sin(p2)])
                wts.append(0.5*wi*(2*np.pi/n2)**2)
    return np.array(pts), np.array(wts)


NOD, WT = s3_nodes()
assert abs(WT.sum() - 2*np.pi**2) < 1e-12


def integrand_z0(y, eps, T, which):
    """sum_m [L_v(f e_k)]_m d_m g for all k at points y (N,4); f, g radial; which = 'f' (f pinned at T) or 'g'."""
    r2 = (y**2).sum(1)
    by = y @ bvec
    v = r2[:, None]*bvec[None, :] - 2*by[:, None]*y
    if which == 'f':
        f, fp = TdD(r2, T), dTdD_dr2(r2, T)
        g, gp = Dfun(r2, eps, T), dD_dr2(r2, eps, T)
    else:
        f, fp = Dfun(r2, eps, T), dD_dr2(r2, eps, T)
        g, gp = TdD(r2, T), dTdD_dr2(r2, T)
    df = 2*y*fp[:, None]          # grad f
    dg = 2*y*gp[:, None]
    vdf = (v*df).sum(1)
    out = np.zeros((y.shape[0], 4))
    for k in range(4):
        # d_m v^k = 2 y_m b_k - 2 b_m y_k - 2 (b.y) delta_mk
        dvk = 2*y*bvec[k] - 2*bvec[None, :]*y[:, k:k+1] - 2*by[:, None]*(np.arange(4) == k)[None, :]
        out[:, k] = vdf*dg[:, k] + f*(dvk*dg).sum(1)
    return out


def TdF_z0_quadrature(eps, T):
    xr, wr = np.polynomial.legendre.leggauss(40)
    lo, hi = np.log(1e-4*np.sqrt(eps)), np.log(40*np.sqrt(T))
    edges = np.linspace(lo, hi, 60)
    tot = np.zeros(4)
    for a, c in zip(edges[:-1], edges[1:]):
        lr = 0.5*(c - a)*xr + 0.5*(c + a); wl = 0.5*(c - a)*wr
        for lri, wli in zip(lr, wl):
            r = np.exp(lri)
            y = r*NOD
            val = integrand_z0(y, eps, T, 'f') + integrand_z0(y, eps, T, 'g')
            tot += wli*r**4*(val*WT[:, None]).sum(0)          # d^4y = r^3 dr dOmega = r^4 dlogr dOmega
    return tot


def closed_form_theory(tau):
    return -(tau - 1)*(3*tau**2 + 8*tau + 3)/(16*np.pi**2*(tau + 1)**3)


r1 = {}
for tau in [2.0, 10.0, 100.0]:
    q = TdF_z0_quadrature(1.0, tau)
    r1['tau=%g' % tau] = {'TdF_k/b_k_quadrature': (q/bvec).tolist(), 'theory_closed_form': closed_form_theory(tau),
                          'max_rel_dev': float(np.max(np.abs(q/bvec - closed_form_theory(tau)))/abs(closed_form_theory(tau)))}
res['1_direct_4d_quadrature_z0'] = r1
print('(1) done %.1fs' % (time.time() - t_start), flush=True)

# ------------------------------------------------------------------ (2) judge's own radial reduction, exact
u, Ts, es = sp.symbols('u T epsilon', positive=True)
Dsym = (sp.exp(-u/(4*Ts)) - sp.exp(-u/(4*es)))/(4*sp.pi**2*u)
TDsym = Ts*sp.diff(Dsym, Ts)
# integrand at z = 0 after the angular average (judge's derivation, k_judgement.tex): with f = f(u), g = g(u), u = |y|^2,
# sum_m [L_v(f e_k)]_m d_m g  ->  b_k [ -u^2 f' g' + 2 u f g' ]   (average over S^3),   d^4y = pi^2 u du.
def radial(f, g):
    return sp.pi**2*u*(-u**2*sp.diff(f, u)*sp.diff(g, u) + 2*u*f*sp.diff(g, u))


dens = radial(TDsym, Dsym) + radial(Dsym, TDsym)
TdF_exact = sp.simplify(sp.integrate(sp.expand(dens), (u, 0, sp.oo), conds='none'))
tau = sp.symbols('tau', positive=True)
TdF_tau = sp.factor(sp.simplify(TdF_exact.subs(Ts, tau*es)))
theory_tau = -(tau - 1)*(3*tau**2 + 8*tau + 3)/(16*sp.pi**2*(tau + 1)**3)
res['2_radial_exact'] = {'judge_closed_form_TdF_k/b_k': str(TdF_tau),
                         'equals_theory_closed_form': bool(sp.simplify(TdF_tau - theory_tau) == 0),
                         'limit_tau_to_infinity': str(sp.limit(TdF_tau, tau, sp.oo)),
                         'value_at_tau=1': str(sp.simplify(TdF_tau.subs(tau, 1)))}
print('(2) done %.1fs' % (time.time() - t_start), flush=True)

# ------------------------------------------------------------------ (3) z != 0, exact Gaussian moments (Isserlis)
xi = sp.symbols('xi1:5', real=True)
y0 = sp.symbols('y01:5', real=True)
za_s = sp.symbols('za1:5', real=True)
zb_s = sp.symbols('zb1:5', real=True)
tf, tg = sp.symbols('t_f t_g', positive=True)
bS = [sp.Rational(7, 10), sp.Rational(-11, 10), sp.Rational(2, 5), sp.Rational(9, 10)]
yv = [y0[i] + xi[i] for i in range(4)]
yy = sum(t**2 for t in yv)
by_ = sum(bS[i]*yv[i] for i in range(4))
vv = [yy*bS[n] - 2*by_*yv[n] for n in range(4)]
# f = G_tf(y - za), g = G_tg(y - zb):  d f = -(y-za)/(2 tf) f,  d g = -(y-zb)/(2 tg) g
dfo = [-(yv[i] - za_s[i])/(2*tf) for i in range(4)]
dgo = [-(yv[i] - zb_s[i])/(2*tg) for i in range(4)]
sig2 = sp.symbols('sigma2', positive=True)


def isserlis(poly):
    """E[poly(xi)] for xi ~ N(0, sigma2 I_4), poly a polynomial in xi (exact)."""
    P = sp.Poly(sp.expand(poly), *xi)
    tot = 0
    for mon, cf in P.terms():
        if any(m % 2 for m in mon):
            continue
        val = cf
        for m in mon:
            # E[xi^m] = sigma2^{m/2} (m-1)!!
            val *= sig2**(m//2)*sp.factorial2(m - 1) if m > 0 else 1
        tot += val
    return sp.expand(tot)


Phi_poly = []
for k in range(4):
    dvk = [2*yv[m]*bS[k] - 2*bS[m]*yv[k] - 2*by_*(1 if m == k else 0) for m in range(4)]
    poly = sum(vv[n]*dfo[n] for n in range(4))*dgo[k] + sum(dvk[m]*dgo[m] for m in range(4))
    Phi_poly.append(isserlis(poly))
print('(3) moments done %.1fs' % (time.time() - t_start), flush=True)
# Gaussian prefactor: f g = (4 pi tf)^-2 (4 pi tg)^-2 exp(-|y-za|^2/4tf - |y-zb|^2/4tg)
#   = (4 pi)^-4 (tf tg)^-2 exp(-|za-zb|^2/(4(tf+tg))) exp(-alpha |y - y0|^2),  alpha = 1/(4tf)+1/(4tg),
#   y0 = (za/tf + zb/tg)/(1/tf + 1/tg),  int d^4y exp(-alpha |y-y0|^2) = (pi/alpha)^2,  sigma2 = 1/(2 alpha).
za_n = [sp.Rational(3, 10), sp.Rational(-1, 5), sp.Rational(1, 2), sp.Rational(1, 10)]
zb_n = [sp.Rational(-2, 5), sp.Rational(1, 4), sp.Rational(3, 20), sp.Rational(-3, 10)]
Phi_f = []
alpha = 1/(4*tf) + 1/(4*tg)
y0e = [(za_n[i]/tf + zb_n[i]/tg)/(1/tf + 1/tg) for i in range(4)]
sub = {sig2: 1/(2*alpha)}
sub.update({y0[i]: y0e[i] for i in range(4)})
sub.update({za_s[i]: za_n[i] for i in range(4)})
sub.update({zb_s[i]: zb_n[i] for i in range(4)})
dz2 = sum((za_n[i] - zb_n[i])**2 for i in range(4))
pref = (4*sp.pi)**-4*(tf*tg)**-2*sp.exp(-dz2/(4*(tf + tg)))*(sp.pi/alpha)**2
for k in range(4):
    Phi_f.append(sp.lambdify((tf, tg), pref*Phi_poly[k].subs(sub), 'mpmath'))
mp.mp.dps = 20
r3 = {}
for k in range(4):
    vals = {}
    for Tv in [1e2, 1e4, 1e6, 1e8]:
        e = mp.mpf(1)
        Tm = mp.mpf(Tv)
        g = lambda w: mp.e**w*Tm*(Phi_f[k](Tm, mp.e**w) + Phi_f[k](mp.e**w, Tm))
        val = mp.quad(g, [mp.log(e)] + [mp.log(mp.mpf(10)**j) for j in range(1, int(round(np.log10(Tv))))] + [mp.log(Tm)])
        vals['T=%g' % Tv] = float(val)
    r3['k=%d' % (k + 1)] = {'TdF_k': vals, 'limit_prediction_-3b_k/(16pi^2)': float(-3*bS[k]/(16*sp.pi**2)),
                            'rel_dev_at_T=1e8': abs(vals['T=1e+08'] + 3*float(bS[k])/(16*np.pi**2))/(3*abs(float(bS[k]))/(16*np.pi**2))}
res['3_z_nonzero_exact_moments_eps=1'] = {'z_a': [str(x) for x in za_n], 'z_b': [str(x) for x in zb_n],
                                          'b': [str(x) for x in bS], 'values': r3}
print('(3) done %.1fs' % (time.time() - t_start), flush=True)

# ------------------------------------------------------------------ (4) orientation and normalisation of CL
res['4_composite_line'] = {
    'CL_k(z_a,z_b)': '-2 (-,-)^{-1} F_k(z_a,z_b)  [2 (-,-)^{-1} from A(x)A, -(-,-)^{-1} d_x D from A^*|_V(x)c, (-,-) of S]',
    'lim_T TdCL/dT': '-2 x (-3/(16 pi^2)) b_k (-,-)^{-1} = +3/(8 pi^2) b_k (-,-)^{-1}',
    'coefficient_of_logT_in_one_loop_part': 'hbar d_{CL} W_0 -> (3/(8 pi^2)) hbar sum_a d_{A|_V=b X_a} d_{c=X^a} W_0 (left derivatives)',
    'engine_units_hbar/(16pi^2)': '16 pi^2 x 3/(8 pi^2) = 6'}
res['seconds'] = round(time.time() - t_start, 1)
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
