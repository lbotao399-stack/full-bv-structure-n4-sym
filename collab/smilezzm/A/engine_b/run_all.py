"""Engine B (IP1): tree and one-loop two-point functions of the scalar bilinears O_ab, per graph class, in the
heat-kernel UV-boundary scheme, by the POSITION-SPACE Schwinger-time representation of whole graphs.

Run:  python3 run_all.py      (writes results.json and gates.json next to this file; about 5 minutes)
Derivations and conventions: notes.txt.
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import sympy as sp
import mpmath as mp

from census import (census_summary, fermion_trace_check, ex_integrand_terms, se_bubble_terms, sun_f)
from gauss import gaussian, vec_deriv
from hepp import Term, asymptotics, finite_eps, fit_log
from oneloop import g1_table, g2_table, check_DK_product
from classes import se_exact, se_numeric, qu_exact, qu_numeric, oc_exact, b_exact

T0 = time.time()
LOG = []


def log(*a):
    s = ' '.join(str(x) for x in a)
    print('[%7.1fs] %s' % (time.time() - T0, s), flush=True)
    LOG.append(s)


def sstr(e):
    if isinstance(e, str):
        return e
    e = sp.sympify(e)
    return str(sp.expand(sp.expand_log(sp.nsimplify(e) if e.is_Rational or e.is_Float else e, force=True)))


def fval(e):
    return float(sp.N(e, 20))


def build_terms(pref, nx, n1, phi, ts, scale):
    """Split the Gaussian-integrated integrand into homogeneous nu-terms h_nu / phi^nu (|x|^2 = 1)."""
    n = len(ts)
    terms = []
    lam = sp.Symbol('lam', positive=True)
    for h in (pref * nx, pref * n1):
        h = sp.factor(sp.together(h * scale))
        if h == 0:
            continue
        r = sp.factor(h.subs({t: lam * t for t in ts}, simultaneous=True) / h)
        k = sp.Poly(sp.numer(r), lam).degree() - sp.Poly(sp.denom(r), lam).degree()
        nu = -k - n
        terms.append(Term(h / phi ** nu, nu, phi, ts))
    return terms


def fit_fixed_a(eps_arr, vals, a):
    """b, c, d of vals - a log eps = b + eps (c log eps + d)."""
    e = np.asarray(eps_arr)
    A = np.stack([np.ones_like(e), e * np.log(e), e], axis=1)
    return np.linalg.lstsq(A, np.asarray(vals) - a * np.log(e), rcond=None)[0]


results = {'engine': 'B',
           'method': ('position-space Schwinger-time representation of whole graphs: every edge '
                      'D_eps(r) = int_eps^inf (4 pi t)^-2 exp(-r^2/4t) dt (and derivatives), internal vertices '
                      'integrated exactly as Gaussians (sympy), overall Schwinger scale integrated exactly '
                      '(lower incomplete gamma), eps -> 0 asymptotics a log eps + b by Hepp sectors with an exact '
                      'subtraction of the logarithmic sectors (b as convergent integrals, Gauss-Legendre); '
                      'SE, QU, B, OC, T in closed form; finite-eps evaluations and least-squares extrapolation as '
                      'confirmation. No momentum-space bubble reduction, no dimensional regularisation, no IR '
                      'regulator.')}
gates = {}

# =====================================================================================  1. census
log('census: Wick enumeration with su(2), su(3)')
cs, patterns = census_summary(groups=(2, 3))
census_table = {}
for (N, cl, kinds), v in sorted(cs.items()):
    census_table.setdefault(cl + ' ' + '+'.join(kinds), {})['su(%d)' % N] = {
        'graphs': v['n'], 'graphs_with_zero_colour': v['zero_colour'],
        'S12': round(float(v['S12']), 12), 'A12': round(float(v['A12']), 12), 'T12': round(float(v['T12']), 12),
        'fit_residual': v['resid']}
log(json.dumps(census_table, indent=None)[:400], '...')
# G8: weights per (g^2 C_A dim g) identical for su(2), su(3); 'other' graphs all of zero colour
g8_ok = True
for key, d in census_table.items():
    a, b = d['su(2)'], d['su(3)']
    for s in ('S12', 'A12', 'T12'):
        if abs(a[s] - b[s]) > 1e-9:
            g8_ok = False
    if key.startswith('other') and (a['graphs_with_zero_colour'] != a['graphs'] or b['graphs_with_zero_colour'] != b['graphs']):
        g8_ok = False
max_A12 = max(abs(d[g]['A12']) for d in census_table.values() for g in d)
max_res = max(d[g]['fit_residual'] for d in census_table.values() for g in d)

err_gamma, sc = fermion_trace_check()
log('gamma-matrix trace check: max error', err_gamma, '; sum_m tr(g^m gb^I g^m gb^J) diag', sc[0, 0])

# EX integrand from the enumerated patterns, and comparison with the current form U.U
t5 = sp.symbols('t1:6', positive=True)
ex_terms, ex_acc = ex_integrand_terms(patterns[(3, 'EX')], t5)
ex_terms2, _ = ex_integrand_terms(patterns[(2, 'EX')], t5)
E_EX = [('x', 0, t5[0]), (0, 'o', t5[1]), ('x', 1, t5[2]), (1, 'o', t5[3]), (0, 1, t5[4])]
pref, nx, n1, phi = gaussian(2, E_EX, ex_terms)
d1, d2, d3, d4 = (vec_deriv(0, 'x', t5[0]), vec_deriv(0, 'o', t5[1]), vec_deriv(1, 'x', t5[2]), vec_deriv(1, 'o', t5[3]))
prefh, nxh, n1h, phih = gaussian(2, E_EX, [(1, [(d1, d3)]), (-1, [(d1, d4)]), (-1, [(d2, d3)]), (1, [(d2, d4)])])
X2 = sp.Symbol('X2', positive=True)
ex_hand_diff = sp.simplify(pref * (X2 * nx + n1) - prefh * (X2 * nxh + n1h))
log('EX integrand (enumeration) - hand form U.U:', ex_hand_diff)

# SE bubble: enumeration vs reduced kernel Sigma_b = g^2 C_A [-2 Lap(D^2) + D K], pointwise after s2<->s3
s4 = sp.symbols('s1:5', positive=True)
se_terms, se_acc = se_bubble_terms(patterns[(3, 'SE_bubble')], s4)
E_SE = [('x', 0, s4[0]), (0, 1, s4[1]), (0, 1, s4[2]), (1, 'o', s4[3])]
p_, nx_, n1_, ph_ = gaussian(2, E_SE, se_terms)
direct = p_ * (X2 * nx_ + n1_) * sp.exp(-X2 * ph_)
p0, nx0, n10, ph0 = gaussian(2, E_SE, [(1, [])])
gfun = p0 * n10 * sp.exp(-X2 * ph0)
red = 2 * (-2 * sp.diff(gfun, s4[0]) - sp.Rational(1, 2) * (sp.diff(gfun, s4[1]) + sp.diff(gfun, s4[2])))
sym = lambda e: (e + e.subs({s4[1]: s4[2], s4[2]: s4[1]}, simultaneous=True)) / 2
ibp_diff = sp.simplify(sym(direct) - sym(red))
log('SE bubble: enumeration - reduced kernel (pointwise, symmetrised):', ibp_diff)

results['census'] = {
    'table_weights_per_g2_CA_dimg': census_table,
    'note': ('weights include the vertex coefficients of -S_E, 1/n! and the Lorentz trace of seagull loops; '
             'S12/A12/T12 = least-squares flavour coefficients from 6 random flavour sets; EX and SE_bubble '
             'sum to 0 here only because their spacetime factors (derivative placements) differ per pattern; '
             'the pattern sums are used to build the EX and SE integrands.'),
    'EX_patterns_S12': {str(k): round(float(v), 12) for k, v in ex_acc.items()},
    'SE_bubble_patterns_S12': {str(k): round(float(v), 12) for k, v in se_acc.items()},
    'EX_integrand_minus_hand_form': str(ex_hand_diff),
    'SE_bubble_minus_reduced_kernel': str(ibp_diff),
    'gamma_trace_max_error': float(err_gamma),
    'max_A12_coefficient': max_A12, 'max_flavour_fit_residual': max_res,
    'vanishing_graphs': ('V3 V3 graphs with a cubic vertex carrying both legs of one operator, and V4 / seagull '
                         'contractions with two legs of one f^{abe} on the same delta: colour zero (counted '
                         'in graphs_with_zero_colour); Yukawa and ghost vertices cannot sit on scalar lines; '
                         '-I^CT cubic/quartic vertices and the 27/16 term enter only at higher order or not at all.')}

# =====================================================================================  2. G1, G2
log('one-loop singular parts (G1, G2)')
g1, g1aux = g1_table()
g2 = g2_table()
log('G1', {k: (str(v['X']), str(v['Y'])) for k, v in g1.items()})
log('G2', {k: {kk: str(vv) for kk, vv in v.items()} for k, v in g2.items()})
log('DK product identity residual:', check_DK_product())
g1_target = {'seagull_AV_loop': (2, 0), 'seagull_Aperp_loop': (sp.Rational(5, 2), 0),
             'bubble_AV_Aperp': (sp.Rational(-1, 4), 1), 'gaugino_loop': (-2, -2), 'total': (sp.Rational(9, 4), -1)}
g1_pass = all(sp.simplify(g1[k]['X'] - v[0]) == 0 and sp.simplify(g1[k]['Y'] - v[1]) == 0 for k, v in g1_target.items())
gates['G1'] = {'value': {k: {'eps^-1': sstr(g1[k]['X']), 'log_eps': sstr(g1[k]['Y'])} for k in g1},
               'target': {k: {'eps^-1': sstr(v[0]), 'log_eps': sstr(v[1])} for k, v in g1_target.items()},
               'pass': bool(g1_pass)}
g2_pass = (g2['e']['alpha'] == 2 and g2['e']['beta'] == -2 and g2['f']['alpha'] == 2 and g2['f']['beta'] == 0
           and g2['total']['alpha'] == 4 and g2['total']['beta'] == -2 and g2['total']['coef_ab_AVAV'] == 0
           and g2['a']['vacuum_coefficient_of_hbar_dimg_ab_over_16pi2_eps'] == 2)
gates['G2'] = {'value': {'total': 'unit log eps [%s O_ab %+d (a,b) K] + 2 hbar dim g (a,b)/(16 pi^2 eps) int omega'
                         % (g2['total']['alpha'], int(g2['total']['beta'])),
                         'bubble_e': [str(g2['e']['alpha']), str(g2['e']['beta'])],
                         'triangle_f': [str(g2['f']['alpha']), str(g2['f']['beta'])],
                         'AV_graphs_c_d': [str(g2['c']['coef_ab_AVAV']), str(g2['d']['coef_ab_AVAV'])],
                         'vacuum_a': str(g2['a']['vacuum_coefficient_of_hbar_dimg_ab_over_16pi2_eps'])},
               'target': {'total': 'unit log eps [4 O_ab - 2 (a,b) K] + 2 hbar dim g (a,b)/(16 pi^2 eps) int omega',
                          'bubble_e': ['2', '-2'], 'triangle_f': ['2', '0'], 'AV_graphs_c_d': ['-2', '2'],
                          'vacuum_a': '2'},
               'pass': bool(g2_pass)}
results['oneloop'] = {'G1_aux': {k: str(v) for k, v in g1aux.items()},
                      'G2_triangle_logs': {'f': str(g2['f']['triangle_log']), 'd': str(g2['d']['triangle_log'])}}

# =====================================================================================  3. classes
classes = {}
lnx = 'log(x^2/ell^2)'
# ---- T
classes['T'] = dict(c_S='0', c_T='0', d_S='0', d_T='0', log_eps_S='0', log_eps_T='0',
                    exact={'value': 'D_eps(x)^2/D_0(x)^2 = (1 - exp(-x^2/4 eps))^2 -> 1 (exponentially)',
                           'note': 'tree normalisation 4 hbar dim(g) S12 D_0(x)^2 reproduced: census weight 1 per dim g S12'})

# ---- SE
log('SE exact')
se = se_exact()
se_num = se_numeric([2.0 ** -k for k in (6, 8, 10, 12, 14, 16, 18)])
eps_se = np.array([r['eps'] for r in se_num])
IJ = np.array([r['IJ'] for r in se_num])
cfit = fit_log(eps_se[1:], IJ[1:], 4)
kappa_free = cfit[1]
# log coefficient fixed to its exact value -1; the corrections of I_J are analytic in eps (notes S.3)
_A = np.stack([np.ones_like(eps_se), eps_se, eps_se ** 2, eps_se ** 3], axis=1)
kappa_num = np.linalg.lstsq(_A[1:], IJ[1:] + np.log(eps_se[1:]), rcond=None)[0][0]
log('SE: c_S =', se['c_S'], '; kappa numeric fit', kappa_num, 'exact', fval(se['kappa']))
classes['SE'] = dict(c_S=sstr(se['c_S']), c_T='0', d_S='8', d_T='0', log_eps_S='0', log_eps_T='0',
                     exact={'c_S': sstr(se['c_S']), 'c_S_float': fval(se['c_S']),
                            'decomposition': {
                                'seagulls_minus_their_mass_counterterms': '0 exactly (D_eps(0) = 1/(16 pi^2 eps) exact)',
                                'Laplacian(D^2)_part_with_log_counterterm (both lines)': '8 kappa + 8 log(x^2), kappa = %s' % sstr(se['kappa']),
                                'D K - mass part (both lines)': '27/2 (= -2 * 2 * 9 * (-3/8) in units; exponentially small corrections)',
                                'log_eps_of_graphs_without_counterterm': '-8 S12', 'log_eps_of_counterterm_vertex': '+8 S12'},
                            'kappa': sstr(se['kappa'])},
                     numeric={'IJ_finite_eps': se_num, 'kappa_from_fit_log_coefficient_fixed': kappa_num,
                              'kappa_from_free_fit': kappa_free,
                              'kappa_fit_note': 'free fit includes the log coefficient as a parameter (ill-conditioned); the fixed fit uses the exact log coefficient -1 and analytic corrections',
                              'kappa_exact': fval(se['kappa']),
                              'method': 'mpmath quad of the exact one-dimensional Schwinger representations at finite eps, least-squares fit -log eps + kappa + eps(...)'})

# ---- QU
log('QU exact')
qu = qu_exact()
qu_num = qu_numeric([2.0 ** -k for k in (8, 10, 12, 14, 16, 18, 20, 22)])
eps_q = np.array([r['eps'] for r in qu_num])
I2 = np.array([r['I2'] for r in qu_num])
cq = fit_log(eps_q, I2, 6)
ex_I2 = 32 * (float(mp.euler) - 3 * math.log(2))
log('QU: I2 fit a, b =', cq[:2], 'exact', -32, ex_I2)
# Hepp-sector route for QU (whole graph, 4 edges)
t4 = sp.symbols('t1:5', positive=True)
pq, nxq, n1q, phq = gaussian(1, [('x', 0, t4[0]), ('x', 0, t4[1]), (0, 'o', t4[2]), (0, 'o', t4[3])], [(1, [])])
qu_terms = build_terms(pq, nxq, n1q, phq, t4, (4 * sp.pi) ** 6)
hq = {}
for N in (32, 48):
    a, br, bg, info = asymptotics(qu_terms, 4, N=N, power=5)
    hq[N] = (a, br + bg * float(mp.euler), br, bg)
log('QU Hepp:', hq, 'exact b', ex_I2)
eps_fe = 2.0 ** -np.arange(24, 45, 2)
fe_q = finite_eps(qu_terms, 4, eps_fe, plain_N=24, power=5)
cqf = fit_log(eps_fe, np.array(fe_q), 4)
cqf_fixed = fit_fixed_a(eps_fe, fe_q, -32.0)
log('QU finite-eps Hepp fit (free a):', cqf[:2], ' (a fixed to -32): b =', cqf_fixed[0])
classes['QU'] = dict(c_S=sstr(qu['c_S']), c_T=sstr(qu['c_T']), d_S=sstr(qu['d_S']), d_T=sstr(qu['d_T']),
                     log_eps_S=sstr(qu['log_eps_S']), log_eps_T=sstr(qu['log_eps_T']),
                     exact={'two_point_weight': 'g^2 C_A dim g (S12 - 2 T12) int d^4z D_eps(x-z)^2 D_eps(z)^2 (census)',
                            'I2': sstr(qu['I2']), 'c_S_float': fval(qu['c_S']), 'c_T_float': fval(qu['c_T'])},
                     numeric={'I2_exact_1d_representation': qu_num,
                              'fit_a_b': [float(cq[0]), float(cq[1])], 'exact_a_b': [-32.0, ex_I2],
                              'hepp_b_formula': {str(N): {'a': v[0], 'b': v[1]} for N, v in hq.items()},
                              'hepp_finite_eps_fit_a_b': [float(cqf[0]), float(cqf[1])],
                              'hepp_finite_eps_fit_b_with_a_fixed': float(cqf_fixed[0]),
                              'method': 'three routes: closed form; finite-eps one-dimensional exact representation + fit; Hepp sectors of the whole 4-edge graph (b-formula and finite eps)'})

# ---- EX
log('EX: Hepp sectors of the whole 5-edge graph')
ex_T = build_terms(pref, nx, n1, phi, t5, 4096 * sp.pi ** 6)
hx = {}
for N in (32, 40):
    a, br, bg, info = asymptotics(ex_T, 5, N=N, power=5)
    hx[N] = dict(a=a, b_rat=br, b_gamma=bg, info={k: (v if not isinstance(v, set) else sorted(v)) for k, v in info.items()})
    log('EX N=%d: a=%.15f b_rat=%.15f b_gamma=%.15f' % (N, a, br, bg), info)
b_rat = hx[40]['b_rat']
b_err = abs(hx[40]['b_rat'] - hx[32]['b_rat'])
a_ex = hx[40]['a']
mp.mp.dps = 30
rel = mp.pslq([mp.mpf(b_rat), 1, mp.log(2), mp.log(3)], tol=1e-9, maxcoeff=10 ** 4, maxsteps=10 ** 6)
log('PSLQ on b_rat with {1, log 2, log 3}:', rel)
pslq_scan = {}
for tol in (1e-9, 1e-10, 1e-11):
    for mc in (10 ** 2, 10 ** 3, 10 ** 4):
        r_ = mp.pslq([mp.mpf(b_rat), 1, mp.log(2), mp.log(3)], tol=tol, maxcoeff=mc, maxsteps=10 ** 6)
        r2_ = mp.pslq([mp.mpf(b_rat) + 32 * mp.euler, 1, mp.euler, mp.log(2), mp.log(3)], tol=tol, maxcoeff=mc, maxsteps=10 ** 6)
        pslq_scan['tol=%g,maxcoeff=%d' % (tol, mc)] = {'b_rational_part': [int(x) for x in r_] if r_ else None,
                                                       'b_with_gamma_free': [int(x) for x in r2_] if r2_ else None}
log('PSLQ scan:', pslq_scan)
rel_full = mp.pslq([mp.mpf(b_rat + hx[40]['b_gamma'] * float(mp.euler)), 1, mp.euler, mp.log(2), mp.log(3)],
                   tol=1e-10, maxcoeff=100, maxsteps=10 ** 6)
log('PSLQ on b with {1, gamma_E, log 2, log 3}:', rel_full)
if rel is not None and rel[0] != 0:
    b_rat_exact = -(sp.Integer(rel[1]) + sp.Integer(rel[2]) * sp.log(2) + sp.Integer(rel[3]) * sp.log(3)) / sp.Integer(rel[0])
else:
    b_rat_exact = None
b_gamma_exact = sp.Integer(int(round(-a_ex)))
# finite-eps confirmation (spec route): two resolutions
fe_ex = {}
for N in (16, 20):
    v = finite_eps(ex_T, 5, eps_fe, plain_N=N, power=5)
    c4 = fit_log(eps_fe, np.array(v), 4)
    c2 = fit_log(eps_fe[4:], np.array(v)[4:], 2)
    cf = fit_fixed_a(eps_fe, v, -32.0)
    fe_ex[N] = dict(values=[float(x) for x in v], fit4_a_b=[float(c4[0]), float(c4[1])],
                    fit2_small_eps_a_b=[float(c2[0]), float(c2[1])], fit_b_with_a_fixed_to_minus32=float(cf[0]))
    log('EX finite-eps N=%d: fit4 a,b =' % N, c4[:2], 'fit2 =', c2[:2], 'b (a fixed) =', cf[0],
        '; b (asymptotic formula) =', b_rat + hx[40]['b_gamma'] * float(mp.euler))
# relative two-point contribution: relative = lambda_hat S12 Etilde/16 = -unit S12 Etilde/8
EG = sp.EulerGamma
if b_rat_exact is not None:
    Etilde0 = b_rat_exact + b_gamma_exact * EG
    cS_EX = sp.nsimplify(sp.expand(-Etilde0 / 8))
else:
    cS_EX = None
cS_EX_num = -(b_rat + hx[40]['b_gamma'] * float(mp.euler)) / 8
classes['EX'] = dict(c_S=sstr(cS_EX) if cS_EX is not None else repr(cS_EX_num), c_T='0',
                     d_S=sstr(sp.Integer(int(round(a_ex))) / 8), d_T='0',
                     log_eps_S=sstr(-sp.Integer(int(round(a_ex))) / 8), log_eps_T='0',
                     numeric={'value': cS_EX_num,
                              'digits': int(-math.log10(max(b_err / 8, 1e-16)) - 1),
                              'Etilde_a': a_ex, 'Etilde_b_rational_part': b_rat, 'Etilde_b_gamma_coefficient': hx[40]['b_gamma'],
                              'resolution_N32_vs_N40_difference_of_b': b_err,
                              'hepp': {str(k): v for k, v in hx.items()},
                              'finite_eps_extrapolation': {str(k): v for k, v in fe_ex.items()},
                              'eps_values': [float(e) for e in eps_fe],
                              'pslq_b_rational_part_basis_1_log2_log3': [int(x) for x in rel] if rel else None,
                              'pslq_full_b_basis_1_gamma_log2_log3': [int(x) for x in rel_full] if rel_full else None,
                              'pslq_scan': pslq_scan,
                              'pslq_note': ('the 4-element relation for the rational part (EulerGamma coefficient fixed '
                                            'exactly to 32 = -a) is stable over tol 1e-9..1e-11 and maxcoeff up to 1e4; '
                                            'with EulerGamma free (5 elements) the same relation appears for maxcoeff 100 and '
                                            'tol <= 1e-10, while larger maxcoeff admits spurious relations at this precision'),
                              'residual_of_identified_b_rational_part': float(b_rat - fval(64 - 64 * sp.log(2) - 32 * sp.log(3))),
                              'identification': ('c_S(EX) = %s (NUMERICAL identification by PSLQ; the EulerGamma '
                                                 'coefficient -4 and the log eps coefficient 4 are exact '
                                                 'consequences of a = -32)' % sstr(cS_EX)),
                              'method': ('Etilde = 4096 pi^6 int_[eps,inf)^5 dt U^-2 [P_a P_b/(4U^2) + 2(t1-t2)(t3-t4)/U] '
                                         'exp(-W/(4U)) (Gaussian over z1, z2 from the enumerated vertices); '
                                         'Hepp sectors (120), b as convergent integrals (Gauss-Legendre, y=u^5, '
                                         'N=32 and N=40); relative contribution = -unit S12 Etilde/8')})

# ---- OC, B
oc = oc_exact(g2)
classes['OC'] = dict(c_S='0', c_T='0', d_S='0', d_T='0', log_eps_S=sstr(oc['log_eps_S']), log_eps_T=sstr(oc['log_eps_T']),
                     exact={'value': '-unit log eps [2 alpha S12 + 4 beta T12] with (alpha, beta) = (4, -2) from G2; '
                                     'constant (vacuum) term of O(eps) contributes 0 to the connected function'})
bb = b_exact()
classes['B'] = dict(c_S=sstr(bb['c_S']), c_T='0', d_S='0', d_T='0', log_eps_S='0', log_eps_T='0',
                    exact={'value': '-F_0 = 4 beta_2 unit x (free scalar action) on each line: D_eps -> D_eps - 4 beta_2 unit D_2eps'})
classes['F'] = dict(c_S='0', c_T='0', d_S='0', d_T='0', log_eps_S='0', log_eps_T='0',
                    exact={'value': 'minimal scheme: finite parts psi(O_a) = 0'})

# ---- totals
b2 = sp.Symbol('beta_2')
tot = {}
for key in ('c_S', 'c_T', 'd_S', 'd_T', 'log_eps_S', 'log_eps_T'):
    tot[key] = sp.nsimplify(sp.expand(sum(sp.sympify(classes[c][key], locals={'beta_2': b2}) for c in classes)))
tot['c_S'] = sp.expand(sp.expand_log(tot['c_S'], force=True))
log('TOTAL:', {k: str(v) for k, v in tot.items()})
target = {'c_S': sp.Rational(27, 2) + 4 * sp.log(sp.Rational(3, 2)) - 8 * b2, 'c_T': 8 * (EG - 3 * sp.log(2)),
          'd_T': sp.Integer(8)}
disc = {k: sp.simplify(sp.expand_log(tot[k] - target[k], force=True)) for k in target}
log('difference to spec targets:', disc)
results['classes'] = classes
results['total'] = {k: sstr(v) for k, v in tot.items()}
results['total']['c_S_float_beta2_0'] = fval(tot['c_S'].subs(b2, 0))
results['total']['c_T_float'] = fval(tot['c_T'])
results['total']['formula'] = ('<O1,O2>(x) = 4 hbar dim(g) D_0(x)^2 [ S12 (1 + unit c_S) + T12 unit (c_T + d_T log(x^2/ell^2)) ], '
                               'd_S = 0, log eps cancelled')
results['total']['difference_to_spec_targets'] = {k: str(v) for k, v in disc.items()}
results['other_structures'] = ('0: the antisymmetric (a1,a2)(b1,b2)-(a1,b2)(b1,a2) has coefficient <= %.1e in every class '
                               '(census flavour fits, residual %.1e); only S12 and T12 occur' % (max_A12, max_res))
results['beta_2_dependence'] = 'c_S contains -8 beta_2 (class B only); c_T, d_T independent of beta_2'

# =====================================================================================  4. gates
gates['G3'] = {'value': {'log_eps_S': results['total']['log_eps_S'], 'log_eps_T': results['total']['log_eps_T']},
               'target': {'log_eps_S': '0', 'log_eps_T': '0'},
               'pass': bool(tot['log_eps_S'] == 0 and tot['log_eps_T'] == 0)}
# G4: total and per class: log x^2 coefficient = -(1/4)(tree insertions of D_1 piece at both ends)
#  D_1 pieces: graph (f) 4*alpha_f O, graph (e) 4(alpha_e O + beta_e (a,b) K), field counting -4*(-4 Y) O (Y from G1)
Ytot = g1['total']['Y']
fc = -4 * Ytot
pred = {'EX': (-2 * g2['f']['alpha'], -4 * g2['f']['beta']),
        'QU': (-2 * g2['e']['alpha'], -4 * g2['e']['beta']),
        'SE': (-2 * (-fc), 0)}
per_class = {}
ok4 = True
for c, (pS, pT) in pred.items():
    vS, vT = sp.sympify(classes[c]['d_S']), sp.sympify(classes[c]['d_T'])
    per_class[c] = {'d_S': str(vS), 'd_T': str(vT), 'predicted_from_D1_piece': [str(pS), str(pT)]}
    ok4 = ok4 and (sp.simplify(vS - pS) == 0) and (sp.simplify(vT - pT) == 0)
D1K = 4 * (g2['total']['alpha'] + 6 * g2['total']['beta'] - fc)
gates['G4'] = {'value': {'d_T': str(tot['d_T']), 'd_S': str(tot['d_S']),
                         'log_x2_coefficient_on_K_per_tree': str(sp.Rational(tot['d_S'] * 12 + tot['d_T'] * 36, 12)),
                         'D_1(K)': '%s unit' % D1K, 'D_1(20prime)': '%s unit' % (4 * (g2['total']['alpha'] - fc)),
                         'per_class': per_class},
               'target': {'d_T': '8', 'd_S': '0', 'log_x2_coefficient_on_K_per_tree': '24', 'D_1(K)': '-48 unit',
                          'D_1(20prime)': '0 unit'},
               'pass': bool(tot['d_T'] == 8 and tot['d_S'] == 0 and D1K == -48 and ok4
                            and 4 * (g2['total']['alpha'] - fc) == 0)}
gamma_cS = sp.expand(tot['c_S']).coeff(EG)
gates['G5'] = {'value': str(gamma_cS), 'target': '0', 'pass': bool(gamma_cS == 0),
               'note': 'EulerGamma coefficients per class: SE +8, QU -4, EX -4 (EX: exact, = -(log eps coefficient))'}
gates['G6'] = {'value': {'max_A12': max_A12, 'max_flavour_fit_residual': max_res}, 'target': 0,
               'pass': bool(max_A12 < 1e-9 and max_res < 1e-9)}


# G7: Gram matrices on the 21-dimensional family with a complex basis f_1..f_6
def gram():
    e = np.eye(6)
    f = []
    for k in range(3):
        f.append((e[2 * k] + 1j * e[2 * k + 1]) / np.sqrt(2))
        f.append((e[2 * k] - 1j * e[2 * k + 1]) / np.sqrt(2))
    ops = [(f[i], f[j]) for i in range(6) for j in range(i, 6)]
    cS = fval(tot['c_S'].subs(b2, 0))
    cT = fval(tot['c_T'])
    n = len(ops)
    Gt = np.zeros((n, n), complex)
    G1 = np.zeros((n, n), complex)
    Gb = np.zeros((n, n), complex)
    for k, (a1, b1) in enumerate(ops):
        for l, (a2, b2_) in enumerate(ops):
            A1, B1 = a1.conj(), b1.conj()
            S = np.dot(A1, a2) * np.dot(B1, b2_) + np.dot(A1, b2_) * np.dot(B1, a2)
            T = np.dot(A1, B1) * np.dot(a2, b2_)
            Gt[k, l] = S
            G1[k, l] = cS * S + cT * T   # at |x| = ell, log term 0
            Gb[k, l] = -8 * S
    herm = max(np.max(np.abs(Gt - Gt.conj().T)), np.max(np.abs(G1 - G1.conj().T)), np.max(np.abs(Gb - Gb.conj().T)))
    ev = np.linalg.eigvalsh(Gt)
    # K = sum_I O_{e_I e_I} = 2 (O_{f1 f2} + O_{f3 f4} + O_{f5 f6}) in this basis
    kvec = np.zeros(n, complex)
    for k, (a1, b1) in enumerate(ops):
        for p in range(3):
            if (a1 is f[2 * p] and b1 is f[2 * p + 1]):
                kvec[k] = 2.0
    # traceless subspace: orthogonal complement of the trace functional  X -> sum_I X_II
    trace_fun = np.array([np.dot(a1, b1) for (a1, b1) in ops])  # (a,b) of the operator
    M = np.array([np.eye(n)[k] for k in range(n)], complex)
    # basis of the 20': vectors v with sum_k v_k (a_k, b_k) = 0
    u, s, vh = np.linalg.svd(trace_fun.reshape(1, -1))
    tl = vh[1:].conj().T  # n x 20
    orth_t = np.max(np.abs(kvec.conj() @ Gt @ tl))
    orth_1 = np.max(np.abs(kvec.conj() @ G1 @ tl))
    orth_b = np.max(np.abs(kvec.conj() @ Gb @ tl))
    return dict(hermiticity_defect=float(herm), tree_eigenvalues_min=float(ev.min()), tree_eigenvalues_max=float(ev.max()),
                K_20prime_orthogonality_tree=float(orth_t), K_20prime_orthogonality_oneloop=float(orth_1),
                K_20prime_orthogonality_beta2=float(orth_b), dimension=n, rank_traceless=tl.shape[1])


gr = gram()
log('G7 Gram:', gr)
gates['G7'] = {'value': gr, 'target': 'Hermitian, Gram_tree > 0, (20prime, K) orthogonal',
               'pass': bool(gr['hermiticity_defect'] < 1e-12 and gr['tree_eigenvalues_min'] > 0
                            and max(gr['K_20prime_orthogonality_tree'], gr['K_20prime_orthogonality_oneloop'],
                                    gr['K_20prime_orthogonality_beta2']) < 1e-12)}
gates['G8'] = {'value': {'weights_su2_equal_su3': g8_ok,
                         'colour_factors': 'tree dim g; every one-loop graph g^2 C_A dim g x rational (no 1/N); '
                                           'fermion loop and SE kernels C_A by the Lie-algebra identities f^{acd}f^{bcd}=N delta'},
               'target': 'colour factor dim(g) (tree) or dim(g) tr(ad ad)/(X,Y) (one loop), identical rationals for su(2), su(3)',
               'pass': bool(g8_ok)}

results['conventions'] = {
    'action': 'S/hbar = -S_E, canonical su(N) fields A = g A^a T_a, (X,Y) = (2 hbar/g^2) X^a Y^a (spec normalisation); '
              'S_E = int 1/4 F^2 + 1/2 (D phi)^2 + 1/2 (d.A)^2 + g^2/4 f f phi^4 + 1/2 lambda gammabar^M D_M lambda (10d MW, Pfaffian)',
    'propagator': 'every edge D_eps(r) = (1 - exp(-r^2/4eps))/(4 pi^2 r^2) = int_eps^inf (4 pi t)^-2 exp(-r^2/4t) dt, delta^ab delta_MN; '
                  'fermion -gamma^m d_m D_eps',
    'vertices_of_-S_E': 'V3 = -g f^{abc} d_m phi^a A^b_m phi^c;  seagull -(g^2/2) f f A A phi phi;  V4 = -(g^2/4) f^{abe} f^{cde} phi^a_I phi^b_J phi^c_I phi^d_J;  '
                        'Yukawa from (1/2) lambda gammabar^I g f phi_I lambda; fermion loop -(1/2) Tr log -> +(1/4) Tr(G V G V)',
    'operators': 'O_ab = a_I b_J phi^c_I phi^c_J (= manuscript O_ab / (2 hbar)); normalised to the tree 4 hbar dim g S12 D_0^2',
    'unit': 'unit = -g^2 N/(32 pi^2) = -lambda_hat/2, lambda_hat = g^2 N/(16 pi^2); C_A = N',
    'symmetry_factors': 'automatic from the Wick enumeration (1/n! for identical vertices, all contractions)',
    'singular_part': 'terms eps^-k and log eps of the plane-wave value (log eps = log(eps/ell^2), ell = 1), no constant subtracted',
    'counterterm_vertices': 'J = I - F_0 - I^CT[eps], O(eps) = O - d_delta I^CT(O): -I^CT quadratic vertex on each regulated scalar line '
                            '(both lines), -d_delta I^CT(O_a) as a tree vertex at x resp. 0 joined by two regulated lines; '
                            '-F_0 = -beta_2 unit int(dA_perp,dA_perp) on each line',
    'reporting': '|x|^2 = ell^2 = 1 by scaling for the constants; log(x^2/ell^2) from the eps/x^2 dependence of each graph',
}
results['log'] = LOG
results['runtime_seconds'] = time.time() - T0

with open(os.path.join(HERE, 'results.json'), 'w') as fh:
    json.dump(results, fh, indent=1, default=str)
with open(os.path.join(HERE, 'gates.json'), 'w') as fh:
    json.dump(gates, fh, indent=1, default=str)
log('gates:', {k: v['pass'] for k, v in gates.items()})
log('done')
