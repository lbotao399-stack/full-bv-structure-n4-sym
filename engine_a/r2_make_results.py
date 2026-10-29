"""Assemble engine_a/results_r2.json (round 2) from the exact runs in logs/."""
import json, itertools
import numpy as np
import sympy as sp
import r2_core as C

L = lambda f: json.load(open('logs/' + f))
run7 = {}
for f in ('r2_T1_test.json', 'r2_T1_cubic.json', 'r2_T1_quartic.json'):
    run7.update(L(f))
run13 = L('r2_T1_seed13.json')
van = L('r2_T1_vanishing.json')
T2a = L('r2_T2_routeA.json'); T2b = L('r2_T2b.json'); W = L('r2_T3_weights.json')
regress = L('r2_regress_quadratic.json')

# ---------------- su(2) adjoint trace identity used to convert quartic colour structures (checked here)
F = C.F.real
ad = np.array([F[a].T for a in range(3)])          # (ad_a)_{cb} = F[a,b,c]
T4 = np.einsum('aij,bjk,ckl,dli->abcd', ad, ad, ad, ad)
dd = np.einsum('ab,cd->abcd', np.eye(3), np.eye(3))
assert np.allclose(T4, dd + np.einsum('ad,bc->abcd', np.eye(3), np.eye(3)))

SET_TXT = {
 'AV,AV,AV': '(A|_V, A|_V, A|_V)', 'AV,AP,AP': '(A|_V, A|_Wperp, A|_Wperp)', 'l,l,AV': '(lambda, lambda, A|_V)',
 'l,l,AP': '(lambda, lambda, A|_Wperp)', 'c,AV,AsV': '(c, A|_V, A^*|_V)', 'c,AP,AsP': '(c, A|_Wperp, A^*|_Wperp)',
 'c,l,ls': '(c, lambda, lambda^*)', 'c,c,cs': '(c, c, c^*)', 'AV,AV,AV,AV': '(A|_V)^4',
 'AV,AV,AP,AP': '(A|_V)^2 (A|_Wperp)^2', 'AP,AP,AP,AP': '(A|_Wperp)^4'}
TREE_TXT = {
 'AV,AV,AV': 'int dvol -1/2 (dA, 1/2[A ^ A]) restricted to (A|_V)^3',
 'AV,AP,AP': 'int dvol -1/2 (dA, 1/2[A ^ A]) restricted to A|_V (A|_Wperp)^2',
 'l,l,AV': 'int dvol 1/2 (lambda, rho([A|_V, lambda]))', 'l,l,AP': 'int dvol 1/2 (lambda, rho([A|_Wperp, lambda]))',
 'c,AV,AsV': 'int dvol ([A|_V, c], A^*|_V)', 'c,AP,AsP': 'int dvol ([A|_Wperp, c], A^*|_Wperp)',
 'c,l,ls': 'int dvol ([c, lambda], lambda^*)', 'c,c,cs': 'int dvol 1/2 ([c, c], c^*)',
 'AV,AV,AV,AV': 'int dvol -1/4 (1/2[A ^ A], 1/2[A ^ A]) restricted to (A|_V)^4',
 'AV,AV,AP,AP': 'int dvol -1/4 (1/2[A ^ A], 1/2[A ^ A]) restricted to (A|_V)^2 (A|_Wperp)^2',
 'AP,AP,AP,AP': 'int dvol -1/4 (1/2[A ^ A], 1/2[A ^ A]) restricted to (A|_Wperp)^4'}
STRUCT_DEF = {
 'cubic A (omega=1)': 'legs i=1,2,3 with polarisations a_i (element of V_10^vee (x) g) and incoming momenta p_i (p_3=-p_1-p_2); '
     '(a_i.a_j)(p.a_k) := ([a_i-colour, a_j-colour]-free contraction) = ((a_k(p)), [.,.]) pattern: the structure "(a1.a2)(p1.a3)" is '
     '(a_3(p_1), [a_1, a_2]) with a_1,a_2 contracted by the metric of the stated subspace and the colour parts combined by (X,[Y,Z]); '
     '"eps(a1,a2,a3,p)" uses the volume form of V.',
 'lambda lambda A (omega=0)': '"(u1,rho(a)[.,u2])" = (u_1, rho([a, u_2])) (ESW pairing Sigma x Sigma^* and Clifford action); '
     '"(u1,rho(a) chirality [.,u2])" = the same with the 4d chirality operator rho(q(dvol)) inserted before u_2.',
 'ghost-sector (omega=0)': '"([gamma,a],alpha)" = ([gamma, a], alpha) with V_10 metric; "([gamma,u],v)" = ([gamma,u],v) with the Sigma x Sigma^* pairing; '
     '"([g1,g2],w)" = ([gamma_1,gamma_2], w).',
 'quartic (omega=0)': 'Colour pairings X=(12)(34), Y=(13)(24), Z=(14)(23) of the four g-factors through (-,-) (g = su(2) in the code), '
     'Lorentz pairings x=(12)(34), y=(13)(24), z=(14)(23) of the V_10^vee factors through the metric of V or W_R^perp; "Xy" = product. '
     'Algebra-independent form: with tr_g(ad ad ad ad)[ijkl] the trace of the product of ad of the colour parts of legs i,j,k,l in this order, '
     'the su(2) identity tr_g(ad_a ad_b ad_c ad_d) = (a,b)(c,d) + (a,d)(b,c) (checked in r2_make_results.py) gives '
     'X = (t[1234] + t[1243] - t[1324])/2, Y = (t[1243] + t[1324] - t[1234])/2, Z = (t[1234] + t[1324] - t[1243])/2; the coefficients '
     'of the three orderings t[...] are algebra-independent (every one-loop colour factor is a trace of products of ad).'}
WEIGHT = {'A': 2, 'l': 4, 'c': 1, 'As': -2, 'ls': -4, 'cs': -1}
LEGT = {'AV,AV,AV': 'A A A', 'AV,AP,AP': 'A A A', 'l,l,AV': 'l l A', 'l,l,AP': 'l l A', 'c,AV,AsV': 'c A As', 'c,AP,AsP': 'c A As',
        'c,l,ls': 'c l ls', 'c,c,cs': 'c c cs', 'AV,AV,AV,AV': 'A A A A', 'AV,AV,AP,AP': 'A A A A', 'AP,AP,AP,AP': 'A A A A'}

def ratio_to_tree(logd, treed):
    keys = set(logd) | set(treed)
    r = None
    for k in keys:
        x = sp.sympify(logd.get(k, '0')); t = sp.sympify(treed.get(k, '0'))
        if t == 0:
            if x != 0: return None
            continue
        q = sp.nsimplify(x / t)
        if r is None: r = q
        elif sp.simplify(q - r) != 0: return None
    return r if r is not None else sp.Integer(0)

def trace_basis(logd):
    out = {}
    lor = sorted(set(k[1:] for k in logd if k[0] in 'XYZ'))
    for l in lor:
        cX, cY, cZ = [sp.sympify(logd.get(c + l, '0')) for c in 'XYZ']
        for nm, v in (('t[1234]', (cX - cY + cZ) / 2), ('t[1243]', (cX + cY - cZ) / 2), ('t[1324]', (-cX + cY + cZ) / 2)):
            if v != 0:
                out[nm + ' x ' + l] = str(sp.nsimplify(v))
    return out

def chain_T1(setname, g, kappa):
    n = {'tadpole': 1, 'bubble': 2, 'triangle': 3, 'box': 4}[g['topology']]
    omega = 1 if setname in ('AV,AV,AV', 'AV,AP,AP') else 0
    m = n - 2
    import math
    return [
     'amplitude: (1/(2n)) sum over orderings of the external insertions of str[ P(k_1)(-1)^{|.|} d^2 I_{v_1} P(k_2)(-1)^{|.|} d^2 I_{v_2} ... ], '
     'n=%d, vertices %s, internal lines %s; P(k) = e^{-eps k^2} k^{-2} x {2 (-,-)^{-1} (x) metric^{-1} (A-A), -i k (-,-)^{-1} (c-A^*|_V), '
     'i rho(k) (-,-)^{-1} (lambda-lambda)} = int_eps^infty dt e^{-t k^2} x numerator; line e carries k + q_e, q_e = -(sum of external momenta entering before it)'
     % (n, ', '.join(g['vertices']), ', '.join(g['internal_lines'])),
     'numerator algebra: N(k) = str[...] (exact Gaussian-integer arithmetic, 10d chiral Clifford blocks and su(2) colour only in code), '
     'polynomial of degree %d in (k, external momenta); only (4^m m!)^{-1} (d_k.d_k)^m N with m = %d enters the log eps coefficient'
     % (2 * n - 4 + omega, m),
     'Schwinger integral: int_{t_e in [eps,infty)} prod_e dt_e (16 pi^2)^{-1} (sum_e t_e)^{-2-m} [(4^m m!)^{-1}(d_k.d_k)^m N](k = -sum_e t_e q_e / sum_e t_e) '
     'x exp(-sum_e t_e q_e^2 + (sum_e t_e q_e)^2 / sum_e t_e)',
     'eps-expansion: {t_e >= eps for all e} = {sum_e t_e >= eps / min_e(t_e/sum t)}; the scale integral gives -log eps + O(1), the remaining '
     'integral over the fractions t_e/sum t (simplex, volume 1/(n-1)! = 1/%d) is %s; Gaussian moments give integer powers of sum t, the degree '
     'count leaves no eps^{-1} and no eps^{-1/2}' % (math.factorial(n - 1), 'that of a linear function = value at t_e/sum t = 1/n' if omega == 1 else 'that of a constant'),
     'local part: %s' % (('(%s) u log eps x [%s]' % (kappa, TREE_TXT[setname])) if kappa is not None else
                        'u log eps x (structures %s)' % json.dumps(g['log_eps']))]

# ---------------- T1
T1 = {'definitions': STRUCT_DEF, 'units': 'u = hbar tr_g(ad_X ad_Y)/((X,Y) 16 pi^2); computed with g = su(2), orthonormal for (-,-), '
      'for which tr_g(ad_X ad_Y)/(X,Y) = -2; every number is an exact rational from exact arithmetic (Gaussian integers, Fractions, '
      'exact Gauss-Jordan over Q(i) on complete invariant bases with zero residual on extra samples); two independent random seeds agree.',
      'sets': {}, 'formula_chains_note': 'per-graph chains in formula_chains.T1'}
chains = {'T1': {}, 'T2': {}}
ws = sp.symbols('beta wV wP wl wc')
eqs = []
for key in ['AV,AV,AV', 'AV,AP,AP', 'l,l,AV', 'l,l,AP', 'c,AV,AsV', 'c,AP,AsP', 'c,l,ls', 'c,c,cs', 'AV,AV,AV,AV', 'AV,AV,AP,AP', 'AP,AP,AP,AP']:
    r = run7[key]; r2 = run13[key]
    # independent-seed agreement per graph
    m7 = {g['id']: g['log_eps'] for g in r['graphs']}; m13 = {g['id']: g['log_eps'] for g in r2['graphs']}
    agree = (set(m7) == set(m13)) and all({k: v for k, v in m7[i].items()} == {k: v for k, v in m13[i].items()} for i in m7)
    graphs = []
    for g in r['graphs']:
        kappa = ratio_to_tree(g['log_eps'], r['tree_in_basis'])
        ent = dict(id=g['id'], topology=g['topology'], vertices=g['vertices'], internal_lines=g['internal_lines'],
                   symmetry_factor=g['symmetry_factor'], graphs_in_class=g['n_graphs_in_class'],
                   log_eps=('%s x tree vertex' % kappa) if kappa is not None else {'structures': g['log_eps']},
                   eps_inv='0', eps_inv_half='0')
        if kappa is None and key in ('AV,AV,AV,AV', 'AV,AV,AP,AP', 'AP,AP,AP,AP'):
            ent['log_eps']['trace_basis'] = trace_basis(g['log_eps'])
        graphs.append(ent)
        chains['T1'][g['id']] = chain_T1(key, g, kappa)
    kap = sp.sympify(r['computed_multiple'])
    T1['sets'][SET_TXT[key]] = dict(
        tree_vertex=TREE_TXT[key], tree_in_basis=r['tree_in_basis'], graphs=graphs,
        total={'structures': r['total_log_eps'], 'multiple_of_tree': str(kap), 'eps_inv': '0', 'eps_inv_half': '0'},
        predicted_from_A=str(r['predicted_from_A']), verdict=r['verdict'], independent_seed_agrees=agree,
        superficial_degree=1 if key in ('AV,AV,AV', 'AV,AP,AP') else 0)
    legs = key.split(',')
    wmap = {'AV': ws[1], 'AP': ws[2], 'l': ws[3], 'c': ws[4], 'AsV': -ws[1], 'AsP': -ws[2], 'ls': -ws[3], 'cs': -ws[4]}
    eqs.append(sp.Eq(ws[0] + sum(wmap[x] for x in legs), kap))
# quadratic sector (round-1 values reproduced exactly by this engine, logs/r2_regress_quadratic.json)
eqs += [sp.Eq(ws[0] + 2 * ws[1], 4), sp.Eq(ws[0] + 2 * ws[2], 4), sp.Eq(ws[0] + 2 * ws[3], 8), sp.Eq(ws[0] + ws[4] - ws[1], -1)]
sol = sp.solve(eqs, ws, dict=True)
T1['vertex_label_key'] = {
  'AAA': 'cubic part of int dvol -1/4 (F_A,F_A), i.e. int dvol -1/2 (dA, 1/2[A ^ A]) [ESW (31)]',
  'AAAA': 'quartic part int dvol -1/4 (1/2[A ^ A], 1/2[A ^ A])',
  'llA': 'cubic part of int dvol 1/2 (lambda, dslash_A lambda), i.e. 1/2 (lambda, rho([A,lambda]))',
  'AcAs': 'int dvol ([A,c], A^*) (cubic part of (d_A c, A^*))', 'cllS': 'int dvol ([c,lambda], lambda^*)', 'cccS': 'int dvol 1/2 ([c,c], c^*)',
  'bracket after a label': 'the external legs attached to that vertex'}
T1['quadratic_regression'] = {k: v['total_units_u'] for k, v in regress.items()}
T1['vanishing_sets_checked'] = van
T1['test_of_A'] = {
  'statement': 'log eps part of I^CT[eps] = u log eps (Q + {I,-})(2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c))',
  'bracket_rederived': 'with <x_i,xi_i> = +int dvol (x_i,xi_i) for (A,A^*),(lambda^*,lambda),(c^*,c) and {f,g} = Delta_0(fg) - Delta_0(f)g - (-1)^{|f|} f Delta_0(g): '
      '{Psi,-} is a derivation with weights ' + json.dumps(W['Psi_weights']) + ' (exact two-mode computation, r2_T3_weights.py); hence '
      '(Q + {I,-})Psi = {Psi, S} = sum over the terms of S of (sum of the weights of their legs) x term; the quadratic part reproduces the '
      'round-1 log counterterm -int(dA|_V,dA|_V) - int(dA|_Wperp,dA|_Wperp) + 4 int(lambda, dslash lambda) - int(dc,A^*).',
  'predicted_multiples': {SET_TXT[k]: str(run7[k]['predicted_from_A']) for k in run7},
  'computed_multiples': {SET_TXT[k]: run7[k]['computed_multiple'] for k in run7},
  'general_fit': 'log eps part = beta S + {Psi\', S} with Psi\' = w_V int(A|_V,A^*|_V) + w_P int(A|_Wperp,A^*|_Wperp) + (lambda, c weights) '
      '[the most general degree -1, weight-4 local functional]; exact solution of all 15 sector equations: ' + str(sol),
  'coupling_direction': 'beta = 0: the coefficient along S itself (= d/ds of S under (-,-) -> (1+s)(-,-), since every term of S is linear in (-,-)) '
      'vanishes at one loop; the full one-loop log counterterm is (Q + {I,-})-exact.'}
T1['verdict'] = ('Assumption (A) holds exactly in every sector with superficial degree >= 0 (quadratic, all cubic, all quartic; the Spin(6)-odd sets '
                 'vanish); computed multiple = predicted multiple for all 11 sets; no eps^{-1} or eps^{-1/2} pieces in cubic/quartic sectors; '
                 'the coupling-direction coefficient is 0 at one loop.')

# ---------------- T2
names = {('sAV', 'AV'): ('([A,c],A^*)[ext c; A^*|_V on special edge, A|_V on P line]', 'AAA[ext A|_V; A|_V on special edge, A|_V on P line]'),
         ('sAV', 'gh'): ('([A,c],A^*)[ext c; A|_V on special edge, A^*|_V on P line (ghost line, A^* end)]',
                         '([A,c],A^*)[ext A|_V; A^*|_V on special edge, c on P line (ghost line, c end)]'),
         ('sAP', 'AP'): ('([A,c],A^*)[ext c; A^*|_Wperp on special edge, A|_Wperp on P line]', 'AAA[ext A|_V; A|_Wperp on special edge, A|_Wperp on P line]'),
         ('sl', 'l'): ('([c,lambda],lambda^*)[ext c; lambda^* on special edge, lambda on P line]', '1/2(lambda, dslash_A lambda)[ext A|_V; lambda on special edge, lambda on P line]'),
         ('sc', 'gh'): ('1/2([c,c],c^*)[ext c; c^* on special edge, c on P line (ghost line, c end)]', '([A,c],A^*)[ext A|_V; c on special edge, A^*|_V on P line (ghost line, A^* end)]')}
SL = {'sAV': 'A|_V - A^*|_V', 'sAP': 'A|_Wperp - A^*|_Wperp', 'sl': 'lambda^* - lambda', 'sc': 'c^* - c'}
PL = {'AV': 'A|_V - A|_V', 'AP': 'A|_Wperp - A|_Wperp', 'gh': 'c - A^*|_V', 'l': 'lambda - lambda'}
t2graphs = []
for g in T2a['graphs']:
    k = (g['special_edge'], g['P_line'])
    gid = 'D(c,A|_V) | special edge (K_eps - K_0) %s ; P(eps,infty) line %s' % (SL[k[0]], PL[k[1]])
    t2graphs.append(dict(id=gid, vertices=list(names[k]), internal_lines=['special edge ' + SL[k[0]], 'P line ' + PL[k[1]]],
                         symmetry_factor='1 (both orders of the two insertions in (1/2) str(K\' S H P S H) summed)',
                         eps_inv=g['eps_inv_coeff_of_i(gamma,iota_p a)'] + ' x u eps^{-1} int(dc, A|_V)', log_eps='0',
                         finite=g['finite_coeff_of_i p^2 (gamma,iota_p a)'] + ' x u int(d^*dc, d^*A|_V)'))
    chains['T2'][gid] = [
      'amplitude: d_P of D^(0) = (1/2)({I,I}_eps - {I,I}_0) = (1/2){I,I}_{K_eps-K_0}; loop part = (1/2) str(K\'(k)(-1)^{|.|} d^2 I_{v_1} P(k-p_1)(-1)^{|.|} d^2 I_{v_2}) '
      '(identity checked by brute force in 0 dimensions, r2_verify_special0d.py); vertices %s ; %s' % names[k],
      'special edge: K\'(k) = (e^{-eps k^2} - 1) x [-sum (x_i (x) xi_i + xi_i (x) x_i) (x) (-,-)^{-1}] = -int_0^eps dt_0 k^2 e^{-t_0 k^2} x [...]; '
      'P line: int_eps^infty dt_1 e^{-t_1 (k-p_1)^2} x numerator',
      'numerator algebra: M(k) = -k^2 str[...], cubic polynomial in (k, p_1) (exact)',
      'Schwinger integral: int_0^eps dt_0 int_eps^infty dt_1 (16 pi^2)^{-1}(t_0+t_1)^{-2}[M + (4(t_0+t_1))^{-1}(d_k.d_k)M](k = t_1 p_1/(t_0+t_1)) e^{-t_0 t_1 p_1^2/(t_0+t_1)}',
      'eps-expansion: rescale t_0, t_1 by eps; with r = t_1/(t_0+t_1) in [1/2,1] and (t_0+t_1)/eps in [1/r, 1/(1-r)]: '
      'eps^{-1}: int_{1/2}^1 dr (2r-1)(d_k.d_k)M(r p_1)/4 ; eps^0: int_{1/2}^1 dr log(r/(1-r)) [M(r p_1) - r(1-r)p_1^2 (d_k.d_k)M(r p_1)/4]; '
      'no log eps (the special edge confines t_0 to [0,eps]); all log 2 terms cancel inside each class',
      'local part: (%s) u eps^{-1} int(dc, A|_V) + (%s) u int(d^*dc, d^*A|_V)' % (g['eps_inv_coeff_of_i(gamma,iota_p a)'], g['finite_coeff_of_i p^2 (gamma,iota_p a)'])]
t2graphs.append(dict(id='D(c,A|_V) | tree vertex -Q I^CT[eps]', vertices=['-Q I^CT[eps]'], internal_lines=[], symmetry_factor='1',
                     eps_inv='+9/2 x u eps^{-1} int(dc, A|_V)  [-(9/4) u eps^{-1} Q int(A|_V,A|_V), Q int(A|_V,A|_V) = -2 int(dc, A|_V)]',
                     log_eps='0  [Q int(dA|_V,dA|_V) = 0, Q int(dc,A^*) = 0 and Q int(lambda,dslash lambda) = Q int(dA|_Wperp,dA|_Wperp) = 0 in this sector]',
                     finite='0'))
chains['T2']['D(c,A|_V) | tree vertex -Q I^CT[eps]'] = [
  'D^(1) = Delta_eps I - Q I^CT[eps] - {I, I^CT[eps]}_eps; Delta_eps I = 0 (every self-contraction of a vertex of I is a trace of ad, = 0); '
  '{I_3, I^CT_2}_eps is cubic; so only -Q I^CT[eps] enters the quadratic sector',
  'Q = {S_free,-} on the two-mode model (c at +p, A|_V at -p): Q int(A|_V,A|_V) = -2 int(dc, A|_V), Q int(d^*A|_V,d^*A|_V) = -2 int(d^*dc, d^*A|_V), '
  'Q int(dA|_V,dA|_V) = 0, Q int(dc,A^*) = 0 (r2_T2b.py)',
  'local part: +(9/2) u eps^{-1} int(dc, A|_V)']
T2 = {
 'sector': 'external legs one gauge ghost c and one A|_V (ghost number 1, quadratic)',
 'other_quadratic_sectors': 'degree-1 quadratic local functionals have leg content (c, A|_V), (c, A|_Wperp) or (c, lambda): (c, A|_Wperp) has no '
     'Spin(6)-invariant (one vector of W_R^perp; also checked: every special-edge bubble numerator vanishes exactly), (c, lambda) has an odd number of '
     'spinor factors (no Spin(4)-invariant). No other leg content has degree 1 with two legs.',
 'graphs': t2graphs,
 'eps_inv_total': 'special-edge graphs -(9/2) u eps^{-1} int(dc, A|_V) + tree vertex (+9/2) u eps^{-1} int(dc, A|_V) = 0',
 'log_total': '0',
 'finite_total': '-(27/8) u int(d^*dc, d^*A|_V)',
 'cross_check_route_b': {'statement': 'Q P(eps,infty) = K_eps gives Q W_1 = (special-edge graphs) at every eps; in this sector Q W_1 = Q(A|_V self-energy) '
                              '(Q annihilates the (c,A^*) self-energy, which is proportional to p). Exact finite parts of every self-energy class '
                              '(per-edge cutoff; log 2 terms appear in the transverse part only):',
                         'self_energy_totals_units_u': T2b['selfenergy_totals'], 'self_energy_classes': T2b['selfenergy_classes'],
                         'Q_maps_two_mode': T2b['Q_two_mode'],
                         'result': 'longitudinal finite coefficient 27/8 -> Q W_1 = -(27/8) u int(d^*dc, d^*A|_V) and eps^{-1}: -(9/2) u eps^{-1} int(dc,A|_V): '
                                   'identical to the special-edge sum (exact agreement, including the cancellation of all log 2 terms)'},
 'obstruction': '-(27/8) u int(d^*dc, d^*A|_V) = (27/16) u Q int(d^*A|_V, d^*A|_V)',
 'q_exact': True,
 'finite_counterterm': 'I^CT[eps] -> I^CT[eps] + (27/16) u int(d^*A|_V, d^*A|_V) (i.e. J = I - I^CT[eps] - (27/16) u int(d^*A|_V,d^*A|_V)); '
     'the obstruction in this sector then vanishes. Equivalently: the renormalised A|_V self-energy carries the finite local longitudinal term '
     '(27/16) u int(d^*A|_V, d^*A|_V), which this counterterm removes. No int(A|_V,A|_V) finite counterterm is needed (eps^{-1} cancels exactly).',
 'uniqueness': 'the only dimension-4 local functionals with leg content (c, A|_V) are multiples of int(d^*dc, d^*A|_V) (int(dc, d^*dA|_V) = 0); '
     'it is Q-exact, so the quadratic obstruction is always removable by a multiple of int(d^*A|_V,d^*A|_V).',
 'cubic_sector': 'OPEN (not computed). Reduction established here: after the quadratic finite counterterm the cubic part of the obstruction is the '
     'eps -> 0 limit of the special-edge 1PI loops with three legs (leg contents (c,A|_V,A|_V), (c,A|_Wperp,A|_Wperp), (c,lambda,lambda), (c,c,A^*|_V)); '
     'special-edge bridges cancel against the bridge part of -{I,I^CT}_eps; s I^CT[eps] has no cubic part (log part: s^2 Psi = 0 by T1; eps^{-1} part: '
     '{I, int(A,A)} = 2 int([A,c],A) = 0). The three-edge special-edge region integrals (t_0 in [0,eps], two t_e in [eps,infty), Taylor order 2 in '
     'external momenta) are not done.',
 'verdict': 'Quadratic (c,A|_V) sector: eps^{-1} cancels exactly against -Q of (9/4)u eps^{-1} int(A|_V,A|_V); no log eps; finite obstruction '
     '-(27/8)u int(d^*dc,d^*A|_V) = (27/16)u Q int(d^*A|_V,d^*A|_V): Q-exact, removed by the finite counterterm (27/16)u int(d^*A|_V,d^*A|_V).'}

# ---------------- T3
sw3 = {'traceless scalar bilinear (traceless part of int omega (iota_{w1}A, iota_{w2}A))': -0,
       'trace scalar bilinear (sum_w int omega (iota_w A, iota_w A))': -24, 'letter iota_w A (w in W_R^perp)': -4,
       'letter A|_V (transverse)': -4, 'letter lambda': -8, 'letter c': -2}
T3 = {'normalisations': 'SW (3): Delta = -E; SW (17): Delta = 2E (superconformal/spec.md); conversion: coefficient(17) = -2 x coefficient(3)',
      'table_units_u_coefficient_multiplying_O': {k: {'SW(3)': v, 'SW(17)': -2 * v} for k, v in sw3.items()},
      'letter_c': 'computed: single letters have no one-loop 1PI graph, so the log-eps coefficient of partial_delta I^CT[eps] is 0; {Psi, c} = +1 x c '
                  '(exact two-mode computation, r2_T3_weights.py); coefficient = 2(0 - 1) = -2 (SW (3)) = +4 (SW (17)). Cross-check: '
                  'the (c,A^*) two-point log coefficient -1 = w_c - w_A with w_A = 2 gives w_c = 1.',
      'Psi_weights': W['Psi_weights'],
      'Konishi': 'trace: -24 u (SW 3) = +48 u (SW 17); with u = -g^2 N/(32 pi^2) (standard su(N) normalisation): 3 g^2 N/(4 pi^2) (SW 3), -3 g^2 N/(2 pi^2) (SW 17).'}

res = {'engine': 'A (round 2): momentum space, per-edge Schwinger parameters t_e in [eps,infty), exact Gaussian loop integral, exact eps-asymptotics, '
                 'exact arithmetic throughout (no floating-point fits)',
       'files': ['engine_a/r2_core.py', 'engine_a/r2_graphs.py', 'engine_a/r2_T1.py', 'engine_a/r2_T1_vanishing.py', 'engine_a/r2_regress_quadratic.py',
                 'engine_a/r2_verify_special0d.py', 'engine_a/r2_T2.py', 'engine_a/r2_T2b.py', 'engine_a/r2_T3_weights.py', 'engine_a/r2_make_results.py'],
       'T1': T1, 'T2': T2, 'T3': T3,
       'formula_chains': {'common_schwinger_formula': 'int d^4k/(2pi)^4 N(k) exp(-sum_e t_e (k+q_e)^2) = (16 pi^2)^{-1}(sum_e t_e)^{-2} sum_m (4^m m!)^{-1}(sum_e t_e)^{-m} '
                                                       '[(d_k.d_k)^m N](k = -sum_e t_e q_e/sum_e t_e) exp(-sum_e t_e q_e^2 + (sum_e t_e q_e)^2/sum_e t_e)',
                          'T1': chains['T1'], 'T2': chains['T2']}}
json.dump(res, open('results_r2.json', 'w'), indent=1)
print('general fit:', sol)
print('seed agreement:', {k: v['independent_seed_agrees'] for k, v in T1['sets'].items()})
