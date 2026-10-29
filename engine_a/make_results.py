"""Assemble engine_a/results.json from the logs (all numbers in units hbar c_adj/(16 pi^2) unless stated)."""
import json
L = lambda f: json.load(open('logs/' + f))
ga2, ga3 = L('graphs_a_su2.json'), None
sa2, sa3 = L('sector_a_su2.json'), L('sector_a_su3.json')
sb2, sb3 = L('sector_b_su2.json'), L('sector_b_su3.json')
fc, bv, bw = L('faces_su2.json'), L('bv_exact.json'), L('bv_ward.json')
VERT = {'AAA+AAA': ['AAA', 'AAA'], 'AcAs+AcAs': ['AcA*', 'AcA*'], 'llA+llA': ['lambda lambda A', 'lambda lambda A'],
        'AAAA': ['AAAA'], 'O+AAAA': ['O', 'AAAA'], 'O+AAA+AAA': ['O', 'AAA', 'AAA']}
graphs = []
for e in ga2['graphs']:
    ext, vt, lines = e['id'].split(':')
    ent = {'id': e['id'], 'sector': e['sector'], 'vertices': VERT[vt], 'internal_lines': lines.split(','),
           'symmetry_factor': '1/2 (n=2 term 1/4 str[(P(-1)^|.|d2I)^2], both cyclic orders)' if e['kind'] == 'bubble' else '1/2 (n=1 term 1/2 str[P(-1)^|.|d2I_4])',
           'eps_inv': e['eps_inv'], 'log_eps': e['log_eps'], 'notes': ''}
    if e['id'] in fc:
        ent['integrand_and_faces'] = fc[e['id']]
    if e['id'] in ('Aperp-Aperp:AcAs+AcAs:c,As', 'As-c:AcAs+AcAs:Aperp,c-As'):
        ent['notes'] = 'vanishes: the W^perp part of A^* has no propagator (P has only the c-A^*|_V component)'
    graphs.append(ent)
for e in sb2['graphs']:
    ext, vt, lines = e['id'].split(':')
    graphs.append({'id': e['id'], 'sector': 'O-insertion, external ' + ext, 'vertices': VERT[vt], 'internal_lines': lines.split(','),
                   'symmetry_factor': '1/2 (delta-linear term 1/2 str[(P(-1)^|.|d2I)^n P(-1)^|.|d2O], both orders for n=2)',
                   'eps_inv': {'O': '0', 'trace': '0', 'AV_trace': '0'}, 'log_eps': e['log_eps'],
                   'notes': 'structures: O = int omega (iota_w1 A, iota_w2 A); trace = (w1,w2) int omega sum_w (iota_w A, iota_w A); AV_trace = (w1,w2) int omega (A_V, A_V). UV-boundary split: each of the 3 (2) edges carries 1/3 (1/2) of the log coefficient (constant numerator).'})
graphs.append({'id': 'vacuum:O:Aperp', 'sector': 'O-insertion, no external field', 'vertices': ['O'], 'internal_lines': ['Aperp'],
               'symmetry_factor': '1/2', 'eps_inv': {'(w1,w2) dim(g) int omega [units hbar/(16 pi^2)]': '2'}, 'log_eps': {}, 'notes': 'one-edge graph: int_eps^infty dt t^-2 = 1/eps; not proportional to c_adj'})
tot = sa2['sectors']
ICT = {}
for sec, Bs in (('AV-AV', ['B1', 'B2', 'B3']), ('Aperp-Aperp', ['B4', 'B5']), ('l-l', ['B6']), ('As-c', ['B7'])):
    for B in Bs:
        ICT[B] = {'eps_inv': tot[sec]['eps_inv_rational'][B], 'log_eps': tot[sec]['log_eps_rational'][B]}
ICT['mixing A_V-A_perp'] = {'eps_inv': '0', 'log_eps': '0'}
res = {
 'conventions': ('Engine A. Euclidean R^4, weight e^{+S/hbar}, S = [ESW (31)] reduced from 10d (coefficient -1/4 of the Hodge pairing (F_A,F_A)), '
   'pairing <x,xi> = +int dvol (x,xi) for (A,A^*),(lambda^*,lambda),(c^*,c) [classical_spec section 3]. Q = {S_free,-}; Q^GF of DESIGN form fixed by '
   'partial_P(Q(y_a y_b)) = partial_{K_eps}(y_a y_b) (C07 Lemma 10.2.1, T=infty): P(eps,infty)(p) = e^{-eps|p|^2}|p|^{-2} x {2 (A x A, all ten directions), -i p (c x A^*|_V and A^*|_V x c), i rho(p) (lambda x lambda)} x (-,-)^{-1}. '
   'One loop: W_1 = sum_n 1/(2n) str[(P (-1)^{|.|} d2I)^n] (verified by brute force in a 0-dim super model, tadpole-free). Each internal edge t_e in [eps,infty); Gaussian loop integral exact; '
   'eps-asymptotics exact by the face (boundary) method. I^CT[eps] = Sing_eps of the one-loop 1PI graphs (I - I^CT the regularized interaction). '
   'Fourier Phi(x) = int d^4p/(2pi)^4 e^{ip.x} Phi(p). Colour: orthonormal basis of g for (-,-) in code only; c_adj = the scalar by which X -> sum_a [t_a,[t^a,X]] acts (t^a dual basis w.r.t. (-,-)); '
   'every result checked for g = su(2) and su(3) (identical rational multiples of c_adj). Units: hbar c_adj/(16 pi^2), separately for eps^{-1} and log eps. Realization: 10d chiral 16x16 gamma blocks (9 real symmetric + one i*identity in W^perp).'),
 'graphs': graphs,
 'ICT_quadratic': ICT,
 'O_eps': {
   'scalar_bilinear': {
     'definition': 'O = int omega (iota_{w1}A, iota_{w2}A); O(eps) = O - partial_delta I^CT[eps]',
     'partial_delta_ICT_log_eps': {'O': '4', '(w1,w2) int omega sum_w (iota_w A, iota_w A)': '-2', '(w1,w2) int omega (A_V,A_V)': '0'},
     'partial_delta_ICT_eps_inv': {'vacuum: (w1,w2) dim(g) int omega, units hbar/(16 pi^2)': '2', 'field-dependent': '0'},
     'traceless_part': 'O(eps) = O - 4 hbar c_adj/(16 pi^2) log(eps) O',
     'trace_part': 'sum_w O_ww (eps) = (1 + 8 hbar c_adj/(16 pi^2) log eps) sum_w O_ww - 12 hbar dim(g)/(16 pi^2 eps) int omega',
     'other_insertions': 'external (A_V,A_perp): 0 (Spin(6) selection); one external field: 0 (colour, f_abb); external lambda lambda, (A^*,c): superficial degree < 0',
     'su3_check': sb3['fits']},
   'single_letters': {'int omega (X, iota_w A)': 'O(eps) = O (no 1PI one-loop graph: a univalent vertex cannot lie on a loop; the reducible graphs are made finite by the quadratic I^CT)',
                      'int omega (X, A_V(v))': 'O(eps) = O (same reason)', 'lambda, c letters': 'O(eps) = O (same reason)'}},
 'dilatation': {
   'definition': 'Delta = -E on R^4, acting on a field of Weyl weight w as -(E+w); weights A 1, lambda 3/2, c 0, A^* 3, lambda^* 5/2, c^* 4',
   'traceless': '0', 'trace': '-24',
   'letters': {'iota_w A (W^perp)': '-4', 'A_V': '-4 (plus longitudinal/other structures: Q O != 0)', 'lambda': '-8', 'c': '-2',
               'note': 'Feynman-gauge (gauge dependent); not gauge invariant letters also receive extra structures from (Q + {I,-})O != 0'},
   'route_i': ('Scaling map (D_lambda phi)(x) = lambda^w phi(lambda x): W(P(eps,infty),F) o D_lambda = W(P(lambda^2 eps,infty), F o D_lambda) (kernel scaling of P). '
     'Hence for O[infty]: (E-generator) O[infty] = O_{omega-dot}[infty] + 2 hbar ( [log-eps coefficient of partial_delta I^CT] + [tree insertion of the log-eps coefficient of I^CT into O] ). '
     'Assumption (A): the log-eps coefficient of I^CT is (Q+{I,-}) of the local functional 2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c) [units hbar c_adj/16pi^2] (checked: unique solution at the quadratic level). '
     'Assumption (B): the Delta-component of the quantum action is corrected by 2 hbar times that functional so that the equivariant QME holds along Delta at one loop. '
     'Tree identity (verified, two-mode model): insertion of Q(Psi) into O = Q(insertion of Psi into O) - {Psi, O}. Result: gamma_O O = 2 hbar ( log-coefficient of partial_delta I^CT - {Psi,O} ). '
     '{Psi,O}: bilinear 4 O, iota_wA 2 O, A_V 2 O, lambda 4 O, c 1 O. Vertex log coefficients: traceless 4, trace -8, letters 0. '
     'Cross-check: the letter values equal the anomalous exponents of the renormalized two-point functions: from the log-eps counterterms -1 (dA,dA) and 4 (lambda, dslash lambda) the A and lambda propagators scale with extra exponents -4 and -8 (scaling dimensions 1 + gamma, 3/2 + gamma).'),
   'route_ii': ('UV-boundary form: d/dlambda P(eps/lambda^2,infty) = 2 eps (Q^GF x 1) K_eps, i.e. the dilatation acts on each internal line as the special edge at t_e = eps; '
     'sum over special edges of a graph = 2 eps d/deps of the graph. Evaluated directly as face integrals (t_e = eps, other t >= eps): '
     'each face of every two-edge self-energy graph carries exactly 1/2 of its log coefficient (integrands symmetric under t1<->t2 after summing both orders), eps^{-1} faces c/8 each cancel against eps d/deps of the eps^{-1} counterterm; '
     'O-insertion triangle: each of the three faces carries 1/3, O+AAAA bubble: each face 1/2. Sum of UV-boundary graphs reproduces 2 x (log coefficients) used in route (i) exactly; with (A),(B) the same gamma follows: traceless 0, trace -24, letters iota_wA -4, lambda -8, c -2.')},
 'selection_rules': [
   'Translations, Spin(4) rotations and Spin(6)=SL(4) R-symmetry: H = Laplacian x id and Q^GF are covariant, so K_eps and P(eps,infty) are invariant; the corresponding components of (d_CE+Q)O[infty] are the classical actions, no one-loop correction; counterterms are invariant (checked: no A_V-A_perp mixing, only Spin(6)-invariant tensors in partial_delta I^CT).',
   'Dilatation: only log-eps counterterms contribute (eps^{-1} terms are scale covariant and drop); one-loop gamma as above.',
   'Superficial degree: 1PI with one O insertion: degree = w_O - sum of external weights; scalar bilinear: two A externals degree 0 (log), one external degree 1 vanishes by colour (f contracted with (-,-)^{-1}), lambda lambda and (A^*,c) externals negative degree.',
   'Q^i_alpha, Qbar (weight +1/2 relative to O): output local functionals of weight w_O + 1/2, R-rep in 6 x 4 (or 4bar), chirality S_+ (S_-), ghost degree 0 coefficient: int omega (iota_w A, Gamma-contraction with lambda) [classical], and possible int omega ([c, .], lambda^*)-type terms; the Feynman-type Q^GF is not supersymmetric, so UV-boundary corrections are allowed only in these structures.',
   'S, Sbar (weight -1/2): structures of weight w_O - 1/2 with one explicit x (from the conformal Killing field); K (weight -1): structures with x^2 / two x; ghost-degree and Lorentz/R selection as above.',
   'Ghost number: I^CT has degree 0: quadratic sector only (A,A), (lambda,lambda), (A^*,c); (lambda^*,lambda), (c^*,c) excluded by degree.']
}
json.dump(res, open('results.json', 'w'), indent=1)
print('written', len(graphs), 'graphs')
