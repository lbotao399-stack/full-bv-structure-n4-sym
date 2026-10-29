"""Assemble engine B results.json and report.md from the computed data files."""
import json
import sympy as sp

tp = json.load(open('two_point_all.json'))
og = json.load(open('o_graphs.json'))
rd = json.load(open('report_data.json'))
rdo = json.load(open('report_data_O.json'))
r2 = json.load(open('route2_numeric.json'))

def clean(d):
    if isinstance(d, str):
        return {} if d.startswith('NONLOCAL:0') else d
    return {k: v for k, v in d.items() if v not in ('0', 0)}

meta = {
 'AV-AV:AAA+AAA:AV,AV': (['AAA', 'AAA'], ['A_V', 'A_V'], '1/2', 'gauge loop'),
 'AV-AV:AAA+AAA:Ap,Ap': (['AAA', 'AAA'], ['A_perp', 'A_perp'], '1/2', 'scalar loop through the A_V (A_perp)^2 cubic term'),
 'AV-AV:AcAs+AcAs:c-As,c-As': (['AcA*', 'AcA*'], ['c-A*_V', 'A*_V-c'], '1', 'ghost loop; both orientations are the two blocks c(x)A*, A*c of P'),
 'AV-AV:lamlamA+lamlamA:lam,lam': (['lam lam A', 'lam lam A'], ['lam', 'lam'], '1/2', 'gaugino loop (lam (x) lam propagator)'),
 'AV-AV:AAAA:AV': (['AAAA'], ['A_V'], '1/2', 'quartic tadpole'),
 'AV-AV:AAAA:Ap': (['AAAA'], ['A_perp'], '1/2', 'quartic tadpole'),
 'Ap-Ap:AAA+AAA:Ap,AV': (['AAA', 'AAA'], ['A_perp', 'A_V'], '1', 'gauge exchange'),
 'Ap-Ap:AcAs+AcAs:c-As,c-As': (['AcA*', 'AcA*'], ['c-A*_perp'], '-', 'identically zero: P has no c (x) A*_perp component'),
 'Ap-Ap:lamlamA+lamlamA:lam,lam': (['lam lam A (Yukawa)', 'lam lam A (Yukawa)'], ['lam', 'lam'], '1/2', 'gaugino loop'),
 'Ap-Ap:AAAA:AV': (['AAAA'], ['A_V'], '1/2', 'quartic tadpole'),
 'Ap-Ap:AAAA:Ap': (['AAAA'], ['A_perp'], '1/2', 'quartic tadpole'),
 'AV-Ap:AAA+AAA:all': (['AAA', 'AAA'], ['all'], '-', 'identically zero (no A_V A_V A_perp vertex, so(6) selection)'),
 'AV-Ap:AAAA:all': (['AAAA'], ['all'], '-', 'identically zero'),
 'lam-lam:lamlamA+lamlamA:lam,AV': (['lam lam A', 'lam lam A'], ['lam', 'A_V'], '1', 'gauge exchange'),
 'lam-lam:lamlamA+lamlamA:lam,Ap': (['lam lam A (Yukawa)', 'lam lam A (Yukawa)'], ['lam', 'A_perp'], '1', 'scalar exchange'),
 'lam-lam:all': (['lam lam A', 'lam lam A'], ['all'], 'sum', 'check: equals the sum of the two lam-lam graphs'),
 'As-c:AcAs+AcAs:AV,c-As': (['AcA* (A* external)', 'AcA* (c external)'], ['A_V', 'c-A*_V'], '1', 'ghost self-energy'),
 'As-c:all': (['AcA*', 'AcA*'], ['all'], 'sum', 'check: equals the As-c graph; A*_perp external legs give zero'),
}
graphs = []
for gid, v in tp.items():
    vs, ls, sym, note = meta[gid]
    graphs.append(dict(id=gid, sector=v['sector'], vertices=vs, internal_lines=ls, symmetry_factor=sym,
                       eps_inv=clean(v['eps_inv']), log_eps=clean(v['log_eps']), notes=note,
                       units='hbar c_adj/(16 pi^2); coefficient of the basis functional in I^CT[eps]'))
ometa = {
 'Ap-Ap+O:AAA+AAA+O:Ap,AV,Ap': (['O', 'AAA', 'AAA'], ['A_perp', 'A_V', 'A_perp'], '1', 'vertex correction by gauge exchange'),
 'Ap-Ap+O:AAAA+O:Ap,Ap': (['O', 'AAAA'], ['A_perp', 'A_perp'], '1/2', 'quartic (A_perp)^4 bubble'),
 'AV-AV+O:AAA+AAA+O:Ap,Ap,Ap': (['O', 'AAA', 'AAA'], ['A_perp', 'A_perp', 'A_perp'], '1', 'mixing into int omega (A_V,A_V)'),
 'AV-AV+O:AAAA+O:Ap,Ap': (['O', 'AAAA'], ['A_perp', 'A_perp'], '1/2', 'mixing into int omega (A_V,A_V)'),
}
for gid, v in og.items():
    if gid.startswith('vac'):
        graphs.append(dict(id='vac:O:Ap', sector='O vacuum', vertices=['O'], internal_lines=['A_perp'], symmetry_factor='1/2 (x2 orderings)',
                           eps_inv={'int omega': '2 dim(g) (w1,w2)  [units hbar/(16 pi^2), not proportional to c_adj]'}, log_eps={},
                           notes='self-contraction of O; P^{AA}(x,x) = 2 (-,-)^{-1} (metric)^{-1} /(16 pi^2 eps)'))
        continue
    vs, ls, sym, note = ometa[gid]
    graphs.append(dict(id=gid, sector=v['sector'], vertices=vs, internal_lines=ls, symmetry_factor=sym,
                       eps_inv={k: x for k, x in v['eps_inv'].items() if x != '0'}, log_eps={k: x for k, x in v['log_eps'].items() if x != '0'},
                       notes=note + '; structures: O = int omega (iota_w1 A, iota_w2 A); trace_mix = (w1,w2) int omega sum_w (iota_w A, iota_w A); AV_mix = (w1,w2) int omega (A_V,A_V)',
                       units='hbar c_adj/(16 pi^2); coefficient in partial_delta I^CT[eps]'))

ICT = {'B1': {'eps_inv': '0', 'log_eps': '-1'}, 'B2': {'eps_inv': '0', 'log_eps': '0'}, 'B3': {'eps_inv': '9/4', 'log_eps': '0'},
       'B4': {'eps_inv': '0', 'log_eps': '-1'}, 'B5': {'eps_inv': '9/4', 'log_eps': '0'}, 'B6': {'eps_inv': '0', 'log_eps': '4'},
       'B7': {'eps_inv': '0', 'log_eps': '-1'},
       'other structures': 'none: the A_V-A_perp, lam-A, (A*_perp, c) and all degree-0 mixed sectors vanish identically; no eps^{-1/2} terms (odd Gaussian moments vanish)'}
# totals check from graphs
tot = {}
for g in graphs:
    if g['sector'] in ('AV-AV', 'Ap-Ap', 'lam-lam', 'As-c') and g['symmetry_factor'] != 'sum':
        for kind in ('eps_inv', 'log_eps'):
            for k, v in g[kind].items():
                tot.setdefault((k, kind), 0)
                tot[(k, kind)] += sp.Rational(v)
for (k, kind), v in tot.items():
    assert str(v) == ICT[k][kind], (k, kind, v, ICT[k][kind])

conventions = (
 "Engine B. Weight e^{+S/hbar}; S = int dvol(-1/4 (F_A,F_A) + 1/2 (lambda, dslash_A lambda) + (d_A c, A*) + ([c,lambda],lambda*) + 1/2([c,c],c*)) "
 "[ESW (31)] dimensionally reduced to Euclidean R^4 (d only along R^4, A in V_10^vee (x) g, A_V its R^4 part, A_perp = W_R^perp part). "
 "Explicit realization (code only): S_+, S_- two-component blocks, sigma = (1, i sigma_k): S_+ -> S_-; W = C^4, W_R^perp realised on Wedge^2 W by "
 "antisymmetric 't Hooft blocks; Gamma^mu = gamma^mu (x) 1, Gamma^{4+i} = gamma_5 (x) gamma_6^i; Sigma = S_+(x)W + S_-(x)W^*. "
 "Contraction rule d_P(f_X f_Y) = P^{XY}; P(eps,infty) = N(q) e^{-eps q^2}/q^2 with N_AA = 2 (-,-)^{-1} (all ten components), "
 "N_{lam lam} = i rho(q)|_{Sigma*->Sigma} read in Sigma(x)Sigma via (-,-), N_{c A*} = N_{A* c} = -i q (-,-)^{-1} (A*_V only), all derived as "
 "-(Pi Hessian of the free action)^{-1} (Feynman gauge for A, ghost sector on Im Q^GF). One-loop formula (0-dim verified): "
 "W_1 = (hbar/2) sum_n (1/n) STr[(P Pi d^2 I)^n], Pi = parity operator. Fourier Phi(x) = int d^4p/(2pi)^4 e^{ip.x} Phi(p). "
 "c_adj = the scalar by which X -> sum_a [t_a,[t^a,X]] acts on simple g (t^a dual to t_a w.r.t. (-,-)); every colour factor reduces to "
 "c_adj (c1,c2) (bubbles: +1, quartic tadpoles: -1 per ([c1,t],[c2,t^]) contraction, checked on su(2), su(3) in non-orthonormal bases "
 "with rescaled pairings). Units: hbar c_adj/(16 pi^2) for eps^{-1} and log eps coefficients. I^CT[eps] = singular part of W_1 "
 "(I - I^CT[eps] is the regularised interaction). Method B: eps d/deps of each graph = -(sum over edges of the graph with that edge pinned at "
 "t_e = eps); pinned graphs expanded in eps (Taylor in external momenta), remaining t_f = eps u_f, exact u-integrals over [1,infty).")

O_eps = {
 'scalar_bilinear': {
   'definition': 'O = int omega (iota_w1 A, iota_w2 A), w1,w2 in W_R^perp; O(eps) = O - partial_delta I^CT[eps]',
   'partial_delta_ICT': {
     'O': {'eps_inv': '0', 'log_eps': '4'},
     '(w1,w2) int omega sum_w (iota_w A, iota_w A)  [orthonormal basis of W_R^perp]': {'eps_inv': '0', 'log_eps': '-2'},
     '(w1,w2) int omega (A_V, A_V)': {'eps_inv': '0', 'log_eps': '0 (graph AV-AV+O:AAA+AAA+O gives +2, AV-AV+O:AAAA+O gives -2)'},
     'constant (vacuum)': {'eps_inv': '2 dim(g) (w1,w2) int omega   [units hbar/(16 pi^2), NOT proportional to c_adj]', 'log_eps': '0'},
   },
   'traceless_symmetric_part': {'partial_delta_ICT_log_eps': '4 O', 'eps_inv': 'none'},
   'trace_part (O_tr = sum_w (iota_w A, iota_w A))': {'partial_delta_ICT_log_eps': '-8 O_tr', 'eps_inv': '2 dim(g) dim(W_R^perp)=12 dim(g) int omega [units hbar/16pi^2]'},
 },
 'single_letters': {
   'int omega (X, iota_w A)': 'partial_delta I^CT[eps] = 0 at one loop (no 1PI graph: a univalent vertex cannot lie on the loop; the tadpole loop vanishes by colour)',
   'int omega (X, A_V(v))': 'partial_delta I^CT[eps] = 0 at one loop (same reason)',
   'int omega (X, (zeta, lambda)), zeta a fixed element of Sigma^*': 'partial_delta I^CT[eps] = 0 at one loop (same reason)',
 },
}

dil = {
 'normalisation': 'Delta = -E on functions, acting on a field of weight w as -(E+w); weights A:1, lambda:3/2, c:0, A*:3, lambda*:5/2, c*:4',
 'traceless': {'gamma': '0', 'kappa_classical_Delta_on_functional': '-8',
               'decomposition': 'gamma = 2*(4) [vertex log] + 2*(-1)*8 [quadratic counterterm B4 attached to both legs by P(0,infty)] - (-8) [equivariant compensation -Q O_1, O_1 from I^(1)(Delta) = -2 hbar R] = 8 - 16 + 8 = 0'},
 'trace': {'gamma': '-24', 'kappa_classical_Delta_on_functional': '-32',
           'decomposition': 'gamma = 2*(4 - 12) + 16*(-1) + 8 = -24; in terms of the standard coupling (hbar c_adj = -g^2 C_A/2 for (-,-) = (2 hbar/g^2) x positive Killing-normalised form) this is 12 g^2 C_A/(16 pi^2) = 3 g^2 C_A/(4 pi^2) [sanity check only]'},
 'letters': {
   'int omega (X, iota_w A)': {'gamma': '-4', 'kappa': '-8', 'formula': 'gamma = -b_4/(-1/4), b_4 = -1 the log eps coefficient of B4, -1/4 the coefficient of B4 in the free action'},
   'int omega (X, A_V(v))': {'gamma_transverse': '-4', 'gamma_longitudinal': '0',
       'other_structures': '-4 hbar c_adj/(16pi^2) int omega (X, (d (d^*d)^{-1} d^* A_V)(v)) i.e. [(d_CE+Q)O[infty]](Delta) = -(E+1)O + 4 u (O - O_longitudinal), u = hbar c_adj/16pi^2; the longitudinal projection is nonlocal and gauge (Q-direction) only'},
   'int omega (X, (zeta, lambda))': {'gamma': '-8', 'formula': 'gamma = -b_6/(1/2), b_6 = 4'},
   'remark': 'letters are gauge dependent (Feynman gauge only here); gamma_letter = -(log eps coefficient)/(free-action coefficient) of the corresponding quadratic basis element',
 },
}

json.dump(dict(conventions=conventions, graphs=graphs, ICT_quadratic=ICT, O_eps=O_eps, dilatation=dil), open('results.json', 'w'), indent=1)
print('results.json written; totals verified against per-graph sums')

# ---- dilatation routes and selection rules ----
res = json.load(open('results.json'))
res['dilatation']['route_i'] = {
 'formula': "Delta.O[infty] = (Delta.O)[infty] - 2 hbar (d_delta I^CT)_log[infty] - 2 hbar sum_k b_k d^cross_{P(0,infty)}(B_k, O)   [from D_lambda^* P(eps,infty) = P(lambda^2 eps,infty), Delta.P = 2 eps L_eps, eps-independence of O[infty]; Delta.B3,5 = 2 B3,5 kill all eps^{-1} terms]; "
            "[(d_CE+Q)O[infty]](Delta) = Delta.O[infty] - Q O_1[infty], O_1 = -2 hbar d^cross_P(R,O), where sum_k b_k B_k = Q R, R = 2 int (A,A*) + (lambda, c parts), forced by the xi-linear equivariant master equation at scale infty: Q I^(1)(Delta) = Delta.I[infty] = -2 hbar Q R",
 'scalar_bilinear': "d^cross(B4,O) = 8 O; -Q O_1 = -8 u O; result: -(E+2)O + u[(-8+16-8) O + 4 (w1,w2) O_tr]  (u = hbar c_adj/16pi^2) => gamma_traceless = 0, gamma_trace = -24",
 'general': "gamma_O = 2 b'_O + 8 b_4 (scalar bilinears), gamma_letter = -b_X/k_X (k_X = coefficient of B_X in the free action)",
}
res['dilatation']['route_ii'] = {
 'method': 'direct UV-boundary graphs: every edge replaced in turn by Delta.P = 2 eps L_eps (pinned at t_e = eps), exact Gaussian factor kept, finite eps (2e-2 ... 6.25e-4 for O graphs, 2e-3 ... 1.25e-4 for self-energy), random rational external momenta/polarisations, eps->0 by fit to 2a/eps - 2b + eps(c log eps + d)',
 'self_energy_Aperp_all_graphs': r2['self_energy_Aperp'],
 'O_insertion_1PI_Aperp': {k: v for k, v in r2['O_vertex_Aperp'].items()},
 'assembly': 'pinned 1PI O-loops -> -2u[4 O - 2(w1,w2) O_tr]; reducible chains O-P-(self-energy - I^CT): pinned loop edges + Delta.(-I^CT) -> -2 b_4 x (B4 kernel) per leg = +8u per leg (2a/eps pieces cancel exactly), pinned tree edge (Sigma - I^CT) 2 eps L_eps -> 0; plus the tree -Q O_1 = -8u O; total identical to route (i): gamma_traceless = 0, gamma_trace = -24',
 'agreement': 'relative 3e-6 (self-energy finite part), 3.7e-6 (O vertex), pole 5e-10',
}
res['selection_rules'] = [
 "P_mu, M_munu, R^i_j (so(6)=sl(4)): H = Laplacian (x) id and Q^GF are Euclidean- and R-covariant, so (X(x)1+1(x)X)P(eps,infty)=0, X.K_eps=0, counterterms X-invariant => [(d_CE+Q)O[infty]](X) = (X.O)[infty] exactly: no one-loop correction for any local O; traceless and trace parts of Sym^2(W_R^perp) cannot mix.",
 "Delta: regulator not covariant (Delta.P = 2 eps L_eps); computed: gamma_traceless = 0, gamma_trace = -24, letters gauge dependent.",
 "K_mu (dimension -1, component of ghost degree |O|): anomalous part fixed by [K_mu,P_nu] = 2M_munu + 2 delta_munu E [SW (4)] to -2 gamma_O int omega x_mu O; a translation-invariant extra term would need a colour-singlet local vector of dimension w_O - 1: none exists for scalar bilinears or letters (the only dimension-1 vector is the adjoint A_V) => output = -2 gamma_O int omega x_mu O + Q-exact.",
 "Q^i_alpha, Qbar_{alphadot i} (dimension +1/2, odd; component degree |O|): allowed outputs = fermionic colour-singlet local functionals of dimension w_O + 1/2 in S_+(x)W (resp. S_-(x)W^*) (x) R-rep of O: for the bilinears int omega (iota_w A, lambda) [classical type], gauge-variant int omega (A_V, rho(.) lambda), BRST-type int omega (c, lambda*) (dimension 0 + 5/2, degree 0); nothing else. Not excluded by covariance (linear supersymmetry does not commute with Q^GF); any such term is supported at t = eps.",
 "S^alpha_i, Sbar^{alphadot i} (dimension -1/2): a translation-invariant output would need a fermionic colour-singlet local operator of dimension w_O - 1/2 = 3/2: none (lambda is adjoint, (c,lambda*) has dimension 5/2) => only x-weighted Q-type structures fixed by {K,Q} ~ S, plus Q-exact terms.",
]
res['checks'] = [
 '10d Clifford, Gamma11, pairing symmetry, ESW Assumption 3.1 (realization.py)',
 '0-dim Grassmann verification of W_1 = (hbar/2) sum (1/n) STr[(P Pi d^2 I)^n] with d_P(f_X f_Y) = P^{XY}, incl. odd backgrounds, quartic and O insertions (toy_verify*.py); 7 other variants fail',
 'propagators derived by inverting the free Hessian (lam: i rho(q)/q^2; A: 2/q^2; c-A*: -i q/q^2)',
 'colour: bubble = +c_adj(c1,c2), quartic tadpole = -c_adj(c1,c2) per ([c1,t],[c2,t^]) contraction, on su(2) and su(3), non-orthonormal bases, rescaled pairings',
 'B2 log coefficient 0 (gauge -1/8 + ghost +1/8): transverse log part',
 '10d covariance of I^CT_quad: B1 and B4 log coefficients equal (-1), B3 and B5 pole coefficients equal (9/4)',
 'Q-exactness of the log part: sum b_k B_k = Q(2 int (A,A*) + ...)',
 'per-graph Feynman-parameter form (B) equals the direct pinned-edge u-integrals (report_data.py)',
 'lam-lam and As-c sector totals reproduced by an all-blocks run',
 'O-insertion: int omega (A_V,A_V) mixing cancels between triangle (+2) and bubble (-2)',
 'route (ii) numerics at finite eps and finite momenta agree with the eps-expansion to <4e-6',
]
res['formula_chain_data'] = {'two_point_words': rd, 'O_words': rdo}
json.dump(res, open('results.json', 'w'), indent=1)
print('results.json updated')
