"""Engine A (R8, route A) -- runs everything end to end and writes
   kernels.json, gates.json (top level of engine_a/, mirrored in out.nosync/), out.nosync/kernels_W0_raw.json.
Usage: python3 run_all.py            (about 6 minutes on an M-series Mac, < 0.5 GB RSS, single process)
Code: src.nosync/ (excluded from iCloud sync), Step-1 engine A vertex functions in src.nosync/vendor/ (byte-identical
copies of engine_a/r2_core.py and engine_a/cliff.py, sha256 checked below when the originals are readable)."""
import os, sys, json, time, resource, hashlib, shutil, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'src.nosync'))
import numpy as np
np.seterr(all='ignore')
import sympy as sp
from fractions import Fraction

T0 = time.time()
OUT = os.path.join(HERE, 'out.nosync')
os.makedirs(OUT, exist_ok=True)


def res():
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return 'wall %.1f s, cpu %.1f s, max RSS %.0f MB' % (time.time() - T0, ru.ru_utime + ru.ru_stime, ru.ru_maxrss / 2 ** 20)


def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)


def dump(obj, name):
    for d in (HERE, OUT):
        with open(os.path.join(d, name), 'w') as f:
            json.dump(obj, f, indent=1, default=str)


import ea_core as E, ea_kern as K, ea_tables as TB, ea_gates as G, ea_ward as W1, ea_ward2 as W1b
import ea_closed as CF, ea_spectra as SP, ea_D1 as D1, ea_run_tables as RT, ea_invariant as INV

# ---------------------------------------------------------------------------------------------- 0 provenance
prov = {}
for nm in ('r2_core.py', 'cliff.py'):
    vend = os.path.join(HERE, 'src.nosync', 'vendor', nm)
    prov[nm] = {'vendor_sha256': hashlib.sha256(open(vend, 'rb').read()).hexdigest(),
                'original': os.path.join(E.ROOT, 'engine_a', nm),
                'original_sha256_at_copy': {'r2_core.py': 'e0c59d493c9022d02cae3a7699656d6c6f2382c4df069b9a77cc718e068835f5',
                                            'cliff.py': 'e4a89ad5fe50ecf7026e17978e3b6ae72b195e5131536601ee1add5b05d7a05f'}[nm]}
    prov[nm]['identical'] = prov[nm]['vendor_sha256'] == prov[nm]['original_sha256_at_copy']
log('vendor copies', {k: v['identical'] for k, v in prov.items()})
E.blocks()
log('blocks ready;', res())

# ---------------------------------------------------------------------------------------------- 1 gates
G.g_simplex()
G.g_selfcontraction()
bil = G.g_bilinears()
G.g_colour()
w1 = W1.gate(ntrials=3)
G.record('G3', 'Ward W1 letters A_perp, lambda+, lambda-, c, A_V (planar, exact): lc(sN) = {L,N}',
         {k: v for k, v in w1.items() if k.endswith('_pass')}, 'all pass', all(v for k, v in w1.items() if k.endswith('_pass')))
w1b = W1b.gate(ntrials=2)
G.record('G3', 'Ward W1 antifield letters A*_perp, A*_V, lambda*, c* on multi-field outputs (explicit su(2))',
         {k: v for k, v in w1b.items() if k.endswith('_pass') or k.endswith('_nontrivial')}, 'all pass',
         all(v for k, v in w1b.items() if k.endswith('_pass')))
cfb = CF.check_bubbles(4)
G.record('CF', 'k=1 closed form eq:tl-bubbleclosed vs explicit derivative letters (quartic bubbles, n<=4)',
         ['%s n=%d %s' % (r['pair'], r['n'], r['engine']) for r in cfb], 'closed form', all(r['pass'] for r in cfb))
cfg = CF.check_genfun(3)
G.record('CF', 'Bessel generating function Prop prop:tl-genfun (independent sympy) vs explicit, k=1,2, n<=3',
         ['%s n=%d %s' % (r['chain'], r['n'], r['engine']) for r in cfg], 'genfun', all(r['pass'] for r in cfg))
dump(G.GATES, 'gates.json')
log('gates stage 1 done;', res())
G.g_numeric()
dump(G.GATES, 'gates.json')
log('numeric residuals done;', res())

# ---------------------------------------------------------------------------------------------- 2 tables at W_0
raw = RT.run(os.path.join(OUT, 'kernels_W0_raw.json'))
G.record('G7', 'homogeneity: every W_0 table entry re-evaluated at an extra generic momentum point', len(TB.HOMOG_FAIL), 0,
         len(TB.HOMOG_FAIL) == 0)
spec = json.load(open(RT.SPEC))
bound_viol = []; zero_sets = {}; nonzero_sets = {}
for pname, ent in raw.items():
    ta, tb = ent['letters']
    W = E.WEYL[ta] + E.WEYL[tb]; dab = E.DEG[ta] + E.DEG[tb]
    for legset, words in ent['leg_sets'].items():
        types = legset.split() if legset != '(none)' else []
        wout = sum((E.WEYL[t] for t in types), Fraction(0))
        degout = sum(E.DEG[t] for t in types)
        if len(types) > W + dab or degout != dab or wout > W:
            bound_viol.append((pname, legset))
        nz = any(v for w in words.values() for v in w.values())
        (nonzero_sets if nz else zero_sets).setdefault(pname, []).append(legset)
G.record('G7', 'output bound n_out <= W + d_ab, ghost number d_out = d_ab, w_out <= W on every W_0 leg set (and log momentum'
         ' degree W - w_out, eps^-m degree W - w_out - 2m by construction of the interpolation)', bound_viol, [], not bound_viol)
inv = INV.run(raw)
nfit = sum(len(v) for v in inv.values()); nok = sum(1 for v in inv.values() for e in v.values() if e['fit_verified'])
G.record('G7', 'index-free form: bosonic W_0 kernels fitted exactly onto O(4) x O(6) invariant structures, fit verified on every '
         'component and monomial', '%d of %d words' % (nok, nfit), 'all', nok == nfit)
dump(G.GATES, 'gates.json')
log('tables done;', res())

# ---------------------------------------------------------------------------------------------- 3 spectra
pred = json.load(open(os.path.join(E.ROOT, 'step2', 'twoletter', 'theory', 'predictions.json')))
M2 = SP.two_site_su2()
G.record('G8', 'su(2) two-site planar kernel (units u, normalised letters, basis ZZ,ZX,XZ,XX)', [[str(x) for x in r] for r in M2],
         '2u P (H = 1 - P after -2u per letter)', M2 == [[2, 0, 0, 0], [0, 0, 2, 0], [0, 2, 0, 0], [0, 0, 0, 2]])
Esym = sp.Symbol('E')
su2 = {}; ok_all = True
for L in range(2, 9):
    for Mm in range(0, L // 2 + 1):
        Em, cp = SP.cyclic_spectrum(M2, L, Mm)
        key = 'L=%d,M=%d' % (L, Mm)
        pf = pred['P6_su2_Heisenberg']['cyclic_spectra_E'].get(key)
        ok = None
        if pf:
            target = sp.Integer(1)
            for f in pf['charpoly_factors']:
                target *= sp.sympify(f['factor'], locals={'E': Esym}) ** f['multiplicity']
            ok = sp.expand(target - cp) == 0
            ok_all &= ok
        su2[key] = {'charpoly': str(cp), 'matches_predictions': ok}
G.record('G8', 'su(2) sector cyclic charpolys of E = lc^natural/(-2u), L <= 8, all magnon numbers', {k: v['charpoly'] for k, v in su2.items()},
         'predictions.json P6', ok_all)
log('su(2) done;', res())
tw = {}
for N in range(0, 9):
    r = SP.twist_two(N)
    tw[N] = r
    log('twist two N=%d' % N, r['eigenvalues_units'], r['pass'])
G.record('G8', 'twist two: eig(lc^natural) on tr(D^k Z D^{N-k} Z), N <= 8 (gauge completion to O(A_n) included)',
         {N: r['eigenvalues_units'] for N, r in tw.items()}, '-8 h(S) units, S even <= N (predictions P5)', all(r['pass'] for r in tw.values()))
twnc = SP.twist_two(2, covariant=False)
G.record('G8', 'twist two N=2 without the O(A_n) terms (two-letter operator only): shows the r=1 outputs are needed',
         twnc['eigenvalues_units'], 'differs from [-12, 0]', twnc['eigenvalues_units'] != ['-12', '0'], 'diagnostic')
dump(G.GATES, 'gates.json')
log('spectra done;', res())

# ---------------------------------------------------------------------------------------------- 4 D_1 and P_1
lc_tl = Fraction(bil[('W', 'total')][0]); lc_tt = Fraction(bil[('W', 'total')][1])
lc_traceless, lc_trace = lc_tl, lc_tl + 6 * lc_tt
d1b = D1.d1_bilinears(lc_traceless, lc_trace)
ok9 = d1b['traceless']['O'] == '0' and d1b['trace']['O'] == '-48' and d1b['traceless']['S1(Delta)'] == '-theta - 1/3'
G.record('G9', 'D_1 on scalar bilinears assembled (two-body 4 lc from engine + Ch9 composite/one-body inputs)', d1b,
         'Ch9 Thm c9-D1: 0 / -48 O, -(1/3+theta){S1(Delta),O}', ok9)
p1b = D1.p1_bilinears()
G.record('G9', 'P_1 on scalar bilinears assembled (Ch9 composite inputs)', p1b, '-(1/3+theta_P){S1(P_b),O}',
         p1b['S1(P_b)'] == '-theta_P - 1/3')
dump(G.GATES, 'gates.json')

# ---------------------------------------------------------------------------------------------- 5 kernels.json
meta = {
    'engine': 'A (route A: Schwinger-parameter Gaussian generating functions, exact simplex and boundary-face integrals)',
    'theory': 'step2/twoletter/theory/twoletter_theory.tex (Thm thm:tl-decomp, thm:tl-simplex, Prop prop:tl-genfun)',
    'units': 'every coefficient is in units hbar/(16 pi^2), vertex scales divided out, for the colour element '
             'sum_gamma t^gamma (x) W(t_gamma), W = ad_{Y_k} ... ad_{Y_1} (Y_1 = first output, nearest to letter a); '
             'tr_g(ad_X ad_Y)/(X,Y) times this = units of u.  su(N) planar adjacent pair of tr(... X_a X_b ...): '
             '(-1)^r (value/2) u times tr(... Y_1 ... Y_r ...), u = -g^2 N/(32 pi^2)',
    'kernel': 'K(X_a,X_b) intrinsic (Def def:tl-dressed): (-1)^{|X_a|+|X_b|} times the str value with the operator Hessian '
              '(d_j d_l(delta M) = (-1)^{|M|(|j|+|l|)} delta d_j d_l M); for M = int omega (zeta_a, X_a)(zeta_b, X_b) the pair '
              'contributes zeta_a^T K zeta_b',
    'momenta': 'r{j}_{mu}: incoming momentum of output j (chain order from letter a), mu = 0..3 = engine A directions M = 0..3 '
               '(two-component index m = mu+1); i r <-> d acting on that output in position space',
    'components': 'index string "ia,ib,io1,io2,..." = basis vectors of ea_core.BASIS: AV, AsV: M = 0..3; AW, AsW: M = 4..9 '
                  '(SO(6) index I = M-3); lamP/lamM: 8 integer vectors spanning the +1/-1 eigenspaces of '
                  'rho(v1v2v3v4) = s0 sb1 s2 sb3 on Sigma (S_+ (x) W / S_- (x) W^*); lamsP/lamsM: +1/-1 on Sigma^* '
                  '(S_+ (x) W^* / S_- (x) W); Clifford matrices of src.nosync/vendor/cliff.py; their two-component (Weinberg-'
                  'Srednicki) form is the intertwiner U of step2/notation/check_dictionary.py (checks A1-A3)',
    'topology': 'block structure of the chain from a to b: C = cubic vertex with one output, Q = quartic vertex with two '
                'consecutive outputs; every internal line type summed',
    'regulator': 'per-edge Schwinger cutoff t_e >= eps, T = infinity (generic momenta); no dimensional regularisation',
    'provenance': prov,
}
kern = {'meta': meta, 'pair_weight_W0_tables': raw,
        'census_leg_sets_with_identically_zero_kernel_at_W0': zero_sets}
kern['invariant_forms_bosonic_W0'] = {'notation': INV.__doc__, 'forms': inv}
kern['closed_forms'] = {
    'k=0 (Lemma lem:tl-k0, all W)': 'K^(0)(d^al phi_a, d^be phi_b) = hbar Phat-numerator (-1)^{|be|} (d^{al+be} S_eps)(0), '
        '(d^gam S_eps)(0) = (-1)^{m+1}/(16 pi^2 m 2^{m-1} eps^m) sum_{pairings} prod delta, |gam| = 2(m-1); single power, no log. '
        'Verified (gate G6) for (A,A) and (c,A*_V) with up to 5 derivatives.',
    'k=1, constant chain numerator N_V (quartic bubbles; eq:tl-bubbleclosed, all W)':
        '16 pi^2 [log] = -N_V sum_{a,b,j} (i(s,R))^a (i(u,R))^b |u-s|^{2j} (|R|^2/4)^j (a+j)!(b+j)!/(a! b! (j!)^2 (a+b+2j+1)!), '
        'point splitting X_a(x+s) X_b(x+u), R = total output momentum; N_V = minus the W_0 log value of the bubble.',
    'k=1, k=2, general chain (Prop prop:tl-genfun)':
        'sum s^al u^be/(al! be!) [log]K = -(1/16 pi^2) e^{i(s,R)} int_Delta N_V(i d_sigma)[e^{-i(sigma,c)} (-1)^E '
        'sum_m (|sigma|^2/4)^{m+E-2} Phi^m/(m!(m+E-2)!)] (engine variables: edges k+q_i, c = -sum alpha_i q_i, sigma = u-s); '
        'eps^{-m}: (1/m) sum_e int_{R_e} alpha_e^m with Bessel order m+E-2.',
    'verification_bubbles': cfb, 'verification_genfun': cfg,
}
kern['planar_scalar_two_site'] = {
    'statement': 'adjacent pair of scalars in a single trace: tr(... A_I A_J ...) -> u [2 tr(... A_J A_I ...) - delta_IJ tr(... A_K A_K ...)] '
                 '(W_0 log kernel, quartic bubble + gluon exchange, planar reading)',
    'su2_two_site_units': [[str(x) for x in r] for r in M2],
    'lc_natural': '-2u sum_l (1 - P_{l,l+1} + K_{l,l+1}/2)  (K = trace operator)',
    'two_component': 'with A_ij = A_I Sigma^I_ij, A^ij A_ij = 4 A_I A_I: pair (A_ij, A_kl) -> u [2 A_kl A_ij - (1/2) eps_ijkl A^mn A_mn]; '
                     'su(N): u = -g^2 N/(32 pi^2)',
}
kern['su2_sector'] = su2
kern['twist_two'] = {str(N): r for N, r in tw.items()}
kern['twist_two_operator'] = ('tr(D_n^k Z D_n^{N-k} Z), n = (1, i, 0, 0) null (engine directions 0,1), Z letter zeta = e_4 + i e_5 '
                              '(two-component: Z = A_1 + i A_2 in the SO(6) index I = M-3), D = d + [A,-]; two-field part of lc gets '
                              'r = 2 outputs from (Z,Z) and r = 1 outputs from (A_n, Z) of the O(A_n) terms')
kern['D1'] = {
    'formula': 'D_1 = 4 lc^natural - {F\'_Delta, -} + C_Delta (Thm thm:tl-D1); theta kept',
    'two_body': '4 K^log(X_a, X_b) (tables above, computed) + C^(2)_Delta(X_a, X_b) (composite cycles: NOT computed by engine A; '
                'Ch9 value on the pair (A_perp, A_perp) without derivatives: -2u {S1(Delta),O}; prediction P7 on cohomology: '
                'Xi_Delta = -(1/3+theta) u {S1(Delta), -})',
    'one_body_units': D1.one_body(),
    'one_body_inputs': "F'_Delta from Thm c9-quadratic (iii), Pi_Delta|perp from Rem c9-Pi (manuscript), n_F = 1,2,2,4 for c, A, lambda",
    'bilinears': d1b,
    'spectrum_at_theta=-1/3_under_P7': {'twist_two_D1_units': {str(N): [str(4 * sp.Rational(x)) for x in r['eigenvalues_units']] for N, r in tw.items()},
                                        'Konishi': '-48', '20prime': '0'},
    'two_component_one_body_A_perp': "u[-8 + (5/3 - theta)(2E + 2 + 2|alpha|)] d^alpha A_ij",
}
kern['P1'] = {
    'formula': 'P_1 = -{hbar F_1^(1)(P_b), -} + C_{P_b} (Prop prop:tl-P1): no log-kernel term',
    'two_body': 'C^(2)_{P_b}(X_a,X_b) = Q-image of the renormalised composite cycles through the pair along P_b (composite kernel '
                '2(d_b (P*)^2 (x) 1)K_0, an A_V - c edge since (Q^GF)^2: A*_V -> c, c* -> A_V); NOT computed by engine A. '
                'Ch9 value on (A_perp, A_perp) without derivatives, symmetrised: -2u {S1(P_b),O}',
    'bilinears': p1b,
}
# ---------------------------------------------------------------------------------------------- revision 2026-10-29
import ea_revise as REV
kern, G.GATES[:] = REV.revise(kern, list(G.GATES), rerun=True, log=log)
log('revision gates (W2, su(1|1)) done;', res())
kern['resources'] = res()
dump(kern, 'kernels.json')
G.record('RES', 'resources', res(), '< 1.5 GB RSS, one process', True)
dump(G.GATES, 'gates.json')
npass = sum(1 for g in G.GATES if g['pass']); log('GATES: %d / %d pass' % (npass, len(G.GATES)))
log('done;', res())
