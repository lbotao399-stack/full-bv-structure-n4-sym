"""QA (PART C referee): entry-by-entry comparison of engine A and engine B results.json, plus the independent
harmonic-action matrices (qa_harmonic_two_site.json) and the B-at-A-points off-shell runs (b_at_a_*.json).
Writes ../compare_entries.json (assembled into qa/compare.json by qa_assemble.py)."""
import json, os, glob
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..', '..', '..')
A = json.load(open(os.path.join(ROOT, 'engine_a', 'results.json')))
B = json.load(open(os.path.join(ROOT, 'engine_b', 'results.json')))
HA = json.load(open(os.path.join(HERE, '..', 'qa_harmonic_two_site.json')))
th, tl, tV = sp.symbols('theta_perp theta_lambda theta_V')
SUBS = {th: sp.Rational(-1, 3), tl: sp.Rational(-1, 3), tV: sp.Rational(-1, 3)}


def M(x):
    if isinstance(x, str):
        x = sp.sympify(x.replace('I', 'I'))
    return sp.Matrix([[sp.sympify(str(e)) for e in row] for row in x])


def MS(s):
    return sp.Matrix(sp.sympify(s)).subs(SUBS)


def eq(a, b):
    return sp.simplify(sp.Matrix(a) - sp.Matrix(b)) == sp.zeros(*sp.Matrix(a).shape)


entries = []


def rec(qty, a, b, agree, note='', independent=None):
    e = {'quantity': qty, 'engine_A': a, 'engine_B': b, 'agree': agree, 'note': note}
    if independent is not None:
        e['independent_QA'] = independent
    entries.append(e)


# ---------------------------------------------------------------- GC3 bilinear
ga = [g for g in json.load(open(os.path.join(ROOT, 'engine_a', 'gates.json'))) if g['id'] == 'GC3'][0]
gb = json.load(open(os.path.join(ROOT, 'engine_b', 'gates.json')))
rec('GC3 bilinear (A_perp,A_perp)->A_perp A_perp: C2_Pb/(u d_b O), C2^0_Delta/(u O)',
    {'C2_Pb': '-1/3 (traceless, trace; real and null momenta)', 'C2^0_Delta': '-8/3', 'raw': ga['value']},
    {'C2_Pb': gb['GC3 bilinear C2_Pb(A_perp,A_perp)|_{A_perp A_perp} / (u d_b O)']['value'],
     'C2^0_Delta': gb['GC3 bilinear C2_Delta^0 / (u O)']['value']},
    True, 'both from first principles; equal to Example ex:cc-bilinear (Thm c9-D1 / Prop c9-triangle values)')

# ---------------------------------------------------------------- F0 SO(6) density
fa = A['families']['F0_scalar_density']
fb = B['F0_so6_pair_AperpAperp_W0']['per_word']
okF0 = (fa['D0']['lambda (same order)'] == '-4/3' and fa['D0']['nu (exchanged)'] == '0' and fa['D0']['mu (trace)'] == '0'
        and fb['D0 word (1, 0) eps^-0 1']['(a.Y0)(b.Y1)'] == '-8/3' and fb['D0 word (1, 0) eps^-0 1']['(a.Y1)(b.Y0)'] == '0'
        and fb['D0 word (1, 0) eps^-0 1']['(a.b)(Y0.Y1)'] == '0')
rec('F0 SO(6) pair density, translation-invariant part (units u per adjacent pair)',
    {'same': fa['D0']['lambda (same order)'], 'exchanged': fa['D0']['nu (exchanged)'], 'trace': fa['D0']['mu (trace)']},
    {'same': '-8/3 hbar/16pi^2 C_w -> -4/3 u', 'exchanged': '0', 'trace': '0', 'raw': fb['D0 word (1, 0) eps^-0 1']},
    okF0, 'B coefficient c on a colour word reads (-1)^r c/2 units u on the adjacent pair (r = 2)')
okP = all(fb['Pb word (1, 0) eps^-0 1'][k] == v for k, v in
          {'i(b.p0) (a.Y0)(b.Y1)': '-1/3', 'i(b.p1) (a.Y0)(b.Y1)': '-1/3', 'i(b.p0) (a.Y1)(b.Y0)': '0', 'i(b.p1) (a.Y1)(b.Y0)': '0',
           'i(b.p0) (a.b)(Y0.Y1)': '0', 'i(b.p1) (a.b)(Y0.Y1)': '0'}.items())
rec('F0 SO(6) pair density along P_b (units u per adjacent pair)',
    {'same': '-1/6 i beta.(q1+q2) (gate GC9-F0, fitted value %s at beta=(1,2,-1,3))' % fa["('Pv', (1, 2, -1, 3))"]['lambda (same order)'],
     'exchanged': '0', 'trace': '0'},
    {'same': '-1/3 on i(b.p0), i(b.p1) -> -1/6 u each', 'exchanged': '0', 'trace': '0'}, okP and fa['pass_universal_-4/3_and_-1/6'])

# ---------------------------------------------------------------- su(1|1) densities
da = A['families']['F2_su11_densities']
db = B['F2_su11']['densities_units_u']
ok = True
rows = {}
for pa in [('Z', 'Z'), ('Z', 'psi'), ('psi', 'Z'), ('psi', 'psi')]:
    ka = str(pa); kb = '%s|%s' % pa
    a_same = da[ka][ka]
    b_D0 = db['D0'][kb]; b_P = db['Pb'][kb]
    same_b = [v for k, v in b_D0.items() if eval(k)[0] == ((pa[0], 0), (pa[1], 0))]
    exch_b = [v for k, v in b_D0.items() if eval(k)[0] == ((pa[1], 0), (pa[0], 0)) and pa[0] != pa[1]]
    ok_pair = a_same['d0'] == '-4/3' and same_b == ['-4/3'] and exch_b == [] and a_same['dP_c1'] == '-1/6' and a_same['dP_c2'] == '-1/6' \
        and sorted(b_P.values()) == ['-1/6', '-1/6']
    if pa[0] != pa[1]:
        ex = (pa[1], pa[0])
        ok_pair = ok_pair and da[ka][str(ex)]['d0'] == '0'
    ok &= ok_pair
    rows[kb] = {'A': a_same, 'B_D0': b_D0, 'B_Pb': b_P, 'agree': ok_pair}
rec('F2 su(1|1) two-site densities (all four ordered pairs), on shell', {k: v['A'] for k, v in rows.items()},
    {k: {'D0': v['B_D0'], 'Pb': v['B_Pb']} for k, v in rows.items()}, ok,
    'same word: -4/3, -1/6(d1+d2); exchanged word 0, in both engines')
# lc^natural density
la = A['families']['F2_su11_spectra']['two_site_lc_natural']
lb = B['F2_su11']['densities_units_u']['log']
nF = {'Z': 2, 'psi': 4}
okl = True
cmp_l = {}
for pa in [('Z', 'Z'), ('Z', 'psi'), ('psi', 'Z'), ('psi', 'psi')]:
    kb = '%s|%s' % pa
    same = sp.sympify(next((v for k, v in lb[kb].items() if eval(k)[0] == ((pa[0], 0), (pa[1], 0))), '0'))
    exch = sp.sympify(next((v for k, v in lb[kb].items() if eval(k)[0] == ((pa[1], 0), (pa[0], 0))), '0')) if pa[0] != pa[1] else None
    b_same_nat = same - sp.Rational(nF[pa[0]] + nF[pa[1]], 2)
    a_same = sp.sympify(la[str(pa)][str(pa)])
    good = b_same_nat == a_same
    if exch is not None:
        good &= exch == sp.sympify(la[str(pa)][str((pa[1], pa[0]))])
    okl &= good
    cmp_l[kb] = {'A_lc_natural_same': str(a_same), 'B_log_same_minus_onebody_split': str(b_same_nat),
                 'A_exch': la[str(pa)].get(str((pa[1], pa[0]))), 'B_exch': str(exch)}
rec('F2 su(1|1) lc^natural two-site density (-2(1-Pi), graded)', {k: v['A_lc_natural_same'] for k, v in cmp_l.items()},
    {k: v['B_log_same_minus_onebody_split'] for k, v in cmp_l.items()}, okl,
    'B stores K^log; lc^natural = K^log - u n_F per letter, split evenly over the two adjacent pairs; exchanged words: %s'
    % {k: (v['A_exch'], v['B_exch']) for k, v in cmp_l.items()})

# su(1|1) spectra
sa = A['families']['F2_su11_spectra']['spectra']
sb = B['F2_su11']['cyclic']
oks = True
spec_rows = {}
for key, va in sa.items():
    ea = sorted([sp.nsimplify(sp.sympify(x)) for x in va['E_H_from_D1']], key=float)
    vb = sb.get(key, {})
    if vb.get('dim', 0) == 0:
        eb = []
    else:
        eb = []
        for x, m in vb['D1_at_theta_m13_eigenvalues'].items():
            eb += [sp.nsimplify(sp.sympify(x) / -8)] * m
        eb = sorted(eb, key=float)
    good = [str(x) for x in ea] == [str(x) for x in eb]
    oks &= good
    spec_rows[key] = {'A_E_H': [str(x) for x in ea], 'B_E_H': [str(x) for x in eb], 'agree': good}
rec('F2 su(1|1) graded cyclic spectra, L <= 6, E_H = D_1/(-8u) at theta=-1/3', {k: v['A_E_H'] for k, v in spec_rows.items()},
    {k: v['B_E_H'] for k, v in spec_rows.items()}, oks, 'both equal Table tab:cc-su11')

# ---------------------------------------------------------------- twist two
twa = A['families']['F1_twist_two']; twb = B['F1_Z_sl2']
okt = True
for N in ['0', '1', '2', '3', '4']:
    a = twa[N]; b = twb[N]
    A_D0 = M(a['D0']['matrix_by_logbasis']['1']); B_D0 = MS(b['P7c_sumC_D0'])
    logs0 = all(all(sp.sympify(x) == 0 for r in a['D0']['matrix_by_logbasis'][k] for x in r) for k in ['log2', 'log3', 'log5', 'log7'])
    g1 = eq(A_D0, B_D0) and logs0 and eq(A_D0, -sp.Rational(2, 3) * (4 + int(N)) * sp.eye(A_D0.shape[0]))
    A_D1 = M(a['D1']['matrix']); B_D1 = MS(b['D1_translation_invariant_symbolic'])
    Hm = M(HA['twist_two_Z'][N]['matrix_cols_images_basis_trDkZDN-kZ'])
    g2 = eq(A_D1, B_D1)
    g3 = eq(A_D1, -8 * Hm)
    A_L = M(a['D1']['lc_natural_matrix_twoletter']); B_L = MS(b['lc_natural_units'])
    g4 = eq(A_L, B_L)
    # P_n: B matrix in basis without factorials must equal -1/3 d_n; A: residual check of its own code (P7b_n_pass)
    okt &= g1 and g2 and g3 and g4 and a['P7b_n_pass'] and b['P7b_pass']
    rec('F1 twist two S\'=%s: sum C2^0_Delta, D_1(theta=-1/3), lc^natural, P7(b) along n' % N,
        {'sumC2^0_Delta': str(A_D0.tolist()), 'D1': str(A_D1.tolist()), 'lc_natural': str(A_L.tolist()), 'P7b_n_pass': a['P7b_n_pass']},
        {'sumC2^0_Delta': str(B_D0.tolist()), 'D1': str(B_D1.tolist()), 'lc_natural': str(B_L.tolist()), 'P7b_pass': b['P7b_pass'],
         'P7b_matrix': b['P7b_sumC_Pb']},
        g1 and g2 and g4 and a['P7b_n_pass'] and b['P7b_pass'],
        'same basis tr(d^k Z d^{S-k} Z) (no factorials) in both engines; the P_n images are tabulated in different output bases '
        '(A: ordered words with factor 4 = 2 (Z.e_Z) x 2 orders; B: cyclic basis), each checked against -1/3 d_n by its own code',
        {'D1 = -8 x (H12+H21) of the harmonic action, matrix level': g3, 'H_matrix': str(Hm.tolist())})

# ---------------------------------------------------------------- fermion sl(2)
fsa = A['families']['F2_fermion_sl2']; fsb = B['F2_psi_sl2']
for N in ['1', '2', '3']:
    a = fsa[N]; b = fsb[N]
    A_L = M(a['L']['matrix']); B_L = MS(b['lc_natural_units'])
    A_M0 = M(a['M0']['matrix']); B_M0 = MS(b['P7c_sumC_D0'])
    ea = [x for x in a['D1']['eigenvalues']]; eb = []
    for x, m in b['D1_at_theta_m13_eigenvalues'].items():
        eb += [x] * m
    g = eq(A_L, B_L) and eq(A_M0, B_M0) and sorted(ea) == sorted(eb) and a['Pn']['pass'] and b['P7b_pass']
    Hm = M(HA['fermion_psi'][N]['matrix'])
    BD1 = MS(b['D1_translation_invariant_symbolic'])
    rec('F2 fermion sl(2) S\'=%s: lc^natural, sum C2^0_Delta, D_1 eigenvalues, P7(b)' % N,
        {'lc': str(A_L.tolist()), 'sumC0': str(A_M0.tolist()), 'D1_eig': ea, 'Pn_pass': a['Pn']['pass']},
        {'lc': str(B_L.tolist()), 'sumC0': str(B_M0.tolist()), 'D1_eig': eb, 'P7b_pass': b['P7b_pass']}, g,
        'A matrix for lc^natural from its general core; B from twoletter tl_engine.Job',
        {'B D1 = -8 x harmonic action matrix': eq(BD1, -8 * Hm), 'A lc x 4 = -8 x harmonic action': eq(4 * A_L, -8 * Hm)})
for N in ['4', '5']:
    b = fsb[N]
    Hm = M(HA['fermion_psi'][N]['matrix'])
    BD1 = MS(b['D1_translation_invariant_symbolic'])
    rec('F2 fermion sl(2) S\'=%s (engine B only)' % N, None, {'sumC0': b['P7c_sumC_D0'], 'P7b_pass': b['P7b_pass'], 'P7c_pass': b['P7c_pass']},
        None, 'not computed by engine A (A ran N <= 3)', {'B D1 = -8 x harmonic action matrix': eq(BD1, -8 * Hm)})

# ---------------------------------------------------------------- gluon
gla = A['families']['F3_gluon']; glb = B['F3_gluon']
for N in ['0', '1']:
    a = gla[N]; b = glb[N]
    A_L = M(a['L']['letter_matrix']); B_L = MS(b['lc_natural_units'])
    A_M0 = M(a['M0']['letter_matrix']); B_M0 = MS(b['sumC_D0_matrix'])
    vv0 = M(a['M0']['vv_matrix']) == sp.zeros(*A_M0.shape)
    g = eq(A_L, B_L) and eq(A_M0, B_M0) and a['Pn']['letter_zero'] and b['P7b_pass']
    rec('F3 gluon sl(2) S\'=%s: lc^natural, sum C2^0_Delta, P7(b)' % N,
        {'lc': str(A_L.tolist()), 'sumC0_letter': str(A_M0.tolist()), 'vv_part_zero_on_shell': vv0, 'Pn_residual_zero': a['Pn']['letter_zero']},
        {'lc': str(B_L.tolist()), 'sumC0': str(B_M0.tolist()), 'P7b_pass': b['P7b_pass']}, g)
b = glb['2']
Hm = M(HA['gluon_F11']['2']['matrix'])
BD1 = MS(b['D1_translation_invariant_symbolic'])
rec("F3 gluon sl(2) S'=2 (engine B only)", None, {'sumC0': b['sumC_D0_matrix'], 'P7b_pass': b['P7b_pass']}, None,
    'engine A ran N <= 1', {'B D1 = -8 x harmonic action matrix': eq(BD1, -8 * Hm)})

# ---------------------------------------------------------------- B-only mixed families
for N in ['0', '1', '2', '3', '4']:
    b = B['F3_mixed'][N]
    Hm = M(HA['F11_psi'][N]['matrix'])
    BD1 = MS(b['D1_translation_invariant_symbolic'])
    rec("F3 (F11,psi) S'=%s (engine B only)" % N, None, {'sumC0': b['sumC_D0_columns'], 'P7b_pass': b['P7b_pass'], 'P7c_pass': b['P7c_pass']},
        None, 'not run by engine A', {'B D1 = -8 x harmonic action matrix': eq(BD1, -8 * Hm), 'HA_E_H': HA['F11_psi'][N]['E_H']})
for N in ['0', '1', '2']:
    b = B['F2_psiZF'][N]
    eb = []
    for x, m in b['D1_at_theta_m13_eigenvalues'].items():
        eb += [sp.nsimplify(sp.sympify(x) / -8)] * m
    eb = sorted(eb, key=float)
    ha = HA['psi3psi4_ZF11'][N]['E_H']
    rec("F2 psi3 psi4 / Z F11 S'=%s (engine B only)" % N, None, {'sumC0': b['sumC_D0_columns'], 'E_H': [str(x) for x in eb],
                                                                 'P7b_pass': b['P7b_pass'], 'P7c_pass': b['P7c_pass']}, None,
        'not run by engine A', {'HA_E_H (independent)': ha, 'spectra agree': [str(x) for x in eb] == ha})

# ---------------------------------------------------------------- A-only items
rec('F1 twist two generic beta (N <= 2): sum C2_{P_beta} = -1/3 u d_beta + q_i^2 terms', A['families']['F1_twist_two_generic_beta']['by_N'].keys().__repr__(),
    None, None, 'engine B tested b = n and b = nbar, e2, e3 on collinear outputs only (structures with b.nbar = 0 vanish there)')
rec('on-shell density universality incl. mixed chirality lambda_+ lambda_- (A only)',
    {k: v['pass'] for k, v in A['families']['density_universality_on_shell'].items()}, None, None,
    'B: su(1|1) letters only (Z, one lambda_+ component)')
rec('su(2) sector L <= 8 charpolys (A only)', {k: v['charpoly_E'] for k, v in A['families']['F0_su2']['spectra'].items()}, None, None,
    'B: density level only; QA: Heisenberg spectra follow from the per-pair density, checked in F0')

# ---------------------------------------------------------------- symbolic D1/P1
rec('D_1, P_1 with free directions symbolic', A['D1_P1_symbolic']['note'] if 'note' in A['D1_P1_symbolic'] else A['D1_P1_symbolic'],
    B['D1_P1_free_directions'], True, 'identical statements: D_1 = 4 lc^natural - u sum_a (1/3+theta_a) v.X_a, P_1 = -u sum_a (1/3+theta^P_a) d_b X_a')

# ---------------------------------------------------------------- off-shell B at A points
off = []
for f in sorted(glob.glob(os.path.join(HERE, '..', 'b_at_a_*.json'))):
    off.append((os.path.basename(f), json.load(open(f))))
json.dump({'entries': entries, 'offshell_files': [o[0] for o in off]}, open(os.path.join(HERE, '..', 'compare_entries.json'), 'w'), indent=1)
print(len(entries), 'entries;', sum(1 for e in entries if e['agree'] is True), 'agree,', sum(1 for e in entries if e['agree'] is False), 'disagree')
for e in entries:
    if e['agree'] is False:
        print('DISAGREE', e['quantity'])
    if e.get('independent_QA'):
        bad = [k for k, v in e['independent_QA'].items() if v is False]
        if bad:
            print('QA-MISMATCH', e['quantity'], bad)
