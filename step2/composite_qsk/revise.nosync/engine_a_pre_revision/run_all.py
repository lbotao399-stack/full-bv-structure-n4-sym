"""Engine A (route A), PART C: composite kernels C^(2)_Delta, C^(2)_{P_b} (special-edge composite cycles), D_1 and P_1 on
the families, Conjecture P7 (b), (c) and the spectra at theta = -1/3.  Runs everything and writes, here and in out.nosync/:
gates.json, results.json, PROGRESS.txt (appended), run_all_log.txt.  Code in src.nosync/ (imports the two-letter engine A
read-only from step2/twoletter/engine_a/src.nosync).  Single process, pure Python (fractions, numpy, sympy, mpmath).
Usage: python3 run_all.py [--quick]"""
import os, sys, json, time, resource, warnings
warnings.filterwarnings('ignore')
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'src.nosync'))
import numpy as np
np.seterr(all='ignore')
import sympy as sp
T0 = time.time()
OUT = os.path.join(HERE, 'out.nosync'); os.makedirs(OUT, exist_ok=True)
QUICK = '--quick' in sys.argv
LOGF = open(os.path.join(OUT, 'run_all_log.txt'), 'w')


def res():
    ru = resource.getrusage(resource.RUSAGE_SELF)
    return 'wall %.1f s, cpu %.1f s, max RSS %.0f MB' % (time.time() - T0, ru.ru_utime + ru.ru_stime, ru.ru_maxrss / 2 ** 20)


def log(*a):
    s = '[%7.1fs] ' % (time.time() - T0) + ' '.join(str(x) for x in a)
    print(s, flush=True); LOGF.write(s + '\n'); LOGF.flush()


def progress(msg):
    line = time.strftime('%Y-%m-%d %H:%M:%S') + ' run_all: ' + msg + '\n'
    for d in (HERE, OUT):
        with open(os.path.join(d, 'PROGRESS.txt'), 'a') as f:
            f.write(line)


def dump(obj, name):
    for d in (HERE, OUT):
        with open(os.path.join(d, name), 'w') as f:
            json.dump(obj, f, indent=1, default=str)


import cq_core as Q, cq_gval as GV, cq_gates as GT, cq_twist as TW, cq_su11 as S11, cq_fsl2 as FS, cq_gluon as GL
import cq_families as FM, cq_tables as TB
import ea_spectra as ES
Q.RULE['koszul'] = 'R2'; GV.INTR['on'] = True
progress('started%s' % (' (quick)' if QUICK else ''))
R = {'meta': {'engine': 'A (route A: Schwinger/simplex exact)', 'theory': 'step2/composite_qsk/composite/theory/composite_theory.tex',
              'spec': 'step2/composite_qsk/composite/theory/engine_spec.json', 'imports_read_only': Q.provenance(),
              'units': 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); kernel values in hbar/(16 pi^2) x colour word, planar factor (-1)^r/2',
              'conventions': {'koszul': 'R2: internal vertex and special edge odd relative to odd outputs only (not to each other)',
                              'intrinsic_pair_sign': '(-1)^{|X_a|+|X_b|} as ea_kern.kernel', 'global_constant': 1,
                              'vertex_vertex_kappa': 'symbolic: C^(2) = letter part + kappa * vertex-vertex part, kappa in {+1,-1,0} (GS passes for each)',
                              'log_basis': ['1', 'log2', 'log3', 'log5', 'log7']}},
     'families': {}}
rng = np.random.default_rng(2026)

# ------------------------------------------------------------------------------------------ 1 gates
log('stage 1: gates'); progress('stage 1 gates')
GT.g_weights(); GT.g_k1(); GT.g_cores(rng); GT.g_bilinear(); GT.g_mirror(rng); GT.g_position(rng)
dump(GT.GATES, 'gates.json'); log('gates stage 1 done;', res())

# ------------------------------------------------------------------------------------------ 2 F0
log('stage 2: F0 scalars'); progress('stage 2 F0')
so6 = FM.so6_density(log=log)
R['families']['F0_scalar_density'] = so6
GT.record('GC9-F0', 'P7(b) planar scalar sector: (A_perp, A_perp) density C_{P_beta} = -1/6 i beta.(q1+q2) (Y_i Y_j), no trace, no exchange (off shell)',
          so6[str(('Pv', (1, 2, -1, 3)))], '-1/6 u d_beta per adjacent pair', so6['pass_universal_-4/3_and_-1/6'])
GT.record('GC10-F0', 'P7(c) planar scalar sector: density C^0 = -4/3 (Y_i Y_j), no trace, no exchange; Xi^0 = 0 per letter at theta = -1/3',
          so6['D0'], '-4/3 u per adjacent pair = -2/3 u N', so6['pass_universal_-4/3_and_-1/6'])
pred = json.load(open(os.path.join(Q.ROOT, 'step2', 'twoletter', 'theory', 'predictions.json')))
M2 = ES.two_site_su2(); Esym = sp.Symbol('E'); su2 = {}; ok = True
for L in range(2, 9):
    for Mm in range(0, L // 2 + 1):
        Em, cp = ES.cyclic_spectrum(M2, L, Mm)          # E = lc^natural/(-2u) = D_1/(-8u) since Xi = 0 (density above)
        pf = pred['P6_su2_Heisenberg']['cyclic_spectra_E'].get('L=%d,M=%d' % (L, Mm))
        good = None
        if pf:
            target = sp.Integer(1)
            for f in pf['charpoly_factors']:
                target *= sp.sympify(f['factor'], locals={'E': Esym}) ** f['multiplicity']
            good = sp.expand(target - cp) == 0; ok &= good
        su2['L=%d,M=%d' % (L, Mm)] = {'charpoly_E': str(cp), 'pass': good}
R['families']['F0_su2'] = {'note': 'D_1 = 4 lc^natural + Xi_Delta with Xi_Delta density 0 (F0_scalar_density) and lc^natural from the two-letter engine; E = D_1/(-8u)', 'spectra': su2}
GT.record('GC11-F0', 'su(2) sector L <= 8: D_1/(-8u) charpolys at theta = -1/3 = Heisenberg (predictions P6)', 'all (L, M)', 'predictions.json P6', ok)
dump(GT.GATES, 'gates.json'); dump(R, 'results.json'); log('F0 done;', res())

# ------------------------------------------------------------------------------------------ 3 F2 su(1|1), fermion sl(2)
log('stage 3: F2'); progress('stage 3 F2')
dens = S11.densities(npts=3, log=log)
R['families']['F2_su11_densities'] = dens
ok = all(v['d0'] == '-4/3' and v['dP_c1'] == '-1/6' and v['dP_c2'] == '-1/6' for p_ in dens.values() for w, v in p_.items() if eval(w)[0] == eval(list(p_.keys())[0])[0] and True and False) if False else None
okd = True; rows = {}
for pr, d in dens.items():
    for w, v in d.items():
        same = (w == pr)
        tgt0 = '-4/3' if same else '0'; tgtc = '-1/6' if same else '0'
        g = v['d0'] == tgt0 and v['dP_c1'] == tgtc and v['dP_c2'] == tgtc and v['dP_fit_exact'] and v['d0_constant']
        okd &= g; rows['%s->%s' % (pr, w)] = [v['d0'], v['dP_c1'], v['dP_c2']]
GT.record('GC9/10-F2a', 'su(1|1) two-site densities on shell: same word d0 = -4/3, dP = -1/6 (q1+q2); exchanged word 0 (all ordered pairs of Z, psi)',
          rows, 'P7 (b), (c) per adjacent pair', okd)
su11 = FM.su11_spectra(Lmax=6, log=log)
R['families']['F2_su11_spectra'] = su11
GT.record('GC11-F2a', 'su(1|1) graded cyclic spectra L <= 6 of D_1/(-8u) at theta = -1/3 = harmonic action (ha_su11.json); lc^natural density -2(1-Pi), Xi^0 density 0',
          {k: v['E_H_from_D1'] for k, v in su11['spectra'].items()}, 'ha_su11.json', su11['pass'])
fs = {}
for N in ([1, 2] if QUICK else [1, 2, 3]):
    fs[N] = FS.run(N, n=(1, -1j, 0, 0), log=log)
    R['families']['F2_fermion_sl2'] = fs; dump(R, 'results.json')
GT.record('GC10-F2b', 'fermion sl(2): sum C^(2),0_Delta = -2/3 (4+N) on the two-field part (on shell)', {N: r['M0']['matrix'] for N, r in fs.items()}, '-2/3(4+N) identity',
          all(sp.Matrix(sp.sympify(r['M0']['matrix'])) == -sp.Rational(2, 3) * (4 + N) * sp.eye(len(r['M0']['matrix'])) for N, r in fs.items()))
GT.record('GC9-F2b', 'fermion sl(2): sum C^(2)_{P_n} = -1/3 u d_n on the two-field part (on shell, exact)', {N: r['Pn'] for N, r in fs.items()}, 'residual 0',
          all(r['Pn']['pass'] for r in fs.values()))
GT.record('GC11-F2b', 'fermion sl(2): D_1 eigenvalues at theta = -1/3 = -32 h(n+1), n odd (log kernel eig -8h(n+1) from the general core)', {N: r['D1'] for N, r in fs.items()}, '-8 u E_H',
          all(r['D1']['pass'] for r in fs.values()))
dump(GT.GATES, 'gates.json'); dump(R, 'results.json'); log('F2 done;', res())

# ------------------------------------------------------------------------------------------ 4 F1 twist two
log('stage 4: F1 twist two'); progress('stage 4 F1')
tw = {}; ok9 = True; ok10 = True; ok11 = True
Ns = [0, 1, 2, 3] if QUICK else [0, 1, 2, 3, 4]
for N in Ns:
    r = TW.run(N, kinds=('D0', 'Pn'), log=log)
    M0 = sp.Matrix([[sp.sympify(x) for x in row] for row in r['D0']['matrix_by_logbasis']['1']])
    logs0 = all(all(sp.sympify(x) == 0 for row in r['D0']['matrix_by_logbasis'][b] for x in row) for b in ['log2', 'log3', 'log5', 'log7'])
    g10 = logs0 and M0 == -sp.Rational(2, 3) * (4 + N) * sp.eye(M0.shape[0])
    # P_n: image of O_k must be -1/3 (a1 + a2) B_k, a_i = i n.q_i; B_k = 4 (a1^k a2^{N-k} + a1^{N-k} a2^k) (zZ.eZ = 2)
    Pn = sp.Matrix([[sp.sympify(x) for x in row] for row in r['Pn']['matrix_by_logbasis']['1']])
    a1, a2 = sp.symbols('a1 a2')
    g9 = True
    for ik, k in enumerate(range(N // 2 + 1)):
        img = sum(Pn[a, ik] * a1 ** a * a2 ** (N + 1 - a) for a in range(N + 2))
        Bk = 4 * (a1 ** k * a2 ** (N - k) + a1 ** (N - k) * a2 ** k) if 2 * k != N else 8 * a1 ** k * a2 ** k
        g9 &= sp.expand(img - sp.Rational(-1, 3) * (a1 + a2) * Bk) == 0
    lc = ES.twist_two(N)
    Lm = sp.Matrix([[sp.sympify(x) for x in row] for row in lc['matrix_units']])
    D1 = 4 * Lm + M0 + sp.Rational(2, 3) * (4 + N) * sp.eye(M0.shape[0])
    ev = sorted([sp.nsimplify(e) for e, m in D1.eigenvals().items() for _ in range(m)], key=float)
    target = sorted([-32 * sum(sp.Rational(1, j) for j in range(1, S + 1)) for S in range(0, N + 1, 2)], key=float)
    g11 = ev == target
    r['D1'] = {'matrix': [[str(D1[i, j]) for j in range(D1.shape[1])] for i in range(D1.shape[0])], 'eigenvalues': [str(x) for x in ev],
               'target_-32h(S)': [str(x) for x in target], 'lc_natural_matrix_twoletter': lc['matrix_units'], 'pass': g11}
    r['P7c_pass'] = bool(g10); r['P7b_n_pass'] = bool(g9)
    ok9 &= g9; ok10 &= g10; ok11 &= g11
    tw[N] = r
    R['families']['F1_twist_two'] = tw; dump(R, 'results.json')
    log('twist two N=%d: P7(c) %s P7(b, n) %s D1 eig %s target %s' % (N, g10, g9, ev, target))
pg = {}
for N in ([1] if QUICK else [1, 2]):
    pg[N] = TW.run(N, kinds=('Pg',), log=log)['Pg']
R['families']['F1_twist_two_generic_beta'] = {'ansatz': '(i beta.q_i)(i n.q1)^a (i n.q2)^{N-a} [rows 0..2N+1], (beta.n)(q1.q2, q1^2, q2^2)(i n.q1)^a (i n.q2)^{N-1-a}', 'by_N': pg}
GT.record('GC10-F1', 'twist two: sum C^(2),0_Delta = -2/3 (4+N) on the two-field part, N <= %d (exact, rational)' % max(Ns), {N: tw[N]['D0']['matrix_by_logbasis']['1'] for N in Ns}, '-2/3(4+N)', ok10)
GT.record('GC9-F1', 'twist two: sum C^(2)_{P_n} = -1/3 u d_n exactly (two-field part), N <= %d; generic beta: = -1/3 u d_beta + (q_i^2 terms, EOM)' % max(Ns),
          {'P_n': ok9, 'generic_beta_matrices': {N: v['matrix_by_logbasis']['1'] for N, v in pg.items()}}, '-1/3 u d_b on cohomology', ok9)
GT.record('GC11-F1', 'twist two: D_1 = 4 lc^natural + sum C^(2)_Delta - sum{F_Delta, X} at theta = -1/3: eigenvalues -32 h(S), S even <= N', {N: tw[N]['D1']['eigenvalues'] for N in Ns}, '-8u E_H, E_H = 4h(S)', ok11)
dump(GT.GATES, 'gates.json'); dump(R, 'results.json'); log('F1 done;', res())

# ------------------------------------------------------------------------------------------ 5 on-shell universality of densities
log('stage 5: universality'); progress('stage 5 universality')
uni = FM.universality(npts=2, log=log)
R['families']['density_universality_on_shell'] = uni
GT.record('GU', 'on-shell planar densities, general polarisations (letter vectors of definite chirality): same word -4/3 / -1/6 i beta.(q1+q2), exchanged word 0, for pairs of A_perp, lambda_+, lambda_-',
          {k: v['pass'] for k, v in uni.items()}, 'all', all(v['pass'] for v in uni.values()))
dump(GT.GATES, 'gates.json'); dump(R, 'results.json'); log('universality done;', res())

# ------------------------------------------------------------------------------------------ 6 F3 gluon
log('stage 6: F3 gluon'); progress('stage 6 F3')
gl = {}
for N in ([0] if QUICK else [0, 1]):
    gl[N] = GL.run(N, log=log)
    R['families']['F3_gluon'] = gl; dump(R, 'results.json')
GT.record('GC10-F3', 'gluon sl(2): sum C^(2),0_Delta = -2/3 (6+N) on the two-field part (letter part; vertex-vertex part on shell)',
          {N: {'letter': r['M0']['letter_matrix'], 'vv': r['M0']['vv_matrix']} for N, r in gl.items()}, '-2/3(6+N)',
          all(sp.Matrix(sp.sympify(r['M0']['letter_matrix'])) == -sp.Rational(2, 3) * (6 + N) * sp.eye(len(r['M0']['letter_matrix'])) for N, r in gl.items()))
GT.record('GC9-F3', 'gluon sl(2): sum C^(2)_{P_n} = -1/3 u d_n on the two-field part (on shell)', {N: {'letter_zero': r['Pn']['letter_zero'], 'vv_zero': r['Pn']['vv_zero']} for N, r in gl.items()},
          'residual 0', all(r['Pn']['letter_zero'] for r in gl.values()))
GT.record('GC11-F3', 'gluon sl(2): D_1 eigenvalues at theta = -1/3 = -32 h(n+2), n even (all kappa)', {N: r['D1'] for N, r in gl.items()}, '-8u E_H',
          all(r['D1']['kappa=1']['pass'] and r['D1']['kappa=-1']['pass'] for r in gl.values()))
GT.record('GW-F3', 'gauge Ward identity: two-field part with one longitudinal leg vanishes on shell (letter and vertex-vertex parts separately)',
          {N: {k: [v['letter_part_zero'], v['vv_part_zero']] for k, v in r['ward'].items()} for N, r in gl.items()}, 'zero',
          all(all(v['letter_part_zero'] and v['vv_part_zero'] for v in r['ward'].values()) for r in gl.values()))
dump(GT.GATES, 'gates.json'); dump(R, 'results.json'); log('F3 done;', res())

# ------------------------------------------------------------------------------------------ 7 tables
log('stage 7: W0 tables'); progress('stage 7 tables')
tab = TB.run(os.path.join(Q.ROOT, 'step2', 'composite_qsk', 'composite', 'theory', 'src.nosync', 'census_sp.json'), npts=1 if QUICK else 2, log=log)
dump(tab, 'kernels_C2.json')
R['tables'] = {'file': 'kernels_C2.json', 'pairs': {p: {ls: v['nonzero'] for ls, v in d['leg_sets'].items()} for p, d in tab.items()}}
dump(R, 'results.json')

# ------------------------------------------------------------------------------------------ 8 D_1, P_1 with symbolic free directions
import cq_symbolic as SY
R['D1_P1_symbolic'] = SY.run()
dump(R, 'results.json')

# ------------------------------------------------------------------------------------------ summary
npass = sum(1 for g in GT.GATES if g['pass']); ntot = len(GT.GATES)
R['summary'] = {'gates_pass': npass, 'gates_total': ntot, 'failed': [g['id'] for g in GT.GATES if not g['pass']], 'resources': res()}
dump(GT.GATES, 'gates.json'); dump(R, 'results.json')
log('ALL DONE: %d/%d gates;' % (npass, ntot), res())
progress('finished: %d/%d gates; %s' % (npass, ntot, res()))
