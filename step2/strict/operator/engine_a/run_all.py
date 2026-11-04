"""Engine A (route A), operator side of the strict Costello route (R10): one script for every task of
step2/strict/operator/theory/engine_spec.json that route A computes.

Run:  python3 run_all.py [--quick] [--only STAGE[,STAGE..]] [--force]
      single process, pure Python (fractions, numpy, sympy, mpmath); code in src.nosync/, caches in out.nosync/;
      resumable: every stage writes out.nosync/stage_<name>.json and every table chunk out.nosync/tables/*.json;
      a watchdog thread aborts above 1.45 GB RSS; prints resources at the end.
Writes results.json, gates.json (here and in out.nosync/), appends to PROGRESS.txt."""
import os, sys, time, json, threading, resource, glob
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'src.nosync')
OUT = os.path.join(HERE, 'out.nosync')
sys.path.insert(0, SRC)
os.makedirs(OUT, exist_ok=True)
import numpy as np
import sympy as sp
import op_common as C

QUICK = '--quick' in sys.argv
FORCE = '--force' in sys.argv
ONLY = None
for i, a in enumerate(sys.argv):
    if a == '--only' and i + 1 < len(sys.argv):
        ONLY = set(sys.argv[i + 1].split(','))
T_START = time.time()
LIMIT_MB = 1450


def log(msg):
    line = '%s %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(line, flush=True)
    with open(os.path.join(OUT, 'run_all_log.txt'), 'a') as f:
        f.write(line + '\n')


def progress(msg):
    with open(os.path.join(HERE, 'PROGRESS.txt'), 'a') as f:
        f.write('%s run_all: %s\n' % (time.strftime('%Y-%m-%d %H:%M:%S'), msg))


def watchdog():
    while True:
        time.sleep(10)
        r = C.rss_mb()
        if r > LIMIT_MB:
            progress('watchdog: RSS %.0f MB above %d MB, aborting' % (r, LIMIT_MB))
            os._exit(3)


threading.Thread(target=watchdog, daemon=True).start()


def jdump(obj, path):
    tmp = path + '.tmp'
    json.dump(obj, open(tmp, 'w'), indent=1, default=str)
    os.replace(tmp, path)


def stage(name, fn):
    path = os.path.join(OUT, 'stage_%s.json' % name)
    if os.path.exists(path) and not FORCE and (ONLY is None or name not in ONLY):
        try:
            d = json.load(open(path))
            log('stage %s: cached' % name)
            return d
        except Exception:
            pass
    if ONLY is not None and name not in ONLY and not os.path.exists(path):
        log('stage %s: skipped (--only)' % name)
        return None
    t0 = time.time()
    log('stage %s: start' % name)
    res, gates = fn()
    d = {'result': res, 'gates': [{'gate': g[0], 'value': g[1], 'target': g[2], 'pass': bool(g[3])} for g in gates],
         'seconds': round(time.time() - t0, 1), 'max_rss_MB': round(C.maxrss_mb(), 1)}
    jdump(d, path)
    msg = 'stage %s: done %.1f s, %s' % (name, time.time() - t0, ', '.join('%s:%s' % (g['gate'][:60], 'pass' if g['pass'] else 'FAIL') for g in d['gates']))
    log(msg); progress(msg[:600])
    return d


# ------------------------------------------------------------------------------------------------ stages
def st_OP0():
    import op_ops as O
    res, gates = O.run(log=lambda s: None)
    res = {k: (v if not isinstance(v, (np.bool_,)) else bool(v)) for k, v in res.items()}
    return res, [(g[0], g[1], True, g[1]) for g in gates]


TABLE_JOBS = [(g, 0, 'all') for g in ('P', 'D0', 'M0', 'R', 'xi', 'xibar', 'eta0', 'etabar0', 'K0')]
TABLE_JOBS_EXT = [(g, 1, 'all') for g in ('P', 'D0', 'R', 'xi', 'xibar', 'M0', 'eta0', 'etabar0', 'K0')]
TABLE_JOBS_EXT2 = [(g, 2, 'field') for g in ('R', 'xi', 'xibar', 'P', 'D0', 'M0')] + [(g, 3, 'field') for g in ('P', 'D0')]
TABLE_JOBS_EXT3 = [(g, 2, 'field') for g in ('eta0', 'etabar0', 'K0')] + [(g, 3, 'field') for g in ('R', 'xi', 'xibar', 'M0')]
if QUICK:
    TABLE_JOBS = [(g, 0, 'field') for g in ('P', 'D0', 'M0', 'R', 'xi')]
    TABLE_JOBS_EXT = []


def st_tables():
    import op_runtables as RT
    RT.run(TABLE_JOBS, log=log)
    return {'jobs': [list(j) for j in TABLE_JOBS]}, []


def st_tables_ext():
    import op_runtables as RT
    RT.run(TABLE_JOBS_EXT, log=log)
    return {'jobs': [list(j) for j in TABLE_JOBS_EXT]}, []


def st_tables_ext2():
    import op_runtables as RT
    RT.run(TABLE_JOBS_EXT2, log=log)
    return {'jobs': [list(j) for j in TABLE_JOBS_EXT2]}, []


def st_tables_ext3():
    import op_runtables as RT
    RT.run([] if QUICK else TABLE_JOBS_EXT3, log=log)
    return {'jobs': [list(j) for j in TABLE_JOBS_EXT3]}, []


def _twist_one(N):
    import op_fam as F, ea_spectra as ES
    r = F.twist(N, log=log)
    lc = ES.twist_two(N)
    Lm = sp.Matrix([[sp.sympify(x) for x in row] for row in lc['matrix_units']])
    M0 = sp.Matrix([[sp.sympify(x) for x in row] for row in r['sum_Upsilon0_Delta_matrix']])
    pred = 4 * Lm + (16 - sp.Rational(2, 3) * (4 + N)) * sp.eye(M0.shape[0])
    r['G-D1sum'] = bool((M0 - pred) == sp.zeros(*M0.shape))
    return r


def mk_twist_ext(name, Ns):
    def fn():
        out = _partial_family(name, _twist_one, [] if QUICK else Ns)
        ok = all(r['pass_spectrum'] and r['G-D1sum'] and r['Upsilon_Pn = -1/3 d_n (two-field part)'] for r in out.values())
        return out, [('G-spec F3 extension N = %s: eigenvalues -32 h(S), G-D1sum, P_n' % Ns, {N: out[N]['D1_eigenvalues'] for N in out}, '-32 h(S)', ok)]
    return fn


def mk_fsl2_ext(name, Ns):
    def fn():
        import op_fam as F
        out = _partial_family(name, lambda N: F.fsl2(N, log=log), [] if QUICK else Ns)
        ok = all(r['pass_spectrum'] for r in out.values())
        return out, [('G-spec F4 extension N = %s: eigenvalues -32 h(n+1)' % Ns, {N: out[N]['D1_eigenvalues'] for N in out}, '-32 h(n+1)', ok)]
    return fn


def st_tables_summary():
    import op_tables as T
    files = sorted(glob.glob(os.path.join(T.TAB, '*.json')))
    by = {}
    for f in files:
        try:
            c = json.load(open(f))
        except Exception:
            continue
        by.setdefault((c['generator'], c['n']), []).append(c)
    summ = {}
    gates = []
    tot = {'rows': 0, 'conj_D<2': [0, 0], 'conj_D>=2': [0, 0], 'polevar': [0, 0], 'F_single': [0, 0, 0]}
    compact = []
    for (g, n), chunks in sorted(by.items()):
        s = T.summarize(g, chunks)
        rows = [r for c in chunks for r in c['rows']]
        s['pairs'] = len(chunks)
        pv = [v for r in rows for v in r.get('pole_variation_check', {}).values()]
        s['pole_variation_checks'] = '%d/%d' % (sum(pv), len(pv))
        fs = [r for r in rows if 'F_x_single_output' in r]
        s['F_single_output_rows'] = len(fs)
        s['F_single_output_nonzero'] = sum(1 for r in fs if r['F_x_single_output'] != '0')
        s['conj_with_F_holds'] = sum(1 for r in fs if r.get('conj_holds_with_F'))
        s['seconds'] = round(sum(c['seconds'] for c in chunks), 1)
        summ['%s_n%d' % (g, n)] = s
        tot['rows'] += len(rows)
        if g in T.EVEN:
            tot['conj_D<2'][0] += s['conj_holds_D<2']; tot['conj_D<2'][1] += s['conj_rows_D<2']
            tot['conj_D>=2'][0] += s['conj_holds_D>=2_without_F']; tot['conj_D>=2'][1] += s['conj_rows_D>=2']
            tot['F_single'][0] += s['F_single_output_rows']; tot['F_single'][1] += s['F_single_output_nonzero']
            tot['F_single'][2] += s['conj_with_F_holds']
        tot['polevar'][0] += sum(pv); tot['polevar'][1] += len(pv)
        for r in rows:
            compact.append({k: r[k] for k in ('generator', 'pair', 'n', 'leg_set', 'ordering', 'D', 'eps0', 'poles',
                                              'placements_eps0', 'line_part_eps0') if k in r} |
                           {k: r[k] for k in ('lc', 'ct', 'conj_holds', 'conj_diff', 'F_x_single_output', 'conj_holds_with_F',
                                              'pole_variation_check') if k in r})
    jdump(compact, os.path.join(OUT, 'tables_compact.json'))
    gates.append(('G-D1pair (Conj op-pair) pair by pair on F8, rows with D < 2: Upsilon^0_Delta = 4 lc + C^0_Delta, '
                  'Upsilon_P = C_P, Upsilon^0_M = C^0_M, Upsilon_R = 0', '%d/%d' % tuple(tot['conj_D<2']), 'all',
                  tot['conj_D<2'][0] == tot['conj_D<2'][1] and tot['conj_D<2'][1] > 0))
    gates.append(('G-D1pair on F8 rows with D >= 2 with Upsilon = [eps^0] D_x (finite special vertices excluded, see notes)',
                  '%d/%d' % tuple(tot['conj_D>=2']), 'all', tot['conj_D>=2'][0] == tot['conj_D>=2'][1]))
    gates.append(('G-pole (Ward form, pair by pair): [eps^-m] of the placements other than the line = -4 m ct_m along Delta, 0 along P_b, M, R',
                  '%d/%d' % tuple(tot['polevar']), 'all', tot['polevar'][0] == tot['polevar'][1]))
    gates.append(('finite special vertex {S(x), ct_m}_{H^m} of Def op-defect (iii) on single-output D >= 2 rows: nonzero rows / rows; '
                  'pair form with F_x added holds on', '%d rows, F_x nonzero on %d, pair form with F_x added holds on %d (exactly the rows with F_x = 0)' % tuple(tot['F_single']),
                  'diagnostic (F_x cancels against the special edge on a tree; not part of Upsilon)', True))
    gates.append(('G-noLog (structural): every [eps^0] and pole from the simplex formula (no log eps by construction, Thm c11-oddsimplex)',
                  'all %d rows' % tot['rows'], 'no log eps', True))
    return {'by_generator_n': summ, 'totals': tot, 'compact_rows_file': 'out.nosync/tables_compact.json',
            'chunk_dir': 'out.nosync/tables/', 'generators': {
                'P': 'P_b, b = %s' % list(map(int, __import__('op_tables').B_DIR)),
                'D0': 'Delta, intrinsic part at x_O = 0 (position part = 2 x_O^mu Upsilon_{P_mu}, G-pos)',
                'M0': 'M_Om, v = Om x, Om = %s, intrinsic part' % __import__('op_tables').OMEGA.astype(int).tolist(),
                'R': 'R in so(W_R^perp), R = %s' % __import__('op_tables').RMAT.astype(int).tolist(),
                'xi': 'xi in S_+ (x) W (engine lam_+ basis, coefficients %s)' % __import__('op_tables').XI_COEF,
                'xibar': 'xibar in S_- (x) W^* (engine lam_- basis, same coefficients)',
                'eta0': 'eta in S_+ (x) W^* (engine lams_+), intrinsic part (position part = x_O^mu Upsilon_{rho(e_mu) eta})',
                'etabar0': 'etabar in S_- (x) W (engine lams_-), intrinsic part',
                'K0': 'K_b intrinsic part at x_O = 0'}}, gates


def st_F0():
    import op_fam as F
    out, gates = F.f0(log=log)
    r7, g7 = F.f0_R7(log=log)
    out['R7_descendants'] = r7
    # cocycle scheme psi on F0 (Def c10-hyp (ii)): psi(O) = kappa (w1,w2) int omega (A_V, A_V), s psi(O) = -2 kappa (w1,w2) (dc, A_V)
    l0 = sp.sympify(out['L_1^(0)(O) / ((w1,w2) int omega (dc, A_V)), psi = 0'][0].replace('log2', 'L2').replace('log3', 'L3').replace(' + -', ' - '),
                    locals={'L2': sp.Symbol('log2'), 'L3': sp.Symbol('log3')})
    out['finite part psi with L_1^(0) = 0 on O'] = 'psi(O) = (%s) u (w1,w2) int omega (A_V, A_V)' % sp.expand(l0 / 2)
    return out, gates + g7


def st_F1():
    import op_fam as F
    st = json.load(open(os.path.join(OUT, 'stage_F0.json')))
    dens = st['result']['so6_density_ordered']
    return F.so6_spectra(dens, Lmax=6 if QUICK else 8, log=log)


def st_F2():
    import op_fam as F
    return F.su11(Lmax=5 if QUICK else 6, log=log)


def st_twist():
    import op_fam as F, ea_spectra as ES
    Ns = [0, 1, 2] if QUICK else [0, 1, 2, 3, 4]
    out = {}
    gates = []
    okS, okP, okD = True, True, True
    for N in Ns:
        r = F.twist(N, log=log)
        # G-D1sum on twist two: sum Upsilon^0_Delta = 4 lc + sum C^0_Delta = 4 lc^nat + 16 - 2/3 (4+N) (two-letter lc^nat, composite
        # engine sum C^0 = -2/3 (4+N)), on the two-field part
        lc = ES.twist_two(N)
        Lm = sp.Matrix([[sp.sympify(x) for x in row] for row in lc['matrix_units']])
        M0 = sp.Matrix([[sp.sympify(x) for x in row] for row in r['sum_Upsilon0_Delta_matrix']])
        pred = 4 * Lm + (16 - sp.Rational(2, 3) * (4 + N)) * sp.eye(M0.shape[0])
        r['G-D1sum: sum Upsilon^0 - (4 lc^nat + 16 - 2/3(4+N))'] = [[str(x) for x in row] for row in (M0 - pred).tolist()]
        okD &= (M0 - pred) == sp.zeros(*M0.shape)
        okS &= r['pass_spectrum']; okP &= r['Upsilon_Pn = -1/3 d_n (two-field part)']
        out[N] = r
    gates.append(('G-spec F3 (twist two, gauge completion): D_1 = sum Upsilon^0_Delta - one-body: eigenvalues -32 h(S), S even <= N, N <= %d' % max(Ns),
                  {N: out[N]['D1_eigenvalues'] for N in Ns}, '-32 h(S)', okS))
    gates.append(('G-D1sum on F3: sum Upsilon^0_Delta = sum (4 lc + C^0_Delta) on the two-field part', okD, True, okD))
    gates.append(('Thm op-a (b) on F3: sum Upsilon_{P_n} = -1/3 d_n (two-field part); P_1 = -(1/3 + theta_P) d_n', okP, True, okP))
    return out, gates


def st_fsl2():
    import op_fam as F
    Ns = [1] if QUICK else [1, 2, 3]
    out = {}
    ok = True
    for N in Ns:
        r = F.fsl2(N, log=log)
        out[N] = r
        ok &= r['pass_spectrum']
    return out, [('G-spec F4 (fermion sl(2), gauge completion, on shell): D_1 eigenvalues -32 h(n+1), n odd <= N, N <= %d' % max(Ns),
                  {N: out[N]['D1_eigenvalues'] for N in Ns}, '-32 h(n+1)', ok)]


def st_pos():
    import op_pos as P
    samples = [('AW|AW', 0, 'AV', 3), ('AV|AV', 0, 'AV AV', 5), ('c|AsV', 0, 'AV', 7), ('AW|lamP', 0, 'AW lamP', 9),
               ('lamP|lamM', 0, 'AV', 11), ('AV|AW', 1, 'AW', 13), ('AW|AW', 1, 'AV AV', 15), ('AV|c', 1, 'AsV c', 17)]
    if QUICK:
        samples = samples[:2]
    rows, ok = P.run(samples, log=log)
    return {'xO': [1, -1, 2, 0], 'rows': rows}, [
        ('G-pos Delta: Upsilon_Delta(x_O) = 2 x_O^mu Upsilon_{P_mu} + Upsilon^0_Delta (graph sums, sample rows)', ok['Delta'], True, ok['Delta']),
        ('G-pos M: Upsilon_M(x_O) = (Om x_O)^mu Upsilon_{P_mu} + Upsilon^0_M', ok['M'], True, ok['M']),
        ('G-pos eta: Upsilon_eta(x_O) = x_O^mu Upsilon_{rho(e_mu)eta} + Upsilon^0_eta', ok['eta'], True, ok['eta']),
        ('G-pos K_b: Upsilon_K(x_O) = v(x_O)^mu Upsilon_{P_mu} - (b,x_O) Upsilon^0_Delta + Upsilon^0_{M_{omega_xO}} + Upsilon^0_K (Killing data at x_O from v(x) directly)',
         ok['K'], True, ok['K'])]


def st_nolog():
    import op_nolog as N
    rows, ok = N.run(npts=1 if QUICK else 2, log=log)
    return {'rows': rows}, [('G-noLog (numerical, E = 2 configurations, cubic T5 at the special edge): log eps coefficient 0, pole and eps^0 = simplex values',
                             [r['numeric_fit'] for r in rows], 'C = 0', ok)]


def st_components():
    """OP6 on F0: every component on the scalar bilinears from the master formula, with the one-body terms symbolic"""
    import op_onebody as OB
    f0 = json.load(open(os.path.join(OUT, 'stage_F0.json')))['result']
    comp = {
        'D_1': {'traceless': f0['D_1 traceless (symbolic theta)'] + '  [theta_perp = -1/3 (action side): 0]',
                'trace': f0['D_1 trace (symbolic theta)'] + '  [theta_perp = -1/3: -48]',
                'position part at x_O': '2 x_O^mu P_1(O)|_{P_mu} (G-pos, graph by graph)'},
        'P_1': f0['P_1(O) / d_b O (symbolic thetaP)'] + ' x d_b O',
        'M_1': 'two-body 0 (Upsilon^0_M = 0 on the two-output leg set); one-body -thetaM_perp M_Om.O; position part (Om x_O)^mu P_1|_{P_mu}',
        'R_1': 'two-body 0; one-body -thetaR_perp R.O (Ob^(1)(R) = 0: the direction is free; R8: removable by psi)',
        'Q_1 (xi in S_+ (x) W)': {'two-body (psi = 0)': f0['Q_1 two-body fit (xi in S_+ (x) W)'],
                                  'traceless / F1': f0['Q_1(O)|traceless / F1 (symbolic theta^xi_perp)P']},
        'Q_1 (xibar in S_- (x) W^*)': {'two-body (psi = 0)': f0['Q_1 two-body fit (xi in S_- (x) W^*)'],
                                       'traceless / F1': f0['Q_1(O)|traceless / F1 (symbolic theta^xi_perp)M']},
        'S_1': 'two-body: Q_1 two-body along the Killing spinor rho(x_O) eta (position part, G-pos graph by graph) + Upsilon^0_eta, '
               'whose only two-letter leg set on (A_perp, A_perp) is a single lam (zero on the colour singlet); one-body '
               '-{hbar F_1^(1)(eta), O}: not provided by the action side (open)',
        'K_1': 'two-body: v(x_O)^mu P_1|_{P_mu} - (b,x_O) D^0_1 + M^0_1|_{omega_xO} + K^0_1 (G-pos); Upsilon^0_K(A_perp,A_perp) on (none): %s, '
               'on A_V: zero on the singlet; one-body -{hbar F_1^(1)(K_b), O}: not provided by the action side (open)' % f0['Upsilon^0_K(A_perp,A_perp) on (none) (vacuum)'],
        'L_1^(0)': {'psi = 0': f0['L_1^(0)(O) / ((w1,w2) int omega (dc, A_V)), psi = 0'][0] + ' x (w1,w2) int omega (dc, A_V)',
                    'cocycle scheme': f0['finite part psi with L_1^(0) = 0 on O']},
        'one-body inputs': OB.__doc__.split('Action-side input')[1][:2500],
    }
    return comp, []


def st_twocomp():
    import qa_twocomp as TC
    f0 = json.load(open(os.path.join(OUT, 'stage_F0.json')))['result']
    fit = f0['Q_1 two-body fit (xi in S_+ (x) W)']
    out = TC.forms(fit['a1'], fit['a2'], fit['a3'], fit['a4'], '79/12 - 3 log(3/2)')
    out['c_Q'] = 'Q_1(O)|_traceless = c_Q {S^(1)(xi), O}, c_Q = 79/12 - 3 log(3/2) - theta^xi_perp (particular solution (log(3/2) - 79/12)(T1+T2) of the action side)'
    out['D_1 scalar density'] = ('D_1 tr(... A_ij A_kl ...) adjacent density = -8u (1 - P + K/2) on the SO(6) indices of 1/4 w^{ij} A_ij '
                                 '(Upsilon^0_Delta ordered: same -8/3, exchanged 16, trace -8 per ordered word, planar factor 1/2; '
                                 'one-body -(22/3 + 2 theta) per A letter, theta = -1/3)')
    out['D_1 bilinears'] = 'D_1 int omega w1^{ij} w2^{kl} (A_ij, A_kl)/16 = 0 (traceless), -48 u int omega (A_ij, A^{ij}) / 16 x (trace part)'
    out['P_1'] = 'P_1(O) = -(1/3 + thetaP_perp) u b^m d_m O'
    out['su(1|1)'] = 'two-site D_1 density -8u (1 - Pi) on {Z = 1/4 w^{ij} A_ij, psi = lam_alpha^i component}, graded permutation Pi'
    out['R7'] = 'D_1( int omega [(1/4 w1^{kl} A_kl, w2_{ji} xi^{alpha i} lam_alpha^j) + (1<->2)] ) = -16 u (1/4 w1^{ij} w2_ij) int omega (A_{ji}, xi^{alpha i} lam_alpha^j)'
    return out, []


def st_GS():
    import op_tables as T
    cen = T.census()
    pairs = ['AV|AW', 'AW|AW', 'AW|lamP', 'AW|lamM', 'lamP|lamM', 'lamP|lamP', 'AV|lamP', 'c|AsV', 'AW|AsV', 'lamP|AsV']
    out = {}
    ok = True
    for g in ('P', 'D0', 'M0', 'xi', 'xibar', 'eta0'):
        rows = T.graded_symmetry(g, cen, pairs if not QUICK else pairs[:3])
        out[g] = rows
        for r in rows:
            a_, b_ = r['pair'].split('|'); k_ = r['odd_outputs']; n_ = len(r['leg_set'].split())
            r['predicted'] = str((-1) ** (C.PAR[a_] * C.PAR[b_] + k_ * (k_ - 1) // 2 + n_))
        ok &= all(r['ratio'] == r['predicted'] for r in rows)
        log('  GS %s: %d rows, ratios %s, rule %s' % (g, len(rows), sorted({r['ratio'] for r in rows}, key=str),
                                                     all(r['ratio'] == r['predicted'] for r in rows)))
    return out, [('graded symmetry (Rem op-pair (1)): Upsilon_x(X_b,X_a; reversed outputs) = (-1)^{|X_a||X_b| + k(k-1)/2 + r} Upsilon_x(X_a,X_b) '
                  '(k odd outputs, r outputs: Koszul sign of the letters, of the odd outputs and of the colour word) on every tested row, every generator',
                  {g: '%d rows' % len(rows) for g, rows in out.items()}, 'rule on all rows', ok)]


def st_tfix():
    """Prop prop:op-tfix: free one-body directions acting as -t_x rho_0(x) on cohomology, t constant on the blocks
    P, K, Delta, so(V), so(W), S+W, S-W*, S+W*, S-W (Sigma = S+W + S-W*, Sigma^* = S+W* + S-W); the arity-two relations
    give t_[x,y] = t_x + t_y for every bracket of psl(4|4) between blocks with rho_0([x,y]) != 0.  Exact linear algebra."""
    names = ['P', 'K', 'D', 'soV', 'soW', 'Qp', 'Qm', 'Sp', 'Sm']     # Qp = S+W (xi), Qm = S-W* (xibar), Sp = S+W* (eta), Sm = S-W (etabar)
    t = {n: sp.Symbol('t_' + n) for n in names}
    brackets = [  # (x, y, [z with nonzero coefficient in [x,y]])
        ('D', 'P', ['P']), ('D', 'K', ['K']), ('D', 'Qp', ['Qp']), ('D', 'Qm', ['Qm']), ('D', 'Sp', ['Sp']), ('D', 'Sm', ['Sm']),
        ('soV', 'P', ['P']), ('soW', 'Qp', ['Qp']), ('soV', 'Qp', ['Qp']), ('soW', 'Sp', ['Sp']),
        ('Qp', 'Qm', ['P']), ('Sp', 'Sm', ['K']),
        ('Qp', 'Sp', ['D', 'soV', 'soW']), ('Qm', 'Sm', ['D', 'soV', 'soW']),
        ('P', 'Sp', ['Qm']), ('P', 'Sm', ['Qp']), ('K', 'Qp', ['Sm']), ('K', 'Qm', ['Sp']), ('K', 'P', ['D', 'soV'])]
    eqs = []
    for x, y, zs in brackets:
        for z in zs:
            eqs.append(sp.Eq(t[z], t[x] + t[y]))
    sol = sp.solve(eqs, list(t.values()), dict=True)[0]
    free = sorted({s for v in sol.values() for s in v.free_symbols} | {t[n] for n in names if t[n] not in sol}, key=str)
    gen = {n: sp.simplify(sol.get(t[n], t[n])) for n in names}
    # compare with a B + (c/2) lambda: B = +1 on S+W (xi^i upper), S-W (etabar^i upper), -1 on S-W*, S+W*; lambda = -1, +1, -2, +2 on Sigma, Sigma^*, P, K
    a, c = sp.symbols('a c')
    Bv = {'P': 0, 'K': 0, 'D': 0, 'soV': 0, 'soW': 0, 'Qp': 1, 'Qm': -1, 'Sp': -1, 'Sm': 1}
    lam = {'P': -2, 'K': 2, 'D': 0, 'soV': 0, 'soW': 0, 'Qp': -1, 'Qm': -1, 'Sp': 1, 'Sm': 1}
    target = {n: a * Bv[n] + c / 2 * lam[n] for n in names}
    # the solution space is 2-dimensional and equals span{B, lambda}
    M = sp.Matrix([[gen[n].coeff(f) for f in free] for n in names])
    T2 = sp.Matrix([[Bv[n], lam[n]] for n in names])
    same_span = len(free) == 2 and M.rank() == 2 and sp.Matrix.hstack(M, T2).rank() == 2
    zero = all(gen[n] == 0 for n in ('D', 'soV', 'soW'))
    out = {'general_solution': {n: str(gen[n]) for n in names}, 'free_parameters': [str(f) for f in free],
           'B': Bv, 'lambda': lam, 'span{B, lambda}': bool(same_span), 't_Delta = t_soV = t_soW = 0': bool(zero),
           'two-component': 'residual freedom multiplies the one-loop actions of Q^i_alpha, Qbar_{alphadot i}, S^alpha_i, Sbar^{alphadot i}, P_m, K_m by 1-(a-c/2), 1-(-a-c/2), 1-(-a+c/2), 1-(a+c/2), 1+c, 1-c (Rem op-tc-main)'}
    return out, [('Prop op-tfix: t_[x,y] = t_x + t_y on the psl(4|4) block brackets => t = a B + (c/2) lambda, t_Delta = t_so(V) = t_so(W) = 0',
                  out['general_solution'], 'a B + (c/2) lambda', same_span and zero)]


def st_onebody():
    """OP5: the one-body terms -{hbar F_1^(1)(x), X} on single letters (units u), free directions symbolic"""
    import op_onebody as OB
    th = OB.thD
    tab = {
        'Delta (at x_O; v. = 2 x.d + 2w)': {
            'A|_perp': '-(22/3 + theta_perp (2w + 2x.d)) X',
            'A|_V': '-(22/3 + theta_V (2w + 2x.d)) X  (+ kappa_V: s-exact, zero on cohomology)',
            'lambda (both chiralities)': '-(47/3 + theta_lambda (2w + 2x.d)) X  (+ kappa_lambda: s-exact)',
            'c': '-(8/3 [flagged; alternative 14/3] + theta_V 2x.d) c',
            'A*': '-(-22/3 + theta (2w + 2x.d)) X*, w = 3', 'lambda*': '-(-47/3 + theta (2w + 2x.d)) X*, w = 5/2',
            'c*': '-(-8/3 [flagged] + theta (2w + 2x.d)) c*, w = 4',
            'action values': 'theta_perp = theta_lambda = theta_V = -1/3 (Ob^(2)(Delta,P_b), Thm at-theta; both action engines)',
            'per letter at the action values, x_O = 0': {'A (w = 1 + k)': str(sp.expand(OB.at_action(OB.delta_letter('AW')) + sp.Rational(2, 3) * sp.Symbol('k'))),
                                                         'lambda (w = 3/2 + k)': str(sp.expand(OB.at_action(OB.delta_letter('lamP')) + sp.Rational(2, 3) * sp.Symbol('k')))}},
        'P_b': {'X in sector d': '-thetaP_d d_b X (Ob^(1)(P_b)|_2 = 0; theta_P free, fixed up to a B + (c/2) lambda by the arity-two relations)'},
        'M_Om': {'A|_V': '-(1/3 spin_Om + thetaM_V M_Om.) A_V  (spin term from both action engines)',
                 'A|_perp': '-thetaM_perp M_Om. A_perp',
                 'lambda': '-(m_lambda spin_Om + thetaM_lambda M_Om.) lambda; m_lambda = 1/3 (engine B: S(M)_spin includes lambda), '
                           'engine A: -1/12 int(rho(Om) lam, dslash lam) without preimage normalisation (FLAGGED, both kept)',
                 'c': '-thetaM_V (Om x).d c'},
        'R': {'X in sector d': '-thetaR_d R.X (Ob^(1)(R) = 0, both engines: free)'},
        'xi (both chiralities)': {'A|_V': '-(p + thetaXi_V) Gamma(xi, lam)|_V-arrow: -(p + thetaXi_V) xi.A_V, p = log(3/2) - 79/12',
                                  'A|_perp': '-(p + thetaXi_perp) xi.A_perp',
                                  'lambda': '-(thetaXi_V xi.|_{rho(dA_V)} + thetaXi_perp xi.|_{rho(dA_perp)}) lambda (T1, T2 carry no lambda*)',
                                  'c': '0', 'outer direction G_B': '-a B: -a on S_+ (x) W, +a on S_- (x) W^* (times the classical action)',
                                  'beta-shifts': '{N_f, S_2(xi)} directions (action side, symbolic)'},
        'eta': {'all letters': 'not provided by the action side (Ob^(1)(eta)|_2 intrinsic only): unknown one-body map F_eta'},
        'K_b': {'all letters': 'not provided by the action side (T7 open): unknown one-body map F_K'},
    }
    return tab, []


def st_gluon():
    import op_fam as F
    Ns = [0] if QUICK else [0, 1]
    out = {}
    ok, okw, okp = True, True, True
    for N in Ns:
        r = F.gluon(N, log=log)
        out[N] = r
        ok &= r['pass_spectrum']; okw &= r['ward_longitudinal_zero']; okp &= r['Pn']['letter_zero']
    return out, [('G-spec F6 (gluon sl(2), gauge completion, on shell): D_1 eigenvalues -32 h(n+2), n even <= N, N <= %d' % max(Ns),
                  {N: out[N]['D1_eigenvalues'] for N in Ns}, '-32 h(n+2)', ok),
                 ('F6 Ward: one longitudinal leg -> two-field part of sum Upsilon (D0, P_n, P_beta) = 0 on shell', okw, True, okw),
                 ('Thm op-a (b) on F6: sum Upsilon_{P_n} = -1/3 d_n (two-field part, on shell)', okp, True, okp)]


def _partial_family(name, fn, Ns):
    """run fn(N) for N in Ns with a per-N partial cache out.nosync/partial_<name>.json (resumable)"""
    path = os.path.join(OUT, 'partial_%s.json' % name)
    done = {}
    if os.path.exists(path):
        try:
            done = {int(k): v for k, v in json.load(open(path))['partial'].items()}
        except Exception:
            done = {}
    for N in Ns:
        if N in done:
            log('  %s N=%d: cached' % (name, N))
            continue
        done[N] = fn(N)
        jdump({'partial': {str(k): v for k, v in sorted(done.items())}}, path)
        progress('%s N=%d done: D1 eigenvalues %s, pass %s' % (name, N, done[N]['D1_eigenvalues'], done[N]['pass_spectrum']))
    return {N: done[N] for N in Ns}


def _desc_gate(out, Ns, nblocks):
    """(R7) along P_n, one-sided, as a matrix identity on the two-field part: D_1^(N+1) o d_n = d_n o D_1^(N), with
    d_n tr(d^j X d^{N-j} Y) = tr(d^{j+1} X d^{N-j} Y) + tr(d^j X d^{N+1-j} Y) on each block of the basis"""
    res = {}
    ok = True
    for N in Ns:
        if N + 1 not in out:
            continue
        Ma = sp.Matrix([[sp.sympify(x) for x in row] for row in out[N]['D1_matrix']])
        Mb = sp.Matrix([[sp.sympify(x) for x in row] for row in out[N + 1]['D1_matrix']])
        na, nbb = N + 1, N + 2
        Dm = sp.zeros(nblocks * nbb, nblocks * na)
        for blk in range(nblocks):
            for j in range(na):
                Dm[blk * nbb + j + 1, blk * na + j] += 1
                Dm[blk * nbb + j, blk * na + j] += 1
        diff = sp.simplify(Mb * Dm - Dm * Ma)
        z = diff == sp.zeros(*diff.shape)
        res['%d->%d' % (N, N + 1)] = bool(z)
        ok &= z
    return res, ok


def st_gluon_ext():
    import op_fam as F
    out = _partial_family('gluon_ext', lambda N: F.gluon(N, log=log), [] if QUICK else [2])
    ok = all(r['pass_spectrum'] for r in out.values())
    okw = all(r['ward_longitudinal_zero'] for r in out.values()); okp = all(r['Pn']['letter_zero'] for r in out.values())
    return out, [('G-spec F6 extension N = 2: D_1 eigenvalues -32 h(n+2), n even <= N', {N: out[N]['D1_eigenvalues'] for N in out}, '-32 h(n+2)', ok),
                 ('F6 Ward N = 2: one longitudinal leg -> two-field part of sum Upsilon = 0 on shell', okw, True, okw),
                 ('Thm op-a (b) on F6, N = 2: sum Upsilon_{P_n} = -1/3 d_n', okp, True, okp)]


def st_f7(Ns=None, name='f7'):
    import op_fam2 as F2
    Ns = Ns or ([0, 1] if QUICK else [0, 1, 2, 3])
    out = _partial_family('f7', lambda N: F2.f7(N, log=log), Ns)
    ok = all(r['pass_spectrum'] for r in out.values()); okf = all(r['fit_exact'] for r in out.values())
    okp = all(r['Pn_letter_zero'] for r in out.values())
    dres, dok = _desc_gate(out, Ns, 1)
    out['R7_P_n_descendants'] = dres
    return out, [('G-spec F7 (F_11, psi; gauge completion, on shell): D_1 eigenvalues = -8u x Prop c11-targets (v), S\' <= %d' % max(Ns),
                  {N: out[N]['D1_eigenvalues'] for N in Ns}, {N: out[N]['target_-8u_x_Prop_c11-targets_(v)'] for N in Ns}, ok and okf),
                 ('Thm op-a (b) on F7: sum Upsilon_{P_n} = -1/3 d_n (two-field part, on shell)', okp, True, okp),
                 ('(R7) along P_n on F7 (one-sided, matrix identity): D_1(d_n O) = d_n D_1(O) on the two-field part', dres, 'all True', dok)]


def st_f5(Ns=None, name='f5'):
    import op_fam2 as F2
    Ns = Ns or ([0, 1] if QUICK else [0, 1, 2])
    out = _partial_family('f5', lambda N: F2.f5(N, log=log), Ns)
    ok = all(r['pass_spectrum'] for r in out.values()); okf = all(r['fit_exact'] for r in out.values())
    okp = all(r['Pn_letter_zero'] for r in out.values())
    dres, dok = _desc_gate(out, Ns, 2)
    out['R7_P_n_descendants'] = dres
    return out, [('G-spec F5 (psi3 psi4 with Z\' F_11; gauge completion, on shell): D_1 eigenvalues = -8u x Prop c11-targets (iv), S\' <= %d' % max(Ns),
                  {N: out[N]['D1_eigenvalues'] for N in Ns}, {N: out[N]['target_-8u_x_Prop_c11-targets_(iv)'] for N in Ns}, ok and okf),
                 ('Thm op-a (b) on F5: sum Upsilon_{P_n} = -1/3 d_n (two-field part, on shell)', okp, True, okp),
                 ('(R7) along P_n on F5 (one-sided, matrix identity): D_1(d_n O) = d_n D_1(O) on the two-field part', dres, 'all True', dok)]


def st_spec_desc():
    """(R7) along P_n, one-sided, spectral form on F3, F4, F6: spec D_1 at S' = N is contained in spec D_1 at S' = N + 1"""
    res = {}
    ok = True
    for fam, stg in (('F3', ['twist', 'twist_ext', 'twist_ext2']), ('F4', ['fsl2', 'fsl2_ext', 'fsl2_ext2']), ('F6', ['gluon', 'gluon_ext'])):
        ev = {}
        for sname in stg:
            p = os.path.join(OUT, 'stage_%s.json' % sname)
            if not os.path.exists(p):
                continue
            for k, v in json.load(open(p))['result'].items():
                if isinstance(v, dict) and 'D1_eigenvalues' in v:
                    ev[int(k)] = [sp.sympify(x) for x in v['D1_eigenvalues']]
        r = {}
        for N in sorted(ev):
            if N + 1 in ev:
                big = list(ev[N + 1])
                inc = True
                for x in ev[N]:
                    hit = [i for i, y in enumerate(big) if sp.simplify(x - y) == 0]
                    if not hit:
                        inc = False; break
                    big.pop(hit[0])
                r['%d->%d' % (N, N + 1)] = inc
                ok &= inc
        res[fam] = r
    return res, [('(R7) along P_n, spectral form on F3, F4, F6: the D_1 spectrum at S\' = N lies in the spectrum at S\' = N + 1 '
                  '(descendants: d_n is injective on the two-field part)', res, 'True at every step', ok)]


def st_twocomp_fam():
    """two-component (Weinberg-Srednicki) form of the family operators and of the reported D_1 / P_1 statements (Def c8-twocomp,
    dictionary orientation), units u = -g^2 N/(32 pi^2)"""
    out = {
        'letters': {'D_n': 'n^m (d_m + [A_m, -]) with n = sigma^m_{1 1dot}: the null direction of the engine (n = (1,-i,0,0) in the engine frame)',
                    'Z': 'Z = iota_w A = (1/4) w^{ij} A_ij, (w,w) = 0 (twist two, su(1|1))',
                    'psi': 'psi = lambda_1^i (alpha = 1 component of lambda_alpha^i, fixed i)',
                    'psi3, psi4': 'psi^(3) = lambda_1^3, psi^(4) = lambda_1^4: the alpha = 1 components with the two SU(4) labels other than those of the su(1|1) pair '
                                  '(engine SO(6) weights (-,+,+), (+,-,+))',
                    "Z'": "Z' = A^{34} = (1/2) eps^{34kl} A_kl: the scalar with the SU(4) weight of lambda^3 lambda^4 (engine w' = e_4 - i e_5, read off from "
                          "the A_perp output of Upsilon^0_Delta(psi3, psi4))",
                    'F_11': 'F_11 = (sigma^m sigmabar^n)_1^gamma eps_{gamma 1} F_mn, F_mn = d_m A_n - d_n A_m + [A_m, A_n] (two-field part f = d_n A_m - d_m A_n '
                            'with m = (0,0,1,i) in the engine frame; three-field part [A_n, A_m])'},
        'F3': 'tr(D_n^k Z D_n^{S-k} Z): D_1 = -32 h(S) u on the primaries, S even',
        'F4': 'tr(D_n^k psi D_n^{S-k} psi): D_1 = -32 h(n+1) u, n odd',
        'F5': "tr(D_n^k lambda_1^3 D_n^{S-k} lambda_1^4) with tr(D_n^k A^{34} D_n^{S-k} F_11): D_1 = -8u x {0, 6, 25/3, 49/5, ...} (Prop c11-targets (iv))",
        'F6': 'tr(D_n^k F_11 D_n^{S-k} F_11): D_1 = -32 h(n+2) u, n even',
        'F7': 'tr(D_n^k F_11 D_n^{S-k} lambda_1^i): D_1 = -8u x {6, 25/3, 49/5, ...} (Prop c11-targets (v))',
        'P_1': 'P_1 along n on F3-F7 (two-field part): sum Upsilon_{P_n} = -(1/3) u n^m d_m O, so P_1 = -(1/3 + theta_P) u n^m d_m O',
        'one-body': 'per letter -(phi_X + 2 theta w_X) u at the action values theta = -1/3, phi_A = 22/3 (A_ij and A_m), phi_lambda = 47/3 (lambda_alpha^i, lambdabar_{alphadot i})'}
    return out, []


def main():
    stages = [('OP0', st_OP0), ('tables', st_tables), ('tables_summary', st_tables_summary), ('F0', st_F0), ('F1', st_F1),
              ('F2', st_F2), ('pos', st_pos), ('nolog', st_nolog), ('GS', st_GS), ('tfix', st_tfix), ('twist', st_twist), ('fsl2', st_fsl2), ('gluon', st_gluon),
              ('onebody', st_onebody), ('components', st_components), ('twocomp', st_twocomp), ('twocomp_fam', st_twocomp_fam),
              ('f7', st_f7), ('f5', st_f5), ('gluon_ext', st_gluon_ext), ('fsl2_ext', mk_fsl2_ext('fsl2_ext', [4])),
              ('twist_ext', mk_twist_ext('twist_ext', [5])), ('spec_desc', st_spec_desc),
              ('tables_ext', st_tables_ext), ('tables_summary_n1', st_tables_summary),
              ('f7_ext', lambda: st_f7([0, 1, 2, 3, 4])), ('f5_ext', lambda: st_f5([0, 1, 2, 3])),
              ('fsl2_ext2', mk_fsl2_ext('fsl2_ext2', [5])), ('twist_ext2', mk_twist_ext('twist_ext2', [6])), ('spec_desc_ext', st_spec_desc),
              ('tables_ext2', st_tables_ext2), ('tables_summary_final', st_tables_summary),
              ('tables_ext3', st_tables_ext3), ('tables_summary_ext3', st_tables_summary)]
    results = {'meta': {'engine': 'A (route A: exact simplex/Schwinger generating functions)', 'side': 'operator',
                        'date': time.strftime('%Y-%m-%d %H:%M'), 'unit': 'u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2)',
                        'spec': 'step2/strict/operator/theory/engine_spec.json', 'theory': 'step2/strict/operator/theory/operator_theory.tex',
                        'vendored_sha256': C.provenance(),
                        'status_labels': 'derived-from-construction (every number from the master formula, Upsilon = [eps^0] of the defect cycles); '
                                         'computed by one engine (route A); comparison with route B open'}, 'stages': {}}
    gates = []
    for name, fn in stages:
        d = stage(name, fn)
        if d is None:
            continue
        results['stages'][name] = d['result']
        for g in d['gates']:
            gates.append(dict(g, stage=name))
        results['meta']['wall_s_this_invocation'] = round(time.time() - T_START, 1)
        results['meta']['stage_seconds_total'] = round(results['meta'].get('stage_seconds_total', 0) + d.get('seconds', 0), 1)
        results['meta']['stage_max_rss_MB'] = max(results['meta'].get('stage_max_rss_MB', 0), d.get('max_rss_MB', 0))
        results['meta']['max_rss_MB'] = round(C.maxrss_mb(), 1)
        for dd in (HERE, OUT):
            jdump(results, os.path.join(dd, 'results.json'))
            jdump({'gates': gates, 'n_pass': sum(1 for g in gates if g['pass']), 'n_total': len(gates)}, os.path.join(dd, 'gates.json'))
    ru = resource.getrusage(resource.RUSAGE_SELF)
    msg = 'finished: %d/%d gates pass; wall %.1f s, cpu %.1f s, max RSS %.0f MB' % (
        sum(1 for g in gates if g['pass']), len(gates), time.time() - T_START, ru.ru_utime + ru.ru_stime, C.maxrss_mb())
    log(msg); progress(msg)


if __name__ == '__main__':
    main()
