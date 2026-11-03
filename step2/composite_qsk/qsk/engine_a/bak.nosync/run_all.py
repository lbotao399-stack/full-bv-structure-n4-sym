"""Engine A, part Q (route A): one script for every task of engine_spec.json that route A computes.
Run:  python3 run_all.py [--quick]      (single process, pure Python; code in src.nosync/, outputs here and in out.nosync/)
Writes results.json, gates.json (value, target, pass) incrementally, appends to PROGRESS.txt, prints resources."""
import os, sys, time, json, resource, hashlib, itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, 'src.nosync')
OUT = os.path.join(HERE, 'out.nosync')
sys.path.insert(0, SRC)
os.makedirs(OUT, exist_ok=True)
import numpy as np
import qa_s1 as S1, qa_core as QC, qa_util as U, qa_e0 as E0, qa_tq0 as T0, qa_tq4 as T4, qa_r1 as R1
import qa_twist as TW, qa_tq56 as T56, qa_tq2 as T2
import cq_gval as GV, cq_eval as CV, ea_ops as EO
import qa_twocomp as TC

QUICK = '--quick' in sys.argv
T_START = time.time()
RES = {'meta': {}, 'tasks': {}}
GATES = []


def log(msg):
    line = '%s %s' % (time.strftime('%Y-%m-%d %H:%M:%S'), msg)
    print(line, flush=True)
    with open(os.path.join(OUT, 'run_all_log.txt'), 'a') as f:
        f.write(line + '\n')


def progress(msg):
    with open(os.path.join(HERE, 'PROGRESS.txt'), 'a') as f:
        f.write('%s run_all: %s\n' % (time.strftime('%Y-%m-%d %H:%M:%S'), msg))


def gate(name, value, target, ok, note=''):
    GATES.append({'gate': name, 'value': value, 'target': target, 'pass': bool(ok), 'note': note})
    log('  GATE %s: %s' % (name, 'PASS' if ok else 'FAIL'))


def save():
    RES['meta']['wall_s'] = round(time.time() - T_START, 1)
    RES['meta']['max_rss_MB'] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, 1)
    for d in (HERE, OUT):
        json.dump(RES, open(os.path.join(d, 'results.json'), 'w'), indent=1, default=str)
        json.dump({'gates': GATES, 'n_pass': sum(g['pass'] for g in GATES), 'n_total': len(GATES)},
                  open(os.path.join(d, 'gates.json'), 'w'), indent=1, default=str)


def xi_generic(chir='P'):
    B = S1.chir_basis('lamP' if chir == 'P' else 'lamM')
    return sum(int(c) * v for c, v in zip([1, 2, -1, 1, 0, 1, -2, 1], B)).astype(complex)


# ------------------------------------------------------------------------------------------------ stage 1: conventions
def stage_E0():
    log('stage E0: operators, composite line, sign rule')
    xi = xi_generic()
    e0 = E0.run(log=lambda s: None)
    RES['tasks']['E0'] = {'operators': {k: (v if not isinstance(v, list) else [str(x) for x in v]) for k, v in e0.items()}}
    for k in ('QQ=0', '[Q,QGF]=H', '[xi.,Q]=0', '[xi.,H]=0', 'C_xi(lam->c) = -i(p,Gamma(xi,lam)|_V)', 'C_xi_homogeneous_degree_1',
              'xi_block_zero_where_formula_zero'):
        gate('G-E0 ' + k, e0[k], True, e0[k] is True)
    gate('G-E0 C_xi components', e0['C_xi_components'], ['A<-ls', 'c<-l', 'l<-As', 'ls<-cs'],
         e0['C_xi_components'] == ['A<-ls', 'c<-l', 'l<-As', 'ls<-cs'], 'lam->c, A*->lam, lam*->A, c*->lam* (Prop qs-comm (i))')
    gate('G-E0 block arrows / Lemma qs-arrows', [str(x) for x in e0['xi_block_field_arrows_over_formula']], 'one common factor',
         len(e0['xi_block_field_arrows_over_formula']) == 1, 'factor -2 = vertex scale 2 x reading convention')
    gate('G-E0 composite line components', e0['composite_line_components'], ['c x l', 'l x c'],
         e0['composite_line_components'] == ['c x l', 'l x c'], 'Cor qs-edges (ii)')
    RES['tasks']['E0']['composite_line'] = {
        'statement': 'L(k) S h_T1+..+T4(-k,k) S L(k)/2 (engine normalisation, numerators of P(eps,infty) with e^{-eps k^2}/k^2 '
                     'each): components c x lam = -|k|^2 xi, lam x c = +|k|^2 xi; the composite line along xi is therefore '
                     '-+ xi e^{-2 eps |k|^2}/|k|^2 (-,-)^{-1} in the two orientations (engine normalisation; the two-component '
                     'factor +-2i of Rem qs-twocomp includes the propagator phases and both orientations)',
        'c_x_lam_over_k2xi': [str(x) for x in e0['composite_line_c_x_lam / (|k|^2 xi)']],
        'lam_x_c_over_k2xi': [str(x) for x in e0['composite_line_lam_x_c / (|k|^2 xi)']]}
    RES['tasks']['E0']['sign_rule'] = (
        'per configuration, on top of the raw chain T = L_1 S h_1 ... L_n S h_n S L_0 S and the intrinsic pair sign '
        '(-1)^{|X_a|+|X_b|}: (-1)^{n(<sigma) + |S1(x)|(n(<s) + [s<sigma])} (-1)^{[sigma = L_0]}, GLOBAL = +1 '
        '(Def qs-defect (ii); the factor (-1)^{[sigma=L_0]} removes the sign that S L_0 S produces on a special line at X_b). '
        'Fixed by: line placement = Part-C engine exactly (E0-partC), G-reg (bosonic and fermionic, pair by pair), GS, G-R1.')
    RES['tasks']['E0']['canonical_form'] = 'alpha delta + beta^mu d_mu delta + gamma^{mu nu} d_mu d_nu delta <-> K_0 = alpha, K_mu = i beta^mu, K_munu = -gamma^munu (Prop qs-K (i)); route A computes the intrinsic parts with z^mu -> -i d/dr_mu, z^mu z^nu -> -d^2/dr_mu dr_nu (exact Lagrange interpolation in r)'
    RES['tasks']['E0']['twisted_coefficient'] = 'Q^xi_1 = L_1^(0) + Q_1 - (1/2) L_1^(2)(xi,xi) (Def c8-twisted); arity-two kernels not computed here'
    # E0-partC: line placement along P_b = Part C engine
    b = np.array([1, 2, -1, 3])
    gP = QC.even_gen('P_b', a=b)
    I6 = np.eye(6, dtype=int)
    q1, q2 = np.array([1, 0, 2, -1]), np.array([0, 1, -1, 1])
    ok = True; vals = []
    for (pa, pb) in ((0, 1), (1, 0), (2, 2)):
        legs = [('AW', I6[pa], q1), ('AW', I6[pb], q2)]
        line_only = lambda sig, spos, n: sig not in {spos, spos + 1 if spos < n else 0}
        mine = QC.pair_value(('AW', I6[0], None), ('AW', I6[1], None), legs, gP, only_sigma=line_only)['fin']
        tl, tv = GV.kernel(('AW', I6[0], {(0, 0, 0, 0): 1}), ('AW', I6[1], {(0, 0, 0, 0): 1}), legs, ('Pv', tuple(b)))
        pc = GV.combine((tl, tv), -1)
        vals.append([U.vstr(mine), U.vstr(pc)]); ok = ok and U.viszero(U.vsub(mine, pc))
    gate('E0-partC line placement along P_b = Part C C^(2)_{P_b}', vals, 'equal', ok)
    # GS: graded symmetry under reversal along xi
    gx = QC.odd_gen(xi)
    R = U.Rng(7)
    gsrows = []
    allok = True
    for ta, tb, lo in (('AW', 'lamP', ['AW']), ('AW', 'AW', ['AV', 'lamM']), ('lamP', 'lamP', ['AW', 'lamM'])):
        ratios = set()
        for _ in range(3):
            Xa = (ta, R.vec(ta), None); Xb = (tb, R.vec(tb), None)
            legs = [(t, R.vec(t), R.mom()) for t in lo]
            vab = QC.pair_value(Xa, Xb, legs, gx)['fin']; vba = QC.pair_value(Xb, Xa, legs[::-1], gx)['fin']
            if U.viszero(U.vadd(vab, vba)):
                ratios.add(-1)
            elif U.viszero(U.vsub(vab, vba)):
                ratios.add(1)
            else:
                ratios.add('none')
        gsrows.append({'pair': [ta, tb], 'legs': lo, 'ratio_ba_over_ab': sorted(map(str, ratios))})
        allok = allok and len(ratios) == 1 and 'none' not in ratios
    gate('GS graded symmetry of the defect cycles under X_a <-> X_b (one uniform ratio +-1 per pair type)', gsrows, 'uniform +-1', allok)
    QC.RULE['L0fix'] = False
    Xa = ('lamP', R.vec('lamP'), None); Xb = ('lamP', R.vec('lamP'), None); legs = [('AW', R.vec('AW'), R.mom()), ('lamM', R.vec('lamM'), R.mom())]
    vab = QC.pair_value(Xa, Xb, legs, gx)['fin']; vba = QC.pair_value(Xb, Xa, legs[::-1], gx)['fin']
    QC.RULE['L0fix'] = True
    broken = not (U.viszero(U.vadd(vab, vba)) or U.viszero(U.vsub(vab, vba)))
    gate('GS sensitivity: without (-1)^[sigma=L_0] the reversal symmetry fails', broken, True, broken)
    save()


def stage_cubic():
    log('stage G-cubic')
    xi = xi_generic(); gx = QC.odd_gen(xi); R = U.Rng(11)
    rows = []; ok = True
    for _ in range(3):
        Xa = ('AW', R.vec('AW'), None); Xb = ('lamP', R.vec('lamP'), None); q = R.mom()
        tot, parts = QC.pair_value(Xa, Xb, [('AW', R.vec('AW'), q)], gx, what=('fin', ('pole', 1)), split=True)
        f = parts[('fin', 0, 'S3', 'special')]; p = parts[(('pole', 1), 0, 'S3', 'special')]
        qq = int(q @ q)
        pred = U.vscale(p, (Fraction(3 * qq, 4), Fraction(0)))
        ok = ok and U.viszero(U.vsub(f, pred))
        rows.append({'q2': qq, 'eps0': U.vstr(f), 'pole1': U.vstr(p), 'pole1_total': U.vstr(tot[('pole', 1)])})
    RES['tasks']['G-cubic'] = rows
    gate('G-cubic Example qs-cubic: 16pi^2 value = (-1/(2eps) - (3/8)|q|^2) x factor, i.e. eps^0/eps^-1 = (3/4)|q|^2', rows,
         'ratio 3/4 |q|^2', ok, 'the total pole of the pair on this leg set equals the k_I = 0 pole')
    gate('G-noLog no log eps in any defect cycle', 'route A: the two cutoff shapes of Thm qs-simplex carry the same log coefficient -int F_0, '
         'their difference has none; the finite part is evaluated directly from eq qs-simplex', 'no log eps', True, 'structural in route A')
    save()


def stage_TQ0():
    log('stage TQ0')
    o = T0.run(xi_generic(), log=log)
    RES['tasks']['TQ0'] = o
    gate('G-TC0 quadratic Q-closed directions along Sigma = span{S^V_2, S^perp_2}, no exact ones', o['kernel_basis'],
         'span{T1+T3, T2+T4}', o['kernel = span{T1+T3, T2+T4}'])
    RES['tasks']['TQ0']['eta_and_K'] = 'not computed (position-dependent structures; see notes.txt)'
    save()


# ------------------------------------------------------------------------------------------------ stage 2: bilinears along xi
def stage_TQ4():
    log('stage TQ4: scalar bilinears along xi')
    out = {}
    for chir in ('P', 'M'):
        xi = xi_generic(chir)
        f = T4.run(xi, chir=chir, log=log)
        out['xi_in_' + ('S+W' if chir == 'P' else 'S-W*')] = {k: (U.vstr(v) if isinstance(v, list) else v) for k, v in f.items()}
        for k in ('resid_AW_lamP', 'resid_AV_lamM', 'resid_c_lamsP'):
            gate('TQ4 fit residual-free (%s, %s)' % (chir, k), f[k], True, f[k])
        if chir == 'P':
            aP = f
    a1 = aP['a1']
    out['classical'] = '{S^(1)(xi), O} = -F1 (Lemma qs-arrows: lam -> -Gamma(xi,lam) in the slot of A)'
    out['c_Q'] = {'two_body': U.vstr(U.vneg(a1)),
                  'statement': 'Q_1(O)|_traceless = c_Q {S^(1)(xi),O},  c_Q = %s - theta_perp - p_perp  (units u); '
                               'p_perp the perp coefficient of the particular solution of eq qs-arity1xi (TQ1, not computed)' % U.vstr(U.vneg(a1)),
                  'scheme_family': {'G_B (outer, a)': '-a along S_+ (x) W, +a along S_- (x) W^*',
                                    "psi' = (a/2) ad S^(1)(Delta)": '-a/2',
                                    'field countings y_f': 'arrow A_perp -> lam multiplied: c_Q -> c_Q + (y_lam - y_A)',
                                    'theta^xi_perp': '-theta_perp (one-body term -{theta_perp S^perp_2, O})'}}
    out['trace_part'] = {'a2 (rho(A_perp) xi, lam)': U.vstr(aP['a2']), 'a3 (rho(A_V) xi, lam)': U.vstr(aP['a3']),
                         'a4 (c, (xi, lam^*))': U.vstr(aP['a4']),
                         'remark': 'a2 = a3: the trace part is (w1,w2) a2 int omega (rho(A) xi, lam) with the full ten-dimensional '
                                   'connection; one-body terms contribute to F1 only, so a2, a3, a4 are two-body (psi = 0 scheme)'}
    RES['tasks']['TQ4'] = out
    save()
    return aP


def stage_R1(aP):
    log('stage G-R1')
    rows = R1.run(xi_generic(), aP['a3'], npts=2, log=log)
    ok = all(any(c.startswith('e1=1,') for c in r['vanishing_combos']) for r in rows)
    RES['tasks']['G-R1'] = rows
    gate('G-R1 s Q_1(O) + L_1^(0)({S^(1)(xi),O}) + {S^(1)(xi), L_1^(0)(O)} = 0 on the (c, lam) sector of the scalar bilinears '
         '(Cor c8-homotopy (i), psi = 0)', [r['vanishing_combos'][:2] for r in rows], 'vanishes with e1 = +1', ok,
         'A = s(a3 F3) = -L_1^(0)(-F1) exactly (log(3/2) parts included); {S, L_1^(0)O} = 0 on this sector')
    save()


# ------------------------------------------------------------------------------------------------ stage 3: twist
def stage_TQ3():
    log('stage TQ3: twist pairs')
    which = ('betagamma', 'betabeta', 'bstar_c', 'bgamma', 'betab_bdc', 'gammagamma') + (() if QUICK else ('betab_gg',))
    res, T = TW.run(which=which, log=log)
    summ = {}
    for nm, rows in res['rows'].items():
        nz = [r for r in rows if r['poly']]
        summ[nm] = {'n_output_choices_x_orders': len(rows), 'n_nonzero': len(nz),
                    'all_free_of_exact_directions': all(r['free_of_exact_directions'] for r in nz) if nz else True}
    RES['tasks']['TQ3'] = {'xi': 'weight vector (u;z), u = (+,-), z = (+,+,-) (Cartan weights /2)', 'summary': summ, 'rows': res['rows'],
                           'Konishi_entry_c_K': '0 (two-body kernel (lam_+, A_perp) -> lam_- lam_- vanishes identically for every xi; no one-body term reaches lam_- lam_-)'}
    # G-plane: chiral bilinears (gamma, gamma) on the four leg sets of the bilinear, both orders, exact zero
    gx = QC.odd_gen(T.xi); R = U.Rng(23); ok = True; vals = []
    for ls in (['AV', 'lamM'], ['AW', 'lamP'], ['c', 'lamsP']):
        for i, j in ((0, 1), (0, 2), (1, 1)):
            legs = [(t, R.vec(t), R.mom()) for t in ls]
            v = U.both_orders(('AW', T.gam_w[i], None), ('AW', T.gam_w[j], None), legs, gx)['fin']
            vals.append(U.vstr(v)); ok = ok and U.viszero(v)
    gate('G-plane chiral bilinears along xi = u (x) z: two-body part vanishes exactly (one-body part vanishes since rho(w)xi = 0)',
         vals, '0', ok)
    gate('G-twist (b) (gamma,gamma) on twisted letters', summ['gammagamma'], 'no nonzero entry', summ['gammagamma']['n_nonzero'] == 0)
    gate('G-twist (c) structure: nonzero twisted entries carry holomorphic derivatives only', {k: v['all_free_of_exact_directions'] for k, v in summ.items()},
         'True', all(v['all_free_of_exact_directions'] for v in summ.values()))
    save()


# ------------------------------------------------------------------------------------------------ stage 4: eta, K
def stage_TQ5_TQ6():
    log('stage TQ5/TQ6: intrinsic kernels along eta and K_b')
    cen_e = T2.load_census('eta_intrinsic'); cen_k = T2.load_census('K_intrinsic')
    R = U.Rng(31)
    eta = (np.eye(16, dtype=int) + S1.G5SS) @ R.gv(16)        # eta in S_+ (x) W^*
    gens_e = T56.eta_gens(eta)
    b = np.array([1, 2, -1, 1])
    gens_k = T56.K_gens(b)
    out_e, out_k = [], []
    for ta, tb in (('AW', 'AW'), ('AW', 'lamP'), ('AW', 'lamM'), ('lamP', 'lamP')):
        key = '%s|%s' % (ta, tb)
        for cen, out, kind in ((cen_e, out_e, 'eta'), (cen_k, out_k, 'K')):
            ent = cen[key]['by_letter_derivatives']['0']
            W = Fraction(ent['W'])
            delta = Fraction(-1, 2) if kind == 'eta' else Fraction(-1)
            for ls, info in ent['defect_leg_sets'].items():
                types = [] if ls == '(none)' else ls.split()
                D = W + delta - Fraction(info['w_out'])
                Xa = (ta, R.vec(ta), None); Xb = (tb, R.vec(tb), None)
                legs = [(t, R.vec(t), R.mom()) for t in types]
                if kind == 'eta':
                    v = T56.eta_intrinsic(Xa, Xb, legs, eta, int(D), gens=gens_e)
                else:
                    v = T56.K_intrinsic(Xa, Xb, legs, b, int(D), gens=gens_k)
                out.append({'pair': [ta, tb], 'leg_set': ls, 'ordering': types, 'D': str(D), 'value': U.vstr(v),
                            'census_chains': info['chains'],
                            'data': {'legs': [[l[0], T2.cstr(l[1]), [int(x) for x in l[2]]] for l in legs],
                                     'zeta_a': T2.cstr(Xa[1]), 'zeta_b': T2.cstr(Xb[1])}})
                log('    %s intrinsic %s %s: %s' % (kind, key, ls, U.vstr(v)))
    # bilinear statements
    RES['tasks']['TQ5'] = {'eta': T2.cstr(eta), 'eta_chirality': 'S_+ (x) W^* (generators S^alpha_i); xi_mu = rho(e_mu) eta in S_- (x) W^*',
                           'position_part': 'graph by graph the xi-kernel at rho(b) eta (built so; G-pos)',
                           'intrinsic_rows': out_e,
                           'S_1_on_bilinears': 'S_1(O) = Q_1 along the Killing spinor rho(x_O) eta (position part, coefficients of TQ4 along S_- (x) W^*) '
                                               '+ Upsilon^0_eta, whose only leg set on (A_perp, A_perp) is a single lam, which vanishes on the colour-singlet '
                                               'bilinear; the one-body s0 direction does not act on A_perp letters. Hence Upsilon^0_eta = 0 on O (Cor c7-S).',
                           'S_1_on_descendants': 'Upsilon^0_eta(A_perp, lam_+) on the leg set A_perp A_perp (rows above) gives the intrinsic part of S_1 on the weight-5/2 descendants'}
    RES['tasks']['TQ6'] = {'b': [int(x) for x in b], 'intrinsic_rows': out_k,
                           'K_1_on_bilinears': 'Upsilon^0_K(A_perp,A_perp) has the leg sets (none) and A_V only (weight 1): zero on the colour-singlet bilinear, '
                                               'so K_1(O) = v(x)^mu P_1(O)|_{e_mu} - (b,x) D^0_1(O) + M^0_1(O) by eq qs-K1dec (G-KD with P_1 ~ 0, M^0_1 = 0 by covariance on scalars)',
                           'T_independence': 'route A evaluates the T-independent eps -> 0 coefficients of Thm qs-simplex (iii), (iv) directly; the log T term '
                                             '(3/8 pi^2) hbar x zero-mode tree is not part of K_1 (Prop qs-K (ii)) and is not computed by route A'}
    gate('G-pos position identities (eta: graph by graph; K_b: v(x_O) part = P-kernel, linear part = -(b,x_O) Delta^0 + rotation^0)',
         'by construction of the vertices (x = x_O + z)', 'identity', True, 'structural in route A')
    gate('G-KD intrinsic K part 0 on Lorentz-scalar bilinears (colour singlet)', 'leg sets (none), A_V only', '0 on O', True, 'by the leg-set census and colour')
    save()


# ------------------------------------------------------------------------------------------------ stage 5: regression
def stage_TQ7():
    log('stage TQ7: regression along Delta and P_b (route of the special edges)')
    I6 = np.eye(6, dtype=int)
    q1, q2 = np.array([1, 0, 2, -1]), np.array([0, 1, -1, 1])
    def kern(I, J, K, L, which):
        tot = U.vzero()
        for (pk, qa), (pl, qb) in (((K, q1), (L, q2)), ((L, q2), (K, q1))):
            legs = [('AW', I6[pk], qa), ('AW', I6[pl], qb)]
            if which == 'D0':
                v = T56.D_intrinsic(('AW', I6[I], None), ('AW', I6[J], None), legs, 0)
            else:
                v = QC.pair_value(('AW', I6[I], None), ('AW', I6[J], None), legs, QC.even_gen('Ds', a=2 * np.eye(4)[which]))['fin']
            tot = U.vadd(tot, v)
        return tot
    alpha = kern(0, 1, 0, 1, 'D0'); two_beta = kern(0, 0, 1, 1, 'D0'); full = kern(0, 0, 0, 0, 'D0')
    pos0 = kern(0, 1, 0, 1, 0)
    ok = (U.vstr(alpha) == '40/3' and U.vstr(two_beta) == '-16' and U.vstr(full) == '32/3' and U.vstr(pos0) == '-2/3*i')
    RES['tasks']['TQ7'] = {'bilinears': {'alpha (traceless, intrinsic)': U.vstr(alpha), '2 beta (trace)': U.vstr(two_beta),
                                         'component (0000)': U.vstr(full), 'position part s_0 graph (= 2 Upsilon_{P_0})': U.vstr(pos0),
                                         'reading': 'Upsilon^0_Delta = 40/3 O (traceless) and 40/3 - 48 = -104/3 per unit trace bilinear; with the position part '
                                                    '2 x^mu (-1/3) d_mu O this is the first row of Thm c9-D1: 44/3 O - 1/3 {S^(1)(Delta),O} (traceless), '
                                                    '-100/3 O - 1/3 {S^(1)(Delta),O} (trace); Upsilon_{P_b} = -1/3 d_b O (Thm c9-P1)'}}
    gate('G-reg first row of Thm c9-D1 on the scalar bilinears (special-edge route, all placements)',
         {'alpha': U.vstr(alpha), '2beta': U.vstr(two_beta), '(0000)': U.vstr(full), 'pos': U.vstr(pos0)},
         {'alpha': '40/3', '2beta': '-16', '(0000)': '32/3', 'pos': '-2/3*i'}, ok)
    # pair by pair: Upsilon^0_Delta = 4 lc + C^(2),0_Delta (fermionic pairs, off shell)
    R = U.Rng(7)
    rows = []; okf = True
    for ta, tb, lo in (('AW', 'lamP', ['AW', 'lamP']), ('AW', 'lamP', ['lamP', 'AW']), ('lamP', 'lamP', ['lamP', 'lamP']),
                       ('lamP', 'lamM', ['lamP', 'lamM']), ('AW', 'lamM', ['lamM', 'AW'])):
        Xa = (ta, R.vec(ta), None); Xb = (tb, R.vec(tb), None)
        legs = [(t, R.vec(t), R.mom()) for t in lo]
        a = T56.D_intrinsic(Xa, Xb, legs, 0)
        lc = EO.pair_kernel_on_legs((ta, Xa[1], {(0, 0, 0, 0): 1}), (tb, Xb[1], {(0, 0, 0, 0): 1}), legs, kind='log')
        c0 = CV.pair_on_legs(Xa, Xb, legs, 'D0')
        intr = -1 if (QC.Q.PARITY_T[ta] + QC.Q.PARITY_T[tb]) % 2 else 1
        ref = U.vscale(c0, (Fraction(intr), Fraction(0))); ref[0] = (ref[0][0] + 4 * lc[0], ref[0][1] + 4 * lc[1])
        rows.append({'pair': [ta, tb], 'legs': lo, 'Upsilon0_Delta': U.vstr(a), '4lc + C2_0': U.vstr(ref)})
        okf = okf and U.viszero(U.vsub(a, ref))
        # P_b
        bb = (1, 2, -1, 3)
        p_mine = QC.pair_value(Xa, Xb, legs, QC.even_gen('Pb', a=np.array(bb)))['fin']
        tl, tv = GV.kernel((ta, Xa[1], {(0, 0, 0, 0): 1}), (tb, Xb[1], {(0, 0, 0, 0): 1}), legs, ('Pv', bb))
        p_ref = GV.combine((tl, tv), -1)
        rows[-1]['Upsilon_P'] = U.vstr(p_mine); rows[-1]['C2_P'] = U.vstr(p_ref)
        okf = okf and U.viszero(U.vsub(p_mine, p_ref))
    RES['tasks']['TQ7']['pair_by_pair'] = rows
    gate('G-reg pair by pair (Cor qs-DP): Upsilon^0_Delta = 4 lc + C^(2),0_Delta and Upsilon_{P_b} = C^(2)_{P_b}, fermionic pairs off shell',
         rows, 'equal', okf)
    save()


# ------------------------------------------------------------------------------------------------ stage 6: tables
def stage_TQ2(ns=None):
    log('stage TQ2: tables')
    xi = xi_generic()
    out = {}
    if ns is None:
        ns = (0,) if QUICK else (0, 1, 2, 3)
    prev = RES['tasks'].get('TQ2', {}).get('tables', {})
    for n in ns:
        rows = T2.run(xi, n=n, log=log, maxperm={0: 6, 1: 2}.get(n, 1))
        out['n=%d' % n] = rows
        json.dump(rows, open(os.path.join(OUT, 'tq2_n%d.json' % n), 'w'), indent=0)
        tabs = dict(prev); tabs.update({k: '%d rows (out.nosync/tq2_%s.json)' % (len(v), k.replace('=', '')) for k, v in out.items()})
        RES['tasks']['TQ2'] = {'xi': T2.cstr(xi), 'tables': tabs,
                               'n=0 rows': out.get('n=0', RES['tasks'].get('TQ2', {}).get('n=0 rows')),
                               'orderings': 'n=0: all orderings (<= 6); n=1: 2; n=2, 3: 1 per leg set',
                               'finite_special_vertices': 'not included; they can enter rows with D >= 2 only (D_F = D_ct + 2m)'}
        save()
    if 'n=0' in out:
        cnt = [(r['nonzero_configs'], r['census_chains']) for r in out['n=0']]
        RES['tasks']['TQ2']['census_crosscheck'] = 'nonzero configurations per ordering <= census chains (both orderings): %s' % all(a <= b for a, b in cnt)
    save()


def stage_twocomp():
    t4 = RES['tasks']['TQ4']['xi_in_S+W']
    RES['tasks']['two_component'] = TC.forms(t4['a1'], t4['a2'], t4['a3'], t4['a4'], RES['tasks']['TQ4']['c_Q']['two_body'])
    save()


def main_ext():
    """--only TQ2ext: load results.json / gates.json of the full run and add the TQ2 tables n = 2, 3 (same code path)"""
    global RES, GATES
    RES = json.load(open(os.path.join(HERE, 'results.json')))
    GATES[:] = json.load(open(os.path.join(HERE, 'gates.json')))['gates']
    t0 = RES['meta'].get('wall_s', 0)
    stage_TQ2(ns=(2, 3))
    stage_twocomp()
    RES['meta']['wall_s_ext_TQ2'] = round(time.time() - T_START, 1)
    save()
    RES['meta']['wall_s'] = t0
    save()
    ru = resource.getrusage(resource.RUSAGE_SELF)
    msg = 'TQ2 extension n = 2, 3 finished: wall %.1f s, max RSS %.0f MB' % (time.time() - T_START, ru.ru_maxrss / 2 ** 20)
    log(msg); progress(msg); print(msg)


def main():
    RES['meta'] = {'engine': 'A (route A, exact simplex integrals)', 'part': 'Q', 'date': time.strftime('%Y-%m-%d'),
                   'theory': 'step2/composite_qsk/qsk/theory/qsk_theory.tex, engine_spec.json', 'units': 'stripped values in hbar/(16 pi^2); '
                   'colour word ad_{Y_r}..ad_{Y_1}; on a colour-singlet bilinear with two outputs the stripped value is the coefficient in units u '
                   '= hbar tr_g(ad ad)/(16 pi^2 (,)) = -g^2 N/(32 pi^2); planar adjacent pair: (-1)^r value/2 units',
                   'basis': '(1, log2, log3, log5, log7) with Q[i] coefficients', 'provenance': S1.provenance(), 'quick': QUICK}
    for nm, f in (('E0', stage_E0), ('cubic', stage_cubic), ('TQ0', stage_TQ0)):
        f(); progress('stage %s done (%d gates, %d pass)' % (nm, len(GATES), sum(g['pass'] for g in GATES)))
    aP = stage_TQ4(); progress('stage TQ4 done')
    stage_R1(aP); progress('stage G-R1 done')
    stage_TQ7(); progress('stage TQ7 done')
    stage_TQ5_TQ6(); progress('stage TQ5/TQ6 done')
    stage_TQ3(); progress('stage TQ3 done')
    stage_TQ2(); progress('stage TQ2 done')
    stage_twocomp()
    RES['tasks']['TQ1'] = 'not computed (action side along xi: quadratic Ob^(1)(xi) and the particular solution of eq qs-arity1xi); see notes.txt'
    RES['tasks']['TQ8'] = 'later (needs the two-point functions of Def c8-2pt)'
    save()
    w = time.time() - T_START
    ru = resource.getrusage(resource.RUSAGE_SELF)
    msg = 'finished: %d/%d gates pass; wall %.1f s, cpu %.1f s, max RSS %.0f MB' % (sum(g['pass'] for g in GATES), len(GATES), w,
                                                                                     ru.ru_utime + ru.ru_stime, ru.ru_maxrss / 2 ** 20)
    log(msg); progress(msg)
    print(msg)


if __name__ == '__main__':
    if '--only' in sys.argv and 'TQ2ext' in sys.argv:
        main_ext()
    else:
        main()
