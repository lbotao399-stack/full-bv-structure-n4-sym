"""Assembly of results.json and gates.json (operator engine B) from the stored outputs in out.nosync/."""
import os, sys, json, time, glob, hashlib
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ob_rows as RW

ROOT = os.path.abspath(os.path.join(HERE, '..'))
OUT = RW.OUT


def load(name):
    fn = os.path.join(OUT, name)
    return json.load(open(fn)) if os.path.exists(fn) else None


def sha(fn):
    return hashlib.sha256(open(fn, 'rb').read()).hexdigest()


def main():
    gates = []
    res = {}
    res['meta'] = dict(
        engine='B (route B: momentum-space numerator algebra, exact Feynman-parameter integrals, numerical eps-extrapolation at finite T)',
        side='operator', date=time.strftime('%Y-%m-%d %H:%M'),
        theory='step2/strict/operator/theory/operator_theory.tex, engine_spec.json',
        unit='u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); pair kernels in units hbar/(16 pi^2) x colour words C_w = (1 (x) ad_{Y_w1}...ad_{Y_wr})(-,-)^{-1}, X_a on the first factor; planar adjacent pair: a C_w -> (-1)^r a/2 units u (outputs in trace order Y_wr...Y_w1; labelled form of two outputs of the same field type divided by 2)',
        regulator='Schwinger cutoff t >= eps on every edge P(eps,T); special edge K_eps - K_0 (t in [0,eps]); T = infinity for the exact coefficients (generic momenta); finite T = 1, 4, 16 along K_b in the numerical gate; no dimensional regularisation; no IR datum other than T',
        vendored={os.path.basename(f): sha(f) for f in sorted(glob.glob(os.path.join(HERE, 'vendor', '*.py')))},
        vendored_note='src.nosync/vendor: verified code of step2/composite_qsk/composite/engine_b, step2/composite_qsk/qsk/engine_b, step2/twoletter/engine_b, Step-1 engine_b; edits: import paths only (ce_core EB1, ce_job TL, qb_core CEB/VEND -> the flat vendor folder); sha256 of the originals in vendor_sha256.txt',
        action_side_inputs='step2/strict/action/engine_a/results.json and engine_b/results.json (both engines agree on every number used: Ob1(Delta)|_2 = u[-11/3 int(dA,dA) + 47/3 int(lam,dslash lam) - 14/3 int(dc,A*)], theta_Delta = -1/3 in every sector, Ob1(P_b)|_2 = 0, Ob1(R) = 0, F_1^(1)(M_w)|_2 = (1/3) u S(M_w)_spin + sector parts, F_1^(1)(xi)|_2 = (L - 79/12)(T1 + T2) + theta^xi directions)')
    import ce_job as CJ
    import qb_core as B
    res['meta']['component_names'] = {slot: {str(c): CJ.T.comp_name(slot, c) for c in B.SLOT[slot][1]}
                                      for slot in ('AV', 'AW', 'lamP', 'lamM', 'c', 'AsV', 'AsW', 'lamsP', 'lamsM', 'cs')}
    res['meta']['ring_names'] = ('o<j>_<c>: polarisation component c of output j; r<j>_<m>: momentum of output j (i r = d on '
                                 'that output); za_<c>, zb_<c>: letter components of X_a, X_b; s<m>, u<m>: point-splitting '
                                 'derivative variables of X_a, X_b (coefficient of s_m <-> d_m X_a); b<m>: P_b and K_b direction; '
                                 'w<mn>: M_w with w_nm = -w_mn; R<MN>: R in so(W_R^perp); xi<j>, et<j>: Sigma, Sigma^* components. '
                                 'Two-component names of the engine components: component_names; engine vector index 0 <-> m = 4 and '
                                 'the S_+ / W index maps: step2/notation/dictionary.json (Def c8-twocomp, Srednicki Wick-rotated)')
    res['findings'] = [
        dict(id='B-F1', severity='major',
             statement='Finite special vertex double count. The S-channel finite special vertex F_x = -sum_m ((-1)^m/m!) {S^(1)(x), ct_m(X_a X_b)}_{H^m} of Def op-defect (iii) (and the I-channel {I, ct^(1)_{x,m}}_{H^m}, eq c10-finitesp in arity zero) is cancelled by the (Br) graphs in which an output of the bare pair cycle is joined to S^(1)(x) (resp. to a vertex of I) by the special edge: (Br, bare cycle) + (special part of -{S^(1)(x), d_delta I^CT[eps](O)}_eps) = (e^{-eps H} - 1) x (renormalised cycle) -> 0. Thm op-master (ii) lets (Br) tend to 0 while (iii) keeps F_x: the pieces are counted once on each side.',
             evidence='(lamP,lamM) -> A_V, n = 0 (the only n = 0 field-pair row with ct_1 != 0): [eps^0] D^special_Delta = 4 lc exactly and [eps^0] D_{P_b} = [eps^0] D_{M_w} = [eps^0] D_R = 0, so Upsilon_x := [eps^0] D_x satisfies Conjecture op-pair and the independent second evaluation (C_x from the composite kernel form of CJob, lc from the arity-zero counterterm) while Upsilon + F_x misses it by exactly F_Delta = 6 H ct_1 and F_{P_b} = H d_b ct_1; with F_x Thm op-a (b) would fail on the cocycle tr(lambda lambda Z) (sum F_{P_b} = d_b H ct_1(O) != 0).',
             consequence='All components below use Upsilon_x = [eps^0] D_x (defect cycles only); F_x is tabulated separately (tables: key F) and the pair gates are reported in both readings.'),
        dict(id='B-F4', severity='info',
             statement='Strict-route family spectra (null-mode defect cycles, no use of 4 lc + C as an input): F3 twist two S\' <= 6 and F4 fermion sl(2) S\' <= 5 have Spec D_1 = -32 h(j) u (scalars, even j) and -32 h(j+1) u (fermions, odd j) on the sl(2) primaries, F2 su(1|1) L <= 6 has D_1 = -8 u sum(1 - Pi) as matrices; on every sector the special-edge class equals 4 lc and the line class equals the P7 (c) target, i.e. Thm op-a holds family by family as an identity between two independent evaluations.'),
        dict(id='B-F2', severity='minor',
             statement='Placements along P_b, M_w (spin part), R at the special edge cancel identically configuration by configuration (the two configurations with S^(1)(x) on either side of the special edge carry opposite Koszul signs and equal matrices since x. commutes with Q^GF and H); along these generators Upsilon_x is the line class alone, i.e. Upsilon_{P_b} = C_{P_b}, Upsilon^0_M = C^0_M and Upsilon_R = 0 pair by pair on every computed row.'),
        dict(id='B-F3', severity='minor',
             statement='The reduced pole identity [eps^-m] D^0_Delta = -4m ct_m, [eps^-m] D_{P,M,R} = 0 (scaling/covariance of the counterterm map) holds on the rows where -ct^(1)_{x,m}(sO), s ct^(1)_{x,m}(O) and the vacuum-pole x spectator terms vanish; rows with sO containing a pair with arity-one poles (e.g. (A_V, lambda), sO ⊃ (dc) lambda) need the full eq c11-pole (open).')]
    # ---------------- OP0
    op0 = load('op0.json')
    if op0:
        res['OP0'] = dict(gates=[dict(id=g['id'], pass_=g['pass_'], value=str(g['value'])[:400]) for g in op0['gates']])
        for g in op0['gates']:
            gates.append(dict(id=g['id'], value=str(g['value'])[:300], target=g['target'], pass_=g['pass_']))
    else:
        res['OP0'] = 'not run'
    # ---------------- OP1/OP3 tables
    for n in (0, 1):
        ts = load(f'tables_summary_n{n}.json')
        if not ts:
            continue
        rows = ts['summary']
        idx = []
        for e in rows:
            idx.append({k: e[k] for k in ('pair', 'legset', 'n', 'outs', 'has', 'D', 'G_noLog', 'G_D1pair', 'G_D1pair_with_theory_F', 'G_pole', 'finite_special_vertex_nonzero')}
                       | {'evaluation': ('symbolic' if e['pol'] == 'sym' else e['pol'][:300]), 'file': f"out.nosync/tables/{e['pair'].replace('|', '_')}__{e['legset'].replace(' ', '_').replace('(', '').replace(')', '')}__n{n}.json"})
        res[f'OP1_tables_n{n}'] = dict(rows=idx, note='kernels (exact strings) per generator (P: b symbolic; D0: Upsilon^0_Delta; M0: w symbolic, w_mn = -w_nm; R: R_MN symbolic; K0: b symbolic) and placement class (line = C_x, special), poles (keys 1, 2), transcendental bases (1, log2, log3); F: theory finite special vertex (not added, finding B-F1); LC: ct_m (keys 1, 2) and lc (key log) of the pair; in the file of each row')
        full = [e for e in rows if all(e['has'].values())]
        nl = sum(1 for e in rows if not e['G_noLog'])
        gates.append(dict(id=f'G-noLog (n={n})', pass_=nl == 0, value=f'{len(rows) - nl}/{len(rows)} rows without log eps', target='no log eps in any defect cycle (structural: bounded special parameter)'))
        for g in ('D0', 'P', 'M0', 'R'):
            ok = [e for e in full if g in e['G_D1pair']]
            bad = [f"{e['pair']} [{e['legset']}]" for e in ok if not e['G_D1pair'][g]]
            gates.append(dict(id=f'G-D1pair {g} (n={n})', pass_=not bad,
                              value=f'{len(ok) - len(bad)}/{len(ok)} complete rows pass' + (f'; failing: {bad}' if bad else ''),
                              target={'D0': 'Upsilon^0_Delta = 4 lc + C^0_Delta', 'P': 'Upsilon_{P_b} = C_{P_b}', 'M0': 'Upsilon^0_{M_w} = C^0_{M_w}', 'R': 'Upsilon_R = 0'}[g] + ' pair by pair (Conj op-pair), Upsilon = [eps^0] D'))
            wF = [e for e in ok if e['finite_special_vertex_nonzero']]
            if wF:
                badF = [f"{e['pair']} [{e['legset']}]" for e in wF if not e['G_D1pair_with_theory_F'].get(g, True)]
                gates.append(dict(id=f'G-D1pair {g} with the theory finite vertex F_x (n={n}, rows with F != 0)', pass_=not badF,
                                  value=f'{len(wF) - len(badF)}/{len(wF)} pass' + (f'; failing: {badF}' if badF else ''),
                                  target='same identity with Upsilon = [eps^0] D + F_x (Def op-defect (iii)); finding B-F1'))
            pl = [e for e in rows if g in e['G_pole']]
            badp = [f"{e['pair']} [{e['legset']}]" for e in pl if not e['G_pole'][g]]
            gates.append(dict(id=f'G-pole reduced form {g} (n={n})', pass_=not badp,
                              value=f'{len(pl) - len(badp)}/{len(pl)} rows' + (f'; not satisfied (omitted -ct^(1)(sO) terms, finding B-F3): {badp}' if badp else ''),
                              target='[eps^-m] D^0_Delta = -4m ct_m; [eps^-m] D_{P,M,R} = 0'))
    # ---------------- F0, F1, conj, pos, OP8
    for name, key in (('f0.json', 'OP6_F0_even'), ('f0_odd.json', 'OP6_F0_odd_and_arity0'), ('fam_f1.json', 'OP6_F1'),
                      ('conj.json', 'OP3_second_evaluation'), ('pos.json', 'OP4_Kb_decomposition'), ('op8_kt.json', 'OP8_Kb_finite_T'),
                      ('r7.json', 'OP7_relations')):
        d = load(name)
        if d is None:
            res[key] = 'not run'
            continue
        res[key] = d.get('results', d)
        for g in d.get('gates', []):
            gates.append(dict(id=g['id'], value=str(g['value'])[:600], target=g['target'], pass_=g['pass_']))
        if name == 'fam_f1.json':
            gates.append(dict(id='G-spec F1 (su(2) sector)', pass_=d['gate_density'],
                              value='planar density of D_1 at theta = -1/3 equals -8 (1 - P) on {ZZ, ZX, XZ, XX}' if d['gate_density'] else 'density differs',
                              target='Spec D_1 = -8 u E_H on tr(Z^{L-M} X^M) (Prop c11-targets, [B03c]); spectra for L <= 8 listed'))
    # ---------------- F2, F3, F4 strict route (null-mode defect cycles)
    for fn, key, fam in ((sorted(glob.glob(os.path.join(OUT, 'sl2_Z_S*.json'))), 'OP6_F3_twist_two', 'F3'),
                         (sorted(glob.glob(os.path.join(OUT, 'sl2_psi_S*.json'))), 'OP6_F4_fermion_sl2', 'F4'),
                         (sorted(glob.glob(os.path.join(OUT, 'su11_L*.json'))), 'OP6_F2_su11', 'F2')):
        if not fn:
            res[key] = 'not run'
            continue
        fn = max(fn, key=lambda f: int(''.join(ch for ch in os.path.basename(f).split('_')[-1] if ch.isdigit())))
        d = json.load(open(fn))
        res[key] = d
        sect = d.get('matrices') or d.get('sectors')
        sect = {k: v for k, v in sect.items() if v.get('dim')}
        spec_key = 'G_spec_pass' if fam != 'F2' else 'G_spec_matrix_identity'
        bad = [k for k, v in sect.items() if not v[spec_key]]
        gates.append(dict(id=f'G-spec {fam} (strict: D_1 from the defect cycles, theta = -1/3)', pass_=not bad,
                          value=f'{len(sect) - len(bad)}/{len(sect)} sectors' + (f'; failing {bad}' if bad else '') + f' ({os.path.basename(fn)})',
                          target=('Spec D_1 = -8 u E_H, E_H = 4 h(j) (scalars) / 4 h(j+1) (fermions) on the sl(2) primaries' if fam != 'F2'
                                  else 'D_1 = -8 u sum_l (1 - Pi_{l,l+1}) (graded permutation) as a matrix on the graded-cyclic basis')))
        for gk, tg in (('G_D1pair_family', 'special-edge class = 4 lc (Conj op-pair, summed on the family)'),
                       ('P7c_line_equals_target', 'line class = sum C^0_Delta = P7 (c) target'),
                       ('P7b_line_equals_minus_third_dn', 'line class of Upsilon_{P_n} = -1/3 d_n (P7 (b))'),
                       ('Upsilon_Pn_special_zero', 'special class of Upsilon_{P_n} = 0 (finding B-F2)')):
            badk = [k for k, v in sect.items() if not v.get(gk, False)]
            gates.append(dict(id=f'{fam}: {gk}', pass_=not badk, value=f'{len(sect) - len(badk)}/{len(sect)} sectors' + (f'; failing {badk}' if badk else ''), target=tg))
        if d.get('n_nonrational'):
            gates.append(dict(id=f'{fam}: transcendental parts', pass_=False, value=str(d['n_nonrational']) + ' nonrational planar entries', target='none'))
    # ---------------- row inventory
    allrows = glob.glob(os.path.join(OUT, 'rows', '*.json'))
    fails = glob.glob(os.path.join(OUT, 'rows', '*.fail*'))
    res['row_inventory'] = dict(rows=len(allrows), failed_markers=sorted(os.path.basename(f) for f in fails))
    res['open'] = [
        'G-pole in the full form of eq c11-pole (terms -ct^(1)_{x,m}(sO), s ct^(1)_{x,m}(O), vacuum-pole x spectator) on the rows where the reduced form does not apply',
        'tables along xi (both chiralities) and eta: rows with <= 2 outputs (xi) / <= 1 output (eta) at n = 0 only; three-output rows of the even generators at fixed numeric polarisations',
        'n = 1 rows only on the family pairs with <= 2 outputs (route B cutoff n <= 1)',
        'families F5-F7 (spectra), F9 (twist), Thm op-c (K = P^dagger, Gram_1), B_T(K_n), Lemma op-Cid (composite loops of I), Prop op-typeI, arity-two components (R2-R6 two-sided): not computed in route B (F2 L <= 6, F3 S\' <= 6, F4 S\' <= 5 are computed, strict route, OP6_F2/F3/F4)',
        'missing pair rows (row_inventory.failed_markers without a row file): even-generator n = 0 three-output EP rows of the lambda pairs (TIMEOUT 2400 s at numeric letters) and n = 1 lambda-pair rows above the memory cap; their leg sets are reported incomplete (no G-D1pair verdict)',
        'memory: two failing row attempts of the 2026-11-03 resume reached ru_maxrss 1501 and 1521 MB before the 1450 MB watchdog (0.3 s polling) stopped them; the queue was restarted at cap 1250 / watchdog 1400 MB; no completed row exceeded 1.1 GB',
        'one-body terms along eta and K_b: no action-side particular solution available; carried as symbols']
    json.dump(res, open(os.path.join(ROOT, 'results.json'), 'w'), indent=1, default=str)
    json.dump(dict(gates=gates, counts=dict(passed=sum(1 for g in gates if g['pass_']), failed=sum(1 for g in gates if not g['pass_']))),
              open(os.path.join(ROOT, 'gates.json'), 'w'), indent=1, default=str)
    import shutil
    fin = os.path.join(OUT, 'final')
    os.makedirs(fin, exist_ok=True)
    for f in ('results.json', 'gates.json'):
        shutil.copy(os.path.join(ROOT, f), os.path.join(fin, f))
    nm = os.path.join(HERE, 'notes_master.txt')
    if os.path.exists(nm):
        shutil.copy(nm, os.path.join(fin, 'notes.txt'))
    print('gates', sum(1 for g in gates if g['pass_']), 'pass,', sum(1 for g in gates if not g['pass_']), 'fail')
    for g in gates:
        if not g['pass_']:
            print('  FAIL', g['id'], str(g['value'])[:200])


if __name__ == '__main__':
    main()
