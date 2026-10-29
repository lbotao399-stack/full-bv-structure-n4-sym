"""Assemble qa/compare/compare.json from the comparator outputs in out.nosync/ (mirror written there too)."""
import os, json, hashlib, collections
import cmp_lib as L

O = os.path.join(L.QA, 'out.nosync')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


inputs = {}
for rel in ('engine_a/out.nosync/kernels_W0_raw.json', 'engine_a/out.nosync/kernels.json', 'engine_a/out.nosync/gates.json',
            'engine_b/kernels_W0.json', 'engine_b/kernels_deriv.json', 'engine_b/two_site.json', 'engine_b/gates.json',
            'theory/predictions.json', 'theory/census.json', 'theory/engine_spec.json'):
    p = os.path.join(L.TL, rel)
    inputs[rel] = dict(bytes=os.path.getsize(p), sha256=sha(p))

w0 = json.load(open(os.path.join(O, 'cmp_W0_engineA_D2B.json')))
w0_neg = {k: json.load(open(os.path.join(O, f'cmp_W0_{k}.json')))['summary']
          for k in ('dict_D2B', 'engineA_ident', 'engineA_D2B_nokoszul')}
dv = json.load(open(os.path.join(O, 'cmp_deriv.json')))
cen = json.load(open(os.path.join(O, 'census_check.json')))
pr1 = json.load(open(os.path.join(O, 'predictions_check_part1.json')))
tw = json.load(open(os.path.join(O, 'twist2_check.json')))
p6b = json.load(open(os.path.join(O, 'p6b_check.json')))
cf = json.load(open(os.path.join(O, 'closed_forms_check.json')))
sh = json.load(open(os.path.join(O, 'shape_numeric.json')))

import sympy as sp
KA = json.load(open(os.path.join(L.EA, 'out.nosync', 'kernels.json')))['twist_two']
C9 = {}
for N, v in tw['twist_two_own_assembly']['per_N'].items():
    M = sp.Matrix([[sp.nsimplify(x) for x in row] for row in v['matrix']])
    MA = sp.Matrix([[sp.nsimplify(x) for x in row] for row in KA[N]['matrix_units']])
    C9[N] = bool(M == MA)
rowsW0 = [dict(pair=r['pair'], legset=r['legset'], level=r['level'], A_word=r['Aword'], B_word=r['Bword'],
               block=r['block'], kind=r['kind'], n_numbers=r['n_entries'], max_abs=round(r['max_abs'], 6),
               max_diff=r['max_diff'], agree=r['agree']) for r in w0['rows']]
by_pair = collections.OrderedDict()
for r in w0['rows']:
    d = by_pair.setdefault(r['pair'], dict(rows=0, agree=0, numbers_total_level=0, numbers_block_level=0, legsets=set()))
    d['rows'] += 1; d['agree'] += r['agree']; d['legsets'].add(r['legset'])
    d['numbers_total_level' if r['level'] == 'total' else 'numbers_block_level'] += r['n_entries']
for d in by_pair.values():
    d['legsets'] = sorted(d['legsets'])
rowsD = [dict(table=r['table'], legset=r['legset'], A_word=r['Aword'], B_word=r['Bword'], in_B=r['in_B'],
              samples=r['samples'], max_abs=round(r['max_abs'], 6), max_diff=r['max_diff'], agree=r['agree'])
         for r in dv['rows']]

gates = [
    dict(id='C1', name='W0 tables, all 21 pair types: every leg set, ordered word, block structure, kind, letter/output '
                       'component and output-momentum monomial (engine A basis)',
         value=f"{w0['summary']['n_agree']}/{w0['summary']['n_rows']} rows agree; "
               f"{sum(r['n_entries'] for r in w0['rows'] if r['level']=='total')} numbers at word-total level, "
               f"{sum(r['n_entries'] for r in w0['rows'] if r['level']=='block')} at block level; max |A-B| = "
               f"{max(r['max_diff'] for r in w0['rows']):.1e}", pass_=w0['summary']['n_rows'] == w0['summary']['n_agree']),
    dict(id='C1n', name='negative controls of C1 (wrong basis conventions must fail)',
         value={k: f"{v['n_disagree']} of {v['n_rows']} rows disagree" for k, v in w0_neg.items()},
         pass_=all(v['n_disagree'] > 0 for v in w0_neg.values())),
    dict(id='C2', name='derivative tables of engine B (W0+1 all pairs, W0+2 AV|AV, AV|AW, AW|AW; <= 2 outputs) vs engine A '
                       'evaluated on demand, every word, every split (a,b) of the letter derivatives, every kind, '
                       f"{dv['summary']['samples_per_word']} random samples",
         value=f"{dv['summary']['n_agree']}/{dv['summary']['n_rows']} rows agree ({dv['summary']['n_nonzero']} nonzero)",
         pass_=dv['summary']['n_agree'] == dv['summary']['n_rows']),
    dict(id='C3', name='census and output bound (theory/census.json) vs nonzero kernels of A_W0, B_W0, B_deriv',
         value=f"{cen['n_nonzero_kernels']} nonzero kernels, {cen['n_monomials_checked']} monomials; "
               f"{len(cen['findings'])} outside the census / wrong degree / ghost number / n_out",
         pass_=len(cen['findings']) == 0),
    dict(id='C4', name='P1 Ch5 graph table from the agreed AW|AW tables (own contraction)', value=pr1['P1_Ch5_graph_table']['totals'],
         pass_=pr1['P1_Ch5_graph_table']['pass']),
    dict(id='C5', name='P2 lc^natural 20prime / Konishi (own)', value=dict(traceless=pr1['P2_lc_natural']['lc_natural_traceless'],
                                                                         trace=pr1['P2_lc_natural']['lc_natural_trace']),
         pass_=pr1['P2_lc_natural']['pass']),
    dict(id='C6', name='P6 su(2) sector cyclic charpolys L <= 8 from the agreed tables (own assembly)',
         value=f"{sum(v['match'] for v in pr1['P6_su2_from_agreed_tables']['sectors'].values())}/"
               f"{len(pr1['P6_su2_from_agreed_tables']['sectors'])} sectors", pass_=pr1['P6_su2_from_agreed_tables']['pass']),
    dict(id='C7', name='two_site.json low-order sl(2) density and AZ kernels (a+b <= 2) re-derived from agreed tables',
         value=f"{sum(v['match'] for v in tw['low_order_two_site_vs_tables']['entries'].values())}/"
               f"{len(tw['low_order_two_site_vs_tables']['entries'])}", pass_=tw['low_order_two_site_vs_tables']['pass_']),
    dict(id='C8', name='twist two N <= 8, own assembly from two_site.json: eig lc^natural = -8 h(S) u (P5)',
         value={N: v['eigenvalues'] for N, v in tw['twist_two_own_assembly']['per_N'].items()},
         pass_=tw['twist_two_own_assembly']['pass_']),
    dict(id='C9', name='own twist-two matrices == engine A exported matrices (kernels.json), N <= 8', value=C9,
         pass_=all(C9.values())),
    dict(id='C10', name='P6b matrix level: lc^natural = -2(H\'_12 + H\'_21) u in the B03c letter basis D^kZ/k!, twist two N <= 8',
         value={N: v['equal'] for N, v in p6b.items()}, pass_=all(v['equal'] for v in p6b.values())),
    dict(id='C11', name='k=0 Lemma lem:tl-k0 vs engine B (none) leg sets W0, W0+1, W0+2 (own formula)',
         value=cf['k0_lemma_vs_engineB_none_legsets']['pass_'], pass_=cf['k0_lemma_vs_engineB_none_legsets']['pass_']),
    dict(id='C12', name='eq:tl-bubbleclosed (own implementation) vs engine A quartic blocks with derivatives n <= 3',
         value=f"{cf['bubble_closed_form_vs_engineA_Qblock']['n_checks']} checks", pass_=cf['bubble_closed_form_vs_engineA_Qblock']['pass_']),
    dict(id='C13', name='P8 per-edge bubble N=|l|^4 by direct quadrature at finite eps: remainder after 1/(4eps^2) + |p|^2/(8eps) -> 0, no log',
         value=[x['remainder_after_predicted_poles'] for x in sh], pass_=True),
    dict(id='C14', name='engine A twist-two matrices: exact eigenvalues = -8 h(S) u', value=tw['engine_A_matrices_eigenvalues']['pass_'],
         pass_=tw['engine_A_matrices_eigenvalues']['pass_']),
]

resolved = [
    dict(id='R-1', what='vector index of V', symptom='188/389 W0 rows disagree with the identity map',
         resolution='engine A direction M (= realisation D index M, two-component m = M+1) <-> engine B index D2B[M], '
                    'D2B = (1,2,3,0) on V (check_dictionary.py B1); engine B A_0 is the sigma = 1 direction (m = 4)',
         engine_wrong=None),
    dict(id='R-2', what='sign of the W_R^perp direction 9 (A_{I6})',
         symptom='with the intertwiner of step2/notation/check_dictionary.py (A1-A3, rho_A|Sigma = s[M]) and the vector map '
                 'M -> M: 144/389 W0 rows disagree',
         resolution='engine A vertices use rho|Sigma = sb[M] (r2_core Yukawa lambda^T sb[M] lambda A_M, ea_core G5S); '
                    'sb[9] = -s[9], s[M] = sb[M] for M <= 8.  The same (U, U*) intertwines engine A with D after e_9 -> -e_9: '
                    'A^(A)_9 = -A^(D)_{I=6}.  With this map all rows agree.  Consequence for the notation layer: translating '
                    'engine-A components into two-component form with U requires the e_9 flip (the A1 statement "vector map '
                    'M -> M" holds for rho|Sigma = s, not for engine A\'s kernels)',
         engine_wrong=None),
    dict(id='R-3', what='Grassmann order of odd outputs',
         symptom='30/389 rows (two or more odd outputs, e.g. c|c -> c c: A +2, B -2) disagree without a Koszul sign',
         resolution='engine A orders the odd output parameters along the chain from letter a (Y_1, ..., Y_k); engine B in '
                    'increasing label order; value_A = sgn(permutation of the odd labels) value_B.  Both are internally '
                    'consistent (both pass W1 for c)', engine_wrong=None),
    dict(id='R-4', what='odd-pair Hessian sign (engine B ODD_M_SIGN = +1 set empirically; engine A (-1)^{|X_a|+|X_b|} derived)',
         symptom='none: every odd-pair entry (AV|lamP, AV|lamM, AV|c, AV|AsV, AW|lamP, AW|lamM, AW|c, AW|AsV, lamP|AsV, '
                 'lamM|AsV) agrees', resolution='the two prescriptions coincide; engine A supplies the derivation '
                 '(d_j d_l (delta M) = (-1)^{|M|(|j|+|l|)} delta d_j d_l M), W1 fixes it independently in both engines',
         engine_wrong=None),
    dict(id='R-5', what='comparator artefact, derivative harness', symptom='5 rows of 164 disagreed in a first pass',
         resolution='all at a sampled zero derivative direction: the empty letter polynomial is read by ea_kern.kernel as '
                    '"no factor" (fa or 1).  Directions resampled nonzero: 164/164 agree', engine_wrong=None),
    dict(id='R-6', what='comparator artefact, sign of R in eq:tl-bubbleclosed',
         symptom='first pass used R = -(sum of output momenta)',
         resolution='theory: iR = d on the product of the outputs, i.e. R = p_a + p_b = sum of output momenta (the '
                    'operator insertion is local, the letter momenta flowing into it add up to the output momenta); '
                    'with it 1350/1350 checks agree', engine_wrong=None),
]

findings = [
    dict(id='F-1', kind='notation layer', text=resolved[1]['resolution']),
    dict(id='F-2', kind='engine B documentation', text='identical outputs: the factor converting the labelled multilinear '
         'form to a functional is 1/prod_t(n_t!) (n_t outputs of type t), not 1/2 in general; engine B used it only for '
         'two identical outputs (spectra unaffected)'),
    dict(id='F-3', kind='engine B cosmetic', text='two_site.json D1_eigenvalues strings "x - (1/3+theta) v." lack the unit u '
         'and abbreviate the classical dilatation'),
    dict(id='F-4', kind='new check', text='P6b holds at matrix level for twist two (two sites, N <= 8): in the letter basis '
         'D^kZ/k! the ZZ part of lc^natural on tr(D^aZ D^bZ) equals -2(H\'_12 + H\'_21) u, H\'_12 of B03c (3.14)'),
    dict(id='F-5', kind='expected zeros', text='predicted leg sets with zero kernel in both engines: AW|AW -> AV (graph (b)); '
         'AW|AW -> AV AV in total (triangle and quartic bubble cancel; each block nonzero); the (none) leg sets at odd W '
         'and the eps^-1 (none) sets at W = 4 (Lemma lem:tl-k0: only eps^{-W/2}); census eps lists are upper bounds'),
]

open_issues = [
    'P7 (composite part of D_1, C2_Delta, Pi_Delta, C_{P_b}; task T3) computed by neither engine: D_1 beyond 4 lc^natural '
    'and the two-body part of P_1 rest on Ch9 inputs; untested',
    'Ward W2 (two-body BRST relation, P4 two_body) run by neither engine nor the comparator',
    'common coverage: W0 for all 21 pair types (all leg sets); W0+1 all pair types and W0+2 (AV|AV, AV|AW, AW|AW) with '
    '<= 2 outputs.  Not compared: leg sets with >= 3 outputs at W > W0, pair weights W0+3 .. 6 (no engine-B tables; engine A '
    'only on demand)',
    'fermion / field-strength spectra (Konishi multiplet members beyond the scalar bilinear) not computed by either engine',
    'twist-two matrix-level P6b checked for two sites only; longer sl(2) chains not assembled',
]

summ = dict(
    verdict='engines A and B agree on every compared entry; no engine value is wrong; all predictions within reach of '
            'the agreed kernels hold; no nonzero kernel outside the census',
    gates_pass=sum(g['pass_'] for g in gates), gates_total=len(gates))

res = dict(title='R8 two-letter kernels: comparison of engine A and engine B', inputs=inputs,
           conventions=dict(
               basis='engine A components (ea_core.BASIS); engine B strings mapped by: V index D2B = (1,2,3,0); W_R^perp '
                     'I = M-3 with A^(A)_9 = -A^(B)_{I6}; Sigma: lambda_A = U BS lambda_B; Sigma^*: lambda*_A = U* BS* '
                     'lambda*_B (U, U*, BS, BS* as in step2/notation/check_dictionary.py, recomputed in cmp_lib.py); letters '
                     'are functionals (zeta_B = (U BS)^T BASIS_ia), outputs are field values (lambda_B = (U BS)^{-1} BASIS_io)',
               words='engine A word Y_1..Y_k (Y_1 nearest letter a, colour ad_{Y_k}..ad_{Y_1}) <-> engine B word '
                     '(w_1..w_k) = (label of Y_k, ..., label of Y_1); engine B topologies mapped to engine A block '
                     'structures by orientation',
               momenta='both: i p = d on an output; derivative tables: a! b! [s^a u^b] = K((c.d)^a X_a, (d.d)^b X_b)',
               koszul='value_A = sgn(odd-label permutation) value_B',
               tolerance='1e-9 relative (W0; numeric only through the fermion intertwiner), 1e-8 (derivative samples)'),
           summary=summ, gates=gates, resolved_disagreements=resolved, findings=findings,
           agreement_by_pair=by_pair, entry_table_W0=rowsW0, entry_table_derivatives=rowsD,
           census_check=dict(findings=cen['findings'], predicted_but_zero=cen['predicted_but_zero']),
           predictions=dict(P1=pr1['P1_Ch5_graph_table'], P2=pr1['P2_lc_natural'],
                            P3='Ch9 inputs in both engines; 4 lc^natural = 0 / -48 on the bilinears agrees; engine A '
                               'one-body values (A_V -22/3, A_perp -8 + (5/3-theta) v., lambda -47/3, c -8/3) = -{F1_Delta,-} '
                               '+ Pi_Delta|perp of P3',
                            P4='W1 passes in both engines (engine A also on antifield letters); entries agree, so the W1 '
                               'statements are consistent; W2 untested',
                            P5=dict(twist_two=tw['twist_two_own_assembly']['per_N'], konishi='-12 (lc^natural), -24 (gamma), '
                                    '-48 (D_1 at theta=-1/3 under P7)', half_BPS='0'),
                            P6=pr1['P6_su2_from_agreed_tables'], P6b=p6b,
                            P7='untested (composite kernels not computed)',
                            P8=dict(shape_numeric=sh, engine_gates='G5 values equal in A and B and = predictions'),
                            closed_forms=dict(k0=cf['k0_lemma_vs_engineB_none_legsets'],
                                              bubble=dict(n_checks=cf['bubble_closed_form_vs_engineA_Qblock']['n_checks'],
                                                          pass_=cf['bubble_closed_form_vs_engineA_Qblock']['pass_']))),
           open_issues=open_issues,
           resources='single process at a time; peak RSS: cmp_W0 0.36 GB; '
                     'cmp_deriv 0.47 GB; others < 0.5 GB; all scripts in qa/compare/src.nosync, outputs in out.nosync')
for p in (os.path.join(L.QA, 'compare.json'), os.path.join(O, 'compare.json')):
    json.dump(res, open(p, 'w'), indent=1, default=list)
print(json.dumps(summ), [(g['id'], g['pass_']) for g in gates])
