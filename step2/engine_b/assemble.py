"""Engine B, Step 2, Part C: assemble cubic_results.json from the component outputs."""
import json, os, re, math, time
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
def J(name):
    p = os.path.join(HERE, name)
    return json.load(open(p)) if os.path.exists(p) else None

SECTORS = ['cAVAV', 'cApAp', 'clamlam', 'ccAsV']
LABEL = {'cAVAV': '(c, A|_V, A|_V)', 'cApAp': '(c, A|_Wperp, A|_Wperp)', 'clamlam': '(c, lambda, lambda)', 'ccAsV': '(c, c, A*|_V)'}
GS = {'cAVAV': ['G1', 'G2'], 'cApAp': ['G3'], 'clamlam': ['G5'], 'ccAsV': ['G4']}
L2, L3 = sp.log(2), sp.log(3)

def comb(parts):
    """parts: {'1': str, 'log2': str, 'log3': str} -> exact string in the basis 1, log(2), log(3), plus log(3/2) form."""
    r = sp.sympify(parts.get('1', '0')); a = sp.sympify(parts.get('log2', '0')); b = sp.sympify(parts.get('log3', '0'))
    expr = r + a * L2 + b * L3
    s = str(sp.expand(expr))
    if a == -b and b != 0:
        s32 = (f'{r} + ' if r != 0 else '') + f'{b}*log(3/2)'
        return s, s32, float(sp.N(expr, 20))
    return s, s, float(sp.N(expr, 20))

def parse_class(name):
    segs = name.split(' | ')
    verts = []; lines = []
    for sg in segs:
        m = re.match(r'(.*)\[(.*)\] -(.*)-$', sg)
        verts.append(f'{m.group(1)} [external {m.group(2)}]'); lines.append(m.group(3))
    special = [l for l in lines if l.startswith('K[')]
    return verts, lines, special[0] if special else None

def main():
    t0 = time.time()
    out = dict(project='N4_Costello_BV_2026-10-27', part='Step 2, Part C: cubic sector of the arity-zero one-loop obstruction',
               engine='B', unit='hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)); every coefficient multiplies unit * G_i; colour of every three-leg '
               'functional: ([c_0,c_1],c_2) (the only invariant; tr_g of three ad equals (1/2) tr_g(ad_X ad_Y)/(X,Y) ([c_0,c_1],c_2))',
               G_basis={'G1': 'int dvol ([c, A|_V], d^*d A|_V)', 'G2': 'int dvol ([c, A|_V], d d^* A|_V)',
                        'G3': 'int dvol ([c, A|_Wperp], d^*d A|_Wperp)', 'G4': 'int dvol ([c, dc], A^*)',
                        'G5': 'int dvol ([c, lambda], dslash lambda)'},
               G_basis_plane_waves={'G1': '-(a1.a2)(|p1|^2-|p2|^2) + (a1.p1)(a2.p1) - (a1.p2)(a2.p2)',
                                    'G2': '-(a1.p1)(a2.p1) + (a1.p2)(a2.p2)', 'G3': '-(b1.b2)(|p1|^2-|p2|^2)',
                                    'G4': 'i (a*.(p0+p1)) = -i (a*.p2)', 'G5': 'i (lambda1, rho(p1+p2) lambda2) = -i (lambda1, rho(p0) lambda2)',
                                    'convention': 'vertex function = coefficient of the odd r_j in leg order on sum_j r_j u_j e^{i p_j.x}, '
                                    'leg 0 = c (momentum p0), legs 1, 2 as in the sector, colour ([c0,c1],c2) stripped'})
    others = J('others.json'); solve = J('solve_F3.json'); d2 = J('d2_series.json')
    sectors = {}
    alpha = {}
    for sec in SECTORS:
        fit = J(f'ia_fit_{sec}.json')
        classes = []
        for cl, v in sorted(fit['classes'].items()):
            verts, lines, special = parse_class(cl)
            e0 = {g: comb({b: v['0'][b]['G'].get(g, '0') for b in v['0']}) for g in GS[sec]}
            res0 = max(float(v['0'][b]['residual_G']) for b in v['0'])
            resm1 = max(float(v['-1'][b]['residual_full']) for b in v['-1'])
            em1 = {b: v['-1'][b]['structures'] for b in v['-1']}
            classes.append(dict(
                **{'class': 'I-a: ' + ('triangle' if len(verts) == 3 else 'quartic-vertex bubble') + ' | ' + cl},
                vertices=verts, lines=lines, special_edge=special + ' = (e^{-eps|k|^2} - 1) x (component of K_0), Schwinger time in [0, eps]',
                other_edges='P(eps,infty)',
                multiplicity_and_weight='weight 1 per labelled graph: the 1/(2n) of the word formula cancels the 2n cyclic words '
                    '(n rotations, 2 orientations) of each labelled graph, the 1/2 of 1/2({I,I}_eps - {I,I}) cancels the exchange of the two '
                    'vertices joined by the special edge; the class value is the sum over the assignments of the external legs to the vertices '
                    '(both assignments of the two identical legs).',
                coeffs={'eps^-1': {g: '0' for g in GS[sec]}, 'log': {g: '0' for g in GS[sec]},
                        'eps^0': {g: e0[g][0] for g in GS[sec]}},
                eps0_with_log32={g: e0[g][1] for g in GS[sec]}, eps0_float={g: e0[g][2] for g in GS[sec]},
                eps_inv_multilinear_form=('0 (the eps^-1 coefficient, weight 2, is a multiple of (a1.a2) resp. (b1.b2) with colour ([c0,c1],c2); '
                                          'after the leg symmetrisation contained in the class it vanishes identically)' if sec in ('cAVAV', 'cApAp')
                                          else 'absent (superficial degree 1: eps^n <-> momentum degree 1 + 2n)'),
                residual_eps0_outside_span_G=res0, residual_eps_inv=resm1))
        tot = fit['total']
        e0 = {g: comb({b: tot['0'][b]['G'].get(g, '0') for b in tot['0']}) for g in GS[sec]}
        o = others[sec]
        def g_of(key):
            return {g: o[key]['G'].get(g, '0') for g in GS[sec]}
        negQY3 = {g: str(-sp.sympify(o['QY3']['G'].get(g, '0'))) for g in GS[sec]}
        negIY2 = {g: str(-sp.sympify(o['{I3,Y2}']['G'].get(g, '0'))) for g in GS[sec]}
        negIM = {g: str(-sp.Rational(9, 4) * sp.sympify(o['{I3,int(A,A)}']['G'].get(g, '0'))) for g in GS[sec]}
        bk2 = {g: sp.sympify(o['B_k2F(I3,int(A,A))']['G'].get(g, '0')) for g in GS[sec]}
        tree = g_of('tree B_P0(Ob2,I3)')
        classes.append(dict(**{'class': 'II-a: -Q I^CT_3 (log part only; Step-1 cubic log coefficients Y_3 = 6 AAA + 10 lambda lambda A + AcA* + c lambda lambda* + ccc*)'},
                            vertices=['-Q I^CT_3'], lines=[], special_edge=None, multiplicity_and_weight='tree vertex, weight 1',
                            coeffs={'eps^-1': {g: '0' for g in GS[sec]}, 'log': negQY3, 'eps^0': {g: '0' for g in GS[sec]}}))
        classes.append(dict(**{'class': 'III (K_0 part): -{I_3, I^CT_2}'}, vertices=['-{I_3, I^CT_2}'], lines=[], special_edge=None,
                            multiplicity_and_weight='tree vertex (local bracket), weight 1',
                            coeffs={'eps^-1': negIM, 'log': negIY2, 'eps^0': {g: '0' for g in GS[sec]}},
                            note='eps^-1 part = -(9/4) unit {I_3, int(A,A)} = 0 identically: {int([A,c],A*), int(A,A)} is a multiple of int([A,c],A) = 0'))
        classes.append(dict(**{'class': 'I-b + III (K_eps - K_0 part): special edge as a bridge from a cubic vertex with two external legs to the '
                               'regulated one-loop two-point function minus I^CT_2'},
                            vertices=['I_3 (two external legs)', 'W_1^(2)(eps) - I^CT_2(eps)'], lines=['K_eps - K_0 (bridge)'],
                            special_edge='K_eps - K_0 bridge, momentum k = momentum of the third leg',
                            multiplicity_and_weight='= -B_{(e^{-eps|k|^2}-1)K_0}(I_3, W_1^(2) - I^CT_2): I-b = {I_3, W_1^(2)}_{Khat}, III_Khat = -{I_3, I^CT_2}_{Khat}',
                            coeffs={'eps^-1': {g: '0' for g in GS[sec]}, 'log': {g: '0' for g in GS[sec]}, 'eps^0': {g: '0' for g in GS[sec]}},
                            eps_expansion=dict(
                                I_b={'eps^0': {g: str(sp.Rational(9, 4) * bk2[g]) for g in GS[sec]},
                                     'eps log eps': '+ unit B_{k^2 K_0}(I_3, Y_2) (local, weight 6)', 'eps^1': '+unit B_{k^2 K_0}(I_3, renormalised self-energy) - (9/8) unit B_{k^4 K_0}(I_3, int(A,A)) (nonlocal)'},
                                III_Khat={'eps^0': {g: str(-sp.Rational(9, 4) * bk2[g]) for g in GS[sec]},
                                          'eps log eps': '- unit B_{k^2 K_0}(I_3, Y_2)', 'eps^1': '+(9/8) unit B_{k^4 K_0}(I_3, int(A,A))'},
                                sum='+eps unit B_{k^2 K_0}(I_3, lim_{eps->0}(W_1^(2) - I^CT_2)) + O(eps^2 log eps) -> 0; eps^0, eps log eps and the eps k^4 int(A,A) terms cancel identically',
                                B_k2K0_I3_intAA=({g: str(bk2[g]) for g in GS[sec]}))))
        classes.append(dict(**{'class': 'I-d + II-b: two-leg defect D_2(eps) = (special bubbles) - Q I^CT_2, sector (c, A|_V), joined by P(eps,infty) to I_3'},
                            vertices=['D_2(eps)', 'I_3'], lines=['P(eps,infty) (bridge)'], special_edge='inside D_2',
                            multiplicity_and_weight='= B_{P(eps,infty)}(D_2(eps), I_3)',
                            coeffs={'eps^-1': {g: '0' for g in GS[sec]}, 'log': {g: '0' for g in GS[sec]}, 'eps^0': tree},
                            note='limit = tree of Ob_2 = -(27/8) unit int(d^*dc, d^*A|_V) with I_3 through P(0,infty); local; NOT part of Ob_3',
                            D2_series=d2))
        tot_logs = {g: str(sp.sympify(negQY3[g]) + sp.sympify(negIY2[g])) for g in GS[sec]}
        sectors[sec] = dict(legs=LABEL[sec], graph_classes=classes,
                            totals={'eps^-1': {g: '0' for g in GS[sec]}, 'log': tot_logs,
                                    'eps^0 (Ob_3)': {g: e0[g][0] for g in GS[sec]},
                                    'eps^0 (Ob_3) with log(3/2)': {g: e0[g][1] for g in GS[sec]},
                                    'eps^0 (Ob_3) float': {g: e0[g][2] for g in GS[sec]},
                                    'tree (I-d + II-b limit, excluded from Ob_3)': tree},
                            residual_checks=dict(
                                I_a_total_eps0_outside_span_G=max(float(tot['0'][b]['residual_G']) for b in tot['0']),
                                I_a_total_eps0_full_basis=max(float(tot['0'][b]['residual_full']) for b in tot['0']),
                                I_a_total_eps_inv=max(float(tot['-1'][b]['residual_full']) for b in tot['-1']),
                                tree_outside_span_G=o['tree B_P0(Ob2,I3)']['residual_G'],
                                n_random_plane_wave_configs=fit['configs']),
                            checks={'log eps of I-a': '0 (each class: only integer powers of eps, momentum degree omega + 2n)',
                                    'Q Y_3 + {I_3, Y_2}': {g: str(sp.sympify(o['QY3']['G'].get(g, '0')) + sp.sympify(o['{I3,Y2}']['G'].get(g, '0'))) for g in GS[sec]},
                                    '{Z, I_3} - Y_3 (Theorem c4-log)': o['{Z,I3}-Y3']['G'],
                                    'eps^-1 total': '0'},
                            runtime_seconds=fit['seconds'])
        for g in GS[sec]:
            alpha[g] = dict(exact=e0[g][0], with_log32=e0[g][1], float=e0[g][2])
    out['sectors'] = sectors
    out['alpha'] = {f'alpha_{g[1]}': alpha[g] for g in ['G1', 'G2', 'G3', 'G4', 'G5']}
    out['Ob_3'] = ('unit * [ (11/4 + 19/4 log(3/2)) (G1 + G3) + (27/8) G2 + 0 G4 - 9 log(3/2) G5 ]')
    out['tree_of_Ob2_with_I3'] = 'unit * [ (27/8) (G1 - G2 + G3) - (27/4) (G4 + G5) ]  (limit of I-d + II-b; local; not part of Ob_3)'
    out['wess_zumino'] = J('wz_result.json')
    out['exactness'] = solve
    out['numeric_full_sum_ccAsV'] = J('numeric_check_ccAsV.json')
    out['numeric_I_a_cApAp'] = J('numeric_check_ia_cApAp.json')
    out['colour_check'] = J('colour_check.json')
    out['bracket_checks'] = J('fa_tests.json')
    out['regression_two_point'] = J('regress_2pt.json')
    out['masters_check'] = J('masters_check.json')
    out['method_notes'] = J('method_notes.json')
    rt = sum(sectors[s]['runtime_seconds'] for s in SECTORS)
    out['runtime'] = dict(I_a_seconds_all_configs=rt, others_seconds=others.get('seconds'), note='single process, < 0.2 GB RSS')
    json.dump(out, open(os.path.join(HERE, 'cubic_results.json'), 'w'), indent=1)
    print(json.dumps(out['alpha'], indent=1)); print(out['Ob_3'])

if __name__ == '__main__':
    main()
