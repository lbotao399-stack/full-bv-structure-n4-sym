"""Assemble results.json, gates.json (complete/) from the fit outputs of both routes, with the entry-by-entry
agreement of route A (engine A code path: sa_core/GSimplex cycles, mode model) and route B (engine B code path:
eb_cycle word formula, eb_qb Q)."""
import json, os
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')
TOP = os.path.dirname(OUT)


def ld(name):
    fn = os.path.join(OUT, name)
    return json.load(open(fn)) if os.path.exists(fn) else None


def agree(a, b):
    if a is None or b is None:
        return 'missing'
    return 'agree' if a == b else 'DIFFER'


results, gates = {}, []


def gate(gid, route, ok, what, value=None):
    gates.append(dict(id=gid, route=route, pass_=bool(ok) if ok is not None else None, what=what, value=value))


# ---------------------------------------------------------------- odd pairs (xi, xi')
for pair, label in (('pm', '(xi_+, xi_-)'), ('pp', '(xi_+, xi_+\')'), ('mm', '(xi_-, xi_-\')')):
    rb, ra = ld(f'routeB_pair_{pair}_fit.json'), ld(f'routeA_pair_{pair}_fit.json')
    ent = {'blocks': {}}
    for bn in (rb or ra or {'blocks': {}})['blocks']:
        eb = rb['blocks'][bn] if rb else None
        ea = ra['blocks'][bn] if ra else None
        thb = eb['residual']['solution_theta'] if eb and eb['residual'] else None
        tha = ea['residual']['solution_theta'] if ea and ea['residual'] else None
        ent['blocks'][bn] = dict(
            routeB=None if not eb else {k: eb[k] for k in ('basis_dim', 'expected', 'closed_dim', 'exact_dim', 'H_minus1_dim', 'ob_in_span', 'pole_sum_zero', 'ob_closed')},
            routeA=None if not ea else {k: ea[k] for k in ('basis_dim', 'expected', 'closed_dim', 'exact_dim', 'H_minus1_dim', 'ob_in_span', 'pole_sum_zero', 'ob_closed')},
            conditions_B=thb, conditions_A=tha, agreement_conditions=agree(thb, tha))
        for r, e in (('B', eb), ('A', ea)):
            if e:
                gate('G-pole2 %s %s' % (label, bn), r, e['pole_sum_zero'], 'sum of the eps^-1 coefficients of the (C2) cycles = 0')
    for r, x in (('B', rb), ('A', ra)):
        if x:
            gate('G-WZ2 %s' % label, r, x['G_WZ2_full']['pass_'], x['G_WZ2_full']['rule'] + ' (corrected signs, finding F1)', dict(QOb2_zero=x['G_WZ2_full']['QOb2_zero']))
    results['Ob2' + label] = ent

for ch in ('p', 'm'):
    rb, ra = ld(f'routeB_xiP_{ch}_fit.json'), ld(f'routeA_xiP_{ch}_fit.json')
    eb, ea = (rb or {}).get('blocks', {}).get('all'), (ra or {}).get('blocks', {}).get('all')
    cb = eb['residual']['solution_theta'] if eb else None
    ca = ea['residual']['solution_theta'] if ea else None
    results['Ob2(xi_%s, P_b)' % ('+' if ch == 'p' else '-')] = dict(
        routeB=None if not eb else {k: eb[k] for k in ('basis_dim', 'closed_dim', 'exact_dim', 'H_dim', 'ob_in_span', 'ob_closed', 'pole_sum_zero')},
        routeA=None if not ea else {k: ea[k] for k in ('basis_dim', 'closed_dim', 'exact_dim', 'H_dim', 'ob_in_span', 'ob_closed', 'pole_sum_zero')},
        Ob2_value_B=None if not eb else eb['ob_coords'], conditions_B=cb, conditions_A=ca, agreement_conditions=agree(cb, ca))
    for r, x in (('B', rb), ('A', ra)):
        if x:
            gate('G-WZ2 (xi_%s,P_b)' % ch, r, x['G_WZ2']['QOb2_zero'], x['G_WZ2']['rule'])
            gate('G-zero (xi_%s,P_b)' % ch, r, bool(x['blocks']['all']['residual']['solution_theta']), '[Ob2(xi,P_b)] = 0 given uniform theta^P ([xi,P_b] = 0)')

for ch in ('p', 'm'):
    b, a = ld(f'B_dxi{ch}_fit.json'), ld(f'A_dxi{ch}_fit.json')
    results['Ob2(Delta, xi_%s)' % ('+' if ch == 'p' else '-')] = dict(
        routeB=None if not b else {k: b[k] for k in ('basis_dim', 'closed', 'exact', 'H_minus1', 'ob_in_span', 'pole_sum_zero', 'G_F1Delta', 'solvable')},
        routeA=None if not a else {k: a[k] for k in ('basis_dim', 'closed', 'exact', 'H_minus1', 'ob_in_span', 'pole_sum_zero', 'G_F1Delta', 'solvable')},
        conditions_B=b and b.get('conditions'), conditions_A=a and a.get('conditions'),
        F1_2_exact_part_B=b and b.get('psi_exact'), F1_2_exact_part_A=a and a.get('psi_exact'),
        agreement=agree(b and (b.get('conditions'), b.get('psi_exact')), a and (a.get('conditions'), a.get('psi_exact'))))
    for r, x in (('B', b), ('A', a)):
        if x:
            gate('G-pole2 (Delta,xi_%s)' % ch, r, x['pole_sum_zero'], 'sum a_1 = 0')
            gate('G-F1Delta', r, x['G_F1Delta'], 'Q F1part(Delta)|2 = Ob1(Delta)|2 (regression of T1)')

for ch in ('p', 'm'):
    b, a = ld(f'B_eta{ch}_fit.json'), ld(f'A_eta{ch}_fit.json')
    keys = ('basis0', 'pole_sum_zero', 'ob_in_span', 'G_eta_closed', 'closed0', 'exact0', 'H0', 'deg_m1_basis', 'H_m1_kernel', 'G_pos', 'G_WZ1', 'G_WZ1_printed_sign', 'R1_closed', 'solvable')
    results['F1(eta_%s)|2' % ('+' if ch == 'p' else '-')] = dict(
        routeB=None if not b else {k: b.get(k) for k in keys}, routeA=None if not a else {k: a.get(k) for k in keys},
        F1_particular_B=b and b.get('F1_particular'), F1_particular_A=a and a.get('F1_particular'),
        free_directions_B=b and b.get('F1_free_directions'), free_directions_A=a and a.get('F1_free_directions'),
        agreement=agree(b and (b.get('F1_particular'), b.get('F1_free_directions')), a and (a.get('F1_particular'), a.get('F1_free_directions'))))
    for r, x in (('B', b), ('A', a)):
        if x:
            gate('G-eta-closed (%s)' % ch, r, x['G_eta_closed'], 'Q S_2(eta) = 0')
            gate('G-pos eta (%s)' % ch, r, x.get('G_pos'), 'position part of Ob1(eta)|2 = Ob1(rho(e_mu) eta)|2 = (79/24 - L/2)(g1+g3) + 27/8 g2 (computed, not structural)')
            gate('G-WZ1 eta (%s)' % ch, r, x.get('G_WZ1'), 'Q Ob1(eta) = -{Ob0, S(eta)} (corrected sign; printed sign fails: %s)' % (not x.get('G_WZ1_printed_sign')))
            gate('G-pole eta (%s), intrinsic part' % ch, r, x['pole_sum_zero'], 'alpha-part of the eps^-1 coefficient vanishes on two legs (position part of the pole not stored)')

for ch in ('p', 'm'):
    b, a = ld(f'B_xieta{ch}_fit.json'), ld(f'A_xieta{ch}_fit.json')
    if b or a:
        keys = ('basis_dim', 'closed', 'exact', 'H_minus1', 'ob_in_span', 'ob_closed', 'pole_sum_zero', 'G_bracket', 'solvable')
        results['Ob2(xi_+, eta_%s)' % ('+' if ch == 'p' else '-')] = dict(
            routeB=None if not b else {k: b.get(k) for k in keys}, routeA=None if not a else {k: a.get(k) for k in keys},
            conditions_B=b and b.get('conditions'), conditions_A=a and a.get('conditions'),
            agreement=agree(b and b.get('conditions'), a and a.get('conditions')))
        for r, x in (('B', b), ('A', a)):
            if x:
                gate('G-bracket (xi,eta_%s)' % ch, r, x.get('G_bracket'), '[xi,eta] formula = exact fit of {S(xi),S(eta)} - Q S^(2) on two configurations')
                gate('G-pole2 (xi,eta_%s), intrinsic part' % ch, r, x.get('pole_sum_zero'), 'alpha-part of the eps^-1 coefficient = 0')

b, a = ld('B_Kb_fit.json'), ld('A_Kb_fit.json')
if b or a:
    keys = ('pole_sum_zero', 'ob_in_span', 'G_K_closed', 'basis0', 'closed0', 'exact0', 'H0', 'deg_m1_basis', 'closed_m1', 'G_WZ1', 'G_WZ1_printed_sign', 'R1_closed', 'solvable')
    results['F1(K_b)|2'] = dict(routeB=None if not b else {k: b.get(k) for k in keys}, routeA=None if not a else {k: a.get(k) for k in keys},
                                F1_particular_B=b and b.get('F1_particular'), F1_particular_A=a and a.get('F1_particular'),
                                free_B=b and b.get('F1_free_directions'), free_A=a and a.get('F1_free_directions'),
                                post_B=ld('B_Kb_post.json'), post_A=ld('A_Kb_post.json'),
                                agreement=agree(b and b.get('F1_particular'), a and a.get('F1_particular')))
    for r, x, pst in (('B', b, ld('B_Kb_post.json')), ('A', a, ld('A_Kb_post.json'))):
        if x:
            gate('G-K-closed', r, x.get('G_K_closed'), 'Q S_2(K_b) = 0 (validates the K_b vertex: weights, spin, L_v)')
            gate('G-WZ1 K_b', r, x.get('G_WZ1'), 'Q Ob1(K_b) = -{Ob0, S(K_b)} (corrected sign; printed sign fails: %s)' % (not x.get('G_WZ1_printed_sign')))
            gate('G-T K_b', r, True, 'structural: Lemma at-special; the eps^0 coefficients are computed as T-independent simplex coefficients (no numerical finite-T run)')
        if pst:
            gate('G-pole K_b', r, pst['G_pole'], 'sum a_1 = (9/4){int(A,A), S(K_b)} (eq at-pole), full canonical data incl. position parts')

gates += [dict(id='GB-S2', route='B', pass_=True, what='Q S^(2)_2 + S_2(P_Gamma|V) - {S_2(xi1),S_2(xi2)} = 0 on two legs, (+,-),(+,+),(-,-),(-,+) (t_s2gate.py); opposite bracket sign fails'),
          dict(id='GA-S2', route='A', pass_=True, what='same identity in the exact mode model (t_gateA.py), (+,-),(+,+),(-,-); opposite sign fails'),
          dict(id='GB-tree', route='B', pass_=True, what='density bracket: M_Sfree = q, M_T1 = engine-B map_T1, {N,Sfree} = QN, C-rule Q{N,F} = -{QN,F} - (-1)^|N|{N,QF} (t_tree*.py)'),
          dict(id='GA-dens/GA-Q/GA-dtree', route='A', pass_=True, what='densities = mode-model values, Q with qmat_A = mode-model Q, density bracket = mode-model bracket (even and odd G) (t_gateA.py, t_dtree.py)'),
          dict(id='G-bracket', route='B', pass_=True, what='[xi,eta] from eq c2-oddodd: S([xi,eta]) = {S(xi),S(eta)} - Q S^(2) fitted exactly onto c_D S(Delta) + S(M_omega) + S(R_r); c_D = (xi,eta)/2, omega = formula, r = -(transcribed formula) (t_bfit*.py)')]
json.dump(results, open(os.path.join(TOP, 'results.json'), 'w'), indent=1)
json.dump(gates, open(os.path.join(TOP, 'gates.json'), 'w'), indent=1)
print(len(results), 'result entries,', len(gates), 'gates;', sum(1 for g in gates if g['pass_'] is False), 'failing')
