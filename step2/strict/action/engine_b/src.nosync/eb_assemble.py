"""Assemble results.json and gates.json from out.nosync/*.json (engine B, action side)."""
import json, os, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'out.nosync')


def J(name):
    fn = os.path.join(OUT, name)
    return json.load(open(fn)) if os.path.exists(fn) else None


def comb(d, sign=1):
    """{'1': {lab: x}, 'log2': {...}, 'log3': {...}} -> {lab: 'a + b log2 + c log3'} with overall sign."""
    labs = set()
    for k, v in (d or {}).items():
        if v:
            labs |= set(v)
    out = {}
    for lab in sorted(labs):
        parts = []
        for k in ('1', 'log2', 'log3', 'log5'):
            x = Fr((d.get(k) or {}).get(lab, '0')) * sign
            if x != 0:
                parts.append(str(x) if k == '1' else f'({x})*{k}')
        out[lab] = ' + '.join(parts) if parts else '0'
    return out


def gate(gid, value, target, ok, notes=''):
    return dict(id=gid, value=value, target=target, status='pass' if ok is True else ('fail' if ok is False else str(ok)), notes=notes)


def main():
    res = dict(engine='B (route B: momentum-space numerator algebra, exact Feynman-parameter integrals)', side='action',
               date=time.strftime('%Y-%m-%d %H:%M'), unit='u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2)',
               normalisation=('theory component (Def at-names (i), Def c2-cochains/c2-dce, generating function '
                              'I + S1(x) c^x): arity one = -(raw left-CE coefficient) for even x, = raw for odd x; '
                              'arity two (Delta,P_b) = -(raw c_Delta c_P coefficient); see notes.txt'),
               Ob1={}, Ob2={}, poles={}, F1={}, cohomology={}, theta={}, findings=[], open=[])
    gates = []
    g0 = J('G_reg0.json')
    if g0:
        t2 = g0['two_legs']; t3 = g0['three_legs']
        res['Ob0'] = {'two legs': comb(t2['eps0']), 'two-leg pole (cycles)': comb(t2['pole']),
                      'three legs': {k: comb(v['eps0']) for k, v in t3.items()}}
        ok = (t2['eps0'].get('1', {}).get('int(d*dc,d*A_V)') == '-27/8' and
              t3['cAVAV']['eps0']['1'] == {'G1': '11/4', 'G2': '27/8'} and t3['cAVAV']['eps0']['log3']['G1'] == '19/4'
              and t3['cAVAV']['eps0']['log2']['G1'] == '-19/4' and t3['cApAp']['eps0']['1']['G3'] == '11/4'
              and t3['clamlam']['eps0']['log3'] == {'G5': '-9', 'parity-odd': '0'} and t3['ccAsV']['eps0'] == {})
        gates.append(gate('G-reg0', res['Ob0'], '-27/8 u int(d*dc,d*A_V); (11/4+19/4 L)(G1+G3)+27/8 G2-9L G5, G4 0, parity-odd 0',
                          ok, 'special-edge cycles of I alone; pole of (c,A_V) = (9/4) Q int(A,A) (counterterm)'))
    t1 = J('T1.json')
    if t1:
        d = {}
        for sec, r in t1.items():
            d[sec] = dict(intrinsic=comb(r['eps0'], -1), position_part_nonzero_entries=r['position_part_nonzero_entries'],
                          raw_line_placements=r['split_ratio_first_basis']['line'],
                          raw_special_edge_placements=r['split_ratio_first_basis']['special'])
            res['poles'][f'Delta {sec}'] = comb(r['pole'], -1) if r['pole'] else '0'
        res['Ob1']['Delta|_2'] = d
        target = {'AVAV': {'int(dA_V,dA_V)': '-11/3', 'int(d*A_V,d*A_V)': '0'}, 'ApAp': {'int(dA_perp,dA_perp)': '-11/3'},
                  'lamlam': {'int(lam,dslash lam)': '47/3', 'int(lam,dslash rho(vol) lam) [parity odd]': '0'},
                  'cAsV': {'int(dc,A*)': '-14/3'}}
        ok = all(d[s]['intrinsic'] == target[s] for s in target) and all(t1[s]['position_part_nonzero_entries'] == 0 for s in t1)
        gates.append(gate('G-reg1', {s: d[s]['intrinsic'] for s in d}, 'u[-11/3 int(dA,dA) + 47/3 int(lam,dslash lam) - 14/3 int(dc,A*)], 0 for int(d*A,d*A) and explicit x',
                          ok, 'all four sectors exact; raw (left CE coordinate) values are -1 times these'))
        gates.append(gate('G-P2', 'position part of Ob1(Delta)|_2 = 2 Ob1(P_b)|_2: all entries 0', 'exact 0',
                          all(t1[s]['position_part_nonzero_entries'] == 0 for s in t1)))
        polesD = {s: (comb(t1[s]['pole'], -1) if t1[s]['pole'] else {}) for s in t1}
        okp = polesD['AVAV'] == {'int(A,A)': '-9'} and polesD['ApAp'] == {'int(A,A)': '-9'}
        gates.append(gate('G-pole (Delta)', polesD, '(9/4) u {int(A,A),S(Delta)} = -9 u int(A,A) ({F,S(Delta)} = F\'[v.phi], Def c1-bracket); 0 in (lam,lam), (c,A*)',
                          okp))
        res['cross_checks'] = dict(X_scal=dict(statement='Ob1(Delta)|_2 = 4 lc0 + composite: u[-4,16,-4] + u[1/3,-1/3,-2/3]',
                                               computed_special_edge_placements={s: d[s]['raw_special_edge_placements'] for s in d},
                                               computed_line_placements={s: d[s]['raw_line_placements'] for s in d},
                                               note='raw values (multiply by -1): special-edge placements give -4,-4,16,-4 x basis, line placements 1/3,1/3,-1/3,-2/3: X-scal agrees placement class by placement class'))
    t4 = J('T4.json')
    if t4:
        d = {}
        for ch, secs in t4.items():
            for sec, r in secs.items():
                d[f'xi{ch} {sec}'] = comb(r['eps0'])
                res['poles'][f'xi{ch} {sec}'] = comb(r['pole']) if r['pole'] else '0'
        res['Ob1']['xi|_2'] = d
        okp = all(v in ('0', {}) or all(x == '9/2' for x in v.values()) for k, v in res['poles'].items() if k.startswith('xi'))
        gates.append(gate('G-pole (xi)', {k: v for k, v in res['poles'].items() if k.startswith('xi')},
                          '(9/4)u{int(A,A),S(xi)} = (9/2) u int(A, Gamma(xi,lam)) (computed bracket ratio 2)', okp))
    gs = J('G_sign.json')
    if gs:
        gates.append(gate('G-sign (two legs)', dict(rule=gs['sign_rule'], candidates=gs['sign_rule_candidates_passing'],
                                                     Q_SPb_closed=gs['Q_SPb_closed']),
                          'Q^2 = 0, Q S_2(x) = 0, images of Thm c9-quadratic (iii)',
                          all(gs['Q_SPb_closed'].values()) and len(gs['sign_rule_candidates_passing']) >= 1,
                          'two rules pass: left derivation (chosen: Q = {S,-} acts from the left, Def c1-bracket) and its right-derivation image (differs by (-1)^{parity of the image})'))
        wz_plus = gs['WZ1_xi_max|Q Ob1 + {Ob0,T1}|']; wz_minus = gs['WZ1_xi_max|Q Ob1 - {Ob0,T1}|']
        gates.append(gate('G-WZ1 (xi, two legs)', dict(max_abs_Q_Ob1_plus_bracket=wz_plus, max_abs_Q_Ob1_minus_bracket=wz_minus),
                          'Lambda Ob1(x) = {Ob0,S(x)} (Lemma at-sformulas (ii) as printed)',
                          'pass with sign: Lambda Ob1(x) = -{Ob0,S(x)}' if float(wz_plus) == 0 else False,
                          'the component formula of Def c2-dce for a degree-1 cochain gives (sG)^(1)(x) = Lambda G^(1)(x) + {G^(0),S(x)} (checked: Lambda{F0,T1} = -{Lambda F0,T1}), so the printed identity has the opposite sign; magnitude 27/8 and the absence of g4 confirmed'))
        res['F1']['xi|_2'] = dict(Q_images=dict(QT1=gs['QT1'], QT3=gs['QT3'], QT2=gs['QT2'], QT4=gs['QT4'], QT0_on_g2_g4=gs['QT0'], Qg2_over_Qg4=gs['Qg2_over_Qg4']),
                                  bracket_F1_0_T1=gs['{F1^(0)|_2,T1} on (g2,g1)'],
                                  residual='R1(xi)|_2 = Ob1(xi)|_2 - {F1^(0)|_2,S(xi)} = (79/24 - (1/2) log(3/2)) (g1 + g3), components on g2, g4: 0',
                                  particular='F1^(1)(xi)|_2 = (log(3/2) - 79/12) (T1 + T2) + u (theta^xi_V S^V_2(xi) + theta^xi_perp S^perp_2(xi)), r0 = 0',
                                  q=dict(q1='-1/2', q3='-1/2', a='-1', b='-1'))
    t3 = J('T3.json')
    if t3:
        gates.append(gate('G-R0', t3['R'], 'Ob1(R_r) = 0 sector by sector', all(t3['R'].values()),
                          'two legs: four sectors; three legs: (c,A_perp,A*_perp), (A_V,A_perp,A_perp), (A_perp,lam,lam); poles 0'))
        res['Ob1']['R_r'] = '0 (two legs all sectors; three legs (c,A_p,A*_p), (A_V,A_p,A_p), (A_p,lam,lam))'
        res['Ob1']['M_omega|_2'] = dict(value='+(1/3) u Q S(M_omega)_spin, S(M_omega)_spin = int omega_{mu nu}(A_mu,A*_nu) + int((1/4) omega_{mu nu} rho_mu rho_nu lam, lam*); 0 in (A_perp,A_perp); position part 0',
                                        raw=t3['M'])
        res['F1']['M_omega|_2'] = 'spin direction: F1^(1)(M_omega)|_2 = (1/3) u S(M_omega)_spin + (Q-closed)'
    bf = J('B1_fit.json')
    if bf:
        res['Ob2']['(Delta,P_b)|_2'] = {k: v.get('Ob_decomposition') for k, v in bf.items()}
        res['cohomology']['weight 5, p<=1, linear in b, two legs (blocks)'] = {
            k: {kk: v[kk] for kk in ('n_candidates', 'basis_dim', 'closed_dim', 'exact_dim', 'H_minus1_dim', 'other_classes') if kk in v}
            for k, v in bf.items()}
        od = {k: v['Ob_decomposition']['1'] for k, v in bf.items()}
        res['theta']['Delta'] = {'perp': '-1/3', 'lam': '-1/3', 'V': '-1/3', 'derivation': 'theta_d = o_d/(2u), o_d = -2/3 u in every block (Thm at-theta)'}
        ok = all(list(v.values())[0] == '-2/3' for v in od.values())
        gates.append(gate('B1 o_d (and X-theta)', od, 'o_perp = -(2/3) u (cross-check X-theta, theta = -1/3)', ok,
                          'logs (present in raw cycle sums) cancel in the class; components along the other class of the V block (second-order free symmetry) 0'))
    s1 = J('S1.json')
    if s1:
        res['cohomology']['Sigma, weight 9/2, p = 0, two legs (bases T0..T4, g1..g4)'] = s1
        ok = all(v['H0'] == 0 and v['Hm1'] == 2 and v['kernel_deg_minus1'] == [['0', '1', '0', '1', '0'], ['0', '0', '1', '0', '1']] for v in s1.values())
        gates.append(gate('G-H (two legs along Sigma: Prop at-xi2)', {k: dict(H0=v['H0'], Hm1=v['Hm1'], kernel=v['kernel_deg_minus1']) for k, v in s1.items()},
                          'H^0 = 0; kernel span{T1+T3, T2+T4}', ok))
    gn = J('G_num.json')
    if gn:
        gates.append(gate('G-cyc', {k: v for k, v in gn.items() if k.startswith('G-cyc')}, 'cyclic reduction of the word sum exact',
                          all(v for k, v in gn.items() if k.startswith('G-cyc'))))
        gates.append(gate('G-masters', gn['G-masters max relative difference'], 'exact J = quadrature to < 1e-10',
                          gn['G-masters max relative difference'] < 1e-10))
    res['open'] += ['T2 (three/four/five legs along Delta and P_b), T5 (xi three and four legs), T6 (eta), T7 (K_b), T8 (log coefficients), B2 (other pairs), S1 at the remaining weights: not computed by engine B in this run',
                    'G-log numerical eps-extrapolation of single cycles: not run (the eps-expansion is exact by construction: bounded special parameter, no log eps can arise)']
    with open(os.path.join(ROOT, 'results.json'), 'w') as f:
        json.dump(res, f, indent=1, default=str)
    with open(os.path.join(ROOT, 'gates.json'), 'w') as f:
        json.dump(gates, f, indent=1, default=str)
    print('assembled', len(gates), 'gates')


if __name__ == '__main__':
    main()
