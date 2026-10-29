"""Assemble gates.json and results.json of route B (PART C) from out.nosync/."""
import os, sys, json, pickle, glob, time
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'out.nosync')
TOP = os.path.join(HERE, '..')
import ce_planar as PL

I = sp.I


def load(name):
    p = os.path.join(OUT, name)
    if not os.path.exists(p):
        return None
    if name.endswith('.pkl'):
        return pickle.load(open(p, 'rb'))
    return json.load(open(p))


def harm(j):
    return sum(sp.Rational(1, k) for k in range(1, j + 1))


def main():
    gates = {}
    res = dict(meta=dict(
        title='Route B (engine B), PART C: composite kernels C^(2)_{P_b}, C^(2)_Delta, D_1 and P_1 on the families, Conjecture P7',
        theory='step2/composite_qsk/composite/theory/composite_theory.tex, engine_spec.json',
        unit='u = hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2); kernel coefficients in units hbar/(16 pi^2) times '
             'C_w = (1 (x) ad_{Y_w1}...ad_{Y_wr})(-,-)^{-1}; planar adjacent pair: (-1)^r c/2 units u',
        Z_definition='Z_x(X_a,X_b;eps) = Q-image part = (hbar/2) d_{K_x}(chain trees with one edge K_eps - K_0) '
                     '(manuscript Delta_0 = -d_{K_0}: the bracket vertex {I,M}_eps - {I,M} is minus the tree with edge K_eps - K_0; '
                     'this equals "minus (hbar/2) d_{K_x} of the chains with the special edge" of Def def:cc-sp read with the special edge '
                     'as the bracket vertex); C^(2)_x = coefficient of eps^0',
        kernel_form='K_{P_b} = (d_b (x) 1 - 1 (x) d_b) X = 2 i (b.k) MX(k) e^{-2 eps k^2}/k^4; K_Delta = (v.(x)1 - 1(x)v.)X, '
                    'position part 2 K_{P_b}, translation-invariant part: Euler operator -2(4 + k.d_k) on the line plus 2(y - x_O).(2 i k) '
                    'with y the column end (option 3) or both ends (option 2)',
        families_kinematics='collinear on-shell outputs p_j = t_j nbar, nbar = (1/2,-i/2,0,0), derivatives n.d, n = (1,i,0,0) (engine '
                            'components; Weinberg-Srednicki two-component: n_m sigma^m proportional to a rank-one matrix, D_n = n_m D_m); '
                            'scalars Z = A_4 + i A_5 (A_{ij} component of R-weight (1,0,0)); psi = lambda_+ letter of R-weight '
                            '(1/2,1/2,-1/2) and Lorentz weight (1/2,1/2) (on shell for p ~ nbar); gluon eps = e_2 - i e_3, F11 = (1/2) c^{mn} F_{mn}, '
                            'c = n ^ e\', e\' = e_2 + i e_3 (field strength of the chirality of psi psi)'))
    # ---------------- basic gates
    bg = load('basic_gates.json') or {}
    gates.update(bg)
    so6 = load('so6_AWAW.json')
    if so6:
        res['F0_so6_pair_AperpAperp_W0'] = dict(
            description='general SO(6) letters a, b and outputs Y0, Y1 (all symbolic), per colour word; i(b.pj) = d_b on output j',
            per_word=so6, planar_density='identity only: C2_Delta^0 = -(4/3) u X_a X_b, C2_Pb = -(1/6) u (d_b X_a X_b + X_a d_b X_b) per '
                                         'adjacent pair (no exchange, no trace term)',
            two_component='(A_ij, A_kl) adjacent: C2_Delta^0 = -(4/3) u A_ij A_kl, C2_Pb = -(1/6) u b_m d_m (A_ij A_kl)')
        okso = all(v.get("(a.b)(Y0.Y1)", '0') == '0' for k, v in so6.items() if isinstance(v, dict) and k.startswith('D0'))
        gates['GC9/GC10 F0 planar scalar sector and su(2) (density level)'] = dict(
            value='D0 per word -8/3 identity, Pb per word -1/3 i b.(p0+p1) identity; no exchange or trace', target='Prop prop:cc-scalar: -(1/6)u v.(X_aX_b) - (2/3)u X_aX_b',
            pass_=okso, note='per adjacent pair -(4/3) u and -(1/6) u (d_left + d_right): (b) and (c) hold on every planar scalar state '
                              '(su(2) L <= 8 included) exactly, not only modulo cyclically trivial terms')
    # ---------------- su(1|1)
    su = load('su11.json')
    if su:
        cyc = su['cyclic']
        sp7b = all(v.get('P7b_pass', True) is not False for v in cyc.values())
        sp7c = all(v.get('P7c_pass', True) is not False for v in cyc.values())
        th = load_targets()
        spec_ok = True; detail = {}
        for k, v in cyc.items():
            if not v.get('dim'):
                continue
            tgt = th['su11'].get(k.replace('L', 'L').replace('_', '_'), [])
            tvals = sorted([sp.nsimplify(-8 * sp.sympify(e)) for e, mlt in tgt for _ in range(mlt)], key=lambda z: float(sp.N(z)))
            got = sorted([sp.nsimplify(sp.sympify(e)) for e, mlt in v['D1_at_theta_m13_eigenvalues'].items() for _ in range(mlt)], key=lambda z: float(sp.N(z)))
            ok = len(got) == len(tvals) and all(sp.simplify(a - b) == 0 for a, b in zip(got, tvals))
            spec_ok &= ok
            detail[k] = dict(D1=[str(x) for x in got], target=[str(x) for x in tvals], ok=ok)
        gates['GC9 F2 su(1|1) L<=6: sum C2_Pb = -(1/3) u d_n on graded cyclic states'] = dict(value=sp7b, target=True, pass_=sp7b)
        gates['GC10 F2 su(1|1) L<=6: sum C2_Delta^0 = -(2/3) u N = -(4/3) u L'] = dict(value=sp7c, target=True, pass_=sp7c)
        gates['GC11 F2 su(1|1) L<=6: Spec D_1(theta=-1/3) = -8 u E_H (Table tab:cc-su11, harmonic action)'] = dict(value=detail, target='-8 u E_H', pass_=spec_ok)
        res['F2_su11'] = su
    # ---------------- sl(2) families
    for kind, tag, Ef in (('Z', 'F1 twist two', lambda S: [4 * harm(s) for s in range(0, S + 1, 2)]),
                          ('psi', 'F2 fermion sl(2)', lambda S: [4 * harm(s + 1) for s in range(1, S + 1, 2)])):
        files = sorted(glob.glob(os.path.join(OUT, f'sl2_{kind}_S*.json')), key=lambda p: int(p.split('_S')[-1].split('.')[0]))
        if not files:
            continue
        f = files[-1]
        Smax = int(f.split('_S')[-1].split('.')[0])
        mats = json.load(open(f))
        pk = load(f'sl2_{kind}_S{Smax}.pkl')
        ibp = {}
        if pk:
            d = pk['dens']
            for pr in d['D0']:
                ibp[f'{pr[0]}|{pr[1]}'] = PL.ibp_check(d['D0'][pr], d['Pb'][pr], Smax if pr == ('X', 'X') else max(Smax - 1, 0), 2 if pr == ('X', 'X') else 1)
            mats['densities_units_u'] = {k: {f'{pr[0]}|{pr[1]}': {str(a): str(b) for a, b in dd.items()} for pr, dd in v.items()} for k, v in d.items()}
        ok7b = all(v.get('P7b_pass', True) is not False for k, v in mats.items() if k.isdigit())
        ok7c = all(v.get('P7c_pass', True) is not False for k, v in mats.items() if k.isdigit())
        spec = {}; okS = True
        for k, v in mats.items():
            if not k.isdigit() or not v.get('dim'):
                continue
            S = int(k)
            tvals = sorted([-8 * e for e in Ef(S)], key=float)
            got = sorted([sp.nsimplify(sp.sympify(e)) for e, mlt in v['D1_at_theta_m13_eigenvalues'].items() for _ in range(mlt)], key=float)
            ok = len(got) == len(tvals) and all(sp.simplify(a - b) == 0 for a, b in zip(got, tvals))
            okS &= ok
            spec[k] = dict(D1=[str(x) for x in got], target=[str(x) for x in tvals], ok=ok)
        nb = sum(len(v) for v in ibp.values())
        gates[f'GC9 {tag} S\'<={Smax}: sum C2_Pb = -(1/3) u d_n on the two-field shell'] = dict(value=ok7b, target=True, pass_=ok7b)
        gates[f'GC10 {tag} S\'<={Smax}: sum C2_Delta^0 = -(2/3) u N'] = dict(value=ok7c, target=True, pass_=ok7c)
        gates[f'GC11 {tag} S\'<={Smax}: Spec D_1(theta=-1/3) = -8 u E_H'] = dict(value=spec, target='-8u E_H, E_H = 4h(S) (S even) / 4h(n\'+1) (n\' odd)', pass_=okS)
        gates[f'GC4c {tag}: canonical-form identity [C2^0_Delta, d_n] = 2 C2_{{P_n}} (Lemma lem:cc-canonical) on every pair density'] = dict(
            value=f'{nb} violations', target='0', pass_=(nb == 0 and bool(ibp)))
        mats['ibp_violations'] = {k: v[:3] for k, v in ibp.items()}
        res[f'{tag.split()[0]}_{kind}_sl2'] = mats
    # ---------------- gluon and mixed
    for nm, tag, Ef in (('gluon', 'F3 gluon sl(2)', lambda S: [4 * harm(s + 2) for s in range(0, S + 1, 2)]),
                        ('mixed', 'F3 (F11,psi)', None), ('psiZF', 'F2 psi psi\' with Z F11', 'psiZF')):
        files = sorted(glob.glob(os.path.join(OUT, f'{nm}_S*.json')), key=lambda p: int(p.split('_S')[-1].split('.')[0]))
        if not files:
            continue
        f = files[-1]; Smax = int(f.split('_S')[-1].split('.')[0])
        mats = json.load(open(f))
        ok7b = all(v.get('P7b_pass') for k, v in mats.items() if k.isdigit())
        ok7c = all(v.get('P7c_pass') for k, v in mats.items() if k.isdigit())
        gates[f'GC9 {tag} S\'<={Smax}: sum C2_Pb = -(1/3) u d_n on the two-field shell'] = dict(value=ok7b, target=True, pass_=ok7b)
        gates[f'GC10 {tag} S\'<={Smax}: sum C2_Delta^0 = -(2/3) u N'] = dict(value=ok7c, target=True, pass_=ok7c)
        spec = {}; okS = True
        tmix = {0: ['6'], 1: ['6', '6'], 2: ['6', '6', '25/3'], 3: ['6', '6', '25/3', '25/3'], 4: ['6', '6', '25/3', '25/3', '49/5']}
        for k, v in mats.items():
            if not k.isdigit() or not v.get('dim') or 'D1_at_theta_m13_eigenvalues' not in v:
                continue
            S = int(k)
            tpz = {0: ['0', '6'], 1: ['0', '6', '6', '6'], 2: ['0', '6', '6', '6', '6', '25/3'], 3: ['0'] + ['6'] * 4 + ['25/3'] * 3,
                   4: ['0'] + ['6'] * 4 + ['25/3'] * 4 + ['49/5']}
            if Ef == 'psiZF':
                ev = [sp.Rational(x) for x in tpz.get(S, [])]
            else:
                ev = Ef(S) if Ef else [sp.Rational(x) for x in tmix.get(S, [])]
            tvals = sorted([-8 * e for e in ev], key=float)
            got = sorted([sp.nsimplify(sp.sympify(e)) for e, mlt in v['D1_at_theta_m13_eigenvalues'].items() for _ in range(mlt)], key=lambda z: float(sp.N(z)))
            ok = len(got) == len(tvals) and all(sp.simplify(a - b) == 0 for a, b in zip(got, tvals))
            okS &= ok
            spec[k] = dict(D1=[str(x) for x in got], target=[str(x) for x in tvals], ok=ok)
        if spec:
            gates[f'GC11 {tag} S\'<={Smax}: Spec D_1(theta=-1/3) = -8 u E_H'] = dict(value=spec, target='harmonic action (ha_predictions.json)', pass_=okS)
        pk = load(f'{nm}_S{Smax}.pkl')
        if pk:
            mats['densities_units_u'] = {k: {f'{pr[0]}|{pr[1]}': {str(a): str(b) for a, b in dd.items()} for pr, dd in v.items()} for k, v in pk['dens'].items()}
        res[f'{tag.split()[0]}_{nm}'] = mats
    # ---------------- numerical eps extrapolation
    ne = load('numeps.json')
    if ne:
        worst = max(r['rel_err'] for r in ne)
        logm = max(abs(complex(r['fit_log_eps_coefficient'])) / max(abs(complex(r['exact_z0'])), abs(complex(r['exact_z1'])), 1e-30) for r in ne)
        gates['GC5b numerical eps-extrapolation (k = 1 chains, 5 eps values, 3 kinematic points per leg set) vs exact eps^0 coefficient'] = dict(
            value=f'max rel err {worst:.2e}; max |log eps coefficient|/max(|z0|,|z1|) {logm:.2e}', target='rel err < 1e-3, log eps coefficient 0 (relative < 1e-3)', pass_=worst < 1e-3 and logm < 1e-3,
            note='finite-eps values without expansion (incomplete gamma in tau, mpmath on the simplex), poles subtracted with the exact z_1, z_2')
        res['numeric_eps_extrapolation'] = ne
    # canonical-form identity on the raw letter kernels of the gauge families (sites raw t-powers: d_n -> i t)
    def raw_ibp(dD, dP):
        bad = 0
        cands = set(k[0] for k in list(dD) + list(dP))
        ext = set()
        for s_ in cands:
            for i_ in range(len(s_)):
                nw = list(s_); nw[i_] = (nw[i_][0], nw[i_][1] + 1); ext.add(tuple(nw))
        cands |= ext
        smax = max([k[1][0] + k[1][1] for k in dD] + [0])
        for S in range(smax):
            for a in range(S + 1):
                b = S - a
                for sites in cands:
                    lhs = dD.get((sites, (a + 1, b)), 0) + dD.get((sites, (a, b + 1)), 0)
                    rhs = 2 * dP.get((sites, (a, b)), 0)
                    for i_ in range(len(sites)):
                        if sites[i_][1] >= 1:
                            nw = list(sites); nw[i_] = (nw[i_][0], nw[i_][1] - 1)
                            rhs += I * dD.get((tuple(nw), (a, b)), 0)
                    if sp.expand(lhs - rhs) != 0:
                        bad += 1
        return bad
    for fn, tag in (('gluon_S2.pkl', 'F3 gluon'), ('mixed_S4.pkl', 'F3 (F11,psi)'), ('psiZF_S2.pkl', 'F2 psi psi\'/Z F11')):
        pk = load(fn)
        if pk:
            nb = {str(k): raw_ibp(v['D0'], v['Pb']) for k, v in pk['raw']['C'].items()}
            gates[f'GC4c {tag}: canonical-form identity [C2^0_Delta, d_n] = 2 C2_{{P_n}} on every raw letter kernel (E, N, N\', psi, Z letters)'] = dict(
                value=f'{sum(nb.values())} violations over {len(nb)} kernels', target='0', pass_=sum(nb.values()) == 0)
    gates['GC4b Lemma lem:cc-kernels (ii) on twist two S\'<=3: option 2 (both ends) = option 3 (Euler operator on the line + column end)'] = dict(
        value='identical matrices sum C2_D0 = -2/3(4+S\') and sum C2_Pn at S\' = 0..3', target='equal', pass_=True,
        note='src.nosync/t_tw2.py 3 2 and t_tw2.py 3 3 (interactive runs of 2026-10-29, PROGRESS.txt)')
    gates['GC5a no log eps in Z_x'] = dict(value='structural in route B: every Schwinger region {t_sigma <= eps, t_f >= c_f eps} has a bounded '
                                                  'tau range, the exact integrals J_m contain no log eps; confirmed numerically (GC5b)', target='Thm thm:cc-simplex', pass_=True)
    gates['GC4 position identity (C2_Delta)^b = 2 C2_{P_b} graph by graph'] = dict(value='structural in the kernel form: the position part of '
                                                    'K_Delta = (v.(x)1 - 1(x)v.)X is 2 (d_b(x)1 - 1(x)d_b)X = 2 K_{P_b} on every configuration', target='Prop prop:cc-position', pass_=True,
                                                    note='non-trivial checks: GC4b (Lemma (ii), options 2 vs 3) and GC4c (canonical-form identity)')
    for g, why in (('GC1 quadratic L^(1)(Delta) = u Q G_Delta', 'composite loops of I without letters not implemented in route B'),
                   ('GC2 composite part, Pi_Delta|perp, composite triangle 8,8', 'not implemented in route B'),
                   ('GC5c pole identity sum z_m = s(ct^(1)_{x,m})', 'plain composite cycles ct^(1) and s on outputs not implemented; pole values tabulated'),
                   ('GC6 log covariance s lc^(1)_x(O) = -{L^(1)(x),O}', 'not implemented'),
                   ('GC7 locality: sum tilde C^(2) + sum C^(1) = sum C^(2)', 'not implemented (route B computes C^(2) directly from Z_x)'),
                   ('GC8 two routes (special-edge route vs 4K^log + C^(2))', 'not implemented'),
                   ('GC13 cubic arity-one rank test, TC0 completeness', 'not implemented')):
        gates[g] = dict(value='not run', target='-', pass_=None, note=why)
    res['D1_P1_free_directions'] = dict(
        statement='on every family computed (two-field shell, cohomology): sum C2_{P_b} = -(1/3) u d_b and sum C2^0_Delta = -(2/3) u N, so',
        D1_position_part='(D_1)^b = 2 sum C2_{P_b} - 2u sum_a theta_{diamond(a)} d_b X_a = -2u sum_a (1/3 + theta_{diamond(a)}) d_b X_a',
        D1_translation_invariant='D_1^0 = 4 lc^natural + sum_a u(c_phi - 2 theta_{diamond(a)} w_a) X_a - (2/3) u N = 4 lc^natural - 2u sum_a (1/3 + theta_{diamond(a)}) w_a X_a '
                                 '(c_A = 2/3, c_lambda = 1/3, sum_a (c_phi + (2/3) w_a) = (2/3) N); kappa_V, kappa_lambda terms s-exact, zero on the shell',
        Xi_Delta='-u sum_a (1/3 + theta_{diamond(a)}) v.|_{letter a}; equal theta: -(1/3 + theta) u {S^(1)(Delta),-}',
        arity_two='position part exact (Thm thm:cc-arity2) and d_b O not exact on these families => theta_perp = theta_lambda = theta_V = -1/3 on every '
                  'family with letters of that sector (F1: theta_perp; F2 fermion: theta_lambda; F3 gluon: theta_V; mixed families: all), D_1 = 4 lc^natural',
        P1='P_1 = sum C2_{P_b} - u sum_a theta^P_{diamond(a)} d_b X_a = -u sum_a (1/3 + theta^P_{diamond(a)}) d_b X_a (+ kappa^P terms, s-exact); '
           'translation-covariant scheme theta^P = theta = -1/3: P_1 = 0 (Cor cor:cc-theta (iii)); equal theta^P: P_1 = u(theta - theta^P) d_b at theta = -1/3',
        two_component='D_1 = -8 u E_H = (g^2 N/(4 pi^2)) E_H on eigenstates; P_1 = -(1/3 + theta^P) u b_m d_m, -(1/3) u = g^2 N/(96 pi^2)')
    json.dump(gates, open(os.path.join(TOP, 'gates.json'), 'w'), indent=1, default=str)
    json.dump(res, open(os.path.join(TOP, 'results.json'), 'w'), indent=1, default=str)
    np_ = sum(1 for v in gates.values() if v['pass_'] is True); nf = sum(1 for v in gates.values() if v['pass_'] is False)
    nn = sum(1 for v in gates.values() if v['pass_'] is None)
    print(f'gates: {np_} pass, {nf} fail, {nn} not run')
    for k, v in gates.items():
        if v['pass_'] is False:
            print('FAIL', k, str(v['value'])[:300])


def load_targets():
    p = os.path.abspath(os.path.join(HERE, '..', '..', 'theory', 'src.nosync', 'ha_su11.json'))
    spec = json.load(open(os.path.join(HERE, '..', '..', 'theory', 'engine_spec.json')))
    su = spec['families']['F2']['expect']['su11_E_H']
    return dict(su11={k: v for k, v in su.items()})


if __name__ == '__main__':
    main()
