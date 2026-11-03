"""Tasks T3, T6, B2 (partial), T2 (partial), T5 (tables), T8 (partial) of engine_spec.json, route A."""
from fractions import Fraction
import numpy as np
import sa_core as SA
import sa_ob as SO
import sa_fit as SF
import sa_mode as MD
from sa_tasks import rngs, gate, two_leg_fit, zero_check, B4, V0, _same
from r2_core import Q2

_MM = {}


def mmodel(p):
    k = tuple(int(round(x.real)) for x in p)
    if k not in _MM:
        _MM[k] = MD.MM(np.array(k))
    return _MM[k]


def mbasis(builder):
    """basis callable from a mode-model functional builder(mm) -> GA (two legs, first leg at mode 0)"""
    cache = {}
    def f(L):
        mm = mmodel(L[0]['p0'])
        key = id(mm)
        if key not in cache:
            cache[key] = builder(mm)
        v = mm.value(cache[key], [(SA.E.SLOT[l['t']], l['pol']) for l in L])
        return (v.re, v.im)
    return f


def Tbasis(xi):
    return {'T0': mbasis(lambda mm: MD.F_T0(mm, xi)), 'T1': mbasis(lambda mm: MD.F_T1(mm, xi)),
            'T2': mbasis(lambda mm: MD.F_T1(mm, xi, range(4, 10))), 'T3': mbasis(lambda mm: MD.F_T3(mm, xi)),
            'T4': mbasis(lambda mm: MD.F_T3(mm, xi, range(4, 10)))}


def fit_T(spec, xi, chir, other=None, sgn=1, seed=41, what='fin'):
    """fit a degree -1, weight 9/2 two-leg density onto T0..T4 along xi of chirality chir"""
    o = 'M' if chir == 'P' else 'P'
    Tb = Tbasis(xi)
    secs = [('AV', 'lams' + chir, ['T0', 'T3']), ('AW', 'lams' + o, ['T4']), ('lam' + chir, 'AsW', ['T2']),
            ('lam' + o, 'AsV', ['T1'])]
    rng, rp, rpol = rngs(seed)
    out = {}
    for t1, t2, names in secs:
        S = []
        for _ in range(len(names) + 2):
            p = rp(); L = [SA.leg(t1, rpol(t1), p), SA.leg(t2, rpol(t2), -p)]
            if what == 'fin':
                d = SO.density(spec, L, other=other)
                a = d['fin']['alpha'].get('(Y0,Y1)', V0())
            else:
                a = plain_density(spec, L)
            S.append(([Tb[n](L) for n in names], SA.vscale(a, sgn)))
        c, _ = SF.fit(S, len(names))
        for n, x in zip(names, c):
            out[n] = SF.vstr(x)
    return out


def plain_density(spec, L):
    """log eps coefficient of the plain arity-one cycles (canonical alpha; spec without position part)"""
    g = SO.gens_of(spec)
    assert len(g) == 1
    r = SA.cycle_sum(L, [g['0']], slotnames=['x'], what=('log',))
    return SO.reduce_colour(r, 2, 'log', 0).get('(Y0,Y1)', V0())


# ================================================================================================ T3
def task_T3():
    res, gates = {}, []
    Rm = np.zeros((6, 6)); Rm[0, 1], Rm[1, 0] = 1, -1; Rm[2, 4], Rm[4, 2] = 2, -2; Rm[3, 5], Rm[5, 3] = -1, 1
    specR = {'kind': 'R', 'R': Rm}
    z2 = zero_check(specR, [('AV', 'AV'), ('AW', 'AW'), ('lamP', 'lamM'), ('c', 'AsV'), ('AW', 'AV'), ('lamP', 'lamP'),
                            ('lamM', 'lamM')], seed=51)
    rng, rp, rpol = rngs(53)
    z3 = {}
    for types in (('c', 'AW', 'AsW'), ('AW', 'lamP', 'lamP'), ('AV', 'AW', 'AW')):
        p0, p1 = rp(), rp(); moms = [p0, p1, -p0 - p1]
        L = [SA.leg(t, rpol(t), m) for t, m in zip(types, moms)]
        d = SO.density(specR, L)
        z3[str(types)] = SF.vstr(d['fin']['alpha'].get('(Y0,[Y1,Y2])', V0()))
    res['Ob1_R'] = {'two legs': z2, 'three legs (alpha, colour (Y0,[Y1,Y2]))': z3, 'r': Rm.tolist()}
    okR = all(v.startswith('alpha=0,') for v in z2.values()) and all(v == '0' for v in z3.values())
    gates.append(gate('G-R0', {'two legs': z2, 'three legs': z3}, 'Ob^(1)(R_r) = 0', okR, 'cycle sums per sector'))
    # M_omega: v = Om x
    Om = np.zeros((4, 4)); Om[0, 1], Om[1, 0] = 1, -1; Om[2, 3], Om[3, 2] = 2, -2; Om[0, 2], Om[2, 0] = -1, 1
    specM = {'kind': 'M', 'Om': Om}
    spin = mbasis(lambda mm: mm.Q(mm.quad('A', 'As', lambda i, j, p1, p2: Om[i, j] if (i < 4 and j < 4) else 0)))
    out = {}
    S = []
    for t1, t2 in (('AV', 'AV'), ('c', 'AsV')):
        for _ in range(3):
            p = rp(); L = [SA.leg(t1, rpol(t1), p), SA.leg(t2, rpol(t2), -p)]
            d = SO.density(specM, L)
            S.append(([spin(L)], d['fin']['alpha'].get('(Y0,Y1)', V0())))
            assert all(SA.viszero(d['fin']['beta'][m].get('(Y0,Y1)', V0())) for m in range(4))
    c, _ = SF.fit(S, 1)
    out['V sector (A_V,A_V)_2 + (c,A*_V)_1'] = {'basis': 'Q int Om_{mu nu}(A_mu, A*_nu)', 'coef': SF.vstr(c[0])}
    out['other sectors'] = zero_check(specM, [('AW', 'AW'), ('AV', 'AW')], seed=57)
    Rs = sum(Om[m, n] * (SF.SS[m] @ SF.SB[n]) for m in range(4) for n in range(4))
    sl = lambda p: sum(p[m] * SF.SB[m] for m in range(4))
    def bR(L):     # int (rho(Om) lam, dslash lam), rho(Om) = sum Om_{mu nu} s[mu] sb[nu] on Sigma
        u, p, w, q = L[0]['pol'], L[0]['p0'], L[1]['pol'], L[1]['p0']
        return 1j * ((Rs @ u) @ (sl(q) @ w)) - 1j * ((Rs @ w) @ (sl(p) @ u))
    S = []
    for _ in range(3):
        p = rp(); L = [SA.leg('lamP', rpol('lamP'), p), SA.leg('lamM', rpol('lamM'), -p)]
        d = SO.density(specM, L)
        S.append(([bR(L), SF.b_ldl(L)], d['fin']['alpha'].get('(Y0,Y1)', V0())))
    c2, _ = SF.fit(S, 2)
    out['lambda sector (lam_+,lam_-)_1'] = {'basis': ['int(rho(Om) lam, dslash lam), rho(Om) = sum Om_{mu nu} s[mu] sb[nu]',
                                                      'int(lam, dslash lam)'], 'coef': [SF.vstr(x) for x in c2]}
    out['position part (beta) on two legs'] = 'zero (= Ob(P_{Om b})|_2 = 0)'
    res['Ob1_M_2'] = out
    res['F1_M_2_spin'] = {'F1(M_Om)|2 contains': '%s u int Om_{mu nu}(A_mu, A*_nu) + sector parts' % SF.vstr(c[0]),
                          'rule': 'Q F1 = Ob1 on (A_V,A_V)_2 + (c,A*_V)_1 (Prop at-two (iii)); v = Om x, v.A = L_v A'}
    okM = all(v.startswith('alpha=0,') for v in out['other sectors'].values())
    gates.append(gate('M_omega two legs', out, 'Ob1(M)|2 in the V sector = k Q int Om(A,A*); 0 elsewhere', okM, ''))
    return res, gates


# ================================================================================================ T6 (eta, intrinsic, two legs)
def task_T6():
    res, gates = {}, []
    for chir in ('P', 'M'):
        o = 'M' if chir == 'P' else 'P'
        eta = SO.eta_generic(chir); spec = {'kind': 'eta', 'eta': eta}
        ss = lambda v: sum(v[m] * SF.SS[m] for m in range(4))
        def bAd(L, eta=eta):        # int (eta, lam) d^*A_V : (-i q.a)(u.eta)
            return -1j * (L[0]['p0'] @ L[0]['pol'][:4]) * (L[1]['pol'] @ eta)
        def bArho(L, eta=eta):      # int (eta, rho(dA_V) lam) = sum_{mu != nu} i q_mu a_nu (s[mu] sb[nu] u).eta
            q, a, u = L[0]['p0'], L[0]['pol'], L[1]['pol']
            return sum(1j * q[m] * a[n] * ((SF.SS[m] @ (SF.SB[n] @ u)) @ eta) for m in range(4) for n in range(4) if m != n)
        def bWrho(L, eta=eta):      # int (eta, rho(dA_perp) lam) = sum_{mu, I} i q_mu a_I (s[mu] sb[I] u).eta
            q, a, u = L[0]['p0'], L[0]['pol'], L[1]['pol']
            return sum(1j * q[m] * a[I] * ((SF.SS[m] @ (SF.SB[I] @ u)) @ eta) for m in range(4) for I in range(4, 10))
        def bc(L, eta=eta):         # int (rho(dc) eta, lam*) = sum i p_mu c (s[mu] eta).s
            p, c, s = L[0]['p0'], L[0]['pol'][0], L[1]['pol']
            return 1j * c * ((ss(p) @ eta) @ s)
        secs = [('AV', 'lam' + chir, ['int(eta,lam) d*A_V', 'int(eta,rho(dA_V)lam)'], [bAd, bArho]),
                ('AW', 'lam' + o, ['int(eta,rho(dA_perp)lam)'], [bWrho]),
                ('c', 'lams' + o, ['int(rho(dc)eta,lam*)'], [bc])]
        out = {}
        for t1, t2, names, bfs in secs:
            rng, rp, rpol = rngs(61)
            S, betanz = [], False
            for _ in range(len(bfs) + 2):
                p = rp(); L = [SA.leg(t1, rpol(t1), p), SA.leg(t2, rpol(t2), -p)]
                d = SO.density(spec, L)
                S.append(([f(L) for f in bfs], d['fin']['alpha'].get('(Y0,Y1)', V0())))
                betanz = betanz or any(not SA.viszero(d['fin']['beta'][m].get('(Y0,Y1)', V0())) for m in range(4))
            try:
                c, _ = SF.fit(S, len(bfs))
                out['(%s,%s)' % (t1, t2)] = {'basis': names, 'coef': [SF.vstr(x) for x in c],
                                             'beta (position part along b = Ob1(rho(b)eta))': 'nonzero' if betanz else 'zero'}
            except ValueError as e:
                out['(%s,%s)' % (t1, t2)] = {'basis': names, 'fit': 'inconsistent: %s' % e,
                                             'table': [[str(SF.q2(v)) for v in s[0]] + [SF.vstr(s[1])] for s in S]}
        res['Ob1_eta%s_2_intrinsic' % ('+' if chir == 'P' else '-')] = out
    res['G-pos'] = ('route A: the position part of every cycle along Delta, M_om, eta is by construction the cycle along the '
                    'derivative of the Killing data (beta^mu = -i G_mu(0), G_mu = cycle with S(d_mu v)); structural, no '
                    'separate number')
    return res, gates


# ================================================================================================ B2 (pairs, partial)
def task_B2():
    res, gates = {}, []
    for chir in ('P', 'M'):
        xi = SO.xi_generic(chir)
        out = fit_T({'kind': 'Delta'}, xi, chir, other={'kind': 'xi', 'xi': xi}, sgn=1)
        closed = _same(out['T1'], out['T3']) and _same(out['T2'], out['T4'])
        res['Ob2_Delta_xi%s_2' % ('+' if chir == 'P' else '-')] = {
            'on T0..T4 (p=0 intrinsic)': out, 'Q-closed (T0 = 0, T1 = T3, T2 = T4)': closed and out['T0'] == '0',
            'class': 'Ob2 = %s S^V_2(xi) + %s S^perp_2(xi) (no exact functionals at degree -1, weight 9/2)' % (out['T1'], out['T2'])}
        gates.append(gate('G-WZ2 (Delta,xi%s), two legs, closedness' % ('+' if chir == 'P' else '-'), out,
                          'Q Ob2(Delta,xi)|2 = -Ob1(xi)|2 + {Ob1(Delta)|2, S(xi)} - ...: class in span{S^V_2, S^perp_2}',
                          closed and out['T0'] == '0', 'normalisation (-1)^{|Delta|} = +1, ghosts theta.. c^xi c^Delta'))
    res['open pairs'] = ['(K_b,P_b\') finite T', '(xi+,xi-) and (xi,xi\'), (xi,eta): S^(2)(x,y) vertex (c2-S2) not built here',
                         '(eta,eta\') finite T', '(xi,P_b)']
    return res, gates


# ================================================================================================ T2 (Delta, P_b: three legs, partial)
def task_T2():
    res, gates = {}, []
    rng, rp, rpol = rngs(71)
    secs = [(('c', 'c', 'cs'), ['cccS']), (('c', 'AV', 'AsV'), ['AcAs']), (('c', 'AW', 'AsW'), ['AcAs']),
            (('c', 'lamP', 'lamsP'), ['cllS']), (('c', 'lamM', 'lamsM'), ['cllS']), (('AV', 'lamP', 'lamM'), ['llA']),
            (('AW', 'lamP', 'lamP'), ['llA']), (('AW', 'lamM', 'lamM'), ['llA']), (('AV', 'AV', 'AV'), ['AAA']),
            (('AV', 'AW', 'AW'), ['AAA'])]
    for nm, spec in (('Delta', {'kind': 'Delta'}), ('P_b', {'kind': 'P', 'b': B4})):
        out = {}
        for types, names in secs:
            S, betanz = [], False
            for _ in range(3):
                p0, p1 = rp(), rp(); moms = [p0, p1, -p0 - p1]
                L = [SA.leg(t, rpol(t), m) for t, m in zip(types, moms)]
                d = SO.density(spec, L)
                S.append(([MD.value3(n, L) for n in names], d['fin']['alpha'].get('(Y0,[Y1,Y2])', V0())))
                betanz = betanz or any(not SA.viszero(d['fin']['beta'][m].get('(Y0,[Y1,Y2])', V0())) for m in range(4))
            key = str(types)
            try:
                c, _ = SF.fit(S, len(names))
                out[key] = {'basis (classical cubic structure, s2_ga)': names, 'coef': [SF.vstr(x) for x in c],
                            'beta': 'nonzero' if betanz else 'zero'}
            except ValueError:
                out[key] = {'basis': names, 'fit': 'not in the span of the classical structure (needs the full cubic basis)',
                            'table (basis value, alpha)': [[str(SF.q2(v)) for v in s[0]] + [SF.vstr(s[1])] for s in S],
                            'beta': 'nonzero' if betanz else 'zero'}
        res['Ob1_%s_3' % nm] = out
    return res, gates


# ================================================================================================ T5 (xi, three legs: tables)
def task_T5():
    res, gates = {}, []
    rng, rp, rpol = rngs(81)
    xi = SO.xi_generic('P'); spec = {'kind': 'xi', 'xi': xi}
    out = {}
    for types in (('c', 'lamP', 'AsW'), ('c', 'lamM', 'AsV'), ('lamP', 'lamP', 'lamP'), ('lamP', 'lamM', 'lamM'),
                  ('AV', 'AV', 'lamM'), ('AV', 'AW', 'lamP'), ('AW', 'AW', 'lamM'), ('c', 'AV', 'lamsP'), ('c', 'AW', 'lamsM')):
        p0, p1 = rp(), rp(); moms = [p0, p1, -p0 - p1]
        pols = [rpol(t) for t in types]
        L = [SA.leg(t, pv, m) for t, pv, m in zip(types, pols, moms)]
        d = SO.density(spec, L, what=('fin', ('pole', 1)))
        out[str(types)] = {'momenta': [m.tolist() for m in moms], 'pols': [[str(x) for x in pv] for pv in pols],
                           'alpha (Y0,[Y1,Y2])': SF.vstr(d['fin']['alpha'].get('(Y0,[Y1,Y2])', V0())),
                           'pole a_1': SF.vstr(d[('pole', 1)]['alpha'].get('(Y0,[Y1,Y2])', V0()))}
    res['Ob1_xi+_3_table'] = {'xi': [str(x) for x in xi], 'rows': out,
                              'note': 'seeded plane-wave values (colour (Y0,[Y1,Y2])); basis decomposition not done'}
    poles0 = all(v['pole a_1'] == '0' for v in out.values())
    gates.append(gate('G-pole (xi, three legs)', 'all zero' if poles0 else 'NONZERO', 'no poles beyond two legs', poles0, ''))
    return res, gates


# ================================================================================================ T8 (log coefficients, partial)
def task_T8():
    res, gates = {}, []
    # along xi: lc1(xi)|2 on T0..T4 and {N, S2(xi)}
    l = {}
    for chir in ('P', 'M'):
        xi = SO.xi_generic(chir)
        lc = fit_T({'kind': 'xi', 'xi': xi}, xi, chir, what='log', seed=91)
        Tb = Tbasis(xi)
        def Nbr(mm, xi=xi):
            N = mm.M.quad('AAs', mm.pairs) * 2 - mm.M.quad('lsl', mm.pairs) * 4 - mm.M.quad('csc', mm.pairs)
            S2 = MD.F_T1(mm, xi) + MD.F_T1(mm, xi, range(4, 10)) + MD.F_T3(mm, xi) + MD.F_T3(mm, xi, range(4, 10))
            return mm.br(N, S2, 1)
        NB = mbasis(Nbr)
        rng, rp, rpol = rngs(93)
        o = 'M' if chir == 'P' else 'P'
        nb = {}
        for t1, t2, names in (('AV', 'lams' + chir, ['T0', 'T3']), ('AW', 'lams' + o, ['T4']), ('lam' + chir, 'AsW', ['T2']),
                              ('lam' + o, 'AsV', ['T1'])):
            S = []
            for _ in range(len(names) + 1):
                p = rp(); L = [SA.leg(t1, rpol(t1), p), SA.leg(t2, rpol(t2), -p)]
                v = NB(L)
                S.append(([Tb[n](L) for n in names], [v] + [(0, 0)] * 4))
            c, _ = SF.fit(S, len(names))
            for n, x in zip(names, c):
                nb[n] = SF.vstr(x)
        tot = {k: SF.vstr([a + b for a, b in zip(SF_p(lc[k]), SF_p(nb[k]))]) for k in lc}
        key = 'xi' + ('+' if chir == 'P' else '-')
        l[key] = tot
        res['lc1_%s_2' % key] = {'lc1 on T0..T4': lc, '{N,S2(xi)} on T0..T4': nb, 'lc1 + {uN,S(xi)}': tot}
    def lval(t):
        ok = _same(t['T1'], t['T3']) and _same(t['T2'], t['T4']) and _same(t['T1'], t['T2']) and t['T0'] == '0'
        return t['T1'] if ok else None
    lp, lm = lval(l['xi+']), lval(l['xi-'])
    res['G-logex'] = {'l_xi+': lp, 'l_xi-': lm, 'rule': 'l_x = c lambda_x + a B_x, lambda_xi = -1, B_xi+ = +1, B_xi- = -1'}
    if lp is not None and lm is not None:
        a = [(x - y) / 2 for x, y in zip(SF_p(lp), SF_p(lm))]; c = [-(x + y) / 2 for x, y in zip(SF_p(lp), SF_p(lm))]
        res['G-logex'].update({'c (inner)': SF.vstr(c), 'a (outer, along [G_B])': SF.vstr(a)})
    gates.append(gate('G-logex (xi, two legs)', res['G-logex'], 'lc1(xi) + {uN,S(xi)} = l_xi u S(xi) (no exact functionals)',
                      lp is not None and lm is not None, 'plain cycles (no special edge), log eps coefficient'))
    return res, gates


def SF_p(s):
    from sa_tasks import SF_parse
    return SF_parse(s)
