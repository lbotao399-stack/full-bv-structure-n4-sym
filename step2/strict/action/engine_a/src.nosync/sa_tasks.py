"""Tasks of engine_spec.json (route A).  Each task returns (results dict, gates list).  Coefficients in units u."""
import itertools
from fractions import Fraction
import numpy as np
import sa_core as SA
import sa_ob as SO
import sa_fit as SF
import sa_mode as MD
from r2_core import Q2

B4 = np.array([1, -1, 2, 1])            # generic translation direction b
V0 = SA.vzero


def rngs(seed):
    rng = np.random.default_rng(seed)
    def rp():
        while True:
            p = rng.integers(-2, 3, size=4)
            if p @ p and not (p == 0).sum() > 2:
                return p
    def rpol(t):
        B = SA.E.BASIS[t]
        while True:
            c = rng.integers(-2, 3, size=len(B))
            if c.any():
                return sum(ci * np.asarray(v) for ci, v in zip(c, B)).astype(complex)
    return rng, rp, rpol


def gate(gid, value, target, ok, notes=''):
    return {'id': gid, 'value': value, 'target': target, 'pass': bool(ok), 'notes': notes}


def two_leg_fit(spec, t1, t2, bfs, n=4, seed=1, other=None, sgn=1, poles=False, pfs=None, mode_bfs=None):
    """fit alpha (p = 0 intrinsic) of the cycle sum along spec (or (spec, other)) on two legs (t1,t2) onto the basis
    functions bfs (callables of the legs); returns coefficient strings, raw coefficients, beta-zero flag, pole fit"""
    rng, rp, rpol = rngs(seed)
    S, P, bz = [], [], True
    what = ('fin', ('pole', 1)) if poles else ('fin',)
    for _ in range(n):
        p = rp()
        L = [SA.leg(t1, rpol(t1), p), SA.leg(t2, rpol(t2), -p)]
        d = SO.density(spec, L, other=other, what=what)
        a = SA.vscale(d['fin']['alpha'].get('(Y0,Y1)', V0()), sgn)
        bz = bz and all(SA.viszero(d['fin']['beta'][m].get('(Y0,Y1)', V0())) for m in range(4))
        S.append(([f(L) for f in bfs], a))
        if poles:
            P.append(([f(L) for f in pfs], SA.vscale(d[('pole', 1)]['alpha'].get('(Y0,Y1)', V0()), sgn)))
    c, rank = SF.fit(S, len(bfs))
    out = {'coef': [SF.vstr(x) for x in c], 'rank': rank, 'beta_zero': bz, 'raw': c}
    if poles:
        cp, _ = SF.fit(P, len(pfs))
        out['pole'] = [SF.vstr(x) for x in cp]; out['pole_raw'] = cp
    return out


def zero_check(spec, pairs, seed=9, other=None):
    rng, rp, rpol = rngs(seed)
    out = {}
    for (t1, t2) in pairs:
        p = rp()
        L = [SA.leg(t1, rpol(t1), p), SA.leg(t2, rpol(t2), -p)]
        d = SO.density(spec, L, other=other)
        tot = d['fin']['alpha'].get('(Y0,Y1)', V0())
        bz = all(SA.viszero(d['fin']['beta'][m].get('(Y0,Y1)', V0())) for m in range(4))
        out['(%s,%s)' % (t1, t2)] = 'alpha=%s, beta=%s' % (SF.vstr(tot), '0' if bz else 'nonzero')
    return out


def q2s(x):
    return str(SF.q2(x))


# ================================================================================================ T0
def task_T0():
    res, gates = {}, []
    # G-reg0, two legs: (c, A_V); pole against -Q I^CT_2 = (9/2) eps^-1 int(dc, A_V)
    b_dcAV = lambda L: 1j * L[0]['pol'][0] * (L[0]['p0'] @ L[1]['pol'][:4])
    f = two_leg_fit({'kind': 'none'}, 'c', 'AV', [SF.b_ddcdA], poles=True, pfs=[b_dcAV]) if False else None
    rng, rp, rpol = rngs(3)
    S, P = [], []
    for _ in range(3):
        p = rp(); L = [SA.leg('c', rpol('c'), p), SA.leg('AV', rpol('AV'), -p)]
        r = SA.cycle_sum(L, [], what=('fin', ('pole', 1)))
        v = SO.reduce_colour(r, 2, 'fin', 0).get('(Y0,Y1)', V0()); vp = SO.reduce_colour(r, 2, ('pole', 1), 0).get('(Y0,Y1)', V0())
        S.append(([SF.b_ddcdA(L)], v)); P.append(([b_dcAV(L)], vp))
    c, _ = SF.fit(S, 1); cp, _ = SF.fit(P, 1)
    res['Ob0_2'] = {'int(d*dc,d*A_V)': SF.vstr(c[0]), 'special-edge pole eps^-1 on int(dc,A_V)': SF.vstr(cp[0])}
    gates.append(gate('G-reg0 (two legs)', SF.vstr(c[0]), '-27/8', SF.vstr(c[0]) == '-27/8',
                      'arity-zero special-edge cycles (c, A_V); pole ' + SF.vstr(cp[0]) + ' = -(9/2), cancelled by -Q I^CT_2'))
    # G-reg0, three legs (cubic obstruction) on G1..G5 (mode-model values, colours e0, e1, e2)
    r3 = reg0_three()
    res['Ob0_3'] = r3['coef']
    tgt = {'G1': '11/4 + (19/4)*log3 + (-19/4)*log2', 'G3': '11/4 + (19/4)*log3 + (-19/4)*log2', 'G2': '27/8', 'G4': '0',
           'G5': '(9)*log2 + (-9)*log3'}
    okk = all(_same(r3['coef'][k], tgt[k]) for k in tgt)
    gates.append(gate('G-reg0 (three legs)', r3['coef'], '(11/4+19/4 L)(G1+G3) + 27/8 G2 - 9 L G5, L = log(3/2)', okk,
                      'cubic special-edge cycles, colour (Y0,[Y1,Y2]) from tr(ad ad ad) = (C/2)(Y_c,[Y_b,Y_a])'))
    # G-sign: derivation rule of Lambda on brackets in the mode model, and Q S_2(x) = 0
    res['G-sign'] = sign_tests()
    gates.append(gate('G-sign', res['G-sign']['summary'], 'component formulas consistent with s^2 = 0',
                      res['G-sign']['ok'], res['G-sign']['note']))
    return res, gates


def _same(a, b):
    def parse(s):
        d = {}
        for part in s.split(' + '):
            part = part.strip()
            if part == '0':
                continue
            if '*' in part and part.split('*')[-1] in ('log2', 'log3', 'log5', 'log7'):
                nm = part.split('*')[-1]; c = part[:-(len(nm) + 1)].strip('()')
            else:
                nm = '1'; c = part
            d[nm] = d.get(nm, Fraction(0)) + Fraction(c)
        return {k: v for k, v in d.items() if v}
    return parse(a) == parse(b)


def reg0_three():
    """arity-zero three-leg cycles on G1..G5 (s2_ga definitions)"""
    rng, rp, rpol = rngs(11)
    secs = {('c', 'AV', 'AV'): ['G1', 'G2'], ('c', 'AW', 'AW'): ['G3'], ('c', 'c', 'AsV'): ['G4'],
            ('c', 'lamP', 'lamM'): ['G5'], ('c', 'lamM', 'lamP'): ['G5']}
    coef = {}
    for types, names in secs.items():
        S = []
        for _ in range(len(names) + 2):
            p0, p1 = rp(), rp()
            moms = [p0, p1, -p0 - p1]
            if not (moms[2] @ moms[2]):
                continue
            L = [SA.leg(t, rpol(t), m) for t, m in zip(types, moms)]
            r = SA.cycle_sum(L, [], what=('fin',))
            v = SO.reduce_colour(r, 3, 'fin', 0).get('(Y0,[Y1,Y2])', V0())
            S.append(([MD.value3(nm, L) for nm in names], v))
        c, _ = SF.fit(S, len(names))
        for nm, x in zip(names, c):
            if nm in coef:
                assert _same(coef[nm], SF.vstr(x)), (nm, coef[nm], SF.vstr(x))
            coef[nm] = SF.vstr(x)
    return {'coef': coef}


def sign_tests():
    p = np.array([1, 2, -1, 1]); mm = MD.MM(p); xi = SO.xi_generic('P')
    T1 = MD.F_T1(mm, xi); T3 = MD.F_T3(mm, xi); T2 = MD.F_T1(mm, xi, range(4, 10)); T4 = MD.F_T3(mm, xi, range(4, 10))
    a = MD.F_named(mm, 'dsAdsA')
    LV = [('c', np.array([1])), ('l', SA.E.BASIS['lamM'][0].astype(complex) + 2 * SA.E.BASIS['lamM'][2])]
    lhs = mm.value(mm.Q(mm.br(a, T1, 0)), LV)
    t1 = mm.value(mm.br(mm.Q(a), T1, 1), LV)
    c07 = (lhs + t1).iszero()               # Lambda{a,b} = -{Lambda a, b} - (-1)^{|a|}{a, Lambda b}
    std = (lhs - t1).iszero()               # Lambda{a,b} = {Lambda a, b} + (-1)^{|a|+1}{a, Lambda b}
    LA = [('A', np.array([1, 0, 2, -1, 0, 0, 0, 0, 0, 0], dtype=complex)), LV[1]]
    LP = [('A', np.array([0, 0, 0, 0, 1, 2, 0, -1, 0, 1], dtype=complex)),
          ('l', SA.E.BASIS['lamP'][1].astype(complex) - SA.E.BASIS['lamP'][4])]
    qsv = mm.value(mm.Q(T1 + T3), LA); qsp = mm.value(mm.Q(T2 + T4), LP)
    return {'Q{a,T1} on (c,lam)': str(lhs), '{Qa,T1}': str(t1), 'C07 rule holds': c07, 'standard rule holds': std,
            'Q(T1+T3)': str(qsv), 'Q(T2+T4)': str(qsp), 'ok': c07 and qsv.iszero() and qsp.iszero(),
            'summary': 'Lambda{a,b} = -{Lambda a,b} - (-1)^{|a|}{a,Lambda b} (C07 bracket of Def c1-bracket); '
                       'Q S^V_2(xi) = Q S^perp_2(xi) = 0',
            'note': 'with this rule s^2 = 0 forces (sG)^(1)(x) = Lambda G^(1)(x) + {G^(0),S(x)} for degree-1 G, i.e. '
                    'Lambda Ob^(1)(x) = -{Ob^(0),S(x)}; Lemma at-sformulas (ii) prints +{G^(0),S(x)} (derived with the '
                    'derivation rule of the standard bracket (-1)^{|f|}{f,g}).  Degree-0 formula (i) is unaffected.'}


# ================================================================================================ T1 (Delta, two legs)
def task_T1():
    res, gates = {}, []
    spec = {'kind': 'Delta'}
    secs = [('AV', 'AV', ['int(dA_V,dA_V)', 'int(d*A_V,d*A_V)'], [SF.b_dAVdAV, SF.b_dsAdsA]),
            ('AW', 'AW', ['int(dA_perp,dA_perp)'], [SF.b_dAWdAW]),
            ('lamP', 'lamM', ['int(lam,dslash lam)'], [SF.b_ldl]),
            ('c', 'AsV', ['int(dc,A*)'], [SF.b_dcAs])]
    out = {}
    allz = True
    for t1, t2, names, bfs in secs:
        pfs = [SF.b_AA] if t1[0] == 'A' and t2[0] == 'A' else [lambda L: 0]
        r = two_leg_fit(spec, t1, t2, bfs, poles=True, pfs=pfs)
        out['(%s,%s)' % (t1, t2)] = {'basis': names, 'coef': r['coef'], 'beta (position part)': 'zero' if r['beta_zero'] else 'NONZERO',
                                     'pole a_1 on int(A,A)': r['pole'] if t1[0] == 'A' else 'n/a (zero)'}
        allz = allz and r['beta_zero']
    res['Ob1_Delta_2'] = out
    res['Ob1_Delta_2_zero_sectors'] = zero_check(spec, [('AV', 'AW'), ('lamP', 'lamP'), ('c', 'AsW'), ('AV', 'lamP'), ('c', 'c')])
    got = [out['(AV,AV)']['coef'], out['(AW,AW)']['coef'], out['(lamP,lamM)']['coef'], out['(c,AsV)']['coef']]
    ok = got == [['-11/3', '0'], ['-11/3'], ['47/3'], ['-14/3']]
    gates.append(gate('G-reg1', {'int(dA,dA)': got[0][0] + ' (A_V), ' + got[1][0] + ' (A_perp)', 'int(d*A_V,d*A_V)': got[0][1],
                                 'int(lam,dslash lam)': got[2][0], 'int(dc,A*)': got[3][0], 'explicit x': 'zero' if allz else 'NONZERO'},
                      'u[-11/3 int(dA,dA) + 47/3 int(lam,dslash lam) - 14/3 int(dc,A*)], 0 for int(d*A_V,d*A_V) and x', ok and allz,
                      'special-edge cycles only; component normalisation Def at-names (i), ghost of Delta on the right'))
    pl = [out['(AV,AV)']['pole a_1 on int(A,A)'], out['(AW,AW)']['pole a_1 on int(A,A)']]
    gates.append(gate('G-pole (Delta, two legs)', pl, '(9/4){int(A,A),S(Delta)} = (9/4) v.int(A,A) = -9 int(A,A)',
                      pl[0][-1] == '-9' and pl[1][-1] == '-9',
                      '{S_X,F} = +D_X F measured in the Part-D mode model (d4_signs), so {int(A,A),S(Delta)} = -4 int(A,A); '
                      'lambda and ghost sectors: pole 0'))
    # P_b two legs (G-P2), direct
    specP = {'kind': 'P', 'b': B4}
    zp = zero_check(specP, [('AV', 'AV'), ('AW', 'AW'), ('lamP', 'lamM'), ('c', 'AsV')], seed=13)
    res['Ob1_P_2'] = zp
    okp = all(v.startswith('alpha=0,') for v in zp.values())
    gates.append(gate('G-P2', zp, 'Ob^(1)(P_b)|_2 = 0', okp, 'direct cycles along P_b; also beta(Delta) = 0 (position part = 2 Ob(P_b))'))
    # X-scal cross-check with the plain-cycle log coefficients (lc^(0) on two legs)
    lc = lc0_two()
    res['lc0_2'] = lc
    gates.append(gate('X-scal (cross-check)', {'lc0': lc, 'Ob1(Delta)|2 - 4 lc0': 'composite part'},
                      'Ob1(Delta)|2 = 4 lc0 + u[1/3, -1/3, -2/3]', lc == {'int(dA,dA)': '-1', 'int(lam,dslash lam)': '4', 'int(dc,A*)': '-1'},
                      'cross-check only (Thm at-second (v)); lc0 from the plain arity-zero cycles (log eps coefficient)'))
    return res, gates


def lc0_two():
    rng, rp, rpol = rngs(17)
    out = {}
    for t1, t2, nm, bf in (('AW', 'AW', 'int(dA,dA)', SF.b_dAWdAW), ('lamP', 'lamM', 'int(lam,dslash lam)', SF.b_ldl),
                           ('c', 'AsV', 'int(dc,A*)', SF.b_dcAs)):
        S = []
        for _ in range(2):
            p = rp(); L = [SA.leg(t1, rpol(t1), p), SA.leg(t2, rpol(t2), -p)]
            r = SA.cycle_sum(L, [], what=('log',))
            S.append(([bf(L)], SO.reduce_colour(r, 2, 'log', 0).get('(Y0,Y1)', V0())))
        c, _ = SF.fit(S, 1)
        out[nm] = SF.vstr(c[0])
    return out


# ================================================================================================ T4 (xi, two legs)
def task_T4():
    res, gates = {}, []
    p = np.array([1, 2, -1, 1]); mm = MD.MM(p)
    for chir in ('P', 'M'):
        o = 'M' if chir == 'P' else 'P'
        xi = SO.xi_generic(chir); spec = {'kind': 'xi', 'xi': xi}
        secs = [('AV', 'lam' + o, ['g1', 'g2'], [lambda L: SF.b_g1(L, xi), lambda L: SF.b_g2(L, xi)], [lambda L: SF.b_GVA(L, xi)]),
                ('AW', 'lam' + chir, ['g3'], [lambda L: SF.b_g3(L, xi)], [lambda L: SF.b_GPA(L, xi)]),
                ('c', 'lams' + chir, ['g4'], [lambda L: SF.b_g4(L, xi)], [lambda L: 0])]
        ob, raw, pole = {}, {}, {}
        for t1, t2, names, bfs, pfs in secs:
            r = two_leg_fit(spec, t1, t2, bfs, poles=True, pfs=pfs, seed=21)
            for nm, c, cr in zip(names, r['coef'], r['raw']):
                ob[nm] = c; raw[nm] = cr
            pole['(%s,%s)' % (t1, t2)] = r['pole']
        key = 'xi_' + ('+' if chir == 'P' else '-')
        res['Ob1_%s_2' % key] = {'o1 (g1)': ob['g1'], 'o2 (g2)': ob['g2'], 'o3 (g3)': ob['g3'], 'o4 (g4)': ob['g4'],
                                 'pole a_1': {'int(Gamma_V,A_V)': pole['(AV,lam%s)' % o][0], 'int(Gamma_perp,A_perp)': pole['(AW,lam%s)' % chir][0]},
                                 'zero sectors': zero_check(spec, [('AV', 'lam' + chir), ('AW', 'lam' + o), ('c', 'lams' + o)], seed=23)}
        # mode model: Q T_i, brackets, WZ
        T = {'T1': MD.F_T1(mm, xi), 'T2': MD.F_T1(mm, xi, range(4, 10)), 'T3': MD.F_T3(mm, xi), 'T4': MD.F_T3(mm, xi, range(4, 10)),
             'T0': MD.F_T0(mm, xi)}
        g = {k: MD.F_g(mm, xi, k) for k in (1, 2, 3, 4)}
        _, _, rpol = rngs(29)
        LV = [('A', rpol('AV')), ('l', rpol('lam' + o))]; LP = [('A', rpol('AW')), ('l', rpol('lam' + chir))]
        Lc = [('c', np.array([1])), ('ls', rpol('lams' + chir))]
        rat = lambda F, G_, L: mm.value(F, L) / mm.value(G_, L)
        q1 = rat(mm.Q(T['T1']), g[1], LV); q3 = rat(mm.Q(T['T2']), g[3], LP)
        q1b = rat(mm.Q(T['T3']), g[1], LV); q3b = rat(mm.Q(T['T4']), g[3], LP)
        a_ = rat(mm.Q(T['T0']), g[2], LV); b_ = rat(mm.Q(T['T0']), g[4], Lc)
        AA = MD.F_named(mm, 'AA')
        bAA = rat(mm.br(AA, T['T1'], 0), MD.F_GA(mm, xi, range(4)), LV)
        dsT1 = rat(mm.br(MD.F_named(mm, 'dsAdsA'), T['T1'], 0), g[2], LV)
        # WZ on (c, lam): Q Ob1 vs {Ob0|2, T1}
        Lcl = [('c', np.array([1])), ('l', rpol('lam' + o))]
        Ob1 = None
        for k, nm in ((1, 'g1'), (2, 'g2'), (3, 'g3'), (4, 'g4')):
            cc = raw[nm]
            assert all(SF.q2(x).iszero() for x in cc[1:]) or nm in ('g1', 'g3')
        o2 = SF.q2(raw['g2'][0]); o4 = SF.q2(raw['g4'][0])
        lhs = mm.value(mm.Q(g[2]), Lcl) * o2 + mm.value(mm.Q(g[4]), Lcl) * o4
        rhs = mm.value(mm.br(MD.F_named(mm, 'ddcdA'), T['T1'], 1), Lcl) * Q2(Fraction(-27, 8))
        wz_c07 = (lhs + rhs).iszero(); wz_print = (lhs - rhs).iszero()
        res['S1_xi_' + key] = {'q1 (QT1 = q1 g1)': str(q1), 'QT3/g1': str(q1b), 'q3 (QT2 = q3 g3)': str(q3), 'QT4/g3': str(q3b),
                               'a (QT0 = a g2 + b g4)': str(a_), 'b': str(b_),
                               '{int(A,A),T1}/int(Gamma_V,A_V)': str(bAA), '{int(d*A_V,d*A_V),T1}/g2': str(dsT1),
                               'H^0(Q) on two legs': '0 (closed: g1, g3, a g2 + b g4 = image)',
                               'kernel of Q at degree -1': 'span{T1+T3, T2+T4}'}
        gates.append(gate('G-WZ1 (xi%s, two legs)' % ('+' if chir == 'P' else '-'),
                          {'Q Ob1|2 on (c,lam)': str(lhs), '{Ob0|2,T1}': str(rhs)},
                          'Lambda Ob1(x) = -{Ob0,S(x)} (C07 bracket; the printed Lemma at-sformulas sign gives +)', wz_c07,
                          'printed-sign version holds: %s' % wz_print))
        pv = res['Ob1_%s_2' % key]['pole a_1']
        gates.append(gate('G-pole (xi%s, two legs)' % ('+' if chir == 'P' else '-'), pv,
                          '(9/4){int(A,A),T1+T2} = (9/2)[int(Gamma_V,A_V) + int(Gamma_perp,A_perp)]',
                          pv['int(Gamma_V,A_V)'] == '9/2' and pv['int(Gamma_perp,A_perp)'] == '9/2' and str(bAA) == '2',
                          '{int(A,A),T1} = %s int(Gamma_V,A_V) (mode model)' % bAA))
        # residual and particular solution: R = Ob1 - {F0|2, T1}, F0|2 = 27/16 u int(d*A_V,d*A_V) + beta (Q-exact)
        r1 = [SF.q2(x) for x in raw['g1']]; r3 = [SF.q2(x) for x in raw['g3']]
        r2 = [SF.q2(x) for x in raw['g2']]; r4 = [SF.q2(x) for x in raw['g4']]
        r2[0] = r2[0] - Q2(Fraction(27, 16)) * dsT1
        # closedness: (r2, r4) = r0 (a, b)
        r0 = r2[0] / a_
        closed = (r4[0] - r0 * b_).iszero() and all(x.iszero() for x in r2[1:] + r4[1:])
        F1 = {'T1': SF.vstr([x / q1 for x in r1]), 'T2': SF.vstr([x / q3 for x in r3]), 'T0': str(r0),
              'free': 'u (theta^xi_V S^V_2(xi) + theta^xi_perp S^perp_2(xi)), S^V_2 = T1+T3, S^perp_2 = T2+T4'}
        res['F1_' + key + '_2'] = {'R1|2 (residual) on g1..g4': [SF.vstr(r1), SF.vstr(r2), SF.vstr(r3), SF.vstr(r4)],
                                   'residual closed': closed, 'particular solution F1(xi)|2': F1,
                                   'beta-directions': 'F1 shifted by {N_f, S_2(xi)} (Prop at-xi2 (vi)), symbolic',
                                   'equivalent form': 'mod free directions: F1 = -2 c (T1+T2) = +2 c (T3+T4), c = 79/24 - 1/2 log(3/2)'}
    return res, gates


# ================================================================================================ B1 (Delta, P_b)
def task_B1():
    res, gates = {}, []
    b = B4
    spec, other = {'kind': 'Delta'}, {'kind': 'P', 'b': b}
    sp = lambda p: sum(p[m] * SF.SB[m] for m in range(4))
    ss = lambda v: sum(v[m] * SF.SS[m] for m in range(4))
    def L01(L):
        return L[0]['pol'], L[0]['p0'], L[1]['pol'], L[1]['p0']
    bas = {
        ('AW', 'AsW'): [('int(d_b A_perp, A*_perp) = S(P_b)_perp', lambda L: 1j * (b @ L[0]['p0']) * (L[0]['pol'][4:] @ L[1]['pol'][4:]))],
        ('AV', 'AsV'): [('int(d_b A_V, A*_V)', lambda L: 1j * (b @ L[0]['p0']) * (L[0]['pol'][:4] @ L[1]['pol'][:4])),
                        ('int(d.A_V, b.A*_V)', lambda L: 1j * (L[0]['p0'] @ L[0]['pol'][:4]) * (b @ L[1]['pol'][:4])),
                        ('int(d_mu (b.A_V), A*_mu)', lambda L: 1j * (L[0]['p0'] @ L[1]['pol'][:4]) * (b @ L[0]['pol'][:4]))],
        ('c', 'cs'): [('int(d_b c, c*)', lambda L: 1j * (b @ L[0]['p0']) * L[0]['pol'][0] * L[1]['pol'][0])],
        ('lamP', 'lamsP'): [('int(d_b lam, lam*)', lambda L: 1j * (b @ L[0]['p0']) * (L[0]['pol'] @ L[1]['pol'])),
                            ('int(rho(b) dslash lam, lam*)', lambda L: 1j * ((ss(b) @ (sp(L[0]['p0']) @ L[0]['pol'])) @ L[1]['pol']))],
        ('lamM', 'lamsM'): None}
    bas[('lamM', 'lamsM')] = bas[('lamP', 'lamsP')]
    fits = {}
    for (t1, t2), bl in bas.items():
        r = two_leg_fit(spec, t1, t2, [f for _, f in bl], other=other, n=len(bl) + 2, seed=31, poles=True, pfs=[f for _, f in bl])
        fits[(t1, t2)] = r
        res.setdefault('Ob2_Delta_P_2', {})['(%s,%s)' % (t1, t2)] = {'basis': [nm for nm, _ in bl], 'coef': r['coef'],
                                                                     'beta (position part)': 'zero' if r['beta_zero'] else 'NONZERO',
                                                                     'poles a_1': r['pole']}
    res['Ob2_Delta_P_2_zero_sectors'] = zero_check(spec, [('AV', 'AsW'), ('lamP', 'lamsM'), ('c', 'AsV'), ('AV', 'AV')], other=other, seed=33)
    # exact directions (mode model): Q int(rho(b) lam*, lam*), Q int(iota_b A_V, c*)
    p = np.array([1, 2, -1, 1]); mm = MD.MM(p)
    sb_ = ss(b)
    Fl = mm.quad('ls', 'ls', lambda i, j, p1, p2: sb_[j, i])
    FV = mm.quad('A', 'cs', lambda i, j, p1, p2: b[i] if i < 4 else 0)
    QFl, QFV = mm.Q(Fl), mm.Q(FV)
    _, rp, rpol = rngs(37)
    def proj(QF, t1, t2, bl, n):
        S = []
        for _ in range(n):
            pol1, pol2 = rpol(t1), rpol(t2)
            L = [SA.leg(t1, pol1, p), SA.leg(t2, pol2, -p)]
            v = mm.value(QF, [(SA.E.SLOT[t1], pol1), (SA.E.SLOT[t2], pol2)])
            S.append(([f(L) for _, f in bl], [(v.re, v.im)] + [(0, 0)] * 4))
        c, _ = SF.fit(S, len(bl))
        return [SF.q2(x[0]) for x in c]
    out = {}
    # perp
    o_perp = [SF.q2(x) for x in fits[('AW', 'AsW')]['raw'][0]]
    out['perp'] = SF.vstr(o_perp)
    # lambda (both chiralities): Ob = c1 B1 + c2 B2 = o_lam (-B1) + y E
    for t1, t2 in (('lamP', 'lamsP'), ('lamM', 'lamsM')):
        e = proj(QFl, t1, t2, bas[(t1, t2)], 4)
        c1 = [SF.q2(x) for x in fits[(t1, t2)]['raw'][0]]; c2 = [SF.q2(x) for x in fits[(t1, t2)]['raw'][1]]
        y = [x / e[1] for x in c2]
        o = [-(a - yy * e[0]) for a, yy in zip(c1, y)]
        out['lambda (%s)' % t1] = SF.vstr(o)
        out['lambda (%s) exact part y' % t1] = SF.vstr(y)
        out['Q int(rho(b)lam*,lam*) on (%s)' % t1] = [str(x) for x in e]
    # V: Ob = cV (B1,B2,B3) on (A_V,A*_V) and cc on (c,c*) = o_V (B1 - Bc) + y E_V
    eV = proj(QFV, 'AV', 'AsV', bas[('AV', 'AsV')], 5); ec = proj(QFV, 'c', 'cs', bas[('c', 'cs')], 2)
    cV = [[SF.q2(x) for x in fits[('AV', 'AsV')]['raw'][k]] for k in range(3)]
    cc = [SF.q2(x) for x in fits[('c', 'cs')]['raw'][0]]
    # unknowns o_V, y per log component: equations B1: cV0 = o + y eV0; B2: cV1 = y eV1; B3: cV2 = y eV2; Bc: cc = -o + y ec
    oV, yV, consistent = [], [], True
    for k in range(5):
        rows = [[1, eV[0]], [0, eV[1]], [0, eV[2]], [-1, ec[0]]]
        rhs = [cV[0][k], cV[1][k], cV[2][k], cc[k]]
        x, rank, free = SF.solve(rows, rhs)
        if x is None:
            consistent = False; x = [Q2(0), Q2(0)]
        oV.append(x[0]); yV.append(x[1])
    out['V'] = SF.vstr(oV); out['V exact part y'] = SF.vstr(yV); out['V consistent (class in span)'] = consistent
    out['Q int(iota_b A_V,c*) on (A_V,A*_V) basis'] = [str(x) for x in eV]; out['... on (c,c*)'] = [str(x) for x in ec]
    res['B1_classes'] = out
    th = {d: SF.vstr([SF.q2(x) / 2 for x in SF_parse(out[k])]) for d, k in (('perp', 'perp'), ('lambda', 'lambda (lamP)'), ('V', 'V'))}
    res['theta_Delta'] = {'theta_' + d: v for d, v in th.items()}
    res['theta_Delta']['rule'] = 'theta^Delta_d = o_d/(2u) (Thm at-theta); o_d = coefficient of [S(P_b)_d] in [Ob^(2)(Delta,P_b)|_2]'
    lamsame = out['lambda (lamP)'] == out['lambda (lamM)']
    gates.append(gate('X-theta (cross-check)', out['perp'], 'o_perp = -2/3', out['perp'] == '-2/3',
                      'cross-check only; theta = o_perp/(2u)'))
    allbeta = all(f['beta_zero'] for f in fits.values())
    poles0 = all(all(x == '0' for x in f['pole']) for f in fits.values())
    gates.append(gate('G-pole (Delta,P_b), arity two', 'all a_1 zero' if poles0 else 'NONZERO', 'sum a_m = 0', poles0, ''))
    gates.append(gate('B1 structure', {'beta (position part) zero': allbeta, 'lambda chiralities equal': lamsame,
                                       'V class in span{S(P_b)_V} + Im Q': consistent},
                      'class of Ob2(Delta,P_b)|2 in span{[S(P_b)_d]}', allbeta and consistent, ''))
    return res, gates


def SF_parse(s):
    """inverse of SF.vstr for real/complex rationals with log terms"""
    out = [Q2(0)] * 5
    if s == '0':
        return out
    names = SF.BASIS_NAMES
    for part in s.split(' + '):
        part = part.strip()
        nm = '1'
        for k in names[1:]:
            if part.endswith('*' + k):
                nm = k; part = part[:-(len(k) + 1)].strip('()')
        idx = names.index(nm)
        if '*i' in part:
            if part.startswith('(') or '+' in part:
                a, bb = part.strip('()').split('+')
                val = Q2(Fraction(a), Fraction(bb.replace('*i', '')))
            else:
                val = Q2(0, Fraction(part.replace('*i', '')))
        else:
            val = Q2(Fraction(part))
        out[idx] = out[idx] + val
    return out
