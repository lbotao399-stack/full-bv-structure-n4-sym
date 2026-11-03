"""Gate G-sign (two-leg part) and the two-leg solution along xi (T4/S2):
(a) the Koszul rule of Q on plane-wave densities: search over (alpha, beta, eta) with Q^2 = 0, Q(T1+T3) = Q(T2+T4) = 0
    and the four images of Theorem c9-quadratic (iii); must be unique;
(b) classical closedness Q S_2(x) = 0 for the sector parts of S(P_b), S(Delta) (translation-invariant part), S(xi);
(c) Q-images along xi: QT1 = q1 g1, QT2 = q3 g3, QT0 = a g2 + b g4, Qg2 = -Qg4;
(d) G-WZ1 along xi (two legs) in the form Lambda Ob^(1)(xi) = -{Ob^(0), S(xi)} (component formula of Def c2-dce for a
    degree-1 cochain) and the particular solution F1^(1)(xi)|_2."""
import json, os, random, itertools
import numpy as np
import sympy as sp
import eb_qb as Q

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')


def run(log=print):
    random.seed(1)
    out = {}
    species = ['c', 'A', 'As', 'lam']
    passing = []
    xi = np.zeros(16, complex); xi[:8] = [random.randint(-2, 2) for _ in range(8)]
    T = Q.xi_terms(xi)
    for a in (0, 1):
        for b in (0, 1):
            for et in itertools.product((1, -1), repeat=4):
                eta = dict(zip(species, et)); eta.update(cs=1, lams=1)
                rule = (a, b, eta)
                QQ = lambda F: Q.apply_Q(F, rule)
                e = [Q.check_identity(QQ(Q.F_AAs_p), Q.dsum((-0.5, Q.F_dAdA_p)), [['A', 'A']]),
                     Q.check_identity(QQ(Q.F_AAs_V), Q.dsum((-0.5, Q.F_dAdA_V), (-1, Q.F_dcAs)), [['A', 'A'], ['c', 'As'], ['As', 'c']]),
                     Q.check_identity(QQ(Q.F_csc), Q.dsum((-1, Q.F_dcAs)), [['c', 'As'], ['As', 'c']]),
                     max(Q.check_identity(QQ(QQ(F)), lambda l: 0j, [s]) for F, s in
                         [(Q.F_AAs_V, ['c', 'As']), (Q.F_AAs_V, ['c', 'A']), (Q.F_csc, ['c', 'A']), (Q.F_lamslam, ['lam', 'lam']),
                          (Q.F_AA, ['c', 'c']), (Q.F_AAs_V, ['A', 'c'])]),
                     Q.check_identity(QQ(Q.dsum((1, T['T1']), (1, T['T3']))), lambda l: 0j, [['lam', 'A'], ['A', 'lam'], ['lam', 'c'], ['c', 'lams']]),
                     Q.check_identity(QQ(Q.dsum((1, T['T2']), (1, T['T4']))), lambda l: 0j, [['lam', 'A'], ['A', 'lam']])]
                if max(e) == 0:
                    passing.append(dict(alpha=a, beta=b, eta=et))
    out['sign_rule_candidates_passing'] = passing
    out['sign_rule_unique'] = len(passing) == 1
    out['sign_rule'] = 'Q acts leg by leg (c->A_V, A->A*, A*->c*, lam->lam*, symbol q(k) of Prop c1-Q) with sign (-1)^{number of odd legs before}; Q int(lam,lam*) = +int(lam,dslash lam), i.e. Q int(lam*,lam) = -int(lam,dslash lam) with (lam*,lam) = -(lam,lam*)'
    log(f'G-sign rule: {passing}')
    # (b) closedness of sector parts (translation part of S(Delta) is S(P_b)-type; weights part checked with x by canon)
    bvec = np.array([1, -2, 2, 1], complex)
    bd = lambda s: np.dot(1j * Q.pv(s[0]), bvec)
    SPV = Q.dsum((1, Q.dens(1, ('A', 'As'), lambda s: bd(s) * np.dot(s[0]['u'][:4], s[1]['u'][:4]))),
                 (-1, Q.dens(1, ('c', 'cs'), lambda s: bd(s) * s[0]['u'][0] * s[1]['u'][0])))
    SPp = Q.dens(1, ('A', 'As'), lambda s: bd(s) * np.dot(s[0]['u'][4:], s[1]['u'][4:]))
    SPl = Q.dens(-1, ('lam', 'lams'), lambda s: bd(s) * (s[0]['u'] @ Q.MP @ s[1]['u']))
    out['Q_SPb_closed'] = {k: Q.check_identity(Q.apply_Q(F), lambda l: 0j, ls) == 0 for k, F, ls in
                           [('V', SPV, [['c', 'As'], ['A', 'A'], ['As', 'c']]), ('perp', SPp, [['A', 'A']]), ('lam', SPl, [['lam', 'lam']])]}
    # (c) images along xi
    g = Q.g_terms(xi)
    LS = [['lam', 'A'], ['A', 'lam']]
    out['QT1'] = [str(x) for x in Q.fit_on(Q.apply_Q(T['T1']), [g['g1'], g['g2']], LS)]
    out['QT3'] = [str(x) for x in Q.fit_on(Q.apply_Q(T['T3']), [g['g1'], g['g2']], LS)]
    out['QT2'] = [str(x) for x in Q.fit_on(Q.apply_Q(T['T2']), [g['g3']], LS)]
    out['QT4'] = [str(x) for x in Q.fit_on(Q.apply_Q(T['T4']), [g['g3']], LS)]
    out['QT0'] = [str(x) for x in Q.fit_on(Q.apply_Q(T['T0']), [g['g2'], g['g4']], LS + [['c', 'lams'], ['lams', 'c']])]
    out['Qg2_over_Qg4'] = [str(x) for x in Q.fit_on(Q.apply_Q(g['g2']), [Q.apply_Q(g['g4'])], [['lam', 'c'], ['c', 'lam']])]
    log(f'Q images along xi: {out}')
    # (d) WZ1 and the particular solution
    Ob0 = Q.dsum((-27 / 8, Q.F_dsdc_dsA))
    lhs = Q.apply_Q(Q.dsum((27 / 8, g['g2'])))          # Q of Ob1(xi)|_2 (g1, g3 closed)
    rhs = Q.apply_vf(Ob0, Q.map_T1(xi))                 # {T1, Ob0} = {Ob0, T1}
    out['WZ1_xi_lhs_plus_rhs'] = Q.check_identity(lhs, Q.dsum((-1, rhs)), [['lam', 'c'], ['c', 'lam']], ntrial=6) if False else None
    out['WZ1_xi_max|Q Ob1 + {Ob0,T1}|'] = Q.check_identity(Q.dsum((1, lhs), (1, rhs)), lambda l: 0j, [['lam', 'c'], ['c', 'lam']], ntrial=6)
    out['WZ1_xi_max|Q Ob1 - {Ob0,T1}|'] = Q.check_identity(Q.dsum((1, lhs), (-1, rhs)), lambda l: 0j, [['lam', 'c'], ['c', 'lam']], ntrial=6)
    F0 = Q.dsum((27 / 16, Q.F_dsA_dsA))
    out['{F1^(0)|_2,T1} on (g2,g1)'] = [str(x) for x in Q.fit_on(Q.apply_vf(F0, Q.map_T1(xi)), [g['g2'], g['g1']], LS)]
    # pole structure {int(A,A), S(xi)} = 2 int(Gamma(xi,lam), A)
    gAG = Q.dens(1, ('lam', 'A'), lambda s: sum((xi @ Q.MR[M] @ s[0]['u']) * s[1]['u'][M] for M in range(10)))
    T12 = Q.dsum((1, T['T1']), (1, T['T2']))
    out['{int(A,A),S(xi)} / int(Gamma,A)'] = [str(x) for x in Q.fit_on(Q.apply_vf(Q.F_AA, lambda leg: (lambda m: m(leg))(lambda lg: None if lg['sp'] != 'lam' else dict(sp='A', u=np.array([xi @ Q.MR[M] @ lg['u'] for M in range(10)]), p=lg['p']))), [gAG], LS)]
    log(f'WZ1 / brackets: {out}')
    return out


if __name__ == '__main__':
    with open(os.path.join(OUT, 'G_sign.json'), 'w') as f:
        json.dump(run(), f, indent=1, default=str)
