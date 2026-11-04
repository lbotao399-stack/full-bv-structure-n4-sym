"""(xi, P_b) on two legs, both routes: class of Ob^(2)(xi,P_b)|_2 in H^{-1}(Q) at weight 11/2 (linear in xi, b),
G-WZ2: Q Ob2 = {Ob1(xi),S(P_b)} = 0 on two legs; residual R2 = Ob2 + {F1(xi),S(P_b)} + {S(xi),F1(P_b)}
((-1)^{|xi|} = -1, [xi,P_b] = 0, S^(2)(xi,P_b) = 0), F1(P_b)|_2 = sum_d theta^P_d S(P_b)_d + exact,
{F1(xi),S(P_b)} = 0 (translation of a translation-invariant functional; computed)."""
import sys, os, json
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import fitlib as FL, bilin as BL, classfit as CF
THP = sp.symbols('thP_V thP_perp thP_lam')


def route_b(chir):
    import cb_common as B, cb_s2 as S, cb_tree as T, cb_xiP as X, cb_pairs as CP
    Q = B.Q
    RL = BL.Real(B.R, B.RB, B.MP)
    raw = json.load(open(os.path.join(B.COMPLETE, 'out.nosync', f'B_xiP_{"p" if chir == "+" else "m"}_raw.json')))
    cvec = lambda a: np.array([complex(x[0], x[1]) for x in a])
    cfg = {}
    for sec, recs in raw.items():
        cfg[sec] = []
        for r in recs:
            legs = [dict(sp=l['sp'], u=cvec(l['u']), p=[Fr(x) for x in l['p']]) for l in r['legs']]
            d = r['raw'].get('k', {})
            ob = FL.clog_parse(d['0']) if '0' in d else FL.lvec_zero()
            pole = FL.clog_parse(d['1']) if '1' in d else FL.lvec_zero()
            cfg[sec].append(dict(legs=legs, ctx=dict(xi=cvec(r['xi']), b=np.array(r['b'], complex)), ob=ob, pole=pole))

    def proj(L, pr):
        def Lp(s, ctx):
            s2 = []
            for leg, rg in zip(s, pr):
                u = np.zeros_like(leg['u']); u[list(rg)] = leg['u'][list(rg)]
                s2.append(dict(sp=leg['sp'], u=u, p=leg['p']))
            return L(s2, ctx)
        return Lp

    def cands_of_sector(sec):
        s0, r0, s1, r1 = X.SECT[chir][sec]
        nm = 2 if s0 == 'A' else 1
        return [BL.Cand(c.label, c.slots, proj(c.L, (tuple(r0), tuple(r1)))) for c in BL.lin_cands(RL, (s0, s1), nm)]
    deg2 = lambda b: BL.lin_cands(RL, ('lam', 'cs'), 0) + BL.lin_cands(RL, ('As', 'lams'), 0)
    dv = lambda cd, c: FL.to_fr(Q.dens(1, cd.slots, lambda s: cd.L(s, c['ctx']))(c['legs']))
    qv = lambda cd, c: FL.to_fr(Q.apply_Q(Q.dens(1, cd.slots, lambda s: cd.L(s, c['ctx'])))(c['legs']))
    import random
    random.seed(17)
    tg = []
    for spcs in (['c', 'lams'], ['lams', 'c'], ['A', 'lam'], ['lam', 'A'], ['A', 'lams'], ['lam', 'As']):
        for _ in range(3):
            legs = [Q.rand_leg(s) for s in spcs]; legs[1]['p'] = [-x for x in legs[0]['p']]
            tg.append(dict(legs=legs, ctx=dict(xi=CP.rspin(CP.SPIN[chir]), b=np.array([random.randint(-2, 2) for _ in range(4)], complex))))

    def named(c):
        xi, b = c['ctx']['xi'], c['ctx']['b']
        s2 = S.S2xi_dens(xi)
        out = {}
        for d in ('V', 'perp', 'lam'):
            out['{S(xi),S(P_b)_%s}' % d] = FL.to_fr(T.bracket(s2, 0, S.translation_dens(b, (d,)))(c['legs']))
        Tx = Q.xi_terms(xi)
        out['{T12(xi),S(P_b)}'] = FL.to_fr(T.bracket(Q.dsum((1, Tx['T1']), (1, Tx['T2'])), 0, S.translation_dens(b))(c['legs']))
        out['{S^V(xi),S(P_b)}'] = FL.to_fr(T.bracket(Q.dsum((1, Tx['T1']), (1, Tx['T3'])), 0, S.translation_dens(b))(c['legs']))
        return out
    return 'routeB', cfg, cands_of_sector, deg2, dv, qv, tg, named, B.log


def route_a(chir):
    import ca_common as A, ca_tree as T, cb_common as B
    Q = B.Q
    RL = BL.Real(A.SB, A.SS, A.MP)
    import route_a as RA
    raw = json.load(open(os.path.join(A.COMPLETE, 'out.nosync', f'A_xiP_{"p" if chir == "+" else "m"}_raw.json')))
    cvec = lambda a: np.array([complex(x[0], x[1]) for x in a])
    cfg = {}
    for sec, recs in raw.items():
        cfg[sec] = []
        for r in recs:
            L = [A.SA.leg(l['t'], cvec(l['pol']), np.array(l['p0'])) for l in r['legs']]
            ob = FL.lvec_scale(RA.v5_to_lvec(r['raw']['0']), (Fr(-1), Fr(0)))
            pole = FL.lvec_scale(RA.v5_to_lvec(r['raw']['1']), (Fr(-1), Fr(0)))
            cfg[sec].append(dict(legsA=L, legs=A.generic_legs(L), ctx=dict(xi=cvec(r['xi']), b=np.array(r['b'], complex)), ob=ob, pole=pole))
    import ca_xiP as XA

    def cands_of_sector(sec):
        t0, t1 = XA.SECT[chir][sec]
        pr = []
        sps = []
        for t in (t0, t1):
            if t in RA.SPRNG:
                sps.append(RA.SPRNG[t][0]); pr.append(tuple(RA.SPRNG[t][1]))
            else:
                sps.append(RA.G5[t][0]); pr.append(RA.chir_proj(t))
        nm = 2 if sps[0] == 'A' else 1
        return [BL.Cand(c.label, c.slots, RA.project(c.L, tuple(pr))) for c in BL.lin_cands(RL, tuple(sps), nm)]
    deg2 = lambda b: BL.lin_cands(RL, ('lam', 'cs'), 0) + BL.lin_cands(RL, ('As', 'lams'), 0)
    dv = lambda cd, c: FL.to_fr(Q.dens(1, cd.slots, lambda s: cd.L(s, c['ctx']))(c['legs']))
    qv = lambda cd, c: FL.to_fr(T.apply_Q(Q.dens(1, cd.slots, lambda s: cd.L(s, c['ctx'])))(c['legs']))
    rp, rpol = XA.CPA.rngs(77)
    tg = []
    for t0, t1 in (('c', 'lamsP'), ('c', 'lamsM'), ('AV', 'lamM'), ('AW', 'lamP'), ('AV', 'lamP'), ('AW', 'lamM')):
        for _ in range(3):
            p = rp(); L = [A.SA.leg(t0, rpol(t0), p), A.SA.leg(t1, rpol(t1), -p)]
            tg.append(dict(legsA=L, legs=A.generic_legs(L), ctx=dict(xi=rpol('lamP' if chir == '+' else 'lamM'), b=rp().astype(complex))))

    def named(c):
        xi, b = c['ctx']['xi'], c['ctx']['b']
        p = np.array([int(round(z.real)) for z in c['legsA'][0]['p0']])
        mm = A.MD.MM(p)
        D = T.mm_deriv(mm, T.mm_S2xi(mm, xi), 0)
        out = {}
        for d in ('V', 'perp', 'lam'):
            out['{S(xi),S(P_b)_%s}' % d] = RA._val(mm, D(T.mm_translation(mm, b, (d,))), c)
        DP = T.mm_deriv(mm, T.mm_translation(mm, b), 1)
        # {T12(xi), S(P_b)} = (-1)^{0}{S(P_b), T12} ... |T12| = 0, |S(P_b)| = 1: {g,f} = (-1)^{|f||g|}{f,g} -> equal
        out['{T12(xi),S(P_b)}'] = RA._val(mm, DP(T.mm_T(mm, xi, 'T1') + T.mm_T(mm, xi, 'T2')), c)
        out['{S^V(xi),S(P_b)}'] = RA._val(mm, DP(T.mm_S2xi(mm, xi, ('V',))), c)
        return out
    return 'routeA', cfg, cands_of_sector, deg2, dv, qv, tg, named, A.log


def residual(obv, nfc, bname):
    # R2 = Ob2 + {F1(xi),S(P_b)} + {S(xi),F1(P_b)}, F1(xi) translation invariant -> {F1(xi),S(P_b)} computed = 0
    R = obv + THP[0] * nfc['{S(xi),S(P_b)_V}'] + THP[1] * nfc['{S(xi),S(P_b)_perp}'] + THP[2] * nfc['{S(xi),S(P_b)_lam}']
    return R


def main(route, chir):
    nm, cfg, cos, deg2, dv, qv, tg, named, log = (route_b if route == 'B' else route_a)(chir)
    blocks = {'all': dict(sectors=list(cfg.keys()))}
    r = CF.analyse(f'{nm} (xi{chir},P_b)', blocks, cfg, cos, deg2, dv, qv, tg, named, residual, THP, log=lambda s: log(s, 'fit_xiP'))
    # WZ2: Q Ob2|2 must vanish ({Ob1(xi),S(P_b)} = 0 on two legs)
    r['G_WZ2'] = dict(QOb2_zero=all(z == '0' for z in r['QOb']), rule='Q Ob2(xi,P_b)|2 = {Ob1(xi)|2, S(P_b)} = 0')
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync', f'{nm}_xiP_{"p" if chir == "+" else "m"}_fit.json')
    json.dump(r, open(out, 'w'), indent=1)
    log(f'{nm} (xi{chir},P_b): G-WZ2 Q Ob2 = 0: {r["G_WZ2"]["QOb2_zero"]}', 'fit_xiP')
    return r


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
