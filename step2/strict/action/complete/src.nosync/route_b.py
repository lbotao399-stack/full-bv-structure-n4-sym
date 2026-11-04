"""route B adapter for pairfit_core (engine B code: cycles eb_cycle, Q eb_qb.apply_Q, bracket cb_tree)."""
import json, os
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import cb_common as B, cb_s2 as S, cb_tree as T, bilin as BL, fitlib as FL, cb_pairs as CP
Q = B.Q
name = 'routeB'
RL = BL.Real(B.R, B.RB, B.MP)
OUT = os.path.join(B.COMPLETE, 'out.nosync')
KIND = {('c', 'cs'): ('ccs', 1), ('A', 'As'): ('AAs', 1), ('lam', 'lams'): ('llS', 1)}


def cands(kind, nmom, proj):
    if kind in ('llS', 'lsls'):
        cs = BL.spinor_cands(RL, kind, nmom) + BL.spinor_cands_b(RL, kind, nmom)
    else:
        cs = BL.boson_cands(RL, kind, nmom)
    if proj is None:
        return cs
    return [BL.Cand(c.label + '|proj', c.slots, project(c.L, proj)) for c in cs]


def project(L, proj):
    def Lp(s, ctx):
        s2 = []
        for leg, rg in zip(s, proj):
            if rg is None:
                s2.append(leg); continue
            u = np.zeros_like(leg['u']); idx = list(rg); u[idx] = leg['u'][idx]
            s2.append(dict(sp=leg['sp'], u=u, p=leg['p']))
        return L(s2, ctx)
    return Lp


def sector_cands(pair, sec):
    s0, r0, s1, r1 = CP.SECTORS[pair][sec]
    k, nm = KIND[(s0, s1)]
    return cands(k, nm, (None if r0 is None else tuple(r0), None if r1 is None else tuple(r1)))


def cvec(a):
    return np.array([complex(x[0], x[1]) for x in a])


def load(pair):
    raw = json.load(open(os.path.join(OUT, f'B_pair_{pair}_raw.json')))
    cfg = {}
    for sec, recs in raw.items():
        Ls = []
        for r in recs:
            legs = [dict(sp=l['sp'], u=cvec(l['u']), p=[Fr(x) for x in l['p']]) for l in r['legs']]
            ob, pole = FL.lvec_zero(), FL.lvec_zero()
            for key, sgn in (('raw2', -1), ('raw1', 1)):
                d = r[key].get('k', {})
                if '0' in d:
                    ob = FL.lvec_add(ob, FL.clog_parse(d['0']), sgn)
                if '1' in d:
                    pole = FL.lvec_add(pole, FL.clog_parse(d['1']), sgn)
            Ls.append(dict(legs=legs, ctx=dict(xi1=cvec(r['xi1']), xi2=cvec(r['xi2'])), ob=ob, pole=pole))
        cfg[sec] = Ls
    return cfg


def dens_of(cd, ctx):
    return Q.dens(1, cd.slots, lambda s: cd.L(s, ctx))


def dens_value(cd, c):
    return FL.to_fr(dens_of(cd, c['ctx'])(c['legs']))


def Q_value(cd, c):
    return FL.to_fr(Q.apply_Q(dens_of(cd, c['ctx']))(c['legs']))


def named(c, pair):
    x1, x2 = c['ctx']['xi1'], c['ctx']['xi2']
    b = S.gam(x1, x2)[:4]
    s2a, s2b = S.S2xi_dens(x1), S.S2xi_dens(x2)
    T1, T2 = Q.xi_terms(x1), Q.xi_terms(x2)
    out = {}
    if pair == 'pm':
        for d in ('V', 'perp', 'lam'):
            out['S(P_G)_' + d] = S.translation_dens(b, (d,))
    out['{S^V(xi1),S(xi2)}'] = T.bracket(Q.dsum((1, T1['T1']), (1, T1['T3'])), 0, s2b)
    out['{S^perp(xi1),S(xi2)}'] = T.bracket(Q.dsum((1, T1['T2']), (1, T1['T4'])), 0, s2b)
    out['{S(xi1),S^V(xi2)}'] = T.bracket(s2a, 0, Q.dsum((1, T2['T1']), (1, T2['T3'])))
    out['{S(xi1),S^perp(xi2)}'] = T.bracket(s2a, 0, Q.dsum((1, T2['T2']), (1, T2['T4'])))
    out['{T12(xi1),S(xi2)}+{S(xi1),T12(xi2)}'] = Q.dsum(
        (1, T.bracket(Q.dsum((1, T1['T1']), (1, T1['T2'])), 0, s2b)),
        (1, T.bracket(s2a, 0, Q.dsum((1, T2['T1']), (1, T2['T2'])))))
    return {k: FL.to_fr(f(c['legs'])) for k, f in out.items()}


def deg0_targets(pair, n=3, seed=7):
    import random
    random.seed(seed)
    c1, c2 = CP.CHIR[pair]
    out = []
    for spcs in (['c', 'As'], ['As', 'c'], ['A', 'A'], ['lam', 'lam'], ['c', 'c'], ['A', 'cs'], ['A', 'lams'], ['lam', 'As']):
        for _ in range(n):
            legs = [Q.rand_leg(s) for s in spcs]
            legs[1]['p'] = [-x for x in legs[0]['p']]
            out.append(dict(legs=legs, ctx=dict(xi1=CP.rspin(CP.SPIN[c1]), xi2=CP.rspin(CP.SPIN[c2]))))
    return out


def Ob1dens(xi, c13, c2):
    g = Q.g_terms(xi)
    return Q.dsum((c13, g['g1']), (c13, g['g3']), (c2, g['g2']))


def wz_rhs(t):
    x1, x2 = t['ctx']['xi1'], t['ctx']['xi2']
    parts = {}
    for nm, (c13, c2v) in (('g13', (1.0, 0.0)), ('g2', (0.0, 1.0))):
        f = Q.dsum((1, T.bracket(S.S2xi_dens(x2), 0, Ob1dens(x1, c13, c2v))),
                   (1, T.bracket(S.S2xi_dens(x1), 0, Ob1dens(x2, c13, c2v))))
        parts[nm] = FL.gq(FL.to_fr(f(t['legs'])))
    ob0 = sp.Rational(-27, 8) * FL.gq(FL.to_fr(T.bracket(S.s2_dens(x1, x2), 0, Q.F_dsdc_dsA)(t['legs'])))
    L = sp.log(3) - sp.log(2)
    return sp.Rational(79, 24) * parts['g13'] + sp.Rational(27, 8) * parts['g2'] - sp.Rational(1, 2) * L * parts['g13'] - ob0
