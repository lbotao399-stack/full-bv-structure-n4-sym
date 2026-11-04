"""route A adapter for pairfit_core (engine A code: cycles sa_core/GSimplex, mode-model Q and bracket sa_mode via
ca_tree.mm_deriv; candidate densities in engine A's realisation, Q on densities with qmat_A, gate GA-Q)."""
import json, os
import numpy as np
import sympy as sp
import ca_common as A, ca_tree as T, ca_pairs as CPA, bilin as BL, fitlib as FL, cb_common as B
from fractions import Fraction as Fr
Q = B.Q
name = 'routeA'
RL = BL.Real(A.SB, A.SS, A.MP)
OUT = os.path.join(A.COMPLETE, 'out.nosync')
SPRNG = {'AV': ('A', range(4)), 'AW': ('A', range(4, 10)), 'AsV': ('As', range(4)), 'AsW': ('As', range(4, 10)),
         'c': ('c', None), 'cs': ('cs', None)}
G5 = {'lamP': ('lam', 'P'), 'lamM': ('lam', 'M'), 'lamsP': ('lams', 'P'), 'lamsM': ('lams', 'M')}
KIND = {('c', 'cs'): ('ccs', 1), ('A', 'As'): ('AAs', 1), ('lam', 'lams'): ('llS', 1)}


def chir_proj(t):
    """projector onto the chirality of letter type t (engine A: lamP = (1 + G5S)/2 etc.)"""
    E = A.SA.E
    if t in ('lamP', 'lamM'):
        G = np.asarray(E.G5S, complex)
    else:
        G = np.asarray(E.G5SS, complex)
    s = 1 if t.endswith('P') else -1
    return 0.5 * (np.eye(16) + s * G)


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
        for leg, P in zip(s, proj):
            if P is None:
                s2.append(leg); continue
            if isinstance(P, np.ndarray):
                u = P @ leg['u']
            else:
                u = np.zeros_like(leg['u']); idx = list(P); u[idx] = leg['u'][idx]
            s2.append(dict(sp=leg['sp'], u=u, p=leg['p']))
        return L(s2, ctx)
    return Lp


def sector_cands(pair, sec):
    t0, t1 = CPA.SECTORS[pair][sec]
    pr = []
    sps = []
    for t in (t0, t1):
        if t in SPRNG:
            sps.append(SPRNG[t][0]); pr.append(None if SPRNG[t][1] is None else tuple(SPRNG[t][1]))
        else:
            sps.append(G5[t][0]); pr.append(chir_proj(t))
    k, nm = KIND[tuple(sps)]
    return cands(k, nm, tuple(pr))


def cvec(a):
    return np.array([complex(x[0], x[1]) for x in a])


def v5_to_lvec(v):
    keys = FL.LOGK
    return {k: (Fr(v[i][0]), Fr(v[i][1])) for i, k in enumerate(keys)}


def load(pair):
    raw = json.load(open(os.path.join(OUT, f'A_pair_{pair}_raw.json')))
    cfg = {}
    for sec, recs in raw.items():
        Ls = []
        for r in recs:
            L = [A.SA.leg(l['t'], cvec(l['pol']), np.array(l['p0'])) for l in r['legs']]
            ob = FL.lvec_add(FL.lvec_scale(v5_to_lvec(r['raw2']['0']), (Fr(-1), Fr(0))), v5_to_lvec(r['raw1']['0']))
            pole = FL.lvec_add(FL.lvec_scale(v5_to_lvec(r['raw2']['1']), (Fr(-1), Fr(0))), v5_to_lvec(r['raw1']['1']))
            Ls.append(dict(legsA=L, legs=A.generic_legs(L), ctx=dict(xi1=cvec(r['xi1']), xi2=cvec(r['xi2'])), ob=ob, pole=pole))
        cfg[sec] = Ls
    return cfg


def dens_of(cd, ctx):
    return Q.dens(1, cd.slots, lambda s: cd.L(s, ctx))


def dens_value(cd, c):
    return FL.to_fr(dens_of(cd, c['ctx'])(c['legs']))


def Q_value(cd, c):
    return FL.to_fr(T.apply_Q(dens_of(cd, c['ctx']))(c['legs']))


def _val(mm, F, c):
    v = mm.value(F, A.mm_legs(c['legsA']))
    return (Fr(v.re), Fr(v.im))


def named(c, pair):
    x1, x2 = c['ctx']['xi1'], c['ctx']['xi2']
    p = np.array([int(round(z.real)) for z in c['legsA'][0]['p0']])
    mm = A.MD.MM(p)
    S2a, S2b = T.mm_S2xi(mm, x1), T.mm_S2xi(mm, x2)
    D1, D2 = T.mm_deriv(mm, S2a, 0), T.mm_deriv(mm, S2b, 0)
    out = {}
    if pair == 'pm':
        b = T.gam(x1, x2)[:4]
        for d in ('V', 'perp', 'lam'):
            out['S(P_G)_' + d] = _val(mm, T.mm_translation(mm, b, (d,)), c)
    SV1, SP1 = T.mm_S2xi(mm, x1, ('V',)), T.mm_S2xi(mm, x1, ('perp',))
    SV2, SP2 = T.mm_S2xi(mm, x2, ('V',)), T.mm_S2xi(mm, x2, ('perp',))
    # {Y, S(xi2)} = {S(xi2), Y} (both even)
    out['{S^V(xi1),S(xi2)}'] = _val(mm, D2(SV1), c)
    out['{S^perp(xi1),S(xi2)}'] = _val(mm, D2(SP1), c)
    out['{S(xi1),S^V(xi2)}'] = _val(mm, D1(SV2), c)
    out['{S(xi1),S^perp(xi2)}'] = _val(mm, D1(SP2), c)
    T12a = T.mm_T(mm, x1, 'T1') + T.mm_T(mm, x1, 'T2')
    T12b = T.mm_T(mm, x2, 'T1') + T.mm_T(mm, x2, 'T2')
    out['{T12(xi1),S(xi2)}+{S(xi1),T12(xi2)}'] = _val(mm, D2(T12a) + D1(T12b), c)
    return out


def wz_targets(pair):
    return deg0_targets(pair, n=1, seed=19)


def deg0_targets(pair, n=3, seed=9):
    rp, rpol = CPA.rngs(seed)
    c1, c2 = CPA.CHIR[pair]
    out = []
    for t0, t1 in (('c', 'AsV'), ('AV', 'AV'), ('AW', 'AW'), ('lamP', 'lamM'), ('AV', 'cs'), ('lamP', 'lamsP'), ('lamM', 'lamsM')):
        for _ in range(n):
            p = rp()
            L = [A.SA.leg(t0, rpol(t0), p), A.SA.leg(t1, rpol(t1), -p)]
            out.append(dict(legsA=L, legs=A.generic_legs(L), ctx=dict(xi1=rpol(c1), xi2=rpol(c2))))
    return out


def wz_rhs(t):
    x1, x2 = t['ctx']['xi1'], t['ctx']['xi2']
    p = np.array([int(round(z.real)) for z in t['legsA'][0]['p0']])
    mm = A.MD.MM(p)
    D1, D2 = T.mm_deriv(mm, T.mm_S2xi(mm, x1), 0), T.mm_deriv(mm, T.mm_S2xi(mm, x2), 0)
    g13a = A.MD.F_g(mm, x1, 1) + A.MD.F_g(mm, x1, 3); g13b = A.MD.F_g(mm, x2, 1) + A.MD.F_g(mm, x2, 3)
    g2a, g2b = A.MD.F_g(mm, x1, 2), A.MD.F_g(mm, x2, 2)
    # {Ob1(xi1), S(xi2)} = {S(xi2), Ob1(xi1)} (S even, Ob1 odd)
    v13 = FL.gq(_val(mm, D2(g13a) + D1(g13b), t)); v2 = FL.gq(_val(mm, D2(g2a) + D1(g2b), t))
    DS2 = T.mm_deriv(mm, T.mm_S2(mm, x1, x2), 0)
    vob0 = sp.Rational(-27, 8) * FL.gq(_val(mm, DS2(A.MD.F_named(mm, 'ddcdA')), t))
    L = sp.log(3) - sp.log(2)
    return sp.Rational(79, 24) * v13 + sp.Rational(27, 8) * v2 - sp.Rational(1, 2) * L * v13 - vob0
