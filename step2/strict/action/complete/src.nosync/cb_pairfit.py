"""Route B: H^{-1}(Q) on the two-leg degree -1, weight 5, translation-invariant functionals bilinear in (xi1, xi2), the
class of Ob^(2)(xi1,xi2)|_2 (raw data of cb_pairs.py), gates G-pole2 and G-WZ2 (corrected signs, finding F1):
   Q Ob2|_2 = Ob1([x,y])|_2 + {Ob1(xi1)|_2, S2(xi2)} + {S2(xi1), Ob1(xi2)|_2} - {Ob0|_2, S^(2)_2}   (odd x, y),
and the residual R2 = Ob2 - F1([x,y]) + {F1(xi1),S(xi2)} + {S(xi1),F1(xi2)} (eq at-eqs, (-1)^{|x|} = -1) with
F1(xi)|_2 = (L - 79/12)(T1 + T2) + theta^xi_V S^V_2(xi) + theta^xi_perp S^perp_2(xi), F1(P_b)|_2 = sum_d theta^P_d
S(P_b)_d (+ exact), Ob1(P_b)|_2 = 0, Ob1(xi)|_2 = (79/24 - L/2)(g1 + g3) + 27/8 g2, Ob0|_2 = -27/8 u int(d*dc,d*A_V)."""
import sys, os, json, random, itertools
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import cb_common as B, cb_s2 as S, cb_tree as T, bilin as BL, fitlib as FL
Q = B.Q
OUT = os.path.join(B.COMPLETE, 'out.nosync')
RL = BL.Real(B.R, B.RB, B.MP)
L32 = {'1': (Fr(0), Fr(0)), 'log2': (Fr(-1), Fr(0)), 'log3': (Fr(1), Fr(0))}

BLOCKS = {
    'mm': {'M': dict(sectors=['AV_Asp', 'Ap_AsV'], deg1=[('AAs', 1)], deg2=[('Acs', 0)], count=4),
           'L': dict(sectors=['l+_ls-', 'l-_ls+'], deg1=[('llS', 1)], deg2=[('lsls', 0)], count=6)},
    'pm': {'V': dict(sectors=['c_cs', 'AV_AsV'], deg1=[('ccs', 1), ('AAs', 1)], deg2=[('Acs', 0)], count=5),
           'P': dict(sectors=['Ap_Asp'], deg1=[('AAs', 1)], deg2=[('Acs', 0)], count=2),
           'L': dict(sectors=['l+_ls+', 'l-_ls-'], deg1=[('llS', 1)], deg2=[('lsls', 0)], count=8)},
    'pp': {'M': dict(sectors=['AV_Asp', 'Ap_AsV'], deg1=[('AAs', 1)], deg2=[('Acs', 0)], count=4),
           'L': dict(sectors=['l+_ls-', 'l-_ls+'], deg1=[('llS', 1)], deg2=[('lsls', 0)], count=6)},
}


def cands_of(kind, nmom, proj=None):
    if kind in ('llS', 'lsls'):
        cs = BL.spinor_cands(RL, kind, nmom) + BL.spinor_cands_b(RL, kind, nmom)
    else:
        cs = BL.boson_cands(RL, kind, nmom)
    if proj is None:
        return cs
    return [BL.Cand(c.label + '|proj' + str(proj), c.slots, project(c.L, proj)) for c in cs]


def project(L, proj):
    """local projection of the slot fields onto component ranges (chirality / V-perp split): proj = (rng0, rng1)."""
    def Lp(s, ctx):
        s2 = []
        for leg, rg in zip(s, proj):
            if rg is None:
                s2.append(leg); continue
            u = np.zeros_like(leg['u']); idx = list(rg); u[idx] = leg['u'][idx]
            s2.append(dict(sp=leg['sp'], u=u, p=leg['p']))
        return L(s2, ctx)
    return Lp


def cvec(a):
    return np.array([complex(x[0], x[1]) for x in a])


def load(pair):
    raw = json.load(open(os.path.join(OUT, f'B_pair_{pair}_raw.json')))
    cfg = {}
    for sec, recs in raw.items():
        L = []
        for r in recs:
            legs = [dict(sp=l['sp'], u=cvec(l['u']), p=[Fr(x) for x in l['p']]) for l in r['legs']]
            ctx = dict(xi1=cvec(r['xi1']), xi2=cvec(r['xi2']))
            ob = FL.lvec_zero(); pole = FL.lvec_zero()
            for key, sgn in (('raw2', -1), ('raw1', 1)):
                d = r[key].get('k', {})
                if '0' in d:
                    ob = FL.lvec_add(ob, FL.clog_parse(d['0']), sgn)
                if '1' in d:
                    pole = FL.lvec_add(pole, FL.clog_parse(d['1']), sgn)
            L.append(dict(legs=legs, ctx=ctx, ob=ob, pole=pole))
        cfg[sec] = L
    return cfg


def dval(fn, legs):
    return FL.to_fr(fn(legs))


def F1dens(xi, f, thV, thP):
    Tx = Q.xi_terms(xi)
    return Q.dsum((f, Tx['T1']), (f, Tx['T2']), (thV, Tx['T1']), (thV, Tx['T3']), (thP, Tx['T2']), (thP, Tx['T4']))


def Ob1dens(xi, c13, c2):
    g = Q.g_terms(xi)
    return Q.dsum((c13, g['g1']), (c13, g['g3']), (c2, g['g2']))


def named_functionals(ctx, pair):
    """tree-level degree -1 functionals at the config (densities)."""
    x1, x2 = ctx['xi1'], ctx['xi2']
    b = S.gam(x1, x2)[:4]
    s2a, s2b = S.S2xi_dens(x1), S.S2xi_dens(x2)
    T1, T2 = Q.xi_terms(x1), Q.xi_terms(x2)
    out = {}
    if pair == 'pm':
        out['S(P_G)_V'] = S.translation_dens(b, ('V',))
        out['S(P_G)_perp'] = S.translation_dens(b, ('perp',))
        out['S(P_G)_lam'] = S.translation_dens(b, ('lam',))
    out['{S^V(xi1),S(xi2)}'] = T.bracket(Q.dsum((1, T1['T1']), (1, T1['T3'])), 0, s2b)
    out['{S^perp(xi1),S(xi2)}'] = T.bracket(Q.dsum((1, T1['T2']), (1, T1['T4'])), 0, s2b)
    out['{S(xi1),S^V(xi2)}'] = T.bracket(s2a, 0, Q.dsum((1, T2['T1']), (1, T2['T3'])))
    out['{S(xi1),S^perp(xi2)}'] = T.bracket(s2a, 0, Q.dsum((1, T2['T2']), (1, T2['T4'])))
    out['{T12(xi1),S(xi2)}+{S(xi1),T12(xi2)}'] = Q.dsum(
        (1, T.bracket(Q.dsum((1, T1['T1']), (1, T1['T2'])), 0, s2b)),
        (1, T.bracket(s2a, 0, Q.dsum((1, T2['T1']), (1, T2['T2'])))))
    return out


def deg0_targets(pair, n=3, seed=7):
    random.seed(seed)
    import cb_pairs as CP
    c1, c2 = CP.CHIR[pair]
    out = []
    for spcs in (['c', 'As'], ['As', 'c'], ['A', 'A'], ['lam', 'lam'], ['c', 'c'], ['A', 'cs'], ['A', 'lams'], ['lam', 'As']):
        for _ in range(n):
            legs = [Q.rand_leg(s) for s in spcs]
            legs[1]['p'] = [-x for x in legs[0]['p']]
            out.append(dict(legs=legs, ctx=dict(xi1=CP.rspin(CP.SPIN[c1]), xi2=CP.rspin(CP.SPIN[c2]))))
    return out


def analyse(pair, log=B.log):
    cfg = load(pair)
    res = {'pair': pair, 'blocks': {}}
    tg = deg0_targets(pair)
    # theta symbols
    th = sp.symbols('thP_V thP_perp thP_lam thx1_V thx1_perp thx2_V thx2_perp')
    QOB = None
    for bname, blk in BLOCKS[pair].items():
        confs = [c for s in blk['sectors'] for c in cfg[s]]
        import cb_pairs as CP
        cands = []
        for sec in blk['sectors']:
            s0, r0, s1, r1 = CP.SECTORS[pair][sec]
            kind = {('c', 'cs'): 'ccs', ('A', 'As'): 'AAs', ('lam', 'lams'): 'llS'}[(s0, s1)]
            nm = [n for (k, n) in blk['deg1'] if k == kind][0]
            cands += cands_of(kind, nm, (None if r0 is None else tuple(r0), None if r1 is None else tuple(r1)))
        rows = [[dval(Q.dens(1, cd.slots, lambda s, cd=cd, ctx=c['ctx']: cd.L(s, ctx)), c['legs']) for c in confs] for cd in cands]
        keep = FL.independent_rows(rows)
        basis = [cands[i] for i in keep]; brows = [rows[i] for i in keep]
        log(f'pair {pair} block {bname}: {len(cands)} candidates, independent {len(basis)} (expected {blk["count"]})', 'pairfitB')
        # closedness
        qrows = []
        for cd in basis:
            qrows.append([FL.gq(dval(Q.apply_Q(Q.dens(1, cd.slots, lambda s, cd=cd, ctx=t['ctx']: cd.L(s, ctx))), t['legs'])) for t in tg])
        Mq = sp.Matrix(qrows).T
        ker = Mq.nullspace()
        # exact image
        c2 = [c for (k, nm) in blk['deg2'] for c in cands_of(k, nm)]
        img = []
        for cd in c2:
            r = [dval(Q.apply_Q(Q.dens(1, cd.slots, lambda s, cd=cd, ctx=c['ctx']: cd.L(s, ctx))), c['legs']) for c in confs]
            if any(z != (0, 0) for z in r):
                x = FL.solve_coords(brows, r)
                assert x is not None, 'Q image outside the degree -1 basis'
                img.append(sp.Matrix(x))
        imM = sp.Matrix.hstack(*img) if img else sp.zeros(len(basis), 0)
        kerM = sp.Matrix.hstack(*ker) if ker else sp.zeros(len(basis), 0)
        Hdim = kerM.rank() - imM.rank()
        # obstruction coordinates
        obc = FL.lsolve_coords(brows, [c['ob'] for c in confs])
        poles = [c['pole'] for c in confs]
        pole_zero = all(all(v == (0, 0) for v in p.values()) for p in poles)
        ob_in_span = all(v is not None for v in obc.values())
        obv = FL.combine(obc) if ob_in_span else None
        # named tree functionals
        nf = {}
        for c in confs:
            for k, f in named_functionals(c['ctx'], pair).items():
                nf.setdefault(k, []).append(dval(f, c['legs']))
        nfc = {}
        for k, r in nf.items():
            x = FL.solve_coords(brows, r)
            nfc[k] = None if x is None else sp.Matrix(x)
        # G-WZ2: Q Ob2 on deg-0 targets vs RHS
        L = sp.log(3) - sp.log(2)
        wz = None
        if False and ob_in_span:
            lhs = sp.zeros(len(tg), 1)
            for i, cd in enumerate(basis):
                lhs += obv[i] * sp.Matrix(qrows[i])
            rr = []
            for t in tg:
                x1, x2 = t['ctx']['xi1'], t['ctx']['xi2']
                parts = {}
                for nm, (c13, c2v) in (('g13', (1.0, 0.0)), ('g2', (0.0, 1.0))):
                    f = Q.dsum((1, T.bracket(S.S2xi_dens(x2), 0, Ob1dens(x1, c13, c2v))),
                               (1, T.bracket(S.S2xi_dens(x1), 0, Ob1dens(x2, c13, c2v))))
                    parts[nm] = FL.gq(dval(f, t['legs']))
                parts['rat'] = sp.Rational(79, 24) * parts['g13'] + sp.Rational(27, 8) * parts['g2']
                parts['L'] = sp.Rational(-1, 2) * parts['g13']
                parts['ob0'] = sp.Rational(-27, 8) * FL.gq(dval(T.bracket(S.s2_dens(x1, x2), 0, Q.F_dsdc_dsA), t['legs']))
                rr.append(parts['rat'] + L * parts['L'] - parts['ob0'])
            wz = sp.simplify(lhs - sp.Matrix(rr))
            wz_ok = all(sp.nsimplify(sp.expand(z)) == 0 for z in wz)
        # residual class conditions
        cond = None
        if ob_in_span and all(v is not None for v in nfc.values()):
            f = L - sp.Rational(79, 12)
            R = obv + f * nfc['{T12(xi1),S(xi2)}+{S(xi1),T12(xi2)}'] \
                + th[3] * nfc['{S^V(xi1),S(xi2)}'] + th[4] * nfc['{S^perp(xi1),S(xi2)}'] \
                + th[5] * nfc['{S(xi1),S^V(xi2)}'] + th[6] * nfc['{S(xi1),S^perp(xi2)}']
            if pair == 'pm':
                R = R - th[0] * nfc['S(P_G)_V'] - th[1] * nfc['S(P_G)_perp'] - th[2] * nfc['S(P_G)_lam']
            # R in Im Q  <=>  exists y: R = imM y
            ys = sp.symbols('y0:%d' % imM.shape[1]) if imM.shape[1] else ()
            eqs = list(R - (imM * sp.Matrix(ys) if ys else sp.zeros(len(basis), 1)))
            sol = sp.solve(eqs, list(ys) + list(th), dict=True)
            cond = dict(equations=[FL.lstr(e) for e in eqs], solution=[{str(k): FL.lstr(v) for k, v in s.items()} for s in sol])
            # also: R closed?
            Rclosed = sp.Matrix.hstack(kerM, R).rank() == kerM.rank() if kerM.shape[1] else all(sp.expand(z) == 0 for z in R)
        ob_class = None
        if ob_in_span:
            # decomposition of Ob2 onto named closed functionals modulo exact
            names = [k for k in nfc if nfc[k] is not None]
            cs = sp.symbols('z0:%d' % len(names)); ys = sp.symbols('w0:%d' % imM.shape[1])
            eqs = list(obv - sum((c * nfc[k] for c, k in zip(cs, names)), sp.zeros(len(basis), 1))
                       - (imM * sp.Matrix(ys) if ys else sp.zeros(len(basis), 1)))
            sol = sp.solve(eqs, list(cs) + list(ys), dict=True)
            ob_class = dict(names=names, solution=[{str(k): FL.lstr(v) for k, v in s.items()} for s in sol])
        if ob_in_span:
            lhs = sp.zeros(len(tg), 1)
            for i, cd in enumerate(basis):
                lhs += obv[i] * sp.Matrix(qrows[i])
            QOB = lhs if QOB is None else QOB + lhs
        res['blocks'][bname] = dict(
            basis=[cd.label for cd in basis], n_candidates=len(cands), basis_dim=len(basis), expected=blk['count'],
            closed_dim=kerM.rank(), exact_dim=imM.rank(), H_minus1_dim=Hdim,
            ob_coords={k: (None if v is None else [FL.lstr(z) for z in v]) for k, v in obc.items()},
            ob_in_span=ob_in_span, pole_sum_zero=pole_zero,
            ob_closed=(None if not ob_in_span else bool(sp.Matrix.hstack(kerM, obv).rank() == kerM.rank()) if kerM.shape[1] else None),
            G_WZ2=(None if wz is None else bool(wz_ok)),
            named_coords={k: (None if v is None else [FL.lstr(z) for z in v]) for k, v in nfc.items()},
            ob_class=ob_class, residual=cond, residual_closed=(None if cond is None else bool(Rclosed)))
        log(f'  block {bname}: closed {kerM.rank()}, exact {imM.rank()}, H^-1 {Hdim}, Ob in span {ob_in_span}, '
            f'poles 0 {pole_zero}, G-WZ2 {res["blocks"][bname]["G_WZ2"]}', 'pairfitB')
    # G-WZ2 on the full functional (all blocks), corrected signs (finding F1)
    L = sp.log(3) - sp.log(2)
    rr = []
    for t in tg:
        x1, x2 = t['ctx']['xi1'], t['ctx']['xi2']
        parts = {}
        for nm, (c13, c2v) in (('g13', (1.0, 0.0)), ('g2', (0.0, 1.0))):
            f = Q.dsum((1, T.bracket(S.S2xi_dens(x2), 0, Ob1dens(x1, c13, c2v))),
                       (1, T.bracket(S.S2xi_dens(x1), 0, Ob1dens(x2, c13, c2v))))
            parts[nm] = FL.gq(dval(f, t['legs']))
        ob0 = sp.Rational(-27, 8) * FL.gq(dval(T.bracket(S.s2_dens(x1, x2), 0, Q.F_dsdc_dsA), t['legs']))
        rr.append(sp.Rational(79, 24) * parts['g13'] + sp.Rational(27, 8) * parts['g2'] - sp.Rational(1, 2) * L * parts['g13'] - ob0)
    wzres = [sp.nsimplify(sp.expand(z)) for z in (QOB - sp.Matrix(rr))]
    res['G_WZ2_full'] = dict(pass_=all(z == 0 for z in wzres), max_rhs=max(abs(complex(sp.N(z))) for z in rr),
                             rule='Q Ob2 = {Ob1(xi1),S(xi2)} + {S(xi1),Ob1(xi2)} - {Ob0,S^(2)} (+ Ob1([x,y]) = 0 on two legs)')
    log(f'pair {pair}: G-WZ2 on the full two-leg functional: {res["G_WZ2_full"]["pass_"]} (|rhs| max {res["G_WZ2_full"]["max_rhs"]:.1f})', 'pairfitB')
    json.dump(res, open(os.path.join(OUT, f'B_pair_{pair}_fit.json'), 'w'), indent=1)
    return res


if __name__ == '__main__':
    for pair in sys.argv[1:] or ['pm', 'pp']:
        r = analyse(pair)
        for bn, b in r['blocks'].items():
            print(bn, json.dumps({k: b[k] for k in ('basis_dim', 'expected', 'closed_dim', 'exact_dim', 'H_minus1_dim', 'ob_in_span', 'pole_sum_zero', 'G_WZ2', 'residual_closed')}))
            print('   ob_class', b['ob_class'])
            print('   residual solution', b['residual'] and b['residual']['solution'])
