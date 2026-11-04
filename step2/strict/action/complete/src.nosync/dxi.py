"""(Delta, xi) on two legs, p <= 1, both routes: Ob^(2)(Delta,xi)|_2 (intrinsic + position parts, defect cycles),
gates G-pole2, G-WZ2 (Q Ob2 = Ob1([Delta,xi]) - {Ob1(Delta),S(xi)} + {S(Delta),Ob1(xi)}, corrected signs), G-F1Delta
(Q F1part(Delta) = Ob1(Delta)|_2), and the residual
   R2 = Ob2 - F1([Delta,xi]) - {F1(Delta),S(xi)} - {S(Delta),F1(xi)},   [Delta, xi] = -xi,
with F1(xi)|_2 = (L - 79/12)(T1 + T2) + th^xi_V S^V_2 + th^xi_perp S^perp_2 and
F1(Delta)|_2 = 22/3 int(A,A*) - 47/3 int(lam*,lam) - 8/3 int(c*,c) + sum_d th^D_d S(Delta)_d + exact (F12);
conditions for [R2] = 0 in H^{-1}(Q) at weight 9/2, p <= 1.   usage: python3 dxi.py A|B +|-"""
import sys, os, json, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import fitlib as FL, bilin as BL, dtree as DT
ROUTE, CH = sys.argv[1], sys.argv[2]
sys.argv = [sys.argv[0], ROUTE, CH]
import eta1 as H      # realisation, xi_terms, g_terms, spec_of, rand_u, legs_at_of, cands machinery
Q, TT, RL, log0 = H.Q, H.TT, H.RL, H.log
log = lambda s: log0('(Delta,xi) ' + s)
OUT = H.OUT
L32 = sp.log(3) - sp.log(2)
if ROUTE == 'B':
    SECM1 = {'+': [('AV', 'ls+'), ('Ap', 'ls-'), ('lam+', 'Asp'), ('lam-', 'AsV')],
             '-': [('AV', 'ls-'), ('Ap', 'ls+'), ('lam-', 'Asp'), ('lam+', 'AsV')]}[CH]
    SECM2 = {'+': [('lam-', 'cs'), ('AsV', 'ls+'), ('Asp', 'ls-')], '-': [('lam+', 'cs'), ('AsV', 'ls-'), ('Asp', 'ls+')]}[CH]
    SEC0 = {'+': [('lam-', 'AV'), ('lam+', 'Ap'), ('c', 'ls+')], '-': [('lam+', 'AV'), ('lam-', 'Ap'), ('c', 'ls-')]}[CH]
    XIT = 'lam+' if CH == '+' else 'lam-'
else:
    SECM1 = {'+': [('AV', 'lamsP'), ('AW', 'lamsM'), ('lamP', 'AsW'), ('lamM', 'AsV')],
             '-': [('AV', 'lamsM'), ('AW', 'lamsP'), ('lamM', 'AsW'), ('lamP', 'AsV')]}[CH]
    SECM2 = {'+': [('lamM', 'cs'), ('AsV', 'lamsP'), ('AsW', 'lamsM')], '-': [('lamP', 'cs'), ('AsV', 'lamsM'), ('AsW', 'lamsP')]}[CH]
    SEC0 = {'+': [('lamM', 'AV'), ('lamP', 'AW'), ('c', 'lamsP')], '-': [('lamP', 'AV'), ('lamM', 'AW'), ('c', 'lamsM')]}[CH]
    XIT = 'lamP' if CH == '+' else 'lamM'


def cands(sector, Dint, Dpos, ckey='xi'):
    t0, t1 = sector
    (s0, P0), (s1, P1) = H.spec_of(t0), H.spec_of(t1)
    slots = (s0, s1); pr = H.projector((P0, P1))
    out = []
    if Dint is not None:
        for c in BL.lin_cands(RL, slots, Dint, vecs_extra=(), ckey=ckey, cspace='S'):
            out.append((c.label, lambda ctx, c=c: DT.PFn(Q.dens(1, slots, lambda s: c.L(pr(s), ctx)), None)))
    if Dpos is not None:
        for c in BL.lin_cands(RL, slots, Dpos, vecs_extra=('e',), ckey=ckey, cspace='S'):
            def mk(ctx, c=c):
                ws = []
                for mu in range(4):
                    e = np.zeros(4, complex); e[mu] = 1
                    cx = dict(ctx); cx['e'] = e
                    ws.append(Q.dens(1, slots, lambda s, cx=cx: c.L(pr(s), cx)))
                return DT.PFn(None, ws)
            out.append(('x.' + c.label, mk))
    return out


def D_of(sector, deg):
    s0 = H.spec_of(sector[0])[0]
    if deg == -1:
        return (1, 2) if s0 == 'A' else (0, 1)
    if deg == -2:
        return (None, 0)
    return (2, 3) if s0 == 'lam' else (2, 3)


def pf(w0=None, wmu=None):
    return DT.PFn(w0, wmu)


def S_Delta(parts=('perp', 'lam', 'V')):
    """S(Delta) = sum_mu x^mu 2[(d_mu A, A*) - (d_mu lam, lam*) - (d_mu c, c*)] + 2(A,A*) - 3(lam,lam*), sector parts."""
    MP_ = H.MP_
    dmu = lambda leg, mu: 1j * complex(leg['p'][mu])
    w0t, wmut = [], [[] for _ in range(4)]
    if 'perp' in parts:
        w0t.append((2, Q.dens(1, ('A', 'As'), lambda s: np.dot(s[0]['u'][4:], s[1]['u'][4:]))))
        for mu in range(4):
            wmut[mu].append((2, Q.dens(1, ('A', 'As'), lambda s, mu=mu: dmu(s[0], mu) * np.dot(s[0]['u'][4:], s[1]['u'][4:]))))
    if 'V' in parts:
        w0t.append((2, Q.dens(1, ('A', 'As'), lambda s: np.dot(s[0]['u'][:4], s[1]['u'][:4]))))
        for mu in range(4):
            wmut[mu].append((2, Q.dens(1, ('A', 'As'), lambda s, mu=mu: dmu(s[0], mu) * np.dot(s[0]['u'][:4], s[1]['u'][:4]))))
            wmut[mu].append((-2, Q.dens(1, ('c', 'cs'), lambda s, mu=mu: dmu(s[0], mu) * s[0]['u'][0] * s[1]['u'][0])))
    if 'lam' in parts:
        w0t.append((-3, Q.dens(1, ('lam', 'lams'), lambda s: s[0]['u'] @ MP_ @ s[1]['u'])))
        for mu in range(4):
            wmut[mu].append((-2, Q.dens(1, ('lam', 'lams'), lambda s, mu=mu: dmu(s[0], mu) * (s[0]['u'] @ MP_ @ s[1]['u']))))
    return pf(Q.dsum(*w0t) if w0t else None, [Q.dsum(*w) for w in wmut])


def F1part_Delta():
    MP_ = H.MP_
    AAs = Q.dens(1, ('A', 'As'), lambda s: np.dot(s[0]['u'], s[1]['u']))
    lsl = Q.dens(1, ('lam', 'lams'), lambda s: s[0]['u'] @ MP_ @ s[1]['u'])      # (lam*,lam) = (lam,lam*)
    csc = Q.dens(1, ('cs', 'c'), lambda s: s[0]['u'][0] * s[1]['u'][0])
    return dict(AAs=AAs, lsl=lsl, csc=csc)


def ob1_Delta_dens():
    MP_ = H.MP_
    MR = [MP_ @ H.R_[M] for M in range(10)]
    pv = lambda l: np.array([complex(x) for x in l['p']])
    def dA(leg):
        u = np.asarray(leg['u'], complex); p = np.zeros(10, complex); p[:4] = pv(leg)
        return 1j * (np.outer(p, u) - np.outer(u, p))
    dAdA = Q.dens(1, ('A', 'A'), lambda s: 0.5 * np.sum(dA(s[0]) * dA(s[1])))
    ldl = Q.dens(1, ('lam', 'lam'), lambda s: sum((s[0]['u'] @ MR[m] @ s[1]['u']) * 1j * pv(s[1])[m] for m in range(4)))
    dcAs = Q.dens(1, ('c', 'As'), lambda s: s[0]['u'][0] * np.dot(1j * pv(s[0]), s[1]['u'][:4]))
    return dict(dAdA=dAdA, ldl=ldl, dcAs=dcAs)


def run_loops(nconf=8, seed=4141):
    fn = os.path.join(OUT, f'{ROUTE}_dxi{"p" if CH == "+" else "m"}_raw.json')
    if os.path.exists(fn):
        return json.load(open(fn))
    rng = np.random.default_rng(seed + (0 if ROUTE == 'A' else 1) + (0 if CH == '+' else 2))
    out = {}
    if ROUTE == 'B':
        E, Y, G, P = H.E, H.Y, H.G, H.P
        tdir = G.tmpl_dir(); tw = G.tmpl_weights(2, 3, 0, 'D'); tq, t3 = G.tmpl_xi()
    for sec in SECM1:
        recs = []
        for ic in range(nconf):
            xi = H.rand_u(XIT, rng)
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            while not any(p):
                p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            u0, u1 = H.rand_u(sec[0], rng), H.rand_u(sec[1], rng)
            la = H.legs_at_of(sec[0], sec[1], u0, u1, p)
            t0 = time.time()
            if ROUTE == 'B':
                def gdir(mu, r):
                    e = np.zeros(4, complex); e[mu] = 1
                    return Y.Gen('Delta', 1, tdir, gvec=e, inj=r)
                al, be, top = P.position_split(la, gdir, Y.Gen('Delta', 1, tw), extra_gens=[Y.Gen('xi', 0, tq, gvec=xi, t3=t3)], gen_first=True)
                def cv(d, m):
                    return FL.clog_parse(d['k'][m].to_json()) if 'k' in d and m in d['k'] else FL.lvec_zero()
                neg = lambda lv: FL.lvec_scale(lv, (Fr(-1), Fr(0)))          # theory = -(raw)
                alpha, pa = neg(cv(al, 0)), neg(cv(al, 1))
                betas = [neg(cv(be.get(mu, {}), 0)) for mu in range(4)]
                topz = all(v.is_zero() for d in top.values() for v in d.values())
            else:
                A, RA = H.A, H.RA
                L = [A.SA.leg(sec[0], u0, np.array([int(x) for x in p])), A.SA.leg(sec[1], u1, -np.array([int(x) for x in p]))]
                d = A.SO.density({'kind': 'Delta'}, L, other={'kind': 'xi', 'xi': xi}, what=('fin', ('pole', 1)))
                alpha = RA.v5_to_lvec(d['fin']['alpha'].get('(Y0,Y1)', A.SA.vzero()))
                pa = RA.v5_to_lvec(d[('pole', 1)]['alpha'].get('(Y0,Y1)', A.SA.vzero()))
                betas = [FL.lvec_scale(RA.v5_to_lvec(d['fin']['beta'][mu].get('(Y0,Y1)', A.SA.vzero())), (Fr(0), Fr(1))) for mu in range(4)]
                topz = True
            js = lambda lv: {k: [str(v[0]), str(v[1])] for k, v in lv.items()}
            recs.append(dict(xi=[[z.real, z.imag] for z in xi], p=[str(x) for x in p], u0=[[z.real, z.imag] for z in u0],
                             u1=[[z.real, z.imag] for z in u1], alpha=js(alpha), pole=js(pa), beta=[js(b) for b in betas],
                             top_zero=topz, sec=round(time.time() - t0, 1)))
            log(f'xi{CH} sector {sec} config {ic}: {recs[-1]["sec"]} s')
        out['%s|%s' % sec] = recs
        json.dump(out, open(fn, 'w'))
    return out


def main():
    raw = run_loops()
    confs = []
    for key, recs in raw.items():
        t0, t1 = key.split('|')
        for r in recs:
            la = H.legs_at_of(t0, t1, H.cvec(r['u0']), H.cvec(r['u1']), [Fr(x) for x in r['p']])
            ob = [H.lv_of(r['alpha'])] + [H.lv_of(b) for b in r['beta']]
            confs.append(dict(la=la, ctx=dict(xi=H.cvec(r['xi'])), ob=ob, pole=H.lv_of(r['pole']), top=r['top_zero']))
    res = dict(route=ROUTE, chirality=CH, pole_sum_zero=all(all(v == (0, 0) for v in c['pole'].values()) for c in confs),
               top_zero=all(c['top'] for c in confs))
    cm1 = [c for s in SECM1 for c in cands(s, *D_of(s, -1))]
    rows = [[z for c in confs for z in DT.canon(mk(c['ctx']), c['la'])] for (_, mk) in cm1]
    keep = FL.independent_rows(rows)
    basis = [cm1[i] for i in keep]; brows = [rows[i] for i in keep]
    log(f'xi{CH}: degree -1 candidates {len(cm1)}, independent {len(basis)}')
    obc = FL.lsolve_coords(brows, [lv for c in confs for lv in c['ob']])
    res['ob_in_span'] = all(v is not None for v in obc.values())
    obv = FL.combine(obc)
    # exact directions
    im, imlab = [], []
    for s in SECM2:
        for (lab, mk) in cands(s, *D_of(s, -2)):
            r = [z for c in confs for z in DT.canon(TT.Q(mk(c['ctx'])), c['la'])]
            if any(z != (0, 0) for z in r):
                x = FL.solve_coords(brows, r)
                assert x is not None
                im.append(sp.Matrix(x)); imlab.append(lab)
    imM = sp.Matrix.hstack(*im) if im else sp.zeros(len(basis), 0)
    # closedness targets (degree 0, p <= 1)
    rng = np.random.default_rng(31)
    tg = []
    for (t0, t1) in SEC0:
        for _ in range(8):
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            tg.append(dict(la=H.legs_at_of(t0, t1, H.rand_u(t0, rng), H.rand_u(t1, rng), p), ctx=dict(xi=H.rand_u(XIT, rng))))
    qrows = [[FL.gq(z) for t in tg for z in DT.canon(TT.Q(mk(t['ctx'])), t['la'])] for (_, mk) in basis]
    kerM = sp.Matrix.hstack(*sp.Matrix(qrows).T.nullspace())
    res.update(basis_dim=len(basis), closed=kerM.rank(), exact=imM.rank(), H_minus1=kerM.rank() - imM.rank())
    # named tree-level functionals on the loop configurations
    named = {}
    Sd = {d: S_Delta((d,)) for d in ('perp', 'lam', 'V')}
    Fp = F1part_Delta()
    for c in confs:
        xi = c['ctx']['xi']
        T = H.xi_terms(xi)
        S2 = Q.dsum(*[(1, T[k]) for k in ('T1', 'T2', 'T3', 'T4')])
        T12 = Q.dsum((1, T['T1']), (1, T['T2']))
        SV, SP = Q.dsum((1, T['T1']), (1, T['T3'])), Q.dsum((1, T['T2']), (1, T['T4']))
        items = {'T12': pf(T12), 'S^V(xi)': pf(SV), 'S^perp(xi)': pf(SP)}
        # {S(xi), F1part(Delta)}: S(xi) even, derivation
        for k, f in Fp.items():
            items['{S(xi),%s}' % k] = DT.PFn(TT.brd(S2, 0, f), None)
        for d, SDd in Sd.items():
            items['{S(xi),S(D)_%s}' % d] = TT.br(S2, 0, SDd)
        # {S(Delta), F1(xi)} = {F1(xi), S(Delta)} (F1(xi) even): derivation by the even F1(xi) pieces on S(Delta)
        SDall = S_Delta()
        items['{T12,S(D)}'] = TT.br(T12, 0, SDall)
        items['{S^V(xi),S(D)}'] = TT.br(SV, 0, SDall)
        items['{S^perp(xi),S(D)}'] = TT.br(SP, 0, SDall)
        for k, fn in items.items():
            named.setdefault(k, []).extend(DT.canon(fn, c['la']))
    nc = {k: FL.solve_coords(brows, v) for k, v in named.items()}
    missing = [k for k, v in nc.items() if v is None]
    res['named_outside_basis'] = missing
    N = {k: sp.Matrix(v) for k, v in nc.items() if v is not None}
    # G-F1Delta: Q F1part(Delta) = Ob1(Delta)|_2 (regression target of T1)
    ob1D = ob1_Delta_dens()
    worst = 0
    rr = np.random.default_rng(3)
    for spcs in (['A', 'A'], ['lam', 'lam'], ['c', 'As'], ['As', 'c']):
        for _ in range(3):
            legs = [Q.rand_leg(s) for s in spcs]; legs[1]['p'] = [-x for x in legs[0]['p']]
            lhs = (22 / 3) * TT.Qd(Fp['AAs'])(legs) - (47 / 3) * TT.Qd(Fp['lsl'])(legs) - (8 / 3) * TT.Qd(Fp['csc'])(legs)
            rhs = (-11 / 3) * ob1D['dAdA'](legs) + (47 / 3) * ob1D['ldl'](legs) - (14 / 3) * ob1D['dcAs'](legs)
            worst = max(worst, abs(lhs - rhs))
    res['G_F1Delta'] = bool(worst < 1e-6)
    thx = sp.symbols('thx_V thx_perp'); thD = sp.symbols('thD_perp thD_lam thD_V')
    f = L32 - sp.Rational(79, 12)
    if not missing and res['ob_in_span']:
        F1xi = f * N['T12'] + thx[0] * N['S^V(xi)'] + thx[1] * N['S^perp(xi)']
        brSF = sp.Rational(22, 3) * N['{S(xi),AAs}'] - sp.Rational(47, 3) * N['{S(xi),lsl}'] - sp.Rational(8, 3) * N['{S(xi),csc}'] \
            + thD[0] * N['{S(xi),S(D)_perp}'] + thD[1] * N['{S(xi),S(D)_lam}'] + thD[2] * N['{S(xi),S(D)_V}']
        brFS = f * N['{T12,S(D)}'] + thx[0] * N['{S^V(xi),S(D)}'] + thx[1] * N['{S^perp(xi),S(D)}']
        # R2 = Ob2 + F1(xi) - {F1(Delta),S(xi)} - {S(Delta),F1(xi)};  {F1(Delta),S(xi)} = {S(xi),F1(Delta)},
        # {S(Delta),F1(xi)} = {F1(xi),S(Delta)}  (graded symmetry, S(xi) and F1(xi) even)
        R = obv + F1xi - brSF - brFS
        ys = sp.symbols('y0:%d' % imM.shape[1])
        sol = FL.lin_solve(list(R - imM * sp.Matrix(ys)), list(ys) + list(thx) + list(thD))
        res['solvable'] = bool(sol)
        if sol:
            res['conditions'] = {str(k): FL.lstr(v) for k, v in sol[0].items() if (k in thx or k in thD) and v != k}
            res['psi_exact'] = {imlab[i]: FL.lstr(sol[0][ys[i]]) for i in range(len(ys)) if sol[0][ys[i]] != 0}
        res['R2_closed_check'] = bool(sp.Matrix.hstack(kerM, R.subs({t: 0 for t in thx + thD})).rank() == kerM.rank())
        # class of Ob2 itself on the named closed directions S^V, S^perp modulo exact
        z = sp.symbols('z0:2')
        s2 = FL.lin_solve(list(obv - z[0] * N['S^V(xi)'] - z[1] * N['S^perp(xi)'] - imM * sp.Matrix(ys)), list(z) + list(ys))
        res['Ob2_class_on_SV_Sperp'] = None if not s2 else {str(k): FL.lstr(s2[0][k]) for k in z}
    log(f'xi{CH}: basis {len(basis)}, closed {res["closed"]}, exact {res["exact"]}, H^-1 {res["H_minus1"]}, Ob in span '
        f'{res["ob_in_span"]}, poles 0 {res["pole_sum_zero"]}, G-F1Delta {res["G_F1Delta"]}, outside {missing}, '
        f'solvable {res.get("solvable")}, conditions {res.get("conditions")}, Ob2 class {res.get("Ob2_class_on_SV_Sperp")}')
    res['ob_coords'] = {k: [FL.lstr(z) for z in v] for k, v in obc.items() if v is not None}
    json.dump(res, open(os.path.join(OUT, f'{ROUTE}_dxi{"p" if CH == "+" else "m"}_fit.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
