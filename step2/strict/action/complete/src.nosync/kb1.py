"""Arity one along K_b on two legs (position degree <= 2), both routes: Ob^(1)(K_b)|_2 by the defect cycles with the
vertex S(K_b) = sum_{r<=s} x^r x^s S_rs + sum_r x^r S_r (momentum injection, second order; canonical data canon2),
v = (x,x) b - 2(b,x) x, div v = -8(b,x), d(v,-) = 4(x,-)^(b,-) (ch02);  x^r x^s part: translations by
c^{rs} = delta_rs b - b_r e_s - b_s e_r (symmetric coefficients f_rs);  x^r part: (L_v A)_V with
(d_n v^m) -> 2 d_rn b_m - 2 b_n d_mr - 2 b_r d_mn, A_perp: -2 b_r, lam: N_r = -4 b_r + sum_n b_n rho_r rho_n.
Gates: G-K-closed (Q S_2(K_b) = 0), G-pole, G-WZ1 (Q Ob1 = -{Ob0, S(K_b)}), G-T (structural: Lemma at-special, the
eps^0 coefficients are T-independent; computed with ordinary Schwinger parameters up to infinity = the T-independent
simplex coefficients), S1-K (H^0, H^{-1} at weight 3, p <= 2), and the solution of Q F1(K_b)|_2 = R1(K_b)|_2 =
Ob1(K_b) - {F1^(0), S(K_b)} with every free direction.   usage: python3 kb1.py A|B"""
import sys, os, json, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import fitlib as FL, bilin as BL, dtree as DT
ROUTE = sys.argv[1]
sys.argv = [sys.argv[0], ROUTE, '+']
import eta1 as H
Q, TT, RL = H.Q, H.TT, H.RL
log = lambda s: H.log('K_b ' + s)
OUT = H.OUT
R_, RB_, MP_ = H.R_, H.RB_, H.MP_
if ROUTE == 'B':
    SEC0 = [('c', 'AsV'), ('AV', 'AV'), ('Ap', 'Ap'), ('lam+', 'lam-')]
    SECM1 = [('c', 'cs'), ('AV', 'AsV'), ('Ap', 'Asp'), ('lam+', 'ls+'), ('lam-', 'ls-')]
    SECM2 = [('AV', 'cs'), ('ls+', 'ls-')]
    TGT1 = [('c', 'AV')]
else:
    SEC0 = [('c', 'AsV'), ('AV', 'AV'), ('AW', 'AW'), ('lamP', 'lamM')]
    SECM1 = [('c', 'cs'), ('AV', 'AsV'), ('AW', 'AsW'), ('lamP', 'lamsP'), ('lamM', 'lamsM')]
    SECM2 = [('AV', 'cs'), ('lamsP', 'lamsM')]
    TGT1 = [('c', 'AV')]
D0 = {('c', 'As'): 0, ('A', 'A'): 1, ('lam', 'lam'): 0, ('c', 'cs'): -1, ('A', 'As'): -1, ('lam', 'lams'): -1,
      ('A', 'cs'): -2, ('lams', 'lams'): -2, ('c', 'A'): 1}


def cands(sector):
    t0, t1 = sector
    (s0, P0), (s1, P1) = H.spec_of(t0), H.spec_of(t1)
    slots = (s0, s1); pr = H.projector((P0, P1))
    d0 = D0[slots]
    out = []
    E4 = np.eye(4, dtype=complex)
    for p in (0, 1, 2):
        D = d0 + p
        if D < 0:
            continue
        ex = ('b',) + (('e',) if p >= 1 else ()) + (('f',) if p == 2 else ())
        for c in BL.vec_cands(RL, slots, D, ex):
            if p == 0:
                out.append((c.label, lambda ctx, c=c: DT.PFn2(Q.dens(1, slots, lambda s: c.L(pr(s), ctx)))))
            elif p == 1:
                def mk(ctx, c=c):
                    ws = []
                    for mu in range(4):
                        cx = dict(ctx); cx['e'] = E4[mu]
                        ws.append(Q.dens(1, slots, lambda s, cx=cx: c.L(pr(s), cx)))
                    return DT.PFn2(None, ws)
                out.append(('x.' + c.label, mk))
            else:
                def mk2(ctx, c=c):
                    wm = {}
                    for (m, n) in DT.PAIRS2:
                        cx1 = dict(ctx); cx1['e'] = E4[m]; cx1['f'] = E4[n]
                        cx2 = dict(ctx); cx2['e'] = E4[n]; cx2['f'] = E4[m]
                        if m == n:
                            wm[(m, n)] = Q.dens(1, slots, lambda s, cx1=cx1: c.L(pr(s), cx1))
                        else:
                            wm[(m, n)] = Q.dens(1, slots, lambda s, cx1=cx1, cx2=cx2: c.L(pr(s), cx1) + c.L(pr(s), cx2))
                    return DT.PFn2(None, None, wm)
                out.append(('xx.' + c.label, mk2))
    return out


def trans(c):
    bd = lambda leg: 1j * np.dot([complex(x) for x in leg['p']], c)
    return Q.dsum((1, Q.dens(1, ('A', 'As'), lambda s: bd(s[0]) * np.dot(s[0]['u'], s[1]['u']))),
                  (1, Q.dens(-1, ('c', 'cs'), lambda s: bd(s[0]) * s[0]['u'][0] * s[1]['u'][0])),
                  (1, Q.dens(-1, ('lam', 'lams'), lambda s: bd(s[0]) * (s[0]['u'] @ MP_ @ s[1]['u']))))


def lin_mats(b, r):
    MV = np.zeros((10, 10), complex)
    for n in range(4):
        for m in range(4):
            MV[n, m] = 2 * (r == n) * b[m] - 2 * b[n] * (m == r) - 2 * b[r] * (m == n)
    for I in range(4, 10):
        MV[I, I] = -2 * b[r]
    N = -4 * b[r] * np.eye(16) + sum(b[n] * (RB_[r] @ R_[n]) for n in range(4))
    return MV, N


def lin_dens(b, r, parts=('V', 'perp', 'lam')):
    MV, N = lin_mats(b, r)
    t = []
    if 'V' in parts:
        t.append((1, Q.dens(1, ('A', 'As'), lambda s: (MV[:4, :4] @ s[0]['u'][:4]) @ s[1]['u'][:4])))
    if 'perp' in parts:
        t.append((1, Q.dens(1, ('A', 'As'), lambda s: (MV[4:, 4:] @ s[0]['u'][4:]) @ s[1]['u'][4:])))
    if 'lam' in parts:
        t.append((-1, Q.dens(1, ('lam', 'lams'), lambda s: (N @ s[0]['u']) @ MP_ @ s[1]['u'])))
    return Q.dsum(*t)


def S_K(b, parts=('V', 'perp', 'lam')):
    """S(K_b) as PFn2; sector parts by species (V: A_V and c, perp: A_perp, lam)."""
    def tr(c):
        bd = lambda leg: 1j * np.dot([complex(x) for x in leg['p']], c)
        t = []
        if 'V' in parts:
            t += [(1, Q.dens(1, ('A', 'As'), lambda s: bd(s[0]) * np.dot(s[0]['u'][:4], s[1]['u'][:4]))),
                  (1, Q.dens(-1, ('c', 'cs'), lambda s: bd(s[0]) * s[0]['u'][0] * s[1]['u'][0]))]
        if 'perp' in parts:
            t += [(1, Q.dens(1, ('A', 'As'), lambda s: bd(s[0]) * np.dot(s[0]['u'][4:], s[1]['u'][4:])))]
        if 'lam' in parts:
            t += [(1, Q.dens(-1, ('lam', 'lams'), lambda s: bd(s[0]) * (s[0]['u'] @ MP_ @ s[1]['u'])))]
        return Q.dsum(*t)
    E4 = np.eye(4)
    wmn = {}
    for (m, n) in DT.PAIRS2:
        c = (m == n) * b - b[m] * E4[n] - b[n] * E4[m]
        wmn[(m, n)] = tr(c if m == n else 2 * c)
    return DT.PFn2(None, [lin_dens(b, r, parts) for r in range(4)], wmn)


# ------------------------------------------------------------------ loops
def letter_values(la, b, nm):
    """pseudo-densities: loop value of the letter (kind, idx) with injection r = -(sum of leg momenta)."""
    raise NotImplementedError


def run_loops(nconf=8, seed=6161):
    fn = os.path.join(OUT, f'{ROUTE}_Kb_raw.json')
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    rng = np.random.default_rng(seed + (0 if ROUTE == 'A' else 1))
    if ROUTE == 'B':
        E, Y, G = H.E, H.Y, H.G
        tdir = G.tmpl_dir()
        import sympy as sy
    for sec in SEC0:
        key = '%s|%s' % sec
        if key in out and len(out[key]) >= nconf:
            continue
        recs = out.get(key, [])
        while len(recs) < nconf:
            b = rng.integers(-2, 3, size=4).astype(complex)
            while not np.any(b):
                b = rng.integers(-2, 3, size=4).astype(complex)
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            while not any(p):
                p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            u0, u1 = H.rand_u(sec[0], rng), H.rand_u(sec[1], rng)
            la = H.legs_at_of(sec[0], sec[1], u0, u1, p)
            t0 = time.time()
            if ROUTE == 'B':
                lt = {}
                for r in range(4):
                    MV, N = lin_mats(b, r)
                    MVs = sy.Matrix(10, 10, [sy.nsimplify(complex(z).real) for z in MV.flatten()])
                    Ns = sy.Matrix(16, 16, [sy.nsimplify(complex(z).real) + sy.I * sy.nsimplify(complex(z).imag) for z in N.flatten()])
                    def LA(s, MVs=MVs):
                        return sum((MVs[P, M] * s[0]['u'][M] * s[1]['u'][P] for P in range(10) for M in range(10) if MVs[P, M] != 0), 0)
                    def Ll(s, Ns=Ns):
                        return ((Ns * sy.Matrix(s[0]['u'])).T * H.B.C.Mp * sy.Matrix(s[1]['u']))[0]
                    H.B.V.TERMS['KA%d' % r] = (sy.Integer(1), ('A', 'As'), LA, ('k', (0, 1)))
                    H.B.V.TERMS['Kl%d' % r] = (sy.Integer(-1), ('lam', 'lams'), Ll, ('k', (0, 1)))
                    lt[r] = E.build_template(['KA%d' % r, 'Kl%d' % r], [], extra_odd=1, name=None)
                E4 = np.eye(4)
                def letter(kind, idx, rr):
                    if kind == 'lin':
                        return Y.Gen('K', 1, lt[idx], inj=rr)
                    m, n = idx
                    c = (m == n) * b - b[m] * E4[n] - b[n] * E4[m]
                    c = c if m == n else 2 * c
                    return Y.Gen('K', 1, tdir, gvec=c / 2, inj=rr)
                def pseudo(kind, idx):
                    def f(legs):
                        rr = [-sum(Fr(l['p'][m]) for l in legs) for m in range(4)]
                        v = Y.evaluate(legs, [letter(kind, idx, rr)], mlist=(0, 1))
                        val = v.get('k', {})
                        f.pole.append(val.get(1, E.CLog()))
                        x = val.get(0, E.CLog())
                        return -complex(x)   # theory = -(raw) for even x; complex for canon2 -> exact via to_fr
                    f.pole = []
                    return f
            else:
                A, RA = H.A, H.RA
                def pseudo(kind, idx):
                    def f(legs):
                        rr = [-sum(Fr(l['p'][m]) for l in legs) for m in range(4)]
                        L = [A.SA.leg(sec[0], legs[0]['u'], np.array([int(x) for x in legs[0]['p']])),
                             A.SA.leg(sec[1], legs[1]['u'], np.array([int(Fr(x)) for x in legs[1]['p']]))]
                        if kind == 'lin':
                            r = idx
                            om = np.array([[2 * (r == m) * b[n] - 2 * b[m] * (n == r) for n in range(4)] for m in range(4)]).real
                            g = A.SA.gen_even('Klin%d' % r, om=om, dv=float(-8 * b[r].real))
                        else:
                            m, n = idx
                            E4 = np.eye(4)
                            c = (m == n) * b - b[m] * E4[n] - b[n] * E4[m]
                            c = c if m == n else 2 * c
                            g = A.SA.gen_even('Kq', a=c.real)
                        res = A.SA.cycle_sum(L, [g], slotnames=['K'], r0=np.array([complex(x) for x in rr]), what=('fin', ('pole', 1)))
                        v0 = A.SO.reduce_colour(res, 2, 'fin', 0).get('(Y0,Y1)', A.SA.vzero())
                        v1 = A.SO.reduce_colour(res, 2, ('pole', 1), 0).get('(Y0,Y1)', A.SA.vzero())
                        f.pole.append(v1)
                        f.logs.append(v0)
                        return complex(float(v0[0][0]), float(v0[0][1]))
                    f.pole, f.logs = [], []
                    return f
            # canonical data of the loop sum: canon2 over pseudo-densities, exact parts from the log-vectors
            res_v = {}
            for comp in ['1', 'log2', 'log3', 'log5', 'log7']:
                res_v[comp] = None
            # evaluate per log component: run canon2 with pseudo-densities that return each component
            vals = {}
            def comp_pseudo(kind, idx):
                base = pseudo(kind, idx)
                cache = {}
                def g(legs):
                    k = tuple(tuple(str(x) for x in l['p']) for l in legs)
                    if k not in cache:
                        if ROUTE == 'B':
                            rr = [-sum(Fr(l['p'][m]) for l in legs) for m in range(4)]
                            v = H.Y.evaluate(legs, [letter(kind, idx, rr)], mlist=(0, 1)).get('k', {})
                            cache[k] = (FL.lvec_scale(FL.clog_parse(v[0].to_json()) if 0 in v else FL.lvec_zero(), (Fr(-1), Fr(0))),
                                        FL.lvec_scale(FL.clog_parse(v[1].to_json()) if 1 in v else FL.lvec_zero(), (Fr(-1), Fr(0))))
                        else:
                            base(legs)
                            cache[k] = (RA.v5_to_lvec(base.logs[-1]), RA.v5_to_lvec(base.pole[-1]))
                    return cache[k]
                return g
            lins = {r: comp_pseudo('lin', r) for r in range(4)}
            quads = {mn: comp_pseudo('quad', mn) for mn in DT.PAIRS2}
            canon_all = {}
            for which in (0, 1):
                for comp in FL.LOGK:
                    for part in (0, 1):
                        sel = lambda g, which=which, comp=comp, part=part: (lambda legs: complex(float(g(legs)[which][comp][0]), float(g(legs)[which][comp][1])) if part == 0 else 0)
                        pass
            # exact canonical data per log component (re, im) from the cached log-vectors
            def comp_fn(g, which, comp):
                return lambda legs: tuple(g(legs)[which][comp])
            data = {}
            for which in (0, 1):
                for comp in FL.LOGK:
                    fn2 = DT.PFn2(None, [comp_fn(lins[r], which, comp) for r in range(4)], {mn: comp_fn(quads[mn], which, comp) for mn in DT.PAIRS2})
                    data[(which, comp)] = DT.canon2(fn2, la)
            canon_ob = [{comp: data[(0, comp)][i] for comp in FL.LOGK} for i in range(15)]
            canon_pole = [{comp: data[(1, comp)][i] for comp in FL.LOGK} for i in range(15)]
            js = lambda lv: {k: [str(v[0]), str(v[1])] for k, v in lv.items()}
            recs.append(dict(b=[int(x.real) for x in b], p=[str(x) for x in p], u0=[[z.real, z.imag] for z in u0],
                             u1=[[z.real, z.imag] for z in u1], ob=[js(x) for x in canon_ob], pole=[js(x) for x in canon_pole],
                             sec=round(time.time() - t0, 1)))
            log(f'sector {sec} config {len(recs) - 1}: {recs[-1]["sec"]} s')
            out[key] = recs
            json.dump(out, open(fn, 'w'))
    return out


def main():
    raw = run_loops()
    confs = []
    for key, recs in raw.items():
        t0, t1 = key.split('|')
        for r in recs:
            la = H.legs_at_of(t0, t1, H.cvec(r['u0']), H.cvec(r['u1']), [Fr(x) for x in r['p']])
            confs.append(dict(la=la, ctx=dict(b=np.array(r['b'], complex)), ob=[H.lv_of(x) for x in r['ob']],
                              pole=[H.lv_of(x) for x in r['pole']]))
    res = dict(route=ROUTE, pole_sum_zero=all(all(v == (0, 0) for x in c['pole'] for v in x.values()) for c in confs))
    t0 = time.time()
    c0 = [c for s in SEC0 for c in cands(s)]
    rows = [[z for c in confs for z in DT.canon2(mk(c['ctx']), c['la'])] for (_, mk) in c0]
    keep = FL.independent_rows(rows)
    basis = [c0[i] for i in keep]; brows = [rows[i] for i in keep]
    log(f'degree-0 candidates {len(c0)}, independent {len(basis)} ({time.time() - t0:.0f} s)')
    obc = FL.lsolve_coords(brows, [lv for c in confs for lv in c['ob']])
    res['ob_in_span'] = all(v is not None for v in obc.values())
    obv = FL.combine(obc) if res['ob_in_span'] else None
    rng = np.random.default_rng(71)
    tg = []
    for (t0_, t1_) in TGT1:
        for _ in range(12):
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            tg.append(dict(la=H.legs_at_of(t0_, t1_, H.rand_u(t0_, rng), H.rand_u(t1_, rng), p), ctx=dict(b=rng.integers(-2, 3, size=4).astype(complex))))
    Qd = lambda fn: DT.tree_apply2(TT.Qd, fn)
    qrows = [[FL.gq(z) for t in tg for z in DT.canon2(Qd(mk(t['ctx'])), t['la'])] for (_, mk) in basis]
    ker = sp.Matrix(qrows).T.nullspace()
    kerM = sp.Matrix.hstack(*ker) if ker else sp.zeros(len(basis), 0)
    # degree -1 basis on random configurations, Q images in degree-0 coordinates
    rngm = np.random.default_rng(73)
    cfm = []
    for (t0_, t1_) in SECM1:
        for _ in range(8):
            p = [Fr(int(x)) for x in rngm.integers(-3, 4, size=4)]
            cfm.append(dict(la=H.legs_at_of(t0_, t1_, H.rand_u(t0_, rngm), H.rand_u(t1_, rngm), p), ctx=dict(b=rngm.integers(-2, 3, size=4).astype(complex))))
    named_m1 = {'S(K)_perp': lambda ctx: S_K(ctx['b'], ('perp',)), 'S(K)_lam': lambda ctx: S_K(ctx['b'], ('lam',)),
                'S(K)_V': lambda ctx: S_K(ctx['b'], ('V',))}
    cm1 = [(k, f) for k, f in named_m1.items()] + [c for s in SECM1 for c in cands(s)]
    rws = [[z for c in cfm for z in DT.canon2(mk(c['ctx']), c['la'])] for (_, mk) in cm1]
    keepm = FL.independent_rows(rws)
    basm = [cm1[i] for i in keepm]
    log(f'degree -1: {len(cm1)} candidates, independent {len(basm)} ({time.time() - t0:.0f} s)')
    im, imlab = [], []
    for (lab, mk) in basm:
        r = [z for c in confs for z in DT.canon2(Qd(mk(c['ctx'])), c['la'])]
        x = FL.solve_coords(brows, r)
        assert x is not None, 'Q image outside the degree-0 basis'
        im.append(sp.Matrix(x)); imlab.append(lab)
    imM = sp.Matrix.hstack(*im)
    # exact directions at degree -1 (Q of degree -2)
    im2 = 0
    for s in SECM2:
        for (lab, mk) in cands(s):
            r = [z for c in cfm for z in DT.canon2(Qd(mk(c['ctx'])), c['la'])]
            if any(z != (0, 0) for z in r):
                im2 += 1
    res.update(basis0=len(basis), closed0=kerM.rank(), exact0=imM.rank(), H0=kerM.rank() - imM.rank(),
               deg_m1_basis=len(basm), closed_m1=len(imM.nullspace()))
    # gates
    QSK = [z for c in confs for z in DT.canon2(Qd(S_K(c['ctx']['b'])), c['la'])]
    res['G_K_closed'] = all(z == (0, 0) for z in QSK)
    F10 = Q.F_dsA_dsA
    brF = [z for c in confs for z in DT.canon2(DT.tree_apply2(lambda d: TT.brd(F10, 0, d), S_K(c['ctx']['b'])), c['la'])]
    xF = FL.solve_coords(brows, brF)
    if res['ob_in_span'] and xF is not None:
        R1 = obv - sp.Rational(27, 16) * sp.Matrix(xF)
        res['R1_closed'] = bool(sp.Matrix.hstack(kerM, R1).rank() == kerM.rank())
        lhs = sp.zeros(len(qrows[0]), 1)
        for i in range(len(basis)):
            lhs += obv[i] * sp.Matrix(qrows[i])
        rhs = sp.Matrix([FL.gq(z) for t in tg for z in DT.canon2(DT.tree_apply2(lambda d: TT.brd(Q.F_dsdc_dsA, 1, d), S_K(t['ctx']['b'])), t['la'])])
        res['G_WZ1'] = all(sp.nsimplify(sp.expand(z)) == 0 for z in (lhs - sp.Rational(27, 8) * rhs))
        res['G_WZ1_printed_sign'] = all(sp.nsimplify(sp.expand(z)) == 0 for z in (lhs + sp.Rational(27, 8) * rhs))
        ys = sp.symbols('f0:%d' % len(im))
        sol = FL.lin_solve(list(imM * sp.Matrix(ys) - R1), list(ys))
        res['solvable'] = bool(sol)
        if sol:
            s0 = sol[0]; free = [y for y in ys if s0[y] == y]
            res['F1_particular'] = {imlab[i]: FL.lstr(s0[ys[i]].subs({f: 0 for f in free})) for i in range(len(ys)) if s0[ys[i]].subs({f: 0 for f in free}) != 0}
            res['F1_free_directions'] = [{imlab[i]: FL.lstr(sp.diff(s0[ys[i]], f)) for i in range(len(ys)) if sp.diff(s0[ys[i]], f) != 0} for f in free]
    log(f'pole 0 {res["pole_sum_zero"]}, Ob in span {res["ob_in_span"]}, G-K-closed {res["G_K_closed"]}, H^0 {res["H0"]}, '
        f'closed degree -1 {res["closed_m1"]}, G-WZ1 {res.get("G_WZ1")} (printed {res.get("G_WZ1_printed_sign")}), '
        f'R1 closed {res.get("R1_closed")}, solvable {res.get("solvable")}')
    res['ob_coords'] = {k: [FL.lstr(z) for z in v] for k, v in obc.items() if v is not None}
    json.dump(res, open(os.path.join(OUT, f'{ROUTE}_Kb_fit.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
