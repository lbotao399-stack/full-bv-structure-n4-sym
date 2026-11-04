"""Arity one along eta in Sigma^* on two legs (p <= 1), both routes: Ob^(1)(eta)|_2 (intrinsic and position parts, defect
cycles), gates G-eta-closed (Q S_2(eta) = 0), G-pos (position part = Ob^(1)(rho(e_mu) eta)|_2, compared with the
computed Ob^(1)(xi)|_2 = (79/24 - L/2)(g1 + g3) + 27/8 g2), G-WZ1 (Q Ob1(eta) = -{Ob0, S(eta)}, corrected sign),
S1-eta (H^0(Q), H^{-1}(Q) at weight 7/2, p <= 1), and the solution of Q F1(eta)|_2 = R1(eta)|_2 = Ob1 - {F1^(0), S(eta)}
with F1^(0)|_2 = 27/16 u int(d*A_V, d*A_V) and every free direction explicit.
usage: python3 eta1.py A|B +|-"""
import sys, os, json, random, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import fitlib as FL, bilin as BL, dtree as DT
ROUTE, CH = sys.argv[1], sys.argv[2]
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')
L32 = sp.log(3) - sp.log(2)
import cb_common as B
Q = B.Q
if ROUTE == 'B':
    E, Y, G, P = B.E, B.Y, B.G, B.P
    R_, RB_, MP_ = B.R, B.RB, B.MP
    TT = DT.Tree(lambda k: Q.qmat([complex(x) for x in k]), DT.k0_from_pairing(B.MP))
    log = lambda s: B.log(s, 'eta1B')
    LAMP, LAMM, LSP, LSM = B.LAMP, B.LAMM, B.LAMS_P, B.LAMS_M
    def P_sig(nm):
        return {'lam+': ('lam', LAMP), 'lam-': ('lam', LAMM), 'ls+': ('lams', LSP), 'ls-': ('lams', LSM)}[nm]
else:
    import ca_common as A
    R_, RB_, MP_ = A.SB, A.SS, A.MP
    TT = DT.Tree(A.qmat_A, DT.k0_from_pairing(A.MP))
    log = lambda s: A.log(s, 'eta1A')
    import route_a as RA
RL = BL.Real(R_, RB_, MP_)


def xi_terms(xi):
    """T0..T4 densities in the route's realisation (eq at-S1list), generic legs."""
    MR = [MP_ @ R_[M] for M in range(10)]
    def rhoF(F):
        v = np.zeros(16, complex)
        for M in range(10):
            for N in range(10):
                if F[M, N] != 0:
                    v += 0.5 * F[M, N] * (RB_[M] @ (R_[N] @ xi))
        return v
    def dA(leg, rng):
        u = np.array([leg['u'][M] if M in rng else 0 for M in range(10)], complex)
        p = np.zeros(10, complex); p[:4] = [complex(x) for x in leg['p']]
        return 1j * (np.outer(p, u) - np.outer(u, p))
    T1 = Q.dens(-1, ('lam', 'As'), lambda s: sum((xi @ MR[M] @ s[0]['u']) * s[1]['u'][M] for M in range(4)))
    T2 = Q.dens(-1, ('lam', 'As'), lambda s: sum((xi @ MR[M] @ s[0]['u']) * s[1]['u'][M] for M in range(4, 10)))
    T3 = Q.dens(0.5, ('A', 'lams'), lambda s: rhoF(dA(s[0], range(4))) @ MP_ @ s[1]['u'])
    T4 = Q.dens(0.5, ('A', 'lams'), lambda s: rhoF(dA(s[0], range(4, 10))) @ MP_ @ s[1]['u'])
    return dict(T1=T1, T2=T2, T3=T3, T4=T4)


def g_terms(xi):
    MR = [MP_ @ R_[M] for M in range(10)]
    pv = lambda l: np.array([complex(x) for x in l['p']])
    gam = lambda l, M: xi @ MR[M] @ l['u']
    def g1(s):
        l, a = s; p = pv(a); return sum(gam(l, n) * ((p @ p) * a['u'][n] - p[n] * (p @ a['u'][:4])) for n in range(4))
    def g2(s):
        l, a = s; p = pv(a); return sum(gam(l, n) * p[n] * (p @ a['u'][:4]) for n in range(4))
    def g3(s):
        l, a = s; p = pv(a); return sum(gam(l, M) * (p @ p) * a['u'][M] for M in range(4, 10))
    def g4(s):
        c, ls = s; p = pv(c); return (p @ p) * c['u'][0] * (xi @ MP_ @ ls['u'])
    return dict(g1=Q.dens(1, ('lam', 'A'), g1), g2=Q.dens(1, ('lam', 'A'), g2), g3=Q.dens(1, ('lam', 'A'), g3),
                g4=Q.dens(1, ('c', 'lams'), g4))


def S_eta(eta, parts=('V', 'perp')):
    """S_2(eta) = sum_mu x^mu S_2(rho(e_mu) eta) + int s0, s0 = -(rho(A_perp) eta, lam*); sector parts S^V, S^perp."""
    wmu = []
    for mu in range(4):
        T = xi_terms(RB_[mu] @ eta)
        terms = []
        if 'V' in parts:
            terms += [(1, T['T1']), (1, T['T3'])]
        if 'perp' in parts:
            terms += [(1, T['T2']), (1, T['T4'])]
        wmu.append(Q.dsum(*terms))
    w0 = None
    if 'perp' in parts:
        w0 = Q.dens(-1, ('A', 'lams'), lambda s: sum(s[0]['u'][I] * ((RB_[I] @ eta) @ MP_ @ s[1]['u']) for I in range(4, 10)))
    return DT.PFn(w0, wmu)


# ------------------------------------------------------------------ sectors
if ROUTE == 'B':
    SEC0 = {'+': [('AV', 'lam+'), ('Ap', 'lam-'), ('c', 'ls-')], '-': [('AV', 'lam-'), ('Ap', 'lam+'), ('c', 'ls+')]}[CH]
    SECM1 = {'+': [('AV', 'ls-'), ('Ap', 'ls+'), ('lam+', 'AsV'), ('lam-', 'Asp')],
             '-': [('AV', 'ls+'), ('Ap', 'ls-'), ('lam-', 'AsV'), ('lam+', 'Asp')]}[CH]
    TGT1 = {'+': [('c', 'lam+')], '-': [('c', 'lam-')]}[CH]
else:
    SEC0 = {'+': [('AV', 'lamP'), ('AW', 'lamM'), ('c', 'lamsM')], '-': [('AV', 'lamM'), ('AW', 'lamP'), ('c', 'lamsP')]}[CH]
    SECM1 = {'+': [('AV', 'lamsM'), ('AW', 'lamsP'), ('lamP', 'AsV'), ('lamM', 'AsW')],
             '-': [('AV', 'lamsP'), ('AW', 'lamsM'), ('lamM', 'AsV'), ('lamP', 'AsW')]}[CH]
    TGT1 = {'+': [('c', 'lamP')], '-': [('c', 'lamM')]}[CH]


def spec_of(t):
    """(species, projector or index range) of a sector letter in the route's conventions."""
    if ROUTE == 'B':
        m = {'AV': ('A', range(4)), 'Ap': ('A', range(4, 10)), 'AsV': ('As', range(4)), 'Asp': ('As', range(4, 10)),
             'c': ('c', None), 'cs': ('cs', None)}
        if t in m:
            return m[t]
        return P_sig(t)
    if t in RA.SPRNG:
        return RA.SPRNG[t]
    return (RA.G5[t][0], RA.chir_proj(t))


def projector(pr):
    def apply(leg, P):
        if P is None:
            return leg
        if isinstance(P, np.ndarray):
            return dict(sp=leg['sp'], u=P @ leg['u'], p=leg['p'])
        u = np.zeros_like(leg['u']); u[list(P)] = leg['u'][list(P)]
        return dict(sp=leg['sp'], u=u, p=leg['p'])
    return lambda s: [apply(l, P) for l, P in zip(s, pr)]


def rand_u(t, rng):
    spc, P = spec_of(t)
    while True:
        u = np.array([rng.integers(-2, 3) for _ in range(DT.DIM[spc])], dtype=complex)
        if isinstance(P, np.ndarray):
            u = P @ u * 2
        elif P is not None:
            w = np.zeros_like(u); w[list(P)] = u[list(P)]; u = w
        if np.any(u):
            return np.round(u)


def rand_spinor(rng, t):
    spc, P = spec_of(t)
    return rand_u(t, rng)


def cands(sector, deg, ckey='eta'):
    """PFn candidates (intrinsic + x^mu parts) of a sector at weight 7/2 with the derivative counts of the census."""
    t0, t1 = sector
    (s0, P0), (s1, P1) = spec_of(t0), spec_of(t1)
    slots = (s0, s1)
    pr = projector((P0, P1))
    cs = 'Ss'
    if deg == 0:
        Dint, Dpos = 1, 2
    else:
        Dint, Dpos = (0, 1) if s0 == 'A' else (None, 0)
    out = []
    if Dint is not None:
        for c in BL.lin_cands(RL, slots, Dint, vecs_extra=(), ckey=ckey, cspace=cs):
            out.append((c.label, slots, lambda ctx, c=c: DT.PFn(Q.dens(1, slots, lambda s: c.L(pr(s), ctx)), None)))
    for c in BL.lin_cands(RL, slots, Dpos, vecs_extra=('e',), ckey=ckey, cspace=cs):
        def mk(ctx, c=c):
            ws = []
            for mu in range(4):
                e = np.zeros(4, complex); e[mu] = 1
                cx = dict(ctx); cx['e'] = e
                ws.append(Q.dens(1, slots, lambda s, cx=cx: c.L(pr(s), cx)))
            return DT.PFn(None, ws)
        out.append(('x.' + c.label, slots, mk))
    return out


def legs_at_of(t0, t1, u0, u1, p):
    (s0, _), (s1, _) = spec_of(t0), spec_of(t1)
    return lambda r: [dict(sp=s0, u=u0, p=list(p)), dict(sp=s1, u=u1, p=[-p[m] - r[m] for m in range(4)])]


# ------------------------------------------------------------------ loops
def run_loops(nconf=8, seed=909):
    fn = os.path.join(OUT, f'{ROUTE}_eta{"p" if CH == "+" else "m"}_ob1_raw.json')
    if os.path.exists(fn):
        return json.load(open(fn))
    rng = np.random.default_rng(seed + (0 if ROUTE == 'A' else 1) + (0 if CH == '+' else 2))
    out = {}
    etat = ('ls+' if CH == '+' else 'ls-') if ROUTE == 'B' else ('lamsP' if CH == '+' else 'lamsM')
    if ROUTE == 'B':
        tq, t3 = G.tmpl_xi()
        def s0tmpl(eta):
            import sympy as sy
            ev = sy.Matrix([sy.nsimplify(complex(z).real) + sy.I * sy.nsimplify(complex(z).imag) for z in eta])
            def L_s0(s):
                a, ls = s
                w = sy.zeros(16, 1)
                for I in range(4, 10):
                    w += a['u'][I] * (B.C.Rbm[I] * ev)
                return (w.T * B.C.Mp * sy.Matrix(ls['u']))[0]
            B.V.TERMS['CS_s0'] = (sy.Integer(-1), ('A', 'lams'), L_s0, ('k', (0, 1)))
            return E.build_template(['CS_s0'], [], extra_odd=0, name=None)
    for sec in SEC0:
        recs = []
        for ic in range(nconf):
            eta = rand_spinor(rng, etat)
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            while not any(p):
                p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            u0, u1 = rand_u(sec[0], rng), rand_u(sec[1], rng)
            la = legs_at_of(sec[0], sec[1], u0, u1, p)
            t0 = time.time()
            if ROUTE == 'B':
                gdir = lambda mu, r, eta=eta: Y.Gen('xi', 0, tq, gvec=RB_[mu] @ eta, t3=t3, inj=r)
                al, be, top = P.position_split(la, gdir, Y.Gen('s0', 0, s0tmpl(eta), gvec=None), stats={})
                def cv(d, m):
                    return FL.clog_parse(d['k'][m].to_json()) if 'k' in d and m in d['k'] else FL.lvec_zero()
                alpha = cv(al, 0); pa = cv(al, 1)
                betas = [cv(be.get(mu, {}), 0) for mu in range(4)]
                topz = all(v.is_zero() for d in top.values() for v in d.values())
            else:
                L = [A.SA.leg(sec[0], u0, np.array([int(x) for x in p])), A.SA.leg(sec[1], u1, -np.array([int(x) for x in p]))]
                d = A.SO.density({'kind': 'eta', 'eta': eta}, L, what=('fin', ('pole', 1)))
                alpha = RA.v5_to_lvec(d['fin']['alpha'].get('(Y0,Y1)', A.SA.vzero()))
                pa = RA.v5_to_lvec(d[('pole', 1)]['alpha'].get('(Y0,Y1)', A.SA.vzero()))
                # engine A: beta^mu = -i G_mu(0); here beta_mu = G_mu(0) = i beta_A
                betas = [FL.lvec_scale(RA.v5_to_lvec(d['fin']['beta'][mu].get('(Y0,Y1)', A.SA.vzero())), (Fr(0), Fr(1))) for mu in range(4)]
                topz = True
            js = lambda lv: {k: [str(v[0]), str(v[1])] for k, v in lv.items()}
            recs.append(dict(eta=[[z.real, z.imag] for z in eta], p=[str(x) for x in p], u0=[[z.real, z.imag] for z in u0],
                             u1=[[z.real, z.imag] for z in u1], alpha=js(alpha), pole=js(pa), beta=[js(b) for b in betas],
                             top_zero=topz, sec=round(time.time() - t0, 1)))
            log(f'eta{CH} sector {sec} config {ic}: {recs[-1]["sec"]} s')
        out['%s|%s' % sec] = recs
        json.dump(out, open(fn, 'w'))
    return out


def lv_of(d):
    return {k: (Fr(v[0]), Fr(v[1])) for k, v in d.items()}


def cvec(a):
    return np.array([complex(x[0], x[1]) for x in a])


def main():
    raw = run_loops()
    confs = []
    for key, recs in raw.items():
        t0, t1 = key.split('|')
        for r in recs:
            la = legs_at_of(t0, t1, cvec(r['u0']), cvec(r['u1']), [Fr(x) for x in r['p']])
            ob = [lv_of(r['alpha'])] + [lv_of(b) for b in r['beta']]
            confs.append(dict(la=la, ctx=dict(eta=cvec(r['eta'])), ob=ob, pole=lv_of(r['pole']), top=r['top_zero']))
    pole_zero = all(all(v == (0, 0) for v in c['pole'].values()) for c in confs)
    # degree-0 basis on the loop configurations
    c0 = [c for s in SEC0 for c in cands(s, 0)]
    def vals(fn_maker, c):
        return DT.canon(fn_maker(c['ctx']), c['la'])
    rows = [[z for c in confs for z in vals(mk, c)] for (_, _, mk) in c0]
    keep = FL.independent_rows(rows)
    basis0 = [c0[i] for i in keep]; brows = [rows[i] for i in keep]
    log(f'eta{CH}: degree-0 candidates {len(c0)}, independent {len(basis0)}')
    # obstruction coordinates (5 canonical numbers per configuration)
    obl = [lv for c in confs for lv in c['ob']]
    obc = FL.lsolve_coords(brows, obl)
    ob_in = all(v is not None for v in obc.values())
    obv = FL.combine(obc) if ob_in else None
    # targets of degree +1 for H^0: (c, lam) with p <= 1
    rng = np.random.default_rng(55)
    etat = ('ls+' if CH == '+' else 'ls-') if ROUTE == 'B' else ('lamsP' if CH == '+' else 'lamsM')
    tg = []
    for (t0, t1) in TGT1:
        for _ in range(10):
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            tg.append(dict(la=legs_at_of(t0, t1, rand_u(t0, rng), rand_u(t1, rng), p), ctx=dict(eta=rand_spinor(rng, etat))))
    qrows0 = [[FL.gq(z) for t in tg for z in DT.canon(TT.Q(mk(t['ctx'])), t['la'])] for (_, _, mk) in basis0]
    ker0 = sp.Matrix(qrows0).T.nullspace()
    # degree -1: named family N (x^mu T_i(rho(e_mu) eta), i = 0..4, s0, u0 = int(rho(A_V) eta, lam*)) and the generic
    # candidates; independence on random degree -1 configurations, Q images in the degree-0 coordinates
    def named_m1(ctx):
        e = ctx['eta']
        out = {}
        for i in ('T0', 'T1', 'T2', 'T3', 'T4'):
            ws = []
            for mu in range(4):
                xi = RB_[mu] @ e
                if i == 'T0':
                    ws.append(Q.dens(1, ('A', 'lams'), lambda s, xi=xi: (-np.dot(1j * np.array([complex(x) for x in s[0]['p']]), s[0]['u'][:4])) * (xi @ MP_ @ s[1]['u'])))
                else:
                    ws.append(xi_terms(xi)[i])
            out['x.' + i + '(rho(e)eta)'] = DT.PFn(None, ws)
        out['s0'] = DT.PFn(Q.dens(-1, ('A', 'lams'), lambda s: sum(s[0]['u'][I] * ((RB_[I] @ e) @ MP_ @ s[1]['u']) for I in range(4, 10))), None)
        out['u0=int(rho(A_V)eta,lam*)'] = DT.PFn(Q.dens(1, ('A', 'lams'), lambda s: sum(s[0]['u'][M] * ((RB_[M] @ e) @ MP_ @ s[1]['u']) for M in range(4))), None)
        return out
    rngm = np.random.default_rng(77)
    cfm = []
    for (t0, t1) in SECM1:
        for _ in range(8):
            p = [Fr(int(x)) for x in rngm.integers(-3, 4, size=4)]
            cfm.append(dict(la=legs_at_of(t0, t1, rand_u(t0, rngm), rand_u(t1, rngm), p), ctx=dict(eta=rand_spinor(rngm, etat))))
    nm1 = list(named_m1(cfm[0]['ctx']).keys())
    allm1 = [(k, None, (lambda ctx, k=k: named_m1(ctx)[k])) for k in nm1] + [(lab, sl, mk) for (lab, sl, mk) in [c for s in SECM1 for c in cands(s, -1)]]
    rws = [[z for c in cfm for z in DT.canon(mk(c['ctx']), c['la'])] for (_, _, mk) in allm1]
    keepm = FL.independent_rows(rws)
    basm = [allm1[i] for i in keepm]
    res_named_indep = [b[0] for b in basm if b[0] in nm1]
    log(f'eta{CH}: degree -1: {len(allm1)} candidates, independent {len(basm)}; named independent: {res_named_indep}')
    im, imlab = [], []
    for (lab, slots, mk) in basm:
        r = [z for c in confs for z in DT.canon(TT.Q(mk(c['ctx'])), c['la'])]
        x = FL.solve_coords(brows, r)
        assert x is not None, 'Q image outside the degree-0 basis'
        im.append(sp.Matrix(x)); imlab.append(lab)
    imM = sp.Matrix.hstack(*im)
    ker0M = sp.Matrix.hstack(*ker0) if ker0 else sp.zeros(len(basis0), 0)
    H0 = ker0M.rank() - imM.rank()
    # closed degree -1 functionals = kernel of the image map on the degree -1 candidates (they map injectively?)
    kerm1 = imM.nullspace()
    # named functionals
    named = {}
    S_e = lambda ctx: S_eta(ctx['eta'])
    F10 = Q.dsum((27 / 16, Q.F_dsA_dsA)) if False else Q.F_dsA_dsA
    for c in confs:
        e = c['ctx']['eta']
        for nm, fn in (('QS(eta)', TT.Q(S_eta(e))),
                       ('{dsAdsA,S(eta)}', TT.br(F10, 0, S_eta(e))),
                       ('{ddcdA,S(eta)}', TT.br(Q.F_dsdc_dsA, 1, S_eta(e)))):
            named.setdefault(nm, []).extend(DT.canon(fn, c['la']))
        # G-pos: position part = Ob1(rho(e_mu) eta): build from g's
        for nm, coef in (('pos_g13', (1, 0)), ('pos_g2', (0, 1))):
            ws = []
            for mu in range(4):
                g = g_terms(RB_[mu] @ e)
                ws.append(Q.dsum((coef[0], g['g1']), (coef[0], g['g3']), (coef[1], g['g2'])))
            named.setdefault(nm, []).extend(DT.canon(DT.PFn(None, ws), c['la']))
    nc = {k: FL.solve_coords(brows, v) for k, v in named.items()}
    QS_zero = all(z == (0, 0) for z in named['QS(eta)'])
    res = dict(route=ROUTE, chirality=CH, n_configs=len(confs), basis0=[b[0] for b in basis0], pole_sum_zero=pole_zero,
               top_zero=all(c['top'] for c in confs), ob_in_span=ob_in, G_eta_closed=QS_zero,
               closed0=ker0M.rank(), exact0=imM.rank(), H0=H0, deg_m1_basis=len(basm), named_m1_independent=res_named_indep, H_m1_kernel=len(kerm1))
    log(f'eta{CH}: Ob in span {ob_in}, poles 0 {pole_zero}, Q S(eta) = 0 {QS_zero}, closed0 {ker0M.rank()}, exact0 '
        f'{imM.rank()}, H^0 {H0}, closed degree -1 {len(kerm1)}')
    if ob_in:
        # G-pos: beta part of Ob equals (79/24 - L/2) pos_g13 + 27/8 pos_g2 at the level of the full functional:
        # compare the position components only: build the candidate x-part and compare the beta-rows
        pos_pred = sp.Rational(79, 24) * sp.Matrix(nc['pos_g13']) - sp.Rational(1, 2) * L32 * sp.Matrix(nc['pos_g13']) \
            + sp.Rational(27, 8) * sp.Matrix(nc['pos_g2']) if nc['pos_g13'] is not None else None
        # beta rows of Ob vs beta rows of the prediction (evaluate both on the configurations)
        Bm = sp.Matrix([[FL.gq(z) for z in r] for r in brows]).T
        ob_vals = Bm * obv; pr_vals = Bm * pos_pred
        nconf = len(confs)
        dpos = [sp.nsimplify(sp.expand(ob_vals[5 * i + 1 + mu] - pr_vals[5 * i + 1 + mu])) for i in range(nconf) for mu in range(4)]
        res['G_pos'] = all(z == 0 for z in dpos)
        # G-WZ1: Q Ob1 = -{Ob0, S(eta)}, Ob0 = -27/8 ddcdA  ->  Q Ob1 = 27/8 {ddcdA, S(eta)} on degree +1 targets
        lhs = sp.zeros(len(qrows0[0]), 1)
        for i in range(len(basis0)):
            lhs += obv[i] * sp.Matrix(qrows0[i])
        rhs = sp.Matrix([FL.gq(z) for t in tg for z in DT.canon(TT.br(Q.F_dsdc_dsA, 1, S_eta(t['ctx']['eta'])), t['la'])])
        res['G_WZ1'] = all(sp.nsimplify(sp.expand(z)) == 0 for z in (lhs - sp.Rational(27, 8) * rhs))
        res['G_WZ1_printed_sign'] = all(sp.nsimplify(sp.expand(z)) == 0 for z in (lhs + sp.Rational(27, 8) * rhs))
        # residual and the solution of Q F1 = R1
        R1 = obv - sp.Rational(27, 16) * sp.Matrix(nc['{dsAdsA,S(eta)}'])
        ys = sp.symbols('f0:%d' % len(im))
        eqs = list(imM * sp.Matrix(ys) - R1)
        sol = FL.lin_solve(eqs, list(ys))
        res['R1_closed'] = bool(sp.Matrix.hstack(ker0M, R1).rank() == ker0M.rank())
        res['solvable'] = bool(sol)
        if sol:
            s0 = sol[0]
            free = [y for y in ys if s0[y] == y]
            res['F1_particular'] = {imlab[i]: FL.lstr(s0[ys[i]].subs({f: 0 for f in free})) for i in range(len(ys)) if s0[ys[i]].subs({f: 0 for f in free}) != 0}
            res['F1_free_directions'] = [{imlab[i]: FL.lstr(sp.diff(s0[ys[i]], f)) for i in range(len(ys)) if sp.diff(s0[ys[i]], f) != 0} for f in free]
            # express S^V(eta), S^perp(eta) in the degree -1 candidates (they must span the free directions)
            SVc, SPc = [], []
            for c in confs:
                SVc += DT.canon(S_eta(c['ctx']['eta'], ('V',)), c['la'])
            res['n_free'] = len(free)
        log(f'eta{CH}: G-pos {res["G_pos"]}, G-WZ1 {res["G_WZ1"]} (printed sign {res["G_WZ1_printed_sign"]}), R1 closed '
            f'{res["R1_closed"]}, solvable {res["solvable"]}, free directions {res.get("n_free")}')
        res['ob_coords'] = {k: [FL.lstr(z) for z in v] for k, v in obc.items()}
    json.dump(res, open(os.path.join(OUT, f'{ROUTE}_eta{"p" if CH == "+" else "m"}_fit.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
