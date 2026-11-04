"""(xi, eta) on two legs, p <= 1, both routes: Ob^(2)(xi,eta)|_2 by the defect cycles of
J = I + c S(xi) + c' S(eta) - c c' S^(2)(xi, rho(x) eta) (eta and S^(2) by momentum injection), theory = -raw(two
letters) + raw(S^(2)); the bracket [xi,eta] = c_D Delta + M_omega + R_r identified at tree level from
S([xi,eta]) = {S(xi),S(eta)} - Q S^(2)(xi,eta) (gate G-bracket: exact fit onto S(Delta), S(M_omega), S(R_r));
residual R2 = Ob2 - F1([xi,eta]) + {F1(xi),S(eta)} + {S(xi),F1(eta)} ({F1^(0),S^(2)} = 0), with
F1(xi)|_2 = (L - 79/12)(T1+T2) + th^xi_V S^V_2 + th^xi_perp S^perp_2, F1(eta)|_2 from eta1.py (+ th^eta_V S^V(eta) +
th^eta_perp S^perp(eta)), F1(Delta)|_2 = F1part(Delta) + sum th^D_d S(Delta)_d, F1(M)|_2 = 1/3 S(M)_spin +
sum th^M_d S(M)_d (engines A, B task T3; finding F8), F1(R)|_2 = sum th^R_d S(R)_d (Ob1(R) = 0).
usage: python3 xieta.py A|B +|- (eta chirality; xi in S_+ (x) W)"""
import sys, os, json, time
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import fitlib as FL, bilin as BL, dtree as DT
ROUTE, ECH = sys.argv[1], sys.argv[2]
sys.argv = [sys.argv[0], ROUTE, ECH]
import eta1 as H
import dxi as DX
Q, TT, RL = H.Q, H.TT, H.RL
log = lambda s: H.log('(xi,eta) ' + s)
OUT = H.OUT
L32 = sp.log(3) - sp.log(2)
R_, RB_, MP_ = H.R_, H.RB_, H.MP_
if ROUTE == 'B':
    XIT, ETT = 'lam+', ('ls+' if ECH == '+' else 'ls-')
    SECM1 = {'+': [('c', 'cs'), ('AV', 'AsV'), ('Ap', 'Asp'), ('lam+', 'ls+'), ('lam-', 'ls-')],
             '-': [('AV', 'Asp'), ('Ap', 'AsV'), ('lam+', 'ls-'), ('lam-', 'ls+')]}[ECH]
    SECM2 = {'+': [('AV', 'cs'), ('ls+', 'ls-')], '-': [('Ap', 'cs'), ('ls+', 'ls+'), ('ls-', 'ls-')]}[ECH]
    SEC0 = {'+': [('c', 'AsV'), ('AV', 'AV'), ('Ap', 'Ap'), ('lam+', 'lam-')], '-': [('c', 'Asp'), ('AV', 'Ap'), ('lam+', 'lam+'), ('lam-', 'lam-')]}[ECH]
else:
    XIT, ETT = 'lamP', ('lamsP' if ECH == '+' else 'lamsM')
    SECM1 = {'+': [('c', 'cs'), ('AV', 'AsV'), ('AW', 'AsW'), ('lamP', 'lamsP'), ('lamM', 'lamsM')],
             '-': [('AV', 'AsW'), ('AW', 'AsV'), ('lamP', 'lamsM'), ('lamM', 'lamsP')]}[ECH]
    SECM2 = {'+': [('AV', 'cs'), ('lamsP', 'lamsM')], '-': [('AW', 'cs'), ('lamsP', 'lamsP'), ('lamsM', 'lamsM')]}[ECH]
    SEC0 = {'+': [('c', 'AsV'), ('AV', 'AV'), ('AW', 'AW'), ('lamP', 'lamM')], '-': [('c', 'AsW'), ('AV', 'AW'), ('lamP', 'lamP'), ('lamM', 'lamM')]}[ECH]


def cands(sector, Dint, Dpos):
    t0, t1 = sector
    (s0, P0), (s1, P1) = H.spec_of(t0), H.spec_of(t1)
    slots = (s0, s1); pr = H.projector((P0, P1))
    out = []
    if Dint is not None:
        for c in BL.bil2_cands(RL, slots, Dint, ()):
            out.append((c.label, lambda ctx, c=c: DT.PFn(Q.dens(1, slots, lambda s: c.L(pr(s), ctx)), None)))
    if Dpos is not None:
        for c in BL.bil2_cands(RL, slots, Dpos, ('e',)):
            def mk(ctx, c=c):
                ws = []
                for mu in range(4):
                    e = np.zeros(4, complex); e[mu] = 1
                    cx = dict(ctx); cx['e'] = e
                    ws.append(Q.dens(1, slots, lambda s, cx=cx: c.L(pr(s), cx)))
                return DT.PFn(None, ws)
            out.append(('x.' + c.label, mk))
    return out


def gam(Q1, Q2):
    return np.array([Q1 @ MP_ @ R_[M] @ Q2 for M in range(10)])


def gam_star(w1, w2):
    return np.array([w1 @ RB_[M].T @ MP_ @ w2 for M in range(10)])


def s2_dens(Q1, Q2):
    G12 = gam(Q1, Q2)
    L = lambda s: 0.25 * (G12 @ gam_star(s[0]['u'], s[1]['u'])) - 0.5 * (Q1 @ MP_ @ s[0]['u']) * (Q2 @ MP_ @ s[1]['u'])
    return Q.dsum((1, Q.dens(1, ('lams', 'lams'), L)), (1, Q.dens(-1, ('A', 'cs'), lambda s: (G12 @ s[0]['u']) * s[1]['u'][0])))


def S2_xieta(xi, eta):
    return DT.PFn(None, [s2_dens(xi, RB_[mu] @ eta) for mu in range(4)])


def trans(b, parts=('V', 'perp', 'lam')):
    bd = lambda leg: 1j * np.dot([complex(x) for x in leg['p']], b)
    t = []
    if 'V' in parts:
        t += [(1, Q.dens(1, ('A', 'As'), lambda s: bd(s[0]) * np.dot(s[0]['u'][:4], s[1]['u'][:4]))),
              (1, Q.dens(-1, ('c', 'cs'), lambda s: bd(s[0]) * s[0]['u'][0] * s[1]['u'][0]))]
    if 'perp' in parts:
        t += [(1, Q.dens(1, ('A', 'As'), lambda s: bd(s[0]) * np.dot(s[0]['u'][4:], s[1]['u'][4:])))]
    if 'lam' in parts:
        t += [(1, Q.dens(-1, ('lam', 'lams'), lambda s: bd(s[0]) * (s[0]['u'] @ MP_ @ s[1]['u'])))]
    return Q.dsum(*t)


def S_M(om, parts=('V', 'perp', 'lam'), spin_only=False):
    """S(M_omega), v = omega x (v^r = om[r,m] x^m): x^m translation by b = om[:, m]; spin: V: om_{rn}(A_r, A*_n),
    lam: +1/4 sum om_{mn} (rho_m rho_n lam, lam*) (-(v.lam, lam*), v.lam spin part -1/4 om rho rho)."""
    w0 = []
    if 'V' in parts:
        w0.append((1, Q.dens(1, ('A', 'As'), lambda s: sum(om[r, n] * s[0]['u'][r] * s[1]['u'][n] for r in range(4) for n in range(4)))))
    if 'lam' in parts:
        Rl = sum(om[m, n] * (RB_[m] @ R_[n]) for m in range(4) for n in range(4))
        w0.append((0.25, Q.dens(1, ('lam', 'lams'), lambda s: (Rl @ s[0]['u']) @ MP_ @ s[1]['u'])))
    wm = None if spin_only else [trans(om[:, m], parts) for m in range(4)]
    return DT.PFn(Q.dsum(*w0) if w0 else None, wm)


def S_R(r, parts=('perp', 'lam')):
    """S(R_r) = int (r A_perp, A*_perp) - int (r.lam, lam*), r.lam = 1/4 sum r_IJ rho_I rho_J lam (6x6 r on 4..9)."""
    w0 = []
    if 'perp' in parts:
        w0.append((1, Q.dens(1, ('A', 'As'), lambda s: sum(r[I, J] * s[0]['u'][4 + J] * s[1]['u'][4 + I] for I in range(6) for J in range(6)))))
    if 'lam' in parts:
        Rl = sum(r[I, J] * (RB_[4 + I] @ R_[4 + J]) for I in range(6) for J in range(6))
        w0.append((-0.25, Q.dens(1, ('lam', 'lams'), lambda s: (Rl @ s[0]['u']) @ MP_ @ s[1]['u'])))
    return DT.PFn(Q.dsum(*w0), None)


OMB = [(i, j) for i in range(4) for j in range(i + 1, 4)]
RBS = [(i, j) for i in range(6) for j in range(i + 1, 6)]


def antis(n, ij):
    m = np.zeros((n, n)); m[ij[0], ij[1]] = 1; m[ij[1], ij[0]] = -1
    return m


def bracket_fit(xi, eta, rng):
    """[xi,eta] = c_D Delta + M_omega + R_r from S_2([xi,eta]) = {S(xi),S(eta)} - Q S^(2) (exact fit on probe legs)."""
    S2x = Q.dsum(*[(1, v) for k, v in H.xi_terms(xi).items()])
    X = DT.psum((1, TT.br(S2x, 0, H.S_eta(eta))), (-1, TT.Q(S2_xieta(xi, eta))))
    basis = [DX.S_Delta()] + [S_M(antis(4, ij)) for ij in OMB] + [S_R(antis(6, ij)) for ij in RBS]
    rows, rhs = [[] for _ in basis], []
    for spc in (('A', 'As'), ('lam', 'lams'), ('c', 'cs'), ('As', 'A'), ('lams', 'lam')):
        for _ in range(4):
            u0 = np.array(rng.integers(-2, 3, size=DT.DIM[spc[0]]), complex); u1 = np.array(rng.integers(-2, 3, size=DT.DIM[spc[1]]), complex)
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            la = lambda r, u0=u0, u1=u1, p=p, spc=spc: [dict(sp=spc[0], u=u0, p=p), dict(sp=spc[1], u=u1, p=[-p[m] - r[m] for m in range(4)])]
            rhs += DT.canon(X, la)
            for i, b in enumerate(basis):
                rows[i] += DT.canon(b, la)
    x = FL.solve_coords(rows, rhs)
    if x is None:
        return None
    cD = x[0]; om = np.zeros((4, 4), complex); rr = np.zeros((6, 6), complex)
    for k, ij in enumerate(OMB):
        om += complex(x[1 + k]) * antis(4, ij)
    for k, ij in enumerate(RBS):
        rr += complex(x[7 + k]) * antis(6, ij)
    return complex(cD), om, rr


def run_loops(nconf=8, seed=5151):
    fn = os.path.join(OUT, f'{ROUTE}_xieta{"p" if ECH == "+" else "m"}_raw.json')
    if os.path.exists(fn):
        return json.load(open(fn))
    rng = np.random.default_rng(seed + (0 if ROUTE == 'A' else 1) + (0 if ECH == '+' else 2))
    out = {}
    if ROUTE == 'B':
        E, Y, G, P = H.E, H.Y, H.G, H.P
        tq, t3 = G.tmpl_xi()
        import cb_s2 as CS
    for sec in SECM1:
        recs = []
        for ic in range(nconf):
            xi, eta = H.rand_u(XIT, rng), H.rand_u(ETT, rng)
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            while not any(p):
                p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            u0, u1 = H.rand_u(sec[0], rng), H.rand_u(sec[1], rng)
            la = H.legs_at_of(sec[0], sec[1], u0, u1, p)
            t0 = time.time()
            if ROUTE == 'B':
                import sympy as sy
                ev = sy.Matrix([sy.nsimplify(complex(z).real) + sy.I * sy.nsimplify(complex(z).imag) for z in eta])
                def L_s0(s):
                    a, ls = s
                    w = sy.zeros(16, 1)
                    for I in range(4, 10):
                        w += a['u'][I] * (B_C.Rbm[I] * ev)
                    return (w.T * B_C.Mp * sy.Matrix(ls['u']))[0]
                B_C = H.B.C
                H.B.V.TERMS['XE_s0'] = (sy.Integer(-1), ('A', 'lams'), L_s0, ('k', (0, 1)))
                ts0 = E.build_template(['XE_s0'], [], extra_odd=0, name=None)
                gx = Y.Gen('xi', 0, tq, gvec=xi, t3=t3)
                gdir = lambda mu, r, eta=eta: Y.Gen('etamu', 0, tq, gvec=RB_[mu] @ eta, t3=t3, inj=r)
                a2, b2, top2 = P.position_split(la, gdir, Y.Gen('s0', 0, ts0), extra_gens=[gx], gen_first=True)
                tS2 = {mu: E.build_template(CS.s2_terms(xi, RB_[mu] @ eta, f'xe{mu}'), [], extra_odd=0, name=None) for mu in range(4)}
                gS2 = lambda mu, r: Y.Gen('S2', 0, tS2[mu], inj=r)
                a1, b1, top1 = P.position_split(la, gS2, None)
                def cv(d, m):
                    return FL.clog_parse(d['k'][m].to_json()) if 'k' in d and m in d['k'] else FL.lvec_zero()
                comb = lambda x2, x1: FL.lvec_add(FL.lvec_scale(x2, (Fr(-1), Fr(0))), x1)
                alpha, pa = comb(cv(a2, 0), cv(a1, 0)), comb(cv(a2, 1), cv(a1, 1))
                betas = [comb(cv(b2.get(mu, {}), 0), cv(b1.get(mu, {}), 0)) for mu in range(4)]
                topz = all(v.is_zero() for t in (top1, top2) for d in t.values() for v in d.values())
            else:
                A, RA = H.A, H.RA
                import ca_s2 as AS
                pa_ = np.array([int(x) for x in p])
                L = [A.SA.leg(sec[0], u0, pa_), A.SA.leg(sec[1], u1, -pa_)]
                d = A.SO.density({'kind': 'eta', 'eta': eta}, L, other={'kind': 'xi', 'xi': xi}, what=('fin', ('pole', 1)))
                neg = lambda v: FL.lvec_scale(RA.v5_to_lvec(v), (Fr(-1), Fr(0)))
                alpha = neg(d['fin']['alpha'].get('(Y0,Y1)', A.SA.vzero()))
                pa = neg(d[('pole', 1)]['alpha'].get('(Y0,Y1)', A.SA.vzero()))
                betas = [FL.lvec_scale(neg(d['fin']['beta'][mu].get('(Y0,Y1)', A.SA.vzero())), (Fr(0), Fr(1))) for mu in range(4)]
                I4 = np.eye(4, dtype=int)
                for mu in range(4):
                    Lm = A.SO._shift_last(L, mu)
                    r = A.SA.cycle_sum(Lm, [], s2=AS.gen_s2(xi, A.SS[mu] @ eta), r1=I4[mu], what=('fin', ('pole', 1)))
                    for w, tgt in (('fin', 'a'), (('pole', 1), 'p')):
                        d1 = A.SO.reduce_colour(r, 2, w, 1).get('(Y0,Y1)', A.SA.vzero())
                        d0 = A.SO.reduce_colour(r, 2, w, 0).get('(Y0,Y1)', A.SA.vzero())
                        # alpha += -i t^1 ; beta_mu = G_mu(0) (here: t^0 coefficient)
                        if tgt == 'a':
                            alpha = FL.lvec_add(alpha, FL.lvec_scale(RA.v5_to_lvec(d1), (Fr(0), Fr(-1))))
                            betas[mu] = FL.lvec_add(betas[mu], RA.v5_to_lvec(d0))
                        else:
                            pa = FL.lvec_add(pa, FL.lvec_scale(RA.v5_to_lvec(d1), (Fr(0), Fr(-1))))
                topz = True
            js = lambda lv: {k: [str(v[0]), str(v[1])] for k, v in lv.items()}
            recs.append(dict(xi=[[z.real, z.imag] for z in xi], eta=[[z.real, z.imag] for z in eta], p=[str(x) for x in p],
                             u0=[[z.real, z.imag] for z in u0], u1=[[z.real, z.imag] for z in u1], alpha=js(alpha),
                             pole=js(pa), beta=[js(b) for b in betas], top_zero=topz, sec=round(time.time() - t0, 1)))
            log(f'eta{ECH} sector {sec} config {ic}: {recs[-1]["sec"]} s')
        out['%s|%s' % sec] = recs
        json.dump(out, open(fn, 'w'))
    return out


def F1eta(eta, part):
    """the eta solution from eta1.py (route's own fit): particular part as PFn"""
    e = eta
    def xT(i, coef):
        ws = []
        for mu in range(4):
            xi = RB_[mu] @ e
            ws.append(H.xi_terms(xi)[i])
        return DT.PFn(None, ws)
    terms = []
    for lab, c in part.items():
        cs = sp.nsimplify(c)
        if lab.startswith('x.T'):
            terms.append((cs, xT(lab[2:4], cs)))
        elif lab == 's0':
            terms.append((cs, DT.PFn(Q.dens(-1, ('A', 'lams'), lambda s: sum(s[0]['u'][I] * ((RB_[I] @ e) @ MP_ @ s[1]['u']) for I in range(4, 10))), None)))
        elif lab.startswith('u0'):
            terms.append((cs, DT.PFn(Q.dens(1, ('A', 'lams'), lambda s: sum(s[0]['u'][M] * ((RB_[M] @ e) @ MP_ @ s[1]['u']) for M in range(4))), None)))
        else:
            raise ValueError(lab)
    return terms


def main():
    raw = run_loops()
    etafit = json.load(open(os.path.join(OUT, f'{ROUTE}_eta{"p" if ECH == "+" else "m"}_fit.json')))
    F1part_eta = etafit['F1_particular']
    confs = []
    for key, recs in raw.items():
        t0, t1 = key.split('|')
        for r in recs:
            la = H.legs_at_of(t0, t1, H.cvec(r['u0']), H.cvec(r['u1']), [Fr(x) for x in r['p']])
            confs.append(dict(la=la, ctx=dict(xi=H.cvec(r['xi']), eta=H.cvec(r['eta'])),
                              ob=[H.lv_of(r['alpha'])] + [H.lv_of(b) for b in r['beta']], pole=H.lv_of(r['pole']), top=r['top_zero']))
    res = dict(route=ROUTE, eta_chirality=ECH, pole_sum_zero=all(all(v == (0, 0) for v in c['pole'].values()) for c in confs),
               top_zero=all(c['top'] for c in confs))
    t0 = time.time()
    cm1 = [c for s in SECM1 for c in cands(s, 0, 1)]
    rows = [[z for c in confs for z in DT.canon(mk(c['ctx']), c['la'])] for (_, mk) in cm1]
    keep = FL.independent_rows(rows)
    basis = [cm1[i] for i in keep]; brows = [rows[i] for i in keep]
    log(f'eta{ECH}: degree -1 candidates {len(cm1)}, independent {len(basis)} ({time.time() - t0:.0f} s)')
    obc = FL.lsolve_coords(brows, [lv for c in confs for lv in c['ob']])
    res['ob_in_span'] = all(v is not None for v in obc.values())
    obv = FL.combine(obc) if res['ob_in_span'] else None
    im, imlab = [], []
    IMV, IMVL = [], []
    for s in SECM2:
        for (lab, mk) in cands(s, None, 0):
            r = [z for c in confs for z in DT.canon(TT.Q(mk(c['ctx'])), c['la'])]
            if any(z != (0, 0) for z in r):
                IMV.append(r); IMVL.append(lab)
                x = FL.solve_coords(brows, r)
                if x is None:
                    continue
                im.append(sp.Matrix(x)); imlab.append(lab)
    imM = sp.Matrix.hstack(*im) if im else sp.zeros(len(basis), 0)
    rng = np.random.default_rng(41)
    tg = []
    for (t0_, t1_) in SEC0:
        for _ in range(6):
            p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
            tg.append(dict(la=H.legs_at_of(t0_, t1_, H.rand_u(t0_, rng), H.rand_u(t1_, rng), p), ctx=dict(xi=H.rand_u(XIT, rng), eta=H.rand_u(ETT, rng))))
    log(f'eta{ECH}: exact images {len(im)} ({time.time() - t0:.0f} s)')
    qrows = [[z for t in tg for z in DT.canon(TT.Q(mk(t['ctx'])), t['la'])] for (_, mk) in basis]
    rkQ = len(FL.independent_rows(qrows))
    rkIm = len(FL.independent_rows([[ (FL.Fr(str(sp.re(v))), FL.Fr(str(sp.im(v)))) for v in col] for col in im])) if im else 0
    res.update(basis_dim=len(basis), closed=len(basis) - rkQ, exact=rkIm, H_minus1=len(basis) - rkQ - rkIm)
    if obv is not None:
        QOB = sp.zeros(len(qrows[0]), 1)
        for i in range(len(basis)):
            QOB += obv[i] * sp.Matrix([FL.gq(z) for z in qrows[i]])
        res['ob_closed'] = all(sp.expand(z) == 0 for z in QOB)
    else:
        res['ob_closed'] = None
    log(f'eta{ECH}: closedness done ({time.time() - t0:.0f} s)')
    # named functionals per configuration
    named = {}
    brk_ok = True
    rngp = np.random.default_rng(9)
    Fp = DX.F1part_Delta()
    for c in confs:
        xi, eta = c['ctx']['xi'], c['ctx']['eta']
        T = H.xi_terms(xi)
        S2x = Q.dsum(*[(1, T[k]) for k in ('T1', 'T2', 'T3', 'T4')])
        SVx, SPx = Q.dsum((1, T['T1']), (1, T['T3'])), Q.dsum((1, T['T2']), (1, T['T4']))
        T12 = Q.dsum((1, T['T1']), (1, T['T2']))
        items = {'{S^V(xi),S(eta)}': TT.br(SVx, 0, H.S_eta(eta)), '{S^perp(xi),S(eta)}': TT.br(SPx, 0, H.S_eta(eta)),
                 '{T12(xi),S(eta)}': TT.br(T12, 0, H.S_eta(eta)),
                 '{S(xi),S^V(eta)}': TT.br(S2x, 0, H.S_eta(eta, ('V',))), '{S(xi),S^perp(eta)}': TT.br(S2x, 0, H.S_eta(eta, ('perp',)))}
        items['{S(xi),F1part(eta)}'] = DT.psum(*[(complex(sp.N(cf)) if False else 1, TT.br(S2x, 0, fn)) for cf, fn in []]) if False else None
        parts = F1eta(eta, F1part_eta)
        cD, om, rr = bracket_formula(xi, eta)
        if len(named.get('S(M)_spin', [])) < 10:      # gate G-bracket on the first two configurations
            bf = bracket_fit(xi, eta, rngp)
            if bf is None or abs(bf[0] - cD) > 1e-9 or np.abs(bf[1] - om).max() > 1e-9 or np.abs(bf[2] - rr).max() > 1e-9:
                brk_ok = False
                break
        for k, fn in items.items():
            if fn is not None:
                named.setdefault(k, []).extend(DT.canon(fn, c['la']))
        # {S(xi), F1part(eta)}: linear in the symbolic coefficients -> store per part
        for i, (cf, fn) in enumerate(parts):
            named.setdefault('{S(xi),eta-part%d}' % i, []).extend(DT.canon(TT.br(S2x, 0, fn), c['la']))
        # F1([xi,eta]) pieces
        for d in ('perp', 'lam', 'V'):
            named.setdefault('cD S(D)_' + d, []).extend([FL._cm(FL.to_fr(cD), z) for z in DT.canon(DX.S_Delta((d,)), c['la'])])
            if d != 'V':
                named.setdefault('S(R)_' + d, []).extend(DT.canon(S_R(rr, (d,)), c['la']))
        for d in ('perp', 'lam', 'V'):
            named.setdefault('S(M)_' + d, []).extend(DT.canon(S_M(om, (d,)), c['la']))
        named.setdefault('S(M)_spin', []).extend(DT.canon(S_M(om, ('V', 'lam'), spin_only=True), c['la']))
        cz = FL.to_fr(cD)
        for kk in ('AAs', 'lsl', 'csc'):
            named.setdefault('cD F1part(D)_' + kk, []).extend([FL._cm(cz, z) for z in DT.canon(DT.PFn(Fp[kk], None), c['la'])])
    res['G_bracket'] = brk_ok
    log(f'eta{ECH}: named functionals done ({time.time() - t0:.0f} s)')
    nc = {k: FL.solve_coords(brows, v) for k, v in named.items()}
    res['named_outside_basis'] = [k for k, v in nc.items() if v is None]
    log(f'eta{ECH}: basis {len(basis)}, closed {res["closed"]}, exact {res["exact"]}, H^-1 {res["H_minus1"]}, Ob in span '
        f'{res["ob_in_span"]}, Ob closed {res["ob_closed"]}, poles 0 {res["pole_sum_zero"]}, G-bracket {brk_ok}, outside {res["named_outside_basis"]}')
    if brk_ok:
        # direct test on the value vectors (configurations x canonical components): [R] = 0 <=> R = Q psi on every
        # configuration; independent of the completeness of the candidate basis
        LS = FL.LSYM
        obvals = sp.Matrix([sum(FL.gq(lv[k]) * LS[k] for k in FL.LOGK) for c in confs for lv in c['ob']])
        V = {k: sp.Matrix([FL.gq(z) for z in v]) for k, v in named.items()}
        thx = sp.symbols('thx_V thx_perp'); the = sp.symbols('the_V the_perp')
        thD = sp.symbols('thD_perp thD_lam thD_V'); thM = sp.symbols('thM_perp thM_lam thM_V'); thR = sp.symbols('thR_perp thR_lam')
        f = L32 - sp.Rational(79, 12)
        etac = [sp.nsimplify(c) for c in F1part_eta.values()]
        F1xi_br = f * V['{T12(xi),S(eta)}'] + thx[0] * V['{S^V(xi),S(eta)}'] + thx[1] * V['{S^perp(xi),S(eta)}']
        F1eta_br = sum((cf * V['{S(xi),eta-part%d}' % i] for i, cf in enumerate(etac)), sp.zeros(len(obvals), 1)) \
            + the[0] * V['{S(xi),S^V(eta)}'] + the[1] * V['{S(xi),S^perp(eta)}']
        F1br = sp.Rational(22, 3) * V['cD F1part(D)_AAs'] - sp.Rational(47, 3) * V['cD F1part(D)_lsl'] - sp.Rational(8, 3) * V['cD F1part(D)_csc'] \
            + thD[0] * V['cD S(D)_perp'] + thD[1] * V['cD S(D)_lam'] + thD[2] * V['cD S(D)_V'] \
            + sp.Rational(1, 3) * V['S(M)_spin'] + thM[0] * V['S(M)_perp'] + thM[1] * V['S(M)_lam'] + thM[2] * V['S(M)_V'] \
            + thR[0] * V['S(R)_perp'] + thR[1] * V['S(R)_lam']
        R = obvals - F1br + F1xi_br + F1eta_br
        keepI = FL.independent_rows(IMV)
        Im = sp.Matrix([[FL.gq(z) for z in IMV[i]] for i in keepI]).T
        ys = sp.symbols('y0:%d' % len(keepI))
        unk = list(thx) + list(the) + list(thD) + list(thM) + list(thR)
        # rational and log parts separately: unknowns are numbers in Q + Q log2 + Q log3
        sub = {}
        for u in list(ys) + unk:
            sub[u] = sp.Symbol(str(u) + '_0') + sp.Symbol(str(u) + '_2') * sp.log(2) + sp.Symbol(str(u) + '_3') * sp.log(3)
        eqs = []
        for e in (R - Im * sp.Matrix(ys)):
            e = sp.expand(e.subs(sub))
            for part in (sp.log(2), sp.log(3)):
                eqs.append(e.coeff(part))
            eqs.append(e.subs({sp.log(2): 0, sp.log(3): 0}))
        allu = [sp.Symbol(str(u) + s_) for u in list(ys) + unk for s_ in ('_0', '_2', '_3')]
        sol = FL.lin_solve(eqs, allu)
        res['solvable'] = bool(sol)
        res['n_value_equations'] = len(R)
        if sol:
            s0 = sol[0]
            res['conditions'] = {}
            for u in unk:
                val = sum(s0[sp.Symbol(str(u) + s_)] * w for s_, w in (('_0', 1), ('_2', sp.log(2)), ('_3', sp.log(3))))
                free_self = all(s0[sp.Symbol(str(u) + s_)] == sp.Symbol(str(u) + s_) for s_ in ('_0', '_2', '_3'))
                if not free_self:
                    res['conditions'][str(u)] = FL.lstr(val)
        log(f'eta{ECH}: solvable {res["solvable"]}, conditions {res.get("conditions")}')
    json.dump(res, open(os.path.join(OUT, f'{ROUTE}_xieta{"p" if ECH == "+" else "m"}_fit.json'), 'w'), indent=1)



def bracket_formula(xi, eta):
    """[xi, eta] from eq c2-oddodd: v(x)^nu = Gamma(xi, rho(x) eta)^nu|_V = M[nu,mu] x^mu, M = (xi,eta) 1 + omega
    (Prop c2-special), rotation of W_perp (w_I, R w_J) = (rho(q(w_I ^ w_J)) xi, eta), q(w^w') = rho_I rho_J (I != J)."""
    M = np.array([[gam(xi, RB_[mu] @ eta)[nu] for mu in range(4)] for nu in range(4)])
    pr = xi @ MP_ @ eta
    om = M - pr * np.eye(4)
    r = np.zeros((6, 6), complex)
    for I in range(6):
        for J in range(6):
            if I != J:
                r[I, J] = -(RB_[4 + I] @ (R_[4 + J] @ xi)) @ MP_ @ eta
    # sign: r of S(R_r) (S_R convention) = minus the matrix elements of eq c2-oddodd as transcribed here; fixed by the
    # exact fit of S([xi,eta]) = {S(xi),S(eta)} - Q S^(2) (bracket_fit, gate G-bracket)
    return pr / 2, om, r


if __name__ == '__main__':
    main()
