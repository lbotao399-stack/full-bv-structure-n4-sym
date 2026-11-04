"""Route-independent analysis of the two-leg arity-two obstruction of odd pairs (xi1, xi2): H^{-1}(Q) on the degree -1
weight-5 block bases, classes of Ob^(2), G-pole2, G-WZ2 (corrected signs, finding F1), the residual
R2 = Ob2 - F1([x,y]) + {F1(xi1),S(xi2)} + {S(xi1),F1(xi2)} (eq at-eqs) and the conditions on the free coefficients
for [R2] = 0 in H^{-1}(Q).  A route supplies: configurations (legs in generic form, ctx = spinors, Ob2, poles),
candidate densities in its realisation, Q of densities, the named tree-level functionals (its own Q and bracket), and
the right side of G-WZ2 on degree-0 targets."""
import sympy as sp
import fitlib as FL

BLOCKS = {
    'pm': {'V': dict(sectors=['c_cs', 'AV_AsV'], deg2=[('Acs', 0)], count=5),
           'P': dict(sectors=['Ap_Asp'], deg2=[('Acs', 0)], count=2),
           'L': dict(sectors=['l+_ls+', 'l-_ls-'], deg2=[('lsls', 0)], count=8)},
    'pp': {'M': dict(sectors=['AV_Asp', 'Ap_AsV'], deg2=[('Acs', 0)], count=4),
           'L': dict(sectors=['l+_ls-', 'l-_ls+'], deg2=[('lsls', 0)], count=6)},
}
BLOCKS['mm'] = BLOCKS['pp']
TH = sp.symbols('thP_V thP_perp thP_lam thx1_V thx1_perp thx2_V thx2_perp')


def analyse(route, pair, log=print):
    cfg = route.load(pair)
    tg = route.deg0_targets(pair)
    tgw = route.wz_targets(pair) if hasattr(route, 'wz_targets') else tg
    res = {'pair': pair, 'route': route.name, 'blocks': {}}
    QOB = None
    L = sp.log(3) - sp.log(2)
    for bname, blk in BLOCKS[pair].items():
        confs = [c for s in blk['sectors'] for c in cfg[s]]
        cands = [cd for s in blk['sectors'] for cd in route.sector_cands(pair, s)]
        rows = [[route.dens_value(cd, c) for c in confs] for cd in cands]
        keep = FL.independent_rows(rows)
        basis = [cands[i] for i in keep]; brows = [rows[i] for i in keep]
        qrows = [[FL.gq(route.Q_value(cd, t)) for t in tg] for cd in basis]
        kerM_list = sp.Matrix(qrows).T.nullspace()
        img = []
        for (k, nm) in blk['deg2']:
            for cd in route.cands(k, nm, None):
                r = [route.Q_value(cd, c) for c in confs]
                if any(z != (0, 0) for z in r):
                    x = FL.solve_coords(brows, r)
                    assert x is not None, 'Q image outside the degree -1 basis'
                    img.append(sp.Matrix(x))
        imM = sp.Matrix.hstack(*img) if img else sp.zeros(len(basis), 0)
        kerM = sp.Matrix.hstack(*kerM_list) if kerM_list else sp.zeros(len(basis), 0)
        obc = FL.lsolve_coords(brows, [c['ob'] for c in confs])
        pole_zero = all(all(v == (0, 0) for v in c['pole'].values()) for c in confs)
        ob_in_span = all(v is not None for v in obc.values())
        obv = FL.combine(obc) if ob_in_span else None
        nf = {}
        for c in confs:
            for k, v in route.named(c, pair).items():
                nf.setdefault(k, []).append(v)
        nfc = {k: (None if FL.solve_coords(brows, r) is None else sp.Matrix(FL.solve_coords(brows, r))) for k, r in nf.items()}
        if ob_in_span:
            lhs = sp.zeros(len(tgw), 1)
            for i, cd in enumerate(basis):
                lhs += obv[i] * sp.Matrix([FL.gq(route.Q_value(cd, t)) for t in tgw])
            QOB = lhs if QOB is None else QOB + lhs
        cond, Rclosed = None, None
        if ob_in_span and all(v is not None for v in nfc.values()):
            f = L - sp.Rational(79, 12)
            R = obv + f * nfc['{T12(xi1),S(xi2)}+{S(xi1),T12(xi2)}'] \
                + TH[3] * nfc['{S^V(xi1),S(xi2)}'] + TH[4] * nfc['{S^perp(xi1),S(xi2)}'] \
                + TH[5] * nfc['{S(xi1),S^V(xi2)}'] + TH[6] * nfc['{S(xi1),S^perp(xi2)}']
            if pair == 'pm':
                R = R - TH[0] * nfc['S(P_G)_V'] - TH[1] * nfc['S(P_G)_perp'] - TH[2] * nfc['S(P_G)_lam']
            ys = sp.symbols('y0:%d' % imM.shape[1]) if imM.shape[1] else ()
            eqs = [sp.expand(e) for e in (R - (imM * sp.Matrix(ys) if ys else sp.zeros(len(basis), 1)))]
            sol = sp.solve(eqs, list(ys) + list(TH), dict=True)
            th_only = [{str(k): FL.lstr(v) for k, v in s.items() if str(k).startswith('th')} for s in sol]
            cond = dict(solution_theta=th_only, solvable=bool(sol))
            Rclosed = bool(sp.Matrix.hstack(kerM, R).rank() == kerM.rank()) if kerM.shape[1] else None
        res['blocks'][bname] = dict(
            basis=[cd.label for cd in basis], basis_dim=len(basis), expected=blk['count'], closed_dim=kerM.rank(),
            exact_dim=imM.rank(), H_minus1_dim=kerM.rank() - imM.rank(),
            ob_coords={k: (None if v is None else [FL.lstr(z) for z in v]) for k, v in obc.items()},
            ob_in_span=ob_in_span, pole_sum_zero=pole_zero,
            ob_closed=(None if not ob_in_span or not kerM.shape[1] else bool(sp.Matrix.hstack(kerM, obv).rank() == kerM.rank())),
            residual=cond, residual_closed=Rclosed)
        log(f'{route.name} pair {pair} block {bname}: basis {len(basis)}/{blk["count"]}, closed {kerM.rank()}, exact '
            f'{imM.rank()}, H^-1 {kerM.rank() - imM.rank()}, Ob in span {ob_in_span}, poles 0 {pole_zero}, '
            f'theta: {cond and cond["solution_theta"]}')
    rr = [route.wz_rhs(t) for t in tgw]
    wzres = [sp.nsimplify(sp.expand(z)) for z in (QOB - sp.Matrix(rr))]
    res['G_WZ2_full'] = dict(pass_=all(z == 0 for z in wzres), QOb2_zero=all(sp.expand(z) == 0 for z in QOB),
                             rule='Q Ob2 = {Ob1(xi1),S(xi2)} + {S(xi1),Ob1(xi2)} - {Ob0,S^(2)} on two legs (Ob1(P)|2 = 0)')
    log(f'{route.name} pair {pair}: G-WZ2 (full functional) {res["G_WZ2_full"]["pass_"]}, Q Ob2|2 = 0: {res["G_WZ2_full"]["QOb2_zero"]}')
    return res
