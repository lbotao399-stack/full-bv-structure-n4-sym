"""Generic exact class analysis on two-leg blocks (any route, any generator pair):
  basis = independent candidates per block (evaluated on the loop configurations), closed subspace (Q on degree +1
  targets), exact image (Q of the lower-degree candidates), coordinates of the obstruction and of named tree-level
  functionals, the residual R(theta) and the conditions for [R] = 0 (R = Q psi, psi explicit)."""
import sympy as sp
import fitlib as FL


def analyse(name, blocks, cfg, cands_of_sector, deg2_cands, dens_value, Q_value, targets, named, residual, symbols,
            log=print):
    out = {'blocks': {}}
    QOB = None
    for bname, blk in blocks.items():
        confs = [c for s in blk['sectors'] for c in cfg[s]]
        cands = [cd for s in blk['sectors'] for cd in cands_of_sector(s)]
        rows = [[dens_value(cd, c) for c in confs] for cd in cands]
        keep = FL.independent_rows(rows)
        basis = [cands[i] for i in keep]; brows = [rows[i] for i in keep]
        qrows = [[FL.gq(Q_value(cd, t)) for t in targets] for cd in basis]
        ker = sp.Matrix(qrows).T.nullspace() if basis else []
        img, imlab = [], []
        for cd in deg2_cands(bname):
            r = [Q_value(cd, c) for c in confs]
            if any(z != (0, 0) for z in r):
                x = FL.solve_coords(brows, r)
                assert x is not None, 'Q image outside the basis'
                img.append(sp.Matrix(x)); imlab.append(cd.label)
        imM = sp.Matrix.hstack(*img) if img else sp.zeros(len(basis), 0)
        kerM = sp.Matrix.hstack(*ker) if ker else sp.zeros(len(basis), 0)
        obc = FL.lsolve_coords(brows, [c['ob'] for c in confs])
        pole_zero = all(all(v == (0, 0) for v in c['pole'].values()) for c in confs)
        ob_in_span = all(v is not None for v in obc.values())
        obv = FL.combine(obc) if ob_in_span else None
        nf = {}
        for c in confs:
            for k, v in named(c).items():
                nf.setdefault(k, []).append(v)
        nfc = {}
        for k, r in nf.items():
            x = FL.solve_coords(brows, r)
            nfc[k] = None if x is None else sp.Matrix(x)
        if ob_in_span and basis:
            lhs = sp.zeros(len(targets), 1)
            for i in range(len(basis)):
                lhs += obv[i] * sp.Matrix(qrows[i])
            QOB = lhs if QOB is None else QOB + lhs
        cond = None
        Rclosed = None
        if ob_in_span and all(v is not None for v in nfc.values()):
            R = residual(obv, nfc, bname)
            ys = sp.symbols('y0:%d' % imM.shape[1]) if imM.shape[1] else ()
            eqs = [sp.expand(e) for e in (R - (imM * sp.Matrix(ys) if ys else sp.zeros(len(basis), 1)))]
            sol = FL.lin_solve(eqs, list(ys) + list(symbols))
            cond = dict(solution_theta=[{str(k): FL.lstr(v) for k, v in s.items() if k in symbols and v != k} for s in sol],
                        psi=[{imlab[i]: FL.lstr(s.get(ys[i], ys[i])) for i in range(len(ys))} for s in sol],
                        solvable=bool(sol))
            Rclosed = bool(sp.Matrix.hstack(kerM, R).rank() == kerM.rank()) if kerM.shape[1] else all(sp.expand(z) == 0 for z in R)
        out['blocks'][bname] = dict(
            basis=[cd.label for cd in basis], basis_dim=len(basis), expected=blk.get('count'), closed_dim=kerM.rank(),
            exact_dim=imM.rank(), H_dim=kerM.rank() - imM.rank(),
            ob_coords={k: (None if v is None else [FL.lstr(z) for z in v]) for k, v in obc.items()},
            ob_in_span=ob_in_span, pole_sum_zero=pole_zero,
            ob_closed=(None if not ob_in_span else (bool(sp.Matrix.hstack(kerM, obv).rank() == kerM.rank()) if kerM.shape[1] else all(z == 0 for z in obv))),
            named_coords={k: (None if v is None else [FL.lstr(z) for z in v]) for k, v in nfc.items()},
            residual=cond, residual_closed=Rclosed)
        log(f'{name} block {bname}: basis {len(basis)}/{blk.get("count")}, closed {kerM.rank()}, exact {imM.rank()}, '
            f'H {kerM.rank() - imM.rank()}, Ob in span {ob_in_span}, Ob closed {out["blocks"][bname]["ob_closed"]}, '
            f'poles 0 {pole_zero}, theta: {cond and cond["solution_theta"]}')
    out['QOb'] = None if QOB is None else [FL.lstr(z) for z in QOB]
    return out
