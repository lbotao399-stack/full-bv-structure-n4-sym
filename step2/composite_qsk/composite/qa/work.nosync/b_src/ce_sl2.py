"""Two-site families with covariant derivatives: tr((n.D)^a X (n.D)^b X), D = d + [A, -], X = Z (F1 twist two) or
X = psi (F2 fermion sl(2)); two-field part on collinear on-shell outputs.
Composite part: pairs (X,X) with two outputs and pairs (A_n, X), (X, A_n) of the covariant derivatives with one output;
log part (twoletter route-B job tl_engine.Job): the same pairs; D_1 = 4 lc^natural + Xi_Delta with the free directions
theta kept symbolic (kappa terms vanish on the two-field shell)."""
import os, sys, json, time, pickle, math
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
import ce_job as CJ
import ce_planar as PL
import ce_su11 as SU
import tl_engine as T

I = sp.I
n, nbar = SU.n, SU.nbar
NV = {0: 1, 1: I}


def letters(kind):
    if kind == 'Z':
        return dict(X=('AW', SU.ZP, SU.ZO, 0))
    zpsi, upsi, onshell, norm = SU.psi_data()
    assert onshell
    return dict(X=('lamP', zpsi, upsi, 1))


def compute(kind, Smax, with_log=True, verbose=True):
    t0 = time.time()
    LX = letters(kind)['X']
    sX, zX, uX, pX = LX
    odd = {'X': pX, 'A': 0}
    dens = {'D0': {}, 'Pb': {}, 'log': {}}
    nonrat = []
    # (X, X) -> X X
    job = CJ.CJob(sX, sX, [sX, sX], pol=[uX, uX], za=zX, zb=zX, deriv=('null', Smax, n, nbar), nullw=False)
    res = job.kernels(mlist=(0,), verbose=verbose)
    d1, r1 = PL.density(job, res['Pb'], ['X', 'X'], [pX, pX], 0, avec=n, smax=Smax)
    DIRS = (('nbar', nbar), ('e2', [0, 0, 1, 0]), ('e3', [0, 0, 0, 1]))
    for nmb, bv in DIRS:
        dX, _ = PL.density(job, res['Pb'], ['X', 'X'], [pX, pX], 0, avec=bv, smax=Smax)
        dens.setdefault('Pb_' + nmb, {})[('X', 'X')] = dX
    d2, r2 = PL.density(job, res['D0'], ['X', 'X'], [pX, pX], 0, smax=Smax)
    dens['Pb'][('X', 'X')] = d1; dens['D0'][('X', 'X')] = d2; nonrat += r1 + r2
    if verbose:
        print(kind, 'XX done', round(time.time() - t0, 1), flush=True)
    # (A_n, X) -> X and (X, A_n) -> X
    for (A, B, sa, sb, za, zb) in (('A', 'X', 'AV', sX, NV, zX), ('X', 'A', sX, 'AV', zX, NV)):
        j1 = CJ.CJob(sa, sb, [sX], pol=[uX], za=za, zb=zb, deriv=('null', max(Smax - 1, 0), n, nbar), nullw=False)
        r1_ = j1.kernels(mlist=(0,), verbose=False)
        e1, q1 = PL.density(j1, r1_['Pb'], ['X'], [pX], 0, avec=n, smax=max(Smax - 1, 0))
        for nmb, bv in (('nbar', nbar), ('e2', [0, 0, 1, 0]), ('e3', [0, 0, 0, 1])):
            eX, _ = PL.density(j1, r1_['Pb'], ['X'], [pX], 0, avec=bv, smax=max(Smax - 1, 0))
            dens.setdefault('Pb_' + nmb, {})[(A, B)] = eX
        e2, q2 = PL.density(j1, r1_['D0'], ['X'], [pX], 0, smax=max(Smax - 1, 0))
        dens['Pb'][(A, B)] = e1; dens['D0'][(A, B)] = e2; nonrat += q1 + q2
        if verbose:
            print(kind, A, B, 'done', round(time.time() - t0, 1), flush=True)
    if with_log:
        tj = T.Job(sX, sX, [sX, sX], pol=[uX, uX], za=zX, zb=zX, deriv=('null', Smax, n, nbar))
        K = tj.kernel(keep_sequences=False)
        dl, _ = PL.density(tj, K['total'], ['X', 'X'], [pX, pX], logkey='log', smax=Smax)
        dens['log'][('X', 'X')] = dl
        for (A, B, sa, sb, za, zb) in (('A', 'X', 'AV', sX, NV, zX), ('X', 'A', sX, 'AV', zX, NV)):
            tj1 = T.Job(sa, sb, [sX], pol=[uX], za=za, zb=zb, deriv=('null', max(Smax - 1, 0), n, nbar))
            K1 = tj1.kernel(keep_sequences=False)
            dl1, _ = PL.density(tj1, K1['total'], ['X'], [pX], logkey='log', smax=max(Smax - 1, 0))
            dens['log'][(A, B)] = dl1
        if verbose:
            print(kind, 'log kernels done', round(time.time() - t0, 1), flush=True)
    return dens, odd, nonrat, time.time() - t0


def three_field_words(a, b):
    """three-field part of tr((n.D)^a X (n.D)^b X): list of (coefficient, word) with sites (letter, nd)."""
    out = []
    for r in range(a):
        s_ = a - 1 - r; cf = math.comb(a, r + 1)
        out.append((cf, (('A', r), ('X', s_), ('X', b))))
        out.append((-cf, (('X', s_), ('A', r), ('X', b))))
    for r in range(b):
        s_ = b - 1 - r; cf = math.comb(b, r + 1)
        out.append((cf, (('X', a), ('A', r), ('X', s_))))
        out.append((-cf, (('X', a), ('X', s_), ('A', r))))
    return out


def two_field_image(a, b, dd, odd):
    """sum over pairs of the operator tr((n.D)^a X (n.D)^b X), projected on two-site words of X."""
    res = {}
    w2 = (('X', a), ('X', b))
    cw, cs = PL.canon(w2, odd)
    if cs == 0:
        return res, False
    for k, v in PL.apply_pairs(w2, {('X', 'X'): dd[('X', 'X')]}, odd).items():
        res[k] = res.get(k, 0) + v
    one = {('A', 'X'): dd[('A', 'X')], ('X', 'A'): dd[('X', 'A')]}
    for cf, w3 in three_field_words(a, b):
        for k, v in PL.apply_pairs(w3, one, odd).items():
            if len(k) == 2:
                res[k] = res.get(k, 0) + cf * v
    # express in terms of the canonical image of the input (sign of the input canonicalisation)
    return {k: v * cs for k, v in res.items() if v != 0}, True


def matrices(kind, dens, odd, Smax):
    thX = sp.Symbol('theta_perp') if kind == 'Z' else sp.Symbol('theta_lambda')
    cphi = sp.Rational(2, 3) if kind == 'Z' else sp.Rational(1, 3)
    wX = 1 if kind == 'Z' else sp.Rational(3, 2)
    nFX = 2 if kind == 'Z' else 4
    out = {}
    for S in range(0, Smax + 1):
        B0 = PL.cyclic_basis({'X': 2}, 2, odd, S)
        if not B0:
            out[S] = dict(dim=0)
            continue
        B1 = PL.cyclic_basis({'X': 2}, 2, odd, S + 1)

        def img(w, tab):
            r, ok = two_field_image(w[0][1], w[1][1], tab, odd)
            return r
        MD = PL.matrix(B0, B0, lambda w: img(w, dens['D0']))
        MP = PL.matrix(B0, B1, lambda w: img(w, dens['Pb']))
        other_zero = all(PL.matrix(B0, B1, lambda w: img(w, dens['Pb_' + nm_])).is_zero_matrix for nm_ in ('nbar', 'e2', 'e3') if ('Pb_' + nm_) in dens)
        TP = PL.matrix(B0, B1, lambda w: PL.deriv_map(w, odd)) * sp.Rational(-1, 3)
        TD = -sp.Rational(2, 3) * (4 + S) * sp.eye(len(B0))
        res = dict(dim=len(B0), basis=[str(w) for w in B0],
                   P7b_sumC_Pb=str(MP.tolist()), P7b_target=str(TP.tolist()), P7b_pass=(MP - TP).is_zero_matrix and other_zero,
                   P7b_other_directions_b_nbar_e2_e3_vanish=other_zero,
                   P7c_sumC_D0=str(MD.tolist()), P7c_target=str(TD.tolist()), P7c_pass=(MD - TD).is_zero_matrix)
        if dens.get('log'):
            Ml = PL.matrix(B0, B0, lambda w: img(w, dens['log']))
            lcn = Ml - 2 * nFX * sp.eye(len(B0))
            # one-body: per letter u(c_phi - 2 theta (w_X + |alpha|)) ; sum over the two letters = 2 c_phi - 2 theta (2 w_X + S)
            ob = 2 * cphi - 2 * thX * (2 * wX + S)
            Xi0 = MD + ob * sp.eye(len(B0))
            D1 = (4 * lcn + Xi0).applyfunc(sp.expand)
            D1c = D1.subs(thX, sp.Rational(-1, 3))
            res.update(lc_natural_units=str(lcn.tolist()), lc_natural_eigenvalues={str(k): v for k, v in lcn.eigenvals().items()},
                       Xi_Delta_translation_invariant=str(Xi0.applyfunc(sp.expand).tolist()),
                       D1_translation_invariant_symbolic=str(D1.tolist()),
                       D1_at_theta_m13_eigenvalues={str(k): v for k, v in D1c.eigenvals().items()})
        out[S] = res
    return out


if __name__ == '__main__':
    kind = sys.argv[1]
    Smax = int(sys.argv[2])
    dens, odd, nonrat, secs = compute(kind, Smax)
    pickle.dump(dict(dens=dens, odd=odd, nonrat=nonrat), open(os.path.join(HERE, '..', 'out.nosync', f'sl2_{kind}_S{Smax}.pkl'), 'wb'))
    print('nonrational entries:', nonrat[:5], len(nonrat))
    mats = matrices(kind, dens, odd, Smax)
    ibp = {f'{pr[0]}|{pr[1]}': PL.ibp_check(dens['D0'][pr], dens['Pb'][pr], Smax if pr == ('X', 'X') else max(Smax - 1, 0), 2 if pr == ('X', 'X') else 1)
           for pr in dens['D0']}
    print('IBP violations', {k: len(v) for k, v in ibp.items()})
    mats['ibp_violations'] = {k: v[:5] for k, v in ibp.items()}
    mats['densities_units_u'] = {k: {f'{pr[0]}|{pr[1]}': {str(a): str(b) for a, b in d.items()} for pr, d in v.items()} for k, v in dens.items()}
    for S, v in mats.items():
        if not isinstance(S, int):
            continue
        print(S, v.get('dim'), 'P7b', v.get('P7b_pass'), 'P7c', v.get('P7c_pass'), 'D1(-1/3)', v.get('D1_at_theta_m13_eigenvalues'),
              'lcn', v.get('lc_natural_eigenvalues'))
        if v.get('dim') and not (v.get('P7b_pass') and v.get('P7c_pass')):
            print('   Pb', v['P7b_sumC_Pb'], 'target', v['P7b_target'])
            print('   D0', v['P7c_sumC_D0'], 'target', v['P7c_target'])
    json.dump({str(k): v for k, v in mats.items()}, open(os.path.join(HERE, '..', 'out.nosync', f'sl2_{kind}_S{Smax}.json'), 'w'), indent=1, default=str)
    print('seconds', round(secs, 1))
