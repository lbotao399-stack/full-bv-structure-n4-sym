"""Family F2 (two gaugino components): tr((n.D)^a psi (n.D)^b psi') and tr((n.D)^a Z (n.D)^b F11), psi, psi' the lambda_+
letters of R-weights (1/2,1/2,-1/2), (1/2,-1/2,1/2) (sum = weight of Z = A_4 + i A_5) and Lorentz weight (1/2,1/2);
the family closes on tr(Z F11).  Two-field part on collinear on-shell legs: sector A (psi at t0, psi' at t1),
sector B (Z at t0, gluon at t1).  Cross-sector mixing through the two-output kernels (psi,psi') -> (Z, A) and
(Z, F) -> (psi, psi')."""
import os, sys, json, time, pickle, math
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
import ce_job as CJ
import ce_planar as PL
import ce_su11 as SU
import ce_gluon as GL
import ce_mixed as MX
import ce_letters as LT

I = sp.I
n, nbar = SU.n, SU.nbar


def psi_pair():
    lets = LT.lam_letters()
    Zm = sp.Matrix.hstack(*[v for v, ch, c in lets]).T
    U = Zm.inv()
    out = {}
    for nm, tgt in (('p', (sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(1, 2))),
                    ('q', (sp.Rational(1, 2), sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2)))):
        k = [i for i, (v, ch, c) in enumerate(lets) if ch == tgt][0]
        z = lets[k][0]; u = U[:, k]
        rho = sum((nbar[m] * LT.R10[m] for m in range(4)), sp.zeros(16, 16))
        assert (rho * u).applyfunc(sp.expand).is_zero_matrix
        out[nm] = ({i: z[i] for i in range(16) if z[i] != 0}, {i: u[i] for i in range(16) if u[i] != 0})
    return out


def compute(Smax, verbose=True):
    t0 = time.time()
    pq = psi_pair()
    LZ = dict(GL.LZ); PRE = dict(GL.PRE)
    LZ['p'] = pq['p'][0]; LZ['q'] = pq['q'][0]; LZ['Z'] = SU.ZP
    for x in ('p', 'q', 'Z'):
        PRE[x] = None
    SL = {'E': 'AV', 'N': 'AV', 'Np': 'AV', 'p': 'lamP', 'q': 'lamP', 'Z': 'AW'}
    OUTS = {'P': ('lamP', pq['p'][1], 1), 'Q': ('lamP', pq['q'][1], 1), 'Zo': ('AW', SU.ZO, 0), 'G': ('AV', GL.EPS, 0)}
    jobs = []
    for (A, B) in [('p', 'q'), ('q', 'p')]:
        jobs.append((A, B, ['P', 'Q'], Smax))
        jobs.append((A, B, ['Zo', 'G'], Smax))
    for (A, B) in [('Z', 'E'), ('Z', 'Np'), ('E', 'Z'), ('Np', 'Z')]:
        jobs.append((A, B, ['P', 'Q'], Smax + 1))
        jobs.append((A, B, ['Zo', 'G'], Smax + 1))
    for (A, B, o) in [('N', 'p', 'P'), ('p', 'N', 'P'), ('N', 'q', 'Q'), ('q', 'N', 'Q'), ('N', 'Z', 'Zo'), ('Z', 'N', 'Zo'),
                      ('N', 'E', 'G'), ('E', 'N', 'G'), ('N', 'Np', 'G'), ('Np', 'N', 'G')]:
        jobs.append((A, B, [o], Smax + 1))
    raw = {'C': {}, 'L': {}}
    for (A, B, kinds, N) in jobs:
        outs = [OUTS[k][0] for k in kinds]; pol = [OUTS[k][1] for k in kinds]; oddo = [OUTS[k][2] for k in kinds]
        job = GL.CJobP(SL[A], SL[B], outs, pol=pol, za=LZ[A], zb=LZ[B], deriv=('null', N, n, nbar), nullw=False,
                       pre={'a': PRE[A], 'b': PRE[B]})
        res = job.kernels(mlist=(0,), verbose=False)
        dP, b1 = MX.raw_density_mixed(job, res['Pb'], kinds, oddo, 0, avec=n, smax=N)
        dD, b2 = MX.raw_density_mixed(job, res['D0'], kinds, oddo, 0, smax=N)
        tj = GL.TJobP(SL[A], SL[B], outs, pol=pol, za=LZ[A], zb=LZ[B], deriv=('null', N, n, nbar), pre={'a': PRE[A], 'b': PRE[B]})
        K = tj.kernel(keep_sequences=False)
        dL, b3 = MX.raw_density_mixed(tj, K['total'], kinds, oddo, logkey='log', smax=N)
        key = (A, B, tuple(kinds))
        raw['C'][key] = dict(Pb=dP, D0=dD); raw['L'][key] = dict(log=dL)
        if verbose:
            print('psiZF', key, 'done', round(time.time() - t0, 1), 's', flush=True)
    dens = {}
    for tab in ('Pb', 'D0', 'log'):
        src = raw['L'] if tab == 'log' else raw['C']
        g = lambda A, B, ks: src.get((A, B, ks), {}).get(tab, {})
        d = {}
        for ks in (('P', 'Q'), ('Zo', 'G')):
            d[('p', 'q', ks)] = g('p', 'q', ks); d[('q', 'p', ks)] = g('q', 'p', ks)
            d[('Z', 'F', ks)] = GL.lin((1, GL.shift(g('Z', 'E', ks), 0, 1)), (-1, g('Z', 'Np', ks)))
            d[('F', 'Z', ks)] = GL.lin((1, GL.shift(g('E', 'Z', ks), 1, 0)), (-1, g('Np', 'Z', ks)))
        for (A, B, o) in [('N', 'p', 'P'), ('p', 'N', 'P'), ('N', 'q', 'Q'), ('q', 'N', 'Q'), ('N', 'Z', 'Zo'), ('Z', 'N', 'Zo'),
                          ('N', 'E', 'G'), ('E', 'N', 'G')]:
            d[(A, B, (o,))] = g(A, B, (o,))
        d[('N', 'F', ('G',))] = GL.lin((1, GL.shift(g('N', 'E', ('G',)), 0, 1)), (-1, g('N', 'Np', ('G',))))
        d[('F', 'N', ('G',))] = GL.lin((1, GL.shift(g('E', 'N', ('G',)), 1, 0)), (-1, g('Np', 'N', ('G',))))
        dens[tab] = d
    return dens, raw, time.time() - t0


ODD = {'p': 1, 'q': 1, 'P': 1, 'Q': 1, 'Z': 0, 'Zo': 0, 'F': 0, 'N': 0, 'E': 0, 'G': 0}


def sval(site, t):
    k, c = site
    if k in ('p', 'q', 'Z'):
        return (I * t) ** c
    if k in ('P', 'Q', 'Zo', 'G'):
        return t ** c
    return GL.site_val(site, t)


def word_val(word, t0, t1, sector):
    """sector A: legs psi (t0), psi' (t1); B: Z (t0), gluon (t1).  Labelled forms; graded for fermions."""
    if len(word) != 2:
        return 0
    if sector == 'A':
        kind0 = ('p', 'P'); kind1 = ('q', 'Q')
    else:
        kind0 = ('Z', 'Zo'); kind1 = ('F', 'E', 'N', 'G')
    s0, s1 = word
    if s0[0] in kind0 and s1[0] in kind1:
        return sp.expand(sval(s0, t0) * sval(s1, t1))
    if s1[0] in kind0 and s0[0] in kind1:
        # word (X1, X0) with X0 on leg 0: labelled form coefficient of r0 r1 of tr(S1 S0): fermions give -1
        sg = -1 if (ODD[s0[0]] and ODD[s1[0]]) else 1
        return sp.expand(sg * sval(s1, t0) * sval(s0, t1))
    return 0


def dens_for(dd, pairs):
    """restrict to letter pairs; apply_pairs keys densities by (letterA, letterB): merge output kinds."""
    out = {}
    for (A, B, ks), d in dd.items():
        if (A, B) in pairs:
            m = out.setdefault((A, B), {})
            for k, v in d.items():
                m[k] = m.get(k, 0) + v
    return out


def three_field_A(a, b):
    out = []
    for r in range(a):
        s_ = a - 1 - r; cf = math.comb(a, r + 1)
        out += [(cf, (('N', r), ('p', s_), ('q', b))), (-cf, (('p', s_), ('N', r), ('q', b)))]
    for r in range(b):
        s_ = b - 1 - r; cf = math.comb(b, r + 1)
        out += [(cf, (('p', a), ('N', r), ('q', s_))), (-cf, (('p', a), ('q', s_), ('N', r)))]
    return out


def three_field_B(a, b):
    out = []
    for r in range(a):
        s_ = a - 1 - r; cf = math.comb(a, r + 1)
        out += [(cf, (('N', r), ('Z', s_), ('F', b))), (-cf, (('Z', s_), ('N', r), ('F', b)))]
    for r in range(b):
        s_ = b - 1 - r; cf = math.comb(b, r + 1)
        out += [(cf, (('Z', a), ('N', r), ('F', s_))), (-cf, (('Z', a), ('F', s_), ('N', r)))]
    for k in range(b + 1):
        cf = math.comb(b, k)
        out += [(cf, (('Z', a), ('N', k), ('E', b - k))), (-cf, (('Z', a), ('E', b - k), ('N', k)))]
    return out


def image(sector_in, a, b, dd, t0, t1):
    """values of the two-field image in sectors A and B."""
    two = dens_for(dd, [('p', 'q'), ('q', 'p'), ('Z', 'F'), ('F', 'Z')])
    one = dens_for(dd, [('N', 'p'), ('p', 'N'), ('N', 'q'), ('q', 'N'), ('N', 'Z'), ('Z', 'N'), ('N', 'F'), ('F', 'N'), ('N', 'E'), ('E', 'N')])
    w2 = (('p', a), ('q', b)) if sector_in == 'A' else (('Z', a), ('F', b))
    tf = three_field_A(a, b) if sector_in == 'A' else three_field_B(a, b)
    vals = {'A': 0, 'B': 0}
    for w, v in PL.apply_pairs(w2, two, ODD).items():
        for sec in ('A', 'B'):
            vals[sec] += v * word_val(w, t0, t1, sec)
    for cf, w3 in tf:
        for w, v in PL.apply_pairs(w3, one, ODD).items():
            if len(w) == 2:
                for sec in ('A', 'B'):
                    vals[sec] += cf * v * word_val(w, t0, t1, sec)
    return {k: sp.expand(v) for k, v in vals.items()}


def matrices(dens, Smax):
    t0, t1 = sp.symbols('t0 t1')
    thl, thp, thV = sp.symbols('theta_lambda theta_perp theta_V')
    out = {}
    for S in range(0, Smax + 1):
        basis = [('A', a, S - a) for a in range(S + 1)] + [('B', a, S - a) for a in range(S + 1)]
        V = {('A', c, d): word_val((('p', c), ('q', d)), t0, t1, 'A') for c in range(S + 2) for d in range(S + 2)}
        V.update({('B', c, d): word_val((('Z', c), ('F', d)), t0, t1, 'B') for c in range(S + 2) for d in range(S + 2)})

        def decompose(vals, SS):
            bas = [('A', c, SS - c) for c in range(SS + 1)] + [('B', c, SS - c) for c in range(SS + 1)]
            cs = sp.symbols(f'x0:{len(bas)}')
            res = []
            sol_all = {}
            eqs = []
            for sec in ('A', 'B'):
                expr = sp.expand(vals[sec] - sum(ci * V[b_] for ci, b_ in zip(cs, bas) if b_[0] == sec))
                if expr != 0:
                    eqs += sp.Poly(expr, t0, t1).coeffs()
            sol = sp.solve(eqs, cs, dict=True) if eqs else [{}]
            if not sol:
                return None
            return [sp.nsimplify(sol[0].get(ci, 0)) for ci in cs]
        cols = {'Pb': [], 'D0': [], 'log': []}
        okP = okD = True
        for (sec, a, b) in basis:
            for tab in cols:
                img = image(sec, a, b, dens[tab], t0, t1)
                if tab == 'Pb':
                    tgt = {'A': 0, 'B': 0}; tgt[sec] = sp.Rational(-1, 3) * I * (t0 + t1) * V[(sec, a, b)]
                    okP &= all(sp.expand(img[x] - tgt[x]) == 0 for x in ('A', 'B'))
                if tab == 'D0':
                    Nn = (4 + S) if sec == 'A' else (5 + S)
                    tgt = {'A': 0, 'B': 0}; tgt[sec] = -sp.Rational(2, 3) * Nn * V[(sec, a, b)]
                    okD &= all(sp.expand(img[x] - tgt[x]) == 0 for x in ('A', 'B'))
                cols[tab].append(decompose(img, S + 1 if tab == 'Pb' else S))
        r = dict(dim=len(basis), basis=[f'{"tr(d^%d psi d^%d psi\')" % (a, b) if s == "A" else "tr(d^%d Z d^%d F11)" % (a, b)}' for s, a, b in basis],
                 P7b_pass=okP, P7c_pass=okD, sumC_D0_columns=str(cols['D0']), sumC_Pb_columns=str(cols['Pb']))
        if all(x is not None for x in cols['D0'] + cols['log']):
            MD = sp.Matrix(cols['D0']).T; ML = sp.Matrix(cols['log']).T
            nF = [8 if s == 'A' else 4 for s, a, b in basis]
            lcn = ML - sp.diag(*nF)
            ob = []
            for s, a, b in basis:
                if s == 'A':
                    ob.append(sp.Rational(2, 3) - 2 * thl * (3 + a + b))
                else:
                    ob.append(sp.Rational(4, 3) - 2 * thp * (1 + a) - 2 * thV * (2 + b))
            Xi0 = MD + sp.diag(*ob)
            D1 = (4 * lcn + Xi0).applyfunc(sp.expand)
            D1c = D1.subs({thl: sp.Rational(-1, 3), thp: sp.Rational(-1, 3), thV: sp.Rational(-1, 3)})
            r.update(lc_natural_units=str(lcn.tolist()), lc_natural_eigenvalues={str(k): v for k, v in lcn.eigenvals().items()},
                     Xi_Delta_translation_invariant=str(Xi0.applyfunc(sp.expand).tolist()), D1_translation_invariant_symbolic=str(D1.tolist()),
                     D1_at_theta_m13_eigenvalues={str(k): v for k, v in D1c.eigenvals().items()})
        out[S] = r
    return out


if __name__ == '__main__':
    Smax = int(sys.argv[1])
    dens, raw, secs = compute(Smax)
    pickle.dump(dict(dens=dens, raw=raw), open(os.path.join(HERE, '..', 'out.nosync', f'psiZF_S{Smax}.pkl'), 'wb'))
    mats = matrices(dens, Smax)
    for S, v in mats.items():
        print(S, v['dim'], 'P7b', v['P7b_pass'], 'P7c', v['P7c_pass'], 'D1(-1/3)', v.get('D1_at_theta_m13_eigenvalues'), 'lcn', v.get('lc_natural_eigenvalues'))
        if not (v['P7b_pass'] and v['P7c_pass']):
            print('   D0', v['sumC_D0_columns'])
    json.dump({str(k): v for k, v in mats.items()}, open(os.path.join(HERE, '..', 'out.nosync', f'psiZF_S{Smax}.json'), 'w'), indent=1, default=str)
    print('seconds', round(secs, 1))
