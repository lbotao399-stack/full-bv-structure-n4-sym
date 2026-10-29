"""Family F3 (F11, psi): tr((n.D)^a F11 (n.D)^b psi), two-field part on collinear on-shell legs: gluon (t0 nbar, eps)
and gaugino (t1 nbar, u_psi).  Letters as in ce_gluon (E, N, N') and ce_sl2 (psi); pairs: (F,psi), (psi,F) two outputs;
(N,F), (F,N), (N,E), (E,N) -> gluon and (N,psi), (psi,N) -> psi one output (spectators E, F, psi; N spectators vanish)."""
import os, sys, json, time, pickle, math
import sympy as sp
from sympy.polys.domains import QQ_I
HERE = os.path.dirname(os.path.abspath(__file__))
import ce_job as CJ
import ce_planar as PL
import ce_su11 as SU
import ce_gluon as GL
import tl_engine as T

I = sp.I
n, nbar = SU.n, SU.nbar


def raw_density_mixed(job, tab, kinds, oddo, which_m=0, avec=None, smax=0, logkey=None):
    """as ce_gluon.raw_density, outputs of kinds (labels) e.g. ['G','P']; sites (kind, e) raw t_label^e; trace-order
    sign of odd outputs; identical field type factor only if the two kinds coincide."""
    nm = job.names
    iss, isu = nm.index('ss'), nm.index('su')
    nout = len(kinds)
    it = [nm.index(f't{j}') for j in range(nout)]
    iw = [nm.index(f'w{j}') for j in range(nout) if f'w{j}' in nm]
    ident = 2 if (nout == 2 and kinds[0] == kinds[1]) else 1
    dens = {}; bad = []
    for w, dm in tab.items():
        items = ([('1', dm.get(logkey))] if dm.get(logkey) else []) if logkey else list(dm.get(which_m, {}).items())
        for basis, P in items:
            if avec is not None:
                import ce_fam as CF
                P = CF.a_subst(job, P, avec)
            for mon, cf in P.items():
                if any(mon[i] for i in iw):
                    bad.append(('w', w)); continue
                a, b = mon[iss], mon[isu]
                if a + b > smax:
                    continue
                v = (sp.Rational(int(cf.x.numerator), int(cf.x.denominator)) + I * sp.Rational(int(cf.y.numerator), int(cf.y.denominator))) * math.factorial(a) * math.factorial(b)
                if basis != '1':
                    v = v * sp.Symbol(basis)          # transcendental part kept as a symbol (log2, log3, ...)
                order = list(reversed(w))
                sg = PL.perm_sign([j for j in order if oddo[j]])
                sites = tuple((kinds[j], mon[it[j]]) for j in order)
                key = (sites, (a, b))
                dens[key] = dens.get(key, 0) + v * sp.Rational((-1) ** len(w), 2) * sg / ident
    return {k: sp.expand(v) for k, v in dens.items() if sp.expand(v) != 0}, bad


def compute(Smax, verbose=True):
    t0 = time.time()
    zpsi, upsi, onshell, norm = SU.psi_data()
    LZ = dict(GL.LZ); PRE = dict(GL.PRE)
    LZ['psi'] = zpsi; PRE['psi'] = None
    SLOTOF = {'E': 'AV', 'N': 'AV', 'Np': 'AV', 'psi': 'lamP'}
    raw = {'C': {}, 'L': {}}
    jobs = []
    for (A, B) in [('E', 'psi'), ('Np', 'psi'), ('psi', 'E'), ('psi', 'Np')]:
        jobs.append((A, B, ['AV', 'lamP'], [GL.EPS, upsi], ['G', 'P'], [0, 1], Smax + 1))
    for (A, B) in [('N', 'E'), ('E', 'N'), ('N', 'Np'), ('Np', 'N')]:
        jobs.append((A, B, ['AV'], [GL.EPS], ['G'], [0], Smax))
    for (A, B) in [('N', 'psi'), ('psi', 'N')]:
        jobs.append((A, B, ['lamP'], [upsi], ['P'], [1], Smax))
    for (A, B, outs, pol, kinds, oddo, N) in jobs:
        job = GL.CJobP(SLOTOF[A], SLOTOF[B], outs, pol=pol, za=LZ[A], zb=LZ[B], deriv=('null', N, n, nbar), nullw=False,
                       pre={'a': PRE[A], 'b': PRE[B]})
        res = job.kernels(mlist=(0,), verbose=False)
        dP, b1 = raw_density_mixed(job, res['Pb'], kinds, oddo, 0, avec=n, smax=N)
        dD, b2 = raw_density_mixed(job, res['D0'], kinds, oddo, 0, smax=N)
        tj = GL.TJobP(SLOTOF[A], SLOTOF[B], outs, pol=pol, za=LZ[A], zb=LZ[B], deriv=('null', N, n, nbar), pre={'a': PRE[A], 'b': PRE[B]})
        K = tj.kernel(keep_sequences=False)
        dL, b3 = raw_density_mixed(tj, K['total'], kinds, oddo, logkey='log', smax=N)
        raw['C'][(A, B)] = dict(Pb=dP, D0=dD, bad=b1 + b2)
        raw['L'][(A, B)] = dict(log=dL, bad=b3)
        if verbose:
            print('mixed pair', A, B, outs, 'done', round(time.time() - t0, 1), 's bad', len(b1 + b2 + b3), flush=True)
    dens = {}
    for key in ('Pb', 'D0', 'log'):
        src = raw['L'] if key == 'log' else raw['C']
        g = lambda A, B: src[(A, B)][key]
        dens[key] = {('F', 'psi'): GL.lin((1, GL.shift(g('E', 'psi'), 1, 0)), (-1, g('Np', 'psi'))),
                     ('psi', 'F'): GL.lin((1, GL.shift(g('psi', 'E'), 0, 1)), (-1, g('psi', 'Np'))),
                     ('N', 'F'): GL.lin((1, GL.shift(g('N', 'E'), 0, 1)), (-1, g('N', 'Np'))),
                     ('F', 'N'): GL.lin((1, GL.shift(g('E', 'N'), 1, 0)), (-1, g('Np', 'N'))),
                     ('N', 'E'): g('N', 'E'), ('E', 'N'): g('E', 'N'),
                     ('N', 'psi'): g('N', 'psi'), ('psi', 'N'): g('psi', 'N')}
    return dens, raw, time.time() - t0


ODD = {'F': 0, 'N': 0, 'E': 0, 'G': 0, 'psi': 1, 'P': 1}


def site_val(site, t):
    k, c = site
    if k == 'psi':
        return (I * t) ** c
    if k == 'P':
        return t ** c
    return GL.site_val(site, t)


def word_val(word, t0, t1):
    """leg 0 gluon (t0), leg 1 gaugino (t1); one site of each type."""
    gl = [s for s in word if s[0] in ('F', 'E', 'N', 'G')]
    fe = [s for s in word if s[0] in ('psi', 'P')]
    if len(gl) != 1 or len(fe) != 1:
        return 0
    return sp.expand(site_val(gl[0], t0) * site_val(fe[0], t1))


def three_field(a, b):
    out = []
    for r in range(a):
        s_ = a - 1 - r; cf = math.comb(a, r + 1)
        out += [(cf, (('N', r), ('F', s_), ('psi', b))), (-cf, (('F', s_), ('N', r), ('psi', b)))]
    for k in range(a + 1):
        cf = math.comb(a, k)
        out += [(cf, (('N', k), ('E', a - k), ('psi', b))), (-cf, (('E', a - k), ('N', k), ('psi', b)))]
    for r in range(b):
        s_ = b - 1 - r; cf = math.comb(b, r + 1)
        out += [(cf, (('F', a), ('N', r), ('psi', s_))), (-cf, (('F', a), ('psi', s_), ('N', r)))]
    return out


def image_val(a, b, dd, t0, t1):
    tot = 0
    two = PL.apply_pairs((('F', a), ('psi', b)), {('F', 'psi'): dd[('F', 'psi')], ('psi', 'F'): dd[('psi', 'F')]}, ODD)
    for w, v in two.items():
        tot += v * word_val(w, t0, t1)
    one = {k: dd[k] for k in (('N', 'F'), ('F', 'N'), ('N', 'E'), ('E', 'N'), ('N', 'psi'), ('psi', 'N'))}
    for cf, w3 in three_field(a, b):
        for w, v in PL.apply_pairs(w3, one, ODD).items():
            if len(w) == 2:
                tot += cf * v * word_val(w, t0, t1)
    return sp.expand(tot)


def matrices(dens, Smax):
    t0, t1 = sp.symbols('t0 t1')
    thV, thl = sp.symbols('theta_V theta_lambda')
    out = {}
    for S in range(0, Smax + 1):
        basis = [(a, S - a) for a in range(S + 1)]

        def decompose(val, SS):
            bas = [(c, SS - c) for c in range(SS + 1)]
            cs = sp.symbols(f'x0:{len(bas)}')
            expr = sp.expand(val - sum(ci * word_val((('F', c), ('psi', d)), t0, t1) for ci, (c, d) in zip(cs, bas)))
            eqs = sp.Poly(expr, t0, t1).coeffs() if expr != 0 else []
            sol = sp.solve(eqs, cs, dict=True) if eqs else [{}]
            if not sol:
                return None
            return [sp.nsimplify(sol[0].get(ci, 0)) for ci in cs]
        okP = okD = True
        MP = []; MD = []; ML = []
        r_trans = []
        for (a, b) in basis:
            V2 = word_val((('F', a), ('psi', b)), t0, t1)
            vP = image_val(a, b, dens['Pb'], t0, t1)
            vD = image_val(a, b, dens['D0'], t0, t1)
            vL = image_val(a, b, dens['log'], t0, t1)
            r_trans.append(all(not str(sym).startswith('log') for x in (vP, vD, vL) for sym in sp.sympify(x).free_symbols))
            okP &= sp.expand(vP - sp.Rational(-1, 3) * I * (t0 + t1) * V2) == 0
            okD &= sp.expand(vD + sp.Rational(2, 3) * (5 + S) * V2) == 0
            MP.append(decompose(vP, S + 1)); MD.append(decompose(vD, S)); ML.append(decompose(vL, S))
        r = dict(dim=len(basis), basis=[f'tr(d^{a}F d^{b}psi)' for a, b in basis], P7b_pass=okP, P7c_pass=okD, logs_cancel_in_images=all(r_trans),
                 sumC_Pb_in_basis=str(MP), sumC_D0_columns=str(MD))
        if all(x is not None for x in MD + ML):
            MDm = sp.Matrix(MD).T; MLm = sp.Matrix(ML).T
            lcn = MLm - 6 * sp.eye(len(basis))
            obs = [sp.Rational(1, 1) - 2 * thV * (2 + a) - 2 * thl * (sp.Rational(3, 2) + b) for (a, b) in basis]
            Xi0 = MDm + sp.diag(*obs)
            D1 = (4 * lcn + Xi0).applyfunc(sp.expand)
            r.update(lc_natural_units=str(lcn.tolist()), lc_natural_eigenvalues={str(k): v for k, v in lcn.eigenvals().items()},
                     Xi_Delta_translation_invariant=str(Xi0.applyfunc(sp.expand).tolist()), D1_translation_invariant_symbolic=str(D1.tolist()),
                     D1_at_theta_m13_eigenvalues={str(k): v for k, v in D1.subs({thV: sp.Rational(-1, 3), thl: sp.Rational(-1, 3)}).eigenvals().items()})
        out[S] = r
    return out


if __name__ == '__main__':
    Smax = int(sys.argv[1])
    dens, raw, secs = compute(Smax)
    pickle.dump(dict(dens=dens, raw=raw), open(os.path.join(HERE, '..', 'out.nosync', f'mixed_S{Smax}.pkl'), 'wb'))
    mats = matrices(dens, Smax)
    for S, v in mats.items():
        print(S, v['dim'], 'P7b', v['P7b_pass'], 'P7c', v['P7c_pass'], 'D1(-1/3)', v.get('D1_at_theta_m13_eigenvalues'), 'lcn', v.get('lc_natural_eigenvalues'))
        if not (v['P7b_pass'] and v['P7c_pass']):
            print('   ', v['sumC_Pb_in_basis'], v['sumC_D0_columns'])
    json.dump({str(k): v for k, v in mats.items()}, open(os.path.join(HERE, '..', 'out.nosync', f'mixed_S{Smax}.json'), 'w'), indent=1, default=str)
    print('seconds', round(secs, 1))
