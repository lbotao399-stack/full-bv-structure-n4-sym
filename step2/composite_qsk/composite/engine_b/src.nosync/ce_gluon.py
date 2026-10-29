"""Family F3 gluon sl(2): tr((n.D)^a F11 (n.D)^b F11), two-field part on collinear on-shell gluons
p_j = t_j nbar, polarisation eps = e_2 - i e_3 (n.eps = nbar.eps = 0, transverse; light-cone gauge).

Letters (A|_V slot): E = e'.A, N = n.A, N' = (e'.d)(n.A), e' = e_2 + i e_3, n = e_0 + i e_1 (engine components);
F11 := (1/2) c^{mn} F_{mn}, c = n ^ e' (weight (1,1) component, the chirality of psi psi):
    F11 = (n.d)(e'.A) - (e'.d)(n.A) + [A_n, A_e'],     (n.d)^s F11|_2 = E_{s+1} - N'_s.
On shell (t nbar, eps): E_k -> 2 (i t)^k, N_k -> 0, N'_k -> 0, F_c -> 2 (i t)^{c+1}.
Pairs: (F,F) two outputs; (N,F), (F,N), (N,E), (E,N) one output (covariant derivatives and [A_n,A_e']); pairs with an
N spectator vanish on these legs.  Composite: route B (ce_job); log: twoletter route-B job (tl_engine.Job), both with
an extra first-order derivative factor (i e'.p) for the letter N'."""
import os, sys, json, time, pickle, math
import sympy as sp
from sympy.polys.domains import QQ_I
HERE = os.path.dirname(os.path.abspath(__file__))
import ce_job as CJ
import ce_planar as PL
import ce_su11 as SU
import tl_engine as T

I = sp.I
n, nbar = SU.n, SU.nbar
ep = [0, 0, 1, I]
EPS = {2: 1, 3: -I}
LZ = {'E': {2: 1, 3: I}, 'N': {0: 1, 1: I}, 'Np': {0: 1, 1: I}}
PRE = {'E': None, 'N': None, 'Np': ep}


def _wrap_dfac(cls):
    class W(cls):
        def __init__(self, *a, pre=None, **k):
            self._pre = pre or {}
            super().__init__(*a, **k)

        def dfac(self, p, which):
            base = super().dfac(p, which)
            v = self._pre.get(which)
            if v is None:
                return base
            R = self.R
            f = R(QQ_I(0, 1)) * sum((R(QQ_I.from_sympy(sp.sympify(v[m]))) * p[m] for m in range(4)), R.zero)
            return base * f
    return W


CJobP = _wrap_dfac(CJ.CJob)
TJobP = _wrap_dfac(T.Job)


def raw_density(job, tab, nout, which_m=0, avec=None, smax=0, logkey=None):
    """{(sites, (a,b)): value}, sites = (('G', e), ...) in trace order, value t^e raw (no i-factors), letters
    a! b! [ss^a su^b]; planar (-1)^r/2, identical field type of the outputs (A,A): 1/2."""
    nm = job.names
    iss, isu = nm.index('ss'), nm.index('su')
    it = [nm.index(f't{j}') for j in range(nout)]
    iw = [nm.index(f'w{j}') for j in range(nout) if f'w{j}' in nm]
    dens = {}
    bad = []
    ident = 2 if nout == 2 else 1
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
                es = [mon[it[j]] for j in range(nout)]
                sites = tuple(('G', es[j]) for j in order)
                key = (sites, (a, b))
                dens[key] = dens.get(key, 0) + v * sp.Rational((-1) ** len(w), 2) / ident
    return {k: sp.expand(v) for k, v in dens.items() if sp.expand(v) != 0}, bad


def shift(d, da, db):
    """relabel the letter derivative orders (a, b) -> (a + da, b + db) : E_{s+1} carries one more n-derivative."""
    out = {}
    for (sites, (a, b)), v in d.items():
        if a - da < 0 or b - db < 0:
            continue
        out[(sites, (a - da, b - db))] = out.get((sites, (a - da, b - db)), 0) + v
    return out


def lin(*terms):
    out = {}
    for c, d in terms:
        for k, v in d.items():
            out[k] = out.get(k, 0) + c * v
    return {k: sp.expand(v) for k, v in out.items() if sp.expand(v) != 0}


def compute(Smax, verbose=True):
    t0 = time.time()
    raw = {}
    for tag in ('C', 'L'):
        raw[tag] = {}
    pairs2 = [('E', 'E'), ('E', 'Np'), ('Np', 'E'), ('Np', 'Np')]
    pairs1 = [('N', 'E'), ('E', 'N'), ('N', 'Np'), ('Np', 'N')]
    for (A, B) in pairs2 + pairs1:
        nout = 2 if (A, B) in pairs2 else 1
        N = Smax + 2 if nout == 2 else Smax + 1     # one output: (E_{b+1}, N_k) with b + 1 + k <= S' + 1
        outs = ['AV'] * nout
        pol = [EPS] * nout
        job = CJobP('AV', 'AV', outs, pol=pol, za=LZ[A], zb=LZ[B], deriv=('null', N, n, nbar), nullw=False,
                    pre={'a': PRE[A], 'b': PRE[B]})
        res = job.kernels(mlist=(0,), verbose=False)
        dP, b1 = raw_density(job, res['Pb'], nout, 0, avec=n, smax=N)
        dD, b2 = raw_density(job, res['D0'], nout, 0, smax=N)
        tj = TJobP('AV', 'AV', outs, pol=pol, za=LZ[A], zb=LZ[B], deriv=('null', N, n, nbar), pre={'a': PRE[A], 'b': PRE[B]})
        K = tj.kernel(keep_sequences=False)
        dL, b3 = raw_density(tj, K['total'], nout, logkey='log', smax=N)
        raw['C'][(A, B)] = dict(Pb=dP, D0=dD, bad=b1 + b2)
        raw['L'][(A, B)] = dict(log=dL, bad=b3)
        if verbose:
            print('gluon pair', A, B, 'nout', nout, 'done', round(time.time() - t0, 1), 's', 'bad', len(b1 + b2 + b3), flush=True)
    # assemble F-letter densities
    dens = {}
    for key in ('Pb', 'D0', 'log'):
        src = raw['L'] if key == 'log' else raw['C']
        g = lambda A, B: src[(A, B)][key]
        dFF = lin((1, shift(g('E', 'E'), 1, 1)), (-1, shift(g('E', 'Np'), 1, 0)), (-1, shift(g('Np', 'E'), 0, 1)), (1, g('Np', 'Np')))
        dNF = lin((1, shift(g('N', 'E'), 0, 1)), (-1, g('N', 'Np')))
        dFN = lin((1, shift(g('E', 'N'), 1, 0)), (-1, g('Np', 'N')))
        dens[key] = {('F', 'F'): dFF, ('N', 'F'): dNF, ('F', 'N'): dFN, ('N', 'E'): g('N', 'E'), ('E', 'N'): g('E', 'N')}
    return dens, raw, time.time() - t0


def site_val(site, t):
    k, c = site
    if k == 'F':
        return 2 * (I * t) ** (c + 1)
    if k == 'E':
        return 2 * (I * t) ** c
    if k == 'N':
        return 0
    if k == 'G':
        return t ** c
    raise KeyError(site)


def word_val(word, t0, t1):
    assert len(word) == 2
    return sp.expand(site_val(word[0], t0) * site_val(word[1], t1) + site_val(word[0], t1) * site_val(word[1], t0))


def three_field(a, b):
    out = []
    for r in range(a):
        s_ = a - 1 - r; cf = math.comb(a, r + 1)
        out += [(cf, (('N', r), ('F', s_), ('F', b))), (-cf, (('F', s_), ('N', r), ('F', b)))]
    for r in range(b):
        s_ = b - 1 - r; cf = math.comb(b, r + 1)
        out += [(cf, (('F', a), ('N', r), ('F', s_))), (-cf, (('F', a), ('F', s_), ('N', r)))]
    for k in range(a + 1):
        cf = math.comb(a, k)
        out += [(cf, (('N', k), ('E', a - k), ('F', b))), (-cf, (('E', a - k), ('N', k), ('F', b)))]
    for k in range(b + 1):
        cf = math.comb(b, k)
        out += [(cf, (('F', a), ('N', k), ('E', b - k))), (-cf, (('F', a), ('E', b - k), ('N', k)))]
    return out


def image_val(a, b, dd, t0, t1):
    odd = {'F': 0, 'N': 0, 'E': 0, 'G': 0}
    tot = 0
    two = PL.apply_pairs((('F', a), ('F', b)), {('F', 'F'): dd[('F', 'F')]}, odd)
    for w, v in two.items():
        tot += v * word_val(w, t0, t1)
    one = {k: dd[k] for k in (('N', 'F'), ('F', 'N'), ('N', 'E'), ('E', 'N'))}
    for cf, w3 in three_field(a, b):
        for w, v in PL.apply_pairs(w3, one, odd).items():
            if len(w) == 2:
                tot += cf * v * word_val(w, t0, t1)
    return sp.expand(tot)


def matrices(dens, Smax):
    t0, t1 = sp.symbols('t0 t1')
    thV = sp.Symbol('theta_V')
    out = {}
    for S in range(0, Smax + 1):
        basis = [(a, S - a) for a in range(S // 2 + 1)]
        res = {}
        # values; express each image in the basis of on-shell values of tr(d^c F d^d F), c + d = S (+1 for Pb)
        def decompose(val, SS):
            bas = [(c, SS - c) for c in range(SS // 2 + 1)]
            cs = sp.symbols(f'x0:{len(bas)}')
            expr = sp.expand(val - sum(ci * word_val((('F', c), ('F', d)), t0, t1) for ci, (c, d) in zip(cs, bas)))
            eqs = sp.Poly(expr, t0, t1).coeffs() if expr != 0 else []
            sol = sp.solve(eqs, cs, dict=True) if eqs else [{}]
            if not sol:
                return None
            return [sp.nsimplify(sol[0].get(ci, 0)) for ci in cs]
        MP = []; MD = []; ML = []
        okP = okD = True
        r_trans = []
        for (a, b) in basis:
            V2 = word_val((('F', a), ('F', b)), t0, t1)
            vP = image_val(a, b, dens['Pb'], t0, t1)
            vD = image_val(a, b, dens['D0'], t0, t1)
            vL = image_val(a, b, dens['log'], t0, t1)
            trans = [x for x in (vP, vD, vL) for sym in sp.sympify(x).free_symbols if str(sym).startswith('log')]
            r_trans.append(len(trans) == 0)
            okP &= sp.expand(vP - sp.Rational(-1, 3) * I * (t0 + t1) * V2) == 0
            okD &= sp.expand(vD + sp.Rational(2, 3) * (6 + S) * V2) == 0
            MP.append(decompose(vP, S + 1)); MD.append(decompose(vD, S)); ML.append(decompose(vL, S))
        MDm = sp.Matrix(MD).T if all(x is not None for x in MD) else None
        MLm = sp.Matrix(ML).T if all(x is not None for x in ML) else None
        r = dict(dim=len(basis), basis=[f'tr(d^{a}F d^{b}F)' for a, b in basis], P7b_pass=okP, P7c_pass=okD, logs_cancel_in_images=all(r_trans),
                 sumC_Pb_in_basis=str(MP), sumC_D0_matrix=str(MDm.tolist() if MDm is not None else None))
        if MLm is not None and MDm is not None:
            lcn = MLm - 4 * sp.eye(len(basis))
            ob = sp.Rational(4, 3) - 2 * thV * (4 + S)
            Xi0 = MDm + ob * sp.eye(len(basis))
            D1 = (4 * lcn + Xi0).applyfunc(sp.expand)
            r.update(lc_natural_units=str(lcn.tolist()), lc_natural_eigenvalues={str(k): v for k, v in lcn.eigenvals().items()},
                     Xi_Delta_translation_invariant=str(Xi0.applyfunc(sp.expand).tolist()), D1_translation_invariant_symbolic=str(D1.tolist()),
                     D1_at_theta_m13_eigenvalues={str(k): v for k, v in D1.subs(thV, sp.Rational(-1, 3)).eigenvals().items()})
        out[S] = r
    return out


if __name__ == '__main__':
    Smax = int(sys.argv[1])
    dens, raw, secs = compute(Smax)
    pickle.dump(dict(dens=dens, raw=raw), open(os.path.join(HERE, '..', 'out.nosync', f'gluon_S{Smax}.pkl'), 'wb'))
    mats = matrices(dens, Smax)
    for S, v in mats.items():
        print(S, v['dim'], 'P7b', v['P7b_pass'], 'P7c', v['P7c_pass'], 'D1(-1/3)', v.get('D1_at_theta_m13_eigenvalues'), 'lcn', v.get('lc_natural_eigenvalues'))
        if not (v['P7b_pass'] and v['P7c_pass']):
            print('   ', v['sumC_Pb_in_basis'], v['sumC_D0_matrix'])
    json.dump({str(k): v for k, v in mats.items()}, open(os.path.join(HERE, '..', 'out.nosync', f'gluon_S{Smax}.json'), 'w'), indent=1, default=str)
    print('seconds', round(secs, 1))
