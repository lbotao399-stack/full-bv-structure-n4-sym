"""
verify.py -- machine verification of the classical psl(4|4)-equivariant BV action
    S_g = S + S^(1) + S^(2)          (S^(k>=3) = 0)
of the reduced 4d N=4 gauge-BV theory (classical/classical_spec.md), in the form [ESW Def 1.27].

Usage:
    python3 verify.py --sections conv,alg,ar1            (fast sections)
    python3 verify.py --sections ar2 --shard 0/2         (expensive sections are sharded; <= 2 processes)
    python3 verify.py --merge                            (writes verify.json from the partial files)

Explicit bases (10d gamma blocks, su(2), spinor components) live in model.py / algebra.py only.
All arithmetic exact (Gaussian rationals); random generators have small Gaussian-integer components
(dense: every basis component generically nonzero), seeds fixed.  A multilinear identity that holds for
random dense arguments holds identically with probability 1 (Schwartz-Zippel); no floating point anywhere.
"""
import sys, os, json, time, random, itertools, argparse, glob
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
from fractions import Fraction as Fr
import engine as E
import model as Mo
import algebra as Al
import tower as Tw

KL, WL, WA, AS = Fr(1, 4), Fr(3, 8), Fr(1, 4), Fr(-1)
CR = sp.Integer(1)
Al.set_cR(CR)

TYPES_EVEN = ['P', 'M', 'D', 'K', 'R']
TYPES_ODD = ['Q', 'S']
TYPES = TYPES_EVEN + TYPES_ODD


class Gen:
    def __init__(self, seed):
        self.rng = random.Random(seed)

    def ri(self, lo=-2, hi=2, nz=True):
        while True:
            r = self.rng.randint(lo, hi)
            if r or not nz:
                return sp.Integer(r)

    def rsp(self):
        return [sp.Integer(self.rng.randint(-2, 2)) + sp.I * self.rng.randint(-1, 1) for _ in range(16)]

    def rom(self):
        om = [[0] * 4 for _ in range(4)]
        for a in range(4):
            for b in range(a + 1, 4):
                r = self.ri()
                om[a][b], om[b][a] = r, -r
        return om

    def rR(self):
        R = {}
        for i in Mo.WP:
            for j in Mo.WP:
                if i < j:
                    r = self.ri()
                    R[(i, j)], R[(j, i)] = r, -r
        return R

    def make(self, t):
        if t == 'P':
            return Al.P([self.ri() for _ in range(4)])
        if t == 'M':
            return Al.Mrot(self.rom())
        if t == 'D':
            return Al.Dil().scale(self.ri())
        if t == 'K':
            return Al.Kb([self.ri() for _ in range(4)])
        if t == 'R':
            return Al.Rsym(self.rR())
        if t == 'Q':
            return Al.Qgen(self.rsp())
        if t == 'S':
            return Al.Sgen(self.rsp())
        raise ValueError(t)


def key_of(x):
    if isinstance(x, Al.Even):
        return ('e',) + tuple(str(c) for c in x.c)
    return ('o',) + tuple(str(c) for c in list(x.Q) + list(x.eta))


class Ctx:
    def __init__(self, kl=KL, wl=WL, wa=WA, aS=AS, s2_wperp=True, s2_on=True):
        self.V = E.Vars()
        self.S = Mo.build_S(self.V)
        self.kl, self.wl, self.wa, self.aS = kl, wl, wa, aS
        self.s2_wperp, self.s2_on = s2_wperp, s2_on
        self.c1, self.c2 = {}, {}
        self.T = Tw.Tower(self.V, self.S, {1: self.S1, 2: self.S2})

    def S1(self, args):
        (x,) = args
        k = key_of(x)
        if k not in self.c1:
            if isinstance(x, Al.Even):
                vf, R = Al.even_to_engine(self.V, x)
                self.c1[k] = Mo.S1_even(self.V, vf, R, self.kl, self.wl, wA=self.wa)
            else:
                Q, eta = Al.odd_to_engine(x)
                self.c1[k] = Mo.S1_odd(self.V, Q, eta, aS=self.aS)
        return self.c1[k]

    def S2(self, args):
        x, y = args
        if not (isinstance(x, Al.Odd) and isinstance(y, Al.Odd)) or not self.s2_on:
            return {}
        k = (key_of(x), key_of(y))
        if k not in self.c2:
            q1, e1 = Al.odd_to_engine(x)
            q2, e2 = Al.odd_to_engine(y)
            self.c2[k] = Mo.S2_odd(self.V, Mo.killing_spinor(self.V, q1, e1), Mo.killing_spinor(self.V, q2, e2),
                                   wperp=self.s2_wperp)
        return self.c2[k]

    def zero(self, P):
        return E.total_derivative_test(self.V, P)[0]


def check_MC(ctx, args):
    t = time.time()
    tot, parts = ctx.T.MC(args, split=True)
    ok, bad = E.total_derivative_test(ctx.V, tot)
    nz = {k: len(v) for k, v in parts.items() if v}
    return {'holds': ok, 'nonzero_terms_before_integration': len(tot), 'contributing_pieces': nz,
            'seconds': round(time.time() - t, 1)}


# =====================================================================================
# conventions
# =====================================================================================
def sec_conv():
    import conventions as Cv
    out = {}
    t = time.time()
    out['ESW_Def1.26_readings_equivalent_to_a_generating_function_model_arity<=3_gl(2|1)'] = Cv.run(3)
    out['seconds_readings'] = round(time.time() - t, 1)
    # ESW Lemma 3.6 item 3 fixes ESW's bracket sign: compare with classical_spec bracket
    ctx = Ctx()
    V, S = ctx.V, ctx.S
    Fd = Mo.Fields(V)
    g = Gen(101)
    G = [Fr(int(g.ri())) for _ in range(4)]
    X = {}
    for mu in range(4):
        for i in range(Mo.DIMG):
            E.padd_into(X, E.pmul(Fd.A(mu, i), Fd.cs(i)), (-G[mu], Fr(0)))
    B = E.bracket(V, S, X)
    GA = lambda j: E.padd(*[E.pscale(Fd.A(mu, j), G[mu]) for mu in range(4)])
    T = {}
    for i in range(Mo.DIMG):
        for mu in range(4):
            Dc = E.padd(Fd.c(i, Mo.e_mu(mu)), Fd.bracket_g(lambda j: Fd.A(mu, j), lambda k: Fd.c(k), i))
            E.padd_into(T, E.pmul(Dc, Fd.cs(i)), (G[mu], Fr(0)))
        for M in range(10):
            D = Fd.bracket_g(lambda j: Fd.A(M, j), GA, i)
            if M < 4:
                E.padd_into(D, E.padd(*[E.pscale(Fd.A(mu, i, Mo.e_mu(M)), G[mu]) for mu in range(4)]))
            E.padd_into(T, E.pmul(D, Fd.As(M, i)), (Fr(-1), Fr(0)))
        E.padd_into(T, E.pmul(Fd.bracket_g(GA, lambda k: Fd.c(k), i), Fd.cs(i)), (Fr(-1), Fr(0)))
        for a in range(16):
            E.padd_into(T, E.pmul(Fd.bracket_g(lambda j: Fd.lam(a, j), GA, i), Fd.lams(a, i)))
    out['ESW_Lemma3.6_item3_{S,-iota_G A c^*}_reproduced_by_classical_spec_bracket'] = ctx.zero(E.psub(B, T))
    out['(control) opposite bracket sign reproduces it'] = ctx.zero(E.padd(B, T))
    # ESW Lemma 3.6 as printed: {S,S2(Q1,Q2)} + S1(Gamma(Q1,Q2)) - {S1(Q1),S1(Q2)} = 0
    Q1, Q2 = g.make('Q'), g.make('Q')
    Bq = E.bracket(V, ctx.S1([Q1]), ctx.S1([Q2]))
    SG = ctx.S1([Al.bracket(Q1, Q2)])
    SS2 = E.bracket(V, S, ctx.S2([Q1, Q2]))
    out['ESW_Lemma3.6_printed: {S,S2}+S1(Gamma)-{S1(Q1),S1(Q2)}=0'] = ctx.zero(E.padd(SS2, SG, E.pscale(Bq, -1)))
    out['(Def1.26 printed, |f|=|g|=0) {S,S2}-S1(Gamma)-{S1(Q1),S1(Q2)}=0'] = ctx.zero(E.padd(SS2, E.pscale(SG, -1), E.pscale(Bq, -1)))
    # graded symmetry of the bracket at density level (used nowhere implicitly; tower.py computes every shuffle)
    v = g.make('M')
    Fv, Fq = ctx.S1([v]), ctx.S1([Q1])
    gs = {}
    for name, (A_, B_) in {'odd,odd': (Fv, ctx.S1([g.make('P')])), 'odd,even': (Fv, Fq), 'even,even': (Fq, ctx.S1([Q2]))}.items():
        pa, pb = E.parity(V, A_), E.parity(V, B_)
        lhs = E.bracket(V, A_, B_)
        rhs = E.pscale(E.bracket(V, B_, A_), (-1) ** (pa * pb))
        gs[name] = ctx.zero(E.psub(lhs, rhs))
    out['bracket_graded_symmetry_{F,G}=(-1)^{|F||G|}{G,F}'] = gs
    out['adopted'] = ('bracket = classical_spec (= ESW, item 3); cochain bracket = ESW Def 1.26 with |g| = degree in C^bullet; '
                      'd_CE = ESW Def 1.26 with |f| = degree in C^bullet and one extra overall sign (ESW Lemma 3.6 item 2)')
    return out


# =====================================================================================
# algebra
# =====================================================================================
def sec_alg():
    import structure as St
    out = {}
    g = Gen(202)
    cR = sp.Symbol('cR')
    Q1, Q2, eta = g.make('Q'), g.make('Q'), g.make('S')
    Al.set_cR(cR)
    # super-Jacobi, x,y odd: [x,[y,z]] = [[x,y],z] - [y,[x,z]]
    J = Al.bracket(Q1, Al.bracket(Q2, eta)) + Al.bracket(Al.bracket(Q1, Q2), eta).scale(-1) + Al.bracket(Q2, Al.bracket(Q1, eta))
    eqs = [sp.expand(c) for c in list(J.Q) + list(J.eta) if sp.expand(c) != 0]
    out['cR_from_super_Jacobi(Q,Q,S)'] = str(sp.solve(eqs, cR))
    Al.set_cR(CR)
    t = time.time()
    bad = []
    for T in itertools.combinations_with_replacement(TYPES, 3):
        x, y, z = [g.make(c) for c in T]
        px, py = Al.parity(x), Al.parity(y)
        Jx = Al.bracket(x, Al.bracket(y, z)) + Al.bracket(Al.bracket(x, y), z).scale(-1) + Al.bracket(y, Al.bracket(x, z)).scale(-(-1) ** (px * py))
        if not Jx.is_zero():
            bad.append(''.join(T))
    out['super_Jacobi_all_84_type_triples'] = {'failures': bad, 'seconds': round(time.time() - t, 1)}

    def is_ckv(v):
        d = v.div4()
        return all(sp.expand(sp.diff(v.c[n], Al.X[m]) + sp.diff(v.c[m], Al.X[n]) - (sp.Rational(1, 2) * d if m == n else 0)) == 0
                   for m in range(4) for n in range(4))
    out['[odd,odd] V-part is a conformal Killing vector'] = {a + b: is_ckv(Al.bracket(g.make(a), g.make(b))) for a, b in [('Q', 'Q'), ('Q', 'S'), ('S', 'S')]}
    out['[odd,odd] has so(W_perp) part'] = {a + b: bool(Al.bracket(g.make(a), g.make(b)).Rmatrix()) for a, b in [('Q', 'Q'), ('Q', 'S'), ('S', 'S')]}
    B = St.basis()
    C = St.Coords(B)
    out['dimensions'] = {'even': C.ne, 'odd': len(C.odd), 'subspaces': {k: len(v) for k, v in B.items()}}
    t = time.time()
    out['SW(17)_block_rule_violations'] = St.block_rule_check(B, C)
    out['seconds_block_rule'] = round(time.time() - t, 1)
    w = {}
    for n in B:
        ws = set()
        for x in B[n]:
            br = Al.bracket(Al.Dil(), x)
            cx, cb = C.coords(x), C.coords(br)
            ws |= set(str(sp.simplify(b / a)) for a, b in zip(cx, cb) if a != 0)
            if any(a == 0 and b != 0 for a, b in zip(cx, cb)):
                ws.add('not an eigenvector')
        w[n] = sorted(ws)
    out['ad_Delta_weights (Delta = 2E)'] = w
    rng = random.Random(203)

    def rcomb(names):
        el = None
        for n in names:
            for x in B[n]:
                r = sp.Integer(rng.randint(-2, 2))
                if r == 0:
                    continue
                el = x.scale(r) if el is None else el + x.scale(r)
        return el
    ev, od = ['P', 'K', 'M', 'Mbar', 'Delta', 'R'], ['Q', 'Qbar', 'Sbar', 'S']
    kf = {}
    kf['even,even'] = str(St.killing(C, rcomb(ev), rcomb(ev)))
    kf['odd,odd'] = str(St.killing(C, rcomb(od), rcomb(od)))
    x = rcomb(['R'])
    kf['R,R (same element)'] = str(St.killing(C, x, x))
    x = rcomb(['M', 'Mbar', 'Delta', 'P', 'K'])
    kf['conformal,conformal (same element)'] = str(St.killing(C, x, x))
    out['Killing_form_str(ad_x ad_y)'] = kf
    return out


# =====================================================================================
# arity 1
# =====================================================================================
def solve_coeffs(ctx, base, pieces):
    """unique c with Euler[{S, base + sum c_k piece_k}] = 0 (exact linear algebra over Q(i))."""
    V, S = ctx.V, ctx.S
    Eb = E.euler_all(V, E.bracket(V, S, base))
    Ep = [E.euler_all(V, E.bracket(V, S, p)) for p in pieces]
    cs = sp.symbols('c0:%d' % len(pieces))
    keys = set(Eb)
    for e in Ep:
        keys |= set(e)
    rows = []
    tosp = lambda z: sp.Rational(z[0].numerator, z[0].denominator) + sp.I * sp.Rational(z[1].numerator, z[1].denominator)
    for k in keys:
        mons = set(Eb.get(k, {}))
        for e in Ep:
            mons |= set(e.get(k, {}))
        for m in mons:
            expr = tosp(Eb.get(k, {}).get(m, E.ZEROC))
            for j, e in enumerate(Ep):
                expr += cs[j] * tosp(e.get(k, {}).get(m, E.ZEROC))
            rows.append(expr)
    A, b = sp.linear_eq_to_matrix(rows, cs)
    sol = sp.linsolve((A, b), cs)
    return [str(s) for s in sol], len(rows)


def sec_ar1():
    out = {}
    g = Gen(303)
    ctx = Ctx()
    V = ctx.V
    der = {}
    # rotation coefficient kl
    v = g.make('M')
    vf, R = Al.even_to_engine(V, v)
    base = Mo.S1_even(V, vf, R, Fr(0), Fr(0))
    der['kl from so(V) rotation'] = solve_coeffs(ctx, base, [E.psub(Mo.S1_even(V, vf, R, Fr(1), Fr(0)), base)])
    v = g.make('R')
    vf, R = Al.even_to_engine(V, v)
    base = Mo.S1_even(V, vf, R, Fr(0), Fr(0))
    der['kl from so(W_perp) rotation'] = solve_coeffs(ctx, base, [E.psub(Mo.S1_even(V, vf, R, Fr(1), Fr(0)), base)])
    v = Al.Dil()
    vf, R = Al.even_to_engine(V, v)
    base = Mo.S1_even(V, vf, R, Fr(0), Fr(0), wA=Fr(0))
    der['(wl, wA) from Delta'] = solve_coeffs(ctx, base, [E.psub(Mo.S1_even(V, vf, R, Fr(0), Fr(1), wA=Fr(0)), base),
                                                          E.psub(Mo.S1_even(V, vf, R, Fr(0), Fr(0), wA=Fr(1)), base)])
    v = g.make('K')
    vf, R = Al.even_to_engine(V, v)
    base = Mo.S1_even(V, vf, R, Fr(0), Fr(0), wA=Fr(0))
    der['(kl, wl, wA) from K'] = solve_coeffs(ctx, base, [E.psub(Mo.S1_even(V, vf, R, Fr(1), Fr(0), wA=Fr(0)), base),
                                                          E.psub(Mo.S1_even(V, vf, R, Fr(0), Fr(1), wA=Fr(0)), base),
                                                          E.psub(Mo.S1_even(V, vf, R, Fr(0), Fr(0), wA=Fr(1)), base)])
    eta = g.make('S')
    q, e = Al.odd_to_engine(eta)
    T1, T2, T3 = Mo.S1_odd(V, q, e, parts=True)
    der['aS from S-supercharge'] = solve_coeffs(ctx, E.padd(T1, T2), [T3])
    out['derived_coefficients'] = der
    inv = {}
    for t in TYPES:
        for trial in range(2):
            x = g.make(t)
            tt = time.time()
            inv['%s#%d' % (t, trial)] = {'holds': ctx.zero(E.bracket(V, ctx.S, ctx.S1([x]))), 'seconds': round(time.time() - tt, 1)}
    out['{S, S^(1)(x)} = 0'] = inv
    neg = {}
    for name, kw, t in [('aS=0', {'aS': Fr(0)}, 'S'), ('wl=1/4 (wrong Weyl weight of lam)', {'wl': Fr(1, 4)}, 'D'),
                        ('kl=-1/4', {'kl': Fr(-1, 4)}, 'M'), ('wA=0 (no Weyl weight on A|_Wperp)', {'wa': Fr(0)}, 'K')]:
        c2 = Ctx(**kw)
        neg[name] = c2.zero(E.bracket(c2.V, c2.S, c2.S1([g.make(t)])))
    out['negative_controls_invariance_holds'] = neg
    return out


# =====================================================================================
# arity 2, 3, 4
# =====================================================================================
def patterns_ar2():
    return [list(p) for p in itertools.combinations_with_replacement(TYPES, 2)] + [['S', 'Q'], ['K', 'P'], ['S', 'K']]


def patterns_ar3():
    pats = []
    for e in TYPES_EVEN:
        for o in [('Q', 'Q'), ('Q', 'S'), ('S', 'S')]:
            pats.append([e, o[0], o[1]])
    pats += [['Q', 'Q', 'Q'], ['Q', 'Q', 'S'], ['Q', 'S', 'S'], ['S', 'S', 'S']]
    pats += [['Q', 'K', 'S'], ['S', 'Q', 'D']]                                   # other orderings
    pats += [['P', 'K', 'Q'], ['K', 'K', 'S'], ['M', 'R', 'D']]                  # <= 1 odd argument (structural)
    return pats


def sec_arity(n, shard, nshard):
    ctx = Ctx()
    pats = patterns_ar2() if n == 2 else patterns_ar3() if n == 3 else [['Q', 'S', 'Q', 'S'], ['S', 'S', 'Q', 'Q'], ['K', 'Q', 'S', 'S']]
    out = {}
    for idx, p in enumerate(pats):
        if idx % nshard != shard:
            continue
        g = Gen(1000 * n + idx)
        args = [g.make(t) for t in p]
        out[''.join(p)] = check_MC(ctx, args)
        print(n, ''.join(p), out[''.join(p)], flush=True)
    if n == 2 and shard == 0:
        neg = {}
        g = Gen(4242)
        Q, eta, Q2 = g.make('Q'), g.make('S'), g.make('Q')
        Al.set_cR(sp.Integer(0))
        c0 = Ctx()
        neg['cR=0 (no so(W_perp) part in [Q,S]): QS closes'] = check_MC(c0, [Q, eta])['holds']
        Al.set_cR(CR)
        c1 = Ctx(s2_wperp=False)
        neg['S2 without W_perp components of Gamma: QQ closes'] = check_MC(c1, [Q, Q2])['holds']
        neg['S2 without W_perp components of Gamma: QS closes'] = check_MC(c1, [Q, eta])['holds']
        c2 = Ctx(s2_on=False)
        neg['S2 = 0: QS closes'] = check_MC(c2, [Q, eta])['holds']
        out['negative_controls'] = neg
        print('neg', neg, flush=True)
    return out


def sec_nonvacuous():
    """all-odd arity-3 patterns cancel already at density level: record that each individual term is nonzero."""
    ctx = Ctx()
    out = {}
    for idx, p in enumerate([['Q', 'Q', 'Q'], ['Q', 'Q', 'S'], ['Q', 'S', 'S'], ['S', 'S', 'S']]):
        g = Gen(3000 + idx)
        x = [g.make(t) for t in p]
        terms = {}
        for i in range(3):
            j, k = [m for m in range(3) if m != i]
            B = E.bracket(ctx.V, ctx.S1([x[i]]), ctx.S2([x[j], x[k]]))
            terms['{S1(x%d),S2(x%d,x%d)}' % (i + 1, j + 1, k + 1)] = {'density_terms': len(B), 'zero_mod_total_derivative': ctx.zero(B)}
        out[''.join(p)] = terms
        print('nv', ''.join(p), terms, flush=True)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--sections', default='')
    ap.add_argument('--shard', default='0/1')
    ap.add_argument('--merge', action='store_true')
    a = ap.parse_args()
    if a.merge:
        res = {}
        for f in sorted(glob.glob(os.path.join(HERE, 'out', 'verify_*.json'))):
            d = json.load(open(f))
            for k, v in d.items():
                res.setdefault(k, {}).update(v) if isinstance(v, dict) else res.__setitem__(k, v)
        json.dump(res, open(os.path.join(HERE, 'verify.json'), 'w'), indent=1, default=str)
        print('wrote verify.json')
        return
    sh, ns = [int(x) for x in a.shard.split('/')]
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    for s in a.sections.split(','):
        t = time.time()
        if s == 'conv':
            r = sec_conv()
        elif s == 'alg':
            r = sec_alg()
        elif s == 'ar1':
            r = sec_ar1()
        elif s == 'nv':
            r = sec_nonvacuous()
        elif s == 'forms':
            import conventions as Cv
            r = {'component_equations_from_adopted_convention (keys: parity pattern of x_1..x_n, e=even, o=odd)': Cv.component_forms()}
        elif s in ('ar2', 'ar3', 'ar4'):
            r = sec_arity(int(s[2]), sh, ns)
        else:
            continue
        r['_seconds_total'] = round(time.time() - t, 1)
        name = {'conv': 'conventions', 'alg': 'algebra', 'ar1': 'arity1_invariance', 'ar2': 'arity2_closure',
                'ar3': 'arity3_jacobi', 'ar4': 'arity4', 'nv': 'arity3_nonvacuity', 'forms': 'component_equations'}[s]
        fn = os.path.join(HERE, 'out', 'verify_%s_%d.json' % (s, sh))
        json.dump({name: r}, open(fn, 'w'), indent=1, default=str)
        print('wrote', fn, flush=True)


if __name__ == '__main__':
    main()
