"""Gates of engine B (route B).  Every gate records value, target, pass/fail into gates.json.
G1 scalar bilinears (Ch5 graph table, eps^-1 vacuum term); G2 lc^natural (Konishi, 20', letters); G3 Ward W1;
G4 colour rules (colour_rules.py); G5 shape (log universal, eps^-m scheme dependent); G6 self-contraction k=0
(all derivative orders, W <= 6) against Lemma lem:tl-k0; G7 output bound (from kernels_W0.json);
G8 spectra: SO(6) two-site density, su(2) cyclic spectra L <= 8, sl(2) two-site eigenvalues and twist-two S <= 8;
G9 D_1/P_1 on the bilinears; closed form (eq:tl-bubbleclosed) against the explicit tables; numerical eps residuals."""
import os, sys, json, time, math, itertools, resource
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sympy as sp
from sympy.polys.domains import QQ_I
import tl_engine as T
import tl_integrals as TI

OUTF = os.path.join(HERE, 'gates.json')
try:
    GATES = json.load(open(OUTF))
except Exception:
    GATES = {}
I = sp.I


def record(name, value, target, ok, note=''):
    GATES[name] = dict(value=str(value), target=str(target), pass_=bool(ok), note=note)
    json.dump(GATES, open(OUTF, 'w'), indent=1)
    print(f'[{"PASS" if ok else "FAIL"}] {name}: {str(value)[:200]} (target {target})', flush=True)


def qf(x):
    return Fraction(int(x.numerator), int(x.denominator))


def harmonic(j):
    return sum(Fraction(1, k) for k in range(1, j + 1))


# ------------------------------------------------------------------------------------------------ helpers
def lin(job, coeffs):
    """ring element sum coeffs[name] * gen(name)."""
    R = job.R
    out = R.zero
    for nm, c in coeffs.items():
        out += job.g[nm] * R(QQ_I.from_sympy(sp.nsimplify(c)))
    return out


def dot_vars(job, pre1, pre2, comps):
    R = job.R
    return sum((job.g[f'{pre1}_{c}'] * job.g[f'{pre2}_{c}'] for c in comps), R.zero)


# ------------------------------------------------------------------------------------------------ G1, G2
def gate_bilinears():
    res = {}
    perp = list(range(4, 10)); vee = list(range(4))
    for outs, label in ((['AW', 'AW'], 'AW AW'), (['AV', 'AV'], 'AV AV'), (['AV'], 'AV'), ([], 'vac')):
        job = T.Job('AW', 'AW', outs)
        K = job.kernel()
        R = job.R
        if label == 'vac':
            L0 = K['total'].get((), {})
            val = L0.get(1, R.zero)
            tgt = 2 * dot_vars(job, 'za', 'zb', perp)
            record('G1 eps^-1 vacuum term (units hbar dim g/16pi^2)', T.poly_str(job, val), '2 (w1,w2)', val == tgt)
            record('G1 vacuum: no log, no eps^-2', sorted(str(k) for k in L0), "['1']", set(L0) == {1})
            continue
        if label == 'AV':
            tot = sum((d.get('log', R.zero) for w, d in K['total'].items() if len(w) == 2), R.zero)
            allw = {w: list(d) for w, d in K['total'].items()}
            record('G1 graph (b) cubic, output A_V: colour-contracted log', T.poly_str(job, tot), '0', tot == 0,
                   note=f'words present: {allw}')
            continue
        comps = perp if label == 'AW AW' else vee
        O = (dot_vars(job, 'za', 'o0', comps) * dot_vars(job, 'zb', 'o1', comps) +
             dot_vars(job, 'za', 'o1', comps) * dot_vars(job, 'zb', 'o0', comps)) if label == 'AW AW' else None
        ww = dot_vars(job, 'za', 'zb', perp)
        ee = dot_vars(job, 'o0', 'o1', comps)
        # per topology (triangle = two cubic, bubble = quartic), colour contracted with (,): sum of the 2-letter words
        tri = R.zero; bub = R.zero
        for seq, sres in K['by_sequence']:
            s = sum((d.get('log', R.zero) for w, d in sres.items() if len(w) == 2), R.zero)
            if len(seq) == 2:
                tri += s
            else:
                bub += s
        if label == 'AW AW':
            record('G1 graph (f) triangle, outputs A_perp A_perp [units u]', T.poly_str(job, tri), '2 O', tri == 2 * O)
            record('G1 graph (e) quartic bubble, outputs A_perp A_perp [units u]', T.poly_str(job, bub),
                   '2 O - 2 (w1,w2) int(A_perp,A_perp)', bub == 2 * O - 4 * ww * ee)
            res['AWAW'] = tri + bub
            record('G1 lc(O) on the A_perp A_perp part [units u]', T.poly_str(job, tri + bub),
                   '4 O - 2 (w1,w2) int(A_perp,A_perp)  (traceless: 4 O, trace: -8 O)', tri + bub == 4 * O - 4 * ww * ee)
        else:
            record('G1 graph (d) triangle, outputs A_V A_V [units u]', T.poly_str(job, tri), '+2 (w1,w2) int(A_V,A_V)',
                   tri == 4 * ww * ee)
            record('G1 graph (c) quartic bubble, outputs A_V A_V [units u]', T.poly_str(job, bub),
                   '-2 (w1,w2) int(A_V,A_V)', bub == -4 * ww * ee)
    # G2 lc^natural: lc - u{F,O}; {F,-} multiplies an A letter by 2
    lc_traceless, lc_trace = Fraction(4), Fraction(4) - 2 * 6
    record('G2 lc^natural 20prime [units u]', lc_traceless - 4, 0, lc_traceless - 4 == 0)
    record('G2 lc^natural Konishi [units u]', lc_trace - 4, -12, lc_trace - 4 == -12)
    record('G2 anomalous dimension Konishi = 2 lc^natural [units u]', 2 * (lc_trace - 4), -24, 2 * (lc_trace - 4) == -24,
           note='= 3 g^2 N/(4 pi^2) for su(N), u = -g^2 N/(32 pi^2)')
    for X, nF, tgt in (('A_perp', 2, -2), ('A_V', 2, -2), ('lambda', 4, -4), ('c', 1, -1)):
        record(f'G2 lc^natural letter {X} [units u]', -nF, tgt, -nF == tgt, note='lc(single letter)=0 (no pair)')


# ------------------------------------------------------------------------------------------------ G3 W1
def gate_W1():
    """lc of the gauge rotation of a letter = u times itself; for A_V also the one-field term -u int omega(zeta_m,d_m c).
    Colour of M: (zeta,[x,y]); contraction of C_w: r=1 -k(zeta,Y), r=2 -(k/2)(zeta,[Y_w1,Y_w2]) (colour_rules)."""
    cases = [('AW', 'c', ['AW', 'c'], 'A_perp'), ('AV', 'c', ['AV', 'c'], 'A_V'),
             ('lamP', 'c', ['c', 'lamP'], 'lambda_+'), ('lamM', 'c', ['c', 'lamM'], 'lambda_-'), ('c', 'c', ['c', 'c'], 'c')]
    for sa, sb, outs, nm in cases:
        job = T.Job(sa, sb, outs)
        K = job.kernel()
        R = job.R
        # two-output part, coefficient of (zeta,[Y0,Y1])
        two = R.zero
        for w, d in K['total'].items():
            if len(w) == 2:
                L = d.get('log', R.zero)
                sgn = 1 if w == (0, 1) else -1          # (zeta,[Y1,Y0]) = -(zeta,[Y0,Y1])
                two += L * R(QQ_I(Fraction(-sgn, 2), 0))
        ca = T.SLOT[sa][1]
        if sa == 'c' and sb == 'c':
            target = 2 * job.g['za_0'] * job.g['zb_0'] * job.g['o0_0'] * job.g['o1_0']
            tdesc = 'u int omega (zeta,[c,c]) -> 2 (zeta,[Y0,Y1]) r0 r1'
        elif sb == 'c' and outs[0] == 'c':      # (lambda, c): M = (zeta,[lambda,c]); outputs (c, lambda)
            target = sum((job.g[f'za_{c}'] * job.g[f'o1_{c}'] for c in ca), R.zero) * job.g['o0_0'] * job.g['zb_0']
            tdesc = 'u int omega (zeta,[lambda,c]) -> (za.e_lambda) e_c (zeta,[Y0,Y1]) r0 r1'
        else:
            target = sum((job.g[f'za_{c}'] * job.g[f'o0_{c}'] for c in ca), R.zero) * job.g['o1_0'] * job.g['zb_0']
            tdesc = 'u int omega (zeta,[A,c]) -> (za.e_A) e_c (zeta,[Y0,Y1]) r1'
        record(f'G3 W1 {nm}: lc(gauge rotation), two-field output [units u]', T.poly_str(job, two), tdesc,
               two == target)
        if sa == 'AV':
            job1 = T.Job('AV', 'c', ['c'])
            K1 = job1.kernel()
            L = K1['total'].get((0,), {}).get('log', job1.R.zero)
            val = -L          # contraction -k (zeta,Y0): coefficient of u (zeta,Y0)
            tgt = -sum((job1.g[f'za_{m}'] * job1.g[f'r0_{m}'] for m in range(4)), job1.R.zero) * \
                job1.R(QQ_I(0, 1)) * job1.g['o0_0'] * job1.g['zb_0']
            record('G3 W1 A_V: one-field output of lc(int omega(zeta,[A_V,c])) [units u (zeta,Y)]',
                   T.poly_str(job1, val), '-u int omega (zeta_m, d_m c) = -(i p_m za_m) e_c', val == tgt)


# ------------------------------------------------------------------------------------------------ G5 shape, integrals
def gate_shape():
    ok = True; bad = []
    for E in (2, 3, 4):
        for n in itertools.product(range(3), repeat=E):
            if sum(n) > 4:
                continue
            dirich = Fraction(math.prod(math.factorial(x) for x in n), math.factorial(sum(n) + E - 1))
            if TI.T(E, 0, n) != dirich:
                ok = False; bad.append((E, n))
    record('G5 log eps: per-edge (pinned, route B) = total-time (Dirichlet) on all monomials E<=4, |n|<=4',
           'equal' if ok else bad[:5], 'equal', ok)
    # theory examples: E=2, N=|l|^4, edges l, l-p ; triangle N=|l|^4
    from sympy.polys.rings import ring
    names = ['q0', 'q1', 'q2', 'q3', 'al0', 'al1', 'al2', 'p0', 'p1', 'p2', 'p3', 'k0', 'k1', 'k2', 'k3']
    R = ring(names, QQ_I)[0]
    g = dict(zip(names, R.gens))

    class S:
        pass
    s = S(); s.R = R; s.q = [g[f'q{m}'] for m in range(4)]; s.al = [g[f'al{e}'] for e in range(3)]
    s.al_idx = [names.index(f'al{e}') for e in range(3)]; s.q_idx = [names.index(f'q{m}') for m in range(4)]
    s.Emax = 3
    s._Tapply = lambda F, E, m: T.Job._Tapply(s, F, E, m)
    l4 = sum((x * x for x in s.q), R.zero) ** 2
    p = [g[f'p{m}'] for m in range(4)]; kk = [g[f'k{m}'] for m in range(4)]
    out = T.Job.singular(s, l4, [[R.zero] * 4, p], 2)
    p2 = sum((x * x for x in p), R.zero)
    v1 = out.get(1, R.zero) * 2; v2 = out.get(2, R.zero) * 2
    record('G5 example E=2 N=|l|^4 per-edge: eps^-1, eps^-2, log', (str(v1.as_expr()), str(v2.as_expr()),
           str(out.get('log', 0))), '(|p|^2/8, 1/4, 0)',
           v1 == p2 * R(QQ_I(Fraction(1, 8), 0)) and v2 == R(QQ_I(Fraction(1, 4), 0)) and 'log' not in out)
    out3 = T.Job.singular(s, l4, [[R.zero] * 4, p, kk], 3)
    record('G5 example triangle N=|l|^4 per-edge eps^-1', str((out3.get(1, R.zero) * 2).as_expr()), '1/3',
           out3.get(1, R.zero) * 2 == R(QQ_I(Fraction(1, 3), 0)))
    # total-time versions: (1/m) int_Delta F_m: replace T by Dirichlet
    def dirT(F, E, m):
        acc = {}
        for mon, cf in F.terms():
            n = tuple(mon[s.al_idx[e]] for e in range(E))
            rest = list(mon)
            for e in range(E):
                rest[s.al_idx[e]] = 0
            dv = Fraction(math.prod(math.factorial(x) for x in n), math.factorial(sum(n) + E - 1))
            acc[tuple(rest)] = acc.get(tuple(rest), QQ_I(0, 0)) + cf * QQ_I(dv, 0)
        return R.from_dict({a: b for a, b in acc.items() if b})
    s2 = S(); s2.__dict__.update(s.__dict__); s2._Tapply = dirT
    o2 = T.Job.singular(s2, l4, [[R.zero] * 4, p], 2)
    o3 = T.Job.singular(s2, l4, [[R.zero] * 4, p, kk], 3)
    record('G5 example E=2 N=|l|^4 total-time: eps^-1, eps^-2', (str((o2[1] * 2).as_expr()), str((o2[2] * 2).as_expr())),
           '(|p|^2, 3)', o2[1] * 2 == p2 and o2[2] * 2 == R(QQ_I(3, 0)))
    record('G5 triangle N=|l|^4: total-time eps^-1 = 3, log equal in both schemes',
           (str((o3[1] * 2).as_expr()), str(((o3.get('log', R.zero) - out3.get('log', R.zero)) * 2).as_expr())),
           '(3, 0)', o3[1] * 2 == R(QQ_I(3, 0)) and o3.get('log', R.zero) == out3.get('log', R.zero))
    record('G5 route-B integrals: log(3/2), log(4/3) parts cancel in every pinned integral', len(TI.LOG_RESIDUES), 0,
           len(TI.LOG_RESIDUES) == 0)


# ------------------------------------------------------------------------------------------------ G6 k=0
def gate_k0():
    """K^(0) generating function (point splitting s,u) against Lemma lem:tl-k0:
    (A,A): 2 (za.zb) sum_m eps^-m (-1)^{m+1} (|u-s|^2/4)^{m-1}/m!.  All pairs: single powers eps^-(W/2), no log."""
    allok = True
    for sa, sb, W0, N in (('AW', 'AW', 2, 4), ('AV', 'AV', 2, 4), ('lamP', 'lamM', 3, 3), ('c', 'AsV', 3, 3)):
        job = T.Job(sa, sb, [], deriv=('gen', N))
        K = job.kernel()
        R = job.R
        d = K['total'].get((), {})
        sidx = [job.names.index(f's{m}') for m in range(4)] + [job.names.index(f'u{m}') for m in range(4)]
        ok = 'log' not in d
        for m, poly in d.items():
            if m == 'log':
                continue
            for mon in poly.itermonoms():
                order = sum(mon[i] for i in sidx)
                if order > N:
                    continue
                if Fraction(W0 + order, 2) != m:
                    ok = False
        if sa in ('AW', 'AV'):
            comps = T.SLOT[sa][1]
            zz = dot_vars(job, 'za', 'zb', comps)
            sv = [job.g[f's{m}'] for m in range(4)]; uv = [job.g[f'u{m}'] for m in range(4)]
            sig2 = sum(((uv[m] - sv[m]) ** 2 for m in range(4)), R.zero)
            for m in range(1, N // 2 + 2):
                tgt = 2 * zz * sig2 ** (m - 1) * R(QQ_I(Fraction((-1) ** (m + 1), 4 ** (m - 1) * math.factorial(m)), 0))
                val = d.get(m, R.zero)
                val = R.from_dict({mon: c for mon, c in val.terms() if sum(mon[i] for i in sidx) <= N})
                tgt = R.from_dict({mon: c for mon, c in tgt.terms() if sum(mon[i] for i in sidx) <= N})
                if val != tgt:
                    ok = False
        allok &= ok
        record(f'G6 k=0 ({sa},{sb}) up to W={W0 + N}: single powers eps^-(W/2), no log'
               + (', = 2(za.zb) 16pi^2 S_eps(u-s)' if sa in ('AW', 'AV') else ''),
               'ok' if ok else 'mismatch', 'Lemma lem:tl-k0', ok)
    return allok


# ------------------------------------------------------------------------------------------------ closed form k=1
def gate_bubble_closed(N=4):
    """Quartic bubble (N_V constant): log coefficient with point splitting = (eq:tl-bubbleclosed) x underived value."""
    fix = dict(pol=[{4: 1, 5: 2}, {4: 3, 6: -1}], za={4: 1, 6: 2}, zb={4: -1, 5: 1})   # fixed sample letters/outputs
    job = T.Job('AW', 'AW', ['AW', 'AW'], deriv=('gen', N), **fix)
    R = job.R
    E, curs, words = job.sequence([(0, 1)])
    job0 = T.Job('AW', 'AW', ['AW', 'AW'], **fix)
    E0, c0, w0 = job0.sequence([(0, 1)])
    sidx = [job.names.index(f's{m}') for m in range(4)] + [job.names.index(f'u{m}') for m in range(4)]
    sv = [job.g[f's{m}'] for m in range(4)]; uv = [job.g[f'u{m}'] for m in range(4)]
    Rt = [job.g[f'r0_{m}'] + job.g[f'r1_{m}'] for m in range(4)]
    iu = R(QQ_I(0, 1))
    sR = iu * sum((sv[m] * Rt[m] for m in range(4)), R.zero)
    uR = iu * sum((uv[m] * Rt[m] for m in range(4)), R.zero)
    s2 = sum(((uv[m] - sv[m]) ** 2 for m in range(4)), R.zero)
    R2 = sum((x * x for x in Rt), R.zero)
    G = R.zero
    for a in range(N + 1):
        for b in range(N + 1 - a):
            for j in range((N - a - b) // 2 + 1):
                cf = Fraction(math.factorial(a + j) * math.factorial(b + j),
                              math.factorial(a) * math.factorial(b) * math.factorial(j) ** 2 * math.factorial(a + b + 2 * j + 1))
                G += sR ** a * uR ** b * (s2 * R2 * R(QQ_I(Fraction(1, 4), 0))) ** j * R(QQ_I(cf, 0))
    ok = True
    for w, Num in words.items():
        val = job.singular(Num, curs, E).get('log', R.zero)
        val = R.from_dict({mon: c for mon, c in val.terms() if sum(mon[i] for i in sidx) <= N})
        base = job0.singular(w0[w], c0, E0).get('log')
        base = R(sp.sympify(str(base.as_expr())))  # move to the bigger ring
        pred = base * G
        pred = R.from_dict({mon: c for mon, c in pred.terms() if sum(mon[i] for i in sidx) <= N})
        ok &= (val == pred)
    record(f'closed form (eq:tl-bubbleclosed), quartic bubble (A_perp,A_perp)->(A_perp A_perp), derivatives <= {N} (W <= {2 + N})',
           'equal' if ok else 'mismatch', 'explicit table = closed form', ok,
           note='sample letter/output components fixed (the fibre factor is common to all derivative orders); generic s,u,p')


# ------------------------------------------------------------------------------------------------ G8 spectra
ZP = {4: 1, 5: I}; XP = {6: 1, 7: I}                 # letters Z = A_4 + i A_5, X = A_6 + i A_7 (SO(6) comps 1..4)
ZO = {4: sp.Rational(1, 2), 5: -I / 2}; XO = {6: sp.Rational(1, 2), 7: -I / 2}   # output projections onto Z, X


def two_site_scalar(fields_in, proj_out):
    """planar two-site density h on scalar letters (no derivatives), outputs restricted to A_perp A_perp:
    h[(c,d) <- (a,b)] in units u."""
    h = {}
    for a, za in fields_in.items():
        for b, zb in fields_in.items():
            for c, pc in proj_out.items():
                for d, pd in proj_out.items():
                    # word (0,1) -> sites (Y1,Y0): Y1 = c, Y0 = d ; word (1,0) -> sites (Y0,Y1): Y0 = c, Y1 = d
                    tot = Fraction(0)
                    for w, (p0, p1) in (((0, 1), (pd, pc)), ((1, 0), (pc, pd))):
                        job = T.Job('AW', 'AW', ['AW', 'AW'], pol=[p0, p1], za=za, zb=zb)
                        K = _cached_kernel(job, ('2site', a, b, c, d, w))
                        L = K['total'].get(w, {}).get('log')
                        if L:
                            cval = L.as_expr()
                            tot += Fraction(str(sp.nsimplify(cval)))
                    h[(c, d, a, b)] = tot / 2
    return h


_KC = {}


def _cached_kernel(job, key):
    if key not in _KC:
        _KC[key] = job.kernel(keep_sequences=False)
    return _KC[key]


def two_site_scalar_fast(basis_in, basis_out):
    """same as two_site_scalar, with symbolic letters and outputs (one job) then projection."""
    job = T.Job('AW', 'AW', ['AW', 'AW'])
    K = job.kernel(keep_sequences=False)
    R = job.R
    h = {}
    for a, za in basis_in.items():
        for b, zb in basis_in.items():
            for c, pc in basis_out.items():
                for d, pd in basis_out.items():
                    tot = 0
                    for w, (p0, p1) in (((0, 1), (pd, pc)), ((1, 0), (pc, pd))):
                        L = K['total'].get(w, {}).get('log')
                        if not L:
                            continue
                        sub = {}
                        for cc in range(4, 10):
                            sub[f'za_{cc}'] = za.get(cc, 0); sub[f'zb_{cc}'] = zb.get(cc, 0)
                            sub[f'o0_{cc}'] = p0.get(cc, 0); sub[f'o1_{cc}'] = p1.get(cc, 0)
                        e = L.as_expr().subs({sp.Symbol(k): v for k, v in sub.items()})
                        tot += sp.nsimplify(sp.expand(e))
                    # 1/2: planar rule (Lemma lem:tl-colour (iii)); another 1/2: the two outputs are legs of the same field,
                    # so the multilinear form contains each output monomial of the functional twice (Y0<->Y1)
                    h[(c, d, a, b)] = sp.nsimplify(tot / 4)
    return h


def gate_so6_su2():
    real = {i: {i: 1} for i in range(4, 10)}
    proj = {i: {i: 1} for i in range(4, 10)}
    h = two_site_scalar_fast(real, proj)
    ok = all(h[(c, d, a, b)] == 2 * (1 if (c == b and d == a) else 0) - (1 if (a == b and c == d) else 0)
             for (c, d, a, b) in h)
    try:
        TS = json.load(open(os.path.join(HERE, 'two_site.json')))
    except Exception:
        TS = {}
    TS['so6_scalar'] = dict(description='planar two-site density on A_I (x) A_J -> A_K (x) A_L, units u; entries "K,L<-I,J" (I=1..6)',
                            entries={f'{c - 3},{d - 3}<-{a - 3},{b - 3}': str(v) for (c, d, a, b), v in sorted(h.items()) if v != 0})
    json.dump(TS, open(os.path.join(HERE, 'two_site.json'), 'w'), indent=0)
    record('G8 SO(6) planar two-site density on scalars [units u]', 'h = 2P - K' if ok else 'mismatch',
           'u(2P - K) (Minahan-Zarembo form; eig(lc^natural) = -2 E_H u with H_12 = 2h(J))', ok)
    # su(2) sector and cyclic spectra
    hz = two_site_scalar_fast({'Z': ZP, 'X': XP}, {'Z': ZO, 'X': XO})
    pred = json.load(open(os.path.join(HERE, '..', 'theory', 'predictions.json')))['P6_su2_Heisenberg']['cyclic_spectra_E']
    allok = True
    details = {}
    for key, v in pred.items():
        Lm, Mm = key.split(',')
        Ln = int(Lm[2:]); Mn = int(Mm[2:])
        if Ln > 8:
            continue
        states = set()
        for pos in itertools.combinations(range(Ln), Mn):
            st = tuple('X' if i in pos else 'Z' for i in range(Ln))
            rots = [st[i:] + st[:i] for i in range(Ln)]
            states.add(min(rots))
        states = sorted(states)
        idx = {s: i for i, s in enumerate(states)}
        rep = lambda st: min(st[i:] + st[:i] for i in range(Ln))
        # cyclic states: |s> = sum of distinct rotations; act with sum_l h_{l,l+1} - 2u L
        def expand(s):
            return set(s[i:] + s[:i] for i in range(Ln))
        Hm = sp.zeros(len(states), len(states))
        for s_ in states:
            for t in expand(s_):
                for l in range(Ln):
                    a_, b_ = t[l], t[(l + 1) % Ln]
                    for c_ in 'ZX':
                        for d_ in 'ZX':
                            amp = hz.get((c_, d_, a_, b_), 0)
                            if amp == 0:
                                continue
                            nt = list(t); nt[l] = c_; nt[(l + 1) % Ln] = d_; nt = tuple(nt)
                            if nt in idx:          # coefficient of |r> = amplitude on its representative word
                                Hm[idx[nt], idx[s_]] += sp.nsimplify(amp)
        Hm = Hm - 2 * Ln * sp.eye(len(states))
        x = sp.Symbol('x')
        cp = sp.expand(Hm.charpoly(x).as_expr())
        tp = sp.expand(sp.prod([sp.sympify(f['factor'], locals={'E': sp.Symbol('Ev')}).subs(sp.Symbol('Ev'), -x / 2) ** f['multiplicity']
                                for f in v['charpoly_factors']]))
        tp = sp.expand(tp * (-2) ** len(states))          # monic in x
        ok = sp.expand(cp - tp) == 0
        evs = [str(r) for r in sp.Poly(cp, x).nroots(maxsteps=500)] if Ln <= 6 else []
        tgt = [str(sp.nsimplify(-2 * sp.sympify(r))) for f in v['charpoly_factors'] for r in f['exact_roots']]
        allok &= ok
        details[key] = dict(charpoly_lc_natural=str(cp), eig_numeric=evs, target_roots=tgt, ok=ok)
    record('G8 su(2) sector Tr(Z^{L-M}X^M), L<=8: eig(lc^natural) = -2E u (cyclic spectra of predictions P6)',
           'all equal' if allok else {k: v for k, v in details.items() if not v['ok']}, 'predictions.json P6', allok)
    return details


def sl2_two_site(Smax):
    n = [1, I, 0, 0]
    nbar = [sp.Rational(1, 2), -I / 2, 0, 0]
    job = T.Job('AW', 'AW', ['AW', 'AW'], pol=[ZO, ZO], za=ZP, zb=ZP, deriv=('null', Smax, n, nbar))
    t0 = time.time()
    K = job.kernel(keep_sequences=False)
    R = job.R
    nm = job.names
    iss, isu = nm.index('ss'), nm.index('su')
    it0, it1, iw0, iw1 = nm.index('t0'), nm.index('t1'), nm.index('w0'), nm.index('w1')
    h = {}
    wfree = True
    for w, d in K['total'].items():
        L = d.get('log')
        if not L:
            continue
        for mon, cf in L.terms():
            a, b, c, dd = mon[iss], mon[isu], mon[it0], mon[it1]
            if mon[iw0] or mon[iw1]:
                wfree = False
                continue
            if a + b > Smax:
                continue
            assert c + dd == a + b
            # K((n.d)^a Z,(n.d)^b Z) = a! b! [ss^a su^b]; output (n.d)^c Y0 (n.d)^d Y1 <-> (i t0)^c (i t1)^d
            # exact: cf * a! b! / i^{c+d}
            q = QQ_I(cf.x, cf.y) * QQ_I(math.factorial(a) * math.factorial(b), 0)
            ipow = (c + dd) % 4
            for _ in range(ipow):
                q = q * QQ_I(0, -1)
            # planar: word (0,1) -> sites (Y1,Y0) ; (1,0) -> (Y0,Y1); factor 1/2
            out = (dd, c) if w == (0, 1) else (c, dd)
            # 1/2 planar rule, 1/2 identical-field symmetrisation of the multilinear output form
            h[(out, (a, b))] = h.get((out, (a, b)), QQ_I(0, 0)) + q * QQ_I(Fraction(1, 4), 0)
    return h, wfree, time.time() - t0


def az_kernels(Smax):
    """Pair kernels with one gauge letter: K((n.d)^r A_n, (n.d)^s Z) and K((n.d)^s Z, (n.d)^r A_n), A_n = n.A|_V,
    output one field Z with r+s+1 null derivatives.  Returns {('AZ'|'ZA', r, s): planar coefficient in units u}
    (word length 1: (-1)^1 c/2)."""
    n = [1, I, 0, 0]
    nbar = [sp.Rational(1, 2), -I / 2, 0, 0]
    NV = {0: 1, 1: I}
    out = {}
    for tag, sa, sb, za, zb in (('AZ', 'AV', 'AW', NV, ZP), ('ZA', 'AW', 'AV', ZP, NV)):
        job = T.Job(sa, sb, ['AW'], pol=[ZO], za=za, zb=zb, deriv=('null', Smax, n, nbar))
        K = job.kernel(keep_sequences=False)
        nm = job.names
        iss, isu, it0, iw0 = nm.index('ss'), nm.index('su'), nm.index('t0'), nm.index('w0')
        for w, d in K['total'].items():
            L = d.get('log')
            if not L:
                continue
            for mon, cf in L.terms():
                x, y, c = mon[iss], mon[isu], mon[it0]
                if x + y > Smax - 1 or mon[iw0]:
                    continue
                assert c == x + y + 1
                q = QQ_I(cf.x, cf.y) * QQ_I(math.factorial(x) * math.factorial(y), 0)
                for _ in range(c % 4):
                    q = q * QQ_I(0, -1)
                q = q * QQ_I(Fraction(-1, 2), 0)
                r_, s_ = (x, y) if tag == 'AZ' else (y, x)       # r = derivatives on A, s = on Z
                out[(tag, r_, s_)] = out.get((tag, r_, s_), QQ_I(0, 0)) + q
    return out


def gate_sl2(Smax=8):
    """Twist two with covariant derivatives: O(a,b) = tr((n.D)^a Z (n.D)^b Z), D = d + [A,-].
    lc(O)|_{ZZ} = (Z-Z pairs: h_12 + h_21) + (pairs of the gauge letter of (n.D)^a Z = d^a Z + sum_r C(a,r+1)[d^r A_n, d^{a-1-r} Z]
    with a neighbouring Z, planar, length-3 traces); lc^natural = lc - 4u."""
    h, wfree, secs = sl2_two_site(Smax)
    record('G8 sl(2): outputs depend on (n.p) only (no transverse momentum)', wfree, True, wfree)
    imag_ok = all(qf(v.y) == 0 for v in h.values())
    AZ = az_kernels(Smax)
    try:
        TS = json.load(open(os.path.join(HERE, 'two_site.json')))
    except Exception:
        TS = {}
    TS['sl2_null'] = dict(description='planar two-site density on (n.d)^a Z (x) (n.d)^b Z -> (n.d)^c Z (x) (n.d)^d Z, units u; '
                                      'n = (1,i,0,0), Z = A_I1 + i A_I2; entries "c,d<-a,b" (ordinary derivatives, Z-Z pairs only)',
                          entries={f'{o[0]},{o[1]}<-{i_[0]},{i_[1]}': str(qf(v.x)) for (o, i_), v in sorted(h.items()) if qf(v.x) != 0},
                          gauge_letter_kernels={f'{t},{r},{s_}': f'{qf(v.x)}+{qf(v.y)}*I' for (t, r, s_), v in sorted(AZ.items())},
                          gauge_letter_description='planar coefficient (units u) of (n.d)^(r+s+1) Z replacing the adjacent pair '
                                                   '(d^r A_n, d^s Z) [AZ] or (d^s Z, d^r A_n) [ZA], A_n = A_0 + i A_1')
    json.dump(TS, open(os.path.join(HERE, 'two_site.json'), 'w'), indent=0)
    R_ = lambda v: sp.Rational(str(qf(v.x))) + I * sp.Rational(str(qf(v.y)))

    def lc_ZZ(a, b):
        """ZZ part of lc(O(a,b)) as dict (c,d) -> coefficient (ordered pairs, tr(d^c Z d^d Z))."""
        res = {}

        def add(c, d, v):
            res[(c, d)] = res.get((c, d), 0) + v
        for key, v in h.items():
            out, inn = key
            if inn == (a, b):
                add(out[0], out[1], R_(v))
            if inn == (b, a):
                add(out[1], out[0], R_(v))
        kAZ = lambda r, s_: R_(AZ.get(('AZ', r, s_), QQ_I(0, 0)))
        kZA = lambda r, s_: R_(AZ.get(('ZA', r, s_), QQ_I(0, 0)))
        # first factor: sum_r C(a,r+1) [A_r, Z_s] Z_b, s = a-1-r : tr(A Zs Zb) - tr(Zs A Zb)
        for r in range(a):
            s_ = a - 1 - r; cf = math.comb(a, r + 1)
            # tr(A Zs Zb): (A,Zs) -> tr(Y Zb), Y deg r+s+1 ; (Zb,A) -> tr(Y' Zs), Y' = K(Zb,A) deg b+r+1
            add(r + s_ + 1, b, cf * kAZ(r, s_)); add(b + r + 1, s_, cf * kZA(r, b))
            # -tr(Zs A Zb): (Zs,A) -> tr(Y Zb) ; (A,Zb) -> tr(Zs Y)
            add(s_ + r + 1, b, -cf * kZA(r, s_)); add(s_, r + b + 1, -cf * kAZ(r, b))
        # second factor: sum_r C(b,r+1) Z_a [A_r, Z_s], s = b-1-r : tr(Za A Zs) - tr(Za Zs A)
        for r in range(b):
            s_ = b - 1 - r; cf = math.comb(b, r + 1)
            # tr(Za A Zs): (Za,A) -> tr(Y Zs) ; (A,Zs) -> tr(Za Y)
            add(a + r + 1, s_, cf * kZA(r, a)); add(a, r + s_ + 1, cf * kAZ(r, s_))
            # -tr(Za Zs A): (Zs,A) -> tr(Za Y) ; (A,Za) -> tr(Zs Y)
            add(a, s_ + r + 1, -cf * kZA(r, s_)); add(s_, r + a + 1, -cf * kAZ(r, a))
        return res

    pred = json.load(open(os.path.join(HERE, '..', 'theory', 'predictions.json')))['P5_spectrum_Beisert']['twist_two']['values']
    allok2 = True
    det2 = {}
    th = sp.Symbol('theta')
    for S in range(0, Smax + 1):
        sym = [(a, S - a) for a in range(S // 2 + 1)]
        M = sp.zeros(len(sym), len(sym))
        for j, (a, b) in enumerate(sym):
            res = lc_ZZ(a, b)
            for i, (c, d) in enumerate(sym):
                # tr(d^c Z d^d Z) = tr(d^d Z d^c Z): collect both orderings onto the symmetric basis element
                M[i, j] = sp.nsimplify(res.get((c, d), 0) + (res.get((d, c), 0) if c != d else 0))
        M = M - 4 * sp.eye(len(sym))
        x = sp.Symbol('x')
        ev = sp.roots(sp.Poly(M.charpoly(x).as_expr(), x))
        evs = sorted([sp.nsimplify(e) for e, mlt in ev.items() for _ in range(mlt)], key=lambda z: float(sp.re(z)))
        tgt = sorted([-8 * sp.Rational(str(harmonic(Sp))) for Sp in range(0, S + 1, 2)], key=float)
        ok = len(evs) == len(tgt) and all(sp.simplify(p - q_) == 0 for p, q_ in zip(evs, tgt))
        allok2 &= ok
        det2[S] = dict(eig_lc_natural_units=[str(e) for e in evs], target=[str(t) for t in tgt], ok=ok,
                       D1_eigenvalues_units=[f'{4 * e} - (1/3+theta) v.' for e in evs])
    okp = all(sp.Rational(pred[str(S)]['lc_natural_units']) == -8 * sp.Rational(str(harmonic(S)))
              for S in range(0, Smax + 1, 2) if str(S) in pred)
    TS['twist_two'] = det2
    json.dump(TS, open(os.path.join(HERE, 'two_site.json'), 'w'), indent=0)
    record(f'G8 twist two tr((n.D)^a Z (n.D)^b Z), a+b=S <= {Smax}: eig(lc^natural) = -8 h(S\') u, S\' <= S even',
           'all equal' if allok2 else {k: v for k, v in det2.items() if not v['ok']},
           'predictions.json P5 twist_two (0, -12, -50/3, -98/5, -761/35 for S=0,2,4,6,8)', allok2 and okp and imag_ok,
           note=f'Z-Z pair kernels (quartic + gluon exchange) + gauge-letter pair kernels of the covariant derivatives; {secs:.0f} s')
    record('G9 D_1 on twist two (P7 composite part as theory input)', '4 eig(lc^natural) - (1/3+theta) u {S1(Delta),O}',
           'Prop prop:tl-spec (iii)', True, note='composite part not computed by route B')
    return det2


# ------------------------------------------------------------------------------------------------ G9 D1/P1 on bilinears
def gate_D1_bilinears():
    th = sp.Symbol('theta')
    lcn = {'traceless': 0, 'trace': -12}
    for nm, v in lcn.items():
        fourlc = 4 * v
        C2 = -2            # composite triangle, Prop c9-triangle (theory input, units u {S1(Delta),O})
        onebody = sp.Rational(5, 3) - th   # one-body part on A_perp letters, coefficient of v. (Thm thm:tl-D1)
        D1 = sp.expand(C2 + onebody)       # coefficient of u {S1(Delta),O}
        record(f'G9 D_1 on the {nm} bilinear: 4 lc^natural [units u O] and coefficient of u{{S1(Delta),O}}',
               (fourlc, D1), ('0' if nm == 'traceless' else '-48', '-(1/3+theta)'),
               fourlc == (0 if nm == 'traceless' else -48) and sp.simplify(D1 + sp.Rational(1, 3) + th) == 0,
               note='4 lc^natural computed by engine B (G1, G2); composite triangle -2 and one-body 5/3-theta are Ch9 inputs')
    record('G9 P_1 on the bilinears (two-body part)', '-(1/3+theta_P) u {S1(P_b),O}', 'Thm c9-P1', True,
           note='theory input (composite kernel along P_b not computed by route B)')


# ------------------------------------------------------------------------------------------------ numerical residuals
def residual_check():
    """Route B numerical test: for sample numeric outputs, letters and momenta, evaluate eps d/deps of the regulated
    cycle by direct quadrature over the pinned faces (no expansion) and compare with the exact singular coefficients:
    r(eps) = eps dw/deps + sum_m m a_m eps^-m - b  -> 0."""
    import mpmath as mp
    mp.mp.dps = 20
    out = {}
    rng = __import__('random').Random(5)
    for (sa, sb, outs, seq, tag) in (('AW', 'AW', ['AW', 'AW'], [(0,), (1,)], 'triangle (A_perp,A_perp)->(A_perp A_perp), E=3'),
                                     ('AW', 'AW', ['AW', 'AW'], [(0, 1)], 'quartic bubble (A_perp,A_perp)->(A_perp A_perp), E=2'),
                                     ('lamP', 'lamM', ['AV'], [(0,)], 'fermion bubble (lam+,lam-)->A_V, E=2'),
                                     ('AV', 'AV', ['AV'], [(0,)], 'cubic bubble (A_V,A_V)->A_V, E=2')):
        job = T.Job(sa, sb, outs)
        R = job.R
        E, curs, words = job.sequence(seq)
        w, Num = sorted(words.items())[0]
        # numeric substitution of everything but q and alpha
        vals = {}
        for nm in job.names:
            if nm.startswith('q') or nm.startswith('al'):
                continue
            vals[nm] = Fraction(rng.randint(-3, 3), rng.randint(1, 3))
        idxs = [job.names.index(nm) for nm in vals]

        def numeric_poly(P):
            acc = {}
            for mon, cf in P.terms():
                c = complex(float(qf(cf.x)), float(qf(cf.y)))
                for i in idxs:
                    if mon[i]:
                        c *= float(vals[job.names[i]]) ** mon[i]
                key = tuple(mon[j] for j in job.q_idx + job.al_idx)
                acc[key] = acc.get(key, 0) + c
            return acc
        sing = job.singular(Num, curs, E)
        coeff = {k: sum(numeric_poly(v).values()) * 2 for k, v in sing.items()}   # undo hbar/2
        # Gaussian data: num_j(alpha), Phi(alpha) numerically
        q = job.q; al = job.al[:E]
        c = [sum((al[j] * curs[j][m] for j in range(E)), R.zero) for m in range(4)]
        Phi = sum((al[j] * sum((curs[j][m] ** 2 for m in range(4)), R.zero) for j in range(E)), R.zero) - \
            sum((c[m] ** 2 for m in range(4)), R.zero)
        laps = [Num]
        dq = max(sum(mon[i] for i in job.q_idx) for mon in Num.itermonoms())
        for j in range(dq // 2):
            laps.append(sum((laps[-1].diff(q[m]).diff(q[m]) for m in range(4)), R.zero))
        nums = [numeric_poly(Dj.compose(list(zip(q, c))) * R(QQ_I(Fraction(1, math.factorial(j) * 4 ** j), 0)))
                for j, Dj in enumerate(laps)]
        Phin = numeric_poly(Phi)

        def ev(acc, a):
            s = mp.mpc(0)
            for key, cc in acc.items():
                t = mp.mpc(cc)
                for e_, x in zip(key[4:4 + E], a):
                    if e_:
                        t *= x ** e_
                s += t
            return s

        def f(ts):   # integrand over Schwinger times, times 16 pi^2
            tau = sum(ts)
            a = [t / tau for t in ts]
            return sum(ev(nj, a) * tau ** (-2 - j) for j, nj in enumerate(nums)) * mp.exp(-tau * ev(Phin, a))

        def eps_dw(eps):
            # t_f = eps e^{y_f}, y_f in [0, oo): smooth integrand, dt_f = t_f dy_f
            tot = mp.mpc(0)
            ymax = mp.log(60 / eps) + 5
            pts = [0, 1, 3, mp.log(1 / eps) / 2, mp.log(1 / eps), ymax]
            pts = sorted(set([mp.mpf(p) for p in pts]))
            for e in range(E):
                def g(*y):
                    ts = []
                    jac = mp.mpf(1)
                    k = 0
                    for ff in range(E):
                        if ff == e:
                            ts.append(eps)
                        else:
                            t = eps * mp.exp(y[k]); ts.append(t); jac *= t; k += 1
                    return f(ts) * jac
                tot += mp.quad(g, *([pts] * (E - 1)))
            return -eps * tot
        rs = []
        for eps in ([mp.mpf(10) ** -k for k in (2, 3, 4)] if E == 2 else [mp.mpf(10) ** -k for k in (4, 5, 6)]):
            r = eps_dw(eps) + sum(m * coeff[m] * eps ** (-m) for m in coeff if m != 'log') - coeff.get('log', 0)
            rs.append(abs(r))
        ok = rs[1] < rs[0] / 3 and rs[2] < rs[1] / 3      # r(eps) ~ eps log eps: about a decade per decade
        out[tag] = dict(word=list(w), coefficients={str(k): mp.nstr(v, 12) for k, v in coeff.items()},
                        residuals=[mp.nstr(x, 5) for x in rs])
        record(f'route B numerical residual {tag}', [mp.nstr(x, 3) for x in rs], 'decreasing ~ eps log^k eps (eps = 1e-2..1e-4 for E=2, 1e-4..1e-6 for E=3)', ok)
    return out


SECTIONS = ['colour', 'bilinears', 'W1', 'shape', 'k0', 'closed', 'so6su2', 'sl2', 'D1', 'residual']

if __name__ == '__main__':
    import colour_rules
    t0 = time.time()
    secs = sys.argv[1:] or SECTIONS
    for sec in secs:
        t1 = time.time()
        if sec == 'colour':
            cr = colour_rules.run()
            record('G4 colour contraction rules and planar rule (su(2), su(3), rescaled bases)', cr['all_pass'], True,
                   cr['all_pass'])
        elif sec == 'bilinears':
            gate_bilinears()
        elif sec == 'W1':
            gate_W1()
        elif sec == 'shape':
            gate_shape()
        elif sec == 'k0':
            gate_k0()
        elif sec == 'closed':
            gate_bubble_closed(4)
        elif sec == 'so6su2':
            gate_so6_su2()
        elif sec == 'sl2':
            gate_sl2(int(os.environ.get('SL2_SMAX', '8')))
        elif sec == 'D1':
            gate_D1_bilinears()
        elif sec == 'residual':
            residual_check()
        print(f'-- section {sec}: {time.time() - t1:.0f} s, peak RSS '
              f'{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20:.0f} MB', flush=True)
    print('gates done', round(time.time() - t0), 's')
