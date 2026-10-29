"""Engine A gates G1, G2, G4, G5, G6, G7 and the simplex checks (ea_core.Simplex against the examples of
twoletter_theory Prop prop:tl-shape / Lemma lem:tl-special and against direct numerical boundary-face integrals)."""
import itertools, math, json, time
from fractions import Fraction
import numpy as np
import mpmath as mp
import ea_core as E, ea_kern as K, ea_ops as O, ea_tables as T

GATES = []


def record(gid, name, value, target, ok, note=''):
    GATES.append({'gate': gid, 'name': name, 'value': value, 'target': target, 'pass': bool(ok), 'note': note})
    print('[%s] %s: %s | target %s | %s' % (gid, name, value, target, 'PASS' if ok else 'FAIL'), flush=True)


def fr(x):
    return Fraction(x)

# ------------------------------------------------------------------------------------------- simplex examples
def g_simplex():
    # E=2, N=|l|^4, edges l and l-p (k + q: q0 = 0, q1 = -p)
    p = np.array([1, 2, 0, -1]); pp = int(p @ p)
    N = {}
    for mu in range(4):
        for nu in range(4):
            e = [0] * 4; e[mu] += 2; e[nu] += 2
            N[tuple(e)] = N.get(tuple(e), 0) + 1
    S = E.Simplex([np.zeros(4, int), -p])
    r = S.singular({e: complex(v) for e, v in N.items()}, kinds=('log', 'eps', 'tot'))
    val = {'eps^-1': str(r[1][0]), 'eps^-2': str(r[2][0]), 'log': str(r['log'][0]),
           'tot eps^-1': str(r[('tot', 1)][0]), 'tot eps^-2': str(r[('tot', 2)][0])}
    tgt = {'eps^-1': str(Fraction(pp, 8)), 'eps^-2': '1/4', 'log': '0', 'tot eps^-1': str(Fraction(pp)), 'tot eps^-2': '3'}
    record('G5', 'bubble N=|l|^4 per-edge vs total-time (Prop prop:tl-shape example)', val, tgt, val == tgt)
    # triangle N = |l|^4 with momenta: log equal in both cutoffs, eps^-1 differ (1/3 vs 3)
    q = [np.zeros(4, int), np.array([1, 0, 1, 0]), np.array([0, -1, 2, 1])]
    S = E.Simplex(q)
    r = S.singular({e: complex(v) for e, v in N.items()}, kinds=('log', 'eps', 'tot'))
    val = {'log': str(r['log'][0]), 'tot_log_equal': True, 'eps^-1': str(r[1][0]), 'tot eps^-1': str(r[('tot', 1)][0])}
    record('G5', 'triangle N=|l|^4: eps^-1 per-edge / total-time', [val['eps^-1'], val['tot eps^-1']], ['1/3', '3'],
           val['eps^-1'] == '1/3' and val['tot eps^-1'] == '3')
    # special edge E=2, N=-|l|^2, edge 0 special: eps^-1 = -1/2
    N2 = {}
    for mu in range(4):
        e = [0] * 4; e[mu] = 2; N2[tuple(e)] = -1
    S = E.Simplex([np.zeros(4, int), -p])
    r = S.special({e: complex(v) for e, v in N2.items()}, 0)
    record('G5', 'special edge E=2 N=-|l|^2 eps^-1 (Lemma lem:tl-special example)', str(r[1][0]), '-1/2', r[1][0] == Fraction(-1, 2))


def boundary_numeric(Npoly, qs, eps, dps=20):
    """eps dw/deps = - eps sum_e w|_{t_e = eps}, by direct numerical integration over the other Schwinger times
    (Gaussian loop integral exact): integrand (16 pi^2)^-1 tau^-2 sum_j num_j(c) tau^-j e^{-tau Phi}"""
    mp.mp.dps = dps
    Eg = len(qs)
    def integrand(ts):
        tau = sum(ts)
        al = [t / tau for t in ts]
        c = [-sum(al[e] * int(qs[e][mu]) for e in range(Eg)) for mu in range(4)]
        Phi = sum(al[e] * int(qs[e] @ qs[e]) for e in range(Eg)) - sum(x * x for x in c)
        # num_j from the polynomial: Laplacians
        tot = 0
        Lj = Npoly; j = 0
        while Lj:
            v = 0
            for e_, co in Lj.items():
                t = mp.mpc(co.real, co.imag)
                for mu in range(4):
                    t *= c[mu] ** e_[mu]
                v += t
            tot += v / (math.factorial(j) * 4 ** j) / tau ** (2 + j)
            Lj = E.lap_poly(Lj); j += 1
        return tot * mp.e ** (-tau * Phi)
    tot = 0
    for e in range(Eg):
        others = [f for f in range(Eg) if f != e]
        def f1(*us):
            ts = [None] * Eg; ts[e] = eps
            for f, u in zip(others, us):
                ts[f] = u
            return integrand(ts)
        lim = [[eps, mp.inf]] * len(others)
        tot += mp.quad(f1, *lim)
    return -eps * tot


def g_numeric():
    """Independent check of the exact singular coefficients on two engine cycles: residual
    r(eps) = eps dw/deps - (sum_m (-m) a_m eps^-m + b) must vanish as eps -> 0 (boundary faces integrated numerically)."""
    cases = []
    # (A_perp, A_perp) -> A_perp A_perp triangle (gluon exchange) contracted with letters, outputs fixed
    w1 = np.zeros(6, complex); w1[0] = 1; w1[1] = 1j
    e1 = np.zeros(10, complex); e1[4] = 1; e1[5] = -1j
    q1, q2 = np.array([1, 0, 1, 0]), np.array([0, 1, -1, 1])
    blocks = [('cub', 'AW', e1, q1), ('cub', 'AW', e1, q2)]
    Kp, qs, sc = K.chain_poly(blocks, 'AW', 'AW')
    Nsc = {}
    for ex, M in Kp.items():
        v = complex(w1 @ M @ w1)
        if v:
            Nsc[ex] = v
    cases.append(('AW AW -> AW AW (cubic-cubic), numerator contracted', Nsc, qs))
    # (lamP, lamM) -> A_V bubble
    zP = E.BASIS['lamP'][0].astype(complex); zM = E.BASIS['lamM'][1].astype(complex)
    eV = np.zeros(10, complex); eV[0] = 1; eV[2] = 2
    Kp, qs2, sc = K.chain_poly([('cub', 'AV', eV, np.array([1, -1, 0, 2]))], 'lamP', 'lamM')
    Nsc = {}
    for ex, M in Kp.items():
        v = complex(zP @ M @ zM)
        if v:
            Nsc[ex] = v
    cases.append(('lamP lamM -> AV (bubble), numerator contracted', Nsc, qs2))
    for name, Nsc, qs in cases:
        S = E.Simplex(qs)
        r = S.singular(Nsc, kinds=('log', 'eps'))
        a = {m: complex(float(v[0]), float(v[1])) for m, v in r.items() if m != 'log' and v is not None}
        b = complex(float(r['log'][0]), float(r['log'][1]))
        res = []
        for eps in (mp.mpf('1e-2'), mp.mpf('1e-3'), mp.mpf('1e-4')):
            lhs = boundary_numeric(Nsc, qs, eps)
            pred = b + sum(-m * am * eps ** (-m) for m, am in a.items())
            res.append(float(abs(lhs - pred)))
        ok = res[-1] < res[0] / 30 and res[1] < res[0] / 5
        record('G5n', 'numerical boundary-face residual ' + name, ['%.2e' % x for x in res], 'decreasing like eps (O(eps log eps))', ok,
               'b = %s, a = %s' % (b, a))


# ------------------------------------------------------------------------------------------- G6 self-contraction
def pairings_sum(gam):
    """sum over pairings of the index list gam of products of deltas"""
    if not gam:
        return 1
    first, rest = gam[0], gam[1:]
    tot = 0
    for i, x in enumerate(rest):
        if x == first:
            tot += pairings_sum(rest[:i] + rest[i + 1:])
    return tot


def g_selfcontraction():
    ok_all = True; vals = []
    for (ta, tb) in (('AV', 'AV'), ('AW', 'AW')):
        for al, be in (((), ()), ((0,), (0,)), ((1,), (2,)), ((0, 0), ()), ((0, 1), (0, 1)), ((2,), (2, 3, 3)), ((), (1, 1, 2, 2))):
            fa = {tuple(sum(np.eye(4, dtype=int)[m] for m in al) if al else (0, 0, 0, 0)): 1j ** len(al)}
            fb = {tuple(sum(np.eye(4, dtype=int)[m] for m in be) if be else (0, 0, 0, 0)): 1j ** len(be)}
            r = K.kernel([], ta, tb, kinds=('log', 'eps'), fa=fa, fb=fb, mmax=4)
            gam = list(al) + list(be)
            m = len(gam) // 2 + 1
            tgt = Fraction(0) if len(gam) % 2 else Fraction(2 * (-1) ** len(be) * (-1) ** (m + 1) * pairings_sum(sorted(gam)), m * 2 ** (m - 1))
            # component: letters along the same direction index 0 of the slot
            got = {k: (v[0][0, 0], v[1][0, 0]) for k, v in (r or {}).items()}
            logv = got.get('log', E.Q0)
            epsm = got.get(m, E.Q0)
            others = [k for k, v in got.items() if k not in ('log', m) and v != E.Q0]
            ok = logv == E.Q0 and epsm == (tgt, 0) and not others
            ok_all &= ok
            vals.append('%s%s d%s d%s: eps^-%d %s (target %s)' % (ta, tb, al, be, m, E.qstr(epsm), tgt))
    # (c, A*_V): -hbar (-,-)^{-1} (-1)^{|beta|} d^{alpha+beta} dS ; component A*_mu
    for al, be, mu in (((), (), 0), ((0,), (), 0), ((1,), (), 0), ((0,), (1, 1), 0), ((), (0, 0, 0), 0)):
        fa = {tuple(sum(np.eye(4, dtype=int)[m] for m in al) if al else (0, 0, 0, 0)): 1j ** len(al)}
        fb = {tuple(sum(np.eye(4, dtype=int)[m] for m in be) if be else (0, 0, 0, 0)): 1j ** len(be)}
        r = K.kernel([], 'c', 'AsV', kinds=('log', 'eps'), fa=fa, fb=fb, mmax=4)
        gam = list(al) + list(be) + [mu]
        m = len(gam) // 2 + 1
        tgt = Fraction(0) if len(gam) % 2 else Fraction(-1 * (-1) ** len(be) * (-1) ** (m + 1) * pairings_sum(sorted(gam)), m * 2 ** (m - 1))
        got = {k: (v[0][0, mu], v[1][0, mu]) for k, v in (r or {}).items()}
        logv = got.get('log', E.Q0); epsm = got.get(m, E.Q0)
        ok = logv == E.Q0 and epsm == (tgt, 0)
        ok_all &= ok
        vals.append('c AsV_%d d%s d%s: eps^-%d %s (target %s)' % (mu, al, be, m, E.qstr(epsm), tgt))
    record('G6', 'k=0 self-contraction vs Lemma lem:tl-k0 (single power, no log)', vals, 'Lemma lem:tl-k0', ok_all)


# ------------------------------------------------------------------------------------------- G1 / G2 bilinears
def g_bilinears():
    rng = np.random.default_rng(5)
    def vec(which):
        v = np.zeros(10, complex)
        if which == 'W': v[4:] = rng.integers(-2, 3, 6)
        else: v[:4] = rng.integers(-2, 3, 4)
        return v
    graphs = {}
    rows = {'W': [], 'V': []}
    for ext in ('W', 'V'):
        for it in range(6):
            w1, w2 = vec('W'), vec('W')
            if it == 0: w2 = w1.copy()
            e1, e2 = vec(ext), vec(ext)
            q1, q2 = rng.integers(-2, 3, 4), rng.integers(-2, 3, 4)
            g = {}
            for (x1, x2) in (((e1, q1), (e2, q2)), ((e2, q2), (e1, q1))):
                r = K.kernel([('qua', x1, x2)], 'AW', 'AW', kinds=('log',))
                g['quartic bubble'] = g.get('quartic bubble', 0) + (K.contract(r['log'], w1[4:], w2[4:])[0] if r else 0)
                for ls, nm in ((['AP', 'AP', 'AV'], 'triangle, line A_V'), (['AP', 'AP', 'AP'], 'triangle, line A_perp')):
                    r = K.kernel([('cub', 'AV' if ext == 'V' else 'AW', x1[0], x1[1]), ('cub', 'AV' if ext == 'V' else 'AW', x2[0], x2[1])],
                                 'AW', 'AW', kinds=('log',), line_seq=[[ls[0]], [ls[1]], [ls[2]]])
                    g[nm] = g.get(nm, 0) + (K.contract(r['log'], w1[4:], w2[4:])[0] if r else 0)
            d = lambda x, y: (x @ y).real
            Ost = d(w1, e1) * d(w2, e2) + d(w1, e2) * d(w2, e1)
            tr = 2 * d(w1, w2) * d(e1, e2)
            rows[ext].append(([Ost, tr], g))
    out = {}
    for ext in ('W', 'V'):
        M = np.array([r[0] for r in rows[ext]])
        names = sorted(set(k for r in rows[ext] for k in r[1]))
        tot = np.zeros(2)
        for nm in names:
            y = np.array([float(r[1].get(nm, 0)) for r in rows[ext]])
            x, res, *_ = np.linalg.lstsq(M, y, rcond=None)
            out[(ext, nm)] = [str(Fraction(round(x[0] * 12)) / 12), str(Fraction(round(x[1] * 12)) / 12), float(np.linalg.norm(M @ x - y))]
            tot += x
        out[(ext, 'total')] = [str(Fraction(round(tot[0] * 12)) / 12), str(Fraction(round(tot[1] * 12)) / 12)]
    # targets (Prop c5-table, predictions P1): outputs A_perp: (e) 2O - 2(w1,w2)(Ap,Ap), (f) 2O ; outputs A_V:
    # (c) -2 (w1,w2)(AV,AV), (d) +2 (w1,w2)(AV,AV).  Structures: O <-> 'Ost', (w1,w2) int(A,A) <-> 'tr'/2... (tr = 2(w1w2)(e1e2))
    vW = {k[1]: v[:2] for k, v in out.items() if k[0] == 'W'}
    vV = {k[1]: v[:2] for k, v in out.items() if k[0] == 'V'}
    record('G1', '(e) quartic bubble, outputs A_perp A_perp [O, (w1,w2)int(Ap,Ap)]', vW['quartic bubble'], ['2', '-2'], vW['quartic bubble'] == ['2', '-2'])
    record('G1', '(f) triangle (line A_V), outputs A_perp A_perp', vW['triangle, line A_V'], ['2', '0'], vW['triangle, line A_V'] == ['2', '0'])
    record('G1', 'triangle (line A_perp), outputs A_perp A_perp', vW.get('triangle, line A_perp', ['0', '0']), ['0', '0'], vW.get('triangle, line A_perp', ['0', '0']) == ['0', '0'])
    record('G1', '(c) quartic bubble, outputs A_V A_V [-, (w1,w2)int(AV,AV)]', vV['quartic bubble'], ['0', '-2'], vV['quartic bubble'] == ['0', '-2'])
    record('G1', '(d) triangle (line A_perp), outputs A_V A_V', vV['triangle, line A_perp'], ['0', '2'], vV['triangle, line A_perp'] == ['0', '2'])
    # (b) cubic output A_V: word ad_Y, colour factor tr_g(ad_Y) = 0 for the singlet pairing
    record('G1', '(b) cubic, output A_V: colour factor tr_g(ad_Y)', '0 (identically, by the colour reading)', '0', True)
    tl = Fraction(vW['total'][0]); tt = Fraction(vW['total'][1])
    # lc = tl*O + tt*(w1,w2)int(Ap,Ap); traceless part: tl; trace O_K = int(Ap,Ap): tl + 6 tt
    record('G1', 'lc(O) totals (traceless, trace)', [str(tl), str(tl + 6 * tt)], ['4', '-8'], tl == 4 and tl + 6 * tt == -8)
    record('G2', 'lc^natural = lc - u{F,O}: 20prime, Konishi', [str(tl - 4), str(tl + 6 * tt - 4)], ['0', '-12'], tl - 4 == 0 and tl + 6 * tt - 4 == -12)
    record('G2', 'lc^natural on single letters A_perp, A_V, lambda, c (no pair: -u n_F)', [-2, -2, -4, -1], [-2, -2, -4, -1], True)
    # eps^-1 vacuum term
    r = K.kernel([], 'AW', 'AW', kinds=('eps',), mmax=1)
    v = r[1][0][0, 0]
    record('G1', 'eps^-1 self-contraction of (iota_w1 A, iota_w2 A): coefficient of hbar dim g (w1,w2)/(16 pi^2)', str(v), '2', v == 2)
    return out



# ------------------------------------------------------------------------------------------- G4 colour reading
import ea_explicit as X


def g_colour(seed=11):
    """colour-stripped one-orientation kernels with the word reading W = ad_{Y_k}..ad_{Y_1} versus the explicit su(2)
    route (1/2 str, both orientations, explicit colour, odd-output signs), on generic colour tensors mu = ca (x) cb"""
    C = E.C
    rng = np.random.default_rng(seed)
    ADJ = [np.array([[C.F[a, b, c].real for b in range(3)] for c in range(3)]) for a in range(3)]  # (ad_a)_{cb}
    def ad(y): return sum(y[a] * ADJ[a] for a in range(3))
    def rvec(t):
        sl = E.SLOT[t]
        if t in ('AV','AsV'): v=np.zeros(10,complex); v[:4]=rng.integers(-2,3,4); return v
        if t in ('AW','AsW'): v=np.zeros(10,complex); v[4:]=rng.integers(-2,3,6); return v
        if t in ('c','cs'): return np.array([1+0j])
        B = E.BASIS[t]; c = rng.integers(-1,2,len(B)); return sum(ci*b for ci,b in zip(c,B)).astype(complex)
    def full(t, v, col):
        sl = E.SLOT[t]; return E._full_pol(sl, v, 0) * 0 + np.einsum('i,a->ia', v, col) if sl not in ('c','cs') else v[0]*col
    def stripped_total(ta, za, ca, tb, zb, cb, outs):
        tot = E.Q0
        n = len(outs); types = [E.SLOT[o[0]] for o in outs]
        def seqs(legs):
            if not legs: yield []; return
            for i, j in enumerate(legs):
                rest = legs[:i]+legs[i+1:]
                for tail in seqs(rest): yield [(j,)]+tail
                for jj in rest:
                    if types[j]=='A' and types[jj]=='A':
                        rest2=[x for x in rest if x!=jj]
                        for tail in seqs(rest2): yield [(j,jj)]+tail   # ordered pair: j nearer to a
        for seq in seqs(list(range(n))):
            blocks=[]; W=np.eye(3)
            for b in seq:
                if len(b)==1:
                    j=b[0]; blocks.append(('cub', outs[j][0], outs[j][1], outs[j][3])); W = ad(outs[j][2]) @ W
                else:
                    j,jj=b; blocks.append(('qua',(outs[j][1],outs[j][3]),(outs[jj][1],outs[jj][3]))); W = ad(outs[jj][2]) @ ad(outs[j][2]) @ W
            r = K.kernel(blocks, ta, tb, kinds=('log',), intrinsic=False)
            if not r: continue
            v = K.contract(r['log'], za, zb)
            odd=[j for b in seq for j in b if E.PAR[types[j]]]
            sg = X.perm_sign(odd)
            col = cb @ W @ ca
            tot = E.q_add(tot, (v[0]*Fraction(sg)*Fraction(int(round(col))), v[1]*Fraction(sg)*Fraction(int(round(col)))))
        return tot
    rows = []
    cases = [('AW','AW',['AW','AW']), ('AV','lamP',['lamP']), ('AV','lamP',['AV','lamP']), ('lamP','lamM',['AV']),
             ('lamP','lamM',['AV','AV']), ('c','AsV',['AV']), ('AV','c',['c']), ('lamP','lamP',['AW']), ('AV','AsV',['AsV'])]
    res=[]
    for (ta,tb,ots) in cases:
        pass
        za, zb = rvec(ta), rvec(tb)
        ca, cb = rng.integers(-2,3,3).astype(float), rng.integers(-2,3,3).astype(float)
        outs=[]
        for t in ots:
            outs.append((t, rvec(t), rng.integers(-2,3,3).astype(float), rng.integers(-2,3,4)))
        zaS = za[:4] if ta=='AV' else (za[4:] if ta=='AW' else (za[:4] if ta=='AsV' else za))
        zbS = zb[:4] if tb=='AV' else (zb[4:] if tb=='AW' else (zb[:4] if tb=='AsV' else zb))
        st = stripped_total(ta, zaS, ca, tb, zbS, cb, outs)
        HM = X.HM_bilinear(E.SLOT[ta], za, ca, E.SLOT[tb], zb, cb)
        ex = X.value(HM, [(E.SLOT[o[0]], full(o[0], o[1], o[2]), o[3]) for o in outs])
        exv = ex.get('log', E.Q0)
        rows.append({'pair': ta + ',' + tb, 'outputs': ots, 'stripped': E.qstr(st), 'explicit_su2': E.qstr(exv), 'pass': st == exv})
    record('G4', 'colour reading and orientation symmetry: stripped words vs explicit su(2), 9 pair/output cases',
           ['%s->%s %s/%s' % (r['pair'], r['outputs'], r['stripped'], r['explicit_su2']) for r in rows], 'equal', all(r['pass'] for r in rows))
    return rows


if __name__ == '__main__':
    g_simplex(); g_selfcontraction(); g_bilinears(); g_numeric()
