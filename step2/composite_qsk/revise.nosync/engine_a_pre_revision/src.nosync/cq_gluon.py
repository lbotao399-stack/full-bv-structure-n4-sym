"""Engine A, part C, family F3 (gluon sl(2)): tr(D_n^k F D_n^{N-k} F), F = n^mu m^nu F_{mu nu} = d_n A_m - d_m A_n + [A_n, A_m]
(n = (1,-i,0,0), m = (0,0,1,i): the two null vectors sharing the undotted index of psi, so F is the component F_11 of
the self-dual field strength up to normalisation), D_n = d_n + [A_n, -].  Expanded to three fields (gauge completion:
the commutator in F and the covariant derivatives).  On-shell gluon legs (null Gaussian q, polarisation e with e.q = 0).
Two-field part in the basis B_j = tr(d^j f d^{N-j} f), f = d_n A_m - d_m A_n, j <= N/2.
Composite kernels as letter part + kappa * vertex-vertex part (special edge between two vertices of I).
Ward gate: the two-field part of sum C^(2)_x(O) with one leg longitudinal (e_1 = q_1) must vanish on shell
(Q-closedness of the quadratic part, the A*|_V c outputs give the free equation of motion on the other leg)."""
import itertools, math, time
from fractions import Fraction
import numpy as np
import sympy as sp
import cq_gops as GO, cq_gval as GV, cq_eval as V, cq_su11 as S11, cq_fsl2 as F2

NV = np.array([1, -1j, 0, 0]); MV = np.array([0, 0, 1, 1j])
LOGSYM = S11.LOGSYM


def poly_nm(j, l):
    """(i n.p)^j (i m.p)^l"""
    a = F2.nullpoly(NV, j); b = F2.nullpoly(MV, l)
    out = {}
    for e1, c1 in a.items():
        for e2, c2 in b.items():
            e = tuple(x + y for x, y in zip(e1, e2))
            out[e] = out.get(e, 0) + c1 * c2
    return {e: complex(round(c.real), round(c.imag)) for e, c in out.items() if abs(c) > 1e-12}


def A_n(j):
    return ('AV', NV.astype(complex), poly_nm(j, 0))


def A_m(j, l=0):
    return ('AV', MV.astype(complex), poly_nm(j, l))


def f_letters(j):
    """d_n^j f as [(coef, letter)]"""
    return [(1, ('AV', MV.astype(complex), poly_nm(j + 1, 0))), (-1, ('AV', NV.astype(complex), poly_nm(j, 1)))]


def DkF(k):
    """D_n^k F to two fields: list of (coef, [letters])"""
    out = [(c, [L]) for c, L in f_letters(k)]
    for j in range(k + 1):                                   # d^k [A_n, A_m]
        c = math.comb(k, j)
        out.append((c, [A_n(j), A_m(k - j)])); out.append((-c, [A_m(k - j), A_n(j)]))
    for i in range(k):                                       # [d^i A_n, d^{k-1-i} f]
        c = math.comb(k, i + 1)
        for cf, L in f_letters(k - 1 - i):
            out.append((c * cf, [A_n(i), L])); out.append((-c * cf, [L, A_n(i)]))
    return out


def op(k, N):
    """tr(D^k F D^{N-k} F) to three fields"""
    a, b = DkF(k), DkF(N - k)
    out = []
    for (c1, w1), (c2, w2) in itertools.product(a, b):
        if len(w1) + len(w2) <= 3:
            out.append((c1 * c2, w1 + w2))
    return out


def two_field(o):
    return [(c, w) for c, w in o if len(w) == 2]


def onshell_leg(q, seed=0):
    """polarisation e with e.q = 0 and f(e, q) != 0 (Gaussian integers)"""
    rng = np.random.default_rng(seed)
    r = np.array([1, 2, -1, 1], dtype=complex)
    for _ in range(50):
        w = rng.integers(-2, 3, 4) + 1j * rng.integers(-2, 3, 4)
        a = complex(np.dot(w, q)); b = complex(np.dot(r, q))
        if abs(b) < 1e-9:
            r = rng.integers(-2, 3, 4).astype(complex); continue
        e = b * w - a * r                                    # e.q = b (w.q) - a (r.q) = 0
        f = 1j * np.dot(NV, q) * np.dot(MV, e) - 1j * np.dot(MV, q) * np.dot(NV, e)
        if abs(f) > 1e-9:
            return np.array([complex(round(z.real), round(z.imag)) for z in e])
    raise RuntimeError('no polarisation')


def kin(npts, seed=2):
    qs = [q for q in S11.NULLS if abs(complex(np.dot(NV, q))) + abs(complex(np.dot(MV, q))) > 1e-9]
    pairs = list(itertools.combinations(range(len(qs)), 2))
    rng = np.random.default_rng(seed); rng.shuffle(pairs)
    out = []
    for i, j in pairs:
        q1, q2 = qs[i], qs[j]
        if abs(complex(np.dot(q1, q2))) < 1e-9:
            continue
        out.append((q1, q2))
        if len(out) == npts:
            break
    return out


def run(N, kinds=('L', 'M0', 'Pn', 'ward'), log=print, npts_extra=2):
    ks = list(range(N // 2 + 1))
    ops = [op(k, N) for k in ks]
    nb = len(ks)
    pts = kin(nb + npts_extra + 2)
    legsets = [[('AV', onshell_leg(q1, 1), q1), ('AV', onshell_leg(q2, 2), q2)] for q1, q2 in pts]
    A = sp.Matrix([[F2.symc(GO.op_graded(two_field(op(j, N)), lg)) for j in ks] for lg in legsets])
    out = {'N': N}
    mats = {}
    caches = {}
    for kind in kinds:
        t0 = time.time()
        if kind == 'ward':
            res = {}
            for kd in ('D0', ('Pv', tuple(NV)), ('Pv', (1, 2, -1, 3))):
                vals = []
                for (q1, q2) in pts[:3]:
                    lg = [('AV', q1.astype(complex), q1), ('AV', onshell_leg(q2, 2), q2)]   # leg 1 longitudinal
                    for o in ops:
                        a, b = GO.planar_graded(o, lg, lambda x, y, w, kd=kd: GV.kernel(x, y, w, kd), cache=caches.setdefault(str(kd), {}))
                        vals.append((F2.symv(a), F2.symv(b)))
                res[str(kd)] = {'letter_part_zero': all(sp.simplify(a) == 0 for a, b in vals),
                                'vv_part_zero': all(sp.simplify(b) == 0 for a, b in vals),
                                'letter': [str(a) for a, b in vals], 'vv': [str(b) for a, b in vals]}
                log('  gluon N=%d ward %s: letter zero %s, vv zero %s' % (N, kd, res[str(kd)]['letter_part_zero'], res[str(kd)]['vv_part_zero']))
            out['ward'] = res
            continue
        Yl, Yv = [], []
        for lg, (q1, q2) in zip(legsets, pts):
            rl, rv = [], []
            for o in ops:
                if kind == 'L':
                    val = GO.planar_graded(o, lg, lambda a, b, w: GV.logkernel(a, b, w), cache=caches.setdefault('L', {}))
                    rl.append(F2.symv(val) - 4 * F2.symc(GO.op_graded(two_field(o), lg))); rv.append(0)
                else:
                    kd = 'D0' if kind == 'M0' else ('Pv', tuple(NV))
                    a, b = GO.planar_graded(o, lg, lambda x, y, w, kd=kd: GV.kernel(x, y, w, kd), cache=caches.setdefault(str(kd), {}))
                    ya, yb = F2.symv(a), F2.symv(b)
                    if kind == 'Pn':
                        tgt = sp.Rational(-1, 3) * F2.symc(1j * np.dot(NV, q1 + q2)) * F2.symc(GO.op_graded(two_field(o), lg))
                        ya = ya - tgt
                    rl.append(sp.expand(ya)); rv.append(sp.expand(yb))
            Yl.append(rl); Yv.append(rv)
        Yl = sp.Matrix(Yl); Yv = sp.Matrix(Yv)
        if kind == 'Pn':
            out[kind] = {'letter_residual_vs_-1/3_i_n.dq': [str(x) for x in Yl[:3, :]], 'vv_part': [str(x) for x in Yv[:3, :]],
                         'letter_zero': all(sp.simplify(x) == 0 for x in Yl), 'vv_zero': all(sp.simplify(x) == 0 for x in Yv),
                         'seconds': round(time.time() - t0, 1)}
            log('  gluon N=%d Pn: letter residual zero %s, vv zero %s (%.1fs)' % (N, out[kind]['letter_zero'], out[kind]['vv_zero'], time.time() - t0))
            continue
        inv = (A.H * A).inv()
        Ml = sp.simplify(inv * (A.H * Yl)); Mv = sp.simplify(inv * (A.H * Yv))
        rl = sp.simplify(A * Ml - Yl); rv = sp.simplify(A * Mv - Yv)
        ok = all(x == 0 for x in rl) and all(x == 0 for x in rv)
        mats[kind] = (Ml, Mv)
        out[kind] = {'letter_matrix': [[str(Ml[i, j]) for j in range(nb)] for i in range(nb)],
                     'vv_matrix': [[str(Mv[i, j]) for j in range(nb)] for i in range(nb)], 'fit_exact': ok,
                     'seconds': round(time.time() - t0, 1)}
        log('  gluon N=%d %s: letter %s vv %s exact %s (%.1fs)' % (N, kind, Ml.tolist(), Mv.tolist(), ok, time.time() - t0))
    if 'L' in mats and 'M0' in mats:
        kap = sp.Symbol('kappa')
        D1 = 4 * mats['L'][0] + mats['M0'][0] + kap * mats['M0'][1] + sp.Rational(2, 3) * (6 + N) * sp.eye(nb)
        target = sorted([-32 * sum(sp.Rational(1, j) for j in range(1, s + 3)) for s in range(0, N + 1, 2)])
        out['D1'] = {'matrix(kappa)': [[str(sp.simplify(D1[i, j])) for j in range(nb)] for i in range(nb)],
                     'target_-32h(n+2)_n_even<=N': [str(x) for x in target]}
        for kv in (1, -1, 0):
            ev = D1.subs(kap, kv).eigenvals()
            got = sorted([sp.nsimplify(e) for e, m in ev.items() for _ in range(m)], key=lambda z: float(sp.re(z)))
            out['D1']['kappa=%d' % kv] = {'eigenvalues': [str(x) for x in got],
                                         'pass': len(got) == len(target) and all(sp.simplify(a - b) == 0 for a, b in zip(got, target))}
        log('  gluon N=%d D1: %s' % (N, {k: v for k, v in out['D1'].items() if k.startswith('kappa=')}))
    return out
