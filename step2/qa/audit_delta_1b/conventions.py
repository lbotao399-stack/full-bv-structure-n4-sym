"""Which conventions for (second-order operator of a kernel, induced derivation of an odd operator)
are compatible with the manuscript's kernel-level data (P = int (QGF x 1) K_t dt untwisted,
(Q x 1 + 1 x Q) P = K_eps - K_T, Delta_t = -d_{K_t}, bracket of def:c3-delta)?

Requirements tested on functionals (C07 6.4 and the BV formalism):
 R1 [d_P, Qf] = d_{K_eps} - d_{K_T}           (C07 Lemma 10.2.1 in operator form; lem:c8-transfer proof)
 R2 Qf^2 = 0,  R3 Qf Delta + Delta Qf = 0,  R4 Qf = {S_free, -} for a quadratic S_free,
 R5 [D_X, d_phi] = -d_{(X x 1 + 1 x X) phi} for even and odd X and kernels phi of both parities
     (Lemma d-calculus (i) of delta_theory.tex, used there for odd X = Q and X = P*M).
Here Qf = sigma * D_Q with sigma = +-1.
"""
import json, os, random, sys, itertools
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp import Alg, mat_mul, mat_add, mat_scale, tensor_act
from model import build

HERE = os.path.dirname(os.path.abspath(__file__))

def test_functionals(alg, rng, n_cubic=600, n_quartic=3, terms=40):
    n = alg.n
    fs = []
    for i in range(n):
        fs.append(alg.x(i))
    for i in range(n):
        for j in range(i, n):
            m, s = alg.sort((i, j))
            if m is not None:
                fs.append({m: Fr(1)})
    for _ in range(n_cubic):
        m, s = alg.sort(tuple(rng.randrange(n) for _ in range(3)))
        if m is not None:
            fs.append({m: Fr(1)})
    for _ in range(n_quartic):
        f = {}
        while len(f) < terms:
            m, s = alg.sort(tuple(rng.randrange(n) for _ in range(rng.choice([3, 4]))))
            if m is None:
                continue
            f[m] = Fr(rng.randint(-5, 5) or 1, rng.randint(1, 4))
        # split into parity-homogeneous parts
        for p in (0, 1):
            g = {mm: c for mm, c in f.items() if alg.parity_of_mono(mm) == p}
            if g:
                fs.append(g)
    return fs

def solve_exact(rows, rhs):
    """solve rows * u = rhs exactly (consistent, possibly overdetermined); None if inconsistent.
    rows: list of dense lists or of dicts {col: coef}; sparse Gauss-Jordan elimination."""
    if not rows:
        return []
    if isinstance(rows[0], dict):
        ncol = 1 + max([max(r) for r in rows if r] + [-1])
        R = [(dict(r), b) for r, b in zip(rows, rhs)]
    else:
        ncol = len(rows[0])
        R = [({c: v for c, v in enumerate(r) if v != 0}, b) for r, b in zip(rows, rhs)]
    pivots = {}   # col -> (row dict normalised, rhs)
    for r, b in R:
        r = dict(r)
        # reduce by existing pivots
        changed = True
        while True:
            hit = [c for c in r if c in pivots]
            if not hit:
                break
            for c in hit:
                if c not in r:
                    continue
                f = r[c]
                pr, pb = pivots[c]
                for cc, vv in pr.items():
                    r[cc] = r.get(cc, 0) - f * vv
                    if r[cc] == 0:
                        del r[cc]
                b = b - f * pb
        if not r:
            if b != 0:
                return None
            continue
        c0 = min(r)
        pv = r[c0]
        r = {cc: vv / pv for cc, vv in r.items()}
        b = b / pv
        # eliminate c0 from existing pivots
        for c, (pr, pb) in list(pivots.items()):
            if c0 in pr:
                f = pr[c0]
                for cc, vv in r.items():
                    pr[cc] = pr.get(cc, 0) - f * vv
                    if pr[cc] == 0:
                        del pr[cc]
                pivots[c] = (pr, pb - f * b)
        pivots[c0] = (r, b)
    u = [Fr(0)] * ncol
    for c, (pr, pb) in pivots.items():
        # free columns set to 0: pr = {c:1, free...}
        u[c] = pb
    return u

def bracket(alg, K, f, g, conv):
    """{f,g}_K = Delta(fg) - Delta(f) g - (-1)^{|f|} f Delta(g), Delta = -d_K (K odd)"""
    Lap = lambda h: alg.scale(-1, alg.dphi(K, h, conv))
    pf = alg.parity(f)
    return alg.add(Lap(alg.mul(f, g)), (-1, alg.mul(Lap(f), g)), (-(-1) ** pf, alg.mul(f, Lap(g))))

def generator_of(alg, K, op_on_linear, target_degree, conv):
    """quadratic S of functional degree target_degree with {S, x^k}_K = op_on_linear(x^k) for all k"""
    n = alg.n
    monos = []
    for i in range(n):
        for j in range(i, n):
            if alg.deg[i] + alg.deg[j] != target_degree:
                continue
            m, s = alg.sort((i, j))
            if m is not None:
                monos.append(m)
    cols = []
    for mono in monos:
        cols.append([bracket(alg, K, {mono: Fr(1)}, alg.x(k), conv) for k in range(n)])
    rows, rhs = [], []
    for k in range(n):
        tgt = op_on_linear(alg.x(k))
        keys = set(tgt)
        for c in cols:
            keys |= set(c[k])
        for key in keys:
            rows.append({ci: c[k][key] for ci, c in enumerate(cols) if key in c[k]}); rhs.append(tgt.get(key, Fr(0)))
    u = solve_exact(rows, rhs)
    if u is None:
        return None
    return alg.clean({mo: v for mo, v in zip(monos, u)})

def run(N=4, seed=3):
    rng = random.Random(seed)
    res = {}
    m = build(N=N)
    alg, par = m['alg'], m['par']
    e, T = m['eps'], m['T']
    P = m['P'](e, T); Ke = m['Kt'](e); KT = m['Kt'](T); K0 = m['K0']
    Q = m['Q']; M = m['M']
    fs = test_functionals(alg, rng)
    res['number of test functionals'] = len(fs)
    for dconv in ('A', 'B'):
        for qconv in ('naive', 'canon'):
            for sig in (1, -1):
                key = 'd-conv %s, D-conv %s, sigma %+d' % (dconv, qconv, sig)
                Qf = lambda f: alg.scale(sig, alg.D(Q, 1, f, qconv))
                r = {}
                # R1
                err = Fr(0)
                for f in fs:
                    lhs = alg.add(alg.dphi(P, Qf(f), dconv), (-1, Qf(alg.dphi(P, f, dconv))))
                    rhs = alg.add(alg.dphi(Ke, f, dconv), (-1, alg.dphi(KT, f, dconv)))
                    err = max(err, alg.dist(lhs, rhs))
                r['R1 [d_P,Qf] = d_Keps - d_KT'] = str(err)
                # R2
                r['R2 Qf^2 = 0'] = str(max(alg.maxabs(Qf(Qf(f))) for f in fs[:400]))
                # R3
                Lap = lambda h: alg.scale(-1, alg.dphi(K0, h, dconv))
                r['R3 Qf Delta + Delta Qf = 0'] = str(max(alg.maxabs(alg.add(Qf(Lap(f)), Lap(Qf(f)))) for f in fs))
                # R4
                Sfree = generator_of(alg, K0, Qf, 0, dconv)
                r['R4 S_free exists'] = Sfree is not None
                res[key] = r
    # R5 Lemma (i), for each convention pair, X in {M (even), Q (odd), P*M (odd)}, phi in {P (even), K_eps (odd)}
    PM = mat_mul(m['Pstar'](e, T), M)
    for dconv in ('A', 'B'):
        for qconv in ('naive', 'canon'):
            key = 'R5 d-conv %s, D-conv %s' % (dconv, qconv)
            r = {}
            for (xn, X, pX) in (('M', M, 0), ('Q', Q, 1), ('P*M', PM, 1)):
                for (pn, phi, pphi) in (('P', P, 0), ('K_eps', Ke, 1)):
                    Xphi = tensor_act(X, pX, phi, par)
                    s = -1 if (pX and pphi) else 1
                    err = Fr(0)
                    for f in fs[:1800]:
                        lhs = alg.add(alg.D(X, pX, alg.dphi(phi, f, dconv), qconv),
                                      (-s, alg.dphi(phi, alg.D(X, pX, f, qconv), dconv)))
                        rhs = alg.scale(-1, alg.dphi(Xphi, f, dconv))
                        err = max(err, alg.dist(lhs, rhs))
                    r['[D_%s, d_%s] + d_{(X x1+1x X)%s}' % (xn, pn, pn)] = str(err)
            res[key] = r
    return res

if __name__ == '__main__':
    out = run()
    with open(os.path.join(HERE, 'conventions_results.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
