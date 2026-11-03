"""QA: canonical twisted letters built directly (basis-free recipe) for engine A's xi (mapped to B's basis) and for
B's xi; beta-beta evaluated with engine B's code (read-only).  Also: A's letters versus the canonical ones."""
import os, sys, json, random
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'qa/src.nosync'))
import qq_ab_compare as AB
import qq_mirror_B as MB
import realization as RZ
R = RZ.R
random.seed(21)
def canon(xi):
    off = 0 if any(xi[i] != 0 for i in range(8)) else 8
    X = sp.Matrix(2, 4, lambda a, A: xi[off + 4 * a + A])
    j = [k for k in range(4) if any(X[a, k] != 0 for a in range(2))][0]
    u = X[:, j]; i = [a for a in range(2) if X[a, j] != 0][0]; z = X[i, :].T / u[i]
    ub = sp.Matrix([sp.conjugate(x) for x in u]); zb = sp.Matrix([sp.conjugate(x) for x in z])
    T = sp.Matrix([list(z)]).nullspace()        # t with t.z = 0
    betas = []
    for t in T:
        v = [0] * 16
        for a in range(2):
            for A in range(4):
                v[off + 4 * a + A] = sp.expand(ub[a] * t[A])
        betas.append(v)
    off2 = 8 - off
    dcs = []
    for s in ((1, 0), (0, 1)):
        v = [0] * 16
        for a in range(2):
            for A in range(4):
                v[off2 + 4 * a + A] = sp.expand(s[a] * zb[A])
        dcs.append(v)
    cols = sp.Matrix.hstack(*[R[4 + K] * sp.Matrix(xi) for K in range(6)])
    N = cols.nullspace()
    gams = [[sp.conjugate(x) for x in n] for n in N]
    return off, betas, dcs, gams
def run(xi, label):
    off, betas, dcs, gams = canon(xi)
    tl = 'lamP' if off == 0 else 'lamM'; to = 'lamM' if off == 0 else 'lamP'
    res = []
    for (i, j) in ((0, 1), (0, 2), (1, 2)):
        Y0 = [sum(random.randint(1, 3) * g[c] for g in gams) for c in range(6)]
        Y1 = [random.randint(1, 3) * dcs[0][c] + random.randint(1, 3) * dcs[1][c] for c in range(16)]
        moms = [[random.randint(-3, 3) for _ in range(4)] for _ in range(2)]
        v = MB.value(tl, tl, MB.vecd(betas[i], range(off, off + 8)), MB.vecd(betas[j], range(off, off + 8)), ['AW', to],
                     [{4 + c: Y0[c] for c in range(6)}, MB.vecd(Y1, range(8 - off, 16 - off))], moms, list(xi))
        print(label, (i, j), v, flush=True)
        res.append({k: str(x) for k, x in v.items()})
    return res
out = {}
xiB = [0] * 16; xiB[3] = 1
out['B_xi canonical'] = run(xiB, 'B xi')
d = json.load(open(os.path.join(QSK, 'qa/out.nosync/twistA_data.json')))
xiA = [sp.nsimplify(x) for x in sp.expand(AB.S * sp.Matrix([AB.cplx(x) for x in d['xi']]))]
out['A_xi canonical'] = run(xiA, 'A xi')
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/canon_B.json'), 'w'), indent=1)
