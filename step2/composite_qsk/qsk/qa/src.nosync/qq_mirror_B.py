"""QA: symmetry tests of engine B's code (read-only import) on the beta-beta defect cycles.
g = rho(e_0) rho(e_4) (rotation by pi in the (e_0, e_4) plane, an element of Spin(4+6) preserving V + W_perp) maps
S_+ (x) W to S_- (x) W^*.  Data transform: xi -> g xi, spinor outputs -> g pol, spinor letters -> g^{-T} zeta, vector
components 0 (V) and 4 (W_perp) change sign (outputs and momenta).  The kernel must be invariant.
Second test: the permutation symmetry of W_perp exchanging the planes (4,7) and (5,8) (with its spinor lift)."""
import os, sys, json, random
import sympy as sp

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_b/src.nosync'))
import qb_job as QJ            # noqa: E402
import realization as RZ       # noqa: E402

LOG = {'1': 1, 'log2': sp.log(2), 'log3': sp.log(3), 'log5': sp.log(5)}
g = sp.simplify(RZ.Rb[0] * RZ.R[4])
assert all(g[i, j] == 0 for i in range(8) for j in range(8)) and all(g[i, j] == 0 for i in range(8, 16) for j in range(8, 16))
ginvT = g.inv().T


def value(ta, tb, za, zb, outs, pols, moms, xi):
    job = QJ.DJob(ta, tb, outs, 'xi', xi, pol=pols, za=za, zb=zb, tag='mir%d' % random.randint(0, 10 ** 8))
    res = job.kernel(mlist=(0,))
    sub = {}
    for j, m in enumerate(moms):
        for mu in range(4):
            sub[sp.Symbol('r%d_%d' % (j, mu))] = m[mu]
    out = {}
    for wd, dm in res['total'].items():
        tot = 0
        for bb, P in dm.get(0, {}).items():
            tot += sp.expand(P.as_expr(*[sp.Symbol(n) for n in job.names]).subs(sub)) * LOG[bb]
        out[str(list(wd))] = sp.nsimplify(sp.expand(tot))
    return out


def vecd(v, rng_):
    return {i: v[i] for i in rng_ if v[i] != 0}


random.seed(11)
xi = sp.zeros(16, 1); xi[3] = 1
# B's own twisted data: beta_1 beta_2 -> gamma^3 (dual vector), dc (comp 11)
zaB = sp.zeros(16, 1); zaB[0] = 1
zbB = sp.zeros(16, 1); zbB[1] = 1
gam = [[RZ.Sig[K][J - 1, 3] for K in range(6)] for J in (1, 2, 3)]
comp = [[RZ.Sig[K][I - 1, J - 1] for K in range(6)] for (I, J) in ((1, 2), (1, 3), (2, 3))]
Linv = sp.Matrix(gam + comp).inv()
polA = [0] * 10
for c in range(6):
    polA[4 + c] = Linv[c, 2]
poll = sp.zeros(16, 1); poll[11] = 1; poll[15] = 2
moms = [[random.randint(-3, 3) for _ in range(4)] for _ in range(2)]
v0 = value('lamP', 'lamP', vecd(zaB, range(8)), vecd(zbB, range(8)), ['AW', 'lamM'],
           [{c: polA[c] for c in range(4, 10)}, vecd(poll, range(8, 16))], moms, [xi[i] for i in range(16)])
print('B original (beta1 beta2 -> gamma3 dc):', v0, flush=True)
# mirrored
xi2 = g * xi
za2 = ginvT * zaB; zb2 = ginvT * zbB
pol2 = list(polA); pol2[4] = -pol2[4]
poll2 = g * poll
moms2 = [[-m[0]] + m[1:] for m in moms]
v1 = value('lamM', 'lamM', vecd(za2, range(8, 16)), vecd(zb2, range(8, 16)), ['AW', 'lamP'],
           [{c: pol2[c] for c in range(4, 10)}, vecd(poll2, range(8))], moms2, [xi2[i] for i in range(16)])
print('B mirrored:', v1, flush=True)
json.dump({'original': {k: str(v) for k, v in v0.items()}, 'mirrored': {k: str(v) for k, v in v1.items()}},
          open(os.path.join(QSK, 'qa/out.nosync/mirror_B.json'), 'w'), indent=1)
