"""(1) vertex tensors (closed forms) vs. Grassmann plane-wave evaluation of the action terms;
(2) calibration on Step-1 (Thm c4-obstruction table): special-edge bubbles with legs (c, A_V):
    eps^-1 coefficient of int(dc,A_V) and eps^0 coefficient of int(d^*dc, d^*A_V) per graph and total.
"""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
from rb_alg import KILLING_RATIO
import rb_func as RF
import rb_loop as RL

here = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(99)
out = {}


def rint(shape, lo=-3, hi=4):
    return rng.integers(lo, hi, size=shape).astype(float)


# (1) vertex validation
FUN = {'AcAst': RF.F_AcAst, 'clamlamst': RF.F_clamlamst, 'ccc': RF.F_ccc, 'lamlamA': RF.F_lamlamA,
       'AAA': RF.F_AAA, 'AAAA': RF.F_I4}
vv = {}
for vt, (slots, Tf) in RL.VERTS.items():
    worst = 0.0; mx = 0.0
    for _ in range(4):
        n = len(slots)
        ps = [rint(4) for _ in range(n - 1)]; ps.append(-sum(ps))
        vecs = [rint(RF.SHAPE[s]) for s in slots]
        legs = [(slots[i], vecs[i], ps[i]) for i in range(n)]
        g = RF.value(FUN[vt], legs)
        T = Tf(ps)
        for v in vecs:
            T = np.tensordot(v.reshape(-1), T, axes=(0, 0))
        worst = max(worst, abs(g - T)); mx = max(mx, abs(g))
    vv[vt] = {'max_dev': float(worst), 'max_val': float(mx)}
    print('vertex', vt, vv[vt], flush=True)
out['vertex_validation'] = vv


# (2) two-point special bubbles, legs (c, A_V)
def F_dcAV(f):
    return RF.int_pair(RF.dexpr(f['c']), RF.AV(f), f, RF.dot)


def classkey(graph, parts):
    verts, kinds = graph['verts'], graph['kinds']
    n = len(verts)
    vs = tuple(sorted(v[0] for v in verts))
    es = []
    for j in range(n):
        t1 = verts[j][4]; t2 = verts[(j + 1) % n][5]
        lab = {('A', 'A'): 'A', ('c', 'Ast'): 'gh', ('Ast', 'c'): 'gh', ('lam', 'lam'): 'lam',
               ('A', 'Ast'): 'A', ('Ast', 'A'): 'A', ('lam', 'lamst'): 'lam', ('lamst', 'lam'): 'lam',
               ('c', 'cst'): 'c', ('cst', 'c'): 'c'}[(t1, t2)]
        if lab == 'A':
            lab += '_' + parts[j]
        es.append(kinds[j] + ':' + lab)
    return vs + tuple(sorted(es))


graphs = list(RL.enumerate_cycles([('c',), ('A',)]))
print('number of 2-point patterns', len(graphs), flush=True)
configs = []
for _ in range(3):
    p = rint(4)
    while not p.any():
        p = rint(4)
    X0 = rint(3); X1 = rint(3); a = rint((10, 3)); a[4:] = 0
    configs.append((p, X0, X1, a))
results = {}
for ci, (p, X0, X1, a) in enumerate(configs):
    legs = [('c', X0, p), ('A', a, -p)]
    e1 = RF.value(F_dcAV, legs); e2 = RF.value(RF.F_E2, legs)
    ext = [[(X0, p)], [(a.reshape(-1), -p)]]
    for g in graphs:
        for parts in RL.a_parts(g):
            fn = lambda l, g=g, parts=parts: RL.graph_N(g, parts, ext, l)
            poly = RL.interpolate(fn)
            if not poly:
                continue
            chk = RL.p_eval_check(poly, fn)
            assert chk < 1e-6, chk
            Ptot = [p, -p]
            cvec = [np.zeros(4), Ptot[1]]
            jx = g['kinds'].index('X')
            (r1, l1), (i1, il1) = RL.coefficient(poly, cvec, 1, jx, 'bub')
            (r3, l3), (i3, il3) = RL.coefficient(poly, cvec, 3, jx, 'bub')
            key = classkey(g, parts)
            acc = results.setdefault(str(key), [[0, 0, 0, 0] for _ in configs])
            w = Fr(1, 2)
            acc[ci][0] += w * (float(r1) + 1j * float(i1)) / KILLING_RATIO   # eps^-1, rational part
            acc[ci][1] += w * (float(l1) + 1j * float(il1)) / KILLING_RATIO  # eps^-1, log2 part
            acc[ci][2] += w * (float(r3) + 1j * float(i3)) / KILLING_RATIO
            acc[ci][3] += w * (float(l3) + 1j * float(il3)) / KILLING_RATIO
    for key in results:
        results[key][ci] = [complex(x) for x in results[key][ci]] + [e1, e2]
cal = {}
tot = [0, 0]
for key, rows in results.items():
    em1 = [r[0] / r[4] for r in rows if abs(r[4]) > 0]
    e0 = [r[2] / r[5] for r in rows if abs(r[5]) > 0]
    lg = max(max(abs(r[1]) for r in rows), max(abs(r[3]) for r in rows))
    if max(abs(x) for x in em1 + e0) < 1e-12:
        continue
    cal[key] = {'eps^-1 / int(dc,A_V)': [str(Fr(x.real).limit_denominator(1000)) for x in em1],
                'eps^0 / int(d*dc,d*A_V)': [str(Fr(x.real).limit_denominator(1000)) for x in e0],
                'imag_parts_max': float(max(abs(x.imag) for x in em1 + e0)), 'log2_parts_max': float(lg)}
    tot[0] += em1[0]; tot[1] += e0[0]
    print(key, cal[key], flush=True)
cal['TOTAL'] = {'eps^-1': str(Fr(tot[0].real).limit_denominator(1000)), 'eps^0': str(Fr(tot[1].real).limit_denominator(1000))}
print('TOTAL', cal['TOTAL'], flush=True)
out['two_point_special_bubbles'] = cal
json.dump(out, open(os.path.join(here, 'rb_run_calib.json'), 'w'), indent=1)
