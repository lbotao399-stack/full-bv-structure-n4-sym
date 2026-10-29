"""I-a classes of the cubic obstruction, own engine (rb_loop), all four sectors.
Triangles: weight 1; quartic-cubic bubbles: weight 1/2 summed over both X positions (|Aut| = 2).
For each class: eps^-1 (omega = 2 only) and eps^0 coefficients, exact in Q + Q log(3/2), fitted to the G basis
on NCONF random integer configurations (residual must vanish), in units hbar tr(ad ad)/(16 pi^2 (X,Y)).
usage: python3 rb_run_cubic.py SECTOR [NCONF]
"""
import sys, os, json, time, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from fractions import Fraction as Fr
from rb_alg import KILLING_RATIO
import rb_func as RF
import rb_loop as RL

here = os.path.dirname(os.path.abspath(__file__))
sector = sys.argv[1]
NCONF = int(sys.argv[2]) if len(sys.argv) > 2 else 4
rng = np.random.default_rng({'cAVAV': 11, 'cll': 22, 'cAWAW': 33, 'ccAs': 44}[sector])


def rint(shape, lo=-3, hi=4):
    return rng.integers(lo, hi, size=shape).astype(float)


def rand_mom():
    while True:
        p1, p2 = rint(4), rint(4)
        p0 = -p1 - p2
        if p0.any() and p1.any() and p2.any():
            return p0, p1, p2


SPEC = {'cAVAV': (('c', 'A', 'A'), 'V', ['G1', 'G2'], 2),
        'cAWAW': (('c', 'A', 'A'), 'W', ['G3'], 2),
        'cll': (('c', 'lam', 'lam'), None, ['G5'], 1),
        'ccAs': (('c', 'c', 'Ast'), 'V', ['G4'], 1)}
slots, part, Gnames, omega = SPEC[sector]
GF = {'G1': RF.F_G1, 'G2': RF.F_G2, 'G3': RF.F_G3, 'G4': RF.F_G4, 'G5': RF.F_G5}


def make_conf():
    ps = rand_mom()
    legs = []
    for s, p in zip(slots, ps):
        v = rint(RF.SHAPE[s])
        if s in ('A', 'Ast'):
            if part == 'V':
                v[4:] = 0
            else:
                v[:4] = 0
        legs.append((s, v, p))
    return legs


EDGE_LAB = {('A', 'A'): 'A', ('c', 'Ast'): 'gh', ('Ast', 'c'): 'gh', ('lam', 'lam'): 'lam',
            ('A', 'Ast'): 'A', ('Ast', 'A'): 'A', ('lam', 'lamst'): 'lam', ('lamst', 'lam'): 'lam',
            ('c', 'cst'): 'c', ('cst', 'c'): 'c'}


def classkey(graph, parts, ext_slots):
    verts, kinds = graph['verts'], graph['kinds']
    n = len(verts)
    vs = []
    for j in range(n):
        vs.append(verts[j][0] + '[' + ','.join(ext_slots[j]) + ']')
    es = []
    for j in range(n):
        t1 = verts[j][4]; t2 = verts[(j + 1) % n][5]
        lab = EDGE_LAB[(t1, t2)]
        if lab == 'A':
            lab += '_' + parts[j]
        es.append(kinds[j] + ':' + lab)
    return ' | '.join(sorted(vs)) + ' || ' + ' '.join(sorted(es))


topologies = {'triangle': [(slots[0],), (slots[1],), (slots[2],)],
              'bubble': [(slots[0],), (slots[1], slots[2])]}
confs = [make_conf() for _ in range(NCONF)]
Gv = [[RF.value(GF[g], L) for g in Gnames] for L in confs]
t0 = time.time()
classes = {}       # key -> per conf: [cm1_rat, cm1_log, c0_rat, c0_log]
exact = {}         # key -> per conf: 4 x (Fr re, Fr im)
npat = {}
for topo, vext in topologies.items():
    graphs = list(RL.enumerate_cycles(vext))
    for ci, L in enumerate(confs):
        if topo == 'triangle':
            ext = [[(L[0][1].reshape(-1), L[0][2])], [(L[1][1].reshape(-1), L[1][2])], [(L[2][1].reshape(-1), L[2][2])]]
        else:
            ext = [[(L[0][1].reshape(-1), L[0][2])], [(L[1][1].reshape(-1), L[1][2]), (L[2][1].reshape(-1), L[2][2])]]
        Ptot = [sum((e[1] for e in ext[j]), np.zeros(4)) for j in range(len(ext))]
        cvec = [np.zeros(4)]
        for j in range(1, len(ext)):
            cvec.append(cvec[-1] + Ptot[j])
        w = Fr(1) if topo == 'triangle' else Fr(1, 2)
        kind = 'tri' if topo == 'triangle' else 'bub'
        for g in graphs:
            for parts in RL.a_parts(g):
                fn = lambda l, g=g, parts=parts: RL.graph_N(g, parts, ext, l)
                poly = RL.interpolate(fn)
                key = topo + ': ' + classkey(g, parts, vext)
                if ci == 0:
                    npat[key] = npat.get(key, 0) + 1
                acc = classes.setdefault(key, [[0j, 0j, 0j, 0j] for _ in confs])
                ex = exact.setdefault(key, [[[Fr(0), Fr(0)] for _ in range(4)] for _ in confs])
                if not poly:
                    continue
                assert RL.p_eval_check(poly, fn) < 1e-6
                jx = g['kinds'].index('X')
                if omega == 2:
                    (r, lg), (ri, li) = RL.coefficient(poly, cvec, 0, jx, kind)
                    acc[ci][0] += complex(float(w * r), float(w * ri)) / KILLING_RATIO
                    acc[ci][1] += complex(float(w * lg), float(w * li)) / KILLING_RATIO
                    for idx, (a_, b_) in [(0, (r, ri)), (1, (lg, li))]:
                        ex[ci][idx][0] += w * a_ / Fr(-2); ex[ci][idx][1] += w * b_ / Fr(-2)
                (r, lg), (ri, li) = RL.coefficient(poly, cvec, omega, jx, kind)
                acc[ci][2] += complex(float(w * r), float(w * ri)) / KILLING_RATIO
                acc[ci][3] += complex(float(w * lg), float(w * li)) / KILLING_RATIO
                for idx, (a_, b_) in [(2, (r, ri)), (3, (lg, li))]:
                    ex[ci][idx][0] += w * a_ / Fr(-2); ex[ci][idx][1] += w * b_ / Fr(-2)
    print(topo, 'done', round(time.time() - t0, 1), 's', flush=True)


def fit(vals):
    A = np.array(Gv, dtype=complex); b = np.array(vals, dtype=complex)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    res = float(np.max(np.abs(A @ sol - b))) if len(b) else 0.0
    out = []
    for s in sol:
        out.append((str(Fr(s.real).limit_denominator(100000)), float(abs(s.imag)),
                    float(abs(s.real - float(Fr(s.real).limit_denominator(100000))))))
    return out, res


GvF = [[(Fr(complex(x).real), Fr(complex(x).imag)) for x in row] for row in Gv]


def exact_check(key_rows, idx, coefs):
    # verify  value = sum_i coef_i G_i  exactly for every configuration (Fractions)
    worst = Fr(0)
    for ci in range(len(confs)):
        re, im = key_rows[ci][idx]
        sr = sum((Fr(coefs[i]) * GvF[ci][i][0] for i in range(len(coefs))), Fr(0))
        si = sum((Fr(coefs[i]) * GvF[ci][i][1] for i in range(len(coefs))), Fr(0))
        worst = max(worst, abs(re - sr), abs(im - si))
    return str(worst)


report = {'sector': sector, 'G': Gnames, 'omega': omega, 'NCONF': NCONF, 'classes': {}}
total = [[0j] * 4 for _ in confs]
for key, rows in sorted(classes.items()):
    allzero = max(abs(x) for r in rows for x in r) < 1e-12
    entry = {'patterns': npat.get(key, 0)}
    if allzero:
        entry['value'] = 'identically 0'
    else:
        for idx, nm in [(0, 'eps^-1 rational'), (1, 'eps^-1 log(3/2)'), (2, 'eps^0 rational'), (3, 'eps^0 log(3/2)')]:
            if omega == 1 and idx < 2:
                continue
            coef, res = fit([r[idx] for r in rows])
            entry[nm] = {'coef_on_' + ','.join(Gnames): [c[0] for c in coef], 'resid': res,
                         'max_imag': max(c[1] for c in coef), 'max_rationalisation_error': max(c[2] for c in coef),
                         'EXACT_max_deviation': exact_check(exact[key], idx, [c[0] for c in coef])}
    report['classes'][key] = entry
    for ci in range(len(confs)):
        for idx in range(4):
            total[ci][idx] += rows[ci][idx]
    print(key, json.dumps(entry), flush=True)
tot = {}
for idx, nm in [(0, 'eps^-1 rational'), (1, 'eps^-1 log(3/2)'), (2, 'eps^0 rational'), (3, 'eps^0 log(3/2)')]:
    if omega == 1 and idx < 2:
        continue
    coef, res = fit([t[idx] for t in total])
    totex = [[[sum((exact[k][ci][idx][0] for k in exact), Fr(0)), sum((exact[k][ci][idx][1] for k in exact), Fr(0))] for idx in range(4)] for ci in range(len(confs))]
    tot[nm] = {'coef_on_' + ','.join(Gnames): [c[0] for c in coef], 'resid': res,
               'EXACT_max_deviation': exact_check(totex, idx, [c[0] for c in coef])}
report['TOTAL_I-a'] = tot
print('TOTAL', json.dumps(tot), flush=True)
report['seconds'] = round(time.time() - t0, 1)
json.dump(report, open(os.path.join(here, 'rb_cubic_%s.json' % sector), 'w'), indent=1)
