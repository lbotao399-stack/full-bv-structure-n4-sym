"""Step 2, Part D4, engine A: one-loop 1PI two-point functions W_1^(2) (eps-expansion incl. finite parts) in the sectors
A|_Wperp, A|_V, lambda, (c, A^*|_V); plane-wave values (coefficient of the odd parameters in leg order), legs
(field 1 at momentum p, field 2 at -p), fitted exactly on the invariant structures.  Units hbar tr(ad ad)/(16 pi^2 (X,Y))."""
import os, sys, json, time
import numpy as np
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import s2_core as C
from r2_core import Q2
import d4_two_point as TP
import d4_pi_sectors as PS
from cliff import s as SIG, sb as SIGB
SB = np.array(SIGB); SS = np.array(SIG)
with np.errstate(all='ignore'):
    CHIS = np.rint((SS[0] @ SB[1] @ SS[2] @ SB[3]).real)   # 4d chirality on Sigma
GR = C.G.real


def structs(sec, a, b, p):
    p = np.asarray(p, float); p2 = float(p @ p)
    cd = PS.cdot
    if sec == 'AP':
        s_ = sum(cd(a[M], b[M]) for M in range(4, 10))
        return [s_, s_ * p2], ['(a,b)', 'p^2 (a,b)']
    if sec == 'AV':
        ab = sum(cd(a[m], b[m]) for m in range(4))
        pa = sum(p[m] * a[m] for m in range(4)); pb = sum(p[m] * b[m] for m in range(4))
        return [ab, ab * p2, cd(pa, pb)], ['(a.b)', 'p^2 (a.b)', '(p.a)(p.b)']
    if sec == 'l':
        rp = np.einsum('m,mAB->AB', p, SB[:4])
        with np.errstate(all='ignore'):
            v1 = np.einsum('Aa,AB,Bb,ab->', a, rp, b, GR); v2 = np.einsum('Aa,AB,Bb,ab->', a, rp @ CHIS, b, GR)
        return [v1, v2], ['(u1, rho(p) u2)', '(u1, rho(p) chi u2)']
    if sec == 'gh':
        pb = sum(p[m] * b[m] for m in range(4))
        return [1j * cd(a, pb)], ['i (gamma, p.b)']
    raise ValueError(sec)


LEGS = {'AP': (('A', 'P'), ('A', 'P')), 'AV': (('A', 'V'), ('A', 'V')), 'l': (('l', None), ('l', None)),
        'gh': (('c', None), ('As', 'V'))}


def run(sec, nsamp=4, p=(1, 2, 0, 1), seed=2):
    rng = np.random.default_rng(seed)
    samples = []
    for it in range(nsamp):
        (t1, w1), (t2, w2) = LEGS[sec]
        a = C.rpol(rng, t1, w1); b = C.rpol(rng, t2, w2)
        val = TP.self_energy([(t1, a), (t2, b)], p)
        samples.append((val, structs(sec, a, b, p)[0]))
    names = structs(sec, a, b, p)[1]
    keys = sorted(set(k for v, _ in samples for k in v))
    out = {'structures': names, 'momentum': list(p), 'coeffs': {}}
    for key in keys:
        rows = [[sp.nsimplify(complex(x).real) + sp.I * sp.nsimplify(complex(x).imag) for x in st] for _, st in samples]
        rhs = [PS.tosym(v.get(key, Q2(0))) for v, _ in samples]
        sol = PS.fit(rows, rhs)
        out['coeffs'][key] = {'c': [str(x) for x in sol[0]], 'exact_fit': bool(sol[1])}
    return out


if __name__ == '__main__':
    res = {}
    for sec in sys.argv[1].split(','):
        t0 = time.time()
        res[sec] = run(sec)
        print(sec, 'time %.1f' % (time.time() - t0), json.dumps(res[sec]), flush=True)
    json.dump(res, open(os.path.join(HERE, 'logs', 'd4_se_%s.json' % C.ALG), 'w'), indent=1)
