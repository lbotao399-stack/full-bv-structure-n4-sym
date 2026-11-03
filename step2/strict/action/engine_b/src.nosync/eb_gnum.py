"""Gates G-cyc (cyclic reduction of the word sum, arity one and two) and G-masters (exact Feynman-parameter integrals
ce_int.J against direct mpmath quadrature)."""
import json, os, random
from fractions import Fraction as Fr
import numpy as np
import eb_core as E, eb_cycle as Y, eb_gens as G, eb_t1
import ce_int as CI

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')


def same(a, b):
    for k in set(a) | set(b):
        for m in set(a.get(k, {})) | set(b.get(k, {})):
            x = a.get(k, {}).get(m, E.CLog()); y = b.get(k, {}).get(m, E.CLog())
            d = CI.LogNum(0)
            if not (x.re - y.re).is_zero() or not (x.im - y.im).is_zero():
                return False
    return True


def run(log=print):
    random.seed(31)
    out = {}
    tdir = G.tmpl_dir(); tw = G.tmpl_weights(2, 3, 0, 'D')
    legs = eb_t1.make_legs('AVAV')([Fr(0)] * 4)
    e = np.zeros(4, complex); e[1] = 1
    gens1 = [Y.Gen('D', 1, tdir, gvec=e, inj=[Fr(0)] * 4)]
    a = Y.evaluate(legs, gens1, cyclic=True); b = Y.evaluate(legs, gens1, cyclic=False)
    out['G-cyc arity one (A_V,A_V), Vdir'] = same(a, b)
    gens2 = [Y.Gen('D', 1, tw), Y.Gen('P', 1, tdir, gvec=np.array([1, -1, 2, 0], complex) / 2)]
    legs2 = eb_t1.make_legs('cAsV')([Fr(0)] * 4)
    a = Y.evaluate(legs2, gens2, cyclic=True); b = Y.evaluate(legs2, gens2, cyclic=False)
    out['G-cyc arity two (c,A*_V), (V0,P_b)'] = same(a, b)
    worst = 0.0
    for (m, aa, nvec) in [(0, 1, (2,)), (1, 0, (3,)), (0, 0, (0, 0)), (0, 1, (2, 0)), (1, 0, (1, 1)), 
                          (1, 1, (0, 2))]:
        ex = float(CI.J(m, aa, nvec, (1,) * len(nvec)))
        nm = float(CI.J_numeric(m, aa, nvec, (1,) * len(nvec), dps=15))
        worst = max(worst, abs(ex - nm) / max(1e-30, abs(ex)))
    out['G-masters max relative difference'] = worst
    log(f'G-num: {out}')
    return out


if __name__ == '__main__':
    with open(os.path.join(OUT, 'G_num.json'), 'w') as f:
        json.dump(run(), f, indent=1)
