"""Form of [C07, 6.4] '[d_phi, Q] = d_{Q phi}' in the two consistent convention systems, for kernels of
both parities and for even and odd operators X:  [d_phi, D_X] - c * d_{(X x 1 + 1 x X)phi} for c = +1, -1
(graded commutators).  Kernels: P (even), K_eps (odd), the composite kernel Psi of the dilatation (odd).
"""
import json, os, random, sys
from fractions import Fraction as Fr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sp import mat_mul, mat_add, tensor_act
from model import build
from conventions import test_functionals
from d1_checks import System

HERE = os.path.dirname(os.path.abspath(__file__))

def main():
    m = build(N=4)
    alg, par = m['alg'], m['par']
    e, T = m['eps'], m['T']
    P = m['P'](e, T); Ps = m['Pstar'](e, T); Ke = m['Kt'](e)
    M, Q = m['M'], m['Q']
    PM = mat_mul(Ps, M)
    Psi2 = tensor_act(PM, 1, P, par)                 # (P*M x 1 + 1 x P*M) P, odd kernel
    fs = test_functionals(alg, random.Random(7), n_cubic=500, n_quartic=3, terms=40)
    out = {}
    for name in ('II', 'I'):
        sysm = System(m, name)
        r = {}
        for (xn, X, pX) in (('M', M, 0), ('Q', Q, 1), ('P*M', PM, 1)):
            for (pn, phi, pphi) in (('P', P, 0), ('K_eps', Ke, 1), ('Psi', Psi2, 1)):
                Xphi = tensor_act(X, pX, phi, par)
                s = -1 if (pX and pphi) else 1
                for c in (1, -1):
                    err = max(alg.dist(alg.add(sysm.dphi(phi, sysm.D(X, pX, f)), (-s, sysm.D(X, pX, sysm.dphi(phi, f)))),
                                       alg.scale(c, sysm.dphi(Xphi, f))) for f in fs)
                    r['[d_%s, D_%s] - (%+d) d_{(X x1 + 1x X)%s}' % (pn, xn, c, pn)] = str(err)
                r['|(X x1 + 1x X)%s| for X=%s' % (pn, xn)] = str(max([abs(v) for v in Xphi.values()] + [Fr(0)]))
        out['system ' + name] = r
    return out

if __name__ == '__main__':
    res = main()
    with open(os.path.join(HERE, 'uniform_c07_results.json'), 'w') as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))
