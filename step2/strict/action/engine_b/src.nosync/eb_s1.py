"""Task S1 (two legs, along Sigma, weight 9/2, p = 0): Q on span{T0..T4} (degree -1) and span{g1..g4} (degree 0),
Prop at-xi2: kernel span{T1+T3, T2+T4}, H^0(Q) = 0.  Ranks over Q(i) on random conserving configurations.
(The weight-5 tables of the pair (Delta,P_b) are in B1_fit.json.)"""
import json, os, random
import numpy as np
import sympy as sp
import eb_qb as Q

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'out.nosync')


def cfgs(legsets, n=5):
    out = []
    for spcs in legsets:
        for _ in range(n):
            legs = [Q.rand_leg(s) for s in spcs]
            legs[-1]['p'] = [-sum(l['p'][m] for l in legs[:-1]) for m in range(4)]
            out.append(legs)
    return out


def mat(fns, cs):
    return sp.Matrix([[sp.nsimplify(complex(f(l)).real) + sp.I * sp.nsimplify(complex(f(l)).imag) for l in cs] for f in fns])


def run(log=print):
    random.seed(9)
    res = {}
    for chir, rng in (('+', range(8)), ('-', range(8, 16))):
        xi = np.zeros(16, complex)
        for c in rng:
            xi[c] = random.randint(-2, 2) or 1
        T = Q.xi_terms(xi); g = Q.g_terms(xi)
        c1 = cfgs([['lam', 'As'], ['As', 'lam'], ['A', 'lams'], ['lams', 'A']])
        c0 = cfgs([['lam', 'A'], ['A', 'lam'], ['c', 'lams'], ['lams', 'c']])
        cm1 = cfgs([['lam', 'c'], ['c', 'lam']])
        Tn = ['T0', 'T1', 'T2', 'T3', 'T4']
        MT = mat([T[k] for k in Tn], c1)
        QT = mat([Q.apply_Q(T[k]) for k in Tn], c0)
        kerT = QT.T.nullspace()
        gn = ['g1', 'g2', 'g3', 'g4']
        Mg = mat([g[k] for k in gn], c0)
        Qg = mat([Q.apply_Q(g[k]) for k in gn], cm1)
        kerg = Qg.T.nullspace()
        res[f'xi{chir}'] = dict(dim_deg_minus1=MT.rank(), dim_deg0=Mg.rank(), rank_Q_deg_minus1=QT.rank(),
                                kernel_deg_minus1=[[str(sp.nsimplify(x)) for x in v] for v in kerT],
                                closed_deg0_dim=len(kerg), image_deg0_dim=QT.rank(),
                                H0=len(kerg) - QT.rank(), Hm1=len(kerT), exact_deg_minus1=0,
                                basis_labels=dict(deg_minus1=Tn, deg0=gn))
        log(f'S1 xi{chir}: {res[f"xi{chir}"]}')
    return res


if __name__ == '__main__':
    with open(os.path.join(OUT, 'S1.json'), 'w') as f:
        json.dump(run(), f, indent=1)
