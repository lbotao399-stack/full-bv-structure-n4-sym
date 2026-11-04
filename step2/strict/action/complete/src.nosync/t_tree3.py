import numpy as np, random, itertools
import cb_common as B, cb_tree as T
Q, C, V = B.Q, B.C, B.V
random.seed(5)
Pi = np.diag([(-1) ** p for p in C.PARITY]).astype(complex)
Sfree = Q.dsum((-0.25, Q.F_dAdA_V), (-0.25, Q.F_dAdA_p), (0.5, Q.F_ldl), (1, Q.F_dcAs))
N = Q.dsum((2, Q.F_AAs_V), (2, Q.F_AAs_p), (-4, Q.F_lamslam), (-1, Q.F_csc))
QN = Q.apply_Q(N)
def mk(variant):
    def Mm(g, p):
        mp = [-x for x in p]
        H = T.hessian(g, p, mp)
        base = {'HT': H.T, 'H': H}[variant[0]]
        M = -T.K0 @ base
        if variant[1] == 'PiL': M = Pi @ M
        if variant[1] == 'PiR': M = M @ Pi
        if variant[1] == 'PiLR': M = Pi @ M @ Pi
        return M
    return Mm
cfgs = []
for spcs in (['A', 'A'], ['lam', 'lam'], ['c', 'As'], ['As', 'c'], ['A', 'c'], ['c', 'A'], ['lam', 'A'], ['As', 'A']):
    for _ in range(2):
        legs = [Q.rand_leg(s) for s in spcs]; legs[1]['p'] = [-x for x in legs[0]['p']]
        cfgs.append(legs)
for variant in itertools.product(('HT', 'H'), ('', 'PiL', 'PiR', 'PiLR')):
    Mm = mk(variant)
    cache = {}
    for (a, b) in ((0, 0), (1, 0), (0, 1), (1, 1)):
        def h(legs):
            tot = 0j
            for j, leg in enumerate(legs):
                key = tuple(float(x) for x in leg['p'])
                if key not in cache: cache[key] = Mm(N, leg['p'])
                img = cache[key] @ T.leg_vec(leg)
                nb = sum(V.PAR[l['sp']] for l in legs[:j]); na = sum(V.PAR[l['sp']] for l in legs[j + 1:])
                for nl in T.split_leg(img, leg['p']):
                    new = list(legs); new[j] = nl
                    sg = (-1) ** (a * nb + b * na) * (-1) ** (V.PAR[leg['sp']] if False else 0)
                    tot += sg * Sfree(new)
            return tot
        for s in (1, -1):
            w = max(abs(QN(l) - s * h(l)) for l in cfgs)
            if w < 1e-9: print('match', variant, (a, b), s)
print('done')
