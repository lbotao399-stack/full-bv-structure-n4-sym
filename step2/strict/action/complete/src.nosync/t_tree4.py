import numpy as np, random
import cb_common as B, cb_tree as T
Q, C, V = B.Q, B.C, B.V
random.seed(7)
Pi = np.diag([(-1) ** p for p in C.PARITY]).astype(complex)
N = Q.dsum((2, Q.F_AAs_V), (2, Q.F_AAs_p), (-4, Q.F_lamslam), (-1, Q.F_csc))
F = Q.dsum((1, Q.F_dcA), (0.5, Q.dens(1, ('c', 'A'), lambda s: s[0]['u'][0] * s[1]['u'][2] * 3)))
def odd_br(g, f, variant):
    cache = {}
    def h(legs):
        tot = 0j
        for j, leg in enumerate(legs):
            key = tuple(float(x) for x in leg['p'])
            if key not in cache:
                M = -T.K0 @ T.hessian(g, leg['p'], [-x for x in leg['p']]).T
                cache[key] = Pi @ M if variant == 'PiL' else M
            img = cache[key] @ T.leg_vec(leg)
            nb = sum(V.PAR[l['sp']] for l in legs[:j]); na = sum(V.PAR[l['sp']] for l in legs[j + 1:])
            sg = 1 if variant == 'PiL' else (-1) ** (nb + na)
            for nl in T.split_leg(img, leg['p']):
                new = list(legs); new[j] = nl
                tot += sg * f(new)
        return tot
    return h
QN = Q.apply_Q(N); QF = Q.apply_Q(F)
rhs = Q.dsum((-1, T.bracket(QN, 0, F)), (1, T.bracket(QF, 0, N)))
cfgs = []
for spcs in (['c', 'A'], ['A', 'c'], ['c', 'As'], ['As', 'c'], ['A', 'A'], ['c', 'c'], ['lam', 'lam']):
    for _ in range(2):
        legs = [Q.rand_leg(s) for s in spcs]; legs[1]['p'] = [-x for x in legs[0]['p']]
        cfgs.append(legs)
for variant in ('plain', 'PiL'):
    lhs = Q.apply_Q(odd_br(N, F, variant))
    print(variant, 'C-rule residual', max(abs(lhs(l) - rhs(l)) for l in cfgs), ' size', max(abs(rhs(l)) for l in cfgs),
          ' {N,F} size', max(abs(odd_br(N, F, variant)(l)) for l in cfgs))
