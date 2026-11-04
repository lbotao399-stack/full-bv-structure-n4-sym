import numpy as np, random
from fractions import Fraction as Fr
import cb_common as B, bilin as BL, fitlib as FL, cb_pairs as CP
Q = B.Q
RL = BL.Real(B.R, B.RB, B.MP)
random.seed(1)
cands = BL.spinor_cands(RL, 'llS', 1) + BL.spinor_cands_b(RL, 'llS', 1)
for sec, (c1, c2, r0, r1) in (('l+ls+ pm', '+', '-', B.LAMP, B.LAMS_P), ('l-ls- pm', '+', '-', B.LAMM, B.LAMS_M)) if False else []:
    pass
def test(c1, c2, r0, r1, n=24):
    confs = []
    for _ in range(n):
        p = [Fr(random.randint(-3, 3)) for _ in range(4)]
        legs = [dict(sp='lam', u=CP.rpol('lam', r0), p=p), dict(sp='lams', u=CP.rpol('lams', r1), p=[-x for x in p])]
        confs.append((legs, dict(xi1=CP.rspin(CP.SPIN[c1]), xi2=CP.rspin(CP.SPIN[c2]))))
    rows = [[FL.to_fr(Q.dens(1, cd.slots, lambda s, cd=cd, ctx=ctx: cd.L(s, ctx))(legs)) for legs, ctx in confs] for cd in cands]
    keep = FL.independent_rows(rows)
    return keep
print('pm l+ls+', test('+', '-', B.LAMP, B.LAMS_P))
print('pm l-ls-', test('+', '-', B.LAMM, B.LAMS_M))
print('pp l+ls-', test('+', '+', B.LAMP, B.LAMS_M))
print('pp l-ls+', test('+', '+', B.LAMM, B.LAMS_P))
