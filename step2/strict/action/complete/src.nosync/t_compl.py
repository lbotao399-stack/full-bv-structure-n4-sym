import sys, time
sys.argv = ['x', 'B', '+']
import numpy as np, xieta as XE, eta1 as H, dtree as DT, fitlib as FL, bilin as BL
from fractions import Fraction as Fr
Q, TT = XE.Q, XE.TT
rng = np.random.default_rng(3)
sec = ('lam+', 'ls+')
confs = []
for _ in range(14):
    p = [Fr(int(x)) for x in rng.integers(-3, 4, size=4)]
    confs.append(dict(la=H.legs_at_of(sec[0], sec[1], H.rand_u(sec[0], rng), H.rand_u(sec[1], rng), p), ctx=dict(xi=H.rand_u('lam+', rng), eta=H.rand_u('ls+', rng))))
def target(c):
    xi, eta = c['ctx']['xi'], c['ctx']['eta']
    S2x = Q.dsum(*[(1, v) for k, v in H.xi_terms(xi).items()])
    return TT.br(S2x, 0, H.S_eta(eta, ('V',)))
tv = [z for c in confs for z in DT.canon(target(c), c['la'])]
for (mk_, mc) in ((3, 6), (4, 8), (5, 8)):
    t0 = time.time()
    (s0, P0), (s1, P1) = H.spec_of(sec[0]), H.spec_of(sec[1]); pr = H.projector((P0, P1)); slots = (s0, s1)
    cs = []
    for D, ex in ((0, ()), (1, ('e',))):
        for c in BL.bil2_cands(XE.RL, slots, D, ex, maxk=mk_, maxc=mc):
            if D == 0:
                cs.append(lambda ctx, c=c: DT.PFn(Q.dens(1, slots, lambda s: c.L(pr(s), ctx)), None))
            else:
                def mk(ctx, c=c):
                    ws = []
                    for mu in range(4):
                        e = np.zeros(4, complex); e[mu] = 1; cx = dict(ctx); cx['e'] = e
                        ws.append(Q.dens(1, slots, lambda s, cx=cx: c.L(pr(s), cx)))
                    return DT.PFn(None, ws)
                cs.append(mk)
    rows = [[z for c in confs for z in DT.canon(m(c['ctx']), c['la'])] for m in cs]
    keep = FL.independent_rows(rows)
    x = FL.solve_coords([rows[i] for i in keep], tv)
    print(mk_, mc, len(cs), 'indep', len(keep), 'target in span', x is not None, round(time.time() - t0))
