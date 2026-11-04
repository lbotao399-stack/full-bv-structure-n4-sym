"""K_b post-checks (route given): G-pole along K_b (sum a_1 = (9/4){int(A,A), S(K_b)}, eq at-pole) on the stored
canonical pole data, and the exactness of the closed degree -1 directions beyond the sector parts S(K_b)_d
(Q images of the degree -2 candidates x x (A_V,c*), x x (lam*,lam*))."""
import sys, json, os
from fractions import Fraction as Fr
import numpy as np, sympy as sp
ROUTE = sys.argv[1]
sys.argv = [sys.argv[0], ROUTE]
import kb1 as K
H, Q, TT, DT, FL = K.H, K.Q, K.TT, K.DT, K.FL
raw = json.load(open(os.path.join(K.OUT, f'{ROUTE}_Kb_raw.json')))
FAA = Q.dens(1, ('A', 'A'), lambda s: np.dot(s[0]['u'], s[1]['u']))     # int (A, A)
ok, worst = True, 0
for key, recs in raw.items():
    t0, t1 = key.split('|')
    for r in recs:
        la = H.legs_at_of(t0, t1, H.cvec(r['u0']), H.cvec(r['u1']), [Fr(x) for x in r['p']])
        pred = DT.canon2(DT.tree_apply2(lambda d: TT.brd(FAA, 0, d), K.S_K(np.array(r['b'], complex))), la)
        pole = [H.lv_of(x) for x in r['pole']]
        for pv, lv in zip(pred, pole):
            tgt = (Fr(9, 4) * pv[0], Fr(9, 4) * pv[1])
            if any(v != (0, 0) for k, v in lv.items() if k != '1') or lv['1'] != tgt:
                ok = False
print('G-pole K_b (sum a_1 = 9/4 {int(A,A), S(K_b)}):', ok)
# exactness of the extra closed directions
fit = json.load(open(os.path.join(K.OUT, f'{ROUTE}_Kb_fit.json')))
rngm = np.random.default_rng(73)
cfm = []
for (a, b) in K.SECM1:
    for _ in range(8):
        p = [Fr(int(x)) for x in rngm.integers(-3, 4, size=4)]
        cfm.append(dict(la=H.legs_at_of(a, b, H.rand_u(a, rngm), H.rand_u(b, rngm), p), ctx=dict(b=rngm.integers(-2, 3, size=4).astype(complex))))
lab2mk = {lab: mk for s in K.SECM1 for (lab, mk) in K.cands(s)}
lab2mk.update({'S(K)_perp': lambda ctx: K.S_K(ctx['b'], ('perp',)), 'S(K)_lam': lambda ctx: K.S_K(ctx['b'], ('lam',)), 'S(K)_V': lambda ctx: K.S_K(ctx['b'], ('V',))})
Qd = lambda fn: DT.tree_apply2(TT.Qd, fn)
img = []
for s in K.SECM2:
    for (lab, mk) in K.cands(s):
        img.append([z for c in cfm for z in DT.canon2(Qd(mk(c['ctx'])), c['la'])])
keep = FL.independent_rows(img)
imrows = [img[i] for i in keep]
out = []
for fd in fit['F1_free_directions']:
    vec = None
    for lab, cf in fd.items():
        cfv = sp.nsimplify(cf)
        vals = [FL.gq(z) for c in cfm for z in DT.canon2(lab2mk[lab](c['ctx']), c['la'])]
        vec = [cfv * v for v in vals] if vec is None else [x + cfv * v for x, v in zip(vec, vals)]
    tv = [(Fr(str(sp.re(v))), Fr(str(sp.im(v)))) for v in vec]
    x = FL.solve_coords(imrows, tv) if imrows else None
    out.append(dict(direction=fd, exact=x is not None))
print('exact images', len(imrows)); [print(o) for o in out]
json.dump(dict(G_pole=ok, free_direction_exactness=out, n_exact_images=len(imrows)), open(os.path.join(K.OUT, f'{ROUTE}_Kb_post.json'), 'w'), indent=1)
vecs = []
for fd in fit['F1_free_directions']:
    vec = None
    for lab, cf in fd.items():
        cfv = sp.nsimplify(cf)
        vals = [FL.gq(z) for c in cfm for z in DT.canon2(lab2mk[lab](c['ctx']), c['la'])]
        vec = [cfv * v for v in vals] if vec is None else [x + cfv * v for x, v in zip(vec, vals)]
    vecs.append([(Fr(str(sp.re(v))), Fr(str(sp.im(v)))) for v in vec])
r_extra = len(FL.independent_rows(vecs[3:])); r_all = len(FL.independent_rows(vecs[3:] + imrows)); r_tot = len(FL.independent_rows(vecs + imrows))
print('rank extra', r_extra, 'rank extra+exact', r_all, 'rank all closed + exact', r_tot, ' => H^{-1} =', r_tot - len(imrows))
d = json.load(open(os.path.join(K.OUT, f'{ROUTE}_Kb_post.json')))
d.update(rank_extra=r_extra, rank_extra_plus_exact=r_all, H_minus1=r_tot - len(imrows))
json.dump(d, open(os.path.join(K.OUT, f'{ROUTE}_Kb_post.json'), 'w'), indent=1)
