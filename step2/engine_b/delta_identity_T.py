"""finite-eps test of the UV-boundary identity for the translation letter, sector (A|_Wperp): Q(class (i)_T)(eps) vs B1_T(eps)."""
import json, sys, time, math
import numpy as np, sympy as sp
import delta_engine as D, delta_canon as C, delta_numeric_i as NI
I = sp.I
u1 = [sp.Integer(x) for x in [0]*4 + [-1, -3, -3, 2, -1, 1]]; u2 = [sp.Integer(x) for x in [0]*4 + [1, 2, -1, -2, 3, -2]]
p1 = [sp.Integer(x) for x in (-1, 1, -3, -2)]
uu = float(sum(a*b for a, b in zip(u1, u2))); pn = np.array([float(x) for x in p1]); p2 = float(pn @ pn)
mu = 0
out = {}
for eps in [float(x) for x in sys.argv[1].split(',')]:
    # class (i), T letter along e_mu, kappa = 0, legs (A at p, A* at -p)
    legs = [dict(sp='A', u=u1, p=list(p1)), dict(sp='As', u=u2, p=[-x for x in p1]), dict(sp='xi', u=[sp.Integer(1)], p=[0]*4)]
    recs = D.enumerate_terms(legs, [('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC), ('SM', (2,), D.T_TERMS(mu))], 'P', [0, 1])
    gi = NI.evaluate_graphs(recs, {'v': lambda e: sp.expand(e)}, eps, n=40)['v']
    # B1, T letter, kappa = 0, legs (A at p, A at -p)
    legsB = [dict(sp='A', u=u1, p=list(p1)), dict(sp='A', u=u2, p=[-x for x in p1]), dict(sp='xi', u=[sp.Integer(1)], p=[0]*4)]
    recsB = []
    for part in ([('cubic', (0,), D.V.CUBIC), ('cubic', (1,), D.V.CUBIC), ('SM', (2,), D.T_TERMS(mu))], [('quartic', (0, 1), D.V.QUARTIC), ('SM', (2,), D.T_TERMS(mu))]):
        recsB += D.enumerate_terms(legsB, part, 'K', [0, 1])
    gb = NI.evaluate_special(recsB, {'v': lambda e: sp.expand(e)}, eps, n=40, n1=20)['v']
    out[str(eps)] = dict(classi_T_over_uu_p=gi / (uu * float(p1[mu])), Q_classi_T=-p2 * gi / (uu * float(p1[mu])),
                         B1_T_over_uu_p=gb / (uu * float(p1[mu])))
    out[str(eps)] = {k: [v.real, v.imag] for k, v in out[str(eps)].items()}
    print(eps, out[str(eps)], flush=True)
json.dump(out, open('delta_identity_T.json', 'w'), indent=1)
