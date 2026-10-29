"""Direct finite-eps evaluation of the order-hbar part of Q W(P(eps,infty), I - I^CT[eps]) in the sectors
(c, A|_Wperp, A|_Wperp) and (c, c, A^*|_V), with Richardson-type extrapolation eps -> 0
(refute_cubic_1; independent code: ft.py, loop.py).

  Q W at a configuration e:  W(e + r Q e) = W(e) + r (QW)(e)      (Q on legs: prop c1-Q; r passes to the front)
  order-hbar W on a leg set: 1PI one-loop (triangles, quartic bubbles)
                             + self-energy subgraph joined to I_3 by P(eps,infty)  [W_1^{1PI}(e_l, chi_jm)]
                             - I^CT_2(e_l, chi_jm) - I^CT_3(e)      [and optionally - F_2(e_l, chi_jm)]
  I^CT_2, I^CT_3: Step-1 counterterms (thm c4-quadratic, thm c4-log); F_2 = (27/16) unit int(d^*A_V, d^*A_V).
"""
import json
import sys
import time
import numpy as np
import ft
import loop

SEED = int(sys.argv[1]) if len(sys.argv) > 1 else 11
NCONF = int(sys.argv[2]) if len(sys.argv) > 2 else 2
rng = np.random.default_rng(SEED)
quad = loop.Quad(nt=int(sys.argv[3]) if len(sys.argv) > 3 else 80)
EPS = 2e-4 * 0.5 ** np.arange(9)


def rand_vec(n):
    return rng.normal(size=n) + 1j * rng.normal(size=n)


def A_leg(p, comps, mono=()):
    f = np.zeros((10, 3), complex)
    for M in comps:
        f[M] = rand_vec(3)
    return ft.Leg('A', f.reshape(30), p, mono)


def As_leg(p, comps, mono):
    f = np.zeros((10, 3), complex)
    for M in comps:
        f[M] = rand_vec(3)
    return ft.Leg('As', f.reshape(30), p, mono)


def fit(vals, eps, with_inv=False):
    cols = [np.ones_like(eps), np.log(eps), eps * np.log(eps), eps, eps ** 2 * np.log(eps), eps ** 2]
    if with_inv:
        cols = [1 / eps] + cols
    B = np.stack(cols, 1)
    sol, *_ = np.linalg.lstsq(B, vals, rcond=None)
    return sol


def cnum(z):
    z = complex(z)
    return [z.real, z.imag]


def run_sector(name, legs, target_mono, Gfun, Gslots, Gmono, claim_total, claim_alpha):
    t0 = time.time()
    rec = {'sector': name}
    # tree level: Q I_3 = 0
    tree = 0
    legsets = []
    for j in range(3):
        QL = ft.Q_leg(legs[j])
        if QL is None:
            continue
        S = list(legs)
        S[j] = QL
        legsets.append((j, S))
        tree += ft.multilinear(ft.CUBIC, S).get(target_mono, 0)
    G = ft.multilinear([(Gslots, Gfun)], legs)[Gmono]
    rec['G_value'] = cnum(G)
    rec['tree_QI3'] = cnum(tree)
    tls = [(j, loop.ThreeLeg(S, rng)) for j, S in legsets]
    rec['poly_fit_max_rel_err'] = max(t.max_fit_err() for _, t in tls)
    comps = {}
    totals, totalsF = [], []
    for e in EPS:
        tot, totF = 0, 0
        for j, tl in tls:
            c = tl.components(e, quad, with_F=True)
            for k, d in c.items():
                v = d.get(target_mono, 0)
                comps.setdefault(f'leg{j+1}:{k}', []).append(v)
                if k != 'tree_F2':
                    tot += v
                totF += v
        totals.append(tot)
        totalsF.append(totF)
    totals = np.array(totals)
    totalsF = np.array(totalsF)
    sol = fit(totals, EPS)
    solF = fit(totalsF, EPS)
    sol_inv = fit(totals, EPS, with_inv=True)
    # stability: fits on subsets
    sub = []
    for drop in range(3):
        idx = np.arange(drop, len(EPS))
        sub.append(fit(totals[idx], EPS[idx])[0] / (ft.UNIT * G))
    rec['raw_total_over_unitG'] = [[float(e), cnum(v / (ft.UNIT * G))] for e, v in zip(EPS, totals)]
    rec['limit_over_unitG'] = cnum(sol[0] / (ft.UNIT * G))
    rec['log_eps_coeff_over_unitG'] = cnum(sol[1] / (ft.UNIT * G))
    rec['eps_inv_coeff_over_unitG(fit incl. 1/eps)'] = cnum(sol_inv[0] / (ft.UNIT * G))
    rec['limit_with_F2_over_unitG'] = cnum(solF[0] / (ft.UNIT * G))
    rec['log_eps_coeff_with_F2_over_unitG'] = cnum(solF[1] / (ft.UNIT * G))
    rec['limit_subset_fits'] = [cnum(x) for x in sub]
    rec['claim_total(alpha+tree)'] = claim_total
    rec['claim_alpha'] = claim_alpha
    # component diagnostics: log eps and eps^-1 and finite part of each piece
    cd = {}
    for k, vs in comps.items():
        s2 = fit(np.array(vs), EPS, with_inv=True)
        cd[k] = {'eps^-1': cnum(s2[0] / (ft.UNIT * G)), 'log': cnum(s2[2] / (ft.UNIT * G)), 'finite': cnum(s2[1] / (ft.UNIT * G))}
    rec['components_over_unitG'] = cd
    rec['seconds'] = time.time() - t0
    return rec


out = {'seed': SEED, 'eps': EPS.tolist(), 'unit': ft.UNIT, 'configs': []}
L32 = np.log(1.5)
claim_W_alpha = 11 / 4 + 19 / 4 * L32
claim_W_total = claim_W_alpha + 27 / 8
claim_S_alpha = 0.0
claim_S_total = -27 / 4
for n in range(NCONF):
    p1 = rng.normal(size=4)
    p2 = rng.normal(size=4)
    p3 = -p1 - p2
    conf = {'p1': p1.tolist(), 'p2': p2.tolist(), 'p3': p3.tolist()}
    legsW = [ft.Leg('c', rand_vec(3), p1, (1,)), A_leg(p2, range(4, 10)), A_leg(p3, range(4, 10))]
    conf['cAWAW'] = run_sector('(c,A_W,A_W)', legsW, (0, 1), ft.T3_G3, ('c', 'A', 'A'), (1,), claim_W_total, claim_W_alpha)
    print(json.dumps({k: conf['cAWAW'][k] for k in ('tree_QI3', 'limit_over_unitG', 'log_eps_coeff_over_unitG', 'limit_with_F2_over_unitG', 'limit_subset_fits', 'poly_fit_max_rel_err', 'seconds')}), flush=True)
    legsS = [ft.Leg('c', rand_vec(3), p1, (1,)), ft.Leg('c', rand_vec(3), p2, (2,)), As_leg(p3, range(4), (3,))]
    conf['ccAsV'] = run_sector('(c,c,A*_V)', legsS, (0, 1, 2, 3), ft.T3_G4, ('c', 'c', 'As'), (1, 2, 3), claim_S_total, claim_S_alpha)
    print(json.dumps({k: conf['ccAsV'][k] for k in ('tree_QI3', 'limit_over_unitG', 'log_eps_coeff_over_unitG', 'limit_with_F2_over_unitG', 'limit_subset_fits', 'poly_fit_max_rel_err', 'seconds')}), flush=True)
    out['configs'].append(conf)
    json.dump(out, open(f'main_cubic_seed{SEED}.json', 'w'), indent=1)
print('claims: (c,A_W,A_W) total', claim_W_total, 'alpha_3', claim_W_alpha, '; (c,c,A*_V) total', claim_S_total, 'alpha_4', claim_S_alpha)
