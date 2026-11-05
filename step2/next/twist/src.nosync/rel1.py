import env, numpy as np, itertools
import freeq as F, qa_twist as TW, ea_core as E, qa_e0 as E0
O = E.OFS; NS = E.NS
T = TW.Twist(); xi = T.xi
kap = F.kappa_S2(xi)
def Mp(p): return F.M_free(xi, np.asarray(p, complex), kap)
def row(slot, vec):
    w = np.zeros(NS, complex); n = len(vec); w[O[slot]:O[slot]+n] = vec; return w
# letters as functionals: gamma: w on A_perp slots 4..9 (gam_w); beta: dual rows on lam; dc: dual rows on lam
rng = np.random.default_rng(5)
P = [rng.integers(-3, 4, 4).astype(complex) for _ in range(3)]
def show(w, p, tol=1e-9):
    v = w @ Mp(p)
    out = {}
    for s in ('c','A','As','l','ls','cs'):
        blk = v[O[s]:O[s]+E.SZ[s]]
        if np.abs(blk).max() > tol: out[s] = blk
    return out
# R3: closedness of gamma, beta, dc letters
for nm, ws in (('gamma', [row('A', np.concatenate([np.zeros(4), g])) for g in T.gam_w]),
               ('beta', [row('l', b) for b in T.beta_z]), ('dc(lambar z)', [row('l', d) for d in T.dc_z])):
    for k, w in enumerate(ws):
        print(nm, k, 'Q w components:', {s: np.round(v, 6).tolist() for s, v in show(w, P[0]).items()})
# R1: Q(iota_n A) for n = d1, d2 (holomorphic) and exact directions
for lab, wv in (('d1', T.nv[(0, 1)]), ('d2', T.nv[(-1, 0)]), ('e(1,0)', T.nv[(1, 0)]), ('e(0,-1)', T.nv[(0, -1)])):
    w = row('A', np.concatenate([wv, np.zeros(6)]))
    for p in P[:2]:
        v = w @ Mp(p)
        cpart = v[O['c']]; lpart = v[O['l']:O['l']+16]
        # decompose lpart on dual rows: coefficients = lpart @ W[k] (dual rows x weight vectors = delta)
        coefs = {k: complex(np.round(lpart @ T.W[k], 9)) for k in T.W if abs(lpart @ T.W[k]) > 1e-9}
        print(lab, 'p', p.real.tolist(), 'c-coef', np.round(cpart, 9), ' (-i p.n =', np.round(-1j * (p @ wv), 9), ') lam coefs', coefs)
