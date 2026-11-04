"""Independent toy checks (no engine code), 30-digit mpmath.
(1) Two-edge cycle with one special edge K_eps - K_0 (factor e^{-eps l^2} - 1) and one edge P(eps, T) (Schwinger time in
    [eps, T]), constant vertex numerator, 4d Gaussian loop integral done in closed form:
       16 pi^2 G(eps, T; q) = int_eps^T ds [ (eps + s)^{-2} exp(-eps s q^2/(eps + s)) - s^{-2} ].
    Prediction at T = oo (by hand, and the qsk closed form B2): 16 pi^2 G = -1/(2 eps) - (3/8) q^2 + O(eps), no log eps.
(2) Finite T: T d/dT G = T[(eps+T)^{-2} e^{-eps T q^2/(eps+T)} - T^{-2}] = O(eps) at fixed T; with a position weight s^k
    on the pinned edge (k = 1, 2 mimic the two r-derivatives of K_b) T d/dT stays O(eps) at fixed T but grows like T^{k-1}.
(3) Bridge cancellation behind finding QA-F1: a divergent subgraph S(eps) = a_1/eps + a_0 + l log eps (here the
    one-letter toy S = G above) on a tree edge carrying (e^{-eps p^2} - 1): the bare bridge tends to -p^2 a_1 != 0, the
    counterterm-joined bridge to +p^2 a_1, their sum to 0."""
import mpmath as mp, json
mp.mp.dps = 30
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/strict/operator'

def G(eps, q2, T=mp.inf, k=0):
    f = lambda s: s**k * ((eps + s)**-2 * mp.e**(-eps * s * q2 / (eps + s)) - s**-2)
    return mp.quad(f, [eps, 10 * eps, 1000 * eps, T]) if T == mp.inf else mp.quad(f, [eps, T])

out = {}
q2 = mp.mpf(3)
rows = []
for e in [mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')]:
    val = G(e, q2)
    rows.append({'eps': str(e), '16pi2 G': mp.nstr(val, 20), 'G + 1/(2eps)': mp.nstr(val + 1 / (2 * e), 15),
                 'target -(3/8) q^2': mp.nstr(-mp.mpf(3) / 8 * q2, 15)})
out['toy1_T_infinity'] = rows
rows = []
for T in [1, 4, 16]:
    for e in [mp.mpf('1e-3'), mp.mpf('1e-4')]:
        T = mp.mpf(T)
        d = {k: mp.nstr(T * T**k * ((e + T)**-2 * mp.e**(-e * T * q2 / (e + T)) - T**-2), 12) for k in (0, 1, 2)}
        rows.append({'T': str(T), 'eps': str(e), 'T dG/dT for position weight s^k, k=0,1,2': d})
out['toy2_finite_T'] = rows
p2 = mp.mpf(2)
rows = []
for e in [mp.mpf('1e-3'), mp.mpf('1e-4'), mp.mpf('1e-5')]:
    S = G(e, q2)
    a1 = -mp.mpf(1) / 2            # pole coefficient of S (exact, from (1))
    bare = (mp.e**(-e * p2) - 1) * S
    ctj = (mp.e**(-e * p2) - 1) * (-a1 / e)
    rows.append({'eps': str(e), 'bare bridge': mp.nstr(bare, 12), 'target -p^2 a_1': mp.nstr(-p2 * a1, 12),
                 'counterterm bridge': mp.nstr(ctj, 12), 'sum': mp.nstr(bare + ctj, 12)})
out['toy3_bridge'] = rows
json.dump(out, open(f'{ROOT}/qa/out.nosync/toy_special_cycle.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
