"""Validation 1: Step-1 arity-zero obstruction, sector (c, A|_V), five special-edge bubbles (ch04; engine A
logs/r2_T2_routeA.json): eps^-1 total -9/2, finite total -27/8 (coefficients of i(gamma, iota_p a)/eps and
i p^2 (gamma, iota_p a), unit hbar tr_g(ad ad)/(16 pi^2 (X,Y)))."""
import json, time, numpy as np
import kn_engine as E
t0 = time.time()
jet = E.jet_trivial()
p = np.array([1, 2, 0, 1]); gam = np.array([1, -2, 1], complex)
a = np.zeros((10, 3), complex); a[:4] = np.array([[1, 0, 2], [-1, 1, 0], [0, 2, 1], [1, 1, -1]])
S = np.einsum('a,m,ma->', gam.real, p, a[:4].real)
bc = E.bg_insertion('c', gam, p, 'c')
bA = E.bg_insertion('A', a, -p, 'A')
cyc = E.enumerate_cycles([bc, bA], jet, ['c'], special=True)
res = {}
epss = [2e-4, 5e-4, 1e-3, 2e-3, 4e-3]
T = 1e7
tot = np.zeros(len(epss), complex)
for perm, Ls, cy in cyc:
    if cy.is_zero():
        continue
    poly = cy.numerator()
    asm = E.Assembled(cy, poly)
    vals = []
    for eps in epss:
        v = E.integrate_cycle(asm, eps, T, per_decade=14)[0] * E.prefactor(cy)
        vals.append(v)
    vals = np.array(vals)
    key = '%s|%s' % (''.join('cA'[i] for i in perm), ','.join(Ls))
    # fit v = A/eps + B + C eps + Dd eps log eps
    M = np.array([[1 / e, 1, e, e * np.log(e)] for e in epss])
    coef = np.linalg.lstsq(M, vals, rcond=None)[0]
    res[key] = {'eps_inv_per_iS_cadj': (coef[0] / (1j * S * (-2))).real, 'fin_per_ip2S_cadj': (coef[1] / (1j * (p @ p) * S * (-2))).real}
    tot += vals
M = np.array([[1 / e, 1, e, e * np.log(e)] for e in epss])
coef = np.linalg.lstsq(M, tot, rcond=None)[0]
res['TOTAL'] = {'eps_inv_per_iS_cadj': (coef[0] / (1j * S * (-2))).real, 'fin_per_ip2S_cadj': (coef[1] / (1j * (p @ p) * S * (-2))).real,
                'imag_parts': [(coef[0] / (1j * S * (-2))).imag, (coef[1] / (1j * (p @ p) * S * (-2))).imag]}
res['expected'] = {'eps_inv_total': -4.5, 'finite_total': -27 / 8, 'per_graph_engineA': 'sAV|AV -9/4,-5/8; sAV|gh 1/8,-1/12; sAP|AP -9/2,-5/4; sl|l 2,-4/3; sc|gh 1/8,-1/12'}
res['time_s'] = time.time() - t0
print(json.dumps(res, indent=1))
json.dump(res, open('kn_test_step1_cA.json', 'w'), indent=1)
