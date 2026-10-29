"""Validation 2: Step-1 A|_Wperp self-energy (flow, no special edge), cross coefficient of W_1 on
A = a e^{ipx} + b e^{-ipx}: expected (manuscript thm:c4-quadratic, engine A Part D) in unit = -2 x (hbar/16pi^2):
Sigma/unit = (a,b) [ (9/2)/eps + p^2 (43/8 - 2 log(eps p^2) - 2 gamma_E - 2 log 2) ]."""
import json, time, numpy as np
import kn_engine as E
t0 = time.time()
jet = E.jet_trivial()
p = np.array([1, 2, 0, 1]); a = np.zeros((10, 3), complex); b = np.zeros((10, 3), complex)
a[4, 0] = 1; a[5, 1] = 2; b[4, 0] = 1; b[5, 1] = -1; b[6, 2] = 3
ab = np.sum(a * b).real
bA = E.bg_insertion('A', a, p, 'a')
bB = E.bg_insertion('A', b, -p, 'b')
b4 = E.bg_insertion('A', a, np.zeros(4), 'a', quartic_with=('A', b, 'b'))
cycles = E.enumerate_cycles([bA, bB], jet, [], special=False) + E.enumerate_cycles([b4], jet, [], special=False)
epss = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3]
T = 1e8
tot = np.zeros(len(epss), complex)
per = {}
for perm, Ls, cy in cycles:
    if cy.is_zero():
        continue
    asm = E.Assembled(cy, cy.numerator())
    vals = np.array([E.integrate_cycle(asm, e, T, per_decade=14)[0] * E.prefactor(cy) for e in epss])
    per['%d|%s' % (cy.n, ','.join(Ls))] = vals
    tot += vals
M = np.array([[1 / e, np.log(e), 1, e, e * np.log(e)] for e in epss])
coef = np.linalg.lstsq(M, tot, rcond=None)[0]
unit = -2.0
p2 = p @ p
exp_c = ab * p2 * (43 / 8 - 2 * np.log(p2) - 2 * np.euler_gamma - 2 * np.log(2))
res = {'eps_inv/unit/(a,b)': (coef[0] / unit / ab).real, 'expected_eps_inv': 4.5,
       'log_eps/unit/(a,b)p^2': (coef[1] / unit / ab / p2).real, 'expected_log_eps': -2.0,
       'finite/unit': (coef[2] / unit).real, 'expected_finite': exp_c,
       'per_graph_last_eps': {k: (v[-1] / unit).real for k, v in per.items()}, 'time_s': time.time() - t0}
print(json.dumps(res, indent=1))
json.dump(res, open('kn_test_step1_AW.json', 'w'), indent=1)
