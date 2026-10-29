"""Validation 3 (S^(1)(v) insertion, composite A|_V--c line): dilatation v = 2E, sector (A|_Wperp field leg at p,
A^*|_Wperp antifield leg), flow 1PI graphs with S^(1)(Delta) on the cycle (engine A Part D 'Pi_int'):
alpha delta(Sigma) + beta^mu d_mu delta(Sigma) with beta/alpha = (10/3 p)/(8/3) = (5/4) p, finite and local at T = infinity.
v = 2E = sum_m 2 y^m e_m  ->  -2i sum_m d/dP_m [graph with u = e_m e^{iP.y}] at P = 0."""
import json, time, numpy as np
import kn_engine as E
t0 = time.time()
jet = E.jet_order1(4)
p1 = np.array([1, 2, 0, 1])
a = np.zeros((10, 3), complex); a[4, 0] = 1            # A|_W field leg, colour 0, direction 4
bst = np.zeros((10, 3), complex); bst[4, 0] = 1        # A^*|_W antifield leg
I4 = np.eye(4)
legA = E.bg_insertion('A', a, p1, 'a')
legS = E.bg_insertion('As', bst, -p1, 'b', mom_lin=-I4)
out = {}
for eps in [1e-3, 3e-4, 1e-4]:
    F = np.zeros((4, jet.nb), complex)
    for m in range(4):
        s1 = E.s1v_insertion(I4[m], 'g', pin_lin=I4)
        cyc = E.enumerate_cycles([legA, legS, s1], jet, ['b', 'g'], special=False)
        for perm, Ls, cy in cyc:
            if cy.is_zero():
                continue
            asm = E.Assembled(cy, cy.numerator())
            F[m] += E.integrate_flow_simplex(asm, eps, 1e10, per_decade=10) * E.prefactor(cy)
    alpha = -2j * sum(F[m][jet.lin[m]] for m in range(4))
    beta = -2j * np.array([F[m][0] for m in range(4)])
    out['eps=%g' % eps] = {'alpha': [alpha.real, alpha.imag], 'beta': [[x.real, x.imag] for x in beta],
                           'beta/(alpha p1)': [(beta[m] / (alpha * p1[m])).real if p1[m] else None for m in range(4)]}
out['expected beta/(alpha p1)'] = 1.25
out['time_s'] = time.time() - t0
print(json.dumps(out, indent=1))
json.dump(out, open('kn_test_dil_piint.json', 'w'), indent=1)
