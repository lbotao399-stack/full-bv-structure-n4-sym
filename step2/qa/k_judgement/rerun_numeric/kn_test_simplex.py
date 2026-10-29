"""consistency: tau-analytic simplex integration vs direct log-grid integration on the A_W self-energy bubbles (n=2)
and a flow triangle with an S^(1)(v) insertion (n=3), moderate range."""
import numpy as np, time, json
import kn_engine as E
out = {}
jet = E.jet_trivial()
p = np.array([1, 2, 0, 1]); a = np.zeros((10, 3), complex); b = np.zeros((10, 3), complex)
a[4, 0] = 1; b[4, 0] = 1
bA = E.bg_insertion('A', a, p, 'a'); bB = E.bg_insertion('A', b, -p, 'b')
for perm, Ls, cy in E.enumerate_cycles([bA, bB], jet, [], special=False):
    if cy.is_zero(): continue
    asm = E.Assembled(cy, cy.numerator())
    t0 = time.time(); d = E.integrate_cycle(asm, 1e-3, 1e3, per_decade=16)[0]; t1 = time.time()
    s = E.integrate_flow_simplex(asm, 1e-3, 1e3, per_decade=16)[0]; t2 = time.time()
    out['bubble ' + ','.join(Ls)] = {'direct': [d.real, d.imag], 'simplex': [s.real, s.imag], 't': [t1 - t0, t2 - t1]}
# triangle with S^(1)(u) (u = e_1), jet order 1
jet1 = E.jet_order1(4); I4 = np.eye(4)
p1 = np.array([1, 2, 0, 1])
aa = np.zeros((10, 3), complex); aa[4, 0] = 1
bs = np.zeros((10, 3), complex); bs[4, 0] = 1
legA = E.bg_insertion('A', aa, p1, 'a'); legS = E.bg_insertion('As', bs, -p1, 'b', mom_lin=-I4)
s1 = E.s1v_insertion(I4[1], 'g', pin_lin=I4)
for perm, Ls, cy in E.enumerate_cycles([legA, legS, s1], jet1, ['b', 'g'], special=False):
    if cy.is_zero(): continue
    asm = E.Assembled(cy, cy.numerator())
    t0 = time.time(); d = E.integrate_cycle(asm, 1e-2, 1e2, per_decade=10); t1 = time.time()
    s = E.integrate_flow_simplex(asm, 1e-2, 1e2, per_decade=10); t2 = time.time()
    out['triangle ' + str(perm) + ','.join(Ls)] = {'direct': [abs(x) for x in d], 'simplex': [abs(x) for x in s],
                                                 'maxrel': float(np.max(np.abs(d - s)) / max(1e-30, np.max(np.abs(d)))), 't': [t1 - t0, t2 - t1]}
print(json.dumps(out, indent=1))
