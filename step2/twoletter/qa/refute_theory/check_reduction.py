"""Referee check of Prop prop:tl-reduction (ii)-(iii) (transversality of homogeneous solutions in an output A|_V).

The homogeneous equation of (i) is
    K_r(Q X_a, X_b) + (-1)^{|X_a|} K_r(X_a, Q X_b) - Q K_r(X_a, X_b) = 0 .
For X_a = A|_V the letter term is Q X_a = -d c: a kernel of the pair (dc, X_b) with the SAME number r of outputs.
(ii) drops this term.  Minimal test: X_a = A_mu|_V, X_b = A_perp,I, outputs (A_rho|_V(r1), A_perp,J(r2)), W = 2,
log part of momentum degree 0; the coupled unknown is K(d_mu c, A_perp,I) with outputs (c(r1), A_perp,J(r2)), degree 1.
Most general O(4) x O(6) covariant forms (per colour word):
    K(A_mu, A_perp,I)|_{A_rho A_perp,J} = t  delta_{mu rho} delta_{IJ}
    K(d_mu c, A_perp,I)|_{c A_perp,J}   = (al r1_mu + be r2_mu) delta_{IJ}
Q on the output A_rho(r1): replaced by -i r1_rho c(r1).  sigma = +-1 collects the Koszul sign of the ghost output.
Then the homogeneous equation reads  -(al r1 + be r2) - sigma(-i t r1) = 0  for all r1, r2.
(ii) asserts that every solution is transverse in the A|_V output (depends on p ^ a): here that means t = 0.
Also the own equation of the pair (c, A_perp) at r = 2 is checked to be empty (letters and outputs Q-inert, no r < 2
kernels), so K(dc, A_perp) is not fixed by lower data.
Engine A's exact invariant form of the actual W_0 kernel K(A|_V, A_perp)|_{A|_V A_perp} is quoted from
engine_a/kernels.json (read-only) to show the non-transverse delta(a,o1) structure in the computed kernel.
"""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
t, al, be = sp.symbols('t al be')
sigma = sp.symbols('sigma')
r1 = sp.symbols('r1_0:4'); r2 = sp.symbols('r2_0:4')
res = {}
for sg in [1, -1]:
    eqs = []
    for mu in range(4):
        expr = -(al * r1[mu] + be * r2[mu]) - sg * (-sp.I * t * r1[mu])
        poly = sp.Poly(sp.expand(expr), *(r1 + r2))
        eqs += list(poly.coeffs())
    sol = sp.solve(eqs, [al, be], dict=True)
    res[str(sg)] = dict(solution=str(sol), t_free=bool(sol and t not in sol[0]))
eng = None
try:
    d = json.load(open(os.path.join(HERE, '..', '..', 'engine_a', 'kernels.json')))
    eng = d['invariant_forms_bosonic_W0']['forms']['AV|AW']
except Exception as e:
    eng = 'unavailable: %r' % e
out = dict(homogeneous_solutions=res,
           conclusion='for either Koszul sign the homogeneous system has the one-parameter family al = +-i t, be = 0 with '
                      't arbitrary: a non-transverse (delta_{mu rho}) dependence on the output A|_V that (ii) excludes',
           engine_A_W0_invariant_forms_AV_AW=eng)
json.dump(out, open(os.path.join(HERE, 'check_reduction.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
