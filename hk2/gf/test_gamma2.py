import sympy as sp, models as Mo, bvjet
from test_equiv import odd, gens_even, r
g2 = sp.Symbol('gamma2')
for frame in ('C','Psi'):
    M = Mo.build_model(True); F = Mo.Fields(M, frame); S0 = Mo.S0(F)
    (a1,e1),(a2,e2) = odd['xi1'], odd['xi2']
    v, rot = Mo.bracket_odd_odd(a1,e1,a2,e2)
    s2 = Mo.S2_odd(F,a1,e1,a2,e2,gamma2=g2)
    lhs = M.plus(M.bracket(S0,s2), Mo.S1_even(F,v=v), Mo.S1_even(F,r=rot),
                 M.scale(M.bracket(Mo.S1_odd(F,a1,e1), Mo.S1_odd(F,a2,e2)), -1))
    eqs=set()
    for a in range(len(M.fields)):
        E = M.simplify(M.EL(lhs,a))
        for k,c in E.items():
            eqs.add(sp.expand(c))
    # coefficients are polynomials in x: split
    E2=set()
    for e in eqs:
        E2 |= set(sp.Poly(e,*bvjet.X).coeffs()) if e.free_symbols - {g2} else {e}
    print(frame, 'gamma2 solving the (xi,xi) arity-two equation:', sp.solve(list(E2), g2))
