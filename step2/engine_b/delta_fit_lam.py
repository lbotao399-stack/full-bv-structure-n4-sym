"""fit of the lambda-sector log coefficient beta^mu of class (i) to the basis
 B1(+-) = p^mu (P+- u_f)^T Mp u_a,  B2(+-) = (Rb[mu] R(p) P+- u_f)^T Mp u_a  (R(p) = sum p_nu R[nu]: Sigma -> Sigma^*)."""
import json, sympy as sp
import delta_engine as D
from realization import R, Rb
d = json.load(open('delta_part1_sing_lam.json'))['lam']
P = list(sp.symbols('P0:4'))
b = sp.symbols('b1p b1m b2p b2m')
eqs = []
for s in d:
    uf = sp.Matrix([sp.sympify(x) for x in s['uf']]); ua = sp.Matrix([sp.sympify(x) for x in s['ua']])
    ufp = sp.Matrix([uf[i] if i < 8 else 0 for i in range(16)]); ufm = uf - ufp
    Rp = sum((P[n] * R[n] for n in range(4)), sp.zeros(16, 16))
    for mu in range(4):
        beta = sp.sympify(s['beta'][mu])
        model = 0
        for (cb1, cb2, u) in ((b[0], b[2], ufp), (b[1], b[3], ufm)):
            model += cb1 * P[mu] * (u.T * D.V.Mp * ua)[0] + cb2 * ((Rb[mu] * Rp * u).T * D.V.Mp * ua)[0]
        diff = sp.expand(beta - model)
        eqs += sp.Poly(diff, *P).coeffs()
    # alpha must vanish
    assert sp.sympify(s['alpha']) == 0
sol = sp.solve(eqs, b, dict=True)
print(sol)
json.dump({'basis': 'beta^mu = b1(+-) p^mu (P+- lambda, lambda*) + b2(+-) (rho(e_mu) rho(p) P+- lambda, lambda*); P+ = S_+ (x) W part',
           'solution': {str(k): str(v) for k, v in sol[0].items()}, 'n_equations': len(eqs)}, open('delta_fit_lam.json', 'w'), indent=1)
