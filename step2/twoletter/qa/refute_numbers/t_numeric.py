"""Direct finite-eps check (no momentum expansion): eps dw/deps = -eps sum_e (integrand at t_e = eps, other t >= eps),
evaluated by mpmath quadrature with the exact Gaussian loop integral and the full e^{-tau Phi};
residual r(eps) = eps w' + sum_m m a_m eps^{-m} - b must tend to 0 (a_m, b from the referee's exact singular parts)."""
import json, itertools
import mpmath as mp
import sympy as sp
from ref_core import *
mp.mp.dps = 18
def face_residuals(N, Rs, sing, epss):
    E = len(Rs)
    s = sp.symbols('t0:%d' % E)
    tau = sum(s)
    c = [-sum(s[j] * Rs[j][mu] for j in range(E)) / tau for mu in range(4)]
    tPhi = sum(s[j] * dot4(Rs[j], Rs[j]) for j in range(E)) - sum((sum(s[j] * Rs[j][mu] for j in range(E))) ** 2 for mu in range(4)) / tau
    integrand = 0
    lapN = sp.expand(N); j = 0
    while lapN != 0:
        integrand += lapN.subs({L[m]: c[m] for m in range(4)}, simultaneous=True) / (sp.factorial(j) * 4 ** j * tau ** (2 + j))
        lapN = sp.expand(laplacian(lapN)); j += 1
    integrand = integrand * sp.exp(-tPhi)        # times 1/(16 pi^2): units
    f = sp.lambdify(s, integrand, 'mpmath')
    out = []
    for eps in epss:
        eps = mp.mpf(eps)
        tot = 0
        for e in range(E):
            others = [k for k in range(E) if k != e]
            def g(*u):
                args = [None] * E; args[e] = eps
                for k, x in zip(others, u): args[k] = x
                return f(*args)
            # logarithmic variables t = eps e^x (multi-scale integrand: scales eps and 1/|R|^2)
            Lg = -mp.log(eps)
            pts = [0, Lg / 3, 2 * Lg / 3, Lg, Lg + 2, Lg + 5, mp.inf]
            if E == 2:
                val = mp.quad(lambda x: g(eps * mp.e ** x) * eps * mp.e ** x, pts)
            else:
                val = mp.quad(lambda x, y: g(eps * mp.e ** x, eps * mp.e ** y) * eps ** 2 * mp.e ** (x + y), pts, pts)
            tot += val
        epsw = -eps * tot
        r = epsw + sum(m * sp.N(v, 30) * eps ** (-m) for m, v in sing.items() if m != 'log') - sp.N(sing.get('log', 0), 30)
        out.append((float(eps), complex(r)))
    return out
rep = {}
r1 = [1, 2, -1, 3]; r2 = [2, -1, 3, 1]
y = [0, 0, 0, 0, 0, 1, 0, 2, 0, -1]; y2 = [0, 0, 0, 0, 1, 0, 1, 0, 3, 0]
# E = 2: quartic bubble, letters d0d1 A_7 | d2 A_5 (W = 5): eps^-1 and log
N, Rs = chain_numerator('A', unit(10, 7), 'A', unit(10, 8), [('AAAA', y, y2)], [[r1[m] + r2[m] for m in range(4)]],
                        fa=deriv_factor((1, 1, 0, 0)), fb=deriv_factor((0, 0, 1, 0)))
Rs = [[sp.Rational(x, 3) for x in R] for R in Rs]  # rescale momenta (|R|^2 ~ 1)
sing = singular_parts(N.subs({}), Rs)
res = face_residuals(N, Rs, sing, ['1e-2', '1e-3', '1e-4'])
print('E=2 quartic bubble W=5: singular', sing); print('  residuals', res)
rep['E2_quartic_W5'] = dict(singular={str(k): str(v) for k, v in sing.items()}, residuals=[(e, str(r)) for e, r in res])
# E = 3: gluon exchange, letters d0 A_7 | A_5 (W = 3), log and eps^-1? (W - w_out - 2 = -1 -> log only)
N, Rs = chain_numerator('A', unit(10, 7), 'A', unit(10, 8), [('AAA', y), ('AAA', y2)], [r1, r2],
                        fa=deriv_factor((2, 0, 0, 0)), fb=deriv_factor((0, 0, 0, 0)))
Rs = [[sp.Rational(x, 3) for x in R] for R in Rs]
sing = singular_parts(N, Rs)
res = face_residuals(N, Rs, sing, ['1e-2', '1e-3', '1e-4'])
print('E=3 gluon exchange W=4: singular', sing); print('  residuals', res)
rep['E3_gluon_W4'] = dict(singular={str(k): str(v) for k, v in sing.items()}, residuals=[(e, str(r)) for e, r in res])
json.dump(rep, open('out_numeric.json', 'w'), indent=1)
