"""Same residual test as t_numeric.py for the E = 3 gluon-exchange chain, double precision (scipy dblquad in logarithmic
Schwinger variables), eps = 1e-2 .. 1e-5: residual r(eps) = eps w' + sum m a_m eps^-m - b must decrease like eps."""
import json
import numpy as np, sympy as sp
from scipy import integrate
from ref_core import *
r1 = [1, 2, -1, 3]; r2 = [2, -1, 3, 1]
y = [0, 0, 0, 0, 0, 1, 0, 2, 0, -1]; y2 = [0, 0, 0, 0, 1, 0, 1, 0, 3, 0]
N, Rs = chain_numerator('A', unit(10, 7), 'A', unit(10, 8), [('AAA', y), ('AAA', y2)], [r1, r2], fa=deriv_factor((2, 0, 0, 0)))
Rs = [[sp.Rational(x, 3) for x in R] for R in Rs]
sing = singular_parts(N, Rs)
E = len(Rs); s = sp.symbols('t0:%d' % E); tau = sum(s)
c = [-sum(s[j] * Rs[j][mu] for j in range(E)) / tau for mu in range(4)]
tPhi = sum(s[j] * dot4(Rs[j], Rs[j]) for j in range(E)) - sum((sum(s[j] * Rs[j][mu] for j in range(E))) ** 2 for mu in range(4)) / tau
integ = 0; lapN = sp.expand(N); j = 0
while lapN != 0:
    integ += lapN.subs({L[m]: c[m] for m in range(4)}, simultaneous=True) / (sp.factorial(j) * 4 ** j * tau ** (2 + j))
    lapN = sp.expand(laplacian(lapN)); j += 1
f = sp.lambdify(s, sp.re(integ) * sp.exp(-tPhi), 'numpy')
fi = sp.lambdify(s, sp.im(integ) * sp.exp(-tPhi), 'numpy')
out = []
for eps in [1e-2, 1e-3, 1e-4, 1e-5, 1e-6]:
    tot = 0
    for e in range(E):
        oth = [k for k in range(E) if k != e]
        def g(xv, yv, ff):
            a = [None] * E; a[e] = eps; a[oth[0]] = eps * np.exp(xv); a[oth[1]] = eps * np.exp(yv)
            return ff(*a) * eps ** 2 * np.exp(xv + yv)
        Lg = -np.log(eps)
        val = 0
        for ff in (f, fi):
            v, err = integrate.dblquad(lambda yv, xv: g(xv, yv, ff), 0, Lg + 9, 0, Lg + 9, epsabs=1e-12, epsrel=1e-12)
            val = val + (v if ff is f else 1j * v)
        tot += val
    r = -eps * tot + sum(m * complex(sp.N(v)) * eps ** (-m) for m, v in sing.items() if m != 'log') - complex(sp.N(sing.get('log', 0)))
    out.append((eps, [r.real, r.imag]))
    print(eps, r)
print('singular parts:', sing)
json.dump(dict(singular={str(k): str(v) for k, v in sing.items()}, residuals=out), open('out_numeric_fast.json', 'w'), indent=1)
