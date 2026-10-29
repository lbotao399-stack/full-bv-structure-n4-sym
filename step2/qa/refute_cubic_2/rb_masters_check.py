"""numerical check of my exact masters J3(a,b,c;m), J2(a,b;m) (rb_loop) by scipy adaptive quadrature
(independent of the closed form: tau_3 integral done by quad as well, nested; double precision)."""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy import integrate
from rb_loop import J3, J2
L32 = np.log(1.5); L2 = np.log(2.0)
out = []
cases3 = [(0, 0, 0, 3), (1, 0, 0, 4), (0, 1, 0, 4), (0, 0, 1, 4), (2, 0, 0, 5), (1, 1, 0, 5), (0, 1, 1, 5),
          (0, 2, 0, 5), (1, 0, 1, 5), (0, 0, 0, 4), (0, 1, 0, 5), (1, 1, 1, 6)]
for (a, b, c, m) in cases3:
    r, l = J3(a, b, c, m)
    exact = float(r) + float(l) * L32
    f = lambda t3, t2, t1: t1 ** a * t2 ** b * t3 ** c / (t1 + t2 + t3) ** m
    num, err = integrate.tplquad(f, 0, 1, 1, np.inf, 1, np.inf, epsabs=1e-13, epsrel=1e-11)
    out.append({'J3': [a, b, c, m], 'exact': [str(r), str(l)], 'exact_float': exact, 'quad': num,
                'rel_diff': abs(num - exact) / abs(exact)})
    print(out[-1], flush=True)
for (a, b, m) in [(0, 0, 2), (1, 0, 3), (0, 1, 3), (1, 1, 4), (2, 0, 4), (0, 0, 3), (2, 1, 5)]:
    r, l = J2(a, b, m)
    exact = float(r) + float(l) * L2
    f = lambda t2, t1: t1 ** a * t2 ** b / (t1 + t2) ** m
    num, err = integrate.dblquad(f, 0, 1, 1, np.inf, epsabs=1e-13, epsrel=1e-12)
    out.append({'J2': [a, b, m], 'exact': [str(r), str(l)], 'exact_float': exact, 'quad': num,
                'rel_diff': abs(num - exact) / abs(exact)})
    print(out[-1], flush=True)
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'rb_masters_check.json'), 'w'), indent=1)
