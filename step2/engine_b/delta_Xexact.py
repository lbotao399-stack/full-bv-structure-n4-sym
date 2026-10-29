"""Is the class-(i) log coefficient X Q-exact?  Y = y1 int(iota_x A|_V, c*) + y2 int(rho(x) lambda*, lambda*)."""
import json, random, sympy as sp
import delta_local as DL
import fa
from realization import Rb
V = DL.V; I = sp.I
P = list(sp.symbols('P0:4')); KAP = list(sp.symbols('K0:4'))
rnd = random.Random(4)
def rp(sp_, comps):
    u = [sp.Integer(0)] * V.DIM[sp_]
    for c in comps: u[c] = sp.Integer(rnd.choice([-3, -2, -1, 1, 2, 3]))
    return u
Y1 = DL.canon('A', 'cs', None, lambda u, w, p, m: -I * (u[m] if m < 4 else 0) * w[0], 'Y1')
Y2 = DL.canon('lams', 'lams', None, lambda u, w, p, m: -I * DL.lam_pair(list(Rb[m] * sp.Matrix(u)), w), 'Y2')
X = DL.X_functional()
out = {}
for nm, F in [('Q Y1', fa.Q(Y1)), ('Q Y2', fa.Q(Y2)), ('X', X)]:
    res = {}
    rnd.seed(4)
    for (s1, s2, c1, c2) in [('A', 'As', DL.AVc, DL.AVc), ('c', 'cs', [0], [0]), ('lam', 'lams', range(16), range(16))]:
        u1 = rp(s1, c1); u2 = rp(s2, c2)
        a, b = DL.to_canon(F, s1, s2, u1, u2, P, KAP)
        res[s1] = [str(a), [str(x) for x in b]]
    out[nm] = res
print(json.dumps(out, indent=1))
json.dump(out, open('delta_Xexact.json', 'w'), indent=1)
