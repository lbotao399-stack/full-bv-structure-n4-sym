"""Exact two-mode check of the weights of {Psi,-}, Psi = 2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c)
(units u), on every letter (fields and antifields), and of Q Psi = quadratic log counterterm; bracket [C07 Lemma 4.3.2]
with classical_spec S3 normalisation; Q = {S_free,-}."""
import json
import sympy as sp
from grass import GA, PAR
import free as F
y = F.y
Sf = F.S_free()
def mode_pair(n1, i1, n2, i2):
    return y(n1, i1, 1) * y(n2, i2, -1) + y(n1, i1, -1) * y(n2, i2, 1)
Psi = sum((mode_pair('A', M, 'As', M) for M in range(10)), GA()) * 2 \
    + sum((mode_pair('ls', a, 'l', a) for a in range(16)), GA()) * (-4) + mode_pair('cs', 0, 'c', 0) * (-1)
out = {}
for nm, (O, par) in {'A|_V': (y('A', 1, 1), 0), 'A|_Wperp': (y('A', 6, 1), 0), 'lambda': (y('l', 3, 1), 1), 'c': (y('c', 0, 1), 1),
                     'A^*|_V': (y('As', 2, 1), 1), 'A^*|_Wperp': (y('As', 7, 1), 1), 'lambda^*': (y('ls', 5, 1), 0),
                     'c^*': (y('cs', 0, 1), 0)}.items():
    br = F.bracket(Psi, O, 1).map(sp.expand)          # {Psi, O}
    (k, c), = O.t.items()
    assert set(br.t) <= {k}
    out[nm] = str(sp.nsimplify(br.t.get(k, 0) / c))
print('{Psi, letter}/letter :', out)
# derivation property check on a cubic monomial: {Psi, x y z} = (w_x + w_y + w_z) x y z
m = y('c', 0, 1) * y('A', 1, -1) * y('As', 1, 1)
br = F.bracket(Psi, m, 1).map(sp.expand)
(k, c), = m.t.items()
print('{Psi, c A A^*}/(c A A^*) =', sp.nsimplify(br.t.get(k, 0) / c), ' other terms:', len(set(br.t) - {k}))
json.dump({'Psi_weights': out}, open('logs/r2_T3_weights.json', 'w'), indent=1)
