"""checks: sign of {S_M,-} (derivation induced by M), Q X = 0 for the class-(i) log coefficient X."""
import json, random, sympy as sp
import delta_local as DL
import fa
V = DL.V
P = list(sp.symbols('P0:4')); KAP = list(sp.symbols('K0:4'))
random.seed(5)
def rp(sp_, comps):
    u = [sp.Integer(0)] * V.DIM[sp_]
    for c in comps: u[c] = sp.Integer(random.choice([-3, -2, -1, 1, 2, 3]))
    return u
out = {}
SM = DL.S_M_functional()
# {S_M, int(A,A)} : expected -4 int(A,A) if {S_M,-} = + (derivation induced by M)
br = fa.bracket(SM, fa.MASS)
u1 = rp('A', range(10)); u2 = rp('A', range(10))
a, b = DL.to_canon(br, 'A', 'A', u1, u2, P, KAP)
ref = 2 * sum(x * y for x, y in zip(u1, u2))
out['{S_M, int(A,A)} / int(A,A)'] = [str(sp.simplify(a / ref)), [str(x) for x in b]]
# {S_M, int(dA,dA)} (weight 4): expected 0
FA = fa.local([(1, 'kinA')], 'kinA', 0)       # -1/4 int (dA,dA)
br = fa.bracket(SM, FA)
a, b = DL.to_canon(br, 'A', 'A', u1, u2, P, KAP)
out['{S_M, -1/4 int(dA,dA)}'] = [str(a), [str(x) for x in b]]
# {S_M, I_3} on (A,A,A): expected 0 (arity-one equation); evaluate with an injection leg is not needed: check two-leg-free 3-leg
# Q X = {S_2, X}
X = DL.X_functional()
QX = fa.Q(X)
res = {}
for (s1, s2, c1, c2) in [('A', 'A', DL.AVc, DL.AVc), ('A', 'A', DL.APc, DL.APc), ('c', 'As', [0], DL.AVc), ('lam', 'lam', range(16), range(16))]:
    for k in range(2):
        u1 = rp(s1, c1); u2 = rp(s2, c2)
        a, b = DL.to_canon(QX, s1, s2, u1, u2, P, KAP)
        res.setdefault(f'{s1}{list(c1)[:1]},{s2}', []).append([str(a), [str(x) for x in b]])
out['QX'] = res
# Q S_M = 0 (consistency of the representation)
QS = fa.Q(SM)
res = {}
for (s1, s2, c1, c2) in [('A', 'A', DL.AVc, DL.AVc), ('A', 'A', DL.APc, DL.APc), ('c', 'As', [0], DL.AVc), ('lam', 'lam', range(16), range(16))]:
    u1 = rp(s1, c1); u2 = rp(s2, c2)
    a, b = DL.to_canon(QS, s1, s2, u1, u2, P, KAP)
    res[f'{s1}{list(c1)[:1]},{s2}'] = [str(a), [str(x) for x in b]]
out['Q S_M'] = res
# Q of the separate pieces of X (to see that the cancellation is nontrivial)
for nm, F in [('X_V', DL.canon('A', 'As', lambda u, w, p: DL.dot(u, w, DL.AVc), lambda u, w, p, m: (u[m] if m < 4 else 0) * DL.dot(p, w, DL.AVc), 'XV2')),
              ('X_c', DL.canon('c', 'cs', None, lambda u, w, p, m: -p[m] * u[0] * w[0], 'Xc2'))]:
    u1 = rp('c', [0]); u2 = rp('As', DL.AVc)
    a, b = DL.to_canon(fa.Q(F), 'c', 'As', u1, u2, P, KAP)
    out[f'Q {nm} on (c, A*_V)'] = [str(a), [str(x) for x in b]]
json.dump(out, open('delta_local_tests.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
