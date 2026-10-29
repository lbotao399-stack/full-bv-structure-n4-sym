"""Round-2 check of the arity-zero cubic-order freedom (design prop:o-arity0, targets F0, engine E3b).

Inputs (all transcribed from step2/theory/cubic_structure.tex, as in the audit script a4):
  Q on the cubic structures (prop:s2c-QH; Q of the full lambda-lambda-A term = -G5, its Yukawa part not hit,
  proof of cor:s2c-QI3), {I_3, B_i} on the beta-directions (prop:s2c-br), I_3 in the basis (cor:s2c-QI3),
  Q of the field-counting functionals on the quadratic beta-directions (engine B checks.Q_basis via a4).
New here: parity bookkeeping.  The parity R of cubic_structure.tex (reflection of one direction of V and one
of W_perp) exchanges the chiralities: H8p <-> H8m, N_lp <-> N_lm, fixes every other basis element.
Output: dimensions of the closed space C, of the span V of the six arity-zero vectors, of C/V, split by
parity, and a basis of C/V represented inside the Part-C family (h7 = h8p = h8m = 0) when possible.
"""
import json, os
import sympy as sp

R = sp.Rational
G = ['G1', 'G2', 'G3', 'G4', 'G5']
Hn = ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8p', 'H8m']
names = ['b1', 'b2', 'b3', 'b4'] + Hn
QH = {'H1': {'G1': -2}, 'H2': {'G3': -1}, 'H3': {'G1': R(-1, 2), 'G4': -1}, 'H4': {'G3': R(-1, 2)},
      'H5': {'G5': 1}, 'H6': {'G4': 2}, 'H7': {'G5': -1}, 'H8p': {}, 'H8m': {}}
IB = {0: {'G1': 2}, 1: {'G3': 2}, 2: {'G5': 2}, 3: {'G1': R(1, 2), 'G3': R(1, 2), 'G5': -1}}
I3 = {'H1': R(-1, 4), 'H2': R(-1, 2), 'H3': 1, 'H4': 1, 'H5': 1, 'H6': R(1, 2), 'H7': 1, 'H8p': 1, 'H8m': 1}
M = sp.zeros(5, 13)
for j in range(4):
    for g, v in IB[j].items():
        M[G.index(g), j] += v
for j, h in enumerate(Hn):
    for g, v in QH[h].items():
        M[G.index(g), 4 + j] += v
# field-counting weights of the cubic monomials (field minus antifield), quadratic Q-images (as in a4)
W = {'V': {'H1': 3, 'H2': 1, 'H7': 1}, 'W': {'H2': 2, 'H8p': 1, 'H8m': 1},
     'lp': {'H7': 1, 'H8p': 2}, 'lm': {'H7': 1, 'H8m': 2}, 'c': {'H3': 1, 'H4': 1, 'H5': 1, 'H6': 1}}
Qquad = {'V': [R(-1, 2), 0, 0, -1], 'W': [0, R(-1, 2), 0, 0], 'lp': [0, 0, R(-1, 2), 0],
         'lm': [0, 0, R(-1, 2), 0], 'c': [0, 0, 0, -1]}
vecs = {}
for f in W:
    t = sp.Symbol('t')
    v = sp.Matrix(Qquad[f] + [t * W[f].get(x, 0) * I3[x] for x in Hn])
    sol = sp.solve([e for e in (M * v) if e != 0], t, dict=True)
    assert len(sol) == 1
    vecs['s(N_' + f + ')'] = v.subs(sol[0])
vecs['S'] = sp.Matrix([R(-1, 4), R(-1, 4), R(1, 2), 1] + [I3[x] for x in Hn])
for k, v in vecs.items():
    assert all(e == 0 for e in M * v), k
# parity
Par = sp.eye(13)
i8p, i8m = 4 + Hn.index('H8p'), 4 + Hn.index('H8m')
Par[i8p, i8p] = 0; Par[i8m, i8m] = 0; Par[i8p, i8m] = 1; Par[i8m, i8p] = 1
C = M.nullspace()
Cm = sp.Matrix.hstack(*C)
Vm = sp.Matrix.hstack(*vecs.values())
def dim(cols):
    return int(sp.Matrix.hstack(*cols).rank()) if cols else 0
even = lambda m: (m + Par * m) / 2
odd = lambda m: (m - Par * m) / 2
Ce = [even(Cm[:, j]) for j in range(Cm.shape[1])]
Co = [odd(Cm[:, j]) for j in range(Cm.shape[1])]
Ve = [even(Vm[:, j]) for j in range(Vm.shape[1])]
Vo = [odd(Vm[:, j]) for j in range(Vm.shape[1])]
# the projections lie in C and V since Par preserves C (check) and permutes the six
assert all(e == 0 for e in M * Par * Cm)
res = {'dim_C': len(C), 'dim_C_even': dim(Ce), 'dim_C_odd': dim(Co), 'rank_six': int(Vm.rank()),
       'dim_V_even': dim(Ve), 'dim_V_odd': dim(Vo)}
res['dim_C_mod_V'] = res['dim_C'] - res['rank_six']
res['dim_C_mod_V_even'] = res['dim_C_even'] - res['dim_V_even']
res['dim_C_mod_V_odd'] = res['dim_C_odd'] - res['dim_V_odd']
# Part-C family PC: h7 = h8p = h8m = 0
Rst = sp.zeros(3, 13)
for i, k in enumerate((10, 11, 12)):
    Rst[i, k] = 1
PC = sp.Matrix.vstack(M, Rst).nullspace()
PCm = sp.Matrix.hstack(*PC)
res['dim_PartC'] = len(PC)
res['dim_PartC_cap_V'] = len(sp.Matrix.hstack(Vm, -PCm).nullspace())
res['dim_(PartC+V)'] = int(sp.Matrix.hstack(Vm, PCm).rank())
res['PartC_plus_V_equals_C'] = res['dim_(PartC+V)'] == res['dim_C']
# representatives of C/V: test the Yukawa and H5+H7 directions and the Part-C basis
def vec(d):
    v = sp.zeros(13, 1)
    for k, c in d.items():
        v[names.index(k)] = c
    return v
tests = {'Yukawa H8p+H8m': vec({'H8p': 1, 'H8m': 1}), 'H8p-H8m': vec({'H8p': 1, 'H8m': -1}),
         'H5+H7': vec({'H5': 1, 'H7': 1}), 'kerQ_PartC_1: H1-4H3-2H6': vec({'H1': 1, 'H3': -4, 'H6': -2}),
         'kerQ_PartC_2: H2-2H4': vec({'H2': 1, 'H4': -2})}
res['tests'] = {}
for k, v in tests.items():
    closed = all(e == 0 for e in M * v)
    inV = int(sp.Matrix.hstack(Vm, v).rank()) == int(Vm.rank())
    res['tests'][k] = {'closed': closed, 'in_span_of_the_six': inV}
# a basis of C/V made of Part-C vectors, and the expression of Yukawa and H5+H7 modulo V in it
basis = []
cur = Vm
for j in range(PCm.shape[1]):
    cand = sp.Matrix.hstack(cur, PCm[:, j])
    if cand.rank() > cur.rank():
        basis.append(PCm[:, j]); cur = cand
res['quotient_basis_from_PartC'] = [{names[i]: str(b[i]) for i in range(13) if b[i] != 0} for b in basis]
Bm = sp.Matrix.hstack(*basis)
for k in ('Yukawa H8p+H8m', 'H5+H7'):
    v = tests[k]
    a = sp.symbols('a0:%d' % Bm.shape[1]); cvar = sp.symbols('c0:%d' % Vm.shape[1])
    expr = v - Bm * sp.Matrix(a) - Vm * sp.Matrix(cvar)
    sol = sp.solve(list(expr), list(a) + list(cvar), dict=True)
    res['tests'][k]['mod_V_in_PartC_basis'] = None if not sol else [str(sol[0].get(x, x)) for x in a]
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cubic_quotient.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
