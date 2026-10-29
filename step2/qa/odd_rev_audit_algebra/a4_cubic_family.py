"""Audit script a4: the arity-zero cubic freedom of Part C (F0 of the targets, design prop:o-arity0, E3b).

Inputs, all from step2/theory/cubic_structure.tex: prop:s2c-QH (Q on H1..H6 and on the fermion cubic term
(1/2) int (lambda, rho([A,lambda])), whose A|_V part alone is hit: Q = -G5; 'the Yukawa part of lambda lambda A
is not hit by Q', proof of cor:s2c-QI3), prop:s2c-br ({I_3,B_i}), cor:s2c-QI3 (I_3 in this basis);
engine B delta_results.json checks.Q_basis (Q of the field-counting functionals).

Additional cubic structures (degree 0, Weyl weight 4) not in span{H1..H6}:
  H7  = (1/2) int (lambda, rho([A|_V, lambda]))        (gauge coupling of the gaugino), Q H7 = -G5
  H8p = Yukawa part with two lambda in S+ (x) W,  H8m = two lambda in S- (x) W*;  Q H8p = Q H8m = 0
  (1/2) int (lambda, rho([A|_Wperp, lambda])) = H8p + H8m.
Checks: kernel of Q on the full cubic space; closed (non-equivariant) cubic-order freedom
  {(beta,h): sum beta_i {I_3,B_i} + Q h = 0};  the arity-zero vectors of s(N_f) (five field countings,
  lambda split by chirality) and of the coupling direction S; their rank and position relative to the
  Part-C family (h7 = h8p = h8m = 0).
"""
import json, os
import sympy as sp

G = ['G1', 'G2', 'G3', 'G4', 'G5']
Hn = ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8p', 'H8m']
Bn = ['B1=int(dA_V,dA_V)', 'B2=int(dA_W,dA_W)', 'B3=int(lam,dslash lam)', 'B4=int(dc,A*)']
R = sp.Rational
QH = {'H1': {'G1': -2}, 'H2': {'G3': -1}, 'H3': {'G1': R(-1, 2), 'G4': -1}, 'H4': {'G3': R(-1, 2)},
      'H5': {'G5': 1}, 'H6': {'G4': 2}, 'H7': {'G5': -1}, 'H8p': {}, 'H8m': {}}
IB = {0: {'G1': 2}, 1: {'G3': 2}, 2: {'G5': 2}, 3: {'G1': R(1, 2), 'G3': R(1, 2), 'G5': -1}}
I3 = {'H1': R(-1, 4), 'H2': R(-1, 2), 'H3': 1, 'H4': 1, 'H5': 1, 'H6': R(1, 2), 'H7': 1, 'H8p': 1, 'H8m': 1}

# matrix of the closedness map (beta, h) -> G
M = sp.zeros(5, 4 + 9)
for j in range(4):
    for g, v in IB[j].items():
        M[G.index(g), j] += v
for j, h in enumerate(Hn):
    for g, v in QH[h].items():
        M[G.index(g), 4 + j] += v
res = {}
QHmat = M[:, 4:]
res['Q_rank_on_H1..H6'] = int(QHmat[:, :6].rank())
res['Q_rank_on_H1..H8'] = int(QHmat.rank())
ker = QHmat.nullspace()
res['ker_Q_on_full_cubic_space_dim'] = len(ker)
res['ker_Q_basis'] = [dict((Hn[i], str(v[i])) for i in range(9) if v[i] != 0) for v in ker]
res['QI3_is_zero'] = all(sum(I3[h] * QH[h].get(g, 0) for h in Hn) == 0 for g in G)

# field-counting weights (field minus antifield of the type) per cubic monomial
W = {'V': {'H1': 3, 'H2': 1, 'H7': 1},
     'W': {'H2': 2, 'H8p': 1, 'H8m': 1},
     'lp': {'H7': 1, 'H8p': 2},
     'lm': {'H7': 1, 'H8m': 2},
     'c': {'H3': 1, 'H4': 1, 'H5': 1, 'H6': 1}}
Qquad = {'V': [R(-1, 2), 0, 0, -1], 'W': [0, R(-1, 2), 0, 0], 'lp': [0, 0, R(-1, 2), 0], 'lm': [0, 0, R(-1, 2), 0],
         'c': [0, 0, 0, -1]}
vecs = {}
tvals = {}
for f in W:
    t = sp.Symbol('t')
    h = [t * W[f].get(x, 0) * I3[x] for x in Hn]
    v = sp.Matrix(Qquad[f] + h)
    eqs = list(M * v)
    sol = sp.solve([e for e in eqs if e != 0], t, dict=True)
    tvals[f] = [str(s) for s in sol]
    if len(sol) == 1:
        vecs['s(N_' + f + ')'] = v.subs(sol[0])
    else:
        vecs['s(N_' + f + ')'] = None
vecs['S (coupling)'] = sp.Matrix([R(-1, 4), R(-1, 4), R(1, 2), 1] + [I3[x] for x in Hn])
res['closedness_fixes_relative_sign_t'] = tvals
res['each_vector_closed'] = {k: (None if v is None else all(e == 0 for e in (M * v))) for k, v in vecs.items()}
res['vectors_(beta1..4|H1..H6,H7,H8p,H8m)'] = {k: (None if v is None else [str(e) for e in v]) for k, v in vecs.items()}
V = sp.Matrix.hstack(*[v for v in vecs.values() if v is not None])
res['rank_of_six_arity_zero_vectors'] = int(V.rank())
free = M.nullspace()
res['dim_closed_cubic_order_freedom_full'] = len(free)
# parity-even subspace h8p = h8m
P = sp.zeros(1, 13); P[0, 4 + 7] = 1; P[0, 4 + 8] = -1
res['dim_closed_freedom_parity_even'] = len(sp.Matrix.vstack(M, P).nullspace())
# Part-C restricted family: h7 = h8p = h8m = 0
Rst = sp.zeros(3, 13)
for i, k in enumerate((10, 11, 12)):
    Rst[i, k] = 1
FC = sp.Matrix.vstack(M, Rst).nullspace()
res['dim_PartC_family_(F3_in_span_H1..H6)'] = len(FC)
res['vectors_with_nonzero_H7_or_H8_components'] = [k for k, v in vecs.items() if v is not None and any(v[i] != 0 for i in (10, 11, 12))]
FCm = sp.Matrix.hstack(*FC)
inter = sp.Matrix.hstack(V, -FCm).nullspace()
res['dim_span(six)_intersect_PartC_family'] = len(inter)
# complement: closed directions not reachable by the six (dimension)
res['dim_closed_freedom_mod_span(six)_full'] = len(free) - int(V.rank())
res['dim_closed_freedom_mod_span(six)_parity_even'] = res['dim_closed_freedom_parity_even'] - int(
    sp.Matrix.hstack(*[v for v in vecs.values() if v is not None]).rank())
# Is the Yukawa direction H8p+H8m (closed: Q-closed, no quadratic part) in span of the six?
y = sp.zeros(13, 1); y[11] = 1; y[12] = 1
res['Yukawa_H8p+H8m_in_span(six)'] = int(sp.Matrix.hstack(V, y).rank()) == int(V.rank())
h57 = sp.zeros(13, 1); h57[4 + 4] = 1; h57[4 + 6] = 1
res['H5+H7_closed'] = all(e == 0 for e in (M * h57))
res['H5+H7_in_span(six)'] = int(sp.Matrix.hstack(V, h57).rank()) == int(V.rank())
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a4_cubic_family.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
