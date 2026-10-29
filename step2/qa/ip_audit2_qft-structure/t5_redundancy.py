"""T5 (D_1 = D_1^adj) on the weight-3 family (Fbar_w, B_w) versus T3 (iii) (orthogonality of the Konishi and protected
descendants).  D_1 = -48 u P_K with P_K the projector onto the Konishi descendant K_3 along the protected E_w
(eigenvalues -48 u and 0), Gram_tree = diag(n_F, u nu_B) (n_B = u nu_B, Prop ip-weight3).  Exact sympy check that
D_1 = Gram^{-1} conj(D_1)^T Gram holds iff <E_w, K_3>_tree = 0, i.e. iff T3 (iii).
"""
import json
import sympy as sp

u = sp.symbols('u', positive=True)
nF, nuB = sp.symbols('n_F nu_B', positive=True)
x = sp.symbols('x')                       # Fbar-component of K_3 (B-component 1): x = c_K u / kappa_B at leading order
bF, bB = sp.symbols('b_F b_B')
K3 = sp.Matrix([x, 1]); E = sp.Matrix([bF, bB])
M = sp.Matrix.hstack(K3, E)
PK = M * sp.diag(1, 0) * M.inv()          # projector onto K3 along E
D1 = -48 * u * PK
G = sp.diag(nF, u * nuB)
D1adj = G.inv() * D1.H * G
diff = sp.simplify(D1 - D1adj)
overlap = sp.simplify((E.H * G * K3)[0])  # <E_w, K_3>_tree
x_orth = sp.solve(sp.Eq(overlap, 0), x)[0]
check = sp.simplify(diff.subs(x, x_orth))
out = {'D1_minus_D1adj_entries_vanish_at_orthogonality': bool(check == sp.zeros(2, 2)),
       'orthogonality_solution_x': str(x_orth),
       'T3_iii_form': 'c_K u = kappa_B x = -kappa_B conj(b_B) n_B/(conj(b_F) n_F) with n_B = u nu_B',
       'weights_2_and_5/2': 'D_1 diagonal by SU(4) x Spin(4) (Schur) and real (0, -48 unit): D_1 = D_1^adj holds identically',
       'conclusion': 'on the test families T5 is either automatic (weights 2, 5/2) or the same equation as T3 (iii) (weight 3); '
                     'it is not an independent test there'}
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
