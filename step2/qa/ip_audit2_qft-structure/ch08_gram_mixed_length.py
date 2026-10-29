"""CH08 cor:c8-gram on a family with two lengths (the T3 family O_zeta, B_w): the expansion 'order hbar^0: Y_0 = X_0^adj,
order hbar: Y_1 = X_1^adj + [Y_0, Gram_tree^{-1} Gram_1]' treats Gram_tree as order hbar^0, but the tree norm of the
length-3 operator is of relative order one unit (n_B = u nu_B, Prop ip-weight3).  Exact sympy check of what the exact
relation conj(X)^T Gram = Gram Y gives, entry (O_zeta, B_w), with X = matrix of -Lambda^(1)(xi), Y = that of
-Lambda^(1)(xi^dagger) on the basis (O_zeta, B_w); u = the unit (expansion parameter).
"""
import json
import sympy as sp

u = sp.symbols('u', positive=True)
nz, nuB = sp.symbols('n_zeta nu_B', positive=True)
kB, cB, U1 = sp.symbols('kappa_B c_Bp U_1')
# X: classical kappa_B (zeta -> B, length increasing) + one-loop c_Bp u ; Y: classical part 0 in (zeta <- B), one-loop u U_1
X = sp.Matrix([[0, 0], [-(kB + u * cB), 0]])
Y = sp.Matrix([[0, u * U1], [0, 0]])
G = sp.diag(nz, u * nuB)                      # Gram_tree on (O_zeta, B_w); Gram_1 omitted (zero in this entry)
rel = sp.expand((X.H * G - G * Y)[0, 1])
exact_U1 = sp.solve(sp.Eq(rel, 0), U1)[0]
# cor:c8-gram read literally: Y_0 = X_0^adj with Gram_tree = G, and Y_1 = X_1^adj + [Y_0, G^{-1} Gram_1] (Gram_1 = 0 here)
X0 = X.subs(u, 0); X1 = sp.diff(X, u)
adj = lambda M: G.inv() * M.H * G
Y0_pred = sp.simplify(adj(X0))
Y1_pred = sp.simplify(adj(X1))
out = {'exact_relation_entry_(zeta,B)': str(sp.simplify(rel)),
       'exact_one_loop_entry_U_1': str(sp.simplify(exact_U1)),
       'cor_c8_gram_Y0_pred_entry_(zeta,B)': str(Y0_pred[0, 1]),
       'cor_c8_gram_Y1_pred_entry_(zeta,B)': str(sp.simplify(Y1_pred[0, 1] * u)),
       'classical_Y0_entry_(zeta,B)': '0 (classical components never decrease the length, Prop ip-orders (ii))',
       'conclusion': ('read literally, cor:c8-gram puts the one-loop quantity -conj(kappa_B) u nu_B/n_zeta into Y_0 '
                      '(classical) and predicts no such term in Y_1; the exact relation puts it in the one-loop entry '
                      '(Prop ip-weight3 (i), Thm ip-graded order one half). Gram_tree has entries of different orders '
                      'in the unit on mixed-length families, so the hbar-expansion of cor:c8-gram must be replaced by the '
                      'graded expansion.')}
ok = sp.simplify(exact_U1.subs(u, 0) + sp.conjugate(kB) * nuB / nz) == 0 and sp.simplify(Y0_pred[0, 1] + sp.conjugate(kB) * u * nuB / nz) == 0
out['checks_pass'] = bool(ok)
json.dump(out, open(__file__.replace('.py', '.json'), 'w'), indent=1)
print(json.dumps(out, indent=1))
