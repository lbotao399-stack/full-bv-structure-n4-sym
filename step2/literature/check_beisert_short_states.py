"""Check of derived item [D5] in n4_untwisted.tex: eq. (4.13)-(4.14) of N. Beisert, hep-th/0310252v2.

H(g) on (O_1, O_2), O_1 = eps^{alpha beta} Tr psi_alpha psi_beta, O_2 = eps^{abc} Tr phi_a phi_b phi_c,
convention H O_i = sum_j M_ij O_j.  Entry M_21 is printed as -9 sqrt2 e^{i b2} e^{-i b2} alpha1^3 g^3;
here M_21 = -9 sqrt2 X alpha1^3 g^3 with X left free, to test which X reproduces (4.14).
"""
import sympy as sp

g, a1, a3, eps, b2 = sp.symbols('g alpha1 alpha3 epsilon beta2', positive=True)
X = sp.symbols('X')
M = sp.Matrix([[3 + eps*g**4, -2*sp.sqrt(2)*sp.exp(sp.I*b2)*a1**3*g**3],
               [-9*sp.sqrt(2)*X*a1**3*g**3, 3 + 6*a1**2*g**2 - 18*a1**4*g**4 + 12*a1*a3*g**4]])
lam = sp.symbols('lam')
p = (M - lam*sp.eye(2)).det()
# perturbative roots: lam = 3 + c2 g^2 + c4 g^4
for lead in [0, 6*a1**2]:
    c4 = sp.symbols('c4')
    trial = 3 + lead*g**2 + c4*g**4
    ser = sp.series(p.subs(lam, trial), g, 0, 7).removeO()
    sol = sp.solve(sp.expand(ser).coeff(g, 6), c4)
    print('branch lead', lead, ': g^4 coefficient =', [sp.simplify(s) for s in sol])
print('(4.14) requires: branch 0 -> eps - 6 alpha1^4 ; branch 6alpha1^2 -> -12 alpha1^4 + 12 alpha1 alpha3')
print('solve for X:', sp.solve(sp.Eq(-2*sp.sqrt(2)*sp.exp(sp.I*b2)*(-9*sp.sqrt(2))*X, 36), X))
# first-order eigen-operators (rows of left eigenvectors in the O-basis)
x2 = sp.simplify(M[0, 1]/(3 - M[1, 1]))          # protected: O_1 + x2 O_2
x1 = sp.simplify(M[1, 0]/(M[1, 1] - M[0, 0]))    # Konishi descendant: O_2 + x1 O_1
print('protected  : O_1 + (', sp.limit(x2/g, g, 0), ') g O_2 + O(g^2)')
print('Konishi    : O_2 + (', sp.limit(x1.subs(X, sp.exp(-sp.I*b2))/g, g, 0), ') g O_1 + O(g^2)')
