"""Check of derived item [D2] in n4_untwisted.tex.

Input: coefficients of (Q_2)^a_alpha and (S_2)^alpha_a in Table 1 of
N. Beisert, hep-th/0310252v2 (common phases e^{+i beta_1}, e^{-i beta_1} stripped).
Structures, in the order of Table 1:
  Q_2: [ {ab->alpha b}+{ba->b alpha},  {a beta->alpha beta}-{beta a->beta alpha},
         {ab->b alpha}+{ba->alpha b},  {a beta->beta alpha}-{beta a->alpha beta} ]
  S_2: [ {alpha b->ab}+{b alpha->ba},  {alpha beta->a beta}-{beta alpha->beta a},
         {alpha b->ba}+{b alpha->ab},  -{alpha beta->beta a}+{beta alpha->a beta} ]
Output: the gamma choices that set Q_2=0 or S_2=0, the other generator at that choice,
and the gamma-independent sums Q_2[i]+S_2[i].
"""
import sympy as sp

a, g1, g2, g3, g4 = sp.symbols('alpha1 gamma1 gamma2 gamma3 gamma4')
I = sp.I
q = sp.Rational(1, 4)
Q2 = [q*a**2 - I*g3/2 + I*g4/2, I*g3/2 + I*g4/2, -q*a**2 - I*g1, q*a**2 + I*g2]
S2 = [q*a**2 + I*g3/2 - I*g4/2, -I*g3/2 - I*g4/2, -q*a**2 + I*g1, q*a**2 - I*g2]
gam = [g1, g2, g3, g4]

solQ = sp.solve(Q2, gam, dict=True)
solS = sp.solve(S2, gam, dict=True)
print('Q2 = 0 at', solQ, '; then S2 coefficients =', [sp.simplify(c.subs(solQ[0])) for c in S2])
print('S2 = 0 at', solS, '; then Q2 coefficients =', [sp.simplify(c.subs(solS[0])) for c in Q2])
print('Q2 = S2 = 0 solvable:', bool(sp.solve(Q2 + S2, gam, dict=True)))
print('gamma-independent sums Q2[i]+S2[i] =', [sp.simplify(x + y) for x, y in zip(Q2, S2)])
