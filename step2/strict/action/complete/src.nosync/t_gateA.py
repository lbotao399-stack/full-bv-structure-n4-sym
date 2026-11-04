"""gates GA-Q (Q on densities with qmat_A = mode-model Q), GA-dens (dens values = mode-model values), GA-S2
(Q S^(2) + S(P_Gamma) - {S(xi1),S(xi2)} = 0 in the mode model)."""
import numpy as np, random
import ca_common as A, ca_tree as T, cb_common as B, bilin as BL
from fractions import Fraction as Fr
Q = B.Q
rng = np.random.default_rng(3)
def rpol(t):
    Bs = A.SA.E.BASIS[t]
    while True:
        c = rng.integers(-2, 3, size=len(Bs))
        if c.any():
            return sum(ci * np.asarray(v) for ci, v in zip(c, Bs)).astype(complex)
def rp():
    while True:
        p = rng.integers(-2, 3, size=4)
        if p @ p:
            return p
xi1, xi2 = rpol('lamP'), rpol('lamM')
# GA-dens: T1 as density in A's realisation vs mode model
T1d = Q.dens(-1, ('lam', 'As'), lambda s: sum((xi1 @ A.SB[M] @ s[0]['u']) * s[1]['u'][M] for M in range(4)))
worst = 0
for _ in range(4):
    p = rp(); L = [A.SA.leg('lamM', rpol('lamM'), p), A.SA.leg('AsV', rpol('AsV'), -p)]
    mm = A.MD.MM(p)
    v1 = T.mm_value(mm, T.mm_T(mm, xi1, 'T1'), L); v2 = T1d(A.generic_legs(L))
    worst = max(worst, abs(complex(float(v1[0]), float(v1[1])) - v2))
print('GA-dens (T1):', worst)
# GA-Q: Q T1 via density vs mode model Q T1 ; also Q of T3
worst = 0
for which, spc in (('T1', ('lamM', 'AV')), ('T3', ('lamsP', 'AV')), ('T3', ('lamsP', 'c'))):
    for _ in range(3):
        p = rp(); L = [A.SA.leg(spc[0], rpol(spc[0]), p), A.SA.leg(spc[1], rpol(spc[1]), -p)]
        mm = A.MD.MM(p)
        F = T.mm_T(mm, xi1, which)
        v1 = T.mm_value(mm, mm.Q(F), L)
        if which == 'T1':
            d = T1d
        else:
            d = Q.dens(0.5, ('A', 'lams'), lambda s: sum(1j * s[0]['p'][M] * s[0]['u'][N] * ((A.SS[M] @ (A.SB[N] @ xi1)) @ s[1]['u']) for M in range(4) for N in range(4) if M != N))
        v2 = T.apply_Q(d)(A.generic_legs(L))
        worst = max(worst, abs(complex(float(v1[0]), float(v1[1])) - v2))
print('GA-Q:', worst)
# GA-S2 in the mode model
for nm, c1, c2 in (('(+,-)', 'lamP', 'lamM'), ('(+,+)', 'lamP', 'lamP'), ('(-,-)', 'lamM', 'lamM')):
    x1, x2 = rpol(c1), rpol(c2)
    p = rp(); mm = A.MD.MM(p)
    S2 = T.mm_S2(mm, x1, x2)
    lhs = mm.Q(S2) + T.mm_translation(mm, T.gam(x1, x2)[:4]) - mm.br(T.mm_S2xi(mm, x1), T.mm_S2xi(mm, x2), 0)
    alt = mm.Q(S2) + T.mm_translation(mm, T.gam(x1, x2)[:4]) + mm.br(T.mm_S2xi(mm, x1), T.mm_S2xi(mm, x2), 0)
    print('GA-S2', nm, 'residual terms', len(lhs.t), ' (opposite sign:', len(alt.t), ')  size Q S2', len(mm.Q(S2).t))
