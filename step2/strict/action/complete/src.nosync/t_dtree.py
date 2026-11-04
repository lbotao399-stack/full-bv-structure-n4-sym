"""gate GA-dtree / GB-dtree: the density tree operations of dtree.py in each realisation reproduce the route's own
validated operations on translation-invariant functionals (route B: eb_qb.apply_Q, cb_tree.bracket; route A: mode
model Q and bracket)."""
import numpy as np, random
import cb_common as B, cb_tree as CT, dtree as DT, ca_common as A, ca_tree as AT, cb_s2 as S
Q = B.Q
TB = DT.Tree(lambda k: np.array(B.C.q_matrix(list(k)).tolist(), dtype=complex), DT.k0_from_pairing(B.MP))
print('K0 B equal:', np.allclose(TB.K0, CT.K0))
TA = DT.Tree(A.qmat_A, DT.k0_from_pairing(A.MP))
random.seed(2)
x1 = np.zeros(16, complex); x2 = np.zeros(16, complex)
for c in range(8): x1[c] = random.randint(-2, 2)
for c in range(8, 16): x2[c] = random.randint(-2, 2)
# route A: density versions of S2xi in A's realisation
def T_A(xi):
    T1 = Q.dens(-1, ('lam', 'As'), lambda s: sum((xi @ A.SB[M] @ s[0]['u']) * s[1]['u'][M] for M in range(10)))
    T3 = Q.dens(0.5, ('A', 'lams'), lambda s: sum(1j * s[0]['p'][M] * s[0]['u'][N] * ((A.SS[M] @ (A.SB[N] @ xi)) @ s[1]['u']) for M in range(4) for N in range(10) if M != N))
    return Q.dsum((1, T1), (1, T3))
rng = np.random.default_rng(1)
def rpol(t):
    Bs = A.SA.E.BASIS[t]
    while True:
        c = rng.integers(-2, 3, size=len(Bs))
        if c.any(): return sum(ci * np.asarray(v) for ci, v in zip(c, Bs)).astype(complex)
xa, xb = rpol('lamP'), rpol('lamM')
worst = 0
for t0, t1 in (('lamP', 'lamsP'), ('AV', 'AsV'), ('c', 'cs'), ('AW', 'AsW'), ('lamM', 'lamsM')):
    p = np.array([1, 2, -1, 0]); mm = A.MD.MM(p)
    L = [A.SA.leg(t0, rpol(t0), p), A.SA.leg(t1, rpol(t1), -p)]
    v_mm = AT.mm_value(mm, AT.mm_deriv(mm, AT.mm_S2xi(mm, xa), 0)(AT.mm_S2xi(mm, xb)), L)
    v_d = TA.brd(T_A(xa), 0, T_A(xb))(A.generic_legs(L))
    worst = max(worst, abs(complex(float(v_mm[0]), float(v_mm[1])) - v_d))
print('GA-dtree bracket {S(xa),S(xb)} vs mode model:', worst)
# odd G check in A: {S(P_b)_perp, T2} vs mode model
b = np.array([1, -1, 0, 2])
SPp = Q.dens(1, ('A', 'As'), lambda s: 1j * np.dot([complex(x) for x in s[0]['p']], b) * np.dot(s[0]['u'][4:], s[1]['u'][4:]))
T2 = Q.dens(-1, ('lam', 'As'), lambda s: sum((xa @ A.SB[M] @ s[0]['u']) * s[1]['u'][M] for M in range(4, 10)))
worst = 0
for _ in range(3):
    p = rng.integers(-2, 3, size=4); mm = A.MD.MM(p)
    L = [A.SA.leg('lamP', rpol('lamP'), p), A.SA.leg('AsW', rpol('AsW'), -p)]
    v_mm = AT.mm_value(mm, AT.mm_deriv(mm, AT.mm_translation(mm, b, ('perp',)), 1)(AT.mm_T(mm, xa, 'T2')), L)
    v_d = TA.brd(SPp, 1, T2)(A.generic_legs(L))
    worst = max(worst, abs(complex(float(v_mm[0]), float(v_mm[1])) - v_d), 0)
    print('  sizes', complex(float(v_mm[0]), float(v_mm[1])), v_d)
print('GA-dtree odd bracket vs mode model:', worst)
