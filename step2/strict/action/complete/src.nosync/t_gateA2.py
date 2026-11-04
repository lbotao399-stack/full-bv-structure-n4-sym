import numpy as np, time
import ca_common as A, ca_tree as T
rng = np.random.default_rng(5)
def rpol(t):
    Bs = A.SA.E.BASIS[t]
    while True:
        c = rng.integers(-2, 3, size=len(Bs))
        if c.any():
            return sum(ci * np.asarray(v) for ci, v in zip(c, Bs)).astype(complex)
p = np.array([1, -1, 2, 0]); mm = A.MD.MM(p)
x1, x2 = rpol('lamP'), rpol('lamM')
t0 = time.time()
D = T.mm_deriv(mm, mm.Sfree, 0)
F = T.mm_S2xi(mm, x2)
print('deriv(Sfree) == Q:', len((D(F) - mm.Q(F)).t) == 0, time.time() - t0, flush=True)
G = T.mm_S2xi(mm, x1)
t0 = time.time(); X1 = T.mm_deriv(mm, G, 0)(F); print('deriv time', time.time() - t0, flush=True)
X2 = mm.br(G, F, 0)
print('deriv(S2xi) == br:', len((X1 - X2).t) == 0, flush=True)
# odd G: translation S(P_b) (odd); {S(P_b), T1} vs br
Pb = T.mm_translation(mm, np.array([1, 2, 0, -1]))
F1 = T.mm_T(mm, x2, 'T1')
Y1 = T.mm_deriv(mm, Pb, 1)(F1); Y2 = mm.br(Pb, F1, 1)
print('deriv(odd) == br:', len((Y1 - Y2).t) == 0, len(Y2.t))
