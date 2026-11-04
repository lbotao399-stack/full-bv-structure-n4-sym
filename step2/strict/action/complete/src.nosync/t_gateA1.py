import numpy as np, time
import ca_common as A, ca_tree as T, cb_common as B
Q = B.Q
rng = np.random.default_rng(3)
def rpol(t):
    Bs = A.SA.E.BASIS[t]
    while True:
        c = rng.integers(-2, 3, size=len(Bs))
        if c.any():
            return sum(ci * np.asarray(v) for ci, v in zip(c, Bs)).astype(complex)
p = np.array([1, -1, 2, 0])
t0 = time.time(); mm = A.MD.MM(p); print('MM', time.time() - t0, flush=True)
xi1 = rpol('lamP')
t0 = time.time(); F = T.mm_T(mm, xi1, 'T1'); print('T1', time.time() - t0, len(F.t), flush=True)
t0 = time.time(); QF = mm.Q(F); print('QT1', time.time() - t0, len(QF.t), flush=True)
t0 = time.time(); S2 = T.mm_S2xi(mm, xi1); print('S2xi', time.time() - t0, len(S2.t), flush=True)
t0 = time.time(); X = mm.br(S2, S2, 0); print('br', time.time() - t0, len(X.t), flush=True)
