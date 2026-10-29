"""D4 engine A: independent check of the structure of Y1a.  Ward identity at the scalar-vector vertex X (the composite
line is longitudinal): beta(Y1a) = lambda_1 Bhat(p2) + lambda_2 p1^2 That(p1,p2) with (units 1/(16 pi^2))
  Bhat_nu(p)      = 16 pi^2 int d^4l/(2pi)^4 l_nu / (l^4 (l+p)^2)                 = -p_nu/p^2,
  That_nu(p1,p2)  = 16 pi^2 int d^4l/(2pi)^4 l_nu / (l^4 (l+p2)^2 (l-p1)^2)
                  = - int_{x+y+z=1} x (y p2 - z p1)_nu / Delta^2,  Delta = y p2^2 + z p1^2 - |y p2 - z p1|^2.
Fit lambda_1, lambda_2 at each kinematic point (4 components, 2 unknowns); constancy across points confirms."""
import os, sys, json
import numpy as np
import mpmath as mp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s2_core as C
import d4_y1a as Y
mp.mp.dps = 15
def That(p1, p2):
    p1 = [mp.mpf(int(x)) for x in p1]; p2 = [mp.mpf(int(x)) for x in p2]
    a11 = sum(x * x for x in p1); a22 = sum(x * x for x in p2); a12 = sum(x * y for x, y in zip(p1, p2))
    out = []
    for nu in range(4):
        def f(x, t):
            y = (1 - x) * t; z = (1 - x) * (1 - t)
            D = y * a22 + z * a11 - (y * y * a22 - 2 * y * z * a12 + z * z * a11)
            return -x * (y * p2[nu] - z * p1[nu]) / D ** 2 * (1 - x)
        out.append(mp.quad(f, [0, 1], [0, 1]))
    return out
res = []
pts = [((1, 2, 0, 1), (0, -1, 1, 2)), ((2, -1, 1, 0), (1, 1, -2, 1)), ((1, 1, 1, -1), (-2, 0, 1, 1))]
if C.ALG == 'gl3': pts = pts[:1]
for (p1, k) in pts:
    p1 = np.array(p1); k = np.array(k); p2 = -k - p1
    w1, w2 = 4, 5
    a = np.zeros((10, C.D)); a[w1] = C.traceless(np.array([1, 2, -1, 0, 1, 0, 2, 0, -1][:C.D], float)); a[w2] = C.traceless(np.array([0, 1, 1, 2, 0, -1, 1, 1, 0][:C.D], float))
    b = np.zeros((10, C.D)); b[w1] = C.traceless(np.array([2, 0, 1, 1, 0, -1, 0, 1, 1][:C.D], float)); b[w2] = C.traceless(np.array([1, -1, 2, 0, 1, 1, -1, 0, 2][:C.D], float))
    r = Y.y1a(p1, k, a, b, w1, w2, dps=12)
    beta = [complex(x) / r['T_O'] for x in r['beta']]
    Bh = [-mp.mpf(int(x)) / int(p2 @ p2) for x in p2]
    Th = That(p1, p2)
    A = np.array([[float(Bh[m]), float(int(p1 @ p1) * Th[m])] for m in range(4)])
    yv = np.array([c.real for c in beta])
    lam, resid, rk, sv = np.linalg.lstsq(A, yv, rcond=None)
    fit = A @ lam
    res.append({'p1': p1.tolist(), 'k': k.tolist(), 'p2': p2.tolist(), 'beta/T_O': [str(c) for c in beta],
                'lambda_1': lam[0], 'lambda_2': lam[1], 'max_residual': float(np.max(np.abs(fit - yv))),
                'alpha': str(r['alpha'])})
    print(json.dumps(res[-1]), flush=True)
json.dump(res, open(os.path.join(HERE, 'logs', 'd4_y1a_ward_%s.json' % C.ALG), 'w'), indent=1)
