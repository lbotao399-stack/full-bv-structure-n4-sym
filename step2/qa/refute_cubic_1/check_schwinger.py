"""Check of the Schwinger-integral evaluation in loop.Quad against brute force:
 (1) Gaussian moment formula at fixed Schwinger times (4d Gauss-Hermite with many nodes);
 (2) full triangle and bubble integrals over [eps, infty)^n by scipy nquad in the original t variables."""
import numpy as np
from scipy import integrate
import ft
import loop

rng = np.random.default_rng(5)
c = rng.normal(size=35) + 1j * rng.normal(size=35)
cl = ft.lap_coeffs(c)
s = [rng.normal(size=4), rng.normal(size=4), np.zeros(4)]

# (1) Gaussian moments
t = np.array([0.37, 0.81, 0.55])
T = t.sum()
xg, wg = np.polynomial.hermite.hermgauss(8)
acc = 0
for i0 in range(8):
    for i1 in range(8):
        for i2 in range(8):
            for i3 in range(8):
                u = np.array([xg[i0], xg[i1], xg[i2], xg[i3]]) / np.sqrt(T)
                w = wg[i0] * wg[i1] * wg[i2] * wg[i3] / T ** 2
                kc = -(t[:, None] * np.array(s)).sum(0) / T
                k = kc + u
                expo = -sum(t[e] * (k + s[e]) @ (k + s[e]) for e in range(3)) + T * (k - kc) @ (k - kc)
                acc += w * loop.polyval(c, k[None, :])[0] * np.exp(expo)
brute = acc / (2 * np.pi) ** 4
x = t / T
kc = -(x[:, None] * np.array(s)).sum(0)
D = x[0] * x[1] * (s[0] - s[1]) @ (s[0] - s[1]) + x[0] * x[2] * (s[0] - s[2]) @ (s[0] - s[2]) + x[1] * x[2] * (s[1] - s[2]) @ (s[1] - s[2])
form = (loop.polyval(c, kc[None, :])[0] + 0.25 * loop.polyval(cl, kc[None, :])[0] / T) * np.exp(-T * D) / (16 * np.pi ** 2 * T ** 2)
print('(1) Gaussian moment formula: rel diff', abs(brute - form) / abs(form))

# (2) triangle over [eps,inf)^3 in t variables
eps = 0.05


def integrand3(a, b, cc):
    tt = np.array([eps + a / (1 - a), eps + b / (1 - b), eps + cc / (1 - cc)])
    jac = 1 / ((1 - a) ** 2 * (1 - b) ** 2 * (1 - cc) ** 2)
    T = tt.sum()
    x = tt / T
    kc = -(x[:, None] * np.array(s)).sum(0)
    D = x[0] * x[1] * (s[0] - s[1]) @ (s[0] - s[1]) + x[0] * x[2] * (s[0] - s[2]) @ (s[0] - s[2]) + x[1] * x[2] * (s[1] - s[2]) @ (s[1] - s[2])
    v = (loop.polyval(c, kc[None, :])[0] + 0.25 * loop.polyval(cl, kc[None, :])[0] / T) * np.exp(-T * D) / (16 * np.pi ** 2 * T ** 2)
    return v * jac


re = integrate.nquad(lambda a, b, cc: integrand3(a, b, cc).real, [[0, 1]] * 3, opts={'epsabs': 1e-11, 'epsrel': 1e-9, 'limit': 100})[0]
im = integrate.nquad(lambda a, b, cc: integrand3(a, b, cc).imag, [[0, 1]] * 3, opts={'epsabs': 1e-11, 'epsrel': 1e-9, 'limit': 100})[0]
q = loop.Quad()
mine = q.integrate(3, np.array(s), c, cl, eps)
print('(2) triangle: nquad', re + 1j * im, ' Quad', mine, ' rel diff', abs(re + 1j * im - mine) / abs(mine))

# bubble
s2 = [s[0], np.zeros(4)]


def integrand2(a, b):
    tt = np.array([eps + a / (1 - a), eps + b / (1 - b)])
    jac = 1 / ((1 - a) ** 2 * (1 - b) ** 2)
    T = tt.sum()
    x = tt / T
    kc = -x[0] * s2[0]
    D = x[0] * x[1] * s2[0] @ s2[0]
    v = (loop.polyval(c, kc[None, :])[0] + 0.25 * loop.polyval(cl, kc[None, :])[0] / T) * np.exp(-T * D) / (16 * np.pi ** 2 * T ** 2)
    return v * jac


re = integrate.nquad(lambda a, b: integrand2(a, b).real, [[0, 1]] * 2, opts={'epsabs': 1e-12, 'epsrel': 1e-11, 'limit': 200})[0]
im = integrate.nquad(lambda a, b: integrand2(a, b).imag, [[0, 1]] * 2, opts={'epsabs': 1e-12, 'epsrel': 1e-11, 'limit': 200})[0]
mine = q.integrate(2, np.array(s2), c, cl, eps)
print('(2) bubble: nquad', re + 1j * im, ' Quad', mine, ' rel diff', abs(re + 1j * im - mine) / abs(mine))
