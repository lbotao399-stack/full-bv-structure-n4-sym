"""Independent audit check (theory round 2): longitudinal part p_mu Pi_munu(p) p_nu of the one-loop A|_V self-energy
of the dimensionally reduced 10d SYM (N=4) in Feynman gauge, every propagator multiplied by exp(-eps k^2)
(= P(eps,infty) of the manuscript), written with standard Euclidean Feynman rules (Faddeev-Popov ghosts,
2 Dirac adjoint fermions), canonical normalisation A = g Ahat.  Gamma_1 = 1/2 int Ahat(-p) Pi(p) Ahat(p).
Prediction from compare_r2.json (T2): Pi_munu contains (9/2) eps^{-1} delta_munu + (27/8) p_mu p_nu (units g^2 C_A/16pi^2),
i.e. p.Pi.p = (9/2) p^2/eps + (27/8) p^4 + O(eps), with no log.
"""
import sympy as sp
x, T, eps, a, P2 = sp.symbols('x T epsilon a p2', positive=True)
lp, l2 = sp.symbols('lp l2')       # l.p and l.l  (k = l - x p)
D = 10                              # internal A components (4 + 6)
nD = 2                              # Dirac adjoint fermions

kk = l2 - 2*x*lp + x**2*P2          # k.k
kp = lp - x*P2                      # k.p
# (k+p) = l + (1-x) p
Kk = l2 + 2*(1-x)*lp + (1-x)**2*P2  # (k+p).(k+p)
Kp = lp + (1-x)*P2                  # (k+p).p
kK = l2 + (1-2*x)*lp - x*(1-x)*P2   # k.(k+p)

# gluon/scalar bubble: Z = p (x) a + b (x) p - s delta, a = k - p, b = 2p + k, s = 2 p.k + p^2
pa = kp - P2; pb = 2*P2 + kp
aa = kk - 2*kp + P2; bb = 4*P2 + 4*kp + kk
s = 2*kp + P2
L = P2*aa + 2*pb*pa + bb*P2 - 2*s*(pa + pb) + s**2*D
gl = -sp.Rational(1, 2)*L
gh = kp*Kp
fe = 4*nD*(2*Kp*kp - P2*kK)
num = sp.expand(gl + gh + fe)

def average(expr):
    """Gaussian average over l with weight exp(-T l^2): <l.p>=0, <(l.p)^2>=p^2/(2T), <l^2>=2/T (4d); degree <= 2 in l"""
    poly = sp.Poly(sp.expand(expr), lp, l2)
    out = 0
    for (i, j), c in poly.terms():
        if (i, j) == (0, 0): out += c
        elif (i, j) == (2, 0): out += c*P2/(2*T)
        elif (i, j) == (0, 1): out += c*2/T
        elif i % 2 == 1: pass
        else: raise ValueError((i, j))
    return sp.expand(out)

avg = average(num)          # = A0(x) + A1(x)/T
A0 = sp.expand(avg.coeff(T, 0)); A1 = sp.expand(avg.coeff(T, -1))
assert sp.simplify(avg - A0 - A1/T) == 0
print('A0 =', sp.factor(A0)); print('A1 =', sp.factor(A1))
# J = (1/16pi^2) int_0^1 dx int_{T0}^inf dT [A0/T + A1/T^2] exp(-T x(1-x) p^2),  T0 = eps/min(x,1-x)
# small-T0 expansions:  int_{T0} dT/T e^{-cT}   = -log(c T0) - EulerGamma + O(T0)
#                       int_{T0} dT/T^2 e^{-cT} = 1/T0 + c log(c T0) + (EulerGamma - 1) c + O(T0 log T0)
c = x*(1-x)*P2
def piece(T0):
    return A0*(-sp.log(c*T0) - sp.EulerGamma) + A1*(1/T0 + c*sp.log(c*T0) + (sp.EulerGamma - 1)*c)
lo = sp.integrate(sp.expand(sp.expand_log(piece(eps/x), force=True)), (x, 0, sp.Rational(1, 2)))
hi = sp.integrate(sp.expand(sp.expand_log(piece(eps/(1-x)), force=True)), (x, sp.Rational(1, 2), 1))
bubbles = sp.simplify(sp.expand(lo + hi))
seagull = (D - 1)*P2/eps         # (D-1) I_1 p^2, I_1 = 1/(16 pi^2 eps)
total = sp.expand(sp.expand_log(bubbles + seagull, force=True))
print('bubbles  (units g^2 C_A/16pi^2):', sp.collect(sp.expand(sp.expand_log(bubbles, force=True)), [eps, sp.log(eps), sp.log(P2)]))
print('total p.Pi.p (units g^2 C_A/16pi^2):', sp.collect(total, [eps, sp.log(eps), sp.log(P2)]))
pred = sp.Rational(9, 2)*P2/eps + sp.Rational(27, 8)*P2**2
print('total - prediction =', sp.simplify(total - pred))
