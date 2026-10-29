"""Exact eps-asymptotics of one-loop Schwinger integrals with per-edge lower limit eps (engine A).
Two-edge graph:  J = int_{t1,t2>=eps} dt1 dt2 (16 pi^2)^{-1} T^{-2} [ g(u2) + c/T ] e^{-s t1 t2/T},  T=t1+t2, u2=t2/T.
Boundary method: dJ/deps = -sum_faces int_{t>=eps} G|_{face} ... = -eps^{-1-a} sum_n (-s eps)^n D_n/n!, with
D_n = int_1^infty [G(1,v)+G(v,1)] (v/(1+v))^n dv  (exact, absolutely convergent), hence
  a=0 (G=T^{-2}g):  log-eps coefficient = -D_0 = -int_0^1 g(u) du
  a=1 (G=c T^{-3}):  eps^{-1} coefficient = D_0 = c/4 ;  log-eps coefficient = s D_1 = s c/6.
One-edge graph:  int_eps^infty dt (16 pi^2)^{-1} t^{-2} h = h/eps.
n-edge graph of superficial degree 0: log-eps coefficient = -int_{simplex} (Laplacian part)/(4) du
(projective-integral form of the face sum; numerator constant => -(Lap/4) * vol(simplex) = -(Lap/4)/(n-1)!)."""
import sympy as sp
from fractions import Fraction

def D_n_check():
    """symbolic re-derivation of the constants 1/4 and 1/6 and of -int_0^1 g"""
    v, u = sp.symbols('v u', positive=True)
    G1 = lambda t1, t2: (t1 + t2) ** -3
    D0 = sp.integrate(G1(1, v) + G1(v, 1), (v, 1, sp.oo))
    D1 = sp.integrate((G1(1, v) + G1(v, 1)) * v / (1 + v), (v, 1, sp.oo))
    # a=0 with g(u2)=u2^m : faces
    res = {}
    for m in range(3):
        G0 = lambda t1, t2: (t1 + t2) ** -2 * (t2 / (t1 + t2)) ** m
        Dm = sp.integrate(sp.simplify(G0(1, v) + G0(v, 1)), (v, 1, sp.oo))
        res[m] = sp.nsimplify(Dm)
    return sp.nsimplify(D0), sp.nsimplify(D1), res

if __name__ == '__main__':
    print(D_n_check())
