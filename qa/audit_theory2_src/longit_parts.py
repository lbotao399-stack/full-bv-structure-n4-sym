"""Per-contribution log(eps) coefficients and finite parts of p.Pi.p (same rules as longit_selfenergy.py):
tests that the log cancellation and the totals depend on the relative signs (gauge/scalar bubble, FP ghost, fermions)."""
import sympy as sp
import longit_selfenergy as M
x, T, eps, P2 = M.x, M.T, M.eps, M.P2
def total_of(expr):
    avg = M.average(expr); A0 = sp.expand(avg.coeff(T, 0)); A1 = sp.expand(avg.coeff(T, -1))
    c = x*(1-x)*P2
    def piece(T0):
        return A0*(-sp.log(c*T0) - sp.EulerGamma) + A1*(1/T0 + c*sp.log(c*T0) + (sp.EulerGamma - 1)*c)
    lo = sp.integrate(sp.expand(sp.expand_log(piece(eps/x), force=True)), (x, 0, sp.Rational(1, 2)))
    hi = sp.integrate(sp.expand(sp.expand_log(piece(eps/(1-x)), force=True)), (x, sp.Rational(1, 2), 1))
    return sp.collect(sp.expand(sp.expand_log(lo + hi, force=True)), [eps, sp.log(eps), sp.log(P2)])
for name, e in [('gauge+scalar bubble', M.gl), ('FP ghost', M.gh), ('2 Dirac fermions', M.fe)]:
    print(name, ':', total_of(e))
print('seagull :', (M.D - 1)*P2/eps)
