"""Engine A: D_1 and P_1 on local operators per Thm thm:tl-D1 / Prop prop:tl-P1 (theta, theta_P kept as parameters).
  D_1 = 4 lc^natural - {F'_Delta, -} + C_Delta ;  two-body: 4 K^log(X_a,X_b) + C^(2)_Delta(X_a,X_b) ;
  one-body: -4u n_F(X) X - {F'_Delta, X} + C^(1)_Delta(X).
  P_1 = -{hbar F_1^(1)(P_b), -} + C_{P_b} ; no log-kernel term ; two-body part C^(2)_{P_b}.
Computed by engine A: K^log (all kernels), lc^natural, the spectra.  Inputs from the manuscript (Ch9) and the theory
(not recomputed here): F'_Delta (Thm c9-quadratic (iii)), Pi_Delta|perp (Rem c9-Pi), the composite triangle along Delta
and along P_b on the scalar bilinears (Prop c9-triangle, Thm c9-P1), and prediction P7 for C_Delta on cohomology."""
from fractions import Fraction
import sympy as sp

th, thP = sp.symbols('theta theta_P')
NF = {'c': 1, 'A_V': 2, 'A_perp': 2, 'lambda': 4, 'A*': -2, 'lambda*': -4, 'c*': -1}
# {F'_Delta, X} in units, from F'_Delta = u[-2/3 int(A,A*) + 1/3 int(lambda*,lambda) + 4/3 int(c*,c)]
#   + theta u int(2E(A_perp) + 2 A_perp, A*_perp), with {int(A,A*),A} = A, {int(lambda*,lambda),lambda} = -lambda,
#   {int(c*,c),c} = -c (Lemma c4-count with the field-counting multipliers 2, 4, 1), {int(v.A_perp,A*_perp),X} = v.X
FPRIME_ON = {'A_V': (-Fraction(2, 3), 0), 'A_perp': (-Fraction(2, 3), th), 'lambda': (-Fraction(1, 3), 0),
             'c': (-Fraction(4, 3), 0)}


def one_body():
    """one-body part of D_1 on a letter X: coefficient of X and of v.X, in units"""
    out = {}
    for X, (a, b) in FPRIME_ON.items():
        cX = -4 * NF[X] - a
        cv = -b
        comp = None
        if X == 'A_perp':
            comp = (-Fraction(2, 3), Fraction(5, 3))      # {Pi_Delta, X} = u(5/3 v.X - 2/3 X)  (Rem c9-Pi)
            cX += comp[0]; cv = cv + comp[1]
        out[X] = {'X': str(sp.nsimplify(cX)), 'v.X': str(sp.nsimplify(cv)),
                  'composite_one_body': 'u(5/3 v.X - 2/3 X) (Rem c9-Pi)' if comp else
                  'C^(1)_Delta non-local (log|p|^2, Rem c9-Pi); only C^(1)+C^(2) is local'}
    return out


def d1_bilinears(lc_traceless, lc_trace):
    """assemble D_1 on the scalar bilinears from the engine's lc and the Ch9 composite data"""
    four_lc = {'traceless': 4 * lc_traceless, 'trace': 4 * lc_trace}
    one_X = -8 * 2                                   # two letters A_perp, -8 u each (incl. -4u n_F and Pi, F')
    c2 = -2                                          # composite triangle: -2 u {S1(Delta),O} (Prop c9-triangle)
    v_coef = sp.nsimplify(c2 + Fraction(5, 3) - th)
    return {k: {'O': str(v + one_X), 'S1(Delta)': str(v_coef)} for k, v in four_lc.items()}


def p1_bilinears():
    c2 = -2                                          # composite triangle along P_b, symmetrised (Thm c9-P1 proof)
    one = Fraction(5, 3) - thP                       # composite two-point along P_b on a leg + (-theta_P) from F_1(P_b)
    return {'S1(P_b)': str(sp.nsimplify(c2 + one))}
