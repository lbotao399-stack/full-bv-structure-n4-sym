"""Engine A, part C: D_1 and P_1 on the families with every free direction symbolic (Def def:cc-free):
  hbar F_1^(1)(Delta) = 4uF + F'_Delta, F'_Delta = F'^tr + u sum_d theta_d S_d + u Q(kappa_V int(iota_x A|_V,c*) + kappa_lambda int(rho(x)lambda*,lambda*)),
  hbar F_1^(1)(P_b)  = u sum_d thetaP_d S_d(P_b) + u Q(kappaP ...),   d in {perp, lambda, V}.
One-body terms on a field letter X = d^alpha phi of sector d (Lemma lem:cc-onebody with sector-wise theta):
  -{F'_Delta, X} = u (c_phi - 2 theta_d w_X) X - 2 u theta_d x^b d_b X,   c_{A|_V} = c_{A_perp} = 2/3, c_lambda = 1/3;
  -{hbar F_1^(1)(P_b), X} = -u thetaP_d d_b X;  the kappa directions are s-exact and give s{G, O}: zero on cohomology.
Two-body terms from the engine: on every family run, sum C^(2),0_Delta = -2/3 u N and sum C^(2)_{P_b} = -1/3 u d_b on
cohomology (the computed matrices), so
  D_1 = 4 lc^natural - u sum_a (theta_{d_a} + 1/3) v.X_a ,     P_1 = -u sum_a (thetaP_{d_a} + 1/3) d_b X_a   (on cohomology),
with v.X = 2(x.d + w_X) X.  For equal theta: D_1 = 4 lc^natural - (1/3 + theta) u v., P_1 = -(1/3 + theta_P) u d_b."""
import sympy as sp

th = {d: sp.Symbol('theta_' + d) for d in ('perp', 'lambda', 'V')}
thP = {d: sp.Symbol('thetaP_' + d) for d in ('perp', 'lambda', 'V')}


def family_entry(name, letters_desc, sectors, wsum_expr, spectrum_note):
    """sectors: list of (sector, total Weyl weight of the two-field letters in it)"""
    D0 = -2 * sum((th[d] + sp.Rational(1, 3)) * w for d, w in sectors)
    Db = -2 * sum((th[d] + sp.Rational(1, 3)) for d, w in sectors if True)
    return {'family': name, 'letters (two-field part)': letters_desc,
            'D_1 translation-invariant part on cohomology': '4 lc^natural + (%s) u' % sp.simplify(D0),
            'D_1 position part on cohomology': 'sum_a -2u (theta_{d_a} + 1/3) d_b X_a  (x^b-linear)',
            'P_1 on cohomology': 'sum_a -u (thetaP_{d_a} + 1/3) d_b X_a',
            'kappa directions': 's-exact (zero on cohomology)',
            'theta = -1/3 (all sectors)': 'Xi_Delta = 0, D_1 = 4 lc^natural: ' + spectrum_note,
            'P_1 = 0 on cohomology iff': ' and '.join('thetaP_%s = -1/3' % d for d in sorted(set(d for d, w in sectors)))}


def run(N_twist=4):
    out = {}
    N = sp.Symbol('N')
    out['F0'] = family_entry('scalar bilinears / planar scalars / su(2)', 'A_perp letters (L of them, no derivatives)', [('perp', sp.Symbol('L'))],
                             None, 'spectra 4 x lc^natural (Heisenberg), bilinears 0 / -48 u')
    out['F1'] = family_entry('twist two', 'd_n^k Z d_n^{N-k} Z (weight 2 + N)', [('perp', 2 + N)], None, '-32 h(S) u, S even <= N')
    out['F2_su11'] = family_entry('su(1|1)', 'Z (weight 1), psi (weight 3/2): L - M and M of them',
                                  [('perp', sp.Symbol('L') - sp.Symbol('M')), ('lambda', sp.Rational(3, 2) * sp.Symbol('M'))], None, '-8 u E_H (ha_su11)')
    out['F2_fermion_sl2'] = family_entry('fermion sl(2)', 'd_n^k psi d_n^{N-k} psi (weight 3 + N)', [('lambda', 3 + N)], None, '-32 h(n+1) u, n odd')
    out['F3_gluon'] = family_entry('gluon sl(2)', 'd_n^k f d_n^{N-k} f, f = d A|_V (weight 4 + N)', [('V', 4 + N)], None, '-32 h(n+2) u, n even')
    out['note'] = ('With sector-dependent theta the two-field part sees only the sectors of its leading letters; the gauge completion '
                   '(A|_V letters, sector V) enters the three-field part, so gauge covariance of D_1 couples theta_V to theta_perp, '
                   'theta_lambda (Conjecture P7 (a), not tested here).  Translation-covariant scheme: thetaP = theta, P_1 = D_1^b / 2.')
    return out
