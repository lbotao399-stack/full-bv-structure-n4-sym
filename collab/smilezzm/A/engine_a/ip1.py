"""Engine A, module 3: the one-loop two-point function of the scalar bilinears (IP1), exact.

Method (momentum space, reduction to one-loop objects):
  * Feynman rules of the dimensionally reduced 10d Yang-Mills theory (coupling g, su(N), Feynman
    gauge, weight e^{-S_E} with S_E = int 1/4 F^a_MN F^a_MN + 1/2 (d_m A^a_m)^2 + fermion), all momenta
    in R^4, metric split delta_mn (+) delta_IJ; Majorana-Weyl gaugino (loops: 1/2 x 16-dim trace).
    Vertex factors (all momenta incoming):
      V3^{abc}_{MNP}(k1,k2,k3) = i g f^{abc} [delta_MN (k1-k2)_P + delta_NP (k2-k3)_M + delta_PM (k3-k1)_N]
      V4^{abcd}_{MNPQ} = -g^2 [ f^{abe}f^{cde}(d_MP d_NQ - d_MQ d_NP) + f^{ace}f^{bde}(d_MN d_PQ - d_MQ d_NP)
                                + f^{ade}f^{bce}(d_MN d_PQ - d_MP d_NQ) ]
    Propagator delta^{ab} delta_MN e^{-eps k^2}/k^2.  Physical normalisation hbar = 1; the operators
    are O_{ab} = (a.A^c)(b.A^c) (sum over colour c).  Only ratios to the tree enter the results.
  * The unit: unit = -g^2 N/(32 pi^2), i.e. g^2 N/(16 pi^2) = -2 unit.
  * The one-loop graphs: the self-energy on a line (exact eps-expansion of the four 1PI graphs through
    masters.py, singular part subtracted = the counterterm -I^CT, the beta_2 insertion), the exchange
    graph (reduced EXACTLY through (k+k').(2q-k-k') = 2q^2 - k^2 - (q-k)^2 - k'^2 - (q-k')^2 + (k-k')^2
    and |k|^2 Dhat(k) = e^{-eps k^2}: a finite five-line graph times 2q^2 (contact term at x != 0) and
    graphs with one edge replaced by the heat kernel K_eps, which localise as eps -> 0 to products of
    one-loop masters plus the localisation constant c_loc = int d^4k/(2pi)^4 Dhat_1(k)^2 (e^{-k^2} - 1)),
    the quartic graph (a product of two bubbles), and the operator-counterterm trees.
  * Position space at x != 0 through the Fourier lemma (fourier.py).
All numbers exact (sympy); beta_2 symbolic.
"""
import os
import json
import sympy as sp
from masters import M, one_edge, gauss_average, bubble_scalar, p2, eps, LE, b2, g_E
import fourier

HERE = os.path.dirname(os.path.abspath(__file__))
unit = sp.Symbol('unit')          # the unit  hbar tr_g(ad ad)/(16 pi^2 (X,Y)) = -g^2 N/(32 pi^2)
LX = fourier.LX                   # log(x^2 / ell^2)
L2, L3 = sp.log(2), sp.log(3)
t1, t2 = sp.symbols('t1 t2', positive=True)
s = sp.Symbol('s', positive=True)

# ---------------------------------------------------------------------------------------------
# 1. The scalar self-energy Sigma(p): sum of the 1PI one-loop graphs with two external A_I legs,
#    in the convention  dressed line = Dhat + Dhat Sigma Dhat + ...,  i.e. Sigma = sum of graphs
#    with the (-S_E) vertex factors, propagators, symmetry factors and (-1) per fermion loop.
#    Colour and flavour: every graph is proportional to delta^{ab} delta_IJ; we quote the
#    coefficient of delta^{ab} delta_IJ, in units of g^2 N (symbol gN2 := g^2 N).
# ---------------------------------------------------------------------------------------------
gN2 = sp.Symbol('g2N', positive=True)
lv = sp.symbols('l1:5')
pv = sp.symbols('pc1:5')   # momentum components (pc_i); p2 = |p|^2 is a separate symbol
lsq = sum(x**2 for x in lv)
psq = sum(x**2 for x in pv)
ldotp = sum(x * y for x, y in zip(lv, pv))


def _finish(bracket):
    """bracket from gauss_average -> exact eps-expansion (with p^2 -> p2 symbol)."""
    expr = sp.expand(bracket.subs({s: t1 + t2}))
    # rotation invariance: the components enter only through |p|^2; eliminate pc4^2 and check
    expr = sp.expand(expr.subs(pv[3]**2, p2 - pv[0]**2 - pv[1]**2 - pv[2]**2))
    assert not (expr.free_symbols & set(pv)), 'non-invariant dependence on the momentum components'
    # make the s-dependence explicit again for bubble_scalar: rewrite in t1, t2 and powers of 1/(t1+t2)
    expr = sp.together(expr)
    num, den = sp.fraction(expr)
    # den is a power of (t1+t2) times a number
    den_poly = sp.Poly(den, t1, t2)
    # determine power of s in den
    k = 0
    d = den
    while True:
        q, r = sp.div(sp.Poly(d, t1, t2), sp.Poly(t1 + t2, t1, t2))
        if r.is_zero and not q.is_zero:
            k += 1
            d = q.as_expr()
        else:
            break
    const = sp.simplify(d)
    bracket_t = sp.expand(num / const) * sp.Symbol('s', positive=True)**(-k)
    return bubble_scalar(bracket_t)


def self_energy_graphs():
    """Returns dict graph -> exact eps-expansion of its contribution to Sigma (coefficient of
    delta^{ab} delta_IJ), as sympy expressions in eps, LE, p2, g2N."""
    out = {}
    # seagull with a 4d gluon loop: V4 contracted on legs 3,4 with delta^{cd} delta_mm (m in 1..4):
    #   -g^2 N delta^{ab} delta_IJ * 2 * 4 ; symmetry factor 1/2 ; tadpole integral 1/(16 pi^2 eps)
    out['seagull_V'] = sp.Rational(1, 2) * (-gN2 * 2 * 4) * one_edge()
    # seagull with a scalar loop (M in 1..6): 2*(6-1)
    out['seagull_perp'] = sp.Rational(1, 2) * (-gN2 * 2 * 5) * one_edge()
    # bubble: one gluon line (momentum p - l), one scalar line (momentum l):
    #   V1 V2 = (ig)^2 f^{adc} f^{bdc} delta_IJ (p+l).(-(p+l)) = + g^2 N delta^{ab} delta_IJ (p+l)^2 ; symmetry 1
    num_b = gN2 * (psq + 2 * ldotp + lsq)
    out['bubble'] = _finish(gauss_average(num_b, lv, pv, t1, t2))
    # gaugino loop: (-1)(1/2) g^2 (colour -N)(-i)^2 tr_16[Gamma^I Gamma.l Gamma^J Gamma.(l-p)] / (l^2 (l-p)^2)
    #   tr_16 = 16 (d^{IK} d^{JL} - d^{IJ} d^{KL} + d^{IL} d^{KJ}) -> -16 delta_IJ l.(l-p)
    #   => Sigma_f = + 8 g^2 N delta^{ab} delta_IJ int l.(l-p) Dhat(l) Dhat(l-p)
    num_f = 8 * gN2 * (lsq - ldotp)
    out['gaugino'] = _finish(gauss_average(num_f, lv, pv, t1, t2))
    return out


def singular_part(expr):
    """Terms with eps^{-1} or LE (log eps), as a polynomial in LE and 1/eps."""
    e = sp.expand(expr)
    sing = e.coeff(eps, -1) / eps + e.coeff(LE, 1) * LE
    return sp.expand(sing)


def self_energy():
    graphs = self_energy_graphs()
    total = sp.expand(sum(graphs.values()))
    sing = singular_part(total)
    ren = sp.expand(total - sing)
    return graphs, total, sing, ren


# ---------------------------------------------------------------------------------------------
# 2. Position space: the dressed line.  Sigma_ren(k) = g2N/(16 pi^2) k^2 [A + B log k^2]  (exactly,
#    as eps -> 0); the line Dhat(k) Sigma_ren(k) Dhat(k) -> (A + B log k^2)/k^2 -> position space.
# ---------------------------------------------------------------------------------------------
def line_correction(sigma_ren):
    """FT[Sigma_ren(k)/k^4](x) divided by D_0(x) = 1/(4 pi^2 x^2): an expression in LX, g2N."""
    e = sp.expand(sigma_ren / p2)
    A = e.coeff(sp.log(p2), 0)
    B = e.coeff(sp.log(p2), 1)
    assert sp.simplify(e - A - B * sp.log(p2)) == 0, 'unexpected p-dependence in Sigma_ren'
    tab = fourier.ft_table()
    D0 = tab['1__over_q2']                 # 1/(4 pi^2) times |x|^{-2}
    FTlog = tab['log_q2__over_q2']         # times |x|^{-2}
    return sp.expand(sp.simplify((A * D0 + B * FTlog) / D0))


# ---------------------------------------------------------------------------------------------
# 3. The masters in position space (from the bubble B(k) and the Fourier lemma):
#    B(k) = int d^4l/(2pi)^4 Dhat(l) Dhat(k-l) = M[0,0,2]/(16 pi^2);  J(x) = int dz D_eps(x-z)^2 D_eps(z)^2
#    = FT[B^2](x);  R(x) = int dz D_eps(x-z)^2 D_eps(z) = FT[B Dhat](x).
# ---------------------------------------------------------------------------------------------
def masters_position():
    tab = fourier.ft_table()
    B = sp.expand(M(0, 0, 2) / (16 * sp.pi**2))          # c0 + c1 (LE + log p2)
    c1 = B.coeff(LE, 1)
    c0 = sp.expand(B - c1 * (LE + sp.log(p2)))
    assert sp.simplify(c0.subs(p2, 1) - c0) == 0
    # B = c0 + c1 (LE + log p2).  B^2 = c0^2 + 2 c0 c1 LE + c1^2 LE^2 + (2 c0 c1 + 2 c1^2 LE) log p2 + c1^2 log^2 p2
    J_over_x4 = (2 * c0 * c1 + 2 * c1**2 * LE) * tab['log_q2__over_1'] + c1**2 * tab['log2_q2__over_1']
    # B Dhat = (c0 + c1 LE)/p2 + c1 log p2 / p2
    R_over_x2 = (c0 + c1 * LE) * tab['1__over_q2'] + c1 * tab['log_q2__over_q2']
    D0_over_x2 = tab['1__over_q2']
    J_rel = sp.expand(sp.simplify(J_over_x4 / D0_over_x2**2))   # J / D_0^2
    R_rel = sp.expand(sp.simplify(R_over_x2 / D0_over_x2))      # R / D_0
    # localisation constant c_loc = int d^4k/(2pi)^4 e^{-2k^2}(e^{-k^2}-1)/k^4 = (1/16pi^2) int_0^inf (e^{-3v}-e^{-2v})/v dv
    v = sp.Symbol('v', positive=True)
    c_loc = sp.integrate((sp.exp(-3 * v) - sp.exp(-2 * v)) / v, (v, 0, sp.oo)) / (16 * sp.pi**2)
    c_loc = sp.simplify(c_loc)
    return {'B_c0': c0, 'B_c1': c1, 'J_rel': J_rel, 'R_rel': R_rel, 'c_loc': c_loc}


# ---------------------------------------------------------------------------------------------
# 4. Flavour/colour tensors (explicit numerical check that only S12 and T12 occur)
# ---------------------------------------------------------------------------------------------
def flavour_structures(a1, b1, a2, b2):
    import numpy as np
    S = (a1 @ a2) * (b1 @ b2) + (a1 @ b2) * (b1 @ a2)
    T = (a1 @ b1) * (a2 @ b2)
    return S, T


def flavour_check(seed=1):
    """Build the explicit flavour tensors of each class by contracting the vertex tensors with
    random complex six-vectors; decompose into S12, T12 and report the residual."""
    import numpy as np
    rng = np.random.default_rng(seed)
    d = np.eye(6)
    res = {}
    for trial in range(3):
        a1, b1, a2, b2 = (rng.normal(size=6) + 1j * rng.normal(size=6) for _ in range(4))
        S, T = flavour_structures(a1, b1, a2, b2)
        # tree: sum over the two pairings
        tree = np.einsum('i,j,k,l,ik,jl->', a1, b1, a2, b2, d, d) + np.einsum('i,j,k,l,il,jk->', a1, b1, a2, b2, d, d)
        # exchange: 1/2 sum over leg assignments of delta_{I1 K1} delta_{I2 K2}: O1 legs (I1,I2), O2 legs (K1,K2)
        O1 = np.einsum('i,j->ij', a1, b1); O1 = O1 + O1.T
        O2 = np.einsum('i,j->ij', a2, b2); O2 = O2 + O2.T
        ex = 0.5 * np.einsum('ij,ij->', O1, O2)
        # quartic: V4 flavour tensor on legs (I1,I2,K1,K2), colour (c,c,d,d):
        #   -g^2 N dim g [2 d_{I1 I2} d_{K1 K2} - d_{I1 K2} d_{I2 K1} - d_{I1 K1} d_{I2 K2}], contracted with a1 b1 a2 b2
        qu = 2 * (a1 @ b1) * (a2 @ b2) - (a1 @ b2) * (b1 @ a2) - (a1 @ a2) * (b1 @ b2)
        # operator counterterm trees: -(4 O_1 - 2 (a1.b1) K) with O_2  and the mirror
        oc = -(4 * tree - 2 * (a1 @ b1) * 2 * (a2 @ b2)) - (4 * tree - 2 * (a2 @ b2) * 2 * (a1 @ b1))
        res['trial%d' % trial] = {
            'tree - S12': abs(tree - S), 'EX - S12': abs(ex - S), 'QU - (2 T12 - S12)': abs(qu - (2 * T - S)),
            'OC - (-8 S12 + 8 T12)': abs(oc - (-8 * S + 8 * T))}
    return res


# ---------------------------------------------------------------------------------------------
# 5. Assembly of the classes relative to the tree  4 hbar dim g S12 D_0^2  (coefficients of S12
#    and T12 in units of unit, plus LX = log(x^2/ell^2) terms and LE = log eps terms)
# ---------------------------------------------------------------------------------------------
def assemble():
    graphs, Sig, Sig_sing, Sig_ren = self_energy()
    pos = masters_position()
    J_rel, R_rel, c_loc = pos['J_rel'], pos['R_rel'], pos['c_loc']
    # conversion g2N -> unit:  g2N = -32 pi^2 unit
    to_unit = {gN2: -32 * sp.pi**2 * unit}

    classes = {}
    # (SE): two lines, each D_0 (1 + line_correction); relative to tree: 2 * line_correction (structure S12)
    lc = line_correction(Sig_ren)
    classes['SE'] = {'S': sp.expand(2 * lc.subs(to_unit)), 'T': sp.Integer(0)}
    # (B): beta_2 insertion: Sigma_F = -4 beta_2 unit k^2  ->  line -4 beta_2 unit D_0 ; two lines
    classes['B'] = {'S': sp.expand(-8 * b2 * unit), 'T': sp.Integer(0)}
    # (EX): -g^2 N dim g S12 FT[ iint N_EX Dhat^5 ]; relative to tree dim g S12 D_0^2:
    #   32 pi^2 unit [ 0 (contact) - 4 (R D_0 + D_0^2 c_loc) + J + 2 D_0^2 c_loc ] / D_0^2
    classes['EX'] = {'S': sp.expand(32 * sp.pi**2 * unit * (J_rel - 4 * R_rel - 2 * c_loc)), 'T': sp.Integer(0)}
    # (QU): -g^2 N dim g (2 T12 - S12) J  ->  32 pi^2 unit (2 T12 - S12) J_rel
    classes['QU'] = {'S': sp.expand(-32 * sp.pi**2 * unit * J_rel), 'T': sp.expand(64 * sp.pi**2 * unit * J_rel)}
    # (OC): -unit LE [ 4 <O1,O2> - 2 (a1.b1) <K,O2> ] and mirror:  <K,O2>_tree = 2 (a2.b2) dim g D_0^2
    classes['OC'] = {'S': sp.expand(-8 * unit * LE), 'T': sp.expand(8 * unit * LE)}
    classes['F'] = {'S': sp.Integer(0), 'T': sp.Integer(0)}

    totS = sp.expand(sum(c['S'] for c in classes.values()))
    totT = sp.expand(sum(c['T'] for c in classes.values()))

    def split(e):
        e = sp.expand(e / unit)
        return {'const': sp.expand(e.coeff(LX, 0).coeff(LE, 0)), 'log_x2': sp.expand(e.coeff(LX, 1)),
                'log_eps': sp.expand(e.coeff(LE, 1))}

    report = {'classes': {k: {'S12': split(v['S']), 'T12': split(v['T'])} for k, v in classes.items()},
              'total': {'S12': split(totS), 'T12': split(totT)},
              'self_energy': {'graphs': {k: sp.expand(v) for k, v in graphs.items()}, 'total': Sig,
                              'singular': Sig_sing, 'renormalised': Sig_ren},
              'masters': pos, 'line_correction': lc}
    return report


def main():
    rep = assemble()
    return rep


if __name__ == '__main__':
    rep = main()
    for k, v in rep['self_energy']['graphs'].items():
        print('Sigma_%s =' % k, v)
    print('Sigma_sing =', rep['self_energy']['singular'])
    print('Sigma_ren  =', rep['self_energy']['renormalised'])
    print('line corr  =', rep['line_correction'])
    print('masters:', rep['masters'])
    for k, v in rep['classes'].items():
        print(k, v)
    print('TOTAL', rep['total'])
