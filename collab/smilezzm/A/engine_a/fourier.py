"""Engine A, module 2: Fourier transforms of log|q|^2 / |q|^{2k} at x != 0 (the position-space
images of the eps -> 0 momentum-space expansions), derived from

    int d^4q/(2pi)^4 e^{i q.x} |q|^{2a}     = 4^a Gamma(2+a) / (pi^2 Gamma(-a))  |x|^{-4-2a},
    int d^4q/(2pi)^4 e^{i q.x} |q|^{2a-2}   = 4^a Gamma(1+a) / (4 pi^2 Gamma(1-a)) |x|^{-2-2a},

valid for x != 0 by continuation in the exponent a at fixed dimension four (no change of dimension;
the formulas are the Hankel transform int_0^inf q^{2a+3} J_1(q r)/(q r) dq / (4 pi^2) evaluated
in the strip where it converges and continued in a).  Expanding in a gives the transforms of
log^j |q|^2 |q|^{-2k}; polynomials in q (contact terms) transform to zero at x != 0.
The expansions are checked numerically by a direct radial (Hankel) quadrature with a Gaussian
damping factor in check_fourier().
"""
import sympy as sp

a = sp.Symbol('a')
x2 = sp.Symbol('x2', positive=True)   # |x|^2
LX = sp.Symbol('LX')                  # log(x2)
g_E = sp.EulerGamma


def _series(expr, order):
    return sp.series(expr, a, 0, order).removeO()


def ft_table():
    """Dictionary of transforms (coefficients of |x|^{-4} resp. |x|^{-2}), as sympy expressions in LX."""
    f0 = 4**a * sp.gamma(2 + a) / (sp.pi**2 * sp.gamma(-a)) * sp.exp(-a * LX)      # times |x|^{-4}
    f1 = 4**a * sp.gamma(1 + a) / (4 * sp.pi**2 * sp.gamma(1 - a)) * sp.exp(-a * LX)  # times |x|^{-2}
    s0 = sp.expand(_series(f0, 3))
    s1 = sp.expand(_series(f1, 3))
    tab = {
        'log_q2__over_1': sp.simplify(s0.coeff(a, 1)),            # FT(log q^2)        = c * |x|^{-4}
        'log2_q2__over_1': sp.simplify(2 * s0.coeff(a, 2)),       # FT(log^2 q^2)      = c * |x|^{-4}
        '1__over_q2': sp.simplify(s1.coeff(a, 0)),                # FT(1/q^2)          = c * |x|^{-2}
        'log_q2__over_q2': sp.simplify(s1.coeff(a, 1)),           # FT(log q^2 / q^2)  = c * |x|^{-2}
        'log2_q2__over_q2': sp.simplify(2 * s1.coeff(a, 2)),      # FT(log^2 q^2/q^2)  = c * |x|^{-2}
    }
    return {k: sp.expand(sp.expand_func(v)) for k, v in tab.items()}


def check_fourier(verbose=True):
    """Direct numerical check at |x| = 1 (so LX = 0) and |x| = 2: the radial transform
    int d^4q/(2pi)^4 e^{iqx} f(|q|) = (1/(4 pi^2 r)) int_0^inf q^2 J_1(q r) f(q) dq, computed with a
    damping factor e^{-eta q^2} and extrapolated eta -> 0 (the damping is the regulator itself)."""
    import mpmath as mp
    mp.mp.dps = 20
    tab = ft_table()
    out = {}
    for r in [mp.mpf(1), mp.mpf(2)]:
        LXv = 2 * mp.log(r)
        for key, f in [('log_q2__over_q2', lambda q: mp.log(q**2) / q**2),
                       ('1__over_q2', lambda q: 1 / q**2),
                       ('log_q2__over_1', lambda q: mp.log(q**2)),
                       ('log2_q2__over_q2', lambda q: mp.log(q**2)**2 / q**2)]:
            vals = []
            for eta in [mp.mpf('1e-3'), mp.mpf('5e-4'), mp.mpf('2.5e-4')]:
                g = lambda q: q**2 * mp.besselj(1, q * r) * f(q) * mp.exp(-eta * q**2) / (4 * mp.pi**2 * r)
                vals.append(mp.quadosc(g, [0, mp.inf], period=2 * mp.pi / r))
            # Richardson in eta (linear): v(eta) ~ v0 + c eta
            v0 = 2 * vals[2] - vals[1]
            pred = tab[key].subs(LX, sp.Float(str(LXv), 20))
            pw = -4 if key.endswith('over_1') else -2
            predv = mp.mpf(str(sp.N(pred, 20))) * r**pw
            out[(key, str(r))] = (str(v0), str(predv))
            if verbose:
                print('FT[%s] at |x|=%s: numeric %s  formula %s' % (key, r, mp.nstr(v0, 10), mp.nstr(predv, 10)))
    return out


if __name__ == '__main__':
    for k, v in ft_table().items():
        print(k, '=', v)
    check_fourier()
