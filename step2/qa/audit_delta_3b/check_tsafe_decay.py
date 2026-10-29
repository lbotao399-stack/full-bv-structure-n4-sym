"""Audit #3: large-|y| order of one insertion of S^(1)(x) on a cycle (ch08 prop:c8-Tsafe and the CH08 item of
delta_theory: 'an insertion of S^(1)(Delta) between an A x A edge and a c x A*|_V edge is of order |y|^{-4} at leading
order and converges only because that leading term is odd in y').

Insertion at y of the term (v.A, A^*) of S^(1)(x) (def:c7-vdot): its A leg is contracted with the A (x) A component of
P(0,infty) from a vertex at z1, scalar factor 1/|y - z1|^2 (prop:c3-P, constants dropped), acted on by v. (on a scalar
component v.d + (1/4) div v); its A^* leg is contracted with the A^*|_V (x) c component from a vertex at z2, factor
d_{y^mu} 1/|y - z2|^2.  y = n/u, expansion in u -> 0 (u = 1/|y| for |n| = 1); leading power and parity in n.
Cases: Delta (v = 2y, 1/4 div v = 2), K_b (v = |y|^2 b - 2(b.y) y, 1/4 div v = -2 b.y), translation (v = b).
"""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
u = sp.symbols('u', positive=True)
Y = sp.symbols('y0:4')
NN = [(sp.Rational(3, 5), sp.Rational(4, 5), 0, 0), (sp.Rational(2, 7), sp.Rational(3, 7), sp.Rational(6, 7), 0)]
ZZ1 = (1, 0, 2, -1); ZZ2 = (0, 1, -1, 3); BV = (1, -2, 3, 1)


def norm2(w):
    return sum(t * t for t in w)


def lead(expr, n):
    e = sp.together(expr.subs({Y[m]: n[m] / u for m in range(4)}))
    s = sp.series(e, u, 0, 9).removeO()
    s = sp.expand(s)
    if s == 0:
        return None, 0
    p = min(sp.Poly(s, u).monoms())[0]
    return p, sp.Poly(s, u).coeff_monomial(u ** p)


def main():
    res = {}
    PA = 1 / norm2([Y[m] - ZZ1[m] for m in range(4)])
    PAc = [sp.diff(1 / norm2([Y[m] - ZZ2[m] for m in range(4)]), Y[mu]) for mu in range(4)]
    yb = sum(Y[m] * BV[m] for m in range(4)); yy = norm2(Y)
    cases = {'Delta': ([2 * Y[m] for m in range(4)], 2),
             'K_b': ([yy * BV[m] - 2 * yb * Y[m] for m in range(4)], -2 * yb),
             'translation': (list(BV), 0)}
    for name, (v, wt) in cases.items():
        Aleg = sum(v[m] * sp.diff(PA, Y[m]) for m in range(4)) + wt * PA
        out = {}
        for i, n in enumerate(NN):
            pa, ca = lead(Aleg, n)
            out['n%d: A-leg factor leading power of |y|' % i] = -pa
            comps = []
            for mu in range(4):
                p, cp = lead(Aleg * PAc[mu], n)
                pm, cm = lead(Aleg * PAc[mu], [-t for t in n])
                comps.append({'power of |y|': (-p if p is not None else None), 'coeff(n)': str(cp), 'coeff(-n)': str(cm)})
            out['n%d: product components' % i] = comps
        res[name] = out
    # the composite line needs the A leg in A|_V: full Lie derivative on the one-form index (rotation part included)
    for name, (v, wt) in cases.items():
        out = {}
        for i, n in enumerate(NN):
            comps = []
            for lam in range(4):
                # (L_v P)_{mu lam} = v.d P delta_{mu lam} + (d_mu v^lam) P, contracted with d_mu 1/|y - z2|^2
                tot = 0
                for mu in range(4):
                    Lv = (sum(v[m] * sp.diff(PA, Y[m]) for m in range(4)) if mu == lam else 0) + sp.diff(v[lam], Y[mu]) * PA
                    tot += Lv * PAc[mu]
                p, cp = lead(tot, n)
                pm, cm = lead(tot, [-t for t in n])
                comps.append({'power of |y|': (-p if p is not None else None), 'coeff(n)': str(cp), 'coeff(-n)': str(cm)})
            out['n%d: A|_V leg, product components (index at z1)' % i] = comps
        res[name + ' (A|_V, Lie derivative)'] = out
    res['reading'] = ('Leading large-|y| power of one insertion between an A x A edge (acted on by v.) and an '
                      'A^*|_V x c edge; prop:c8-Tsafe states -6 (constant fields), -5 (Delta), -4 (K).  The composite '
                      'line has its A leg in A|_V, where the Lie derivative carries the rotation part of dv.')
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(HERE, 'check_tsafe_decay.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
