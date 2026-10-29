"""Regression of the round-2 exact engine against the round-1 quadratic counterterms (exact rationals)."""
import json, sys
from fractions import Fraction
import numpy as np
import r2_core as C
from r2_core import Q2
from r2_graphs import GraphSum
rng = np.random.default_rng(101)

def rpol(t, which=None, lo=-1, hi=2):
    x = rng.integers(lo, hi, size=C.SHAPE[t]).astype(complex)
    if which == 'V': x[4:] = 0
    if which == 'P': x[:4] = 0
    return x

def solve_exact(rows, rhs):
    """exact least squares via normal equations over Q(i) (Fractions): rows list of lists of Q2, rhs list of Q2"""
    import sympy as sp
    M = sp.Matrix([[sp.Rational(r.re.numerator, r.re.denominator) + sp.I * sp.Rational(r.im.numerator, r.im.denominator) for r in row] for row in rows])
    b = sp.Matrix([sp.Rational(r.re.numerator, r.re.denominator) + sp.I * sp.Rational(r.im.numerator, r.im.denominator) for r in rhs])
    if all(x == 0 for x in b):
        return [sp.Integer(0)] * M.shape[1], True
    sol, params = M.gauss_jordan_solve(b)
    if params.shape[0] > 0:
        sol = sol.subs({p_: 0 for p_ in params})
    resid = (M * sol - b).applyfunc(sp.simplify)
    return [sp.nsimplify(x) for x in sol], all(r == 0 for r in resid)

def q2(z):
    return Q2.of(z)

def sector(name, t1, w1, t2, w2, structs, omega=2, nsamp=4):
    rows = {}; rhs = {}
    for it in range(nsamp):
        p = rng.integers(-2, 3, size=4) * 2
        if not p.any(): p[0] = 2
        a1 = rpol(t1, w1); a2 = rpol(t2, w2)
        ext = [(t1, w1, a1, p), (t2, w2, a2, -p)]
        G = GraphSum(ext)
        res = G.coefficients(omega)
        S = [Q2.of(sv) for sv in structs(a1, a2, p)]
        for key, v in res.items():
            for part in ('log', 'epsinv'):
                rows.setdefault((key, part), []).append(S)
                rhs.setdefault((key, part), []).append(v[part])
    out = {}
    tot = {}
    for (key, part), rr in rows.items():
        sol, ok = solve_exact(rr, rhs[(key, part)])
        out[(str(key), part)] = ([str(x / C.C_ADJ) for x in sol], ok)
        tot.setdefault(part, [0] * len(sol))
        tot[part] = [a + x / C.C_ADJ for a, x in zip(tot[part], sol)]
    return out, {k: [str(x) for x in v] for k, v in tot.items()}

if __name__ == '__main__' and len(sys.argv) == 1:
    def sAV(a1, a2, p):
        dV = np.einsum('ma,ma->', a1[:4], a2[:4])
        return [p @ p * dV, np.einsum('m,ma->a', p, a1[:4]) @ np.einsum('m,ma->a', p, a2[:4]), dV]
    out, tot = sector('AV-AV', 'A', 'V', 'A', 'V', sAV)
    for k, v in out.items(): print(k, v)
    print('TOTAL AV-AV [p^2(a1.a2), (p.a1)(p.a2), (a1.a2)] in units u:', tot)

def run_all():
    from sector_a import Bcross
    import check_vertices as cv
    cv.d = 3; cv.SH = {'c': (3,), 'A': (10, 3), 'l': (16, 3), 'As': (10, 3), 'ls': (16, 3), 'cs': (3,)}
    SECT = {'AV-AV': (('A', 'V'), ('A', 'V'), ['B1', 'B2', 'B3']),
            'AP-AP': (('A', 'P'), ('A', 'P'), ['B4', 'B5']),
            'AV-AP': (('A', 'V'), ('A', 'P'), ['B1', 'B4']),
            'l-l': (('l', None), ('l', None), ['B6']),
            'As-c': (('As', None), ('c', None), ['B7'])}
    summary = {}
    for nm, ((t1, w1), (t2, w2), Bn) in SECT.items():
        def structs(a1, a2, p, Bn=Bn, t1=t1, t2=t2):
            b1 = (t1, a1, p.astype(float)); b2 = (t2, a2, (-p).astype(float))
            return [Bcross(b, b1, b2, p.astype(float)) for b in Bn]
        out, tot = sector(nm, t1, w1, t2, w2, structs, nsamp=len(Bn) + 2)
        summary[nm] = {'B': Bn, 'total_units_u': tot,
                       'graphs': {k[0] + '|' + k[1]: v for k, v in out.items()}}
        print(nm, Bn, tot, flush=True)
    return summary

if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'all':
    S = run_all()
    json.dump(S, open('logs/r2_regress_quadratic.json', 'w'), indent=1)
