"""Engine B, action side: candidate local functionals (plane-wave values in the convention of vendor/vertices.py, colour
structure stripped) and exact fitting of cycle sums (CLog values) onto them.

Candidate values are exact Gaussian rationals; a cycle value is a CLog (complex in Q + sum_p Q log p).  The fit solves,
for every log-basis element ('1', 'log2', 'log3', ...), the exact rational linear system over the configurations and
records the residual on the remaining configurations (must vanish)."""
from fractions import Fraction as Fr
import sympy as sp
import numpy as np
import eb_core as E
import ce_core as C

V = E.V
I = sp.I
Mp = C.Mp
Rm = C.Rm
Rbm = C.Rbm
MR = [Mp * Rm[M] for M in range(10)]
VOL = Rbm[0] * Rm[1] * Rbm[2] * Rm[3]          # rho(e0 e1 e2 e3) on Sigma


def sleg(leg):
    """numeric leg -> sympy leg for vertices.py."""
    u = [sp.nsimplify(complex(x).real) + I * sp.nsimplify(complex(x).imag) for x in leg['u']]
    p = [sp.Rational(x.numerator, x.denominator) if isinstance(x, Fr) else sp.nsimplify(x) for x in leg['p']]
    return dict(sp=leg['sp'], u=u, p=p)


def to_gq(z):
    z = sp.nsimplify(sp.expand(z))
    re, im = z.as_real_imag()
    return (Fr(str(sp.Rational(re))), Fr(str(sp.Rational(im))))


_REG = set()
_TAG = [0]


def newtag():
    _TAG[0] += 1
    return f't{_TAG[0]}'



def reg(name, coef, slots, fn, colour=('k', (0, 1))):
    if name not in _REG:
        V.TERMS[name] = (sp.nsimplify(coef), slots, fn, colour)
        _REG.add(name)
    return name


def value(terms, legs, key='k'):
    """plane-wave value (Gaussian rational pair) of sum(terms) on numeric legs, colour key component."""
    sl = [sleg(l) for l in legs]
    tot = V.total_vertex(terms, sl)
    return to_gq(tot.get(key, 0))


# ------------------------------------------------------------------ common quadratic functionals
def _dA(x, M, N):
    v = 0
    if M < 4:
        v += I * x['p'][M] * x['u'][N]
    if N < 4:
        v -= I * x['p'][N] * x['u'][M]
    return v


def L_dAdA(rng):
    def f(s):
        a, b = s
        aa = dict(u=[a['u'][M] if M in rng else 0 for M in range(10)], p=a['p'])
        bb = dict(u=[b['u'][M] if M in rng else 0 for M in range(10)], p=b['p'])
        return sum(_dA(aa, M, N) * _dA(bb, M, N) for M in range(10) for N in range(M + 1, 10))
    return f


def L_dsAdsA(s):
    a, b = s
    return (-sum(I * a['p'][m] * a['u'][m] for m in range(4))) * (-sum(I * b['p'][m] * b['u'][m] for m in range(4)))


def L_lamdlam(s):
    a, b = s
    return sum((sp.Matrix(a['u']).T * MR[m] * sp.Matrix(b['u']))[0] * I * b['p'][m] for m in range(4))


def L_lamdlam_odd(s):
    a, b = s
    vb = VOL * sp.Matrix(b['u'])
    return sum((sp.Matrix(a['u']).T * MR[m] * vb)[0] * I * b['p'][m] for m in range(4))


def L_dcAs(s):
    c, As = s
    return sum(I * c['p'][m] * c['u'][0] * As['u'][m] for m in range(4))


def L_AA(rng):
    def f(s):
        a, b = s
        return sum(a['u'][M] * b['u'][M] for M in rng)
    return f


BASIS2 = {}


def basis_delta2():
    """intrinsic two-leg structures of degree 0 at weight 4 (sector -> list of (label, terms))."""
    if BASIS2:
        return BASIS2
    BASIS2['AVAV'] = [('int(dA_V,dA_V)', [reg('F_dAdA_V', 1, ('A', 'A'), L_dAdA(range(4)))]),
                      ('int(d*A_V,d*A_V)', [reg('F_dsds', 1, ('A', 'A'), L_dsAdsA)])]
    BASIS2['ApAp'] = [('int(dA_perp,dA_perp)', [reg('F_dAdA_p', 1, ('A', 'A'), L_dAdA(range(4, 10)))])]
    BASIS2['lamlam'] = [('int(lam,dslash lam)', [reg('F_ldl', 1, ('lam', 'lam'), L_lamdlam)]),
                        ('int(lam,dslash rho(vol) lam) [parity odd]', [reg('F_ldl_odd', 1, ('lam', 'lam'), L_lamdlam_odd)])]
    BASIS2['cAsV'] = [('int(dc,A*)', [reg('F_dcAs', 1, ('c', 'As'), L_dcAs)])]
    BASIS2['AA_pole'] = [('int(A,A)', [reg('F_AA', 1, ('A', 'A'), L_AA(range(10)))])]
    return BASIS2


# ------------------------------------------------------------------ exact fitting
def gq_mul(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def solve_exact(A, b):
    """least-squares-free exact solve: A (rows x cols, Gaussian rational pairs as sympy), b list; returns (x, residual
    rows) using sympy over QQ_I; requires consistency (returns None if inconsistent)."""
    M = sp.Matrix([[sp.Rational(str(a[0])) + I * sp.Rational(str(a[1])) for a in row] for row in A])
    bv = sp.Matrix([sp.Rational(str(x[0])) + I * sp.Rational(str(x[1])) for x in b])
    try:
        sol, params = M.gauss_jordan_solve(bv)
    except ValueError:
        return None
    if params.shape[0]:
        sol = sol.subs({p: 0 for p in params})
    return [sp.nsimplify(sp.expand(x)) for x in sol]


def fit(values, cand_values, labels):
    """values: list of CLog (one per configuration); cand_values: list (per config) of list (per candidate) of Gaussian
    rational pairs.  Returns {log basis: {label: sympy number}} and the consistency flag per log basis element."""
    keys = set()
    for v in values:
        keys |= set(v.re.items().keys()) | set(v.im.items().keys())
    out = {}
    ok = True
    for k in sorted(keys):
        b = []
        for v in values:
            re = v.re.items().get(k, Fr(0))
            im = v.im.items().get(k, Fr(0))
            b.append((re, im))
        sol = solve_exact(cand_values, b)
        if sol is None:
            ok = False
            out[k] = None
        else:
            out[k] = {lab: x for lab, x in zip(labels, sol)}
    return out, ok
