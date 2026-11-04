"""Exact linear algebra for the completion fits: values are dicts {log-basis key: (re, im) Fractions} ('1', 'log2',
'log3', 'log5', 'log7'); candidate values are Gaussian rationals.  Classes in H^{-1}(Q) on a block: closed subspace,
exact image, coordinates."""
from fractions import Fraction as Fr
import sympy as sp
LOGK = ['1', 'log2', 'log3', 'log5', 'log7']


def to_fr(z, bits=24):
    S = 2 ** bits
    re, im = complex(z).real * S, complex(z).imag * S
    ire, iim = round(re), round(im)
    if abs(re - ire) > 1e-5 or abs(im - iim) > 1e-5:
        raise ValueError(f'non-dyadic value {z!r}')
    return (Fr(ire, S), Fr(iim, S))


def gq(z):
    return sp.Rational(z[0].numerator, z[0].denominator) + sp.I * sp.Rational(z[1].numerator, z[1].denominator)


def clog_parse(d):
    """CLog json {'re': {k: str}, 'im': {k: str}} -> {k: (Fr, Fr)}"""
    out = {}
    for k in LOGK:
        out[k] = (Fr(d['re'].get(k, '0')), Fr(d['im'].get(k, '0')))
    return out


def lvec_add(a, b, s=1):
    return {k: (a[k][0] + s * b[k][0], a[k][1] + s * b[k][1]) for k in LOGK}


def lvec_scale(a, z):
    """multiply log-vector by Gaussian rational z=(re,im)"""
    x, y = z
    return {k: (a[k][0] * x - a[k][1] * y, a[k][0] * y + a[k][1] * x) for k in LOGK}


def lvec_zero():
    return {k: (Fr(0), Fr(0)) for k in LOGK}


def rat(z):
    return {k: ((z[0], z[1]) if k == '1' else (Fr(0), Fr(0))) for k in LOGK}


PMOD = 1000000009
IMOD = pow(13, (PMOD - 1) // 4, PMOD)
assert IMOD * IMOD % PMOD == PMOD - 1


def _modz(z):
    def fm(x):
        x = Fr(x)
        return x.numerator % PMOD * pow(x.denominator % PMOD, PMOD - 2, PMOD) % PMOD
    return (fm(z[0]) + IMOD * fm(z[1])) % PMOD


def independent_rows(rows):
    """greedy independent subset of rows (lists of Gaussian rational pairs): incremental elimination over
    Z/p with i -> sqrt(-1) mod p (p = 10^9 + 9); the exact solves downstream re-verify every fit."""
    keep, piv = [], []          # piv: list of (pivot column, reduced row)
    for i, r in enumerate(rows):
        v = [_modz(z) for z in r]
        for (c, pr) in piv:
            if v[c]:
                f = v[c]
                v = [(x - f * y) % PMOD for x, y in zip(v, pr)]
        nz = next((j for j, x in enumerate(v) if x), None)
        if nz is not None:
            inv = pow(v[nz], PMOD - 2, PMOD)
            v = [x * inv % PMOD for x in v]
            piv.append((nz, v)); keep.append(i)
    return keep


def solve_coords(basis_rows, target):
    """target = sum_k x_k basis_rows[k] exactly (None if not in span): pivot rows chosen modulo p, square exact solve
    with sympy, then exact verification on every row."""
    nb = len(basis_rows)
    if nb == 0:
        return [] if all(z == (0, 0) for z in target) else None
    ncol = len(target)
    cols = list(zip(*basis_rows))          # ncol rows of length nb
    sel = independent_rows([list(c) for c in cols])
    if len(sel) < nb:
        return _solve_coords_sympy(basis_rows, target)
    xs = gauss_solve([[(Fr(z[0]), Fr(z[1])) for z in cols[i]] for i in sel], [(Fr(target[i][0]), Fr(target[i][1])) for i in sel])
    if xs is None:
        return None
    x = [sp.Rational(a.numerator, a.denominator) + sp.I * sp.Rational(b.numerator, b.denominator) for a, b in xs]
    for i in range(ncol):
        re = sum(c[0] * z[0] - c[1] * z[1] for c, z in zip(xs, cols[i]))
        im = sum(c[0] * z[1] + c[1] * z[0] for c, z in zip(xs, cols[i]))
        if (re, im) != (Fr(target[i][0]), Fr(target[i][1])):
            return None
    return x


def _solve_coords_sympy(basis_rows, target):
    A = sp.Matrix([[gq(z) for z in r] for r in basis_rows]).T
    b = sp.Matrix([gq(z) for z in target])
    try:
        sol, params = A.gauss_jordan_solve(b)
    except ValueError:
        return None
    if params.shape[0]:
        sol = sol.subs({p: 0 for p in params})
    return [sp.nsimplify(sp.expand(x)) for x in sol]


def lsolve_coords(basis_rows, ltarget):
    """ltarget: list over configs of log-vectors -> {k: coords or None}"""
    out = {}
    for k in LOGK:
        t = [lv[k] for lv in ltarget]
        if all(z == (0, 0) for z in t):
            out[k] = [sp.Integer(0)] * len(basis_rows)
            continue
        out[k] = solve_coords(basis_rows, t)
    return out


LSYM = {'1': sp.Integer(1), 'log2': sp.log(2), 'log3': sp.log(3), 'log5': sp.log(5), 'log7': sp.log(7)}


def combine(lc):
    """{k: coords} -> sympy vector with logs."""
    n = len(next(iter(lc.values())))
    return sp.Matrix([sum(lc[k][i] * LSYM[k] for k in LOGK if lc[k] is not None) for i in range(n)])


def lstr(x):
    return str(sp.nsimplify(sp.expand(x)))


def lin_solve(eqs, unknowns):
    """linear system -> list with one dict {unknown: value} (free unknowns map to themselves), or [] if inconsistent."""
    eqs = [sp.expand(e) for e in eqs if sp.expand(e) != 0]
    if not eqs:
        return [{u: u for u in unknowns}]
    s = sp.linsolve(eqs, list(unknowns))
    if not s:
        return []
    t = list(s)[0]
    return [{u: sp.nsimplify(sp.expand(v)) for u, v in zip(unknowns, t)}]



def _cm(a, b):
    return (a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0])


def _cinv(a):
    n = a[0] * a[0] + a[1] * a[1]
    return (a[0] / n, -a[1] / n)


def gauss_solve(A, b):
    """exact solve of the square system A x = b over Q(i) (pairs of Fractions)."""
    n = len(A)
    M = [list(A[i]) + [b[i]] for i in range(n)]
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != (0, 0)), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        iv = _cinv(M[c][c])
        M[c] = [_cm(iv, z) for z in M[c]]
        for i in range(n):
            if i != c and M[i][c] != (0, 0):
                f = M[i][c]
                M[i] = [(z[0] - w[0], z[1] - w[1]) for z, w in zip(M[i], [_cm(f, y) for y in M[c]])]
    return [M[i][n] for i in range(n)]
