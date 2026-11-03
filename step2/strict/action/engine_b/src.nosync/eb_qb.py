"""Engine B, action side: Q on two-leg local functionals, canonical densities and exact linear algebra.

A local functional (possibly with explicit position, degree <= 1) is given by densities omega_0 and omega_mu (the
functional int(x^mu omega_mu + omega_0)); a density is a function of plane-wave legs (species, numeric polarisation,
momentum, no momentum conservation imposed) returning the vertex value in the convention of vendor/vertices.py
(coefficient of the odd leg parameters in leg order; colour structure (c0,c1) stripped).  Canonical data (Def
c11-canonical): last leg momentum -sum p_j - r,
    alpha = omega_0 + sum_mu (-i) d/dt omega_mu(r = t e_mu)|_0,   beta_mu = omega_mu(r = 0),
the same rule as the injected-momentum evaluation of the cycles (eb_pos).
Q acts on densities by the field substitutions of Prop c1-Q (c -> A_V, A -> A*, A* -> c*, lam -> lam*; symbol
ce_core.q_matrix, Prop c1-Q), leg by leg, with a Koszul sign rule fixed by the gate G-sign (Q^2 = 0, Q S_2(x) = 0,
and the four images of Theorem c9-quadratic (iii)); see find_sign_rule()."""
import itertools, random, math, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)
np_err = None
from fractions import Fraction as Fr
import numpy as np
import sympy as sp
import eb_core as E
import ce_core as C
import eb_pos as P

V = E.V
OFF, DIM, PAR = E.OFF, V.DIM, V.PAR
QTARGET = {'c': 'A', 'A': 'As', 'As': 'cs', 'lam': 'lams'}
_KS = sp.symbols('qk0:4')
_qfun = sp.lambdify([_KS], C.q_matrix(list(_KS)), 'numpy')
R = [np.array(C.Rm[M].tolist(), dtype=complex) for M in range(10)]
RB = [np.array(C.Rbm[M].tolist(), dtype=complex) for M in range(10)]
MP = np.array(C.Mp.tolist(), dtype=complex)
MR = [MP @ R[M] for M in range(10)]


def qmat(k):
    return np.array(_qfun([complex(x) for x in k]), dtype=complex)


def qleg(leg):
    """Q-image of a plane-wave leg (field of species t -> species QTARGET[t])."""
    t = leg['sp']
    if t not in QTARGET:
        return None
    s = QTARGET[t]
    v = np.zeros(E.NT, complex)
    v[OFF[t]:OFF[t] + DIM[t]] = leg['u']
    w = qmat(leg['p']) @ v
    u = w[OFF[s]:OFF[s] + DIM[s]].copy()
    if not np.any(u != 0):
        return None
    return dict(sp=s, u=u, p=leg['p'])


SIGN_RULE = (1, 0, dict(c=1, A=1, As=1, lam=1, cs=1, lams=1))   # fixed by gate G-sign (find in eb_gsign.py)


def apply_Q(dens, rule=None):
    """density -> density of Q(functional): sum over legs j whose species is hit, sign
    (-1)^{alpha n_odd(<j) + beta n_odd(>j)} eta[species of leg j]."""
    rule = rule or SIGN_RULE

    def f(legs):
        a, b, eta = rule
        tot = 0j
        # the Q-image functional has legs whose preimage... we evaluate QF on 'legs' by replacing a leg chi_j of
        # species t by its Q-image, i.e. dens(... qleg(chi_j) ...)
        for j, leg in enumerate(legs):
            ql = qleg(leg)
            if ql is None:
                continue
            nb = sum(PAR[l['sp']] for l in legs[:j])
            na = sum(PAR[l['sp']] for l in legs[j + 1:])
            sg = (-1) ** (a * nb + b * na) * eta[leg['sp']]
            new = list(legs)
            new[j] = ql
            v = dens(new)
            if v != 0:
                tot += sg * v
        return tot
    return f


# ------------------------------------------------------------------ exact conversion and canonical data
def to_fr(z, bits=20):
    S = 2 ** bits
    re, im = z.real * S, z.imag * S
    ire, iim = round(re), round(im)
    if abs(re - ire) > 1e-6 or abs(im - iim) > 1e-6:
        raise ValueError(f'non-dyadic density value {z!r}')
    return (Fr(ire, S), Fr(iim, S))


TV = (0, 1, -1, 2, -2)
WD = P.lagrange_deriv0(TV)


def canon(fn, legs_at):
    """fn = (omega0, [omega_mu x4] or None).  Returns (alpha, [beta_mu]) as Gaussian rational pairs."""
    w0, wmu = fn
    zero = [Fr(0)] * 4
    ar, ai = Fr(0), Fr(0)
    if w0 is not None:
        a = to_fr(w0(legs_at(zero)))
        ar += a[0]; ai += a[1]
    betas = [(Fr(0), Fr(0))] * 4
    if wmu is not None:
        betas = []
        for mu in range(4):
            if wmu[mu] is None:
                betas.append((Fr(0), Fr(0)))
                continue
            for t, w in zip(TV, WD):
                r = [Fr(t) if m == mu else Fr(0) for m in range(4)]
                v = to_fr(wmu[mu](legs_at(r)))
                # -i * w * v
                ar += w * v[1]; ai += -w * v[0]
                if t == 0:
                    bz = v
            betas.append(bz)
    return (ar, ai), betas


def Qfn(fn):
    w0, wmu = fn
    return (apply_Q(w0) if w0 is not None else None,
            [apply_Q(w) if w is not None else None for w in wmu] if wmu is not None else None)


# ------------------------------------------------------------------ modular linear algebra (ranks) over Q(i)
PRIME = 2 ** 61 - 1 if False else 1000000009     # p = 1 mod 4
IMOD = pow(13, (1000000009 - 1) // 4, 1000000009)
assert IMOD * IMOD % 1000000009 == 1000000009 - 1


def gq_mod(z, p=1000000009):
    re, im = z
    def fm(x):
        return x.numerator % p * pow(x.denominator % p, p - 2, p) % p
    return (fm(Fr(re)) + IMOD * fm(Fr(im))) % p


def rank_mod(rows, p=1000000009):
    M = [list(r) for r in rows]
    rk = 0
    ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = None
        for i in range(rk, len(M)):
            if M[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], p - 2, p)
        M[rk] = [x * inv % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def nullspace_mod(rows, p=1000000009):
    """basis of {x : sum_i x_i rows[i] = 0} (left kernel) mod p."""
    n = len(rows)
    ncol = len(rows[0])
    M = [list(rows[i]) + [1 if j == i else 0 for j in range(n)] for i in range(n)]
    rk = 0
    for c in range(ncol):
        piv = None
        for i in range(rk, n):
            if M[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        inv = pow(M[rk][c], p - 2, p)
        M[rk] = [x * inv % p for x in M[rk]]
        for i in range(n):
            if i != rk and M[i][c] % p:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return [row[ncol:] for row in M[rk:]]


# ------------------------------------------------------------------ numpy vertex functions (vertices.vertex convention)
def _perm_sign(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def dens(coef, slots, L):
    """density of the term coef * L(slot legs) with colour (c0,c1) [symmetric], vertices.vertex convention."""
    n = len(slots)

    def f(legs):
        if len(legs) != n:
            return 0j
        tot = 0j
        for sigma in itertools.permutations(range(n)):
            if any(legs[sigma[s]]['sp'] != slots[s] for s in range(n)):
                continue
            odd = [sigma[s] for s in range(n) if PAR[slots[s]]]
            tot += coef * _perm_sign(odd) * L([legs[sigma[s]] for s in range(n)])
        return tot
    return f


def dsum(*terms):
    def f(legs):
        return sum((c * d(legs) for c, d in terms), 0j)
    return f


def pv(leg):
    return np.array([complex(x) for x in leg['p']])


def dA(leg, rng):
    """F_MN = i(p_M u_N - p_N u_M) with u restricted to rng, p only along R^4."""
    u = np.array([leg['u'][M] if M in rng else 0 for M in range(10)], complex)
    p = np.zeros(10, complex)
    p[:4] = pv(leg)
    return 1j * (np.outer(p, u) - np.outer(u, p))


# basic functionals of Theorem c9-quadratic (iii) and the regression
F_AAs_V = dens(1, ('A', 'As'), lambda s: np.dot(s[0]['u'][:4], s[1]['u'][:4]))
F_AAs_p = dens(1, ('A', 'As'), lambda s: np.dot(s[0]['u'][4:], s[1]['u'][4:]))
F_dAdA_V = dens(1, ('A', 'A'), lambda s: 0.5 * np.sum(dA(s[0], range(4)) * dA(s[1], range(4))))
F_dAdA_p = dens(1, ('A', 'A'), lambda s: 0.5 * np.sum(dA(s[0], range(4, 10)) * dA(s[1], range(4, 10))))
F_dcAs = dens(1, ('c', 'As'), lambda s: s[0]['u'][0] * np.dot(1j * pv(s[0]), s[1]['u'][:4]))
F_csc = dens(1, ('cs', 'c'), lambda s: s[0]['u'][0] * s[1]['u'][0])
F_lamslam = dens(1, ('lam', 'lams'), lambda s: s[0]['u'] @ MP @ s[1]['u'])
F_ldl = dens(1, ('lam', 'lam'), lambda s: sum((s[0]['u'] @ MR[m] @ s[1]['u']) * 1j * pv(s[1])[m] for m in range(4)))
F_AA = dens(1, ('A', 'A'), lambda s: np.dot(s[0]['u'], s[1]['u']))
F_dcA = dens(1, ('c', 'A'), lambda s: s[0]['u'][0] * np.dot(1j * pv(s[0]), s[1]['u'][:4]))


def rand_leg(spc, rng=None, lo=-2, hi=2):
    u = np.zeros(DIM[spc], complex)
    idx = range(DIM[spc]) if rng is None else rng
    for c in idx:
        u[c] = random.randint(lo, hi)
    return dict(sp=spc, u=u, p=[Fr(random.randint(-3, 3)) for _ in range(4)])


def check_identity(lhs, rhs, legsets, ntrial=4):
    """max |lhs - rhs| over random legs of the given species lists."""
    worst = 0.0
    for spcs in legsets:
        for _ in range(ntrial):
            legs = [rand_leg(s) for s in spcs]
            legs[-1]['p'] = [-sum(l['p'][m] for l in legs[:-1]) for m in range(4)]
            worst = max(worst, abs(lhs(legs) - rhs(legs)))
    return worst


def rhoF(F, xi):
    """rho(F) xi = (1/2) sum_MN F_MN rho_M rho_N xi (Sigma -> Sigma)."""
    v = np.zeros(16, complex)
    for M in range(10):
        for N in range(10):
            if F[M, N] != 0:
                v += 0.5 * F[M, N] * (RB[M] @ (R[N] @ xi))
    return v


def xi_terms(xi):
    """numpy densities of T0..T4 (eq at-S1list, Prop at-xi2) at the numeric xi in Sigma."""
    T1 = dens(-1, ('lam', 'As'), lambda s: sum((xi @ MR[M] @ s[0]['u']) * s[1]['u'][M] for M in range(4)))
    T2 = dens(-1, ('lam', 'As'), lambda s: sum((xi @ MR[M] @ s[0]['u']) * s[1]['u'][M] for M in range(4, 10)))
    T3 = dens(0.5, ('A', 'lams'), lambda s: rhoF(dA(s[0], range(4)), xi) @ MP @ s[1]['u'])
    T4 = dens(0.5, ('A', 'lams'), lambda s: rhoF(dA(s[0], range(4, 10)), xi) @ MP @ s[1]['u'])
    # T0 = int (d*A_V, (xi, lam*)), d*A = -d.A
    T0 = dens(1, ('A', 'lams'), lambda s: (-np.dot(1j * pv(s[0]), s[0]['u'][:4])) * (xi @ MP @ s[1]['u']))
    return dict(T0=T0, T1=T1, T2=T2, T3=T3, T4=T4)


def g_terms(xi):
    """numpy densities of g1..g4 (Prop at-xi2) at the numeric xi."""
    def gam(l, M):
        return xi @ MR[M] @ l['u']
    def g1(s):
        l, a = s; p = pv(a); p2 = p @ p; pa = p @ a['u'][:4]
        return sum(gam(l, n) * (p2 * a['u'][n] - p[n] * pa) for n in range(4))
    def g2(s):
        l, a = s; p = pv(a); pa = p @ a['u'][:4]
        return sum(gam(l, n) * p[n] * pa for n in range(4))
    def g3(s):
        l, a = s; p = pv(a); p2 = p @ p
        return sum(gam(l, M) * p2 * a['u'][M] for M in range(4, 10))
    def g4(s):
        c, ls = s; p = pv(c)
        return (p @ p) * c['u'][0] * (xi @ MP @ ls['u'])
    return dict(g1=dens(1, ('lam', 'A'), g1), g2=dens(1, ('lam', 'A'), g2), g3=dens(1, ('lam', 'A'), g3),
                g4=dens(1, ('c', 'lams'), g4))


def fit_on(target, basis, legsets, ntrial=6):
    """exact least-squares-free fit of a translation-invariant density onto basis densities on random conserving
    two-leg configurations; returns (coeffs as sympy, residual max)."""
    rows, rhs = [], []
    cfgs = []
    for spcs in legsets:
        for _ in range(ntrial):
            legs = [rand_leg(s) for s in spcs]
            legs[-1]['p'] = [-sum(l['p'][m] for l in legs[:-1]) for m in range(4)]
            cfgs.append(legs)
    for legs in cfgs:
        rows.append([to_fr(f(legs)) for f in basis])
        rhs.append(to_fr(target(legs)))
    M = sp.Matrix([[sp.Rational(str(a)) + sp.I * sp.Rational(str(b)) for a, b in r] for r in rows])
    bv = sp.Matrix([sp.Rational(str(a)) + sp.I * sp.Rational(str(b)) for a, b in rhs])
    sol, params = M.gauss_jordan_solve(bv)
    sol = sol.subs({p: 0 for p in params})
    return [sp.nsimplify(x) for x in sol]


def apply_vf(densf, legmap, rule=None):
    """odd derivation {G,-} acting from the left by the slot substitution legmap (leg -> image leg or None), with the
    Koszul rule of apply_Q (used for G = T1: lam-leg -> A_V leg with polarisation dl_A* T1 = +Gamma(xi,lam)|_V)."""
    rule = rule or SIGN_RULE

    def f(legs):
        a, b, eta = rule
        tot = 0j
        for j, leg in enumerate(legs):
            ql = legmap(leg)
            if ql is None:
                continue
            nb = sum(PAR[l['sp']] for l in legs[:j])
            na = sum(PAR[l['sp']] for l in legs[j + 1:])
            new = list(legs)
            new[j] = ql
            v = densf(new)
            if v != 0:
                tot += (-1) ** (a * nb + b * na) * v
        return tot
    return f


def map_T1(xi):
    def m(leg):
        if leg['sp'] != 'lam':
            return None
        u = np.zeros(10, complex)
        for M in range(4):
            u[M] = xi @ MR[M] @ leg['u']
        return dict(sp='A', u=u, p=leg['p'])
    return m


F_dsdc_dsA = dens(1, ('c', 'A'), lambda s: (pv(s[0]) @ pv(s[0])) * s[0]['u'][0] * (-(1j * pv(s[1])) @ s[1]['u'][:4]))
F_dsA_dsA = dens(1, ('A', 'A'), lambda s: (-(1j * pv(s[0])) @ s[0]['u'][:4]) * (-(1j * pv(s[1])) @ s[1]['u'][:4]))
