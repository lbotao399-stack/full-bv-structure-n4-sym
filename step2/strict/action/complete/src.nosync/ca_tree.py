"""Route A tree level: Q on numpy densities with engine A's symbol (qmat_A, validated against the mode-model Q by
gate GA-Q), mode-model functionals for the named quantities (T's, g's, S^(2), translations) and their brackets."""
import numpy as np
import ca_common as A
import cb_common as B          # only eb_qb.dens / dsum combinatorics (realisation independent)
Q = B.Q
OFFS = {'A': 0, 'lam': 10, 'c': 26, 'As': 27, 'lams': 37, 'cs': 53}
DIMS = {'A': 10, 'lam': 16, 'c': 1, 'As': 10, 'lams': 16, 'cs': 1}
PAR = {'A': 0, 'lam': 1, 'c': 1, 'As': 1, 'lams': 0, 'cs': 0}
QT = {'c': 'A', 'A': 'As', 'As': 'cs', 'lam': 'lams'}


def qleg(leg):
    t = leg['sp']
    if t not in QT:
        return None
    s = QT[t]
    v = np.zeros(54, complex)
    v[OFFS[t]:OFFS[t] + DIMS[t]] = leg['u']
    w = A.qmat_A(leg['p']) @ v
    u = w[OFFS[s]:OFFS[s] + DIMS[s]].copy()
    if not np.any(u):
        return None
    return dict(sp=s, u=u, p=leg['p'])


def apply_Q(dens):
    def f(legs):
        tot = 0j
        for j, leg in enumerate(legs):
            ql = qleg(leg)
            if ql is None:
                continue
            nb = sum(PAR[l['sp']] for l in legs[:j])
            new = list(legs); new[j] = ql
            tot += (-1) ** nb * dens(new)
        return tot
    return f


# ------------------------------------------------------------------ mode-model functionals (route A, exact GA)
def gam(Q1, Q2):
    return np.array([Q1 @ A.SB[M] @ Q2 for M in range(10)])


def gam_star(w1, w2):
    return np.array([w1 @ A.SS[M].T @ w2 for M in range(10)])


def mm_S2(mm, Q1, Q2, parts=('ll', 'Ac')):
    G12 = gam(Q1, Q2)
    F = None
    if 'll' in parts:
        Gs = np.array([[gam_star(np.eye(16)[i], np.eye(16)[j]) @ G12 for j in range(16)] for i in range(16)])
        F = mm.quad('ls', 'ls', lambda i, j, p1, p2: 0.25 * Gs[i, j] - 0.5 * Q1[i] * Q2[j])
    if 'Ac' in parts:
        F2 = mm.quad('A', 'cs', lambda i, j, p1, p2: -G12[i])
        F = F2 if F is None else F + F2
    return F


def mm_S2xi(mm, xi, sect=('V', 'perp')):
    F = None
    for s in sect:
        rng = range(4) if s == 'V' else range(4, 10)
        X = A.MD.F_T1(mm, xi, rng) + A.MD.F_T3(mm, xi, rng)
        F = X if F is None else F + X
    return F


def mm_T(mm, xi, which):
    return {'T1': lambda: A.MD.F_T1(mm, xi, range(4)), 'T2': lambda: A.MD.F_T1(mm, xi, range(4, 10)),
            'T3': lambda: A.MD.F_T3(mm, xi, range(4)), 'T4': lambda: A.MD.F_T3(mm, xi, range(4, 10)),
            'T0': lambda: A.MD.F_T0(mm, xi)}[which]()


def mm_translation(mm, b, sectors=('V', 'perp', 'lam')):
    bd = lambda p1: 1j * complex(np.dot(b, p1[:4]))
    F = None
    def add(X):
        nonlocal F
        F = X if F is None else F + X
    if 'V' in sectors:
        add(mm.quad('A', 'As', lambda i, j, p1, p2: bd(p1) if (i == j and i < 4) else 0))
        add(mm.quad('c', 'cs', lambda i, j, p1, p2: -bd(p1)))
    if 'perp' in sectors:
        add(mm.quad('A', 'As', lambda i, j, p1, p2: bd(p1) if (i == j and i >= 4) else 0))
    if 'lam' in sectors:
        add(mm.quad('l', 'ls', lambda i, j, p1, p2: -bd(p1) if i == j else 0))
    return F


def mm_value(mm, F, L):
    v = mm.value(F, A.mm_legs(L))
    return (v.re, v.im)


def mm_deriv(mm, G, gpar):
    """F -> {G, F} for a quadratic mode-model functional G of parity gpar, as the left derivation fixed by the values
    {G, y_v} = mm.br(G, y_v) on the mode variables (Koszul sign (-1)^{n_odd before} for even G, none for odd G);
    validated against mm.br and mm.Q (gate GA-deriv)."""
    from grass import GA, PAR
    from r2_core import Q2
    cache = {}

    def gy(v):
        if v not in cache:
            cache[v] = mm.br(G, GA.var(v), gpar)
        return cache[v]

    def D(F):
        out = {}
        for key, c in F.t.items():
            nodd = 0
            for i, v in enumerate(key):
                gv = gy(v)
                if gv.t:
                    pre = GA({key[:i]: Q2(1)}); post = GA({key[i + 1:]: Q2(1)})
                    sg = (-1) ** nodd if gpar == 0 else 1
                    term = pre * gv * post * (c * sg)
                    for k2, c2 in term.t.items():
                        out[k2] = out[k2] + c2 if k2 in out else c2
                nodd += PAR[v]
        return GA({k: v for k, v in out.items() if not v.iszero()})
    return D
