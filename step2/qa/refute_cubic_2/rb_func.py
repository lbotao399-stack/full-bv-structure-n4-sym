"""Local functionals on plane waves (multilinear part, (2 pi)^4 delta removed), Q by substitution,
odd bracket from Delta_0 = sum d_x d_xi (derived from K_0 = -delta x inverse pairing, Ch.1 Def c1-bracket).
Everything written from Ch.1 (action, Q table) only.  Refuter #2, resumed run.
"""
import numpy as np
from rb_alg import (EPS, gbr, gadd, gmul, gscale, gcoef, canon, Amat, Bmat, p10, SLOTS, SHAPE, PARITY,
                    d_scalar, dstar_oneform, dstard, ddstar, dA2form, dslash)

RGEN = 0          # generator of the Q-parameter r (sorts first)
THETA = 999       # generator for internal odd leg in brackets


# ---------------------------------------------------------------- expression algebra
# an expression = list of terms (gens, vec, p, legs)
def proj(expr, part):
    out = []
    for g, v, p, L in expr:
        w = np.array(v, dtype=complex)
        if part == 'V':
            w[4:] = 0
        elif part == 'W':
            w[:4] = 0
        out.append((g, w, p, L))
    return out


def ex_map(expr, fn):
    return [(g, fn(v, p), p, L) for g, v, p, L in expr]


def pairs(e1, e2):
    for t1 in e1:
        for t2 in e2:
            if t1[3] & t2[3]:
                continue
            yield t1, t2


def triples(e1, e2, e3):
    for t1 in e1:
        for t2 in e2:
            if t1[3] & t2[3]:
                continue
            for t3 in e3:
                if (t1[3] | t2[3]) & t3[3]:
                    continue
                yield t1, t2, t3


def grass(gens, val):
    sg, key = canon(gens)
    if sg == 0 or abs(val) < 1e-14:
        return {}
    return {key: sg * val}


def total_legs(fields):
    L = frozenset()
    for s in SLOTS:
        for t in fields[s]:
            L = L | t[3]
    return L


def _full(Ls, fields):
    return Ls == total_legs(fields)


# ---------------------------------------------------------------- bilinear / trilinear integrals
def int_pair(e1, e2, fields, contract):
    """int (e1, e2) with contraction function contract(v1, v2) -> number; multilinear part."""
    out = {}
    allL = total_legs(fields)
    for t1, t2 in pairs(e1, e2):
        if (t1[3] | t2[3]) != allL:
            continue
        out = gadd(out, grass(t1[0] + t2[0], contract(t1[1], t2[1])))
    return out


def int_triple(e1, e2, e3, fields, contract):
    out = {}
    allL = total_legs(fields)
    for t1, t2, t3 in triples(e1, e2, e3):
        if (t1[3] | t2[3] | t3[3]) != allL:
            continue
        out = gadd(out, grass(t1[0] + t2[0] + t3[0], contract(t1[1], t2[1], t3[1])))
    return out


def int_quad(e1, e2, e3, e4, fields, contract):
    out = {}
    allL = total_legs(fields)
    for t1 in e1:
        for t2 in e2:
            if t1[3] & t2[3]:
                continue
            for t3 in e3:
                if (t1[3] | t2[3]) & t3[3]:
                    continue
                for t4 in e4:
                    if (t1[3] | t2[3] | t3[3]) & t4[3]:
                        continue
                    if (t1[3] | t2[3] | t3[3] | t4[3]) != allL:
                        continue
                    out = gadd(out, grass(t1[0] + t2[0] + t3[0] + t4[0],
                                          contract(t1[1], t2[1], t3[1], t4[1])))
    return out


# contractions
def dot(v1, v2):
    return np.sum(np.asarray(v1) * np.asarray(v2))


def two_form_pair(w1, w2):
    """(omega, eta) = 1/2 sum_{MN} omega_MN . eta_MN (Hodge normalisation of Def c1-action)."""
    return 0.5 * np.sum(w1 * w2)


# derived expressions
def dexpr(expr):           # d on scalar ghost -> one-form
    return ex_map(expr, d_scalar)


def dA(expr):             # 2-form
    return ex_map(expr, dA2form)


def wedge_br(e1, e2):
    """[X ^ Y]_{MN} = [X_M, Y_N] - [X_N, Y_M] for one-form expressions."""
    out = []
    for t1, t2 in pairs(e1, e2):
        X, Y = t1[1], t2[1]
        br = np.einsum('abc,ma,nb->mnc', EPS, X, Y)
        out.append((t1[0] + t2[0], br - np.transpose(br, (1, 0, 2)), t1[2] + t2[2], t1[3] | t2[3]))
    return out


def br_oneform_scalar(e1, e2):
    """[A, c]_M = [A_M, c]."""
    out = []
    for t1, t2 in pairs(e1, e2):
        out.append((t1[0] + t2[0], np.einsum('abc,ma,b->mc', EPS, t1[1], t2[1]), t1[2] + t2[2], t1[3] | t2[3]))
    return out


def br_scalar_oneform(e1, e2):
    """[c, A]_M = [c, A_M]."""
    out = []
    for t1, t2 in pairs(e1, e2):
        out.append((t1[0] + t2[0], np.einsum('abc,a,mb->mc', EPS, t1[1], t2[1]), t1[2] + t2[2], t1[3] | t2[3]))
    return out


def br_scalar_scalar(e1, e2):
    out = []
    for t1, t2 in pairs(e1, e2):
        out.append((t1[0] + t2[0], np.einsum('abc,a,b->c', EPS, t1[1], t2[1]), t1[2] + t2[2], t1[3] | t2[3]))
    return out


def br_scalar_spinor(e1, e2):
    """[c, lam]_s = [c, lam_s]."""
    out = []
    for t1, t2 in pairs(e1, e2):
        out.append((t1[0] + t2[0], np.einsum('abc,a,sb->sc', EPS, t1[1], t2[1]), t1[2] + t2[2], t1[3] | t2[3]))
    return out


def rho_br_A_lam(eA, elam):
    """rho([A, lam]) = sum_M A(e_M) [A_M, lam]  (Sigma^* valued)."""
    out = []
    for t1, t2 in pairs(eA, elam):
        v = np.einsum('Mst,abc,Ma,tb->sc', np.array([Amat(np.eye(10)[M]) for M in range(10)]), EPS, t1[1], t2[1])
        out.append((t1[0] + t2[0], v, t1[2] + t2[2], t1[3] | t2[3]))
    return out


# ---------------------------------------------------------------- the functionals
def F_free(f):
    A, lam, c, Ast = f['A'], f['lam'], f['c'], f['Ast']
    t1 = gscale(int_pair(dA(A), dA(A), f, two_form_pair), -0.25)
    t2 = gscale(int_pair(lam, ex_map(lam, dslash), f, dot), 0.5)
    t3 = int_pair(dexpr(c), Ast, f, dot)
    return gadd(t1, t2, t3)


def F_AAA(f):
    A = f['A']
    return gscale(int_pair(dA(A), wedge_br(A, A), f, two_form_pair), -0.25)


def F_lamlamA(f):
    lam, A = f['lam'], f['A']
    return gscale(int_pair(lam, rho_br_A_lam(A, lam), f, dot), 0.5)


def F_AcAst(f):
    return int_pair(br_oneform_scalar(f['A'], f['c']), f['Ast'], f, dot)


def F_clamlamst(f):
    return int_pair(br_scalar_spinor(f['c'], f['lam']), f['lamst'], f, dot)


def F_ccc(f):
    return gscale(int_pair(br_scalar_scalar(f['c'], f['c']), f['cst'], f, dot), 0.5)


def F_I3(f):
    return gadd(F_AAA(f), F_lamlamA(f), F_AcAst(f), F_clamlamst(f), F_ccc(f))


def F_I4(f):
    A = f['A']
    w = wedge_br(A, A)
    w = [(g, 0.5 * v, p, L) for g, v, p, L in w]
    return gscale(int_pair(w, w, f, two_form_pair), -0.25)


def F_S(f):
    return gadd(F_free(f), F_I3(f), F_I4(f))


AV = lambda f: proj(f['A'], 'V')
AW = lambda f: proj(f['A'], 'W')


def F_H1(f):
    return int_pair(dA(AV(f)), wedge_br(AV(f), AV(f)), f, two_form_pair)


def F_H2(f):
    # sum_{mu,I} (d_mu A_I, [A_mu, A_I]):  dA_W as (10,10,3) with only (mu,I) entries; [A_V,A_W]_{mu I}
    out = {}
    allL = total_legs(f)
    for t1 in AW(f):
        dW = np.zeros((10, 10, 3), dtype=complex)          # (mu, I, colour)
        dW[:4, 4:, :] = 1j * np.einsum('m,ia->mia', t1[2], t1[1][4:, :])
        for t2, t3 in pairs(AV(f), AW(f)):
            if t1[3] & (t2[3] | t3[3]):
                continue
            if (t1[3] | t2[3] | t3[3]) != allL:
                continue
            brVW = np.einsum('abc,ma,nb->mnc', EPS, t2[1], t3[1])     # [A_mu, A_I]
            val = np.sum(dW[:4, 4:, :] * brVW[:4, 4:, :])
            out = gadd(out, grass(t1[0] + t2[0] + t3[0], val))
    return out


def F_H3(f):
    return int_pair(br_oneform_scalar(AV(f), f['c']), proj(f['Ast'], 'V'), f, dot)


def F_H4(f):
    return int_pair(br_oneform_scalar(AW(f), f['c']), proj(f['Ast'], 'W'), f, dot)


def F_H5(f):
    return int_pair(br_scalar_spinor(f['c'], f['lam']), f['lamst'], f, dot)


def F_H6(f):
    return int_pair(br_scalar_scalar(f['c'], f['c']), f['cst'], f, dot)


def F_YV(f):
    return int_pair(f['lam'], rho_br_A_lam(AV(f), f['lam']), f, dot)


def F_YW(f):
    return int_pair(f['lam'], rho_br_A_lam(AW(f), f['lam']), f, dot)


def F_B1(f):
    return int_pair(dA(AV(f)), dA(AV(f)), f, two_form_pair)


def F_B2(f):
    return int_pair(dA(AW(f)), dA(AW(f)), f, two_form_pair)


def F_B3(f):
    return int_pair(f['lam'], ex_map(f['lam'], dslash), f, dot)


def F_B4(f):
    return int_pair(dexpr(f['c']), f['Ast'], f, dot)


def F_D1(f):
    return int_pair(ex_map(AV(f), dstar_oneform), ex_map(AV(f), dstar_oneform), f, dot)


def F_AA(f):
    return int_pair(f['A'], f['A'], f, dot)


def F_AH(f):
    return int_pair(f['A'], ex_map(f['A'], lambda v, p: np.dot(p, p) * v), f, dot)


def F_G1(f):
    return int_pair(br_scalar_oneform(f['c'], AV(f)), ex_map(AV(f), dstard), f, dot)


def F_G2(f):
    return int_pair(br_scalar_oneform(f['c'], AV(f)), ex_map(AV(f), ddstar), f, dot)


def F_G3(f):
    return int_pair(br_scalar_oneform(f['c'], AW(f)), ex_map(AW(f), dstard), f, dot)


def F_G4(f):
    return int_pair(br_scalar_oneform(f['c'], dexpr(f['c'])), f['Ast'], f, dot)


def F_G5(f):
    return int_pair(br_scalar_spinor(f['c'], f['lam']), ex_map(f['lam'], dslash), f, dot)


def F_E2(f):      # int (d^*d c, d^* A|_V)
    return int_pair(ex_map(f['c'], lambda v, p: np.dot(p, p) * v), ex_map(AV(f), dstar_oneform), f, dot)


# ---------------------------------------------------------------- configurations
def make_fields(legs):
    """legs: list of (slot, vec, p). Odd legs get generator j+1 (leg order)."""
    f = {s: [] for s in SLOTS}
    for j, (slot, vec, p) in enumerate(legs):
        g = (j + 1,) if PARITY[slot] == 1 else ()
        f[slot].append((g, np.array(vec, dtype=complex), np.array(p, dtype=float), frozenset([j])))
    return f


def leg_gens(legs):
    return tuple(j + 1 for j, (slot, v, p) in enumerate(legs) if PARITY[slot] == 1)


def value(F, legs):
    """plane-wave value: coefficient of the odd-leg generators in leg order."""
    return gcoef(F(make_fields(legs)), leg_gens(legs))


# ---------------------------------------------------------------- Q
def Q_fields(f):
    """fields + r Q(fields), r = generator RGEN placed on the left."""
    g = {s: list(f[s]) for s in SLOTS}
    for (gg, v, p, L) in f['c']:
        g['A'].append(((RGEN,) + gg, -d_scalar(v, p), p, L))            # A_V -> A_V - r dc
    for (gg, v, p, L) in f['A']:
        g['Ast'].append(((RGEN,) + gg, -0.5 * dstard(v, p), p, L))      # A* -> A* - (r/2) d*d A
    for (gg, v, p, L) in f['Ast']:
        g['cst'].append(((RGEN,) + gg, dstar_oneform(v, p), p, L))      # c* -> c* + r d^*(A*|_V)
    for (gg, v, p, L) in f['lam']:
        g['lamst'].append(((RGEN,) + gg, dslash(v, p), p, L))           # lam* -> lam* + r dslash lam
    return g


def Qvalue(F, legs):
    """plane-wave value of QF on legs: F(e + rQe) = F(e) + r QF(e)."""
    f = make_fields(legs)
    G = F(Q_fields(f))
    return gcoef(G, (RGEN,) + leg_gens(legs))


# ---------------------------------------------------------------- bracket
PAIRS = [('A', 'Ast'), ('lamst', 'lam'), ('cst', 'c')]   # (even member x, odd partner xi)


def _basis(slot):
    n = int(np.prod(SHAPE[slot]))
    for i in range(n):
        v = np.zeros(n, dtype=complex)
        v[i] = 1
        yield i, v.reshape(SHAPE[slot])


def _deriv(F, sub_legs, slot, k, vec):
    """left derivative of F along the internal mode 'vec e^{ikx}' in 'slot', with external legs sub_legs
    (a list of (j, slot, vec, p)); returns a Grassmann element in the leg generators."""
    f = {s: [] for s in SLOTS}
    for (j, sl, v, p) in sub_legs:
        gg = (j + 1,) if PARITY[sl] == 1 else ()
        f[sl].append((gg, np.array(v, dtype=complex), np.array(p, dtype=float), frozenset([j])))
    odd = PARITY[slot] == 1
    f[slot].append(((THETA,) if odd else (), np.array(vec, dtype=complex), np.array(k, dtype=float),
                    frozenset(['int'])))
    G = F(f)
    if not odd:
        return G
    out = {}
    for key, val in G.items():
        if THETA in key:
            rest = tuple(x for x in key if x != THETA)
            pos = key.index(THETA)             # sorted: THETA is last
            out = gadd(out, {rest: val * (-1) ** pos})
    return out


def bracket_value(F, parF, nF, Gf, nG, legs):
    """plane-wave value of {F, G} on 'legs'; F has nF legs, G has nG legs (as polynomials)."""
    import itertools
    m = len(legs)
    idx = list(range(m))
    total = {}
    for SF in itertools.combinations(idx, nF - 1):
        SG = [j for j in idx if j not in SF]
        if len(SG) != nG - 1:
            continue
        subF = [(j,) + tuple(legs[j]) for j in SF]
        subG = [(j,) + tuple(legs[j]) for j in SG]
        PF = sum((np.array(legs[j][2], dtype=float) for j in SF), np.zeros(4))
        k = -PF
        for (x, xi) in PAIRS:
            for i, e in _basis(x):
                # (d_xi F)(d_x G) + (-1)^{|F|} (d_x F)(d_xi G)
                a = gmul(_deriv(F, subF, xi, k, e), _deriv(Gf, subG, x, -k, e))
                b = gmul(_deriv(F, subF, x, k, e), _deriv(Gf, subG, xi, -k, e))
                total = gadd(total, a, gscale(b, (-1) ** parF))
    return gcoef(total, leg_gens(legs))
