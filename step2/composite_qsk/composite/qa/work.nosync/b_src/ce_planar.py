"""Planar multi-site reading of pair kernels (composite C^(2)_x and logarithmic K^log) for single traces of letters
with null derivatives, on collinear on-shell outputs p_j = t_j nbar (w = 0).

A pair kernel term c hbar/(16 pi^2) C_w on the adjacent pair (X_l X_{l+1}) of a single trace is replaced by
(-1)^r (c/2) Y_wr ... Y_w1 (units u); odd outputs: the multilinear form has the odd parameters in front in label order,
the trace order Y_wr ... Y_w1 costs the sign of that permutation of the odd outputs; identical outputs (same letter):
division by the symmetry factor 2.  Sites: (letter, number of n-derivatives); output (n.d)^c Y_j <-> (i t_j)^c."""
import math
import itertools
import sympy as sp
from sympy.polys.domains import QQ_I
import ce_fam as CF

I = sp.I


def perm_sign(seq):
    s = 1
    seq = list(seq)
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


def density(job, tab, letters_out, odd_out, which_m=0, avec=None, smax=0, logkey=None):
    """tab: {word: {m: {basis: poly}}} (composite) or, if logkey, {word: {'log': poly}} (twoletter job).
    letters_out: output letter names by label; odd_out: parity by label.
    Returns {(sites_out tuple, (a, b)): value} with sites_out = ((letter, c), ...) in trace order."""
    dens = {}
    nonrational = []
    nout = len(letters_out)
    ident = 2 if (nout == 2 and letters_out[0] == letters_out[1]) else 1
    for w, dm in tab.items():
        if logkey is not None:
            items = [('1', dm.get(logkey))] if dm.get(logkey) else []
        else:
            items = list(dm.get(which_m, {}).items())
        for basis, P in items:
            if avec is not None:
                P = CF.a_subst(job, P, avec)
            ent, wfree = CF.null_entries(job, P, nout, smax)
            for (cs, ins), v in ent.items():
                v = sp.nsimplify(v)
                if v == 0:
                    continue
                if basis != '1':
                    nonrational.append((w, basis, cs, ins, v))
                    continue
                order = list(reversed(w))                       # trace order of output labels
                sg = perm_sign([j for j in order if odd_out[j]])
                sites = tuple((letters_out[j], cs[j]) for j in order)
                key = (sites, ins)
                dens[key] = dens.get(key, 0) + v * sp.Rational((-1) ** len(w), 2) * sg / ident
    return {k: v for k, v in dens.items() if v != 0}, nonrational


# ------------------------------------------------------------------ graded cyclic words
def canon(word, odd):
    """word: tuple of sites (letter, nderiv); odd: dict letter -> parity.  Returns (canonical rotation, sign) with the
    graded cyclicity sign tr(X_1 ... X_L) = (-1)^{|X_L|(|X_1|+...+|X_{L-1}|)} tr(X_L X_1 ... X_{L-1})."""
    L = len(word)
    best = None
    w = list(word)
    sign = 1
    cands = []
    for r in range(L):
        cands.append((tuple(w), sign))
        last = w[-1]
        pl = odd[last[0]]
        prest = sum(odd[x[0]] for x in w[:-1])
        sign *= (-1) ** (pl * prest)
        w = [last] + w[:-1]
    m = min(c[0] for c in cands)
    sg = [s for c, s in cands if c == m]
    # a word equal to a rotation of itself with opposite sign vanishes
    if any(s != sg[0] for s in sg):
        return m, 0
    return m, sg[0]


def apply_pairs(word, dens, odd):
    """sum over the L cyclically adjacent pairs (rotated so that the pair is (X_l, X_{l+1})) of the planar density.
    word: tuple of sites (letter, nd); dens keyed by (sites_out, (a, b)) for letters given by dens_letters:
    dens is a dict (letterA, letterB) -> density.  Returns {canonical word: coefficient}."""
    L = len(word)
    out = {}
    w = list(word)
    sign = 1
    for r in range(L):
        # pair (w[0], w[1]) at the front of the rotated word
        X, Y = w[0], w[1]
        dd = dens.get((X[0], Y[0]), {})
        for (sites, ins), v in dd.items():
            if ins != (X[1], Y[1]):
                continue
            nw = tuple(sites) + tuple(w[2:])
            cw, cs = canon(nw, odd)
            if cs == 0:
                continue
            out[cw] = out.get(cw, 0) + v * sign * cs
        # rotate: move last to the front
        last = w[-1]
        sign *= (-1) ** (odd[last[0]] * sum(odd[x[0]] for x in w[:-1]))
        w = [last] + w[:-1]
    return {k: v for k, v in out.items() if v != 0}


def cyclic_basis(letters_multiset, L, odd, nderiv_total=0):
    """canonical graded-cyclic words of length L with the given letter content (dict letter -> count) and a total
    number of derivatives distributed over the sites; zero words removed."""
    lets = []
    for k, c in letters_multiset.items():
        lets += [k] * c
    words = set()
    for perm in set(itertools.permutations(lets)):
        for ds in itertools.product(range(nderiv_total + 1), repeat=L):
            if sum(ds) != nderiv_total:
                continue
            wd = tuple((perm[i], ds[i]) for i in range(L))
            cw, cs = canon(wd, odd)
            if cs != 0:
                words.add(cw)
    return sorted(words)


def matrix(basis_in, basis_out, fn):
    """fn(word) -> {canonical word: coeff}; matrix with columns = basis_in."""
    M = sp.zeros(len(basis_out), len(basis_in))
    idx = {w: i for i, w in enumerate(basis_out)}
    for j, w in enumerate(basis_in):
        for k, v in fn(w).items():
            if k not in idx:
                raise KeyError(f'output word {k} outside basis')
            M[idx[k], j] += v
    return M.applyfunc(sp.nsimplify)


def deriv_map(word, odd):
    """d_n on a word: sum over sites of (nd + 1) (no extra factors: (n.d) on the field)."""
    out = {}
    for i in range(len(word)):
        nw = list(word)
        nw[i] = (nw[i][0], nw[i][1] + 1)
        cw, cs = canon(tuple(nw), odd)
        if cs:
            out[cw] = out.get(cw, 0) + cs
    return out


def ibp_check(dD, dP, smax_in, nsites_out):
    """canonical-form identity [M^0, d_n] = M^n = 2 C_{P_n} (Lemma lem:cc-canonical) on a planar density:
    D0(., (a+1,b)) + D0(., (a,b+1)) = d_n(outputs of D0(., (a,b))) + 2 Pb(., (a,b)), for a + b < smax_in.
    Densities keyed by (sites_out, (a, b)), sites (letter, c).  Returns list of violations."""
    import itertools as _it
    keys_out = set(k[0] for k in list(dD) + list(dP))
    bad = []
    ins_all = set(k[1] for k in list(dD) + list(dP))
    for S in range(0, smax_in):
        for a in range(S + 1):
            b = S - a
            # candidate output site tuples: any appearing with one more derivative
            cands = set()
            for (sites, ins) in list(dD) + list(dP):
                if ins in ((a + 1, b), (a, b + 1), (a, b)):
                    cands.add(sites)
                    for i in range(len(sites)):
                        nw = list(sites); nw[i] = (nw[i][0], nw[i][1] + 1); cands.add(tuple(nw))
            for sites in cands:
                lhs = dD.get((sites, (a + 1, b)), 0) + dD.get((sites, (a, b + 1)), 0)
                rhs = 2 * dP.get((sites, (a, b)), 0)
                for i in range(len(sites)):
                    if sites[i][1] >= 1:
                        nw = list(sites); nw[i] = (nw[i][0], nw[i][1] - 1)
                        rhs += dD.get((tuple(nw), (a, b)), 0)
                if sp.simplify(lhs - rhs) != 0:
                    bad.append(((a, b), sites, str(lhs), str(rhs)))
    return bad
