"""
conventions.py -- combinatorial lock of the [ESW Def 1.26] sign conventions.

Compares, symbolically and for every parity pattern up to arity 4, the component form of
    d_CE S_g + 1/2 {S_g, S_g}                                   [ESW Def 1.26, 1.27(c)]
with the generating-function (GF) model on O(g[1]) (x) O_loc:
    ghosts c^a of parity |t_a|+1, d_CE c^c = -1/2 sum_ab (-1)^{|t_a| p_b} f^c_ab c^a c^b  (odd left derivation),
    {c^I A, c^J B} = (-1)^{|c^I|+|c^J|+|c^J||A|} c^I c^J {A,B}   (bracket of [C07 Lemma 4.3.2] with ghosts passive;
                                                                  this is what engine.bracket does),
    F = S + sum_n (1/n!) sum_a kappa(a) c^{a_1}...c^{a_n} S^(n)(t_{a_1},...,t_{a_n}).
Values S^(k)(x_1..x_k) are formal symbols with parity k + sum|x_i| (S_g of total degree 0 [ESW Def 1.27]),
antisymmetric in the [ESW Def 1.26] sense; {A,B} formal, bilinear, {B,A} = (-1)^{|A||B|}{A,B} [ESW Def 1.15].
Test algebra: gl(2|1) (elementary supermatrices), which has brackets of every parity type.

Readings of the symbols |f|, |g| in Def 1.26:  'map' = degree of the multilinear map (= k for S^(k)),
'tot' = degree in C^bullet(g, O_loc) (= 0 for S_g).
"""
import itertools
from fractions import Fraction as Fr
from collections import defaultdict

# ---------------- gl(2|1) ----------------
IDX = [0, 1, 2]
PI = {0: 0, 1: 0, 2: 1}
BASIS = [(i, j) for i in IDX for j in IDX]
NB = len(BASIS)
PAR = [(PI[i] + PI[j]) % 2 for (i, j) in BASIS]
BI = {b: n for n, b in enumerate(BASIS)}


def lie(a, b):
    """[e_ij, e_kl] = d_jk e_il - (-1)^{|e_ij||e_kl|} d_li e_kj ; returns dict basis->coef."""
    (i, j), (k, l) = BASIS[a], BASIS[b]
    out = defaultdict(Fr)
    if j == k:
        out[BI[(i, l)]] += 1
    if l == i:
        out[BI[(k, j)]] -= (-1) ** (PAR[a] * PAR[b])
    return {k_: v for k_, v in out.items() if v}


def sgn_perm(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]:
                s = -s
    return s


def koszul(args, perm):
    """Koszul sign (parities PAR) of reordering args into [args[perm[0]], args[perm[1]], ...]."""
    s = 0
    for i in range(len(perm)):
        for j in range(i + 1, len(perm)):
            if perm[i] > perm[j]:
                s += PAR[args[perm[i]]] * PAR[args[perm[j]]]
    return -1 if s % 2 else 1


# ---------------- formal values ----------------
def canon_val(k, args):
    """S^(k)(args) -> (sign, key) with args sorted, ESW antisymmetry (-1)^{|x||y|+1} per adjacent swap."""
    a = list(args)
    sign = 1
    n = len(a)
    for i in range(n):
        for j in range(n - 1 - i):
            if a[j] > a[j + 1]:
                sign *= (-1) ** (PAR[a[j]] * PAR[a[j + 1]] + 1)
                a[j], a[j + 1] = a[j + 1], a[j]
    for j in range(n - 1):
        if a[j] == a[j + 1] and PAR[a[j]] == 0:
            return 0, None
    return sign, ('v', k, tuple(a))


def val_parity(key):
    _, k, a = key
    return (k + sum(PAR[x] for x in a)) % 2


def canon_br(s1, k1, s2, k2):
    if k1 is None or k2 is None:
        return 0, None
    if k1 <= k2:
        return s1 * s2, ('b', k1, k2)
    return s1 * s2 * (-1) ** (val_parity(k1) * val_parity(k2)), ('b', k2, k1)


def add(D, key, c):
    if key is None or c == 0:
        return
    D[key] += c
    if D[key] == 0:
        del D[key]


def val_lin(k, args):
    """S^(k) evaluated on args where args entries are dicts (linear combos) -> dict key->coef."""
    out = defaultdict(Fr)
    lists = [list(a.items()) for a in args]
    for combo in itertools.product(*lists):
        c = Fr(1)
        idx = []
        for (b, cb) in combo:
            c *= cb
            idx.append(b)
        s, key = canon_val(k, idx)
        add(out, key, c * s)
    return out


# ---------------- ESW side ----------------
def esw_dCE(k, args, read_f):
    """(d_CE S^(k))(x_1..x_n), n = k+1 [ESW Def 1.26]."""
    n = len(args)
    fdeg = k if read_f.startswith('map') else 0
    if read_f.endswith('flip'):
        fdeg += 1
    out = defaultdict(Fr)
    for i in range(n):
        for j in range(i + 1, n):
            e = PAR[args[i]] * sum(PAR[args[p]] for p in range(i))
            e += PAR[args[j]] * sum(PAR[args[p]] for p in range(j) if p != i)
            e += (i + 1) + (j + 1) + fdeg
            sign = -1 if e % 2 else 1
            br = lie(args[i], args[j])
            rest = [{args[p]: Fr(1)} for p in range(n) if p not in (i, j)]
            for key, c in val_lin(k, [br] + rest).items():
                add(out, key, sign * c)
    return out


def esw_br(k, l, args, read_g):
    """{S^(k), S^(l)}(x_1..x_{k+l}) [ESW Def 1.26]."""
    n = len(args)
    gdeg = l if read_g.startswith('map') else 0
    out = defaultdict(Fr)
    for left in itertools.combinations(range(n), k):
        right = [p for p in range(n) if p not in left]
        perm = list(left) + right
        e1 = gdeg * k + sum(PAR[args[p]] for p in left) * (l + gdeg)
        sign = sgn_perm(perm) * koszul(args, perm) * (-1 if e1 % 2 else 1)
        s1, k1 = canon_val(k, [args[p] for p in left])
        s2, k2 = canon_val(l, [args[p] for p in right])
        s, key = canon_br(s1, k1, s2, k2)
        add(out, key, sign * s)
    return out


def esw_MC(args, read_f, read_g):
    n = len(args)
    out = defaultdict(Fr)
    for key, c in esw_dCE(n - 1, args, read_f).items():
        add(out, key, c)
    for k in range(0, n + 1):
        for key, c in esw_br(k, n - k, args, read_g).items():
            add(out, key, Fr(1, 2) * c)
    return out


# ---------------- GF side ----------------
GP = [(PAR[a] + 1) % 2 for a in range(NB)]   # ghost parities


def gmono_canon(seq):
    """ghost monomial c^{seq[0]} c^{seq[1]}... -> (sign, sorted tuple) or (0, None)."""
    a = list(seq)
    sign = 1
    for i in range(len(a)):
        for j in range(len(a) - 1 - i):
            if a[j] > a[j + 1]:
                sign *= (-1) ** (GP[a[j]] * GP[a[j + 1]])
                a[j], a[j + 1] = a[j + 1], a[j]
    for j in range(len(a) - 1):
        if a[j] == a[j + 1] and GP[a[j]] == 1:
            return 0, None
    return sign, tuple(a)


def kappa(seq, n, pi_mode, s_n):
    """dictionary sign: (-1)^{sum_i p_{a_i} (pi + sum_{j<i} |t_{a_j}|)} * s_n."""
    pi = 0 if pi_mode == 0 else n
    e = 0
    for i, a in enumerate(seq):
        e += GP[a] * (pi + sum(PAR[seq[j]] for j in range(i)))
    return s_n * (-1 if e % 2 else 1)


def gf_F(nmax, pi_mode, svec):
    """F_n terms: dict (ghost tuple sorted, value key) -> coef, for n = 1..nmax. (S itself handled separately.)"""
    F = defaultdict(Fr)
    from math import factorial
    for n in range(1, nmax + 1):
        for seq in itertools.product(range(NB), repeat=n):
            sg, gm = gmono_canon(seq)
            if sg == 0:
                continue
            sv, key = canon_val(n, seq)
            if sv == 0:
                continue
            c = Fr(kappa(seq, n, pi_mode, svec[n]) * sg * sv, factorial(n))
            add(F, (gm, key), c)
    return F


def gf_MC(nmax, pi_mode, svec):
    """coefficient of each sorted ghost monomial of length n in d_CE F + 1/2 {F,F}, F = S + sum F_n."""
    F = gf_F(nmax, pi_mode, svec)
    Sterm = {((), ('v', 0, ())): Fr(1)}
    Fall = defaultdict(Fr)
    for k, v in F.items():
        Fall[k] += v
    for k, v in Sterm.items():
        Fall[k] += v
    out = defaultdict(Fr)
    # d_CE on ghosts
    for (gm, key), c in F.items():
        acc = 0
        for i, a in enumerate(gm):
            sgn_pass = (-1) ** (sum(GP[gm[j]] for j in range(i)) % 2)
            # d c^a = -1/2 sum_{b,d} (-1)^{|t_b| p_d} f^a_{bd} c^b c^d
            for b in range(NB):
                for d in range(NB):
                    f = lie(b, d).get(a, 0)
                    if f == 0:
                        continue
                    coef = Fr(-1, 2) * f * (-1) ** (PAR[b] * GP[d]) * sgn_pass
                    seq = list(gm[:i]) + [b, d] + list(gm[i + 1:])
                    sg, g2 = gmono_canon(seq)
                    if sg == 0:
                        continue
                    add(out, (g2, key), c * coef * sg)
    # 1/2 {F,F}
    items = list(Fall.items())
    for (g1, k1), c1 in items:
        for (g2, k2), c2 in items:
            if len(g1) + len(g2) > nmax + 1:
                continue
            pA = val_parity(k1)
            pc1 = sum(GP[x] for x in g1) % 2
            pc2 = sum(GP[x] for x in g2) % 2
            sign = (-1) ** ((pc1 + pc2 + pc2 * pA) % 2)
            sg, gm = gmono_canon(list(g1) + list(g2))
            if sg == 0:
                continue
            s, bk = canon_br(1, k1, 1, k2)
            add(out, (gm, bk), Fr(1, 2) * c1 * c2 * sign * sg * s)
    # group by ghost monomial
    res = defaultdict(lambda: defaultdict(Fr))
    for (gm, key), c in out.items():
        res[gm][key] += c
    return res


def compare(nmax, read_f, read_g, pi_mode, svec):
    gf = gf_MC(nmax, pi_mode, svec)
    ok = True
    ratios = {}
    for n in range(2, nmax + 1):
        for args in itertools.combinations_with_replacement(range(NB), n):
            sg, gm = gmono_canon(list(args))
            if sg == 0:
                continue
            e = {k: v for k, v in esw_MC(list(args), read_f, read_g).items() if v}
            g = {k: v for k, v in gf.get(gm, {}).items() if v}
            if not e and not g:
                continue
            if set(e) != set(g):
                ok = False
                break
            rs = set(g[k] / e[k] for k in e)
            if len(rs) != 1:
                ok = False
                break
            ratios[args] = rs.pop()
        if not ok:
            break
    return ok, ratios


def run(nmax=3):
    results = {}
    for read_f in ('map', 'tot', 'map_flip', 'tot_flip'):
        for read_g in ('map', 'tot'):
            found = []
            for pi_mode in (0, 1):
              for s1 in (1, -1):
                for s2 in (1, -1):
                    for s3 in (1, -1):
                        svec = {1: s1, 2: s2, 3: s3, 4: 1}
                        ok, ratios = compare(nmax, read_f, read_g, pi_mode, svec)
                        if ok:
                            found.append({'pi_mode': pi_mode, 's1': s1, 's2': s2, 's3': s3,
                                          'ratio_set': sorted(set(str(r) for r in ratios.values()))})
            results['|f|=%s,|g|=%s' % (read_f, read_g)] = found
    return results


if __name__ == '__main__':
    import json
    print(json.dumps(run(3), indent=1))


# ---------------------------------------------------------------------------------------------
# explicit component equations for abstract arguments (adopted convention: d_CE flipped, |f|=|g|=0)
# ---------------------------------------------------------------------------------------------
def component_form(pars, s_nonzero=lambda k, ps: True):
    """pars: tuple of parities of x_1..x_n.  Returns {term string: coefficient} of
    d_CE S_g + 1/2 {S_g,S_g} at (x_1..x_n); brackets [x_i,x_j] kept symbolic."""
    n = len(pars)
    names = ['x%d' % (i + 1) for i in range(n)]
    P = dict(zip(names, pars))

    def canon(k, args):
        a = list(args)
        sign = 1
        for i in range(len(a)):
            for j in range(len(a) - 1 - i):
                if a[j] > a[j + 1]:
                    sign *= (-1) ** (P[a[j]] * P[a[j + 1]] + 1)
                    a[j], a[j + 1] = a[j + 1], a[j]
        return sign, 'S%d(%s)' % (k, ','.join(a)) if k else 'S'

    def vpar(k, args):
        return (k + sum(P[a] for a in args)) % 2

    out = defaultdict(Fr)
    # d_CE S^(n-1) with extra overall sign
    for i in range(n):
        for j in range(i + 1, n):
            e = pars[i] * sum(pars[p] for p in range(i)) + pars[j] * sum(pars[p] for p in range(j) if p != i) + (i + 1) + (j + 1) + 1
            br = '[%s,%s]' % (names[i], names[j])
            P[br] = (pars[i] + pars[j]) % 2
            rest = [names[p] for p in range(n) if p not in (i, j)]
            s, key = canon(n - 1, [br] + rest)
            out[key] += (-1) ** (e % 2) * s
    # 1/2 sum_{k+l=n} {S^(k), S^(l)}
    for k in range(0, n + 1):
        l = n - k
        for left in itertools.combinations(range(n), k):
            right = [p for p in range(n) if p not in left]
            perm = list(left) + right
            e1 = sum(pars[p] for p in left) * l
            ks = 0
            for a in range(len(perm)):
                for b in range(a + 1, len(perm)):
                    if perm[a] > perm[b]:
                        ks += pars[perm[a]] * pars[perm[b]]
            sign = sgn_perm(perm) * (-1) ** (ks % 2) * (-1) ** (e1 % 2)
            la, ra = [names[p] for p in left], [names[p] for p in right]
            s1, A = canon(k, la)
            s2, B = canon(l, ra)
            pa, pb = vpar(k, la), vpar(l, ra)
            if (k, A) > (l, B):
                A, B = B, A
                sign *= (-1) ** (pa * pb)
            out['{%s,%s}' % (A, B)] += Fr(1, 2) * sign * s1 * s2
    return {k: str(v) for k, v in out.items() if v}


def component_forms():
    res = {}
    for pars in [(0, 0), (0, 1), (1, 0), (1, 1), (0, 1, 1), (1, 0, 1), (1, 1, 1), (0, 0, 1), (1, 1, 1, 1)]:
        res[''.join('e' if p == 0 else 'o' for p in pars)] = component_form(pars)
    return res
