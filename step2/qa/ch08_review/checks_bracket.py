#!/usr/bin/env python3
"""Finite-dimensional model of the bracket of [ch01, def:c1-bracket]:
{f,g} = Delta(fg) - Delta(f) g - (-1)^{|f|} f Delta(g),  Delta an odd second-order operator, Delta^2 = 0.

Checks (exact rational arithmetic):
 B1 graded symmetry {g,f} = (-1)^{|f||g|}{f,g}
 B2 ad_{{f,g}} = (-1)^{|f|+1} [ad_f, ad_g]   (graded commutator, parity of ad_f = |f|+1)
 B3 consequence used in Cor c8-homotopy (ii): with (-1)^{|x1|}{S1(x1),S1(x2)} = -S1([x1,x2]) - {S,S2(x1,x2)}
    the commutator of ad_{S1(x1)}, ad_{S1(x2)} equals -ad_{S1([x1,x2])} - ad_{{S,S2}} for all parities
 B4 [ad_S, ad_phi] = -ad_{{S,phi}} for even S and every phi (Prop c8-laws (a))
 B5 the ad_{{S,S2}} term is the graded commutator [ad_S, ad_{S2}] (exact term of Cor c8-homotopy)
Variables: even x_i and odd t_i paired by Delta = sum_i d/dx_i d/dt_i (left derivatives).
"""
import itertools, json, os, random
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
NV = 3          # number of pairs
random.seed(5)


# polynomial: dict {(xexp tuple, odd tuple sorted): Fraction}
def mono_mul(m1, m2):
    (e1, o1), (e2, o2) = m1, m2
    if set(o1) & set(o2):
        return None, 0
    seq = list(o1) + list(o2)
    # sign of sorting permutation
    sign = 1
    arr = seq[:]
    for i in range(len(arr)):
        for j in range(len(arr) - 1 - i):
            if arr[j] > arr[j + 1]:
                arr[j], arr[j + 1] = arr[j + 1], arr[j]
                sign = -sign
    e = tuple(a + b for a, b in zip(e1, e2))
    return (e, tuple(arr)), sign


def mul(f, g):
    out = {}
    for m1, c1 in f.items():
        for m2, c2 in g.items():
            m, s = mono_mul(m1, m2)
            if m is None:
                continue
            out[m] = out.get(m, 0) + s * c1 * c2
    return {k: v for k, v in out.items() if v != 0}


def add(*fs, coeffs=None):
    out = {}
    for i, f in enumerate(fs):
        c = 1 if coeffs is None else coeffs[i]
        for k, v in f.items():
            out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v != 0}


def scal(c, f):
    return {k: c * v for k, v in f.items() if c * v != 0}


def dx(f, i):
    out = {}
    for (e, o), c in f.items():
        if e[i] == 0:
            continue
        e2 = list(e); e2[i] -= 1
        k = (tuple(e2), o)
        out[k] = out.get(k, 0) + c * e[i]
    return {k: v for k, v in out.items() if v != 0}


def dt(f, i):
    out = {}
    for (e, o), c in f.items():
        if i not in o:
            continue
        r = o.index(i)
        o2 = tuple(a for a in o if a != i)
        k = (e, o2)
        out[k] = out.get(k, 0) + c * (-1) ** r
    return {k: v for k, v in out.items() if v != 0}


def Delta(f):
    return add(*[dx(dt(f, i), i) for i in range(NV)])


def parity(f):
    ps = {len(o) % 2 for (e, o) in f}
    assert len(ps) <= 1, 'inhomogeneous'
    return ps.pop() if ps else 0


def br(f, g):
    pf = parity(f)
    return add(Delta(mul(f, g)), mul(Delta(f), g), mul(f, Delta(g)), coeffs=[1, -1, -(-1) ** pf])


def rand_poly(par, deg=3, terms=4):
    f = {}
    while len(f) < terms:
        e = tuple(random.randint(0, 2) for _ in range(NV))
        k = random.choice([kk for kk in range(NV + 1) if kk % 2 == par])
        o = tuple(sorted(random.sample(range(NV), k)))
        f[(e, o)] = f.get((e, o), 0) + Fraction(random.randint(-3, 3) or 1)
    return {k: v for k, v in f.items() if v != 0}


def eq(f, g):
    return add(f, g, coeffs=[1, -1]) == {}


res = {}
# B0: Delta^2 = 0
ok = all(Delta(Delta(rand_poly(p))) == {} for p in (0, 1) for _ in range(5))
res['B0_Delta_squared_zero'] = ok
# B1, B2
ok1 = ok2 = True
alt_fail = 0
for pf, pg, ph in itertools.product((0, 1), repeat=3):
    for _ in range(3):
        f, g, hh = rand_poly(pf), rand_poly(pg), rand_poly(ph)
        ok1 &= eq(br(g, f), scal((-1) ** (pf * pg), br(f, g)))
        lhs = br(br(f, g), hh)
        comm = add(br(f, br(g, hh)), br(g, br(f, hh)), coeffs=[1, -(-1) ** ((pf + 1) * (pg + 1))])
        ok2 &= eq(lhs, scal((-1) ** (pf + 1), comm))
        # alternative without the sign (would be the naive Jacobi ad_{{f,g}} = [ad_f,ad_g])
        if not eq(lhs, comm):
            alt_fail += 1
res['B1_graded_symmetry'] = ok1
res['B2_ad_{f,g}=(-1)^{|f|+1}[ad_f,ad_g]'] = ok2
res['B2b_naive_ad_{f,g}=[ad_f,ad_g]_fails_in_cases'] = alt_fail
# B3: for every parity of x1, x2 (|S1(x)| = 1+|x|): [ad_{f},ad_{g}] with f=S1(x1), g=S1(x2) equals
# (-1)^{|x1|} ad_{{f,g}}; inserting the arity-two equation gives the displayed identity with no residual sign
ok3 = True
for px1, px2 in itertools.product((0, 1), repeat=2):
    pf, pg = (1 + px1) % 2, (1 + px2) % 2
    for _ in range(3):
        f, g, hh = rand_poly(pf), rand_poly(pg), rand_poly(random.randint(0, 1))
        comm = add(br(f, br(g, hh)), br(g, br(f, hh)), coeffs=[1, -(-1) ** ((pf + 1) * (pg + 1))])
        ok3 &= eq(comm, scal((-1) ** px1, br(br(f, g), hh)))
res['B3_[ad_S1(x1),ad_S1(x2)]=(-1)^{|x1|}ad_{S1,S1}_all_parities'] = ok3
# B4: [ad_S, ad_phi] = -ad_{{S,phi}} for even S, phi of either parity
ok4 = True
for pphi in (0, 1):
    for _ in range(3):
        S, phi, hh = rand_poly(0), rand_poly(pphi), rand_poly(random.randint(0, 1))
        comm = add(br(S, br(phi, hh)), br(phi, br(S, hh)), coeffs=[1, -(-1) ** ((0 + 1) * (pphi + 1))])
        ok4 &= eq(comm, scal(-1, br(br(S, phi), hh)))
res['B4_[ad_S,ad_phi]=-ad_{S,phi}_(S_even)'] = ok4
# B5: ad_{{S,S2}} = -[ad_S, ad_S2] for even S and even S2 (S2(x1,x2) on odd x's is even)
ok5 = True
for _ in range(3):
    S, S2, hh = rand_poly(0), rand_poly(0), rand_poly(random.randint(0, 1))
    comm = add(br(S, br(S2, hh)), br(S2, br(S, hh)), coeffs=[1, 1])   # graded commutator of two odd operators
    ok5 &= eq(br(br(S, S2), hh), scal(-1, comm))
res['B5_ad_{S,S2}=-[ad_S,ad_S2]'] = ok5
for k, v in res.items():
    print(k, v)
json.dump(res, open(os.path.join(HERE, 'checks_bracket.json'), 'w'), indent=1)
