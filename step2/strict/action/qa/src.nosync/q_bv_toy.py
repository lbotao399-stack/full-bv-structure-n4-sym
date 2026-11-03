"""QA third method for the sign finding: a finite-dimensional BV algebra (even fields x1,x2 with odd antifields
s1,s2; odd field t with even antifield y), Delta = sum d_x d_s + d_t d_y (left derivatives), the bracket of
Def c1-bracket {f,g} = Delta(fg) - Delta(f) g - (-1)^{|f|} f Delta(g).  Checks on random polynomials:
 (a) graded symmetry {g,f} = (-1)^{|f||g|}{f,g};
 (b) derivation rule of L = {S,-}, S even:  L{a,b} = -{La,b} - (-1)^{|a|}{a,Lb}  (C rule)
     against the standard rule L{a,b} = {La,b} + (-1)^{|a|+1}{a,Lb} (used in the proof of Lemma at-sformulas);
 (c) the arity <= 2 component formulas of s on degree-1 cochains derived from Def c2-dce:
     (sG)^(1)(x) = L G1(x) + {G0,S(x)},
     (sG)^(2)(x,y) = L G2 - G1([x,y]) - {S(x),G1(y)} + (-1)^{|x|}{G1(x),S(y)} + {G0,S2(x,y)},
     tested by s(sF) = 0 for random degree-0 F, with a toy action of the two-dimensional Lie algebra [e1,e2] = e2
     by Hamiltonian vector fields (S = 0 background, S2 = 0), against the printed formulas of the Lemma."""
import random, itertools
from fractions import Fraction as Fr
EVEN = ['x1', 'x2', 'y']; ODD = ['s1', 's2', 't']
PAIRS = [('x1', 's1'), ('x2', 's2'), ('t', 'y')]
# element: dict {(evexp tuple, oddset tuple sorted): coeff}
def mono_mul(m1, m2):
    e = tuple(a + b for a, b in zip(m1[0], m2[0]))
    o1, o2 = list(m1[1]), list(m2[1])
    if set(o1) & set(o2):
        return None, 0
    seq = o1 + o2
    sgn = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                sgn = -sgn
    return (e, tuple(sorted(seq))), sgn
def mul(f, g):
    h = {}
    for m1, c1 in f.items():
        for m2, c2 in g.items():
            m, s = mono_mul(m1, m2)
            if s:
                h[m] = h.get(m, 0) + s * c1 * c2
    return {m: c for m, c in h.items() if c}
def add(*fs):
    h = {}
    for f in fs:
        for m, c in f.items():
            h[m] = h.get(m, 0) + c
    return {m: c for m, c in h.items() if c}
def scal(a, f):
    return {m: a * c for m, c in f.items() if a * c}
def par(f):
    ps = {len(m[1]) % 2 for m in f}
    assert len(ps) <= 1, 'inhomogeneous'
    return ps.pop() if ps else 0
def d_even(f, v):
    i = EVEN.index(v); h = {}
    for (e, o), c in f.items():
        if e[i]:
            e2 = list(e); e2[i] -= 1
            h[(tuple(e2), o)] = h.get((tuple(e2), o), 0) + c * e[i]
    return {m: c for m, c in h.items() if c}
def d_odd(f, v):  # left derivative
    h = {}
    for (e, o), c in f.items():
        if v in o:
            k = o.index(v); s = (-1) ** k
            o2 = tuple(x for x in o if x != v)
            h[(e, o2)] = h.get((e, o2), 0) + s * c
    return {m: c for m, c in h.items() if c}
def d(f, v):
    return d_even(f, v) if v in EVEN else d_odd(f, v)
def Delta(f):
    return add(*[d(d(f, b), a) for a, b in PAIRS])
def br(f, g):
    return add(Delta(mul(f, g)), scal(-1, mul(Delta(f), g)), scal(-(-1) ** par(f), mul(f, Delta(g))))
def var(v):
    if v in EVEN:
        e = [0, 0, 0]; e[EVEN.index(v)] = 1; return {(tuple(e), ()): 1}
    return {((0, 0, 0), (v,)): 1}
def rand_poly(parity, nterm=5):
    f = {}
    while len(f) < nterm:
        e = tuple(random.randint(0, 2) for _ in EVEN)
        k = random.randint(0, 3)
        o = tuple(sorted(random.sample(ODD, k)))
        if k % 2 != parity:
            continue
        f[(e, o)] = f.get((e, o), 0) + random.randint(-3, 3)
    return {m: c for m, c in f.items() if c} or rand_poly(parity, nterm)
def eq(f, g):
    return add(f, scal(-1, g)) == {}
random.seed(5)
okA = okC = True; okStd = True
for trial in range(40):
    a = rand_poly(random.randint(0, 1)); b = rand_poly(random.randint(0, 1)); S = rand_poly(0)
    okA &= eq(br(b, a), scal((-1) ** (par(a) * par(b)), br(a, b)))
    L = lambda f: br(S, f)
    lhs = L(br(a, b))
    okC &= eq(lhs, add(scal(-1, br(L(a), b)), scal(-(-1) ** par(a), br(a, L(b)))))
    okStd &= eq(lhs, add(br(L(a), b), scal((-1) ** (par(a) + 1), br(a, L(b)))))
print('(a) graded symmetry:', okA)
print('(b) C rule  L{a,b} = -{La,b} - (-1)^|a|{a,Lb}:', okC, '; standard rule:', okStd)
# (c) toy action: S = 0, S1(e) = odd Hamiltonian functionals; [e1,e2] = e2.  Choose S1(e1) = x2 s2, S1(e2) = s2
X1, X2, S1v, S2v = var('x1'), var('x2'), var('s1'), var('s2')
S1 = {1: mul(X2, S2v), 2: S2v}
br12 = {(1, 2): (2, 1), (2, 1): (2, -1), (1, 1): (None, 0), (2, 2): (None, 0)}
def bracket_lin(f, x, y):            # f([x,y]) for f linear on g
    z, c = br12[(x, y)]
    return scal(c, f(z)) if z else {}
# arity-two MC component with S = 0, S2 = 0: S1([x,y]) + {S1(x),S1(y)} = 0 (all |x| = 0)?  fix the sign of S1 if needed
mc = add(bracket_lin(lambda z: S1[z], 1, 2), br(S1[1], S1[2]))
if mc:
    S1 = {k: scal(-1, v) for k, v in S1.items()}
    mc = add(bracket_lin(lambda z: S1[z], 1, 2), br(S1[1], S1[2]))
print('toy MC arity two holds:', mc == {})
ok_mine = ok_printed = True
for trial in range(20):
    F0 = rand_poly(0); F1 = {1: rand_poly(1), 2: rand_poly(1)}
    # (sF) for degree 0 (Lemma at-sformulas (i), unaffected): L = 0 here
    G0 = {}
    G1 = {z: br(F0, S1[z]) for z in (1, 2)}
    def G2(x, y):
        return add(bracket_lin(lambda z: F1[z], x, y), br(S1[x], F1[y]), br(F1[x], S1[y]))
    for (x, y) in [(1, 2), (2, 1)]:
        mine = add(scal(-1, bracket_lin(lambda z: G1[z], x, y)), scal(-1, br(S1[x], G1[y])), br(G1[x], S1[y]))
        printed = add(scal(-1, bracket_lin(lambda z: G1[z], x, y)), scal(-1, br(G1[x], S1[y])), scal(-1, br(S1[x], G1[y])))
        ok_mine &= (mine == {}); ok_printed &= (printed == {})
print('(c) s(sF) = 0 at arity two with the Def c2-dce formula:', ok_mine, '; with the printed Lemma at-sformulas (ii):', ok_printed)
