# Audit #2: finite-dimensional check of the sign identities behind design Prop prop:o-laws (a), (d) and the twisted law.
# Graded-commutative polynomial algebra; BV operator Delta of degree +1, Delta^2 = 0; bracket exactly as ch01
# def:c1-bracket: {f,g} = Delta(fg) - Delta(f) g - (-1)^{|f|} f Delta(g).
# Checks (exact rationals, random homogeneous elements):
#  B1  {g,f} = (-1)^{|f||g|} {f,g}                                  (ch01 graded symmetry)
#  B2  for Q a degree-1 derivation of the product with [Q,Delta]=0:  Q{f,g} = -{Qf,g} - (-1)^{|f|}{f,Qg}
#  B3  for even I of degree 0:  {I,{f,g}} = -{{I,f},g} - (-1)^{|f|}{f,{I,g}}
#  B4  s = Q + {I,-}, phi odd (degree -1), O even (degree 0):  s{phi,O} - {phi,sO} = -{s phi, O}
#      i.e. [s, ad_phi] = -ad_{s phi}: law (a) "L_1 -> L_1 + [s,ad_phi]" equals law (d) "L_1 -> L_1 - {G,-}" at G = s phi.
#  B5  for odd phi, odd X: tree-level sign pattern used in the proof of (a) (lem c6-bvexact applied to X = sO) is the
#      same identity B4 with O replaced by sO (checked as s{phi,X} = -{s phi,X} - {phi,sX}).
#  B6  twisted coupling law: with formal symbols, r_xi(S - S2) = S + S2/2 = S^xi - S1 + S2 and N_R S^xi = -S1 + S2.
import random, json
from fractions import Fraction as Fr
random.seed(7)
# variables: (name, degree); parity = degree mod 2. Pairs (x, x*) with deg x* = -1 - deg x.
VARS = [('a1', 0), ('a1*', -1), ('a2', 0), ('a2*', -1), ('c', 1), ('c*', -2), ('l', 0), ('l*', -1)]
# give 'l' odd intrinsic parity (like lambda): parity = (degree + intrinsic) mod 2
INTR = {'l': 1, 'l*': 1}
NV = len(VARS)
par = [ (VARS[i][1] + INTR.get(VARS[i][0], 0)) % 2 for i in range(NV)]
deg = [VARS[i][1] for i in range(NV)]
PAIRS = [(0, 1), (2, 3), (4, 5), (6, 7)]
ODD = [i for i in range(NV) if par[i] == 1]; EVEN = [i for i in range(NV) if par[i] == 0]

def canon(odd_list):
    lst = list(odd_list); sign = 1
    for i in range(len(lst)):
        for j in range(len(lst)-1-i):
            if lst[j] > lst[j+1]:
                lst[j], lst[j+1] = lst[j+1], lst[j]; sign = -sign
    for i in range(len(lst)-1):
        if lst[i] == lst[i+1]: return 0, None
    return sign, tuple(lst)
def mono_mul(m1, m2):
    e1, o1 = m1; e2, o2 = m2
    s, o = canon(o1+o2)
    if s == 0: return 0, None
    return s, (tuple(a+b for a, b in zip(e1, e2)), o)
def mul(f, g):
    out = {}
    for m1, c1 in f.items():
        for m2, c2 in g.items():
            s, m = mono_mul(m1, m2)
            if s: out[m] = out.get(m, 0) + s*c1*c2
    return {m: c for m, c in out.items() if c != 0}
def add(*fs, coeffs=None):
    out = {}
    coeffs = coeffs or [1]*len(fs)
    for f, k in zip(fs, coeffs):
        for m, c in f.items(): out[m] = out.get(m, 0) + k*c
    return {m: c for m, c in out.items() if c != 0}
def scal(k, f): return {m: k*c for m, c in f.items() if k*c != 0}
def var(i):
    e = [0]*len(EVEN); o = ()
    if par[i] == 0: e[EVEN.index(i)] = 1
    else: o = (i,)
    return {(tuple(e), o): Fr(1)}
ONE = {(tuple([0]*len(EVEN)), ()): Fr(1)}
def dleft(i, f):
    out = {}
    for (e, o), c in f.items():
        if par[i] == 0:
            k = EVEN.index(i)
            if e[k] == 0: continue
            e2 = list(e); e2[k] -= 1
            m = (tuple(e2), o); out[m] = out.get(m, 0) + c*e[k]
        else:
            if i not in o: continue
            j = o.index(i); s = (-1)**j
            m = (e, o[:j]+o[j+1:]); out[m] = out.get(m, 0) + s*c
    return {m: c for m, c in out.items() if c != 0}
def Delta(f):
    # sum over pairs d_x d_{x*}; degree: -deg x - deg x* = +1
    return add(*[dleft(i, dleft(j, f)) for i, j in PAIRS])
def mdeg(m):
    e, o = m
    return sum(deg[EVEN[k]]*e[k] for k in range(len(EVEN))) + sum(deg[i] for i in o)
def degof(f):
    ds = {mdeg(m) for m in f}
    assert len(ds) <= 1, ds
    return ds.pop() if ds else 0
def pof(f):
    ps = {(mdeg(m) + sum(INTR.get(VARS[i][0], 0) for i in m[1]) + sum(INTR.get(VARS[EVEN[k]][0], 0)*m[0][k] for k in range(len(EVEN)))) % 2 for m in f}
    assert len(ps) <= 1
    return ps.pop() if ps else 0
def br(f, g):
    # parity used in the sign: |f| in ch01 is the degree; with intrinsic parity the Koszul sign uses total parity
    return add(Delta(mul(f, g)), scal(-1, mul(Delta(f), g)), scal(-(-1)**pof(f), mul(f, Delta(g))))
# random homogeneous polynomial of given degree, up to 3 factors
def monos_of_degree(d, maxlen=3, parity=None, minlen=1):
    res = []
    def rec(start, cur, cd):
        if len(cur) >= minlen and cd == d and (parity is None or (sum(par[i] for i in cur) % 2) == parity): res.append(list(cur))
        if len(cur) == maxlen: return
        for i in range(start, NV):
            if par[i] == 1 and i in cur: continue
            rec(i, cur+[i], cd+deg[i])
    rec(0, [], 0)
    return res
def rand_elem(d, nterms=6, parity=None, maxlen=3, minlen=1):
    if parity is None: parity = d % 2
    ms = monos_of_degree(d, maxlen, parity, minlen); f = {}
    for _ in range(nterms):
        lst = random.choice(ms); m = ONE
        for i in lst: m = mul(m, var(i))
        f = add(f, scal(Fr(random.randint(-5, 5)), m))
    return f
def eq(f, g): return add(f, scal(-1, g)) == {}

# Q: derivation induced by a quadratic S2 through the bracket, Q = {S2,-}; check [Q,Delta]=0 and Q^2 (not needed)
S2 = rand_elem(0, 6, 0, maxlen=2, minlen=2)
Qf = lambda f: br(S2, f)
I = rand_elem(0, 8, 0, maxlen=3, minlen=3)
res = {}
ok = True
for trial in range(6):
    f = rand_elem(random.choice([-1, 0, 1])); g = rand_elem(random.choice([-1, 0, 1]))
    ok &= eq(br(g, f), scal((-1)**(pof(f)*pof(g)), br(f, g)))
res['B1_graded_symmetry'] = ok
ok2 = True; ok3 = True
for trial in range(6):
    f = rand_elem(random.choice([-1, 0])); g = rand_elem(random.choice([-1, 0]))
    lhs = Qf(br(f, g)); rhs = add(scal(-1, br(Qf(f), g)), scal(-(-1)**pof(f), br(f, Qf(g))))
    ok2 &= eq(lhs, rhs)
    lhs = br(I, br(f, g)); rhs = add(scal(-1, br(br(I, f), g)), scal(-(-1)**pof(f), br(f, br(I, g))))
    ok3 &= eq(lhs, rhs)
res['B2_Q_derivation_rule'] = ok2
res['B3_I_derivation_rule'] = ok3
s = lambda f: add(Qf(f), br(I, f))
ok4 = True; ok5 = True; nonzero = False
for trial in range(6):
    phi = rand_elem(-1); O = rand_elem(0)
    lhs = add(s(br(phi, O)), scal(-1, br(phi, s(O))))
    rhs = scal(-1, br(s(phi), O))
    ok4 &= eq(lhs, rhs); nonzero |= (rhs != {})
    X = rand_elem(1)  # odd
    lhs = s(br(phi, X)); rhs = add(scal(-1, br(s(phi), X)), scal(-1, br(phi, s(X))))
    ok5 &= eq(lhs, rhs)
res['B4_[s,ad_phi]=-ad_(s phi)'] = ok4
res['B4_nontrivial'] = nonzero
res['B5_odd_X_rule'] = ok5
res['consequence'] = 'law (a) L1 -> L1 + [s,ad_phi] = L1 - {s phi,-}: identical to law (d) with G = s phi, a = 1'
# B6: formal check of the twisted coupling law
import sympy as sp
S, S1, S2s, c = sp.symbols('S S1 S2 c')
a = lambda k: sp.Integer(-1)**(k*(k-1)//2)/sp.factorial(k)
rG = S*a(0) + 0*a(1) + (-S2s)*a(2)
Sxi = S + S1 + a(2)*S2s
NR = {S: 0, S1: -S1, S2s: -2*S2s}   # N_R on S, S1(xi), S2(xi,xi): from arity-2 and arity-3 components with [R,xi]=xi
NR_Sxi = sp.expand(Sxi.subs(NR, simultaneous=True) - 0) if False else sp.expand(0 + (-S1) + a(2)*(-2*S2s))
res['B6'] = {'r_xi(S - S2)': str(sp.expand(rG)), 'S^xi - S1 + S2': str(sp.expand(Sxi - S1 + S2s)),
             'N_R S^xi': str(NR_Sxi), 'r_xi(G) == S^xi + N_R S^xi': bool(sp.expand(rG - (Sxi + NR_Sxi)) == 0)}
# a_k a_l = a_{k+l} binom(k+l,k) (-1)^{kl}
res['B7_decalage_coefficients'] = all(sp.simplify(a(k)*a(l) - a(k+l)*sp.binomial(k+l, k)*(-1)**(k*l)) == 0 for k in range(7) for l in range(7))
# MC family mu S + S1 + mu^-1 S2: arity-n terms and their coefficients
mu = sp.symbols('mu')
cs = {0: mu, 1: 1, 2: 1/mu}
ar2 = [cs[0]*cs[2], cs[1], cs[1]*cs[1]]      # {S,S2}, S1([.,.]), {S1,S1}
ar3 = [cs[2], cs[1]*cs[2]]                   # S2([.,.],.), {S1,S2}
res['B8_family_uniform'] = {'arity2': [str(x) for x in ar2], 'arity3': [str(sp.simplify(x)) for x in ar3],
                            'd/dmu at 1': 'S - S2'}
print(json.dumps(res, indent=1))
json.dump(res, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/design_audit_2/bv_sign_toy.json', 'w'), indent=1)
