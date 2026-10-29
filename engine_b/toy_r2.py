"""0-dim checks for round 2 (exact Grassmann algebra).
(1) Word expansion with stripped odd parameters: the coefficient of a product of distinct external symbols in
    W_1 = [hbar log(e^{hbar d_P} e^{I/hbar})]_{hbar^1} equals
    sum_{letter partitions} sum_{words} (1/(2n)) sign(odd reordering) STr[P Pi h_{w1} ... P Pi h_{wn}],
    h = c-number Hessian with the external odd symbols stripped to the front (in leg order inside a letter).
(2) Odd special edge: with delta odd and K an odd kernel (field-antifield type, K^{XY} = K^{YX}),
    d/d delta of W_1[P + delta K] (d_{P+delta K} = d_P + delta d_K) equals the STr formula (hbar/2) sum (1/n) STr[(Pi (P + delta K) H)^n]
    (delta in front; the parity operator acts on the ROW index of the kernel, which for K differs by -1 from acting on the row of H)."""
import itertools, random, sympy as sp
from grassmann import G

random.seed(17)
hb, g, h = sp.symbols('hbar g h')
xs = sp.symbols('x0 x1')
internal = [('e', xs[0]), ('e', xs[1]), ('o', 0), ('o', 1)]
par = [0, 0, 1, 1]
n = len(internal)
def var(v): return G.scalar(v[1]) if v[0] == 'e' else G.odd(v[1])
def D(F, v): return F.d_even(v[1]) if v[0] == 'e' else F.d_odd(v[1])
def R(): return sp.Integer(random.choice([-3, -2, -1, 1, 2, 3]))
P = [[sp.Integer(0)] * n for _ in range(n)]
for i in range(n):
    for j in range(i, n):
        if par[i] != par[j]: continue
        if par[i] == 0: P[i][j] = P[j][i] = R()
        elif i != j:
            P[i][j] = R(); P[j][i] = -P[i][j]
K = [[sp.Integer(0)] * n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if par[i] != par[j] and i < j:
            K[i][j] = K[j][i] = R()
X = [var(v) for v in internal]
# externals: legs 0 (even e0), 1 (odd, id 10), 2 (odd, id 11), 3 (even e1)
e0, e1 = sp.symbols('e0 e1')
EXT = [('e', e0), ('o', 10), ('o', 11), ('e', e1)]
def ext(k): return G.scalar(EXT[k][1]) if EXT[k][0] == 'e' else G.odd(EXT[k][1])
epar = [0, 1, 1, 0]
cub = {}; qua = {}
I3 = G(); I4 = G()
for k in range(4):
    for i in range(n):
        for j in range(i, n):
            if (par[i] + par[j] + epar[k]) % 2: continue
            if i == j and par[i] == 1: continue
            c = R(); I3 = I3 + ext(k) * X[i] * X[j] * c
for (k, l) in itertools.combinations(range(4), 2):
    for i in range(n):
        for j in range(i, n):
            if (par[i] + par[j] + epar[k] + epar[l]) % 2: continue
            if i == j and par[i] == 1: continue
            c = R(); I4 = I4 + ext(k) * ext(l) * X[i] * X[j] * c
I = I3 * g + I4 * h

def trunc(expr):
    p = sp.Poly(sp.expand(expr), g, h)
    out = 0
    for (a, b), cf in p.terms():
        if a <= 3 and b <= 1 and a + 2 * b <= 3:
            out += cf * g**a * h**b
    return sp.expand(out)
def T(F): return F.map(trunc)
def dK(F, Kmat):
    out = G()
    for i in range(n):
        for j in range(n):
            if Kmat[i][j] != 0:
                out = out + D(D(F, internal[i]), internal[j]) * (Kmat[i][j] / 2)
    return out
zero = {xs[0]: 0, xs[1]: 0}
def W1_exact(useK):
    Ex = G.scalar(1); term = G.scalar(1)
    for k in range(1, 5):
        term = T(term * I * (1 / hb) * sp.Rational(1, k)); Ex = Ex + term
    res = G(); cur = Ex
    for m in range(0, 5):
        res = res + cur * (hb**m / sp.factorial(m))
        nxt = T(dK(cur, P))
        if useK:
            nxt = nxt + G.odd(99) * T(dK(cur, K))
        cur = nxt
    res = T(res)
    Y = T(res - G.scalar(1))
    logv = Y + T(Y * Y) * sp.Rational(-1, 2) + T(T(Y * Y) * Y) * sp.Rational(1, 3) + T(T(T(Y * Y) * Y) * Y) * sp.Rational(-1, 4)
    W = T(logv * hb).subs(zero).drop_odd([0, 1])
    return W.map(lambda v: sp.expand(v).coeff(hb, 1))

def hess(F):
    return [[D(D(F, internal[k]), internal[i]).subs(zero).drop_odd([0, 1]) for k in range(n)] for i in range(n)]
def strip(Gel, odd_ids):
    """coefficient c-number of the product of the odd ids (in the given order) in front."""
    key = tuple(sorted(odd_ids))
    c = Gel.t.get(key, 0)
    # Gel = c * (sorted product); want Gel = c' * (product in given order): c' = c * sign(given -> sorted)
    s = 1; seq = list(odd_ids)
    for a in range(len(seq)):
        for b in range(a + 1, len(seq)):
            if seq[a] > seq[b]: s = -s
    return c * s
def letter(kind, legs):
    # c-number Hessian of the part of I linear in each listed external
    F = I3.map(lambda v: v) if kind == 'cubic' else I4
    H = hess(F)
    odd_ids = [EXT[k][1] for k in legs if EXT[k][0] == 'o']
    evens = [EXT[k][1] for k in legs if EXT[k][0] == 'e']
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            c = strip(H[i][k], odd_ids)
            c = sp.expand(c)
            for e in evens:
                c = sp.expand(c).coeff(e, 1)
            for e in (e0, e1):
                if e not in evens:
                    c = sp.expand(c).subs(e, 0)
            out[i][k] = c
    return out
def STr(mats):
    M = mats[0]
    for Xm in mats[1:]:
        M = [[sp.expand(sum(M[i][j] * Xm[j][k] for j in range(n))) for k in range(n)] for i in range(n)]
    return sp.expand(sum((-1) ** par[i] * M[i][i] for i in range(n)))
def PiH(Hm): return [[(-1) ** par[i] * Hm[i][k] for k in range(n)] for i in range(n)]
def word_formula(legs, order_g, order_h, Kmode=False):
    tot = 0
    parts = [[('cubic', (k,)) for k in legs]] if order_h == 0 else \
            [[('quartic', pr)] + [('cubic', (k,)) for k in legs if k not in pr] for pr in itertools.combinations(legs, 2)]
    for part in parts:
        nn = len(part)
        for perm in itertools.permutations(range(nn)):
            word = [part[i] for i in perm]
            odd_seq = [EXT[k][1] for (_, ids) in word for k in ids if EXT[k][0] == 'o']
            s = 1
            for a in range(len(odd_seq)):
                for b in range(a + 1, len(odd_seq)):
                    if odd_seq[a] > odd_seq[b]: s = -s
            hs = [PiH(letter(kd, ids)) for kd, ids in word]
            if not Kmode:
                tot += sp.Rational(1, 2 * nn) * s * STr([x for pair in zip([P] * nn, hs) for x in pair])
            else:
                # derivative in delta: K in each position; delta moved to the front past the odd r's before it
                for pos in range(nn):
                    props = [K if j == pos else P for j in range(nn)]
                    odd_before = sum(1 for (_, ids) in word[:pos] for k in ids if EXT[k][0] == 'o')
                    tot += -sp.Rational(1, 2 * nn) * s * (-1) ** odd_before * STr([x for pair in zip(props, hs) for x in pair])   # (Pi K)_{ij} = -(K Pi)_{ij}: K joins opposite parities
    return sp.expand(tot)
def coeff_exact(W, legs, order_g, order_h, delta=False):
    odd_ids = ([99] if delta else []) + [EXT[k][1] for k in legs if EXT[k][0] == 'o']
    c = strip(W, odd_ids)
    c = sp.expand(c).coeff(g, order_g).coeff(h, order_h)
    for k in legs:
        if EXT[k][0] == 'e': c = sp.expand(c).coeff(EXT[k][1], 1)
    for e in (e0, e1):
        c = sp.expand(c).subs(e, 0)
    return sp.expand(c)

W = W1_exact(False)
ok = True
for legs, og, oh in [((0, 1, 2), 3, 0), ((0, 1, 2), 1, 1), ((1, 2, 3), 3, 0), ((1, 2, 3), 1, 1), ((0, 1), 2, 0), ((1, 2), 2, 0), ((0, 3, 1), 3, 0)]:
    legs = tuple(sorted(legs))
    a = coeff_exact(W, legs, og, oh); b = word_formula(legs, og, oh)
    print('legs', legs, 'g^%d h^%d' % (og, oh), 'exact', a, 'words', b, a == b); ok &= (a == b)
WK = W1_exact(True)
for legs, og, oh in [((0, 1), 2, 0), ((1, 2), 2, 0), ((0, 3), 2, 0), ((0, 1, 2), 3, 0), ((0, 1, 2), 1, 1)]:
    a = coeff_exact(WK, legs, og, oh, delta=True); b = word_formula(legs, og, oh, Kmode=True)
    print('K-edge legs', legs, 'g^%d h^%d' % (og, oh), 'exact', a, 'words', b, a == b); ok &= (a == b)
print('ALL OK' if ok else 'MISMATCH')
