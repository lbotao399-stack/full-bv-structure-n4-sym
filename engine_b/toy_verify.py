"""0-dim verification of the one-loop formula used by engine B.
Truth: W = hbar log( e^{hbar d_P} e^{I/hbar} ),  d_P = 1/2 sum_ij P^{ij} D_j D_i (D_i applied first; left derivatives),
so that d_P(f_X f_Y) = P^{XY} (contraction of an adjacent ordered pair).
Candidate: W_1 = (hbar/2) sum_n (1/n) STr[(P S H)^n],  H_ij = D_i D_j I (D_j first),
STr X = sum_i (-1)^{|i|} X_ii, S = diag((-1)^{|j|}) or 1.  Internal variables propagate, external ones are background."""
import itertools, random, sympy as sp
from grassmann import G

random.seed(7)
hb, g, h, d = sp.symbols('hbar g h d')
x = sp.symbols('x0 x1'); e0 = sp.Symbol('e0')
ODD_INT = [0, 1, 2]; ODD_EXT = [10, 11]
internal = [('e', x[0]), ('e', x[1]), ('o', 0), ('o', 1), ('o', 2)]
par = [0, 0, 1, 1, 1]

def var(v):
    return G.scalar(v[1]) if v[0] == 'e' else G.odd(v[1])
def D(F, v):
    return F.d_even(v[1]) if v[0] == 'e' else F.d_odd(v[1])
def R():
    return sp.Integer(random.choice([-3, -2, -1, 1, 2, 3]))

n = len(internal)
P = [[sp.Integer(0)] * n for _ in range(n)]
for i in range(n):
    for j in range(i, n):
        if par[i] != par[j]:
            continue
        if par[i] == 0:
            P[i][j] = P[j][i] = R()
        elif i != j:
            P[i][j] = R(); P[j][i] = -P[i][j]

X = [var(v) for v in internal]
eta = [G.odd(10), G.odd(11)]
E0 = G.scalar(e0)
I3 = G()
for i in range(n):
    for j in range(i, n):
        if par[i] == par[j]:
            I3 = I3 + E0 * X[i] * X[j] * R()
        else:
            for k in range(2):
                I3 = I3 + eta[k] * X[i] * X[j] * R()
I4 = G()
for i in range(n):
    for j in range(i, n):
        if par[i] == par[j]:
            I4 = I4 + E0 * E0 * X[i] * X[j] * R() + eta[0] * eta[1] * X[i] * X[j] * R()
        else:
            I4 = I4 + E0 * eta[0] * X[i] * X[j] * R() + eta[1] * E0 * X[i] * X[j] * R()
O = G()
for i in range(n):
    for j in range(i, n):
        if par[i] == par[j]:
            O = O + X[i] * X[j] * R()
I = I3 * g + I4 * h + O * d

def trunc(expr):
    p = sp.Poly(sp.expand(expr), g, h, d)
    out = 0
    for (a, b, c), cf in p.terms():
        if a <= 2 and b <= 1 and c <= 1:
            out += cf * g**a * h**b * d**c
    return sp.expand(out)
def T(F):
    return F.map(trunc)

def dP(F):
    out = G()
    for i in range(n):
        for j in range(n):
            if P[i][j] != 0:
                out = out + D(D(F, internal[i]), internal[j]) * (P[i][j] / 2)
    return out

# exp(I/hbar) truncated
E = G.scalar(1); term = G.scalar(1)
for k in range(1, 5):
    term = T(term * I * (1 / hb) * sp.Rational(1, k))
    E = E + term
# e^{hbar dP}
res = G(); cur = E
for m in range(0, 4):
    res = res + cur * (hb**m / sp.factorial(m))
    cur = T(dP(cur))
res = T(res)
Y = T(res - G.scalar(1))
logv = Y + T(Y * Y) * sp.Rational(-1, 2) + T(T(Y * Y) * Y) * sp.Rational(1, 3)
W = T(logv * hb)
zero = {x[0]: 0, x[1]: 0}
W = W.subs(zero).drop_odd(ODD_INT)
def coeff_hbar(F, k):
    return F.map(lambda v: sp.expand(v).coeff(hb, k))
W1 = coeff_hbar(W, 1)

def candidate(S_par, order, strace):
    Hm = [[None] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if order == 'ij':
                Hm[i][j] = D(D(I, internal[j]), internal[i])
            else:
                Hm[i][j] = D(D(I, internal[i]), internal[j])
            Hm[i][j] = Hm[i][j].subs(zero).drop_odd(ODD_INT)
    M = [[G() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for k in range(n):
            acc = G()
            for j in range(n):
                if P[i][j] != 0:
                    s = (-1) ** par[j] if S_par else 1
                    acc = acc + Hm[j][k] * (P[i][j] * s)
            M[i][k] = T(acc)
    tot = G(); Mp = M
    for nn in range(1, 4):
        tr = G()
        for i in range(n):
            tr = tr + Mp[i][i] * ((-1) ** par[i] if strace else 1)
        tot = tot + tr * (hb / 2 / nn)
        Mp2 = [[G() for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for k in range(n):
                acc = G()
                for j in range(n):
                    acc = acc + Mp[i][j] * M[j][k]
                Mp2[i][k] = T(acc)
        Mp = Mp2
    return coeff_hbar(T(tot), 1)

if __name__ == '__main__':
    for S_par in (True, False):
        for order in ('ij', 'ji'):
            for strace in (True, False):
                C = candidate(S_par, order, strace)
                diff = W1 - C
                ok = diff.iszero()
                print(f"S_parity={S_par} H_order={order} supertrace={strace}: match={ok}")
    print("nonzero check, W1 has", len(W1.t), "odd-sectors; sample:", list(W1.t.items())[:2])

def check_contraction_rule():
    # d_P(X_i X_j) must equal P^{ij} (agreed edge rule: adjacent ordered pair -> P in the same order)
    bad = []
    for i in range(n):
        for j in range(n):
            v = dP(X[i] * X[j])
            if not (v - G.scalar(P[i][j])).iszero():
                bad.append((i, j))
    return bad
