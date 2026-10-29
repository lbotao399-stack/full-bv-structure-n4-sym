"""0-dim check of the one-loop formula  W_1 = sum_n 1/(2n) str[(P S d2I)^n]  (S = parity operator,
d2I_{bc} = d_b d_c I with left derivatives, P_{ab} = partial_P(y_a y_b)), against brute-force
hbar log(e^{hbar partial_P} e^{I/hbar})."""
import random, sympy as sp
from grass import GA, PAR, newvar, dP
random.seed(3)
hb, g = sp.symbols('hbar g')
par = [0, 0, 1, 1, 1, 1]
V = [newvar(i, p) for i, p in enumerate(par)]
def rnd():
    return sp.Rational(random.randint(-5, 5), random.randint(1, 3))
P = {}
for a in V:
    for b in V:
        if a <= b and PAR[a] == PAR[b]:
            if PAR[a] == 1 and a == b:
                continue
            x = rnd(); P[(a, b)] = x
            P[(b, a)] = x * (-1) ** (PAR[a] * PAR[b])
def monos(deg):
    import itertools
    out = []
    for comb in itertools.combinations_with_replacement(V, deg):
        odd = [v for v in comb if PAR[v]]
        if len(set(odd)) < len(odd) or len(odd) % 2:
            continue
        out.append(comb)
    return out
def rpoly(deg, n):
    F = GA()
    ms = monos(deg); random.shuffle(ms)
    for m in ms[:n]:
        t = GA.const(rnd())
        for v in m:
            t = t * GA.var(v)
        F = F + t
    return F
def rpoly_notad(deg, n):
    F = GA()
    ms = monos(deg); random.shuffle(ms); k = 0
    for m in ms:
        t = GA.const(rnd())
        for v in m:
            t = t * GA.var(v)
        if dP(t, P, V).t:   # drop monomials with a self-contraction (tadpoles)
            continue
        F = F + t; k += 1
        if k == n: break
    return F
I3 = rpoly_notad(3, 14); I4 = rpoly(4, 12)
print('tadpole check', dP(I3, P, V).t)
# brute force
E0 = GA.const(1)
E1 = I3 * (1 / hb)
E2 = I4 * (1 / hb) + I3 * I3 * (1 / (2 * hb ** 2))
def expdP(F, m):
    out = F; term = F
    for k in range(1, m + 1):
        term = dP(term, P, V) * (hb / k)
        out = out + term
    return out
F1 = expdP(E1, 2); F2 = expdP(E2, 4)
W2 = (F2 - F1 * F1 * sp.Rational(1, 2)) * hb
Wq = W2.degree_part(2).map(lambda c: sp.expand(c))
W1q = Wq.map(lambda c: sp.expand(c).coeff(hb, 1))
# formula
def hess(I):
    return {(b, c): I.dleft(c).dleft(b) for b in V for c in V}
def PSH(H):
    M = {}
    for a in V:
        for c in V:
            s = GA()
            for b in V:
                if P.get((a, b), 0) != 0:
                    s = s + H[(b, c)] * (P[(a, b)] * (-1) ** PAR[b])
            M[(a, c)] = s
    return M
def mm(X, Y):
    return {(a, c): sum((X[(a, b)] * Y[(b, c)] for b in V), GA()) for a in V for c in V}
def strace(X):
    return sum((X[(a, a)] * (-1) ** PAR[a] for a in V), GA())
X3 = PSH(hess(I3)); X4 = PSH(hess(I4))
form = (strace(X4) * sp.Rational(1, 2) + strace(mm(X3, X3)) * sp.Rational(1, 4)).degree_part(2)
diff = (W1q - form).map(lambda c: sp.simplify(c))
print('brute terms', len(W1q.t), 'formula terms', len(form.t), 'diff terms', len(diff.t))
# variants for comparison: plain trace, no parity
def trace(X):
    return sum((X[(a, a)] for a in V), GA())
alt = (trace(X4) * sp.Rational(1, 2) + trace(mm(X3, X3)) * sp.Rational(1, 4)).degree_part(2)
print('plain-trace variant diff terms', len((W1q - alt).map(sp.simplify).t))
