"""0-dim brute force: sign of the special-edge one-loop bubble.
D0 = (1/2){I3,I3}_K, {f,g}_K = Delta_K(fg) - Delta_K(f) g - (-1)^{|f|} f Delta_K(g), Delta_K = -d_K.
Compare d_P D0 (one-loop, degree 2) with c * str(K S H P S H) (H = Hessian of I3)."""
import sys, random, itertools
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import sympy as sp
from grass import GA, PAR, newvar, dP
random.seed(7)
par = [0, 0, 0, 1, 1, 1, 0, 1, 1, 0]
V = [newvar(500 + i, p) for i, p in enumerate(par)]
def rnd(): return sp.Rational(random.randint(-4, 4), random.randint(1, 2))
P = {}; K = {}
PPAIRS = [(0, 1), (2, 6), (3, 4), (7, 8), (5, 7), (9, 9)]
for a in V:
    for b in V:
        if (V.index(a), V.index(b)) in PPAIRS:
            x = rnd(); P[(a, b)] = x; P[(b, a)] = x * (-1) ** (PAR[a] * PAR[b])
        if a < b and PAR[a] != PAR[b]:
            x = rnd(); K[(a, b)] = x; K[(b, a)] = x
def monos(deg):
    out = []
    for comb in itertools.combinations_with_replacement(V, deg):
        odd = [v for v in comb if PAR[v]]
        if len(set(odd)) < len(odd) or len(odd) % 2: continue
        out.append(comb)
    return out
def rpoly_notad(deg, n):
    F = GA(); ms = monos(deg); random.shuffle(ms); k = 0
    for m in ms:
        t = GA.const(rnd())
        for v in m: t = t * GA.var(v)
        if dP(t, P, V).t or dP(t, K, V).t: continue
        F = F + t; k += 1
        if k == n: break
    return F
I3 = rpoly_notad(3, 25)
def dK(F): return dP(F, K, V)
def brK(f, g, pf):
    return dK(f * g) * (-1) - dK(f) * (-1) * g - (f * (dK(g) * (-1))) * ((-1) ** pf)
D0 = brK(I3, I3, 0) * sp.Rational(1, 2)
X = dP(D0, P, V).degree_part(2).map(sp.expand)
def hess(I): return {(b, c): I.dleft(c).dleft(b) for b in V for c in V}
def KSH(M_, H):
    out = {}
    for a in V:
        for c in V:
            s_ = GA()
            for b in V:
                if M_.get((a, b), 0) != 0:
                    s_ = s_ + H[(b, c)] * (M_[(a, b)] * (-1) ** PAR[b])
            out[(a, c)] = s_
    return out
def mm(X_, Y_): return {(a, c): sum((X_[(a, b)] * Y_[(b, c)] for b in V), GA()) for a in V for c in V}
def strace(X_): return sum((X_[(a, a)] * (-1) ** PAR[a] for a in V), GA())
H = hess(I3)
T = strace(mm(KSH(K, H), KSH(P, H))).degree_part(2).map(sp.expand)
c = sp.symbols('c')
keys = set(X.t) | set(T.t)
sol = sp.solve([sp.expand(X.t.get(k, 0) - c * T.t.get(k, 0)) for k in keys], [c], dict=True)
print('terms', len(X.t), len(T.t), 'd_P D0 = c str(KSHPSH):', sol)
# theta identity: W(P + theta K) with theta odd: brute force one-loop degree-2 part, theta-linear
TH = newvar(599, 1)
PK = dict(P)
for (a, b), v in K.items():
    PK[(a, b)] = PK.get((a, b), 0) + GA.var(TH) * v   # entries with odd theta (left)
# dP with GA-valued kernel entries: implement directly
def dPg(F, M_):
    r = GA()
    for a in V:
        Fa = F.dleft(a)
        if not Fa.t: continue
        for b in V:
            m = M_.get((a, b), 0)
            if isinstance(m, GA):
                r = r + m * Fa.dleft(b) * sp.Rational(1, 2)
            elif m != 0:
                r = r + Fa.dleft(b) * (m / 2)
    return r
hb = sp.symbols('hbar')
# one loop, quadratic in background: from e^{hbar d} e^{I3/hbar} : terms I3^2/(2 hbar^2) with two contractions -> hbar^0 ... W = hbar log
E2 = I3 * I3 * sp.Rational(1, 2)
Z = dPg(dPg(E2, PK), PK) * sp.Rational(1, 2)      # hbar^2/2! d^2 (I3^2/(2hbar^2)) -> hbar^0 ; W gets factor hbar -> order hbar
Zd = GA({k: sp.expand(v) for k, v in Z.t.items() if sum(1 for x in k if x != TH) == 2 and TH in k})
# remove disconnected: (d I3)^2 terms vanish (tadpole-free for P and K) ; theta-linear part:
Wth = Zd                       # coefficient of hbar^1 in W, theta-linear, degree 2
# -d_theta (left): remove TH
def dth(F): return F.dleft(TH)
lhs = dth(Wth).map(sp.expand)
sol2 = sp.solve([sp.expand(-lhs.t.get(k, 0) - c * T.t.get(k, 0)) for k in set(lhs.t) | set(T.t)], [c], dict=True)
print('-d_theta W(P+theta K)|_1loop = c str(KSHPSH):', sol2)
sol3 = sp.solve([sp.expand(-lhs.t.get(k, 0) - c * X.t.get(k, 0)) for k in set(lhs.t) | set(X.t)], [c], dict=True)
print('-d_theta W(P+theta K) = c * d_P D0 :', sol3)
