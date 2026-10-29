"""0-dim brute-force check of the one-loop special-edge formula used for the QME defect graphs (engine A, round 2).
Brute force: X = d_P( (1/2){I,I}_K ), {f,g}_K = Delta_K(fg) - Delta_K(f) g - (-1)^{|f|} f Delta_K(g), Delta_K = -d_K,
d_K = (1/2) sum K_ab d_b d_a  [C07 S6.4, S10.1; CG2 S7.2.8 product rule].  K pairs even with odd coordinates (field-antifield),
P pairs equal parities; I cubic, tadpole-free for P.  Candidate: X_quadratic = c * str(K S H P S H) with H_bc = d_b d_c I."""
import random, itertools
import sympy as sp
from grass import GA, PAR, newvar, dP
random.seed(11)
par = [0, 0, 0, 1, 1, 1, 0, 1, 1, 0]
V = [newvar(200 + i, p) for i, p in enumerate(par)]
def rnd():
    return sp.Rational(random.randint(-4, 4), random.randint(1, 2))
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
        if len(set(odd)) < len(odd) or len(odd) % 2:
            continue
        out.append(comb)
    return out
def rpoly_notad(deg, n):
    F = GA(); ms = monos(deg); random.shuffle(ms); k = 0
    for m in ms:
        t = GA.const(rnd())
        for v in m:
            t = t * GA.var(v)
        if dP(t, P, V).t:
            continue
        F = F + t; k += 1
        if k == n: break
    return F
I3 = rpoly_notad(3, 30)
print('I3 terms', len(I3.t), 'P entries', len(P), 'K entries', len(K))
def DeltaK(Fn):
    return dP(Fn, K, V) * (-1)
def brK(f, g, pf):
    return DeltaK(f * g) - DeltaK(f) * g - (f * DeltaK(g)) * ((-1) ** pf)
X = dP(brK(I3, I3, 0) * sp.Rational(1, 2), P, V).degree_part(2).map(sp.expand)
def hess(I):
    return {(b, c): I.dleft(c).dleft(b) for b in V for c in V}
def KSH(Kmat, H):
    M = {}
    for a in V:
        for c in V:
            s_ = GA()
            for b in V:
                if Kmat.get((a, b), 0) != 0:
                    s_ = s_ + H[(b, c)] * (Kmat[(a, b)] * (-1) ** PAR[b])
            M[(a, c)] = s_
    return M
def mm(X_, Y_):
    return {(a, c): sum((X_[(a, b)] * Y_[(b, c)] for b in V), GA()) for a in V for c in V}
def strace(X_):
    return sum((X_[(a, a)] * (-1) ** PAR[a] for a in V), GA())
H = hess(I3)
T1 = strace(mm(KSH(K, H), KSH(P, H))).degree_part(2)
T2 = strace(mm(KSH(P, H), KSH(K, H))).degree_part(2)
# solve X = c1 T1 + c2 T2
c1, c2 = sp.symbols('c1 c2')
eqs = []
keys = set(X.t) | set(T1.t) | set(T2.t)
for k in keys:
    eqs.append(sp.expand(X.t.get(k, 0) - c1 * T1.t.get(k, 0) - c2 * T2.t.get(k, 0)))
sol = sp.solve(eqs, [c1, c2], dict=True)
print('terms', len(X.t), len(T1.t), len(T2.t))
print('T1 == T2 :', all(sp.expand(T1.t.get(k, 0) - T2.t.get(k, 0)) == 0 for k in keys))
print('solution X = c1 str(K S H P S H) + c2 str(P S H K S H):', sol)
