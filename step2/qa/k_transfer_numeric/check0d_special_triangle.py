"""0-dim brute-force check (Grassmann algebra of engine_a/grass.py) of the one-loop formulas used in this folder.
Model: field coordinates V (mixed parity), P pairs equal parities (graded symmetric), K pairs opposite parities
(field-antifield), Delta_K = -d_K, {f,g}_K = Delta_K(fg) - Delta_K(f) g - (-1)^{|f|} f Delta_K(g)  [C07 S6.4, S10.1;
CG2 S7.2.8].  J = I3 + J2 with J2 = gamma * (odd quadratic), gamma an odd constant (the ghost of the even generator
K_v; no edge touches it).  I3 and J2 are tadpole free for P.
Brute force: the one-loop part of d/dsigma W(P, J + sigma*D) with D a special-edge vertex, linear in the extra vertex Y:
   L1 = [d_P^2(D Y) - (d_P^2 D) Y]/2 = X(d_P D, Y) + (1/2) XX(D, Y),
X(f,g) = d_P(fg) - (d_P f) g - f d_P g (one cross contraction; X(d_P D, Y) = the graph with Y on a tree leg),
(1/2) XX(D,Y) = the graphs with the special edge and two P edges on the cycle (1PI) or the special edge as a bridge.
Candidates (S = parity, H = Hessian with left derivatives):
   D-class  D = (1/2){I3,I3}_K, Y = J2 :  (1/2) str(K S H3 P S H2 P S H3)
   B-class  D = {I3,J2}_K,      Y = I3 :  (1/2) [str(K S H2 P S H3 P S H3) + str(K S H3 P S H3 P S H2)] + bridge
   bridge  = c * {W1q, J2}_K  with W1q = (1/4) str(P S H3 P S H3) (the one-loop self-energy) -- coefficient c found.
Flow (no special edge): (1/6) sum_orderings str(P S Ha P S Hb P S Hc) = brute one-loop part of W(P, I3 + J2)."""
import sys, random, itertools, json
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import sympy as sp
from grass import GA, PAR, newvar, dP
random.seed(20261028)
par = [0, 0, 0, 1, 1, 1, 0, 1, 1, 0]
V = [newvar(300 + i, p) for i, p in enumerate(par)]
GAM = newvar(399, 1)                       # odd constant (ghost of K_v)
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
        if len(set(odd)) < len(odd):
            continue
        out.append(comb)
    return out
def rpoly_notad(deg, n, parity):
    F = GA(); ms = [m for m in monos(deg) if sum(PAR[v] for v in m) % 2 == parity]; random.shuffle(ms); k = 0
    for m in ms:
        t = GA.const(rnd())
        for v in m:
            t = t * GA.var(v)
        if dP(t, P, V).t:
            continue
        F = F + t; k += 1
        if k == n: break
    return F
I3 = rpoly_notad(3, 18, 0)
J2 = GA.var(GAM) * rpoly_notad(2, 8, 1)       # even overall
def dK(F): return dP(F, K, V)
def brK(f, g, pf):
    DK = lambda F: dK(F) * (-1)
    return DK(f * g) - DK(f) * g - (f * DK(g)) * ((-1) ** pf)
def X(f, g):
    return dP(f * g, P, V) - dP(f, P, V) * g - f * dP(g, P, V)
def deg2(F): return F.degree_part(2).map(sp.expand) if False else GA({k: sp.expand(c) for k, c in F.t.items() if sum(1 for v in k if v != GAM) == 2 and sp.expand(c) != 0})
def L1(D, Y):
    return (dP(dP(D * Y, P, V), P, V) - dP(dP(D, P, V), P, V) * Y) * sp.Rational(1, 2)
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
def same(A_, B_):
    keys = set(A_.t) | set(B_.t)
    return all(sp.simplify(A_.t.get(k, 0) - B_.t.get(k, 0)) == 0 for k in keys)
H3 = hess(I3); H2 = hess(J2)
out = {}
# ---- D-class
D0 = brK(I3, I3, 0) * sp.Rational(1, 2)
LD = deg2(L1(D0, J2)); PRD = deg2(X(dP(D0, P, V), J2))
onepiD = deg2(LD - PRD)
candD = deg2(strace(mm(mm(KSH(K, H3), KSH(P, H2)), KSH(P, H3))) * sp.Rational(1, 2))
out['D_class_1PI_equals_half_str(K S H3 P S H2 P S H3)'] = same(onepiD, candD)
out['D_class_terms'] = [len(onepiD.t), len(candD.t)]
# ---- B-class
D1 = brK(I3, J2, 0)
LB = deg2(L1(D1, I3)); PRB = deg2(X(dP(D1, P, V), I3))
onepiB = deg2(LB - PRB)
candB = deg2((strace(mm(mm(KSH(K, H2), KSH(P, H3)), KSH(P, H3))) + strace(mm(mm(KSH(K, H3), KSH(P, H3)), KSH(P, H2)))) * sp.Rational(1, 2))
W1q = deg2(strace(mm(KSH(P, H3), KSH(P, H3))) * sp.Rational(1, 4))
bridge = deg2(brK(W1q, J2, 0))
cB = sp.symbols('cB')
resid = (onepiB - candB - bridge * cB)
eqs = [sp.expand(c) for c in resid.t.values()]
sol = sp.solve(eqs, [cB], dict=True)
out['B_class_terms'] = [len(onepiB.t), len(candB.t), len(bridge.t)]
out['B_class: 1PI+bridge = half str(..) + c*{W1q,J2}_K, c ='] = str(sol)
# ---- the one-loop D-vertex alone: (1/2) str(K S H3 P S H3)  (engine A check, repeated)
out['bubble_D0_equals_half_str(KSH3PSH3)'] = same(deg2(dP(D0, P, V)), deg2(strace(mm(KSH(K, H3), KSH(P, H3))) * sp.Rational(1, 2)))
# ---- flow (no special edge): one-loop W(P, I3+J2), part linear in J2 and quadratic in I3
hb = sp.symbols('hbar')
E3 = (I3 * I3 * J2) * sp.Rational(1, 2)        # coefficient of hbar^-3 in e^{J/hbar} at order I3^2 J2
# one-loop, connected, field degree 2: hbar^-3 * hbar^3 d_P^3/3! ... W = hbar log; linear in J2, I3^2:
# connected part of d_P^3(I3^2 J2/2)/6 minus disconnected pieces (I3 tadpole free, J2 tadpole free)
Z3 = dP(dP(dP(E3, P, V), P, V), P, V) * sp.Rational(1, 6)
# disconnected at degree 2 with 3 contractions: [I3 I3 loop (2 contractions) -> deg 2] * J2 untouched -> deg 4; none at deg2
# except (I3-I3 vacuum? no) ; tree-level pieces need fewer contractions. Subtract nothing; check against formula
flow_brute = deg2(Z3)
Hs = [H3, H3, H2]
form = GA()
for perm in set(itertools.permutations([0, 1, 2])):
    Ms = [KSH(P, Hs[i]) for i in perm]
    form = form + strace(mm(mm(Ms[0], Ms[1]), Ms[2]))
form = deg2(form * sp.Rational(1, 6) * sp.Rational(1, 2))   # H = H3 + H3' : the two I3 copies are identical -> 1/2! from e^J
# note: the expansion of (H3 + H2)^3 linear in H2 contains each ordering of (H3,H3,H2) with multiplicity 1 for the
# 3 distinct cyclic positions of H2: sum over the 3 distinct sequences, times 1/(2*3).
form2 = GA()
for seq in [(H2, H3, H3), (H3, H2, H3), (H3, H3, H2)]:
    Ms = [KSH(P, h) for h in seq]
    form2 = form2 + strace(mm(mm(Ms[0], Ms[1]), Ms[2]))
form2 = deg2(form2 * sp.Rational(1, 6))
out['flow_triangle_brute_equals_(1/6)sum_seq str(PSH PSH PSH)'] = same(flow_brute, form2)
out['flow_terms'] = [len(flow_brute.t), len(form2.t)]
print(json.dumps(out, indent=1))
json.dump(out, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_transfer_numeric/check0d_special_triangle.json', 'w'), indent=1)
# ---- refinement: separate the tadpole / tree-leg topologies that the brute force contains
def deg1(F): return GA({k: sp.expand(c) for k, c in F.t.items() if sum(1 for v in k if v != GAM) == 1 and sp.expand(c) != 0})
T1 = deg1(strace(mm(KSH(P, H3), KSH(P, H2))) * sp.Rational(1, 2))     # one-loop 1-pt function with I3 and J2
a_, b_ = sp.symbols('a_ b_')
tadK = deg2(brK(T1, I3, 1))      # T1 is odd (gamma), bracket parity of first argument
residD = onepiD - candD * a_ - tadK * b_
solD = sp.solve([sp.expand(c) for c in residD.t.values()], [a_, b_], dict=True)
out['D_class: 1PI_brute = a*half str(KSH3PSH2PSH3) + b*{T1,I3}_K (special tadpole)'] = str(solD)
residF = flow_brute - form2 * a_ - deg2(X(W1q, J2)) * b_
c_ = sp.symbols('c_')
residF = residF - deg2(X(T1, I3)) * c_
solF = sp.solve([sp.expand(c) for c in residF.t.values()], [a_, b_, c_], dict=True)
out['flow: brute = a*(1/6)sum str(PSH)^3 + b*X(W1q,J2) + c*X(T1,I3)'] = str(solF)
print(json.dumps(out, indent=1))
json.dump(out, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_transfer_numeric/check0d_special_triangle.json', 'w'), indent=1)
