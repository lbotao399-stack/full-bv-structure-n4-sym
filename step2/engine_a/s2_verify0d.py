"""Step 2, engine A: 0-dim brute-force validation of the one-loop special-edge graph formulas.

Facts used (exact identities, derived in s2_notes.txt):
  D_0 part of Q W(P,J) = hbar Delta_K W + (1/2){W,W}_K   (K = K_eps - K_0, valid when partial_K I = 0),
so at order hbar:  hbar Delta_K W_0 + {W_0, W_1}_K.
W_0 = tree generating function, dW_0/ds = (1/2) C_P(W_0,W_0)  =>  W_0 = I + (1/2)C_P(I,I) + (1/4)C_P(I,C_P(I,I)) + ...
C_X(f,g) = d_X(fg) - d_X(f) g - f d_X(g),  d_X = (1/2) sum X_ab d_b d_a (left derivatives, [C07 S6.4]),
{f,g}_X = -C_X(f,g)  (Delta_X = -d_X, product rule of [CG2 S7.2.8]).
Claim tested here (cubic background degree, cubic I_3 + quartic I_4, d_K I = 0):
  Delta_K W_0 |_deg3  =  (1/2) str[K S H3 P S H3 P S H3] + (1/2)(str[K S H3 P S H4] + str[K S H4 P S H3])      (I-a)
                        + C_P( (1/2) str[K S H3 P S H3] |_deg2 , I_3 )                                            (I-d)
and, in degree 2:  Delta_K( (1/2) C_P(I3,I3) )|_deg2 = d_P((1/2){I3,I3}_K)|_loop = (1/2) str[K S H3 P S H3].
S = parity operator, H = left Hessian (l-slot differentiated first), K odd field-antifield kernel, P even kernel."""
import sys, random, itertools
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import sympy as sp
from grass import GA, PAR, newvar, dP

def run(seed, ptadfree=False):
    random.seed(seed)
    base = 500 + 100 * seed
    par = [0, 0, 1, 1, 0, 1, 1, 0, 1, 0]
    V = [newvar(base + i, p) for i, p in enumerate(par)]
    def rnd():
        return sp.Rational(random.randint(-4, 4), random.randint(1, 3))
    P = {}; K = {}
    for i, a in enumerate(V):
        for j, b in enumerate(V):
            if j < i:
                continue
            if PAR[a] == PAR[b] and random.random() < 0.5:
                x = rnd()
                if a == b and PAR[a]:
                    continue
                P[(a, b)] = x; P[(b, a)] = x * (-1) ** (PAR[a] * PAR[b])
            if PAR[a] != PAR[b] and random.random() < 0.6:
                x = rnd(); K[(a, b)] = x; K[(b, a)] = x
    def monos(deg):
        out = []
        for comb in itertools.combinations_with_replacement(V, deg):
            odd = [v for v in comb if PAR[v]]
            if len(set(odd)) < len(odd) or len(odd) % 2:
                continue
            out.append(comb)
        return out
    def rpoly(deg, n):
        F = GA(); ms = monos(deg); random.shuffle(ms); k = 0
        for m in ms:
            t = GA.const(rnd())
            for v in m:
                t = t * GA.var(v)
            if not t.t or dP(t, K, V).t:          # K-tadpole free monomials only
                continue
            if ptadfree and deg == 3 and dP(t, P, V).t:   # optionally also P-tadpole free (as in the field theory)
                continue
            F = F + t; k += 1
            if k == n:
                break
        return F
    I3 = rpoly(3, 25); I4 = rpoly(4, 25)
    assert not dP(I3, K, V).t and not dP(I4, K, V).t
    def C(X, f, g):
        return dP(f * g, X, V) - dP(f, X, V) * g - f * dP(g, X, V)
    def DeltaK(F):
        return dP(F, K, V) * (-1)
    def brK(f, g, pf):
        return DeltaK(f * g) - DeltaK(f) * g - (f * DeltaK(g)) * ((-1) ** pf)
    W0_5 = C(P, I3, I4) + C(P, I3, C(P, I3, I3)) * sp.Rational(1, 4)
    LHS3 = DeltaK(W0_5).degree_part(3).map(sp.expand)
    W0_4 = C(P, I3, I3) * sp.Rational(1, 2)
    LHS2 = DeltaK(W0_4).degree_part(2).map(sp.expand)
    LHS2b = dP(brK(I3, I3, 0) * sp.Rational(1, 2), P, V).degree_part(2).map(sp.expand)
    def hess(I):
        return {(b, c): I.dleft(c).dleft(b) for b in V for c in V}
    def XSH(X, H):
        M = {}
        for a in V:
            for c in V:
                s_ = GA()
                for b in V:
                    if X.get((a, b), 0) != 0:
                        s_ = s_ + H[(b, c)] * (X[(a, b)] * (-1) ** PAR[b])
                M[(a, c)] = s_
        return M
    def mm(X_, Y_):
        return {(a, c): sum((X_[(a, b)] * Y_[(b, c)] for b in V), GA()) for a in V for c in V}
    def strace(X_):
        return sum((X_[(a, a)] * (-1) ** PAR[a] for a in V), GA())
    H3 = hess(I3); H4 = hess(I4)
    KH3, PH3, KH4, PH4 = XSH(K, H3), XSH(P, H3), XSH(K, H4), XSH(P, H4)
    tri = strace(mm(mm(KH3, PH3), PH3)).degree_part(3)
    bub = (strace(mm(KH3, PH4)) + strace(mm(KH4, PH3))).degree_part(3)
    sb2 = (strace(mm(KH3, PH3)) * sp.Rational(1, 2)).degree_part(2)
    Id = C(P, sb2, I3).degree_part(3)
    RHS3 = (tri * sp.Rational(1, 2) + bub * sp.Rational(1, 2) + Id).map(sp.expand)
    def diff(A, B):
        keys = set(A.t) | set(B.t)
        return sum(1 for k in keys if sp.expand(A.t.get(k, 0) - B.t.get(k, 0)) != 0), len(keys)
    out = {'seed': seed, 'P_tadpole_free_I3': ptadfree, 'terms_LHS3': len(LHS3.t), 'mismatch_deg3': diff(LHS3, RHS3),
           'mismatch_deg2_DeltaK_W0': diff(LHS2, sb2), 'mismatch_deg2_dP_bracket': diff(LHS2b, sb2),
           'terms_triangle': len(tri.t), 'terms_quartic_bubble': len(bub.t), 'terms_Id': len(Id.t)}
    # sign sensitivity: the wrong overall sign of the I-a part must fail
    wrong = (tri * sp.Rational(-1, 2) + bub * sp.Rational(1, 2) + Id).map(sp.expand)
    out['control_wrong_sign_mismatch'] = diff(LHS3, wrong)
    return out

if __name__ == '__main__':
    for sd in (1, 2, 3):
        print(run(sd), flush=True)
    # d_P((1/2){I,I}_K)|_deg2 contains, besides the loop, P-tadpoles on one vertex (absent in the field theory: tr ad = 0);
    # with d_P I_3 = 0 it equals (1/2) str[K S H3 P S H3] as well:
    for sd in (4, 5):
        print(run(sd, ptadfree=True), flush=True)
