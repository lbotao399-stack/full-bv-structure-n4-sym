"""Tree-level Q^xi_0 on the twisted letters, quadratic parts, engine realisation:
 Q0 lam = s'([c,lam] + 1/2 rho(F_A) xi),  F_MN = [A_M, A_N] on W_R^perp (rho(F) xi = sum_{M<N} F_MN s[M] sb[N] xi);
 Q0 lam* ⊃ s_*( dslash_A lam ) with the Yukawa sum_{M in W} sb[M][A_M, lam]  (covariant completion of op_Q)."""
import env, numpy as np, itertools, json, sympy as sp
import qa_twist as TW, qa_s1 as S1, sector as S
T = TW.Twist(); xi = T.xi
SB, SS = S1.SB, S1.SS
X = S1._X(xi)
go = [np.concatenate([np.zeros(4), g]) for g in T.gam_out]       # 10-vectors: field directions of gamma(j)
def gq(z): return S.gq(z)
res = {}
# c_{i;jk}: coefficient of [gamma(j), gamma(k)] (j<k) in 1/2 rho(F)xi projected on beta(i)
for i in range(3):
    b = T.beta_z[i]
    for j, k in ((0, 1), (0, 2), (1, 2)):
        tot = 0
        for M in range(4, 10):
            for N in range(4, 10):
                if M < N:
                    coef = b @ X[M, N]
                    tot += 0.5 * coef * (go[j][M] * go[k][N] - go[k][M] * go[j][N])
        if abs(tot) > 1e-12: res['c_%d;%d%d' % (i, j, k)] = str(gq(tot))
# Yukawa projection: v = W[(u';z)] (the lam* part of b); coefficient of [gamma(j), beta(k)]
v = T.W[(-1, 1) + T.z]
lamdir = [T.W[T.u + t] for t in T.taus]                              # field directions of beta(k)
for j in range(3):
    for k in range(3):
        tot = sum(go[j][M] * (v @ SB[M] @ lamdir[k]) for M in range(4, 10))
        if abs(tot) > 1e-12: res['d_%d%d' % (j, k)] = str(gq(tot))
print(json.dumps(res, indent=1)); json.dump(res, open('q0tree.json', 'w'), indent=1)
