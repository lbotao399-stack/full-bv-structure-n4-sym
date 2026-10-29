"""Check of Lemma lem:tl-colour (ii): colour action of the kernel element C_w = (1 (x) w)(-,-)^{-1},
w = ad_{Y_1} ... ad_{Y_r}, on a single trace, for gl(N) and su(N), N = 2, 3, 4, r = 0..3, with random integer
matrices (exact integer arithmetic; kappa = 1).
  adjacent:  sum_A tr(M1 T_A w(T^A)) = (-1)^r N tr(M1 Y_r..Y_1) + sum_{S != empty} (-1)^{r-|S|} tr(Y_S^->) tr(M1 Y_{S^c}^<-)
  general :  sum_A tr(M1 T_A M2 w(T^A)) = sum_S (-1)^{r-|S|} tr(M2 Y_S^->) tr(Y_{S^c}^<- M1)
  su(N), r = 0: adjacent (N - 1/N) tr M1, general tr M1 tr M2 - (1/N) tr(M1 M2).
Writes checks_colour_results.json."""
import itertools, json, os
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261029)


def E(N, i, j):
    m = np.zeros((N, N), dtype=object); m[i, j] = 1; return m


def mat(N):
    return np.array(rng.integers(-3, 4, size=(N, N)), dtype=object)


def tr(m):
    return sum(m[i, i] for i in range(m.shape[0]))


def mul(*ms):
    out = ms[0]
    for m in ms[1:]:
        out = out.dot(m)
    return out


def word(Ys, T):
    out = T
    for Y in reversed(Ys):
        out = Y.dot(out) - out.dot(Y)
    return out


def ordered(Ys, S, rev=False):
    N = Ys[0].shape[0] if Ys else None
    idx = sorted(S, reverse=rev)
    out = np.identity(N, dtype=object)
    for i in idx:
        out = out.dot(Ys[i])
    return out


def basis(N, alg):
    if alg == 'gl':
        B = [E(N, i, j) for i in range(N) for j in range(N)]
        D = [E(N, j, i) for i in range(N) for j in range(N)]
        return B, D
    # su(N): E_ij (i != j) and H_k = E_kk - E_{k+1,k+1}; dual basis w.r.t. tr via the inverse Gram matrix (Fractions)
    B = [E(N, i, j) for i in range(N) for j in range(N) if i != j]
    for k in range(N - 1):
        B.append(E(N, k, k) - E(N, k + 1, k + 1))
    G = [[Fraction(tr(x.dot(y))) for y in B] for x in B]
    n = len(B)
    # invert G exactly
    A = [row[:] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(G)]
    for c in range(n):
        p = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[p] = A[p], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [x - f * y for x, y in zip(A[r], A[c])]
    Ginv = [row[n:] for row in A]
    D = []
    for a in range(n):
        m = np.zeros((N, N), dtype=object)
        for b in range(n):
            m = m + Ginv[a][b] * B[b]
        D.append(m)
    return B, D


def run():
    res = []
    ok_all = True
    for N in (2, 3, 4):
        for alg in ('gl', 'su'):
            B, D = basis(N, alg)
            for r in range(0, 4):
                for trial in range(2):
                    M1, M2 = mat(N), mat(N)
                    Ys = [mat(N) for _ in range(r)]
                    adj = sum(tr(mul(M1, Ta, word(Ys, Td))) for Ta, Td in zip(B, D))
                    gen = sum(tr(mul(M1, Ta, M2, word(Ys, Td))) for Ta, Td in zip(B, D))
                    if r == 0 and alg == 'su':
                        adj_f = Fraction(N) * tr(M1) - Fraction(1, N) * tr(M1)
                        gen_f = tr(M1) * tr(M2) - Fraction(1, N) * tr(mul(M1, M2))
                    else:
                        I = np.identity(N, dtype=object)
                        adj_f = 0; gen_f = 0
                        for size in range(r + 1):
                            for S in itertools.combinations(range(r), size):
                                Sc = [i for i in range(r) if i not in S]
                                sg = (-1) ** (r - len(S))
                                YS = ordered(Ys, S) if r else I
                                YSc = ordered(Ys, Sc, rev=True) if r else I
                                adj_f += sg * tr(YS) * tr(mul(M1, YSc))
                                gen_f += sg * tr(mul(M2, YS)) * tr(mul(YSc, M1))
                        if r >= 1:
                            # display form: S = empty term is (-1)^r N tr(M1 Y_r..Y_1)
                            pass
                    ok = (Fraction(adj) == Fraction(adj_f)) and (Fraction(gen) == Fraction(gen_f))
                    ok_all &= ok
                    res.append(dict(N=N, alg=alg, r=r, trial=trial, ok=bool(ok)))
    # the S = empty term of the adjacent formula is (-1)^r N tr(M1 Y_r ... Y_1): check separately
    N = 3; Ys = [mat(N) for _ in range(3)]; M1 = mat(N)
    I = np.identity(N, dtype=object)
    sgn_ok = ((-1) ** 3 * tr(I) * tr(mul(M1, Ys[2], Ys[1], Ys[0])) ==
              (-1) ** 3 * tr(I) * tr(mul(M1, ordered(Ys, [0, 1, 2], rev=True))))
    out = dict(description=__doc__, cases=res, empty_set_term_ok=bool(sgn_ok), all_pass=bool(ok_all and sgn_ok))
    with open(os.path.join(HERE, 'checks_colour_results.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print('colour identities:', 'ALL PASS' if out['all_pass'] else 'FAIL', '(%d cases)' % len(res))


if __name__ == '__main__':
    run()
