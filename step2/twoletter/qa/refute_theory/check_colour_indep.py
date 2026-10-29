"""Referee check of Lemma tl-colour (ii), (iii) with random matrices, plus the cases the lemma does not state:
 (a) the cyclically adjacent pair (X_n, X_1) of a single trace (pair order reversed relative to the trace order),
 (b) a pair across two traces of a multi-trace operator,
 (c) su(N) r = 0 general position.
Colour element of a kernel: C_w = sum_A T_A (x) w(T^A), w = ad_{Y_1} ... ad_{Y_r} (word as in Lemma tl-colour (i)),
first factor in the slot of X_a, second in the slot of X_b; kappa = 1.
Output: check_colour_indep.json
"""
import json, os, itertools
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(20261029)


def basis_gl(N):
    E = []
    for i in range(N):
        for j in range(N):
            m = np.zeros((N, N), complex); m[i, j] = 1; E.append(m)
    dual = [m.T.copy() for m in E]          # tr(E_ij E_ji) = 1
    return E, dual


def basis_su(N):
    E, dual = basis_gl(N)
    # project onto traceless: sum_A T_A (x) T^A |_su = sum_gl - (1/N) 1 (x) 1
    return E, dual


def ad_word(Ys, T):
    out = T
    for Y in reversed(Ys):                  # w = ad_{Y_1} ... ad_{Y_r}: apply ad_{Y_r} first
        out = Y @ out - out @ Y
    return out


def Cw_apply(f, Ys, N, su):
    """sum_A f(T_A, w(T^A)) with the su(N) correction"""
    E, dual = basis_gl(N)
    s = sum(f(E[A], ad_word(Ys, dual[A])) for A in range(len(E)))
    if su:
        I = np.eye(N)
        s = s - (1.0 / N) * f(I, ad_word(Ys, I))
    return s


def rnd(N):
    return rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))


def Sexp(Ys):
    """w(T) = sum_S (-1)^{r-|S|} Y_S^-> T Y_{S^c}^<-  (lemma's expansion) as list of (sign, left, right)"""
    r = len(Ys)
    N = Ys[0].shape[0] if r else None
    terms = []
    for mask in itertools.product([0, 1], repeat=r):
        S = [i for i in range(r) if mask[i]]
        Sc = [i for i in range(r) if not mask[i]]
        L = np.eye(N, dtype=complex)
        for i in S:
            L = L @ Ys[i]
        R = np.eye(N, dtype=complex)
        for i in reversed(Sc):
            R = R @ Ys[i]
        terms.append(((-1) ** (r - len(S)), L, R))
    return terms


res = []
for N in [2, 3, 4]:
    for su in [False, True]:
        for r in [0, 1, 2, 3]:
            Ys = [rnd(N) for _ in range(r)]
            if su:
                Ys = [Y - np.trace(Y) / N * np.eye(N) for Y in Ys]
            M1, M2 = rnd(N), rnd(N)
            # adjacent: tr(M1 Xa Xb)
            lhs = Cw_apply(lambda a, b: np.trace(M1 @ a @ b), Ys, N, su)
            if r == 0:
                rhs_lemma = (N - (1.0 / N if su else 0)) * np.trace(M1)
            else:
                rhs_lemma = sum(sg * np.trace(L) * np.trace(M1 @ R) for sg, L, R in Sexp(Ys))
            # general position: tr(M1 Xa M2 Xb)
            lhs2 = Cw_apply(lambda a, b: np.trace(M1 @ a @ M2 @ b), Ys, N, su)
            if r == 0:
                rhs2 = np.trace(M1) * np.trace(M2) - ((1.0 / N) * np.trace(M1 @ M2) if su else 0)
            else:
                rhs2 = sum(sg * np.trace(M2 @ L) * np.trace(R @ M1) for sg, L, R in Sexp(Ys))
            # (a) cyclic: pair (X_b, X_a) appearing as tr(M1 Xb Xa)  [X_a first in the operator's pair order but
            #     second in the trace]: not covered by the lemma's formula; compute and compare with the naive
            #     reading "pair replaced by Y_r ... Y_1" (with N factor)
            lhs3 = Cw_apply(lambda a, b: np.trace(M1 @ b @ a), Ys, N, su)
            if r >= 1:
                Yrev = np.eye(N, dtype=complex)
                for Y in reversed(Ys):
                    Yrev = Yrev @ Y
                Yfwd = np.eye(N, dtype=complex)
                for Y in Ys:
                    Yfwd = Yfwd @ Y
                planar3_rev = N * (-1) ** r * np.trace(M1 @ Yrev)
                planar3_fwd = N * np.trace(M1 @ Yfwd)
                # leading-N part of lhs3: the S = full set term: tr(Y^->)... compute exactly via expansion
                lead3 = sum(sg * np.trace(R) * np.trace(M1 @ L) for sg, L, R in Sexp(Ys))
            else:
                planar3_rev = planar3_fwd = lead3 = None
            # (b) cross-trace pair: tr(M1 Xa) tr(M2 Xb)
            lhs4 = Cw_apply(lambda a, b: np.trace(M1 @ a) * np.trace(M2 @ b), Ys, N, su)
            if r >= 1:
                rhs4 = sum(sg * np.trace(M1 @ R @ M2 @ L) for sg, L, R in Sexp(Ys))
            else:
                rhs4 = np.trace(M1 @ M2) - ((1.0 / N) * np.trace(M1) * np.trace(M2) if su else 0)
            res.append(dict(N=N, su=su, r=r,
                            adjacent_ok=bool(abs(lhs - rhs_lemma) < 1e-8 * (1 + abs(lhs))),
                            general_ok=bool(abs(lhs2 - rhs2) < 1e-8 * (1 + abs(lhs2))),
                            reversed_pair_equals_lemma_with_S_and_Sc_swapped=(None if r == 0 else bool(abs(lhs3 - lead3) < 1e-8 * (1 + abs(lhs3)))),
                            cross_trace_formula_ok=bool(abs(lhs4 - rhs4) < 1e-8 * (1 + abs(lhs4)))))

summary = dict(adjacent=all(x['adjacent_ok'] for x in res), general=all(x['general_ok'] for x in res),
               reversed=all(x['reversed_pair_equals_lemma_with_S_and_Sc_swapped'] for x in res if x['r'] > 0),
               cross_trace=all(x['cross_trace_formula_ok'] for x in res))
# (iii) normalisation: tr_g(ad_X ad_Y) = 2N tr(XY) on su(N)
iii = []
for N in [2, 3, 4]:
    X, Y = rnd(N), rnd(N)
    X -= np.trace(X) / N * np.eye(N); Y -= np.trace(Y) / N * np.eye(N)
    E, dual = basis_gl(N)
    t = sum(np.trace(dual[A] @ (X @ (Y @ E[A] - E[A] @ Y) - (Y @ E[A] - E[A] @ Y) @ X)) for A in range(len(E)))
    iii.append(dict(N=N, ratio=complex(t / np.trace(X @ Y)).real))
out = dict(summary=summary, rows=res, iii_trg_over_tr=iii,
           note='cross-trace pairs merge the two traces: tr(M1 Xa) tr(M2 Xb) -> sum_S (-1)^{r-|S|} tr(M1 Y_{S^c}^<- M2 Y_S^->); '
                'this case (multi-trace operators) is not stated in Lemma tl-colour (ii).')
json.dump(out, open(os.path.join(HERE, 'check_colour_indep.json'), 'w'), indent=1)
print(json.dumps(dict(summary=summary, iii=iii), indent=1))
