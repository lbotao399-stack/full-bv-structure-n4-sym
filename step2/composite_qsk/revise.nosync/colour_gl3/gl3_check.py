"""GC12 (gl(3) realisation): the planar adjacent rule of Lemma c10-colour (iv) for the colour words of the composite
kernels, C_w = (1 (x) ad_{Y_wr}...ad_{Y_w1})(-,-)^{-1}, contracted into tr(M1 X_a X_b).
For gl(N) with (X,Y) = kappa tr(XY) the Casimir is sum_A t_A (x) t^A = kappa^{-1} P (no 1/N term), and
   sum_A tr(M1 t_A w(t^A)),  w = ad_{Y_r}...ad_{Y_1},
equals kappa^{-1} sum_{S subset {1..r}} (-1)^{r-|S|} tr(Y_S ordered) tr(M1 Y_{S^c} reversed): the N-part is
(-1)^r N kappa^{-1} tr(M1 Y_1...Y_r) (S empty), the planar adjacent rule; checked for r = 1, 2, 3 with random complex matrices,
orthonormal and rescaled (non-orthonormal) bases."""
import itertools, json, os
import numpy as np

rng = np.random.default_rng(7)


def gl_basis(N, rescale):
    B = []
    for i in range(N):
        for j in range(N):
            M = np.zeros((N, N), complex); M[i, j] = 1; B.append(M)
    if rescale:
        A = rng.normal(size=(len(B), len(B))) + 3 * np.eye(len(B))
        B = [sum(A[a, b] * B[b] for b in range(len(B))) for a in range(len(B))]
    return B


def check(N, rescale, kappa=-1.5):
    B = gl_basis(N, rescale); n = len(B)
    G = np.array([[kappa * np.trace(B[a] @ B[b]) for b in range(n)] for a in range(n)])
    Gi = np.linalg.inv(G)
    dual = [sum(Gi[a, b] * B[b] for b in range(n)) for a in range(n)]
    rnd = lambda: rng.normal(size=(N, N)) + 1j * rng.normal(size=(N, N))
    out = {}
    M1 = rnd()
    for r in (1, 2, 3):
        Ys = [rnd() for _ in range(r)]
        def w(T):
            for Y in Ys:
                T = Y @ T - T @ Y
            return T
        full = sum(np.trace(M1 @ B[a] @ w(dual[a])) for a in range(n))
        # expansion of ad_{Y_r}...ad_{Y_1}(T) = sum_S (-1)^{r-|S|} Y_S^{left} T Y_{S^c}^{right}
        tot = 0
        for S in itertools.product((0, 1), repeat=r):
            left = np.eye(N, dtype=complex); right = np.eye(N, dtype=complex)
            for k, Y in enumerate(Ys):          # apply in order Y_1 first: T -> Y T - T Y
                if S[k]:
                    left = Y @ left
                else:
                    right = right @ Y
            sign = (-1) ** (r - sum(S))
            # sum_A tr(M1 t_A L t^A R) = kappa^{-1} tr(L) tr(M1 R)  ...  wait order: t_A (x) t^A = kappa^{-1} P
            tot += sign * np.trace(left) * np.trace(M1 @ right) / kappa
        out['r=%d exact gl identity (Casimir = P/kappa, no 1/N term)' % r] = bool(np.isclose(full, tot))
    return out


if __name__ == '__main__':
    res = {}
    for resc in (False, True):
        res['gl(3)%s' % (' rescaled' if resc else '')] = check(3, resc)
    res['all_pass'] = all(v for r in res.values() if isinstance(r, dict) for v in r.values())
    print(json.dumps(res, indent=1))
    json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gl3_check.json'), 'w'), indent=1)
