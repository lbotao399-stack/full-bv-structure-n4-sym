"""Colour contractions of the word form C_w = (1 (x) ad_{Y_w1}...ad_{Y_wr})(-,-)^{-1} used by the gates, and their
numerical verification (gate G4) in su(2), su(3) (orthonormal and rescaled bases) and the planar rule in gl(3)/su(3).
Rules (k := tr(ad_X ad_Y)/(X,Y), so hbar k/(16 pi^2) is the unit u):
  pairing m(x,y) = (x,y):            m(C_w) = tr_g(w):  r=0 dim g, r=1 0, r=2 k (Y1,Y2), r=3 (k/2)(Y1,[Y2,Y3])
  rotation m(x,y) = (zeta,[x,y]):    m(C_w) = -tr_g(ad_zeta w): r=1 -k (zeta,Y1), r=2 -(k/2)(zeta,[Y1,Y2])
  planar (su(N), (X,Y) = kappa tr XY): tr(.. X_a X_b ..) -> N kappa^{-1} (-1)^r tr(.. Y_r..Y_1 ..) + multi-traces,
  so a kernel c hbar/(16 pi^2) C_w contributes (-1)^r c/2 units (Lemma lem:tl-colour (iii))."""
import itertools
import numpy as np


def su_basis(N, rng, rescale=False):
    """basis of su(N) (anti-hermitian) as matrices; pairing (X,Y) = kappa tr(XY), kappa = -3/2 (rescaled) or -1."""
    B = []
    for i in range(N):
        for j in range(i + 1, N):
            M = np.zeros((N, N), complex); M[i, j] = 1; M[j, i] = -1; B.append(M)
            M = np.zeros((N, N), complex); M[i, j] = 1j; M[j, i] = 1j; B.append(M)
    for d in range(1, N):
        M = np.zeros((N, N), complex)
        for i in range(d):
            M[i, i] = 1j
        M[d, d] = -1j * d
        B.append(M)
    if rescale:   # deliberately non-orthonormal basis
        A = rng.normal(size=(len(B), len(B))) + 3 * np.eye(len(B))
        B = [sum(A[a, b] * B[b] for b in range(len(B))) for a in range(len(B))]
    return B


def check(N, rescale, rng):
    B = su_basis(N, rng, rescale)
    n = len(B)
    kappa = -1.5 if rescale else -1.0
    pair = lambda X, Y: kappa * np.trace(X @ Y)
    G = np.array([[pair(X, Y) for Y in B] for X in B])
    Ginv = np.linalg.inv(G)
    dual = [sum(Ginv[a, b] * B[b] for b in range(n)) for a in range(n)]   # (t_a, t^b) = delta

    def coords(X):
        return np.array([pair(X, d) for d in dual])
    br = lambda X, Y: X @ Y - Y @ X

    def adm(Y):  # matrix of ad_Y in basis B: column b = coords of [Y, t_b]
        return np.array([coords(br(Y, B[b])) for b in range(n)]).T
    rnd = lambda: sum((rng.normal() * B[a]) for a in range(n))
    k = np.trace(adm(B[0]) @ adm(B[1])) / pair(B[0], B[1]) if abs(pair(B[0], B[1])) > 1e-9 else \
        np.trace(adm(B[0]) @ adm(B[0])) / pair(B[0], B[0])
    Y1, Y2, Y3, Z = rnd(), rnd(), rnd(), rnd()
    out = {}
    # C_w contracted with the pairing: sum_A (t_A, w t^A) = tr(w)
    w2 = adm(Y1) @ adm(Y2)
    # direct: sum_A (t_A, w t^A) with w t^A expressed through the matrix of w on coordinates
    def contract_pair(W):
        s = 0
        for a in range(n):
            wt = sum(W[c, :] @ coords(dual[a]) * B[c] for c in range(n))
            s += pair(B[a], wt)
        return s

    def contract_rot(W):
        s = 0
        for a in range(n):
            wt = sum(W[c, :] @ coords(dual[a]) * B[c] for c in range(n))
            s += pair(Z, br(B[a], wt))
        return s
    out['pair_r2'] = bool(np.isclose(contract_pair(w2), k * pair(Y1, Y2)))
    w3 = adm(Y1) @ adm(Y2) @ adm(Y3)
    out['pair_r3'] = bool(np.isclose(contract_pair(w3), k / 2 * pair(Y1, br(Y2, Y3))))
    out['pair_r1'] = bool(np.isclose(contract_pair(adm(Y1)), 0))
    out['rot_r1'] = bool(np.isclose(contract_rot(adm(Y1)), -k * pair(Z, Y1)))
    out['rot_r2'] = bool(np.isclose(contract_rot(w2), -k / 2 * pair(Z, br(Y1, Y2))))
    # planar rule: tr(M1 X_a X_b) -> kappa^{-1} sum_A tr(M1 t_A w(t^A)); compare N-part with (-1)^r N tr(M1 Y_r..Y_1)
    M1 = rnd() @ rnd()
    full = sum(np.trace(M1 @ B[a] @ (sum(w2[c, :] @ coords(dual[a]) * B[c] for c in range(n)))) for a in range(n))
    # su(N): w(T) = [Y1,[Y2,T]] ; exact gl-formula pieces
    # careful: sum over S subsets of {1,2} with (-1)^{r-|S|} tr(Y_S->) tr(M1 Y_{S^c}<-)
    S_terms = (N * np.trace(M1 @ Y2 @ Y1) - np.trace(Y1) * np.trace(M1 @ Y2) - np.trace(Y2) * np.trace(M1 @ Y1)
               + np.trace(Y1 @ Y2) * np.trace(M1)) / kappa
    # the B basis is a basis of su(N) and (t_A,t^B) = delta for kappa tr: sum_A t_A (x) t^A = kappa^{-1}(P - 1/N)
    out['planar_adjacent'] = bool(np.isclose(full, S_terms))
    out['k_over_2N'] = complex(k * kappa / (2 * N))
    return out


def run():
    rng = np.random.default_rng(11)
    res = {}
    for N in (2, 3):
        for resc in (False, True):
            res[f'su({N}){" rescaled" if resc else ""}'] = check(N, resc, rng)
    ok = all(v for r in res.values() for kk, v in r.items() if isinstance(v, bool))
    return dict(cases=res, all_pass=ok)


if __name__ == '__main__':
    print(run())
