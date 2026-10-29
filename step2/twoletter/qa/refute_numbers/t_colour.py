"""Planar colour reading (Lemma lem:tl-colour (ii),(iii)) by explicit gl(N)/su(N) matrices: replacing X_a (x) X_b in
tr(M1 X_a X_b) by C_w = sum_A T_A (x) w(T^A), w = ad_{Y_r} ... ad_{Y_1} (Y_1 nearest a, engine/referee convention):
leading N-term must be (-1)^r N tr(M1 Y_1 ... Y_r)."""
import numpy as np, json
rng = np.random.default_rng(3)
out = {}
for N in (3, 4, 5):
    # orthonormal basis of gl(N) for tr: E_ij, dual E_ji
    basis = [(np.eye(N)[:, [i]] @ np.eye(N)[[j], :], np.eye(N)[:, [j]] @ np.eye(N)[[i], :]) for i in range(N) for j in range(N)]
    M1 = rng.normal(size=(N, N)); Ys = [rng.normal(size=(N, N)) for _ in range(3)]
    for r in (1, 2, 3):
        Y = Ys[:r]
        def w(T):
            for y in Y:            # apply ad_{Y_1} first, ad_{Y_r} last
                T = y @ T - T @ y
            return T
        lhs = sum(np.trace(M1 @ Ta @ w(Tb)) for Ta, Tb in basis)
        prod = np.eye(N)
        for y in Y:
            prod = prod @ y
        lead = (-1) ** r * N * np.trace(M1 @ prod)
        # full gl(N) formula: sum over subsets S (Lemma, relabelled): compute exactly as check of the remainder being multitrace
        out['N=%d r=%d' % (N, r)] = dict(total=float(lhs), leading=float(lead))
# exact gl(N) check of the remainder for r = 2: lhs - lead = sum_{S != empty} ... compute directly
ok = True
for N in (3, 4):
    M1 = rng.normal(size=(N, N)); Y1, Y2 = rng.normal(size=(N, N)), rng.normal(size=(N, N))
    basis = [(np.eye(N)[:, [i]] @ np.eye(N)[[j], :], np.eye(N)[:, [j]] @ np.eye(N)[[i], :]) for i in range(N) for j in range(N)]
    w = lambda T: Y2 @ (Y1 @ T - T @ Y1) - (Y1 @ T - T @ Y1) @ Y2
    lhs = sum(np.trace(M1 @ Ta @ w(Tb)) for Ta, Tb in basis)
    # w(T) = Y2 Y1 T - Y2 T Y1 - Y1 T Y2 + T Y1 Y2 ; adjacent contraction sum_A tr(M1 T_A P T^A Q) = tr(P) tr(M1 Q)
    pred = np.trace(Y2 @ Y1) * np.trace(M1) - np.trace(Y2) * np.trace(M1 @ Y1) - np.trace(Y1) * np.trace(M1 @ Y2) + N * np.trace(M1 @ Y1 @ Y2)
    ok &= abs(lhs - pred) < 1e-9
print('gl(N) r=2 adjacent formula, leading term +N tr(M1 Y1 Y2):', ok)
json.dump(dict(r2_formula_ok=bool(ok)), open('out_colour.json', 'w'))
