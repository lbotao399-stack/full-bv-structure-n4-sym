"""P6b (conditional, matrix level): own twist-two matrix (check_twist2.py, ZZ part of lc^natural on tr(D^aZ D^bZ)) in the
letter basis of B03c, D^kZ/k! <-> (a^dag)^k|0>, two-site cyclic states  vs  -2(H'_12 + H'_21) = -4 H'_12 on the
site-symmetric states, H'_12 from predictions.json P6b (B03c (3.14))."""
import os, json, math
import sympy as sp
import cmp_lib as L
r = json.load(open(os.path.join(L.QA, 'out.nosync', 'twist2_check.json')))['twist_two_own_assembly']['per_N']
P = json.load(open(os.path.join(L.TL, 'theory', 'predictions.json')))['P6b_sl2_two_site']['levels']
out = {}
for N in range(0, 9):
    M = sp.Matrix([[sp.nsimplify(x) for x in row] for row in r[str(N)]['matrix']])   # basis g_a = tr(d^a Z d^{N-a} Z), a <= N/2
    basis = list(range(0, N // 2 + 1))
    # g_a = a!(N-a)! tr(e_a e_{N-a});  cyclic state tr(e_a e_b) <-> f_a = e_a(x)e_b + e_b(x)e_a (a != b), e_a(x)e_a (a = b)
    S = sp.diag(*[sp.factorial(a) * sp.factorial(N - a) for a in basis])          # g = S f
    Mf = S * M * S.inv()                                                           # matrix on f basis
    hh = lambda k: sum(sp.Rational(1, j) for j in range(1, k + 1))
    H = sp.Matrix(N + 1, N + 1, lambda k, kp: hh(k) + hh(N - k) if k == kp else -sp.Rational(1, abs(k - kp)))  # B03c (3.14)
    if str(N) in P:
        assert H == sp.Matrix([[sp.nsimplify(x) for x in row] for row in P[str(N)]['matrix']])
    # -4 H'_12 on symmetric states f_a
    pred = sp.zeros(len(basis), len(basis))
    for j, a in enumerate(basis):
        vec = sp.zeros(N + 1, 1)
        vec[a] += 1; vec[N - a] += 1          # f_a = e_a(x)e_{N-a} + e_{N-a}(x)e_a  (= 2 e(x)e for a = N/2)
        img = -4 * H * vec
        for i, b in enumerate(basis):
            pred[i, j] = img[b] if b != N - b else img[b] / 2
    out[N] = dict(own=[[str(x) for x in row] for row in Mf.tolist()], pred=[[str(x) for x in row] for row in pred.tolist()],
                  equal=bool(sp.simplify(Mf - pred) == sp.zeros(*pred.shape)))
    print(N, out[N]['equal'], Mf.tolist(), pred.tolist())
json.dump(out, open(os.path.join(L.QA, 'out.nosync', 'p6b_check.json'), 'w'), indent=1)
