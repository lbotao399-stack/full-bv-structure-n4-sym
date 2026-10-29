"""(1) P6b at matrix level from the referee's own twist-two kernels: on v_k = O_k/(k!(N-k)!),
    lc^natural v_k = -4 sum_k' H'[k',k] v_k',  H'[k',k] = (h(k)+h(N-k)) delta - 1/|k-k'|  (B03c (3.14)).
(2) twist two without the O(A_n) terms of the covariant derivatives (compare report side remark).
(3) su(2) Heisenberg spectra L <= 8 from the referee-verified two-site density u[2P - K] vs predictions.json P6."""
import json, itertools
import sympy as sp
import numpy as np
T = json.load(open('out_twist2.json'))
KZZ = {eval(k): {eval(a): sp.sympify(b) for a, b in v.items()} for k, v in T['KZZ'].items()}
KAZ = {eval(k): {int(a): sp.sympify(b) for a, b in v.items()} for k, v in T['KAZ'].items()}
KZA = {eval(k): {int(a): sp.sympify(b) for a, b in v.items()} for k, v in T['KZA'].items()}
h = lambda n: sum(sp.Rational(1, j) for j in range(1, n + 1))
def binom(a, b): return sp.binomial(a, b)
def lc_nat_ZZ(N, k, with_A=True):
    """lc^natural(O_k)|_{ZZ}: dict (c,d) [c<=d] -> coefficient (units u), O_k = tr(D^k Z D^{N-k} Z)"""
    out = {}
    def add(c, d, v):
        key = (min(c, d), max(c, d)); out[key] = out.get(key, 0) + v
    # ZZ part, L = 2: both adjacencies -> full kernel value
    for (c, d), v in KZZ.get((k, N - k), {}).items():
        add(c, d, v)
    add(k, N - k, -4)
    if with_A:
        # three-letter terms: D^m Z ∋ sum_{i<=j<m} C(j,i)[A_i, Z_{m-1-i}]
        def terms(mp):
            for j in range(mp):
                for i in range(j + 1):
                    yield binom(j, i), i, mp - 1 - i
        # tr([A_i,Z_m] Z_s) = tr(A_i Z_m Z_s) - tr(Z_m A_i Z_s); planar r = 1 kernels: (-1)^1 c/2 adjacent, +c/2 wrap
        for (first, other) in [(k, N - k), (N - k, k)]:
            for c0, i, m in terms(first):
                s = other
                # tr(A_i Z_m Z_s): adjacent (A_i,Z_m) -> -K/2 ; wrap pair (A_i, Z_s) [X_1, X_L] -> +K/2 reversed
                for o, v in KAZ.get((i, m), {}).items(): add(o, s, -c0 * v / 2)
                for o, v in KAZ.get((i, s), {}).items(): add(m, o, c0 * v / 2)
                # - tr(Z_m A_i Z_s): adjacent (Z_m, A_i) -> -K/2 ; adjacent (A_i, Z_s) -> -K/2
                for o, v in KZA.get((m, i), {}).items(): add(o, s, c0 * v / 2)
                for o, v in KAZ.get((i, s), {}).items(): add(m, o, c0 * v / 2)
    return out
res = {}
for N in range(1, 9):
    ks = list(range(N // 2 + 1))
    # basis O_k <-> ZZ part tr(Z_k Z_{N-k}) (k <= N/2)
    M = sp.zeros(len(ks), len(ks))
    for r, k in enumerate(ks):
        for (c, d), v in lc_nat_ZZ(N, k).items():
            M[r, ks.index(c)] += v
    # B03c prediction on v_k = O_k/(k!(N-k)!): lc v_k = -4 sum_k' H'[k',k] v_k'
    Hp = sp.zeros(N + 1, N + 1)
    for k in range(N + 1):
        for kp in range(N + 1):
            Hp[kp, k] = (h(k) + h(N - k)) if k == kp else -sp.Rational(1, abs(k - kp))
    ok = True
    for r, k in enumerate(ks):
        # lc(O_k) = k!(N-k)! lc(v_k) = k!(N-k)! (-4) sum_k' H'[k',k] O_k'/(k'!(N-k')!)
        pred = {}
        for kp in range(N + 1):
            key = min(kp, N - kp)
            pred[key] = pred.get(key, 0) + sp.factorial(k) * sp.factorial(N - k) * (-4) * Hp[kp, k] / (sp.factorial(kp) * sp.factorial(N - kp))
        for c in ks:
            # O_c has ZZ part tr(Z_c Z_{N-c}); for c = N/2 the two orderings coincide
            ok &= sp.simplify(M[r, ks.index(c)] - pred.get(c, 0)) == 0
    # without A terms
    M0 = sp.zeros(len(ks), len(ks))
    for r, k in enumerate(ks):
        for (c, d), v in lc_nat_ZZ(N, k, with_A=False).items():
            M0[r, ks.index(c)] += v
    ev0 = sorted(sum([[x] * m for x, m in M0.eigenvals().items()], []))
    res[N] = dict(P6b_matrix_level=bool(ok), eig_without_A_terms=[str(x) for x in ev0])
    print('N =', N, ' P6b matrix level:', ok, '  eig without O(A_n) terms:', ev0)
# su(2) sector from the two-site density u[2P - K] (K = 0 on Z, X): lc^natural = -2u sum (1-P)
P = json.load(open('../../theory/predictions.json'))['P6_su2_Heisenberg']['cyclic_spectra_E']
okall = True
for L in range(2, 9):
    for Mg in range(0, L // 2 + 1):
        states = [s for s in itertools.product([0, 1], repeat=L) if sum(s) == Mg]
        # cyclic classes
        cls = {}
        for s in states:
            rep_ = min(tuple(s[i:] + s[:i]) for i in range(L))
            cls.setdefault(rep_, []).append(s)
        reps = sorted(cls)
        idx = {r: i for i, r in enumerate(reps)}
        H = sp.zeros(len(reps), len(reps))
        for r in reps:
            # apply sum_l (1 - P_{l,l+1}) to the cyclic sum of r (states of the class with multiplicity)
            for l in range(L):
                s = list(r); t = list(r); t[l], t[(l + 1) % L] = t[(l + 1) % L], t[l]
                rt = min(tuple(t[i:] + t[:i]) for i in range(L))
                H[idx[r], idx[r]] += 1
                H[idx[r], idx[rt]] -= 1
        # cyclic invariance: operators tr(...) with symmetry factors -- use the trace basis; eigenvalues are basis independent
        Hn = np.array(H.T.tolist(), dtype=float)
        ev = sorted(np.round(np.linalg.eigvals(Hn).real, 9))
        key = 'L=%d,M=%d' % (L, Mg)
        if key in P:
            pr = sorted(np.round(P[key]['roots'], 9))
            m = len(ev) == len(pr) and np.allclose(ev, pr, atol=1e-7)
            okall &= m
            if not m: print('su(2) mismatch', key, ev, pr)
print('su(2) cyclic spectra L<=8 (H from u[2P-K]) equal predictions P6:', okall)
res['su2_ok'] = bool(okall)
json.dump(res, open('out_p6b_su2.json', 'w'), indent=1)
