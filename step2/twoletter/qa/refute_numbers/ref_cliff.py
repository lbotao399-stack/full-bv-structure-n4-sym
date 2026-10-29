"""Referee's own 10d chiral Clifford realisation (independent construction): 9 mutually anticommuting real symmetric
16x16 matrices g_1..g_9 squaring to 1, built from 4-fold Pauli products by a search in REVERSED order (so the basis
differs from engine A's), sig[M] = g_M (M = 0..8), sig[9] = i*1 ; sigb[M] = g_M, sigb[9] = -i*1.
sig[M]: Sigma -> Sigma^* (symmetric bilinear on Sigma), sigb[M]: Sigma^* -> Sigma.  Checked: Clifford relations on both
sides, symmetry, the 10d Fierz identity rho(Gamma(Q,Q))Q = 0."""
import itertools
import numpy as np

P = {'1': np.eye(2, dtype=int), 'X': np.array([[0, 1], [1, 0]]), 'Z': np.array([[1, 0], [0, -1]]),
     'E': np.array([[0, 1], [-1, 0]])}


def kron(ws):
    out = np.array([[1]])
    for w in ws:
        out = np.kron(out, P[w])
    return out


def build():
    cands = []
    for w in itertools.product('EZX1', repeat=4):          # reversed alphabet order
        if w.count('E') % 2 == 0 and w != ('1',) * 4:
            cands.append(kron(w))
    cands = cands[::-1]
    best = []

    def rec(ch, st):
        nonlocal best
        if len(ch) > len(best):
            best = list(ch)
        if len(best) == 9:
            return True
        for i in range(st, len(cands)):
            if all(not np.any(cands[i] @ c + c @ cands[i]) for c in ch):
                if rec(ch + [cands[i]], i + 1):
                    return True
        return False
    rec([], 0)
    assert len(best) == 9
    g = best
    sig = [m.astype(complex) for m in g] + [1j * np.eye(16)]
    sigb = [m.astype(complex) for m in g] + [-1j * np.eye(16)]
    return sig, sigb


def checks(sig, sigb):
    ok = {}
    ok['cliff_Sigma'] = all(np.allclose(sigb[M] @ sig[N] + sigb[N] @ sig[M], 2 * (M == N) * np.eye(16))
                            for M in range(10) for N in range(10))
    ok['cliff_Sigmastar'] = all(np.allclose(sig[M] @ sigb[N] + sig[N] @ sigb[M], 2 * (M == N) * np.eye(16))
                                for M in range(10) for N in range(10))
    ok['symmetric'] = all(np.allclose(sig[M], sig[M].T) and np.allclose(sigb[M], sigb[M].T) for M in range(10))
    rng = np.random.default_rng(7)
    Q = rng.normal(size=16) + 1j * rng.normal(size=16)
    G = np.array([Q @ sig[M] @ Q for M in range(10)])
    ok['fierz'] = bool(np.allclose(sum(G[M] * sig[M] for M in range(10)) @ Q, 0))
    return ok


SIG, SIGB = build()
if __name__ == '__main__':
    print(checks(SIG, SIGB))
