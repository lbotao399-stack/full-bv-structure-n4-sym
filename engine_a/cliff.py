"""10d chiral Clifford blocks (engine A).  rho(v)|_Sigma = sb[M] v^M : Sigma -> Sigma^*,
rho(v)|_Sigma* = s[M] v^M : Sigma^* -> Sigma.  Pairing (u, w) = u^T w for u in Sigma, w in Sigma^*.
Directions M = 0..3 : V (R^4); M = 4..9 : W_R^perp.  The single non-real generator is M = 9."""
import itertools
import numpy as np

X = np.array([[0, 1], [1, 0]]); Z = np.array([[1, 0], [0, -1]]); E = np.array([[0, 1], [-1, 0]]); I2 = np.eye(2, dtype=int)
PAULI = {'1': I2, 'X': X, 'Z': Z, 'E': E}

def kron(*ms):
    out = np.array([[1]])
    for m in ms:
        out = np.kron(out, m)
    return out

def _find9():
    cands = []
    for w in itertools.product('1XZE', repeat=4):
        if w.count('E') % 2 == 0 and w != ('1',) * 4:
            cands.append(kron(*[PAULI[c] for c in w]))
    def anti(a, b):
        return not np.any(a @ b + b @ a)
    best = []
    def rec(chosen, start):
        nonlocal best
        if len(chosen) > len(best):
            best = list(chosen)
        if len(best) == 9:
            return True
        for i in range(start, len(cands)):
            if all(anti(cands[i], c) for c in chosen):
                if rec(chosen + [cands[i]], i + 1):
                    return True
        return False
    rec([], 0)
    assert len(best) == 9
    return best

_g9 = _find9()
s = [m.astype(complex) for m in _g9] + [1j * np.eye(16)]
sb = [m.astype(complex) for m in _g9] + [-1j * np.eye(16)]
DIM10 = 10
VDIR = [0, 1, 2, 3]
PERP = [4, 5, 6, 7, 8, 9]

def checks():
    out = {}
    ok1 = all(np.allclose(s[M] @ sb[N] + s[N] @ sb[M], 2 * (M == N) * np.eye(16)) for M in range(10) for N in range(10))
    ok2 = all(np.allclose(sb[M] @ s[N] + sb[N] @ s[M], 2 * (M == N) * np.eye(16)) for M in range(10) for N in range(10))
    sym = all(np.allclose(s[M], s[M].T) and np.allclose(sb[M], sb[M].T) for M in range(10))
    real4 = all(np.allclose(sb[M].imag, 0) for M in VDIR)
    # [ESW Assumption 3.1] rho(Gamma(Q1,Q2))Q3 + cyclic = 0 for Q_i in Sigma, Gamma(Q1,Q2)^M = Q1^T sb[M] Q2
    rng = np.random.default_rng(1)
    Q = [rng.normal(size=16) for _ in range(3)]
    def G(a, b):
        return np.array([a @ sb[M] @ b for M in range(10)])
    def rhoS(v, q):
        return sum(v[M] * sb[M] for M in range(10)) @ q
    a31 = rhoS(G(Q[0], Q[1]), Q[2]) + rhoS(G(Q[1], Q[2]), Q[0]) + rhoS(G(Q[2], Q[0]), Q[1])
    out.update(clifford_Sigma=ok1, clifford_Sigmastar=ok2, symmetric=sym, dslash_real=real4,
               ESW_assumption_3_1=bool(np.allclose(a31, 0)))
    g5 = s[0] @ sb[1] @ s[2] @ sb[3]
    out['4d_chirality_eigs_on_Sigma'] = sorted(np.round(np.linalg.eigvals(g5).real, 6).tolist()).count(1.0)
    return out

if __name__ == '__main__':
    print(checks())
