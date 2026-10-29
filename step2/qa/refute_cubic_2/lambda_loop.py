"""Refuter #2 (Part C): gaugino-loop triangle  cll*[c] K[lam] llA[A] lam llA[A] lam  in (c,A_V,A_V) and
(c,A_W,A_W).  Own derivation.

 D = {([c,lambda],lambda^*), 1/2(lambda,rho([A,lambda]))}_K = B(kappa[c0,lam_0], lam_2, A_2)
     (left derivative along delta lambda = kappa[c,lambda]; the antisymmetry of the numeric tensor B in the
      two lambda slots gives coefficient 1), B(Q1 X, Q2 Y, v Z) = (Q1, rho(v) Q2)(X,[Z,Y]).
 V1 multilinear = B(lam_x, lam_y, A_1).
 Contractions (lam_0,lam_x), (lam_2,lam_y) with P_{lam (x) lam}(p) = (,)^{-1} (x) i rho(p)|_{Sigma^*->Sigma} s(p);
 Koszul: c0 lam_0 lam_2 lam_x lam_y -> c0 (lam_0 lam_x)(lam_2 lam_y): sign -1.
 Gamma = - r0 int kappa(p+p0) (i)^2 s(p) s(p1-p) tr_Sigma(rho(p) rho(a2) rho(p1-p) rho(a1)) C,
 C = sum_{a,b} ([X0,T_a],[X2,T_b])(T_a,[X1,T_b]);   kappa = -int_0^eps |q|^2 e^{-t1 q^2}
   => Gamma = - r0 int dt T^-2 ... |p+p0|^2 tr(...) C  ;  /unit = 1/2.
 tr_Sigma of four Clifford factors on a 16-dim semi-spin module of so(10):
   16[(a,b)(c,d) - (a,c)(b,d) + (a,d)(b,c)]   (checked numerically below on an explicit realisation).
"""
import sys, os, random, json, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
from ghost_classes import (u, t1, t2, t3, w, s, I, cross, integrate_lines, cmul, S, rnd, TB, solve2, G1, G2)
from more_classes import G3, run, summarise

DIM = 10


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def emb(p):
    return list(p) + [0] * (DIM - len(p))


def tr4(a, b, c, d):
    return 16 * (dot(a, b) * dot(c, d) - dot(a, c) * dot(b, d) + dot(a, d) * dot(b, c))


def lam_graph(p0, p1, a1, X0, X1, p2, a2, X2):
    sp0 = [s * x for x in p0]; sp1 = [s * x for x in p1]
    p = list(u)
    q = [p[i] + sp0[i] for i in range(4)]
    pp = [sp1[i] - p[i] for i in range(4)]
    C = sum(dot(cross(X0, TB[a]), cross(X2, TB[b])) * dot(TB[a], cross(X1, TB[b])) for a in range(3) for b in range(3))
    F = sp.expand(dot(q, q) * tr4(emb(p), emb(a2), emb(pp), emb(a1)) * C)
    res = integrate_lines(F, [(t1, sp0), (t2, [0] * 4), (t3, [-x for x in sp1])], 3, orders=(0, 1, 2))
    return {k: cmul((Fr(1, 2), Fr(0)), v) for k, v in res.items()}


def check_trace():
    """explicit so(10) realisation: 9 real symmetric anticommuting 16x16 + i*1; verify the 4-trace."""
    Id = np.eye(2); X = np.array([[0, 1], [1, 0]]); Z = np.array([[1, 0], [0, -1]]); E = np.array([[0, 1], [-1, 0]])
    mats = {'I': Id, 'X': X, 'Z': Z, 'E': E}
    cands = []
    for w3 in itertools.product('IXZE', repeat=3):
        if w3.count('E') % 2 == 1:
            M = np.kron(np.kron(mats[w3[0]], mats[w3[1]]), mats[w3[2]])
            cands.append(M)
    # greedy/backtracking for 7 mutually anticommuting
    def anti(A, B):
        return np.allclose(A @ B + B @ A, 0)
    best = []
    def bt(start, cur):
        nonlocal best
        if len(cur) > len(best):
            best = list(cur)
        if len(best) == 7:
            return True
        for i in range(start, len(cands)):
            if all(anti(cands[i], c) for c in cur):
                cur.append(cands[i])
                if bt(i + 1, cur):
                    return True
                cur.pop()
        return False
    bt(0, [])
    assert len(best) == 7
    G = []
    Z8 = np.zeros((8, 8))
    for Ei in best:
        G.append(np.block([[Z8, Ei], [Ei.T, Z8]]))
    G.append(np.block([[Z8, np.eye(8)], [np.eye(8), Z8]]))
    G.append(np.block([[np.eye(8), Z8], [Z8, -np.eye(8)]]))
    sig = [g.astype(complex) for g in G] + [1j * np.eye(16)]
    sigb = [g.astype(complex) for g in G] + [-1j * np.eye(16)]
    for M in range(10):
        assert np.allclose(sig[M], sig[M].T) and np.allclose(sigb[M], sigb[M].T)
        for N in range(10):
            assert np.allclose(sigb[M] @ sig[N] + sigb[N] @ sig[M], 2 * (M == N) * np.eye(16))
    rng = np.random.default_rng(3)
    worst = 0
    for _ in range(20):
        vs = [rng.integers(-3, 4, 10) for _ in range(4)]
        S_ = lambda v: sum(v[M] * sig[M] for M in range(10))
        Sb = lambda v: sum(v[M] * sigb[M] for M in range(10))
        a, b, c, d = vs
        lhs = np.trace(Sb(a) @ S_(b) @ Sb(c) @ S_(d))
        rhs = tr4(list(a), list(b), list(c), list(d))
        worst = max(worst, abs(lhs - rhs))
    return worst


if __name__ == '__main__':
    out = {}
    out['trace_formula_max_abs_dev'] = float(check_trace())
    print('trace check', out['trace_formula_max_abs_dev'], flush=True)
    rows = run(lam_graph, 5, [G1, G2], extW=False, seed=41)
    out['lambda_triangle_(c,A_V,A_V) claimed (10/9 + 8/3 L) G1 + 4/3 G2'] = summarise(rows, 2)
    print(json.dumps(out), flush=True)
    rows = run(lam_graph, 5, [G3], extW=True, seed=43)
    out['lambda_triangle_(c,A_W,A_W) claimed (7/3 + 4 L) G3'] = summarise(rows, 1)
    print(json.dumps(out), flush=True)
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lambda_loop.json'), 'w'), indent=1)
