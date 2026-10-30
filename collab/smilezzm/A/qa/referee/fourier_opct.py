"""Referee check (4) and graph-combinatorics checks.

(4a) Fourier lemma at x != 0 by an independent route: 4d Hankel transform
       FT[f](x) = 1/(4 pi^2 |x|) int_0^inf q^2 f(q) J_1(q|x|) dq
     with EXPONENTIAL (Abel) damping e^{-eta q} (the note used Gaussian damping), eta -> 0 by polynomial
     extrapolation (the damped integrals are analytic in eta near 0), at |x| = 1 and |x| = 2.
(4b) Operator counterterm of ch05 thm:c5-Oeps / prop:c5-table, graphs (e) and (f), from the referee's own vertex
     functions (feynman_rules.py: polarisation of -S_E) with explicit su(2), su(3) structure constants, random complex
     six-vectors a, b; log coefficient int d^4l/(2pi)^4 |l|^2 D(l)^3 -> -(1/16pi^2) log eps (and B -> same).
(4c) Exchange (EX) and quartic (QU) classes: vertex/flavour/colour factors and symmetry factors from the vertex
     functions of O and of -S_E, relative to the tree, at random momenta.
(4d) Tree Gram matrix of the 21 bilinears in the basis e_I e_J (I <= J): eigenvalues.
Usage: python3 fourier_opct.py
"""
import itertools
import json
import os
import sys

import mpmath as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from feynman_rules import su_structure, V3, V4, unit_field, scalar_field, vec, V3_tensor, V4_tensor  # noqa: E402

mp.mp.dps = 22
PART = sys.argv[1] if len(sys.argv) > 1 else 'all'
out = {}
rng = np.random.default_rng(2026)

# ---------------- (4a) Fourier lemma ----------------------------------------
g, L2 = mp.euler, mp.log(2)
entries = {
    'log q^2': (lambda q: mp.log(q * q), lambda x: -1 / (mp.pi**2 * x**4)),
    'log^2 q^2': (lambda q: mp.log(q * q)**2, lambda x: 2 / (mp.pi**2 * x**4) * (mp.log(x * x) + 2 * g - 1 - 2 * L2)),
    '1/q^2': (lambda q: 1 / (q * q), lambda x: 1 / (4 * mp.pi**2 * x**2)),
    'log q^2/q^2': (lambda q: mp.log(q * q) / (q * q), lambda x: -(mp.log(x * x) + 2 * g - 2 * L2) / (4 * mp.pi**2 * x**2)),
}
ETAS = [mp.mpf(e) for e in ('0.4', '0.32', '0.25', '0.2', '0.16', '0.125', '0.1', '0.08')]


def damped(fq, x, eta):
    zeros = [mp.besseljzero(1, n) / x for n in range(1, 4000, 3)]
    pts = [mp.mpf(0)] + [z for z in zeros if eta * z < 64]
    tot = mp.mpf(0)
    for a, b in zip(pts[:-1], pts[1:]):
        tot += mp.quad(lambda q: q * q * fq(q) * mp.besselj(1, q * x) * mp.e**(-eta * q), [a, b])
    return tot / (4 * mp.pi**2 * x)


four = {}
for name, (fq, claim) in (entries.items() if PART in ('all', 'fourier', 'fourier3') else []):
    for x in ([mp.mpf(3)] if PART == 'fourier3' else [mp.mpf(1), mp.mpf(2)]):
        vals = [damped(fq, x, e) for e in ETAS]
        # Neville (polynomial) extrapolation to eta = 0, and with one point less as a stability estimate
        def neville(xs, ys, x0=0):
            P = list(ys)
            n = len(xs)
            for k in range(1, n):
                for i in range(n - k):
                    P[i] = ((x0 - xs[i + k]) * P[i] + (xs[i] - x0) * P[i + 1]) / (xs[i] - xs[i + k])
            return P[0]
        e_all = neville(ETAS, vals)
        e_less = neville(ETAS[1:], vals[1:])
        c = claim(x)
        four['%s @|x|=%s' % (name, mp.nstr(x, 2))] = dict(extrapolated=mp.nstr(e_all, 14), extrap_one_point_less=mp.nstr(e_less, 14),
                                                          claimed=mp.nstr(c, 14), rel_diff=mp.nstr((e_all - c) / c, 3))
        print(name, x, four['%s @|x|=%s' % (name, mp.nstr(x, 2))], flush=True)
out['fourier'] = four

# ---------------- (4b) operator counterterm G2 -------------------------------
def O_vertex(X, Y, a, b):
    """vertex function of O_ab = 2 hbar int a_I b_J A^c_I A^c_J (hbar = 1), symmetrised in its two legs."""
    aX, bX = X[:, 4:] @ a, X[:, 4:] @ b
    aY, bY = Y[:, 4:] @ a, Y[:, 4:] @ b
    return 2 * (np.dot(aX, bY) + np.dot(aY, bX))


G2 = {}
comb = {}
for N in ((2, 3) if PART in ('all', 'algebra') else ()):
    f = su_structure(N)
    dimG = f.shape[0]
    D = dimG * 10
    sidx = [c * 10 + M for c in range(dimG) for M in range(4, 10)]       # scalar components in the flattened index
    # validation of the fast tensors against the polarisation definitions (random fields, random momenta)
    k0, k1 = vec(rng.normal(size=4)), vec(rng.normal(size=4)); k2 = -(k0 + k1)
    As = [rng.normal(size=(dimG, 10)) + 1j * rng.normal(size=(dimG, 10)) for _ in range(4)]
    T3 = V3_tensor([k0, k1, k2], f)
    slow3 = V3([(As[0], k0), (As[1], k1), (As[2], k2)], f)
    fast3 = np.einsum('xyz,x,y,z->', T3, *(A.reshape(-1) for A in As[:3]))
    Ss = [A.copy() for A in As]
    for A in Ss:
        A[:, :4] = 0
    T4 = V4_tensor(f, list(range(4, 10)))
    slow4 = V4(Ss, f)
    fast4 = np.einsum('wxyz,w,x,y,z->', T4, *(A[:, 4:].reshape(-1) for A in Ss))
    valid = dict(V3=float(abs(slow3 - fast3)), V4_scalar=float(abs(slow4 - fast4)))
    a = rng.normal(size=6) + 1j * rng.normal(size=6)
    b = rng.normal(size=6) + 1j * rng.normal(size=6)
    ab = np.dot(a, b)
    ns = dimG * 6

    def Omat(a, b):
        """vertex function of O_ab on two scalar legs, flattened (colour, flavour)."""
        m = 2 * (np.outer(a, b) + np.outer(b, a))
        return np.kron(np.eye(dimG), m)
    Ov = Omat(a, b)
    Kv = np.kron(np.eye(dimG), 4 * np.eye(6))
    # (e) bubble: (1/2) sum_{X,Y} O(X,Y) V4(X,Y,E1,E2)  -> matrix in (E1,E2)
    bub = 0.5 * np.einsum('xy,xyuv->uv', Ov, T4)
    # (f) triangle at zero external momenta, l averaged over the 4 unit directions; X,Y scalar, Z all components
    tri = 0
    for mu in range(4):
        lv = vec(np.eye(4)[mu])
        Va = V3_tensor([lv, 0 * lv, -lv], f)[np.ix_(sidx, sidx, range(D))]     # (X, E1, Z)
        Vb = V3_tensor([-lv, 0 * lv, lv], f)[np.ix_(sidx, sidx, range(D))]     # (Y, E2, Z)
        tri = tri + np.einsum('xy,xuz,yvz->uv', Ov, Va, Vb)
    tri = tri / 4.0
    Mfit = np.array([Ov.reshape(-1), (ab * Kv).reshape(-1)]).T
    solE, *_ = np.linalg.lstsq(Mfit, bub.reshape(-1), rcond=None)
    solF, *_ = np.linalg.lstsq(Mfit, tri.reshape(-1), rcond=None)
    G2['su%d' % N] = dict(validation=valid,
                          bubble_e=dict(O=str(complex(2 * solE[0] / N)), abK=str(complex(2 * solE[1] / N)),
                                        max_residual=float(np.abs(Mfit @ solE - bub.reshape(-1)).max())),
                          triangle_f=dict(O=str(complex(2 * solF[0] / N)), abK=str(complex(2 * solF[1] / N)),
                                          max_residual=float(np.abs(Mfit @ solF - tri.reshape(-1)).max())),
                          claimed=dict(bubble_e='2 O - 2 (a,b) K', triangle_f='2 O', total='4 O - 2 (a,b) K'))
    print('G2', N, G2['su%d' % N], flush=True)
    # ---- (4c) EX and QU relative to the tree ----
    a1, b1, a2, b2 = [rng.normal(size=6) + 1j * rng.normal(size=6) for _ in range(4)]
    S12 = np.dot(a1, a2) * np.dot(b1, b2) + np.dot(a1, b2) * np.dot(b1, a2)
    T12 = np.dot(a1, b1) * np.dot(a2, b2)
    O1, O2 = Omat(a1, b1), Omat(a2, b2)
    tree = 0.5 * np.sum(O1 * O2)                                   # |Aut| = 2
    qu = 0.25 * np.einsum('xy,zw,xyzw->', O1, O2, T4)              # |Aut| = 2 x 2
    k, kp, q = (vec(rng.normal(size=4)) for _ in range(3))
    VA = V3_tensor([k, -kp, kp - k], f)[np.ix_(sidx, sidx, range(D))]
    VB = V3_tensor([q - k, -(q - kp), k - kp], f)[np.ix_(sidx, sidx, range(D))]
    ex = 0.5 * np.einsum('ab,cd,acz,bdz->', O1, O2, VA, VB)        # |Aut| = 2
    num = np.dot(k + kp, 2 * q - k - kp)
    comb['su%d' % N] = dict(tree_over_4dimS12=str(complex(tree / (4 * dimG * S12))),
                            QU_over_tree=str(complex(qu / tree)), QU_claimed=str(complex(-N * (2 * T12 - S12) / S12)),
                            EX_over_tree_over_num=str(complex(ex / tree / num)), EX_claimed=-N)
    print('comb', N, comb['su%d' % N], flush=True)
out['G2'] = G2
out['combinatorics'] = comb

# ---------------- (4d) tree Gram matrix in the basis e_I e_J ------------------
pairs = [(I, J) for I in range(6) for J in range(I, 6)]
E = np.eye(6)
Gt = np.array([[np.dot(E[I1], E[I2]) * np.dot(E[J1], E[J2]) + np.dot(E[I1], E[J2]) * np.dot(E[J1], E[I2])
                for (I2, J2) in pairs] for (I1, J1) in pairs])
ev = sorted(np.linalg.eigvalsh(Gt).round(12).tolist())
Kvec = np.array([1.0 if I == J else 0.0 for (I, J) in pairs])
out['gram_tree'] = dict(basis='O_{e_I e_J}, I<=J (21)', eigenvalues_in_units_of_4hbar_dimg_over_16pi4=ev,
                        note_claim='1 (twenty times), 2 (once, on K/sqrt6)',
                        K_norm=float(Kvec @ Gt @ Kvec), Kover_sqrt6_norm=float(Kvec @ Gt @ Kvec / 6))
print(out['gram_tree'], flush=True)
with open(os.path.join(HERE, 'fourier_opct_%s.json' % PART), 'w') as fh:
    json.dump(out, fh, indent=1, default=str)
