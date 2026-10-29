"""Step 2, Part D4, engine A: sign and normalisation constants of the bracket conventions (ch01, ch03), measured EXACTLY
in the plane-wave mode model of s2_ga (two modes p, -p; su(2) or gl(3) colour).

Test operator: a constant even operator X on the fields that preserves the pairing,
    X = x_A on the slot of A, -x_A on A^*,  x_l on lambda, -x_l on lambda^*,  x_c on c, -x_c on c^*,
(<X e0, e1> + <e0, X e1> = 0), with the functional built exactly as S^(1)(Delta) is built from v. (def:c7-vdot):
    S_X = int (X A, A^*) - int (X lambda, lambda^*) - int (X c, c^*).
For X = M = v. (v = 2E) this is S_M = S^(1)(Delta) (ex:c2-dil); the E-part of M is a linear operator of the same kind,
and every statement below is linear in X with no X-dependent parity, so the constants measured here hold for M.

Measured (D_X = derivation of functionals induced by X, i.e. F(e) -> d/ds F(e + s X e)):
  sigma : {S_X, F} = sigma * D_X F               (F = int(A,A), int(lambda, dslash lambda), int(dc,A^*), cubic I)
  s     : Q int(Y A|_P, A^*|_P) = s * int(Y A|_P, -(1/2) d^*d A|_P)           (Y constant on A|_P)
  c1    : Q C_P(S_X, O) = c1 * D_X O,  O = int (iota_w1 A, iota_w2 A),  P = P(0,infty) (Pnum/p^2)
  c2    : C_P(S_X, int(A|_P, A|_P)) = c2 * int( X (|p|^{-2} A|_P), A^*|_P )       (A|_P self-energy-type leg)
  c2gh  : C_P(S_X, W_gh) for a ghost two-point W_gh = int (dc, A^*|_V): the (A|_V, A^*|_V) and (c, c^*) outputs
  c2l   : C_P(S_X, int(lambda, dslash lambda))
and the derived ratio kappa := c1 / s used in the arity-one master equation (see d4_notes.txt)."""
import os, sys, json, itertools
from fractions import Fraction
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
import s2_core as C
from r2_core import Q2
import s2_ga as GAm
from grass import GA, PAR

q2 = GAm.q2


def quadP(M, name, pairs, comps):
    """restricted bilinears: 'AAs' = int(A_c, A^*_c), 'AA' = int(A_c, A_c), for components M in comps"""
    acc = {}
    G = C.G; D = C.D
    for (m1, m2) in pairs:
        if (M.moms[m1] + M.moms[m2]).any():
            continue
        for a in range(D):
            for b in range(D):
                g = G[a, b]
                if g == 0:
                    continue
                for c in comps:
                    t2 = 'As' if name == 'AAs' else 'A'
                    M.add(acc, q2(g), [M.y('A', (c, a), m1), M.y(t2, (c, b), m2)])
    return M.ga(acc)


def build(p):
    M = GAm.Model([p, -p])
    pm = [(0, 1), (1, 0)]
    pairs = [(0, 1), (1, 0)]
    return M, pm, pairs


def S_X(M, pairs, xA, xl, xc, comps=range(10)):
    out = GA()
    if xA:
        out = out + quadP(M, 'AAs', pairs, comps) * Fraction(xA)
    if xl:
        out = out - M.quad('lsl', pairs) * Fraction(xl)
    if xc:
        out = out - M.quad('csc', pairs) * Fraction(xc)
    return out


def Pkernel(M, scale_inv_p2=True):
    p = M.moms[0]
    p2 = int(p @ p)
    return M.kernel('P', 0, 1, Q2(Fraction(1, p2)) if scale_inv_p2 else Q2(1))


def ratio(A, B):
    """A = r * B exactly (GA); return r or None"""
    if not A.t and not B.t:
        return 'both zero'
    if not B.t:
        return None
    r = None
    for k in set(A.t) | set(B.t):
        a = A.t.get(k, Q2(0)); b = B.t.get(k, Q2(0))
        if b.iszero():
            if not a.iszero():
                return None
            continue
        x = a / b
        if r is None:
            r = x
        elif not (x - r).iszero():
            return None
    return r


def run():
    out = {'algebra': C.ALG}
    p = np.array([1, 2, -1, 1])
    M, pm, pairs = build(p)
    Sf = M.S_free(pairs)
    xA, xl, xc = 3, 5, 7
    SX = S_X(M, pairs, xA, xl, xc)
    # --- sigma: {S_X, F} = sigma D_X F
    tests = {}
    F_AA = M.quad('AA', pairs); tests['int(A,A)'] = (F_AA, 2 * xA)
    F_ldl = M.quad('ldl', pairs); tests['int(lambda,dslash lambda)'] = (F_ldl, 2 * xl)
    F_dcAs = M.quad('dcAs', pairs); tests['int(dc,A^*)'] = (F_dcAs, xc - xA)
    F_dAdA = M.quad('dAdA', pairs); tests['int(dA,dA)'] = (F_dAdA, 2 * xA)
    sig = {}
    for name, (F, w) in tests.items():
        br = M.bracket0(SX, F, 1, pm)
        sig[name] = str(ratio(br, F * Fraction(w)))
    out['sigma_tests'] = sig
    # --- s: Q int(Y A|_P, A^*|_P) = s int(Y A|_P, -(1/2) d^*d A|_P); the latter = -(1/2) int(dA|_P, dA|_P)
    PC = range(4, 10)
    T = quadP(M, 'AAs', pairs, PC)
    QT = M.Q_apply(T, Sf, pm)
    # int(A|_P, -(1/2) d^*d A|_P) in modes: -(1/2) p^2 (A(p), A(-p)) summed over both mode assignments
    ref = quadP(M, 'AA', pairs, PC) * Fraction(-int(p @ p), 2)
    out['s'] = str(ratio(QT, ref))
    # --- c1 with O = int (iota_w1 A, iota_w2 A), w1 = e_4, w2 = e_5 (traceless) and w1 = w2 = e_4
    Pk = Pkernel(M)
    for lab, (w1, w2) in {'traceless(e4,e5)': (4, 5), 'diag(e4,e4)': (4, 4)}.items():
        acc = {}
        for (m1, m2) in pairs:
            for a in range(C.D):
                for b in range(C.D):
                    g = C.G[a, b]
                    if g == 0:
                        continue
                    M.add(acc, q2(g), [M.y('A', (w1, a), m1), M.y('A', (w2, b), m2)])
        O = M.ga(acc)
        SXP = S_X(M, pairs, xA, 0, 0, comps=PC)
        tree = M.C(Pk, SXP, O)
        Qt = M.Q_apply(tree, Sf, pm)
        DXO = O * Fraction(2 * xA)
        out['c1_' + lab] = str(ratio(Qt, DXO))
        # also the reversed order C_P(O, S_X)
        out['c1_reversed_order_' + lab] = str(ratio(M.Q_apply(M.C(Pk, O, SXP), Sf, pm), DXO))
    # --- c2: C_P(S_X, int(A|_P,A|_P)) = c2 * int(X |p|^{-2} A|_P, A^*|_P)
    W = quadP(M, 'AA', pairs, PC)
    SXP = S_X(M, pairs, xA, 0, 0, comps=PC)
    t2 = M.C(Pk, SXP, W)
    ref2 = quadP(M, 'AAs', pairs, PC) * Fraction(xA, int(p @ p))
    out['c2_AP'] = str(ratio(t2, ref2))
    # Q of it, symmetrised (functional of A only) vs D_X W = 2 xA W
    out['Q C_P(S_X, int(A_P,A_P)) / D_X int(A_P,A_P)'] = str(ratio(M.Q_apply(t2, Sf, pm), W * Fraction(2 * xA)))
    # A|_V: C_P(S_X, int(A|_V,A|_V))
    VC = range(0, 4)
    WV = quadP(M, 'AA', pairs, VC)
    SXV = S_X(M, pairs, xA, 0, 0, comps=VC)
    out['c2_AV'] = str(ratio(M.C(Pk, SXV, WV), quadP(M, 'AAs', pairs, VC) * Fraction(xA, int(p @ p))))
    # ghost two-point W_gh = int(dc, A^*|_V); C_P(S_X, W_gh) splits into the (A_V,A^*_V) output (S_X's A^* contracted with c)
    # and the (c, c^*) output (S_X's c contracted with A^*)
    Wgh = M.quad('dcAs', pairs)
    SXA = S_X(M, pairs, xA, 0, 0, comps=VC)
    SXc = S_X(M, pairs, 0, 0, xc)
    gA = M.C(Pk, SXA, Wgh)
    gc = M.C(Pk, SXc, Wgh)
    # references: int(X (P_{A* c} p-structure) ...) -- report the explicit plane-wave coefficients instead of a ratio
    def planewave(F, legs):
        rs = GAm.rvars(len(legs))
        L = []
        for (t, pol, mo), r in zip(legs, rs):
            L.append((t, pol, mo, r if PAR_T[t] else None))
        Fs = M.subst(F, L, [0, 1])
        odd = [r for (t, pol, mo), r in zip(legs, rs) if PAR_T[t]]
        return GAm.coeff(Fs, odd) if odd else Fs.t.get((), Q2(0))
    PAR_T = C.PARITY
    rng = np.random.default_rng(3)
    aV = np.zeros((10, C.D)); aV[:4] = rng.integers(-2, 3, size=(4, C.D))
    bV = np.zeros((10, C.D)); bV[:4] = rng.integers(-2, 3, size=(4, C.D))
    gam = rng.integers(-2, 3, size=C.D).astype(float); w = rng.integers(-2, 3, size=C.D).astype(float)
    pA = (aV @ np.eye(C.D)); dotp = sum(p[m] * (aV[m] @ C.G.real @ bV[m]) for m in range(4))
    vA = planewave(gA, [('A', aV, 0), ('As', bV, 1)])
    vc = planewave(gc, [('c', gam, 0), ('cs', w, 1)])
    p2 = int(p @ p)
    # expected structures: (A,A^*) output proportional to (p.a)(p.b)/p^2 * x ... report raw numbers with the invariants
    pa = np.array([p @ aV[:4, c] for c in range(C.D)]); pb = np.array([p @ bV[:4, c] for c in range(C.D)])
    inv_papb = float(pa @ C.G.real @ pb)
    out['C_P(S_X(A_V part), int(dc,A^*)) on (A_V(p)=a, A^*_V(-p)=b)'] = {'value': str(vA), 'p.a G p.b': inv_papb,
                                                                       'p^2': p2, 'x_A': xA}
    gw = float(gam @ C.G.real @ w)
    out['C_P(S_X(c part), int(dc,A^*)) on (c(p)=gamma, c^*(-p)=w)'] = {'value': str(vc), '(gamma,w)': gw, 'p^2': p2,
                                                                      'x_c': xc}
    # lambda: C_P(S_X, int(lambda, dslash lambda))
    Wl = M.quad('ldl', pairs)
    SXl = S_X(M, pairs, 0, xl, 0)
    gl = M.C(Pk, SXl, Wl)
    out['C_P(S_X(lambda part), int(lambda,dslash lambda)) / int(lambda,lambda^*) x_l'] = str(
        ratio(gl, M.quad('lsl', pairs) * Fraction(xl)))
    return out


if __name__ == '__main__':
    res = run()
    print(json.dumps(res, indent=1))
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    json.dump(res, open(os.path.join(HERE, 'logs', 'd4_signs_%s.json' % C.ALG), 'w'), indent=1)
