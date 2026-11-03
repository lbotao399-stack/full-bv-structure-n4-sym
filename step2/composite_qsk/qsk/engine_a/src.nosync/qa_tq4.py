"""Engine A, part Q, task TQ4: the scalar bilinears along xi in S_+ (x) W (and S_- (x) W^* by the chirality swap).
O = int omega (iota_w1 A, iota_w2 A).  Two-body part = Upsilon_xi(A_perp, A_perp) on the leg sets of Table
tab:qs-fieldpairs (lamM; AV lamM; AW lamP; c lamsP), fitted onto the four functionals of Cor c7-Q:
  F1 = int omega [(iota_w1 A, (rho(w2) xi, lam)) + (iota_w2 A, (rho(w1) xi, lam))]      (classical type)
  F2 = (w1,w2) int omega (rho(A|_perp) xi, lam),  F3 = (w1,w2) int omega (rho(A|_V) xi, lam),
  F4 = (w1,w2) int omega (c, (xi, lam^*)).
Classical transformation (Cor c7-Q; Lemma qs-arrows with the sign corrected in the revision, A -> +Gamma(xi,lam)): {S^(1)(xi), O} = +F1."""
from fractions import Fraction
import numpy as np
import qa_s1 as S1, qa_core as QC, qa_util as U

SB = S1.SB


def rho_w_xi(w, xi):
    return sum(w[i] * (SB[4 + i] @ xi) for i in range(6))


def F_on_legs(kind, w1, w2, xi, legs):
    """the functional F_kind on the (unordered) legs, as an exact Gaussian rational"""
    tw = complex(np.dot(w1, w2))
    tot = 0j
    d = {l[0]: l for l in legs}
    if kind == 1:
        a = d['AW'][1]; l = d['lamP'][1]
        tot = np.dot(w1, a) * np.dot(rho_w_xi(w2, xi), l) + np.dot(w2, a) * np.dot(rho_w_xi(w1, xi), l)
    elif kind == 2:
        a = d['AW'][1]; l = d['lamP'][1]
        tot = tw * sum(a[i] * np.dot(SB[4 + i] @ xi, l) for i in range(6))
    elif kind == 3:
        a = d['AV'][1]; l = d['lamM'][1]
        tot = tw * sum(a[m] * np.dot(SB[m] @ xi, l) for m in range(4))
    elif kind == 4:
        c = d['c'][1][0]; ls = d['lamsP'][1]
        tot = tw * c * np.dot(xi, ls)
    return U.gq(tot)


def run(xi, npts=3, seed=41, log=print, chir='P'):
    """chir 'P': xi in S_+ (x) W (letters/legs as stated); 'M': xi in S_- (x) W^* (lamP <-> lamM, lamsP <-> lamsM)"""
    R = U.Rng(seed)
    gx = QC.odd_gen(xi)
    sw = (lambda t: t) if chir == 'P' else (lambda t: {'lamP': 'lamM', 'lamM': 'lamP', 'lamsP': 'lamsM', 'lamsM': 'lamsP'}.get(t, t))
    res = {}
    # traceless points (w1.w2 = 0) for a1; general points for a2..a4; extra points for the residual check
    def pts(traceless):
        while True:
            w1 = R.gv(6); w2 = R.gv(6)
            if traceless:
                # make w2 orthogonal to w1 (bilinear dot): w2 -> (w1.w1) w2 - (w1.w2) w1
                w2 = np.dot(w1, w1) * w2 - np.dot(w1, w2) * w1
                if np.dot(w1, w1) == 0:
                    continue
            else:
                if np.dot(w1, w2) == 0:
                    continue
            return w1, w2
    fits = {}
    for ls, kinds in ((['AW', 'lamP'], (1, 2)), (['AV', 'lamM'], (3,)), (['c', 'lamsP'], (4,))):
        rows = []
        for ip in range(npts + len(kinds)):
            w1, w2 = pts(traceless=(1 in kinds and ip == 0))
            legs = [(sw(t), R.vec(sw(t)), R.mom()) for t in ls]
            Xa = ('AW', w1, None); Xb = ('AW', w2, None)
            v = U.both_orders(Xa, Xb, legs, gx)['fin']
            legs_std = [(t, l[1], l[2]) for t, l in zip(ls, legs)]
            Fs = [F_on_legs(k, w1, w2, xi if chir == 'P' else xi, legs_std) for k in kinds]
            rows.append((v, Fs))
        # solve for the coefficients (exact): use the first len(kinds) independent rows
        coef = solve(rows[:len(kinds)], len(kinds))
        resid = []
        for v, Fs in rows[len(kinds):]:
            pred = U.vzero()
            for c, F in zip(coef, Fs):
                pred = U.vadd(pred, U.vscale(c, F))
            resid.append(U.viszero(U.vsub(v, pred)))
        for k, c in zip(kinds, coef):
            fits['a%d' % k] = c
        fits['resid_' + '_'.join(ls)] = all(resid)
        log('  TQ4 leg set %s: %s  residual-free on %d extra points: %s' % (ls, ['a%d = %s' % (k, U.vstr(c)) for k, c in zip(kinds, coef)], len(resid), all(resid)))
    return fits


def solve(rows, n):
    """rows: [(value, [F_1..F_n] Gaussian rationals)], n unknowns (values) -- exact Gaussian elimination per basis slot"""
    import sympy as sp
    A = sp.Matrix(n, n, lambda i, j: sp.Rational(rows[i][1][j][0]) + sp.I * sp.Rational(rows[i][1][j][1]))
    Ainv = A.inv()
    out = []
    for j in range(n):
        acc = U.vzero()
        for i in range(n):
            z = sp.nsimplify(Ainv[j, i])
            zr, zi = sp.re(z), sp.im(z)
            acc = U.vadd(acc, U.vscale(rows[i][0], (Fraction(int(sp.fraction(zr)[0]), int(sp.fraction(zr)[1])), Fraction(int(sp.fraction(zi)[0]), int(sp.fraction(zi)[1])))))
        out.append(acc)
    return out
