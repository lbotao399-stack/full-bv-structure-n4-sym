# k_tr_defectlog.py -- the coefficient of log T in the defect graphs along K_b at a FIXED UV boundary eps
# (k_transfer.tex, Theorem thm:ktr-defectlog).
#
# Statement checked.  At one loop, as T -> infinity at fixed eps, the coefficient of log T of the defect flow along K_b
# is (3/(8 pi^2)) hbar Z_b(special-edge tree flow at eps), and by the tree-level UV-boundary identity together with
# Z_b Q = -Q Z_b (the constant modes are Q-closed), this equals -(3/(8 pi^2)) hbar Q Y(eps), with Y(eps) the zero-mode
# tree built with the regulated edges P(eps, infinity) (scalar factor e^{-eps|k|^2}/|k|^2).
#   (1) Action, two-point sectors: Y_W(eps) has one internal line; both terms of Q Y_W(eps) carry the same line momentum,
#       so Q Y_W(eps) = e^{-eps|p|^2} Q Y_W(0) = 0 identically in eps (no log T in the two-point defect at any eps).
#   (2) Scalar bilinear O = int omega (iota_w1 A, iota_w2 A), output sector (A|_Wperp, A|_Wperp): Q Y_O(eps) is
#       a1 r(p) r(q) + b1 r(q)^2 + a2 r(q) r(p) + b2 r(p)^2,  r(s) = e^{-eps |s|^2},
#       with a1 + a2 + b1 + b2 = 0 (Q Y_O(0) = 0) but a nonzero coefficient of eps: the defect of the bilinear carries
#       eps log T at fixed eps, which vanishes as eps -> 0 at fixed T.
# Realisation and conventions: cubic_core.py (su(2) orthonormal, explicit 16-dim spinors), the same functions as
# k_tr_zeromode.py (copied, so that importing does not rerun that script).
import sys
import json
import random
from fractions import Fraction as Fr
import mpmath as mp

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory')
import cubic_core as cc  # noqa: E402

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory/ktr_out/k_tr_defectlog.json'
GQ, ZERO = cc.GQ, cc.ZERO
T_ = [[Fr(1 if b == a else 0) for b in range(3)] for a in range(3)]
CAS = Fr(-2)
b = [Fr(2, 3), Fr(-1, 5), Fr(1, 4), Fr(3, 7)]


def vec10(entries):
    out = [Fr(0)]*10
    for i, x in entries.items():
        out[i] = Fr(x)
    return out


def neg(k):
    return tuple(-x for x in k)


def ksq(k):
    return sum(x*x for x in k)


def bdot(k):
    return sum(b[i]*k[i] for i in range(4))


def mode(slot, u, X, k):
    u = [x if isinstance(x, GQ) else GQ(x) for x in u]
    g = cc.new_gen(cc.SPAR[slot], 'ext')
    return g, cc.Mode(cc.gen_elem(g), slot, u, X, k)


def phi_internal(k):
    out = []
    for i in range(4, 10):
        for a in range(3):
            g = cc.new_gen(0, 'int')
            out.append((g, (i, a), cc.Mode(cc.gen_elem(g), 'A', [GQ(x) for x in vec10({i: 1})], T_[a], k, internal=True)))
    return out


def one_leg(F, exts, ints):
    val = cc.evaluate(F, exts, [m for (_, _, m) in ints], maxint=1)
    out = {}
    for (g, lab, _) in ints:
        d = cc.dleft(val, g)
        d = {m: c for m, c in d.items() if all(cc.KIND[h] != 'int' for h in m)}
        if d:
            out[lab] = d
    return out


def two_legs(F, exts, ints1, ints2):
    val = cc.evaluate(F, exts + [m for (_, _, m) in ints2], [m for (_, _, m) in ints1], maxint=2)
    out = {}
    for (g1, l1, _) in ints1:
        d1 = cc.dleft(val, g1)
        if not d1:
            continue
        for (g2, l2, _) in ints2:
            d = cc.dleft(d1, g2)
            d = {m: c for m, c in d.items() if all(cc.KIND[h] != 'int' for h in m)}
            if d:
                out[(l1, l2)] = d
    return out


def line(da, db, k):
    out = {}
    for lab, x in da.items():
        if lab in db:
            out = cc.aadd(out, cc.ascale(cc.amul(x, db[lab], maxint=0), Fr(2)/ksq(k)))
    return out


def Qimg(u, k):
    return [Fr(-1, 2)*ksq(k)*x for x in u]


# ------------------------------------------------------------------ (1) action, sector (A_Wperp, A_Wperp)
def YW_phi_val(u, X, k, us, Xs):
    tot = ZERO
    for a in range(3):
        gA, mA = mode('A', vec10({mu: b[mu] for mu in range(4)}), T_[a], (0, 0, 0, 0))
        gc, mc = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))
        g1, m1 = mode('A', u, X, k)
        g2, m2 = mode('As', us, Xs, neg(k))
        da = one_leg(cc.I3, [mA, m1], phi_internal(neg(k)))
        dg = one_leg(cc.I3, [mc, m2], phi_internal(k))
        d = cc.dleft(cc.dleft(line(da, dg, k), gc), gA)
        tot = tot + cc.ext_coeff(d, [g1, g2])
    return tot


res = {}
p = (Fr(1, 2), Fr(-1, 3), Fr(2, 5), Fr(1, 7))
q = (Fr(-2, 7), Fr(3, 5), Fr(-1, 4), Fr(1, 3))
uW = vec10({4: Fr(1), 6: Fr(-2, 3), 9: Fr(1, 2)})
uW2 = vec10({4: Fr(1, 3), 5: Fr(1, 3), 7: Fr(2), 9: Fr(-1, 2)})
X1 = [Fr(1, 2), Fr(-1), Fr(1, 3)]
X3 = [Fr(1), Fr(1, 5), Fr(-1, 2)]
t1 = YW_phi_val(uW, X1, p, Qimg(uW2, neg(p)), X3)
t2 = YW_phi_val(uW2, X3, neg(p), Qimg(uW, p), X1)
res['1_action_two_point_(A_Wperp,A_Wperp)'] = {
    'term_1_Y_W(phi1(p); Q phi2(-p))': str(t1), 'term_2_Y_W(phi2(-p); Q phi1(p))': str(t2),
    'internal_line_momentum_both_terms': '|p|^2 = %s (term 1), %s (term 2)' % (ksq(p), ksq(neg(p))),
    'Q_Y_W(eps)': 'e^{-eps|p|^2} (term_1 + term_2) = %s' % str(t1 + t2),
    'identically_zero_in_eps': (t1 + t2).iszero() and ksq(p) == ksq(neg(p))}

# ------------------------------------------------------------------ (2) scalar bilinear, sector (A_Wperp, A_Wperp)
w1 = vec10({4: Fr(1), 5: Fr(1, 2), 7: Fr(-1)})
w2 = vec10({5: Fr(2), 6: Fr(-1, 3), 9: Fr(1)})


def O_two(l1, l2, W1, W2):
    (i, a), (j, c) = l1, l2
    return (W1[i]*W2[j] + W1[j]*W2[i])*(1 if a == c else 0)


def O_one(u, X, l, W1, W2):
    (i, a) = l
    return (sum(W1[k]*u[k] for k in range(4, 10))*W2[i] + W1[i]*sum(W2[k]*u[k] for k in range(4, 10)))*X[a]


def YO_parts(phi, phis, W1, W2):
    """(alpha, beta): the two zero-mode trees of O on (phi(k1), phi^*(k2)), unregulated edges 2/|k|^2.
    alpha: lines of momenta k1, k2;  beta: two lines of momentum k2 (k_transfer.tex, Prop prop:ktr-values (iii))."""
    (u, X, k1), (us, Xs, k2) = phi, phis
    al, be = ZERO, ZERO
    for a in range(3):
        gA, mA = mode('A', vec10({mu: b[mu] for mu in range(4)}), T_[a], (0, 0, 0, 0))
        gc, mc = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))
        g1, m1 = mode('A', u, X, k1)
        g2, m2 = mode('As', us, Xs, k2)
        da = one_leg(cc.I3, [mA, m1], phi_internal(neg(k1)))
        dg = one_leg(cc.I3, [mc, m2], phi_internal(neg(k2)))
        val = {}
        for la, xa in da.items():
            for lg, xg in dg.items():
                o = O_two(la, lg, W1, W2)
                if o:
                    val = cc.aadd(val, cc.ascale(cc.amul(xa, xg, maxint=0), o*Fr(2)/ksq(k1)*Fr(2)/ksq(k2)))
        al = al + cc.ext_coeff(cc.dleft(cc.dleft(val, gc), gA), [g1, g2])
        kl = neg(k2)
        dg2 = one_leg(cc.I3, [mc, m2], phi_internal(kl))
        dm = two_legs(cc.I3, [mA], phi_internal(neg(k2)), phi_internal(neg(kl)))
        val = {}
        for (l1, l2), xm in dm.items():
            if l2 not in dg2:
                continue
            o = O_one(u, X, l1, W1, W2)
            if o:
                val = cc.aadd(val, cc.ascale(cc.amul(xm, dg2[l2], maxint=0), o*Fr(2)/ksq(k2)*Fr(2)/ksq(kl)))
        val = cc.amul(cc.gen_elem(g1), val, maxint=0)
        be = be + cc.ext_coeff(cc.dleft(cc.dleft(val, gc), gA), [g1, g2])
    return al, be


def gq_sym(x):
    return complex(float(x.re), float(x.im))


def Sval(W1, W2, U1, U2):
    d = lambda A_, B_: sum(A_[k]*B_[k] for k in range(4, 10))
    return d(W1, U1)*d(W2, U2) + d(W2, U1)*d(W1, U2)


def bilinear_case(W1, W2, tag, P, Qm, U1, U2, Xa, Xb):
    a1, b1 = YO_parts((U1, Xa, P), (Qimg(U2, Qm), Xb, Qm), W1, W2)      # phi(p), Q phi'(q)
    a2, b2 = YO_parts((U2, Xb, Qm), (Qimg(U1, P), Xa, P), W1, W2)      # phi'(q), Q phi(p)
    # closed forms (k_transfer.tex, Theorem thm:ktr-defectlog): a1 = 2 i C (b,p) S/|p|^2, a2 = 2 i C (b,q) S/|q|^2,
    # b1 = -a2, b2 = -a1, S = (w1,u)(w2,u') + (w2,u)(w1,u'), colour (X1, X3)
    S = Sval(W1, W2, U1, U2)*cc.cp(Xa, Xb)
    a1_cf = GQ(0, 2)*CAS*bdot(P)/ksq(P)*S
    a2_cf = GQ(0, 2)*CAS*bdot(Qm)/ksq(Qm)*S
    closed = {'a1_closed_form_ok': (a1 - a1_cf).iszero(), 'a2_closed_form_ok': (a2 - a2_cf).iszero(),
              'b1=-a2': (b1 + a2).iszero(), 'b2=-a1': (b2 + a1).iszero()}
    tot0 = a1 + a2 + b1 + b2
    p2, q2 = ksq(P), ksq(Qm)
    # first order in eps: -eps [ (a1+a2)(p2+q2) + 2 b1 q2 + 2 b2 p2 ]
    c1 = (a1 + a2)*(p2 + q2) + b1*(2*q2) + b2*(2*p2)
    mp.mp.dps = 30
    vals = {}
    for e in ['1e-2', '1e-4', '1e-6']:
        ee = mp.mpf(e)
        rp, rq = mp.e**(-ee*p2), mp.e**(-ee*q2)
        v = gq_sym(a1 + a2)*rp*rq + gq_sym(b1)*rq**2 + gq_sym(b2)*rp**2
        vals['eps=' + e] = {'Q_Y_O(eps)': mp.nstr(v, 12), 'Q_Y_O(eps)/eps': mp.nstr(v/ee, 12)}
    return {'tag': tag, 'a1': str(a1), 'b1': str(b1), 'a2': str(a2), 'b2': str(b2), 'closed_forms': closed,
            'Q_Y_O(0)=a1+a2+b1+b2': str(tot0), 'Q_Y_O(0)_zero': tot0.iszero(),
            'Q_Y_O(eps)=(r_q-r_p)(a1 r_p - a2 r_q), first order eps(|p|^2-|q|^2)(a1-a2)': str((a1 - a2)*(ksq(P) - ksq(Qm))),
            'coefficient_of_-eps_in_Q_Y_O(eps)': str(c1), 'eps_coefficient_nonzero': not c1.iszero(),
            'predicted_Q_Y_O(eps)/eps_at_eps->0': str(-c1), 'numerics': vals}


cases = []
cases.append(bilinear_case(w1, w2, 'generic (w1,w2) (traceless + trace)', p, q, uW, uW2, X1, X3))
# traceless: w1 orthogonal to w2 within W_R^perp, both nonzero
w1t = vec10({4: Fr(1), 5: Fr(1)})
w2t = vec10({4: Fr(1), 5: Fr(-1), 6: Fr(2)})
cases.append(bilinear_case(w1t, w2t, 'traceless (w1,w2)=0', p, q, uW, uW2, X1, X3))
# trace: sum_w (iota_w A, iota_w A) = sum over an orthonormal basis; by linearity sum of the diagonal cases w1 = w2 = e_I
def trace_case(P, Qm, U1, U2, Xa, Xb):
    A1 = B1 = A2 = B2 = ZERO
    for I in range(4, 10):
        e = vec10({I: 1})
        a1, b1 = YO_parts((U1, Xa, P), (Qimg(U2, Qm), Xb, Qm), e, e)
        a2, b2 = YO_parts((U2, Xb, Qm), (Qimg(U1, P), Xa, P), e, e)
        A1, B1, A2, B2 = A1 + a1, B1 + b1, A2 + a2, B2 + b2
    p2, q2 = ksq(P), ksq(Qm)
    c1 = (A1 + A2)*(p2 + q2) + B1*(2*q2) + B2*(2*p2)
    return {'tag': 'trace sum_w (iota_w A, iota_w A) (with the 1/2 of w1=w2 not divided out)', 'a1': str(A1), 'b1': str(B1),
            'a2': str(A2), 'b2': str(B2), 'Q_Y_O(0)': str(A1 + A2 + B1 + B2), 'Q_Y_O(0)_zero': (A1 + A2 + B1 + B2).iszero(),
            'coefficient_of_-eps': str(c1), 'eps_coefficient_nonzero': not c1.iszero()}


cases.append(trace_case(p, q, uW, uW2, X1, X3))
# a second kinematic point for the generic case
random.seed(5)
p_b = tuple(Fr(random.randint(-5, 5), random.randint(1, 4)) for _ in range(4))
q_b = tuple(Fr(random.randint(-5, 5), random.randint(1, 4)) for _ in range(4))
cases.append(bilinear_case(w1, w2, 'generic, second kinematic point', p_b, q_b, uW, uW2, X1, X3))
res['2_scalar_bilinear_(A_Wperp,A_Wperp)'] = cases
res['conclusion'] = ('Two-point sectors of the action: Q Y_W(eps) = 0 identically, so the two-point defect along K_b has no '
                     'log T at any eps.  Scalar bilinear: Q Y_O(0) = 0 (gauge invariance) and Q Y_O(eps) = O(eps) with a '
                     'nonzero eps coefficient, so the defect of the bilinear grows like eps log T at fixed eps and the '
                     'limit eps -> 0 at fixed T removes it.')
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
