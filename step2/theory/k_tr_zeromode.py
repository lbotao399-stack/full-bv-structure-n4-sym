# k_tr_zeromode.py -- zero-mode trees: the coefficients of log T in the flows along K_v (arity one) and along
# (eta1, eta2) (arity two), exact, in the realisation of cubic_core.py (su(2) orthonormal, explicit 16-dim spinors).
# k_transfer.tex, Theorem thm:ktr-flowlog and Propositions prop:ktr-values, prop:ktr-consistency.
#
# Statement checked here.  At one loop, the coefficient of log T in W(P(eps,T), I_psl)^(1)(K_v) and in
# O[T]^(1)(K_v) (eps -> 0 first) is hbar * 3/(8 pi^2) * Y, with
#     Y_W = sum_a d_{A|_V = b (x) X_a} d_{c = X^a} W(P(0,infinity), I)_tree,
#     Y_O = sum_a d_{A|_V = b (x) X_a} d_{c = X^a} O[infinity]_tree,
# (left derivatives along the constant modes; {X_a} a basis of g, {X^a} the dual basis for (-,-)); the coefficient of
# log T at arity two along (eta1,eta2) has the same form with the modes (c = X^a, c = (v_1, x) X_a).
# Values are reported as the coefficient of kappa_1 kappa_2 when the functional is evaluated on
# kappa_1 * (mode 1) + kappa_2 * (mode 2), (2 pi)^4 delta stripped.
# Checks: (1) the two-vertex tree of W by cubic_core.tree against the hand contraction used for three-vertex trees;
# (2) Y_W in the sectors (A|_Wperp, A^*|_Wperp), (A|_V, A^*|_V), (lambda, lambda^*), (c, c^*): nonzero, with closed forms;
# (3) Q Y_W = 0 and Q Y_O = 0 (equivariant consistency at large T);
# (4) Y_W = sigma * Q Z with Z the arity-two zero-mode tree (v_1 = b), sigma exact;
# (5) Y_O for a scalar bilinear: nonzero (two independent momentum structures);
# (6) the defect counterpart (special edge in place of the internal line) in the sector (A|_Wperp, A|_Wperp):
#     proportional to e^{-eps|p|^2} - 1 and zero after the sum over the two legs.
import sys
import json
import itertools
import random
from fractions import Fraction as Fr

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory')
import cubic_core as cc  # noqa: E402

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory/ktr_out/k_tr_zeromode.json'
GQ, ZERO, ONE, IU = cc.GQ, cc.ZERO, cc.ONE, cc.IU
res = {}

T_ = [[Fr(1 if b == a else 0) for b in range(3)] for a in range(3)]      # orthonormal basis of su(2), X^a = X_a
CAS = Fr(-2)                                                              # tr_g(ad_X ad_Y)/(X,Y) for this basis
b = [Fr(2, 3), Fr(-1, 5), Fr(1, 4), Fr(3, 7)]                             # v_1 = b
p = (Fr(1, 2), Fr(-1, 3), Fr(2, 5), Fr(1, 7))
q = (Fr(-2, 7), Fr(3, 5), Fr(-1, 4), Fr(1, 3))


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


def mode(slot, u, X, k, parity=None):
    u = [x if isinstance(x, GQ) else GQ(x) for x in u]
    g = cc.new_gen(cc.SPAR[slot] if parity is None else parity, 'ext')
    return g, cc.Mode(cc.gen_elem(g), slot, u, X, k)


def coeff(val, gens):
    return cc.ext_coeff(val, gens)


# ------------------------------------------------------------------ (1)+(2) two-vertex trees by cubic_core.tree
def YW_sector(make_ext):
    """sum_a d_A d_c of (1/2) tree(I3, I3) on the external modes from make_ext(); returns GQ coefficient."""
    tot = ZERO
    for a in range(3):
        gA, mA = mode('A', vec10({mu: b[mu] for mu in range(4)}), T_[a], (0, 0, 0, 0))
        gc, mc = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))
        exts, gens, kset = make_ext()
        cache = cc.make_internal(kset)
        dPt = cc.dP_terms(kset, cache[1])
        val = cc.tree(cc.I3, 0, cc.I3, [mA, mc] + exts, kset, cache, dPt)
        val = cc.ascale(val, Fr(1, 2))
        d = cc.dleft(cc.dleft(val, gc), gA)
        tot = tot + coeff(d, gens)
    return tot


uW = vec10({4: Fr(1), 6: Fr(-2, 3), 9: Fr(1, 2)})
uWs = vec10({4: Fr(3, 4), 6: Fr(1, 5), 8: Fr(-1)})
X1 = [Fr(1, 2), Fr(-1), Fr(1, 3)]
X2 = [Fr(-2, 3), Fr(1, 4), Fr(1)]


def ext_scalar():
    g1, m1 = mode('A', uW, X1, p)
    g2, m2 = mode('As', uWs, X2, neg(p))
    return [m1, m2], [g1, g2], [p, neg(p)]


YW_phi = YW_sector(ext_scalar)
dotWW = sum(uW[i]*uWs[i] for i in range(4, 10))
pred_phi = GQ(0, -2)*bdot(p)/ksq(p)*dotWW*CAS*cc.cp(X1, X2)          # hand formula, k_transfer.tex (cubic vertex -1/2 (d phi,[A,phi]))
res['YW_(A_Wperp,A*_Wperp)'] = {'value': str(YW_phi), 'hand_formula_-2i(b.p)/|p|^2 (u,u*) tr(adad)/(,) (X1,X2)': str(pred_phi),
                                'agree': (YW_phi - pred_phi).iszero()}

# vector sector (A|_V, A^*|_V)
uV = vec10({0: Fr(1, 3), 1: Fr(-1), 2: Fr(1, 2), 3: Fr(2)})
uVs = vec10({0: Fr(-1, 2), 1: Fr(2, 3), 2: Fr(1), 3: Fr(-1, 4)})


def ext_vector():
    g1, m1 = mode('A', uV, X1, p)
    g2, m2 = mode('As', uVs, X2, neg(p))
    return [m1, m2], [g1, g2], [p, neg(p)]


YW_vec = YW_sector(ext_vector)
# structure test: express as alpha (b.p)/p^2 (u.u*) + beta (b.u)(p.u*)/p^2 + gamma (b.u*)(p.u)/p^2 + delta (b.p)(p.u)(p.u*)/p^4 ;
# the coefficients follow from the hand contraction in k_transfer.tex; here only the value and its colour factor are recorded
def YW_vec_cfg(uv, uvs, pp_):
    def mk():
        g1, m1 = mode('A', uv, X1, pp_)
        g2, m2 = mode('As', uvs, X2, neg(pp_))
        return [m1, m2], [g1, g2], [pp_, neg(pp_)]
    return YW_sector(mk)


def vstructs(uv, uvs, pp_):
    du = lambda a_, c_: sum(a_[i]*c_[i] for i in range(4))
    pre = GQ(0, 1)*CAS*cc.cp(X1, X2)/ksq(pp_)
    return [pre*bdot(pp_)*du(uv, uvs), pre*du(b, uv)*du(pp_, uvs), pre*du(pp_, uv)*du(b, uvs),
            pre*bdot(pp_)*du(pp_, uv)*du(pp_, uvs)/ksq(pp_)]


random.seed(11)
cfgs = []
for _ in range(6):
    uv = vec10({mu: Fr(random.randint(-4, 4), random.randint(1, 3)) for mu in range(4)})
    uvs = vec10({mu: Fr(random.randint(-4, 4), random.randint(1, 3)) for mu in range(4)})
    pp_ = tuple(Fr(random.randint(-5, 5), random.randint(1, 4)) for _ in range(4))
    cfgs.append((uv, uvs, pp_))
import sympy as _sp
rows = [[_sp.nsimplify(x.re) + _sp.I*_sp.nsimplify(x.im) for x in vstructs(*c)] for c in cfgs]
rhs_ = [YW_vec_cfg(*c) for c in cfgs]
rhs_s = [_sp.nsimplify(x.re) + _sp.I*_sp.nsimplify(x.im) for x in rhs_]
coef_ = _sp.Matrix(rows[:4]).LUsolve(_sp.Matrix(rhs_s[:4]))
check_ = [_sp.simplify(sum(coef_[j]*rows[k][j] for j in range(4)) - rhs_s[k]) for k in (4, 5)]
res['YW_(A_V,A*_V)'] = {'value_at_default_config': str(YW_vec),
                        'basis_(times i tr(adad)/(,) (X1,X2)/|p|^2)': ['(b.p)(u.u*)', '(b.u)(p.u*)', '(p.u)(b.u*)', '(b.p)(p.u)(p.u*)/|p|^2'],
                        'coefficients_exact_from_4_configs': [str(c_) for c_ in coef_],
                        'verification_on_2_further_configs_residuals': [str(c_) for c_ in check_]}

# ghost sector (c, c^*)
def ext_ghost():
    g1, m1 = mode('c', [Fr(1)], X1, p)
    g2, m2 = mode('cs', [Fr(1)], X2, neg(p))
    return [m1, m2], [g1, g2], [p, neg(p)]


YW_gh = YW_sector(ext_ghost)
pred_gh_candidate = bdot(p)/ksq(p)*CAS*cc.cp(X1, X2)
res['YW_(c,c*)'] = {'value': str(YW_gh), 'ratio_to_(b.p)/|p|^2 tr(adad)/(,) (X1,X2)': str(YW_gh/GQ(pred_gh_candidate))}

# spinor sector (lambda, lambda^*)
random.seed(7)
ul = [Fr(random.randint(-3, 3), random.randint(1, 4)) for _ in range(16)]
uls = [Fr(random.randint(-3, 3), random.randint(1, 4)) for _ in range(16)]


def ext_spinor():
    g1, m1 = mode('lam', ul, X1, p)
    g2, m2 = mode('lams', uls, X2, neg(p))
    return [m1, m2], [g1, g2], [p, neg(p)]


YW_lam = YW_sector(ext_spinor)
# compare with the two spinor structures (lambda^*, rho(b) rho(p) lambda)/p^2 and (b.p)(lambda^*,lambda)/p^2
rb = cc.rho_vec(b, ul)                        # Sigma -> Sigma^*  (rho(b) on Sigma)
rpb = cc.rho_vec([IU*p[mu] for mu in range(4)], rb, bar=True)   # rho(ip) rho(b) u : Sigma -> Sigma
s1 = cc.spin_pair(rpb, uls)
rp = cc.rho_vec([IU*p[mu] for mu in range(4)], ul)
rbp = cc.rho_vec(b, rp, bar=True)
s2 = cc.spin_pair(rbp, uls)
s0 = cc.spin_pair(ul, uls)*bdot(p)
CX = CAS*cc.cp(X1, X2)/ksq(p)
res['YW_(lambda,lambda*)'] = {'value': str(YW_lam), 'closed_form_(rho(ip)rho(b)u,u*) tr(adad)/(,) (X1,X2)/|p|^2': str(s1*CX),
                              'agree_closed_form': (YW_lam - s1*CX).iszero(),
                              'structures(times tr(adad)/(,)(X1,X2)/p^2)': {'(rho(ip)rho(b)u,u*)': str(s1), '(rho(b)rho(ip)u,u*)': str(s2),
                                                                           '(b.p)(u,u*)': str(s0)}}

# ------------------------------------------------------------------ hand contraction (phi-lines only)
def phi_internal(k, tag):
    """internal modes of the slot A, W^perp components, all colours, momentum k."""
    out = []
    for i in range(4, 10):
        for a in range(3):
            g = cc.new_gen(0, 'int')
            out.append((g, (i, a), cc.Mode(cc.gen_elem(g), 'A', [GQ(x) for x in vec10({i: 1})], T_[a], k, internal=True)))
    return out


def one_leg(F, exts, ints):
    """d/d(internal) of F evaluated on exts + one internal leg: dict label -> algebra element."""
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


def line(da, db, k, prop=None):
    """sum over labels of da[l] * P * db[l], P = 2/|k|^2 (A_Wperp (x) A_Wperp, orthonormal colour)."""
    pr = Fr(2)/ksq(k) if prop is None else prop
    out = {}
    for lab, x in da.items():
        if lab in db:
            out = cc.aadd(out, cc.ascale(cc.amul(x, db[lab], maxint=0), pr))
    return out


# check of the hand contraction against cubic_core.tree in the scalar sector
tot = ZERO
for a in range(3):
    gA, mA = mode('A', vec10({mu: b[mu] for mu in range(4)}), T_[a], (0, 0, 0, 0))
    gc, mc = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))
    g1, m1 = mode('A', uW, X1, p)
    g2, m2 = mode('As', uWs, X2, neg(p))
    ints_a = phi_internal(neg(p), 'a')       # leg at the A-stub vertex: momentum -p
    ints_g = phi_internal(p, 'g')
    da = one_leg(cc.I3, [mA, m1], ints_a)
    dg = one_leg(cc.I3, [mc, m2], ints_g)
    val = line(da, dg, p)
    d = cc.dleft(cc.dleft(val, gc), gA)
    tot = tot + coeff(d, [g1, g2])
res['check_hand_contraction_vs_cubic_core_tree'] = {'hand': str(tot), 'tree': str(YW_phi), 'agree': (tot - YW_phi).iszero()}

# ------------------------------------------------------------------ (5) scalar bilinear O = int omega (iota_w1 A, iota_w2 A)
w1 = vec10({4: Fr(1), 5: Fr(1, 2), 7: Fr(-1)})
w2 = vec10({5: Fr(2), 6: Fr(-1, 3), 9: Fr(1)})


def O_two(l1, l2):
    """second derivative of O along the internal basis modes l = (i, a): (w1_i w2_j + w1_j w2_i)(T_a, T_b) (hat omega = 1)."""
    (i, a), (j, c) = l1, l2
    return (w1[i]*w2[j] + w1[j]*w2[i])*(1 if a == c else 0)


def O_one(u, X, l):
    (i, a) = l
    return (sum(w1[k]*u[k] for k in range(4, 10))*w2[i] + w1[i]*sum(w2[k]*u[k] for k in range(4, 10)))*X[a]


def YO(phi_mode_data, phis_mode_data, Q_on=None):
    """sum_a d_A d_c O_tree in the sector (phi(p1), phi^*(q1)), omega carrying -(p1+q1).
    phi_mode_data = (u, X, k), phis_mode_data = (u*, X*, k*).  Returns (alpha, beta) GQ values of the two trees."""
    (u, X, k1), (us, Xs, k2) = phi_mode_data, phis_mode_data
    alpha_t, beta_t = ZERO, ZERO
    for a in range(3):
        gA, mA = mode('A', vec10({mu: b[mu] for mu in range(4)}), T_[a], (0, 0, 0, 0))
        gc, mc = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))
        g1, m1 = mode('A', u, X, k1)
        g2, m2 = mode('As', us, Xs, k2)
        # (alpha): V(A-stub, phi_ext, line1) - O(line1, line2) - V(c-stub, line2, phi*_ext)
        ia = phi_internal(neg(k1), 'a'); ig = phi_internal(neg(k2), 'g')
        da = one_leg(cc.I3, [mA, m1], ia); dg = one_leg(cc.I3, [mc, m2], ig)
        val = {}
        for la, xa in da.items():
            for lg, xg in dg.items():
                o = O_two(la, lg)
                if o == 0:
                    continue
                val = cc.aadd(val, cc.ascale(cc.amul(xa, xg, maxint=0), o*Fr(2)/ksq(k1)*Fr(2)/ksq(k2)))
        dd = cc.dleft(cc.dleft(val, gc), gA)
        alpha_t = alpha_t + coeff(dd, [g1, g2])
        # (beta): O(phi_ext, line1) - V(line1, A-stub, line2) - V(line2, c-stub, phi*_ext); line momenta: k2-side
        kl = neg(k2)                                   # momentum of the phi line at the c-stub vertex
        ig2 = phi_internal(kl, 'g2')                    # at V_g
        im2 = phi_internal(neg(kl), 'm2')               # at V_a, towards V_g
        im1 = phi_internal(neg(k2), "m1")               # at V_a, towards O: 0 + k2 + (-k2) = 0
        dg2 = one_leg(cc.I3, [mc, m2], ig2)
        dm = two_legs(cc.I3, [mA], im1, im2)
        val = {}
        for (l1, l2), xm in dm.items():
            if l2 not in dg2:
                continue
            o = O_one(u, X, l1)
            if o == 0:
                continue
            t = cc.amul(xm, dg2[l2], maxint=0)
            val = cc.aadd(val, cc.ascale(t, o*Fr(2)/ksq(k2)*Fr(2)/ksq(kl)))
        # O carries no generator of its own here: multiply by the phi_ext generator to keep the extraction uniform
        val = cc.amul(cc.gen_elem(g1), val, maxint=0)
        dd = cc.dleft(cc.dleft(val, gc), gA)
        beta_t = beta_t + coeff(dd, [g1, g2])
    return alpha_t, beta_t


al, be = YO((uW, X1, p), (uWs, X2, q))
# hand prediction (k_transfer.tex): alpha = -64 i C (b.p) <O-structure> /(p^2 q^2) * ..., beta = +64 i C (b.q)/q^4 * ...
Ostruct_pu = sum(w1[k]*uW[k] for k in range(4, 10))*sum(w2[k]*uWs[k] for k in range(4, 10)) + \
    sum(w2[k]*uW[k] for k in range(4, 10))*sum(w1[k]*uWs[k] for k in range(4, 10))
pa = GQ(0, -4)*CAS*bdot(p)/(ksq(p)*ksq(q))*Ostruct_pu*cc.cp(X1, X2)
pb = GQ(0, 4)*CAS*bdot(q)/(ksq(q)*ksq(q))*Ostruct_pu*cc.cp(X1, X2)
res['YO_scalar_bilinear'] = {'alpha_tree(V-O-V)': str(al), 'beta_tree(O-V-V)': str(be), 'total': str(al + be),
                             'hand_alpha=-4 i C (b.p)/(p^2 q^2) S': str(pa), 'hand_beta=4 i C (b.q)/q^4 S': str(pb),
                             'agree_alpha': (al - pa).iszero(), 'agree_beta': (be - pb).iszero(),
                             'total_nonzero': not (al + be).iszero()}

# ------------------------------------------------------------------ (3) Q Y_O = 0 and Q Y_W = 0
# Q on functionals: (Q F)(e) = d/dr F(e + r Q e); on an A|_Wperp mode, Q e = -1/2 |k|^2 e in the slot A^*|_Wperp.
uW2 = vec10({4: Fr(1, 3), 5: Fr(1, 3), 7: Fr(2), 9: Fr(-1, 2)})   # not orthogonal to uW
X3 = [Fr(1), Fr(1, 5), Fr(-1, 2)]
Qimg = lambda u, k: [Fr(-1, 2)*ksq(k)*x for x in u]
a1, b1 = YO((uW, X1, p), (Qimg(uW2, q), X3, q))
a2, b2 = YO((uW2, X3, q), (Qimg(uW, p), X1, p))
QYO = a1 + b1 + a2 + b2
res['Q_YO'] = {'value': str(QYO), 'zero': QYO.iszero(),
               'alpha_part': str(a1 + a2), 'beta_part': str(b1 + b2)}
# Q Y_W in (A_Wperp, A_Wperp): Y_W(phi1; Q phi2) + Y_W(phi2; Q phi1)
def YW_phi_val(u, X, k, us, Xs):
    tot = ZERO
    for a in range(3):
        gA, mA = mode('A', vec10({mu: b[mu] for mu in range(4)}), T_[a], (0, 0, 0, 0))
        gc, mc = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))
        g1, m1 = mode('A', u, X, k)
        g2, m2 = mode('As', us, Xs, neg(k))
        da = one_leg(cc.I3, [mA, m1], phi_internal(neg(k), 'a'))
        dg = one_leg(cc.I3, [mc, m2], phi_internal(k, 'g'))
        d = cc.dleft(cc.dleft(line(da, dg, k), gc), gA)
        tot = tot + coeff(d, [g1, g2])
    return tot


QYW = YW_phi_val(uW, X1, p, Qimg(uW2, neg(p)), X3) + YW_phi_val(uW2, X3, neg(p), Qimg(uW, p), X1)
res['Q_YW_(A_Wperp,A_Wperp)'] = {'value': str(QYW), 'zero': QYW.iszero(),
                                 'single_terms_nonzero': [str(YW_phi_val(uW, X1, p, Qimg(uW2, neg(p)), X3)), str(YW_phi_val(uW2, X3, neg(p), Qimg(uW, p), X1))]}

# ------------------------------------------------------------------ (4) arity two: Z = sum_a d_{c=X^a} d_{c=(b,x) X_a} W_tree
# In Z the two c-stubs sit on two vertices ([A,c],A^*) joined by a phi line; the linear stub (b,x) X_a at the vertex
# with phi*_2 acts as -i (b . d_k) on the plane wave e^{i(k,x)} at k = 0; the only k-dependence is the line
# propagator 2/|m(k)|^2 with m(k) = -(k + k_2) at that vertex, so -i b.d_k (2/|m|^2) = -i * 4 (b.(k+k_2))/|k+k_2|^4 * (-1)...
# computed exactly: d/dk_mu |k+k2|^-2 = -2 (k+k2)_mu |k+k2|^-4.
def Z_val(us1, Xs1, k1, us2, Xs2, k2):
    """coefficient of kappa_1 kappa_2 in Z(kappa_1 phi*_1 + kappa_2 phi*_2), with phi*_i = (us_i, Xs_i, k_i), k1 + k2 = 0."""
    tot = ZERO
    for a in range(3):
        for lin_on in (1, 2):
            gc0, mc0 = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))          # constant stub X^a
            gcl, mcl = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))          # linear stub (b,x) X_a, momentum derivative below
            g1, m1 = mode('As', us1, Xs1, k1)
            g2, m2 = mode('As', us2, Xs2, k2)
            if lin_on == 2:
                ext0, extl, kl = [mc0, m1], [mcl, m2], k2
            else:
                ext0, extl, kl = [mc0, m2], [mcl, m1], k1
            # line momentum at the linear-stub vertex: -kl (stub momentum 0); derivative of the propagator:
            # -i (b . d_k) [2/|k + kl|^2] at k = 0 = -i * 2 * (-2) (b . kl)/|kl|^4 = 4 i (b.kl)/|kl|^4
            prop = GQ(0, 4)*bdot(kl)/(ksq(kl)**2)
            d0 = one_leg(cc.I3, ext0, phi_internal(kl, 'z0'))
            dl = one_leg(cc.I3, extl, phi_internal(neg(kl), 'zl'))
            val = {}
            for lab, x in dl.items():
                if lab in d0:
                    val = cc.aadd(val, cc.amul(x, d0[lab], maxint=0))
            val = {m: c*prop for m, c in val.items()}
            d = cc.dleft(cc.dleft(val, gcl), gc0)
            tot = tot + coeff(d, [g1, g2])
    return tot


# Q Z in the sector (phi, phi^*): Q acts on one phi^* argument: QZ(phi; phi^*) = d/dr Z(r Q phi + phi^*)
# with Z bilinear in two odd arguments: coefficient of kappa_phi kappa_s in d/dr Z(r kappa_phi Qphi + kappa_s phi^*)
# = Zc(Qphi, phi^*) [order (r kappa_phi) then kappa_s] -> d/dr removes r from the front: + Zc.
# The momentum derivative at fixed external momenta also differentiates the overall delta-function; summing the two
# assignments of the odd stubs, the delta'-pieces add up to minus one copy of the translation-invariant kernel
# b.(y_lin - y_const) P(y_1 - y_2), and each naive assignment equals one copy (k_transfer.tex, proof of
# Prop prop:ktr-consistency (iii)).  Hence Z = (1/2) * (naive sum computed by Z_val).
QZ = Z_val(Qimg(uW, p), X1, p, uWs, X2, neg(p))*Fr(1, 2)
YW_same = YW_phi_val(uW, X1, p, uWs, X2)
res['arity_two_Z_and_QZ'] = {'QZ_(phi(p),phi*(-p))': str(QZ), 'YW_(phi(p),phi*(-p))': str(YW_same),
                             'ratio_YW/QZ': str(YW_same/QZ) if not QZ.iszero() else 'QZ=0'}

# ------------------------------------------------------------------ (6) defect counterpart in the sector (A_Wperp, A_Wperp)
# special edge between A^*|_Wperp (at the c-stub vertex) and A|_Wperp (at the A-stub vertex): K_eps - K_0 in place of the
# A (x) A line; its scalar factor is (e^{-eps|p|^2} - 1) times the component of K_0 (Prop c4-defect1).  The tree value is
# (e^{-eps|p|^2}-1) * Ytilde(p) with Ytilde the contraction of the two vertices through K_0; Ytilde(p) + Ytilde(-p) = 0.
def Ytilde(u1, Xa_, k, u2, Xb_):
    tot = ZERO
    for a in range(3):
        gA, mA = mode('A', vec10({mu: b[mu] for mu in range(4)}), T_[a], (0, 0, 0, 0))
        gc, mc = mode('c', [Fr(1)], T_[a], (0, 0, 0, 0))
        g1, m1 = mode('A', u1, Xa_, k)
        g2, m2 = mode('A', u2, Xb_, neg(k))
        # at the c-stub vertex the special-edge end is the antifield slot A^*|_Wperp at momentum k
        ints = []
        for i in range(4, 10):
            for c_ in range(3):
                g = cc.new_gen(1, 'int')
                ints.append((g, (i, c_), cc.Mode(cc.gen_elem(g), 'As', [GQ(x) for x in vec10({i: 1})], T_[c_], k, internal=True)))
        dg = one_leg(cc.I3, [mc, m2], ints)
        da = one_leg(cc.I3, [mA, m1], phi_internal(neg(k), 'a'))
        val = {}
        for lab, x in da.items():
            if lab in dg:
                val = cc.aadd(val, cc.amul(x, dg[lab], maxint=0))      # K_0 component: pairing 1 (up to the global sign of K_0)
        d = cc.dleft(cc.dleft(val, gc), gA)
        tot = tot + coeff(d, [g1, g2])
    return tot


yt1 = Ytilde(uW, X1, p, uW2, X3)
yt2 = Ytilde(uW2, X3, neg(p), uW, X1)
res['defect_counterpart_(A_Wperp,A_Wperp)'] = {'Ytilde(phi1(p),phi2(-p))': str(yt1), 'Ytilde(phi2(-p),phi1(p))': str(yt2),
                                               'sum_over_legs': str(yt1 + yt2), 'sum_zero': (yt1 + yt2).iszero(),
                                               'full_value': '(e^{-eps|p|^2}-1) * sum_over_legs = 0'}
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
