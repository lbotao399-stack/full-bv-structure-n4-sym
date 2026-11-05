"""
The non-minimal sector alone (frame C): the trivial pair (cbar, b) with antifields, and the reduced pair
(cbar, cbar^*), their equivariant extensions under the even part of psl(4|4), the arity-two term, and the
elimination of b along the gauge fermion kappa int (cbar, b).
cbar, b have Weyl weight 2; v.f = v(f) + (1/2) div v f on both.
"""
import sympy as sp
import bvjet
from test_equiv import gens_even

X = bvjet.X


def model(with_b=True):
    M = bvjet.Model()
    M.add('cb', 1)
    M.add('cbs', 0)
    M.pair('cb', 'cbs', -1)
    if with_b:
        M.add('b', 0)
        M.add('bs', 1)
        M.pair('b', 'bs', -1)
    return M


def vdot(M, v, P, w=2):
    out = {}
    for m in range(4):
        if v[m] != 0:
            out = M.plus(out, M.mul(M.const(v[m]), M.Dtot(P, m)))
    dv = sp.expand(sum(sp.diff(v[m], X[m]) for m in range(4)))
    if dv != 0:
        out = M.plus(out, M.mul(M.const(sp.Rational(w, 4) * dv), P))
    return out


def ckv_bracket(v, w):
    return [sp.expand(sum(v[n] * sp.diff(w[m], X[n]) - w[n] * sp.diff(v[m], X[n]) for n in range(4)))
            for m in range(4)]


def run():
    res = {}
    vs = gens_even
    # ---------------- reduced pair
    M = model(False)
    cb, cbs = M.jet('cb'), M.jet('cbs')
    Snm = M.scale(M.mul(cbs, cbs), -sp.Rational(1, 4))
    S1 = lambda v: M.scale(M.mul(vdot(M, v, cb), cbs), -1)
    S2 = lambda u, g=-1: M.scale(M.mul(cb, vdot(M, u, cb, 0)), g)
    res['red_ar1_{Snm,S1(v)}=0'] = all(M.is_zero_functional(M.bracket(Snm, S1(vs[k]))) for k in vs)
    res['red_ar2_even-even'] = all(M.is_zero_functional(M.plus(M.bracket(S1(vs[a]), S1(vs[b])), S1(ckv_bracket(vs[a], vs[b]))))
                                   for a in vs for b in vs)
    res['red_ar2_odd-odd_{Snm,S2(u)}+S1(u)=0_gamma2=-1'] = all(M.is_zero_functional(M.plus(M.bracket(Snm, S2(vs[k])), S1(vs[k]))) for k in vs)
    res['NEG_red_gamma2=+1_fails'] = not M.is_zero_functional(M.plus(M.bracket(Snm, S2(vs['K'], 1)), S1(vs['K'])))
    # arity three with one even and two odd arguments, reduced to: {S1(v), S2(u)} + S2([v,u]) = 0
    res['red_ar3_{S1(v),S2(u)}+S2([v,u])=0'] = all(M.is_zero_functional(M.plus(M.bracket(S1(vs[a]), S2(vs[b])), S2(ckv_bracket(vs[a], vs[b]))))
                                                 for a in vs for b in vs)
    res['red_ar4_{S2,S2}=0'] = M.is_zero_functional(M.bracket(S2(vs['K']), S2(vs['P'])))
    # ---------------- trivial pair (cbar, b), frame C:  S_tp = int (b, cbar^*)
    M = model(True)
    cb, cbs, bb, bs = M.jet('cb'), M.jet('cbs'), M.jet('b'), M.jet('bs')
    Stp = M.mul(bb, cbs)
    T1 = lambda v: M.plus(M.scale(M.mul(vdot(M, v, cb), cbs), -1), M.mul(vdot(M, v, bb), bs))
    gb = sp.Symbol('gamma_b')
    T2 = lambda u, g=gb: M.scale(M.mul(vdot(M, u, cb), bs), g)
    res['tp_ar1'] = all(M.is_zero_functional(M.bracket(Stp, T1(vs[k]))) for k in vs)
    res['tp_ar2_even-even'] = all(M.is_zero_functional(M.plus(M.bracket(T1(vs[a]), T1(vs[b])), T1(ckv_bracket(vs[a], vs[b]))))
                                  for a in vs for b in vs)
    lhs = M.plus(M.bracket(Stp, T2(vs['K'])), T1(vs['K']))
    eqs = set()
    for a in range(len(M.fields)):
        for k, c in M.simplify(M.EL(lhs, a)).items():
            eqs |= set(sp.Poly(sp.expand(c), *X).coeffs())
    sol = sp.solve(list(eqs), gb, dict=True)
    res['tp_gamma_b'] = str(sol)
    g = sol[0][gb]
    res['tp_ar2_odd-odd_all_u'] = all(M.is_zero_functional(M.plus(M.bracket(Stp, T2(vs[k], g)), T1(vs[k]))) for k in vs)
    res['tp_ar3'] = all(M.is_zero_functional(M.plus(M.bracket(T1(vs[a]), T2(vs[b], g)), T2(ckv_bracket(vs[a], vs[b]), g)))
                        for a in vs for b in vs)
    # elimination of b along Psi_1 = kappa (cbar, b): b^* -> b^* + kappa cbar, cbar^* -> cbar^* - kappa b, then b^* = 0
    kap = -1
    sub_cbs = M.plus(cbs, M.scale(bb, -kap))
    sub_bs = M.scale(cb, kap)                     # b^* + kappa cbar at b^* = 0
    Sk = M.mul(bb, sub_cbs)
    # E_b(Sk) = 0 :  cbar^* - 2 kappa b = 0  ->  b = cbar^*/(2 kappa)
    Eb = M.simplify(M.EL(Sk, M.index['b']))
    res['tp_Eb_of_S_tp_Psi1'] = str(Eb)
    bsol = M.scale(cbs, sp.Rational(1, 2 * kap))
    Sred = M.plus(M.mul(bsol, cbs), M.scale(M.mul(bsol, bsol), -kap))
    Mr = model(False)
    res['tp_reduced_action_equals_-1/4(cbar*,cbar*)'] = M.is_zero_functional(M.plus(Sred, M.scale(M.mul(cbs, cbs), sp.Rational(1, 4))))
    ok1 = True
    for k in vs:
        v = vs[k]
        T1k = M.plus(M.scale(M.mul(vdot(M, v, cb), sub_cbs), -1), M.mul(vdot(M, v, bb), sub_bs))
        # b-dependence is a total derivative: the remainder is -(v.cbar, cbar^*)
        if not M.is_zero_functional(M.plus(T1k, M.scale(M.mul(vdot(M, v, cb), cbs), 1))):
            ok1 = False
    res['tp_S1_after_Psi1_at_bs=0_equals_-(v.cbar,cbar*)_mod_total_derivatives'] = ok1
    T2k = M.scale(M.mul(vdot(M, vs['K'], cb), sub_bs), g)
    red2 = M.scale(M.mul(cb, vdot(M, vs['K'], cb, 0)), -1)
    res['tp_S2_after_Psi1_equals_reduced_cbar_cbar_vertex'] = M.is_zero_functional(M.plus(T2k, M.scale(red2, -1)))
    return res


if __name__ == '__main__':
    for k, v in run().items():
        print(k, ':', v)
