"""Audit 2b: independent recomputation of Q Pi_Psi in the (c, A*_V) and (A_V, A_V) output sectors (Part D, thm:d-D2 (ii)).

Q(class (i)) = (hbar/2) d_Psi(Q W_0) = -hbar x [special-edge composite graphs]: the composite line (S^(1)(Delta) with
edges P_AA on its A leg and P_{A* c} on its A*|_V leg) closes a loop on the tree I_3 -(K_eps - K_0)- I_3.
Complete list of such graphs with the given external legs (one external leg on each cubic vertex):
  (c, A*_V) output:  (a) V1 = ([A_V(Psi), c_ext], A*_X), V2 = ([A_X, c(Psi)], A*_ext), K-hat (A* -- A)
                     (b) V1 = ([A_V(Psi), c_X], A*_ext), V2 = (1/2)([c(Psi), c_ext], c*_X), K-hat (c -- c*)
  (A_V, A_V) output: V1 = AAA(A_ext, A_X, A_V(Psi)), V2 = ([A_ext, c(Psi)], A*_X), K-hat (A -- A*), both assignments.
Predictions (engine B Ob - 4 C_log, and Q of engine A's class (i), see q_of_engineA_classi.json):
  Q Pi_Psi|(c,A*) = -(2/3) int(dc, A*),  Q Pi_Psi|(A_V,A_V) = +(1/3) int(dA_V, dA_V)   (units).
Hence Ob(Delta)|(c,A*) = -4 - 2/3 = -14/3 and Ob(Delta)|(A_V,A_V) = -4 + 1/3 = -11/3 (coefficient of int(dA_V,dA_V)).
Writes vc_sectors_from_scratch.json.  Single process, sympy.
"""
import itertools, json, os, sys
import sympy as sp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from perp_sector_from_scratch import (eps, u, t, sk, m, kap, PV, PSQ, NV, DIM, COL, trunc, exp_shift, inv_shift, levi,
                                      unitvec, V_AAA, V_AcAs, koszul, coeff_t_sk, to_units_triangle, TMAX)

HERE = os.path.dirname(os.path.abspath(__file__))
Khat = lambda q: 1 - sp.exp(-eps * u) * exp_shift(q)          # -(e^{-eps|m+q|^2} - 1), argument m + q
PAA = 2 * sp.exp(-eps * u) / u                                # P_AA at momentum -m

def PAsc_m(kv):
    """P_{A* c} at the A* momentum k6 = -m - kv: -i k6_nu e^{-eps k6^2}/k6^2 = i (m+kv)_nu e^{..}/|m+kv|^2"""
    E = exp_shift(kv); Iv = inv_shift(kv)
    return [trunc(sp.I * (m[nu] + kv[nu]) * E * Iv) * sp.exp(-eps * u) / u for nu in range(4)]

def V_ccc(c1, c2, cs):
    """(1/2)([c,c],c*) vertex function, legs (c, c, c*), coefficient of r1 r2: ([T1,T2],T3)"""
    return levi(c1, c2, cs)

def finish(G_w, G_nu, order):
    """alpha = G_w + sum_nu (-i d/dkappa_nu G_nu)|0, beta^nu = -i G_nu|0, at order p^order, in units, with the -hbar"""
    aw = to_units_triangle(coeff_t_sk(G_w, order, 0))
    aE = 0
    for nu in range(4):
        aE += -sp.I * to_units_triangle(sp.diff(coeff_t_sk(G_nu[nu], order, 1), kap[nu]))
    be = [-sp.I * to_units_triangle(coeff_t_sk(G_nu[nu], order, 0)) for nu in range(4)]
    return sp.simplify(aw + aE), [sp.simplify(b) for b in be]

def cAs_graphs(nu0):
    """external c (colour 0, momentum p, leg X1) and A*_V (component nu0, colour 0, momentum p2 = -p - kappa, leg X2)"""
    p = PV * t; kv = kap * sk; p2 = -p - kv
    w = unitvec(nu0)
    G_w = 0; G_nu = [0, 0, 0, 0]
    PAsc = PAsc_m(kv)
    # ---- graph (a): V1(A_Psi k2, c_ext p, A*_X k1), V2(A_X k3, c_Psi k4, A*_ext p2), S(A k5 = m, A* k6 = -m - kv)
    k2 = -m; k1 = -p - k2                                     # k1 = m - p
    par = [0, 1, 1, 0, 1, 1, 0, 1]
    names = ['A_Psi', 'c_ext', 'As_X', 'A_X', 'c_Psi', 'As_ext', 'A_S', 'As_S']
    order_ = [names.index(x) for x in ['As_X', 'A_X', 'A_Psi', 'A_S', 'As_S', 'c_Psi', 'c_ext', 'As_ext']]
    sgn = koszul(par, order_)
    KH = Khat(-p)                                             # K-hat at k1 = m - p (A* first slot; symmetric fibre)
    for M, cb, cd in itertools.product(range(NV), range(COL), range(COL)):
        eA = unitvec(M)
        v1 = V_AcAs(eA, cb, 0, eA, cd)                        # ([A_Psi(b), c_ext(0)], A*_X(d)), A*_X component M
        v2 = V_AcAs(eA, cd, cb, w, 0)                         # ([A_X(d), c_Psi(cb)], A*_ext(0)): needs M = nu0
        if v1 == 0 or v2 == 0:
            continue
        base = trunc(sgn * v1 * v2 * KH * PAA * PAsc[M])
        G_w += 2 * base
        for nu in range(4):
            G_nu[nu] += trunc(2 * sp.I * m[nu] * base)
    # ---- graph (b): V1(A_Psi k2, c_X k1, A*_ext p2), V2(c_Psi k4, c_ext p, c*_X k3), S(A m, A* -m - kv)
    k2 = -m; k1 = -k2 - p2                                    # k1 = m + p + kv
    par = [0, 1, 1, 1, 1, 0, 0, 1]
    names = ['A_Psi', 'c_X', 'As_ext', 'c_Psi', 'c_ext', 'cs_X', 'A_S', 'As_S']
    order_ = [names.index(x) for x in ['c_X', 'cs_X', 'A_Psi', 'A_S', 'As_S', 'c_Psi', 'c_ext', 'As_ext']]
    sgn = koszul(par, order_)
    KH = Khat(p + kv)                                         # K-hat (c x c*) at k1 = m + p + kv
    for M, cb, cd in itertools.product(range(NV), range(COL), range(COL)):
        eA = unitvec(M)
        v1 = V_AcAs(eA, cb, cd, w, 0)                         # ([A_Psi(b), c_X(d)], A*_ext(0)): M = nu0
        v2 = V_ccc(cb, 0, cd)                                 # ([c_Psi(cb), c_ext(0)], c*_X(d))
        if v1 == 0 or v2 == 0:
            continue
        base = trunc(sgn * v1 * v2 * KH * PAA * PAsc[M])
        G_w += 2 * base
        for nu in range(4):
            G_nu[nu] += trunc(2 * sp.I * m[nu] * base)
    out = {}
    for order in range(0, 3):
        al, be = finish(sp.expand(G_w), [sp.expand(g) for g in G_nu], order)
        out[order] = (al, be)
    return out

def VV_graphs(a, b):
    """external A_V legs: X1 (component a, momentum p), X2 (component b, momentum p2 = -p - kappa), colours 0;
    both assignments of (X1, X2) to (V1 = AAA, V2 = AcA*)."""
    p = PV * t; kv = kap * sk; p2 = -p - kv
    G_w = 0; G_nu = [0, 0, 0, 0]
    PAsc = PAsc_m(kv)
    for asg in (1, 2):
        (q1, e1), (q2, e2) = ((p, unitvec(a)), (p2, unitvec(b))) if asg == 1 else ((p2, unitvec(b)), (p, unitvec(a)))
        k2 = -m; k1 = -q1 - k2                                # V1: A_ext q1, A_X k1, A_Psi k2 ; K-hat at k1 = m - q1
        par = [0, 0, 0, 0, 1, 1, 0, 1]
        names = ['X_V1', 'A_X', 'A_Psi', 'X_V2', 'c_Psi', 'As_X', 'A_S', 'As_S']
        order_ = [names.index(x) for x in ['A_X', 'As_X', 'A_Psi', 'A_S', 'As_S', 'c_Psi', 'X_V1', 'X_V2']]
        sgn = koszul(par, order_)
        KH = Khat(-q1)
        for X, M, ca, cb in itertools.product(range(NV), range(NV), range(COL), range(COL)):
            eX = unitvec(X); eM = unitvec(M)
            v1 = V_AAA(e1, list(q1), 0, eX, list(k1), ca, eM, list(k2), cb)
            if v1 == 0:
                continue
            v2 = V_AcAs(e2, 0, cb, eX, ca)                    # ([A_ext(0), c_Psi(cb)], A*_X(ca)), A*_X component X
            if v2 == 0:
                continue
            base = trunc(sgn * v1 * v2 * KH * PAA * PAsc[M])
            G_w += 2 * base
            for nu in range(4):
                G_nu[nu] += trunc(2 * sp.I * m[nu] * base)
    out = {}
    for order in range(0, 4):
        al, be = finish(sp.expand(G_w), [sp.expand(g) for g in G_nu], order)
        out[order] = (al, be)
    return out

if __name__ == '__main__':
    res = {'(c,A*_V)': {}, '(A_V,A_V)': {}}
    pv = [1, 2, 0, -1]
    for nu0 in range(4):
        r = cAs_graphs(nu0)
        tgt = -sp.Rational(2, 3) * sp.I * pv[nu0]           # -(2/3) int(dc,A*) on (X_c at p, X_A* comp nu0): -(2/3) i p_nu0
        res['(c,A*_V)']['A* component %d' % nu0] = {
            'order p^%d' % k: {'alpha': str(r[k][0]), 'beta': [str(x) for x in r[k][1]]} for k in r}
        res['(c,A*_V)']['A* component %d' % nu0]['alpha(order p^1) / [-(2/3) i p_nu0] (expected 1, or 0 = 0)'] = \
            str(sp.simplify(r[1][0] / tgt)) if tgt != 0 else 'target 0, alpha = ' + str(r[1][0])
    for (a, b) in [(2, 2), (0, 1), (0, 0), (1, 3)]:
        r = VV_graphs(a, b)
        tgt = 2 * sp.Rational(1, 3) * (PSQ * (1 if a == b else 0) - pv[a] * pv[b])
        res['(A_V,A_V)']['(%d,%d)' % (a, b)] = {
            'order p^%d' % k: {'alpha': str(r[k][0]), 'beta': [str(x) for x in r[k][1]]} for k in r}
        res['(A_V,A_V)']['(%d,%d)' % (a, b)]['alpha(order p^2) - target (1/3) int(dA_V,dA_V)'] = str(sp.simplify(r[2][0] - tgt))
        print('done', a, b, flush=True)
    with open(os.path.join(HERE, 'vc_sectors_from_scratch.json'), 'w') as f:
        json.dump(res, f, indent=1)
    print(json.dumps(res, indent=1))
