# j_zeromode_trees_engineB.py -- judge's independent evaluation of the zero-mode trees of W along K_b,
#   Y_W = sum_a d_{A|_V = b (x) X_a} d_{c = X^a} W(P(0,infinity), I)_tree      (left derivatives),
# in the sectors (A|_Wperp, A^*|_Wperp), (A|_V, A^*|_V), (c, c^*), with the Step-1 vertex functions of ENGINE B
# (engine_b/vertices.py: independent realisation, su(2) in a deliberately non-orthonormal basis with a rescaled pairing).
# Neither the theory agent's cubic_core.py nor the numerics agent's kn engine (engine A core) is used.
#
# Conventions (ch03 Def c3-P / Prop c3-P): the value of d_P on the product of the coordinate functionals of a slot at x
# and a slot at y, in this order, in momentum space with e^{ikx} at the first slot:
#     A (x) A        : 2 (-,-)^{-1} delta_MN / k^2
#     c (x) A^*|_V   : -(-,-)^{-1} i k_m / k^2
#     A^*|_V (x) c   : -(-,-)^{-1} i k_m / k^2
# W = hbar log(e^{hbar d_P} e^{I/hbar}): every two-vertex tree with distinct legs has coefficient 1 (vertex functions are
# the full multilinear coefficients, engine B: coefficient of the odd parameters in leg order).
# Grassmann bookkeeping: explicit (odd parameters ordered; internal pair moved to adjacency in the order (X-leg, Y-leg)
# and replaced by the d_P value; the left derivative d_{c-mode} moves the c-mode parameter to the front).
# Prediction (theory, k_transfer.tex Prop prop:ktr-values):
#     (A_W,A*_W): -2 i C (b.p)/|p|^2 (u,u*)(X1,X2);   (A_V,A*_V): i C/|p|^2 [-2(b.p)(u.u*) + (b.u)(p.u*)](X1,X2);
#     (c,c*):     i C (b.p)/|p|^2 (X1,X2);     C = tr_g(ad_X ad_Y)/(X,Y) = c_adj of engine B.
# Then the coefficient of log T of the 1PI flow is (3/(8 pi^2)) hbar Y_W (judge: j_composite_kernel.py), i.e. in units
# hbar C/(16 pi^2): 6 Y_W / C, evaluated at the kinematics of the numerics agent (p1=(1,2,0,1), b=(1,0,-1,2), its
# polarisations) and compared with its stored slopes c_f (its X0 = minus the coefficient of delta).
import sys
import json
import itertools
import sympy as sp

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_b')
import vertices as VB          # noqa: E402
from realization import su2_data   # noqa: E402

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_judgement/indep/j_zeromode_trees_engineB.json'
L = su2_data()
n = L['n']; ft = L['ft']; kinv = L['kinv']; kap = L['kap']; C = L['cadj']
PAR = VB.PAR
I = sp.I


def unitc(a):
    return [sp.Integer(1) if i == a else sp.Integer(0) for i in range(n)]


def dualc(a):
    return [kinv[a, e] for e in range(n)]


def ftil(c0, c1, c2):
    return sum(c0[i]*c1[j]*c2[k]*ft[i, j, k] for i in range(n) for j in range(n) for k in range(n) if ft[i, j, k] != 0)


def pairc(c0, c1):
    return sum(c0[i]*kap[i, j]*c1[j] for i in range(n) for j in range(n))


def vertex_value(legs):
    """legs: list of dict(sp, u, p, col) -> full vertex value (colour contracted), cubic terms only"""
    tot = 0
    for t in VB.CUBIC:
        for key, val in VB.vertex(t, [dict(sp=l['sp'], u=l['u'], p=l['p']) for l in legs]).items():
            assert key == 'f'
            tot += val*ftil(legs[0]['col'], legs[1]['col'], legs[2]['col'])
    return sp.expand(tot)


def prop_value(spX, spY, k, mX, mY):
    """d_P value on (coordinate functional of slot spX, component mX, at x with e^{ikx}) (x) (spY, mY at y)"""
    k2 = sum(x**2 for x in k)
    if spX == 'A' and spY == 'A':
        return 2*(1 if mX == mY else 0)/k2
    if spX == 'c' and spY == 'As':
        return -I*k[mY]/k2 if mY < 4 else 0
    if spX == 'As' and spY == 'c':
        return -I*k[mX]/k2 if mX < 4 else 0
    return 0


def comps(spc):
    if spc == 'A':
        return list(range(10))
    if spc == 'As':
        return list(range(4))          # only A^*|_V is internal (c (x) A^*|_V)
    if spc == 'c':
        return [0]
    raise ValueError(spc)


def basis_u(spc, m):
    d = VB.DIM[spc]
    return [sp.Integer(1) if i == m else sp.Integer(0) for i in range(d)]


def sign_move_to_front(order, x):
    i = order.index(x)
    s = (-1)**i
    return s, [x] + order[:i] + order[i + 1:]


def tree_sector(ext1, ext2, bvec, p):
    """ext1 = (sp, u, col) at momentum p, parameter r1; ext2 = (sp, u, col) at -p, parameter r2.
    Returns Y_W = coefficient of r1 r2 (in this order) of sum_a d_{kappa_A} d_{kappa_c} W_tree."""
    zero = [0, 0, 0, 0]
    mp_ = [-x for x in p]
    total = 0
    detail = {}
    for a in range(n):
        modeA = dict(sp='A', u=[bvec[i] if i < 4 else 0 for i in range(10)], p=zero, col=unitc(a), name='kA')
        modec = dict(sp='c', u=[1], p=zero, col=dualc(a), name='kc')
        e1 = dict(sp=ext1[0], u=ext1[1], p=list(p), col=ext1[2], name='r1')
        e2 = dict(sp=ext2[0], u=ext2[1], p=mp_, col=ext2[2], name='r2')
        # partitions: (ext1 + one mode | ext2 + other mode)
        for (m1, m2) in ((modeA, modec), (modec, modeA)):
            kX = mp_                              # incoming momentum of the internal leg at X  (p + 0 + kX = 0)
            kY = list(p)
            for spX, spY in (('A', 'A'), ('As', 'c'), ('c', 'As')):
                for mX in comps(spX):
                    for mY in comps(spY):
                        pv = prop_value(spX, spY, kX, mX, mY)
                        if pv == 0:
                            continue
                        for cX in range(n):
                            for cY in range(n):
                                if kinv[cX, cY] == 0:
                                    continue
                                iX = dict(sp=spX, u=basis_u(spX, mX), p=kX, col=unitc(cX), name='iX')
                                iY = dict(sp=spY, u=basis_u(spY, mY), p=kY, col=unitc(cY), name='iY')
                                LX = [e1, m1, iX]
                                LY = [iY, m2, e2]
                                vX = vertex_value(LX)
                                if vX == 0:
                                    continue
                                vY = vertex_value(LY)
                                if vY == 0:
                                    continue
                                order = [l['name'] for l in LX + LY if PAR[l['sp']]]
                                sgn = 1
                                if PAR[spX]:
                                    s1, order = sign_move_to_front(order, 'iX')
                                    # bring iY right after iX
                                    rest = order[1:]
                                    s2, rest = sign_move_to_front(rest, 'iY')
                                    order = ['iX'] + rest
                                    sgn *= s1*s2
                                    order = order[2:]
                                # left derivative along the odd c-mode
                                s3, order = sign_move_to_front(order, 'kc')
                                order = order[1:]
                                sgn *= s3
                                # remaining odd parameters must be the odd ones of (r1, r2) in this order
                                target = [x for x, l in (('r1', e1), ('r2', e2)) if PAR[l['sp']]]
                                assert sorted(order) == sorted(target)
                                if order != target:
                                    sgn *= -1
                                val = sgn*vX*pv*kinv[cX, cY]*vY
                                total += val
                                key = '%s|%s%s' % ('A-at-X' if m1 is modeA else 'c-at-X', spX, spY)
                                detail[key] = detail.get(key, 0) + val
    return sp.simplify(total), {k: str(sp.simplify(v)) for k, v in detail.items()}


res = {'lie_algebra': 'su(2), engine B non-orthonormal basis, pairing scale 3/2', 'C = c_adj': str(C)}
dot = lambda x, y: sum(a*c for a, c in zip(x, y))
# generic rational kinematics
p = [sp.Rational(1, 2), sp.Rational(-1, 3), sp.Rational(2, 5), sp.Rational(1, 7)]
b = [sp.Rational(2, 3), sp.Rational(-1, 5), sp.Rational(1, 4), sp.Rational(3, 7)]
X1, X2 = unitc(0), unitc(2)
p2 = dot(p, p)
# (A_W, A*_W)
uW = [0, 0, 0, 0, sp.Integer(1), 0, sp.Rational(-2, 3), 0, 0, sp.Rational(1, 2)]
uWs = [0, 0, 0, 0, sp.Rational(3, 4), 0, sp.Rational(1, 5), 0, -1, 0]
Y, det = tree_sector(('A', uW, X1), ('As', uWs, X2), b, p)
pred = -2*I*C*dot(b, p)/p2*dot(uW, uWs)*pairc(X1, X2)
res['(A_W,A*_W)'] = {'Y_W': str(Y), 'prediction': str(sp.nsimplify(pred)), 'ratio': str(sp.nsimplify(sp.simplify(Y/pred))), 'trees': det}
# (A_V, A*_V)
uV = [sp.Rational(1, 3), -1, sp.Rational(1, 2), 2, 0, 0, 0, 0, 0, 0]
uVs = [sp.Rational(-1, 2), sp.Rational(2, 3), 1, sp.Rational(-1, 4), 0, 0, 0, 0, 0, 0]
Y, det = tree_sector(('A', uV, X1), ('As', uVs, X2), b, p)
pred = I*C/p2*(-2*dot(b, p)*dot(uV, uVs) + dot(b, uV)*dot(p, uVs))*pairc(X1, X2)
res['(A_V,A*_V)'] = {'Y_W': str(Y), 'prediction': str(sp.nsimplify(pred)), 'ratio': str(sp.nsimplify(sp.simplify(Y/pred))), 'trees': det}
# (c, c*)
Y, det = tree_sector(('c', [1], X1), ('cs', [1], X2), b, p)
pred = I*C*dot(b, p)/p2*pairc(X1, X2)
res['(c,c*)'] = {'Y_W': str(Y), 'prediction': str(sp.nsimplify(pred)), 'ratio': str(sp.nsimplify(sp.simplify(Y/pred))), 'trees': det}

# ------------------------------------------------------------------ comparison with the numerics agent's slopes c_f
import numpy as np
P1 = [1, 2, 0, 1]; BV = [1, 0, -1, 2]
rng3 = np.random.default_rng(3).integers(-2, 3, size=4); rng4 = np.random.default_rng(4).integers(-2, 3, size=4)
cmp = {}
pp = sum(x*x for x in P1); bp = sum(x*y for x, y in zip(BV, P1))
# log T coefficient of the graph (coefficient of W) in units hbar C/(16 pi^2):  6 * (Y_W / C)
cmp['(A_W,A*_W)'] = {'judge_6Y_W/C': str(6*(-2*I*sp.Rational(bp, pp))), 'numerics_graph=-X0_slope': str(-6*I),
                     'numerics_c_f(X0)': '6i'}
uu = [int(x) for x in rng3]; us = [int(x) for x in rng4]
vv = 6*I*sp.Rational(-2*bp*sum(x*y for x, y in zip(uu, us)) + sum(x*y for x, y in zip(BV, uu))*sum(x*y for x, y in zip(P1, us)), pp)
cmp['(A_V,A*_V)'] = {'u': uu, 'u*': us, 'judge_6Y_W/C': str(vv), 'numerics_graph=-X0_slope': str(46*I), 'numerics_c_f(X0)': '-46i'}
cmp['(c,c*)'] = {'judge_6Y_W/C': str(6*I*sp.Rational(bp, pp)), 'numerics_graph=-X0_slope': str(3*I), 'numerics_c_f(X0)': '-3i'}
res['comparison_with_numerics_agent_flow_slopes_(units)'] = cmp
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
