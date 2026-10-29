# T1 audit: log eps coefficients (multiples of the tree vertex, unit hbar tr_g(ad ad)/((X,Y)16pi^2))
# of the (A|_V)^3 vertex and of the (lambda, lambda, A|_Wperp) vertex, plus the quadratic regressions
# (A|_V two-point and lambda two-point) that fix the normalisation.  Fresh script; imports only core.py.
import json
import sys
import numpy as np
from fractions import Fraction
import core as C

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
RANGE = int(sys.argv[2]) if len(sys.argv) > 2 else 3
np.seterr(all="ignore")


def rand_pol_V():
    a = np.zeros((C.ND, C.NC))
    a[:C.NV, :] = rng.integers(-RANGE, RANGE + 1, size=(C.NV, C.NC))
    return a.reshape(-1).astype(complex)


def rand_pol_W():
    a = np.zeros((C.ND, C.NC))
    a[C.NV:, :] = rng.integers(-RANGE, RANGE + 1, size=(C.ND - C.NV, C.NC))
    return a.reshape(-1).astype(complex)


def rand_mom():
    return np.array(rng.integers(-RANGE, RANGE + 1, size=4), dtype=float)


def unit(mu):
    v = np.zeros(4)
    v[mu] = 1
    return v


# ---------------------------------------------------------------- bosonic vertex as PM
_g3cache = {}


def G3(moms):
    key = tuple(tuple(m) for m in moms)
    if key not in _g3cache:
        _g3cache[key] = C.gamma3([np.asarray(m) for m in moms])
    return _g3cache[key]


def vert_boson_pm(a_bg, p_bg, P):
    """v(k + tau P) for background (a_bg, p_bg): Gamma3(p_bg; -(k+P+p_bg); k+P), split linear in k (tau^0)
    plus constant (tau^1)."""
    out = {}
    z = np.zeros(4)
    for mu in range(4):
        e = unit(mu)
        L = np.einsum('i,iwu->wu', a_bg, G3([z, -e, e]))
        if np.any(L != 0):
            out[(0, 0, C.unitvec(mu))] = L
    Cc = np.einsum('i,iwu->wu', a_bg, G3([p_bg, -(P + p_bg), P]))
    if np.any(Cc != 0):
        out[(1, 0, C.Z4)] = Cc
    return out


def quartic_pm(aj, ak):
    G4 = C.gamma4()
    M = np.einsum('i,j,ijwu->wu', aj, ak, G4)
    return {(0, 0, C.Z4): M.astype(complex)}


def chain(factors, taumax):
    P = factors[0]
    for f in factors[1:]:
        P = C.pm_mul(P, f, taumax)
    return P


def proj(indices, n):
    Pm = np.zeros((n, n), dtype=complex)
    for i in indices:
        Pm[i, i] = 1
    return Pm


VIDX = [C.bidx(M, a) for M in range(C.NV) for a in range(C.NC)]
WIDX = [C.bidx(M, a) for M in range(C.NV, C.ND) for a in range(C.NC)]
PV = proj(VIDX, C.NB)
PW = proj(WIDX, C.NB)


def triangle(kind, vs, ps, order, taumax):
    """Tr(G V_a G V_b G V_c) integrand PM for ordering (a,b,c):
    tr[G(q) v_a(q+p_b+p_c) G(q+p_b+p_c) v_b(q+p_c) G(q+p_c) v_c(q)]."""
    a, b, c = order
    pb, pc = ps[b], ps[c]
    if kind == 'B':
        fa = vert_boson_pm(vs[a], ps[a], pb + pc)
        fb = vert_boson_pm(vs[b], ps[b], pc)
        fc = vert_boson_pm(vs[c], ps[c], np.zeros(4))
    else:
        fa, fb, fc = (C.const_pm(vs[a]), C.const_pm(vs[b]), C.const_pm(vs[c]))
    return chain([C.prop_pm(kind, np.zeros(4), taumax), fa, C.prop_pm(kind, pb + pc, taumax), fb,
                  C.prop_pm(kind, pc, taumax), fc], taumax)


def bubble_q(vi, pi, P_jk, vjk, taumax):
    """Tr(G V_i G V_jk): tr[G(q) v_i(q + p_j + p_k) G(q + p_j + p_k) v_jk]."""
    fi = vert_boson_pm(vi, pi, P_jk)
    return chain([C.prop_pm('B', np.zeros(4), taumax), fi, C.prop_pm('B', P_jk, taumax), vjk], taumax)


def run_AVAVAV():
    p1, p2 = rand_mom(), rand_mom()
    p3 = -p1 - p2
    ps = [p1, p2, p3]
    a = [rand_pol_V(), rand_pol_V(), rand_pol_V()]
    tree = np.einsum('i,j,k,ijk->', a[0], a[1], a[2], G3(ps))
    tree_c = C.cfrac(tree)
    res = {}
    # bosonic triangles: +1/2 (T123 + T132)
    TB = [triangle('B', a, ps, o, 1) for o in [(0, 1, 2), (0, 2, 1)]]
    for name, Pr in [('triangle A|_V loop', PV), ('triangle A|_Wperp loop', PW)]:
        s = C.cadd(*[C.log_avg(T, 1, Pr) for T in TB])
        res[name] = C.cscale(s, Fraction(1, 2))
    # bosonic bubbles with quartic vertex: +1/2 sum_i Tr(G V_i G V_jk)
    bubs = []
    for i, (j, k) in [(0, (1, 2)), (1, (0, 2)), (2, (0, 1))]:
        vjk = quartic_pm(a[j], a[k])
        bubs.append(bubble_q(a[i], ps[i], ps[j] + ps[k], vjk, 1))
    for name, Pr in [('bubble A|_V loop (quartic)', PV), ('bubble A|_Wperp loop (quartic)', PW)]:
        s = C.cadd(*[C.log_avg(T, 1, Pr) for T in bubs])
        res[name] = C.cscale(s, Fraction(1, 2))
    # gaugino triangle: -1/2 (T123 + T132)
    Y = [C.Y_of(x) for x in a]
    TF = [triangle('F', Y, ps, o, 1) for o in [(0, 1, 2), (0, 2, 1)]]
    res['triangle lambda loop'] = C.cscale(C.cadd(*[C.log_avg(T, 1) for T in TF]), Fraction(-1, 2))
    # ghost triangle (c - A^*|_V lines): -1/2 (T123 + T132)
    Vg = [C.ghost_V(x) for x in a]
    TG = [triangle('G', Vg, ps, o, 1) for o in [(0, 1, 2), (0, 2, 1)]]
    res['triangle c-A^* loop'] = C.cscale(C.cadd(*[C.log_avg(T, 1) for T in TG]), Fraction(-1, 2))
    # number = <f>/(2 tree)   (C_log = -<f>/(16pi^2), unit = -2/(16pi^2))
    out = {}
    for k, v in res.items():
        out[k] = C.cdiv(v, C.cscale(tree_c, 2))
    tot = C.cadd(*res.values())
    out['TOTAL'] = C.cdiv(tot, C.cscale(tree_c, 2))
    return out, {'p1': p1.tolist(), 'p2': p2.tolist(), 'tree': C.cstr(tree_c)}


def run_AV2():
    """A|_V two-point log coefficient as multiple of the tree -1/4(dA,dA) (transverse pols)."""
    p = rand_mom()
    while np.dot(p, p) == 0:
        p = rand_mom()
    ps = [p, -p]
    # transverse polarisations
    def tpol():
        a = np.zeros((C.ND, C.NC))
        for col in range(C.NC):
            v = rng.integers(-RANGE, RANGE + 1, size=4).astype(float)
            v = v * np.dot(p, p) - p * np.dot(p, v)
            a[:C.NV, col] = v
        return a.reshape(-1).astype(complex)
    a = [tpol(), tpol()]
    tree = C.tree_quadratic_A(a[0], a[1], p)
    tree_c = C.cfrac(tree)
    # Tr(G V1 G V2) = tr[G(q) v1(q+p2) G(q+p2) v2(q)]
    def bub(kind, v1, v2, taumax=2):
        if kind == 'B':
            f1 = vert_boson_pm(v1, ps[0], ps[1])
            f2 = vert_boson_pm(v2, ps[1], np.zeros(4))
        else:
            f1, f2 = C.const_pm(v1), C.const_pm(v2)
        return chain([C.prop_pm(kind, np.zeros(4), taumax), f1, C.prop_pm(kind, ps[1], taumax), f2], taumax)
    res = {}
    TB = bub('B', a[0], a[1])
    res['A|_V loop'] = C.cscale(C.log_avg(TB, 2, PV), Fraction(1, 2))
    res['A|_Wperp loop'] = C.cscale(C.log_avg(TB, 2, PW), Fraction(1, 2))
    res['lambda loop'] = C.cscale(C.log_avg(bub('F', C.Y_of(a[0]), C.Y_of(a[1])), 2), Fraction(-1, 2))
    res['ghost loop'] = C.cscale(C.log_avg(bub('G', C.ghost_V(a[0]), C.ghost_V(a[1])), 2), Fraction(-1, 2))
    out = {k: C.cdiv(v, C.cscale(tree_c, 2)) for k, v in res.items()}
    out['TOTAL'] = C.cdiv(C.cadd(*res.values()), C.cscale(tree_c, 2))
    # longitudinal log check: pols along p
    al = []
    for _ in range(2):
        b = np.zeros((C.ND, C.NC))
        for col in range(C.NC):
            b[:C.NV, col] = p * rng.integers(-RANGE, RANGE + 1)
        al.append(b.reshape(-1).astype(complex))
    TBl = bub('B', al[0], al[1])
    lres = [C.cscale(C.log_avg(TBl, 2, PV), Fraction(1, 2)), C.cscale(C.log_avg(TBl, 2, PW), Fraction(1, 2)),
            C.cscale(C.log_avg(bub('F', C.Y_of(al[0]), C.Y_of(al[1])), 2), Fraction(-1, 2)),
            C.cscale(C.log_avg(bub('G', C.ghost_V(al[0]), C.ghost_V(al[1])), 2), Fraction(-1, 2))]
    out['longitudinal log (absolute, must be 0)'] = C.cadd(*lres)
    out['longitudinal log per loop A_V,A_W,lambda,ghost'] = [C.cstr(x) for x in lres]
    return out


def bilin_avg(P, tau_target, u1, u2):
    re = Fraction(0)
    im = Fraction(0)
    for (t, m, mono), M in P.items():
        if t != tau_target or sum(mono) - 2 * m != -4:
            continue
        w = C.ang_avg(mono)
        if w == 0:
            continue
        r, i = C.cfrac(u1 @ M @ u2)
        re += w * r
        im += w * i
    return (re, im)


def rand_spinor():
    return rng.integers(-RANGE, RANGE + 1, size=C.NF).astype(complex)


def run_lamlam():
    """lambda two-point: Sigma0 = int_k sum_A Y_A G_F(k) Y_A 2/(k - p2)^2 vs tree K = i p2.sigmabar."""
    p2 = rand_mom()
    u1, u2 = rand_spinor(), rand_spinor()
    Kt = np.kron(np.eye(C.NC), 1j * sum(p2[mu] * C.SB[mu] for mu in range(4)))
    tree_c = C.cfrac(u1 @ Kt @ u2)
    res = {}
    for name, Ms in [('internal A|_V', range(C.NV)), ('internal A|_Wperp', range(C.NV, C.ND))]:
        tot = (Fraction(0), Fraction(0))
        for M in Ms:
            for aa in range(C.NC):
                YAm = C.YA(M, aa)
                D = C.pm_scale(C.inv_sq_shift(-p2, 1, C.NF), 2.0)
                P = chain([C.const_pm(YAm), C.prop_pm('F', np.zeros(4), 1), C.const_pm(YAm), D], 1)
                tot = C.cadd(tot, bilin_avg(P, 1, u1, u2))
        res[name] = tot
    out = {k: C.cdiv(v, C.cscale(tree_c, 2)) for k, v in res.items()}
    out['TOTAL'] = C.cdiv(C.cadd(*res.values()), C.cscale(tree_c, 2))
    return out


def run_lamlamAW():
    """(lambda, lambda, A|_Wperp): tau^0 log coefficient, three classes."""
    a3 = rand_pol_W()
    u1, u2 = rand_spinor(), rand_spinor()
    YW = C.Y_of(a3)
    tree_c = C.cfrac(u1 @ YW @ u2)
    gF = C.hatG_fermion_lin()
    # class (ii)/(iii): sum_A Y_A G_F Y_W G_F Y_A D ; integrand tau^0 = sum_A Y_A hatG(q) Y_W hatG(q) Y_A 2/q^6
    res = {}
    for name, Ms in [('external A|_Wperp on lambda line, internal A|_V', range(C.NV)),
                     ('external A|_Wperp on lambda line, internal A|_Wperp', range(C.NV, C.ND))]:
        tot = (Fraction(0), Fraction(0))
        for M in Ms:
            for aa in range(C.NC):
                YAm = C.YA(M, aa)
                D = C.pm_scale(C.inv_sq_shift(np.zeros(4), 0, C.NF), 2.0)
                P = chain([C.const_pm(YAm), C.prop_pm('F', np.zeros(4), 0), C.const_pm(YW),
                           C.prop_pm('F', np.zeros(4), 0), C.const_pm(YAm), D], 0)
                tot = C.cadd(tot, bilin_avg(P, 0, u1, u2))
        res[name] = tot
    # class (i): sum_{AB} Y_A G_F(-q) Y_B [D v3(q) D]_{AB}; tau^0: v3 = Gamma3(0, -q, q) contracted with a3
    z = np.zeros(4)
    Lmu = [np.einsum('i,iwu->wu', a3, G3([z, -unit(mu), unit(mu)])) for mu in range(4)]
    YAll = [C.YA(M, aa) for M in range(C.ND) for aa in range(C.NC)]  # index A = M*3 + a = bidx
    tot = (Fraction(0), Fraction(0))
    # integrand = sum_{mu,nu} q_mu q_nu/q^6 * 4 * sum_AB (L_mu)_AB Y_A (-g_nu) Y_B ; <q_mu q_nu> = delta/4
    acc = np.zeros((C.NF, C.NF), dtype=complex)
    for mu in range(4):
        Lm = Lmu[mu]
        for A in range(C.NB):
            for B in range(C.NB):
                if Lm[A, B] != 0:
                    acc += Lm[A, B] * (YAll[A] @ (-gF[mu]) @ YAll[B])
    # <f> = 4 * (1/4) * u1 acc u2
    r, i = C.cfrac(u1 @ acc @ u2)
    res['external A|_Wperp on boson line (A|_V and A|_Wperp internal)'] = (r, i)
    out = {k: C.cdiv(v, C.cscale(tree_c, 2)) for k, v in res.items()}
    out['TOTAL'] = C.cdiv(C.cadd(*res.values()), C.cscale(tree_c, 2))
    # cross-check class (ii)+(iii) by the closed tau^0 formula (1/2) sum_{A,mu} Y_A g_mu Y_W g_mu Y_A
    chk = np.zeros((C.NF, C.NF), dtype=complex)
    for A in range(C.NB):
        for mu in range(4):
            chk += YAll[A] @ gF[mu] @ YW @ gF[mu] @ YAll[A]
    r2, i2 = C.cfrac(u1 @ chk @ u2)
    out['check (ii)+(iii) closed formula'] = C.cdiv(C.cscale((r2, i2), Fraction(1, 2)), C.cscale(tree_c, 2))
    return out


def s(d):
    return {k: (C.cstr(v) if isinstance(v, tuple) else v) for k, v in d.items()}


if __name__ == '__main__':
    report = {'clifford_ok': C.check_clifford()}
    report['AVAVAV'], report['AVAVAV_config'] = run_AVAVAV()
    report['AVAVAV'] = s(report['AVAVAV'])
    report['AV_two_point'] = s(run_AV2())
    report['lambda_two_point'] = s(run_lamlam())
    report['lam_lam_AWperp'] = s(run_lamlamAW())
    print(json.dumps(report, indent=1))
