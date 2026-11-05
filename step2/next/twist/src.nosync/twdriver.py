"""Twisted pair entries of Q^xi_1 = sum over arities k = 0, 1, 2 of the two-body defect-cycle kernels, with inputs on
closed representatives of the holomorphic letters and outputs projected letter by letter onto H(Q^xi_free) (p).
Output functional of one chain order (Y1 near X_a, Y2 near X_b) = value x F[Y1, Y2] (engine convention: value =
coefficient of theta_1 theta_2 for chains visiting the outputs in leg order, colour word (1 x ad_Y2 ad_Y1)(-,-)^-1)."""
import env, itertools, json, time
import numpy as np, sympy as sp
from fractions import Fraction
import sector as S, qa_twist as TW, qa_core as QC, qa_util as U, qa2, ea_core as E
T = S.T
GX = QC.odd_gen(T.xi); G0 = QC.zero_gen()
LT_OF = {}
def legtype(ci):
    slot, key = S.COMPS[ci][0], S.COMPS[ci][1]
    if slot == 'A':  return 'AV' if key[1] < 2 else 'AW'
    if slot == 'As': return 'AsV' if key[1] < 2 else 'AsW'
    if slot == 'c':  return 'c'
    if slot == 'cs': return 'cs'
    k = key[1]
    plus = (k[0] == -k[1])          # S_+ in Sigma
    if slot == 'l':  return 'lamP' if plus else 'lamM'
    if slot == 'ls': return 'lamsP' if plus else 'lamsM'   # functional (W[k], lam*) pairs with W[k] in Sigma
def legpol(ci):
    """polarisation (letter row) on which the component functional is 1 and the other basis functionals 0"""
    slot, key, vec = S.COMPS[ci][0], S.COMPS[ci][1], S.COMPS[ci][2]
    if slot in ('A', 'As'):
        pl, s = key[1], key[2]
        nb = S.nul10(pl, -s)
        v = nb / (vec @ nb)
        return v[:4] if pl < 2 else v[4:]
    if slot == 'l':  return T.W[key[1]]
    if slot == 'ls': return T.dual[key[1]]
    return np.array([1 + 0j])
def weyl(b):  ci, m = b; return S.COMPS[ci][4] + len(m)
def deg(b):   ci, m = b; return S.COMPS[ci][5]
# ---------------------------------------------------------------- projections
def projection(chi, t, target=None, maxdeg=3):
    """returns (basis, pi) with pi[i] = coordinate of p(basis_i) along the class of `target` (a dict name->coef of a
    closed element) or along the unique class if dim H = 1.  Complement of Z chosen among basis vectors."""
    R = S.cohomology(chi, t, maxdeg)
    B0, Z, Bim = R['B0'], R['Z'], R['Bim']
    n = len(B0)
    Zm = sp.Matrix.hstack(*Z) if Z else sp.zeros(n, 0)
    Bm = sp.Matrix.hstack(*Bim) if Bim else sp.zeros(n, 0)
    # H basis: complement of B in Z
    hz = []
    cur = Bm
    for z in Z:
        test = sp.Matrix.hstack(cur, z) if cur.shape[1] else z
        if test.rank() > (cur.rank() if cur.shape[1] else 0):
            hz.append(z); cur = test
    assert len(hz) >= 1, 'empty cohomology'
    # complement C of Z among standard basis vectors
    comp = []
    cur = Zm
    for i in range(n):
        e = sp.zeros(n, 1); e[i] = 1
        test = sp.Matrix.hstack(cur, e) if cur.shape[1] else e
        if test.rank() > (cur.rank() if cur.shape[1] else 0):
            comp.append(e); cur = test
    full = sp.Matrix.hstack(*( ([Bm] if Bm.shape[1] else []) + hz + comp))
    # express each basis vector: e_i = B a + H h + C c ; pi = h-coordinates
    if target is None:
        assert len(hz) == 1, ('dim H > 1, give target', len(hz))
        tvec = hz[0]
    else:
        tvec = sp.Matrix([target.get(S.name(b), 0) for b in B0])
    sol = full.inv() if full.shape[0] == full.shape[1] else None
    assert sol is not None, ('basis not square', full.shape)
    nb = Bm.shape[1]
    H = len(hz)
    # class of target in hz coordinates
    tc = sol * tvec
    th = [sp.nsimplify(tc[nb + j]) for j in range(H)]
    assert all(sp.simplify(tc[i]) == 0 for i in range(nb + H, n)), 'target not closed'
    jt = [j for j in range(H) if th[j] != 0]
    assert len(jt) == 1 or H == 1
    pi = []
    for i in range(n):
        e = sp.zeros(n, 1); e[i] = 1
        c = sol * e
        hcoord = [c[nb + j] for j in range(H)]
        # coordinate along the target class (H = 1, or the target is one of the hz)
        j0 = jt[0]
        val = sp.nsimplify(hcoord[j0] / th[j0])
        if H > 1:
            for j in range(H):
                if j != j0 and hcoord[j] != 0:
                    pass   # other classes: recorded separately by the caller if needed
        pi.append(val)
    return B0, pi, R
# ---------------------------------------------------------------- engine evaluation
def eval_pair(La, Lb, cu, cv, k, D, log=print):
    """kernel of the input reps (La, Lb) (lists of (coef, zeta, deriv, ltype)) on output legs of components cu, cv at
    arity k, both chain orders; returns {order: {(mon_u, mon_v): value-vector}} with monomials in the y of each leg"""
    gen = {0: G0, 1: GX, 2: GX}[k]
    tu, tv = legtype(cu), legtype(cv)
    pu, pv = legpol(cu), legpol(cv)
    wl = list(S.YV)
    res = {}
    for order in ((0, 1), (1, 0)):
        def fun(x, order=order):
            q = []
            for j in range(2):
                qq = np.zeros(4, complex)
                for a, w in enumerate(wl):
                    qq = qq + x[4 * j + a] * S.NV[w]
                q.append(qq)
            legs = [(tu, pu, q[0]), (tv, pv, q[1])]
            legs = [legs[i] for i in order]
            tot = U.vzero()
            for ca, za, fa, ta in La:
                for cb, zb, fb, tb in Lb:
                    if k < 2:
                        v = QC.pair_value((ta, za, fa), (tb, zb, fb), legs, gen)['fin']
                    else:
                        v = qa2.pair_value2((ta, za, fa), (tb, zb, fb), legs, gen)['fin']
                    c = complex(ca) * complex(cb)
                    tot = U.vadd(tot, TW.vscale_c(v, c))
            return tot
        poly = TW.interpolate(fun, D, 8)
        out = {}
        for e, v in poly.items():
            mu = tuple(sorted(sum([[wl[a]] * e[a] for a in range(4)], [])))
            mv = tuple(sorted(sum([[wl[a]] * e[4 + a] for a in range(4)], [])))
            out[(mu, mv)] = v
        res[order] = out
    return res
def vec_to_sym(v):
    """5 Gaussian rationals on (1, log2, log3, log5, log7) -> sympy"""
    L2, L3, L5, L7 = sp.log(2), sp.log(3), sp.log(5), sp.log(7)
    basis = [1, L2, L3, L5, L7]
    return sp.nsimplify(sum((sp.Rational(a.numerator, a.denominator) + sp.I * sp.Rational(b.numerator, b.denominator)) * basis[i]
                            for i, (a, b) in enumerate(v)))

def assemble(La, Lb, w_in, d_in, secU, secV, ks=(0, 1, 2), log=print, cache=None):
    """secU = (chi, t, target-or-None).  Returns {'UX': coefficient of F[U near X_a, V near X_b],
    'VX': coefficient of F[V near X_a, U near X_b]} (sympy), with per-arity parts."""
    BU, piU, _ = projection(*secU)
    BV, piV, _ = projection(*secV)
    groups = {}
    for i, u in enumerate(BU):
        if piU[i] == 0: continue
        for j, v in enumerate(BV):
            if piV[j] == 0: continue
            k2 = 2 * (weyl(u) + weyl(v) - w_in)
            if k2.denominator != 1: continue
            k = int(k2)
            if k not in ks or deg(u) + deg(v) != d_in + 1 - k: continue
            Dg = len(u[1]) + len(v[1])
            groups.setdefault((u[0], v[0], k, Dg), []).append((i, j))
    out = {'UX': 0, 'VX': 0, 'parts': {}}
    for (cu, cv, k, Dg), ijs in sorted(groups.items(), key=lambda x: x[0][2]):
        t0 = time.time()
        res = eval_pair(La, Lb, cu, cv, k, Dg)
        for (i, j) in ijs:
            mu, mv = tuple(sorted(BU[i][1])), tuple(sorted(BV[j][1]))
            w = piU[i] * piV[j]
            a = vec_to_sym(res[(0, 1)].get((mu, mv), U.vzero()))
            b = vec_to_sym(res[(1, 0)].get((mu, mv), U.vzero()))
            out['UX'] += w * a; out['VX'] += w * b
            key = 'k%d:%s|%s' % (k, S.name(BU[i]), S.name(BV[j]))
            out['parts'][key] = (sp.nsimplify(w), sp.nsimplify(a), sp.nsimplify(b))
        log('   arity %d legs %s x %s (deg %d): %d monomial pairs, %.1fs' % (k, S.COMPS[cu][1], S.COMPS[cv][1], Dg, len(ijs), time.time() - t0))
    out['UX'] = sp.nsimplify(sp.expand(out['UX'])); out['VX'] = sp.nsimplify(sp.expand(out['VX']))
    return out
# ---------------------------------------------------------------- letters and sectors
from fractions import Fraction as _Fr
def _h(k): return tuple(_Fr(x, 2) for x in k)
GCH = [S.wt5((0, 0), tuple(_Fr(a + b, 2) for a, b in zip(t, T.z))) for t in T.taus]
BCH = [tuple(-x for x in _h(T.u + t)) for t in T.taus]
D1, D2 = (0, 1), (-1, 0)
def vsum(*cs): return tuple(sum(x) for x in zip(*cs))
CHI_b = S.wt5((-1, 1))
UPZ = (-1, 1) + T.z     # key of W[(u'; z)]
def L_c(): return [(1, np.array([1 + 0j]), None, 'c')]
def L_gamma(i): return [(1, T.gam_w[i], None, 'AW')]
def L_beta(i): return [(1, T.beta_z[i], None, 'lamP')]
def L_b(which=('F', 'ls')):
    out = []
    if 'F' in which: out += [(0.5 * c, z, f, 'AV') for (c, z, f) in TW.b_letters(T)]
    if 'ls' in which: out += [(1, T.W[UPZ], None, 'lamsP')]
    return out
def sec_dc(d):   return (S.wt5(d), 3, None)
def sec_dgamma(i, d): return (vsum(GCH[i], S.wt5(d)), 4, None)
def sec_dbeta(i, d):  return (vsum(BCH[i], S.wt5(d)), 5, None)
def sec_db(d):   return (vsum(CHI_b, S.wt5(d)), 6, None)

# ---------------------------------------------------------------- per-term arity (inputs with mixed weights)
WEYL_T = dict(E.WEYL); DEG_T = dict(E.DEG)
def term_wd(term):
    c, z, f, lt = term
    nd = 0 if f is None else max(sum(e) for e in f)
    return WEYL_T[lt] + nd, DEG_T[lt]
def eval_terms(ta, tb, cu, cv, k, Dg):
    return eval_pair([ta], [tb], cu, cv, k, Dg)
def assemble2(La, Lb, secU, secV, ks=(0, 1, 2), log=print):
    BU, piU, _ = projection(*secU)
    BV, piV, _ = projection(*secV)
    out = {'UX': 0, 'VX': 0, 'parts': {}}
    for ia, ta in enumerate(La):
        wa, da = term_wd(ta)
        for ib, tb in enumerate(Lb):
            wb, db = term_wd(tb)
            groups = {}
            for i, u in enumerate(BU):
                if piU[i] == 0: continue
                for j, v in enumerate(BV):
                    if piV[j] == 0: continue
                    k2 = 2 * (weyl(u) + weyl(v) - wa - wb)
                    if k2.denominator != 1: continue
                    k = int(k2)
                    if k not in ks or deg(u) + deg(v) != da + db + 1 - k: continue
                    groups.setdefault((u[0], v[0], k, len(u[1]) + len(v[1])), []).append((i, j))
            for (cu, cv, k, Dg), ijs in sorted(groups.items(), key=lambda x: x[0][2]):
                t0 = time.time()
                res = eval_pair([ta], [tb], cu, cv, k, Dg)
                for (i, j) in ijs:
                    mu, mv = tuple(sorted(BU[i][1])), tuple(sorted(BV[j][1]))
                    w = piU[i] * piV[j]
                    a = vec_to_sym(res[(0, 1)].get((mu, mv), U.vzero()))
                    b = vec_to_sym(res[(1, 0)].get((mu, mv), U.vzero()))
                    out['UX'] += w * a; out['VX'] += w * b
                    key = 'k%d:t%d,%d:%s|%s' % (k, ia, ib, S.name(BU[i]), S.name(BV[j]))
                    out['parts'][key] = (sp.nsimplify(w), sp.nsimplify(a), sp.nsimplify(b))
                log('   terms (%d,%d) arity %d legs %s x %s (deg %d): %.1fs' % (ia, ib, k, S.COMPS[cu][1], S.COMPS[cv][1], Dg, time.time() - t0))
    out['UX'] = sp.nsimplify(sp.expand(out['UX'])); out['VX'] = sp.nsimplify(sp.expand(out['VX']))
    return out

# ---------------------------------------------------------------- targeted interpolation (any degree, few variables)
def eval_pair_targeted(ta, tb, cu, cv, k, monos, log=None):
    """monos: list of (mu, mv) monomials needed (y-multisets for leg u and leg v); homogeneous of one total degree.
    Restricts the leg momenta to the null directions occurring in monos and solves exactly for all monomials of that
    degree in those variables.  Returns {order: {(mu, mv): value}}."""
    gen = {0: G0, 1: GX, 2: GX}[k]
    tu, tv = legtype(cu), legtype(cv)
    pu, pv = legpol(cu), legpol(cv)
    varsU = sorted(set(v for mu, mv in monos for v in mu)); varsV = sorted(set(v for mu, mv in monos for v in mv))
    Dg = len(monos[0][0]) + len(monos[0][1])
    assert all(len(a) + len(b) == Dg for a, b in monos)
    nv = len(varsU) + len(varsV)
    allm = [m for m in itertools.product(range(Dg + 1), repeat=nv) if sum(m) == Dg]
    rng = np.random.default_rng(11)
    pts = []
    if nv == 0:
        pts = [()]
    while nv > 0 and len(pts) < len(allm) + 2:
        p = tuple(int(x) for x in rng.integers(-3, 4, nv))
        if p not in pts: pts.append(p)
    A = sp.Matrix([[sp.prod([sp.Integer(p[i]) ** m[i] for i in range(nv)]) for m in allm] for p in pts])
    res = {}
    for order in ((0, 1), (1, 0)):
        vals = []
        for p in pts:
            qU = sum((p[i] * S.NV[varsU[i]] for i in range(len(varsU))), np.zeros(4, complex))
            qV = sum((p[len(varsU) + i] * S.NV[varsV[i]] for i in range(len(varsV))), np.zeros(4, complex))
            legs = [(tu, pu, qU), (tv, pv, qV)]
            legs = [legs[i] for i in order]
            if k < 2:
                v = QC.pair_value((ta[3], ta[1], ta[2]), (tb[3], tb[1], tb[2]), legs, gen)['fin']
            else:
                v = qa2.pair_value2((ta[3], ta[1], ta[2]), (tb[3], tb[1], tb[2]), legs, gen)['fin']
            c = complex(ta[0]) * complex(tb[0])
            vals.append(TW.vscale_c(v, c))
        # solve A x = vals for each of the 5 basis components (real and imaginary parts)
        sol = {}
        for comp in range(5):
            bvec = sp.Matrix([sp.Rational(v[comp][0].numerator, v[comp][0].denominator) + sp.I * sp.Rational(v[comp][1].numerator, v[comp][1].denominator) for v in vals])
            x = (A.T * A).LUsolve(A.T * bvec)
            assert sp.simplify(A * x - bvec) == sp.zeros(len(pts), 1), 'not homogeneous of degree %d' % Dg
            for mi, m in enumerate(allm):
                if x[mi] != 0: sol.setdefault(m, [0] * 5)[comp] = x[mi]
        out = {}
        for m, comps in sol.items():
            mu = tuple(sorted(sum([[varsU[i]] * m[i] for i in range(len(varsU))], [])))
            mv = tuple(sorted(sum([[varsV[i]] * m[len(varsU) + i] for i in range(len(varsV))], [])))
            out[(mu, mv)] = comps
        res[order] = out
    return res
L2_, L3_, L5_, L7_ = sp.log(2), sp.log(3), sp.log(5), sp.log(7)
def comps_to_sym(cs): return sp.nsimplify(sum(c * b for c, b in zip(cs, [1, L2_, L3_, L5_, L7_])))
def projection_rev(chi, t, target=None, maxdeg=3):
    """projection with the complement of Z chosen from the END of the basis (alternative p, homotopy-data check)"""
    R = S.cohomology(chi, t, maxdeg)
    B0 = R['B0']; n = len(B0)
    order = list(range(n))[::-1]
    perm = sp.Matrix([[1 if order[j] == i else 0 for j in range(n)] for i in range(n)])
    Z = R['Z']; Bim = R['Bim']
    Zm = sp.Matrix.hstack(*Z); Bm = sp.Matrix.hstack(*Bim) if Bim else sp.zeros(n, 0)
    hz = [z for z in Z if (sp.Matrix.hstack(Bm, z).rank() if Bm.shape[1] else z.rank()) > (Bm.rank() if Bm.shape[1] else 0)][:1]
    comp = []; cur = Zm
    for i in order:
        e = sp.zeros(n, 1); e[i] = 1
        test = sp.Matrix.hstack(cur, e)
        if test.rank() > cur.rank(): comp.append(e); cur = test
    full = sp.Matrix.hstack(*(([Bm] if Bm.shape[1] else []) + hz + comp))
    sol = full.inv(); nb = Bm.shape[1]
    pi = [sp.nsimplify((sol * sp.Matrix([1 if j == i else 0 for j in range(n)]))[nb]) for i in range(n)]
    return B0, pi, R
def assemble3(La, Lb, secU, secV, ks=(0, 1, 2), log=print, rev=(False, False)):
    BU, piU, _ = (projection_rev if rev[0] else projection)(*secU)
    BV, piV, _ = (projection_rev if rev[1] else projection)(*secV)
    out = {'UX': 0, 'VX': 0, 'parts': {}}
    for ia, ta in enumerate(La):
        wa, da = term_wd(ta)
        for ib, tb in enumerate(Lb):
            wb, db = term_wd(tb)
            groups = {}
            for i, u in enumerate(BU):
                if piU[i] == 0: continue
                for j, v in enumerate(BV):
                    if piV[j] == 0: continue
                    k2 = 2 * (weyl(u) + weyl(v) - wa - wb)
                    if k2.denominator != 1: continue
                    k = int(k2)
                    if k not in ks or deg(u) + deg(v) != da + db + 1 - k: continue
                    groups.setdefault((u[0], v[0], k, len(u[1]) + len(v[1])), []).append((i, j))
            for (cu, cv, k, Dg), ijs in sorted(groups.items(), key=lambda x: x[0][2]):
                t0 = time.time()
                monos = [(tuple(sorted(BU[i][1])), tuple(sorted(BV[j][1]))) for (i, j) in ijs]
                res = eval_pair_targeted(ta, tb, cu, cv, k, monos)
                for (i, j), (mu, mv) in zip(ijs, monos):
                    w = piU[i] * piV[j]
                    a = comps_to_sym(res[(0, 1)].get((mu, mv), [0] * 5)); b = comps_to_sym(res[(1, 0)].get((mu, mv), [0] * 5))
                    out['UX'] += w * a; out['VX'] += w * b
                    out['parts']['k%d:t%d,%d:%s|%s' % (k, ia, ib, S.name(BU[i]), S.name(BV[j]))] = (sp.nsimplify(w), sp.nsimplify(a), sp.nsimplify(b))
                log('   terms (%d,%d) arity %d legs %s x %s (deg %d, %d pairs): %.1fs' % (ia, ib, k, S.COMPS[cu][1], S.COMPS[cv][1], Dg, len(ijs), time.time() - t0))
    out['UX'] = sp.nsimplify(sp.expand(out['UX'])); out['VX'] = sp.nsimplify(sp.expand(out['VX']))
    return out
