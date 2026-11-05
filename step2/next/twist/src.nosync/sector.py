"""Single-letter Q^xi_free cohomology in a charge sector, exact (sympy over Q(i)).
Functional basis element: (slot, comp-vector key, y-monomial). y_v = null coordinates of the plane-wave momentum,
p = sum_v y_v n_v; derivative along n_w <-> 2 i y_{-w}.  Charges (functional weights): iota_n A, <n,A*>: weight(n);
lam: dual row of W[k]: -k; lam*: (W[k], lam*): +k; c, c*: 0; y_v: -v.  Twisted degree t = d + 2 w (Weyl)."""
import env, itertools, numpy as np, sympy as sp
from fractions import Fraction
import freeq as F, qa_twist as TW, ea_core as E, qa_e0 as E0
O = E.OFS; NS = E.NS
T = TW.Twist(); XI = T.xi
KAP = F.kappa_S2(XI)
V4W = [(1, 0), (-1, 0), (0, 1), (0, -1)]           # so4 vector weights (planes 0,1)
def wt5(so4=(0, 0), so6=(0, 0, 0)):
    return tuple(Fraction(x) for x in tuple(so4) + tuple(so6))
QXI = wt5((Fraction(1, 2), Fraction(-1, 2)), (Fraction(1, 2), Fraction(1, 2), Fraction(-1, 2)))
def nul10(plane, s):
    a, b = TW.CART[plane]; v = np.zeros(10, complex); v[a] = 1; v[b] = 1j * s; return v
# components with charges
COMPS = []   # (slot, key, vector(len SZ), charge, weyl, deg)
for pl in range(5):
    for s in (1, -1):
        w = [0] * 5; w[pl] = s
        COMPS.append(('A', ('A', pl, s), nul10(pl, s), wt5(w[:2], w[2:]), Fraction(1), 0))
        COMPS.append(('As', ('As', pl, s), nul10(pl, s), wt5(w[:2], w[2:]), Fraction(3), -1))
for k, vec in T.W.items():
    h = tuple(Fraction(x, 2) for x in k)
    COMPS.append(('l', ('l', k), T.dual[k], tuple(-x for x in h), Fraction(3, 2), 0))
    COMPS.append(('ls', ('ls', k), vec, h, Fraction(5, 2), -1))
COMPS.append(('c', ('c',), np.array([1 + 0j]), wt5(), Fraction(0), 1))
COMPS.append(('cs', ('cs',), np.array([1 + 0j]), wt5(), Fraction(4), -2))
YV = V4W
Y = {v: sp.Symbol('y_%d_%d' % (v[0] + 1, v[1] + 1)) for v in YV}
def ycharge(v): return wt5((-v[0], -v[1]))
NV = {v: nul10(0 if v[0] else 1, v[0] if v[0] else v[1])[:4] for v in YV}
def gq(z, den=1 << 16):
    z = complex(z); fr = Fraction(z.real).limit_denominator(den); fi = Fraction(z.imag).limit_denominator(den)
    assert abs(float(fr) - z.real) < 1e-9 and abs(float(fi) - z.imag) < 1e-9, z
    return sp.Rational(fr.numerator, fr.denominator) + sp.I * sp.Rational(fi.numerator, fi.denominator)
# symbolic M(p) in y: M is polynomial of degree <= 2 in p; get it from evaluations
def M_sym():
    def Mnum(p): return F.M_free(XI, np.asarray(p, complex), KAP)
    M0 = Mnum(np.zeros(4))
    lin = [Mnum(np.eye(4)[m]) for m in range(4)]
    quad = {}
    for m in range(4):
        for n in range(m, 4):
            pp = np.eye(4)[m] + np.eye(4)[n]
            quad[(m, n)] = Mnum(pp)
    # M(p) = M0 + sum p_m L_m + sum_{m<=n} p_m p_n Q_mn ; solve: Mnum(e_m) = M0 + L_m + Q_mm ; Mnum(2 e_m) = M0 + 2L_m + 4Q_mm
    L = {}; Qd = {}
    for m in range(4):
        M2 = Mnum(2 * np.eye(4)[m])
        Qd[m] = (M2 - 2 * lin[m] + M0) / 2
        L[m] = lin[m] - M0 - Qd[m]
    Qo = {}
    for m in range(4):
        for n in range(m + 1, 4):
            Qo[(m, n)] = quad[(m, n)] - M0 - L[m] - L[n] - Qd[m] - Qd[n]
    # check at a random point
    rng = np.random.default_rng(2); p = rng.integers(-3, 4, 4).astype(complex)
    rec = M0 + sum(p[m] * L[m] for m in range(4)) + sum(p[m] ** 2 * Qd[m] for m in range(4)) + sum(p[m] * p[n] * Qo[(m, n)] for (m, n) in Qo)
    assert np.allclose(rec, Mnum(p))
    pc = [sum(Y[v] * complex(NV[v][m]) for v in YV) for m in range(4)]
    pc = [sp.nsimplify(0) + sum(Y[v] * gq(NV[v][m]) for v in YV) for m in range(4)]
    return M0, L, Qd, Qo, pc
M0, L, Qd, Qo, PC = M_sym()
def apply_M(fvec):
    """fvec: dict comp_index(0..53) -> sympy poly; returns w M as dict"""
    out = {}
    def add(j, val):
        if val != 0: out[j] = sp.expand(out.get(j, 0) + val)
    for i, f in fvec.items():
        for j in range(NS):
            coef = gq(M0[i, j])
            for m in range(4):
                if abs(L[m][i, j]) > 1e-12: coef += gq(L[m][i, j]) * PC[m]
                if abs(Qd[m][i, j]) > 1e-12: coef += gq(Qd[m][i, j]) * PC[m] ** 2
            for (m, n), Qm in Qo.items():
                if abs(Qm[i, j]) > 1e-12: coef += gq(Qm[i, j]) * PC[m] * PC[n]
            if coef != 0: add(j, f * coef)
    return out
def comp_functional(c):
    slot, key, vec = c[0], c[1], c[2]
    w = {}
    for a, x in enumerate(vec):
        if abs(x) > 1e-12: w[O[slot] + a] = gq(x)
    return w
# decompose a raw functional (dict index -> poly) into the component basis: for each slot, the component functionals
# are dual/primal vectors; raw row r on a slot = sum_c alpha_c comp_vec_c  (solve)
SLOTCOMPS = {}
for ci, c in enumerate(COMPS):
    SLOTCOMPS.setdefault(c[0], []).append(ci)
def decompose(raw):
    res = {}
    for slot, cis in SLOTCOMPS.items():
        idx = list(range(O[slot], O[slot] + E.SZ[slot]))
        if not any(i in raw for i in idx): continue
        B = sp.Matrix([[gq(COMPS[ci][2][a]) if a < len(COMPS[ci][2]) else 0 for ci in cis] for a in range(E.SZ[slot])])
        if slot == 'A' or slot == 'As':
            B = sp.Matrix([[gq(COMPS[ci][2][a]) for ci in cis] for a in range(10)])
        rvec = sp.Matrix([raw.get(i, 0) for i in idx])
        sol = B.solve(rvec)   # alpha_c (polynomials)
        for k, ci in enumerate(cis):
            if sp.expand(sol[k]) != 0: res[ci] = sp.expand(sol[k])
    return res

def monomials(maxdeg):
    out = [()]
    for d in range(1, maxdeg + 1):
        out += list(itertools.combinations_with_replacement(YV, d))
    return out
def mcharge(m):
    ch = wt5()
    for v in m: ch = tuple(a + b for a, b in zip(ch, ycharge(v)))
    return ch
def in_lattice(diff):
    """diff = k q_xi for integer k ?"""
    k = diff[0] / QXI[0]
    if k.denominator != 1: return None
    return int(k) if all(diff[i] == k * QXI[i] for i in range(5)) else None
def sector_basis(chi, t, maxdeg=3):
    B = []
    for ci, c in enumerate(COMPS):
        for m in monomials(maxdeg):
            w = c[4] + len(m); d = c[5]
            if d + 2 * w != t: continue
            ch = tuple(a + b for a, b in zip(c[3], mcharge(m)))
            k = in_lattice(tuple(a - b for a, b in zip(ch, chi)))
            if k is None: continue
            B.append((ci, m))
    return B
def mono_expr(m):
    e = sp.Integer(1)
    for v in m: e *= Y[v]
    return e
def to_basis_vec(dec, basis):
    """dec: dict ci -> poly in y ; express as vector on basis (ci, monomial)"""
    vec = [0] * len(basis)
    pos = {b: i for i, b in enumerate(basis)}
    for ci, poly in dec.items():
        P = sp.Poly(sp.expand(poly), *[Y[v] for v in YV]) if poly != 0 else None
        if P is None: continue
        for mon, cf in P.terms():
            m = tuple(sorted(sum([[YV[i]] * e for i, e in enumerate(mon)], [])))
            m = tuple(sorted(m))
            key = None
            for b in basis:
                if b[0] == ci and tuple(sorted(b[1])) == m: key = b; break
            if key is None: raise KeyError(('not in sector', COMPS[ci][1], m, cf))
            vec[pos[key]] += cf
    return vec
def Qmatrix(basis_src, basis_tgt):
    cols = []
    for (ci, m) in basis_src:
        f = {k: v * mono_expr(m) for k, v in comp_functional(COMPS[ci]).items()}
        r = decompose(apply_M(f))
        cols.append(to_basis_vec(r, basis_tgt))
    return sp.Matrix(cols).T if cols else sp.zeros(len(basis_tgt), 0)
def cohomology(chi, t, maxdeg=3):
    Bm, B0, Bp = sector_basis(chi, t - 1, maxdeg), sector_basis(chi, t, maxdeg), sector_basis(chi, t + 1, maxdeg + 1)
    D0 = Qmatrix(B0, Bp)          # t -> t+1
    Dm = Qmatrix(Bm, B0)          # t-1 -> t
    Z = D0.nullspace()
    Bim = Dm.columnspace()
    return dict(Bm=Bm, B0=B0, Bp=Bp, D0=D0, Dm=Dm, Z=Z, Bim=Bim)
def name(b):
    ci, m = b
    return '%s%s' % (COMPS[ci][1], ''.join('*y%s' % (v,) for v in m))
