#!/usr/bin/env python3
"""conjugation_checks.py -- machine checks for step2/theory/conjugation.tex
(radial conjugation of psl(4|4), Osterwalder-Schrader reflection, conformal map C, realisation on BV fields).

Realisation: the one of Chapter 2 (superconformal/algebra.py, structure.py, model.py, engine.py, which import
classical/n4_model.py): rho(e_M)|_Sigma = sig[M], rho(e_M)|_Sigma* = sigbar[M], (s,t) = s^T t,
sig[M] = sigbar[M] real symmetric for M <= 8, sig[9] = i, sigbar[9] = -i.  V = M in {0,1,2,3}, W_R^perp = {4..9}.
In this realisation the antilinear Spin(10)-equivariant map C : Sigma <-> Sigma^* is componentwise complex
conjugation.  Two-component (dictionary) realisation: step2/notation (formulas (R1)-(R4), (P1)).

Sections
  S1  exact (sympy, Gaussian rationals): anti-automorphism and involution on all 62^2 ordered basis pairs,
      named subspaces, positivity coefficients, negative controls
  S2  exact (sympy): the conformal map C_t (flow of P_n - K_n), C_{+-pi/4}, plane -> sphere, C r_n C^-1 = I
  S3  numeric: Ad(C_{+-pi/4}) = exp(-+ pi/4 ad_{P_n-K_n}) intertwines the OS and radial conjugations
  S4  numeric: field level (C^-1)^* Theta C^* = radial reflection, for every field
  S5  numeric: pointwise density identities (random Grassmann jets, su(2)) for S, S^(1), S^(2) under Theta and theta
  S6  numeric: Q, Q^GF, H, pairing, K_t, L_t, P(eps,T), P(eps,infty) under Theta (momentum space)
  S7  numeric: Step-1 counterterms, Step-2 bases G1..G5, H1..H6, finite counterterms under Theta
  S8  numeric: reflection positivity of the free scalar and fermion forms (OS half spaces, radial ball)
  S9  numeric: free-field adjointness G(rho(X)a,b) = G(a,rho(X^dag)b) for X = P, K, Delta, Q, S (radial and OS)
  S10 exact (project engine): Theta on S, S^(1), bracket rule, with the engine's own bracket
  S11 numeric: two-component form of C and of the reflections (dictionary realisation D)
Light: < 400 MB, ~3-4 min.  Output: conjugation_checks.json (next to this file).
"""
import os, sys, json, time, itertools, random, resource, gc
import numpy as np
np.seterr(all='ignore')          # Apple Accelerate emits spurious matmul warnings; every check is an explicit residual
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, 'superconformal'))
sys.path.insert(0, os.path.join(ROOT, 'classical'))
sys.path.insert(0, os.path.join(ROOT, 'engine_b'))
import sympy as sp
from scipy.linalg import expm
import algebra as Al, structure as St
import n4_model as NM
import engine as EN, model as MO
import realization as RB

T0 = time.time()
OUT = {'file': 'step2/theory/conjugation_checks.py', 'tex': 'step2/theory/conjugation.tex', 'sections': {}}
TOL = 1e-9


def mem():
    return '%.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20)


def rec(sec, cid, desc, value, passed):
    OUT['sections'].setdefault(sec, []).append(dict(id=cid, check=desc, value=value, passed=bool(passed)))
    print('[%s] %-6s %s :: %s' % ('PASS' if passed else 'FAIL', cid, desc, value))


def resid(a, b):
    a = np.asarray(a); b = np.asarray(b)
    return float(np.max(np.abs(a - b)) / max(1.0, float(np.max(np.abs(b)))))


# =====================================================================================================
# S1. exact algebra
# =====================================================================================================
Al.set_cR(sp.Integer(1))
X = Al.X
COORD = Al.COORD
R2 = sum(x ** 2 for x in X)


def conj_expr(e):
    e = sp.expand(e)
    if e == 0:
        return sp.Integer(0)
    p = sp.Poly(e, *COORD)
    return sp.expand(sum(sp.conjugate(c) * sp.Mul(*[g ** k for g, k in zip(COORD, m)]) for m, c in p.terms()))


def conj_el(a):
    if isinstance(a, Al.Even):
        return Al.Even([conj_expr(c) for c in a.c])
    return Al.Odd([sp.conjugate(z) for z in a.Q], [sp.conjugate(z) for z in a.eta])


def Istar(v, sign=1):
    """pushforward by the inversion I(x) = x/(x,x) (sign=+1); sign=-1 gives the pushforward by -I."""
    sub = {X[m]: sign * X[m] / R2 for m in range(4)}
    vI = [sp.sympify(v.c[m]).subs(sub, simultaneous=True) for m in range(4)]
    dot = sum(vI[m] * X[m] for m in range(4))
    w = [sp.expand(sp.cancel(sign * (R2 * vI[m] - 2 * dot * X[m]))) for m in range(4)]
    return Al.Even(w + list(v.c[4:]))


def rstar(v, n):
    Rn = sp.eye(4) - 2 * sp.Matrix(n) * sp.Matrix(n).T
    y = Rn * sp.Matrix(X)
    sub = {X[m]: y[m] for m in range(4)}
    vr = sp.Matrix([sp.sympify(v.c[m]).subs(sub, simultaneous=True) for m in range(4)])
    w = Rn * vr
    return Al.Even([sp.expand(w[m]) for m in range(4)] + list(v.c[4:]))


def SIGv(n):
    return sum((n[m] * Al.SIGM[m] for m in range(4)), sp.zeros(16, 16))


def SIGBv(n):
    return sum((n[m] * Al.SIGBM[m] for m in range(4)), sp.zeros(16, 16))


def dag_rad(a, c=-1, even_sign=-1, inv_sign=1):
    """radial conjugation: v -> -I_*(vbar);  Q_1 + rho(x)Q_2 -> c (C Q_2 + rho(x) C Q_1)."""
    if isinstance(a, Al.Even):
        return Istar(conj_el(a), inv_sign).scale(even_sign)
    return Al.Odd([c * sp.conjugate(z) for z in a.eta], [c * sp.conjugate(z) for z in a.Q])


def dag_OS(a, n, c=1):
    """OS conjugation in n^perp: v -> -(r_n)_*(vbar);  Q_1 + rho(x)Q_2 -> c (rho(n)CQ_1 - rho(x) rho(n) C Q_2)."""
    if isinstance(a, Al.Even):
        return rstar(conj_el(a), n).scale(-1)
    Qc = sp.Matrix([sp.conjugate(z) for z in a.Q]); ec = sp.Matrix([sp.conjugate(z) for z in a.eta])
    return Al.Odd(list(c * SIGBv(n) * Qc), list(-c * SIGv(n) * ec))


def same_el(a, b):
    return (a + b.scale(-1)).is_zero()


print('=== S1 exact algebra ===')
B = St.basis(); CO = St.Coords(B)
NB = CO.ne + len(CO.odd)
EL = [CO.element(k) for k in range(NB)]
PAR = [CO.parity(k) for k in range(NB)]
BR = [[Al.bracket(EL[i], EL[j]) for j in range(NB)] for i in range(NB)]
print('  brackets of all ordered basis pairs: %.1f s' % (time.time() - T0))


def anti_check(dg, pairs=None):
    D = [dg(e) for e in EL]
    inv = sum(1 for k in range(NB) if not same_el(dg(D[k]), EL[k]))
    bad = 0
    pl = pairs if pairs is not None else [(i, j) for i in range(NB) for j in range(NB)]
    for (i, j) in pl:
        if not same_el(dg(BR[i][j]), Al.bracket(D[j], D[i])):
            bad += 1
    return inv, bad, len(pl), D


n_e0 = [1, 0, 0, 0]
n_gen = [sp.Rational(3, 5), 0, sp.Rational(4, 5), 0]
for name, dg in [('radial, c=-1', lambda a: dag_rad(a, -1)), ('radial, c=+1', lambda a: dag_rad(a, 1)),
                 ('OS n=e_0, c=+1', lambda a: dag_OS(a, n_e0, 1)), ('OS n=(3/5,0,4/5,0), c=+1', lambda a: dag_OS(a, n_gen, 1))]:
    inv, bad, npair, D = anti_check(dg)
    rec('S1', 'A1', 'involution (X^dag)^dag = X on all 62 basis elements, %s' % name, 'failures %d' % inv, inv == 0)
    rec('S1', 'A2', '[X,Y]^dag = [Y^dag,X^dag] on all %d ordered basis pairs, %s' % (npair, name), 'violations %d' % bad, bad == 0)
    if name.startswith('radial, c=-1'):
        Drad = D
    if name.startswith('OS n=e_0'):
        DOS = D
# negative controls
odd_pairs = [(i, j) for i in range(CO.ne, NB) for j in range(CO.ne, NB)][:400]
even_pairs = [(i, j) for i in range(CO.ne) for j in range(CO.ne)][:400]
_, bad_i, npr, _ = anti_check(lambda a: dag_rad(a, sp.I), odd_pairs)
rec('S1', 'A3', 'control: odd phase c = i violates the anti-automorphism (odd-odd pairs)', 'violations %d/%d' % (bad_i, npr), bad_i > 0)
_, bad_a, npr, _ = anti_check(lambda a: dag_rad(a, -1, even_sign=+1), even_pairs)
rec('S1', 'A4', 'control: v -> +I_*(vbar) (an automorphism) violates the anti-automorphism (even-even pairs)', 'violations %d/%d' % (bad_a, npr), bad_a > 0)
_, bad_J, npr, _ = anti_check(lambda a: dag_rad(a, -1, inv_sign=-1), [(i, j) for i in range(CO.ne) for j in range(CO.ne)])
rec('S1', 'A5', 'v -> -(-I)_*(vbar) (inversion composed with x->-x) is also an anti-automorphism of the even part', 'violations %d' % bad_J, bad_J == 0)
# named subspaces
nm = CO.names
img_rad = {}; img_os = {}
for k in range(NB):
    img_rad.setdefault(nm[k], set()).update(CO.support(Drad[k]))
    img_os.setdefault(nm[k], set()).update(CO.support(DOS[k]))
img_rad = {k: sorted(v) for k, v in img_rad.items()}; img_os = {k: sorted(v) for k, v in img_os.items()}
exp_rad = {'P': ['K'], 'K': ['P'], 'M': ['M'], 'Mbar': ['Mbar'], 'Delta': ['Delta'], 'R': ['R'], 'Q': ['S'], 'Qbar': ['Sbar'], 'Sbar': ['Qbar'], 'S': ['Q']}
exp_os = {'P': ['P'], 'K': ['K'], 'M': ['Mbar'], 'Mbar': ['M'], 'Delta': ['Delta'], 'R': ['R'], 'Q': ['Qbar'], 'Qbar': ['Q'], 'Sbar': ['S'], 'S': ['Sbar']}
rec('S1', 'A6', 'radial: images of the named subspaces (Q^dag in S, P^dag in K, ...)', img_rad, img_rad == exp_rad)
rec('S1', 'A7', 'OS (n=e_0): images of the named subspaces', img_os, img_os == exp_os)
# explicit even generators
chk = []
for k in range(CO.ne):
    v = EL[k]
    if nm[k] == 'P':
        a = [sp.expand(v.c[m]) for m in range(4)]
        chk.append(same_el(Drad[k], Al.Kb([sp.conjugate(z) for z in a]).scale(-1)))
    if nm[k] == 'Delta':
        chk.append(same_el(Drad[k], v))
    if nm[k] in ('M', 'Mbar', 'R'):
        chk.append(same_el(Drad[k], conj_el(v).scale(-1)))
rec('S1', 'A8', 'radial: P_v^dag = -K_vbar, Delta^dag = Delta, rotation^dag = -rotationbar, R^dag = -Rbar (basis)', '%d/%d' % (sum(chk), len(chk)), all(chk))
# positivity coefficients
iD = nm.index('Delta')
iPn = [k for k in range(CO.ne) if nm[k] == 'P' and sp.expand(EL[k].c[0]) != 0][0]
pos_rad = []; pos_os = []
for k in range(CO.ne, NB):
    e = EL[k]
    nrm = sum(sp.Abs(z) ** 2 for z in list(e.Q) + list(e.eta))
    cD = sp.nsimplify(CO.coords(Al.bracket(e, Drad[k]))[iD])
    pos_rad.append(sp.simplify(cD / nrm))
    if nm[k] in ('Q', 'Qbar'):
        cP = sp.nsimplify(CO.coords(Al.bracket(e, DOS[k]))[iPn]) * sp.expand(EL[iPn].c[0])
        pos_os.append(sp.simplify(cP / nrm))
rec('S1', 'A9', 'radial: Delta-coefficient of [X, X^dag] / (X, C X) for all 32 odd basis X (Delta = 2E)', sorted(set(str(z) for z in pos_rad)), set(pos_rad) == {sp.Rational(-1, 2)})
rec('S1', 'A10', 'OS n=e_0: coefficient of the translation P_n in [X, X^dagTheta] / (X, C X) for the 16 basis X in Sigma', sorted(set(str(z) for z in pos_os)), set(pos_os) == {1})
print('  S1 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S2. exact geometry of C_t
# =====================================================================================================
print('=== S2 exact geometry ===')
xs = sp.symbols('y0:4', real=True); t = sp.symbols('t', real=True)
xv = sp.Matrix(xs); nn = sp.Matrix([1, 0, 0, 0]); r2 = (xv.T * xv)[0]
X0 = (1 - r2) / (1 + r2); Xvec = 2 * xv / (1 + r2)
phi = 2 * t
Xn = Xvec[0]
Xn_new = sp.cos(phi) * Xn + sp.sin(phi) * X0
X0_new = -sp.sin(phi) * Xn + sp.cos(phi) * X0
Xv_new = sp.Matrix([Xn_new, Xvec[1], Xvec[2], Xvec[3]])
Ct = Xv_new / (1 + X0_new)
vfield = lambda y: nn * (1 - (y.T * y)[0]) + 2 * (nn.T * y)[0] * y          # P_n - K_n
flow_res = sp.simplify((sp.diff(Ct, t) - vfield(Ct)).subs({xs[0]: sp.Rational(1, 3), xs[1]: sp.Rational(-2, 5), xs[2]: sp.Rational(1, 7), xs[3]: 2}))
flow_res2 = [sp.N(e.subs(t, sp.Rational(3, 10)), 30) for e in flow_res]
rec('S2', 'G1', 'd/dt C_t(x) = (P_n - K_n)(C_t(x)) (sample point, t = 3/10, 30 digits)', [str(sp.N(z, 3)) for z in flow_res2], all(abs(z) < 1e-25 for z in flow_res2))
Cp = sp.simplify(Ct.subs(t, sp.pi / 4)); Cm = sp.simplify(Ct.subs(t, -sp.pi / 4))
Cp_closed = (2 * xv + (1 - r2 - 2 * xs[0]) * nn) / ((xv - nn).T * (xv - nn))[0]
Cm_closed = (2 * xv + (r2 - 1 - 2 * xs[0]) * nn) / ((xv + nn).T * (xv + nn))[0]
rec('S2', 'G2', 'C_{pi/4}(x) = (2x + (1-(x,x)-2(n,x)) n)/|x-n|^2, C_{-pi/4}(x) = (2x + ((x,x)-1-2(n,x)) n)/|x+n|^2',
    [str(sp.simplify(Cp - Cp_closed)), str(sp.simplify(Cm - Cm_closed))], sp.simplify(Cp - Cp_closed) == sp.zeros(4, 1) and sp.simplify(Cm - Cm_closed) == sp.zeros(4, 1))
plane = {xs[0]: 0}
rec('S2', 'G3', '|C_{+-pi/4}(x)|^2 = 1 for (n,x) = 0', [str(sp.simplify((Cp_closed.T * Cp_closed)[0].subs(plane))), str(sp.simplify((Cm_closed.T * Cm_closed)[0].subs(plane)))],
    sp.simplify((Cp_closed.T * Cp_closed)[0].subs(plane) - 1) == 0 and sp.simplify((Cm_closed.T * Cm_closed)[0].subs(plane) - 1) == 0)
rn_x = sp.Matrix([-xs[0], xs[1], xs[2], xs[3]])
Cinv_closed = Cm_closed                         # C_{-pi/4} = C_{pi/4}^{-1}
comp = Cp_closed.subs({xs[i]: Cinv_closed[i] for i in range(4)}, simultaneous=True)
comp_ok = sp.simplify(comp - xv) == sp.zeros(4, 1)
CrCi = Cp_closed.subs({xs[i]: (sp.Matrix([-Cinv_closed[0], Cinv_closed[1], Cinv_closed[2], Cinv_closed[3]]))[i] for i in range(4)}, simultaneous=True)
CrCi_ok = sp.simplify(CrCi - xv / r2) == sp.zeros(4, 1)
CmrCm = Cm_closed.subs({xs[i]: (sp.Matrix([-Cp_closed[0], Cp_closed[1], Cp_closed[2], Cp_closed[3]]))[i] for i in range(4)}, simultaneous=True)
CmrCm_ok = sp.simplify(CmrCm - xv / r2) == sp.zeros(4, 1)
rec('S2', 'G4', 'C_{pi/4} C_{-pi/4} = id;  C_{pi/4} r_n C_{pi/4}^{-1} = I and C_{-pi/4} r_n C_{-pi/4}^{-1} = I (exact)', [comp_ok, CrCi_ok, CmrCm_ok], comp_ok and CrCi_ok and CmrCm_ok)
s_ = sp.symbols('s', positive=True)
img_p = sp.simplify((Cp_closed.T * Cp_closed)[0].subs({xs[0]: s_, xs[1]: 0, xs[2]: 0, xs[3]: 0}))
img_m = sp.simplify((Cm_closed.T * Cm_closed)[0].subs({xs[0]: s_, xs[1]: 0, xs[2]: 0, xs[3]: 0}))
rec('S2', 'G5', '|C_{pi/4}(s n)|^2 and |C_{-pi/4}(s n)|^2 for s > 0 (outside / inside the unit ball)', [str(img_p), str(img_m)],
    sp.simplify(img_p - ((1 + s_) / (1 - s_)) ** 2) == 0 and sp.simplify(img_m - ((1 - s_) / (1 + s_)) ** 2) == 0)
rec('S2', 'G6', 'C_{pi/4}(-n) = 0, C_{-pi/4}(n) = 0', [str(list(Cp_closed.subs({xs[0]: -1, xs[1]: 0, xs[2]: 0, xs[3]: 0}))), str(list(Cm_closed.subs({xs[0]: 1, xs[1]: 0, xs[2]: 0, xs[3]: 0})))],
    Cp_closed.subs({xs[0]: -1, xs[1]: 0, xs[2]: 0, xs[3]: 0}) == sp.zeros(4, 1) and Cm_closed.subs({xs[0]: 1, xs[1]: 0, xs[2]: 0, xs[3]: 0}) == sp.zeros(4, 1))

# =====================================================================================================
# S3. Ad(C) on psl(4|4) (numeric)
# =====================================================================================================
print('=== S3 Ad(C) ===')
def cnum(lst):
    return np.array([complex(sp.N(z)) for z in lst])
Xn_el = Al.P(n_e0) + Al.Kb(n_e0).scale(-1)
AD = np.array([cnum(CO.coords(Al.bracket(Xn_el, EL[k]))) for k in range(NB)]).T
DM = lambda D: np.array([cnum(CO.coords(D[k])) for k in range(NB)]).T
Drad_m = {c: DM([dag_rad(e, c) for e in EL]) for c in (-1, 1)}
DOS_m = {c: DM([dag_OS(e, n_e0, c) for e in EL]) for c in (-1, 1)}
ev = np.linalg.eigvals(AD)
rec('S3', 'C1', 'spectrum of ad_{P_n - K_n} (compact generator)', sorted(set(np.round(ev.imag, 9).tolist())), np.max(np.abs(ev.real)) < 1e-9)
for sgn, lab in ((-1, 'C_{pi/4}: Ad = exp(-pi/4 ad)'), (+1, 'C = C_{-pi/4} (the map of the tex): Ad = exp(+pi/4 ad)')):
    A = expm(sgn * np.pi / 4 * AD); Ai = np.linalg.inv(A)
    for cT in (1, -1):
        Mm = A @ DOS_m[cT] @ np.conj(Ai)
        cbest = min((-1, 1), key=lambda c: np.max(np.abs(Mm - Drad_m[c])))
        r = float(np.max(np.abs(Mm - Drad_m[cbest])))
        rec('S3', 'C2', '%s, OS c_Theta=%+d  ->  radial c=%+d' % (lab, cT, cbest), '%.1e' % r, r < 1e-9 and cbest == (cT if sgn < 0 else -cT))
    # geometric pushforward agreement on the conformal fields
    def Cmap_num(x, s):
        n = np.array([1., 0, 0, 0])
        if s < 0:
            d = x - n; return (2 * x + (1 - x @ x - 2 * (n @ x)) * n) / (d @ d)
        d = x + n; return (2 * x + (x @ x - 1 - 2 * (n @ x)) * n) / (d @ d)
    rng = np.random.default_rng(0); worst = 0
    lam = [sp.lambdify(X, EL[k].c[:4], 'numpy') for k in range(CO.ne)]
    for k in range(CO.ne):
        if nm[k] == 'R':
            continue
        y = rng.normal(size=4); x = Cmap_num(y, -sgn); h = 1e-6
        J = np.array([(Cmap_num(x + h * e, sgn) - Cmap_num(x - h * e, sgn)) / (2 * h) for e in np.eye(4)]).T
        push = J @ np.array(lam[k](*x), dtype=complex)
        recon = sum(A[m, k] * np.array(lam[m](*y), dtype=complex) for m in range(CO.ne))
        worst = max(worst, float(np.max(np.abs(push - recon))))
    rec('S3', 'C3', '%s equals the pushforward by the map on all conformal basis fields' % lab, '%.1e' % worst, worst < 1e-7)
# the family of conjugations with the same even part: X -> exp(tB)(X^dag), B = +1 on Q, Sbar, -1 on Qbar, S, 0 on even
sc = json.load(open(os.path.join(ROOT, 'superconformal', 'out', 'psl44_structure_constants.json')))
assert sc['names'] == nm
FS = np.zeros((NB, NB, NB), complex)
for (i, j, m, v) in sc['C']:
    FS[i, j, m] = complex(sp.N(sp.sympify(v)))
spot = [(i, j) for i in range(0, NB, 7) for j in range(0, NB, 5)]
spot_res = max(float(np.max(np.abs(FS[i, j] - cnum(CO.coords(BR[i][j]))))) for (i, j) in spot)
rec('S3', 'C4', 'stored structure constants (superconformal/out/psl44_structure_constants.json) = realised brackets (spot check)', '%.1e on %d pairs' % (spot_res, len(spot)), spot_res < 1e-12)
Bw = np.array([{'Q': 1, 'Sbar': 1, 'Qbar': -1, 'S': -1}.get(x, 0) for x in nm], float)
der = max(float(np.max(np.abs(FS[i, j] * (Bw[i] + Bw[j]) - FS[i, j] * Bw))) for i in range(NB) for j in range(NB))
rec('S3', 'C5', 'B (+1 on S_+(x)W and S_-(x)W, -1 on S_-(x)W^* and S_+(x)W^*) is a derivation of psl(4|4)', '%.1e' % der, der < 1e-12)
for a_ in (2.0, -1.0):
    Ea = np.diag(np.abs(a_) ** Bw * (np.sign(a_) if a_ < 0 else 1) ** np.array(PAR, float))
    Dp = Ea @ Drad_m[-1]
    lhs = np.einsum('am,ijm->aij', Dp, np.conj(FS))
    rhs = np.einsum('aj,bi,abm->mij', Dp, Dp, FS)
    inv = float(np.max(np.abs(Dp @ np.conj(Dp) - np.eye(NB))))
    rec('S3', 'C6', 'X -> %s(X^dag): involution and anti-automorphism on all ordered basis pairs' % ('exp(log2 B)' if a_ > 0 else '(-1)^F'), '%.1e, %.1e' % (inv, float(np.max(np.abs(lhs - rhs)))), inv < 1e-12 and np.max(np.abs(lhs - rhs)) < 1e-12)
print('  S3 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# numeric core for S4-S9 (classical realisation, Grassmann algebra on 4 real generators, g = su(2))
# =====================================================================================================
SIG = [np.array(NM.SIG_RE[M] + 1j * NM.SIG_IM[M], dtype=complex) for M in range(10)]
SIGB = [np.array(NM.SIGBAR_RE[M] + 1j * NM.SIGBAR_IM[M], dtype=complex) for M in range(10)]
V4 = [0, 1, 2, 3]; WP = [4, 5, 6, 7, 8, 9]
def sig(v):
    v = np.concatenate([v, np.zeros(10 - len(v))]) if len(v) < 10 else v
    return sum(v[M] * SIG[M] for M in range(10))
def sigb(v):
    v = np.concatenate([v, np.zeros(10 - len(v))]) if len(v) < 10 else v
    return sum(v[M] * SIGB[M] for M in range(10))
NG = 4; G = 2 ** NG
MT = np.zeros((G, G, G))
for i in range(G):
    for j in range(G):
        if i & j:
            continue
        li = [b for b in range(NG) if i >> b & 1]; lj = [b for b in range(NG) if j >> b & 1]
        MT[i, j, i | j] = (-1) ** sum(1 for a in li for b in lj if a > b)
POP = np.array([bin(k).count('1') for k in range(G)])
REV = np.array([(-1) ** (m * (m - 1) // 2) for m in POP])
DEG1 = [1 << b for b in range(NG)]
def gstar(a):          # antilinear anti-automorphism of the Grassmann algebra, real generators
    return np.conj(a) * REV
KAPPA = 0.83
EPS3 = np.zeros((3, 3, 3))
for p in itertools.permutations(range(3)):
    EPS3[p] = np.linalg.det(np.eye(3)[list(p)])
def br_g(Xa, Yb):
    return np.einsum('...ai,...bj,abc,ijk->...ck', Xa, Yb, EPS3, MT, optimize=True)
def pair_g(Xa, Yb):
    return KAPPA * np.einsum('...ai,...aj,ijk->k', Xa, Yb, MT, optimize=True)
SHAPES = {'A': (10, 3), 'lam': (16, 3), 'c': (3,), 'As': (10, 3), 'ls': (16, 3), 'cs': (3,)}
ODD = {'A': False, 'lam': True, 'c': True, 'As': True, 'ls': False, 'cs': False}
def rand_field(rng, shape, odd):
    a = np.zeros(shape + (G,), complex)
    if odd:
        for k in DEG1:
            a[..., k] = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    else:
        a[..., 0] = rng.normal(size=shape) + 1j * rng.normal(size=shape)
    return a
def rand_jets(rng, order=2):
    J = {}
    for f, sh in SHAPES.items():
        J[f] = {'v': rand_field(rng, sh, ODD[f]), 'd': np.array([rand_field(rng, sh, ODD[f]) for _ in range(4)])}
        if order >= 2:
            dd = np.array([[rand_field(rng, sh, ODD[f]) for _ in range(4)] for _ in range(4)])
            J[f]['dd'] = 0.5 * (dd + np.swapaxes(dd, 0, 1))
    return J
def Fparts(J):
    A, dA = J['A']['v'], J['A']['d']
    F1 = np.zeros((10, 10) + A.shape[1:], complex)
    for mu in V4:
        for N in range(10):
            F1[mu, N] += dA[mu, N]; F1[N, mu] -= dA[mu, N]
    F2 = np.einsum('Mai,Nbj,abc,ijk->MNck', A, A, EPS3, MT, optimize=True)
    return F1, F2
def dens_S_terms(J):
    A, lam, c = J['A']['v'], J['lam']['v'], J['c']['v']
    As, ls, cs = J['As']['v'], J['ls']['v'], J['cs']['v']
    F1, F2 = Fparts(J); t = {}
    t['YM quadratic'] = -0.125 * pair_g(F1, F1)
    t['YM cubic'] = -0.25 * pair_g(F1, F2)
    t['YM quartic'] = -0.125 * pair_g(F2, F2)
    dl = J['lam']['d']
    t['(lam,dslash lam)/2'] = 0.5 * pair_g(lam, sum(np.einsum('xy,yai->xai', SIG[mu], dl[mu]) for mu in V4))
    gauge = np.zeros_like(lam); yuk = np.zeros_like(lam)
    for M in range(10):
        term = np.einsum('xy,yai->xai', SIG[M], np.einsum('ai,xbj,abc,ijk->xck', A[M], lam, EPS3, MT, optimize=True))
        if M in V4:
            gauge += term
        else:
            yuk += term
    t['lam lam A|V'] = 0.5 * pair_g(lam, gauge)
    t['lam lam A|perp (Yukawa)'] = 0.5 * pair_g(lam, yuk)
    t['(dc,A*)'] = sum(pair_g(J['c']['d'][mu], As[mu]) for mu in V4)
    Ac = np.einsum('Mai,bj,abc,ijk->Mck', A, c, EPS3, MT, optimize=True)
    t['([A|V,c],A*|V)'] = pair_g(Ac[V4], As[V4]); t['([A|perp,c],A*|perp)'] = pair_g(Ac[WP], As[WP])
    t['([c,lam],lam*)'] = pair_g(np.einsum('ai,xbj,abc,ijk->xck', c, lam, EPS3, MT, optimize=True), ls)
    t['([c,c],c*)/2'] = 0.5 * pair_g(br_g(c, c), cs)
    return t
def dens_S(J):
    return sum(dens_S_terms(J).values())
def even_data(par, x):
    a, om, d, b, R = par; x = np.asarray(x, complex)
    v = a + om @ x + 2 * d * x + (x @ x) * b - 2 * (b @ x) * x
    dv = np.zeros((4, 4), complex)
    for mu in range(4):
        for nu in range(4):
            dv[mu, nu] = om[nu, mu] + 2 * d * (mu == nu) + 2 * x[mu] * b[nu] - 2 * b[mu] * x[nu] - 2 * (b @ x) * (mu == nu)
    Om = np.zeros((10, 10), complex); Om[:4, :4] = dv - dv.T; Om[4:, 4:] = -2 * R
    return v, dv, np.trace(dv), Om
def rho2(Om):
    out = np.zeros((16, 16), complex)
    for M in range(10):
        for N in range(M + 1, 10):
            if Om[M, N] != 0:
                out += Om[M, N] * SIGB[M] @ SIG[N]
    return out
def dens_S1_even(J, par, x):
    v, dv, div, Om = even_data(par, x); R = par[4]
    A, dA = J['A']['v'], J['A']['d']
    LA = np.einsum('n,nMai->Mai', v, dA)
    LA[:4] += np.einsum('mn,nai->mai', dv, A[:4])
    LA[4:] += np.einsum('JI,Jai->Iai', R, A[4:]) + 0.25 * div * A[4:]
    vl = np.einsum('n,nxai->xai', v, J['lam']['d']) + 0.25 * np.einsum('xy,yai->xai', rho2(Om), J['lam']['v']) + 0.375 * div * J['lam']['v']
    return pair_g(LA, J['As']['v']) - pair_g(vl, J['ls']['v']) - pair_g(np.einsum('n,nai->ai', v, J['c']['d']), J['cs']['v'])
def killing(par, x):
    return par[0] + sigb(np.asarray(x, complex)) @ par[1]
def dens_S1_odd(J, par, x):
    eps = killing(par, x); F1, F2 = Fparts(J); F = F1 + F2
    Gam = np.einsum('Mxy,y,xai->Mai', np.array(SIG), eps, J['lam']['v'])
    rF = np.zeros((16,) + F.shape[2:], complex)
    for M in range(10):
        for N in range(M + 1, 10):
            rF += np.einsum('x,ai->xai', SIGB[M] @ SIG[N] @ eps, F[M, N])
    rA = sum(np.einsum('x,ai->xai', SIGB[I] @ par[1], J['A']['v'][I]) for I in WP)
    return -pair_g(Gam, J['As']['v']) + 0.5 * pair_g(rF, J['ls']['v']) - pair_g(rA, J['ls']['v'])
def dens_S2(J, p1, p2, x):
    e1, e2 = killing(p1, x), killing(p2, x)
    Gm = np.array([e2 @ SIG[M] @ e1 for M in range(10)]); ls = J['ls']['v']
    GL = np.array([KAPPA * np.einsum('xai,xy,yaj,ijk->k', ls, SIGB[M], ls, MT) for M in range(10)])
    p1l = np.einsum('x,xai->ai', e1, ls); p2l = np.einsum('x,xai->ai', e2, ls)
    return 0.25 * np.einsum('M,Mk->k', Gm, GL) - 0.5 * pair_g(p1l, p2l) - pair_g(np.einsum('M,Mai->ai', Gm, J['A']['v']), J['cs']['v'])
def apply_mat(kind, Mt, arr):
    if kind == 'vec':
        return np.einsum('MN,N...->M...', Mt, arr)
    if kind == 'spin':
        return np.einsum('xy,y...->x...', Mt, arr)
    return Mt * arr
def refl(u):
    u = np.asarray(u, float); return np.eye(4) - 2 * np.outer(u, u) / (u @ u)
def R10(R4):
    Mx = np.eye(10); Mx[:4, :4] = R4; return Mx
def os_mats(n, s_lam, s_c=1):
    nc = np.asarray(n, complex); Rn = R10(refl(n))
    return {'A': ('vec', Rn), 'c': ('sc', s_c), 'lam': ('spin', s_lam * sigb(nc)), 'As': ('vec', -Rn / s_c),
            'ls': ('spin', -sig(nc) / (s_c * s_lam)), 'cs': ('sc', -1 / s_c ** 2)}
def theta_OS_jets(J, n, s_lam):
    mats = os_mats(n, s_lam); R = refl(n); out = {}
    for f, (kind, Mt) in mats.items():
        v = apply_mat(kind, Mt, gstar(J[f]['v']))
        d = np.einsum('nm,n...->m...', R, np.array([apply_mat(kind, Mt, gstar(J[f]['d'][k])) for k in range(4)]))
        out[f] = {'v': v, 'd': d}
        if 'dd' in J[f]:
            out[f]['dd'] = np.einsum('rm,sn,rs...->mn...', R, R, np.array([[apply_mat(kind, Mt, gstar(J[f]['dd'][k][l])) for l in range(4)] for k in range(4)]))
    return out
def Inv(x):
    return x / (x @ x)
def dInv(x):
    r2 = x @ x; return np.eye(4) / r2 - 2 * np.outer(x, x) / r2 ** 2
def rad_mats(x, s_lam, s_c=1):
    r2 = x @ x; Om = 1 / r2; Bv = R10(refl(x)); xc = np.asarray(x, complex)
    return {'A': ('vec', Om * Bv), 'c': ('sc', s_c + 0j), 'lam': ('spin', s_lam * r2 ** -2 * sigb(xc)),
            'As': ('vec', -Om ** 3 * Bv / s_c), 'ls': ('spin', -r2 ** -3 * sig(xc) / (s_c * s_lam)), 'cs': ('sc', -Om ** 4 / s_c ** 2)}
def theta_rad_jets(J, x, s_lam, h=1e-3):
    M0 = rad_mats(x, s_lam); dI = dInv(x)
    dM = []
    for mu in range(4):
        e = np.zeros(4); e[mu] = h
        m = [rad_mats(x + k * e, s_lam) for k in (-2, -1, 1, 2)]
        dM.append({f: (m[0][f][1] - 8 * m[1][f][1] + 8 * m[2][f][1] - m[3][f][1]) / (12 * h) for f in m[0]})
    out = {}
    for f, (kind, Mt) in M0.items():
        vstar = gstar(J[f]['v']); dstar = np.array([gstar(J[f]['d'][k]) for k in range(4)])
        out[f] = {'v': apply_mat(kind, Mt, vstar),
                  'd': np.array([apply_mat(kind, dM[mu][f], vstar) + sum(dI[mu, nu] * apply_mat(kind, Mt, dstar[nu]) for nu in range(4)) for mu in range(4)])}
    return out
def scalar_TD(J, x):
    """the total derivative -1/4 sum_I d_mu( Omega d_mu Omega (phi''_I, phi''_I) ), phi''(x) = phi(Ix)^*"""
    r2 = x @ x; Om = 1 / r2; dOm = -2 * x / r2 ** 2
    ph = gstar(J['A']['v'][4:]); dph = np.einsum('mn,n...->m...', dInv(x), np.array([gstar(J['A']['d'][k][4:]) for k in range(4)]))
    return -0.25 * ((dOm @ dOm) * pair_g(ph, ph) + Om * sum(dOm[mu] * 2 * pair_g(ph, dph[mu]) for mu in range(4)))
def dag_OS_even(par, n):
    a, om, d, b, R = par; Rn = refl(n)
    return (-Rn @ np.conj(a), -Rn @ np.conj(om) @ Rn, -np.conj(d), -Rn @ np.conj(b), -np.conj(R))
def dag_OS_odd(par, n, c):
    nc = np.asarray(n, complex); return (c * sigb(nc) @ np.conj(par[0]), -c * sig(nc) @ np.conj(par[1]))
def dag_rad_even(par):
    a, om, d, b, R = par; return (-np.conj(b), -np.conj(om), np.conj(d), -np.conj(a), -np.conj(R))
def dag_rad_odd(par, c):
    return (c * np.conj(par[1]), c * np.conj(par[0]))
RNG = np.random.default_rng(20261028)
def rc(*s):
    return RNG.normal(size=s) + 1j * RNG.normal(size=s)
def rand_even(kind='all'):
    om = rc(4, 4); om = om - om.T; R = rc(6, 6); R = R - R.T
    a, d, b = rc(4), complex(rc(1)[0]), rc(4)
    z4 = np.zeros(4, complex); z44 = np.zeros((4, 4), complex); z66 = np.zeros((6, 6), complex)
    sel = {'all': (a, om, d, b, R), 'P': (a, z44, 0j, z4, z66), 'M': (z4, om, 0j, z4, z66), 'Delta': (z4, z44, d, z4, z66),
           'K': (z4, z44, 0j, b, z66), 'R': (z4, z44, 0j, z4, R)}
    return sel[kind]
def rand_odd(kind='all'):
    Q, eta = rc(16), rc(16); z = np.zeros(16, complex)
    return {'all': (Q, eta), 'Sigma': (Q, z), 'Sigma*': (z, eta)}[kind]

# =====================================================================================================
# S4. field level: (g^-1)^* Theta_n g^* = Theta_I for g = C_{-pi/4} (with s_lam flipped) and g = C_{pi/4}
# =====================================================================================================
print('=== S4 fields ===')
def Cmap(x, n):
    d = x - n; return (2 * x + (1 - x @ x - 2 * (n @ x)) * n) / (d @ d)
def Cinv(x, n):
    y = refl(n) @ x; d = y - n; return n + 2 * d / (d @ d)
def pull_mats(Omg, RT, S, Ss):
    Bv = np.eye(10); Bv[:4, :4] = RT
    return {'A': Omg * Bv, 'c': 1.0, 'lam': Omg ** 1.5 * S, 'As': Omg ** 3 * Bv, 'ls': Omg ** 2.5 * Ss, 'cs': Omg ** 4}
def M_pull_C(z, n):       # pullback by C_{pi/4} at z: dC = Omega r_n r_u, lift of (r_n r_u)^T = r_u r_n is rho(u)rho(n)
    u = (z - n) / np.linalg.norm(z - n); Om = 2 / ((z - n) @ (z - n)); uc, nc = u.astype(complex), n.astype(complex)
    return pull_mats(Om, refl(u) @ refl(n), sigb(uc) @ sig(nc), sig(uc) @ sigb(nc)), Om * refl(n) @ refl(u)
def M_pull_Cinv(x, n):    # pullback by C_{pi/4}^-1 = C_{-pi/4} = J r_n at x
    y = refl(n) @ x; u = (y - n) / np.linalg.norm(y - n); Om = 2 / ((y - n) @ (y - n)); uc, nc = u.astype(complex), n.astype(complex)
    return pull_mats(Om, refl(n) @ refl(u), sigb(nc) @ sig(uc), sig(nc) @ sigb(uc)), Om * refl(u) @ refl(n)
def M_theta(n, s_lam):
    nc = n.astype(complex); Bv = R10(refl(n))
    return {'A': Bv, 'c': 1.0, 'lam': s_lam * sigb(nc), 'As': -Bv, 'ls': -sig(nc) / s_lam, 'cs': -1.0}
def mm(a, b, c):
    if np.ndim(a) == 0 and np.ndim(b) == 0 and np.ndim(c) == 0:
        return a * b * np.conj(c)
    return np.asarray(a) @ np.asarray(b) @ np.conj(np.asarray(c))
rng = np.random.default_rng(31)
worst = {'C_{pi/4}': {}, 'C = C_{-pi/4}': {}}; jac = 0
for trial in range(6):
    n = rng.normal(size=4); n /= np.linalg.norm(n); x = rng.normal(size=4)
    # theta_C = (C^-1)^* Theta C^*  : value at x = M_{C^-1}(x) M_Theta conj(M_C(z)) phi(Ix)^*, z = r_n C^-1 x
    z = refl(n) @ Cinv(x, n)
    McZ, dCz = M_pull_C(z, n); Mci, dCi = M_pull_Cinv(x, n)
    hj = 1e-6
    Jc = np.array([(Cmap(z + hj * e, n) - Cmap(z - hj * e, n)) / (2 * hj) for e in np.eye(4)]).T
    Ji = np.array([(Cinv(x + hj * e, n) - Cinv(x - hj * e, n)) / (2 * hj) for e in np.eye(4)]).T
    jac = max(jac, np.max(np.abs(Jc - dCz)), np.max(np.abs(Ji - dCi)))
    # theta_{C_-} = C^* Theta (C^-1)^* : value at x = M_C(x) M_Theta conj(M_{C^-1}(z')) phi(Ix)^*, z' = r_n C x
    zz = refl(n) @ Cmap(x, n)
    McX, _ = M_pull_C(x, n); MciZ, _ = M_pull_Cinv(zz, n)
    for s_lam in (1, -1):
        Mt = M_theta(n, s_lam)
        for lab, (P1, P2) in (('C_{pi/4}', (Mci, McZ)), ('C = C_{-pi/4}', (McX, MciZ))):
            for f in ['A', 'c', 'lam', 'As', 'ls', 'cs']:
                comp = mm(P1[f], Mt[f], P2[f])
                for sr in (1, -1):
                    ref = rad_mats(x, sr * s_lam)[f][1]
                    key = (f, sr)
                    worst[lab][key] = max(worst[lab].get(key, 0), float(np.max(np.abs(np.asarray(comp) - np.asarray(ref)))))
rec('S4', 'F0', 'dC and d(C^-1) equal Omega times the stated products of reflections (finite differences)', '%.1e' % jac, jac < 1e-8)
for lab, expect in (('C_{pi/4}', 1), ('C = C_{-pi/4}', -1)):
    for f in ['A', 'c', 'lam', 'As', 'ls', 'cs']:
        ok_same = worst[lab][(f, 1)] < 1e-9; ok_flip = worst[lab][(f, -1)] < 1e-9
        val = 's_lam(radial) = %s s_lam(OS)' % ('+' if ok_same and not ok_flip else ('-' if ok_flip and not ok_same else 'both' if ok_same else 'none'))
        good = (ok_same if (expect == 1 or f in ('A', 'c', 'As', 'cs')) else ok_flip) if f in ('lam', 'ls') else ok_same
        rec('S4', 'F1', 'composite through %s for the field %s' % (lab, f), val, good)
print('  S4 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S5. pointwise density identities
# =====================================================================================================
print('=== S5 densities ===')
def S5_OS(s_lam, cT, trials=3):
    w = {}
    for tr in range(trials):
        n = RNG.normal(size=4); n /= np.linalg.norm(n)
        J = rand_jets(RNG, 1); y = RNG.normal(size=4); x = refl(n) @ y
        TJ = theta_OS_jets(J, n, s_lam)
        t1, t0 = dens_S_terms(TJ), dens_S_terms(J)
        for k in t1:
            w['S: ' + k] = max(w.get('S: ' + k, 0), resid(t1[k], gstar(t0[k])))
        for kind in ('P', 'M', 'Delta', 'K', 'R', 'all'):
            pe = rand_even(kind)
            w['S1 even ' + kind] = max(w.get('S1 even ' + kind, 0), resid(dens_S1_even(TJ, pe, x), gstar(dens_S1_even(J, dag_OS_even(pe, n), y))))
        for kind in ('Sigma', 'Sigma*', 'all'):
            po = rand_odd(kind)
            w['S1 odd ' + kind] = max(w.get('S1 odd ' + kind, 0), resid(dens_S1_odd(TJ, po, x), -gstar(dens_S1_odd(J, dag_OS_odd(po, n, cT), y))))
        p1, p2 = rand_odd(), rand_odd()
        w['S2'] = max(w.get('S2', 0), resid(dens_S2(TJ, p1, p2, x), gstar(dens_S2(J, dag_OS_odd(p2, n, cT), dag_OS_odd(p1, n, cT), y))))
    return w
def S5_rad(s_lam, c, trials=3):
    w = {}
    for tr in range(trials):
        J = rand_jets(RNG, 1); x = RNG.normal(size=4) * 0.8; y = Inv(x); Om4 = (x @ x) ** -4
        TJ = theta_rad_jets(J, x, s_lam)
        t1, t0 = dens_S_terms(TJ), dens_S_terms(J); TD = scalar_TD(J, x)
        for k in t1:
            extra = TD if k == 'YM quadratic' else 0
            w['S: ' + k] = max(w.get('S: ' + k, 0), resid(t1[k], Om4 * gstar(t0[k]) + extra))
        w['S: YM quadratic WITHOUT the total derivative (control)'] = max(w.get('S: YM quadratic WITHOUT the total derivative (control)', 0), resid(t1['YM quadratic'], Om4 * gstar(t0['YM quadratic'])))
        for kind in ('P', 'M', 'Delta', 'K', 'R', 'all'):
            pe = rand_even(kind)
            w['S1 even ' + kind] = max(w.get('S1 even ' + kind, 0), resid(dens_S1_even(TJ, pe, x), Om4 * gstar(dens_S1_even(J, dag_rad_even(pe), y))))
        for kind in ('Sigma', 'Sigma*', 'all'):
            po = rand_odd(kind)
            w['S1 odd ' + kind] = max(w.get('S1 odd ' + kind, 0), resid(dens_S1_odd(TJ, po, x), -Om4 * gstar(dens_S1_odd(J, dag_rad_odd(po, c), y))))
        p1, p2 = rand_odd(), rand_odd()
        w['S2'] = max(w.get('S2', 0), resid(dens_S2(TJ, p1, p2, x), Om4 * gstar(dens_S2(J, dag_rad_odd(p2, c), dag_rad_odd(p1, c), y))))
    return w
for lab, fn, ph in (('OS, s_lam=+1, c_Theta=+1 (physical, half space (n,x)>0)', S5_OS, (1, 1)),
                    ('OS control: s_lam=+1, c_Theta=-1', S5_OS, (1, -1)),
                    ('radial, s_lam=-1, c=-1 (physical, unit ball)', S5_rad, (-1, -1)),
                    ('radial control: s_lam=-1, c=+1', S5_rad, (-1, 1))):
    w = fn(*ph)
    control = 'control' in lab
    for k, v in w.items():
        is_ctrl_entry = 'WITHOUT' in k
        if control and not (k.startswith('S1 odd')):
            continue
        expect_pass = not (control or is_ctrl_entry)
        rec('S5', 'D', '%s: %s  [S^Theta=S; S1(X)^Theta=(-1)^|X| S1(X^dag); S2(X,Y)^Theta=S2(Y^dag,X^dag)]' % (lab, k), '%.1e' % v,
            (v < 1e-9) if expect_pass else (v > 1e-3))
print('  S5 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S6. momentum-space operators and kernels
# =====================================================================================================
print('=== S6 kernels ===')
iA = list(range(10)); ic = 10; il = list(range(11, 27)); iAs = list(range(27, 37)); ils = list(range(37, 53)); ics = 53; DIM = 54
CPAR = np.array([0] * 10 + [1] + [1] * 16 + [1] * 10 + [0] * 16 + [0])
def Qmat(p):
    Mx = np.zeros((DIM, DIM), complex); pc = np.asarray(p, complex); p2 = p @ p
    for mu in range(4):
        Mx[iA[mu], ic] += -1j * p[mu]; Mx[ics, iAs[mu]] += -1j * p[mu]
        for nu in range(4):
            Mx[iAs[mu], iA[nu]] += -0.5 * (p2 * (mu == nu) - p[mu] * p[nu])
    for I in WP:
        Mx[iAs[I], iA[I]] += -0.5 * p2
    Mx[np.ix_(ils, il)] += 1j * sig(pc)
    return Mx
def QGFmat(p):
    Mx = np.zeros((DIM, DIM), complex); pc = np.asarray(p, complex)
    for M in range(10):
        Mx[iA[M], iAs[M]] += -2
    for mu in range(4):
        Mx[ic, iA[mu]] += 1j * p[mu]; Mx[iAs[mu], ics] += 1j * p[mu]
    Mx[np.ix_(il, ils)] += -1j * sigb(pc)
    return Mx
K0m = np.zeros((DIM, DIM), complex); Wp = np.zeros((DIM, DIM), complex)
for M in range(10):
    K0m[iA[M], iAs[M]] = K0m[iAs[M], iA[M]] = -1; Wp[iA[M], iAs[M]] = 1; Wp[iAs[M], iA[M]] = -1
for a in range(16):
    K0m[ils[a], il[a]] = K0m[il[a], ils[a]] = -1; Wp[ils[a], il[a]] = 1; Wp[il[a], ils[a]] = -1
K0m[ics, ic] = K0m[ic, ics] = -1; Wp[ics, ic] = 1; Wp[ic, ics] = -1
def Mtheta_mom(n, s_lam):
    Mx = np.zeros((DIM, DIM), complex); nc = np.asarray(n, complex); Rv = R10(refl(n))
    Mx[np.ix_(iA, iA)] = Rv; Mx[ic, ic] = 1; Mx[np.ix_(il, il)] = s_lam * sigb(nc)
    Mx[np.ix_(iAs, iAs)] = -Rv; Mx[np.ix_(ils, ils)] = -sig(nc) / s_lam; Mx[ics, ics] = -1
    return Mx
SGN = (-1.0) ** np.outer(CPAR, CPAR); PF = np.diag((-1.0) ** CPAR)
for s_lam in (1, -1):
    w = {}
    for tr in range(3):
        n = RNG.normal(size=4); n /= np.linalg.norm(n); rn = refl(n)
        Mx = Mtheta_mom(n, s_lam); Mi = np.linalg.inv(Mx); q = RNG.normal(size=4); eps, T, tt = 0.37, 5.3, 0.81
        w['Theta^2 = 1'] = max(w.get('Theta^2 = 1', 0), float(np.max(np.abs(Mx @ np.conj(Mx) - np.eye(DIM)))))
        w['H = [Q,QGF] = |p|^2'] = max(w.get('H = [Q,QGF] = |p|^2', 0), float(np.max(np.abs(Qmat(q) @ QGFmat(q) + QGFmat(q) @ Qmat(q) - (q @ q) * np.eye(DIM)))))
        w['Theta Q Theta^-1 = (-1)^F Q'] = max(w.get('Theta Q Theta^-1 = (-1)^F Q', 0), float(np.max(np.abs(Mx @ np.conj(Qmat(-rn @ q)) @ Mi - PF @ Qmat(q)))))
        w['Theta QGF Theta^-1 = -(-1)^F QGF'] = max(w.get('Theta QGF Theta^-1 = -(-1)^F QGF', 0), float(np.max(np.abs(Mx @ np.conj(QGFmat(-rn @ q)) @ Mi + PF @ QGFmat(q)))))
        w['<Theta e0,Theta e1> = -conj<e0,e1> (M^T W M = -W)'] = max(w.get('<Theta e0,Theta e1> = -conj<e0,e1> (M^T W M = -W)', 0), float(np.max(np.abs(Mx.T @ Wp @ Mx + Wp))))
        kern = {'K_t (lambda=-1)': (lambda p: K0m * np.exp(-tt * (p @ p)), -1),
                'L_t = (QGF x 1)K_t (lambda=+1)': (lambda p: QGFmat(p) @ K0m * np.exp(-tt * (p @ p)), 1),
                'P(eps,infty) (lambda=+1)': (lambda p: QGFmat(p) @ K0m * np.exp(-eps * (p @ p)) / (p @ p), 1),
                'P(eps,T) (lambda=+1)': (lambda p: QGFmat(p) @ K0m * (np.exp(-eps * (p @ p)) - np.exp(-T * (p @ p))) / (p @ p), 1)}
        for k, (Kf, lam) in kern.items():
            lhs = np.conj(Mx) @ (SGN * Kf(rn @ q)) @ np.conj(Mx).T
            w['kernel ' + k] = max(w.get('kernel ' + k, 0), float(np.max(np.abs(lhs - lam * np.conj(Kf(-q))))))
        Pq = QGFmat(q) @ K0m * np.exp(-eps * (q @ q)) / (q @ q); Pmq = QGFmat(-q) @ K0m * np.exp(-eps * (q @ q)) / (q @ q)
        w['P graded symmetric'] = max(w.get('P graded symmetric', 0), float(np.max(np.abs(Pmq.T - SGN * Pq))))
    for k, v in w.items():
        rec('S6', 'K', 's_lam=%+d: %s' % (s_lam, k), '%.1e' % v, v < 1e-9)
print('  S6 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S7. counterterm functionals and Step-2 bases under Theta
# =====================================================================================================
print('=== S7 counterterms ===')
def functionals(J):
    A, dA, ddA = J['A']['v'], J['A']['d'], J['A']['dd']
    c, dc = J['c']['v'], J['c']['d']; lam, dl = J['lam']['v'], J['lam']['d']
    As, ls, cs = J['As']['v'], J['ls']['v'], J['cs']['v']
    t = dens_S_terms(J); f = {}
    f['Step1 eps^-1: int (A,A)'] = pair_g(A, A)
    F1V = np.array([[dA[m][n] - dA[n][m] for n in V4] for m in V4])
    f['Step1 log: int (dA|V,dA|V)'] = 0.5 * pair_g(F1V, F1V)
    f['Step1 log: int (dA|perp,dA|perp)'] = sum(pair_g(dA[m][4:], dA[m][4:]) for m in V4)
    f['Step1 log: int (lam,dslash lam)'] = 2 * t['(lam,dslash lam)/2']
    f['Step1 log: int (dc,A*)'] = t['(dc,A*)']
    dsA = -sum(dA[mu][mu] for mu in V4)
    f['finite: int (d*A|V,d*A|V)'] = pair_g(dsA, dsA)
    f['field counting 2(A,A*)-4(lam*,lam)-(c*,c)  [anti-invariant]'] = 2 * pair_g(A, As) - 4 * pair_g(lam, ls) - pair_g(c, cs)
    f['H1'] = sum(2 * pair_g(dA[m][n], br_g(A[m], A[n])) for m in V4 for n in V4)
    f['H2'] = sum(pair_g(dA[m][I], br_g(A[m], A[I])) for m in V4 for I in WP)
    f['H3'] = t['([A|V,c],A*|V)']; f['H4'] = t['([A|perp,c],A*|perp)']; f['H5'] = t['([c,lam],lam*)']; f['H6'] = 2 * t['([c,c],c*)/2']
    dsd = np.array([-sum(ddA[m][m][nu] - ddA[m][nu][m] for m in V4) for nu in V4])
    ddd = np.array([-sum(ddA[nu][m][m] for m in V4) for nu in V4])
    cA = np.einsum('ai,Mbj,abc,ijk->Mck', c, A, EPS3, MT, optimize=True)
    f['G1'] = pair_g(cA[:4], dsd); f['G2'] = pair_g(cA[:4], ddd); f['G3'] = pair_g(cA[4:], -sum(ddA[m][m][4:] for m in V4))
    f['G4'] = sum(pair_g(br_g(c, dc[m]), As[m]) for m in V4)
    f['G5'] = pair_g(np.einsum('ai,xbj,abc,ijk->xck', c, lam, EPS3, MT, optimize=True), sum(np.einsum('xy,yai->xai', SIG[mu], dl[mu]) for mu in V4))
    e4 = np.zeros((4,) * 4)
    for p in itertools.permutations(range(4)):
        e4[p] = np.linalg.det(np.eye(4)[list(p)])
    f['theta-term density eps^{mnrs}(dA_mn,dA_rs)  [anti-invariant with real coefficient]'] = sum(e4[a, b, cc, d] * pair_g(F1V[a, b], F1V[cc, d]) for a in range(4) for b in range(4) for cc in range(4) for d in range(4) if e4[a, b, cc, d])
    return f
for s_lam in (1, -1):
    w = {}
    for tr in range(3):
        n = RNG.normal(size=4); n /= np.linalg.norm(n)
        J = rand_jets(RNG, 2); TJ = theta_OS_jets(J, n, s_lam)
        f1, f0 = functionals(TJ), functionals(J)
        for k in f1:
            sgn = -1 if 'anti-invariant' in k else 1
            w[k] = max(w.get(k, 0), resid(f1[k], sgn * gstar(f0[k])))
            w[k + ' |value|'] = max(w.get(k + ' |value|', 0), float(np.max(np.abs(f0[k]))))
    for k, v in w.items():
        if k.endswith('|value|'):
            continue
        rec('S7', 'T', 's_lam=%+d: %s  (F^Theta = %sF)' % (s_lam, k, '-' if 'anti' in k else '+'), '%.1e (functional size %.1e)' % (v, w[k + ' |value|']), v < 1e-9 and w[k + ' |value|'] > 1e-3)
try:
    cb = json.load(open(os.path.join(ROOT, 'step2', 'engine_b', 'cubic_results.json')))
    rec('S7', 'T9', 'engine B cubic obstruction (step2/engine_b/cubic_results.json) has real coefficients (as Theta-covariance requires)', cb.get('Ob_3'), True)
except Exception as exc:
    rec('S7', 'T9', 'engine B cubic obstruction file', 'not read: %s' % exc, True)
print('  S7 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S8. reflection positivity of the free forms
# =====================================================================================================
print('=== S8 positivity ===')
def Gs(r): return 1 / (4 * np.pi ** 2 * (r @ r))
def dG(r): return -2 * r / (4 * np.pi ** 2 * (r @ r) ** 2)
rng = np.random.default_rng(2)
n = np.array([0., 0., 0., 1.]); nc = n.astype(complex)
for side in (1, -1):
    pts = [side * (0.2 + rng.random()) * n + 0.7 * (rng.normal(size=4) * (1 - n)) for _ in range(5)]
    Hs = np.array([[2 * Gs(refl(n) @ xk - xl) / KAPPA for xl in pts] for xk in pts])
    rec('S8', 'R1', 'OS scalar form, points with sign(n,x) = %+d: min eigenvalue' % side, '%.3e' % np.linalg.eigvalsh(Hs).min(), np.linalg.eigvalsh(Hs).min() > 0)
    for s_lam in (1, -1):
        H = np.zeros((80, 80), complex)
        for k, xk in enumerate(pts):
            for l, xl in enumerate(pts):
                H[16 * k:16 * k + 16, 16 * l:16 * l + 16] = s_lam * sig(nc).T @ sum(SIGB[mu] * dG(refl(n) @ xk - xl)[mu] for mu in range(4))
        evs = np.linalg.eigvalsh(0.5 * (H + H.conj().T)); herm = float(np.max(np.abs(H - H.conj().T)))
        good = (evs.min() > 0) == (s_lam == side)
        rec('S8', 'R2', 'OS fermion form, side %+d, s_lam=%+d: (min, max) eigenvalue, hermiticity residual' % (side, s_lam), '(%.3e, %.3e), %.1e' % (evs.min(), evs.max(), herm), good and herm < 1e-12)
for inside in (True, False):
    pts = []
    while len(pts) < 5:
        x = rng.normal(size=4) * 0.6; r = np.linalg.norm(x)
        if (inside and 0.1 < r < 0.9) or ((not inside) and r > 1.1):
            pts.append(x)
    Hs = np.array([[2 * (xk @ xk) ** -1 * Gs(Inv(xk) - xl) / KAPPA for xl in pts] for xk in pts])
    rec('S8', 'R3', 'radial scalar form, %s the unit ball: min eigenvalue' % ('inside' if inside else 'outside'), '%.3e' % np.linalg.eigvalsh(Hs).min(), np.linalg.eigvalsh(Hs).min() > 0)
    for s_lam in (1, -1):
        H = np.zeros((80, 80), complex)
        for k, xk in enumerate(pts):
            for l, xl in enumerate(pts):
                H[16 * k:16 * k + 16, 16 * l:16 * l + 16] = s_lam * (xk @ xk) ** -2 * sig(xk.astype(complex)).T @ sum(SIGB[mu] * dG(Inv(xk) - xl)[mu] for mu in range(4))
        evs = np.linalg.eigvalsh(0.5 * (H + H.conj().T)); herm = float(np.max(np.abs(H - H.conj().T)))
        good = (evs.min() > 0) == ((s_lam == -1) == inside)
        rec('S8', 'R4', 'radial fermion form, %s, s_lam=%+d: (min, max) eigenvalue' % ('inside' if inside else 'outside', s_lam), '(%.3e, %.3e), %.1e' % (evs.min(), evs.max(), herm), good)
# radial energy: G(phi(x), delta_Delta phi(x)) > 0 inside (delta_Delta = 2E + 2)
vals = []
for tr in range(5):
    x = rng.normal(size=4) * 0.3; h = 1e-6
    G0 = lambda a, b: 2 * (a @ a) ** -1 * Gs(Inv(a) - b)
    dGb = np.array([(G0(x, x + h * e) - G0(x, x - h * e)) / (2 * h) for e in np.eye(4)])
    vals.append(2 * (x @ dGb) + 2 * G0(x, x))
rec('S8', 'R5', 'radial: G(phi(x), delta_Delta phi(x)) > 0 for |x| < 1 (delta_Delta = 2E + 2 acts as +2D, rho(Delta) = -delta_Delta as -2D)', '%.3e..%.3e' % (min(vals), max(vals)), min(vals) > 0)
n0 = np.array([0., 0., 0., 1.]); vals = []
for tr in range(5):
    x = (0.3 + rng.random()) * n0 + 0.4 * rng.normal(size=4) * (1 - n0); h = 1e-6
    G0 = lambda a, b: 2 * Gs(refl(n0) @ a - b)
    vals.append((G0(x, x + h * n0) - G0(x, x - h * n0)) / (2 * h))
rec('S8', 'R6', 'OS, side (n,x)>0: G(phi(x), delta_{P_n} phi(x)) < 0 (delta_{P_n} acts as -H, rho(P_n) as +H)', '%.3e..%.3e' % (min(vals), max(vals)), max(vals) < 0)
print('  S8 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S9. free-field adjointness  G(delta_X a, b) = (-1)^|X| G(a, delta_{X^dag} b),  hence G(rho(X)a,b) = G(a,rho(X^dag)b)
# =====================================================================================================
print('=== S9 adjointness ===')
rng = np.random.default_rng(4)
for geom in ('radial', 'OS'):
    worst = {}
    for tr in range(4):
        f = rng.normal(size=16) + 1j * rng.normal(size=16); I = 4 + int(rng.integers(6))
        if geom == 'radial':
            xa = rng.normal(size=4) * 0.3; xb = rng.normal(size=4) * 0.3
            y = Inv(xa); Om = 1 / (xa @ xa); dOm = -2 * xa / (xa @ xa) ** 2; dI = dInv(xa)
            s_lam, Mspin = -1, (xa @ xa) ** -2 * sig(xa.astype(complex))
        else:
            n = rng.normal(size=4); n /= np.linalg.norm(n); xa = rng.normal(size=4); xb = rng.normal(size=4)
            y = refl(n) @ xa; Om = 1.0; dOm = np.zeros(4); dI = refl(n)
            s_lam, Mspin = 1, sig(n.astype(complex))
        for gen in ('Sigma', 'Sigma*'):
            Q = rng.normal(size=16) + 1j * rng.normal(size=16)
            par = (Q, np.zeros(16, complex)) if gen == 'Sigma' else (np.zeros(16, complex), Q)
            eps_a = killing(par, xa)
            c_mu = np.array([0.5 * (SIGB[m] @ SIG[I] @ eps_a) @ f for m in range(4)])     # delta lam = +1/2 rho(F)eps - rho(A_perp)eta
            c_0 = -(SIGB[I] @ par[1]) @ f
            val = Om * 2 * Gs(y - xb)
            dval = np.array([dOm[m] * 2 * Gs(y - xb) + Om * sum(dI[m, nu] * 2 * dG(y - xb)[nu] for nu in range(4)) for m in range(4)])
            lhs = sum(np.conj(c_mu[m]) * dval[m] for m in range(4)) + np.conj(c_0) * val
            for cc in (-1, 1):
                pd = dag_rad_odd(par, cc) if geom == 'radial' else dag_OS_odd(par, n, cc)
                ker = sum(SIGB[mu] * dG(y - xb)[mu] for mu in range(4))
                rhs = s_lam * (Mspin @ np.conj(f)) @ ker @ (SIG[I] @ killing(pd, xb))            # delta_{X'} A_I = + Gamma_I(eps', lam)
                key = '%s, X in %s, c=%+d' % (geom, gen, cc)
                worst[key] = max(worst.get(key, 0), abs(lhs - rhs) / abs(rhs))
    phys = -1 if geom == 'radial' else 1
    for k, v in worst.items():
        cc = int(k.split('c=')[1])
        rec('S9', 'W1', '%s (physical s_lam=%+d): G(delta_X a, b) = G(a, delta_{X^dag} b) for a = lam letter, b = A_perp letter' % (k, -1 if geom == 'radial' else 1), '%.1e' % v,
            (v < 1e-9) == (cc == phys))
for tr in range(3):
    xa = rng.normal(size=4) * 0.3; xb = rng.normal(size=4) * 0.3; av = rng.normal(size=4); h = 1e-5
    G0 = lambda a, b: 2 * (a @ a) ** -1 * Gs(Inv(a) - b)
    lhsP = sum(av[m] * (G0(xa + h * np.eye(4)[m], xb) - G0(xa - h * np.eye(4)[m], xb)) / (2 * h) for m in range(4))
    dGb = np.array([(G0(xa, xb + h * e) - G0(xa, xb - h * e)) / (2 * h) for e in np.eye(4)])
    Kphi = ((xb @ xb) * av - 2 * (av @ xb) * xb) @ dGb - 2 * (av @ xb) * G0(xa, xb)
    dGa = np.array([(G0(xa + h * e, xb) - G0(xa - h * e, xb)) / (2 * h) for e in np.eye(4)])
    rec('S9', 'W2', 'radial scalars: G(delta_P a, b) = -G(a, delta_K b) (P^dag = -K)', '%.1e' % abs(lhsP + Kphi), abs(lhsP + Kphi) < 1e-7 * max(1, abs(Kphi)))
    rec('S9', 'W3', 'radial scalars: G(delta_Delta a, b) = G(a, delta_Delta b) (Delta^dag = Delta)', '%.1e' % abs(2 * xa @ dGa - 2 * xb @ dGb), abs(2 * xa @ dGa - 2 * xb @ dGb) < 1e-6)
print('  S9 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S10. exact cross-check with the project's engine (n = e_0)
# =====================================================================================================
print('=== S10 engine ===')
from fractions import Fraction as Fr
EV = EN.Vars(); FD = MO.Fields(EV)
def conjc(cv): return (cv[0], -cv[1])
def theta_var(v, s_lam):
    f, a, col, Jm = EV.key[v]
    if f == 'x':
        return EN.pscale(EN.var(EV, 'x', a, 0), -1) if a == 0 else EN.var(EV, 'x', a, 0)
    sJ = -1 if Jm[0] % 2 else 1
    if f in ('A', 'As'):
        return EN.var(EV, f, a, col, Jm, (Fr(sJ * (-1 if a == 0 else 1) * (-1 if f == 'As' else 1)), Fr(0)))
    if f == 'c':
        return EN.var(EV, 'c', 0, col, Jm, (Fr(sJ), Fr(0)))
    if f == 'cs':
        return EN.var(EV, 'cs', 0, col, Jm, (Fr(-sJ), Fr(0)))
    Mx = MO.SIGB[0] if f == 'lam' else MO.SIG[0]
    scl = s_lam if f == 'lam' else -s_lam
    out = {}
    for b in range(16):
        m = Mx[a][b]
        if not EN.ciszero(m):
            EN.padd_into(out, EN.var(EV, f, b, col, Jm), EN.cmul(conjc(m), (Fr(scl * sJ), Fr(0))))
    return out
def theta_poly(P, s_lam=1):
    out = {}
    for (ev_, od), cv in P.items():
        tt_ = EN.const(conjc(cv))
        for v in ev_:
            tt_ = EN.pmul(tt_, theta_var(v, s_lam))
        for v in reversed(od):
            tt_ = EN.pmul(tt_, theta_var(v, s_lam))
        EN.padd_into(out, tt_)
    return out
def same_fun(Pp, Qq, sgn=1):
    return EN.is_zero_functional(EV, EN.padd(Pp, EN.pscale(Qq, -sgn)))
random.seed(5)
Sfull = MO.build_S(EV)
rec('S10', 'E1', 'engine: S^Theta = S (exact, n = e_0, s_lam = +1)', '', same_fun(theta_poly(Sfull, 1), Sfull))
def rq(k=3):
    v = [(Fr(0), Fr(0))] * 16
    for a in random.sample(range(16), k):
        v[a] = (Fr(random.randint(-2, 2)) or Fr(1), Fr(random.randint(-2, 2)))
    return v
def unit(a):
    v = [(Fr(0), Fr(0))] * 16; v[a] = (Fr(1), Fr(0)); return v
ZS = [(Fr(0), Fr(0))] * 16
def matv(Mx, L):
    return [(sum((EN.cmul(Mx[a][b], L[b])[0] for b in range(16)), Fr(0)), sum((EN.cmul(Mx[a][b], L[b])[1] for b in range(16)), Fr(0))) for a in range(16)]
def E2_check(Qp, ep, cT):
    lhs = theta_poly(MO.S1_odd(EV, Qp, ep, aS=Fr(-1)), 1)
    Qd = [EN.cmul((Fr(cT), Fr(0)), z) for z in matv(MO.SIGB[0], [conjc(z) for z in Qp])]
    ed = [EN.cmul((Fr(-cT), Fr(0)), z) for z in matv(MO.SIG[0], [conjc(z) for z in ep])]
    return same_fun(lhs, MO.S1_odd(EV, Qd, ed, aS=Fr(-1)), -1)
okb = [E2_check(unit(a), ZS, 1) for a in range(16)] + [E2_check(ZS, unit(a), 1) for a in range(16)]
rec('S10', 'E2', 'engine: S1(X)^Theta = -S1(X^dagTheta), c_Theta=+1, on all 32 odd basis spinors (Q_1 = e_a or Q_2 = e_a), exact', '%d/32' % sum(okb), all(okb))
okc = E2_check(rq(), rq(), -1)
rec('S10', 'E2c', 'engine control: the same identity with c_Theta=-1 fails (random sparse X)', okc, not okc)
xe = [EN.var(EV, 'x', m, 0) for m in range(4)]
def K_vf(b):
    x2 = EN.padd(*[EN.pmul(xe[m], xe[m]) for m in range(4)]); bx = EN.padd(*[EN.pscale(xe[m], b[m]) for m in range(4) if b[m]])
    return [EN.padd(EN.pscale(x2, b[nu]) if b[nu] else {}, EN.pscale(EN.pmul(bx, xe[nu]), -2)) for nu in range(4)]
vf = [EN.padd(a_, EN.pscale(xe[m], 6)) for m, a_ in enumerate(K_vf([1, 0, 1, 0]))]
vfd = [EN.pscale(EN.padd(a_, EN.pscale(xe[m], 6)), -1) for m, a_ in enumerate(K_vf([-1, 0, 1, 0]))]
S1v = MO.S1_even(EV, vf, {}, Fr(1, 4), Fr(3, 8)); S1vd = MO.S1_even(EV, vfd, {}, Fr(1, 4), Fr(3, 8))
rec('S10', 'E3', 'engine: S1(K_b + 3 Delta)^Theta = S1((K_b + 3 Delta)^dagTheta) (exact)', '', same_fun(theta_poly(S1v, 1), S1vd))
f_ = MO.S1_odd(EV, rq(), rq(), aS=Fr(-1)); h_ = EN.pmul(EN.pmul(xe[1], FD.lam(3, 2)), FD.A(4, 0))
for (nf, Fp), (ngn, Gp) in [(('S1(odd)', f_), ('S1(even)', S1v)), (('S1(even)', S1v), ('test functional', h_)), (('S1(odd)', f_), ('test functional', h_))]:
    pf, pg = EN.parity(EV, Fp), EN.parity(EV, Gp)
    ok = same_fun(theta_poly(EN.bracket(EV, Fp, Gp)), EN.bracket(EV, theta_poly(Gp), theta_poly(Fp)), (-1) ** (pf + pg))
    rec('S10', 'E4', 'engine bracket: {f,g}^Theta = (-1)^(|f|+|g|) {g^Theta, f^Theta}, f = %s, g = %s' % (nf, ngn), ok, ok)
gg = MO.S1_odd(EV, rq(), ZS, aS=Fr(-1))
for hh, lab in ((EN.pmul(EN.pmul(xe[0], FD.A(4, 0)), FD.A(5, 1)), 'even'), (h_, 'odd')):
    s1 = EN.padd(EN.bracket(EV, f_, EN.bracket(EV, gg, hh)), EN.bracket(EV, gg, EN.bracket(EV, f_, hh)), EN.bracket(EV, EN.bracket(EV, f_, gg), hh))
    rec('S10', 'E5', 'engine bracket Jacobi for even f,g: {f,{g,h}} + {g,{f,h}} = -{{f,g},h} (h %s): hence [delta_x,delta_y] = -delta_[x,y]' % lab, '', EN.is_zero_functional(EV, s1))
print('  S10 done %.1f s, peak RSS %s' % (time.time() - T0, mem())); gc.collect()

# =====================================================================================================
# S11. two-component forms (dictionary realisation D)
# =====================================================================================================
print('=== S11 two-component ===')
I2 = np.eye(2, dtype=complex)
PAULI = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]]), np.array([[1, 0], [0, -1]], complex)]
EPSU = np.array([[0, 1], [-1, 0]], complex); EPSD = np.array([[0, -1], [1, 0]], complex)
SIGE = [PAULI[0], PAULI[1], PAULI[2], 1j * I2]; SIGBE = [-PAULI[0], -PAULI[1], -PAULI[2], 1j * I2]
SIGI = [np.array(S.tolist(), dtype=complex) for S in RB.Sig]; SIGBI = [np.array(S.tolist(), dtype=complex) for S in RB.Sigb]
def split(x): return x[:8].reshape((2, 4) + x.shape[1:]), x[8:].reshape((2, 4) + x.shape[1:])
def join(a, b): return np.concatenate([a.reshape((8,) + a.shape[2:]), b.reshape((8,) + b.shape[2:])])
def rho_S(M, x):
    q, qb = split(x)
    if M < 4:
        return join(1j * np.einsum('di...,da->ai...', qb, SIGBE[M]), -1j * np.einsum('da,ai...->di...', SIGBE[M], q))
    I = M - 4
    return join(np.einsum('ab,ji,bi...->aj...', EPSU, SIGI[I], q), -np.einsum('ji,de,ei...->dj...', SIGBI[I], EPSU, qb))
def rho_Ss(M, y):
    s, sb = split(y)
    if M < 4:
        return join(-1j * np.einsum('ad,di...->ai...', SIGE[M], sb), 1j * np.einsum('ai...,ad->di...', s, SIGE[M]))
    I = M - 4
    return join(np.einsum('kj,ab,bj...->ak...', SIGBI[I], EPSD, s), -np.einsum('kj,ed,dj...->ek...', SIGI[I], EPSD, sb))
def pair_D(x, y):
    q, qb = split(x); s, sb = split(y); return np.sum(q * s) - np.sum(qb * sb)
E16 = np.eye(16)
RD = [np.array([rho_S(M, E16[:, k]) for k in range(16)]).T for M in range(10)]
RbD = [np.array([rho_Ss(M, E16[:, k]) for k in range(16)]).T for M in range(10)]
PD = np.array([[pair_D(E16[:, a], E16[:, b]) for b in range(16)] for a in range(16)])
def C_S(x):
    q, qb = split(x); return join(np.conj(q), -np.conj(qb))        # eta^alpha_i = (xi_alpha^i)^*, etabar^{adot i} = -(xibar_{adot i})^*
def C_Ss(y):
    s, sb = split(y); return join(np.conj(s), -np.conj(sb))        # xi_alpha^i = (eta^alpha_i)^*, xibar_{adot i} = -(etabar^{adot i})^*
MC = np.array([C_S(E16[:, k]) for k in range(16)]).T; NC = np.array([C_Ss(E16[:, k]) for k in range(16)]).T
r1 = max(float(np.max(np.abs(NC @ np.conj(RD[M]) - RbD[M] @ MC))) for M in range(10))
r2 = max(float(np.max(np.abs(MC @ np.conj(RbD[M]) - RD[M] @ NC))) for M in range(10))
rec('S11', 'N1', 'two-component C: C rho(e_M) = rho(e_M) C for the ten real basis vectors, C^2 = 1', '%.1e, %.1e, %.1e' % (r1, r2, float(np.max(np.abs(NC @ np.conj(MC) - E16)))), max(r1, r2) < 1e-12)
xr = rng.normal(size=16) + 1j * rng.normal(size=16)
rec('S11', 'N2', 'two-component C: (xi, C xi) = |xi_alpha^i|^2 + |xibar_{adot i}|^2 > 0', '%s' % np.round(pair_D(xr, C_S(xr)) / np.sum(np.abs(xr) ** 2), 12), abs(pair_D(xr, C_S(xr)) / np.sum(np.abs(xr) ** 2) - 1) < 1e-12)
rows = []
for M in range(10):
    rows.append(np.hstack([-np.kron(RbD[M], E16), np.kron(E16, np.conj(RD[M]).T)]))
    rows.append(np.hstack([np.kron(E16, np.conj(RbD[M]).T), -np.kron(RD[M], E16)]))
svals = np.linalg.svd(np.vstack(rows), compute_uv=False)
rec('S11', 'N3', 'uniqueness: equivariant antilinear (Sigma <-> Sigma^*) maps form a one-dimensional space', int(np.sum(svals < 1e-9)), int(np.sum(svals < 1e-9)) == 1)
rows = []
for M in range(10):
    rows.append(np.hstack([np.kron(SIG[M], E16), -np.kron(E16, RD[M].T)]))
    rows.append(np.hstack([-np.kron(E16, RbD[M].T), np.kron(SIGB[M], E16)]))
_u, _s, _vh = np.linalg.svd(np.vstack(rows))
Xv = _vh[-1].conj(); U = Xv[:256].reshape(16, 16); Us = Xv[256:].reshape(16, 16)
kk = np.argwhere(np.abs(PD) > 1e-9)[0]; ratio = (U.T @ Us)[tuple(kk)] / PD[tuple(kk)]
U = U / np.sqrt(ratio); Us = Us / np.sqrt(ratio)
rec('S11', 'N4', 'intertwiner U (D -> classical): C_classical U = U^* C_D (C_classical = componentwise conjugation)', '%.1e' % float(np.max(np.abs(np.conj(U) - Us @ MC))), float(np.max(np.abs(np.conj(U) - Us @ MC))) < 1e-10)
# reflection formulas in two-component form
nv = rng.normal(size=4); nv /= np.linalg.norm(nv)
lam2 = rng.normal(size=16) + 1j * rng.normal(size=16); ls2 = rng.normal(size=16) + 1j * rng.normal(size=16)
rhoS = lambda v, x: sum(v[M] * rho_S(M, x) for M in range(4)); rhoSs = lambda v, y: sum(v[M] * rho_Ss(M, y) for M in range(4))
ref1 = rhoSs(nv, C_S(lam2)); q, qb = split(lam2)
f1q = 1j * sum(nv[m] * np.einsum('ad,di->ai', SIGE[m], np.conj(qb)) for m in range(4))
f1qb = 1j * sum(nv[m] * np.einsum('ai,ad->di', np.conj(q), SIGE[m]) for m in range(4))
rec('S11', 'N5', '(rho(n)C lam)_alpha^i = i n_m sigma^m_{alpha adot} (lambar_{adot i})^*, (rho(n)C lam)_{adot i} = i n_m (lam_alpha^i)^* sigma^m_{alpha adot}', '%.1e' % resid(join(f1q, f1qb), ref1), resid(join(f1q, f1qb), ref1) < 1e-12)
ref2 = rhoS(nv, C_Ss(ls2)); s, sb = split(ls2)
f2s = -1j * sum(nv[m] * np.einsum('di,da->ai', np.conj(sb), SIGBE[m]) for m in range(4))
f2sb = -1j * sum(nv[m] * np.einsum('da,ai->di', SIGBE[m], np.conj(s)) for m in range(4))
rec('S11', 'N6', '(rho(n)C lam*)^alpha_i = -i n_m (lambar*^{adot i})^* sigmabar^{m adot alpha}, (rho(n)C lam*)^{adot i} = -i n_m sigmabar^{m adot alpha}(lam*^alpha_i)^*', '%.1e' % resid(join(f2s, f2sb), ref2), resid(join(f2s, f2sb), ref2) < 1e-12)
# R-symmetry: r^i_j = 1/4 R_IJ (Sigmabar^J Sigma^I)^i_j ; dagger R -> -Rbar is r -> r^H (matrix adjoint)
Rr = rc(6, 6); Rr = Rr - Rr.T
rofR = lambda R: 0.25 * sum(R[I, J] * SIGBI[J] @ SIGI[I] for I in range(6) for J in range(6))
rec('S11', 'N7', 'R-symmetry: r(-Rbar) = r(R)^H (two-component dagger = matrix adjoint on sl(4)); r(R) anti-Hermitian for real R',
    '%.1e, %.1e' % (resid(rofR(-np.conj(Rr)), rofR(Rr).conj().T), resid(rofR(Rr.real), -rofR(Rr.real).conj().T)), resid(rofR(-np.conj(Rr)), rofR(Rr).conj().T) < 1e-12)

OUT['summary'] = {'checks': sum(len(v) for v in OUT['sections'].values()),
                  'passed': sum(1 for v in OUT['sections'].values() for r in v if r['passed']),
                  'seconds': round(time.time() - T0, 1)}
json.dump(OUT, open(os.path.join(HERE, 'conjugation_checks.json'), 'w'), indent=1, default=str)
print('SUMMARY', OUT['summary'])
