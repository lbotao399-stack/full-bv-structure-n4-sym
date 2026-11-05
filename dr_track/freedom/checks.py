#!/usr/bin/env python3
"""
checks.py -- sympy / scipy checks for freedom.tex (one-loop superconformal components, N=4 SYM).
Writes checks.json next to this file.  Light: no process above a few hundred MB, runtime ~ 1-2 min.

Blocks (labels used in freedom.tex):
  C1  diagonal derivations of psl(4|4) on the blocks = span{lambda, B}; annihilator (7-dim) and the five project relations
  C2  field-counting exact direction: shift of t_{xi+}+t_{xi-}-t_P and of c_Q(xi+)+c_Q(xi-)-p_1
  C3  operator-side combinations from eq. (c12-t) and Thm c12-F0; inner/outer invariance
  C4  old parametrisation c_Q = -2L - theta_perp - p_perp  vs  new c_Q = 79/12 - 3L - t_xi
  C5  coboundary model on one multiplet: shift of the combination = f_O - f_QO - f_QbarO + f_dO
  C6  Mellin-0 profile functionals: c_E[chi], single-theta universality, k-scaling for chi_k = exp(-u^k), simplex formula
  C7  numerical second-profile test: two model defect cycles with external momentum, heat kernel vs chi_2
  C8  gl(1|1) models: lifting class invariance, kappa*beta = gamma, rigidity of the typical module, twisted d_1
  C9  data: Konishi normalisation, rationality of the computed D_1 spectra, Beisert su(2|3) product c9*c10
"""
import json, math, os, itertools
import sympy as sp
import numpy as np
from scipy import integrate

OUT = {}
L = sp.log(sp.Rational(3, 2))


def record(key, ok, **data):
    OUT[key] = {"pass": bool(ok), **{k: (str(v) if not isinstance(v, (int, float, bool, list, dict)) else v)
                                     for k, v in data.items()}}
    print(("PASS " if ok else "FAIL ") + key)


# ----------------------------------------------------------------------------------------------- C1
blocks = ["Delta", "P", "M", "R", "xi+", "xi-", "eta+", "eta-", "K"]
lam = {"Delta": 0, "P": -2, "M": 0, "R": 0, "xi+": -1, "xi-": -1, "eta+": 1, "eta-": 1, "K": 2}
Bc = {"Delta": 0, "P": 0, "M": 0, "R": 0, "xi+": 1, "xi-": -1, "eta+": -1, "eta-": 1, "K": 0}
t = {b: sp.Symbol("t_" + b) for b in blocks}
# nonzero brackets between blocks (Prop c2-special, Thm c12-fix table): (x, y) -> list of image blocks
brackets = [("Delta", "P", ["P"]), ("Delta", "xi+", ["xi+"]), ("Delta", "xi-", ["xi-"]), ("Delta", "eta+", ["eta+"]),
            ("Delta", "eta-", ["eta-"]), ("Delta", "K", ["K"]), ("M", "P", ["P"]), ("M", "xi+", ["xi+"]),
            ("R", "xi+", ["xi+"]), ("xi+", "xi-", ["P"]), ("eta+", "eta-", ["K"]),
            ("xi+", "eta+", ["Delta", "M", "R"]), ("xi-", "eta-", ["Delta", "M", "R"]),
            ("K", "P", ["Delta", "M"]), ("P", "eta+", ["xi-"]), ("P", "eta-", ["xi+"]),
            ("K", "xi+", ["eta-"]), ("K", "xi-", ["eta+"])]
eqs = []
for x, y, imgs in brackets:
    for z in imgs:
        eqs.append(t[z] - t[x] - t[y])
sol = sp.linsolve(eqs, [t[b] for b in blocks])
sol = list(sol)[0]
free = sorted(set().union(*[s.free_symbols for s in sol]), key=str)
dim_der = len(free)
lamvec = sp.Matrix([lam[b] for b in blocks])
Bvec = sp.Matrix([Bc[b] for b in blocks])
# check lam and B solve the derivation equations
ok_lamB = all(sp.simplify(e.subs({t[b]: lam[b] for b in blocks})) == 0 for e in eqs) and \
          all(sp.simplify(e.subs({t[b]: Bc[b] for b in blocks})) == 0 for e in eqs)
# general solution lies in span{lam, B}
Msol = sp.Matrix(sol)
span_ok = sp.Matrix.hstack(lamvec, Bvec, *[Msol.diff(f) for f in free]).rank() == 2
record("C1a_diagonal_derivations_dim2_span_lambda_B", dim_der == 2 and ok_lamB and span_ok,
       dim=dim_der, general_solution=str(sol))
# annihilator of span{lam, B}
A = sp.Matrix.vstack(lamvec.T, Bvec.T)
ann = A.nullspace()
rels = {
    "t_Delta": {"Delta": 1},
    "t_M - t_Delta": {"M": 1, "Delta": -1},
    "t_R - t_Delta": {"R": 1, "Delta": -1},
    "t_xi+ + t_xi- - t_P": {"xi+": 1, "xi-": 1, "P": -1},
    "t_xi+ + t_eta+ - t_Delta": {"xi+": 1, "eta+": 1, "Delta": -1},
    "t_xi- + t_eta- - t_Delta": {"xi-": 1, "eta-": 1, "Delta": -1},
    "t_eta+ + t_eta- - t_K": {"eta+": 1, "eta-": 1, "K": -1},
}
vecs = []
ok_ann = True
for name, c in rels.items():
    v = sp.Matrix([c.get(b, 0) for b in blocks])
    vecs.append(v)
    ok_ann &= (lamvec.dot(v) == 0 and Bvec.dot(v) == 0)
rank7 = sp.Matrix.hstack(*vecs).rank()
record("C1b_annihilator_dim7_contains_project_relations", len(ann) == 7 and ok_ann and rank7 == 7,
       annihilator_dim=len(ann), rank_of_7_relations=rank7,
       note="t_xi- - t_eta+ is annihilated by B but not by lambda (shift -2c)")
record("C1c_txi-_minus_teta+_not_inner_invariant",
       sp.Matrix([1 if b == "xi-" else (-1 if b == "eta+" else 0) for b in blocks]).dot(lamvec) == -2,
       lambda_pairing=-2)

# ----------------------------------------------------------------------------------------------- C2
# F^(1)(xi)|_2 = p (T1+T2) + t (T1+T2+T3+T4); field counting y: arrows f->f' multiplied by (y_f' - y_f)
# T1: A_V -> lambda, T2: A_perp -> lambda, T3: lambda -> A_V, T4: lambda -> A_perp (Prop c11-freedom (i))
a, yl, yA, p0, t0 = sp.symbols("a y_lambda y_A p0 t0")
T = sp.symbols("T1:5")
F1 = sp.expand(p0 * (T[0] + T[1]) + t0 * (T[0] + T[1] + T[2] + T[3]))
shift = a * ((yl - yA) * T[0] + (yl - yA) * T[1] + (yA - yl) * T[2] + (yA - yl) * T[3])
F1n = sp.expand(F1 + shift)
# read new (p,t): t = coeff of T3, p + t = coeff of T1
tn = F1n.coeff(T[2])
pn = sp.expand(F1n.coeff(T[0]) - tn)
ok_t = sp.simplify(tn - (t0 - a * (yl - yA))) == 0 and sp.simplify(pn - (p0 + 2 * a * (yl - yA))) == 0
# combination t_xi+ + t_xi- - t_P : both chiralities shift equally, P unchanged
comb_shift = sp.simplify(2 * (tn - t0) - 0)
# operator side: c_Q = (two-body -2L) - (coefficient of T2 in F^(1)(xi))  (one-body -{T2,O} = -B1)
cQ_old = -2 * L - F1.coeff(T[1])
cQ_new = -2 * L - F1n.coeff(T[1])
dcQ = sp.simplify(cQ_new - cQ_old)
record("C2_field_counting_moves_claimed_invariants",
       ok_t and sp.simplify(comb_shift + 2 * a * (yl - yA)) == 0 and sp.simplify(dcQ + a * (yl - yA)) == 0,
       t_xi_shift=str(sp.simplify(tn - t0)), particular_part_shift=str(sp.simplify(pn - p0)),
       action_combination_shift=str(comb_shift), cQ_shift=str(dcQ),
       operator_combination_shift=str(sp.simplify(2 * dcQ)))

# ----------------------------------------------------------------------------------------------- C3
td, tP, txp, txm, tep = sp.symbols("t_Delta t_P t_xip t_xim t_etap")
cQ = lambda tx: sp.Rational(79, 12) - 3 * L - tx
cS = lambda te: sp.Rational(79, 12) - 3 * L - te
p1 = -(sp.Rational(1, 3) + tP)
D1tl = -4 * td - sp.Rational(4, 3)
rel1 = {txm: sp.Rational(27, 4) + 2 * L + tP - txp}           # t_xi+ + t_xi- - t_P = 27/4 + 2L
rel2 = {tep: sp.Rational(27, 4) + 2 * L + td - txp}           # t_xi+ + t_eta+ - t_Delta = 27/4 + 2L
comb1 = sp.simplify((cQ(txp) + cQ(txm) - p1).subs(rel1))
comb2 = sp.simplify((cQ(txp) + cS(tep) - D1tl / 4).subs(rel2))
target = sp.Rational(27, 4) - 8 * L
c_, a_ = sp.symbols("c a_out")
# inner: t_z -> t_z + c lambda_z ; outer: t_xi+ -> +a, t_xi- -> -a, t_eta+ -> -a
inner = {txp: txp - c_, txm: txm - c_, tP: tP - 2 * c_, tep: tep + c_, td: td}
outer = {txp: txp + a_, txm: txm - a_, tep: tep - a_}
e1 = cQ(txp) + cQ(txm) - p1
e2 = cQ(txp) + cS(tep) - D1tl / 4
ok_inv = all(sp.simplify(e.subs(s, simultaneous=True) - e) == 0 for e in (e1, e2) for s in (inner, outer))
D1tr = -4 * td - sp.Rational(148, 3)
record("C3_operator_combinations_and_inner_outer_invariance",
       sp.simplify(comb1 - target) == 0 and sp.simplify(comb2 - target) == 0 and ok_inv and
       D1tl.subs(td, -sp.Rational(1, 3)) == 0 and D1tr.subs(td, -sp.Rational(1, 3)) == -48,
       comb1=str(comb1), comb2=str(comb2), D1_traceless=str(D1tl.subs(td, -sp.Rational(1, 3))),
       D1_trace=str(D1tr.subs(td, -sp.Rational(1, 3))))

# ----------------------------------------------------------------------------------------------- C4
thp, pp = sp.symbols("theta_perp p_perp")
cQ_oldpar = -2 * L - thp - pp
ok4 = sp.simplify(cQ_oldpar.subs({pp: L - sp.Rational(79, 12), thp: txp}) - cQ(txp)) == 0
record("C4_old_vs_new_parametrisation_of_cQ", ok4, p_perp="L - 79/12", theta_perp="t_xi")

# ----------------------------------------------------------------------------------------------- C5
fO, fQ, fQb, fdO = sp.symbols("f_O f_QO f_QbarO f_dO")
# coboundary delta = [rho_0, phi], phi = f on each conformal primary/descendant line used:
dcQp = fO - fQ          # c_Q(xi+) shift
dcQm = fO - fQb         # c_Q(xi-) shift
dp1 = fO - fdO          # p_1 shift
shift5 = sp.simplify(dcQp + dcQm - dp1)
# purely weight-dependent phi = f(w): f_QO = f_QbarO = f(5/2), f_dO = f(3)
f2, f52, f3 = sp.symbols("f2 f52 f3")
shift5w = sp.simplify(shift5.subs({fO: f2, fQ: f52, fQb: f52, fdO: f3}))
# translation covariant phi: f_dO = f_O
shift5t = sp.simplify(shift5.subs({fdO: fO}))
record("C5_coboundary_shift_of_combination",
       shift5 == fO - fQ - fQb + fdO and shift5w == f2 - 2 * f52 + f3 and shift5t == 2 * fO - fQ - fQb,
       general=str(shift5), weight_function=str(shift5w), translation_covariant=str(shift5t))

# ----------------------------------------------------------------------------------------------- C6
u = sp.Symbol("u", positive=True)


def mellin0(expr):
    """int_0^oo du/u expr, numeric (mpmath)."""
    f = sp.lambdify(u, expr / u, "mpmath")
    import mpmath
    mpmath.mp.dps = 30
    return mpmath.quad(f, [0, 1, 10, mpmath.inf])


import mpmath
mpmath.mp.dps = 30
res6 = {}
ok6 = True
profiles = {"heat": sp.exp(-u), "chi2": sp.exp(-u ** 2), "PV": 1 / (1 + u), "chi3": sp.exp(-u ** 3)}
pred_cE = {"heat": lambda E: math.log(E / (E + 1)), "chi2": lambda E: 0.5 * math.log(E / (E + 1)),
           "chi3": lambda E: math.log(E / (E + 1)) / 3, "PV": lambda E: -1.0 / E}
for name, ch in profiles.items():
    for E in (1, 2, 3):
        val = float(mellin0(ch ** E * (ch - 1)))
        pr = pred_cE[name](E)
        res6[f"c_{E}[{name}]"] = [val, pr]
        ok6 &= abs(val - pr) < 1e-12
record("C6a_profile_functional_c_E", ok6, values=res6)
# single-theta universality: int du/u (theta chi) chi^a = -1/(a+1)
res6b = {}
ok6b = True
for name, ch in profiles.items():
    th = u * sp.diff(ch, u)
    for aa in (0, 1, 2):
        val = float(mellin0(th * ch ** aa))
        res6b[f"{name},a={aa}"] = val
        ok6b &= abs(val + 1.0 / (aa + 1)) < 1e-12
record("C6b_single_theta_functionals_universal", ok6b, values=res6b)
# multi-theta: (theta chi)^2 chi^a for chi_k equals k * heat-kernel value; heat value rational 1/(2+a)^2
res6c = {}
ok6c = True
for k in (2, 3):
    ch = sp.exp(-u ** k)
    th = u * sp.diff(ch, u)
    for aa in (0, 1):
        val = float(mellin0(th ** 2 * ch ** aa))
        heat = 1.0 / (2 + aa) ** 2
        res6c[f"k={k},a={aa}"] = [val, k * heat]
        ok6c &= abs(val - k * heat) < 1e-12
record("C6c_two_theta_functionals_scale_like_k", ok6c, values=res6c)
# simplex formula of Lemma c12-special at constant F_0: int_{alpha_sigma smallest} log(min_f alpha_f/alpha_sigma) = log(n/(n-1))/(n-1)!
s2 = mpmath.quad(lambda x: mpmath.log((1 - x) / x), [0, 0.5])
s3 = mpmath.quad(lambda x1, x2: (mpmath.log(min(x2, 1 - x1 - x2) / x1) if (x1 <= x2 and x1 <= 1 - x1 - x2) else 0),
                 [0, 1.0 / 3, 0.5, 1], [0, 0.5, 1]) if False else None
# region-split exact 2-simplex integral by symbolic integration
x1, x2 = sp.symbols("x1 x2", positive=True)
# alpha_1 smallest, alpha_2 <= alpha_3 (factor 2 by symmetry): x1 in [0,1/3], x2 in [x1, (1-x1)/2]
I3 = 2 * sp.integrate(sp.integrate(sp.log(x2 / x1), (x2, x1, (1 - x1) / 2)), (x1, 0, sp.Rational(1, 3)))
I3 = sp.nsimplify(sp.simplify(I3))
ok6d = abs(float(s2) - math.log(2)) < 1e-15 and sp.simplify(I3 - sp.log(sp.Rational(3, 2)) / 2) == 0
record("C6d_simplex_formula_nt2_nt3", ok6d, nt2=float(s2), nt3=str(I3),
       relation="(n_t-1)! * simplex value = -c_{n_t-1}[exp(-u)]")

# ----------------------------------------------------------------------------------------------- C7
# J0(q) = (4/pi) int_0^oo l^3 dl int_0^pi sin^2 th dth (chi(l^2)-1) chi(|l+q|^2) / (l^2 |l+q|^2)  -> c_1[chi] as q->0
# J2(q) = (4/pi) int ... (chi(l^2)-1) chi(|l+q|^2)/|l+q|^2 = A + B q^2 + O(q^4),
#    A = int du (chi-1) chi,  B = (1/2) int du u chi''(u) (chi(u)-1)   (q = sqrt(eps) p, scaled form)
def make_chi(name):
    if name == "heat":
        return lambda x: math.exp(-x)
    if name == "chi2":
        return lambda x: math.exp(-x * x)


def J(name, q, kind):
    chi = make_chi(name)

    def inner(l):
        def ang(th):
            k2 = l * l + 2 * l * q * math.cos(th) + q * q
            if k2 <= 0:
                return 0.0
            base = (chi(l * l) - 1.0) * chi(k2) / k2
            if kind == 0:
                base /= (l * l)
            return math.sin(th) ** 2 * base
        v, _ = integrate.quad(ang, 0, math.pi, limit=200, epsabs=1e-13, epsrel=1e-12)
        return l ** 3 * v

    pts = [q] if q > 0 else None
    v1, _ = integrate.quad(inner, 0, 2.0, points=pts, limit=400, epsabs=1e-13, epsrel=1e-12)
    v2, _ = integrate.quad(inner, 2.0, 12.0, limit=400, epsabs=1e-13, epsrel=1e-12)
    return 4.0 / math.pi * (v1 + v2)


res7 = {}
ok7 = True
predA = {"heat": -0.5, "chi2": -0.5 * math.sqrt(math.pi) * (1 - 1 / math.sqrt(2))}  # int (chi-1)chi du
predB = {"heat": -3.0 / 8.0, "chi2": -0.5}
predJ0 = {"heat": -math.log(2.0), "chi2": -0.5 * math.log(2.0)}
for name in ("heat", "chi2"):
    qs = [0.08, 0.12, 0.16, 0.2]
    vals2 = [J(name, q, 2) for q in qs]
    # fit A + B q^2 + C q^4
    M = np.array([[1, q ** 2, q ** 4] for q in qs])
    coef, *_ = np.linalg.lstsq(M, np.array(vals2), rcond=None)
    j0 = [J(name, q, 0) for q in (0.05, 0.1)]
    # J0 -> c_1 with O(q^2 log q) corrections: Richardson-type estimate
    res7[name] = {"A_fit": coef[0], "A_pred": predA[name], "B_fit": coef[1], "B_pred": predB[name],
                  "J0(q=0.05)": j0[0], "J0(q=0.1)": j0[1], "J0_pred(q->0)": predJ0[name]}
    ok7 &= abs(coef[0] - predA[name]) < 2e-6 and abs(coef[1] - predB[name]) < 2e-3 and abs(j0[0] - predJ0[name]) < 2e-2
record("C7_second_profile_numerical_test", ok7, results=res7,
       statement="D=0 log part scales 1/k (-log2 -> -log2/2); D=2 rational part moves (-3/8 -> -1/2)")

# ----------------------------------------------------------------------------------------------- C8
# gl(1|1): N even, Z even central, Q, S odd; [N,Q]=Q, [N,S]=-S, {Q,S}=Z.  Matrices on a Z2-graded space.
def scom(X, Y, px, py):
    return X * Y - (-1) ** (px * py) * Y * X


gens = ["N", "Z", "Q", "S"]
par = {"N": 0, "Z": 0, "Q": 1, "S": 1}


def br(x, y):
    """structure constants: returns dict gen->coeff of [x,y]"""
    table = {("N", "Q"): {"Q": 1}, ("Q", "N"): {"Q": -1}, ("N", "S"): {"S": -1}, ("S", "N"): {"S": 1},
             ("Q", "S"): {"Z": 1}, ("S", "Q"): {"Z": 1}}
    return table.get((x, y), {})


def is_rep(rho):
    ok = True
    for x in gens:
        for y in gens:
            lhs = scom(rho[x], rho[y], par[x], par[y])
            rhs = sp.zeros(*lhs.shape)
            for z, cz in br(x, y).items():
                rhs += cz * rho[z]
            ok &= (lhs - rhs).applyfunc(sp.simplify) == sp.zeros(*lhs.shape)
    return ok


def cocycle_defect(rho0, dl):
    """delta([x,y]) - [delta(x),rho0(y)] - [rho0(x),delta(y)] for all pairs"""
    bad = []
    for x in gens:
        for y in gens:
            lhs = sp.zeros(*rho0["N"].shape)
            for z, cz in br(x, y).items():
                lhs += cz * dl[z]
            rhs = scom(dl[x], rho0[y], par[x], par[y]) + scom(rho0[x], dl[y], par[x], par[y])
            if (lhs - rhs).applyfunc(sp.simplify) != sp.zeros(*lhs.shape):
                bad.append((x, y))
    return bad


n_ = sp.Symbol("n")
gam, kap, h = sp.symbols("gamma kappa hbar")
# atypical extension E: basis (v0 even, v1 odd), N=diag(n, n+1), Q v0 = v1 (beta), S = 0, Z = 0
rho0 = {"N": sp.diag(n_, n_ + 1), "Z": sp.zeros(2), "Q": sp.Matrix([[0, 0], [1, 0]]), "S": sp.zeros(2)}
# one-loop: delta(S) v1 = kappa v0, delta(Z) = gamma Id, delta(Q) = cq (v0 -> v1)
cq = sp.Symbol("c_q")
dl = {"N": sp.zeros(2), "Z": gam * sp.eye(2), "Q": cq * sp.Matrix([[0, 0], [1, 0]]),
      "S": kap * sp.Matrix([[0, 1], [0, 0]])}
bad = cocycle_defect(rho0, dl)
# cocycle iff kappa = gamma (beta = 1)
bad_sub = cocycle_defect(rho0, {k: v.subs(kap, gam) for k, v in dl.items()})
# coboundary phi = diag(f0, f1): delta -> delta + [rho0, phi]
f0, f1 = sp.symbols("f0 f1")
phi = sp.diag(f0, f1)
cob = {x: scom(rho0[x], phi, par[x], 0) for x in gens}
shifts = {x: cob[x] for x in gens}
ok8a = is_rep(rho0) and bad != [] and bad_sub == [] and shifts["S"] == sp.zeros(2) and shifts["Z"] == sp.zeros(2) \
       and shifts["Q"] == (f0 - f1) * sp.Matrix([[0, 0], [1, 0]])
record("C8a_lifting_class_invariant_and_kappa_beta_equals_gamma", ok8a,
       cocycle_condition="kappa*beta = gamma (beta=1)", coboundary_moves="only c_q (one-loop correction of beta)")
# rigidity of the typical module V_{n,z}, z != 0: relative cocycles (delta(N)=delta(Z)=0) are coboundaries
z = sp.Symbol("z", nonzero=True)
rhoT = {"N": sp.diag(n_, n_ + 1), "Z": z * sp.eye(2), "Q": sp.Matrix([[0, 0], [1, 0]]),
        "S": z * sp.Matrix([[0, 1], [0, 0]])}
qa, sa = sp.symbols("q_a s_a")
dlT = {"N": sp.zeros(2), "Z": sp.zeros(2), "Q": qa * sp.Matrix([[0, 0], [1, 0]]), "S": sa * sp.Matrix([[0, 1], [0, 0]])}
badT = cocycle_defect(rhoT, dlT)
# condition for cocycle: solve
conds = []
for x in gens:
    for y in gens:
        lhs = sp.zeros(2)
        for zz, cz in br(x, y).items():
            lhs += cz * dlT[zz]
        rhs = scom(dlT[x], rhoT[y], par[x], par[y]) + scom(rhoT[x], dlT[y], par[x], par[y])
        conds += list((lhs - rhs).applyfunc(sp.expand))
conds = [c for c in conds if c != 0]
solT = sp.solve(conds, [qa, sa], dict=True)
# coboundaries of diag(f0,f1): Q shifts by (f0-f1), S shifts by -(f0-f1) z ... check that every relative cocycle is one
cobT = {x: scom(rhoT[x], phi, par[x], 0) for x in gens}
is_cob = True
for s_ in solT:
    dQ = dlT["Q"].subs(s_)
    dS = dlT["S"].subs(s_)
    eqf = list((dQ - cobT["Q"]).applyfunc(sp.expand)) + list((dS - cobT["S"]).applyfunc(sp.expand))
    eqf = [e for e in eqf if e != 0]
    is_cob &= len(sp.solve(eqf, [f0, f1], dict=True)) > 0 or eqf == []
# dimension deformation of the typical module: delta(Z) = Id, delta(S) = (v1 -> v0), others 0 is a cocycle
dlTd = {"N": sp.zeros(2), "Z": sp.eye(2), "Q": sp.zeros(2), "S": sp.Matrix([[0, 1], [0, 0]])}
ok8b = is_rep(rhoT) and is_cob and cocycle_defect(rhoT, dlTd) == []
record("C8b_rigidity_typical_module_and_dimension_direction", ok8b, relative_cocycles=str(solT))
# atypical 1-dim module (BPS analogue): delta(Z)=gamma forces gamma=0
rho1 = {"N": sp.Matrix([[n_]]), "Z": sp.zeros(1), "Q": sp.zeros(1), "S": sp.zeros(1)}
dl1 = {"N": sp.zeros(1), "Z": sp.Matrix([[gam]]), "Q": sp.zeros(1), "S": sp.zeros(1)}
bad1 = cocycle_defect(rho1, dl1)
record("C8c_protection_atypical_simple_module", ("Q", "S") in bad1 and cocycle_defect(rho1, {k: v.subs(gam, 0) for k, v in dl1.items()}) == [],
       statement="delta(Z)=gamma Id on a 1-dim atypical module is a cocycle iff gamma=0")
# twisted d_1: M = C_n(even) + C_{n+1}(odd), rho0 odd = 0; delta(Q)= kQ (v0->v1), delta(S) = kS (v1->v0)
kQ, kS = sp.symbols("k_Q k_S")
rhoM = {"N": sp.diag(n_, n_ + 1), "Z": sp.zeros(2), "Q": sp.zeros(2), "S": sp.zeros(2)}
dlM = {"N": sp.zeros(2), "Z": sp.zeros(2), "Q": kQ * sp.Matrix([[0, 0], [1, 0]]), "S": kS * sp.Matrix([[0, 1], [0, 0]])}
okM1 = cocycle_defect(rhoM, dlM) == []
# second order: {rho(Q), rho(S)} = rho(Z) at order hbar^2 gives delta_2(Z) = kQ kS Id
rhoQ = h * dlM["Q"]
rhoS = h * dlM["S"]
anti = (rhoQ * rhoS + rhoS * rhoQ).applyfunc(sp.expand)
okM2 = anti == h ** 2 * kQ * kS * sp.eye(2)
# d_1 on H(rho0(Q)) = M is delta(Q): coboundaries [rho0(Q), phi] = 0 since rho0(Q) = 0
okM3 = scom(rhoM["Q"], phi, 1, 0) == sp.zeros(2)
record("C8d_twisted_d1_canonical_and_two_loop_dimension", okM1 and okM2 and okM3,
       d1="v0 -> k_Q v1 (canonical)", two_loop="delta_2(Z) = k_Q k_S")

# ----------------------------------------------------------------------------------------------- C9
g2N = sp.Symbol("g2N", positive=True)
unit = -g2N / (32 * sp.pi ** 2)
D1K = -48 * unit
gamK = 3 * g2N / (4 * sp.pi ** 2)
ok9a = sp.simplify(D1K - 2 * gamK) == 0
spectra = {"traceless": 0, "trace": -48, "twist2": [-32 * sp.harmonic(j) for j in (0, 2, 4, 6)],
           "su2_heisenberg_unit": -8}
ok9b = all(sp.nsimplify(v).is_rational for v in [0, -48] + spectra["twist2"])
al1, b1, b2 = sp.symbols("alpha1 beta1 beta2", real=True)
c9 = al1 * sp.exp(sp.I * (b1 + b2)) / sp.sqrt(2)
c10 = al1 * sp.exp(-sp.I * (b1 + b2)) / sp.sqrt(2)
ok9c = sp.simplify(c9 * c10 - al1 ** 2 / 2) == 0
record("C9_data_Konishi_rationality_Beisert_product", ok9a and ok9b and ok9c,
       Konishi="D_1 = -48 u = 3 g^2 N/(2 pi^2) = 2 gamma_K", twist2_values=str(spectra["twist2"]),
       beisert="c9 c10 = alpha1^2/2 = c1/2 independent of beta1, beta2 [hep-th/0310252 (3.36)-(3.38)]")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "checks.json"), "w") as fh:
    json.dump(OUT, fh, indent=1, default=str)
print("ALL PASS" if all(v["pass"] for v in OUT.values()) else "SOME FAIL")
