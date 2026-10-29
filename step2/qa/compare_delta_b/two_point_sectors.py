"""compare_delta_b / two_point_sectors.py

Exact (sympy) Q-images of engine A's per-graph arity-one two-point functionals along Delta, in all four sectors,
assembled independently of both engines' own assembly code, and compared with
  - the theory note (delta_theory.tex thm:d-D2, ex:d-quadratic; delta_predictions.json):
        Q(class ii) -> 4 C_log  = unit[-4 int(dA,dA) + 16 int(lambda,dslash lambda) - 4 int(dc,A*)],
        Q Pi_Psi    = unit[(1/3) int(dA,dA) - (1/3) int(lambda,dslash lambda) - (2/3) int(dc,A*)],
        coefficient of int(d*A_V,d*A_V) in Ob(Delta) = 0,
        translation-invariant part of Pi_Psi = -(2/3) int(A,A*) + (1/3) int(lambda*,lambda) + (4/3) int(c*,c),
        Pi_Psi local (def:d-PiPsi);
  - engine B: Ob(Delta) = unit[-(11/3) int(dA,dA) + (47/3) int(lambda,dslash lambda) - (14/3) int(dc,A*)],
              F(Delta) = unit[(22/3) int(A,A*) - (47/3) int(lambda*,lambda) - (8/3) int(c*,c)].
Canonical form (both engines): a two-leg functional int int X(p1) Y(p2) [alpha(p1) delta + beta^mu(p1) d_{p2 mu} delta],
first leg momentum p = p1 fixed (here total momentum 0).
Q maps (conventions.py, from the ch01 bracket and the action): Q A* = -(1/2) d*dA, Q A = -dc, Q c* = d*A*,
Q lambda* = dslash lambda, with the Koszul sign (-1)^{|X|} when Q passes the first leg.
Identical output legs are (graded) symmetrised: G_sym = (1/2)(G +- G_swap), swap rule
alpha''(p) = alpha(-p)^T - (div beta^T)(-p), beta''(p) = beta(-p)^T.
Output: two_point_sectors.json
"""
import json, os
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
p = sp.symbols('p0:4', real=True)
P4 = sp.Matrix(p)
pp = sum(c * c for c in p)
le, gE, L2c, L3c = sp.symbols('logeps gamma_E log2 log3')
lp = sp.log(pp)
L2 = lp + gE + L2c
L3 = lp + gE + L3c
L3e = L3 + le
I = sp.I
R = sp.Rational


def vecsym(name):
    return sp.Matrix(sp.symbols(name + '0:4'))


a, b, ap = vecsym('a'), vecsym('b'), vecsym('c')        # A_V at p1, A*_V at p2, second A_V
gam, w = sp.symbols('gam w')                           # c polarisation, c* polarisation
dot = lambda x, y: sum(x[i] * y[i] for i in range(4))


def subs_vec(expr, v, newv):
    return sp.sympify(expr).subs({v[i]: newv[i] for i in range(4)}, simultaneous=True)


def div(beta):
    return sum(sp.diff(beta[i], p[i]) for i in range(4))


def at(expr, q):
    """evaluate a function of p at p -> q (q a list of expressions)"""
    return sp.sympify(expr).subs({p[i]: q[i] for i in range(4)}, simultaneous=True)


# ----------------------------------------------------------------------------------------------- scalar/vector engine
def Q_second_leg(alpha, beta, vsym, qfun, sX):
    """second leg (polarisation vector vsym, linear) -> q(p2) z ; returns canonical (alpha_G, beta_G) at p2* = -p.
    qfun(p2list) returns the replacement vector for vsym (linear in the new polarisation)."""
    p2s = [-x for x in p]
    t = sp.symbols('t0:4')
    rep_star = qfun(p2s)
    aG = subs_vec(alpha, vsym, rep_star)
    bG = [subs_vec(beta[m], vsym, rep_star) for m in range(4)]
    # - sum_mu beta^mu(v = d_{p2 mu} q(p2) z) at p2*
    corr = 0
    for mu in range(4):
        rep_t = qfun(list(t))
        drep = [sp.diff(rep_t[i], t[mu]) for i in range(4)]
        drep = [d.subs({t[j]: p2s[j] for j in range(4)}, simultaneous=True) for d in drep]
        corr += subs_vec(beta[mu], vsym, sp.Matrix(drep))
    return sp.expand(sX * (aG - corr)), [sp.expand(sX * x) for x in bG]


def Q_first_leg(alpha, beta, usym, rvec):
    """first leg polarisation usym -> rvec (function of p1 = p): no conversion"""
    return sp.expand(subs_vec(alpha, usym, rvec)), [sp.expand(subs_vec(x, usym, rvec)) for x in beta]


def sym_bosonic(alpha, beta, u, up):
    """symmetrise over (p, u) <-> (p_d = -p, up)"""
    sw = {}
    for i in range(4):
        sw[u[i]] = up[i]; sw[up[i]] = u[i]
    mp_ = [-x for x in p]
    a_sw = at(sp.sympify(alpha).subs(sw, simultaneous=True), mp_)
    bsw = [sp.sympify(x).subs(sw, simultaneous=True) for x in beta]
    divb = at(div(bsw), mp_)
    b_sw = [at(x, mp_) for x in bsw]
    aS = sp.expand((alpha + a_sw - divb) / 2)
    bS = [sp.expand((beta[m] + b_sw[m]) / 2) for m in range(4)]
    return aS, bS


def clean(e):
    e = sp.expand(e)
    e = e.subs(sp.log(pp), sp.Symbol('LOGPP'))
    return sp.simplify(sp.expand(e))


def fit(alpha, basis):
    """alpha = sum c_i basis_i ; returns dict and residual"""
    cs = sp.symbols('x0:%d' % len(basis))
    eqexpr = sp.expand(sp.together(alpha - sum(c * bb for c, bb in zip(cs, basis))))
    num = sp.numer(sp.together(eqexpr))
    num = sp.expand(num.subs(sp.log(pp), sp.Symbol('LOGPP')))
    allv = [sp.Symbol('LOGPP'), le, gE, L2c, L3c] + list(p) + list(a) + list(b) + list(ap) + [gam, w]
    poly = sp.Poly(num, *[v for v in allv if num.has(v)]) if num != 0 else None
    if poly is None:
        return {}, '0'
    eqs = poly.coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    if not sol:
        return None, str(sp.simplify(eqexpr))
    s = sol[0]
    return {str(c): str(sp.nsimplify(s.get(c, 0))) for c in cs}, '0'


out = {}

# ------------------------------------------------------------------------------------------------ sector A|_Wperp
# 1-dim flavour: (a,b) = a0 b0 ; second A polarisation c0
ab = a[0] * b[0]
data_P = {
    'Pi_int': (R(8, 3) * ab, [R(10, 3) * p[m] * ab for m in range(4)]),
    'ext_leg': (2 * (R(43, 4) - 4 * L2) * ab, [2 * (R(43, 4) - 4 * L2) * p[m] * ab for m in range(4)]),
    'S_M|_perp': (2 * ab, [2 * p[m] * ab for m in range(4)]),
    'int(A,A*)': (ab, [0] * 4),
    'int(EA,A*)': (0, [p[m] * ab for m in range(4)]),
}
qP = lambda p2: sp.Matrix([-R(1, 2) * sum(x * x for x in p2) * ap[0], 0, 0, 0])
basis_P = [pp * a[0] * ap[0]]          # int(dA_perp, dA_perp)
resP = {}
for name, (al, be) in data_P.items():
    aG, bG = Q_second_leg(al, be, b, qP, +1)
    aS, bS = sym_bosonic(aG, bG, a, ap)
    f, resid = fit(aS, basis_P)
    resP[name] = {'Q_sym in units of int(dA_perp,dA_perp)': f.get('x0') if f else None,
                  'beta_sym': [str(clean(x)) for x in bS], 'residual': resid}
out['A_perp'] = resP

# ------------------------------------------------------------------------------------------------ sector A|_V and c
abV = dot(a, b); pa = dot(P4, a); pb = dot(P4, b)
Lo = pa * pb / pp
data_V = {
    'Pi_int': ((abV) * (R(13, 6) + L3e) + 2 * Lo,
               [R(10, 3) * p[m] * abV + (L3e - R(1, 2)) * a[m] * pb - 2 * p[m] * Lo for m in range(4)]),
    'ext_leg_1 (S_M A*_V external, A_V leg on W_1^AV)': (
        2 * (abV * (R(43, 4) - 4 * L2) + Lo * (4 * L2 - 4)),
        [2 * (abV * (R(43, 4) - 4 * L2) + Lo * (4 * L2 - 4)) * p[m] for m in range(4)]),
    'ext_leg_2 (S_M A_V external, A*_V leg on the ghost self-energy)': (
        (2 + 2 * L2) * Lo, [(2 * L2 - 2) * p[m] * Lo for m in range(4)]),
    'F0_graph (S_M - P - (-F_0), F_0 = (27/16) int(d*A_V,d*A_V))': (
        -R(27, 2) * Lo, [-R(27, 2) * p[m] * Lo for m in range(4)]),
    'S_M|_V': (2 * abV, [2 * p[m] * abV for m in range(4)]),
    'int(A_V,A*_V)': (abV, [0] * 4),
    'X_V (log eps coefficient, engine A)': (abV, [a[m] * pb for m in range(4)]),
}
gw = gam * w
data_c = {
    'Pi_int': (-R(8, 3) * gw, [(-R(5, 6) - L3e) * p[m] * gw for m in range(4)]),
    'ext_leg': (0 * gw, [(2 - 2 * L2) * p[m] * gw for m in range(4)]),
    'S_M|_c': (0 * gw, [-2 * p[m] * gw for m in range(4)]),
    'int(c*,c)': (gw, [0] * 4),
    'X_c (log eps coefficient, engine A)': (0 * gw, [-p[m] * gw for m in range(4)]),
}
# (A_V,A_V) output: Q on the A*_V leg: b -> -(1/2)(|p2|^2 c - p2 (p2.c))
qV = lambda p2: sp.Matrix([-R(1, 2) * (sum(x * x for x in p2) * ap[i] - p2[i] * sum(p2[j] * ap[j] for j in range(4))) for i in range(4)])
basis_VV = [pp * dot(a, ap) - dot(P4, a) * dot(P4, ap), dot(P4, a) * dot(P4, ap)]   # int(dA,dA), int(d*A,d*A)
# (c,A*) output: from V: first leg A_V -> -dc : a -> -i p gam ; from c: second leg c* -> d*A* = -i p2.A*, sign (-1)^{|c|} = -1
rV = sp.Matrix([-I * p[i] * gam for i in range(4)])


def qc_fun(p2):
    # replacement for the scalar w: -i p2.b ; returned as a 4-vector in slot 0 for the generic routine
    return sp.Matrix([-I * sum(p2[j] * b[j] for j in range(4)), 0, 0, 0])


basis_cA = [I * pb * gam]          # int(dc, A*)
resV = {}; cA_parts = {}
for name, (al, be) in data_V.items():
    aG, bG = Q_second_leg(al, be, b, qV, +1)
    aS, bS = sym_bosonic(aG, bG, a, ap)
    f, resid = fit(aS, basis_VV)
    resV[name] = {'Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]': [f.get('x0'), f.get('x1')] if f else None,
                  'beta_sym_zero': all(clean(x) == 0 for x in bS), 'residual': resid}
    aC, bC = Q_first_leg(al, be, a, rV)
    cA_parts['V:' + name] = (aC, bC)
wvec = sp.Matrix([w, 0, 0, 0])
for name, (al, be) in data_c.items():
    # treat w as slot 0 of a vector for the generic second-leg routine
    aG, bG = Q_second_leg(al, be, wvec, qc_fun, -1)
    cA_parts['c:' + name] = (aG, bG)
res_cA = {}
for name, (aC, bC) in cA_parts.items():
    f, resid = fit(aC, basis_cA)
    res_cA[name] = {'alpha in units of int(dc,A*)': f.get('x0') if f else None, 'residual': resid,
                    'beta': [str(clean(x)) for x in bC]}
out['A_V_(A_V,A_V)_output'] = resV
out['(c,A*)_output_per_piece'] = res_cA


def total_cA(names):
    aT = sum(cA_parts[n][0] for n in names)
    bT = [sum(cA_parts[n][1][m] for n in names) for m in range(4)]
    f, resid = fit(aT, basis_cA)
    return {'alpha in units of int(dc,A*)': f.get('x0') if f else None, 'residual': resid,
            'beta': [str(clean(x)) for x in bT], 'beta_zero': all(clean(x) == 0 for x in bT)}


out['(c,A*)_output_totals'] = {
    'Q(X_V + X_c) (log eps coefficient; must be 0)': total_cA(['V:X_V (log eps coefficient, engine A)', 'c:X_c (log eps coefficient, engine A)']),
    'Q(S_M|_V + S_M|_c) (must be 0)': total_cA(['V:S_M|_V', 'c:S_M|_c']),
    'Q(int(A_V,A*_V)) (expect -1)': total_cA(['V:int(A_V,A*_V)']),
    'Q(int(c*,c)) (expect -1)': total_cA(['c:int(c*,c)']),
    'Q(Pi_int V + c)': total_cA(['V:Pi_int', 'c:Pi_int']),
    'Q(class ii) without F_0': total_cA(['V:ext_leg_1 (S_M A*_V external, A_V leg on W_1^AV)',
                                         'V:ext_leg_2 (S_M A_V external, A*_V leg on the ghost self-energy)', 'c:ext_leg']),
    'Q(class ii) with F_0': total_cA(['V:ext_leg_1 (S_M A*_V external, A_V leg on W_1^AV)',
                                      'V:ext_leg_2 (S_M A_V external, A*_V leg on the ghost self-energy)', 'c:ext_leg',
                                      'V:F0_graph (S_M - P - (-F_0), F_0 = (27/16) int(d*A_V,d*A_V))']),
}

# ------------------------------------------------------------------------------------------------ sector lambda
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]]); one = sp.eye(2)
sig = [s1, s2, s3, I * one]                 # sigma^m_{alpha alphadot}
sigb = [-s1, -s2, -s3, I * one]             # sigmabar^{m alphadot alpha} (row alphadot, column alpha)


def vs(v, mats):
    return sum((v[m] * mats[m] for m in range(4)), sp.zeros(2, 2))


def rho_S_to_Ss(v):
    """Sigma (xi_alpha, xibar_alphadot) -> Sigma* (eta^alpha, etabar^alphadot), ch08 def:c8-twocomp"""
    M = sp.zeros(4, 4)
    sb = vs(v, sigb)
    # (rho xi)^alpha = i v_m xibar_alphadot sigmabar^{m alphadot alpha}
    for al in range(2):
        for ad in range(2):
            M[al, 2 + ad] = I * sb[ad, al]
    # (rho xi)^alphadot = -i v_m sigmabar^{m alphadot alpha} xi_alpha
    for ad in range(2):
        for al in range(2):
            M[2 + ad, al] = -I * sb[ad, al]
    return M


def rho_Ss_to_S(v):
    M = sp.zeros(4, 4)
    s = vs(v, sig)
    # (rho eta)_alpha = -i v_m sigma^m_{alpha alphadot} etabar^alphadot
    for al in range(2):
        for ad in range(2):
            M[al, 2 + ad] = -I * s[al, ad]
    # (rho eta)_alphadot = i v_m eta^alpha sigma^m_{alpha alphadot}
    for ad in range(2):
        for al in range(2):
            M[2 + ad, al] = I * s[al, ad]
    return M


Pm = sp.diag(1, 1, -1, -1)                     # (xi, eta) = xi_alpha eta^alpha - xibar eta bar
e = [sp.Matrix([1 if j == i else 0 for j in range(4)]) for i in range(4)]
vv = sp.Matrix(sp.symbols('v0:4')); ww = sp.Matrix(sp.symbols('w0:4'))
cliff = {
    'rho rho on Sigma = +|v|^2': sp.simplify(rho_Ss_to_S(vv) * rho_S_to_Ss(vv) - dot(vv, vv) * sp.eye(4)) == sp.zeros(4, 4),
    'rho rho on Sigma = -|v|^2': sp.simplify(rho_Ss_to_S(vv) * rho_S_to_Ss(vv) + dot(vv, vv) * sp.eye(4)) == sp.zeros(4, 4),
    'P rho_{Sigma->Sigma*} symmetric': sp.simplify(Pm * rho_S_to_Ss(vv) - (Pm * rho_S_to_Ss(vv)).T) == sp.zeros(4, 4),
    'rho_{Sigma*->Sigma}^T P symmetric': sp.simplify(rho_Ss_to_S(vv).T * Pm - (rho_Ss_to_S(vv).T * Pm).T) == sp.zeros(4, 4),
    'rho(v)rho(w)+rho(w)rho(v) on Sigma* = 2(v,w)': sp.simplify(rho_S_to_Ss(vv) * rho_Ss_to_S(ww) + rho_S_to_Ss(ww) * rho_Ss_to_S(vv) - 2 * dot(vv, ww) * sp.eye(4)) == sp.zeros(4, 4),
}
out['lambda_realisation_checks'] = {k: bool(v) for k, v in cliff.items()}

# structures as bilinear-form matrices B with value u^T B v (u in Sigma: lambda at p1; v in Sigma*: lambda* at p2)
S1 = Pm
S2 = [Pm * rho_S_to_Ss(e[mu]) * rho_Ss_to_S(P4) for mu in range(4)]        # (u, rho(e^mu) rho(p) v)
phi = 4 * L2 - 4
lam_data = {
    'Pi_int': (-R(8, 3) * S1, [(-R(1, 3) - 2 * L3e) * p[mu] * S1 + (L3e - R(1, 2)) * S2[mu] for mu in range(4)]),
    'ext_leg': (6 * phi * S1, [4 * phi * p[mu] * S1 for mu in range(4)]),
    'S_M|_lambda': (-3 * S1, [-2 * p[mu] * S1 for mu in range(4)]),
    'int(lambda*,lambda)': (S1, [sp.zeros(4, 4)] * 4),
    'X_lambda (log eps coefficient, engine A: -(u, rho(p) rho(e^mu) v))': (
        sp.zeros(4, 4), [-Pm * rho_S_to_Ss(P4) * rho_Ss_to_S(e[mu]) for mu in range(4)]),
    'log eps coefficient of Pi_int (= -2 p^mu (u,v) + (u, rho(e^mu) rho(p) v))': (
        sp.zeros(4, 4), [-2 * p[mu] * S1 + S2[mu] for mu in range(4)]),
}


def mat_at(M, q):
    return M.subs({p[i]: q[i] for i in range(4)}, simultaneous=True)


def Q_lambda(alpha, beta):
    """v -> i rho_{Sigma->Sigma*}(p2) w, Koszul sign -1 (lambda odd); canonical at p2* = -p"""
    p2s = [-x for x in p]
    aG = -(alpha * I * rho_S_to_Ss(p2s) - sum((beta[mu] * I * rho_S_to_Ss(e[mu]) for mu in range(4)), sp.zeros(4, 4)))
    bG = [-(beta[mu] * I * rho_S_to_Ss(p2s)) for mu in range(4)]
    return aG.applyfunc(sp.expand), [x.applyfunc(sp.expand) for x in bG]


def gsym(alpha, beta):
    mp_ = [-x for x in p]
    divbT = sum((sp.diff(beta[mu].T, p[mu]) for mu in range(4)), sp.zeros(4, 4))
    a_sw = mat_at(alpha.T, mp_) - mat_at(divbT, mp_)
    b_sw = [mat_at(beta[mu].T, mp_) for mu in range(4)]
    aS = ((alpha - a_sw) / 2).applyfunc(sp.expand)
    bS = [((beta[mu] - b_sw[mu]) / 2).applyfunc(sp.expand) for mu in range(4)]
    return aS, bS


target = (-I * Pm * rho_S_to_Ss(P4)).applyfunc(sp.expand)   # graded-symmetric alpha of int(lambda, dslash lambda)
# check the target's own graded symmetrisation
t_a, t_b = gsym((I * Pm * rho_S_to_Ss([-x for x in p])).applyfunc(sp.expand), [sp.zeros(4, 4)] * 4)
out['lambda_target_check'] = {'gsym(int(lambda,dslash lambda)) alpha == -i P rho(p)': sp.simplify(t_a - target) == sp.zeros(4, 4)}
resL = {}
for name, (al, be) in lam_data.items():
    aG, bG = Q_lambda(al, be)
    aS, bS = gsym(aG, bG)
    # aS = x * target ?
    xs = sp.Symbol('x')
    diffm = (aS - xs * target).applyfunc(lambda z: sp.expand(z.subs(sp.log(pp), sp.Symbol('LOGPP'))))
    eqs = []
    for z in diffm:
        if z != 0:
            num = sp.numer(sp.together(z))
            vars_ = [v for v in [sp.Symbol('LOGPP'), le, gE, L2c, L3c] + list(p) if num.has(v)]
            eqs += sp.Poly(sp.expand(num), *vars_).coeffs() if vars_ else [num]
    sol = sp.solve(eqs, [xs], dict=True) if eqs else [{xs: 0}]
    bzero = all(sp.simplify(z.subs(sp.log(pp), sp.Symbol('LOGPP'))) == 0 for x in bS for z in x)
    resL[name] = {'Q_gsym in units of int(lambda,dslash lambda)': (str(sp.nsimplify(sol[0][xs])) if sol else 'NO FIT'),
                  'beta_gsym_zero': bzero}
out['lambda_(lambda,lambda)_output'] = resL

# ------------------------------------------------------------------------------------------------ assembly
def num(s):
    return sp.nsimplify(s) if s is not None else None


QPi = {'dAdA_perp': num(resP['Pi_int']['Q_sym in units of int(dA_perp,dA_perp)']),
       'dAdA_V': num(resV['Pi_int']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][0]),
       'd*Ad*A_V': num(resV['Pi_int']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][1]),
       'lamdlam': num(resL['Pi_int']['Q_gsym in units of int(lambda,dslash lambda)']),
       'dcAs': num(out['(c,A*)_output_totals']['Q(Pi_int V + c)']['alpha in units of int(dc,A*)'])}
Qii = {'dAdA_perp': num(resP['ext_leg']['Q_sym in units of int(dA_perp,dA_perp)']),
       'dAdA_V': num(resV['ext_leg_1 (S_M A*_V external, A_V leg on W_1^AV)']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][0])
       + num(resV['ext_leg_2 (S_M A_V external, A*_V leg on the ghost self-energy)']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][0])
       + num(resV['F0_graph (S_M - P - (-F_0), F_0 = (27/16) int(d*A_V,d*A_V))']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][0]),
       'd*Ad*A_V': num(resV['ext_leg_1 (S_M A*_V external, A_V leg on W_1^AV)']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][1])
       + num(resV['ext_leg_2 (S_M A_V external, A*_V leg on the ghost self-energy)']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][1])
       + num(resV['F0_graph (S_M - P - (-F_0), F_0 = (27/16) int(d*A_V,d*A_V))']['Q_sym on (A_V,A_V): [int(dA,dA), int(d*A,d*A)]'][1]),
       'lamdlam': num(resL['ext_leg']['Q_gsym in units of int(lambda,dslash lambda)']),
       'dcAs': num(out['(c,A*)_output_totals']['Q(class ii) with F_0']['alpha in units of int(dc,A*)'])}
theory_QPi = {'dAdA_perp': R(1, 3), 'dAdA_V': R(1, 3), 'd*Ad*A_V': 0, 'lamdlam': -R(1, 3), 'dcAs': -R(2, 3)}
theory_4Clog = {'dAdA_perp': -4, 'dAdA_V': -4, 'd*Ad*A_V': 0, 'lamdlam': 16, 'dcAs': -4}
B_Ob = {'dAdA_perp': -R(11, 3), 'dAdA_V': -R(11, 3), 'd*Ad*A_V': 0, 'lamdlam': R(47, 3), 'dcAs': -R(14, 3)}
cmp = {}
for kk in QPi:
    ObA = QPi[kk] + Qii[kk]
    cmp[kk] = {'Q Pi_int (engine A data)': str(QPi[kk]), 'theory Q Pi_Psi': str(theory_QPi[kk]),
               'Q class ii (engine A data, F_0 incl.)': str(Qii[kk]), 'theory 4 C_log': str(theory_4Clog[kk]),
               'Ob from engine A data': str(ObA), 'engine B Ob': str(B_Ob[kk]),
               'agree_theory': (QPi[kk] == theory_QPi[kk]) and (Qii[kk] == theory_4Clog[kk]),
               'agree_engineB': ObA == B_Ob[kk]}
out['assembly_Ob_Delta'] = cmp

# translation-invariant representative T = t_A int(A,A*) + t_l int(lambda*,lambda) + t_c int(c*,c) with Q T = Q Pi_Psi
# Q int(A_perp,A*_perp) = -1/2 dAdA ; Q int(A_V,A*_V) = -1/2 dAdA_V - dcAs ; Q int(lambda*,lambda) = -lamdlam ; Q int(c*,c) = -dcAs
tA, tl, tc = sp.symbols('tA tl tc')
solT = sp.solve([-tA / 2 - QPi['dAdA_V'], -tl - QPi['lamdlam'], -tA - tc - QPi['dcAs']], [tA, tl, tc], dict=True)[0]
out['translation_invariant_representative_of_Pi_Psi'] = {
    'from engine A Pi_int (via Q)': {str(k_): str(v_) for k_, v_ in solT.items()},
    'perp sector t_A from Q': str(-2 * QPi['dAdA_perp']),
    'theory prediction': {'tA': '-2/3', 'tl': '1/3', 'tc': '4/3'},
    'agree': solT[tA] == -R(2, 3) and solT[tl] == R(1, 3) and solT[tc] == R(4, 3) and -2 * QPi['dAdA_perp'] == -R(2, 3)}
# F(Delta) translation-invariant: Q F = Ob ; F = 4 FC + T ; FC = 2 int(A,A*) - 4 int(lambda*,lambda) - int(c*,c)
Fa, Fl, Fc = 8 + solT[tA], -16 + solT[tl], -4 + solT[tc]
out['F_Delta_translation_invariant'] = {
    'from engine A data: 4 FC + T': {'int(A,A*)': str(Fa), 'int(lambda*,lambda)': str(Fl), 'int(c*,c)': str(Fc)},
    'engine B': {'int(A,A*)': '22/3', 'int(lambda*,lambda)': '-47/3', 'int(c*,c)': '-8/3'},
    'agree': (Fa, Fl, Fc) == (R(22, 3), -R(47, 3), -R(8, 3)),
    'as (11/3) FC + rest': {'rest int(lambda*,lambda)': str(Fl - R(11, 3) * (-4)), 'rest int(c*,c)': str(Fc - R(11, 3) * (-1))}}

# locality of the renormalised class (i) (Pi_int with log eps removed)
loc = {}
for nm, (al, be) in [('A_perp', data_P['Pi_int']), ('A_V', data_V['Pi_int']), ('c', data_c['Pi_int'])]:
    alr = sp.expand(al.subs(le, 0)); ber = [sp.expand(x.subs(le, 0)) for x in be]
    nonloc = any(x.has(sp.log(pp)) for x in [alr] + ber) or any(sp.denom(sp.together(x)).has(p[0]) for x in [alr] + ber)
    loc[nm] = {'nonlocal (log|p|^2 or 1/|p|^2 present)': bool(nonloc)}
al, be = lam_data['Pi_int']
loc['lambda'] = {'nonlocal (log|p|^2 or 1/|p|^2 present)': any(z.subs(le, 0).has(sp.log(pp)) for x in be for z in x)}
out['locality_of_renormalised_Pi_int'] = loc

# engine A internal: A_V total alpha (a.b) constant
const_ab = sp.Rational(43, 2) + sp.Rational(13, 6)
out['engineA_AV_total_alpha_ab_constant'] = {'reported': '70/3', 'sum of reported pieces 43/2 + 13/6': str(const_ab)}

json.dump(out, open(os.path.join(HERE, 'two_point_sectors.json'), 'w'), indent=1, default=str)
print(json.dumps({k: out[k] for k in ['lambda_realisation_checks', 'lambda_target_check', 'assembly_Ob_Delta',
                                      'translation_invariant_representative_of_Pi_Psi', 'F_Delta_translation_invariant',
                                      'locality_of_renormalised_Pi_int', 'engineA_AV_total_alpha_ab_constant']}, indent=1, default=str))
print(json.dumps(out['(c,A*)_output_totals'], indent=1, default=str))
