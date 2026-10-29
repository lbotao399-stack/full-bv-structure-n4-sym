# k_tr_relations.py -- exact checks of the algebraic statements of k_transfer.tex, Section "The definition of K_1 and
# its consistency relations".
#
# (1) Brackets of vector fields (ch02 Def c2-bracket, [v,v'] = v(v') - v'(v)):
#     [K_b, P_c] = (b,c) 2E + w_bc,  w_bc(x) = 2(b,x)c - 2(c,x)b;  [Delta, K_b] = 2 K_b;  [Delta, P_c] = -2 P_c;
#     [K_b, K_c] = 0;  (1/4) div is a 1-cocycle: (1/4)div[v,v'] = v((1/4)div v') - v'((1/4)div v).
# (2) Operator model.  On the local functionals int omega E[f] (E a linear differential expression with polynomial
#     coefficients in a scalar density f of Weyl weight w), every operator considered acts by substitution
#     f -> L f with L a first-order differential operator: the classical delta_v = {S^(1)(v), -} by
#     L_v = v.d + (w/4) div v (ch07 Def c7-vdot, Lemma c7-weights), and the one-loop maps of the ansatz by
#     L = alpha v.d + (kappa/4) div v (ch07 Cor c7-K; the density of K_1(int omega E f) is E applied to that of
#     K_1(int omega f), by integration by parts in omega).  Composition reverses the order: A1 A2 <-> L2 L1, so
#     [A1, A2] <-> -[L1, L2].  Check: [delta_v1, delta_v2] = -delta_[v1,v2]  (Cor c8-homotopy (ii) at order hbar^0).
# (3) The order-hbar components of L o L = 0 (Cor c8-homotopy (ii)) at (K_b,P_c), (P_c,Delta), (K_b,Delta), (K_b,K_c)
#     with the ansatz  P_1 = p delta_P,  D_1 = 2 gamma + d delta_Delta,  M_1 = m delta_rot,
#     K_1 = alpha (v.d) + (kappa/4) div v  (scalar densities; no second-derivative term, ch07 Cor c7-K):
#     solve for alpha, kappa, d, m.  Expected: d = 0, m = 0, alpha = -p, kappa = gamma - p w,
#     i.e. K_1 = (1/4) div v gamma - p delta_K.
# (4) The scheme change psi' = (a/2) delta_Delta (finite part of the operators, Prop c8-laws (b)) shifts
#     L_1 by [s, psi']: P_1 by -a delta_P, K_1 by +a delta_K, D_1 and M_1 by 0.
# (5) Two-component (Euclidean sigma of Def c8-twocomp): the bracket of two special supercharges, (B1) of the
#     dictionary with eps_a(x) = rho(x) eta_a (K1), equals -K_{v1} with v1^m = -i(eta2 sigma^m etabar1 + eta1 sigma^m etabar2)
#     = Gamma(eta1,eta2)|_V^m (Def c8-twocomp); the map (eta1,eta2) -> v1 is onto V.
# (6) Scale invariance of the composite kernel F_k(z_a, z_b; eps, T) (k_tr_kernel.py): T dF/dT at (lambda z, lambda^2 eps,
#     lambda^2 T) equals that at (z, eps, T).
import json
import itertools
import numpy as np
import sympy as sp

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/k_judgement/rerun_theory/k_tr_relations.json'
res = {}
x = sp.symbols('x1:5', real=True)
bS = sp.symbols('b1:5', real=True)
cS = sp.symbols('c1:5', real=True)
dot = lambda a, c: sum(a[i]*c[i] for i in range(4))
xx = dot(x, x)
K = [xx*bS[n] - 2*dot(bS, x)*x[n] for n in range(4)]
Kc = [xx*cS[n] - 2*dot(cS, x)*x[n] for n in range(4)]
P = list(cS)
Dl = [2*x[n] for n in range(4)]
wbc = [2*dot(bS, x)*cS[n] - 2*dot(cS, x)*bS[n] for n in range(4)]


def vbr(v, w):
    return [sp.expand(sum(v[m]*sp.diff(w[n], x[m]) - w[m]*sp.diff(v[n], x[m]) for m in range(4))) for n in range(4)]


def div(v):
    return sp.expand(sum(sp.diff(v[n], x[n]) for n in range(4)))


r1 = {}
r1['[K_b,P_c] - ((b,c) 2E + w_bc)'] = [str(sp.simplify(a - (dot(bS, cS)*Dl[n] + wbc[n]))) for n, a in enumerate(vbr(K, P))]
r1['[Delta,K_b] - 2K_b'] = [str(sp.simplify(a - 2*K[n])) for n, a in enumerate(vbr(Dl, K))]
r1['[Delta,P_c] + 2P_c'] = [str(sp.simplify(a + 2*P[n])) for n, a in enumerate(vbr(Dl, P))]
r1['[K_b,K_c]'] = [str(a) for a in vbr(K, Kc)]
r1['div K_b + 8(b,x)'] = str(sp.simplify(div(K) + 8*dot(bS, x)))
r1['div w_bc'] = str(div(wbc))
cocyc = {}
for (n1, v1), (n2, v2) in itertools.combinations([('K_b', K), ('P_c', P), ('Delta', Dl), ('w_bc', wbc), ('K_c', Kc)], 2):
    lhs = sp.Rational(1, 4)*div(vbr(v1, v2))
    rhs = sum(v1[m]*sp.diff(sp.Rational(1, 4)*div(v2), x[m]) - v2[m]*sp.diff(sp.Rational(1, 4)*div(v1), x[m]) for m in range(4))
    cocyc[n1 + ',' + n2] = str(sp.simplify(lhs - rhs))
r1['(1/4)div cocycle residuals'] = cocyc
res['1_brackets'] = r1

# ------------------------------------------------------------------ (2) operator model
f = sp.Function('f')(*x)
w = sp.symbols('w')


def Lop(alpha, v, kappa):
    """L g = alpha v.d g + (kappa/4) div v g"""
    dv = div(v)
    return lambda g: alpha*sum(v[m]*sp.diff(g, x[m]) for m in range(4)) + kappa*sp.Rational(1, 4)*dv*g


def Lcl(v):
    return Lop(1, v, w)


def comm_ops(L1, L2):
    """[A1, A2] in the model: f -> -(L1 L2 - L2 L1) f"""
    return sp.expand(-(L1(L2(f)) - L2(L1(f))))


r2 = {}
for (n1, v1), (n2, v2) in itertools.combinations([('K_b', K), ('P_c', P), ('Delta', Dl), ('w_bc', wbc)], 2):
    lhs = comm_ops(Lcl(v1), Lcl(v2))
    rhs = sp.expand(-Lcl(vbr(v1, v2))(f))
    r2['[delta_%s, delta_%s] + delta_[%s,%s]' % (n1, n2, n1, n2)] = str(sp.simplify(lhs - rhs))
res['2_classical_commutators_residuals'] = r2

# ------------------------------------------------------------------ (3) order-hbar relations with the rescaling ansatz
al, ka, d, m, p, ga = sp.symbols('alpha kappa d m p gamma')
L_K1 = Lop(al, K, ka)
L_P1 = Lop(p, P, 0)
L_D1 = lambda g: 2*ga*g + d*Lcl(Dl)(g)
L_M1 = Lop(m, wbc, 0)
L_dP, L_dK, L_dD = Lcl(P), Lcl(K), Lcl(Dl)


def as_op(L):
    return sp.expand(L(f))


# (K_b,P_c): [K_1, d_P] + [d_K, P_1] + (b,c) D_1 + M_1(w_bc) = 0 on cohomology
eqKP = sp.expand(comm_ops(L_K1, L_dP) + comm_ops(L_dK, L_P1) + dot(bS, cS)*as_op(L_D1) + as_op(L_M1))
# (P_c,Delta): [P_1, d_D] + [d_P, D_1] + P_1([P,Delta]) = 0,  [P_c,Delta] = 2 P_c
eqPD = sp.expand(comm_ops(L_P1, L_dD) + comm_ops(L_dP, L_D1) + 2*as_op(L_P1))
# (K_b,Delta): [K_1, d_D] + [d_K, D_1] + K_1([K_b,Delta]) = 0,  [K_b,Delta] = -2 K_b
eqKD = sp.expand(comm_ops(L_K1, L_dD) + comm_ops(L_dK, L_D1) - 2*as_op(L_K1))


def coeff_eqs(expr):
    """coefficients of f, d_i f, d_i d_j f and of the monomials in x, b, c."""
    derivs = [f] + [sp.diff(f, x[i]) for i in range(4)] + [sp.diff(f, x[i], x[j]) for i in range(4) for j in range(i, 4)]
    reps = {dv: sp.Symbol('F%d' % k) for k, dv in enumerate(derivs)}
    e = expr
    for dv in reversed(derivs):
        e = e.subs(dv, reps[dv])
    poly = sp.Poly(sp.expand(e), *(list(reps.values()) + list(x) + list(bS) + list(cS)))
    return [c_ for c_ in poly.coeffs()]


eqs = coeff_eqs(eqKP) + coeff_eqs(eqPD) + coeff_eqs(eqKD)
sol = sp.solve(eqs, [al, ka, d, m], dict=True)
res['3_rescaling_relations'] = {'unknowns': 'K_1 = alpha (v.d) + (kappa/4) div v;  D_1 = 2 gamma + d delta_Delta;  M_1 = m delta_rot;  P_1 = p delta_P',
                                'solution': [{str(k_): str(v_) for k_, v_ in s.items()} for s in sol],
                                'expected': {'alpha': '-p', 'kappa': 'gamma - p*w', 'd': '0', 'm': '0'}}
# with d left free, the (K,P) relation alone: alpha = d - p, kappa = gamma + (d - p) w, m = d
solKP = sp.solve(coeff_eqs(eqKP), [al, ka, m], dict=True)
res['3_KP_alone'] = [{str(k_): str(sp.simplify(v_)) for k_, v_ in s.items()} for s in solKP]
solPD = sp.solve(coeff_eqs(eqPD), [d], dict=True)
res['3_PDelta_alone'] = [{str(k_): str(v_) for k_, v_ in s.items()} for s in solPD]
# (K_b,K_c) with the solution: [K_1(K_b), d_Kc] + [d_Kb, K_1(K_c)] = 0
if sol:
    s0 = sol[0]
    L_K1b = Lop(s0[al], K, s0[ka]); L_K1c = Lop(s0[al], Kc, s0[ka])
    eqKK = sp.expand(comm_ops(L_K1b, Lcl(Kc)) + comm_ops(Lcl(K), L_K1c))
    res['3_KK_residual_with_solution'] = str(sp.simplify(eqKK))

# ------------------------------------------------------------------ (4) scheme change psi' = (a/2) delta_Delta
a = sp.symbols('a')
L_psi = lambda g: sp.Rational(1, 2)*a*Lcl(Dl)(g)
shift = {}
for name, v in [('P_c', P), ('K_b', K), ('Delta', Dl), ('w_bc', wbc)]:
    sh = comm_ops(Lcl(v), L_psi)                      # [delta_v, psi'] = arity-one part of [s, psi'] on these operators
    cl = as_op(Lcl(v))
    rat = sp.simplify(sh/cl) if sp.simplify(cl) != 0 else None
    shift[name] = {'[delta_v, psi\']': str(sp.factor(sh)), 'ratio_to_delta_v': str(rat)}
res['4_scheme_change_(a/2)delta_Delta'] = shift

# ------------------------------------------------------------------ (5) two-component bracket of two special supercharges
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]], complex); s3 = np.array([[1, 0], [0, -1]], complex)
sig = [s1, s2, s3, 1j*np.eye(2)]
sigb = [-s1, -s2, -s3, 1j*np.eye(2)]
rng = np.random.default_rng(1028)


def killing(eta, etab, X):
    """(K1): eps_al^i(x) = -i x_m sig^m_{al dal} etabar^{dal i},  epsbar_{dal i}(x) = i x_m eta^al_i sig^m_{al dal}."""
    ep = -1j*sum(X[mm]*sig[mm] for mm in range(4)) @ etab           # [al, i]
    epb = 1j*np.einsum('ai,ad->di', eta, sum(X[mm]*sig[mm] for mm in range(4)))   # [dal, i]
    return ep, epb


def B1_Vpart(e1, eb1, e2, eb2):
    """V part of (B1): i (epsbar1 sigbar^m eps2 + epsbar2 sigbar^m eps1)"""
    return np.array([1j*(np.einsum('di,da,ai->', eb1, sigb[mm], e2) + np.einsum('di,da,ai->', eb2, sigb[mm], e1)) for mm in range(4)])


eta1 = rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4)); etab1 = rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4))
eta2 = rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4)); etab2 = rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4))
v1 = np.array([-1j*(np.einsum('ai,ad,di->', eta2, sig[mm], etab1) + np.einsum('ai,ad,di->', eta1, sig[mm], etab2)) for mm in range(4)])
devs = []
for _ in range(5):
    X = rng.normal(size=4)
    e1, eb1 = killing(eta1, etab1, X); e2, eb2 = killing(eta2, etab2, X)
    br = B1_Vpart(e1, eb1, e2, eb2)
    Kv = (X @ X)*v1 - 2*(v1 @ X)*X
    devs.append(float(np.abs(br + Kv).max()))
# rank of (eta1, eta2) -> v1
M = []
for idx in range(16):
    E = np.zeros(16, complex); E[idx] = 1
    e2t, eb2t = E[:8].reshape(2, 4), E[8:].reshape(2, 4)
    M.append([-1j*(np.einsum('ai,ad,di->', e2t, sig[mm], etab1) + np.einsum('ai,ad,di->', eta1, sig[mm], eb2t)) for mm in range(4)])
rank = int(np.linalg.matrix_rank(np.array(M)))
res['5_two_component_eta_eta_bracket'] = {'max|[eta1,eta2]|_V + K_{v1}| over 5 random x': max(devs),
                                          'rank of eta2 -> Gamma(eta1,eta2)|_V (eta1 generic)': rank}

# ------------------------------------------------------------------ (6) scale invariance of the composite kernel
gh_x, gh_w = np.polynomial.hermite_e.hermegauss(3)
GX = np.array(np.meshgrid(gh_x, gh_x, gh_x, gh_x, indexing='ij')).reshape(4, -1).T
GW = np.prod(np.array(np.meshgrid(gh_w, gh_w, gh_w, gh_w, indexing='ij')).reshape(4, -1), 0)/(2*np.pi)**2
bn = np.array([0.7, -1.1, 0.4, 0.9])


def pair_integral(k, tfv, tgv, zan, zbn):
    alpha = 1/(4*tfv) + 1/(4*tgv)
    y0 = (zan/tfv + zbn/tgv)/(1/tfv + 1/tgv)
    beta = ((zan - zbn)**2).sum()/(4*(tfv + tgv))
    s = np.sqrt(1/(2*alpha))
    Yp = y0[None, :] + s*GX
    yy_ = (Yp**2).sum(1); by_ = Yp @ bn
    vv = yy_[:, None]*bn[None, :] - 2*by_[:, None]*Yp
    df = -(Yp - zan[None, :])/(2*tfv); dg = -(Yp - zbn[None, :])/(2*tgv)
    vdf = (vv*df).sum(1)
    dvk = 2*Yp*bn[k] - 2*bn[None, :]*Yp[:, k:k+1] - 2*by_[:, None]*(np.arange(4) == k)[None, :]
    poly = vdf*dg[:, k] + (dvk*dg).sum(1)
    pref = (4*np.pi*tfv)**-2*(4*np.pi*tgv)**-2*(np.pi/alpha)**2*np.exp(-beta)
    return pref*(poly*GW).sum()


from scipy import integrate as si  # noqa: E402


def TdF(k, eps, T, za, zb):
    g = lambda u: np.exp(u)*T*(pair_integral(k, T, np.exp(u), za, zb) + pair_integral(k, np.exp(u), T, za, zb))
    return si.quad(g, np.log(eps), np.log(T), limit=400, epsabs=0, epsrel=1e-11)[0]


za = np.array([0.3, -0.2, 0.5, 0.1]); zb = np.array([-0.4, 0.25, 0.15, -0.3])
sc = {}
for lam in [2.0, 0.5]:
    a0 = TdF(0, 1.0, 50.0, za, zb); a1 = TdF(0, lam**2, lam**2*50.0, lam*za, lam*zb)
    sc['lambda=%g' % lam] = {'TdF(z,1,50)': a0, 'TdF(lambda z, lambda^2, 50 lambda^2)': a1, 'rel_dev': abs(a1 - a0)/abs(a0)}
res['6_composite_kernel_scale_invariance'] = sc
json.dump(res, open(OUT, 'w'), indent=1)
print(json.dumps(res, indent=1))
