"""Audit script a6: two-component realisation of ch08 Def c8-twocomp / Prop c8-twocomp / Lemma c8-cc and the
two-component formulas of the revised odd design that rely on it.  Exact arithmetic with sympy matrices.

Checks: Euclidean sigma relations; Fierz identity used in design prop:o-konishi; 't Hooft identity;
Clifford relation of rho on V (+) W_perp from the component formulas; symmetry of (rho(v)z1,z2);
Gamma components (eq. for Gamma|_V, Gamma_kl, Gamma^kl) against (rho(v)z1,z2) = (v,Gamma(z1,z2));
rank-one criterion and the alternative form of Gamma_kl(xi,xi) in design lem:o-rankone; Killing spinor
eq:c8-killing = xi + rho(x) eta; rho(e1e2e3e4) = +1 on S+ (x) W; the R-rotation r = 4/3(w n - 1/4) of
design rem:o-alpha with (B3) of the dictionary; conjugation of Lemma c8-cc (rho-compatibility,
positivity, Gamma compatibility); B and ad_Delta versus the radial and reflection conjugations.
"""
import json, os, itertools, random
import sympy as sp

I = sp.I
s1 = sp.Matrix([[0, 1], [1, 0]]); s2 = sp.Matrix([[0, -I], [I, 0]]); s3 = sp.Matrix([[1, 0], [0, -1]])
one = sp.eye(2)
sig = [s1, s2, s3, I * one]            # sigma^m_{alpha alphadot}
sigb = [-s1, -s2, -s3, I * one]        # sigmabar^{m alphadot alpha}
eps_up = sp.Matrix([[0, 1], [-1, 0]])  # eps^{12} = +1
eps_dn = sp.Matrix([[0, -1], [1, 0]])  # eps_{12} = -1
res = {}

# 1. sigma relations
ok = all((sig[m] * sigb[n] + sig[n] * sigb[m] - (-2 if m == n else 0) * one) == sp.zeros(2) for m in range(4) for n in range(4))
ok2 = all((sigb[m] * sig[n] + sigb[n] * sig[m] - (-2 if m == n else 0) * one) == sp.zeros(2) for m in range(4) for n in range(4))
res['sigma sigmabar + sigma sigmabar = -2 delta'] = ok
res['sigmabar sigma + sigmabar sigma = -2 delta'] = ok2
# 2. Fierz sum_m sigbar^{m ad a} sigbar^{m bd b} = -2 eps^{ab} eps^{ad bd}
fz = True
for ad, a, bd, b in itertools.product(range(2), repeat=4):
    lhs = sum(sigb[m][ad, a] * sigb[m][bd, b] for m in range(4))
    rhs = -2 * eps_up[a, b] * eps_up[ad, bd]
    if sp.simplify(lhs - rhs) != 0:
        fz = False
res['Fierz sum_m sigbar^m sigbar^m = -2 eps eps (design prop:o-konishi)'] = fz

# 3. 't Hooft blocks
def levi(*idx):
    return sp.LeviCivita(*idx)
Sig = []
for a in range(3):
    M = sp.zeros(4)
    for b in range(3):
        for c in range(3):
            M[b, c] = levi(a, b, c)
        M[b, 3] = 1 if a == b else 0
        M[3, b] = -1 if a == b else 0
    Sig.append(M)
for a in range(3):
    M = sp.zeros(4)
    for b in range(3):
        for c in range(3):
            M[b, c] = I * levi(a, b, c)
        M[b, 3] = -I if a == b else 0
        M[3, b] = I if a == b else 0
    Sig.append(M)
th = all(sp.simplify(sum(Sig[J][i, j] * Sig[J][k, l] for J in range(6)) - 2 * levi(i, j, k, l)) == 0 for i, j, k, l in itertools.product(range(4), repeat=4))
res["'t Hooft sum_I Sigma^I_ij Sigma^I_kl = 2 eps_ijkl"] = th

# spinor spaces: xi = (xi_{alpha}^i [2x4], xib_{alphadot i} [2x4]); eta = (eta^alpha_i [2x4], etab^{alphadot i} [2x4])
def lower_pair(wI):
    wl = sp.zeros(4)
    for J in range(6):
        wl += wI[J] * Sig[J]
    wu = sp.zeros(4)
    for i, j in itertools.product(range(4), repeat=2):
        wu[i, j] = sp.Rational(1, 2) * sum(levi(i, j, k, l) * wl[k, l] for k in range(4) for l in range(4))
    return wl, wu


def rho_xi(v, wI, xi):
    """rho(v + w) xi in Sigma^*: returns (eta^alpha_i, etab^{alphadot i})"""
    x, xb = xi
    eta = sp.zeros(2, 4); etab = sp.zeros(2, 4)
    wl, wu = lower_pair(wI)
    for al, i in itertools.product(range(2), range(4)):
        # (rho(v)xi)^alpha_i = i v_m xib_{ad i} sigbar^{m ad alpha}
        eta[al, i] += I * sum(v[m] * xb[ad, i] * sigb[m][ad, al] for m in range(4) for ad in range(2))
        # (rho(w)xi)^alpha_j = eps^{alpha beta} w_{ji} xi_beta^i
        eta[al, i] += sum(eps_up[al, be] * wl[i, k] * x[be, k] for be in range(2) for k in range(4))
    for ad, i in itertools.product(range(2), range(4)):
        # (rho(v)xi)^{ad i} = -i v_m sigbar^{m ad alpha} xi_alpha^i
        etab[ad, i] += -I * sum(v[m] * sigb[m][ad, al] * x[al, i] for m in range(4) for al in range(2))
        # (rho(w)xi)^{ad j} = eps^{ad bd} w^{ji} xib_{bd i}
        etab[ad, i] += sum(eps_up[ad, bd] * wu[i, k] * xb[bd, k] for bd in range(2) for k in range(4))
    return (eta, etab)


def rho_eta(v, wI, et):
    """rho(v + w) eta in Sigma: returns (xi_alpha^i, xib_{alphadot i})"""
    e, eb = et
    x = sp.zeros(2, 4); xb = sp.zeros(2, 4)
    wl, wu = lower_pair(wI)
    for al, i in itertools.product(range(2), range(4)):
        # (rho(v)eta)_alpha^i = -i v_m sigma^m_{alpha ad} etab^{ad i}
        x[al, i] += -I * sum(v[m] * sig[m][al, ad] * eb[ad, i] for m in range(4) for ad in range(2))
        # (rho(w)eta)_alpha^k = -w^{kj} eps_{alpha beta} eta^beta_j
        x[al, i] += -sum(wu[i, j] * eps_dn[al, be] * e[be, j] for j in range(4) for be in range(2))
    for ad, i in itertools.product(range(2), range(4)):
        # (rho(v)eta)_{ad i} = i v_m eta^alpha_i sigma^m_{alpha ad}
        xb[ad, i] += I * sum(v[m] * e[al, i] * sig[m][al, ad] for m in range(4) for al in range(2))
        # (rho(w)eta)_{bd k} = -w_{kj} eps_{bd ad} etab^{ad j}
        xb[ad, i] += -sum(wl[i, j] * eps_dn[ad, bd] * eb[bd, j] for j in range(4) for bd in range(2))
    return (x, xb)


def pair(xi, et):
    x, xb = xi; e, eb = et
    return sum(x[al, i] * e[al, i] for al in range(2) for i in range(4)) - sum(xb[ad, i] * eb[ad, i] for ad in range(2) for i in range(4))


rng = random.Random(3)
def rint():
    return sp.Integer(rng.randint(-3, 3)) + I * sp.Integer(rng.randint(-2, 2))
def rxi():
    return (sp.Matrix(2, 4, lambda a, b: rint()), sp.Matrix(2, 4, lambda a, b: rint()))
def rvec(n):
    return [sp.Integer(rng.randint(-3, 3)) for _ in range(n)]

# 4. Clifford: rho(u)rho(u) = (u,u) with (u,u) = sum v_m^2 + sum w_I^2
cl = True
for _ in range(3):
    v = rvec(4); wI = rvec(6); xi = rxi()
    back = rho_eta(v, wI, rho_xi(v, wI, xi))
    nrm = sum(t * t for t in v) + sum(t * t for t in wI)
    if any(sp.expand(back[k] - nrm * xi[k]) != sp.zeros(2, 4) for k in range(2)):
        cl = False
res['Clifford rho(v+w)^2 = |v|^2 + |w|^2 on Sigma'] = cl
# 5. symmetry (rho(u) z1, z2) = (z1, rho(u) z2) for z1 in Sigma, z2 in Sigma (both sides in Sigma^* pairing)
sy = True
for _ in range(3):
    v = rvec(4); wI = rvec(6); z1 = rxi(); z2 = rxi()
    a = pair(z2, rho_xi(v, wI, z1)); b = pair(z1, rho_xi(v, wI, z2))
    if sp.expand(a - b) != 0:
        sy = False
res['(rho(u)z1,z2) = (z1,rho(u)z2) on Sigma x Sigma'] = sy

# 6. Gamma formulas vs (rho(v) z1, z2) = (v, Gamma(z1,z2)), (Gamma,u) = Gamma^m u_m + 1/4 Gamma_ij u^ij
def Gamma_xixi(z1, z2):
    x1, xb1 = z1; x2, xb2 = z2
    GV = [I * sum(xb2[ad, i] * sigb[m][ad, al] * x1[al, i] + xb1[ad, i] * sigb[m][ad, al] * x2[al, i] for ad in range(2) for al in range(2) for i in range(4)) for m in range(4)]
    Gl = sp.zeros(4)
    for k, l in itertools.product(range(4), repeat=2):
        Gl[k, l] = 2 * sum(levi(j, i, k, l) * x1[al, j] * eps_up[al, be] * x2[be, i] for j in range(4) for i in range(4) for al in range(2) for be in range(2)) \
            - 2 * sum(xb1[ad, k] * eps_up[ad, bd] * xb2[bd, l] - xb1[ad, l] * eps_up[ad, bd] * xb2[bd, k] for ad in range(2) for bd in range(2))
    return GV, Gl
gm = True
for _ in range(3):
    v = rvec(4); wI = rvec(6); z1 = rxi(); z2 = rxi()
    lhs = pair(z2, rho_xi(v, wI, z1))
    GV, Gl = Gamma_xixi(z1, z2)
    wl, wu = lower_pair(wI)
    rhs = sum(GV[m] * v[m] for m in range(4)) + sp.Rational(1, 4) * sum(Gl[i, j] * wu[i, j] for i in range(4) for j in range(4))
    if sp.expand(lhs - rhs) != 0:
        gm = False
res['Gamma(xi1,xi2) components reproduce (rho(u)xi1,xi2) = (u,Gamma)'] = gm

# 7. rank one: Gamma(u (x) w, u (x) w) = 0 in all components; rank two nonzero
u = sp.Matrix([2, 1 + I]); w = sp.Matrix([1, 0, 2, -I])
xi1 = (u * w.T, sp.zeros(2, 4))
GV, Gl = Gamma_xixi(xi1, xi1)
res['rank-one xi: Gamma(xi,xi) = 0 (all ten)'] = all(sp.expand(t) == 0 for t in GV) and sp.expand(Gl) == sp.zeros(4)
u2 = sp.Matrix([1, -1]); w2 = sp.Matrix([0, 1, 1, 0])
xi2 = (u * w.T + u2 * w2.T, sp.zeros(2, 4))
GV2, Gl2 = Gamma_xixi(xi2, xi2)
res['rank-two xi: Gamma_kl(xi,xi) != 0'] = sp.expand(Gl2) != sp.zeros(4)
# alternative form in design lem:o-rankone: Gamma_kl(xi,xi) = -2 eps_jikl xi^{alpha j} xi_alpha^i
alt = sp.zeros(4)
for k, l in itertools.product(range(4), repeat=2):
    alt[k, l] = -2 * sum(levi(j, i, k, l) * eps_up[al, be] * xi2[0][be, j] * xi2[0][al, i] for j in range(4) for i in range(4) for al in range(2) for be in range(2))
res['design lem:o-rankone form of Gamma_kl equals ch08 form'] = sp.simplify(alt - Gl2) == sp.zeros(4)

# 8. Killing spinor eq:c8-killing equals xi + rho(x) eta
ks = True
for _ in range(2):
    x = rvec(4); xi = rxi(); et = rxi()
    rx = rho_eta(x, [0] * 6, et)
    eps_c = (xi[0] - I * sp.Matrix(2, 4, lambda al, i: sum(x[m] * sig[m][al, ad] * et[1][ad, i] for m in range(4) for ad in range(2))),
             xi[1] + I * sp.Matrix(2, 4, lambda ad, i: sum(x[m] * et[0][al, i] * sig[m][al, ad] for m in range(4) for al in range(2))))
    if any(sp.expand(xi[k] + rx[k] - eps_c[k]) != sp.zeros(2, 4) for k in range(2)):
        ks = False
res['eq:c8-killing = xi + rho(x) eta'] = ks

# 9. rho(e1 e2 e3 e4) on S+ (x) W  (product of the four V-directions, Sigma -> Sigma* -> Sigma ...)
E = [[1 if j == m else 0 for j in range(4)] for m in range(4)]
xi = (sp.Matrix(2, 4, lambda a, b: rint()), sp.zeros(2, 4))
t = rho_xi(E[3], [0] * 6, xi); t = rho_eta(E[2], [0] * 6, t); t = rho_xi(E[1], [0] * 6, t); t = rho_eta(E[0], [0] * 6, t)
res['rho(e1e2e3e4) = +1 on S+ (x) W'] = sp.expand(t[0] - xi[0]) == sp.zeros(2, 4) and t[1] == sp.zeros(2, 4)

# 10. R-rotation r = 4/3 (w n - 1/4) with n.w = 1: r w = w, tr r = 0
wv = sp.Matrix([1, 2, 0, -1]); nv = sp.Matrix([1, 0, 0, 0])
r = sp.Rational(4, 3) * (wv * nv.T - sp.Rational(1, 4) * sp.eye(4))
res['r^i_j = 4/3(w^i n_j - delta/4): r w = w and traceless'] = (r * wv - wv == sp.zeros(4, 1)) and r.trace() == 0

# 11. conjugation Lemma c8-cc: bar xi = (eta = conj(xi), etab = -conj(xib))
def bar_xi(xi):
    return (xi[0].conjugate(), -xi[1].conjugate())
def bar_eta(et):
    # inverse: xi = conj(eta), xib = -conj(etab)  (involutive)
    return (et[0].conjugate(), -et[1].conjugate())
cc = True
for _ in range(3):
    v = rvec(4); wI = rvec(6); xi = rxi()
    lhs = bar_eta(rho_xi(v, wI, xi))       # conj(rho(u) xi), an element of Sigma
    rhs = rho_eta(v, wI, bar_xi(xi))        # rho(u) bar xi (u real)
    if any(sp.expand(lhs[k] - rhs[k]) != sp.zeros(2, 4) for k in range(2)):
        cc = False
res['conj(rho(u)xi) = rho(u) conj(xi) for real u'] = cc
xi = rxi()
res['(xi, bar xi) > 0'] = bool(sp.re(sp.expand(pair(xi, bar_xi(xi)))) > 0 and sp.im(sp.expand(pair(xi, bar_xi(xi)))) == 0)

# 12. B and ad_Delta versus conjugations:
# radial: X -> X^dagger maps Sigma <-> Sigma^* (S+W <-> S+W*, S-W* <-> S-W)
Bv = {'S+W': 1, 'S-W*': -1, 'S-W': 1, 'S+W*': -1}
eD = {'S+W': -1, 'S-W*': -1, 'S-W': 1, 'S+W*': 1}
rad = {'S+W': 'S+W*', 'S-W*': 'S-W', 'S+W*': 'S+W', 'S-W': 'S-W*'}
# reflection dagger_n: xi -> rho(n) bar(xi): S+W -> (bar) S+W* -> (rho(n)) S-W*; S-W* -> S-W -> S+W
refl = {'S+W': 'S-W*', 'S-W*': 'S+W', 'S+W*': 'S-W', 'S-W': 'S+W*'}
res['B anticommutes with radial dagger'] = all(Bv[rad[k]] == -Bv[k] for k in Bv)
res['B anticommutes with reflection dagger_n'] = all(Bv[refl[k]] == -Bv[k] for k in Bv)
res['ad_Delta anticommutes with radial dagger'] = all(eD[rad[k]] == -eD[k] for k in eD)
res['ad_Delta commutes with reflection dagger_n'] = all(eD[refl[k]] == eD[k] for k in eD)

out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'a6_twocomp.json')
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res, indent=1))
