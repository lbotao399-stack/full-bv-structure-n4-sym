"""Independent numerical check of the two-component formulas of ch09 (Remarks rem:c9-twocomp3, rem:c9-twocompD)
against Definition def:c8-twocomp (eq:c8-pairing and the rho(w) formulas), plus the symmetry facts used in
Proposition prop:c9-basis and the colour identity of Lemma lem:c9-legs.

All components are treated as commuting complex numbers; the Grassmann ordering of the field functionals is
handled analytically in the findings (one sign per swap of two odd components).
numpy only, < 50 MB.
"""
import json, itertools
import numpy as np

rng = np.random.default_rng(20261029)
I2 = np.eye(2)
s1 = np.array([[0, 1], [1, 0]], complex); s2 = np.array([[0, -1j], [1j, 0]]); s3 = np.array([[1, 0], [0, -1]], complex)
sig = [s1, s2, s3, 1j * I2]            # sigma^m_{alpha alphadot}
sigb = [-s1, -s2, -s3, 1j * I2]        # sigmabar^{m alphadot alpha}  (index order [alphadot, alpha])
eps_up = np.array([[0, 1], [-1, 0]], complex)   # eps^{12}=+1
eps_dn = -eps_up                                  # eps_{12}=-1

# 't Hooft blocks Sigma^I_ij
def levi3(a, b, c):
    return np.linalg.det(np.eye(3)[[a, b, c]])
Sig = np.zeros((6, 4, 4), complex)
for a in range(3):
    for b in range(3):
        for c in range(3):
            Sig[a, b, c] = levi3(a, b, c)
            Sig[a + 3, b, c] = 1j * levi3(a, b, c)
    for b in range(3):
        Sig[a, b, 3] = (a == b); Sig[a, 3, b] = -(a == b)
        Sig[a + 3, b, 3] = -1j * (a == b); Sig[a + 3, 3, b] = 1j * (a == b)
eps4 = np.zeros((4, 4, 4, 4))
for p in itertools.permutations(range(4)):
    eps4[p] = np.linalg.det(np.eye(4)[list(p)])
res = {}
res['tHooft_sum_SigSig_eq_2eps'] = float(np.abs(np.einsum('Iij,Ikl->ijkl', Sig, Sig) - 2 * eps4).max())

# Sigma element: xi[alpha,i] (xi_alpha^i), xib[alphadot,i] (xibar_{alphadot i})
# Sigma* element: eta[alpha,i] (eta^alpha_i), etab[alphadot,i] (etabar^{alphadot i})
def pair(xi, xib, eta, etab):            # (xi,eta) = xi_a^i eta^a_i - xib_{ad i} etab^{ad i}
    return np.einsum('ai,ai->', xi, eta) - np.einsum('ai,ai->', xib, etab)

def rho_v_on_Sigma(v, xi, xib):          # -> Sigma*: (rho xi)^a_i = i v_m xib_{ad i} sigb^{m ad a};  (rho xi)^{ad i} = -i v_m sigb^{m ad a} xi_a^i
    eta = 1j * sum(v[m] * np.einsum('di,da->ai', xib, sigb[m]) for m in range(4))
    etab = -1j * sum(v[m] * np.einsum('da,ai->di', sigb[m], xi) for m in range(4))
    return eta, etab

def rho_v_on_Sigmastar(v, eta, etab):    # -> Sigma: (rho eta)_a^i = -i v_m sig^m_{a ad} etab^{ad i}; (rho eta)_{ad i} = i v_m eta^a_i sig^m_{a ad}
    xi = -1j * sum(v[m] * np.einsum('ad,di->ai', sig[m], etab) for m in range(4))
    xib = 1j * sum(v[m] * np.einsum('ai,ad->di', eta, sig[m]) for m in range(4))
    return xi, xib

def wlow(w):  # w_ij = w_I Sigma^I_ij
    return np.einsum('I,Iij->ij', w, Sig)
def wup(w):   # w^ij = 1/2 eps^{ijkl} w_kl
    return 0.5 * np.einsum('ijkl,kl->ij', eps4, wlow(w))

def rho_w_on_Sigma(w, xi, xib):  # (rho w xi)^a_j = eps^{ab} w_ji xi_b^i ; (rho w xi)^{ad j} = eps^{ad bd} w^{ji} xib_{bd i}
    eta = np.einsum('ab,ji,bi->aj', eps_up, wlow(w), xi)
    etab = np.einsum('ab,ji,bi->aj', eps_up, wup(w), xib)
    return eta, etab

def rho_w_on_Sigmastar(w, eta, etab):  # (rho w eta)_a^k = -w^{kj} eps_{ab} eta^b_j ; (rho w eta)_{bd k} = -w_{kj} eps_{bd ad} etab^{ad j}
    xi = -np.einsum('kj,ab,bj->ak', wup(w), eps_dn, eta)
    xib = -np.einsum('kj,ba,aj->bk', wlow(w), eps_dn, etab)
    return xi, xib

def rnd(*shape):
    return rng.normal(size=shape) + 1j * rng.normal(size=shape)

# --- Clifford relation and symmetry on V_10 = V + W_R^perp (real vectors) ---
def rho10_Sigma(u, xi, xib):
    e1, e2 = rho_v_on_Sigma(u[:4], xi, xib); f1, f2 = rho_w_on_Sigma(u[4:], xi, xib)
    return e1 + f1, e2 + f2
def rho10_Sigmastar(u, eta, etab):
    e1, e2 = rho_v_on_Sigmastar(u[:4], eta, etab); f1, f2 = rho_w_on_Sigmastar(u[4:], eta, etab)
    return e1 + f1, e2 + f2
cl = 0.0; sym = 0.0
for _ in range(5):
    u = rng.normal(size=10); u2 = rng.normal(size=10)
    xi, xib = rnd(2, 4), rnd(2, 4)
    a = rho10_Sigmastar(u, *rho10_Sigma(u2, xi, xib)); b = rho10_Sigmastar(u2, *rho10_Sigma(u, xi, xib))
    cl = max(cl, np.abs(a[0] + b[0] - 2 * u.dot(u2) * xi).max(), np.abs(a[1] + b[1] - 2 * u.dot(u2) * xib).max())
    z1, z2 = (rnd(2, 4), rnd(2, 4)), (rnd(2, 4), rnd(2, 4))
    r1 = rho10_Sigma(u, *z1); r2 = rho10_Sigma(u, *z2)
    sym = max(sym, abs(pair(*z1, *r2) - pair(*z2, *r1)))
res['Clifford_rho(u)rho(u2)+rho(u2)rho(u)=2(u,u2)_on_Sigma'] = float(cl)
res['(z1,rho(u)z2)=(z2,rho(u)z1)_commuting_components'] = float(sym)

# --- chirality: rho(e1 e2 e3 e4) on S_+ (x) W ---
def rho_vol_Sigma(xi, xib):
    E = np.eye(4)
    # rho(e1)rho(e2)rho(e3)rho(e4) acting on Sigma: alternate Sigma -> Sigma* -> Sigma
    z = (xi, xib)
    z = rho_v_on_Sigma(E[3], *z); z = rho_v_on_Sigmastar(E[2], *z); z = rho_v_on_Sigma(E[1], *z); z = rho_v_on_Sigmastar(E[0], *z)
    return z
xi = rnd(2, 4); z = rho_vol_Sigma(xi, 0 * xi)
res['rho(e1e2e3e4)_on_S+W_eigenvalue'] = complex(np.round((z[0] / xi).mean(), 12)).real
xib = rnd(2, 4); z = rho_vol_Sigma(0 * xib, xib)
res['rho(e1e2e3e4)_on_S-W*_eigenvalue'] = complex(np.round((z[1] / xib).mean(), 12)).real

# --- parity of the bilinears used in prop:c9-basis (commuting components) ---
q = rng.normal(size=4); w = rng.normal(size=6)
z1, z2 = (rnd(2, 4), rnd(2, 4)), (rnd(2, 4), rnd(2, 4))
def B(zA, zB, op):
    return pair(*zA, *op(zB))
vV = lambda z: rho_v_on_Sigma(q, *z)
vVvol = lambda z: rho_v_on_Sigma(q, *rho_vol_Sigma(*z))
vW = lambda z: rho_w_on_Sigma(w, *z)
vWvol = lambda z: rho_w_on_Sigma(w, *rho_vol_Sigma(*z))
for name, op in [('rho(q)', vV), ('rho(q)rho(vol)', vVvol), ('rho(w)', vW), ('rho(w)rho(vol)', vWvol)]:
    a, b = B(z1, z2, op), B(z2, z1, op)
    res['(z1,%s z2) vs (z2,%s z1)' % (name, name)] = 'symmetric' if abs(a - b) < 1e-10 else ('antisymmetric' if abs(a + b) < 1e-10 else 'neither')

# --- the two-component identities of rem:c9-twocomp3 (commuting components) ---
v = rng.normal(size=4)
z1, z2 = (rnd(2, 4), rnd(2, 4)), (rnd(2, 4), rnd(2, 4))
lhs = pair(*z1, *rho_v_on_Sigma(v, *z2))
rhs = 1j * sum(v[m] * (np.einsum('ai,da,di->', z1[0], sigb[m], z2[1]) + np.einsum('di,da,ai->', z1[1], sigb[m], z2[0])) for m in range(4))
res['(z1,rho(v)z2) = i v_m (z1 sigb z2bar + z1bar sigb z2)'] = float(abs(lhs - rhs))
lhs = pair(*z1, *rho_w_on_Sigma(w, *z2))
rhs = np.einsum('aj,ab,ji,bi->', z1[0], eps_up, wlow(w), z2[0]) - np.einsum('aj,ab,ji,bi->', z1[1], eps_up, wup(w), z2[1])
res['(z1,rho(w)z2) = z1 eps w_low z2 - z1bar eps w_up z2bar'] = float(abs(lhs - rhs))

# parity-odd structure: (zeta1, rho(e^m) D_m rho(vol) lambda) = i{ zeta1bar sigb D lambda - zeta1 sigb D lambdabar }
D = [(rnd(2, 4), rnd(2, 4)) for m in range(4)]   # D_m lambda components
zeta1 = (rnd(2, 4), rnd(2, 4)); E = np.eye(4)
lhs = sum(pair(*zeta1, *rho_v_on_Sigma(E[m], *rho_vol_Sigma(*D[m]))) for m in range(4))
rhs = 1j * sum(np.einsum('di,da,ai->', zeta1[1], sigb[m], D[m][0]) - np.einsum('ai,da,di->', zeta1[0], sigb[m], D[m][1]) for m in range(4))
res['parity-odd structure = i{([c,lb],sigb d l) - ([c,l],sigb d lb)} (orientation e1e2e3e4)'] = float(abs(lhs - rhs))

# --- rem:c9-twocompD ---
x = rng.normal(size=4)
lam_s = (rnd(2, 4), rnd(2, 4))  # lambda^* in Sigma*: (lam*^a_i, lamb*^{ad i})
# (rho(x) dslash lambda, lambda^*) with dslash lambda = sum_m rho(e^m) D_m
dsl = [rho_v_on_Sigma(E[m], *D[m]) for m in range(4)]
dsl = (sum(d[0] for d in dsl), sum(d[1] for d in dsl))
rx = rho_v_on_Sigmastar(x, *dsl)
lhs = pair(*rx, *lam_s)
sgsb = lambda n, m: sig[n] @ sigb[m]        # (sigma^n sigmabar^m)_a^b
sbsg = lambda m, n: sigb[m] @ sig[n]        # (sigmabar^m sigma^n)^{bd}_{ad}
rhs = 0
for n in range(4):
    for m in range(4):
        rhs += -x[n] * np.einsum('ab,bi,ai->', sgsb(n, m), D[m][0], lam_s[0])
        rhs += x[n] * np.einsum('bi,ba,ai->', D[m][1], sbsg(m, n), lam_s[1])
res['(rho(x)dslash lambda, lambda*) two-component formula'] = float(abs(lhs - rhs))
# (rho(x) lambda^*, lambda^*) = -2i x_m lam*^a_i sig^m_{a ad} lamb*^{ad i}   (lambda^* components even)
lhs = pair(*rho_v_on_Sigmastar(x, *lam_s), *lam_s)
rhs = -2j * sum(x[m] * np.einsum('ai,ad,di->', lam_s[0], sig[m], lam_s[1]) for m in range(4))
res['(rho(x)lambda*,lambda*) = -2i x_m lam* sig^m lamb*'] = float(abs(lhs - rhs))

# --- Q(int(rho(x)lambda*,lambda*)) = 2 int(rho(x) dslash lambda, lambda*): needs (rho(x) a, b) = (rho(x) b, a) on Sigma* ---
a_, b_ = (rnd(2, 4), rnd(2, 4)), (rnd(2, 4), rnd(2, 4))
res['(rho(x)a,b) symmetric on Sigma* (commuting)'] = float(abs(pair(*rho_v_on_Sigmastar(x, *a_), *b_) - pair(*rho_v_on_Sigmastar(x, *b_), *a_)))

# --- colour identity of lem:c9-legs on su(3) ---
gm = np.zeros((8, 3, 3), complex)
gm[0][0, 1] = gm[0][1, 0] = 1; gm[1][0, 1] = -1j; gm[1][1, 0] = 1j; gm[2][0, 0] = 1; gm[2][1, 1] = -1
gm[3][0, 2] = gm[3][2, 0] = 1; gm[4][0, 2] = -1j; gm[4][2, 0] = 1j; gm[5][1, 2] = gm[5][2, 1] = 1
gm[6][1, 2] = -1j; gm[6][2, 1] = 1j; gm[7] = np.diag([1, 1, -2]) / np.sqrt(3)
T = -0.5j * gm   # anti-Hermitian basis of su(3)
f = np.zeros((8, 8, 8))
for a in range(8):
    for b in range(8):
        C = T[a] @ T[b] - T[b] @ T[a]
        for c in range(8):
            f[a, b, c] = np.real(-2 * np.trace(C @ T[c]))   # [T_a,T_b] = f_abc T_c with tr(T_aT_b) = -1/2 delta
ad = lambda X: np.einsum('a,abc->cb', X, f)   # (ad_X)_{cb}: ad_X T_b = X_a f_abc T_c
pairing = lambda X, Y: -2 * np.real(np.trace(np.einsum('a,aij->ij', X, T) @ np.einsum('a,aij->ij', Y, T)))  # (X,Y) = X.Y
br = lambda X, Y: np.einsum('a,b,abc->c', X, Y, f)
X, Y, Z = rng.normal(size=8), rng.normal(size=8), rng.normal(size=8)
ratio = np.trace(ad(X) @ ad(Y)) / pairing(X, Y)
res['colour: tr(adX adY adZ) - (1/2) ratio (X,[Y,Z])'] = float(abs(np.trace(ad(X) @ ad(Y) @ ad(Z)) - 0.5 * ratio * pairing(X, br(Y, Z))))
res['colour: ratio tr(adX adY)/(X,Y) (su(3), (X,Y)=X.Y)'] = float(ratio)

print(json.dumps(res, indent=1))
json.dump(res, open(__file__.replace('.py', '_results.json'), 'w'), indent=1)
