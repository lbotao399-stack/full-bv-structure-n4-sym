# Audit #2: large-distance behaviour of a PAIR of floating odd insertions joined by a lambda-lambda edge P(eps,infty)
# (a configuration not covered by design Thm thm:o-ir (iv), which lists single floating S^(1)(eta) (|y|^-5), pairs
# joined by a special edge (|y|^-6), and S^(2) (never floating)).
#
# Floating odd vertex (design Lemma lem:o-legs, two-component):
#   -i( lambar_{adot i} sigbar^{m adot a} eps_a^i(y) + epsbar_{adot i}(y) sigbar^{m adot a} lam_a^i , A*_m ),
#   eps_a^i(y) = xi_a^i - i y_m sig^m_{a adot} etabar^{adot i},  epsbar_{adot i}(y) = xibar_{adot i} + i y_m eta^a_i sig^m_{a adot}.
# Chain on a cycle:  c(z) --[A*|_V c]-- (vertex 1 at y1) --[lam lam]-- (vertex 2 at y2) --[A*|_V c]-- c(z').
# Fixed-orientation part with chirality (lam at y1, lambar at y2), lower-order constants dropped:
#   C1[z,z'] = int int d_mD(y1-z) y1_a d_kD(y1-y2) y2_b d_nD(y2-z') M^{amknb},
#   M^{amknb} = sum_i eta1^a_i (sig^a sigbar^m sig^k sigbar^n sig^b) etabar2^{i}.
# Wick signs + graded symmetry of P(lam lam) + symmetry of Gamma give the physical sum over both orientations and both
# chiralities (derivation in the audit verdict):  Total = C1(eta1,etabar2)[z,z'] - C1(eta1,etabar2)[z',z] + (1<->2).
# Momentum space (hhat(p)=int e^{-ipy}h, D^=1/p^2): with f(z), g(z') the two attachments of the rest of the cycle,
#   C1 -> int dp/(2pi)^4 i * I_fg(p),
#   I_fg(p) = [d_{q_a}(q_m fhat(q)/q^2)]_{q=-p} (p_k/p^2) [d_{p_b}(p_n ghat(p)/p^2)] M^{amknb}.
# fhat(q)=f0+f1.q+..., ghat(p)=g0+g1.p+... (entire: compactly supported rest).  IR region = p -> 0.
import numpy as np, json

rng = np.random.default_rng(20261028)
s1 = np.array([[0,1],[1,0]], complex); s2 = np.array([[0,-1j],[1j,0]], complex); s3 = np.array([[1,0],[0,-1]], complex)
I2 = np.eye(2, dtype=complex)
sig = [s1, s2, s3, 1j*I2]          # sig^m_{a adot}, Euclidean (dictionary)
sigb = [-s1, -s2, -s3, 1j*I2]      # sigbar^{m adot a}
chk = max(np.abs(sig[m]@sigb[n] + sig[n]@sigb[m] + 2*(m == n)*I2).max() for m in range(4) for n in range(4))

def sphere_nodes(ns=24, nphi=24):
    # S^3 in C^2: p = (cos al e^{i ph1}, sin al e^{i ph2}); normalised measure = ds dph1 dph2/(4pi^2), s = sin^2 al
    xs, ws = np.polynomial.legendre.leggauss(ns); s = 0.5*(xs+1); w = 0.5*ws
    ph = 2*np.pi*np.arange(nphi)/nphi
    S, P1, P2 = np.meshgrid(s, ph, ph, indexing='ij'); W = np.meshgrid(w, ph, ph, indexing='ij')[0]/nphi**2
    ca, sa = np.sqrt(1-S), np.sqrt(S)
    pts = np.stack([ca*np.cos(P1), ca*np.sin(P1), sa*np.cos(P2), sa*np.sin(P2)], -1).reshape(-1, 4)
    return pts, W.reshape(-1)
NODES, WTS = sphere_nodes()
assert abs(WTS.sum()-1) < 1e-13

def mats(p, f0, f1, g0, g1, fgrad=None):
    """left (2x2 alpha->alpha), mid (alpha->alphadot), right (alphadot->alphadot) at momenta p (N,4).
    Linear Taylor data (f0,f1),(g0,g1) unless fgrad given: fgrad=(F(q),dF(q),G(p),dG(p)) callables."""
    N = p.shape[0]; pp = (p**2).sum(1)
    q = -p; qq = pp
    if fgrad is None:
        Fq = f0 + q@f1; dFq = np.tile(f1, (N, 1)); Gp = g0 + p@g1; dGp = np.tile(g1, (N, 1))
    else:
        Fq, dFq, Gp, dGp = fgrad[0](q), fgrad[1](q), fgrad[2](p), fgrad[3](p)
    L = np.zeros((N, 2, 2), complex); Mid = np.zeros((N, 2, 2), complex); R = np.zeros((N, 2, 2), complex)
    for a in range(4):
        for m in range(4):
            la = ((a == m)/qq - 2*q[:, a]*q[:, m]/qq**2)*Fq + q[:, m]/qq*dFq[:, a]
            L += la[:, None, None]*(sig[a]@sigb[m])[None]
            rb = ((a == m)/pp - 2*p[:, a]*p[:, m]/pp**2)*Gp + p[:, m]/pp*dGp[:, a]   # here a=b index, m=n index
            R += rb[:, None, None]*(sigb[m]@sig[a])[None]
    for k in range(4):
        Mid += (p[:, k]/pp)[:, None, None]*sig[k][None]
    return L, Mid, R

def I_sym(p, e1, eb2, *args, **kw):
    L, Mid, R = mats(p, *args, **kw)
    Mt = L@Mid@R
    return np.einsum('ai,nab,bi->n', e1, Mt, eb2)   # e1[alpha,i], eb2[alphadot,i]

def bil(e1, eb2, m):  # eta1 sig^m etabar2 summed over i
    return np.einsum('ai,ab,bi->', e1, sig[m], eb2)

res = {'euclid_clifford_check_max_residual': float(chk)}
# random spinor parameters eta=(eta^a_i, etabar^{adot i}) for eta1, eta2
def rnd(): return rng.normal(size=(2, 4)) + 1j*rng.normal(size=(2, 4))
e1, eb1, e2, eb2 = rnd(), rnd(), rnd(), rnd()

# (1) degree -5 part of a single orientation: f0=g0=1, f1=g1=0. pointwise nonzero, odd -> zero sphere average
z4 = np.zeros(4)
deg5 = I_sym(NODES, e1, eb2, 1.0, z4, 1.0, z4)
pred5 = 4*np.einsum('nm,m->n', NODES, np.array([bil(e1, eb2, m) for m in range(4)]))   # hand result 4 (eta1 pslash etabar2)/p^6 at |p|=1
res['deg-5_single_orientation'] = {'max_abs_on_S3': float(np.abs(deg5).max()), 'sphere_avg_abs': float(abs((deg5*WTS).sum())),
                                   'max_dev_from_hand_formula_4(eta1 pslash etabar2)': float(np.abs(deg5-pred5).max())}
# (2) degree -4 pieces: B_c (f0=1, g1=e_c) and C_c (f1=e_c, g0=1), all others 0
B = np.zeros((NODES.shape[0], 4), complex); C = np.zeros_like(B)
for c in range(4):
    ec = np.eye(4)[c]
    B[:, c] = I_sym(NODES, e1, eb2, 1.0, z4, 0.0, ec)
    C[:, c] = I_sym(NODES, e1, eb2, 0.0, ec, 1.0, z4)
avgBC = ((B-C)*WTS[:, None]).sum(0)
pred = 6*np.array([bil(e1, eb2, c) for c in range(4)])
res['deg-4_antisymmetrised_orientation_sum'] = {
    'statement': 'I_fg - I_gf = (f0 g1 - f1 g0)_c (B_c - C_c) + O(p^-3); degree -5 cancels identically',
    'sphere_avg_(B-C)_c': [complex(x).__repr__() for x in avgBC],
    'hand_prediction_6(eta1 sig^c etabar2)': [complex(x).__repr__() for x in pred],
    'max_rel_dev': float(np.abs(avgBC-pred).max()/np.abs(pred).max()),
    'pointwise_hand_formula_8(eta1 pslash etabar2)p_c + 4(eta1 sig^c etabar2)_max_dev': float(np.abs(
        (B-C) - (8*np.einsum('nm,m->n', NODES, np.array([bil(e1, eb2, m) for m in range(4)]))[:, None]*NODES
                 + 4*np.array([bil(e1, eb2, c) for c in range(4)])[None, :])).max())}
# (3) symmetrised over eta1<->eta2: log coefficient proportional to v1^c = -i(eta2 sig^c etabar1 + eta1 sig^c etabar2)
v1 = np.array([-1j*(bil(e2, eb1, c)+bil(e1, eb2, c)) for c in range(4)])
def BC_total(e1_, eb1_, e2_, eb2_):
    out = np.zeros((NODES.shape[0], 4), complex)
    for (ea, ebb) in [(e1_, eb2_), (e2_, eb1_)]:
        for c in range(4):
            ec = np.eye(4)[c]
            out[:, c] += I_sym(NODES, ea, ebb, 1.0, z4, 0.0, ec) - I_sym(NODES, ea, ebb, 0.0, ec, 1.0, z4)
    return out
T = BC_total(e1, eb1, e2, eb2)
avgT = (T*WTS[:, None]).sum(0)
res['symmetrised_total'] = {'sphere_avg': [complex(x).__repr__() for x in avgT],
                            'equals_6i*v1': bool(np.allclose(avgT, 6j*v1, rtol=1e-10, atol=1e-10)),
                            'v1': [complex(x).__repr__() for x in v1]}
# (4) pairs with v1 = Gamma(eta1,eta2)|_V = 0: degree -4 part vanishes pointwise
# solve for eta2 in the kernel of eta2 -> v1 (4 linear conditions on 16 unknowns)
basis = []
for idx in range(16):
    E = np.zeros(16, complex); E[idx] = 1
    e2t, eb2t = E[:8].reshape(2, 4), E[8:].reshape(2, 4)
    basis.append([-1j*(bil(e2t, eb1, c)+bil(e1, eb2t, c)) for c in range(4)])
Amat = np.array(basis).T  # 4 x 16
_, sv, Vh = np.linalg.svd(Amat)
ker = Vh[4:].conj().T  # 16 x 12
coef = rng.normal(size=12)+1j*rng.normal(size=12); vec = ker@coef
e2k, eb2k = vec[:8].reshape(2, 4), vec[8:].reshape(2, 4)
v1k = np.array([-1j*(bil(e2k, eb1, c)+bil(e1, eb2k, c)) for c in range(4)])
Tk = BC_total(e1, eb1, e2k, eb2k)
res['v1_zero_pairs'] = {'|v1|': float(np.abs(v1k).max()), 'max_pointwise_deg-4': float(np.abs(Tk).max()),
                        'compare_generic_pair_max_pointwise': float(np.abs(T).max())}
# (5) radial check with the exact (non-Taylor) rest: fhat(q)=exp(-|q-a|^2), ghat(p)=exp(-|p-c|^2)
aa = np.array([0.3, -0.7, 0.2, 0.5]); cc = np.array([-0.4, 0.1, 0.6, -0.2])
Ff = lambda q: np.exp(-((q-aa)**2).sum(1)); dFf = lambda q: -2*(q-aa)*Ff(q)[:, None]
Gg = lambda p: np.exp(-((p-cc)**2).sum(1)); dGg = lambda p: -2*(p-cc)*Gg(p)[:, None]
def total_exact(p):
    tot = 0
    for (ea, ebb) in [(e1, eb2), (e2, eb1)]:
        tot = tot + I_sym(p, ea, ebb, None, None, None, None, fgrad=(Ff, dFf, Gg, dGg)) \
                  - I_sym(p, ea, ebb, None, None, None, None, fgrad=(Gg, dGg, Ff, dFf))
    return tot
xr, wr = np.polynomial.legendre.leggauss(40)
shells = []
for rho in [1e-1, 1e-2, 1e-3, 1e-4]:
    # int_{rho<|p|<2rho} d^4p = 2pi^2 int r^3 dr <.>_S3 ; use log-radius t in [log rho, log 2rho]
    t = np.log(rho) + (xr+1)/2*np.log(2); wt = wr/2*np.log(2)
    val = 0
    for ti, wi in zip(t, wt):
        r = np.exp(ti); vals = total_exact(r*NODES)
        val += wi*r**4*(vals*WTS).sum()
    shells.append(2*np.pi**2*val)
f0, g0 = Ff(np.zeros((1, 4)))[0], Gg(np.zeros((1, 4)))[0]
f1, g1 = dFf(np.zeros((1, 4)))[0], dGg(np.zeros((1, 4)))[0]
pred_shell = 2*np.pi**2*np.log(2)*(6j*v1)@(f0*g1 - f1*g0)
res['radial_shells_exact_rest'] = {'shell_integrals_rho=1e-1..1e-4': [complex(x).__repr__() for x in shells],
                                   'prediction_2pi^2 log2 * 6i v1.(f0 g1 - f1 g0)': complex(pred_shell).__repr__(),
                                   'rel_dev_at_1e-4': float(abs(shells[-1]-pred_shell)/abs(pred_shell))}
# (6) control: (xi, eta) pair, xi constant (undotted xi_a^i) at y1, eta linear at y2 (epsbar from eta^a_i)
#   structure eta sig^b sigbar^n sig^k sigbar^m xi ; left factor without derivative.  Degree -4 single orientation,
#   cancels in the orientation-antisymmetrised sum -> convergent (degree -3).
def I_xieta(p, et, xi, f0_, f1_, g0_, g1_):
    N = p.shape[0]; pp = (p**2).sum(1); q = -p
    Fq = f0_ + q@f1_; Gp = g0_ + p@g1_
    Lm = np.zeros((N, 2, 2), complex)
    for m in range(4):
        Lm += (q[:, m]/pp*Fq)[:, None, None]*sigb[m][None]       # placed last in the product
    Mid = sum((p[:, k]/pp)[:, None, None]*sig[k][None] for k in range(4))
    Rb = np.zeros((N, 2, 2), complex)
    for b_ in range(4):
        for n in range(4):
            rb = ((b_ == n)/pp - 2*p[:, b_]*p[:, n]/pp**2)*Gp + p[:, n]/pp*g1_[b_]
            Rb += rb[:, None, None]*(sig[b_]@sigb[n])[None]
    Mt = Rb@Mid@Lm
    return np.einsum('ai,nab,bi->n', et, Mt, xi)
xi = rnd(); et = rnd()
single = I_xieta(NODES, et, xi, 1.0, z4, 1.0, z4)
antis = single - I_xieta(NODES, et, xi, 1.0, z4, 1.0, z4)
res['control_xi_eta'] = {'single_orientation_deg-4_sphere_avg': complex((single*WTS).sum()).__repr__(),
                         'pairing_eta.xi': complex(np.einsum('ai,ai->', et, xi)).__repr__(),
                         'antisymmetrised_deg-4_max': float(np.abs(antis).max())}
print(json.dumps(res, indent=1))
json.dump(res, open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/design_audit_2/eta_pair_ir.json', 'w'), indent=1)
