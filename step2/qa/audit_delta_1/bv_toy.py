# Auditor #1, Part D, D1 (a),(b),(c): independent finite-dimensional BV check.
# Conventions taken from the manuscript only:
#   * odd pairing (Def c1-pairing): <even member, odd member> = +fibre pairing, <odd, even> = -<even, odd>;
#     slots: A (even) / A^* (odd), lambda^* (even) / lambda (odd), c^* (even) / c (odd).
#   * convolution [C07 §6.3]: K*e = (-1)^{|e|} sum K' <K'', e>;  K_0 = kernel of the identity (K_0* = id).
#   * partial_phi = 1/2 sum partial_{phi''} partial_{phi'} [C07 §6.4];  Delta_K = -partial_K;
#     {f,g}_K = Delta_K(fg) - Delta_K(f) g - (-1)^{|f|} f Delta_K(g)  (Def c1-bracket, Def c3-delta).
#   * derivation induced by an operator N [C07 §6.2]: D_N = sum_j (N x)^j d/dx_j (left derivatives, coefficient left).
#   * (1 (x) A)(e_i (x) e_j) = (-1)^{|A| |e_i|} e_i (x) A e_j  (Koszul).
#   * W(P,J) = hbar log(e^{hbar partial_P} e^{J/hbar})  (Def c3-W).
# The toy has the slot pattern of the N=4 BV complex (one 'momentum mode' replaced by small matrices):
#   A: 2 even, A^*: 2 odd, c: 1 odd, c^*: 1 even, lambda: 2 odd, lambda^*: 2 even;  xi (odd) = Delta-ghost parameter.
import json, sys, random, os
from fractions import Fraction as Fr
from grassmann import Alg

random.seed(20261028)
NF = 10
# index: 0 xi | 1 A1 2 A2 | 3 As1 4 As2 | 5 c | 6 cs | 7 l1 8 l2 | 9 ls1 10 ls2
PAR = [1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0]
A_, AS_, C_, CS_, L_, LS_ = [1, 2], [3, 4], [5], [6], [7, 8], [9, 10]
FIELDS = list(range(1, NF + 1))
NAMES = {0: 'xi', 1: 'A1', 2: 'A2', 3: 'As1', 4: 'As2', 5: 'c', 6: 'cs', 7: 'l1', 8: 'l2', 9: 'ls1', 10: 'ls2'}
R = range(NF + 1)

def zmat():
    return [[Fr(0)] * (NF + 1) for _ in R]

def mm(X, Y):
    Z = zmat()
    for i in R:
        for k in R:
            if X[i][k]:
                for j in R:
                    if Y[k][j]:
                        Z[i][j] += X[i][k] * Y[k][j]
    return Z

def madd(*Xs):
    Z = zmat()
    for X in Xs:
        for i in R:
            for j in R:
                Z[i][j] += X[i][j]
    return Z

def msc(X, c):
    return [[X[i][j] * c for j in R] for i in R]

def mpar(X):
    ps = set()
    for i in R:
        for j in R:
            if X[i][j]:
                ps.add((PAR[i] + PAR[j]) % 2)
    assert len(ps) <= 1, "inhomogeneous operator"
    return ps.pop() if ps else 0

def gcomm(X, Y):
    s = (-1) ** (mpar(X) * mpar(Y))
    return madd(mm(X, Y), msc(mm(Y, X), -s))

def iszero(X):
    return all(X[i][j] == 0 for i in R for j in R)

# ---------------- pairing, K_0, convolution ----------------
OMEGA = zmat()
for y, eta in zip(A_ + LS_ + CS_, AS_ + L_ + C_):
    OMEGA[y][eta] = Fr(1)
    OMEGA[eta][y] = Fr(-1)

def conv(phi):
    """matrix of phi* : (phi*)[i][l] = (-1)^{p_l} sum_j phi^{ij} omega[j][l]"""
    Z = zmat()
    for (i, j), c in phi.items():
        for l in R:
            if OMEGA[j][l]:
                Z[i][l] += c * OMEGA[j][l] * (-1) ** PAR[l]
    return Z

# K_0: solve sum_j k^{ij} omega[j][l] = (-1)^{p_l} delta_il  (omega is a signed permutation here)
K0 = {}
for i in FIELDS:
    for l in FIELDS:
        if i == l:
            # find j with omega[j][l] != 0
            for j in FIELDS:
                if OMEGA[j][l]:
                    K0[(i, j)] = Fr((-1) ** PAR[l]) / OMEGA[j][l]

def ker_L(Aop, phi):
    out = {}
    for (k, j), c in phi.items():
        for i in R:
            if Aop[i][k]:
                out[(i, j)] = out.get((i, j), 0) + Aop[i][k] * c
    return {k: v for k, v in out.items() if v}

def ker_R(Aop, phi):
    pa = mpar(Aop)
    out = {}
    for (k, l), c in phi.items():
        for j in R:
            if Aop[j][l]:
                out[(k, j)] = out.get((k, j), 0) + (-1) ** (pa * PAR[k]) * Aop[j][l] * c
    return {k: v for k, v in out.items() if v}

def ker_sym_action(Aop, phi):
    a, b = ker_L(Aop, phi), ker_R(Aop, phi)
    out = dict(a)
    for k, v in b.items():
        out[k] = out.get(k, 0) + v
    return {k: v for k, v in out.items() if v}

def ker_add(*phis, coeffs=None):
    out = {}
    coeffs = coeffs or [1] * len(phis)
    for phi, c in zip(phis, coeffs):
        for k, v in phi.items():
            out[k] = out.get(k, 0) + c * v
    return {k: v for k, v in out.items() if v}

def ker_par(phi):
    ps = {(PAR[i] + PAR[j]) % 2 for (i, j) in phi}
    assert len(ps) <= 1
    return ps.pop() if ps else 0

def graded_symmetric(phi):
    for (i, j), v in phi.items():
        if phi.get((j, i), 0) != (-1) ** (PAR[i] * PAR[j]) * v:
            return False
    return True

def pairing(Aop, k, l):
    """<A e_k, e_l>"""
    return sum(Aop[j][k] * OMEGA[j][l] for j in R)

def pairing2(k, Aop, l):
    """<e_k, A e_l>"""
    return sum(OMEGA[k][j] * Aop[j][l] for j in R)

ALG = Alg(PAR, gmax=2)

def rand_poly(alg, degree, nterms, parity=0, vars_=FIELDS, rng=random):
    p = {}
    tries = 0
    while len(p) < nterms and tries < 1000:
        tries += 1
        mono = [0] * (alg.NV + 3)
        vs = [rng.choice(vars_) for _ in range(degree)]
        ok = True
        for v in vs:
            if PAR[v] and mono[v]:
                ok = False
                break
            mono[v] += 1
        if not ok:
            continue
        if sum(mono[v] for v in range(alg.NV) if PAR[v]) % 2 != parity:
            continue
        p[tuple(mono)] = Fr(rng.choice([-3, -2, -1, 1, 2, 3]))
    return p

def lin(alg, coeffs):
    return alg.add(*[alg.var(k, c) for k, c in coeffs.items() if c])

def op_from_hamiltonian(alg, h, K=K0):
    """matrix N with {h, x_j}_K = (N x)^j"""
    N = zmat()
    for j in FIELDS:
        b = bracket(alg, h, alg.var(j), K)
        for key, v in b.items():
            idx = [i for i in range(alg.NV) if key[i]]
            assert len(idx) == 1 and key[idx[0]] == 1 and key[alg.H] == 0, (key, v)
            N[j][idx[0]] += v
    return N

def Delta(alg, f, K):
    return alg.scale(alg.dker(f, K), -1)

def bracket(alg, f, g, K=K0):
    pf = alg.parity(f)
    return alg.add(Delta(alg, alg.mul(f, g), K), alg.scale(alg.mul(Delta(alg, f, K), g), -1),
                   alg.scale(alg.mul(f, Delta(alg, g, K)), -((-1) ** pf)))

def bracket_phi(alg, f, g, phi):
    """{f,g}_phi := partial_phi(fg) - (partial_phi f) g - (-1)^{|f|} f partial_phi g"""
    pf = alg.parity(f)
    return alg.add(alg.dker(alg.mul(f, g), phi), alg.scale(alg.mul(alg.dker(f, phi), g), -1),
                   alg.scale(alg.mul(f, alg.dker(g, phi)), -((-1) ** pf)))

REPORT = {"conventions": "see header of bv_toy.py", "checks": {}}

def record(name, ok, note=""):
    REPORT["checks"][name] = {"pass": bool(ok), "note": note}
    print(("PASS " if ok else "FAIL ") + name + ("  -- " + note if note else ""))
    sys.stdout.flush()

# ===================== T1: conventions reproduce Prop c1-Q =====================
# free action (left entry first): -1/4 a^T G a + 1/2 l^T Gamma l + (Dc, A^*)
Dv = [Fr(0), Fr(1), Fr(2)]           # D: c -> A  (column (1,2))
# G symmetric with G D = 0:  G = (2,-1)^T (2,-1)
Gm = {(1, 1): Fr(4), (1, 2): Fr(-2), (2, 1): Fr(-2), (2, 2): Fr(1)}
gam = Fr(3)                           # Gamma = [[0,gam],[-gam,0]] on (l1,l2)
a = ALG
S0 = {}
for (i, j), v in Gm.items():
    S0 = a.add(S0, a.scale(a.mul(a.var(i), a.var(j)), Fr(-1, 4) * v))
S0 = a.add(S0, a.scale(a.mul(a.var(7), a.var(8)), Fr(1, 2) * gam), a.scale(a.mul(a.var(8), a.var(7)), Fr(1, 2) * (-gam)))
for ia, ias in zip(A_, AS_):
    S0 = a.add(S0, a.scale(a.mul(a.var(5), a.var(ias)), Dv[ia]))   # (Dc)^i As_i ;  c * As ordering, left entry first
# Q from the table of Prop c1-Q
Qm = zmat()
for ia in A_:
    Qm[ia][5] = -Dv[ia]                       # c -> -dc
for (i, j), v in Gm.items():
    Qm[AS_[i - 1]][j] = Fr(-1, 2) * v         # A -> -1/2 d^*d A  (into A^*)
for ia, ias in zip(A_, AS_):
    Qm[6][ias] = Dv[ia]                       # A^* -> d^*(A^*)  (into c^*)
Qm[9][7], Qm[9][8], Qm[10][7], Qm[10][8] = Fr(0), gam, -gam, Fr(0)   # lambda -> Gamma lambda (into lambda^*)
assert mpar(Qm) == 1
Qderived = op_from_hamiltonian(a, S0)
record("T1a {1/2<e,Qe>,-} = D_Q with Q = table of Prop c1-Q (c->-dc, A->-1/2 d*d A, A*->d*A*, lambda->dslash lambda)",
       Qderived == Qm)
record("T1b Q^2 = 0", iszero(mm(Qm, Qm)))
okQ = all(pairing(Qm, k, l) + (-1) ** PAR[k] * pairing2(k, Qm, l) == 0 for k in FIELDS for l in FIELDS)
record("T1c <Qe0,e1> + (-1)^{|e0|}<e0,Qe1> = 0", okQ)
record("T1d (Q(x)1 + 1(x)Q) K_0 = 0", not ker_sym_action(Qm, K0))
K0conv = conv(K0)
record("T1e K_0* = id (C07 convolution convention) and K_0 = -delta x (inverse fibre pairing) on A(x)A*, A*(x)A, ...",
       all(K0conv[i][j] == (1 if i == j and i > 0 else 0) for i in R for j in R) and
       all(K0[(y, e)] == -1 and K0[(e, y)] == -1 for y, e in zip(A_ + LS_ + CS_, AS_ + L_ + C_)) and len(K0) == 10)
f_rand = a.add(rand_poly(a, 3, 8, 0), rand_poly(a, 2, 5, 1), rand_poly(a, 4, 6, 1))
record("T1f {S0,f} = D_Q f on a random polynomial", a.eq(bracket(a, S0, f_rand), a.vf(f_rand, Qm)))

# ===================== T0: sign of {S^(1)(Delta), -} =====================
# M block diagonal on the field slots, antifield blocks fixed by preservation of <,> (Def c7-vdot)
def rnd():
    return Fr(random.choice([-3, -2, -1, 1, 2, 3]), random.choice([1, 2]))
Mm = zmat()
for blk, ablk in [(A_, AS_), (L_, LS_), (C_, CS_)]:
    for i in blk:
        for j in blk:
            Mm[i][j] = rnd()
for blk, ablk in [(A_, AS_), (L_, LS_), (C_, CS_)]:
    for ii, i in enumerate(ablk):
        for jj, j in enumerate(ablk):
            Mm[i][j] = -Mm[blk[jj]][blk[ii]]       # M on antifields = - transpose
SM = {}
for blk, ablk, sgn in [(A_, AS_, 1), (L_, LS_, -1), (C_, CS_, -1)]:
    for ii, i in enumerate(blk):
        Mx_i = lin(a, {j: Mm[i][j] for j in blk})
        SM = a.add(SM, a.scale(a.mul(Mx_i, a.var(ablk[ii])), sgn))
# S_M = (M A, A^*) - (M lambda, lambda^*) - (M c, c^*)   (left entry first; fibre pairing = identity)
record("T0a M preserves <,>:  <Me0,e1> + <e0,Me1> = 0", all(pairing(Mm, k, l) + pairing2(k, Mm, l) == 0 for k in FIELDS for l in FIELDS))
record("T0b (M(x)1 + 1(x)M) K_0 = 0", not ker_sym_action(Mm, K0))
record("T0c S_M odd", a.parity(SM) == 1)
Mderived = op_from_hamiltonian(a, SM)
record("T0d {S_M, -} = +D_M  (sign +1) for S_M = (MA,A*) - (M lambda,lambda*) - (Mc,c*)", Mderived == Mm)
record("T0e {S_M,f} = D_M f on a random polynomial", a.eq(bracket(a, SM, f_rand), a.vf(f_rand, Mm)))
record("T0f (control) {S_M,f} != -D_M f", not a.eq(bracket(a, SM, f_rand), a.scale(a.vf(f_rand, Mm), -1)))

# ===================== T2: kernel / operator identities =====================
def rand_sym_kernel(parity, pattern=None, rng=random):
    phi = {}
    pairs = pattern if pattern is not None else [(i, j) for i in FIELDS for j in FIELDS if i <= j and (PAR[i] + PAR[j]) % 2 == parity]
    for (i, j) in pairs:
        if rng.random() < 0.7:
            v = Fr(rng.choice([-3, -2, -1, 1, 2, 3]), rng.choice([1, 2, 3]))
            phi[(i, j)] = v
            if i != j:
                phi[(j, i)] = (-1) ** (PAR[i] * PAR[j]) * v
            elif PAR[i] == 1:
                phi.pop((i, j))          # odd-odd diagonal must vanish for graded symmetry
    return phi

Pgen = rand_sym_kernel(0)
assert graded_symmetric(Pgen)
phi_odd = rand_sym_kernel(1)
assert graded_symmetric(phi_odd)
for nm, phi in [("even", Pgen), ("odd", phi_odd)]:
    lhs = a.add(a.dker(a.vf(f_rand, Qm), phi), a.scale(a.vf(a.dker(f_rand, phi), Qm), -((-1) ** (1 * ker_par(phi)))))
    rhs = a.dker(f_rand, ker_sym_action(Qm, phi))
    record("T2a [partial_phi, D_Q] = partial_{(Q(x)1+1(x)Q)phi}  (phi %s)" % nm, a.eq(lhs, rhs))
    lhs = a.sub(a.vf(a.dker(f_rand, phi), Mm), a.dker(a.vf(f_rand, Mm), phi))
    rhs = a.scale(a.dker(f_rand, ker_sym_action(Mm, phi)), -1)
    record("T2b [D_M, partial_phi] = -partial_{(M(x)1+1(x)M)phi}  (phi %s)" % nm, a.eq(lhs, rhs))
    # (Q phi)* = [Q, phi*]  graded, operator parity = kernel parity + 1
    lhsm = conv(ker_sym_action(Qm, phi))
    rhsm = gcomm(Qm, conv(phi))
    record("T2c ((Q(x)1+1(x)Q)phi)* = [Q, phi*] (graded; phi* has parity |phi|+1)  (phi %s)" % nm, lhsm == rhsm)
    record("T2d ((M(x)1+1(x)M)phi)* = [M, phi*]  (phi %s)" % nm, conv(ker_sym_action(Mm, phi)) == gcomm(Mm, conv(phi)))
record("T2e P* is odd for even P (so [Q,P*] = QP* + P*Q)", mpar(conv(Pgen)) == 1)

# ===================== flow helpers =====================
def flow(alg, P, J):
    Z = alg.exp_dker(alg.exp(alg.scale(J, 1, h=-1)), P)
    Y = alg.sub(Z, alg.const(1))
    W = alg.scale(alg.log1p(Y), 1, h=1)
    Zinv = alg.inv1p(Y)
    return Z, Zinv, W

def Wprime(alg, P, J, Z, Zinv, X):
    """W'_J(X) = Z^{-1} e^{hbar partial_P}(X e^{J/hbar})"""
    return alg.mul(Zinv, alg.exp_dker(alg.mul(X, alg.exp(alg.scale(J, 1, h=-1))), P))

def make_J(alg, rng=random):
    J3 = rand_poly(alg, 3, 7, 0, rng=rng)
    J4 = rand_poly(alg, 4, 5, 0, rng=rng)
    J2 = rand_poly(alg, 2, 4, 0, rng=rng)
    return alg.add(alg.scale(J3, 1, g=1), alg.scale(J4, 1, g=2), alg.scale(J2, 1, h=1, g=1))

# N=4 slot pattern of P(eps,infty) (Prop c3-P): A(x)A, c(x)A*, A*(x)c, lambda(x)lambda only
PATTERN_N4 = [(i, j) for i in A_ for j in A_ if i <= j] + [(5, j) for j in AS_] + [(7, 8)]
Pn4 = rand_sym_kernel(0, pattern=PATTERN_N4)
for k in list(Pn4):
    pass
assert graded_symmetric(Pn4)

# ===================== T3: D1(a) =====================
def check_a(P, label, SMv=SM, Mv=Mm):
    J = make_J(a)
    Z, Zinv, W = flow(a, P, J)
    assert a.min_hbar(W) >= 0
    Jx = a.add(J, a.mul(a.var(0), SMv))
    Zx, Zxinv, Wx = flow(a, P, Jx)
    X_lhs = a.d(Wx, 0)                                   # left xi-derivative
    X_def = Wprime(a, P, J, Z, Zinv, SMv)
    record("T3.0 d/dxi W(P,J+xi S_M) = Z^{-1} e^{hbar d_P}(S_M e^{J/hbar})  [%s]" % label, a.eq(X_lhs, X_def))
    # conjugation terminates
    c1 = a.sub(a.dker(a.mul(SMv, f_rand), P), a.mul(SMv, a.dker(f_rand, P)))     # [d_P, S_M] f
    N = mm(conv(P), Mv)
    c1_pred = a.add(a.mul(a.dker(SMv, P), f_rand), a.vf(f_rand, N))
    record("T3.1 [partial_P, S_M] = (partial_P S_M) + D_N,  N = P* o M  [%s]" % label, a.eq(c1, c1_pred))
    phi = ker_sym_action(N, P)
    c2 = a.sub(a.dker(a.vf(f_rand, N), P), a.vf(a.dker(f_rand, P), N))          # [d_P, D_N] (D_N odd, d_P even)
    record("T3.2 [partial_P, D_N] = partial_phi,  phi = (N(x)1+1(x)N)P;  [partial_P, partial_phi] = 0  [%s]" % label,
           a.eq(c2, a.dker(f_rand, phi)) and a.eq(a.dker(a.dker(f_rand, phi), P), a.dker(a.dker(f_rand, P), phi)))
    # closed form
    rhs = a.add(SMv, a.scale(a.dker(SMv, P), 1, h=1), a.vf(W, N), a.scale(a.dker(W, phi), Fr(1, 2), h=1),
                a.scale(bracket_phi(a, W, W, phi), Fr(1, 4)))
    record("T3.3 D1(a): W'_J(S_M) = S_M + hbar d_P S_M + D_{P*M} W + (hbar/2) d_phi W + 1/4 {W,W}_phi  [%s]" % label,
           a.eq(X_lhs, rhs))
    rhs2 = a.add(SMv, a.scale(a.dker(SMv, P), 1, h=1), a.vf(W, N), a.scale(a.mul(Zinv, a.dker(Z, phi)), Fr(1, 2), h=2))
    record("T3.4 same with (hbar^2/2) Z^{-1} d_phi Z  [%s]" % label, a.eq(X_lhs, rhs2))
    # controls: wrong signs / normalisations must fail
    bad1 = a.add(SMv, a.scale(a.dker(SMv, P), 1, h=1), a.scale(a.vf(W, N), -1), a.scale(a.dker(W, phi), Fr(1, 2), h=1),
                 a.scale(bracket_phi(a, W, W, phi), Fr(1, 4)))
    bad2 = a.add(SMv, a.scale(a.dker(SMv, P), 1, h=1), a.vf(W, N), a.scale(a.dker(W, phi), 1, h=1),
                 a.scale(bracket_phi(a, W, W, phi), Fr(1, 2)))
    record("T3.5 controls fail: -D_{P*M}W, and phi-terms without the 1/2", (not a.eq(X_lhs, bad1)) and (not a.eq(X_lhs, bad2)))
    # kernel of phi vs P* M P*
    phis = conv(phi)
    PMP = mm(mm(conv(P), Mv), conv(P))
    rel = None
    for c in [Fr(2), Fr(1), Fr(-1), Fr(-2), Fr(0)]:
        if phis == msc(PMP, c):
            rel = c
    record("T3.6 phi* = %s * (P* M P*)  [%s]" % (rel, label), rel is not None, "relation found: phi* = %s P*MP*" % rel)
    return rel, a.dker(SMv, P)

rel_gen, dPS_gen = check_a(Pgen, "generic P")
rel_n4, dPS_n4 = check_a(Pn4, "N=4 slot pattern P")
record("T3.7 constant hbar d_P S_M = 0 for the N=4 slot pattern (no field-antifield components in P)", not dPS_n4,
       "generic P gives %s" % ("0" if not dPS_gen else "nonzero"))

# ===================== T4: D1(b) =====================
hq = rand_poly(a, 2, 6, 0)
SMb = bracket(a, S0, hq)            # odd quadratic, Q-closed
Mb = op_from_hamiltonian(a, SMb)
record("T4.0 M_b: even, preserves <,>, [Q,M_b] = 0, D_Q S_Mb = 0",
       mpar(Mb) == 0 and not ker_sym_action(Mb, K0) and iszero(gcomm(Qm, Mb)) and not a.vf(SMb, Qm))
for P, label in [(Pgen, "generic P"), (Pn4, "N=4 slot pattern P")]:
    K = ker_sym_action(Qm, P)                          # K := Q P   (= K_eps in N=4)
    record("T4.1 (QP)* = QP* + P*Q  [%s]" % label, conv(K) == madd(mm(Qm, conv(P)), mm(conv(P), Qm)))
    N = mm(conv(P), Mb)
    record("T4.2 [Q, P*M] = (QP)* M  (= e^{-eps H} M)  [%s]" % label, gcomm(Qm, N) == mm(conv(K), Mb))
    phi = ker_sym_action(N, P)
    Qphi = ker_sym_action(Qm, phi)
    pred = ker_add(ker_sym_action(gcomm(Qm, N), P), ker_sym_action(N, K), coeffs=[1, -1])
    record("T4.3 (Q(x)1+1(x)Q)phi = ([Q,N](x)1+1(x)[Q,N])P - (N(x)1+1(x)N) K  [%s]" % label, Qphi == pred or not ker_add(Qphi, pred, coeffs=[1, -1]))
    J = make_J(a)
    Z, Zinv, W = flow(a, P, J)
    X = Wprime(a, P, J, Z, Zinv, SMb)
    QX = a.vf(X, Qm)
    QW = a.vf(W, Qm)
    O = a.add(a.vf(J, Qm), a.scale(bracket(a, J, J, K), Fr(1, 2)), a.scale(Delta(a, J, K), 1, h=1))
    record("T4.4 Q W(P,J) = W'_J(Q J + 1/2{J,J}_K + hbar Delta_K J)  (UV-boundary identity)  [%s]" % label,
           a.eq(QW, Wprime(a, P, J, Z, Zinv, O)))
    eJ = a.exp(a.scale(J, 1, h=-1))
    two = a.mul(Zinv, a.exp_dker(a.mul(a.mul(SMb, O), eJ), P))
    WO = Wprime(a, P, J, Z, Zinv, O)
    Cconn = a.scale(a.sub(two, a.mul(X, WO)), 1, h=-1)
    record("T4.5 connected two-insertion term C = hbar^{-1}[Z^{-1}e^{hbar d_P}(S_M O e^{J/hbar}) - W'(S_M) W'(O)] has no negative hbar power [%s]" % label,
           (a.min_hbar(Cconn) is None) or a.min_hbar(Cconn) >= 0)
    bSJ = bracket(a, SMb, J, K)
    DelS = Delta(a, SMb, K)
    rhs_b1 = a.add(Wprime(a, P, J, Z, Zinv, bSJ), a.scale(DelS, 1, h=1), a.scale(Cconn, -1))
    record("T4.6 D1(b) UV-boundary form: Q W'_J(S_M) = W'_J({S_M,J}_K) + hbar Delta_K S_M - C  [%s]" % label, a.eq(QX, rhs_b1))
    record("T4.7 {S_M, F}_K = D_{K* M} F  [%s]" % label, a.eq(bracket(a, SMb, f_rand, K), a.vf(f_rand, mm(conv(K), Mb))))
    # term-by-term form from D1(a)
    F = a.mul(Zinv, a.dker(Z, phi))
    QF = a.add(a.mul(Zinv, a.dker(Z, Qphi)), a.scale(a.mul(QW, F), -1, h=-1),
               a.scale(a.mul(Zinv, a.dker(a.mul(QW, Z), phi)), -1, h=-1))
    record("T4.8 Q(Z^{-1} d_phi Z) = Z^{-1} d_{Q phi} Z - hbar^{-1}(QW) Z^{-1}d_phi Z - hbar^{-1} Z^{-1} d_phi((QW)Z)  [%s]" % label,
           a.eq(a.vf(F, Qm), QF))
    rhs_b2 = a.add(a.vf(W, mm(conv(K), Mb)), a.scale(a.vf(QW, N), -1), a.scale(QF, Fr(1, 2), h=2))
    record("T4.9 D1(b) term-by-term: Q W'_J(S_M) = D_{K*M} W - D_{P*M}(QW) + (hbar^2/2) Q(Z^{-1}d_phi Z)  [%s]" % label, a.eq(QX, rhs_b2))
    record("T4.10 Delta_K S_M is a constant (field-independent)  [%s]" % label,
           all(sum(k[i] for i in range(a.NV)) == 0 for k in DelS))
    bad = a.add(Wprime(a, P, J, Z, Zinv, bSJ), a.scale(DelS, -1, h=1), a.scale(Cconn, -1))
    record("T4.11 control: the sign of hbar Delta_K S_M matters (generic data)  [%s]" % label, (not DelS) or not a.eq(QX, bad))

# ===================== T5: D1(c) =====================
J = make_J(a)
for P, label in [(Pgen, "generic P"), (Pn4, "N=4 slot pattern P")]:
    Z, Zinv, W = flow(a, P, J)
    phiM = ker_sym_action(Mm, P)
    lhs = a.vf(W, Mm)
    rhs = a.sub(Wprime(a, P, J, Z, Zinv, a.vf(J, Mm)), a.scale(a.mul(Zinv, a.dker(Z, phiM)), 1, h=2))
    record("T5.1 D_M W(P,J) = W'_J(D_M J) - d/ds W(P + s (M(x)1+1(x)M)P, J)|_0  [%s]" % label, a.eq(lhs, rhs))
    record("T5.2 D_M W = {S_M, W} (bracket at scale 0) [%s]" % label, a.eq(lhs, bracket(a, SM, W)))

# scaling family: M diagonal with integer weights, P(eps) built from eigen-kernels so that (M(x)1+1(x)M)P = -4 eps dP/deps
wts = {1: 2, 2: -2, 3: -2, 4: 2, 5: 4, 6: -4, 7: 2, 8: 6, 9: -2, 10: -6}
Md = zmat()
for i, m in wts.items():
    Md[i][i] = Fr(m)
record("T5.3 diagonal M preserves <,>", not ker_sym_action(Md, K0))
Peps = {}
rng2 = random.Random(7)
for i in FIELDS:
    for j in FIELDS:
        if i <= j and PAR[i] == PAR[j] and (wts[i] + wts[j]) % 4 == 0 and not (i == j and PAR[i]):
            if rng2.random() < 0.8:
                v = Fr(rng2.choice([-3, -2, -1, 1, 2, 3]), rng2.choice([1, 2]))
                e = -(wts[i] + wts[j]) // 4
                Peps[(i, j)] = (v, e)
                if i != j:
                    Peps[(j, i)] = ((-1) ** (PAR[i] * PAR[j]) * v, e)

def dker_eps(alg, p, phi):
    out = []
    for (i, j), (c, e) in phi.items():
        out.append(alg.scale(alg.d(alg.d(p, i), j), Fr(c, 2), e=e))
    return alg.add(*out)

def exp_dker_eps(alg, p, phi):
    r = dict(p)
    term = dict(p)
    n = 0
    while term:
        n += 1
        term = alg.scale(dker_eps(alg, term, phi), Fr(1, n), h=1)
        r = alg.add(r, term)
    return r

# check (M(x)1+1(x)M)P(eps) = -4 eps dP/deps
lhsk = {}
for (i, j), (c, e) in Peps.items():
    lhsk[(i, j)] = (c * (wts[i] + wts[j]), e)
okk = all(lhsk[k][0] == -4 * Peps[k][1] * Peps[k][0] for k in Peps)
record("T5.4 toy family: (M(x)1+1(x)M)P(eps) = -4 eps dP/deps", okk)
eJ = a.exp(a.scale(J, 1, h=-1))
Z = exp_dker_eps(a, eJ, Peps)
Y = a.sub(Z, a.const(1))
W = a.scale(a.log1p(Y), 1, h=1)
Zinv = a.inv1p(Y)
lhs = a.vf(W, Md)
WpDMJ = a.mul(Zinv, exp_dker_eps(a, a.mul(a.vf(J, Md), eJ), Peps))
rhs = a.add(WpDMJ, a.scale(a.eps_deriv(W), 4))
record("T5.5 D1(c): D_M W(P(eps),J) = W'_J(D_M J) + 4 eps d/deps|_J W(P(eps),J)", a.eq(lhs, rhs))
rhs_bad = a.add(WpDMJ, a.scale(a.eps_deriv(W), -4))
record("T5.6 control: with -4 eps d/deps it fails", not a.eq(lhs, rhs_bad))

json.dump(REPORT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "bv_toy_results.json"), "w"), indent=1)
nf = sum(1 for v in REPORT["checks"].values() if not v["pass"])
print("failures:", nf)

# ===================== T6: further forms used in reports =====================
REPORT2 = {}
def record2(name, ok, note=""):
    record(name, ok, note)

for P, label in [(Pgen, "generic P"), (Pn4, "N=4 slot pattern P")]:
    J = make_J(a)
    Z, Zinv, W = flow(a, P, J)
    N = mm(conv(P), Mm)
    phi = ker_sym_action(N, P)
    kappa = {k_: v / 2 for k_, v in phi.items()}
    record2("T6.1 kappa := kernel of P*MP* is phi/2 and graded symmetric  [%s]" % label,
            conv(kappa) == mm(mm(conv(P), Mm), conv(P)) and graded_symmetric(kappa))
    X = Wprime(a, P, J, Z, Zinv, SM)
    rhs = a.add(SM, a.vf(W, N), a.scale(a.dker(W, kappa), 1, h=1), a.scale(bracket_phi(a, W, W, kappa), Fr(1, 2)))
    record2("T6.2 D1(a) final form: W'_J(S_M) = S_M + D_{P*M}W + hbar d_kappa W + 1/2 {W,W}_kappa  [%s]" % label, a.eq(X, rhs))
    record2("T6.3 hbar d_P S_M = 0 for every even P (S_M odd quadratic pairs opposite parities)  [%s]" % label, not a.dker(SM, P))

# equivariant defect along the xi direction (M_b, Q S_Mb = 0)
for P, label in [(Pgen, "generic P"), (Pn4, "N=4 slot pattern P")]:
    K = ker_sym_action(Qm, P)
    J = make_J(a)
    def Odef(JJ):
        return a.add(a.vf(JJ, Qm), a.scale(bracket(a, JJ, JJ, K), Fr(1, 2)), a.scale(Delta(a, JJ, K), 1, h=1))
    Jx = a.add(J, a.mul(a.var(0), SMb))
    Ox = Odef(Jx)
    pred = a.sub(Odef(J), a.mul(a.var(0), a.add(bracket(a, SMb, J, K), a.scale(Delta(a, SMb, K), 1, h=1))))
    record2("T6.4 O_K(J + xi S_M) = O_K(J) - xi({S_M,J}_K + hbar Delta_K S_M)   (xi left; Q S_M = 0)  [%s]" % label, a.eq(Ox, pred))
    Z, Zinv, W = flow(a, P, J)
    Zx, Zxinv, Wx = flow(a, P, Jx)
    X = Wprime(a, P, J, Z, Zinv, SMb)
    record2("T6.5 Q W(P,J+xi S_M) = Q W - xi Q W'_J(S_M)  [%s]" % label, a.eq(a.vf(Wx, Qm), a.sub(a.vf(W, Qm), a.mul(a.var(0), a.vf(X, Qm)))))
    # right-derivative control: the right xi-derivative of W(P,J+xi S_M) is -W'_J(S_M)
    WxR = {}
    for k_, v in Wx.items():
        if k_[0]:
            # move xi to the right: xi * m = (-1)^{|m|} m * xi; right derivative picks m with that sign
            pm = sum(k_[i] for i in range(1, a.NV) if PAR[i]) % 2
            kk = list(k_); kk[0] = 0
            WxR[tuple(kk)] = WxR.get(tuple(kk), 0) + v * (-1) ** pm
    WxR = {k_: v for k_, v in WxR.items() if v}
    record2("T6.6 right xi-derivative gives -W'_J(S_M) (the left derivative is the one with leading term +S_M)  [%s]" % label,
            a.eq(WxR, a.scale(X, -1)))
    O = Odef(J)
    eJ = a.exp(a.scale(J, 1, h=-1))
    two = a.mul(Zinv, a.exp_dker(a.mul(a.mul(SMb, O), eJ), P))
    C = a.scale(a.sub(two, a.mul(X, Wprime(a, P, J, Z, Zinv, O))), 1, h=-1)
    record2("T6.7 the two-insertion term W''_J(S_M, O_eps(J)) is nonzero for generic J (cannot be dropped unless O_eps(J)=0)  [%s]" % label, bool(C))

json.dump(REPORT, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "bv_toy_results.json"), "w"), indent=1)
nf = sum(1 for v in REPORT["checks"].values() if not v["pass"])
print("failures (final):", nf)
