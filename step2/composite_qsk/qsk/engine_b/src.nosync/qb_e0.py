"""Task E0 (route B): conventions in code, checked exactly on plane waves.

(a) arrows of xi. (Lemma qs-arrows) and the commutators of Prop qs-comm; operator <-> functional through the
    relation q(k) = Pi K_0 H_free(-k,k) verified on the free action (Prop c1-Q, eq c1-free): a quadratic functional F
    with F = (1/2)<e, X e> has X(k) = Pi K_0 H_F(-k,k);
(b) the composite line along xi, 2 (P* xi. P* (x) 1) K_0, as a multiple of xi e^{-2 eps p^2}/p^2;
(d) canonical form (definition, recorded); (e) twisted coefficient (definition, recorded); (f) colour of T5 (engine).
Position-dependent operators (eta., v.) are handled as normal-ordered polynomials sum_a x^a A_a(k) (A_a applied first).
"""
import json, time, itertools, math
import sympy as sp
import qb_core as B

C = B.C; V = B.V
I = sp.I
NT = B.NT; OFF = B.OFF
k = list(sp.symbols('k0:4'))
mk = [-x for x in k]
INT6 = C.INT6
PAIRS = [(a, b) for a in INT6 for b in INT6]
K0 = C.k0_matrix()
PI = C.PI
q = C.q_matrix(k)
qgf = C.qgf_matrix(k)
k2 = sum(x * x for x in k)


def hess(terms, kr, kc):
    bl = V.hessian_blocks(terms, [], list(kr), list(kc), PAIRS)
    M = sp.zeros(NT, NT)
    for (a, b), d in bl.items():
        for key, Mx in d.items():
            for i in range(Mx.rows):
                for j in range(Mx.cols):
                    if Mx[i, j] != 0:
                        M[OFF[a] + i, OFF[b] + j] += Mx[i, j]
    return M.applyfunc(sp.expand)


def op_from_functional(terms):
    return (PI * K0 * hess(terms, mk, k)).applyfunc(sp.expand)


def is_zero(M):
    return all(sp.expand(e) == 0 for e in M)


def blocks(M):
    """set of (target slot, source slot) with a nonzero entry."""
    return sorted(set((B.slot_of(i), B.slot_of(j)) for i in range(NT) for j in range(NT) if M[i, j] != 0))


def graded_sym(E, Em):
    return all(sp.expand(E[i, j] - (-1) ** (B.PARITY[i] * B.PARITY[j]) * Em[j, i]) == 0 for i in range(NT) for j in range(NT))


def homogeneous_degree(M, d):
    for e in M:
        if e != 0:
            P = sp.Poly(e, *k)
            if any(sum(mon) != d for mon in P.monoms()):
                return False
    return True


# ------------------------------------------------------------------ polynomial-coefficient operators
class PolyOp:
    """sum_a x^a A_a(k): apply the constant-coefficient operator with symbol A_a, then multiply by x^a."""

    def __init__(self, d=None):
        self.d = {e: M for e, M in (d or {}).items()}

    @staticmethod
    def const(M):
        return PolyOp({(0, 0, 0, 0): M})

    @staticmethod
    def mult(poly):
        """multiplication by a polynomial in x (sympy expression in xs)."""
        P = sp.Poly(sp.expand(poly), *XS)
        return PolyOp({tuple(m): c * sp.eye(NT) for m, c in P.terms()})

    def __add__(self, o):
        d = dict(self.d)
        for e, M in o.d.items():
            d[e] = (d[e] + M) if e in d else M
        return PolyOp(d)

    def scale(self, c):
        return PolyOp({e: (c * M) for e, M in self.d.items()})

    def __sub__(self, o):
        return self + o.scale(-1)

    def __mul__(self, o):
        out = {}
        for a, A in self.d.items():
            for b, Bm in o.d.items():
                for g in itertools.product(*[range(x + 1) for x in b]):
                    cf = 1
                    for bi, gi in zip(b, g):
                        cf *= math.comb(bi, gi)
                    DA = A
                    for mu, gm in enumerate(g):
                        for _ in range(gm):
                            DA = DA.applyfunc(lambda e, mu=mu: -I * sp.diff(e, k[mu]))
                    if DA.is_zero_matrix:
                        continue
                    e = tuple(a[i] + b[i] - g[i] for i in range(4))
                    term = (cf * DA * Bm)
                    out[e] = (out[e] + term) if e in out else term
        return PolyOp({e: M.applyfunc(sp.expand) for e, M in out.items()})

    def is_zero(self):
        return all(sp.expand(e) == 0 for M in self.d.values() for e in M)

    def gcomm(self, o, pa, pb):
        return self * o - (o * self).scale((-1) ** (pa * pb))


XS = list(sp.symbols('x0:4'))


def spin_gen_Sigma(mu, nu):
    return sp.Rational(1, 2) * (B.Rbm[mu] * B.Rm[nu]) if mu != nu else sp.zeros(16, 16)


def spin_gen_Sigmas(mu, nu):
    return sp.Rational(1, 2) * (B.Rm[mu] * B.Rbm[nu]) if mu != nu else sp.zeros(16, 16)


def rot_gen(mu, nu):
    """Sigma^{mu nu} on the 54 slot components: vector rule on A_V and A*_V, spinor rules on lam (Sigma) and lams
    (Sigma^*), zero on A_perp, A*_perp, c, c*."""
    M = sp.zeros(NT, NT)
    for base in (OFF['A'], OFF['As']):
        # (Sigma A)_m = delta_{m mu} A_nu - delta_{m nu} A_mu
        M[base + mu, base + nu] += 1
        M[base + nu, base + mu] -= 1
    S1 = spin_gen_Sigma(mu, nu); S2 = spin_gen_Sigmas(mu, nu)
    for i in range(16):
        for j in range(16):
            M[OFF['lam'] + i, OFF['lam'] + j] += S1[i, j]
            M[OFF['lams'] + i, OFF['lams'] + j] += S2[i, j]
    return M


WEIGHT = {'A': 1, 'lam': sp.Rational(3, 2), 'c': 0, 'As': 3, 'lams': sp.Rational(5, 2), 'cs': 4}
WDIAG = sp.diag(*[WEIGHT[B.SPEC_OF[i]] for i in range(NT)])


def v_operator(v, omega, divv):
    """v. = v^mu d_mu + (1/2) omega_{mu nu} Sigma^{mu nu} + (w/4) div v (Prop qs-comm (iii)); v, omega, divv polynomials
    in x."""
    op = PolyOp({})
    for mu in range(4):
        op = op + PolyOp.mult(v[mu]) * PolyOp.const(I * k[mu] * sp.eye(NT))
    for mu in range(4):
        for nu in range(4):
            if omega[mu][nu] != 0:
                op = op + PolyOp.mult(sp.Rational(1, 2) * omega[mu][nu]) * PolyOp.const(rot_gen(mu, nu))
    op = op + PolyOp.mult(divv) * PolyOp.const(WDIAG / 4)
    return op


def run():
    t0 = time.time()
    res = {}
    gates = []

    def gate(gid, ok, value, target, note=''):
        gates.append(dict(id=gid, value=value, target=target, pass_=bool(ok), note=note))

    # ---------------- conversion functional <-> operator on the free action
    Xfree = op_from_functional(['kinA', 'kinlam', 'kingh'])
    ok = is_zero(Xfree - q)
    gate('E0-conv', ok, 'Pi K_0 H_free(-k,k) - q(k) = 0' if ok else 'mismatch', 'q(k) = Pi K_0 H_free(-k,k)',
         'free action of eq c1-free (kinA, kinlam, kingh of vertices.py); fixes the reading F = (1/2)<e,Xe>')

    # ---------------- chirality / orientation
    E = sp.eye(16)
    g = B.Rbm[1] * B.Rm[2] * B.Rbm[3] * B.Rm[0]          # rho(e1 e2 e3 e4) on Sigma, engine index 0 = m = 4
    ev_P = all(sp.simplify((g * E[:, j] - E[:, j]).norm()) == 0 for j in range(8))
    ev_M = all(sp.simplify((g * E[:, j] + E[:, j]).norm()) == 0 for j in range(8, 16))
    gate('E0-orient', ev_P and ev_M, f'rho(e1e2e3e4) = +1 on comps 0..7: {ev_P}, -1 on 8..15: {ev_M}',
         'lamP = S_+ (x) W = +1 eigenspace for the orientation dx1dx2dx3dx4 of the dictionary (engine m=0 <-> 4)')

    # ---------------- xi. for xi in S_+ (x) W (symbolic) and in S_- (x) W^*
    xiP = B.XI[:8] + [0] * 8
    xiM = [0] * 8 + B.XI[8:]
    out_xi = {}
    for lab, xi in (('S_+ (x) W', xiP), ('S_- (x) W^*', xiM)):
        quad, cub = B.register_xi_terms(xi, tag='e0' + ('p' if lab.startswith('S_+') else 'm'))
        Xs = op_from_functional(quad)
        Xf = B.xi_op_fields(k, xi)
        # field arrows of the functional-derived operator must be those of Lemma qs-arrows
        field_rows = [i for i in range(NT) if B.SPEC_OF[i] in ('A', 'lam', 'c')]
        Xs_f = sp.zeros(NT, NT)
        for i in field_rows:
            for j in range(NT):
                if B.SPEC_OF[j] in ('A', 'lam', 'c'):
                    Xs_f[i, j] = Xs[i, j]
        ok_arrows = is_zero(Xs_f + Xf)          # functional-derived xi. = minus the arrows printed in Lemma qs-arrows
        # antifield arrows from preservation of the pairing
        Xcomp, _ = B.complete_antifield_arrows(-Xf, k, [('lams', 'As'), ('As', 'lams')], name='u' + lab[2])
        ok_af = is_zero(Xcomp - Xs)
        Xm = Xs.subs({k[m]: mk[m] for m in range(4)}, simultaneous=True)
        pres = is_zero(Xs * K0 + PI * K0 * Xm.T)
        comQ = is_zero(Xs * q + q * Xs)
        comH = is_zero(Xs * (k2 * sp.eye(NT)) - (k2 * sp.eye(NT)) * Xs)
        Cxi = (Xs * qgf + qgf * Xs).applyfunc(sp.expand)
        bl = blocks(Cxi)
        hom = homogeneous_degree(Cxi, 1)
        # lam -> c entry: -i (p, Gamma(xi, lam)|_V)
        lc_ok = True
        for j in range(16):
            lam = [1 if i == j else 0 for i in range(16)]
            G = B.Gam(xi, lam)
            target = sp.expand(I * sum(k[m] * G[m] for m in range(4)))
            if sp.expand(Cxi[OFF['c'], OFF['lam'] + j] - target) != 0:
                lc_ok = False
        Cm = Cxi.subs({k[m]: mk[m] for m in range(4)}, simultaneous=True)
        E1 = (Cxi * K0).applyfunc(sp.expand); E1m = (Cm * K0).applyfunc(sp.expand)
        gsym = graded_sym(E1, E1m)
        # composite line 2 (P* xi. P* (x) 1) K_0 = 2 qgf X qgf k0 e^{-2 eps k^2}/k^4
        Xline = (2 * qgf * Xs * qgf * K0).applyfunc(sp.expand)
        lbl = blocks(Xline)
        # c (x) lam entries / k^2
        mult = None
        cl = [sp.cancel(Xline[OFF['c'], OFF['lam'] + j] / k2) for j in range(16)]
        lcv = [sp.cancel(Xline[OFF['lam'] + j, OFF['c']] / k2) for j in range(16)]
        # compare with xi components (as an element of the lam fibre): entry = mu * xi_j
        ratios = set()
        for j in range(16):
            if xi[j] != 0:
                ratios.add(sp.simplify(cl[j] / xi[j]))
            elif cl[j] != 0:
                ratios.add('nonzero off-chirality')
        ratios2 = set()
        for j in range(16):
            if xi[j] != 0:
                ratios2.add(sp.simplify(lcv[j] / xi[j]))
        # identities of Cor qs-edges (iii) on symbols: P* xi. P* = C_xi Q^GF g^2 - xi. (Q^GF)^2 g^2, [Q, P* xi. P*] = E C_xi g
        lhs = (qgf * Xs * qgf).applyfunc(sp.expand)
        rhs = (Cxi * qgf - Xs * qgf * qgf).applyfunc(sp.expand)
        iii_a = is_zero(lhs - rhs)
        # [Q, qgf X qgf] (graded: both odd -> commutator of odd Q with even composite): Q M - M Q, should = k^2 C_xi
        iii_b = is_zero(q * lhs + lhs * q - k2 * Cxi)
        out_xi[lab] = dict(arrows_from_functional=ok_arrows, antifield_arrows_from_pairing=ok_af,
                           preserves_pairing=pres, comm_Q=comQ, comm_H=comH, C_xi_blocks=[f'{s}<-{t}' for (s, t) in bl],
                           C_xi_order_one=hom, C_xi_lam_to_c=lc_ok, C_xi_K0_graded_symmetric=gsym,
                           line_blocks=[f'{s}<-{t}' for (s, t) in lbl], line_c_lam_over_k2_per_xi=sorted(map(str, ratios)),
                           line_lam_c_over_k2_per_xi=sorted(map(str, ratios2)), cor_edges_iii_split=iii_a,
                           cor_edges_iii_comm=iii_b)
        tagc = 'p' if lab.startswith('S_+') else 'm'
        gate(f'E0-xi-arrows-{tagc}', ok_arrows and ok_af, dict(field_arrows_equal_minus_Lemma=ok_arrows, antifield_arrows_by_pairing=ok_af),
             'xi. induced by S_2(xi) in the convention of Prop c1-Q / Prop c7-susy (F = (1/2)<e,Xe>): lam -> +Gamma(xi,lam), A -> +(1/2) rho(dA) xi; antifield arrows from the pairing',
             'FINDING: Lemma qs-arrows prints lam -> -Gamma(xi,lam), A -> -(1/2) rho(dA) xi, i.e. the reading int(Xe,e*) <-> X, which is minus the convention of Prop c1-Q (int(dc,A*) <-> Q c = -dc); all statements of Prop qs-comm (i) hold for both signs, C_xi and the composite line change sign')
        gate(f'E0-xi-comm-{tagc}', pres and comQ and comH, f'pairing {pres}, [xi.,Q]=0 {comQ}, [xi.,H]=0 {comH}',
             'Prop qs-comm (i)')
        want = {('c', 'lamP'), ('c', 'lamM'), ('lamP', 'AsV'), ('lamM', 'AsV'), ('lamP', 'AsW'), ('lamM', 'AsW'),
                ('AV', 'lamsP'), ('AV', 'lamsM'), ('AW', 'lamsP'), ('AW', 'lamsM'), ('lamsP', 'cs'), ('lamsM', 'cs')}
        okC = set(bl) <= want and hom and lc_ok and gsym and not any(t in ('AV', 'AW', 'c') for (s, t) in bl)
        gate(f'E0-Cxi-{tagc}', okC, dict(blocks=out_xi[lab]['C_xi_blocks'], order_one=hom, lam_to_c=lc_ok, graded_sym=gsym),
             'components lam->c (= +i(p,Gamma(xi,lam)|_V) with xi. of Prop c7-susy), A*->lam, lam*->A, c*->lam*, order one, zero on A, c; (C_xi (x) 1)K_0 graded symmetric',
             'refinement: only A*|_V -> lam and lam* -> A|_V occur (the A*|_perp -> lam and lam* -> A|_perp components cancel between xi.Q^GF and Q^GF xi.)')
        gate(f'E0-line-{tagc}', iii_a and iii_b and set(lbl) <= {('c', 'lamP'), ('c', 'lamM'), ('lamP', 'c'), ('lamM', 'c')}
             and len(ratios) == 1 and 'nonzero off-chirality' not in ratios,
             dict(blocks=out_xi[lab]['line_blocks'], c_lam=sorted(map(str, ratios)), lam_c=sorted(map(str, ratios2))),
             'Cor qs-edges (ii),(iii): composite line in c (x) lam, lam (x) c; multiple of xi e^{-2eps p^2}/p^2')
    res['xi'] = out_xi

    # ---------------- eta: eta. = x^mu xi_mu. + eta^0., xi_mu = rho(e_mu) eta
    out_eta = {}
    for lab, eta in (('S_+ (x) W^*', [0] * 8 + B.ETA[8:]), ('S_- (x) W', B.ETA[:8] + [0] * 8)):
        tg = 'ep' if lab.startswith('S_+') else 'em'
        ev = sp.Matrix(eta)
        xis = [list(B.Rbm[mu] * ev) for mu in range(4)]
        Xmu = []
        for mu in range(4):
            quad, cub = B.register_xi_terms(xis[mu], tag=f'{tg}{mu}')
            Xmu.append(op_from_functional(quad))
        s0 = B.register_s0_term(eta, tag=tg)
        X0 = op_from_functional(s0)
        bl0 = blocks(X0)
        # pairing-preserving completion: (x^mu X_mu (x) 1 + 1 (x) x^mu X_mu) K_0 = -i Pi sum_mu d_{k_mu}(K_0 X_mu(-k)^T);
        # the operator induced by S_2(eta) has the field arrows x^mu X_mu (normal ordered) and antifield arrows fixed by
        # the pairing: add Y (constant, antifield rows) with Y K_0 + Pi K_0 Y(-k)^T = +i Pi sum_mu d_mu(K_0 X_mu(-k)^T)
        RHS = sp.zeros(NT, NT)
        for mu in range(4):
            Xm_ = Xmu[mu].subs({k[m]: mk[m] for m in range(4)}, simultaneous=True)
            RHS += (K0 * Xm_.T).applyfunc(lambda e, mu=mu: sp.expand(sp.diff(e, k[mu])))
        RHS = (I * PI * RHS).applyfunc(sp.expand)
        afrows = [i for i in range(NT) if B.SPEC_OF[i] in ('As', 'lams', 'cs')]
        Rrows = sp.zeros(NT, NT)
        for i in afrows:
            for j in range(NT):
                Rrows[i, j] = RHS[i, j]
        Y = (Rrows * K0.inv()).applyfunc(sp.expand)
        Ym = Y.subs({k[m]: mk[m] for m in range(4)}, simultaneous=True)
        okY = is_zero(Y * K0 + PI * K0 * Ym.T - RHS)
        Yblocks = [f'{s_}<-{t_}' for (s_, t_) in blocks(Y)]
        lhs = sp.zeros(NT, NT)
        for mu in range(4):
            lhs += B.symbol_derivative(q, k, mu) * Xmu[mu]
        lhs_naive = (lhs + X0 * q + q * X0).applyfunc(sp.expand)
        okQ_naive = is_zero(lhs_naive)
        lhs += (X0 + Y) * q + q * (X0 + Y)
        okQ = is_zero(lhs)
        # without s0 the identity fails (s0 is needed)
        lhs2 = sp.zeros(NT, NT)
        for mu in range(4):
            lhs2 += B.symbol_derivative(q, k, mu) * Xmu[mu]
        needs_s0 = not is_zero(lhs2)
        # [eta., H] = 2 sum xi_mu. d_mu : -sigma^mu(H) X_mu = 2 i k_mu X_mu
        H = k2 * sp.eye(NT)
        lhsH = sp.zeros(NT, NT)
        for mu in range(4):
            lhsH += -B.symbol_derivative(H, k, mu) * Xmu[mu]
        lhsH += (X0 + Y) * H - H * (X0 + Y)
        rhsH = sum((2 * I * k[mu] * Xmu[mu] for mu in range(4)), sp.zeros(NT, NT))
        okH = is_zero(lhsH - rhsH)
        C0 = sp.zeros(NT, NT)
        for mu in range(4):
            C0 += B.symbol_derivative(qgf, k, mu) * Xmu[mu]
        C0 += (X0 + Y) * qgf + qgf * (X0 + Y)
        C0 = C0.applyfunc(sp.expand)
        out_eta[lab] = dict(s0_blocks=[f'{s}<-{t}' for (s, t) in bl0], comm_Q_zero=okQ, comm_Q_zero_without_Y=okQ_naive,
                            Y_blocks=Yblocks, Y_solves_pairing=okY, s0_needed=needs_s0,
                            comm_H=okH, C0_blocks=[f'{s}<-{t}' for (s, t) in blocks(C0)],
                            C0_order_le_one=all(sp.Poly(e, *k).total_degree() <= 1 for e in C0 if e != 0))
        gate(f'E0-eta-{tg}', okQ and okY and needs_s0 and okH and set(bl0) <= {('lamP', 'AW'), ('lamM', 'AW'), ('AsW', 'lamsP'), ('AsW', 'lamsM')},
             dict(comm_Q=okQ, comm_Q_without_Y=okQ_naive, Y_blocks=Yblocks, s0_needed=needs_s0, comm_H=okH, s0_blocks=out_eta[lab]['s0_blocks']),
             '[eta.,Q] = 0 and [eta.,H] = 2 xi_mu. d_mu for the operator induced by S_2(eta) = int x^mu S_2(xi_mu) + int s_0; eta^0. arrows A_perp -> lam, lam* -> A*_perp',
             'FINDING: eta. is x^mu xi_mu. (normal ordered) + eta^0. + Y, with Y the constant antifield arrows lam* -> A*, A* -> lam* produced by the integration by parts of x^mu on the antifield arrows (pairing preservation); with Y, [eta.,Q] = sum sigma^mu(Q) xi_mu. + [eta^0.+Y, Q] = 0 exactly; without Y the identity fails. The constant part C^0_eta of [eta.,Q^GF] includes Y.')
    res['eta'] = out_eta

    # ---------------- even generators: Delta (v = 2x) and K_b
    bsym = list(sp.symbols('b0:4'))
    x = XS
    bx = sum(bsym[m] * x[m] for m in range(4)); xx = sum(x[m] ** 2 for m in range(4))
    # spinor generators: equivariance and pairing invariance
    eq_ok = True; pair_ok = True
    for mu in range(4):
        for nu in range(4):
            if mu == nu:
                continue
            S1 = spin_gen_Sigma(mu, nu); S2 = spin_gen_Sigmas(mu, nu)
            for lam_ in range(10):
                vec = [0] * 10
                if lam_ == nu:
                    vec[mu] += 1
                if lam_ == mu:
                    vec[nu] -= 1
                rv = sum((vec[M] * B.Rm[M] for M in range(10)), sp.zeros(16, 16))
                if not (S2 * B.Rm[lam_] - B.Rm[lam_] * S1 - rv).is_zero_matrix:
                    eq_ok = False
            if not (S1.T * B.Mp + B.Mp * S2).is_zero_matrix:
                pair_ok = False
    gate('E0-spin', eq_ok and pair_ok, dict(equivariance=eq_ok, pairing=pair_ok),
         'Sigma^{mu nu} = (1/2) rho_mu rho_nu on Sigma and Sigma^*: rho(Sigma v) = [Sigma, rho(v)], (Sigma s, s*) + (s, Sigma s*) = 0')
    Qop = PolyOp.const(q); Gop = PolyOp.const(qgf); Hop = PolyOp.const(k2 * sp.eye(NT))
    # Delta: v = 2x, omega = 0, div v = 8
    vD = v_operator([2 * x[m] for m in range(4)], [[0] * 4 for _ in range(4)], 8)
    okD1 = vD.gcomm(Qop, 0, 1).is_zero()
    okD2 = (vD.gcomm(Gop, 0, 1) + Gop.scale(4)).is_zero()
    okD3 = (vD.gcomm(Hop, 0, 0) + Hop.scale(4)).is_zero()
    gate('E0-Delta', okD1 and okD2 and okD3, dict(comm_Q=okD1, comm_QGF_eq_m4QGF=okD2, comm_H_eq_m4H=okD3),
         '[v.,Q] = 0, [v.,Q^GF] = -4 Q^GF, [v.,H] = -4H for v = 2E (Prop c7-dil)')
    # K_b
    vK = [xx * bsym[m] - 2 * bx * x[m] for m in range(4)]
    om = [[2 * (x[m] * bsym[n] - x[n] * bsym[m]) for n in range(4)] for m in range(4)]
    vKop = v_operator(vK, om, -8 * bx)
    okK1 = vKop.gcomm(Qop, 0, 1).is_zero()
    comG = vKop.gcomm(Gop, 0, 1)
    rest = comG - (PolyOp.mult(4 * bx) * Gop)
    const_only = all(sum(e) == 0 for e, M in rest.d.items() if not M.applyfunc(sp.expand).is_zero_matrix)
    CK = rest.d.get((0, 0, 0, 0), sp.zeros(NT, NT)).applyfunc(sp.expand)
    CK_order0 = all(sp.Poly(e, *k).total_degree() == 0 for e in CK if e != 0) if not CK.is_zero_matrix else True
    CK_AV = all(sp.expand(CK[OFF['c'], OFF['A'] + m] + 4 * bsym[m]) == 0 for m in range(4))
    CK_AsW = all(CK[i, j] == 0 for j in range(OFF['As'] + 4, OFF['As'] + 10) for i in range(NT))
    # C_K = 2 sigma^rho(Q^GF)(w b_rho - b_nu Sigma^{rho nu})
    CKf = sp.zeros(NT, NT)
    for rho_ in range(4):
        sg = B.symbol_derivative(qgf, k, rho_)
        Mx = WDIAG * bsym[rho_]
        for nu in range(4):
            Mx = Mx - bsym[nu] * rot_gen(rho_, nu)
        CKf += 2 * sg * Mx
    CK_formula = is_zero(CK - CKf)
    comH = vKop.gcomm(Hop, 0, 0)
    rhsH = PolyOp.mult(4 * bx) * Hop
    Wm1 = (WDIAG - sp.eye(NT))
    Hextra = sp.zeros(NT, NT)
    for rho_ in range(4):
        Hextra += -4 * Wm1 * bsym[rho_] * I * k[rho_]
        for nu in range(4):
            Hextra += 4 * bsym[nu] * rot_gen(rho_, nu) * I * k[rho_]
    okKH = (comH - rhsH - PolyOp.const(Hextra)).is_zero()
    gate('E0-K', okK1 and const_only and CK_order0 and CK_AV and CK_AsW and CK_formula and okKH,
         dict(comm_Q=okK1, remainder_constant=const_only, C_K_order0=CK_order0, C_K_AV_eq_m4b=CK_AV, C_K_AsW_zero=CK_AsW,
              C_K_formula=CK_formula, comm_H=okKH, C_K_blocks=[f'{s}<-{t}' for (s, t) in blocks(CK)]),
         '[v.,Q]=0; [v.,Q^GF] = 4(b,x)Q^GF + C_K, C_K = 2 sigma^rho(Q^GF)(w b_rho - b_nu Sigma^{rho nu}), C_K A_V = -4(b,A_V), C_K = 0 on A*_perp; [v.,H] = 4(b,x)H - 4(w-1)(b,d) + 4 b_nu Sigma^{rho nu} d_rho')
    res['K_b'] = dict(C_K_blocks=[f'{s}<-{t}' for (s, t) in blocks(CK)])
    res['recorded'] = {
        'canonical_form': 'two-leg outputs on plane waves alpha delta + beta^mu d_mu delta + gamma^{mu nu} d_mu d_nu delta <-> int omega (x^mu x^nu K_munu + x^mu K_mu + K_0): K_0 = alpha, K_mu = i beta^mu, K_munu = -gamma^{munu}; route B computes position factors as -i d/dr (one) and -d^2/dr dr (two) of a momentum r injected at the S vertex and absorbed by omega',
        'twisted_coefficient': 'Q^xi_1 = L_1^(0) + Q_1 - (1/2) L_1^(2)(xi,xi) (Def c8-twisted, Prop qs-twist (i)); used as stated',
        'colour_T5': 'T5 = (1/4) sum_{MN} ([A_M,A_N], (rho_M rho_N xi, lam*)) carries the colour key f of vertices.py and acts on the cycle as ad_Y of its output Y, like a cubic vertex of I (Rem qs-colour)'}
    res['seconds'] = round(time.time() - t0, 1)
    return res, gates


if __name__ == '__main__':
    import sys
    res, gates = run()
    out = sys.argv[1] if len(sys.argv) > 1 else 'e0_conventions.json'
    json.dump(dict(results=res, gates=gates), open(out, 'w'), indent=1, default=str)
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL', str(g['value'])[:200])
    print('seconds', res['seconds'])
