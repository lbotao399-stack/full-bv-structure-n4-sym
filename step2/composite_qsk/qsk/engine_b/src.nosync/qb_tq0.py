"""Task TQ0 (route B): the quadratic Q-closed directions of degree -1 linear in the Killing data (Prop qs-freexi (i)).

Method: a quadratic functional F = (1/2)<e, X e> corresponds to the operator X = Pi K_0 H_F(-k,k) (E0-conv); F is
Q-closed iff the graded commutator [X, Q] vanishes (X odd: anticommutator).  Position-dependent structures
int x^mu F(xi_mu) are handled as the normal-ordered field arrows x^mu X^f_mu plus the constant antifield arrows Y fixed
by pairing preservation (E0-eta); then [O, Q] = x^mu [X_mu, Q] + sum_mu sigma^mu(Q) X_mu + [Y, Q] and both the x-linear
and the constant parts must vanish.
Along Sigma (xi in S_+ (x) W, symbolic): T1, T2, T3, T4 of Table qs-terms and T5' = int((d^*A|_V) xi, lam*).
Along Sigma^* (eta in S_+ (x) W^*, symbolic): the position structures x^mu T_i(rho(e_mu) eta), i = 1..5, and the two
intrinsic structures of weight 7/2: s_0 = -(rho(A|_perp) eta, lam*) and (rho(A|_V) eta, lam*).
Exactness: a Q-exact quadratic functional of degree -1 would be Q of a degree -2 functional of field weight 9/2
(resp. 7/2): the degree -2 quadratic contents (A*,A*), (A*,lam*), (lam*,lam*), (c*,A), (c*,lam) have weight >= 5;
none exists (counted below)."""
import json, time
import sympy as sp
import qb_core as B
import qb_e0 as E

C = B.C; V = B.V
I = sp.I
NT = B.NT; OFF = B.OFF
k = E.k; mk = E.mk
K0 = E.K0; PI = E.PI; q = E.q


def register_structs(xi, tag):
    xv = sp.Matrix(xi)
    MRxi = [(xv.T * B.MR[M]) for M in range(10)]

    def L_T12(Ms):
        def L(s):
            l, As = s
            u = sp.Matrix(l['u'])
            return sum(((MRxi[M] * u)[0] * As['u'][M] for M in Ms), 0)
        return L

    def L_T34(Ms):
        def L(s):
            a, ls = s
            a2 = dict(u=[a['u'][M] if M in Ms else 0 for M in range(10)], p=a['p'])
            w = B.rho_two(B._dA(a2), xi)
            return (w.T * B.Mp * sp.Matrix(ls['u']))[0]
        return L

    def L_T5p(s):
        a, ls = s
        dsA = -sum((I * a['p'][m] * a['u'][m] for m in range(4)), 0)      # d^* A|_V = -d_m A_m
        return dsA * (xv.T * B.Mp * sp.Matrix(ls['u']))[0]
    V4 = list(range(4)); P6 = list(range(4, 10))
    V.TERMS['T1_' + tag] = (sp.Integer(-1), ('lam', 'As'), L_T12(V4), ('k', (0, 1)))
    V.TERMS['T2_' + tag] = (sp.Integer(-1), ('lam', 'As'), L_T12(P6), ('k', (0, 1)))
    V.TERMS['T3_' + tag] = (sp.Rational(1, 2), ('A', 'lams'), L_T34(V4), ('k', (0, 1)))
    V.TERMS['T4_' + tag] = (sp.Rational(1, 2), ('A', 'lams'), L_T34(P6), ('k', (0, 1)))
    V.TERMS['T5p_' + tag] = (sp.Integer(1), ('A', 'lams'), L_T5p, ('k', (0, 1)))
    return [nm + '_' + tag for nm in ('T1', 'T2', 'T3', 'T4', 'T5p')]


def register_eta_intrinsic(eta, tag):
    ev = sp.Matrix(eta)

    def L_AV(s):
        a, ls = s
        w = sp.zeros(16, 1)
        for M in range(4):
            if a['u'][M] != 0:
                w += a['u'][M] * (B.Rbm[M] * ev)
        return (w.T * B.Mp * sp.Matrix(ls['u']))[0]
    V.TERMS['iAV_' + tag] = (sp.Integer(1), ('A', 'lams'), L_AV, ('k', (0, 1)))
    s0 = B.register_s0_term(eta, tag=tag)
    return s0 + ['iAV_' + tag]


def nullspace_of(maps, coeff_syms):
    """maps: list of matrices (one per structure); returns a basis of {a : sum a_i M_i = 0}."""
    a = sp.symbols(f'a0:{len(maps)}')
    tot = sp.zeros(NT, NT)
    for ai, M in zip(a, maps):
        tot += ai * M
    eqs = []
    for e in tot:
        e = sp.expand(e)
        if e != 0:
            syms = sorted(e.free_symbols - set(a), key=lambda s: s.name)
            if syms:
                eqs += sp.Poly(e, *syms).coeffs()
            else:
                eqs.append(e)
    Mx = sp.Matrix([[sp.expand(eq).coeff(ai) for ai in a] for eq in eqs]) if eqs else sp.zeros(1, len(a))
    ns = Mx.nullspace()
    return [list(v) for v in ns], Mx.rank()


def run():
    t0 = time.time()
    res = {}; gates = []
    # ---------------- along Sigma
    xi = B.XI[:8] + [0] * 8
    names = register_structs(xi, 'tq0')
    ops = [E.op_from_functional([nm]) for nm in names]
    comms = [(X * q + q * X).applyfunc(sp.expand) for X in ops]
    ns, rk = nullspace_of(comms, B.XI[:8])
    lab = ['T1', 'T2', 'T3', 'T4', "T5'"]
    basis = [{lab[i]: str(v[i]) for i in range(5) if v[i] != 0} for v in ns]
    # expected span{T1+T3, T2+T4}
    expect = [sp.Matrix([1, 0, 1, 0, 0]), sp.Matrix([0, 1, 0, 1, 0])]
    M = sp.Matrix.hstack(*[sp.Matrix(v) for v in ns]) if ns else sp.zeros(5, 0)
    ok = len(ns) == 2 and sp.Matrix.hstack(M, *expect).rank() == 2
    # sign check: T1 - T3 not closed
    res['Sigma'] = dict(structures=lab, rank_of_Q_map=rk, closed_basis=basis, expected='S^V_2 = T1 + T3, S^perp_2 = T2 + T4')
    gates.append(dict(id='G-TC0-Sigma', value=basis, target='kernel of Q on {T1..T4, T5\'} = span{T1+T3, T2+T4}; no exact (no degree -2 structure of field weight 9/2)', pass_=ok))
    # ---------------- along Sigma^* (eta in S_+ (x) W^*), position degree <= 1
    eta = [0] * 8 + B.ETA[8:]
    ev = sp.Matrix(eta)
    xis = [list(B.Rbm[mu] * ev) for mu in range(4)]
    pos_ops = []      # per structure i: (X_mu list)
    for i in range(5):
        Xmu = []
        for mu in range(4):
            nm = register_structs(xis[mu], f'tq0e{mu}')
            Xmu.append(E.op_from_functional([nm[i]]))
        pos_ops.append(Xmu)
    intr = register_eta_intrinsic(eta, 'tq0ei')
    intr_ops = [E.op_from_functional([nm]) for nm in intr]
    afrows = [i for i in range(NT) if B.SPEC_OF[i] in ('As', 'lams', 'cs')]
    K0inv = K0.inv()
    xparts = []; cparts = []
    for Xmu in pos_ops:
        RHS = sp.zeros(NT, NT)
        for mu in range(4):
            Xm_ = Xmu[mu].subs({k[m]: mk[m] for m in range(4)}, simultaneous=True)
            RHS += (K0 * Xm_.T).applyfunc(lambda e, mu=mu: sp.expand(sp.diff(e, k[mu])))
        RHS = (I * PI * RHS).applyfunc(sp.expand)
        Rr = sp.zeros(NT, NT)
        for i in afrows:
            for j in range(NT):
                Rr[i, j] = RHS[i, j]
        Y = (Rr * K0inv).applyfunc(sp.expand)
        cpart = (Y * q + q * Y).applyfunc(sp.expand)
        for mu in range(4):
            cpart += B.symbol_derivative(q, k, mu) * Xmu[mu]
        xp = [(Xmu[mu] * q + q * Xmu[mu]).applyfunc(sp.expand) for mu in range(4)]
        xparts.append(xp); cparts.append(cpart.applyfunc(sp.expand))
    for X in intr_ops:
        xparts.append([sp.zeros(NT, NT)] * 4)
        cparts.append((X * q + q * X).applyfunc(sp.expand))
    # stack the x^mu parts (4 blocks) and the constant part
    big = []
    for xp, cp in zip(xparts, cparts):
        big.append(sp.Matrix.vstack(*(xp + [cp])))
    a = sp.symbols(f'a0:{len(big)}')
    tot = sp.zeros(5 * NT, NT)
    for ai, Mx in zip(a, big):
        tot += ai * Mx
    eqs = []
    for e in tot:
        e = sp.expand(e)
        if e != 0:
            syms = sorted(e.free_symbols - set(a), key=lambda s: s.name)
            eqs += sp.Poly(e, *syms).coeffs() if syms else [e]
    Mx = sp.Matrix([[sp.expand(eq).coeff(ai) for ai in a] for eq in eqs])
    nse = Mx.nullspace()
    labe = ['x.T1', 'x.T2', 'x.T3', 'x.T4', "x.T5'", 's0', '(rho(A_V) eta, lam*)']
    basise = [{labe[i]: str(v[i]) for i in range(len(labe)) if v[i] != 0} for v in nse]
    expect_e = [sp.Matrix([1, 0, 1, 0, 0, 0, 0]), sp.Matrix([0, 1, 0, 1, 0, 1, 0])]
    Me = sp.Matrix.hstack(*[sp.Matrix(v) for v in nse]) if nse else sp.zeros(7, 0)
    oke = len(nse) == 2 and sp.Matrix.hstack(Me, *expect_e).rank() == 2
    res['Sigma*'] = dict(structures=labe, closed_basis=basise,
                         expected='S^V_2(eta) = x.(T1+T3), S^perp_2(eta) = x.(T2+T4) + s0: the intrinsic direction s_0 is tied to the perp position part; no closed intrinsic direction alone')
    gates.append(dict(id='G-TC0-Sigma*', value=basise, target='closed span{x.(T1+T3), x.(T2+T4)+s0} (position degree <= 1)', pass_=oke))
    res['exactness_count'] = 'degree -2 quadratic contents: (A*,A*) w=6, (A*,lam*) 11/2, (lam*,lam*) 5, (c*,A) 5, (c*,lam) 11/2, (c*,c*) 8: all > 9/2 (and > 7/2): no exact direction'
    res['K_b'] = 'not run in route B (position degree <= 2; open)'
    res['seconds'] = round(time.time() - t0, 1)
    return res, gates


if __name__ == '__main__':
    import sys
    res, gates = run()
    json.dump(dict(results=res, gates=gates), open(sys.argv[1] if len(sys.argv) > 1 else 'tq0_directions.json', 'w'), indent=1, default=str)
    for g in gates:
        print(g['id'], 'PASS' if g['pass_'] else 'FAIL', g['value'])
    print('seconds', res['seconds'])
