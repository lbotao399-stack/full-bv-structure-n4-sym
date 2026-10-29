"""Part (c): arity two along two special supercharges eta1, eta2 in Sigma^* with Gamma(eta1,eta2)|_V != 0.
Killing spinors eps_a(x) = rho(x) eta_a = sum_m x^m s[m] eta_a (ch02 def:c2-odd; engine A: rho(v)|_{Sigma^*} = s[M] v^M).
Plane-wave decomposition eps_a = zeta e^{i P_a.x} with zeta = s[m_a] eta_a, x^m -> -i d/dP_{a,m} (hyper-dual jets,
P_1 = h1 e_{m1}, P_2 = h2 e_{m2}).  Generating function with even ghosts c1, c2 (odd generators):
   J = I + c1 S^(1)(eps_1) + c2 S^(1)(eps_2) - c1 c2 S^(2)(eps_1, eps_2)
(the relative sign -1 reproduces the arity-two component of the master equation, ch02 lem:c2-components, for odd
arguments); the c1 c2 coefficient is minus the arity-two component, a common factor of all graphs.
S^(1)(eps) (ch02 eq:c2-S1Q, eq:c2-S1S, general spinor field eps):
   - (Gamma(eps,lambda), A^*) + 1/2 (rho(F_A) eps, lambda^*) - 1/4 (rho(A|_W) dslash eps, lambda^*),
   rho(q(k^a)) = rho(k) rho(a) - (k,a), dslash eps = i rho(P) zeta e^{iPx} in Sigma^*; for eps = rho(x) eta the last
   term is -(rho(A|_W) eta, lambda^*).
S^(2)(eps1,eps2) = 1/4 (Gamma(eps1,eps2), Gamma(lambda^*,lambda^*)) - 1/2 (eps1,lambda^*)(eps2,lambda^*)
                   - (iota_{Gamma(eps1,eps2)} A, c^*)   (ch02 eq:c2-S2, pointwise).
All Hessians are scaled by 4 (integer entries)."""
import sys, os, json, time, itertools
import numpy as np
import kn_engine as E

C = E.C
D = E.D
NT = E.NT
s = C.SARR
sb = C.SBARR
F = C.F                                   # f_abc = eps_abc
OUT = os.path.dirname(os.path.abspath(__file__))
I4 = np.eye(4)
ETA1 = np.zeros(16); ETA1[0] = 1; ETA1[5] = 1
ETA2 = np.zeros(16); ETA2[0] = 1; ETA2[5] = -1


def sV(k):
    return np.einsum('m,mab->ab', np.asarray(k[:4], complex), s[:4])


def sbV(k):
    return np.einsum('m,mab->ab', np.asarray(k[:4], complex), sb[:4])


def E_quad(zeta):
    zeta = np.asarray(zeta, complex)
    zsb = np.array([zeta @ sb[N] for N in range(10)])          # (N, beta): (zeta^T sb[N])_beta

    def fn(kj, kl, P):
        H = np.zeros((NT, NT), complex)
        for a in range(D):
            # Term 1: -(Gamma(eps,lambda), A^*):  H[lambda, A^*] = +zsb, H[A^*, lambda] = -zsb
            for N in range(10):
                iS = C.OFF['As'] + N * D + a
                for b in range(16):
                    iL = C.OFF['l'] + b * D + a
                    v = 4 * zsb[N, b]
                    if v != 0:
                        H[iL, iS] += v
                        H[iS, iL] += -v
            # Term 2: 1/2 (rho(dA) eps, lambda^*), A-leg momentum k: (i/2)[s(k) sb[M] zeta - k_M zeta]
            for M in range(10):
                iA = C.OFF['A'] + M * D + a
                vj = 2j * (sV(kj) @ sb[M] @ zeta - (kj[M] if M < 4 else 0) * zeta)
                vl = 2j * (sV(kl) @ sb[M] @ zeta - (kl[M] if M < 4 else 0) * zeta)
                v3 = (-1j * (s[M] @ sbV(P) @ zeta)) if M >= 4 else np.zeros(16)
                for b in range(16):
                    isl = C.OFF['ls'] + b * D + a
                    H[iA, isl] += vj[b] + v3[b]
                    H[isl, iA] += vl[b] + v3[b]
        return H
    return fn


def E_cubic_bgA(zeta, alpha):
    """cubic part 1/4 sum f^{abc} A^a_M A^b_N (s[M] sb[N] zeta).lambda^{*c} with one A = alpha (bg):
    H[A_(L,d), lambda*_(beta,c)] = 1/4 sum_{M,a} alpha^a_M f^{adc} ((s[M]sb[L] - s[L]sb[M]) zeta)_beta"""
    zeta = np.asarray(zeta, complex); alpha = np.asarray(alpha, complex)

    def fn(kj, kl, P):
        H = np.zeros((NT, NT), complex)
        for L in range(10):
            for d in range(D):
                iA = C.OFF['A'] + L * D + d
                for c in range(D):
                    vec = np.zeros(16, complex)
                    for M in range(10):
                        for a in range(D):
                            if alpha[M, a] == 0 or F[a, d, c] == 0:
                                continue
                            vec += alpha[M, a] * F[a, d, c] * ((s[M] @ sb[L] - s[L] @ sb[M]) @ zeta)
                    for b in range(16):
                        isl = C.OFF['ls'] + b * D + c
                        H[iA, isl] += vec[b]          # 4 * 1/4
                        H[isl, iA] += vec[b]
        return H
    return fn


def E_cubic_bgls(zeta, w):
    """cubic part with lambda^* = w (bg): H[A_(M,a), A_(N,b)] = 1/4 sum_c f^{abc} ((s[M]sb[N]-s[N]sb[M]) zeta).w^c"""
    zeta = np.asarray(zeta, complex); w = np.asarray(w, complex)

    def fn(kj, kl, P):
        H = np.zeros((NT, NT), complex)
        for M in range(10):
            for N in range(10):
                X = (s[M] @ sb[N] - s[N] @ sb[M]) @ zeta
                for a in range(D):
                    for b in range(D):
                        v = 0
                        for c in range(D):
                            if F[a, b, c] != 0:
                                v += F[a, b, c] * (X @ w[:, c])
                        if v != 0:
                            H[C.OFF['A'] + M * D + a, C.OFF['A'] + N * D + b] += v
        return H
    return fn


def S2_hess(z1, z2):
    """Hessian of -S^(2)(eps1, eps2) (scaled by 4)"""
    z1 = np.asarray(z1, complex); z2 = np.asarray(z2, complex)
    G = np.array([z1 @ sb[M] @ z2 for M in range(10)])

    def fn(kj, kl, P):
        H = np.zeros((NT, NT), complex)
        Mss = 0.5 * np.einsum('M,Mab->ab', G, s) - 0.5 * (np.outer(z1, z2) + np.outer(z2, z1))
        for c in range(D):
            for b1 in range(16):
                for b2 in range(16):
                    H[C.OFF['ls'] + b1 * D + c, C.OFF['ls'] + b2 * D + c] += -4 * Mss[b1, b2]
            for M in range(10):
                iA = C.OFF['A'] + M * D + c; ics = C.OFF['cs'] + c
                H[iA, ics] += -4 * (-G[M]); H[ics, iA] += -4 * (-G[M])
        return H
    return fn


SECTORS = {
    'AW_AsW': (('A', 'AW', 1), ('As', 'AW', 2)),
    'AV_AsV': (('A', 'AV', 3), ('As', 'AV', 4)),
    'l_ls': (('l', 'l', 5), ('ls', 'l', 6)),
    'c_cs': (('c', 'c', 7), ('cs', 'c', 8)),
}
P1 = np.array([1, 2, 0, 1])


def pol(kind, seed):
    rng = np.random.default_rng(seed)
    if kind == 'AW':
        x = np.zeros((10, 3), complex); x[4, 0] = 1; return x
    if kind == 'AV':
        x = np.zeros((10, 3), complex); x[:4, 0] = rng.integers(-2, 3, size=4); return x
    if kind == 'l':
        x = np.zeros((16, 3), complex); x[:, 0] = rng.integers(-1, 2, size=16); return x
    if kind == 'c':
        x = np.zeros(3, complex); x[0] = 1; return x
    raise ValueError(kind)


def lin(d1=None, d2=None):
    m = np.zeros((2, 4))
    if d1 is not None: m[0] = d1
    if d2 is not None: m[1] = d2
    return m


def build(sector, m1, m2, jet):
    (t1, k1, s1), (t2, k2, s2) = SECTORS[sector]
    pol1, pol2 = pol(k1, s1), pol(k2, s2)
    z1 = s[m1] @ ETA1
    z2 = s[m2] @ ETA2
    e1, e2 = I4[m1], I4[m2]
    leg1 = E.bg_insertion(t1, pol1, P1, 'L1', mom_lin=lin())
    leg2 = E.bg_insertion(t2, pol2, -P1, 'L2', mom_lin=lin(-e1, -e2))
    canon = [x for x, t in (('L1', t1), ('L2', t2)) if E.PAR[t]]
    ins = {}
    ins['E1'] = E.Insertion(E_quad(z1), 4.0, ['c1'], [], 'E1', None, lin(e1, None))
    ins['E2'] = E.Insertion(E_quad(z2), 4.0, ['c2'], [], 'E2', None, lin(None, e2))
    ins['S2'] = E.Insertion(S2_hess(z1, z2), 4.0, ['c1', 'c2'], [], 'S2', None, lin(e1, e2))
    sets = [['bg1', 'bg2', 'E1', 'E2'], ['bg1', 'bg2', 'S2']]
    for a, za, ea, tag in (('1', z1, e1, 'c1'), ('2', z2, e2, 'c2')):
        for j, (tj, pj, legmom_c, legmom_l) in enumerate(((t1, pol1, P1, lin()), (t2, pol2, -P1, lin(-e1, -e2)))):
            name = 'E%sb%d' % (a, j + 1)
            pin_l = legmom_l + (lin(ea, None) if a == '1' else lin(None, ea))
            if tj == 'A':
                ins[name] = E.Insertion(E_cubic_bgA(za, pj), 4.0, [tag, 'L%d' % (j + 1)], [], name, legmom_c, pin_l)
            elif tj == 'ls':
                ins[name] = E.Insertion(E_cubic_bgls(za, pj), 4.0, [tag, 'L%d' % (j + 1)], [], name, legmom_c, pin_l)
            else:
                continue
            other_bg = 'bg2' if j == 0 else 'bg1'
            other_E = 'E2' if a == '1' else 'E1'
            sets.append([name, other_bg, other_E])
    # pairs of cubic E insertions
    for x, y in (('E1b1', 'E2b2'), ('E1b2', 'E2b1')):
        if x in ins and y in ins:
            sets.append([x, y])
    ins['bg1'] = leg1; ins['bg2'] = leg2
    cycles = []
    for st in sets:
        for perm, Ls, cy in E.enumerate_cycles([ins[n] for n in st], jet, canon, special=True):
            if cy.is_zero():
                continue
            names = [I.name for I in cy.ins]
            n = cy.n
            fl = {}
            for nm in ('E1', 'E2'):
                if nm in names:
                    i = names.index(nm)
                    fl[nm] = (cy.Ls[i] in E.PLINES) and (cy.Ls[(i + 1) % n] in E.PLINES)
            chain = False
            if fl.get('E1') and fl.get('E2'):
                i1, i2 = names.index('E1'), names.index('E2')
                if (i2 - i1) % n == 1 and cy.Ls[i2] == 'l':
                    chain = True
                if (i1 - i2) % n == 1 and cy.Ls[i1] == 'l':
                    chain = True
            cls = 'chain_eta_eta' if chain else ('floating_eta' if any(fl.values()) else 'anchored')
            cycles.append(('%s|%s' % ('-'.join(names), ','.join(cy.Ls)), cls, cy))
    return cycles, canon


def main():
    sector = sys.argv[1]
    tag = sys.argv[2]
    mlist = [tuple(int(c) for c in x) for x in os.environ.get('KN_M', '00,01,02,03,10,11,12,13,20,21,22,23,30,31,32,33').split(',')]
    epss = [float(x) for x in os.environ.get('KN_EPS', '1e-1,3e-2,1e-2,3e-3,1e-3').split(',')]
    Ts = [float(x) for x in os.environ.get('KN_T', '1e1,1e2,1e3,1e4,1e5').split(',')]
    vals_pts = [tuple(float(y) for y in x.split(':')) for x in os.environ.get('KN_VALS', '1e-2:1e2,3e-3:1e2,1e-3:1e2,1e-2:1e3,3e-3:1e3,1e-3:1e3').split(',')]
    pd = int(os.environ.get('KN_PD', '8'))
    jet = E.jet_hyperdual()
    t0 = time.time()
    res = {'sector': sector, 'p1': P1.tolist(), 'eta1': ETA1.tolist(), 'eta2': ETA2.tolist(),
           'Gamma_eta1_eta2_V': [float((ETA1 @ s[M] @ ETA2).real) for M in range(4)],
           'units': 'hbar/(16 pi^2) (unit = -2 x this); jets [1, h1, h2, h1 h2]; c1c2 coefficient = -(arity-two component)',
           'per_m': {}}
    for (m1, m2) in mlist:
        cycles, canon = build(sector, m1, m2, jet)
        entry = {'n_cycles': len(cycles), 'classes': {}, 'TdT': {}, 'values': {}}
        for (lab, cls, cy) in cycles:
            entry['classes'][lab] = cls
            asm = E.Assembled(cy, cy.numerator())
            pf = E.prefactor(cy)
            for eps in epss:
                for T in Ts:
                    if T <= 10 * eps:
                        continue
                    v = E.integrate_cycle(asm, eps, T, per_decade=pd, mode='TdT') * pf
                    entry['TdT'].setdefault('%g|%g' % (eps, T), {})[lab] = [[z.real, z.imag] for z in v]
            for (eps, T) in vals_pts:
                if cy.n <= 3 or os.environ.get('KN_FULL4', '1') == '1':
                    v = E.integrate_cycle(asm, eps, T, per_decade=pd) * pf
                    entry['values'].setdefault('%g|%g' % (eps, T), {})[lab] = [[z.real, z.imag] for z in v]
        res['per_m']['%d%d' % (m1, m2)] = entry
        print('m=%d%d: %d cycles, %.1fs' % (m1, m2, len(cycles), time.time() - t0), flush=True)
        json.dump(res, open(os.path.join(OUT, 'raw_eta_%s_%s.json' % (sector, tag)), 'w'))
    res['time_s'] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, 'raw_eta_%s_%s.json' % (sector, tag)), 'w'))


if __name__ == '__main__':
    main()
