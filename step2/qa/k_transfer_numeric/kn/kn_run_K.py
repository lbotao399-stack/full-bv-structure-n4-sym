"""Part (a) / (b): arity one along the special conformal field K_b, v(y) = |y|^2 b - 2 (b.y) y, quadratic sector.
(a) mode 'defect': 1PI one-loop graphs with one special edge K_eps - K_0, one vertex S^(1)(K_b), vertices of I,
    edges P(eps,T); external legs (A|_W,A|_W), (A|_V,A|_V), (lambda,lambda), (c,A^*).
(b) mode 'flow': 1PI one-loop graphs of W(P(eps,T), I + S^(1)(K_b)) (no special edge); external legs
    (A|_W,A^*|_W), (A|_V,A^*|_V), (lambda,lambda^*), (c,c^*).
K_b is reached from the plane-wave insertions S^(1)(e_m e^{iP.y}) (leg 1 at momentum p1 fixed, leg 2 at -p1-P):
    graph_K = - sum_{m,a,c} C^m_{ac} d_a d_c [delta(Sigma+P) f_m(p1,P)]_{P=0},  C^m_{ac} = delta_ac b_m - b_a delta_cm - b_c delta_am,
    = -( X2_{ac} d_a d_c delta(Sigma) + X1_a d_a delta(Sigma) + X0 delta(Sigma) ),
    X2_{ac} = sum_m C^m_{ac} f_m(0), X1_a = 2 sum_{m,c} C^m_{ac} d_c f_m(0), X0 = sum_{m,a,c} C^m_{ac} d_a d_c f_m(0).
Units hbar/(16 pi^2); divide by -2 for the unit hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) (su(2), orthonormal basis).
Usage: python3 kn_run_K.py <mode> <sector> <tag>"""
import sys, os, json, time, pickle
import numpy as np
import kn_engine as E

OUT = os.path.dirname(os.path.abspath(__file__))
I4 = np.eye(4)
P1 = np.array([int(x) for x in os.environ.get('KN_P1', '1,2,0,1').split(',')])
BVEC = np.array([1, 0, -1, 2])


def Cm(b):
    C_ = np.zeros((4, 4, 4))        # C[m, a, c]
    for m in range(4):
        for a in range(4):
            for c in range(4):
                C_[m, a, c] = (a == c) * b[m] - b[a] * (c == m) - b[c] * (a == m)
    return C_


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


SECTORS = {
    'defect': {
        'AW_AW': (('A', 'AW', 1), ('A', 'AW', 2)),
        'AV_AV': (('A', 'AV', 3), ('A', 'AV', 4)),
        'l_l': (('l', 'l', 5), ('l', 'l', 6)),
        'c_As': (('c', 'c', 7), ('As', 'AV', 8)),
    },
    'flow': {
        'AW_AsW': (('A', 'AW', 1), ('As', 'AW', 2)),
        'AV_AsV': (('A', 'AV', 3), ('As', 'AV', 4)),
        'l_ls': (('l', 'l', 5), ('ls', 'l', 6)),
        'c_cs': (('c', 'c', 7), ('cs', 'c', 8)),
    },
}


def build(mode, sector):
    (t1, k1, s1), (t2, k2, s2) = SECTORS[mode][sector]
    pol1, pol2 = pol(k1, s1), pol(k2, s2)
    jet = E.jet_order2(4)
    leg1 = E.bg_insertion(t1, pol1, P1, 'L1')
    leg2 = E.bg_insertion(t2, pol2, -P1, 'L2', mom_lin=-I4)
    canon = [x for x, t in (('L1', t1), ('L2', t2)) if E.PAR[t]] + ['g']
    special = (mode == 'defect')
    cyc_by_m = []
    for m in range(4):
        S = E.s1v_insertion(I4[m], 'g', pin_lin=I4)
        sets = [[leg1, leg2, S]]
        if t1 == 'A' and t2 == 'A':
            q = E.bg_insertion('A', pol1, P1, 'L1', mom_lin=None, quartic_with=('A', pol2, 'L2'))
            q.pin_const = np.zeros(4); q.pin_lin = -I4         # total injected momentum p1 + p2 = -P
            sets.append([q, S])
        cl = []
        for st in sets:
            for perm, Ls, cy in E.enumerate_cycles(st, jet, canon, special=special):
                if cy.is_zero():
                    continue
                names = [I.name for I in cy.ins]
                iS = names.index('S1v')
                n = cy.n
                Lin = cy.Ls[iS]                 # line entering S (j-leg)
                Lout = cy.Ls[(iS + 1) % n]      # line leaving S (l-leg side): line iS+1 enters the next insertion
                floating = (Lin in E.PLINES) and (Lout in E.PLINES)
                cls = 'D_composite' if (floating and special) else ('flow_composite' if floating else 'B_anchored')
                label = '%s|%s' % ('-'.join(names), ','.join(cy.Ls))
                cl.append((label, cls, cy))
        cyc_by_m.append(cl)
    return jet, cyc_by_m, canon


def kb_combine(Fm, jet, b):
    """Fm: (4 [u = e_m], nb) jets -> X0, X1 (4), X2 (4,4)"""
    C_ = Cm(b)
    f0 = Fm[:, 0]
    d1 = np.array([[Fm[m, jet.lin[a]] for a in range(4)] for m in range(4)])     # (m, a)
    d2 = np.zeros((4, 4, 4), complex)
    for m in range(4):
        for a in range(4):
            for c in range(4):
                e = tuple((1 if k == a else 0) + (1 if k == c else 0) for k in range(4))
                v = Fm[m, jet.index[e]]
                d2[m, a, c] = 2 * v if a == c else v
    X2 = np.einsum('mac,m->ac', C_, f0)
    X1 = 2 * np.einsum('mac,mc->a', C_, d1)
    X0 = np.einsum('mac,mac->', C_, d2)
    return X0, X1, X2


def main():
    mode, sector, tag = sys.argv[1], sys.argv[2], sys.argv[3]
    epss = [float(x) for x in os.environ.get('KN_EPS', '1e-1,3e-2,1e-2,3e-3,1e-3,3e-4,1e-4').split(',')]
    Ts = [float(x) for x in os.environ.get('KN_T', '1e1,1e2,1e3,1e4,1e5,1e6').split(',')]
    pd = int(os.environ.get('KN_PD', '10'))
    t0 = time.time()
    jet, cyc_by_m, canon = build(mode, sector)
    print('built', [len(c) for c in cyc_by_m], 'cycles per m; %.1fs' % (time.time() - t0), flush=True)
    # numerators
    asm_by_m = []
    for m in range(4):
        al = []
        for (label, cls, cy) in cyc_by_m[m]:
            al.append((label, cls, cy, E.Assembled(cy, cy.numerator())))
        asm_by_m.append(al)
    print('numerators done %.1fs' % (time.time() - t0), flush=True)
    labels = sorted({lab for al in asm_by_m for (lab, _, _, _) in al})
    classes = {lab: cls for al in asm_by_m for (lab, cls, _, _) in al}
    res = {'mode': mode, 'sector': sector, 'p1': P1.tolist(), 'b': BVEC.tolist(), 'canon_odd': canon,
           'per_decade': pd, 'units': 'hbar/(16 pi^2); unit = -2 x this', 'labels': labels, 'classes': classes,
           'values': {}, 'TdT': {}}
    for eps in epss:
        for T in Ts:
            if T <= eps * 10:
                continue
            key = '%g|%g' % (eps, T)
            per = {}
            perT = {}
            for m in range(4):
                for (label, cls, cy, asm) in asm_by_m[m]:
                    pf = E.prefactor(cy)
                    if mode == 'flow':
                        v = E.integrate_flow_simplex(asm, eps, T, per_decade=pd) * pf
                    else:
                        v = E.integrate_cycle(asm, eps, T, per_decade=pd) * pf
                    vt = E.integrate_cycle(asm, eps, T, per_decade=pd, mode='TdT') * pf
                    per.setdefault(label, np.zeros((4, jet.nb), complex))[m] += v
                    perT.setdefault(label, np.zeros((4, jet.nb), complex))[m] += vt
            out = {}
            outT = {}
            for label in per:
                X0, X1, X2 = kb_combine(per[label], jet, BVEC)
                out[label] = {'X0': [X0.real, X0.imag], 'X1': [[x.real, x.imag] for x in X1],
                              'X2': [[[x.real, x.imag] for x in row] for row in X2]}
                Y0, Y1, Y2 = kb_combine(perT[label], jet, BVEC)
                outT[label] = {'X0': [Y0.real, Y0.imag], 'X1': [[x.real, x.imag] for x in Y1]}
            res['values'][key] = out
            res['TdT'][key] = outT
            print(key, 'done %.1fs' % (time.time() - t0), flush=True)
        json.dump(res, open(os.path.join(OUT, 'raw_%s_%s_%s.json' % (mode, sector, tag)), 'w'))
    res['time_s'] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, 'raw_%s_%s_%s.json' % (mode, sector, tag)), 'w'))
    print('total %.1fs' % (time.time() - t0))


if __name__ == '__main__':
    main()
