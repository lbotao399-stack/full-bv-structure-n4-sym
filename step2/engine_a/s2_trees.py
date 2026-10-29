"""Step 2, engine A, Part C: the non-I-a contributions to the cubic sector of Q W(P(eps,infty), I - I^CT[eps]) at one
loop, evaluated exactly in the mode model (s2_ga), sector by sector, and fitted on the structure bases of s2_Ia.

Contributions (D(eps) = 1/2({I,I}_eps-{I,I}) + hbar Delta_eps I - Q I^CT - {I,I^CT}_eps, Delta_eps I = 0):
  II-a     = -Q I^CT_3                       = -unit log(eps) Q Y_3,       Y_3 = cubic log coefficient (Step 1 multiples)
  III(K_0) = -{I_3, I^CT_2}_0 = C_{K_0}(I_3, I^CT_2)
           = unit [ (9/4) eps^{-1} C_{K0}(I_3, int(A,A)) + log(eps) C_{K0}(I_3, Y_2) ],
             Y_2 = -int(dA,dA) + 4 int(lambda,dslash lambda) - int(dc,A^*)
  I-b      = {I_3, W_1^(2)(eps)}_{K'} = -C_{K'}(I_3, W_1^(2)(eps)),   K'(p) = (e^{-eps p^2} - 1) Khat,
  III(K')  = -{I_3, I^CT_2}_{K'} = C_{K'}(I_3, I^CT_2):
             I-b + III(K') = -C_{K'}(I_3, W_1^(2)(eps) - I^CT_2) = (e^{-eps p^2}-1) x (renormalised self-energy) -> 0;
             eps-expansion through eps^0:  I-b|_{eps^0} = + C_{p^2 Khat}(I_3, w_{-1}),  III(K')|_{eps^0} = - same,
             w_{-1} = eps^{-1} coefficient of W_1^(2) (recomputed here from the self-energy graphs), no eps^{-1}, no log.
  I-d      = C_P(sb(eps), I_3), sb(eps) = special bubbles (2 legs) = sb_{-1}/eps + sb_0 + O(eps) (Step 1, s2_regress2pt),
  II-b     = C_P(-Q I^CT_2, I_3), -Q I^CT_2 = -sb_{-1}/eps;  P(eps) = (1/p^2 - eps + ...) Pnum:
             I-d + II-b = C_{Pnum/p^2}(sb_0, I_3) + O(eps) = tree of Ob_2 with I_3 (Ob_2 = sb_0).
All plane-wave values in units hbar tr_g(ad_X ad_Y)/(16 pi^2 (X,Y)) times the indicated function of eps."""
import json, sys, time, itertools
from fractions import Fraction
import numpy as np
import sympy as sp
import s2_core as C
from r2_core import Q2
import s2_ga as GAm
import s2_Ia as IA
from grass import GA

MULT = json.load(open('/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a/logs/r2_T1_seed13.json'))
Y3MULT = {'AAA': int(MULT['AV,AV,AV']['computed_multiple']), 'llA': int(MULT['l,l,AV']['computed_multiple']),
          'AcAs': int(MULT['c,AV,AsV']['computed_multiple']), 'cllS': int(MULT['c,l,ls']['computed_multiple']),
          'cccS': int(MULT['c,c,cs']['computed_multiple'])}
assert int(MULT['AV,AP,AP']['computed_multiple']) == Y3MULT['AAA'] and int(MULT['l,l,AP']['computed_multiple']) == Y3MULT['llA']
assert int(MULT['c,AP,AsP']['computed_multiple']) == Y3MULT['AcAs']
MASS = Fraction(9, 4)                      # Theorem c4-quadratic (recomputed below from the graphs)
SB_M1 = Fraction(-9, 2)                    # special bubbles: eps^{-1} coefficient of int(dc,A|_V)  (s2_regress2pt)
SB_0 = Fraction(-27, 8)                    # special bubbles: eps^0 coefficient of int(d^*dc,d^*A|_V)
CUBIC = ['AAA', 'llA', 'AcAs', 'cllS', 'cccS']


def selfenergy_mass(rng):
    """eps^{-1} coefficient of the regulated one-loop A|_V and A|_Wperp two-point functions, from the graphs
    (engine_a/r2_graphs, su(2)); returned as the multiple of int(A,A) (plane wave 2 (A(p),A(-p)))"""
    sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/engine_a')
    import r2_core as RC
    from r2_graphs import GraphSum
    out = {}
    for w in ('V', 'P'):
        vals = []
        for it in range(2):
            p = rng.integers(-2, 3, size=4) * 2
            if not p.any(): p[0] = 2
            a1 = np.zeros((10, 3), complex); a2 = np.zeros((10, 3), complex)
            sl = slice(0, 4) if w == 'V' else slice(4, 10)
            a1[sl] = rng.integers(-2, 3, size=a1[sl].shape); a2[sl] = rng.integers(-2, 3, size=a2[sl].shape)
            G = GraphSum([('A', w, a1, p), ('A', w, a2, -p)])
            res = G.coefficients(2)
            tot = Q2(0)
            for v in res.values():
                tot = tot + v['epsinv']
            dot = int(round(np.einsum('ma,ma->', a1.real, a2.real)))
            vals.append(tot * Fraction(1, RC.C_ADJ * 2 * dot))
        assert (vals[0] - vals[1]).iszero()
        out[w] = vals[0]
    return out


def run_sector(name, rng, nsamp_extra=2, verbose=True):
    cfg = IA.SECTORS[name]
    legs_t = cfg['legs']
    t0 = time.time()
    rows = {-1: [], 0: []}
    data = []
    nst = None
    while True:
        while True:
            ps = [rng.integers(-2, 3, size=4) for _ in range(2)]
            ps = [-(ps[0] + ps[1])] + ps
            allm = ps + [-p for p in ps]
            if all(p.any() for p in ps) and len({tuple(p) for p in allm}) == 6:
                break
        L = [(t, w, C.rpol(rng, t, w), p) for (t, w), p in zip(legs_t, ps)]
        st = {e: cfg['st'](L, e) for e in (-1, 0)}
        if nst is None:
            nst = len(st[0][0])
        M = GAm.Model(allm)
        rs = GAm.rvars(3)
        legs = [(t, pol, i, rs[i] if C.PARITY[t] else None) for i, (t, w, pol, p) in enumerate(L)]
        odd_r = [rs[i] for i, (t, w, pol, p) in enumerate(L) if C.PARITY[t]]
        perms = list(itertools.permutations([0, 1, 2]))
        neg = {0: 3, 1: 4, 2: 5}
        pm_all = [(m, neg[m]) for m in range(3)] + [(neg[m], m) for m in range(3)]
        Sfree = M.S_free(pm_all)
        T = {nm: M.cubic(nm, perms) for nm in CUBIC}
        I3 = T['AAA'] + T['llA'] + T['AcAs'] + T['cllS'] + T['cccS']
        Y3 = GA()
        for nm in CUBIC:
            Y3 = Y3 + T[nm] * Y3MULT[nm]
        val = {}
        def full(F):
            return GAm.coeff(M.subst(F, legs, None), odd_r)
        # (a) tree vertex: mode model vs canonical vertex functions of s2_core
        Vc = C.V([(t, pol[None], p) for (t, w, pol, p) in L])
        tree_canon = Q2(0) if isinstance(Vc, int) else Q2.of(complex(np.asarray(Vc).reshape(-1)[0])) * Fraction(1, 2)
        val['tree_vertex_mode_model_minus_canonical'] = full(I3) - tree_canon
        # (b) G_i in the mode model vs the plane-wave formulas of s2_Ia
        Gh = cfg['G'](L)
        for gn, gv in Gh.items():
            val['G_mode_model_minus_formula_' + gn] = full(M.cubic(gn, perms)) - Q2.of(complex(gv))
        # (c) II-a = -unit log(eps) (Q Y_3)
        QY3 = Q2(0)
        for j in range(3):
            others = [lg for i, lg in enumerate(legs) if i != j]
            Fj = M.subst(Y3, others, None)
            QF = M.Q_apply(Fj, Sfree, pm_all, only_modes={j})
            QY3 = QY3 + GAm.coeff(M.subst(QF, [legs[j]], None), odd_r)
        val['IIa_log'] = QY3 * (-1)
        # (d)-(f): contractions with a two-point functional on leg l
        Y2parts = {'dAdA': -1, 'ldl': 4, 'dcAs': -1}
        acc = {k: Q2(0) for k in ('IIIK0_log', 'IIIK0_epsinv', 'Ib_eps0', 'IIIKp_eps0', 'Id_epsinv_nonlocal', 'IIb_epsinv_nonlocal',
                                  'Id_eps0_local', 'IIb_eps0_local', 'tree_Ob2_I3', 'Y3_minus_bracket_I3_Y')}
        for l in range(3):
            i, j = [x for x in range(3) if x != l]
            mneg = neg[l]
            pl = L[l][3]; p2 = int(pl @ pl)
            F = M.subst(I3, [legs[i], legs[j]], None)          # internal leg at mode l
            def two(Fq):
                return M.subst(Fq, [legs[l]], None)            # external leg l, internal at -p_l
            pairs = [(l, mneg), (mneg, l)]
            K0 = M.kernel('K', l, mneg)
            Kp2 = M.kernel('K', l, mneg, scalar=Q2(p2))
            Pd = M.kernel('P', l, mneg, scalar=Q2(Fraction(1, p2)))
            P1 = M.kernel('P', l, mneg)
            Y2 = GA()
            for nm, c_ in Y2parts.items():
                Y2 = Y2 + M.quad(nm, pairs) * c_
            AA = M.quad('AA', pairs)
            def cv(X, G2):
                return GAm.coeff(M.C(X, F, two(G2)), odd_r)
            acc['IIIK0_log'] = acc['IIIK0_log'] + cv(K0, Y2)
            acc['IIIK0_epsinv'] = acc['IIIK0_epsinv'] + cv(K0, AA) * MASS
            acc['Ib_eps0'] = acc['Ib_eps0'] + cv(Kp2, AA) * MASS_G[legs_key(L[l])]
            acc['IIIKp_eps0'] = acc['IIIKp_eps0'] - cv(Kp2, AA) * MASS
            dcAV = M.quad('dcAV', pairs); ddc = M.quad('ddcdA', pairs)
            acc['Id_epsinv_nonlocal'] = acc['Id_epsinv_nonlocal'] + cv(Pd, dcAV) * SB_M1
            acc['IIb_epsinv_nonlocal'] = acc['IIb_epsinv_nonlocal'] - cv(Pd, dcAV) * SB_M1
            acc['Id_eps0_local'] = acc['Id_eps0_local'] - cv(P1, dcAV) * SB_M1
            acc['IIb_eps0_local'] = acc['IIb_eps0_local'] + cv(P1, dcAV) * SB_M1
            acc['tree_Ob2_I3'] = acc['tree_Ob2_I3'] + cv(Pd, ddc) * SB_0
            # {I_3, Y} = -C_{K0}(I_3, Y), Y = 2 int(A,A^*) - 4 int(lambda^*,lambda) - int(c^*,c)
            Y = M.quad('AAs', pairs) * 2 + M.quad('lsl', pairs) * (-4) + M.quad('csc', pairs) * (-1)
            acc['Y3_minus_bracket_I3_Y'] = acc['Y3_minus_bracket_I3_Y'] + cv(K0, Y)
        acc['Y3_minus_bracket_I3_Y'] = full(Y3) + acc['Y3_minus_bracket_I3_Y']
        # -Q I^CT_2 (mass part) = (9/2) eps^{-1} int(dc, A|_V): two-leg configuration c at p_0 (mode 0), A|_V at -p_0 (mode 3)
        rc = GAm.rvars(1)[0]
        g0 = C.rpol(rng, 'c'); a0 = C.rpol(rng, 'A', 'V')
        legs2 = [('c', g0, 0, rc), ('A', a0, 3, None)]
        AA2 = M.quad('AA', [(0, 3), (3, 0)])
        QAA = Q2(0)
        for jj in range(2):
            Fj = M.subst(AA2, [legs2[1 - jj]], None)
            QAA = QAA + GAm.coeff(M.subst(M.Q_apply(Fj, Sfree, pm_all, only_modes={legs2[jj][2]}), [legs2[jj]], None), [rc])
        dcv = GAm.coeff(M.subst(M.quad('dcAV', [(0, 3), (3, 0)]), legs2, None), [rc])
        val['-Q(9/4 int(A,A)) - (9/2) int(dc,A|_V)'] = QAA * (-MASS) - dcv * Fraction(9, 2)
        val['int(dc,A|_V) two-leg value nonzero'] = Q2(0) if not dcv.iszero() else Q2(1)
        val.update(acc)
        val['Id_eps0_total'] = acc['tree_Ob2_I3'] + acc['Id_eps0_local']
        val['IIa+IIIK0_log'] = val['IIa_log'] + acc['IIIK0_log']
        val['Ib+IIIKp_eps0'] = acc['Ib_eps0'] + acc['IIIKp_eps0']
        val['Id+IIb_epsinv'] = acc['Id_epsinv_nonlocal'] + acc['IIb_epsinv_nonlocal']
        val['Id+IIb_eps0_minus_tree'] = val['Id_eps0_total'] + acc['IIb_eps0_local'] - acc['tree_Ob2_I3']
        data.append(dict(st=st, val=val))
        if verbose:
            print(name, C.ALG, 'sample', len(data), 'time %.1f' % (time.time() - t0),
                  {k: str(v) for k, v in val.items() if k in ('tree_vertex_mode_model_minus_canonical', 'IIa+IIIK0_log',
                                                          'Ib+IIIKp_eps0', 'Id+IIb_epsinv', 'Id+IIb_eps0_minus_tree',
                                                          'Y3_minus_bracket_I3_Y')}, flush=True)
        if len(data) >= nst + nsamp_extra:
            break
    names0 = data[0]['st'][0][1]; namesm = data[0]['st'][-1][1]
    out = {'sector': name, 'algebra': C.ALG, 'samples': len(data)}
    exact_zero = ['-Q(9/4 int(A,A)) - (9/2) int(dc,A|_V)', 'tree_vertex_mode_model_minus_canonical', 'IIa+IIIK0_log', 'Ib+IIIKp_eps0', 'Id+IIb_epsinv',
                  'Id+IIb_eps0_minus_tree', 'Y3_minus_bracket_I3_Y'] + [k for k in data[0]['val'] if k.startswith('G_mode')]
    out['exact_zero_checks'] = {k: all(d['val'][k].iszero() for d in data) for k in exact_zero}
    out['IIIK0_epsinv_identically_zero'] = all(d['val']['IIIK0_epsinv'].iszero() for d in data)
    fits = {}
    for k in ('IIa_log', 'IIIK0_log', 'Ib_eps0', 'IIIKp_eps0', 'Id_eps0_local', 'IIb_eps0_local', 'tree_Ob2_I3', 'Id_eps0_total'):
        rows = [[IA.toS(Q2.of(x)) for x in d['st'][0][0]] for d in data]
        sol, ok = IA.solve_exact(rows, [IA.toS(d['val'][k]) for d in data])
        fits[k] = {'structures': {nm: str(x) for nm, x in zip(names0, sol) if x != 0}, 'zero_residual': ok}
    out['fits_structure_basis'] = fits
    out['seconds'] = round(time.time() - t0, 1)
    return out


MASS_G = {}


def legs_key(leg):
    return (leg[0], leg[1])


if __name__ == '__main__':
    rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 99)
    if C.ALG == 'su2':
        mg = selfenergy_mass(np.random.default_rng(3))
    else:
        mg = {'V': Q2(MASS), 'P': Q2(MASS)}
    print('eps^{-1} coefficient of the regulated A two-point function (multiple of int(A,A)):', mg, flush=True)
    # I-b uses the mass of the self-energy carried by the leg on the two-point side: A|_V or A|_Wperp
    MASS_G.update({('A', 'V'): mg['V'], ('A', 'P'): mg['P'], ('c', None): Q2(MASS), ('l', None): Q2(MASS), ('As', 'V'): Q2(MASS)})
    which = sys.argv[1].split(';') if len(sys.argv) > 1 and sys.argv[1] != 'all' else list(IA.SECTORS)
    res = {'selfenergy_eps^-1_from_graphs': {k: str(v) for k, v in mg.items()}}
    for nm in which:
        r = run_sector(nm, rng)
        res[nm] = r
        print(json.dumps(r, indent=1), flush=True)
        json.dump(res, open('logs/s2_trees_%s.json' % C.ALG, 'w'), indent=1)
