#!/usr/bin/env python3
"""QA, Step 2 Part C: plug-back verification of the finite counterterm F_3 in the theory realisation.

Uses the exact plane-wave realisation step2/theory/cubic_core.py (su(2), Q(i) arithmetic) on FRESH random plane waves
(seeds disjoint from those used by the theory agent) and checks, sector by sector,

  (E)  Q F_3 + {I_3, F_2} - sum_i alpha_i G_i = 0,

with the agreed alpha_i = r_i + l_i log(3/2) (rational part and log(3/2) part checked separately: all maps are linear),
F_2 = (27/16) D1 + beta_1 B1 + beta_2 B2 + beta_3 B3 + beta_4 B4, F_3 = particular + beta-directions + kernel directions,
for beta = kappa = 0 and for random rational beta, kappa.  Sectors: the four target sectors and the degree-1 three-leg
sectors without invariant ((c,A|_V,A|_Wperp), (c,A,lambda), (c,c,lambda^*), (c,c,A^*|_Wperp)), where every term must vanish.

Further checks in the same realisation:
  (WZ) Q Ob_3 + {I_3, Ob_2} = 0 on the degree-2 sector (c,c,A|_V), Ob_2 = -(27/8) int(d^*dc, d^*A|_V);
  (Y)  II-a log part  = -Q Y_3, Y_3 = 6 AAA + 10 lambda lambda A + AcA^* + c lambda lambda^* + ccc^* (action vertices),
       III(K_0) log part = -{I_3, Y_2}, Y_2 = -B1 - B2 + 4 B3 - B4 (Step-1 log coefficients), compared with the engines;
  (M)  III(K_eps-K_0) eps^0 part = (9/4) {I_3, int(A, H A)}, compared with the engines.
Writes qa/plugback_theory.json.
"""
import itertools
import json
import random
import sys
import time
from fractions import Fraction as Fr

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/theory')
import cubic_core as C  # noqa: E402
from cubic_core import GQ, ZERO, ONE, Mode  # noqa: E402

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/plugback_theory.json'

SECTORS = {
    'c_AV_AV': [('c', None), ('A', 'V'), ('A', 'V')],
    'c_AW_AW': [('c', None), ('A', 'W'), ('A', 'W')],
    'c_lam_lam': [('c', None), ('lam', None), ('lam', None)],
    'c_c_AsV': [('c', None), ('c', None), ('As', 'V')],
    # degree-1 three-leg sectors without invariant (everything must vanish)
    'c_AV_AW': [('c', None), ('A', 'V'), ('A', 'W')],
    'c_A_lam': [('c', None), ('A', 'all'), ('lam', None)],
    'c_c_lams': [('c', None), ('c', None), ('lams', None)],
    'c_c_AsW': [('c', None), ('c', None), ('As', 'W')],
    # degree-2 sector for the Wess-Zumino check
    'deg2_c_c_AV': [('c', None), ('c', None), ('A', 'V')],
}
SECTOR_G = {'c_AV_AV': ['G1', 'G2'], 'c_AW_AW': ['G3'], 'c_lam_lam': ['G5'], 'c_c_AsV': ['G4']}


def rint(rng, a=-5, b=5):
    return rng.randint(a, b)


def sample(sector, seed):
    rng = random.Random(seed)
    legs = SECTORS[sector]
    n = len(legs)
    while True:
        ks = [tuple(Fr(rint(rng)) for _ in range(4)) for _ in range(n - 1)]
        ks.append(tuple(-sum((k[i] for k in ks), Fr(0)) for i in range(4)))
        sums, ok = set(), True
        for r in range(1, n):
            for sub in itertools.combinations(range(n), r):
                s = tuple(sum((ks[j][i] for j in sub), Fr(0)) for i in range(4))
                if C.iszerok(s):
                    ok = False
                sums.add(s)
        if ok:
            break
    modes, gens = [], []
    for j, (s, pol) in enumerate(legs):
        g = C.new_gen(C.SPAR[s], 'ext')
        gens.append(g)
        if s in ('c', 'cs'):
            u = [ONE]
        elif s in ('A', 'As'):
            idx = C.V4 if pol == 'V' else (C.WP if pol == 'W' else range(10))
            u = [GQ(rint(rng, -3, 3)) if M in idx else ZERO for M in range(10)]
        else:
            u = [GQ(rint(rng, -3, 3)) for _ in range(C.SDIM[s])]
        X = [Fr(rint(rng, -3, 3)) for _ in range(C.DIMG)]
        modes.append(Mode(C.gen_elem(g), s, u, X, ks[j]))
    return modes, gens, sorted(sums)


def val(elem, gens):
    return C.ext_coeff(elem, gens)


def lin(coeffs, values):
    tot = ZERO
    for k, c in coeffs.items():
        if c != 0:
            tot = tot + values[k] * Fr(c)
    return tot


def fmt(x):
    if x.im == 0:
        return str(x.re)
    return f'{x.re}+{x.im}*i'


# agreed alpha (rational part, log(3/2) part)
ALPHA_R = {'G1': Fr(11, 4), 'G2': Fr(27, 8), 'G3': Fr(11, 4), 'G4': Fr(0), 'G5': Fr(0)}
ALPHA_L = {'G1': Fr(19, 4), 'G2': Fr(0), 'G3': Fr(19, 4), 'G4': Fr(0), 'G5': Fr(-9)}


def F3_coeffs(beta, kappa):
    """F_3 = (alpha4/4 - alpha1/2) H1 - alpha4 H3 - 2 alpha3 H4 + alpha5 H5
             + beta1 H1 + 2 beta2 H2 - 2 beta3 H5 + beta4 (H3 + H4 + H5 + H6/2)
             + kappa1 (H1 - 4 H3 - 2 H6) + kappa2 (H2 - 2 H4); returns (rational part, log(3/2) part)."""
    def part(a):
        return {'H1': a['G4'] / 4 - a['G1'] / 2, 'H2': Fr(0), 'H3': -a['G4'], 'H4': -2 * a['G3'], 'H5': a['G5'],
                'H6': Fr(0)}
    hr, hl = part(ALPHA_R), part(ALPHA_L)
    b1, b2, b3, b4 = beta
    k1, k2 = kappa
    add = {'H1': b1 + k1, 'H2': 2 * b2 + k2, 'H3': b4 - 4 * k1, 'H4': b4 - 2 * k2, 'H5': -2 * b3 + b4,
           'H6': b4 / 2 - 2 * k1}
    return {k: hr[k] + add[k] for k in hr}, hl


def main():
    t0 = time.time()
    out = {'realisation': 'step2/theory/cubic_core.py (su(2), exact over Q(i)); fresh seeds 900000+',
           'F3_particular_beta0_kappa0': {
               'rational_part': {k: str(v) for k, v in F3_coeffs((0, 0, 0, 0), (0, 0))[0].items()},
               'log(3/2)_part': {k: str(v) for k, v in F3_coeffs((0, 0, 0, 0), (0, 0))[1].items()}},
           'equation_E': {}, 'WZ': {}, 'Y_log_checks': {}, 'III_Kprime_check': {}}
    rngp = random.Random(424242)
    param_sets = [((0, 0, 0, 0), (0, 0))] + [
        (tuple(Fr(rngp.randint(-9, 9), rngp.randint(1, 5)) for _ in range(4)),
         tuple(Fr(rngp.randint(-9, 9), rngp.randint(1, 5)) for _ in range(2))) for _ in range(3)]
    nsamp = 5
    for sec in list(SECTOR_G) + ['c_AV_AW', 'c_A_lam', 'c_c_lams', 'c_c_AsW']:
        recs = []
        for s in range(nsamp):
            modes, gens, kset = sample(sec, 900000 + 1000 * len(sec) + 17 * s + 7 * sum(map(ord, sec)))
            cache = C.make_internal(kset)
            vG = {g: val(C.evaluate(C.G[g], modes), gens) for g in C.G}
            vQH = {h: val(C.applyQ(C.H[h], modes), gens) for h in C.H}
            vBr = {x: val(C.bracket(C.I3, 0, C.Q2[x], modes, kset, cache), gens)
                   for x in ('D1', 'B1', 'B2', 'B3', 'B4', 'massH')}
            vQI3 = {t: val(C.applyQ(C.I3_terms[t], modes), gens) for t in C.I3_terms}
            rec = {'seed_index': s}
            # equation (E) for every parameter set
            ok_all = True
            for beta, kappa in param_sets:
                hr, hl = F3_coeffs(beta, kappa)
                f2 = {'D1': Fr(27, 16), 'B1': beta[0], 'B2': beta[1], 'B3': beta[2], 'B4': beta[3]}
                lhs_r = lin(hr, vQH) + lin(f2, vBr)
                rhs_r = lin(ALPHA_R, vG)
                lhs_l = lin(hl, vQH)
                rhs_l = lin(ALPHA_L, vG)
                ok = (lhs_r - rhs_r).iszero() and (lhs_l - rhs_l).iszero()
                ok_all = ok_all and ok
                if beta == (0, 0, 0, 0) and kappa == (0, 0):
                    rec['beta0_kappa0'] = {'LHS_rational': fmt(lhs_r), 'RHS_rational': fmt(rhs_r),
                                           'LHS_log32': fmt(lhs_l), 'RHS_log32': fmt(rhs_l), 'holds': ok}
            rec['all_parameter_sets_hold'] = ok_all
            # representative with F_3 = 0: F_2 = (27/16) D1 + (alpha_1/2)(B1 + B2) + (alpha_5/2) B3  (alpha_1 = alpha_3)
            z_r = vBr['D1'] * Fr(27, 16) + (vBr['B1'] + vBr['B2']) * Fr(11, 8)
            z_l = (vBr['B1'] + vBr['B2']) * Fr(19, 8) + vBr['B3'] * Fr(-9, 2)
            rec['F3_zero_representative_holds'] = (z_r - lin(ALPHA_R, vG)).iszero() and (z_l - lin(ALPHA_L, vG)).iszero()
            rec['all_G_values'] = {g: fmt(v) for g, v in vG.items()}
            # (Y) log checks
            qy3 = (vQI3['gauge'] * 6 + vQI3['llA'] * 10 + vQI3['AcAs'] + vQI3['clls'] + vQI3['cccs'])
            IIa_log = -qy3
            IIIK0_log = -(lin({'B1': -1, 'B2': -1, 'B3': 4, 'B4': -1}, vBr))
            rec['IIa_log_value'] = fmt(IIa_log)
            rec['IIIK0_log_value'] = fmt(IIIK0_log)
            rec['IIIKprime_eps0_value'] = fmt(vBr['massH'] * Fr(9, 4))
            rec['QI3_total'] = fmt(sum(vQI3.values(), ZERO))
            rec['_vG'] = vG
            rec['_IIa'] = IIa_log
            rec['_IIIK0'] = IIIK0_log
            rec['_IIIKp'] = vBr['massH'] * Fr(9, 4)
            recs.append(rec)
        # decompose the Y / III(K') values on the sector's G basis (exact solve over the samples)
        dec = {}
        if sec in SECTOR_G:
            basis = SECTOR_G[sec]
            for key in ('_IIa', '_IIIK0', '_IIIKp'):
                rows = [[r['_vG'][g] for g in basis] for r in recs]
                rhs = [r[key] for r in recs]
                x, okc, rank = solve(rows, rhs)
                dec[key[1:]] = {**{g: fmt(v) for g, v in zip(basis, x)}, '_residual_zero': okc, '_rank': rank}
        for r in recs:
            for k in ('_vG', '_IIa', '_IIIK0', '_IIIKp'):
                r.pop(k)
        out['equation_E'][sec] = {'samples': recs, 'all_hold': all(r['all_parameter_sets_hold'] for r in recs),
                                  'F3_zero_representative_all_hold': all(r['F3_zero_representative_holds'] for r in recs)}
        if dec:
            out['Y_log_checks'][sec] = {k: dec[k] for k in ('IIa', 'IIIK0')}
            out['III_Kprime_check'][sec] = dec['IIIKp']
        print(sec, out['equation_E'][sec]['all_hold'], json.dumps(dec), flush=True)
    # (WZ) on (c, c, A|_V)
    wz = []
    for s in range(3):
        modes, gens, kset = sample('deg2_c_c_AV', 910000 + 31 * s)
        cache = C.make_internal(kset)
        qg = {g: val(C.applyQ(C.G[g], modes), gens) for g in C.G}
        b_ob2 = val(C.bracket(C.I3, 0, C.Q2['Ob2'], modes, kset, cache), gens) * Fr(-27, 8)
        tot_r = lin(ALPHA_R, qg) + b_ob2
        tot_l = lin(ALPHA_L, qg)
        wz.append({'QG': {g: fmt(v) for g, v in qg.items()}, '{I3,Ob2}': fmt(b_ob2),
                   'Q Ob3 + {I3,Ob2} (rational part)': fmt(tot_r), 'Q Ob3 (log(3/2) part)': fmt(tot_l),
                   'zero': tot_r.iszero() and tot_l.iszero()})
    out['WZ'] = {'samples': wz, 'all_zero': all(w['zero'] for w in wz)}
    print('WZ', out['WZ']['all_zero'], flush=True)
    out['parameter_sets'] = [{'beta': [str(b) for b in be], 'kappa': [str(k) for k in ka]} for be, ka in param_sets]
    out['seconds'] = round(time.time() - t0, 1)
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print('done', out['seconds'])


def solve(rows, rhs):
    n = len(rows[0])
    A = [list(r) + [b] for r, b in zip(rows, rhs)]
    piv, r = [], 0
    for col in range(n):
        p = next((i for i in range(r, len(A)) if not A[i][col].iszero()), None)
        if p is None:
            continue
        A[r], A[p] = A[p], A[r]
        pv = A[r][col]
        A[r] = [x / pv for x in A[r]]
        for i in range(len(A)):
            if i != r and not A[i][col].iszero():
                f = A[i][col]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(col)
        r += 1
    ok = all(A[i][n].iszero() for i in range(r, len(A)))
    x = [ZERO] * n
    for i, col in enumerate(piv):
        x[col] = A[i][n]
    return x, ok, len(piv)


if __name__ == '__main__':
    main()
