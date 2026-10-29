#!/usr/bin/env python3
"""QA, Step 2 Part C: second plug-back of F_3, in engine B's realisation (step2/engine_b/fa.py: Q = {S_2,-}, local
bracket, vertex functions with colour ([c0,c1],c2) stripped), on configurations not used by engine B (seeds 40-47).

Checks  Q F_3 + {I_3, F_2} - sum_i alpha_i G_i = 0  in the four sectors, with
  F_2 = (27/16) D + beta_1 B1 + beta_2 B2 + beta_3 B3 + beta_4 B4,
  F_3 = (alpha4/4 - alpha1/2) H1 - alpha4 H3 - 2 alpha3 H4 + alpha5 H5
        + beta1 H1 + 2 beta2 H2 - 2 beta3 H5 + beta4 (H3 + H4 + H5 + H6/2) + kappa1 (H1 - 4 H3 - 2 H6) + kappa2 (H2 - 2 H4),
for beta = kappa = 0 and random rational beta, kappa; alpha = r + l log(3/2), both parts separately.
H_j and the quadratic functionals are realised exactly as in engine_b/solve_F3.py (normalisations stated there; they
agree with the theory realisation: same Q H_j and {I_3, X} matrices).  Writes qa/plugback_engine_b.json."""
import json
import random
import sys
import time

sys.path.insert(0, '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/engine_b')
import sympy as sp  # noqa: E402
import fa  # noqa: E402
import sectors3 as S  # noqa: E402
import vertices as V  # noqa: E402

OUT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27/step2/qa/plugback_engine_b.json'
I = sp.I
dot4 = V.dot4
fa.add_term('H1', 1, ('A', 'A', 'A'), lambda s: 2 * sum(I * s[0]['p'][m] * s[0]['u'][n] * s[1]['u'][m] * s[2]['u'][n] for m in range(4) for n in range(4)), ('f', (0, 1, 2)))
fa.add_term('H2', 1, ('A', 'A', 'A'), lambda s: sum(I * s[0]['p'][m] * s[0]['u'][w] * s[1]['u'][m] * s[2]['u'][w] for m in range(4) for w in range(4, 10)), ('f', (0, 1, 2)))
fa.add_term('H3', 1, ('A', 'c', 'As'), lambda s: sum(s[0]['u'][m] * s[2]['u'][m] for m in range(4)) * s[1]['u'][0], ('f', (0, 1, 2)))
fa.add_term('H4', 1, ('A', 'c', 'As'), lambda s: sum(s[0]['u'][w] * s[2]['u'][w] for w in range(4, 10)) * s[1]['u'][0], ('f', (0, 1, 2)))
fa.add_term('H5', 1, ('c', 'lam', 'lams'), lambda s: s[0]['u'][0] * (sp.Matrix(s[1]['u']).T * V.Mp * sp.Matrix(s[2]['u']))[0], ('f', (0, 1, 2)))
fa.add_term('H6', 1, ('c', 'c', 'cs'), lambda s: s[0]['u'][0] * s[1]['u'][0] * s[2]['u'][0], ('f', (0, 1, 2)))


def dA(x, M, N):
    v = 0
    if M < 4:
        v += I * x['p'][M] * x['u'][N]
    if N < 4:
        v -= I * x['p'][N] * x['u'][M]
    return v


fa.add_term('B1', 1, ('A', 'A'), lambda s: sum(dA(s[0], m, n) * dA(s[1], m, n) for m in range(4) for n in range(m + 1, 4)), ('k', (0, 1)))
fa.add_term('B2', 1, ('A', 'A'), lambda s: sum(dA(s[0], m, w) * dA(s[1], m, w) for m in range(4) for w in range(4, 10)), ('k', (0, 1)))
fa.add_term('dsAdsA', 1, ('A', 'A'), lambda s: (I * dot4(s[0]['p'], s[0]['u'])) * (I * dot4(s[1]['p'], s[1]['u'])), ('k', (0, 1)))
H = {h: fa.local([(1, h)], h, 0) for h in ('H1', 'H2', 'H3', 'H4', 'H5', 'H6')}
F2 = {'D': fa.local([(1, 'dsAdsA')], 'int(d*A,d*A)', 0), 'B1': fa.local([(1, 'B1')], 'int(dA_V,dA_V)', 0),
      'B2': fa.local([(1, 'B2')], 'int(dA_perp,dA_perp)', 0), 'B3': fa.local([(2, 'kinlam')], 'int(lam,dslash lam)', 0),
      'B4': fa.local([(1, 'kingh')], 'int(dc,A*)', 0)}
QH = {h: fa.Q(F) for h, F in H.items()}
BR = {k: fa.bracket(fa.I3, F) for k, F in F2.items()}

R = sp.Rational
ALPHA_R = {'G1': R(11, 4), 'G2': R(27, 8), 'G3': R(11, 4), 'G4': R(0), 'G5': R(0)}
ALPHA_L = {'G1': R(19, 4), 'G2': R(0), 'G3': R(19, 4), 'G4': R(0), 'G5': R(-9)}


def F3(beta, kappa):
    def part(a):
        return {'H1': a['G4'] / 4 - a['G1'] / 2, 'H2': 0, 'H3': -a['G4'], 'H4': -2 * a['G3'], 'H5': a['G5'], 'H6': 0}
    hr, hl = part(ALPHA_R), part(ALPHA_L)
    b1, b2, b3, b4 = beta
    k1, k2 = kappa
    add = {'H1': b1 + k1, 'H2': 2 * b2 + k2, 'H3': b4 - 4 * k1, 'H4': b4 - 2 * k2, 'H5': -2 * b3 + b4,
           'H6': b4 / 2 - 2 * k1}
    return {k: hr[k] + add[k] for k in hr}, hl


def fval(F, legs):
    d = F(legs)
    assert set(d) <= {'f'}, d
    return sp.expand(d.get('f', 0))


def main():
    t0 = time.time()
    rng = random.Random(777)
    params = [((0, 0, 0, 0), (0, 0))] + [
        (tuple(R(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(4)),
         tuple(R(rng.randint(-9, 9), rng.randint(1, 5)) for _ in range(2))) for _ in range(3)]
    out = {'realisation': 'step2/engine_b/fa.py (vertex functions, colour ([c0,c1],c2) stripped); configs seeds 40-45',
           'sectors': {}}
    for sec in S.SECTORS:
        recs = []
        for seed in range(40, 46):
            legs = S.config(sec, seed)
            vG = {g: S.G_value(g, legs) for g in ('G1', 'G2', 'G3', 'G4', 'G5') if g in S.SECTOR_G[sec]}
            vQH = {h: fval(F, legs) for h, F in QH.items()}
            vBR = {k: fval(F, legs) for k, F in BR.items()}
            ok_all = True
            first = None
            for beta, kappa in params:
                hr, hl = F3(beta, kappa)
                f2 = {'D': R(27, 16), 'B1': beta[0], 'B2': beta[1], 'B3': beta[2], 'B4': beta[3]}
                lhs_r = sp.expand(sum(hr[h] * vQH[h] for h in hr) + sum(f2[k] * vBR[k] for k in f2))
                lhs_l = sp.expand(sum(hl[h] * vQH[h] for h in hl))
                rhs_r = sp.expand(sum(ALPHA_R[g] * v for g, v in vG.items()))
                rhs_l = sp.expand(sum(ALPHA_L[g] * v for g, v in vG.items()))
                ok = sp.expand(lhs_r - rhs_r) == 0 and sp.expand(lhs_l - rhs_l) == 0
                ok_all = ok_all and ok
                if first is None:
                    first = {'LHS_rational': str(lhs_r), 'RHS_rational': str(rhs_r), 'LHS_log32': str(lhs_l),
                             'RHS_log32': str(rhs_l), 'holds': ok}
            z_r = sp.expand(R(27, 16) * vBR['D'] + R(11, 8) * (vBR['B1'] + vBR['B2']) - sum(ALPHA_R[g] * v for g, v in vG.items()))
            z_l = sp.expand(R(19, 8) * (vBR['B1'] + vBR['B2']) - R(9, 2) * vBR['B3'] - sum(ALPHA_L[g] * v for g, v in vG.items()))
            recs.append({'F3_zero_representative_holds': z_r == 0 and z_l == 0, 'seed': seed, 'G_values': {g: str(v) for g, v in vG.items()}, 'beta0_kappa0': first,
                         'all_parameter_sets_hold': ok_all})
        out['sectors'][sec] = {'samples': recs, 'all_hold': all(r['all_parameter_sets_hold'] for r in recs),
                               'F3_zero_representative_all_hold': all(r['F3_zero_representative_holds'] for r in recs)}
        print(sec, out['sectors'][sec]['all_hold'], flush=True)
    out['parameter_sets'] = [{'beta': [str(b) for b in be], 'kappa': [str(k) for k in ka]} for be, ka in params]
    out['seconds'] = round(time.time() - t0, 1)
    with open(OUT, 'w') as fh:
        json.dump(out, fh, indent=1)
    print('done', out['seconds'])


if __name__ == '__main__':
    main()
