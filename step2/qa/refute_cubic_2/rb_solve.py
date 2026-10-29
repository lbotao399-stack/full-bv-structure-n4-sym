"""Solve Ob_3 = Q F_3 + {I_3, F_2} with MY maps (rb_plugback_maps.json); check the claimed family, kernel,
F_3 = 0 representative; direct plane-wave plug-back of the family (rational and log(3/2) parts separately)
in all four target sectors and four no-invariant sectors; Wess-Zumino on (c,c,A_V) plane waves.
Also: the extended ansatz F_3 in span{H1..H6, Y_V, Y_W} (Y_V = int(lam,rho([A_V,lam])), Y_W likewise).
"""
import sys, os, json, warnings
warnings.filterwarnings('ignore')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import sympy as sp
from fractions import Fraction as Fr
from rb_func import *

here = os.path.dirname(os.path.abspath(__file__))
maps = json.load(open(os.path.join(here, 'rb_plugback_maps.json')))
Gs = ['G1', 'G2', 'G3', 'G4', 'G5']
sec_of = {'G1': ('cAVAV', 0), 'G2': ('cAVAV', 1), 'G3': ('cAWAW', 0), 'G4': ('ccAs', 0), 'G5': ('cll', 0)}


def vec(entry):
    v = []
    for g in Gs:
        s, i = sec_of[g]
        v.append(sp.Rational(entry[s]['coef'][i]))
    return sp.Matrix(v)


QH = {h: vec(maps['QH'][h]) for h in ['H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'YV', 'YW']}
BR = {b: vec(maps['brackets_I3'][b]) for b in ['D1', 'B1', 'B2', 'B3', 'B4', 'AA', 'AH']}
out = {}
out['QH_matrix_rows_G1_G5'] = [[str(QH[h][i]) for h in ['H1', 'H2', 'H3', 'H4', 'H5', 'H6']] for i in range(5)]
out['Q_YV'] = [str(x) for x in QH['YV']]
out['Q_YW'] = [str(x) for x in QH['YW']]
out['brackets'] = {b: [str(x) for x in BR[b]] for b in BR}

a1, a2, a3, a4, a5 = sp.symbols('alpha1:6')
b1, b2, b3, b4 = sp.symbols('beta1:5')
h = sp.symbols('h1:7')
Ob3 = sp.Matrix([a1, a2, a3, a4, a5])
QF3 = sum((h[j] * QH['H%d' % (j + 1)] for j in range(6)), sp.zeros(5, 1))
F2br = sp.Rational(27, 16) * BR['D1'] + b1 * BR['B1'] + b2 * BR['B2'] + b3 * BR['B3'] + b4 * BR['B4']
eqs = list(QF3 + F2br - Ob3)
out['row_equations'] = [str(sp.expand(e)) for e in eqs]
# solvability: G2 row
out['G2_row'] = str(sp.expand(eqs[1]))
sol = sp.solve([e for k, e in enumerate(eqs) if k != 1], [h[0], h[2], h[3], h[4]], dict=True)
out['solution_free_h2_h6'] = {str(k): str(v) for k, v in sol[0].items()}
# claimed family
k1, k2 = sp.symbols('kappa1 kappa2')
claimed = {h[0]: -a1 / 2 + a4 / 4 + b1 + k1, h[1]: 2 * b2 + k2, h[2]: -a4 + b4 - 4 * k1,
           h[3]: -2 * a3 + b4 - 2 * k2, h[4]: a5 - 2 * b3 + b4, h[5]: b4 / 2 - 2 * k1}
res = [sp.simplify(e.subs(claimed).subs(a2, sp.Rational(27, 8))) for e in eqs]
out['claimed_family_residuals'] = [str(r) for r in res]
# kernel of Q on span H
Qm = sp.Matrix.hstack(*[QH['H%d' % (j + 1)] for j in range(6)])
out['rank_Q_spanH'] = int(Qm.rank())
out['kernel_Q_spanH'] = [str(list(v.T)) for v in Qm.nullspace()]
k_claim = [sp.Matrix([1, 0, -4, 0, 0, -2]), sp.Matrix([0, 1, 0, -2, 0, 0])]
out['claimed_kernel_in_kernel'] = [str(list((Qm * v).T)) for v in k_claim]
# image of (F3, beta) map
Mfull = sp.Matrix.hstack(Qm, BR['B1'], BR['B2'], BR['B3'], BR['B4'])
out['rank_(QH,beta_brackets)'] = int(Mfull.rank())
out['G2_in_image'] = bool(sp.Matrix.hstack(Mfull, sp.Matrix([0, 1, 0, 0, 0])).rank() == Mfull.rank())
# extended ansatz with Y_V, Y_W
Qx = sp.Matrix.hstack(Qm, QH['YV'], QH['YW'])
out['rank_Q_extended'] = int(Qx.rank())
out['kernel_Q_extended_basis_H1..H6,YV,YW'] = [str(list(v.T)) for v in Qx.nullspace()]
# F3 = 0 representative
L = sp.log(sp.Rational(3, 2))
alpha_claim = {a1: sp.Rational(11, 4) + sp.Rational(19, 4) * L, a2: sp.Rational(27, 8),
               a3: sp.Rational(11, 4) + sp.Rational(19, 4) * L, a4: 0, a5: -9 * L}
rep = {b1: alpha_claim[a1] / 2, b2: alpha_claim[a1] / 2, b3: alpha_claim[a5] / 2, b4: 0, k1: 0, k2: -alpha_claim[a3]}
hs = [sp.simplify(claimed[h[j]].subs(alpha_claim).subs(rep)) for j in range(6)]
out['F3_zero_rep_h'] = [str(x) for x in hs]
res0 = [sp.simplify(e.subs({h[j]: hs[j] for j in range(6)}).subs(alpha_claim).subs(rep)) for e in eqs]
out['F3_zero_rep_residuals'] = [str(r) for r in res0]
part = {k: sp.simplify(v.subs(alpha_claim).subs({b1: 0, b2: 0, b3: 0, b4: 0, k1: 0, k2: 0})) for k, v in claimed.items()}
out['particular_beta0'] = {str(k): str(v) for k, v in part.items()}
print(json.dumps(out, indent=1), flush=True)

# ---------------------------------------------------------------- direct plane-wave plug-back
rng = np.random.default_rng(4242)


def rint(shape, lo=-3, hi=4):
    return rng.integers(lo, hi, size=shape).astype(float)


def rand_mom(n):
    ps = [rint(4) for _ in range(n - 1)]
    ps.append(-sum(ps))
    return ps


def leg(slot, p, part=None):
    v = rint(SHAPE[slot])
    if slot in ('A', 'Ast') and part == 'V':
        v[4:] = 0
    if slot in ('A', 'Ast') and part == 'W':
        v[:4] = 0
    return (slot, v, p)


def sector_legs(name):
    p = rand_mom(3)
    return {'cAVAV': lambda: [leg('c', p[0]), leg('A', p[1], 'V'), leg('A', p[2], 'V')],
            'cAWAW': lambda: [leg('c', p[0]), leg('A', p[1], 'W'), leg('A', p[2], 'W')],
            'cll': lambda: [leg('c', p[0]), leg('lam', p[1]), leg('lam', p[2])],
            'ccAs': lambda: [leg('c', p[0]), leg('c', p[1]), leg('Ast', p[2], 'V')],
            'ccAV': lambda: [leg('c', p[0]), leg('c', p[1]), leg('A', p[2], 'V')],
            'cAVAW': lambda: [leg('c', p[0]), leg('A', p[1], 'V'), leg('A', p[2], 'W')],
            'cAl': lambda: [leg('c', p[0]), leg('A', p[1]), leg('lam', p[2])],
            'cclst': lambda: [leg('c', p[0]), leg('c', p[1]), leg('lamst', p[2])],
            'ccAsW': lambda: [leg('c', p[0]), leg('c', p[1]), leg('Ast', p[2], 'W')]}[name]()


HF = [F_H1, F_H2, F_H3, F_H4, F_H5, F_H6]
GF = [F_G1, F_G2, F_G3, F_G4, F_G5]
alpha_rat = [sp.Rational(11, 4), sp.Rational(27, 8), sp.Rational(11, 4), 0, 0]
alpha_log = [sp.Rational(19, 4), 0, sp.Rational(19, 4), 0, -9]
param_sets = [(0, 0, 0, 0, 0, 0), (sp.Rational(3, 2), -1, sp.Rational(1, 4), 2, -3, sp.Rational(5, 2)),
              (-2, sp.Rational(7, 3), 1, sp.Rational(-1, 2), sp.Rational(1, 3), 4)]
plug = []
for part_name, alph, d1coef in [('rational', alpha_rat, sp.Rational(27, 16)), ('log(3/2)', alpha_log, 0)]:
    for ps in param_sets:
        bb = dict(zip([b1, b2, b3, b4, k1, k2], ps))
        sub = dict(zip([a1, a2, a3, a4, a5], alph))
        hv = [float(sp.nsimplify(claimed[h[j]].subs(sub).subs(bb))) for j in range(6)]
        # for the log part the G2 row needs alpha2 = 0 and no D1 (solvable: alpha2_log = 0)
        Fb = [float(d1coef), float(bb[b1]), float(bb[b2]), float(bb[b3]), float(bb[b4])]
        worst = 0.0; scale = 0.0
        for sec in ['cAVAV', 'cAWAW', 'cll', 'ccAs', 'cAVAW', 'cAl', 'cclst', 'ccAsW']:
            for _ in range(3):
                Lg = sector_legs(sec)
                lhs = sum(float(alph[i]) * value(GF[i], Lg) for i in range(5))
                rhs = sum(hv[j] * Qvalue(HF[j], Lg) for j in range(6))
                for c, F in zip(Fb, [F_D1, F_B1, F_B2, F_B3, F_B4]):
                    if c != 0:
                        rhs += c * bracket_value(F_I3, 0, 3, F, 2, Lg)
                worst = max(worst, abs(lhs - rhs)); scale = max(scale, abs(lhs))
        plug.append({'part': part_name, 'params(beta1..4,kappa1,kappa2)': [str(x) for x in ps],
                     'h': [str(sp.nsimplify(claimed[h[j]].subs(sub).subs(bb))) for j in range(6)],
                     'max|Ob3 - QF3 - {I3,F2}|': worst, 'max|Ob3|': scale})
        print(plug[-1], flush=True)
out['plane_wave_plugback'] = plug

# Wess-Zumino on (c,c,A_V): Q Ob_3 + {I_3, Ob_2}, Ob_2 = -27/8 int(d*dc, d*A_V)
wz = []
for alph, tag in [(alpha_rat, 'rational part (alpha2 = 27/8)'), (alpha_log, 'log part (alpha2 = 0, Ob_2 absent)')]:
    worst = 0.0; scale = 0.0
    for _ in range(4):
        Lg = sector_legs('ccAV')
        q = sum(float(alph[i]) * Qvalue(GF[i], Lg) for i in range(5))
        b = bracket_value(F_I3, 0, 3, F_E2, 2, Lg) * (-27.0 / 8.0 if tag.startswith('rational') else 0.0)
        worst = max(worst, abs(q + b)); scale = max(scale, abs(q))
    wz.append({'case': tag, 'max|Q Ob3 + {I3,Ob2}|': worst, 'max|Q Ob3|': scale})
    print(wz[-1], flush=True)
# with a wrong alpha2 the identity fails (sensitivity)
Lg = sector_legs('ccAV')
q = Qvalue(F_G2, Lg) * (27 / 8 + 1) - 27 / 8 * bracket_value(F_I3, 0, 3, F_E2, 2, Lg)
wz.append({'case': 'alpha2 = 27/8 + 1 (control)', '|Q Ob3 + {I3,Ob2}|': abs(q)})
out['wess_zumino'] = wz
json.dump(out, open(os.path.join(here, 'rb_solve.json'), 'w'), indent=1, default=str)
