"""Calibration of the independent realisation on Step-1 two-point data (manuscript thm c4-quadratic,
thm c4-obstruction): singular parts of the 1PI two-point functions and the (c, A|_V) obstruction."""
import json
import sys
import numpy as np
import ft
import loop

rng = np.random.default_rng(2026)
quad = loop.Quad()
EPS = np.array([4e-3, 2e-3, 1e-3, 5e-4, 2.5e-4, 1.25e-4, 6.25e-5, 3.125e-5])


def fit_sing(vals, eps):
    # vals = a/eps + b log eps + c + d eps log eps + e eps + f eps^2 log eps + g eps^2
    B = np.stack([1 / eps, np.log(eps), np.ones_like(eps), eps * np.log(eps), eps, eps ** 2 * np.log(eps), eps ** 2], 1)
    sol, *_ = np.linalg.lstsq(B, vals, rcond=None)
    return sol


def rand_vec(n):
    return rng.normal(size=n) + 1j * rng.normal(size=n)


p = np.array([0.7, -0.4, 0.9, 0.3])
out = {}

# ---- A|_V, A|_V
def Aleg(p, comps, mono=()):
    f = np.zeros((10, 3), complex)
    for M in comps:
        f[M] = rand_vec(3)
    return ft.Leg('A', f.reshape(30), p, mono)


def run_pair(L1, L2, name, targets):
    ol = loop.OneLoop([L1, L2], rng)
    vals = {}
    for e in EPS:
        for m, v in ol.value(e, quad).items():
            vals.setdefault(m, []).append(v)
    res = {}
    for m, vs in vals.items():
        sol = fit_sing(np.array(vs), EPS)
        res[str(m)] = {'eps^-1': sol[0], 'log': sol[1]}
    tgt = {}
    for key, (T2f, slots) in targets.items():
        tgt[key] = loop.local(T2f, slots, [L1, L2])
    out[name] = {'fit_err': ol.max_fit_err(), 'fit': {k: {kk: [vv.real, vv.imag] for kk, vv in v.items()} for k, v in res.items()},
                 'targets': {k: {str(m): [v.real, v.imag] for m, v in t.items()} for k, t in tgt.items()}}
    return res, tgt


comps_V = range(4)
comps_W = range(4, 10)
L1 = Aleg(p, comps_V)
L2 = Aleg(-p, comps_V)
res, tgt = run_pair(L1, L2, 'AV_AV', {'AA': (ft.T2_AA, ('A', 'A')), 'dAdA': (ft.T2_dAdA, ('A', 'A')), 'dsAdsA': (ft.T2_dsAdsA, ('A', 'A'))})
m = ()
# decompose the log coefficient on dAdA and dsAdsA using a second polarisation pair would be needed; here
# test the claimed total: eps^-1 -> 9/4 unit (A,A), log -> -unit (dA,dA)
print('AV: eps^-1 / (unit*AA) =', res[str(m)]['eps^-1'] / (ft.UNIT * tgt['AA'][m]),
      ' log / (unit*dAdA) =', res[str(m)]['log'] / (ft.UNIT * tgt['dAdA'][m]))
L1 = Aleg(p, comps_W)
L2 = Aleg(-p, comps_W)
res, tgt = run_pair(L1, L2, 'AW_AW', {'AA': (ft.T2_AA, ('A', 'A')), 'dAdA': (ft.T2_dAdA, ('A', 'A'))})
print('AW: eps^-1 / (unit*AA) =', res[str(m)]['eps^-1'] / (ft.UNIT * tgt['AA'][m]),
      ' log / (unit*dAdA) =', res[str(m)]['log'] / (ft.UNIT * tgt['dAdA'][m]))
# mixed V-W: should vanish
L1 = Aleg(p, comps_V)
L2 = Aleg(-p, comps_W)
res, tgt = run_pair(L1, L2, 'AV_AW', {})
print('AV-AW mixing:', {k: (abs(v['eps^-1']), abs(v['log'])) for k, v in res.items()})

# ---- lambda lambda: legs r1 lambda(p), r2 lambda(-p)
L1 = ft.Leg('l', rand_vec(48), p, (1,))
L2 = ft.Leg('l', rand_vec(48), -p, (2,))
res, tgt = run_pair(L1, L2, 'l_l', {'ldl': (ft.T2_ldl, ('l', 'l'))})
m = (1, 2)
print('lam: eps^-1 =', abs(res[str(m)]['eps^-1']), ' log / (unit*ldl) =', res[str(m)]['log'] / (ft.UNIT * tgt['ldl'][m]))

# ---- c A^*|_V : legs r1 A^*(p), r2 c(-p)
f = np.zeros((10, 3), complex)
for M in range(4):
    f[M] = rand_vec(3)
L1 = ft.Leg('As', f.reshape(30), p, (1,))
L2 = ft.Leg('c', rand_vec(3), -p, (2,))
res, tgt = run_pair(L1, L2, 'As_c', {'dcAs': (ft.T2_dcAs, ('c', 'As'))})
print('cA*: eps^-1 =', abs(res[str(m)]['eps^-1']), ' log / (unit*dcAs) =', res[str(m)]['log'] / (ft.UNIT * tgt['dcAs'][m]))

# ---- obstruction (c, A|_V): lim Q [W_1^{1PI} - I^CT_2]  ->  -27/8 unit int(d^*dc, d^*A|_V)
c1 = ft.Leg('c', rand_vec(3), p, (1,))
A2 = Aleg(-p, comps_V)
legsets = [[ft.Q_leg(c1), A2], [c1, ft.Q_leg(A2)]]
ols = [loop.OneLoop(ls, rng) for ls in legsets]
tot = []
for e in EPS:
    acc = {}
    for ls, ol in zip(legsets, ols):
        loop.addto(acc, ol.value(e, quad))
        loop.addto(acc, loop.ICT2(ls, e), -1.0)
    tot.append(acc.get((0, 1), 0))
tot = np.array(tot)
sol = fit_sing(tot, EPS)
G = loop.local(ft.T2_ddc_dsA, ('c', 'A'), [c1, A2])[(1,)]
print('Ob2 (c,A_V): eps^-1 coeff', abs(sol[0]), ' log coeff', abs(sol[1]), ' limit/(unit*int(d*dc,d*A)) =', sol[2] / (ft.UNIT * G))
out['Ob2'] = {'eps^-1': abs(sol[0]), 'log': abs(sol[1]), 'ratio': [(sol[2] / (ft.UNIT * G)).real, (sol[2] / (ft.UNIT * G)).imag],
              'raw': [[float(e), complex(v).real, complex(v).imag] for e, v in zip(EPS, tot)]}
json.dump(out, open('test_2pt.json', 'w'), indent=1, default=lambda x: [complex(x).real, complex(x).imag])
