"""QA: engine A's code (read-only) on A's twisted beta0 beta1 -> gamma2 dc0 at real and at complex output momenta;
compare with engine B's polynomial (twist_B_on_A2.json) evaluated at the same momenta."""
import os, sys, json
import numpy as np
import sympy as sp
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'engine_a/src.nosync'))
import qa_twist as TW, qa_core as QC, qa_util as U
T = TW.Twist(); gx = QC.odd_gen(T.xi)
Bp = json.load(open(os.path.join(QSK, 'qa/out.nosync/twist_B_on_A2.json')))['beta0 beta1 gamma2']['[0, 1]']
def Bval(q0, q1):
    sub = {sp.Symbol('r0_%d' % m): sp.nsimplify(complex(q0[m]).real) + sp.I * sp.nsimplify(complex(q0[m]).imag) for m in range(4)}
    sub.update({sp.Symbol('r1_%d' % m): sp.nsimplify(complex(q1[m]).real) + sp.I * sp.nsimplify(complex(q1[m]).imag) for m in range(4)})
    L = {'1': 1, 'log2': sp.log(2), 'log3': sp.log(3)}
    return sp.nsimplify(sp.expand(sum(sp.sympify(e).subs(sub) * L[b] for b, e in Bp.items())))
rng = np.random.default_rng(1)
out = []
for kind in ('real', 'real', 'complex_null', 'complex_generic'):
    if kind == 'real':
        q = [rng.integers(-3, 4, 4).astype(complex) for _ in range(2)]
    elif kind == 'complex_null':
        q = [sum(rng.integers(-3, 4) * T.nv[w] for w in T.hol + T.exa) for _ in range(2)]
    else:
        q = [rng.integers(-3, 4, 4) + 1j * rng.integers(-3, 4, 4) for _ in range(2)]
    legs = [('AW', T.gam_out[2], q[0]), ('lamM', T.dc_out[0], q[1])]
    v = QC.pair_value(('lamP', T.beta_z[0], None), ('lamP', T.beta_z[1], None), legs, gx, what=('fin',))
    a = U.vstr(v['fin'])
    print(kind, 'q=', [list(np.round(x, 3)) for x in q], 'A:', a, '| B poly (B basis, S0 scale):', Bval(q[0], q[1]), flush=True)
    out.append({'kind': kind, 'A': a, 'B': str(Bval(q[0], q[1]))})
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/twistA_realmom.json'), 'w'), indent=1)
