import time, sys
import sympy as sp
from sympy.polys.domains import QQ_I
import ce_job as CJ
I = sp.I
w1 = {4: 1, 5: 2, 6: -1, 7: 0, 8: 3, 9: 1}
w2 = {4: 2, 5: -1, 6: 1, 7: 1, 8: 0, 9: -2}
for opt in (3, 2):
    t0 = time.time()
    job = CJ.CJob('AW', 'AW', ['AW', 'AW'], pol='sym', za=w1, zb=w2)
    res = job.kernels(mlist=(0,), option=opt)
    R = job.R; g = job.g
    Y0 = [g[f'o0_{c}'] for c in range(4, 10)]; Y1 = [g[f'o1_{c}'] for c in range(4, 10)]
    W1 = [R(QQ_I.from_sympy(sp.nsimplify(w1[c]))) for c in range(4, 10)]; W2 = [R(QQ_I.from_sympy(sp.nsimplify(w2[c]))) for c in range(4, 10)]
    dot = lambda x, y: sum((a * b for a, b in zip(x, y)), R.zero)
    form = dot(W1, Y0) * dot(W2, Y1) + dot(W1, Y1) * dot(W2, Y0)
    ap = sum((g[f'a{m}'] * (g[f'r0_{m}'] + g[f'r1_{m}']) for m in range(4)), R.zero)
    totP = sum((dm.get(0, {}).get('1', R.zero) for dm in res['Pb'].values()), R.zero)
    totD = sum((dm.get(0, {}).get('1', R.zero) for dm in res['D0'].values()), R.zero)
    others = [(w, m, b) for tab in ('Pb', 'D0') for w, dm in res[tab].items() for m, d in dm.items() for b in d if b != '1']
    print('option', opt, 'Pb/form:', sp.simplify(totP.as_expr() / (ap * form * R(QQ_I(0,1))).as_expr()),
          'D0/form:', sp.simplify(totD.as_expr() / form.as_expr()), 'non-rational parts:', others, 'configs', res['configs'], 'time %.1f' % (time.time() - t0))
