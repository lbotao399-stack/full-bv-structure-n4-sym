import time, sys
import sympy as sp
from sympy.polys.domains import QQ_I
import ce_job as CJ
I = sp.I
w1 = {4: 1, 5: 2, 6: -1, 7: 0, 8: 3, 9: 1}
w2 = {4: 2, 5: -1, 6: 1, 7: 1, 8: 0, 9: -2}
t0 = time.time()
job = CJ.CJob('AW', 'AW', ['AW', 'AW'], pol='sym', za=w1, zb=w2)
res = job.kernels(mlist=(0, 1), want=('Pb', 'D'), verbose=True)
print('configs', res['configs'], 'sec', round(res['seconds'], 1))
R = job.R; g = job.g
def show(tab):
    for w, dm in tab.items():
        for m, d in dm.items():
            for b, P in d.items():
                print('  word', w, 'eps^-%d' % m, b, ':', P.as_expr() if len(P) < 30 else ('%d terms' % len(P)))
print('Pb:'); show(res['Pb'])
print('D0:'); show(res['D0'])
# target: -(1/3) i a.(p0+p1) [(w1.Y0)(w2.Y1)+(w1.Y1)(w2.Y0)] summed over words (colour pairing) in units u
Y0 = [g[f'o0_{c}'] for c in range(4, 10)]; Y1 = [g[f'o1_{c}'] for c in range(4, 10)]
W1 = [R(QQ_I.from_sympy(sp.nsimplify(w1[c]))) for c in range(4, 10)]; W2 = [R(QQ_I.from_sympy(sp.nsimplify(w2[c]))) for c in range(4, 10)]
dot = lambda x, y: sum((a * b for a, b in zip(x, y)), R.zero)
form = dot(W1, Y0) * dot(W2, Y1) + dot(W1, Y1) * dot(W2, Y0)
ap = sum((g[f'a{m}'] * (g[f'r0_{m}'] + g[f'r1_{m}']) for m in range(4)), R.zero)
tot = R.zero
for w, dm in res['Pb'].items():
    tot += dm.get(0, {}).get('1', R.zero)
logs = {b: sum((dm.get(0, {}).get(b, R.zero) for dm in res['Pb'].values()), R.zero) for b in ('log2', 'log3')}
print('Pb total (rational part), words summed:', (tot).as_expr())
print('ratio to i a.(p0+p1) form:')
tgt = ap * form * R(QQ_I(0, 1))
import sympy
print(sympy.simplify(tot.as_expr() / tgt.as_expr()))
print('log parts', {b: v.as_expr() for b, v in logs.items()})
tot0 = sum((dm.get(0, {}).get('1', R.zero) for dm in res['D0'].values()), R.zero)
print('D0 total rational part / form:', sympy.simplify(tot0.as_expr() / form.as_expr()))
print('D0 logs', {b: sum((dm.get(0, {}).get(b, R.zero) for dm in res['D0'].values()), R.zero).as_expr() for b in ('log2', 'log3')})
print('time', time.time() - t0)
