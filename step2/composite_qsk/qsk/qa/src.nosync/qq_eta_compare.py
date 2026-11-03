"""QA: intrinsic eta kernels Upsilon^0_eta, engine A (TQ5 rows, generic eta in A's S_+ (x) W^*) versus engine B's code
(read-only: qb_jobinj.IJob with the momentum h injected at the vertex, as in qb_etak.eta_intrinsic) run on A's data
mapped by the intertwiner (eta in Sigma^*: eta_B = SST0 eta_A).  Upsilon^0_eta = sum_mu (-i d/dh_mu) D[S^(1)(rho(e_mu)
eta)] + D[s_0].  c-power: A = C2^{e/2} B(S0 data), e = n_out_lam - n_let_lam - 1."""
import os, sys, json, time, resource
import sympy as sp
from sympy.polys.domains import QQ_I
ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'qa/src.nosync'))
import qq_ab_compare as AB          # noqa: E402
QJ = AB.QJ
import qb_jobinj as IJ              # noqa: E402
import qb_core as BC                # noqa: E402
A = json.load(open(os.path.join(QSK, 'engine_a/results.json')))['tasks']['TQ5']
etaA = sp.Matrix([AB.cplx(x) for x in A['eta']])
etaB = sp.expand(AB.SST0 * etaA)
print('eta_B components:', [i for i in range(16) if etaB[i] != 0])
LOG = AB.LOG
out = []
t0 = time.time()
for r in A['intrinsic_rows']:
    if len(r['ordering']) != 1:
        continue
    ta, tb = r['pair']; d = r['data']
    tBa, za = AB.to_B(ta, d['zeta_a'], 'letter'); tBb, zb = AB.to_B(tb, d['zeta_b'], 'letter')
    to, po = AB.to_B(d['legs'][0][0], d['legs'][0][1], 'out')
    mom = d['legs'][0][2]
    lam = ('lamP', 'lamM')
    e = (1 if d['legs'][0][0] in lam else 0) - sum(t in lam for t in (ta, tb)) - 1
    assert e % 2 == 0
    fac = AB.C2 ** (e // 2)
    sub = {sp.Symbol('r0_%d' % m): int(mom[m]) for m in range(4)}
    tot = {}
    for mu in range(4):
        ximu = list(BC.Rbm[mu] * etaB)
        job = IJ.IJob(tBa, tBb, [to], 'xi', ximu, hmax=1, pol=[po], za=za, zb=zb, tag='qaeta%d' % mu)
        res = job.kernel(mlist=(0,))
        ex = [0, 0, 0, 0]; ex[mu] = 1
        for w, dm in res['total'].items():
            for bb, P in dm.get(0, {}).items():
                c = job.hcoeff(P, tuple(ex)) * QQ_I(0, -1)
                if c:
                    tot[bb] = tot.get(bb, 0) + sp.expand(c.as_expr(*[sp.Symbol(n) for n in job.names]).subs(sub))
    js = QJ.DJob(tBa, tBb, [to], 's0', list(etaB), pol=[po], za=za, zb=zb, tag='qaetas0')
    rs = js.kernel(mlist=(0,))
    for w, dm in rs['total'].items():
        for bb, P in dm.get(0, {}).items():
            tot[bb] = tot.get(bb, 0) + sp.expand(P.as_expr(*[sp.Symbol(n) for n in js.names]).subs(sub))
    bval = sp.expand(fac * sum(v * LOG[k] for k, v in tot.items()))
    aval = AB.Aval(r['value'])
    rec = {'A_pair': r['pair'], 'A_leg': r['leg_set'], 'A_value': str(aval), 'B_value_at_A_data': str(sp.nsimplify(bval)),
           'agree': sp.expand(bval - aval) == 0, 'agree_with_minus': sp.expand(bval + aval) == 0,
           'seconds': round(time.time() - t0, 1), 'rss_MB': round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, 1)}
    print(rec, flush=True)
    out.append(rec)
    json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/eta_compare.json'), 'w'), indent=1)
