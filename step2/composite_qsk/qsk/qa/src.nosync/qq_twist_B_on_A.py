"""QA: run engine B's code (read-only import) on engine A's twisted beta-beta data mapped by the intertwiner S.
Also report where A's twisted letters sit in B's basis (to compare with B's own twisted letters)."""
import os, sys, json
import sympy as sp

ROOT = '/Users/libotao/Desktop/N4_Costello_BV_2026-10-27'
QSK = os.path.join(ROOT, 'step2/composite_qsk/qsk')
sys.path.insert(0, os.path.join(QSK, 'qa/src.nosync'))
import qq_ab_compare as AB     # noqa: E402  (imports engine B read-only)
QJ = AB.QJ

d = json.load(open(os.path.join(QSK, 'qa/out.nosync/twistA_data.json')))
xiA = sp.Matrix([AB.cplx(x) for x in d['xi']])
xiB = sp.expand(AB.S * xiA)
print('xi_A in B basis (S0):', [sp.nsimplify(x) for x in xiB])
gvec = [sp.nsimplify(xiB[i]) for i in range(16)]
out = []
for bi in range(3):
    tB, zb = AB.to_B('lamP', d['beta_z'][bi], 'letter')
    print('beta_%d (A) as B letter:' % bi, tB, zb)
for k in range(3):
    print('gamma_out %d' % k, d['gam_out'][k])
for a in range(2):
    print('dc_out %d (A) in B:' % a, AB.to_B('lamM', d['dc_out'][a], 'out'))
tBa, za = AB.to_B('lamP', d['beta_z'][0], 'letter')
tBb, zb = AB.to_B('lamP', d['beta_z'][1], 'letter')
for rec in d['A_betabeta_values'][:6]:
    k, a = rec['k'], rec['a']
    tt0, p0 = AB.to_B('AW', d['gam_out'][k], 'out')
    tt1, p1 = AB.to_B('lamM', d['dc_out'][a], 'out')
    job = QJ.DJob(tBa, tBb, [tt0, tt1], 'xi', gvec, pol=[p0, p1], za=za, zb=zb, tag='qatw')
    res = job.kernel(mlist=(0,))
    subsd = {}
    for j in range(2):
        q = [AB.cplx(x) for x in rec['q'][j]]
        for mu in range(4):
            subsd[sp.Symbol(f'r{j}_{mu}')] = q[mu]
    e = 1 + 1 - 2   # n_out_lam - n_let_lam + xi
    fac = AB.C2 ** (e // 2)
    bw = {}
    for wd, dm in res['total'].items():
        tot = 0
        for bb, P in dm.get(0, {}).items():
            tot += sp.expand(P.as_expr(*[sp.Symbol(n) for n in job.names]).subs(subsd)) * AB.LOG.get(bb, 1)
        bw[str(list(wd))] = str(sp.nsimplify(sp.expand(fac * tot)))
    print('k', k, 'a', a, 'A:', rec['A_value'], 'B on A data:', bw, flush=True)
    out.append({'k': k, 'a': a, 'A': rec['A_value'], 'B_on_A_data': bw})
json.dump(out, open(os.path.join(QSK, 'qa/out.nosync/twist_B_on_A.json'), 'w'), indent=1)
