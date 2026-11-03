# quick test of DJob on the G-cubic chain: pair (AW, lamP), outputs [AW], S = T5 at the special edge
import time, sympy as sp
import qb_core as B, qb_job as QJ
xi = B.XI[:8] + [0]*8
t0=time.time()
job = QJ.DJob('AW','lamP',['AW'],'xi',xi, terms='cubic')
res = job.kernel(mlist=(0,1), by_chain=True, verbose=True)
print('configs', res['configs'], 'sec', round(res['seconds'],1))
for seq, sig, cres in res['chains']:
    for w, dm in cres.items():
        for m,d in dm.items():
            for b,P in d.items():
                print(seq, sig, w, m, b, len(P.terms()))
tab = QJ.table_str(job, res['total'])
for kk,v in tab.items(): print(kk, v[:300])
