import mpmath as mp, qb_core as B, qb_job as QJ, qb_numeps as NE
mp.mp.dps=30
xi=B.XI[:8]+[0]*8
job=QJ.DJob('AW','lamP',['AW'],'xi',xi,terms='cubic',tag='res')
vals=NE.rnd_vals(job,7)
seq=job.sequences()[0]
types,curs,Num=NE.config_data(job,seq,0,vals)
ex=NE.exact_coeffs(job,types,curs,Num)
print('exact', ex)
p2=sum(float(vals[f'r0_{m}'])**2 for m in range(4)); print('|p|^2', p2)
for k in range(6,22,2):
    e=mp.mpf(2)**(-k)
    v=NE.numeric_value(job,types,curs,Num,e)
    r=v-ex[1]/e-ex[0]
    print(k, mp.nstr(r,8), mp.nstr(r/e,8), mp.nstr(r/(e*mp.log(e)),8))
