import mpmath as mp, qb_core as B, qb_job as QJ, qb_numeps as NE, qb_numT as NT
mp.mp.dps=30
xi=B.XI[:8]+[0]*8
job=QJ.DJob('AW','lamP',['AW'],'xi',xi,terms='cubic',tag='rT')
vals=NE.rnd_vals(job,7)
seq=job.sequences()[0]
types,curs,Num=NE.config_data(job,seq,0,vals)
for T in (1,4):
    print('T',T)
    for k in range(8,24,3):
        e=mp.mpf(2)**(-k)
        d=NT.value_T(job,types,curs,Num,e,mp.mpf(T))-NE.numeric_value(job,types,curs,Num,e)
        print('  ',k, mp.nstr(d,10), mp.nstr(d/e,10), mp.nstr(d/(e*mp.log(e)),10))
