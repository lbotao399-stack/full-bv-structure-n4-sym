import time, json
import qa_s1 as S1, qa_tq2 as T2, qa_util as U
t0=time.time()
B=S1.chir_basis('lamP')
xi=sum(int(c)*v for c,v in zip([1,2,-1,1,0,1,-2,1],B)).astype(complex)
rows=T2.run(xi,n=0,log=lambda s: print(s, round(time.time()-t0,1), flush=True))
json.dump(rows, open('../out.nosync/tq2_n0_test.json','w'), indent=0)
for r in rows[:12]: print(r['pair'], r['leg_set'], r['ordering'], r['eps0'], r['poles'], r['nonzero_configs'], r['census_chains'])
print('total time', time.time()-t0, len(rows))
