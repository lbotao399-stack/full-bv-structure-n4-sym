import time, json
import qa_twist as TW
t0=time.time()
res,T=TW.run(which=('gammagamma','bgamma','betab_bdc','betab_gg'), log=lambda s: print(s, round(time.time()-t0,1), flush=True))
json.dump(res, open('../out.nosync/tq3_rest.json','w'), indent=1, default=str)
print('done', time.time()-t0)
