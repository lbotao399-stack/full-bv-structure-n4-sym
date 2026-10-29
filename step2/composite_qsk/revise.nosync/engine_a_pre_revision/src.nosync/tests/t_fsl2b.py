import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cq_core, cq_fsl2 as F, json
for N in [1,2,3]:
    t0=time.time()
    r=F.run(N, n=(1,-1j,0,0))
    print(json.dumps(r)[:600]); print('%.1fs'%(time.time()-t0))
