import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import cq_core, cq_gates as GT
rng=np.random.default_rng(1)
t0=time.time()
#GT.g_weights(); print('%.1fs'%(time.time()-t0))
GT.g_k1(); print('%.1fs'%(time.time()-t0))
GT.g_cores(rng); print('%.1fs'%(time.time()-t0))
GT.g_bilinear(); print('%.1fs'%(time.time()-t0))
GT.g_position(rng); print('%.1fs'%(time.time()-t0))
for g in GT.GATES:
    if not g['pass']: print(json.dumps(g, default=str)[:1500])
