import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cq_core, cq_fsl2 as F
for n in [(1,1j,0,0),(1,-1j,0,0),(0,0,1,1j),(0,0,1,-1j)]:
    t0=time.time()
    try:
        r=F.run(1, n=n, kinds=('L',))
    except Exception as e:
        print(n, 'ERR', e)
    print(n, '%.1fs'%(time.time()-t0))
