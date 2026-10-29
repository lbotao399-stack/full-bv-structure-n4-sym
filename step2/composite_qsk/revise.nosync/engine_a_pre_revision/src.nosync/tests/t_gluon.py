import sys, os, time, json
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cq_core, cq_gluon as GL
t0=time.time()
r=GL.run(int(sys.argv[1]), kinds=tuple(sys.argv[2].split(',')))
print(json.dumps({k:v for k,v in r.items() if k!='ward'}, indent=0)[:2500]); print('%.1fs'%(time.time()-t0))
