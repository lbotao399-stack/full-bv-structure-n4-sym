import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cq_su11 as S, cq_gval as GV
print('nulls ok', S.check_nulls())
t0=time.time()
GV.INTR['on'] = (sys.argv[1]=='on') if len(sys.argv)>1 else True
S.densities(npts=3)
print('%.1fs'%(time.time()-t0))
