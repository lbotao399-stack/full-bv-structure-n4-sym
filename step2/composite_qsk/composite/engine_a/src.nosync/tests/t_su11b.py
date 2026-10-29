import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import cq_core as Q, cq_su11 as S, cq_gval as GV
Q.RULE['koszul']=sys.argv[1]; GV.INTR['on']=(sys.argv[2]=='on')
S.densities(npts=3)
