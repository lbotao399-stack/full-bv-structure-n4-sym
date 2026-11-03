import time, numpy as np
from fractions import Fraction
t0=time.time()
import qa_s1 as S1, qa_core as QC
import cq_eval as CV
b = np.array([1,2,-1,3])
gP = QC.even_gen('P_b', a=b)
w1 = np.eye(6,dtype=int)[0]; w2 = np.eye(6,dtype=int)[1]
q1 = np.array([1,0,2,-1]); q2 = np.array([0,1,-1,1])
for (pa,pb) in [((0,1)),((1,0)),((0,0)),((2,2))]:
  legs = [('AW', np.eye(6,dtype=int)[pa], q1), ('AW', np.eye(6,dtype=int)[pb], q2)]
  Xa=('AW',w1,None); Xb=('AW',w2,None)
  # my engine, line placement only
  def line_only(sig, spos, n):
      adj = {spos, spos+1 if spos < n else 0}
      return sig not in adj
  mine = QC.pair_value(Xa,Xb,legs,gP,only_sigma=line_only)['fin']
  allp, parts = QC.pair_value(Xa,Xb,legs,gP,split=True)
  # Part C
  pc = CV.pair_on_legs(Xa,Xb,legs,('Pv',tuple(b))) if False else None
  import cq_gval as GV
  tl,tv = GV.kernel(('AW',w1,{(0,0,0,0):1}),('AW',w2,{(0,0,0,0):1}),[(l,p,q) for l,p,q in legs],('Pv',tuple(b)))
  pcv = GV.combine((tl,tv),-1)
  print((pa,pb),'mine line:',CV.vstr(mine),'| partC:',CV.vstr(pcv),'| all:',CV.vstr(allp['fin']))
  for k,v in sorted(parts.items(), key=str):
      if not CV.viszero(v): print('   ',k,CV.vstr(v))
print('time',time.time()-t0, QC.STATS)
