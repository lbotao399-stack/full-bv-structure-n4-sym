import time, numpy as np
import qa_s1 as S1, qa_tq4 as T4, qa_util as U
t0=time.time()
R=U.Rng(5)
B=S1.chir_basis('lamP')
xi=sum(int(c)*v for c,v in zip([1,2,-1,1,0,1,-2,1],B)).astype(complex)
f=T4.run(xi)
print({k:(U.vstr(v) if isinstance(v,list) else v) for k,v in f.items()}, time.time()-t0)
