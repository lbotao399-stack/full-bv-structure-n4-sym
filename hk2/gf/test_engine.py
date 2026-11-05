import time, sympy as sp
import bvjet, models as Mo, freecx as Fx
X = bvjet.X
t0=time.time()
M = Mo.build_model(True)
for frame in ('min','C','Psi'):
    F = Mo.Fields(M, frame)
    S = Mo.S0(F)
    t1=time.time()
    SS = M.bracket(S,S)
    print(frame, 'terms S0', len(S), ' {S0,S0}=0:', M.is_zero_functional(SS), ' %.1fs'%(time.time()-t1))
